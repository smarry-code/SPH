// test_e2e_bc1.cpp
// End-to-end тест на кейсе BC1 (удар по карбиду бора).
// ТЗ: contact/04-e2e-test-brief.md. Автоматизирует шаг 6 плана
// (contact/01-plan.md) — интеграционную проверку всего конвейера на реальной
// осесимметричной геометрии с последующей сверкой с эталонной кривой.
//
// ВАЖНО ПРО РАЗМЕРНОСТЬ (решение E1, ТЗ §2): кейс требует dim = 2, поэтому
// цель собирается с -DSPH_DIM=2 и компилирует СВОИ копии src/*.cpp. Линковка
// с объектниками цели SPH (dim = 1) невозможна по ODR: layout ParticleSystem
// зависит от dim. См. docs/04-explanation/decisions/adr-002-dim-constexpr.md.
//
// ЧТО ПРОВЕРЯЕТСЯ:
//   1. Геометрия строится и даёт ожидаемое число частиц (ТЗ §1.2, §4.4).
//   2. Конвейер проходит до t_total, шаг по времени не коллапсирует (§4.4).
//   3. Ни одного NaN/Inf в x, vx, rho, p на каждом сохранённом кадре; rho > 0,
//      mass > 0 на всех кадрах (§4.4) — кадры читаются из sph_output_*.txt,
//      которые пишет сам конвейер.
//   4. Кривая vx[0][moni] сверяется с эталоном (§4.3): интерполяция на общую
//      сетку, метрики rel_l2 и max_abs_dev.
//   5. Известное ограничение D2=A: контакт не видит пару B4C<->B4C (§6).
//   6. Отсутствие побочных эффектов в рабочем каталоге (§7).
//
// РЕЖИМЫ (§5):
//   без аргументов  — smoke (delta_x = 8e-4, t_total = 3e-7), цель < 30 с;
//   --full          — полный кейс (delta_x = 2e-4, t_total = 3e-6) + контактные
//                     сценарии §6;
//   --emit-reference— напечатать эталон текущего режима в stdout (см. шапку
//                     contact/e2e/bc1_reference.txt);
//   --verbose       — не глушить вывод конвейера.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"
#include "ext_force.hpp"

static_assert(dim == 2, "test_e2e_bc1 требует сборки с -DSPH_DIM=2 (решение E1, ТЗ §2)");
static_assert(axis_sim, "test_e2e_bc1 требует axis_sim = true (осесимметричная постановка BC1)");

#ifndef SPH_TEST_DATA_DIR
// Задаётся из CMake (target_compile_definitions). Фолбэк нужен только для
// сборки вне CMake и заведомо неверен — тест сообщит об отсутствии эталона.
#define SPH_TEST_DATA_DIR "contact/e2e"
#endif

namespace fs = std::filesystem;

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    if (cond)
    {
        std::printf("  ok   : %s\n", what.c_str());
    }
    else
    {
        std::printf("  FAIL : %s\n", what.c_str());
        ++failures;
    }
}

void note(const std::string& what)
{
    std::printf("  info : %s\n", what.c_str());
}

// ===========================================================================
// 1. Геометрия кейса BC1
// ===========================================================================
// Перенесена из закомментированного блока src/input.cpp:490-596 (внутри
// shock_tube). Здесь — ЕДИНСТВЕННЫЙ исполняемый экземпляр этой постановки;
// src/input.cpp намеренно не изменяется (ТЗ §1.4). Если постановка станет
// частью продового кода — переносить отсюда, а блок в input.cpp удалить.
//
// Исправления дефектов оригинала (ТЗ §1.3):
//   - добавлен двухпроходный подсчёт + sys.resize() ДО записи (в оригинале
//     запись в sys.x[0][ntotal] шла в контейнер нулевого размера — UB);
//   - накопление `r += delta_x` / `z += delta_x` в условии цикла заменено на
//     целочисленные счётчики: число частиц больше не зависит от накопленной
//     ошибки округления double на границе диапазона;
//   - `moni_particle = 12973` (глобал, которого больше нет) заменён на
//     геометрическое правило, воспроизводящее тот же индекс при delta_x = 2e-4
//     и осмысленное при любом другом шаге (см. pick_moni_particle);
//   - sys.axis заполняется после resize (std::vector<char>, bool-семантика).

struct LayerSpec
{
    double r_hi;      // верхняя граница r, включительно
    double z_lo;      // нижняя граница z, включительно
    double z_hi;      // верхняя граница z, включительно
    int    itype;
    double rho;
    double vz;        // vx[1], м/с
    const char* name;
};

struct Bc1Geometry
{
    size_t ntotal = 0;
    size_t moni_particle = 0;
    double delta_x = 0.0;
    std::array<size_t, 4> layer_begin{};
    std::array<size_t, 4> layer_count{};
    std::array<size_t, 4> layer_nr{};
    std::array<size_t, 4> layer_nz{};
};

// Число узлов сетки lo, lo+dx, ... <= hi.
// Относительный допуск 1e-9 нужен там, где (hi-lo)/dx — целое в точной
// арифметике (слой 4 по z даёт ровно 126.0): без него узел на границе то
// попадает, то не попадает в диапазон в зависимости от округления.
size_t grid_count(double lo, double hi, double dx)
{
    if (hi < lo) return 0;
    const double n = std::floor((hi - lo) / dx + 1.0 + 1e-9);
    return n <= 0.0 ? 0u : static_cast<size_t>(n);
}

std::array<LayerSpec, 4> bc1_layers(double dx)
{
    // Опорные плоскости (ТЗ §1.1) — как в оригинальном блоке input.cpp.
    const double Z1 = 0.006;
    const double Z2 = Z1 + 0.003917;   // 0.009917
    const double Z3 = Z2 + 0.009033;   // 0.018950
    const double Z4 = Z3 + 0.0254;     // 0.044350
    return std::array<LayerSpec, 4>{{
        // r_hi,     z_lo,          z_hi,            itype, rho,  vz,     name
        {0.04375,  0.0,           0.006,             8,   626.0, 2210.0, "ударник (полиуретан 626)"},
        {0.0346,   Z1 + dx,       Z2,               10,  2510.0, 2210.0, "керамика летящая (B4C)"},
        {0.0346,   Z2 + 4.0 * dx, Z3 + 4.0 * dx,    10,  2510.0,    0.0, "керамика покоящаяся (B4C)"},
        {0.0254,   Z3 + 5.0 * dx, Z4 + 4.0 * dx,     5,  2638.0,    0.0, "окно (LiF)"},
    }};
}

// Точка съёма эталонной кривой.
// В оригинале был литерал moni_particle = 12973. Разбор индекса при
// delta_x = 2e-4: слой 3 начинается с 10076, n_z = 46, значит
// 12973 = 10076 + 62*46 + 45 — это последний узел по z (тыльная поверхность
// покоящейся керамики, граница с окном) на радиусе r = 62.5*dx = 0.0125 м.
// Правило ниже воспроизводит именно эту точку и остаётся осмысленным при
// любом delta_x (иначе в smoke-режиме индекс 12973 просто вне диапазона).
constexpr double kMoniRadius = 0.0125;

size_t pick_moni_particle(const Bc1Geometry& g)
{
    const size_t n_r = g.layer_nr[2];
    const size_t n_z = g.layer_nz[2];
    size_t best_ir = 0;
    double best_err = 1e300;
    for (size_t ir = 0; ir < n_r; ++ir)
    {
        const double r = (0.5 + static_cast<double>(ir)) * g.delta_x;
        const double err = std::fabs(r - kMoniRadius);
        if (err < best_err) { best_err = err; best_ir = ir; }
    }
    return g.layer_begin[2] + best_ir * n_z + (n_z - 1);
}

Bc1Geometry build_bc1_geometry(ParticleSystem& sys, double delta_x)
{
    const auto layers = bc1_layers(delta_x);

    Bc1Geometry g;
    g.delta_x = delta_x;

    // --- проход 1: только подсчёт (даёт размер ДО любой записи) -------------
    size_t total = 0;
    for (size_t L = 0; L < layers.size(); ++L)
    {
        const size_t n_r = grid_count(0.5 * delta_x, layers[L].r_hi, delta_x);
        const size_t n_z = grid_count(layers[L].z_lo, layers[L].z_hi, delta_x);
        g.layer_nr[L] = n_r;
        g.layer_nz[L] = n_z;
        g.layer_begin[L] = total;
        g.layer_count[L] = n_r * n_z;
        total += n_r * n_z;
    }
    g.ntotal = total;

    // --- resize ДО заполнения (исправление дефекта §1.3) --------------------
    sys.resize(total);

    // --- проход 2: заполнение ----------------------------------------------
    size_t k = 0;
    for (size_t L = 0; L < layers.size(); ++L)
    {
        const LayerSpec& s = layers[L];
        for (size_t ir = 0; ir < g.layer_nr[L]; ++ir)
        {
            const double r = (0.5 + static_cast<double>(ir)) * delta_x;
            for (size_t iz = 0; iz < g.layer_nz[L]; ++iz)
            {
                const double z = s.z_lo + static_cast<double>(iz) * delta_x;
                sys.axis[k]  = (std::fabs(r) <= 0.5 * delta_x) ? 1 : 0;
                sys.hsml[k]  = 2.0 * delta_x;
                sys.u[k]     = 0.0;
                sys.p[k]     = 0.0;
                sys.x[0][k]  = r;
                sys.x[1][k]  = z;
                sys.itype[k] = s.itype;
                sys.rho[k]   = s.rho;
                sys.rho_0[k] = s.rho;
                sys.vx[0][k] = 0.0;
                sys.vx[1][k] = s.vz;
                // Плоская формула массы — из раскомментированной строки
                // оригинала. Два осесимметричных варианта там закомментированы
                // и намеренно не переносятся (менять постановку — не задача теста).
                sys.mass[k]  = s.rho * delta_x * delta_x;
                ++k;
            }
        }
    }

    g.moni_particle = pick_moni_particle(g);
    return g;
}

// ===========================================================================
// 2. Кривая и компаратор (ТЗ §4.3)
// ===========================================================================

struct Curve
{
    std::vector<double> t;
    std::vector<double> v;
    std::vector<int>    itype;
    bool empty() const { return t.empty(); }
};

// Чтение файла формата log.txt / эталона: три поля «время vx[0] itype»,
// строки с '#' и пустые — комментарии. Возможен UTF-8 BOM в начале (его пишет
// main.cpp в реальный log.txt).
Curve read_curve(const fs::path& p, std::string& err)
{
    Curve c;
    std::ifstream f(p);
    if (!f)
    {
        err = "не удалось открыть " + p.string();
        return c;
    }
    std::string line;
    while (std::getline(f, line))
    {
        size_t b = 0;
        if (line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF &&
            static_cast<unsigned char>(line[1]) == 0xBB && static_cast<unsigned char>(line[2]) == 0xBF)
        {
            b = 3;
        }
        while (b < line.size() && (line[b] == ' ' || line[b] == '\t' || line[b] == '\r')) ++b;
        if (b >= line.size() || line[b] == '#') continue;
        std::istringstream is(line.substr(b));
        double t = 0.0, v = 0.0;
        int it = 0;
        if (!(is >> t >> v >> it)) continue;
        c.t.push_back(t);
        c.v.push_back(v);
        c.itype.push_back(it);
    }
    return c;
}

// Линейная интерполяция кривой в точке tq (кривая упорядочена по t).
double interp(const Curve& c, double tq)
{
    const size_t n = c.t.size();
    if (n == 0) return 0.0;
    if (tq <= c.t.front()) return c.v.front();
    if (tq >= c.t.back())  return c.v.back();
    const size_t i = static_cast<size_t>(
        std::lower_bound(c.t.begin(), c.t.end(), tq) - c.t.begin());
    if (i == 0) return c.v.front();
    const double t0 = c.t[i - 1], t1 = c.t[i];
    const double den = t1 - t0;
    if (den <= 0.0) return c.v[i];
    const double a = (tq - t0) / den;
    return c.v[i - 1] + a * (c.v[i] - c.v[i - 1]);
}

struct CompareResult
{
    bool   comparable = false;
    double coverage = 0.0;
    double rel_l2 = 0.0;
    double max_abs_dev = 0.0;
    size_t worst_node = 0;
    double worst_t = 0.0;
    double worst_calc = 0.0;
    double worst_ref = 0.0;
    double ref_amplitude = 0.0;   // max|ref| на общей сетке
    std::string reason;
};

constexpr size_t kNodes      = 200;    // §4.3: 200 узлов общей сетки
constexpr double kMinCover   = 0.8;    // §4.3: покрытие не меньше 80 %
constexpr double kTolRelL2   = 1e-3;   // §4.3, регрессионный режим
constexpr double kTolMaxDev  = 1.0;    // м/с

CompareResult compare_curves(const Curve& calc, const Curve& ref)
{
    CompareResult r;
    if (calc.t.size() < 2 || ref.t.size() < 2)
    {
        r.reason = "одна из кривых содержит меньше двух точек";
        return r;
    }
    const double lo = std::max(calc.t.front(), ref.t.front());
    const double hi = std::min(calc.t.back(),  ref.t.back());
    const double calc_span = calc.t.back() - calc.t.front();
    if (hi <= lo || calc_span <= 0.0)
    {
        r.reason = "диапазоны времени расчёта и эталона не пересекаются";
        return r;
    }
    r.coverage = (hi - lo) / calc_span;
    if (r.coverage < kMinCover)
    {
        std::ostringstream os;
        os << "покрытие эталоном диапазона расчёта " << (r.coverage * 100.0)
           << " % < " << (kMinCover * 100.0) << " % — постановки не совпадают";
        r.reason = os.str();
        return r;
    }

    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < kNodes; ++i)
    {
        const double a = static_cast<double>(i) / static_cast<double>(kNodes - 1);
        const double t = lo + a * (hi - lo);
        const double vc = interp(calc, t);
        const double vr = interp(ref, t);
        const double d  = std::fabs(vc - vr);
        r.ref_amplitude = std::max(r.ref_amplitude, std::fabs(vr));
        if (d > r.max_abs_dev)
        {
            r.max_abs_dev = d;
            r.worst_node  = i;
            r.worst_t     = t;
            r.worst_calc  = vc;
            r.worst_ref   = vr;
        }
        num += (vc - vr) * (vc - vr);
        den += vr * vr;
    }
    // Норма эталона может быть нулевой (кривая начинается с нулей): тогда
    // относительная метрика не определена, и решает только max_abs_dev.
    r.rel_l2 = (den > 0.0) ? std::sqrt(num / den) : (num > 0.0 ? 1e300 : 0.0);
    r.comparable = true;
    return r;
}

// ===========================================================================
// 3. Изоляция побочных эффектов (ТЗ §7)
// ===========================================================================
// time_integration() пишет log.txt в ТЕКУЩИЙ каталог жёстко зашитым именем,
// output() — sph_output_*.txt и sph_gp.csv туда же. Прогон в каталоге
// разработчика затёр бы лог реального расчёта. Поэтому каждый прогон уходит
// в собственный временный каталог, а после — проверка, что рабочий каталог
// не изменился (see check_no_side_effects).

class ScopedWorkDir
{
public:
    explicit ScopedWorkDir(const std::string& tag)
        : prev_(fs::current_path())
    {
        std::ostringstream os;
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        os << "sph_e2e_bc1_" << tag << "_" << static_cast<long long>(stamp)
           << "_" << counter_++;
        dir_ = fs::temp_directory_path() / os.str();
        fs::remove_all(dir_);
        fs::create_directories(dir_);
        fs::current_path(dir_);
    }
    ~ScopedWorkDir()
    {
        std::error_code ec;
        fs::current_path(prev_, ec);
        if (!keep_) fs::remove_all(dir_, ec);
    }
    ScopedWorkDir(const ScopedWorkDir&) = delete;
    ScopedWorkDir& operator=(const ScopedWorkDir&) = delete;

    const fs::path& dir() const { return dir_; }
    void keep() { keep_ = true; }

private:
    fs::path prev_;
    fs::path dir_;
    bool keep_ = false;
    static int counter_;
};
int ScopedWorkDir::counter_ = 0;

// Глушение вывода конвейера: single_step/output/time_integration печатают
// в std::wcout/std::cout десятки строк на кадр. Без этого (а) отчёт теста
// нечитаем, (б) --emit-reference не может печатать эталон в stdout.
class ScopedSilence
{
public:
    explicit ScopedSilence(bool active) : active_(active)
    {
        if (!active_) return;
        old_w_ = std::wcout.rdbuf(sink_w_.rdbuf());
        old_n_ = std::cout.rdbuf(sink_n_.rdbuf());
    }
    ~ScopedSilence()
    {
        if (!active_) return;
        std::wcout.rdbuf(old_w_);
        std::cout.rdbuf(old_n_);
    }
    ScopedSilence(const ScopedSilence&) = delete;
    ScopedSilence& operator=(const ScopedSilence&) = delete;

private:
    bool active_;
    std::wostringstream sink_w_;
    std::ostringstream  sink_n_;
    std::wstreambuf* old_w_ = nullptr;
    std::streambuf*  old_n_ = nullptr;
};

// ===========================================================================
// 4. Проверка кадров (ТЗ §4.4)
// ===========================================================================
// Кадры пишет сам конвейер (output.cpp) в sph_output_<ts>.txt. Колонки:
//   0:N 1:itype 2:mass 3..3+dim-1:x 3+dim..3+2dim-1:vx 3+2dim:rho 3+2dim+1:p ...
struct FrameCheck
{
    size_t frames = 0;
    size_t rows = 0;
    size_t bad_finite = 0;
    size_t bad_rho = 0;
    size_t bad_mass = 0;
    std::string first_problem;
};

FrameCheck check_frames(const fs::path& dir)
{
    FrameCheck fc;
    std::vector<fs::path> files;
    for (const auto& e : fs::directory_iterator(dir))
    {
        const std::string n = e.path().filename().string();
        if (n.rfind("sph_output_", 0) == 0 && e.path().extension() == ".txt")
        {
            files.push_back(e.path());
        }
    }
    std::sort(files.begin(), files.end());
    fc.frames = files.size();

    constexpr size_t kIdxMass = 2;
    const size_t kIdxRho = 3 + 2 * dim;

    std::vector<double> vals;
    for (const auto& f : files)
    {
        std::ifstream in(f);
        std::string line;
        while (std::getline(in, line))
        {
            if (line.empty()) continue;
            // Строка заголовка начинается с 'N'; данные — с цифры.
            const char c0 = line[0];
            if (!(c0 >= '0' && c0 <= '9')) continue;
            vals.clear();
            std::istringstream is(line);
            double v;
            while (is >> v) vals.push_back(v);
            if (vals.size() < kIdxRho + 2) continue;
            ++fc.rows;
            for (double x : vals)
            {
                if (!std::isfinite(x))
                {
                    ++fc.bad_finite;
                    if (fc.first_problem.empty())
                    {
                        fc.first_problem = "не-конечное значение в " + f.filename().string();
                    }
                    break;
                }
            }
            if (!(vals[kIdxMass] > 0.0))
            {
                ++fc.bad_mass;
                if (fc.first_problem.empty())
                {
                    fc.first_problem = "mass <= 0 в " + f.filename().string();
                }
            }
            if (!(vals[kIdxRho] > 0.0))
            {
                ++fc.bad_rho;
                if (fc.first_problem.empty())
                {
                    fc.first_problem = "rho <= 0 в " + f.filename().string();
                }
            }
        }
    }
    return fc;
}

// ===========================================================================
// 5. Прогон одного случая
// ===========================================================================

struct RunSettings
{
    const char* name = "smoke";
    double delta_x = 8e-4;
    double t_total = 3e-7;
    int    contact_force_type = 0;
    bool   silence = true;
    bool   check_frames_flag = true;
    bool   keep_workdir = false;   // --keep-workdir: не удалять каталог прогона
    // Бюджет времени на один прогон, с. Сторож (run_guard.hpp) обрывает расчёт
    // и говорит, чего именно не хватило — времени или продвижения.
    double wall_limit_s = 600.0;   // 10 минут
};

struct RunResult
{
    bool ok = false;
    std::string error;
    Bc1Geometry geom;
    Curve curve;
    FrameCheck frames;
    double t_reached = 0.0;
    double dt_last = 0.0;
    size_t steps = 0;
    // Исход сторожа: None — расчёт дошёл до t_total сам.
    RunGuard::Reason guard_reason = RunGuard::Reason::None;
    std::string guard_message;
    std::string workdir;           // каталог прогона (существует только при --keep-workdir)
};

// Конфигурация расчёта (ТЗ §3). Источник значений:
//   «Керамика/Карбид бора/BC1/Initial_Conditions.txt».
// Файл в кодировке cp1251 и НЕ парсится намеренно: разбор зависел бы от
// рабочего каталога и от кодировки — обе вещи ломают тест непредсказуемо.
void apply_bc1_config(SimulationConfig& cfg, const RunSettings& rs)
{
    cfg.apply_defaults();
    cfg.pa_sph  = 2;
    cfg.nnps    = 2;
    cfg.sle     = 1;
    cfg.skf     = 1;
    cfg.CFL     = 0.1;
    cfg.epsilon = 0.25;
    cfg.alpha   = 1.0;
    cfg.beta    = 1.0;
    cfg.etq     = 0.01;
    cfg.save_t  = 2.5e-8;
    cfg.print_step = 50;
    cfg.t_total = rs.t_total;
    cfg.contact_force_type = rs.contact_force_type;
}

RunResult run_case(const RunSettings& rs)
{
    RunResult res;
    ScopedWorkDir wd(rs.name);
    if (rs.keep_workdir) wd.keep();
    res.workdir = wd.dir().string();

    SimulationConfig cfg;
    apply_bc1_config(cfg, rs);

    ParticleSystem sys;
    res.geom = build_bc1_geometry(sys, rs.delta_x);
    const size_t ntotal = res.geom.ntotal;
    cfg.moni_particle = res.geom.moni_particle;

    // Запас под виртуальные частицы — как в main.cpp.
    cfg.maxn = static_cast<size_t>(static_cast<double>(ntotal) * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    if (sys.size() < cfg.maxn)
    {
        const size_t prev = sys.size();
        sys.resize(cfg.maxn);
        for (size_t i = prev; i < sys.size(); ++i) sys.inzializ_index(i);
    }

    // Начальный шаг. main.cpp ставит dt = 0.005 из ветки shocktube — это
    // одномерная ударная труба, к BC1 отношения не имеющая; для BC1 такой шаг
    // бессмысленен (в contact_force_type = 3 он вошёл бы в импульс первого
    // шага). Берём тот же критерий Куранта, которым конвейер пересчитывает
    // шаг дальше. Решение сверх E1-E3, зафиксировано здесь.
    double dt = compute_cfl_timestep(cfg, sys, ntotal);
    double time = 0.0;

    // Сторож: бюджет времени + детектор застоя. Потолок «шагов до финиша»
    // 2e6 подобран так, чтобы не задевать здоровый расчёт (BC1 при dx = 2e-4
    // финиширует за ~1.4e3 шагов), но ловить коллапс: при dt ~ 1e-12 остаток
    // оценивается в ~2.6e6 шагов. Терпение 200 шагов — чтобы одиночный
    // просадочный шаг не считался застоем.
    RunGuard guard(rs.wall_limit_s, 2.0e6, 200);

    try
    {
        ScopedSilence silence(rs.silence);
        time_integration(cfg, sys, ntotal, dt, time, rs.delta_x, &guard);
    }
    catch (const std::exception& e)
    {
        // Даже при обрыве log.txt уже содержит все выполненные шаги —
        // без этого «прогон не состоялся» не говорит, где именно он умер.
        std::string e_err;
        const Curve partial = read_curve(wd.dir() / "log.txt", e_err);
        std::ostringstream os;
        os << "исключение в конвейере: " << e.what();
        if (!partial.empty())
        {
            os << " [оборвано на шаге " << partial.t.size()
               << ", t = " << partial.t.back() << " с из " << rs.t_total
               << " с, то есть " << (100.0 * partial.t.back() / rs.t_total) << " % расчёта]";
        }
        else
        {
            os << " [до первой записи в log.txt]";
        }
        res.error = os.str();
        res.curve = partial;
        res.steps = partial.t.size();
        res.t_reached = partial.empty() ? 0.0 : partial.t.back();
        return res;
    }
    res.t_reached = time;
    res.dt_last = dt;
    res.guard_reason = guard.reason();
    if (guard.triggered())
    {
        res.guard_message = guard.message();
    }

    std::string err;
    res.curve = read_curve(wd.dir() / "log.txt", err);
    if (!err.empty())
    {
        res.error = err;
        return res;
    }
    // В log.txt ровно одна строка на шаг интегрирования.
    res.steps = res.curve.t.size();

    if (rs.check_frames_flag)
    {
        res.frames = check_frames(wd.dir());
    }
    res.ok = true;
    return res;
}

// ===========================================================================
// 6. Проверки
// ===========================================================================

void report_geometry(const Bc1Geometry& g, bool full)
{
    const auto layers = bc1_layers(g.delta_x);
    for (size_t L = 0; L < 4; ++L)
    {
        std::ostringstream os;
        os << "слой " << (L + 1) << " (" << layers[L].name << "): "
           << g.layer_nr[L] << " x " << g.layer_nz[L] << " = " << g.layer_count[L]
           << ", [" << g.layer_begin[L] << ", " << (g.layer_begin[L] + g.layer_count[L]) << ")";
        note(os.str());
    }
    {
        std::ostringstream os;
        os << "ntotal = " << g.ntotal << ", moni_particle = " << g.moni_particle;
        note(os.str());
    }

    if (full)
    {
        // §1.2: ожидаемое число частиц при delta_x = 2e-4.
        check(g.layer_count[0] == 6789,  "§1.2 слой 1 = 6789 частиц");
        check(g.layer_count[1] == 3287,  "§1.2 слой 2 = 3287 частиц");
        check(g.layer_count[2] == 7958,  "§1.2 слой 3 = 7958 частиц");
        check(g.layer_count[3] == 16129, "§1.2 слой 4 = 16129 частиц");
        check(g.ntotal == 34163,         "§1.2 ntotal = 34163");
        // Целочисленный счётчик (§1.3) обязан дать тот же moni, что литерал
        // 12973 в оригинальном блоке.
        check(g.moni_particle == 12973,  "§1.3 moni_particle = 12973 (как в оригинале)");
    }
    else
    {
        // Регрессионный замок smoke-геометрии: значения получены этим же
        // правилом при delta_x = 8e-4 и зафиксированы, чтобы изменение
        // геометрии не прошло незамеченным.
        check(g.layer_count[0] == 440, "smoke слой 1 = 440 частиц");
        check(g.layer_count[1] == 172, "smoke слой 2 = 172 частицы");
        check(g.layer_count[2] == 516, "smoke слой 3 = 516 частиц");
        check(g.layer_count[3] == 992, "smoke слой 4 = 992 частиц");
        check(g.ntotal == 2120,        "smoke ntotal = 2120");
        check(g.moni_particle == 803,  "smoke moni_particle = 803");
    }

    // moni обязан попасть в слой 3 (покоящаяся керамика) — §1.2.
    const bool in_layer3 = g.moni_particle >= g.layer_begin[2] &&
                           g.moni_particle <  g.layer_begin[2] + g.layer_count[2];
    check(in_layer3, "§1.2 moni_particle лежит в слое 3 (покоящаяся керамика)");
}

void check_run_health(const RunResult& r, const RunSettings& rs)
{
    {
        std::ostringstream os;
        os << "шагов = " << r.steps << ", достигнутое время = " << r.t_reached
           << " (t_total = " << rs.t_total << "), последний dt = " << r.dt_last;
        note(os.str());
    }
    // Сторож времени/застоя. Проверяется ПЕРВЫМ: если он сработал, все
    // последующие проверки говорят об оборванном прогоне, и без этой строки
    // причина обрыва из отчёта не видна.
    if (r.guard_reason == RunGuard::Reason::WallClock)
    {
        check(false, "сторож: расчёт не уложился в бюджет времени (шаг здоров — "
                     "не хватило именно времени)\n" + r.guard_message);
    }
    else if (r.guard_reason == RunGuard::Reason::Stalled)
    {
        check(false, "сторож: расчёт перестал продвигаться (коллапс шага — "
                     "увеличение бюджета не поможет)\n" + r.guard_message);
    }
    else
    {
        note("сторож не срабатывал: расчёт завершился сам, в пределах бюджета "
             + std::to_string(static_cast<long>(rs.wall_limit_s)) + " с");
    }

    // §4.4: расчёт дошёл до t_total, а не оборвался раньше.
    check(r.t_reached >= rs.t_total, "§4.4 расчёт дошёл до t_total");
    check(r.steps > 0, "§4.4 выполнен хотя бы один шаг");

    // §4.4: шаг не сколлапсировал.
    // Опорным берётся ПЕРВЫЙ принятый шаг (t[0] в log.txt), а не последний:
    // при коллапсе последний dt сам по себе крошечный, оценка t_total/dt_last
    // раздувается, и потолок перестаёт что-либо ограничивать.
    // Две проверки: (а) число шагов не больше 10x номинальной оценки,
    // (б) минимальный dt за прогон не упал более чем в 100 раз относительно
    //     первого — это и есть прямой признак коллапса.
    if (r.curve.t.size() >= 2)
    {
        const double dt_first = r.curve.t.front();
        double dt_min = dt_first;
        size_t dt_min_step = 0;
        for (size_t i = 1; i < r.curve.t.size(); ++i)
        {
            const double d = r.curve.t[i] - r.curve.t[i - 1];
            if (d < dt_min) { dt_min = d; dt_min_step = i; }
        }
        {
            std::ostringstream os;
            os << "первый dt = " << dt_first << ", минимальный dt = " << dt_min
               << " (шаг " << dt_min_step << "), отношение = "
               << (dt_min > 0.0 ? dt_first / dt_min : 0.0);
            note(os.str());
        }
        if (dt_first > 0.0)
        {
            const double expected = rs.t_total / dt_first;
            const double cap = 10.0 * expected;
            std::ostringstream os;
            os << "§4.4 число шагов " << r.steps << " < потолка "
               << static_cast<size_t>(cap) << " (10x от оценки " << expected
               << " по первому dt)";
            check(static_cast<double>(r.steps) < cap, os.str());

            constexpr double kMaxDtDrop = 100.0;
            std::ostringstream os2;
            os2 << "§4.4 шаг не сколлапсировал: dt_first/dt_min = "
                << (dt_min > 0.0 ? dt_first / dt_min : 1e300)
                << " < " << kMaxDtDrop;
            check(dt_min > 0.0 && dt_first / dt_min < kMaxDtDrop, os2.str());
        }
        else
        {
            check(false, "§4.4 первый dt должен быть положительным");
        }
    }
    else
    {
        check(false, "§4.4 в log.txt меньше двух шагов — оценить коллапс шага нечем");
    }

    // §4.4: кадры без NaN/Inf, rho > 0, mass > 0.
    {
        std::ostringstream os;
        os << "кадров = " << r.frames.frames << ", строк проверено = " << r.frames.rows;
        note(os.str());
    }
    check(r.frames.frames > 0, "§4.4 конвейер сохранил хотя бы один кадр");
    check(r.frames.rows > 0, "§4.4 кадры содержат данные частиц");
    check(r.frames.bad_finite == 0,
          "§4.4 нет NaN/Inf в кадрах" +
              (r.frames.first_problem.empty() ? std::string() : " (" + r.frames.first_problem + ")"));
    check(r.frames.bad_rho == 0, "§4.4 rho > 0 на всех кадрах");
    check(r.frames.bad_mass == 0, "§4.4 mass > 0 на всех кадрах");

    // Амплитуда снятой кривой. Если она нулевая — волна до монитора не дошла,
    // и сверка с эталоном (§4.3) вырождается в сравнение нулей. Для smoke это
    // штатно: зазор между слоями 2 и 3 равен 4*delta_x, время его закрытия
    // 4*delta_x/2210 (при delta_x = 8e-4 это 1.4e-6 с, то есть больше
    // smoke-t_total = 3e-7 с).
    if (!r.curve.empty())
    {
        double amp = 0.0;
        for (double v : r.curve.v) amp = std::max(amp, std::fabs(v));
        std::ostringstream os;
        os << "max|vx[0][moni]| по расчёту = " << amp << " м/с"
           << "; оценка времени смыкания слоёв 2/3 = "
           << (4.0 * rs.delta_x / 2210.0) << " с";
        note(os.str());
    }

    // Кривая снята с той частицы, на которую рассчитан эталон.
    if (!r.curve.empty())
    {
        const bool same_itype = std::all_of(r.curve.itype.begin(), r.curve.itype.end(),
                                            [](int t) { return t == 10; });
        check(same_itype, "itype частицы-монитора всё время равен 10 (B4C)");
    }
}

// --- §6: известное ограничение D2=A (B4C <-> B4C контакту не виден) --------
// Тест-документация, а НЕ проверка исправления. Когда будет реализовано
// D2=B (идентификатор тела в ParticleSystem), ожидания ниже инвертируются:
// пара «одинаковый материал, разные тела» обязана стать контактной.
void test_known_limitation_b4c_b4c(const Bc1Geometry& g)
{
    std::printf("\n[§6] Известное ограничение D2=A: контакт не различает два тела одного материала\n");

    // (а) доказательство на самой геометрии: слои 2 и 3 целиком itype = 10,
    //     поэтому критерий itype[i] != itype[j] не может сработать ни на одной
    //     паре между ними — независимо от их положений на любом шаге.
    const auto layers = bc1_layers(g.delta_x);
    check(layers[1].itype == 10 && layers[2].itype == 10,
          "§6 слои 2 и 3 — оба itype = 10, критерий D2=A между ними неприменим");

    // (б) операционное подтверждение на ext_force: две перекрывающиеся
    //     частицы B4C дают ноль контактных пар, те же частицы с разными
    //     материалами — одну.
    for (int type = 1; type <= 3; ++type)
    {
        SimulationConfig cfg;
        cfg.apply_defaults();
        cfg.contact_force_type = type;
        cfg.CFL = 0.1;

        ParticleSystem sys;
        sys.resize(2);
        const double h = 1e-3;
        for (size_t i = 0; i < 2; ++i)
        {
            sys.hsml[i] = h;
            sys.mass[i] = 1e-3;
            sys.rho[i] = 2510.0;
            sys.c[i] = 1e4;
            sys.x[0][i] = 0.0;
            sys.x[1][i] = 0.0;
            sys.vx[0][i] = 0.0;
            sys.vx[1][i] = 0.0;
        }
        // Перекрытие: расстояние 0.5*h при среднем радиусе h -> p = 0.5*h > 0.
        sys.x[1][1] = 0.5 * h;
        sys.vx[1][0] = 100.0;   // сближение

        const std::vector<size_t> pi{0}, pj{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        WorkBuffers wb;

        // Пара B4C <-> B4C: оба itype = 10.
        sys.itype[0] = 10; sys.itype[1] = 10;
        ext_force(cfg, sys, 1e-9, 2, 1, pi, pj, dvxdt, wb);
        const size_t pairs_same = wb.contact_pair_count;

        // Та же геометрия, но разные материалы (ударник 8 <-> керамика 10).
        sys.itype[0] = 8; sys.itype[1] = 10;
        ext_force(cfg, sys, 1e-9, 2, 1, pi, pj, dvxdt, wb);
        const size_t pairs_diff = wb.contact_pair_count;

        std::ostringstream a, b;
        a << "§6 contact_force_type=" << type << ": пар B4C<->B4C = 0 (ОГРАНИЧЕНИЕ D2=A)";
        b << "§6 contact_force_type=" << type << ": пар 8<->10 = 1 (контакт работает между материалами)";
        check(pairs_same == 0, a.str());
        check(pairs_diff == 1, b.str());
    }
    note("§6 переход на идентификатор тела (D2=B) — отдельная задача, в объём ТЗ не входит");
}

// --- §7: отсутствие побочных эффектов в рабочем каталоге -------------------
// Снимок берётся ТОЛЬКО по именам, которые пишет конвейер (log.txt,
// sph_output_*.txt, sph_gp.csv, force_timing.csv). Сравнивать весь каталог
// нельзя: при `ctest -j` соседние тесты создают в том же каталоге свои файлы
// (test_particle_system пишет __ps_test_out.txt), и проверка стала бы флакающей.
bool is_pipeline_artifact(const std::string& n)
{
    return n == "log.txt" || n == "sph_gp.csv" || n == "force_timing.csv" ||
           (n.rfind("sph_output_", 0) == 0);
}

struct DirSnapshot
{
    // имя -> размер; отсутствующий файл в отображении просто отсутствует
    std::vector<std::pair<std::string, std::uintmax_t>> artifacts;
};

DirSnapshot snapshot_dir(const fs::path& d)
{
    DirSnapshot s;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(d, ec))
    {
        const std::string n = e.path().filename().string();
        if (!is_pipeline_artifact(n)) continue;
        std::error_code ec2;
        s.artifacts.emplace_back(n, fs::file_size(e.path(), ec2));
    }
    std::sort(s.artifacts.begin(), s.artifacts.end());
    return s;
}

void check_no_side_effects(const fs::path& cwd, const DirSnapshot& before)
{
    std::printf("\n[§7] Побочные эффекты в рабочем каталоге\n");
    const DirSnapshot after = snapshot_dir(cwd);
    {
        std::ostringstream os;
        os << "артефактов конвейера в рабочем каталоге: было " << before.artifacts.size()
           << ", стало " << after.artifacts.size();
        note(os.str());
    }
    check(before.artifacts == after.artifacts,
          "§7 log.txt / sph_output_*.txt / sph_gp.csv в рабочем каталоге не созданы и не изменены");
    if (before.artifacts != after.artifacts)
    {
        for (const auto& a : after.artifacts)
        {
            const auto it = std::find_if(before.artifacts.begin(), before.artifacts.end(),
                                         [&](const auto& b) { return b == a; });
            if (it == before.artifacts.end())
            {
                std::printf("         появился/изменился: %s (%llu байт)\n",
                            a.first.c_str(), static_cast<unsigned long long>(a.second));
            }
        }
    }
}

// ===========================================================================
// 7. Эталон
// ===========================================================================

fs::path reference_path(bool full)
{
    return fs::path(SPH_TEST_DATA_DIR) / (full ? "bc1_reference.txt" : "bc1_reference_smoke.txt");
}

void emit_reference(const RunResult& r, const RunSettings& rs, bool full)
{
    std::printf("# Эталон для end-to-end теста BC1 (удар по карбиду бора)\n");
    std::printf("#\n");
    std::printf("# СГЕНЕРИРОВАНО: test_e2e_bc1 %s--emit-reference\n", full ? "--full " : "");
    std::printf("# Формат: <время, с>  <vx[0][moni_particle], м/с>  <itype[moni_particle]>\n");
    std::printf("# (ровно формат time_integration.cpp:304-305 и BC1/log_*.txt)\n");
    std::printf("#\n");
    std::printf("# ПАРАМЕТРЫ ПОСТАНОВКИ\n");
    std::printf("#   режим:                 %s\n", full ? "full" : "smoke");
    std::printf("#   delta_x:               %g\n", rs.delta_x);
    std::printf("#   ntotal:                %zu\n", r.geom.ntotal);
    std::printf("#   moni_particle:         %zu\n", r.geom.moni_particle);
    std::printf("#   t_total:               %g\n", rs.t_total);
    std::printf("#   CFL:                   0.1\n");
    std::printf("#   nnps / pa_sph / skf:   2 / 2 / 1\n");
    std::printf("#   contact_force_type:    %d\n", rs.contact_force_type);
    std::printf("#   dim / axis_sim:        %zu / %s\n", dim, axis_sim ? "true" : "false");
    std::printf("#   шагов:                 %zu\n", r.steps);
    std::printf("#   commit:                [git rev-parse --short HEAD — вписать вручную]\n");
    std::printf("#   дата:                  [ГГГГ-ММ-ДД — вписать вручную]\n");
    std::printf("#   сгенерировал:          [кто — вписать вручную]\n");
    std::printf("#\n");
    std::printf("# ДОПУСКИ СРАВНЕНИЯ (регрессионный режим, решение E2=A)\n");
    std::printf("#   rel_l2       < %g\n", kTolRelL2);
    std::printf("#   max_abs_dev  < %g м/с\n", kTolMaxDev);
    std::printf("#\n");
    std::printf("# ДАННЫЕ (время  vx[0]  itype)\n");
    for (size_t i = 0; i < r.curve.t.size(); ++i)
    {
        std::printf("%.17g %.17g %d\n", r.curve.t[i], r.curve.v[i], r.curve.itype[i]);
    }
}

void check_against_reference(const RunResult& r, bool full)
{
    const fs::path ref_path = reference_path(full);
    std::printf("\n[§4] Сверка с эталоном: %s\n", ref_path.string().c_str());

    if (!fs::exists(ref_path))
    {
        std::printf("  FAIL : эталонный файл отсутствует: %s\n", ref_path.string().c_str());
        ++failures;
        return;
    }
    std::string err;
    const Curve ref = read_curve(ref_path, err);
    if (!err.empty())
    {
        std::printf("  FAIL : %s\n", err.c_str());
        ++failures;
        return;
    }
    if (ref.t.size() < 2)
    {
        std::printf(
            "  FAIL : эталон не заполнен — в %s нет строк данных (только шапка).\n"
            "         Это ОЖИДАЕМОЕ состояние до того, как владелец признает расчёт\n"
            "         корректным (ТЗ §8 п.5): эталон фиксирует текущее поведение, и\n"
            "         если оно неверно, тест закрепит ошибку.\n"
            "         Как заполнить:\n"
            "           1) убедиться, что расчёт физически корректен (проверяет владелец);\n"
            "           2) ./test_e2e_bc1 %s--emit-reference > %s\n"
            "           3) вписать в шапку commit/дату/автора и просмотреть кривую глазами.\n",
            ref_path.string().c_str(), full ? "--full " : "", ref_path.string().c_str());
        ++failures;
        return;
    }

    const CompareResult c = compare_curves(r.curve, ref);
    if (!c.comparable)
    {
        std::printf("  FAIL : сравнение невозможно: %s\n", c.reason.c_str());
        std::printf("         расчёт: %zu точек, t в [%g, %g]\n",
                    r.curve.t.size(), r.curve.t.front(), r.curve.t.back());
        std::printf("         эталон: %zu точек, t в [%g, %g]\n",
                    ref.t.size(), ref.t.front(), ref.t.back());
        ++failures;
        return;
    }

    std::printf("  info : покрытие = %.1f %%, узлов сетки = %zu\n", c.coverage * 100.0, kNodes);
    if (c.ref_amplitude == 0.0)
    {
        // Вырожденный эталон: на всём диапазоне vx[0][moni] тождественно нулевой.
        // Сравнение формально проходит, но НИЧЕГО не проверяет. Для smoke-режима
        // это штатная ситуация: при t_total = 3e-7 волна до монитора не доходит
        // (ТЗ §5 требует от smoke только «конвейер не падает, нет NaN»).
        std::printf("  ВНИМАНИЕ: эталон тождественно нулевой на всём диапазоне —\n"
                    "         сверка кривой вырождена и регрессию по величине не ловит.\n"
                    "         Смысл имеют только проверки §4.4 и §1.2 выше.\n");
    }
    std::printf("  info : rel_l2 = %.6g (допуск < %g), max_abs_dev = %.6g м/с (допуск < %g)\n",
                c.rel_l2, kTolRelL2, c.max_abs_dev, kTolMaxDev);
    const bool ok_l2 = c.rel_l2 < kTolRelL2;
    const bool ok_dev = c.max_abs_dev < kTolMaxDev;
    if (!ok_l2 || !ok_dev)
    {
        // §7: без этой распечатки диагностировать расхождение невозможно.
        std::printf("  FAIL : расхождение с эталоном\n");
        std::printf("         rel_l2      = %.9g  (допуск < %g) %s\n",
                    c.rel_l2, kTolRelL2, ok_l2 ? "ok" : "ПРЕВЫШЕН");
        std::printf("         max_abs_dev = %.9g  (допуск < %g) %s\n",
                    c.max_abs_dev, kTolMaxDev, ok_dev ? "ok" : "ПРЕВЫШЕН");
        std::printf("         узел максимального расхождения: #%zu, t = %.9g с\n",
                    c.worst_node, c.worst_t);
        std::printf("         расчёт = %.9g м/с, эталон = %.9g м/с, разность = %.9g м/с\n",
                    c.worst_calc, c.worst_ref, c.worst_calc - c.worst_ref);
        ++failures;
    }
    else
    {
        check(true, "§4.3 кривая совпадает с эталоном в пределах допусков");
    }
}

// --- §6, таблица сценариев: устойчивость при contact_force_type = 1..3 -----
void test_contact_scenarios(const RunSettings& base)
{
    for (int type = 1; type <= 3; ++type)
    {
        std::printf("\n[§6] Сценарий contact_force_type = %d (сверяется только устойчивость)\n", type);
        RunSettings rs = base;
        rs.name = "contact";
        rs.contact_force_type = type;
        const RunResult r = run_case(rs);
        if (!r.ok)
        {
            std::printf("  FAIL : прогон не состоялся: %s\n", r.error.c_str());
            ++failures;
            continue;
        }
        check_run_health(r, rs);
        // Эталон для контактных сценариев не требуется (ТЗ §6).
    }
}

void print_usage()
{
    std::printf(
        "Использование: test_e2e_bc1 [--full] [--emit-reference] [--verbose]\n"
        "  (без аргументов)  smoke-режим: delta_x = 8e-4, t_total = 3e-7\n"
        "  --full            полный кейс: delta_x = 2e-4, t_total = 3e-6 + сценарии §6\n"
        "  --emit-reference  напечатать эталон текущего режима в stdout\n"
        "  --verbose         не глушить вывод конвейера\n"
        "  --keep-workdir    не удалять временный каталог прогона (log.txt, кадры)\n"
        "  --time-limit=СЕК  бюджет времени на прогон (по умолчанию 600 = 10 мин;\n"
        "                    0 — без ограничения). При исчерпании расчёт обрывается\n"
        "                    с разбором: не хватило времени или коллапс шага\n");
}

} // namespace

int main(int argc, char** argv)
{
    bool full = false;
    bool emit = false;
    bool verbose = false;
    bool keep_workdir = false;
    double wall_limit_s = 600.0;   // 10 минут на прогон
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--full") full = true;
        else if (a == "--emit-reference") emit = true;
        else if (a == "--verbose") verbose = true;
        else if (a == "--keep-workdir") keep_workdir = true;
        else if (a.rfind("--time-limit=", 0) == 0)
        {
            try
            {
                wall_limit_s = std::stod(a.substr(std::strlen("--time-limit=")));
            }
            catch (const std::exception&)
            {
                std::fprintf(stderr, "Неверное значение --time-limit: %s\n", a.c_str());
                return 2;
            }
            if (wall_limit_s < 0.0) wall_limit_s = 0.0;   // 0 — без ограничения
        }
        else if (a == "-h" || a == "--help") { print_usage(); return 0; }
        else
        {
            std::fprintf(stderr, "Неизвестный аргумент: %s\n", a.c_str());
            print_usage();
            return 2;
        }
    }

    const fs::path cwd = fs::current_path();
    const DirSnapshot before = snapshot_dir(cwd);

    RunSettings rs;
    rs.name = full ? "full" : "smoke";
    rs.delta_x = full ? 2e-4 : 8e-4;
    rs.t_total = full ? 3e-6 : 3e-7;
    rs.silence = !verbose;
    rs.keep_workdir = keep_workdir;
    rs.wall_limit_s = wall_limit_s;

    if (emit)
    {
        // Эталон печатается в stdout — диагностика уходит в stderr, иначе
        // `... --emit-reference > bc1_reference.txt` соберёт мусор.
        std::fprintf(stderr, "test_e2e_bc1: генерация эталона, режим %s...\n", rs.name);
        const RunResult r = run_case(rs);
        if (!r.ok)
        {
            std::fprintf(stderr, "test_e2e_bc1: прогон не состоялся: %s\n", r.error.c_str());
            return 1;
        }
        emit_reference(r, rs, full);
        std::fprintf(stderr, "test_e2e_bc1: напечатано %zu строк данных\n", r.curve.t.size());
        return 0;
    }

    std::printf("=== e2e BC1, режим %s (delta_x = %g, t_total = %g) ===\n",
                rs.name, rs.delta_x, rs.t_total);

    std::printf("\n[§1] Геометрия\n");
    const RunResult r = run_case(rs);
    if (keep_workdir)
    {
        std::printf("  info : каталог прогона сохранён: %s\n", r.workdir.c_str());
    }
    if (!r.ok)
    {
        std::printf("  FAIL : прогон не состоялся: %s\n", r.error.c_str());
        if (!keep_workdir)
        {
            std::printf("         повторить с --keep-workdir, чтобы сохранить log.txt и кадры прогона\n");
        }
        // Частичная кривая всё равно есть — покажем, вёл ли себя шаг по времени
        // нормально до обрыва (коллапс шага и обрыв по плотности — разные болезни).
        if (r.curve.t.size() >= 2)
        {
            const double dt_first = r.curve.t.front();
            double dt_min = dt_first;
            for (size_t i = 1; i < r.curve.t.size(); ++i)
            {
                dt_min = std::min(dt_min, r.curve.t[i] - r.curve.t[i - 1]);
            }
            std::printf("         до обрыва: первый dt = %g, минимальный dt = %g, отношение = %g\n",
                        dt_first, dt_min, dt_min > 0.0 ? dt_first / dt_min : 0.0);
        }
        ++failures;
    }
    else
    {
        report_geometry(r.geom, full);
        std::printf("\n[§4.4] Здоровье расчёта\n");
        check_run_health(r, rs);
        check_against_reference(r, full);
    }

    if (r.ok)
    {
        test_known_limitation_b4c_b4c(r.geom);
    }

    if (full && r.ok)
    {
        test_contact_scenarios(rs);
    }

    check_no_side_effects(cwd, before);

    std::printf("\n=== e2e BC1: %s (провалов: %d) ===\n",
                failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", failures);
    return failures == 0 ? 0 : 1;
}
