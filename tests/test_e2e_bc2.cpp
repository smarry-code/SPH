// test_e2e_bc2.cpp
// End-to-end тест на кейсе BC2 (удар по карбиду бора).
// ТЗ: notes/contact/04-e2e-test-brief.md. Автоматизирует шаг 6 плана
// (notes/contact/01-plan.md) — интеграционную проверку всего конвейера на реальной
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
//   4. Кривая vx[1][moni] (осевая) сверяется с эталоном по ФОРМЕ (§4.3):
//      признаки (пик, фронт 10/50/90 %, спад до 70 %), совмещение по t50
//      переднего фронта, безразмерные отношения расчёт/эталон и грубые ворота.
//      Сверка в общем АБСОЛЮТНОМ времени неприменима: начало отсчёта времени
//      эталона не привязано к моменту удара (см. §2 файла и
//      notes/contact/08-phase-shift.md). Сам компаратор проверяется
//      test_comparator_selfcheck() в обоих режимах.
//   5. Известное ограничение D2=A: контакт не видит пару B4C<->B4C (§6).
//   6. Отсутствие побочных эффектов в рабочем каталоге (§7).
//
// ОДИН ПРОГОН — ОДИН ТИП КОНТАКТА (--contact-type=0..3).
// Раньше один запуск делал базовый прогон и следом три сценария §6, и обрыв
// любого из них останавливал разбор остальных, а в отчёте ctest всё это было
// одним тестом: по «провален» нельзя было сказать, что именно сломалось.
// Теперь тип выбирается аргументом, а перебор делает НАБОР тестов — по одному
// на тип. Каждый падает и читается независимо.
//   contact_force_type = 0 — базовая постановка, только она сверяется с
//     приборным эталоном (§4.3);
//   1..3 — контактные сценарии §6, с них ТЗ требует лишь устойчивости.
//
// СВОДНЫЙ ГРАФИК. Четыре типа — четыре процесса, поэтому в одном из них все
// кривые не собрать. Каждый прогон выгружает свою кривую (--emit-curve), а
// --plot-from рисует их вместе, уже без расчёта. Готовая цель: bc2_plot.
//
// РЕЖИМЫ (§5):
//   без аргументов  — smoke (delta_x = 8e-4, t_total = 3e-7), цель < 30 с;
//   --full          — полный кейс (delta_x = 2e-4, t_total из конфига);
//   --emit-reference— напечатать эталон текущего режима в stdout (см. шапку
//                     tests/data/bc2_reference.txt);
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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "input.hpp"
#include "time_integration.hpp"
#include "ext_force.hpp"
#include "contact_force.hpp"   // contact::Off и имена типов контактной силы
#include "int_force.hpp"   // urs() — диспетчер УРС для init_eos_state

#include "bc2_geometry.hpp"

static_assert(dim == 2, "test_e2e_bc2 требует сборки с -DSPH_DIM=2 (решение E1, ТЗ §2)");
static_assert(axis_sim, "test_e2e_bc2 требует axis_sim = true (осесимметричная постановка BC2)");

#ifndef SPH_CASES_DIR
// Задаётся из CMake (функция sph_case). Фолбэк -- запуск из корня репозитория.
#define SPH_CASES_DIR "cases"
#endif

namespace fs = std::filesystem;

namespace
{

int failures = 0;

// Возвращает cond — чтобы вызывающий мог доложить подробности провала
// сразу после строки FAIL (см. test_comparator_selfcheck).
bool check(bool cond, const std::string& what)
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
    return cond;
}

void note(const std::string& what)
{
    std::printf("  info : %s\n", what.c_str());
}

// ===========================================================================
// 1. Геометрия кейса BC2
// ===========================================================================
// Построение вынесено в tests/bc2_geometry.hpp: его использует и регрессионный
// тест (tests/test_regression.cpp). Здесь — только привязка имён.
using bc2::LayerSpec;
using bc2::Bc2Geometry;
using bc2::bc2_layers;
using bc2::build_bc2_geometry;

// Начальные p и c через штатный диспетчер УРС.
// Без этого sys.c остаётся нулевым после inzializ_index, и самый первый
// compute_cfl_timestep считает шаг как CFL*h/|v| — без скорости звука, то есть
// примерно впятеро завышенным. Конвейер заполняет c только внутри int_force,
// то есть уже ПОСЛЕ того, как этим шагом сделан первый полный шаг.
// Диспетчер берётся общий с конвейером (int_force.hpp), чтобы соответствие
// itype -> УРС не разъезжалось с продовым кодом.
void init_eos_state(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal)
{
    for (size_t i = 0; i < ntotal; ++i)
    {
        urs(cfg.mats, sys, i);
    }
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

// Чтение файла формата log.txt / эталона: три поля «время vx[1] itype»,
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

// ---------------------------------------------------------------------------
// ПОЧЕМУ СРАВНИВАЕТСЯ ФОРМА, А НЕ ЗНАЧЕНИЯ В ОБЩЕМ АБСОЛЮТНОМ ВРЕМЕНИ
// ---------------------------------------------------------------------------
// Эталон — оцифровка графика из статьи-источника (Grady & Moody, SAND96-0551,
// опыт BC2 = оригинальный тест CE6). Статья прямо оговаривает, что начало
// отсчёта времени на этом графике НЕ привязано к моменту удара:
//
//   «the transit times of the first arrival was not measured in these tests,
//    and the profiles are not necessarily aligned in time with expected
//    elastic precursor velocities»
//
// Приборная запись начинается в 1,51 мкс от произвольного триггера. Поэтому
// прежний компаратор (интерполяция обеих кривых на 200 узлов ОБЩЕГО
// АБСОЛЮТНОГО времени с метриками rel_l2 / max_abs_dev) измерял бессмысленную
// величину: любой сдвиг начала отсчёта эталона трактовался как ошибка расчёта.
// Наблюдавшийся «фазовый сдвиг» ~1,07 мкс — артефакт той методики сверки,
// а не дефект расчёта (разбор — notes/contact/08-phase-shift.md,
// notes/contact/09-article-conformance.md).
//
// Новая схема: из каждой кривой независимо извлекаются признаки формы (пик,
// уровни 10/50/90 % на переднем фронте, спад до 70 % пика), кривые совмещаются
// по t50 переднего фронта, и сравниваются БЕЗРАЗМЕРНЫЕ отношения признаков.
// Величина совмещения печатается как справочная и в вердикт не входит.

// Уровни в долях пика, по которым берутся признаки формы.
constexpr double kLevelLo   = 0.10;   // низ переднего фронта
constexpr double kLevelMid  = 0.50;   // опорная точка совмещения
constexpr double kLevelHi   = 0.90;   // верх переднего фронта
constexpr double kLevelFall = 0.70;   // уровень, до которого меряется спад

struct CurveFeatures
{
    bool        valid = false;
    std::string reason;
    double peak   = 0.0;   // максимум v, м/с
    double t_peak = 0.0;   // момент пика, с
    double t10    = 0.0;   // пересечение 10 % пика на переднем фронте, с
    double t50    = 0.0;   // 50 % — опорная точка совмещения, с
    double t90    = 0.0;   // 90 %, с
    double rise   = 0.0;   // длительность фронта 10-90 %, с
    double fall   = 0.0;   // от пика до спада на 70 % пика, с
    bool   fall_valid = false;   // спад до 70 % может не попасть в окно записи
};

// Момент, когда кривая, идя вверх, ПОСЛЕДНИЙ раз перед пиком пересекает level.
// Поиск от пика назад: на реальном сигнале уровень пересекается не один раз
// (шум, упругий предвестник), и содержательно последнее пересечение перед пиком.
bool front_crossing(const Curve& c, size_t i_peak, double level, double& t_out)
{
    for (size_t i = i_peak; i-- > 0;)
    {
        if (c.v[i] <= level && c.v[i + 1] >= level)
        {
            const double den = c.v[i + 1] - c.v[i];
            t_out = (den > 0.0)
                        ? c.t[i] + ((level - c.v[i]) / den) * (c.t[i + 1] - c.t[i])
                        : c.t[i + 1];
            return true;
        }
    }
    return false;
}

// Первое пересечение уровня вниз ПОСЛЕ пика.
bool fall_crossing(const Curve& c, size_t i_peak, double level, double& t_out)
{
    for (size_t i = i_peak; i + 1 < c.v.size(); ++i)
    {
        if (c.v[i] >= level && c.v[i + 1] <= level)
        {
            const double den = c.v[i] - c.v[i + 1];
            t_out = (den > 0.0)
                        ? c.t[i] + ((c.v[i] - level) / den) * (c.t[i + 1] - c.t[i])
                        : c.t[i + 1];
            return true;
        }
    }
    return false;
}

CurveFeatures extract_features(const Curve& c)
{
    CurveFeatures f;
    if (c.t.size() < 2)
    {
        f.reason = "кривая содержит меньше двух точек";
        return f;
    }
    for (size_t i = 0; i < c.t.size(); ++i)
    {
        if (!std::isfinite(c.t[i]) || !std::isfinite(c.v[i]))
        {
            std::ostringstream os;
            os << "в кривой есть неконечное значение (точка #" << i << ")";
            f.reason = os.str();
            return f;
        }
    }
    size_t i_peak = 0;
    for (size_t i = 1; i < c.v.size(); ++i)
    {
        if (c.v[i] > c.v[i_peak]) i_peak = i;
    }
    f.peak   = c.v[i_peak];
    f.t_peak = c.t[i_peak];
    if (!(f.peak > 0.0))
    {
        f.reason = "у кривой нет положительного пика (сигнал не пришёл?)";
        return f;
    }
    if (i_peak == 0)
    {
        f.reason = "пик пришёлся на первую точку — переднего фронта в записи нет";
        return f;
    }
    if (!front_crossing(c, i_peak, kLevelLo * f.peak, f.t10) ||
        !front_crossing(c, i_peak, kLevelMid * f.peak, f.t50) ||
        !front_crossing(c, i_peak, kLevelHi * f.peak, f.t90))
    {
        f.reason = "передний фронт не прослеживается: кривая начинается выше "
                   "10 % пика или нарастание не пересекает уровни 10/50/90 %";
        return f;
    }
    f.rise = f.t90 - f.t10;
    if (!(f.rise > 0.0))
    {
        f.reason = "длительность фронта 10-90 % не положительна "
                   "(уровни пересечены не по порядку)";
        return f;
    }
    // Спад может не попасть в окно записи — это не повод отказываться
    // от сравнения амплитуды и фронта, поэтому признак необязателен.
    double t_fall = 0.0;
    if (fall_crossing(c, i_peak, kLevelFall * f.peak, t_fall))
    {
        f.fall       = t_fall - f.t_peak;
        f.fall_valid = (f.fall > 0.0);
    }
    f.valid = true;
    return f;
}

// ГРУБЫЕ ворота. Тест обязан падать только на том, что означает настоящую
// поломку (кривая вдвое не та), а не на количественном несовпадении с
// приборным профилем: расхождение по амплитуде объяснено недостаточным
// разрешением (notes/contact/08-phase-shift.md, §4), и ужесточать пороги следует
// по мере роста разрешения, а не подгонять их под текущий результат.
//
// Отношение амплитуд: порог ГРУБЫЙ. Текущий расчёт даёт ~0,66 при
// delta_x = 2e-4 м; ворота ловят деградацию вдвое. Ужесточать по мере роста
// разрешения — см. notes/contact/08-phase-shift.md, §4.1-4.3 и §6, п. 1.
constexpr double kAmpRatioLo  = 0.5;
constexpr double kAmpRatioHi  = 2.0;
// Отношение длительностей фронта 10-90 %: порог ГРУБЫЙ. Фронт расчёта
// размазан численной диссипацией примерно вдвое сильнее приборного; ворота
// ловят только кратное расхождение. Ужесточать по мере роста разрешения —
// см. notes/contact/08-phase-shift.md, §4.1.
constexpr double kRiseRatioLo = 1.0 / 3.0;
constexpr double kRiseRatioHi = 3.0;

struct ShapeMatch
{
    bool          comparable = false;
    std::string   reason;
    CurveFeatures calc;
    CurveFeatures ref;
    // Справочная величина: на сколько сдвинут расчёт при совмещении по t50.
    // НЕ ошибка расчёта — см. цитату из статьи выше.
    double shift            = 0.0;   // t50_ref - t50_calc, с
    double amp_ratio        = 0.0;   // расчёт / эталон
    double rise_ratio       = 0.0;   // расчёт / эталон
    double fall_ratio       = 0.0;   // расчёт / эталон
    bool   fall_ratio_valid = false;
    double peak_offset      = 0.0;   // (t_peak_calc + shift) - t_peak_ref, с
};

ShapeMatch compare_curves(const Curve& calc, const Curve& ref)
{
    ShapeMatch m;
    m.calc = extract_features(calc);
    m.ref  = extract_features(ref);
    if (!m.ref.valid)
    {
        m.reason = "эталон: " + m.ref.reason;
        return m;
    }
    if (!m.calc.valid)
    {
        m.reason = "расчёт: " + m.calc.reason;
        return m;
    }
    m.shift       = m.ref.t50 - m.calc.t50;
    m.amp_ratio   = m.calc.peak / m.ref.peak;
    m.rise_ratio  = m.calc.rise / m.ref.rise;
    m.peak_offset = (m.calc.t_peak + m.shift) - m.ref.t_peak;
    if (m.calc.fall_valid && m.ref.fall_valid)
    {
        m.fall_ratio       = m.calc.fall / m.ref.fall;
        m.fall_ratio_valid = true;
    }
    m.comparable = true;
    return m;
}

bool amp_gate_ok(const ShapeMatch& m)
{
    return m.comparable && m.amp_ratio >= kAmpRatioLo && m.amp_ratio <= kAmpRatioHi;
}

bool rise_gate_ok(const ShapeMatch& m)
{
    return m.comparable && m.rise_ratio >= kRiseRatioLo && m.rise_ratio <= kRiseRatioHi;
}

// Метрики печатаются ВСЕГДА, независимо от вердикта: без них расхождение
// не диагностируется, а «зелёный» прогон не проверить глазами.
void print_features(const char* who, const CurveFeatures& f)
{
    constexpr double us = 1e6;
    if (!f.valid)
    {
        std::printf("  info : признаки %s: не извлечены (%s)\n", who, f.reason.c_str());
        return;
    }
    std::printf("  info : признаки %s: пик %.1f м/с при %.4f мкс; "
                "t10/t50/t90 = %.4f/%.4f/%.4f мкс; фронт 10-90 %% = %.4f мкс; ",
                who, f.peak, f.t_peak * us, f.t10 * us, f.t50 * us, f.t90 * us,
                f.rise * us);
    if (f.fall_valid) std::printf("спад до 70 %% = %.4f мкс\n", f.fall * us);
    else              std::printf("спад до 70 %% = н/д (не попал в окно записи)\n");
}

void print_shape_report(const ShapeMatch& m)
{
    constexpr double us = 1e6;
    print_features("эталона", m.ref);
    print_features("расчёта", m.calc);
    if (!m.comparable) return;

    std::printf("  info : совмещение по t50 переднего фронта: расчёт сдвинут на %+.4f мкс\n",
                m.shift * us);
    std::printf("         ВНИМАНИЕ: сдвиг — СПРАВОЧНАЯ величина, он НЕ является ошибкой\n"
                "         расчёта. Начало отсчёта времени на графике статьи не привязано\n"
                "         к моменту удара: «the transit times of the first arrival was not\n"
                "         measured in these tests, and the profiles are not necessarily\n"
                "         aligned in time with expected elastic precursor velocities»\n"
                "         (Grady & Moody, SAND96-0551). Сравнивается ФОРМА кривой после\n"
                "         совмещения; разбор — notes/contact/08-phase-shift.md.\n");
    std::printf("  info : отношения расчёт/эталон: амплитуда = %.3f, фронт 10-90 %% = %.3f, ",
                m.amp_ratio, m.rise_ratio);
    if (m.fall_ratio_valid) std::printf("спад до 70 %% = %.3f\n", m.fall_ratio);
    else                    std::printf("спад до 70 %% = н/д\n");
    std::printf("  info : смещение пика после совмещения = %+.4f мкс\n", m.peak_offset * us);
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
        os << "sph_e2e_bc2_" << tag << "_" << static_cast<long long>(stamp)
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
    Bc2Geometry geom;
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
//   «cases/bc2/source/initial-conditions.txt» (архив кейса, cases/bc2/README.md).
// Архивный файл НЕ парсится намеренно: тест читает свою копию из tests/data/,
// чтобы не зависеть от архива и рабочего каталога.
fs::path config_path()
{
    return fs::path(SPH_CASES_DIR) / "bc2" / "regression" / "initial_conditions.txt";
}

// t_total физического кейса — из файла конфигурации, а не из литерала в коде:
// файл остаётся единственным источником истины по параметрам кейса.
double case_t_total()
{
    SimulationConfig probe;
    probe.apply_defaults();
    const fs::path p = config_path();
    if (!fs::exists(p))
    {
        throw std::runtime_error("не найден конфиг теста: " + p.string());
    }
    parse_sph_config(probe, p.string());
    return probe.t_total;
}

void apply_bc2_config(SimulationConfig& cfg, const RunSettings& rs)
{
    cfg.apply_defaults();

    // Конфигурация читается из файла теми же средствами, что и в боевом
    // запуске (parse_sph_config), — так тест покрывает реальный путь загрузки
    // параметров, а не обходит его присваиванием полей.
    const fs::path cfg_path = config_path();
    if (!fs::exists(cfg_path))
    {
        // parse_sph_config при неоткрывшемся файле молча оставляет умолчания:
        // без этой проверки тест считал бы «не тот» кейс и выглядел рабочим.
        throw std::runtime_error("не найден конфиг теста: " + cfg_path.string());
    }
    parse_sph_config(cfg, cfg_path.string());

    // Ниже — только то, что параметром физического кейса не является:
    //   t_total в smoke  — укорочение ради быстрого прогона в CI (ТЗ §5);
    //   contact_force_type — перебирается по сценариям §6.
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
    apply_bc2_config(cfg, rs);

    ParticleSystem sys;
    res.geom = build_bc2_geometry(sys, rs.delta_x);
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
    // одномерная ударная труба, к BC2 отношения не имеющая; для BC2 такой шаг
    // бессмысленен (в contact_force_type = 3 он вошёл бы в импульс первого
    // шага). Берём тот же критерий Куранта, которым конвейер пересчитывает
    // шаг дальше. Решение сверх E1-E3, зафиксировано здесь.
    // Порядок важен: c заполняется УРС, и только потом считается шаг —
    // иначе критерий Куранта остался бы без скорости звука (см. init_eos_state).
    init_eos_state(cfg, sys, ntotal);
    double dt = compute_cfl_timestep(cfg, sys, ntotal);
    double time = 0.0;

    // Сторож: бюджет времени + детектор застоя. Потолок «шагов до финиша»
    // 2e6 подобран так, чтобы не задевать здоровый расчёт (BC2 при dx = 2e-4
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

void report_geometry(const Bc2Geometry& g, bool full)
{
    const auto layers = bc2_layers();
    for (size_t L = 0; L < 4; ++L)
    {
        std::ostringstream os;
        os << "слой " << (L + 1) << " (" << layers[L].name << "): "
           << g.layer_nr[L] << " x " << g.layer_nz[L] << " = " << g.layer_count[L]
           << ", [" << g.layer_begin[L] << ", " << (g.layer_begin[L] + g.layer_count[L]) << ")"
           << ", dr = " << (g.layer_dr[L] * 1e3) << " мм, dz = " << (g.layer_dz[L] * 1e3) << " мм"
           << ", толщина = " << (static_cast<double>(g.layer_nz[L]) * g.layer_dz[L] * 1e3) << " мм";
        note(os.str());
    }
    {
        std::ostringstream os;
        os << "ntotal = " << g.ntotal << ", moni_particle = " << g.moni_particle
           << ", r(moni) = " << (0.5 * g.layer_dr[2] * 1e3) << " мм (ось)";
        note(os.str());
    }

    // Ключевая величина постановки: толщина ЛЕТЯЩЕЙ пластины B4C задаёт
    // длительность нагружающего импульса (статья: 3.917 мм). Узловая схема
    // до приведения к статье давала здесь 3.600 мм.
    {
        const double h2_mm = static_cast<double>(g.layer_nz[1]) * g.layer_dz[1] * 1e3;
        std::ostringstream os;
        os << "толщина слоя 2 (летящая B4C) = " << std::setprecision(10) << h2_mm
           << " мм (статья: 3.917 мм)";
        check(std::fabs(h2_mm - 3.917) < 1e-9, os.str());
    }

    if (full)
    {
        // §1.2: ожидаемое число частиц при delta_x = 2e-4.
        // Числа пересчитаны после перехода на узлы по центрам ячеек и снятия
        // зазоров между слоями (прежние: 6789 / 3287 / 7958 / 16129,
        // ntotal = 34163, moni = 12973 — они отвечали постановке с зазорами
        // и обрезанной летящей пластиной).
        check(g.layer_count[0] == 6570,  "§1.2 слой 1 = 6570 частиц");
        check(g.layer_count[1] == 3460,  "§1.2 слой 2 = 3460 частиц");
        check(g.layer_count[2] == 7785,  "§1.2 слой 3 = 7785 частиц");
        check(g.layer_count[3] == 16129, "§1.2 слой 4 = 16129 частиц");
        check(g.ntotal == 33944,         "§1.2 ntotal = 33944");
        check(g.moni_particle == 10074,  "§1.3 moni_particle = 10074 (ось, тыл слоя 3)");
    }
    else
    {
        // Регрессионный замок smoke-геометрии: значения получены этим же
        // правилом при delta_x = 8e-4 и зафиксированы, чтобы изменение
        // геометрии не прошло незамеченным.
        check(g.layer_count[0] == 440, "smoke слой 1 = 440 частиц");
        check(g.layer_count[1] == 215, "smoke слой 2 = 215 частиц");
        check(g.layer_count[2] == 473, "smoke слой 3 = 473 частицы");
        check(g.layer_count[3] == 1024, "smoke слой 4 = 1024 частицы");
        check(g.ntotal == 2152,        "smoke ntotal = 2152");
        check(g.moni_particle == 665,  "smoke moni_particle = 665");
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
    // штатно: сборка контактная, удар происходит при t = 0, но волне нужно
    // пройти покоящуюся керамику (9.033 мм) до монитора на её тыльной
    // поверхности. При продольной скорости звука в B4C ~1.4e4 м/с это ~6.4e-7 с,
    // то есть больше smoke-t_total = 3e-7 с.
    if (!r.curve.empty())
    {
        double amp = 0.0;
        for (double v : r.curve.v) amp = std::max(amp, std::fabs(v));
        // Продольная скорость звука в B4C по упругим константам Eagle Picher
        // (materials.hpp): sqrt((K + 4G/3)/rho) = sqrt((230.5e9 + 4*199.4e9/3)/2517).
        constexpr double kCl_B4C = 14043.0;   // м/с
        constexpr double kCeramicThickness = 0.009033;   // м, покоящаяся керамика
        std::ostringstream os;
        os << "max|vx[1][moni]| по расчёту = " << amp << " м/с"
           << "; оценка времени прихода волны к монитору = "
           << (kCeramicThickness / kCl_B4C) << " с";
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
// Проверяется ограничение ТОГО типа контакта, за который отвечает этот прогон.
// Прежде функция перебирала типы 1..3 сама; теперь перебор делает набор тестов,
// по одному прогону на тип, и дублировать его внутри незачем.
void test_known_limitation_b4c_b4c(const Bc2Geometry& g, int contact_type)
{
    std::printf("\n[§6] Известное ограничение D2=A: контакт не различает два тела одного материала\n");

    // (а) доказательство на самой геометрии: слои 2 и 3 целиком itype = 10,
    //     поэтому критерий itype[i] != itype[j] не может сработать ни на одной
    //     паре между ними — независимо от их положений на любом шаге.
    const auto layers = bc2_layers();
    check(layers[1].itype == 10 && layers[2].itype == 10,
          "§6 слои 2 и 3 — оба itype = 10, критерий D2=A между ними неприменим");

    // (б) операционное подтверждение на ext_force: две перекрывающиеся
    //     частицы B4C дают ноль контактных пар, те же частицы с разными
    //     материалами — одну. При contact_force_type = 0 контакта нет вовсе,
    //     и проверять на нём нечего.
    if (contact_type == contact::Off)
    {
        note("§6 contact_force_type = 0: контакт выключен, операционная проверка неприменима");
    }
    else
    {
        const int type = contact_type;
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
            sys.rho[i] = 2517.0;   // B4C по статье (см. bc2_layers)
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
        ext_force(cfg, sys, 1e-9, 2, 2, 1, pi, pj, dvxdt, wb);
        const size_t pairs_same = wb.contact_pair_count;

        // Та же геометрия, но разные материалы (ударник 8 <-> керамика 10).
        sys.itype[0] = 8; sys.itype[1] = 10;
        ext_force(cfg, sys, 1e-9, 2, 2, 1, pi, pj, dvxdt, wb);
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

// Шаг решётки режима. Отдельной функцией, потому что его спрашивают из двух
// мест: настройка прогона и сборка сводного графика, где расчёта нет, а
// геометрия всё равно нужна для подписи. Разъехаться они не должны.
double case_delta_x(bool full)
{
    return full ? 2e-4 : 8e-4;
}

fs::path reference_path(bool full)
{
    return fs::path(SPH_CASES_DIR) / "bc2" / "regression" / (full ? "reference.txt" : "reference_smoke.txt");
}

void emit_reference(const RunResult& r, const RunSettings& rs, bool full)
{
    std::printf("# Эталон для end-to-end теста BC2 (удар по карбиду бора)\n");
    std::printf("#\n");
    std::printf("# СГЕНЕРИРОВАНО: test_e2e_bc2 %s--emit-reference\n", full ? "--full " : "");
    std::printf("# Формат: <время, с>  <vx[1][moni_particle], м/с>  <itype[moni_particle]>\n");
    std::printf("# (ровно формат time_integration.cpp и BC2/log_*.txt)\n");
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
    std::printf("# РЕЖИМ СВЕРКИ: сравнение ФОРМЫ после совмещения по t50 переднего фронта\n");
    std::printf("#   Начало отсчёта времени эталона не привязано к моменту удара\n");
    std::printf("#   (Grady & Moody, SAND96-0551), поэтому сверка в общем абсолютном\n");
    std::printf("#   времени неприменима. Сравниваются безразмерные отношения\n");
    std::printf("#   расчёт/эталон: амплитуда, фронт 10-90 %%, спад до 70 %% пика.\n");
    std::printf("#   ГРУБЫЕ ворота (notes/contact/08-phase-shift.md):\n");
    std::printf("#     отношение амплитуд          в [%g, %g]\n", kAmpRatioLo, kAmpRatioHi);
    std::printf("#     отношение фронтов 10-90 %%   в [%g, %g]\n", kRiseRatioLo, kRiseRatioHi);
    std::printf("#   Величина совмещения печатается справочно и в вердикт не входит.\n");
    std::printf("#\n");
    std::printf("# ДАННЫЕ (время  vx[1]  itype)\n");
    for (size_t i = 0; i < r.curve.t.size(); ++i)
    {
        std::printf("%.17g %.17g %d\n", r.curve.t[i], r.curve.v[i], r.curve.itype[i]);
    }
}

// ===========================================================================
// 7a. График «расчёт против эталона» (--emit-plot)
// ===========================================================================
//
// Рисуется из тех же данных, что уже лежат в памяти теста: r.curve — кривая
// прогона, эталон — reference_path(full). Отдельного разбора файлов и внешнего
// скрипта не требуется, поэтому график не может разойтись с тем, что проверил
// компаратор.
//
// Включается ТОЛЬКО флагом --emit-plot=<путь>. По умолчанию выключено
// намеренно: тест следит за отсутствием побочных эффектов в рабочем каталоге
// (§7, check_no_side_effects), обязан работать на checkout только для чтения и
// запускаться параллельно. Запись в дерево исходников — осознанное действие
// владельца, а не следствие прогона тестов.

// Число в русской типографике: запятая как десятичный разделитель и неразрывные
// пробелы по три разряда. График — часть русскоязычного отчёта, и «0.2 мм»
// рядом с «0,2 мм» в тексте выглядит как разные величины.
std::string ru_num(double value, int precision)
{
    std::ostringstream s;
    s.setf(std::ios::fixed);
    s.precision(precision);
    s << value;
    std::string out = s.str();

    const size_t dot = out.find('.');
    std::string head = (dot == std::string::npos) ? out : out.substr(0, dot);
    const std::string tail = (dot == std::string::npos) ? std::string() : out.substr(dot + 1);

    const bool negative = !head.empty() && head[0] == '-';
    if (negative) head.erase(0, 1);
    for (size_t pos = head.size(); pos > 3;)
    {
        pos -= 3;
        head.insert(pos, " ");   // неразрывный пробел
    }
    if (negative) head.insert(0, "-");

    return tail.empty() ? head : head + "," + tail;
}

// Согласование существительного с числом: 1 частица, 2 частицы, 5 частиц.
// Нужно потому, что ntotal и число шагов меняются от прогона к прогону, и
// зашитая форма («34 163 частицы») при другом delta_x станет неверной.
const char* ru_plural(size_t n, const char* one, const char* few, const char* many)
{
    const size_t mod100 = n % 100;
    if (mod100 >= 11 && mod100 <= 14) return many;
    switch (n % 10)
    {
        case 1:          return one;
        case 2: case 3:
        case 4:          return few;
        default:         return many;
    }
}

// Прореживание: на N интервалов, в каждом сохраняются первая точка и точки
// минимума/максимума. Полка и пики выживают, файл не раздувается.
Curve thin_curve(const Curve& c, size_t buckets = 500)
{
    if (c.t.size() <= buckets) return c;
    Curve out;
    const double t0 = c.t.front(), t1 = c.t.back();
    const double step = (t1 - t0) / static_cast<double>(buckets);
    if (!(step > 0.0)) return c;

    size_t k = 0;
    double edge = t0 + step;
    while (k < c.t.size())
    {
        size_t lo = k, hi = k, first = k;
        while (k < c.t.size() && c.t[k] <= edge)
        {
            if (c.v[k] < c.v[lo]) lo = k;
            if (c.v[k] > c.v[hi]) hi = k;
            ++k;
        }
        if (k == first)   // ни одна точка не попала в интервал — двигаем границу
        {
            edge += step;
            continue;
        }
        std::array<size_t, 3> pick{first, lo, hi};
        std::sort(pick.begin(), pick.end());
        size_t prev = static_cast<size_t>(-1);
        for (size_t idx : pick)
        {
            if (idx == prev) continue;
            prev = idx;
            out.t.push_back(c.t[idx]);
            out.v.push_back(c.v[idx]);
            out.itype.push_back(c.itype[idx]);
        }
        edge += step;
    }
    // Последняя точка кривой должна попасть на график всегда.
    if (!out.t.empty() && out.t.back() != c.t.back())
    {
        out.t.push_back(c.t.back());
        out.v.push_back(c.v.back());
        out.itype.push_back(c.itype.back());
    }
    return out;
}

// Серия графика — кривая одного сценария контакта.
//
// Сценарии разнесены по ОТДЕЛЬНЫМ прогонам теста (по одному на
// contact_force_type), поэтому в одном процессе все четыре кривые не живут.
// Каждый прогон выгружает свою кривую в файл (--emit-curve), а рисование —
// отдельный дешёвый шаг без расчёта (--plot-from). Связь «график не расходится
// с метриками отчёта» сохраняется: в файл пишется ровно тот объект Curve,
// который проверил компаратор, без пересчёта и без повторного чтения log.txt.
struct PlotSeries
{
    Curve  curve;
    int    type    = 0;      // contact_force_type
    bool   aborted = false;  // прогон не дошёл до t_total
    size_t steps   = 0;
};

// Имя файла выгрузки внутри каталога обмена. Фиксировано, потому что по нему
// --plot-from собирает серии, не получая списка файлов.
fs::path curve_dump_path(const fs::path& dir, int type)
{
    return dir / ("bc2_curve_" + std::to_string(type) + ".txt");
}

// Формат — тот же, что у эталона и log.txt: <время> <скорость> <itype>.
// Признаки прогона идут заголовочными комментариями, которые read_curve
// пропускает, а load_series читает.
bool write_curve_dump(const PlotSeries& s, const fs::path& path, std::string& err)
{
    std::error_code ec;
    if (path.has_parent_path()) fs::create_directories(path.parent_path(), ec);
    std::ofstream f(path, std::ios::trunc);
    if (!f)
    {
        err = "не удалось открыть " + path.string();
        return false;
    }
    f << "# Кривая BC2 для сводного графика (см. --plot-from)\n"
      << "# contact_force_type = " << s.type << "\n"
      << "# steps = " << s.steps << "\n"
      << "# aborted = " << (s.aborted ? 1 : 0) << "\n"
      << "# Формат: <время, с>  <vx[1][moni_particle], м/с>  <itype[moni_particle]>\n";
    f.precision(12);
    for (size_t i = 0; i < s.curve.t.size(); ++i)
    {
        f << s.curve.t[i] << ' ' << s.curve.v[i] << ' ' << s.curve.itype[i] << '\n';
    }
    f.close();
    if (!f)
    {
        err = "ошибка записи " + path.string();
        return false;
    }
    return true;
}

// Разбор одного числа из строки вида "# ключ = значение".
bool header_value(const std::string& line, const char* key, long long& out)
{
    const std::string pat = std::string("# ") + key + " = ";
    if (line.rfind(pat, 0) != 0) return false;
    try { out = std::stoll(line.substr(pat.size())); }
    catch (const std::exception&) { return false; }
    return true;
}

bool load_series(const fs::path& path, PlotSeries& out, std::string& err)
{
    out.curve = read_curve(path, err);
    if (!err.empty()) return false;
    if (out.curve.empty())
    {
        err = "в " + path.string() + " нет точек кривой";
        return false;
    }
    std::ifstream f(path);
    std::string line;
    long long v = 0;
    while (std::getline(f, line))
    {
        if (!line.empty() && line[0] != '#') break;
        if (header_value(line, "contact_force_type", v)) out.type = static_cast<int>(v);
        else if (header_value(line, "steps", v))         out.steps = static_cast<size_t>(v);
        else if (header_value(line, "aborted", v))       out.aborted = (v != 0);
    }
    return true;
}

// Названия типов для легенды. Короткие: легенда идёт в одну строку.
const char* contact_type_name(int type)
{
    switch (type)
    {
        case 0:  return "выключен";
        case 1:  return "Пинбол";
        case 2:  return "Свигл";
        case 3:  return "импульс";
        default: return "?";
    }
}

// Геометрия передаётся отдельно от прогона: в режиме --plot-from расчёта нет
// вовсе, а подпись графика всё равно обязана назвать разрешение и число частиц.
void emit_plot(const Bc2Geometry& geom, const std::vector<PlotSeries>& series,
               bool full, const fs::path& out_path)
{
    std::printf("\n[график] %s\n", out_path.string().c_str());

    bool any_curve = false;
    for (const PlotSeries& s : series) any_curve = any_curve || !s.curve.empty();
    if (!any_curve)
    {
        std::printf("  info : кривых расчёта нет, график не построен\n");
        return;
    }

    std::string err;
    Curve ref;
    const fs::path ref_path = reference_path(full);
    if (fs::exists(ref_path))
    {
        ref = read_curve(ref_path, err);
        if (!err.empty())
        {
            std::printf("  info : эталон не прочитан (%s), рисуется только расчёт\n", err.c_str());
            ref = Curve{};
        }
    }
    else
    {
        std::printf("  info : эталона нет (%s), рисуется только расчёт\n", ref_path.string().c_str());
    }

    std::vector<PlotSeries> plot = series;
    for (PlotSeries& s : plot) s.curve = thin_curve(s.curve);

    // Стиль серии. Одного цвета мало: типы 0 и 1 на этом кейсе совпадают до
    // четвёртой значащей цифры, и сплошная линия, нарисованная позже, закрывает
    // раннюю целиком — легенда обещала бы четыре кривые, а на поле их было бы
    // три. Поэтому к цвету добавлены штрих и толщина: нижняя сплошная
    // просматривается сквозь пунктир верхней. Побочно это даёт второй признак
    // помимо цвета, что нужно и при цветовой слепоте.
    // Ключ — ТИП, а не место в списке: если какого-то типа нет (прогон не
    // делали), остальные не должны сменить штрих. Ровно та же причина, по
    // которой цвет берётся из типа.
    struct SeriesStyle { double width; const char* dash; };
    const auto style_of = [](int type) -> SeriesStyle
    {
        switch (type & 3)
        {
            case 0:  return {3.0, ""};
            case 1:  return {2.4, "9 6"};
            case 2:  return {2.0, ""};
            default: return {2.2, "2.5 4"};
        }
    };

    // Пределы осей. Потолок по Y — ближайшее сверху кратное 300, чтобы линии
    // сетки попадали на круглые значения и кривая не упиралась в рамку.
    double vmax = 0.0, tmax = 0.0;
    for (const PlotSeries& s : plot)
    {
        for (size_t i = 0; i < s.curve.t.size(); ++i)
        {
            vmax = std::max(vmax, s.curve.v[i]);
            tmax = std::max(tmax, s.curve.t[i]);
        }
    }
    for (size_t i = 0; i < ref.t.size();  ++i) { vmax = std::max(vmax, ref.v[i]);  tmax = std::max(tmax, ref.t[i]);  }
    const double kGrid = 300.0;
    double ymax = std::ceil(std::max(vmax, 1.0) * 1.08 / kGrid) * kGrid;
    const double xmax_us = std::max(tmax * 1e6, 1e-9) * 1.02;

    // Поле справа расширено под прямые подписи кривых: в светлой теме аква и
    // жёлтый идут ниже 3:1 к фону, и подпись у конца кривой — обязательная
    // страховка, а не украшение. Сверху добавлено место под строку легенды.
    constexpr double W = 880.0, H = 470.0;
    constexpr double ML = 68.0, MR = 104.0, MT = 84.0, MB = 56.0;
    const double PW = W - ML - MR, PH = H - MT - MB;

    const auto sx = [&](double t_us) { return ML + t_us / xmax_us * PW; };
    const auto sy = [&](double v) { return MT + PH - std::min(std::max(v, 0.0), ymax) / ymax * PH; };

    std::ostringstream os;
    os.setf(std::ios::fixed);
    const auto path_of = [&](const Curve& c)
    {
        std::ostringstream d;
        d.setf(std::ios::fixed);
        d.precision(1);
        for (size_t i = 0; i < c.t.size(); ++i)
        {
            d << (i == 0 ? 'M' : 'L') << sx(c.t[i] * 1e6) << ',' << sy(c.v[i]) << ' ';
        }
        return d.str();
    };

    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 880 470\" width=\"880\" height=\"470\" "
          "role=\"img\" aria-label=\"Осевая скорость частицы-монитора BC2: расчёт при "
          "contact_force_type 0-3 и приборный эталон\">\n";
    // Палитра задаётся токенами на :root-уровне svg и переопределяется для
    // тёмной темы, поэтому график читается в обеих.
    os << "<style>\n"
          "svg { --surface-1:#fcfcfb; --text-primary:#0b0b0b; --text-secondary:#52514e;\n"
          "      --grid:#e4e3df; --axis:#b9b8b2; --ref:#0b0b0b;\n"
          "      --s0:#2a78d6; --s1:#eb6834; --s2:#1baf7a; --s3:#eda100;\n"
          "      font-family:-apple-system,BlinkMacSystemFont,\"Segoe UI\",Roboto,sans-serif; }\n"
          "@media (prefers-color-scheme: dark) {\n"
          "  svg { --surface-1:#1a1a19; --text-primary:#ffffff; --text-secondary:#c3c2b7;\n"
          "        --grid:#333330; --axis:#5a5955; --ref:#ffffff;\n"
          "        --s0:#3987e5; --s1:#d95926; --s2:#199e70; --s3:#c98500; }\n"
          "}\n"
          ".ttl { fill:var(--text-primary); font-size:15px; font-weight:600; }\n"
          ".sub { fill:var(--text-secondary); font-size:11.5px; }\n"
          ".axl { fill:var(--text-secondary); font-size:11px; }\n"
          ".tik { fill:var(--text-secondary); font-size:10.5px; }\n"
          ".lbl { fill:var(--text-primary); font-size:11px; font-weight:600; }\n"
          ".leg { fill:var(--text-secondary); font-size:11px; }\n"
          ".ln  { fill:none; stroke-width:2; stroke-linejoin:round; stroke-linecap:round; }\n"
          "</style>\n";
    os.precision(1);
    os << "<rect width=\"880\" height=\"470\" fill=\"var(--surface-1)\"/>\n";
    os << "<text class=\"ttl\" x=\"" << ML << "\" y=\"24\">Осевая скорость частицы-монитора, кейс BC2</text>\n";
    os << "<text class=\"sub\" x=\"" << ML << "\" y=\"41\">"
       << "Тыльная поверхность неподвижной пластины B4C, r = "
       << ru_num(0.5 * geom.layer_dr[2] * 1e3, 2) << " мм (ось). Расчёт: delta_x = "
       << ru_num(geom.delta_x * 1e3, 1) << " мм, "
       << ru_num(static_cast<double>(geom.ntotal), 0) << " "
       << ru_plural(geom.ntotal, "частица", "частицы", "частиц")
       << ". Крестик — прогон оборван до t_total.</text>\n";

    // Легенда — строкой над полем, колонки с постоянным шагом. Текст носит
    // текстовые токены, цвет несёт только штрих слева от подписи, поэтому
    // идентичность не держится на одном цвете.
    {
        constexpr double kLegY = 62.0, kLegPitch = 132.0, kLegRule = 22.0, kLegGap = 6.0;
        // Штрих и толщина в легенде — ровно те же, что у кривой на поле:
        // иначе легенда опознаёт серию не по тому признаку, по которому её
        // приходится узнавать на графике.
        const auto legend_entry = [&](double lx, const std::string& stroke,
                                      const SeriesStyle& st, bool faded,
                                      const std::string& text)
        {
            os << "<line x1=\"" << lx << "\" y1=\"" << kLegY
               << "\" x2=\"" << (lx + kLegRule) << "\" y2=\"" << kLegY
               << "\" stroke=\"" << stroke << "\" stroke-width=\"" << st.width << "\"";
            if (st.dash[0] != '\0') os << " stroke-dasharray=\"" << st.dash << "\"";
            if (faded) os << " opacity=\"0.85\"";
            os << "/>\n";
            os << "<text class=\"leg\" x=\"" << (lx + kLegRule + kLegGap)
               << "\" y=\"" << (kLegY + 3.8) << "\">" << text << "</text>\n";
        };
        double lx = ML;
        for (size_t i = 0; i < plot.size(); ++i)
        {
            std::ostringstream t;
            t << "тип " << plot[i].type << " — " << contact_type_name(plot[i].type);
            legend_entry(lx, "var(--s" + std::to_string(plot[i].type) + ")",
                         style_of(plot[i].type), false, t.str());
            lx += kLegPitch;
        }
        if (!ref.empty()) legend_entry(lx, "var(--ref)", SeriesStyle{2.0, "6 4"}, true,
                                       "эталон (прибор)");
    }

    for (double v = 0.0; v <= ymax + 1e-9; v += kGrid)
    {
        const double y = sy(v);
        os << "<line x1=\"" << ML << "\" y1=\"" << y << "\" x2=\"" << (ML + PW) << "\" y2=\"" << y
           << "\" stroke=\"var(--grid)\" stroke-width=\"1\"/>\n";
        os << "<text class=\"tik\" x=\"" << (ML - 10) << "\" y=\"" << (y + 3.6)
           << "\" text-anchor=\"end\">" << ru_num(v, 0) << "</text>\n";
    }
    for (double t = 0.0; t <= xmax_us + 1e-9; t += 0.5)
    {
        os << "<text class=\"tik\" x=\"" << sx(t) << "\" y=\"" << (MT + PH + 18)
           << "\" text-anchor=\"middle\">" << ru_num(t, 1) << "</text>\n";
    }
    os << "<line x1=\"" << ML << "\" y1=\"" << (MT + PH) << "\" x2=\"" << (ML + PW) << "\" y2=\"" << (MT + PH)
       << "\" stroke=\"var(--axis)\" stroke-width=\"1\"/>\n";
    os << "<text class=\"axl\" x=\"" << (ML + PW / 2) << "\" y=\"" << (H - 14)
       << "\" text-anchor=\"middle\">время, мкс</text>\n";
    os << "<text class=\"axl\" x=\"18\" y=\"" << (MT + PH / 2) << "\" text-anchor=\"middle\" "
       << "transform=\"rotate(-90 18 " << (MT + PH / 2) << ")\">осевая скорость, м/с</text>\n";

    if (!ref.empty())
    {
        os << "<path class=\"ln\" d=\"" << path_of(ref)
           << "\" stroke=\"var(--ref)\" stroke-dasharray=\"6 4\" opacity=\"0.85\"/>\n";
    }
    // Кривые расчёта: стиль — из style_of, тот же, что в легенде.
    for (size_t i = 0; i < plot.size(); ++i)
    {
        if (plot[i].curve.empty()) continue;
        const SeriesStyle st = style_of(plot[i].type);
        os << "<path class=\"ln\" d=\"" << path_of(plot[i].curve)
           << "\" stroke=\"var(--s" << plot[i].type << ")\" stroke-width=\"" << st.width << "\"";
        if (st.dash[0] != '\0') os << " stroke-dasharray=\"" << st.dash << "\"";
        os << "/>\n";
    }

    // Прямые подписи у концов кривых с разведением по вертикали: цвет в
    // светлой теме ниже 3:1 у двух слотов, и подпись — обязательная страховка.
    struct EndLabel
    {
        double y_end;     // конец кривой, куда ставится крестик обрыва
        double y_label;   // та же точка после разведения подписей
        double x_end;
        int    type;
        bool   aborted;
    };
    std::vector<EndLabel> ends;
    for (const PlotSeries& s : plot)
    {
        if (s.curve.empty()) continue;
        const double ye = sy(s.curve.v.back());
        ends.push_back({ye, ye, sx(s.curve.t.back() * 1e6), s.type, s.aborted});
    }
    std::sort(ends.begin(), ends.end(),
              [](const EndLabel& a, const EndLabel& b) { return a.y_label < b.y_label; });
    // Разведение по вертикали в два прохода. Один проход вниз недостаточен:
    // когда кривые кончаются рядом у нижнего края (smoke-режим — все четыре
    // у нуля), подписи уезжают за рамку на подписи оси времени. Обратный
    // проход возвращает пачку внутрь поля.
    constexpr double kLabelMin = 14.0;
    const double y_lo = MT + 6.0, y_hi = MT + PH - 6.0;
    for (size_t i = 1; i < ends.size(); ++i)
    {
        if (ends[i].y_label - ends[i - 1].y_label < kLabelMin)
        {
            ends[i].y_label = ends[i - 1].y_label + kLabelMin;
        }
    }
    if (!ends.empty() && ends.back().y_label > y_hi)
    {
        ends.back().y_label = y_hi;
        for (size_t i = ends.size() - 1; i > 0; --i)
        {
            if (ends[i].y_label - ends[i - 1].y_label < kLabelMin)
            {
                ends[i - 1].y_label = ends[i].y_label - kLabelMin;
            }
        }
        for (EndLabel& e : ends) e.y_label = std::max(e.y_label, y_lo);
    }
    for (const EndLabel& e : ends)
    {
        // Крестик в точке обрыва: обрыв — свойство данных и обязан быть виден
        // на самом графике, а не только в подписи сбоку.
        if (e.aborted)
        {
            os << "<path d=\"M" << (e.x_end - 4) << "," << (e.y_end - 4)
               << "L" << (e.x_end + 4) << "," << (e.y_end + 4)
               << "M" << (e.x_end + 4) << "," << (e.y_end - 4)
               << "L" << (e.x_end - 4) << "," << (e.y_end + 4)
               << "\" stroke=\"var(--s" << e.type << ")\" stroke-width=\"2\" fill=\"none\"/>\n";
        }
        os << "<line x1=\"" << (ML + PW + 6) << "\" y1=\"" << e.y_label
           << "\" x2=\"" << (ML + PW + 24) << "\" y2=\"" << e.y_label
           << "\" stroke=\"var(--s" << e.type << ")\" stroke-width=\"2.5\"/>\n";
        os << "<text class=\"lbl\" x=\"" << (ML + PW + 30) << "\" y=\"" << (e.y_label + 3.8)
           << "\">тип " << e.type << (e.aborted ? " ✕" : "") << "</text>\n";
    }
    os << "</svg>\n";

    std::error_code ec;
    if (out_path.has_parent_path()) fs::create_directories(out_path.parent_path(), ec);
    std::ofstream f(out_path, std::ios::binary | std::ios::trunc);
    if (!f)
    {
        std::printf("  FAIL : не удалось открыть файл для записи графика\n");
        ++failures;
        return;
    }
    const std::string svg = os.str();
    f << svg;
    f.close();
    if (!f)
    {
        std::printf("  FAIL : ошибка записи графика\n");
        ++failures;
        return;
    }
    size_t pts = 0;
    for (const PlotSeries& s : plot) pts += s.curve.t.size();
    std::printf("  ok   : график записан (%zu %s, %zu точек расчёта, %zu эталона, %zu байт)\n",
                plot.size(), ru_plural(plot.size(), "кривая", "кривые", "кривых"),
                pts, ref.t.size(), svg.size());
}

// ===========================================================================
// 7b. Самопроверка компаратора
// ===========================================================================
// Требование отдельное от «тест не падает»: компаратор — измерительный
// инструмент, и он сам обязан быть проверен. Проверки работают ТОЛЬКО
// с данными эталона (прогон расчёта не нужен), поэтому выполняются в обоих
// режимах и стоят почти ничего.
//
// Проверяется ровно то, на чём стоит новая схема сверки:
//   * тождество            — кривая против самой себя даёт отношения 1,0 и
//                            нулевое смещение;
//   * чувствительность к амплитуде — множитель 3 обязан провалить ворота;
//   * НЕчувствительность к сдвигу времени — сдвиг на 5 мкс не меняет ничего
//     (суть схемы: начало отсчёта эталона произвольно, см. цитату в §2);
//   * чувствительность к растяжению — растяжение времени в 5 раз обязано
//     провалить ворота по фронту.
void test_comparator_selfcheck()
{
    std::printf("\n[§4.3] Самопроверка компаратора (на данных эталона)\n");

    const fs::path ref_path = reference_path(true);
    if (!fs::exists(ref_path))
    {
        check(false, "самопроверка компаратора: эталонный файл отсутствует: " +
                         ref_path.string());
        return;
    }
    std::string err;
    const Curve ref = read_curve(ref_path, err);
    if (!err.empty() || ref.t.size() < 2)
    {
        check(false, "самопроверка компаратора: эталон не прочитан" +
                         (err.empty() ? std::string() : (": " + err)));
        return;
    }

    const CurveFeatures f = extract_features(ref);
    if (!check(f.valid, "признаки эталона извлекаются"))
    {
        std::printf("         причина: %s\n", f.reason.c_str());
        return;
    }
    print_features("эталона", f);

    // Допуски самопроверки — чисто арифметические: сдвиг и растяжение времени
    // выполняются в double, поэтому точное равенство ожидать нельзя нигде,
    // кроме тождества (там входные данные побитово те же).
    constexpr double kTolRatio = 1e-9;    // безразмерные отношения
    constexpr double kTolTime  = 1e-12;   // с (1e-6 мкс)

    // 1. Тождество.
    {
        const ShapeMatch m = compare_curves(ref, ref);
        check(m.comparable && m.amp_ratio == 1.0 && m.rise_ratio == 1.0 &&
                  m.fall_ratio_valid && m.fall_ratio == 1.0 &&
                  m.shift == 0.0 && m.peak_offset == 0.0,
              "тождество: эталон против самого себя даёт отношения ровно 1,0 "
              "и нулевое смещение");
        if (m.comparable)
        {
            std::printf("         амплитуда = %.17g, фронт = %.17g, спад = %.17g, "
                        "сдвиг = %.17g с, смещение пика = %.17g с\n",
                        m.amp_ratio, m.rise_ratio, m.fall_ratio, m.shift, m.peak_offset);
        }
    }

    // 2. Чувствительность к амплитуде: x3 обязано провалить ворота.
    {
        Curve scaled = ref;
        for (double& v : scaled.v) v *= 3.0;
        const ShapeMatch m = compare_curves(scaled, ref);
        check(m.comparable && !amp_gate_ok(m),
              "чувствительность к амплитуде: кривая x3 проваливает ворота по амплитуде");
        // Масштаб по величине не должен затрагивать временные признаки:
        // уровни берутся в долях пика, поэтому фронт обязан остаться прежним.
        check(m.comparable && rise_gate_ok(m) &&
                  std::fabs(m.rise_ratio - 1.0) <= kTolRatio,
              "избирательность: множитель по амплитуде не меняет длительность фронта");
        std::printf("         амплитуда = %.6f (ворота [%.3f, %.3f]), фронт = %.6f\n",
                    m.amp_ratio, kAmpRatioLo, kAmpRatioHi, m.rise_ratio);
    }

    // 3. НЕчувствительность к сдвигу времени — суть новой схемы сверки.
    {
        constexpr double kShift = 5e-6;   // произвольные 5 мкс
        Curve shifted = ref;
        for (double& t : shifted.t) t += kShift;
        const ShapeMatch m = compare_curves(shifted, ref);
        check(m.comparable &&
                  std::fabs(m.amp_ratio - 1.0) <= kTolRatio &&
                  std::fabs(m.rise_ratio - 1.0) <= kTolRatio &&
                  m.fall_ratio_valid && std::fabs(m.fall_ratio - 1.0) <= kTolRatio &&
                  std::fabs(m.peak_offset) <= kTolTime &&
                  amp_gate_ok(m) && rise_gate_ok(m),
              "нечувствительность к сдвигу времени: кривая, сдвинутая на 5 мкс, "
              "даёт те же отношения");
        if (m.comparable)
        {
            std::printf("         амплитуда = %.12f, фронт = %.12f, спад = %.12f, "
                        "смещение пика = %.3g с; справочный сдвиг = %.6f мкс\n",
                        m.amp_ratio, m.rise_ratio, m.fall_ratio, m.peak_offset,
                        m.shift * 1e6);
        }
    }

    // 4. Чувствительность к растяжению времени: x5 обязано провалить фронт.
    {
        constexpr double kStretch = 5.0;
        Curve stretched = ref;
        const double t0 = ref.t.front();
        for (double& t : stretched.t) t = t0 + kStretch * (t - t0);
        const ShapeMatch m = compare_curves(stretched, ref);
        check(m.comparable && !rise_gate_ok(m),
              "чувствительность к растяжению: кривая, растянутая по времени в 5 раз, "
              "проваливает ворота по фронту");
        // Растяжение по времени не трогает величины, значит амплитудные ворота
        // обязаны остаться пройденными: ворота различают род дефекта.
        check(m.comparable && amp_gate_ok(m) &&
                  std::fabs(m.amp_ratio - 1.0) <= kTolRatio,
              "избирательность: растяжение по времени не меняет отношение амплитуд");
        std::printf("         фронт = %.6f (ворота [%.3f, %.3f]), амплитуда = %.6f\n",
                    m.rise_ratio, kRiseRatioLo, kRiseRatioHi, m.amp_ratio);
    }
}

void check_against_reference(const RunResult& r, bool full)
{
    // Эталон — приборная запись, начинающаяся в 1.51 мкс. Smoke считает до
    // 0.3 мкс, поэтому его окно с эталоном не пересекается ни при каких
    // условиях: сверять нечего, и падать тут не на чем. Назначение smoke —
    // живучесть конвейера (ТЗ §5), она проверена выше в §4.4.
    if (!full)
    {
        std::printf("\n[§4] Сверка с эталоном: пропущена\n");
        std::printf("  info : smoke считает до t_total = 3e-7 с, волна за это время до монитора\n"
                    "         не доходит: у кривой нет ни переднего фронта, ни пика, то есть\n"
                    "         признаков формы, по которым идёт сверка, не существует. Сверка\n"
                    "         формы выполняется в режиме --full; сам компаратор проверен выше\n"
                    "         самопроверкой, она работает в обоих режимах.\n");
        return;
    }

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
            "           2) ./test_e2e_bc2 %s--emit-reference > %s\n"
            "           3) вписать в шапку commit/дату/автора и просмотреть кривую глазами.\n",
            ref_path.string().c_str(), full ? "--full " : "", ref_path.string().c_str());
        ++failures;
        return;
    }

    const ShapeMatch m = compare_curves(r.curve, ref);
    // Метрики печатаются всегда — и при провале, и при успехе.
    print_shape_report(m);
    if (!m.comparable)
    {
        std::printf("  FAIL : сравнение формы невозможно: %s\n", m.reason.c_str());
        std::printf("         расчёт: %zu точек, t в [%g, %g]\n",
                    r.curve.t.size(), r.curve.t.front(), r.curve.t.back());
        std::printf("         эталон: %zu точек, t в [%g, %g]\n",
                    ref.t.size(), ref.t.front(), ref.t.back());
        ++failures;
        return;
    }

    const bool ok_amp  = amp_gate_ok(m);
    const bool ok_rise = rise_gate_ok(m);
    std::printf("  info : ворота (ГРУБЫЕ, см. notes/contact/08-phase-shift.md): "
                "амплитуда в [%.3f, %.3f] — %s; фронт в [%.3f, %.3f] — %s\n",
                kAmpRatioLo, kAmpRatioHi, ok_amp ? "ok" : "НАРУШЕНЫ",
                kRiseRatioLo, kRiseRatioHi, ok_rise ? "ok" : "НАРУШЕНЫ");
    if (!ok_amp)
    {
        std::printf("  FAIL : отношение амплитуд %.3f вне грубых ворот [%.3f, %.3f] —\n"
                    "         расчёт отличается от приборного профиля кратно, это поломка,\n"
                    "         а не количественная неточность разрешения\n",
                    m.amp_ratio, kAmpRatioLo, kAmpRatioHi);
        ++failures;
    }
    if (!ok_rise)
    {
        std::printf("  FAIL : отношение длительностей фронта 10-90 %% %.3f вне грубых ворот\n"
                    "         [%.3f, %.3f] — форма фронта расходится кратно\n",
                    m.rise_ratio, kRiseRatioLo, kRiseRatioHi);
        ++failures;
    }
    if (ok_amp && ok_rise)
    {
        check(true, "§4.3 форма кривой согласуется с эталоном в пределах грубых ворот");
    }
}

// --- Сводный график: сборка серий из выгрузок отдельных прогонов -----------
// Режим --plot-from: расчёта нет, читаются файлы bc2_curve_<тип>.txt, которые
// оставили прогоны с --emit-curve. Отсутствующий тип — не ошибка сборки
// графика, а сообщение: прогон по этому типу просто ещё не делали.
void emit_plot_from_dir(const fs::path& dir, bool full, const fs::path& out_path)
{
    std::printf("\n[график] сборка из выгрузок: %s\n", dir.string().c_str());

    std::vector<PlotSeries> series;
    for (int type = 0; type <= 3; ++type)
    {
        const fs::path p = curve_dump_path(dir, type);
        if (!fs::exists(p))
        {
            std::printf("  info : нет выгрузки для типа %d (%s) — кривая не рисуется\n",
                        type, p.filename().string().c_str());
            continue;
        }
        PlotSeries s;
        std::string err;
        if (!load_series(p, s, err))
        {
            std::printf("  FAIL : выгрузка типа %d не прочитана: %s\n", type, err.c_str());
            ++failures;
            continue;
        }
        // Тип берётся из ИМЕНИ файла: цвет серии обязан отвечать типу, под
        // которым её положили, даже если заголовок внутри расходится.
        s.type = type;
        std::printf("  info : тип %d — %zu точек, %zu шагов%s\n",
                    type, s.curve.t.size(), s.steps, s.aborted ? ", прогон оборван" : "");
        series.push_back(s);
    }

    if (series.empty())
    {
        std::printf("  FAIL : ни одной выгрузки не найдено, график не построен\n");
        ++failures;
        return;
    }

    // Геометрия нужна только для подписи (разрешение, число частиц): строится
    // без расчёта, это заполнение массивов.
    ParticleSystem sys;
    const Bc2Geometry geom = build_bc2_geometry(sys, case_delta_x(full));
    emit_plot(geom, series, full, out_path);
}

void print_usage()
{
    std::printf(
        "Использование: test_e2e_bc2 [--full] [--contact-type=N] [прочие ключи]\n"
        "  (без аргументов)  smoke-режим: delta_x = 8e-4, t_total = 3e-7\n"
        "  --full            полный кейс: delta_x = 2e-4, t_total из конфига\n"
        "  --contact-type=N  contact_force_type прогона, N = 0..3 (по умолчанию 0).\n"
        "                    Один прогон — один тип: сценарии §6 разнесены по\n"
        "                    отдельным тестам, чтобы падение одного типа не\n"
        "                    скрывало результат остальных. Сверка с эталоном\n"
        "                    (§4.3) делается только при N = 0, для остальных\n"
        "                    типов ТЗ §6 требует лишь устойчивости.\n"
        "  --emit-reference  напечатать эталон текущего режима в stdout\n"
        "  --verbose         не глушить вывод конвейера\n"
        "  --keep-workdir    не удалять временный каталог прогона (log.txt, кадры)\n"
        "  --emit-curve=<f>  выгрузить кривую прогона в файл f для сводного\n"
        "                    графика (формат эталона: время, скорость, itype)\n"
        "  --plot-from=<d>   НЕ считать, а собрать сводный график из выгрузок\n"
        "                    bc2_curve_0..3.txt в каталоге d; требует --emit-plot\n"
        "  --emit-plot=<f>   записать SVG «расчёт против эталона» в файл f.\n"
        "                    По умолчанию выключено: тест не оставляет побочных\n"
        "                    эффектов (§7) и должен работать на checkout только\n"
        "                    для чтения. Готовая цель сборки: cmake --build build\n"
        "                    --target bc2_plot\n"
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
    fs::path plot_path;            // пуст — график не строится (см. emit_plot)
    fs::path curve_path;           // пуст — кривая не выгружается
    fs::path plot_from_dir;        // непуст — режим сборки графика без расчёта
    int  contact_type = 0;         // за какой contact_force_type отвечает прогон
    double wall_limit_s = 600.0;   // 10 минут на прогон
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--full") full = true;
        else if (a == "--emit-reference") emit = true;
        else if (a == "--verbose") verbose = true;
        else if (a == "--keep-workdir") keep_workdir = true;
        else if (a.rfind("--contact-type=", 0) == 0)
        {
            const std::string v = a.substr(std::strlen("--contact-type="));
            try { contact_type = std::stoi(v); }
            catch (const std::exception&)
            {
                std::fprintf(stderr, "Неверное значение --contact-type: %s\n", a.c_str());
                return 2;
            }
            if (contact_type < 0 || contact_type > 3)
            {
                std::fprintf(stderr, "--contact-type должен быть 0..3, получено %d\n", contact_type);
                return 2;
            }
        }
        else if (a.rfind("--emit-curve=", 0) == 0)
        {
            curve_path = a.substr(std::strlen("--emit-curve="));
            if (curve_path.empty())
            {
                std::fprintf(stderr, "--emit-curve= требует путь к файлу\n");
                return 2;
            }
        }
        else if (a.rfind("--plot-from=", 0) == 0)
        {
            plot_from_dir = a.substr(std::strlen("--plot-from="));
            if (plot_from_dir.empty())
            {
                std::fprintf(stderr, "--plot-from= требует путь к каталогу\n");
                return 2;
            }
        }
        else if (a.rfind("--emit-plot=", 0) == 0)
        {
            plot_path = a.substr(std::strlen("--emit-plot="));
            if (plot_path.empty())
            {
                std::fprintf(stderr, "--emit-plot= требует путь к файлу\n");
                return 2;
            }
        }
        else if (a == "--emit-plot")
        {
            std::fprintf(stderr, "--emit-plot требует путь: --emit-plot=<файл.svg>\n");
            return 2;
        }
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

    // Режим сборки графика: расчёта нет, поэтому он идёт до всего остального
    // и не трогает ни RunSettings, ни сторожа.
    if (!plot_from_dir.empty())
    {
        if (plot_path.empty())
        {
            std::fprintf(stderr, "--plot-from требует --emit-plot=<файл.svg>\n");
            return 2;
        }
        std::printf("=== e2e BC2: сводный график, режим %s ===\n", full ? "full" : "smoke");
        emit_plot_from_dir(plot_from_dir, full, plot_path);
        std::printf("\n=== e2e BC2: %s (провалов: %d) ===\n",
                    failures == 0 ? "график построен" : "ГРАФИК НЕ ПОСТРОЕН", failures);
        return failures == 0 ? 0 : 1;
    }

    const fs::path cwd = fs::current_path();
    const DirSnapshot before = snapshot_dir(cwd);

    RunSettings rs;
    rs.name = full ? "full" : "smoke";
    rs.contact_force_type = contact_type;
    rs.delta_x = case_delta_x(full);
    // Полный режим берёт t_total из tests/data/bc2_initial_conditions.txt.
    // Smoke укорачивает расчёт: это решение обвязки (быстрый прогон в CI),
    // а не параметр кейса, поэтому оно и живёт в коде, а не в конфиге.
    rs.t_total = full ? case_t_total() : 3e-7;
    rs.silence = !verbose;
    rs.keep_workdir = keep_workdir;
    rs.wall_limit_s = wall_limit_s;

    if (emit)
    {
        // Эталон печатается в stdout — диагностика уходит в stderr, иначе
        // `... --emit-reference > bc2_reference.txt` соберёт мусор.
        std::fprintf(stderr, "test_e2e_bc2: генерация эталона, режим %s...\n", rs.name);
        const RunResult r = run_case(rs);
        if (!r.ok)
        {
            std::fprintf(stderr, "test_e2e_bc2: прогон не состоялся: %s\n", r.error.c_str());
            return 1;
        }
        emit_reference(r, rs, full);
        std::fprintf(stderr, "test_e2e_bc2: напечатано %zu строк данных\n", r.curve.t.size());
        return 0;
    }

    // Сверка с приборным эталоном осмысленна только для базовой постановки:
    // ТЗ §6 требует от контактных сценариев одной устойчивости, эталона для
    // них нет и быть не может.
    const bool is_base = (contact_type == contact::Off);

    std::printf("=== e2e BC2, режим %s, contact_force_type = %d (delta_x = %g, t_total = %g) ===\n",
                rs.name, contact_type, rs.delta_x, rs.t_total);

    // Самопроверка компаратора идёт ДО прогона: она не зависит от расчёта
    // (работает на данных эталона) и должна выполняться даже если прогон упал.
    // Нужна только там, где компаратор применяется, то есть в базовом прогоне.
    if (is_base) test_comparator_selfcheck();

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
        if (is_base)
        {
            check_against_reference(r, full);
        }
        else
        {
            note("§6 сверка с эталоном для контактного сценария не требуется — "
                 "проверена только устойчивость");
        }
    }

    if (r.ok)
    {
        test_known_limitation_b4c_b4c(r.geom, contact_type);
    }

    // Кривая выгружается и при обрыве: частичная история — это и есть то,
    // что сводный график обязан показать, иначе упавший тип просто исчезнет.
    const PlotSeries series{r.curve, contact_type, !r.ok, r.steps};
    if (!curve_path.empty())
    {
        std::printf("\n[выгрузка] %s\n", curve_path.string().c_str());
        if (r.curve.empty())
        {
            std::printf("  FAIL : кривая пуста, выгружать нечего\n");
            ++failures;
        }
        else
        {
            std::string err;
            if (write_curve_dump(series, curve_path, err))
            {
                std::printf("  ok   : выгружено %zu точек (тип %d%s)\n",
                            r.curve.t.size(), contact_type, r.ok ? "" : ", прогон оборван");
            }
            else
            {
                std::printf("  FAIL : %s\n", err.c_str());
                ++failures;
            }
        }
    }

    // Одиночный --emit-plot рисует ровно одну кривую — кривую ЭТОГО прогона.
    // Сводный график по всем типам строится режимом --plot-from.
    if (!plot_path.empty())
    {
        if (r.curve.empty())
        {
            std::printf("\n[график] %s\n", plot_path.string().c_str());
            std::printf("  FAIL : кривая прогона пуста, рисовать нечего\n");
            ++failures;
        }
        else
        {
            emit_plot(r.geom, std::vector<PlotSeries>{series}, full, plot_path);
        }
    }

    check_no_side_effects(cwd, before);

    std::printf("\n=== e2e BC2: %s (провалов: %d) ===\n",
                failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", failures);
    return failures == 0 ? 0 : 1;
}
