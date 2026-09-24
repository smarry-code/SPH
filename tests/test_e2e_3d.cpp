// test_e2e_3d.cpp
// Сквозной тест конвейера в полном 3D на уменьшенном кейсе BC2.
//
// ЗАЧЕМ. Сквозные тесты BC2/BC6 и регрессия считают только dim = 2 с
// осевой симметрией. Ветка dim = 3 проходит через тот же конвейер другими
// путями: шесть компонент девиатора, три поворота, трёхмерная сетка поиска
// соседей, контакт без зеркальных частиц. Нити (fiber/PLAN.md) живут только в
// 3D, поэтому Solid должен быть проверен именно здесь, до того как Ф4
// перенесёт его на новую схему RHS.
//
// КЕЙС — tests/bc2_mini3d_geometry.hpp: материалы, толщины, плотности и
// скорости BC2, радиусы уменьшены до R_керамики = 4 мм. С приборным эталоном
// его не сверить (боковая разгрузка доходит до оси почти сразу), поэтому
// проверяется то, что следует из законов сохранения и симметрии:
//
//   §1 геометрия: число частиц, отображения симметрии — инволюции;
//   §2 здоровье: сторож, коллапс шага, NaN/Inf, ρ > 0, m > 0 на всех кадрах;
//   §3 сохранение: масса — точно; импульс — поперечный остаётся нулевым,
//      осевой не меняется (парные силы антисимметричны, внешних сил нет);
//      полная энергия m(u + v²/2) не растёт. Сохраняться она не обязана:
//      XSPH-поправка av прибавляется к СКОРОСТИ (time_integration), то есть
//      сглаживает поле скоростей без нагрева, и на этом кейсе за 1 мкс уносит
//      около четверти энергии. Поэтому отдельный прогон --no-xsph (epsilon = 0)
//      проверяет, что без неё энергия сохраняется;
//   §4 симметрия решётки: x → −x, y → −y, x ↔ y. Начальные данные ей
//      симметричны, значит и решение обязано быть симметричным — до уровня
//      округления. Перепутанные в 3D-ветке компоненты (xz вместо yz, знак при
//      повороте) эту симметрию ломают;
//   §5 монитор: тишина до прихода волны, фронт приходит в своё окно.
//
// ОДИН ПРОГОН — ОДИН ТИП КОНТАКТА, как в BC2/BC6.
//
// РЕЖИМЫ:
//   без аргументов — smoke: до kSmokeTTotal = 1 мкс (фронт успевает пройти
//                    монитор при любом типе контакта; останов по времени, а не
//                    по числу шагов: тип 2 режет шаг, и за 100 шагов доходит
//                    лишь до 0.53 мкс);
//   --full         — до t_total из конфига BC2 (3.5 мкс);
//   --no-xsph      — epsilon = 0: проверка сохранения полной энергии.
//
// ВАЖНО ПРО РАЗМЕРНОСТЬ: цель собирается с -DSPH_DIM=3 и компилирует свои
// копии src/*.cpp (решение E1, ADR-002).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "input.hpp"
#include "int_force.hpp"
#include "time_integration.hpp"
#include "run_guard.hpp"

#include "e2e_harness.hpp"
#include "e2e_curve.hpp"
#include "e2e_frames.hpp"
#include "bc2_mini3d_geometry.hpp"

static_assert(dim == 3, "test_e2e_3d требует сборки с -DSPH_DIM=3");

#ifndef SPH_TEST_DATA_DIR
#define SPH_TEST_DATA_DIR "tests/data"
#endif

namespace fs = std::filesystem;

using e2e::check;
using e2e::note;
using e2e::Curve;

namespace
{

// ===========================================================================
// 1. Параметры
// ===========================================================================

constexpr double kSmokeTTotal = 1.0e-6;

// Контрольное число геометрии при Δx = 0.8 мм: замок от незаметного
// изменения постановки (960 + 400 + 880 + 1408).
constexpr size_t kSmokeNtotal = 3648;

fs::path config_path() { return fs::path(SPH_TEST_DATA_DIR) / "bc2_initial_conditions.txt"; }

struct RunSettings
{
    const char* name = "smoke";
    bool   full = false;
    int    contact_force_type = 0;
    bool   no_xsph = false;
    bool   silence = true;
    bool   keep_workdir = false;
    double wall_limit_s = 600.0;
};

struct Totals
{
    double mass = 0.0;
    std::array<double, 3> momentum{};
    double kinetic = 0.0;
    double internal = 0.0;
    double energy() const { return kinetic + internal; }
};

Totals totals(const ParticleSystem& sys, size_t n)
{
    Totals t;
    for (size_t i = 0; i < n; ++i)
    {
        double v2 = 0.0;
        for (size_t d = 0; d < 3; ++d)
        {
            t.momentum[d] += sys.mass[i] * sys.vx[d][i];
            v2 += sys.vx[d][i] * sys.vx[d][i];
        }
        t.mass += sys.mass[i];
        t.kinetic += 0.5 * sys.mass[i] * v2;
        t.internal += sys.mass[i] * sys.u[i];
    }
    return t;
}

// ===========================================================================
// 2. Симметрия
// ===========================================================================
// Поля частицы и их поведение при симметриях решётки. Знак компоненты при
// отражении x → −x меняют ровно те, у которых индекс x встречается нечётное
// число раз: x0, v0, sxy, sxz. Перестановка x ↔ y меняет местами пары
// (x0, x1), (v0, v1), (sxx, syy), (sxz, syz).

enum F : size_t
{
    fX0, fX1, fX2, fV0, fV1, fV2, fRho, fP, fU, fC, fH, fEpl, fD,
    fSxx, fSyy, fSzz, fSxy, fSyz, fSxz, fCount
};

constexpr std::array<const char*, fCount> kName = {
    "x0", "x1", "x2", "v0", "v1", "v2", "rho", "p", "u", "c", "hsml", "e_pl", "D",
    "sxx", "syy", "szz", "sxy", "syz", "sxz"};

// Пол масштаба поля — как в регрессии (tests/regression_snapshot.hpp).
constexpr std::array<double, fCount> kFloor = {
    1e-3, 1e-3, 1e-3, 1.0, 1.0, 1.0, 1.0, 1e6, 1.0, 1.0, 1e-6, 1e-2, 1.0,
    1e6, 1e6, 1e6, 1e6, 1e6, 1e6};

using Fields = std::array<std::vector<double>, fCount>;

Fields take_fields(const ParticleSystem& sys, size_t n)
{
    Fields f;
    auto take = [&](F k, const std::vector<double>& src) {
        f[k].assign(src.begin(), src.begin() + static_cast<std::ptrdiff_t>(n));
    };
    take(fX0, sys.x[0]); take(fX1, sys.x[1]); take(fX2, sys.x[2]);
    take(fV0, sys.vx[0]); take(fV1, sys.vx[1]); take(fV2, sys.vx[2]);
    take(fRho, sys.rho); take(fP, sys.p); take(fU, sys.u); take(fC, sys.c);
    take(fH, sys.hsml); take(fEpl, sys.e_pl); take(fD, sys.D);
    take(fSxx, sys.sxx); take(fSyy, sys.syy); take(fSzz, sys.szz);
    take(fSxy, sys.sxy); take(fSyz, sys.syz); take(fSxz, sys.sxz);
    return f;
}

struct Symmetry
{
    const char* name;
    const std::vector<size_t>* image;
    std::array<size_t, fCount> source{};   // поле образа, с которым сравнивать
    std::array<double, fCount> sign{};
};

Symmetry make_symmetry(const char* name, const std::vector<size_t>& image, int kind)
{
    Symmetry s{name, &image, {}, {}};
    for (size_t k = 0; k < fCount; ++k) { s.source[k] = k; s.sign[k] = 1.0; }
    if (kind == 0)        // x → −x
    {
        for (F k : {fX0, fV0, fSxy, fSxz}) s.sign[k] = -1.0;
    }
    else if (kind == 1)   // y → −y
    {
        for (F k : {fX1, fV1, fSxy, fSyz}) s.sign[k] = -1.0;
    }
    else                  // x ↔ y
    {
        for (auto [a, b] : {std::pair{fX0, fX1}, std::pair{fV0, fV1}, std::pair{fSxx, fSyy},
                            std::pair{fSxz, fSyz}})
        {
            s.source[a] = b;
            s.source[b] = a;
        }
    }
    return s;
}

struct SymmetryResult
{
    double worst = 0.0;      // max по полям и частицам |f − T f| / scale
    size_t worst_field = 0;
    size_t worst_particle = 0;
};

SymmetryResult measure_symmetry(const Fields& f, const Symmetry& s)
{
    SymmetryResult r;
    const size_t n = f[0].size();
    for (size_t k = 0; k < fCount; ++k)
    {
        double max_abs = 0.0;
        for (double v : f[k]) max_abs = std::max(max_abs, std::fabs(v));
        const double scale = std::max(max_abs, kFloor[k]);
        for (size_t i = 0; i < n; ++i)
        {
            const size_t j = (*s.image)[i];
            const double d = std::fabs(f[k][i] - s.sign[k] * f[s.source[k]][j]) / scale;
            if (!(d <= r.worst))
            {
                r.worst = d;
                r.worst_field = k;
                r.worst_particle = i;
            }
        }
    }
    return r;
}

// ===========================================================================
// 3. Прогон
// ===========================================================================

struct RunResult
{
    bool ok = false;
    std::string error;
    bc2m3::Geometry geom;
    Totals before, after;
    Fields fields;
    Curve curve;
    e2e::FrameScan scan;
    double t_reached = 0.0;
    double t_total = 0.0;
    size_t steps = 0;
    RunGuard::Reason guard_reason = RunGuard::Reason::None;
    std::string guard_message;
    std::string workdir;
};

RunResult run_case(const RunSettings& rs)
{
    RunResult res;
    e2e::ScopedWorkDir wd("sph_e2e_3d", std::string(rs.name) + (rs.no_xsph ? "_noxsph" : "") +
                                            "_c" + std::to_string(rs.contact_force_type));
    if (rs.keep_workdir) wd.keep();
    res.workdir = wd.dir().string();

    SimulationConfig cfg;
    ParticleSystem sys;
    double dt = 0.0;
    try
    {
        cfg.apply_defaults();
        if (!fs::exists(config_path()))
        {
            throw std::runtime_error("не найден конфиг теста: " + config_path().string());
        }
        parse_sph_config(cfg, config_path().string());
        cfg.contact_force_type = rs.contact_force_type;
        if (!rs.full) cfg.t_total = kSmokeTTotal;
        if (rs.no_xsph) cfg.epsilon = 0.0;

        res.geom = bc2m3::build_geometry(sys, bc2m3::kDxSmoke);
        const size_t n = res.geom.ntotal;
        cfg.moni_particle = res.geom.moni;
        cfg.maxn = static_cast<size_t>(static_cast<double>(n) * 1.4);
        cfg.max_interaction = 100 * cfg.maxn;
        const size_t prev = sys.size();
        sys.resize(cfg.maxn);
        for (size_t i = prev; i < sys.size(); ++i) sys.inzializ_index(i);
        for (size_t i = 0; i < n; ++i) urs(cfg.mats, sys, i);
        dt = compute_cfl_timestep(cfg, sys, n);
    }
    catch (const std::exception& e)
    {
        res.error = std::string("подготовка прогона не удалась: ") + e.what();
        return res;
    }
    const size_t n = res.geom.ntotal;
    res.t_total = cfg.t_total;
    res.before = totals(sys, n);

    RunGuard guard(rs.wall_limit_s, 2.0e6, 200);

    double time = 0.0;
    try
    {
        e2e::ScopedSilence silence(rs.silence);
        time_integration(cfg, sys, n, dt, time, bc2m3::kDxSmoke, &guard);
    }
    catch (const std::exception& e)
    {
        res.error = std::string("исключение в конвейере: ") + e.what() + " (на шаге " +
                    std::to_string(guard.steps()) + ")";
        return res;
    }
    res.t_reached = time;
    res.guard_reason = guard.reason();
    if (guard.failed()) res.guard_message = guard.message();
    res.after = totals(sys, n);
    res.fields = take_fields(sys, n);

    std::string err;
    res.curve = e2e::read_curve(wd.dir() / "log.txt", err);
    if (!err.empty())
    {
        res.error = err;
        return res;
    }
    res.steps = res.curve.t.size();

    e2e::FrameProbe probe;
    probe.particles = {res.geom.moni};
    probe.damage_itype = 10;
    res.scan = e2e::scan_frames(wd.dir(), probe);

    res.ok = true;
    return res;
}

// ===========================================================================
// 4. Проверки
// ===========================================================================

void check_geometry(const bc2m3::Geometry& g)
{
    std::printf("\n[§1] Геометрия (уменьшенный BC2, Δx = %g м)\n", g.delta_x);
    const auto layers = bc2::bc2_layers();
    for (size_t L = 0; L < 4; ++L)
    {
        std::printf("  info : слой %zu (%s): R = %.3f мм, %zu столбцов x %zu слоёв = %zu, "
                    "dz = %.4f мм\n",
                    L + 1, layers[L].name, g.layer_radius[L] * 1e3, g.layer_columns[L],
                    g.layer_nz[L], g.layer_count[L], g.layer_dz[L] * 1e3);
    }
    note("ntotal = " + std::to_string(g.ntotal) + ", монитор = " + std::to_string(g.moni));

    bool involutions = true, fixed_free = true;
    for (size_t i = 0; i < g.ntotal; ++i)
    {
        for (const auto* m : {&g.mirror_x, &g.mirror_y})
        {
            involutions = involutions && (*m)[(*m)[i]] == i;
            fixed_free = fixed_free && (*m)[i] != i;   // на решётке со сдвигом ½ нет неподвижных
        }
        involutions = involutions && g.swap_xy[g.swap_xy[i]] == i;
    }
    check(involutions, "отображения симметрии — инволюции (образ образа = сама частица)");
    check(fixed_free, "отражения не имеют неподвижных частиц: ось проходит между столбцами");
    check(g.ntotal == kSmokeNtotal,
          "ntotal = " + std::to_string(kSmokeNtotal) + " (замок постановки)");
}

void check_health(const RunResult& r, const RunSettings& rs)
{
    std::printf("\n[§2] Здоровье расчёта\n");
    {
        std::ostringstream os;
        os << "шагов = " << r.steps << ", t = " << r.t_reached << " (t_total = " << r.t_total
           << ")";
        note(os.str());
    }
    if (r.guard_reason == RunGuard::Reason::WallClock ||
        r.guard_reason == RunGuard::Reason::Stalled)
    {
        check(false, "сторож оборвал расчёт:\n" + r.guard_message);
    }
    check(r.t_reached >= r.t_total, "расчёт дошёл до t_total");
    if (r.curve.t.size() >= 2)
    {
        const double dt_first = r.curve.t.front();
        double dt_min = dt_first;
        for (size_t i = 1; i < r.curve.t.size(); ++i)
            dt_min = std::min(dt_min, r.curve.t[i] - r.curve.t[i - 1]);
        std::ostringstream os;
        os << "шаг не сколлапсировал: dt_first/dt_min = " << dt_first / dt_min << " < 100";
        check(dt_min > 0.0 && dt_first / dt_min < 100.0, os.str());
    }
    const e2e::FrameCheck& fc = r.scan.health;
    note("кадров = " + std::to_string(fc.frames) + ", строк проверено = " + std::to_string(fc.rows));
    check(fc.frames > 0 && fc.rows > 0 && fc.short_rows == 0,
          "кадры записаны и содержат все колонки 3D");
    check(fc.bad_finite == 0, "нет NaN/Inf в кадрах" +
                                  (fc.first_problem.empty() ? "" : " (" + fc.first_problem + ")"));
    check(fc.bad_rho == 0, "rho > 0 на всех кадрах");
    check(fc.bad_mass == 0, "mass > 0 на всех кадрах");
    {
        std::ostringstream os;
        os << "макс D по керамике за прогон = " << r.scan.max_damage << " (частица "
           << r.scan.max_damage_index << ", r = " << r.scan.max_damage_r * 1e3
           << " мм, z = " << r.scan.max_damage_z * 1e3 << " мм)";
        note(os.str());
    }
}

// Пороги §3 и §4. Измерено на smoke (2026-09-25):
//   импульс: поперечный ~1e-17, осевой ~1e-14 от |P_z(0)| — уровень
//     округления; порог 1e-9 оставляет пять порядков запаса и ловит любую
//     несимметричную парную силу;
//   энергия без XSPH: +3.4e-3 за 100 шагов (схема слегка разогревает);
//     порог 1 % — втрое выше;
//   симметрия: отражения ~3e-14; порог 1e-9.
constexpr double kMomentumTol  = 1e-9;   // |ΔP| / |P_z(0)|
constexpr double kEnergyGrowth = 0.01;   // ΔE / E(0) сверху: энергия не создаётся
constexpr double kEnergyTol    = 0.01;   // |ΔE| / E(0) при --no-xsph
constexpr double kSymmetryTol  = 1e-9;   // max |f − T f| / scale

void check_conservation(const RunResult& r, const RunSettings& rs)
{
    std::printf("\n[§3] Законы сохранения\n");
    check(r.after.mass == r.before.mass, "масса сохраняется точно");
    const double pz0 = std::fabs(r.before.momentum[2]);
    const double dpx = std::fabs(r.after.momentum[0]) / pz0;
    const double dpy = std::fabs(r.after.momentum[1]) / pz0;
    const double dpz = std::fabs(r.after.momentum[2] - r.before.momentum[2]) / pz0;
    {
        std::ostringstream os;
        os << "импульс: P(0) = (" << r.before.momentum[0] << ", " << r.before.momentum[1] << ", "
           << r.before.momentum[2] << "), P = (" << r.after.momentum[0] << ", "
           << r.after.momentum[1] << ", " << r.after.momentum[2] << ") кг·м/с";
        note(os.str());
    }
    char buf[160];
    std::snprintf(buf, sizeof(buf), "поперечный импульс нулевой: |P_x|, |P_y| / |P_z(0)| = "
                                    "%.2e, %.2e <= %.0e", dpx, dpy, kMomentumTol);
    check(dpx <= kMomentumTol && dpy <= kMomentumTol, buf);
    std::snprintf(buf, sizeof(buf), "осевой импульс сохраняется: |ΔP_z| / |P_z(0)| = %.2e <= %.0e",
                  dpz, kMomentumTol);
    check(dpz <= kMomentumTol, buf);

    const double e0 = r.before.energy();
    const double de = (r.after.energy() - e0) / e0;
    std::snprintf(buf, sizeof(buf),
                  "энергия: E(0) = %.6e Дж, кинетическая %.6e -> %.6e, внутренняя %.6e -> %.6e",
                  e0, r.before.kinetic, r.after.kinetic, r.before.internal, r.after.internal);
    note(buf);
    if (rs.no_xsph)
    {
        std::snprintf(buf, sizeof(buf),
                      "без XSPH полная энергия сохраняется: |ΔE/E(0)| = %.3e <= %.2f",
                      std::fabs(de), kEnergyTol);
        check(std::fabs(de) <= kEnergyTol, buf);
    }
    else
    {
        std::snprintf(buf, sizeof(buf),
                      "полная энергия не создаётся: ΔE/E(0) = %+.3e <= %+.2f "
                      "(убыль — работа XSPH-сглаживания скорости, см. шапку)",
                      de, kEnergyGrowth);
        check(de <= kEnergyGrowth, buf);
    }
}

void check_symmetry(const RunResult& r)
{
    std::printf("\n[§4] Симметрия решётки\n");
    const std::array<Symmetry, 3> syms = {make_symmetry("x → −x", r.geom.mirror_x, 0),
                                          make_symmetry("y → −y", r.geom.mirror_y, 1),
                                          make_symmetry("x ↔ y", r.geom.swap_xy, 2)};
    for (const Symmetry& s : syms)
    {
        const SymmetryResult m = measure_symmetry(r.fields, s);
        const size_t i = m.worst_particle;
        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "%s: худшее расхождение %.2e (поле %s, частица %zu, x = (%.2f %.2f %.2f) мм) "
                      "<= %.0e",
                      s.name, m.worst, kName[m.worst_field], i, r.fields[fX0][i] * 1e3,
                      r.fields[fX1][i] * 1e3, r.fields[fX2][i] * 1e3, kSymmetryTol);
        check(m.worst <= kSymmetryTol, buf);
    }
}

void check_monitor(const RunResult& r)
{
    std::printf("\n[§5] Монитор: тыльная грань керамики, приосевой столбец\n");
    const Curve& c = r.curve;
    // Фронт разбирается только в окне первой волны (t <= kSmokeTTotal = 1 мкс).
    // Дальше монитор разгоняется всей сборкой: в полном прогоне к 3.5 мкс
    // v_z = 468 м/с и пик приходится на последний шаг, так что «50 % пика»
    // по всей кривой относилось бы уже не к фронту.
    size_t n_win = 0;
    while (n_win < c.t.size() && c.t[n_win] <= kSmokeTTotal) ++n_win;
    double vmax = 0.0, t_first = -1.0, t_peak = 0.0;
    size_t i_peak = 0;
    for (size_t i = 0; i < n_win; ++i)
    {
        if (c.v[i] > vmax) { vmax = c.v[i]; t_peak = c.t[i]; i_peak = i; }
        if (t_first < 0.0 && std::fabs(c.v[i]) > 1.0) t_first = c.t[i];
    }
    double t50 = -1.0;
    for (size_t i = 0; i <= i_peak && i < n_win; ++i)
    {
        if (c.v[i] >= 0.5 * vmax) { t50 = c.t[i]; break; }
    }
    // Пик первой волны должен лежать внутри окна, а не на его краю: иначе
    // фронт ещё не прошёл монитор, и t50 ничего не значит.
    const bool peak_inside = n_win > 0 && i_peak + 1 < n_win;
    // Упругая волна в B4C (c_L ≈ 14 км/с) проходит 9.03 мм покоящейся
    // керамики за ~0.65 мкс. На этой сетке фронт сильно размазан: слои
    // взаимодействуют через носитель ядра (2h = 3.2 мм) уже с t = 0, и
    // уровень 50 % приходит раньше, около 0.51 мкс; первые 1 м/с — около
    // 0.21 мкс. Окна ниже ловят кратную ошибку скорости волны, а не
    // точность разрешения.
    constexpr double kQuietUntil = 0.15e-6;
    constexpr double kT50Lo = 0.35e-6, kT50Hi = 0.75e-6;
    char buf[220];
    std::snprintf(buf, sizeof(buf),
                  "монитор (окно 1 мкс): пик v_z = %.3f м/с при %.4f мкс; первое |v_z| > 1 м/с при %.4f мкс; "
                  "50 %% пика при %.4f мкс",
                  vmax, t_peak * 1e6, t_first * 1e6, t50 * 1e6);
    note(buf);
    std::snprintf(buf, sizeof(buf), "тишина до %.2f мкс: |v_z| <= 1 м/с", kQuietUntil * 1e6);
    check(t_first < 0.0 || t_first >= kQuietUntil, buf);
    check(vmax > 10.0, "волна дошла до монитора: пик v_z > 10 м/с");
    std::snprintf(buf, sizeof(buf), "фронт (50 %% пика) приходит в окне [%.2f; %.2f] мкс",
                  kT50Lo * 1e6, kT50Hi * 1e6);
    check(t50 >= kT50Lo && t50 <= kT50Hi && peak_inside, buf);
    bool itype_ok = !c.itype.empty();
    for (int it : c.itype) itype_ok = itype_ok && it == 10;
    check(itype_ok, "монитор всё время частица B4C (itype = 10)");
}

void print_usage()
{
    std::printf(
        "Использование: test_e2e_3d [--full] [--contact-type=N] [ключи]\n"
        "  (без аргументов)  smoke: до t = %g с\n"
        "  --full            до t_total из tests/data/bc2_initial_conditions.txt\n"
        "  --contact-type=N  contact_force_type, N = 0..3 (по умолчанию 0)\n"
        "  --no-xsph         epsilon = 0: проверка сохранения полной энергии\n"
        "  --verbose         не глушить вывод конвейера\n"
        "  --keep-workdir    не удалять временный каталог прогона\n"
        "  --time-limit=СЕК  бюджет времени на прогон (по умолчанию 600)\n",
        kSmokeTTotal);
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    RunSettings rs;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        try
        {
            if (a == "--full") { rs.full = true; rs.name = "full"; }
            else if (a == "--verbose") rs.silence = false;
            else if (a == "--no-xsph") rs.no_xsph = true;
            else if (a == "--keep-workdir") rs.keep_workdir = true;
            else if (a.rfind("--contact-type=", 0) == 0)
            {
                rs.contact_force_type = std::stoi(a.substr(std::strlen("--contact-type=")));
                if (rs.contact_force_type < 0 || rs.contact_force_type > 3)
                    throw std::invalid_argument("N = 0..3");
            }
            else if (a.rfind("--time-limit=", 0) == 0)
            {
                rs.wall_limit_s = std::max(0.0, std::stod(a.substr(std::strlen("--time-limit="))));
            }
            else if (a == "-h" || a == "--help") { print_usage(); return 0; }
            else throw std::invalid_argument("неизвестный аргумент");
        }
        catch (const std::exception& e)
        {
            std::fprintf(stderr, "Неверный аргумент %s: %s\n", a.c_str(), e.what());
            print_usage();
            return 2;
        }
    }

    const fs::path cwd = fs::current_path();
    const e2e::DirSnapshot before = e2e::snapshot_dir(cwd);

    std::printf("=== e2e 3D (уменьшенный BC2), режим %s%s, contact_force_type = %d ===\n",
                rs.name, rs.no_xsph ? " без XSPH" : "", rs.contact_force_type);

    const RunResult r = run_case(rs);
    if (rs.keep_workdir) note("каталог прогона сохранён: " + r.workdir);
    if (!r.ok)
    {
        std::printf("\n[прогон] FAIL : прогон не состоялся: %s\n", r.error.c_str());
        ++e2e::failures;
    }
    else
    {
        check_geometry(r.geom);
        check_health(r, rs);
        check_conservation(r, rs);
        check_symmetry(r);
        check_monitor(r);
    }

    e2e::check_no_side_effects(cwd, before);

    std::printf("\n=== e2e 3D: %s (провалов: %d) ===\n",
                e2e::failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
