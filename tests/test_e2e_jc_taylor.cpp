// test_e2e_jc_taylor.cpp
// E4 (тест Тейлора) и E5 (сходимость по разрешению E4) для материала
// Джонсона-Кука STEEL 4340 (itype 12). План: fiber/jc-implementation-plan.md
// §5.3 E4/E5. Постановка, обоснование стенки и геометрия — tests/jc_taylor.hpp.
//
// РЕЖИМЫ (один исполняемый файл)
//   --case=taylor        (по умолчанию) E4: один прогон на --nr (12 по умолчанию)
//   --case=convergence   E5: прогоны на n_r и 2*n_r (12 и 24), сравнение L_final
//   --case=smoke         короткий прогон (kSmokeSteps шагов, n_r = 6): только
//                        живучесть и здравый смысл, физических ворот нет
//   --nr=N / --dx=МЕТРЫ  число ячеек по радиусу (dx = R0/N); для convergence —
//                        грубый уровень, тонкий вдвое мельче
//   --outdir=ДИР         куда писать кадры и сводки (по умолчанию ./jc_taylor_out)
//   --t-max=СЕК          потолок времени расчёта (по умолчанию 150 мкс)
//   --time-limit=СЕК     бюджет стенных часов на ОДИН прогон (0 — без ограничения)
//   --verbose, --keep-workdir, -h
//
// ОСТАНОВ. Прогон идёт, пока кинетическая энергия половины не упадёт ниже 1 %
// начальной (штатный останов через RunGuard::set_step_limit из наблюдателя
// шага — без исключений), либо до t_max. Если t_max исчерпан раньше, это FAIL
// «стержень не остановился».
//
// ВОРОТА E4 (без внешних данных)
//   1. L_final/L0 в (0.60, 0.95)
//   2. радиус торца не убывает дольше допуска (1 % R0 от бегущего максимума)
//   3. максимум e_pl не дальше 0.05*L0 (1.27 мм) от плоскости удара (z = 0)
//   4. баланс энергии |residual|/E0 < 2 % (сборка SPH_ENERGY_LOG; иначе пропуск
//      с сообщением; на 2026-11 НЕ выполняется, ~20 % — см. g_energy_tol)
//   + здоровье: все величины конечны, ρ > 0, половина не пересекает z = 0,
//     стержень остановился.
//   Сверка с Autodyn: файл tests/data/jc/taylor_4340_autodyn.txt (формат —
//   tests/data/jc/README.md). Есть — ворота 5 % по L_final и R_max; нет —
//   печатается «нет эталона Autodyn» и идут только качественные ворота.
//
// ВОРОТА E5: |L_final(h) − L_final(h/2)| / L_final(h/2) < 3 %.
//
// ВЫХОД (в --outdir): jc_taylor_<тег>_summary.txt (key=value, тот же формат,
// что у файла Autodyn), jc_taylor_<тег>_history.txt (t, KE, L, R_face, e_pl_max),
// jc_taylor_<тег>_frame_<k>.txt (половина z > 0: r z v_r v_z rho p u T e_pl D).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "e2e_harness.hpp"
#include "int_force.hpp"   // urs() — УРС для начального c (шаг по Курану)
#include "jc_taylor.hpp"
#include "run_guard.hpp"
#include "time_integration.hpp"

#ifndef SPH_TEST_DATA_DIR
#define SPH_TEST_DATA_DIR "tests/data"
#endif

static_assert(dim == 2, "тест Тейлора требует сборки с SPH_DIM=2");
static_assert(axis_sim, "тест Тейлора требует осесимметрии (axis_sim)");

namespace fs = std::filesystem;
using e2e::check;
using e2e::note;
using e2e::warn;

namespace
{

// Ворота баланса энергии. По плану E4 — 2 %. ИЗВЕСТНАЯ ПРОБЛЕМА (см. отчёт):
// в осесимметричной постановке с плоской массой rho0*dx^2 невязка лога энергии
// около 20 % и от разрешения не зависит (22/20/18 % при n_r = 6/12/24), поэтому
// ворота на 2 % не проходят. Ключ --energy-warn-only печатает невязку как WARN,
// не роняя тест; ключ --energy-tol задаёт допуск (доля E_mech0).
double g_energy_tol = jct::kEnergyTol;
bool   g_energy_warn_only = false;

// ===========================================================================
// Один прогон
// ===========================================================================
struct Sample
{
    double t = 0.0;
    jct::Measure m;
};

struct RunSettings
{
    std::string tag;
    size_t n_r = jct::kNrE4;
    double t_max = jct::kTMax;
    bool smoke = false;
    bool silence = true;
    bool keep_workdir = false;
    double wall_limit_s = 0.0;
    double epsilon = -1.0;     // XSPH; < 0 — по умолчанию кейса
    fs::path outdir;
};

struct RunResult
{
    bool ok = false;
    std::string error;
    jct::Geometry geom;
    std::vector<Sample> hist;
    jct::Measure fin;          // состояние в момент останова
    double t_end = 0.0;
    double ke0 = 0.0;
    bool stopped = false;      // сработал критерий KE
    size_t steps = 0;
    double wall_s = 0.0;
    double dt_first = 0.0, dt_last = 0.0;
    double z_min_all = 1e30;   // наименьший z половины за прогон
    bool finite = true;
    bool rho_ok = true;
    double R_max = 0.0;        // максимум радиуса торца за прогон
    double R_running_dip = 0.0;// наибольшее убывание R_face от бегущего максимума
    // энергобаланс
    bool have_energy = false;
    double energy_e0 = 0.0, energy_max_res = 0.0;
    size_t energy_rows = 0;
    RunGuard::Reason guard_reason = RunGuard::Reason::None;
    std::string guard_message;
};

void write_frame(const fs::path& path, const ParticleSystem& sys, const jct::Geometry& g, double t)
{
    std::ofstream out(path);
    out << "# t = " << t << " c; половина z > 0; колонки: r z v_r v_z rho p u T e_pl D\n";
    out.precision(9);
    for (size_t k = 0; k < g.nhalf; ++k)
    {
        out << sys.x[0][k] << ' ' << sys.x[1][k] << ' ' << sys.vx[0][k] << ' ' << sys.vx[1][k] << ' '
            << sys.rho[k] << ' ' << sys.p[k] << ' ' << sys.u[k] << ' ' << sys.T[k] << ' '
            << sys.e_pl[k] << ' ' << sys.D[k] << '\n';
    }
}

#ifdef SPH_ENERGY_LOG
// Чтение energy_log.csv: E_mech0 (первая строка) и max|residual|.
void read_energy_csv(const fs::path& p, RunResult& r)
{
    std::ifstream in(p);
    if (!in) return;
    std::string header;
    std::getline(in, header);
    std::vector<std::string> names;
    {
        std::stringstream ss(header);
        std::string c;
        while (std::getline(ss, c, ',')) names.push_back(c);
    }
    const auto col = [&](const char* n) -> int {
        for (size_t i = 0; i < names.size(); ++i)
            if (names[i] == n) return static_cast<int>(i);
        return -1;
    };
    const int c_e = col("E_mech"), c_res = col("residual");
    if (c_e < 0 || c_res < 0) return;
    std::string line;
    bool first = true;
    while (std::getline(in, line))
    {
        std::vector<double> v;
        std::stringstream ss(line);
        std::string c;
        while (std::getline(ss, c, ',')) v.push_back(std::stod(c));
        if (v.size() <= static_cast<size_t>(std::max(c_e, c_res))) continue;
        if (first) { r.energy_e0 = v[static_cast<size_t>(c_e)]; first = false; }
        r.energy_max_res = std::max(r.energy_max_res, std::fabs(v[static_cast<size_t>(c_res)]));
        ++r.energy_rows;
    }
    r.have_energy = r.energy_rows > 0 && r.energy_e0 > 0.0;
}
#endif

RunResult run_case(const RunSettings& rs)
{
    RunResult res;
    std::error_code ec;
    fs::create_directories(rs.outdir, ec);

    e2e::ScopedWorkDir wd("sph_e2e_jc_taylor", rs.tag);
    if (rs.keep_workdir) wd.keep();

    SimulationConfig cfg;
    ParticleSystem sys;
    double dt = 0.0;
    try
    {
        jct::apply_config(cfg, rs.t_max);
        if (rs.epsilon >= 0.0) cfg.epsilon = rs.epsilon;
#ifdef SPH_ENERGY_LOG
        cfg.energy_log_every = 100;
        cfg.energy_log_file = "energy_log.csv";
#endif
        res.geom = jct::build_geometry(sys, rs.n_r);
        const size_t ntotal = res.geom.ntotal;
        cfg.moni_particle = res.geom.index(0, res.geom.n_z - 1);   // тыльный торец на оси
        cfg.maxn = static_cast<size_t>(static_cast<double>(ntotal) * 1.4);
        cfg.max_interaction = 100 * cfg.maxn;
        if (sys.size() < cfg.maxn)
        {
            const size_t prev = sys.size();
            sys.resize(cfg.maxn);
            for (size_t i = prev; i < sys.size(); ++i) sys.inzializ_index(i);
        }
        for (size_t i = 0; i < ntotal; ++i) urs(cfg.mats, sys, i);
        dt = compute_cfl_timestep(cfg, sys, ntotal);
        res.dt_first = dt;
    }
    catch (const std::exception& e)
    {
        res.error = std::string("подготовка прогона не удалась: ") + e.what();
        return res;
    }

    const jct::Geometry& g = res.geom;
    res.ke0 = jct::measure_half(sys, g).ke;
    const double ke_stop = jct::kStopKeFrac * res.ke0;

    RunGuard guard(rs.wall_limit_s, 2.0e6, 200);
    if (rs.smoke) guard.set_step_limit(jct::kSmokeSteps);

    double next_frame_t = 0.0;
    const double frame_dt = 10.0e-6;
    size_t frame_no = 0;
    double R_run_max = 0.0;
    size_t obs_calls = 0;
    bool stop_requested = false;

    // Наблюдатель вызывается ПОСЛЕ шага (RunGuard::should_stop) с новым time.
    // Состояние частиц только читается. Останов — штатный StepLimit.
    guard.set_step_observer([&](double dt_step, double t) {
        res.dt_last = dt_step;
        if (obs_calls++ % jct::kProbeEvery != 0) return;
        const jct::Measure m = jct::measure_half(sys, g);
        res.hist.push_back({t, m});
        res.z_min_all = std::min(res.z_min_all, m.z_min);
        res.R_max = std::max(res.R_max, m.R_face);
        R_run_max = std::max(R_run_max, m.R_face);
        res.R_running_dip = std::max(res.R_running_dip, R_run_max - m.R_face);
        if (!std::isfinite(m.ke) || !std::isfinite(m.L) || !std::isfinite(m.R_face)) res.finite = false;
        if (!rs.smoke && t >= next_frame_t)
        {
            std::ostringstream fn;
            fn << "jc_taylor_" << rs.tag << "_frame_" << frame_no++ << ".txt";
            write_frame(rs.outdir / fn.str(), sys, g, t);
            next_frame_t += frame_dt;
        }
        if (!rs.smoke && !stop_requested && m.ke < ke_stop)
        {
            stop_requested = true;
            guard.set_step_limit(guard.steps() + 1);   // сработает в этом же should_stop
        }
    });

    try
    {
        e2e::ScopedSilence silence(rs.silence);
        time_integration(cfg, sys, g.ntotal, dt, res.t_end, g.dx, &guard);
    }
    catch (const std::exception& e)
    {
        res.error = std::string("исключение в конвейере: ") + e.what();
        res.steps = guard.steps();
        res.wall_s = guard.elapsed_s();
        return res;
    }
    res.steps = guard.steps();
    res.wall_s = guard.elapsed_s();
    res.guard_reason = guard.reason();
    if (guard.failed()) res.guard_message = guard.message();
    res.stopped = stop_requested;

    // Итоговое состояние.
    res.fin = jct::measure_half(sys, g);
    res.z_min_all = std::min(res.z_min_all, res.fin.z_min);
    res.R_max = std::max(res.R_max, res.fin.R_face);
    for (size_t k = 0; k < g.nhalf; ++k)
    {
        if (!(sys.rho[k] > 0.0)) res.rho_ok = false;
        if (!std::isfinite(sys.x[0][k]) || !std::isfinite(sys.x[1][k]) || !std::isfinite(sys.e_pl[k]) ||
            !std::isfinite(sys.vx[1][k]) || !std::isfinite(sys.p[k]))
            res.finite = false;
    }
    if (!rs.smoke)
    {
        std::ostringstream fn;
        fn << "jc_taylor_" << rs.tag << "_frame_final.txt";
        write_frame(rs.outdir / fn.str(), sys, g, res.t_end);
    }
#ifdef SPH_ENERGY_LOG
    read_energy_csv(wd.dir() / "energy_log.csv", res);
#endif
    res.ok = true;
    return res;
}

// ===========================================================================
// Сводка и файлы
// ===========================================================================
void write_summary(const RunSettings& rs, const RunResult& r)
{
    const fs::path hp = rs.outdir / ("jc_taylor_" + rs.tag + "_history.txt");
    {
        std::ofstream h(hp);
        h << "# t[с] KE[Дж] L[м] R_face[м] e_pl_max z_min[м] p_min_face[Па] T_max[К]\n";
        h.precision(9);
        for (const auto& s : r.hist)
            h << s.t << ' ' << s.m.ke << ' ' << s.m.L << ' ' << s.m.R_face << ' ' << s.m.e_pl_max << ' '
              << s.m.z_min << ' ' << s.m.p_min_face << ' ' << s.m.T_max << '\n';
    }
    const fs::path sp = rs.outdir / ("jc_taylor_" + rs.tag + "_summary.txt");
    std::ofstream s(sp);
    s.precision(9);
    s << "# Тест Тейлора STEEL 4340, прогон " << rs.tag << " (формат key=value, единицы СИ)\n";
    s << "L_final=" << r.fin.L << "\n";
    s << "R_max=" << r.R_max << "\n";
    s << "e_pl_axis=" << r.fin.e_pl_axis << "\n";
    s << "e_pl_max=" << r.fin.e_pl_max << "\n";
    s << "e_pl_max_r=" << r.fin.e_pl_max_r << "\n";
    s << "e_pl_max_z=" << r.fin.e_pl_max_z << "\n";
    s << "T_max=" << r.fin.T_max << "\n";
    s << "t_end=" << r.t_end << "\n";
    s << "dx=" << r.geom.dx << "\n";
    s << "n_r=" << r.geom.n_r << "\n";
    s << "n_z=" << r.geom.n_z << "\n";
    s << "particles=" << r.geom.ntotal << "\n";
    s << "steps=" << r.steps << "\n";
    s << "wall_s=" << r.wall_s << "\n";
    if (r.have_energy) s << "energy_residual_rel=" << r.energy_max_res / r.energy_e0 << "\n";
    note("сводка: " + sp.string());
}

void report(const RunSettings& rs, const RunResult& r)
{
    const jct::Geometry& g = r.geom;
    std::printf("\n[%s] n_r = %zu (dx = %.5f мм), n_z = %zu на половину, частиц %zu (обе половины)\n",
                rs.tag.c_str(), g.n_r, g.dx * 1e3, g.n_z, g.ntotal);
    std::printf("  info : шагов %zu, t_end = %.3f мкс, dt: первый %.3e, последний %.3e с, стенные часы %.1f с "
                "(%.2f мс/шаг, %.3g мкс/частицу-шаг)\n",
                r.steps, r.t_end * 1e6, r.dt_first, r.dt_last, r.wall_s,
                r.steps ? 1e3 * r.wall_s / static_cast<double>(r.steps) : 0.0,
                r.steps ? 1e6 * r.wall_s / (static_cast<double>(r.steps) * static_cast<double>(g.ntotal)) : 0.0);
    std::printf("  info : L_final = %.4f мм (L_final/L0 = %.4f), R_max = %.4f мм (R_max/R0 = %.4f)\n",
                r.fin.L * 1e3, r.fin.L / jct::kL0, r.R_max * 1e3, r.R_max / jct::kR0);
    std::printf("  info : e_pl: ось у торца %.4f, максимум %.4f при r = %.3f мм, z = %.3f мм; T_max = %.0f К\n",
                r.fin.e_pl_axis, r.fin.e_pl_max, r.fin.e_pl_max_r * 1e3, r.fin.e_pl_max_z * 1e3, r.fin.T_max);
    std::printf("  info : KE_конец/KE_0 = %.4f; минимум давления в приторцевом слое %.3g Па; "
                "наименьший z половины %.3e м; наибольшее убывание R_face %.3e м\n",
                r.ke0 > 0 ? r.fin.ke / r.ke0 : 0.0, r.fin.p_min_face, r.z_min_all, r.R_running_dip);
}

void check_health(const RunSettings& rs, const RunResult& r)
{
    std::printf("\n[здоровье %s]\n", rs.tag.c_str());
    check(r.ok, r.ok ? "прогон состоялся" : "прогон не состоялся: " + r.error);
    if (!r.ok) return;
    if (r.guard_reason == RunGuard::Reason::WallClock || r.guard_reason == RunGuard::Reason::Stalled)
        check(false, "сторож: " + r.guard_message);
    check(r.finite, "все наблюдаемые величины конечны");
    check(r.rho_ok, "плотность всех частиц положительна");
    check(r.z_min_all > 0.0, "половина z > 0 не пересекает плоскость удара (min z > 0)");
    if (!rs.smoke)
    {
        std::ostringstream os;
        os << "стержень остановился до t_max: KE < " << 100.0 * jct::kStopKeFrac << " % начальной, t = "
           << r.t_end * 1e6 << " мкс (t_max = " << rs.t_max * 1e6 << " мкс)";
        check(r.stopped, os.str());
    }
}

// ---------------------------------------------------------------------------
// Эталон Autodyn
// ---------------------------------------------------------------------------
std::map<std::string, double> read_kv(const fs::path& p)
{
    std::map<std::string, double> kv;
    std::ifstream in(p);
    std::string line;
    while (std::getline(in, line))
    {
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        const auto trim = [](std::string s) {
            const size_t a = s.find_first_not_of(" \t\r");
            const size_t b = s.find_last_not_of(" \t\r");
            return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
        };
        const std::string k = trim(line.substr(0, eq)), v = trim(line.substr(eq + 1));
        if (k.empty() || v.empty()) continue;
        try { kv[k] = std::stod(v); } catch (const std::exception&) {}
    }
    return kv;
}

void check_autodyn(const RunResult& r)
{
    std::printf("\n[сверка с Autodyn]\n");
    const fs::path p = fs::path(SPH_TEST_DATA_DIR) / "jc" / "taylor_4340_autodyn.txt";
    if (!fs::exists(p))
    {
        note("нет эталона Autodyn (" + p.string() + "): проверяются только качественные ворота");
        return;
    }
    const auto kv = read_kv(p);
    const auto rel_gate = [&](const char* key, double sph) {
        const auto it = kv.find(key);
        if (it == kv.end() || !(it->second > 0.0))
        {
            warn(std::string("в эталоне нет ключа ") + key + ": сверка пропущена");
            return;
        }
        const double rel = std::fabs(sph - it->second) / it->second;
        std::ostringstream os;
        os << key << ": SPH " << sph << " м, Autodyn " << it->second << " м, отклонение " << 100.0 * rel
           << " % (ворота " << 100.0 * jct::kAutodynTol << " %)";
        check(rel < jct::kAutodynTol, os.str());
    };
    rel_gate("L_final", r.fin.L);
    rel_gate("R_max", r.R_max);
    const auto it = kv.find("e_pl_axis");
    if (it != kv.end())
    {
        std::ostringstream os;
        os << "e_pl на оси у торца: SPH " << r.fin.e_pl_axis << ", Autodyn " << it->second
           << " (информационно, без ворот: зависит от того, какая ячейка/частица считается «у торца»)";
        note(os.str());
    }
}

void check_e4(const RunSettings& rs, const RunResult& r)
{
    check_health(rs, r);
    if (!r.ok) return;
    report(rs, r);
    if (rs.smoke) return;

    std::printf("\n[ворота E4]\n");
    {
        const double ratio = r.fin.L / jct::kL0;
        std::ostringstream os;
        os << "L_final/L0 = " << ratio << " в (" << jct::kLfLo << ", " << jct::kLfHi << ")";
        check(ratio > jct::kLfLo && ratio < jct::kLfHi, os.str());
    }
    {
        // Допуск на шум/упругую отдачу: радиус торца может убыть не более чем на
        // 1 % R0 от бегущего максимума (после остановки упругая часть радиальной
        // деформации ~σ/E*R ~ 0.5 % R0 частично возвращается).
        std::ostringstream os;
        os << "радиус торца растёт монотонно: наибольшее убывание от бегущего максимума "
           << r.R_running_dip * 1e3 << " мм <= допуск " << jct::kRadiusTol * jct::kR0 * 1e3 << " мм ("
           << 100.0 * jct::kRadiusTol << " % R0)";
        check(r.R_running_dip <= jct::kRadiusTol * jct::kR0, os.str());
        check(r.R_max > 1.05 * jct::kR0, "торец действительно расширился: R_max > 1.05 R0");
    }
    {
        std::ostringstream os;
        const double zlim = jct::kEplFaceL * jct::kL0;
        os << "максимум e_pl (" << r.fin.e_pl_max << ") у ударного торца: z = " << r.fin.e_pl_max_z * 1e3
           << " мм (" << r.fin.e_pl_max_z / (1.3 * r.geom.dx) << " h) <= " << zlim * 1e3 << " мм ("
           << 100.0 * jct::kEplFaceL << " % L0)";
        check(r.fin.e_pl_max_z <= zlim, os.str());
        check(r.fin.e_pl_max > 0.3, "пластическая деформация у торца развита: e_pl_max > 0.3");
    }
    if (r.have_energy)
    {
        const double rel = r.energy_max_res / r.energy_e0;
        std::ostringstream os;
        os << "баланс энергии: max|residual|/E_mech0 = " << 100.0 * rel << " % (< " << 100.0 * g_energy_tol
           << " %), строк лога " << r.energy_rows;
        if (g_energy_warn_only)
        {
            if (rel < g_energy_tol) check(true, os.str());
            else warn(os.str() + " — НЕ УЛОЖИЛСЯ, ворота ослаблены ключом --energy-warn-only");
        }
        else
        {
            check(rel < g_energy_tol, os.str());
        }
    }
    else
    {
#ifdef SPH_ENERGY_LOG
        check(false, "energy_log.csv не прочитан или пуст");
#else
        note("баланс энергии пропущен: сборка без SPH_ENERGY_LOG (цель test_e2e_jc_taylor_elog его проверяет)");
#endif
    }
    check_autodyn(r);
}

void check_smoke(const RunSettings& rs, const RunResult& r)
{
    check_health(rs, r);
    if (!r.ok) return;
    report(rs, r);
    std::printf("\n[smoke: здравый смысл]\n");
    check(r.steps >= jct::kSmokeSteps / 2, "выполнено заметное число шагов");
    check(r.fin.ke < r.ke0, "кинетическая энергия половины убывает после удара");
    check(r.fin.ke > 0.0, "кинетическая энергия положительна");
    check(r.fin.L <= jct::kL0 * (1.0 + 1e-9), "длина половины не растёт");
    check(r.hist.size() >= 2, "наблюдатель снял несколько отсчётов");
}

} // namespace

// ===========================================================================
void print_usage()
{
    std::printf(
        "Использование: test_e2e_jc_taylor [--case=taylor|convergence|smoke] [ключи]\n"
        "  --case=taylor       E4, один прогон (по умолчанию)\n"
        "  --case=convergence  E5, два прогона: n_r и 2*n_r\n"
        "  --case=smoke        короткий прогон на %zu шагов\n"
        "  --nr=N | --dx=М     ячеек по радиусу (taylor: %zu; convergence: 12 и 24; smoke: 6)\n"
        "                      для convergence — грубый уровень, тонкий вдвое мельче\n"
        "  --outdir=ДИР        кадры и сводки (по умолчанию ./jc_taylor_out)\n"
        "  --t-max=СЕК         потолок времени (по умолчанию %.0f мкс)\n"
        "  --time-limit=СЕК    бюджет стенных часов на один прогон (0 — без ограничения)\n"
        "  --xsph=E            константа XSPH epsilon (диагностика; по умолчанию 0.25)\n"
        "  --energy-tol=X      допуск баланса энергии, доля E_mech0 (по умолчанию 0.02)\n"
        "  --energy-warn-only  невязка энергии выше допуска — WARN, а не FAIL\n"
        "  --verbose --keep-workdir -h\n",
        jct::kSmokeSteps, jct::kNrE4, jct::kTMax * 1e6);
}

int main(int argc, char** argv)
{
    std::string mode = "taylor";
    size_t n_r = 0;
    double t_max = jct::kTMax, wall_limit = 0.0, epsilon = -1.0;
    bool verbose = false, keep = false;
    fs::path outdir = fs::current_path() / "jc_taylor_out";

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        const auto val = [&](const char* key) { return a.substr(std::strlen(key)); };
        try
        {
            if (a.rfind("--case=", 0) == 0) mode = val("--case=");
            else if (a.rfind("--nr=", 0) == 0) n_r = static_cast<size_t>(std::stoul(val("--nr=")));
            else if (a.rfind("--dx=", 0) == 0)
            {
                const double dx = std::stod(val("--dx="));
                if (!(dx > 0.0)) throw std::invalid_argument("dx");
                n_r = jct::nr_from_dx(dx);
            }
            else if (a.rfind("--outdir=", 0) == 0) outdir = fs::absolute(val("--outdir="));
            else if (a.rfind("--t-max=", 0) == 0) t_max = std::stod(val("--t-max="));
            else if (a.rfind("--time-limit=", 0) == 0) wall_limit = std::stod(val("--time-limit="));
            else if (a.rfind("--xsph=", 0) == 0) epsilon = std::stod(val("--xsph="));
            else if (a.rfind("--energy-tol=", 0) == 0) g_energy_tol = std::stod(val("--energy-tol="));
            else if (a == "--energy-warn-only") g_energy_warn_only = true;
            else if (a == "--verbose") verbose = true;
            else if (a == "--keep-workdir") keep = true;
            else if (a == "-h" || a == "--help") { print_usage(); return 0; }
            else { std::fprintf(stderr, "Неизвестный ключ: %s\n", a.c_str()); print_usage(); return 2; }
        }
        catch (const std::exception&)
        {
            std::fprintf(stderr, "Неверное значение: %s\n", a.c_str());
            return 2;
        }
    }
    if (mode != "taylor" && mode != "convergence" && mode != "smoke")
    {
        std::fprintf(stderr, "--case должен быть taylor, convergence или smoke\n");
        return 2;
    }
    if (n_r == 0) n_r = (mode == "smoke") ? jct::kNrSmoke : (mode == "convergence") ? jct::kNrConv : jct::kNrE4;

    const fs::path cwd = fs::current_path();
    const e2e::DirSnapshot before = e2e::snapshot_dir(cwd);

    std::printf("E2E JC Тейлор: режим %s, SPH_DIM=%d\n", mode.c_str(), static_cast<int>(dim));
#ifdef SPH_ENERGY_LOG
    note("сборка с SPH_ENERGY_LOG: баланс энергии проверяется");
#else
    note("сборка без SPH_ENERGY_LOG: ворота баланса энергии пропускаются");
#endif
    {
        bool ok = true;
        const jct::Geometry g = jct::plan_geometry(n_r, &ok);
        if (!ok)
            warn("L0/dx не целое: длина половины " + std::to_string(g.L0() * 1e3) +
                 " мм вместо 25.4 мм; берите n_r кратное 3");
    }

    RunSettings rs;
    rs.t_max = t_max;
    rs.silence = !verbose;
    rs.keep_workdir = keep;
    rs.wall_limit_s = wall_limit;
    rs.epsilon = epsilon;
    rs.outdir = outdir;

    if (mode == "smoke")
    {
        rs.tag = "smoke";
        rs.n_r = n_r;
        rs.smoke = true;
        const RunResult r = run_case(rs);
        check_smoke(rs, r);
        if (r.ok) write_summary(rs, r);
    }
    else if (mode == "taylor")
    {
        rs.tag = "nr" + std::to_string(n_r);
        rs.n_r = n_r;
        const RunResult r = run_case(rs);
        check_e4(rs, r);
        if (r.ok) write_summary(rs, r);
    }
    else
    {
        // E5: h и h/2.
        rs.tag = "nr" + std::to_string(n_r);
        rs.n_r = n_r;
        const RunResult rc = run_case(rs);
        check_health(rs, rc);
        if (rc.ok) { report(rs, rc); write_summary(rs, rc); }

        RunSettings rf = rs;
        rf.n_r = 2 * n_r;
        rf.tag = "nr" + std::to_string(2 * n_r);
        const RunResult rfine = run_case(rf);
        check_health(rf, rfine);
        if (rfine.ok) { report(rf, rfine); write_summary(rf, rfine); }

        std::printf("\n[ворота E5: сходимость по разрешению]\n");
        if (rc.ok && rfine.ok)
        {
            const double rel = std::fabs(rc.fin.L - rfine.fin.L) / rfine.fin.L;
            std::ostringstream os;
            os << "L_final(h = " << rc.geom.dx * 1e3 << " мм) = " << rc.fin.L * 1e3 << " мм, L_final(h/2) = "
               << rfine.fin.L * 1e3 << " мм, отличие " << 100.0 * rel << " % (< " << 100.0 * jct::kConvTol << " %)";
            check(rel < jct::kConvTol, os.str());
            // справочно: R_max и e_pl
            std::ostringstream o2;
            o2 << "R_max: " << rc.R_max * 1e3 << " -> " << rfine.R_max * 1e3 << " мм; e_pl_max: "
               << rc.fin.e_pl_max << " -> " << rfine.fin.e_pl_max;
            note(o2.str());
        }
        else
        {
            check(false, "для сравнения нужны оба прогона");
        }
    }

    e2e::check_no_side_effects(cwd, before);
    std::printf("\n%s (провалов: %d)\n", e2e::failures == 0 ? "ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ" : "ЕСТЬ ПРОВАЛЫ", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
