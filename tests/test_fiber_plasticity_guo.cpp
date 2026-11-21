// tests/test_fiber_plasticity_guo.cpp
// Этап 5 дизайна пластичности нити (notes/fiber-model/fiber-plasticity-design.md, раздел 10,
// п. 5), тест Т16: интеграционный прогон сцены Guo 2019 (cases/guo2019/regression/guo2019_scene.hpp)
// с упругой и с упруго-пластичной нитью. ЭТО НЕ КАЛИБРОВКА И НЕ СРАВНЕНИЕ С
// ОПЫТОМ: параметры пластики иллюстративные (заданы ОТНОСИТЕЛЬНО E и eps_break
// материала сцены), экспериментальных данных Guo для пластической модели нет;
// проверяется только качественное «упругая против пластичной».
//
// Пластика задаётся директивой .fib (PLASTIC 0 SY=... H=...) -- прогон заодно
// e2e-проверка разбора директив: из минимального .fib загружаются материалы и
// их неупругие поля переносятся в model.materials сцены (сам .fib сцены пишет
// генератор плетения и пластики не знает).
//
// Режимы запуска:
//   test_fiber_plasticity_guo                        -- smoke для ctest (1T, уровень smoke,
//                                                       V0 = 100 м/с, 20 мкс): жёсткие проверки
//   test_fiber_plasticity_guo --level=scaled --config=1T --v0=... --t-max=... --out=<dir>
//                                                    -- длинный прогон: таблица в <out>/summary.txt,
//                                                       проверки не выполняются (только отчёт)
// Опции: --sy-frac=<доля E*eps_break> (0.2), --h-frac=<доля E> (0.1), --time-limit=<с>,
//   --eps-break-p=<eps_break пластичной>, --eps-break-e=<eps_break упругой> (контроль), --only=elastic|plastic (один прогон).
#include "guo2019_monitor.hpp"
#include "guo2019_scene.hpp"

#include "e2e_harness.hpp"

#include "fiber/fiber_input.hpp"
#include "fiber/fiber_plasticity.hpp"
#include "run_guard.hpp"
#include "time_integration.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using e2e::check;
using e2e::note;
using namespace guo2019;

namespace
{

std::string fmtg(double v, int prec = 6)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", prec, v);
    return buf;
}

// Наблюдатель нити: пиковые деформация/напряжение связей, накопленные энергии
// (с учётом порванных связей: их последнее известное значение остаётся в сумме).
class FiberObserver
{
public:
    FiberObserver(const fiber::Model& m, std::size_t every) : m_(m), every_(every ? every : 1) {}

    void on_step(double /*dt*/, double /*time*/)
    {
        if (++n_ % every_ != 0) return;
        sample();
    }
    void sample()
    {
        const fiber::FiberState& fb = m_.fiber;
        for (const fiber::FiberBond& b : m_.topology.bonds)
        {
            const fiber::FiberMaterial& mat = m_.materials[static_cast<std::size_t>(fb.material_id[b.i])];
            double l2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double dx = fb.x[d][b.j] - fb.x[d][b.i];
                l2 += dx * dx;
            }
            const double eps = std::sqrt(l2) / b.l0 - 1.0;
            const double sig = mat.E * (eps - b.st.eps_p);
            eps_max_ = std::max(eps_max_, eps);
            epsp_max_ = std::max(epsp_max_, b.st.eps_p);
            if (eps > b.st.eps_p) sig_max_ = std::max(sig_max_, sig);
            const double V = interactions::fiber_cross_section_area(mat) * b.l0;
            auto& rec = bonds_[{b.i, b.j}];
            rec.first = b.st.w_p * V;
            rec.second = b.st.w_v * V;
        }
        double wp = 0.0, wv = 0.0;
        for (const auto& kv : bonds_)
        {
            wp += kv.second.first;
            wv += kv.second.second;
        }
        wp_ = wp;
        wv_ = wv;
    }
    double eps_max() const { return eps_max_; }
    double epsp_max() const { return epsp_max_; }
    double sig_max() const { return sig_max_; }
    double wp() const { return wp_; }
    double wv() const { return wv_; }

private:
    const fiber::Model& m_;
    std::size_t every_, n_ = 0;
    double eps_max_ = 0.0, epsp_max_ = 0.0, sig_max_ = 0.0, wp_ = 0.0, wv_ = 0.0;
    std::map<std::pair<std::size_t, std::size_t>, std::pair<double, double>> bonds_;
};

struct RunResult
{
    bool ok = false;
    std::string error;
    std::size_t steps = 0, n_bonds_0 = 0, n_bonds_end = 0;
    double t_end = 0.0, wall = 0.0;
    double v_z_end = 0.0, z_c_end = 0.0, V_R = 0.0;
    double ke_ball_0 = 0.0, ke_ball_end = 0.0;
    double eps_max = 0.0, epsp_max = 0.0, sig_max = 0.0, wp = 0.0, wv = 0.0;
    double F_f_z_max = 0.0, F_s_z_max = 0.0, w_f_max = 0.0;
    Outcome outcome = Outcome::Undetermined;
    double SY = 0.0, H = 0.0, E = 0.0, eps_break = 0.0;
    bool finite = true;
};

// Включает пластику материалов сцены через директиву .fib (e2e парсера).
void apply_plasticity_via_fib(GuoScene& sc, const fs::path& workdir, double sy_frac, double h_frac, double eps_break_plastic,
                              RunResult& r)
{
    const fiber::FiberMaterial& m0 = sc.model.materials.at(0);
    r.E = m0.E;
    r.eps_break = m0.eps_break;
    r.SY = sy_frac * m0.E * m0.eps_break;
    r.H = h_frac * m0.E;
    char buf[512];
    const fs::path p = workdir / "guo2019_plastic_material.fib";
    {
        std::ofstream f(p);
        std::snprintf(buf, sizeof buf, "MATERIAL 0 %.17g %.17g %.17g %.17g %.17g\n", m0.E, m0.nu, m0.d, m0.rho, m0.eps_break);
        f << buf;
        std::snprintf(buf, sizeof buf, "PLASTIC 0 SY=%.17g H=%.17g\n", r.SY, r.H);
        f << buf;
        // Нить с двумя бусинами -- парсер требует хотя бы одну FIBER.
        f << "FIBER 0 2\n0 0 0\n1e-3 0 0\n";
    }
    fiber::Model tmp;
    fiber::load_fiber_file(p.string(), tmp);
    if (tmp.materials.empty() || !tmp.materials.at(0).plastic.enabled())
        throw std::runtime_error("директива PLASTIC не включила пластику материала");
    for (fiber::FiberMaterial& m : sc.model.materials)
    {
        m.plastic = tmp.materials[0].plastic;
        m.rate = tmp.materials[0].rate;
        m.visco = tmp.materials[0].visco;
        m.failure = tmp.materials[0].failure;
        if (eps_break_plastic > 0.0) m.eps_break = eps_break_plastic;
    }
    if (eps_break_plastic > 0.0)
    {
        // l_break хранится в связи (создаётся из eps_break при загрузке .fib): пересчитываем так же, (1+eps_break)*l0.
        for (fiber::FiberBond& b : sc.model.topology.bonds) b.l_break = (1.0 + eps_break_plastic) * b.l0;
        r.eps_break = eps_break_plastic;
    }
}

double g_eps_break_elastic = 0.0; // --eps-break-e: eps_break упругой нити (0 -- как у материала сцены); контроль для --eps-break-p
double g_eps_break_plastic = 0.0; // --eps-break-p: eps_break пластичной нити (0 -- как у материала сцены); SY/H считаются от исходного eps_break

RunResult run_case(Config config, Level level, double V0, bool plastic, double t_max, double sy_frac, double h_frac,
                   double time_limit, const fs::path& workdir, std::size_t obs_every, bool stop_on_outcome)
{
    RunResult r;
    try
    {
        GuoScene sc;
        {
            e2e::ScopedSilence silence(true);
            sc = build_guo_scene(config, level, V0, workdir, kMuDefault, true);
        }
        if (plastic) apply_plasticity_via_fib(sc, workdir, sy_frac, h_frac, g_eps_break_plastic, r);
        else
        {
            r.E = sc.model.materials.at(0).E;
            r.eps_break = sc.model.materials.at(0).eps_break;
            if (g_eps_break_elastic > 0.0)
            {
                for (fiber::FiberMaterial& m : sc.model.materials) m.eps_break = g_eps_break_elastic;
                for (fiber::FiberBond& b : sc.model.topology.bonds) b.l_break = (1.0 + g_eps_break_elastic) * b.l0;
                r.eps_break = g_eps_break_elastic;
            }
        }
        sc.cfg.t_total = t_max;
        sc.cfg.save_t = 1.0e30;
        sc.cfg.print_step = 1000000000;
        double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal), time = 0.0;
        RunGuard guard(time_limit, 2.0e6, 200);
        MonitorOptions mo;
        mo.csv_path = (workdir / "monitor.csv").string();
        mo.guard = &guard;
        mo.n_mon = obs_every * 5;
        mo.n_diag = 50;
        mo.t_max = t_max;
        mo.stop_on_outcome = stop_on_outcome;
        GuoShotMonitor monitor(sc, mo);
        FiberObserver obs(sc.model, obs_every);
        guard.set_step_observer([&](double d, double t) {
            monitor.on_step(d, t);
            obs.on_step(d, t);
        });
        const auto w0 = std::chrono::steady_clock::now();
        {
            e2e::ScopedSilence silence(true);
            time_integration(sc.cfg, sc.model, sc.ntotal, dt, time, sc.spec.dx_s, &guard);
        }
        r.wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
        monitor.finalize();
        obs.sample();
        const MonitorStats& st = monitor.stats();
        r.steps = st.steps;
        r.t_end = st.t_end;
        r.n_bonds_0 = st.n_bonds_0;
        r.n_bonds_end = st.n_bonds_end;
        r.v_z_end = st.v_z_end;
        r.z_c_end = st.z_c_end;
        r.V_R = st.V_R;
        r.outcome = st.outcome;
        r.ke_ball_0 = st.KE_s_0;
        r.ke_ball_end = st.KE_s_end;
        r.F_f_z_max = st.F_f_z_max;
        r.F_s_z_max = st.F_s_z_max;
        r.w_f_max = st.w_f_max;
        r.eps_max = obs.eps_max();
        r.epsp_max = obs.epsp_max();
        r.sig_max = obs.sig_max();
        r.wp = obs.wp();
        r.wv = obs.wv();
        // Конечность состояния нити и шара.
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < sc.model.fiber.size(); ++i)
                if (!std::isfinite(sc.model.fiber.x[d][i]) || !std::isfinite(sc.model.fiber.vx[d][i])) r.finite = false;
        for (const fiber::FiberBond& b : sc.model.topology.bonds)
            if (!std::isfinite(b.st.eps_p) || !std::isfinite(b.st.w_p) || !std::isfinite(b.st.w_v)) r.finite = false;
        r.ok = true;
    }
    catch (const std::exception& ex)
    {
        r.error = ex.what();
    }
    return r;
}

void print_row(const char* tag, const RunResult& r)
{
    std::printf("  %-9s: шагов %zu, t=%s с, wall %.1f с; v_z шара %s м/с, z_c %s м; KE шара %s -> %s Дж; связей %zu -> %zu "
                "(разрывов %zu); eps_max %s, sigma_ep_max %s Па, eps_p_max %s; W_p %s Дж, W_v %s Дж; F_clamp_f_z_max %s Н; "
                "w_f_max %s м; исход %s\n",
                tag, r.steps, fmtg(r.t_end, 4).c_str(), r.wall, fmtg(r.v_z_end, 6).c_str(), fmtg(r.z_c_end, 6).c_str(),
                fmtg(r.ke_ball_0, 5).c_str(), fmtg(r.ke_ball_end, 5).c_str(), r.n_bonds_0, r.n_bonds_end,
                r.n_bonds_0 - r.n_bonds_end, fmtg(r.eps_max, 5).c_str(), fmtg(r.sig_max, 5).c_str(),
                fmtg(r.epsp_max, 5).c_str(), fmtg(r.wp, 5).c_str(), fmtg(r.wv, 5).c_str(), fmtg(r.F_f_z_max, 5).c_str(),
                fmtg(r.w_f_max, 5).c_str(), outcome_name(r.outcome));
}

void write_summary(const fs::path& out, const char* cfg, const char* lvl, double V0, double t_max, const RunResult& e,
                   const RunResult& p)
{
    std::ofstream f(out);
    f << "config=" << cfg << " level=" << lvl << " V0=" << V0 << " t_max=" << t_max << " SY=" << p.SY << " H=" << p.H
      << " E=" << p.E << " eps_break=" << p.eps_break << "\n";
    f << "metric elastic plastic\n";
    auto row = [&](const char* n, double a, double b) { f << n << " " << fmtg(a, 8) << " " << fmtg(b, 8) << "\n"; };
    row("steps", static_cast<double>(e.steps), static_cast<double>(p.steps));
    row("t_end_s", e.t_end, p.t_end);
    row("wall_s", e.wall, p.wall);
    row("v_z_end_mps", e.v_z_end, p.v_z_end);
    row("V_R_mps", e.V_R, p.V_R);
    row("z_c_end_m", e.z_c_end, p.z_c_end);
    row("ke_ball_0_J", e.ke_ball_0, p.ke_ball_0);
    row("ke_ball_end_J", e.ke_ball_end, p.ke_ball_end);
    row("ke_ball_lost_J", e.ke_ball_0 - e.ke_ball_end, p.ke_ball_0 - p.ke_ball_end);
    row("bonds_broken", static_cast<double>(e.n_bonds_0 - e.n_bonds_end), static_cast<double>(p.n_bonds_0 - p.n_bonds_end));
    row("eps_max", e.eps_max, p.eps_max);
    row("sigma_ep_max_Pa", e.sig_max, p.sig_max);
    row("eps_p_max", e.epsp_max, p.epsp_max);
    row("W_p_J", e.wp, p.wp);
    row("W_v_J", e.wv, p.wv);
    row("F_clamp_fiber_z_max_N", e.F_f_z_max, p.F_f_z_max);
    row("F_clamp_solid_z_max_N", e.F_s_z_max, p.F_s_z_max);
    row("w_f_max_m", e.w_f_max, p.w_f_max);
    f << "outcome " << outcome_name(e.outcome) << " " << outcome_name(p.outcome) << "\n";
}

bool parse_arg(const std::string& a, const char* key, std::string& out)
{
    const std::string k = std::string("--") + key + "=";
    if (a.rfind(k, 0) != 0) return false;
    out = a.substr(k.size());
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    Config config = Config::T1;
    Level level = Level::Smoke;
    double V0 = kSmokeV0, t_max = kSmokeTTotal, sy_frac = 0.2, h_frac = 0.1, time_limit = 0.0;
    std::string out, only;
    bool long_run = false;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        std::string v;
        if (parse_arg(a, "config", v))
        {
            if (!parse_config(v, config)) return 2;
            long_run = true;
        }
        else if (parse_arg(a, "level", v))
        {
            if (!parse_level(v, level)) return 2;
            long_run = true;
        }
        else if (parse_arg(a, "v0", v)) V0 = std::stod(v), long_run = true;
        else if (parse_arg(a, "t-max", v)) t_max = std::stod(v), long_run = true;
        else if (parse_arg(a, "sy-frac", v)) sy_frac = std::stod(v);
        else if (parse_arg(a, "h-frac", v)) h_frac = std::stod(v);
        else if (parse_arg(a, "eps-break-p", v)) g_eps_break_plastic = std::stod(v);
        else if (parse_arg(a, "eps-break-e", v)) g_eps_break_elastic = std::stod(v);
        else if (parse_arg(a, "only", v)) only = v; // elastic | plastic: один прогон из пары (длинный режим)
        else if (parse_arg(a, "time-limit", v)) time_limit = std::stod(v);
        else if (parse_arg(a, "out", v)) out = v, long_run = true;
        else
        {
            std::fprintf(stderr, "Неизвестный аргумент: %s\n", a.c_str());
            return 2;
        }
    }

    std::printf("=== Т16 Guo 2019 %s, уровень %s: упругая нить против упруго-пластичной (SY = %.3g·E·eps_break, H = %.3g·E; "
                "параметры ИЛЛЮСТРАТИВНЫЕ, не откалиброваны), V0 = %.1f м/с, t_max = %s с ===\n",
                config_name(config), level_name(level), sy_frac, h_frac, V0, fmtg(t_max, 4).c_str());

    const std::size_t obs_every = long_run ? 5 : 2;
    RunResult re, rp;
    if (long_run)
    {
        const fs::path dir = out.empty() ? fs::path("fiber_plasticity_guo") : fs::path(out);
        fs::create_directories(dir / "elastic");
        fs::create_directories(dir / "plastic");
        const bool do_e = only.empty() || only == "elastic", do_p = only.empty() || only == "plastic";
        if (do_e)
        {
            std::printf("[упругая]\n");
            re = run_case(config, level, V0, false, t_max, sy_frac, h_frac, time_limit, dir / "elastic", obs_every, true);
            if (re.ok) print_row("упругая", re);
            else std::printf("  ОШИБКА: %s\n", re.error.c_str());
        }
        if (do_p)
        {
            std::printf("[пластичная]\n");
            rp = run_case(config, level, V0, true, t_max, sy_frac, h_frac, time_limit, dir / "plastic", obs_every, true);
            if (rp.ok) print_row("пластичная", rp);
            else std::printf("  ОШИБКА: %s\n", rp.error.c_str());
        }
        if (!only.empty()) return ((!do_e || re.ok) && (!do_p || rp.ok)) ? 0 : 1;
        if (re.ok && rp.ok) write_summary(dir / "summary.txt", config_name(config), level_name(level), V0, t_max, re, rp);
        return (re.ok && rp.ok) ? 0 : 1;
    }

    // ---------------------------------------------------------------- smoke (ctest)
    e2e::ScopedWorkDir wd_e("sph_fpl_guo", "elastic");
    e2e::ScopedWorkDir wd_p("sph_fpl_guo", "plastic");
    std::printf("[упругая]\n");
    re = run_case(config, level, V0, false, t_max, sy_frac, h_frac, time_limit, wd_e.dir(), obs_every, false);
    std::printf("[пластичная]\n");
    rp = run_case(config, level, V0, true, t_max, sy_frac, h_frac, time_limit, wd_p.dir(), obs_every, false);
    check(re.ok, "упругий прогон без исключений" + (re.ok ? std::string() : ": " + re.error));
    check(rp.ok, "пластичный прогон без исключений (директива PLASTIC разобрана)" + (rp.ok ? std::string() : ": " + rp.error));
    if (!re.ok || !rp.ok)
    {
        std::printf("\n=== Т16: ПРОВАЛЕН ===\n");
        return 1;
    }
    print_row("упругая", re);
    print_row("пластичная", rp);

    check(re.finite && rp.finite, "состояние нити конечно (нет NaN/inf)");
    check(re.steps == rp.steps || std::fabs(static_cast<double>(re.steps) - static_cast<double>(rp.steps)) < 0.2 * static_cast<double>(re.steps),
          "число шагов сопоставимо (dt не обрушился из-за пластики)");
    check(re.epsp_max == 0.0 && re.wp == 0.0 && re.wv == 0.0, "упругая нить: eps_p = 0, W_p = W_v = 0");
    check(rp.epsp_max > 0.0, "пластичная нить: eps_p > 0 (проход Postprocessing работает на сцене Guo), eps_p_max = " + fmtg(rp.epsp_max, 4));
    check(rp.wp + rp.wv > 0.0, "пластичная нить: W_p + W_v > 0, = " + fmtg(rp.wp + rp.wv, 4) + " Дж");
    check(rp.sig_max <= re.sig_max * (1.0 + 1.0e-9), "пиковое σ_ep нити пластичной <= упругой: " + fmtg(rp.sig_max, 4) + " <= " + fmtg(re.sig_max, 4) + " Па");
    check(rp.sig_max < 0.9 * re.sig_max, "пиковое σ_ep нити пластичной заметно ниже упругой (< 0.9)");
    const double fmax = 2.0 * V0; // |v_z| шара не больше 2 V0
    check(std::fabs(rp.v_z_end) < fmax && std::fabs(re.v_z_end) < fmax, "скорости шара ограничены (ничего не взорвалось)");
    check(rp.ke_ball_end <= rp.ke_ball_0 * (1.0 + 1e-9), "кинетическая энергия шара не растёт (пластичная)");
    const long long broken_e = static_cast<long long>(re.n_bonds_0) - static_cast<long long>(re.n_bonds_end);
    const long long broken_p = static_cast<long long>(rp.n_bonds_0) - static_cast<long long>(rp.n_bonds_end);
    check(broken_p <= broken_e + 5 || static_cast<double>(rp.n_bonds_end) > 0.9 * static_cast<double>(rp.n_bonds_0),
          "число разрывов не взорвалось (пластичная " + std::to_string(broken_p) + ", упругая " + std::to_string(broken_e) + ")");
    {
        // Линейное упрочнение, без RATE/VISCO: пиковое σ_ep связи = SY + H·ε_p той же связи; максимум по связям и моментам
        // наблюдения совпадает с SY + H·ε_p_max (σ растёт только при пластическом нагружении, где σ = SY + H·ε_p).
        const double sig_pred = rp.SY + rp.H * rp.epsp_max;
        const double rel = std::fabs(rp.sig_max - sig_pred) / sig_pred;
        check(rel < 1.0e-6, "σ_max согласовано с SY + H·ε_p_max: " + fmtg(rp.sig_max, 6) + " vs " + fmtg(sig_pred, 6) + " Па (отн. " +
                                fmtg(rel, 3) + ")");
    }
    const double dv = rp.v_z_end - re.v_z_end;
    std::printf("  справка : Δv_z (пластичная - упругая) = %s м/с; E_kin шара потеряно: упругая %s Дж, пластичная %s Дж\n",
                fmtg(dv, 4).c_str(), fmtg(re.ke_ball_0 - re.ke_ball_end, 4).c_str(),
                fmtg(rp.ke_ball_0 - rp.ke_ball_end, 4).c_str());

    std::printf("\n=== Т16 Guo: %s (провалов: %d) ===\n", e2e::failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
