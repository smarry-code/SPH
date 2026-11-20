#pragma once
// tests/clamp_e2e_common.hpp
// Общая обвязка e2e-тестов закрепления (tests/test_e2e_clamp.cpp): опции
// прогона, запуск боевого time_integration(cfg, fiber::Model&, ...) с
// наблюдателем шага, монитор §1.2 и общие проверки §1.3 постановки
// (notes/fabrics/clamp-tests.md; сцены ES*/EF* -- tests/clamp_e2e_solid.hpp и
// tests/clamp_e2e_fiber.hpp, числа кейсов -- tests/clamp_case.hpp).
//
// Почему наблюдатель, а не прогон по одному шагу. Каждый вызов
// time_integration создаёт новый решатель, и первый шаг PredictorCorrector --
// полушаг (v += dt/2*a, notes/fiber-model/fiber-problem-statement.md, «ODE solver»). Прогон по одному шагу
// сделал бы полушагом КАЖДЫЙ шаг. Поэтому сцена идёт одним непрерывным
// прогоном, а монитор снимается в RunGuard::set_step_observer после
// завершения шага (advance, постобработка, time += dt), когда
// model.clamp_report уже относится к этому шагу.
//
// Баланс импульса на шаге n (ТЗ §10 F4/S4): P^{n+1} - P^n = h*(R_solid + R_fiber),
// h = dt/2 на первом шаге PC и dt иначе. Хранимая скорость PC -- база
// v^{n+1} = v^n + h*a(предиктор), поэтому соотношение точное и на шагах n > 1.

#include "fiber/model.hpp"
#include "input.hpp"
#include "interactions/clamp.hpp"
#include "interactions/fiber_bond.hpp"
#include "particle_system.hpp"
#include "run_guard.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"

#include "e2e_harness.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#ifndef SPH_TEST_DATA_DIR
// Запуск из корня репозитория без CMake: путь относительно текущего каталога.
#define SPH_TEST_DATA_DIR "tests/data"
#endif

namespace clamp_e2e
{

using e2e::check;
using e2e::note;
using e2e::warn;

// ===========================================================================
// Опции прогона (разбор -- в main, tests/test_e2e_clamp.cpp)
// ===========================================================================
struct Options
{
    std::string case_name;   // ES1, EF1, ...
    int solver = 0;          // 0 = PredictorCorrector, 1 = явный Эйлер
    bool full = false;       // --full: полный режим (метка slow)
    std::string variant;     // --variant=...: вариант сцены (pull, low-diss, b, ...)
    bool keep_workdir = false;
    double time_limit = 600; // бюджет на ОДИН прогон time_integration, с
    std::string monitor_csv; // --monitor=ФАЙЛ: записать монитор §1.2 в CSV
    // --frames-dir=DIR: сохранить кадры для визуализации (tools/clamp_animate.py):
    // sph_output_*.txt, fiber_output_*.txt, clamp_mask.txt, reaction.csv --
    // в DIR/<тег прогона>. --frames=N -- число кадров на прогон (по t_end).
    std::string frames_dir;
    std::size_t frames = 80;
};

inline const char* solver_name(int s) { return s == 0 ? "pc" : "euler"; }

inline std::string fmt(double v, int prec = 4)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*e", prec, v);
    return buf;
}

// ===========================================================================
// Конфиг решателя: боевой parse_sph_config() по tests/data/clamp_initial_conditions.txt
// ===========================================================================
inline std::filesystem::path config_path()
{
    return std::filesystem::path(SPH_TEST_DATA_DIR) / "clamp_initial_conditions.txt";
}

// Конфиг с разобранным файлом. Бросает, если файла нет: parse_sph_config при
// неоткрывшемся файле молча оставляет умолчания, и тест прошёл бы на чужих
// параметрах.
inline SimulationConfig load_config(int solver, std::size_t maxn)
{
    const auto p = config_path();
    if (!std::filesystem::exists(p)) throw std::runtime_error("нет файла конфигурации " + p.string());
    SimulationConfig cfg;
    cfg.apply_defaults();
    {
        e2e::ScopedSilence silence(true);
        parse_sph_config(cfg, p.string());
    }
    cfg.ode_solver = solver;
    cfg.maxn = std::max<std::size_t>(maxn, 1);
    cfg.max_interaction = 100 * cfg.maxn;
    cfg.moni_particle = 0;
    cfg.save_t = 1.0e30;            // кадр только на первом шаге (во временный каталог)
    cfg.print_step = 1000000000;
    return cfg;
}

// ===========================================================================
// Величины состояния
// ===========================================================================
inline void momentum_solid(const ParticleSystem& s, std::size_t n, double (&P)[3])
{
    P[0] = P[1] = P[2] = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < dim; ++d) P[d] += s.mass[i] * s.vx[d][i];
}

inline void momentum_fiber(const fiber::FiberState& f, double (&P)[3])
{
    P[0] = P[1] = P[2] = 0.0;
    for (std::size_t k = 0; k < f.size(); ++k)
        for (std::size_t d = 0; d < dim; ++d) P[d] += f.mass[k] * f.vx[d][k];
}

inline double kinetic_solid(const ParticleSystem& s, std::size_t n)
{
    double e = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < dim; ++d) e += 0.5 * s.mass[i] * s.vx[d][i] * s.vx[d][i];
    return e;
}

inline double kinetic_fiber(const fiber::FiberState& f)
{
    double e = 0.0;
    for (std::size_t k = 0; k < f.size(); ++k)
        for (std::size_t d = 0; d < dim; ++d) e += 0.5 * f.mass[k] * f.vx[d][k] * f.vx[d][k];
    return e;
}

// Внутренняя энергия Solid sum m*u (§1.2).
inline double internal_solid(const ParticleSystem& s, std::size_t n)
{
    double e = 0.0;
    for (std::size_t i = 0; i < n; ++i) e += s.mass[i] * s.u[i];
    return e;
}

// Упругая энергия связей sum 1/2 (EA/l0)(l - l0)^2 [l > l0] (связь работает только на растяжение);
// для неупругих материалов -- с вычетом eps_p и пружиной Максвелла (см. ниже).
inline double bond_energy(const fiber::Model& m)
{
    double e = 0.0;
    for (const auto& b : m.topology.bonds)
    {
        const fiber::FiberMaterial& mat = m.materials[static_cast<std::size_t>(fiber::bond_material_id(b, m.fiber))];
        double l2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double r = m.fiber.x[d][b.j] - m.fiber.x[d][b.i];
            l2 += r * r;
        }
        const double l = std::sqrt(l2);
        if (mat.inelastic())
        {
            // Упруго-пластичная связь: запасённая энергия 0.5·E·(ε-ε_p)²·A·l0 (ε > ε_p) плюс энергия пружины
            // ветви Максвелла 0.5·q²/E1·A·l0. Пластическая работа W_p и вязкая W_v в U_b НЕ входят (см. fiber_system).
            const double A = interactions::fiber_cross_section_area(mat);
            const double eps = (l - b.l0) / b.l0;
            if (eps > b.st.eps_p) e += 0.5 * mat.E * (eps - b.st.eps_p) * (eps - b.st.eps_p) * A * b.l0;
            if (mat.visco.E1 > 0.0) e += 0.5 * b.st.q * b.st.q / mat.visco.E1 * A * b.l0;
        }
        else if (l > b.l0)
        {
            const double k = mat.E * interactions::fiber_cross_section_area(mat) / b.l0;
            e += 0.5 * k * (l - b.l0) * (l - b.l0);
        }
    }
    return e;
}

// Живучесть (§1.3 п. 1): конечные x, v, rho, p, s; rho > 0, m > 0.
inline bool state_alive(const fiber::Model& m, std::size_t ntotal, std::string* why = nullptr)
{
    const ParticleSystem& s = m.solid;
    for (std::size_t i = 0; i < ntotal; ++i)
    {
        bool ok = s.mass[i] > 0.0 && s.rho[i] > 0.0 && std::isfinite(s.rho[i]) && std::isfinite(s.p[i]) &&
                  std::isfinite(s.u[i]) && std::isfinite(s.sxx[i]) && std::isfinite(s.syy[i]) &&
                  std::isfinite(s.szz[i]) && std::isfinite(s.sxy[i]) && std::isfinite(s.syz[i]) &&
                  std::isfinite(s.sxz[i]);
        for (std::size_t d = 0; d < dim; ++d) ok = ok && std::isfinite(s.x[d][i]) && std::isfinite(s.vx[d][i]);
        if (!ok)
        {
            if (why) *why = "Solid частица " + std::to_string(i);
            return false;
        }
    }
    for (std::size_t k = 0; k < m.fiber.size(); ++k)
    {
        bool ok = m.fiber.mass[k] > 0.0;
        for (std::size_t d = 0; d < dim; ++d) ok = ok && std::isfinite(m.fiber.x[d][k]) && std::isfinite(m.fiber.vx[d][k]);
        if (!ok)
        {
            if (why) *why = "бусина " + std::to_string(k);
            return false;
        }
    }
    return true;
}

// ===========================================================================
// Монитор §1.2
// ===========================================================================
struct StepRecord
{
    std::size_t step = 0;
    double t = 0.0, dt = 0.0, h = 0.0;   // h -- коэффициент приращения скорости на шаге
    double Rs[3] = {}, Rf[3] = {};        // model.clamp_report.solid_force / fiber_force
    double Ps[3] = {}, Pf[3] = {};        // импульс Solid и Fiber ПОСЛЕ шага
    double bal[3] = {};                   // (P^{n+1} - P^n) - h*(Rs + Rf)
    double Ekin = 0.0, Eint = 0.0, Ebond = 0.0;
    double clamp_dx_max = 0.0;            // max |x - x0| по зажатым компонентам
    double clamp_v_max = 0.0;             // max |v| по зажатым компонентам
    double W_clamp = 0.0;                 // sum m v (v - v_prev)/dt по зажатым компонентам
};

class ClampMonitor
{
public:
    // Снимает x0 зажатых компонент и P^0. Вызывать ПОСЛЕ заполнения сцены и
    // model.clamps, ДО прогона.
    ClampMonitor(const fiber::Model& m, std::size_t ntotal, int solver)
        : m_(m), ntotal_(ntotal), solver_(solver)
    {
        collect(m.clamps.solid_mask, ntotal, /*fiber=*/false);
        collect(m.clamps.fiber_mask, m.fiber.size(), /*fiber=*/true);
        totals(P_prev_);
        for (const Comp& c : comps_) v_prev_.push_back(vel(c));
        E0_ = energy_kin() + internal_solid(m.solid, ntotal_) + bond_energy(m);
    }

    // Наблюдатель шага (RunGuard::set_step_observer).
    void on_step(double dt, double time)
    {
        StepRecord r;
        r.step = records_.size() + 1;
        r.t = time;
        r.dt = dt;
        r.h = (solver_ == 0 && records_.empty()) ? 0.5 * dt : dt;
        double Ps[3], Pf[3];
        momentum_solid(m_.solid, ntotal_, Ps);
        momentum_fiber(m_.fiber, Pf);
        for (std::size_t d = 0; d < 3; ++d)
        {
            r.Rs[d] = m_.clamp_report.solid_force[d];
            r.Rf[d] = m_.clamp_report.fiber_force[d];
            r.Ps[d] = Ps[d];
            r.Pf[d] = Pf[d];
            const double P = Ps[d] + Pf[d];
            r.bal[d] = (P - P_prev_[d]) - r.h * (r.Rs[d] + r.Rf[d]);
            P_prev_[d] = P;
            Pmax_ = std::max(Pmax_, std::fabs(P));
        }
        r.Ekin = energy_kin();
        r.Eint = internal_solid(m_.solid, ntotal_);
        r.Ebond = bond_energy(m_);
        for (std::size_t k = 0; k < comps_.size(); ++k)
        {
            const Comp& c = comps_[k];
            const double x = pos(c), v = vel(c);
            if (!(x == c.x0) || !(v == 0.0)) ++stillness_violations_;
            r.clamp_dx_max = std::max(r.clamp_dx_max, std::fabs(x - c.x0));
            r.clamp_v_max = std::max(r.clamp_v_max, std::fabs(v));
            // §1.3 п. 4: sum m v a по зажатым. Сама постановка отмечает, что это следствие
            // неподвижности (v == 0): столбец не независим от проверки выше и ловит только
            // нарушение v == 0.0, а не работу сброса «после шага».
            r.W_clamp += mass(c) * v * (v - v_prev_[k]) / dt;
            v_prev_[k] = v;
        }
        if (!(r.W_clamp == 0.0)) ++work_violations_;
        records_.push_back(r);
    }

    std::function<void(double, double)> observer()
    {
        return [this](double dt, double t) { on_step(dt, t); };
    }

    const std::vector<StepRecord>& records() const { return records_; }
    std::size_t clamped_components() const { return comps_.size(); }
    std::size_t stillness_violations() const { return stillness_violations_; }
    std::size_t work_violations() const { return work_violations_; }
    double initial_energy() const { return E0_; }
    double max_abs_momentum() const { return Pmax_; }

    // max_n |bal^n| (евклидова норма по осям).
    double max_balance_residual() const
    {
        double m = 0.0;
        for (const StepRecord& r : records_)
            m = std::max(m, std::sqrt(r.bal[0] * r.bal[0] + r.bal[1] * r.bal[1] + r.bal[2] * r.bal[2]));
        return m;
    }

    void write_csv(const std::string& path) const
    {
        std::ofstream f(path);
        f << "step,t,dt,h,Rs_x,Rs_y,Rs_z,Rf_x,Rf_y,Rf_z,Ps_x,Ps_y,Ps_z,Pf_x,Pf_y,Pf_z,bal_x,bal_y,bal_z,"
             "Ekin,Eint,Ebond,clamp_dx_max,clamp_v_max,W_clamp\n";
        f.precision(17);
        for (const StepRecord& r : records_)
        {
            f << r.step << ',' << r.t << ',' << r.dt << ',' << r.h;
            for (double v : r.Rs) f << ',' << v;
            for (double v : r.Rf) f << ',' << v;
            for (double v : r.Ps) f << ',' << v;
            for (double v : r.Pf) f << ',' << v;
            for (double v : r.bal) f << ',' << v;
            f << ',' << r.Ekin << ',' << r.Eint << ',' << r.Ebond << ',' << r.clamp_dx_max << ',' << r.clamp_v_max
              << ',' << r.W_clamp << '\n';
        }
    }

private:
    struct Comp
    {
        bool fiber;
        std::size_t i, d;
        double x0;
    };

    void collect(const std::vector<std::uint8_t>& mask, std::size_t n, bool fiber)
    {
        if (mask.empty()) return;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t d = 0; d < dim; ++d)
                if (mask[i] & (1u << d))
                {
                    Comp c{fiber, i, d, 0.0};
                    c.x0 = pos(c);
                    comps_.push_back(c);
                }
    }
    double pos(const Comp& c) const { return c.fiber ? m_.fiber.x[c.d][c.i] : m_.solid.x[c.d][c.i]; }
    double vel(const Comp& c) const { return c.fiber ? m_.fiber.vx[c.d][c.i] : m_.solid.vx[c.d][c.i]; }
    double mass(const Comp& c) const { return c.fiber ? m_.fiber.mass[c.i] : m_.solid.mass[c.i]; }
    double energy_kin() const { return kinetic_solid(m_.solid, ntotal_) + kinetic_fiber(m_.fiber); }
    void totals(double (&P)[3]) const
    {
        double Ps[3], Pf[3];
        momentum_solid(m_.solid, ntotal_, Ps);
        momentum_fiber(m_.fiber, Pf);
        for (std::size_t d = 0; d < 3; ++d) P[d] = Ps[d] + Pf[d];
    }

    const fiber::Model& m_;
    std::size_t ntotal_;
    int solver_;
    std::vector<Comp> comps_;
    std::vector<double> v_prev_;
    std::vector<StepRecord> records_;
    double P_prev_[3] = {};
    double Pmax_ = 0.0;
    double E0_ = 0.0;
    std::size_t stillness_violations_ = 0;
    std::size_t work_violations_ = 0;
};

// ===========================================================================
// Кадры для визуализации (--frames-dir)
// ===========================================================================
// clamp_mask.txt: строки "solid <i> <mask>" и "fiber <k> <mask>" только для
// зажатых частиц (маска -- биты компонент, bit d = компонента d).
inline void write_clamp_masks(const fiber::Model& m, std::size_t ntotal, const std::string& path)
{
    std::ofstream f(path);
    f << "# kind index mask (bit d = компонента d зажата)\n";
    for (std::size_t i = 0; i < std::min(ntotal, m.clamps.solid_mask.size()); ++i)
        if (m.clamps.solid_mask[i] != 0) f << "solid " << i << ' ' << int(m.clamps.solid_mask[i]) << '\n';
    for (std::size_t k = 0; k < m.clamps.fiber_mask.size(); ++k)
        if (m.clamps.fiber_mask[k] != 0) f << "fiber " << k << ' ' << int(m.clamps.fiber_mask[k]) << '\n';
}

// Копирует кадры и маски из рабочего каталога прогона в frames_dir/<tag>
// (при повторе тега -- <tag>_2, _3, ...) и пишет reaction.csv.
inline void save_run_frames(const Options& opt, const std::string& tag, const std::filesystem::path& workdir,
                            const std::vector<std::array<double, 8>>& reaction)
{
    namespace fs = std::filesystem;
    const fs::path base = fs::path(opt.frames_dir) / (opt.case_name + "_" + solver_name(opt.solver) +
                                                      (opt.variant.empty() ? "" : "_" + opt.variant) + "_" + tag);
    fs::path dst = base;
    for (int k = 2; fs::exists(dst); ++k) dst = base.string() + "_" + std::to_string(k);
    fs::create_directories(dst);
    std::size_t nframes = 0;
    for (const auto& e : fs::directory_iterator(workdir))
    {
        const std::string n = e.path().filename().string();
        const bool frame = n.rfind("sph_output_", 0) == 0 || n.rfind("fiber_output_", 0) == 0;
        if (frame || n == "clamp_mask.txt")
        {
            fs::copy_file(e.path(), dst / n, fs::copy_options::overwrite_existing);
            if (n.rfind("sph_output_", 0) == 0) ++nframes;
        }
    }
    std::ofstream f(dst / "reaction.csv");
    f << "step,t,Rs_x,Rs_y,Rs_z,Rf_x,Rf_y,Rf_z\n";
    f.precision(17);
    for (const auto& r : reaction)
    {
        f << static_cast<long long>(r[0]);
        for (std::size_t c = 1; c < r.size(); ++c) f << ',' << r[c];
        f << '\n';
    }
    note("кадры для визуализации (" + std::to_string(nframes) + "): " + dst.string());
}

// ===========================================================================
// Прогон через боевой цикл
// ===========================================================================
struct RunOutcome
{
    bool threw = false;
    std::string what;
    RunGuard::Reason reason = RunGuard::Reason::None;
    std::size_t steps = 0;
    double dt_first = 0.0, dt_min = 0.0;
    double wall_s = 0.0;
    double time = 0.0;
    double t_end = 0.0;
};

// Один непрерывный прогон до t_end (или max_steps шагов, если t_end <= 0).
// Рабочий каталог -- временный (ScopedWorkDir): цикл пишет log.txt и кадры в
// текущий каталог. observer вызывается после каждого шага.
inline RunOutcome run_case(SimulationConfig cfg, fiber::Model& m, std::size_t ntotal, double& dt, double& time,
                           double t_end, std::size_t max_steps, const Options& opt, const std::string& tag,
                           std::function<void(double, double)> observer)
{
    RunOutcome out;
    out.t_end = t_end;
    cfg.t_total = t_end > 0.0 ? t_end : 1.0e30;
    const bool save_frames = !opt.frames_dir.empty();
    // Реакция по шагам для визуализации: снимается своим наблюдателем поверх
    // переданного (сцена про него не знает).
    std::vector<std::array<double, 8>> reaction; // step, t, Rs[3], Rf[3]
    if (save_frames)
    {
        if (t_end > 0.0 && opt.frames > 0) cfg.save_t = t_end / static_cast<double>(opt.frames);
        else if (max_steps > 0 && opt.frames > 0 && dt > 0.0)
            cfg.save_t = dt * static_cast<double>(std::max<std::size_t>(max_steps / opt.frames, 1));
        observer = [inner = std::move(observer), &reaction, &m](double step_dt, double t)
        {
            if (inner) inner(step_dt, t);
            const auto& r = m.clamp_report;
            reaction.push_back({static_cast<double>(reaction.size() + 1), t, r.solid_force[0], r.solid_force[1],
                                r.solid_force[2], r.fiber_force[0], r.fiber_force[1], r.fiber_force[2]});
        };
    }
    RunGuard guard(opt.time_limit, 0.0, 0);
    if (max_steps > 0) guard.set_step_limit(max_steps);
    guard.set_step_observer(std::move(observer));
    e2e::ScopedWorkDir wd("sph_e2e_clamp", tag);
    if (opt.keep_workdir)
    {
        wd.keep();
        note("рабочий каталог сохранён: " + wd.dir().string());
    }
    if (save_frames) write_clamp_masks(m, ntotal, "clamp_mask.txt");
    const auto t0 = std::chrono::steady_clock::now();
    {
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(cfg, m, ntotal, dt, time, /*delta_x=*/0.0, &guard);
        }
        catch (const std::exception& ex)
        {
            out.threw = true;
            out.what = ex.what();
        }
    }
    out.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out.reason = guard.reason();
    out.steps = guard.steps();
    out.dt_first = guard.dt_first();
    out.dt_min = guard.dt_min();
    out.time = time;
    if (save_frames) save_run_frames(opt, tag, wd.dir(), reaction);
    return out;
}

// ===========================================================================
// Общие проверки §1.3
// ===========================================================================
// balance_tol -- абсолютный допуск на max_n |P^{n+1} - P^n - h R| (§1.3 п. 3):
// для сцен только с нитями 1e-12*max|P|, для сцен с Solid -- от невязки той
// же сцены без зажима (выбирает вызывающий и пишет обоснование в clamp_case.hpp).
inline void common_checks(const std::string& tag, const fiber::Model& m, std::size_t ntotal, const ClampMonitor& mon,
                          const RunOutcome& run, double balance_tol)
{
    std::printf("\n[§1.3] Общие проверки%s\n", tag.c_str());
    check(!run.threw, "нет исключений" + (run.threw ? ": " + run.what : std::string()));
    check(!(run.reason == RunGuard::Reason::WallClock || run.reason == RunGuard::Reason::Stalled),
          "сторож не оборвал расчёт (шагов " + std::to_string(run.steps) + ", " + fmt(run.wall_s, 2) + " с)");
    if (run.t_end > 0.0)
        check(run.time >= run.t_end, "расчёт дошёл до t_end = " + fmt(run.t_end) + " (t = " + fmt(run.time) + ")");
    check(run.dt_first > 0.0 && run.dt_min > 0.1 * run.dt_first,
          "dt не схлопнулся (dt_min/dt_first = " + fmt(run.dt_first > 0.0 ? run.dt_min / run.dt_first : 0.0, 3) + ")");
    std::string why;
    check(state_alive(m, ntotal, &why), "живучесть: нет NaN/Inf в x, v, rho, p, s; rho > 0, m > 0" +
                                            (why.empty() ? std::string() : " (" + why + ")"));
    check(mon.records().size() == run.steps, "монитор снят на каждом шаге (" + std::to_string(mon.records().size()) +
                                                 " записей)");
    check(mon.clamped_components() > 0 && mon.stillness_violations() == 0,
          "неподвижность побитово: x == x0 и v == 0.0 у " + std::to_string(mon.clamped_components()) +
              " зажатых компонент на каждом шаге (нарушений " + std::to_string(mon.stillness_violations()) + ")");
    const double res = mon.max_balance_residual();
    check(res <= balance_tol, "баланс импульса по шагам: max |dP - h R| = " + fmt(res) + " <= " + fmt(balance_tol) +
                                  " (max|P| = " + fmt(mon.max_abs_momentum()) + ")");
    check(mon.work_violations() == 0, "нет работы зажима: W_clamp == 0.0 на всём прогоне");
}

// Запись монитора, если просили --monitor.
inline void maybe_write_monitor(const Options& opt, const ClampMonitor& mon, const std::string& suffix = "")
{
    if (opt.monitor_csv.empty()) return;
    const std::string p = opt.monitor_csv + suffix;
    mon.write_csv(p);
    note("монитор записан: " + p);
}

} // namespace clamp_e2e
