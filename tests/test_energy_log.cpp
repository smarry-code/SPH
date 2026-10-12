// tests/test_energy_log.cpp
// Лог механической энергии (include/energy_log.hpp, сборка с SPH_ENERGY_LOG).
//
// Сцены -- smoke Guo 2019 (tests/guo2019_scene.hpp): шар, Al-пластина, ткань с
// трением mu = 0.3, зажимы, контакт Solid-Solid типа 1, XSPH. На smoke-скорости
// связи не рвутся, поэтому для канала D_break -- отдельный прогон 1T с порогом
// разрыва 0.5 %.
//
//   E1 формат: заголовок, строка step=0, далее строки ровно на шагах k, 2k, ...;
//   E2 запасы последней строки совпадают с независимым подсчётом по состоянию
//      после прогона (KE Solid/Fiber, sum m*u, энергия связей) -- в прогоне E4;
//   E3 баланс: все значения конечны, каналы потерь >= 0, |residual| мал
//      относительно E_mech0;
//   E4 лог только читает состояние: прогон с логом и без (energy_log_every = 0)
//      даёт побитово одинаковые x, v, u нитей и Solid.
//
// Цель собирается с -DSPH_ENERGY_LOG независимо от опции CMake.

#include "clamp_e2e_common.hpp"
#include "e2e_harness.hpp"
#include "guo2019_scene.hpp"
#include "run_guard.hpp"
#include "time_integration.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef SPH_ENERGY_LOG
#error "test_energy_log собирается только с -DSPH_ENERGY_LOG"
#endif

static_assert(dim == 3, "test_energy_log требует сборки с SPH_DIM=3");

namespace fs = std::filesystem;
using e2e::check;
using e2e::note;
using namespace guo2019;
using namespace clamp_e2e;

namespace
{

constexpr std::size_t kEvery = 20;
// Допуски невязки баланса относительно E_mech0. residual -- то, что не вошло ни
// в один канал (ошибка схемы PC, h-update, зажимы). Факт: 1T ~0.02 % E0, A_1T
// (удар шара о пластину с пластикой) ~0.5 % E0. Каналы трения, контакта и XSPH --
// работа на средней скорости шага; при прежнем разбиении «вклад после остальных
// сил» невязка A_1T была ~6 % E0, а на полном выстреле 1T трение превышало E0.
constexpr double kResidualTolFabric = 0.005;
constexpr double kResidualTolPlate = 0.02;

struct Row
{
    std::vector<double> v;
};

struct Csv
{
    std::string header;
    std::vector<Row> rows;
};

Csv read_csv(const fs::path& p)
{
    Csv c;
    std::ifstream in(p);
    std::getline(in, c.header);
    std::string line;
    while (std::getline(in, line))
    {
        Row r;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) r.v.push_back(std::stod(cell));
        c.rows.push_back(r);
    }
    return c;
}

// Номера колонок -- по заголовку energy_log.hpp.
enum Col
{
    kStep, kTime, kDt, kKeS, kKeF, kUInt, kUShear, kUBond, kUHff, kUHfs, kEMech,
    kDFric, kDXsph, kDCss, kDBrk, kSUcl, kQav, kWpl, kResid, kNCols
};

const char* kHeader = "step,time,dt,KE_solid,KE_fiber,U_int,U_shear_el,U_bond,U_hertz_ff,U_hertz_fs,"
                      "E_mech,D_friction,D_xsph,D_contact_ss,D_break,S_u_clamp,Q_av,W_plastic,residual";

bool close_rel(double a, double b, double tol = 1e-12)
{
    return std::fabs(a - b) <= tol * std::max({std::fabs(a), std::fabs(b), 1e-300});
}

std::string g(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.6g", v);
    return buf;
}

GuoScene make_smoke_scene(Config config, const fs::path& dir)
{
    GuoScene sc;
    {
        e2e::ScopedSilence silence(true);
        sc = build_guo_scene(config, Level::Smoke, kSmokeV0, dir, kMuDefault, true);
    }
    sc.cfg.t_total = kSmokeTTotal;
    sc.cfg.save_t = 1.0; // кадры не нужны
    sc.cfg.print_step = 1000000000;
    return sc;
}

// E1-E3: полный smoke-прогон конфигурации с логом каждые kEvery шагов.
// eps_break > 0 -- порог разрыва всех связей ткани (1+eps_break)*l0 вместо
// материального: на smoke-скорости связи рвутся, и работает канал D_break.
void test_balance(Config config, double time_limit, double residual_tol, double eps_break = 0.0)
{
    std::printf("\n[E1-E3] Лог энергии, Guo 2019 smoke %s%s, k = %zu\n", config_name(config),
                eps_break > 0.0 ? (", eps_break = " + g(eps_break)).c_str() : "", kEvery);
    e2e::ScopedWorkDir wd("sph_energy_log", std::string(config_name(config)) + (eps_break > 0.0 ? "_brk" : ""));
    GuoScene sc = make_smoke_scene(config, wd.dir());
    if (eps_break > 0.0)
        for (auto& b : sc.model.topology.bonds) b.l_break = (1.0 + eps_break) * b.l0;
    const std::size_t bonds0 = sc.model.topology.bonds.size();
    sc.cfg.energy_log_every = kEvery;
    sc.cfg.energy_log_file = (wd.dir() / "energy_log.csv").string();

    double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal), time = 0.0;
    RunGuard guard(time_limit, 2.0e6, 200);
    bool threw = false;
    std::string what;
    {
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(sc.cfg, sc.model, sc.ntotal, dt, time, sc.spec.dx_s, &guard);
        }
        catch (const std::exception& ex)
        {
            threw = true;
            what = ex.what();
        }
    }
    if (!check(!threw && !guard.failed(), "прогон без исключений и аварий сторожа" + (threw ? ": " + what : guard.message())))
        return;
    const std::size_t steps = guard.steps();
    note("шагов " + std::to_string(steps) + ", t = " + g(time));

    // --- E1 формат ---
    const Csv csv = read_csv(wd.dir() / "energy_log.csv");
    check(csv.header == kHeader, "заголовок CSV");
    check(csv.rows.size() == 1 + steps / kEvery,
          "строк " + std::to_string(csv.rows.size()) + " == 1 + steps/k = " + std::to_string(1 + steps / kEvery));
    bool cols_ok = true, steps_ok = true, finite = true;
    for (std::size_t r = 0; r < csv.rows.size(); ++r)
    {
        const auto& v = csv.rows[r].v;
        if (v.size() != kNCols) { cols_ok = false; continue; }
        if (static_cast<std::size_t>(v[kStep]) != r * kEvery) steps_ok = false;
        for (double x : v) finite = finite && std::isfinite(x);
    }
    check(cols_ok, "в каждой строке " + std::to_string(kNCols) + " колонок");
    check(steps_ok, "строки на шагах 0, k, 2k, ...");
    check(finite, "все значения конечны");
    if (!cols_ok || csv.rows.size() < 2) return;

    const auto& first = csv.rows.front().v;
    const auto& last = csv.rows.back().v;
    check(close_rel(first[kKeS], 0.5 * sc.ball_mass * kSmokeV0 * kSmokeV0, 1e-9),
          "KE_solid(0) = 1/2 m_шара V0^2");
    check(last[kUShear] <= last[kUInt] * (1.0 + 1e-12) + 1e-300, "U_shear_el не превышает U_int");

    // --- E3 баланс ---
    const double e0 = first[kEMech];
    std::printf("  info : E_mech0 = %s; конец: KE_s %s KE_f %s U_int %s (сдвиг %s) U_b %s U_hff %s U_hfs %s\n",
                g(e0).c_str(), g(last[kKeS]).c_str(), g(last[kKeF]).c_str(), g(last[kUInt]).c_str(),
                g(last[kUShear]).c_str(), g(last[kUBond]).c_str(), g(last[kUHff]).c_str(), g(last[kUHfs]).c_str());
    std::printf("  info : каналы: D_fric %s D_xsph %s D_css %s D_brk %s S_ucl %s | Q_av %s W_pl %s | residual %s (%s E0)\n",
                g(last[kDFric]).c_str(), g(last[kDXsph]).c_str(), g(last[kDCss]).c_str(), g(last[kDBrk]).c_str(),
                g(last[kSUcl]).c_str(), g(last[kQav]).c_str(), g(last[kWpl]).c_str(), g(last[kResid]).c_str(),
                g(last[kResid] / e0).c_str());
    double max_res = 0.0;
    for (const auto& r : csv.rows) max_res = std::max(max_res, std::fabs(r.v[kResid]));
    check(first[kResid] == 0.0, "residual(0) = 0");
    check(last[kDFric] >= 0.0, "D_friction >= 0");
    check(last[kDCss] >= 0.0, "D_contact_ss >= 0");
    if (eps_break > 0.0)
    {
        const std::size_t broken = bonds0 - sc.model.topology.bonds.size();
        note("порвано связей " + std::to_string(broken));
        check(broken > 0 && last[kDBrk] > 0.0, "есть разрывы и D_break > 0");
    }
    else
    {
        check(last[kDBrk] >= 0.0, "D_break >= 0");
    }
    check(last[kSUcl] >= 0.0, "S_u_clamp >= 0");
    check(last[kQav] >= 0.0, "Q_av >= 0");
    check(last[kDXsph] > 0.0, "D_xsph > 0 (epsilon > 0)");
    check(max_res <= residual_tol * e0,
          "max|residual| = " + g(max_res / e0) + " E0 <= " + g(residual_tol) + " E0");
}

// E4: лог не возмущает траекторию.
struct Snapshot
{
    std::vector<double> s, f;
};

Snapshot snap(const GuoScene& sc)
{
    Snapshot o;
    const auto& sy = sc.model.solid;
    const auto& fb = sc.model.fiber;
    for (std::size_t i = 0; i < sc.ntotal; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d) { o.s.push_back(sy.x[d][i]); o.s.push_back(sy.vx[d][i]); }
        o.s.push_back(sy.u[i]);
        o.s.push_back(sy.sxx[i]);
    }
    for (std::size_t k = 0; k < fb.size(); ++k)
        for (std::size_t d = 0; d < dim; ++d) { o.f.push_back(fb.x[d][k]); o.f.push_back(fb.vx[d][k]); }
    return o;
}

Snapshot run_steps(Config config, std::size_t every, std::size_t steps, const char* tag)
{
    e2e::ScopedWorkDir wd("sph_energy_log", std::string(config_name(config)) + tag);
    GuoScene sc = make_smoke_scene(config, wd.dir());
    sc.cfg.energy_log_every = every;
    sc.cfg.energy_log_file = (wd.dir() / "energy_log.csv").string();
    double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal), time = 0.0;
    RunGuard guard(100.0, 2.0e6, 200);
    guard.set_step_limit(steps);
    {
        e2e::ScopedSilence silence(true);
        time_integration(sc.cfg, sc.model, sc.ntotal, dt, time, sc.spec.dx_s, &guard);
    }
    check(fs::exists(wd.dir() / "energy_log.csv") == (every > 0),
          std::string("energy_log.csv ") + (every > 0 ? "создан" : "не создан") + " (" + tag + ")");
    if (every > 0 && guard.steps() % every == 0)
    {
        // E2: запасы последней строки (шаг steps) против независимого подсчёта.
        const Csv csv = read_csv(wd.dir() / "energy_log.csv");
        const auto& last = csv.rows.back().v;
        const fiber::Model& m = sc.model;
        check(static_cast<std::size_t>(last[kStep]) == guard.steps(), "E2: последняя строка -- последний шаг");
        check(close_rel(last[kKeS], kinetic_solid(m.solid, sc.ntotal)), "E2: KE_solid = sum 1/2 m v^2 Solid");
        check(close_rel(last[kKeF], kinetic_fiber(m.fiber)), "E2: KE_fiber = sum 1/2 m v^2 бусин");
        check(close_rel(last[kUInt], internal_solid(m.solid, sc.ntotal)), "E2: U_int = sum m u");
        check(close_rel(last[kUBond], bond_energy(m)), "E2: U_bond = sum 1/2 EA/l0 (l-l0)^2");
    }
    return snap(sc);
}

void test_no_perturbation(Config config)
{
    std::printf("\n[E2, E4] Лог не возмущает траекторию; сверка запасов (%s, 300 шагов, 1 поток)\n",
                config_name(config));
    // Один поток: свёртка вкладов контакта по буферам потоков не детерминирована
    // (порядок суммирования зависит от расписания), и два прогона БЕЗ лога на
    // нескольких потоках уже расходятся в последних битах.
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    const Snapshot a = run_steps(config, 0, 300, "_off");
    const Snapshot b = run_steps(config, 10, 300, "_on");
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    check(a.s == b.s, "Solid: x, v, u, sxx побитово совпадают");
    check(a.f == b.f, "Fiber: x, v побитово совпадают");
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    double time_limit = 200.0;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a.rfind("--time-limit=", 0) == 0) time_limit = std::stod(a.substr(std::strlen("--time-limit=")));
    }
    test_balance(Config::A_T1, time_limit, kResidualTolPlate);
    test_balance(Config::T1, time_limit, kResidualTolFabric);
    test_balance(Config::T1, time_limit, kResidualTolFabric, 0.005);
    test_no_perturbation(Config::A_T1);
    std::printf("\n%s: провалов %d\n", e2e::failures == 0 ? "OK" : "FAIL", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
