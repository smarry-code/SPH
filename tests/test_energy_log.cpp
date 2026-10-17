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
//   E5 пластичная нить (PLASTIC [+ VISCO], упрочнение линейное): max|residual| <= доли
//      (W_p + W_v) (0.1; 0.35 для Максвелла + разрывы), W_p_fiber > 0 и совпадает с W_p
//      монитора Guo (guo2019_monitor.hpp) -- точно без разрывов; с разрывами лог больше
//      монитора ровно на приращение W_p связей на их последнем шаге (монитор их уже не
//      видит). Контроли: без W_p + W_v и (при разрывах) без D_break_contact невязка
//      выходит за допуск. Столбцы D_break_contact (потенциал Герца, исчезающий вместе с
//      сегментом порванной связи) и clip_loss (отсечка sigma >= 0) закрывают учёт;
//   E4 лог только читает состояние: прогон с логом и без (energy_log_every = 0)
//      даёт побитово одинаковые x, v, u нитей и Solid.
//
// Цель собирается с -DSPH_ENERGY_LOG независимо от опции CMake.

#include "clamp_e2e_common.hpp"
#include "e2e_harness.hpp"
#include "guo2019_monitor.hpp"
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
// Пластичная нить (E5): допуск невязки от собственной работы нити W_p + W_v, а не от E0.
// Факт (1T, plastic): W_p + W_v ~ 0.7-0.8 % E0, max|residual| ~ 0.02-0.04 % E0; при 0.1 (W_p + W_v)
// запас ~2x, а выпадение W_p (вклад 0.4-0.8 % E0) этот допуск нарушает на порядок.
constexpr double kResidualWFrac = 0.1;
// Максвелл + Кельвин + разрывы (eps_break = 0.005, 8 порванных связей): W_p + W_v = 0.12 % E0 (связи рвутся рано,
// нагружены мало), max|residual| = 0.020 % E0 = 0.17 (W_p + W_v) -- остаток O(dt) (сходится ~ dt, см. дизайн,
// раздел 9, риск 4); до учёта D_break_contact было 0.077 % E0. Допуск 0.35 (W_p + W_v) = 0.043 % E0 (запас ~2x).
constexpr double kResidualWFracBrk = 0.35;
// Упругая нить с разрывами (eps_break = 0.005): факт 0.026 % E0 (до D_break_contact 0.072 %), допуск 0.06 % E0.
constexpr double kResidualTolFabricBrk = 6.0e-4;

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
    kDFric, kDXsph, kDCss, kDBrk, kSUcl, kQav, kWpl, kResid, kWpF, kWvF, kBrkC, kClip, kDCf, kNCols
};

const char* kHeader = "step,time,dt,KE_solid,KE_fiber,U_int,U_shear_el,U_bond,U_hertz_ff,U_hertz_fs,"
                      "E_mech,D_friction,D_xsph,D_contact_ss,D_break,S_u_clamp,Q_av,W_plastic,residual,"
                      "W_p_fiber,W_v_fiber,D_break_contact,clip_loss,D_contact_fiber";

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

// E5: независимая запасённая энергия неупругих связей, константы -- здесь (не из energy_log.hpp):
//   U = 1/2 E (eps - eps_p)^2 A l0 при eps > eps_p  +  1/2 q^2 / E1 A l0 (ветвь Максвелла),
//   A = pi d^2 / 4. Связи сумм по топологии текущего состояния.
double bond_energy_inelastic_indep(const fiber::Model& m)
{
    const double pi = std::acos(-1.0);
    double e = 0.0;
    for (const auto& b : m.topology.bonds)
    {
        const fiber::FiberMaterial& mat = m.materials[static_cast<std::size_t>(m.fiber.material_id[b.i])];
        const double A = pi * mat.d * mat.d / 4.0;
        double l2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double r = m.fiber.x[d][b.j] - m.fiber.x[d][b.i];
            l2 += r * r;
        }
        const double eps = (std::sqrt(l2) - b.l0) / b.l0;
        if (eps > b.st.eps_p) e += 0.5 * mat.E * (eps - b.st.eps_p) * (eps - b.st.eps_p) * A * b.l0;
        if (mat.visco.E1 > 0.0) e += 0.5 * b.st.q * b.st.q / mat.visco.E1 * A * b.l0;
    }
    return e;
}

// E1-E3: полный smoke-прогон конфигурации с логом каждые kEvery шагов.
// eps_break > 0 -- порог разрыва всех связей ткани (1+eps_break)*l0 вместо
// материального: на smoke-скорости связи рвутся, и работает канал D_break.
struct PlasticSetup
{
    bool on = false;
    double sy_frac = 0.02; // sy = sy_frac * E * eps_break материала сцены
    double h_frac = 0.1;   // H = h_frac * E
    bool visco = false;    // Максвелл (E1 = 0.2 E, tau = 2 мкс) + Кельвин (eta = 2 Па*с*1e3)
};

struct BalanceResult
{
    double wp_log = 0.0, wp_mon = 0.0, wv_log = 0.0, resid_rel = 0.0;
};

BalanceResult test_balance(Config config, double time_limit, double residual_tol, double eps_break = 0.0,
                           const PlasticSetup& pl = PlasticSetup(), int contact_type = 0)
{
    BalanceResult res;
    std::printf("\n[E1-E3] Лог энергии, Guo 2019 smoke %s%s%s, k = %zu\n", config_name(config),
                eps_break > 0.0 ? (", eps_break = " + g(eps_break)).c_str() : "",
                contact_type != 0 ? (", закон контакта нитей " + std::to_string(contact_type)).c_str() : "", kEvery);
    e2e::ScopedWorkDir wd("sph_energy_log", std::string(config_name(config)) + (eps_break > 0.0 ? "_brk" : "") + (eps_break >= 1.0 ? "_nb" : "") +
                                                    (pl.on ? (pl.visco ? "_plv" : "_pl") : "") +
                                                    (contact_type != 0 ? "_ct" + std::to_string(contact_type) : ""));
    GuoScene sc = make_smoke_scene(config, wd.dir());
    // Закон нормальной силы контактов нитей (оба контакта): 0 -- Герц (запас U_hertz_*),
    // 1-3 -- работа контакта в канале D_contact_fiber.
    sc.cfg.fiber_solid_contact_type = contact_type;
    sc.cfg.fiber_fiber_contact_type = contact_type;
    if (pl.on)
    {
        for (fiber::FiberMaterial& m : sc.model.materials)
        {
            m.plastic.type = fiber::HardeningType::Linear;
            m.plastic.sy = pl.sy_frac * m.E * m.eps_break;
            m.plastic.H = pl.h_frac * m.E;
            if (pl.visco)
            {
                m.visco.E1 = 0.2 * m.E;
                m.visco.tau = 2.0e-6;
                m.visco.eta = 2.0e3;
            }
        }
    }
    if (eps_break > 0.0)
        for (auto& b : sc.model.topology.bonds) b.l_break = (1.0 + eps_break) * b.l0;
    const std::size_t bonds0 = sc.model.topology.bonds.size();
    sc.cfg.energy_log_every = kEvery;
    sc.cfg.energy_log_file = (wd.dir() / "energy_log.csv").string();

    double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal), time = 0.0;
    RunGuard guard(time_limit, 2.0e6, 200);
    // Монитор Guo (W_p по всем связям на каждом шаге) -- независимая сверка W_p_fiber.
    MonitorOptions mo;
    mo.csv_path = (wd.dir() / "monitor.csv").string();
    mo.guard = &guard;
    mo.n_mon = 1000000;
    mo.n_diag = 1000000;
    mo.t_max = kSmokeTTotal;
    mo.stop_on_outcome = false;
    GuoShotMonitor monitor(sc, mo);
    // W_p монитора фиксируется на шагах, кратных kEvery (строки лога): тот же момент, что и в CSV.
    std::size_t obs_n = 0;
    std::vector<double> wp_mon_at; // [r] -- W_p монитора на шаге r*kEvery
    guard.set_step_observer([&](double d, double t) {
        monitor.on_step(d, t);
        if (++obs_n % kEvery == 0) wp_mon_at.push_back(monitor.stats().W_p_total);
    });
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
        return res;
    monitor.finalize();
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
    if (!cols_ok || csv.rows.size() < 2) return res;

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
    if (contact_type == 0)
    {
        check(last[kDCf] == 0.0, "закон Герц: D_contact_fiber = 0");
    }
    else
    {
        std::printf("  info : закон %d: D_contact_fiber %s (%s E0), U_hertz_ff %s, U_hertz_fs %s\n", contact_type,
                    g(last[kDCf]).c_str(), g(last[kDCf] / e0).c_str(), g(last[kUHff]).c_str(), g(last[kUHfs]).c_str());
        check(last[kDCf] > 0.0, "закон 1-3: D_contact_fiber > 0 (сила только на сближении -- рассеяние)");
        check(last[kUHff] == 0.0 && last[kUHfs] == 0.0, "закон 1-3: потенциала Герца в запасе нет");
    }
    if (eps_break > 0.0 && eps_break < 1.0)
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
    const std::size_t broken_pl = bonds0 - sc.model.topology.bonds.size();
    // Без разрывов пластичной нити допуск E5 от W_p + W_v (ниже); с разрывами -- residual_tol (0.5 % E0).
    if (!pl.on)
    {
        const double tol = eps_break > 0.0 && eps_break < 1.0 ? kResidualTolFabricBrk : residual_tol;
        check(max_res <= tol * e0, "max|residual| = " + g(max_res / e0) + " E0 <= " + g(tol) + " E0");
    }
    else if (broken_pl > 0)
        check(max_res <= residual_tol * e0,
              "max|residual| = " + g(max_res / e0) + " E0 <= " + g(residual_tol) + " E0");
    res.resid_rel = last[kResid] / e0;
    res.wp_log = last[kWpF];
    res.wv_log = last[kWvF];
    // строка r лога (r >= 1) -- шаг r*kEvery; wp_mon_at[r-1] -- монитор на том же шаге
    const std::size_t last_row = csv.rows.size() - 1;
    res.wp_mon = (last_row >= 1 && last_row <= wp_mon_at.size()) ? wp_mon_at[last_row - 1] : -1.0;
    check(last[kWpF] >= 0.0 && last[kWvF] >= 0.0, "W_p_fiber, W_v_fiber >= 0");
    if (pl.on)
    {
        std::printf("  info : W_p_fiber %s W_v_fiber %s | W_p монитора %s | eps_p_max %s | residual/E0 %s\n",
                    g(last[kWpF]).c_str(), g(last[kWvF]).c_str(), g(res.wp_mon).c_str(),
                    g(monitor.stats().eps_p_max).c_str(), g(res.resid_rel).c_str());
        check(last[kWpF] > 0.0, "W_p_fiber > 0 (пластика сработала)");
        if (broken_pl == 0)
            check(close_rel(last[kWpF], res.wp_mon, 1e-6), "W_p_fiber совпадает с W_p монитора (1e-6 отн.)");
        else
        {
            // Монитор видит связи после шага и хранит W_p порвавшихся по ключу (i,j), но приращение W_p связи на
            // шаге разрыва он не видит (связь уже удалена), а лог берёт его до удаления. Проверено численно: разность
            // лог - монитор равна сумме приращений W_p порвавшихся связей на их последних шагах с точностью 1e-6
            // отн. (Максвелл + разрывы: 0.416 мДж против 0.340 мДж, разность 0.0757 мДж = 22 %; связи рвутся на
            // стадии активного течения). Отсюда лог >= монитор; верхняя граница 5 % -- только без Максвелла.
            note("порвано связей " + std::to_string(broken_pl) + ": W_p лога больше монитора на W_p разорванных связей");
            check(last[kWpF] >= res.wp_mon * (1.0 - 1e-9) && last[kDBrk] > 0.0, "W_p_fiber >= W_p монитора (разрывы)");
            if (!pl.visco)
                check(last[kWpF] <= res.wp_mon * 1.05, "W_p_fiber - W_p монитора <= 5 % (разрывы, без Максвелла)");
        }
        if (pl.visco) check(last[kWvF] > 0.0, "W_v_fiber > 0 (Максвелл/Кельвин)");

        // E5: допуск невязки пластичной нити -- 0.1 (W_p + W_v), а не 0.5 % E0.
        const double w_sum = last[kWpF] + last[kWvF];
        const double tol_w = (pl.visco && broken_pl > 0 ? kResidualWFracBrk : kResidualWFrac) * w_sum;
        std::printf("  info : E5 max|residual| = %s E0; допуск %s (W_p+W_v) = %s E0; W_p+W_v = %s E0\n",
                    g(max_res / e0).c_str(), g(tol_w / w_sum).c_str(), g(tol_w / e0).c_str(), g(w_sum / e0).c_str());
        check(max_res <= tol_w, "E5: max|residual| <= " + g(pl.visco && broken_pl > 0 ? kResidualWFracBrk : kResidualWFrac) +
                                    " (W_p + W_v)");
        // Контроли: без W_p + W_v (а при разрывах -- и без D_break_contact) невязка должна выходить за тот же допуск
        // (иначе тест слеп к выпадению соответствующего канала). Максимум по строкам лога.
        double max_noW = 0.0, max_noBc = 0.0;
        for (const auto& r : csv.rows)
        {
            max_noW = std::max(max_noW, std::fabs(r.v[kResid] - r.v[kWpF] - r.v[kWvF]));
            max_noBc = std::max(max_noBc, std::fabs(r.v[kResid] - r.v[kBrkC]));
        }
        check(max_noW > tol_w, "E5 контроль: max|residual - W_p - W_v| = " + g(max_noW / e0) + " E0 выходит за допуск");
        if (broken_pl > 0 && pl.visco)
        {
            std::printf("  info : E5 D_break_contact = %s E0 (D_break = %s E0), clip_loss = %s E0\n", g(last[kBrkC] / e0).c_str(),
                        g(last[kDBrk] / e0).c_str(), g(last[kClip] / e0).c_str());
            check(max_noBc > tol_w, "E5 контроль: max|residual - D_break_contact| = " + g(max_noBc / e0) + " E0 выходит за допуск");
        }
        // Отсечка sigma >= 0: без Максвелла/Кельвина sigma_модели > 0 при eps > eps_p, потерь нет; с ними --
        // малая поправка (<= 5 % W_v), не зависящая от dt.
        if (!pl.visco) check(last[kClip] == 0.0, "clip_loss = 0 без Максвелла/Кельвина");
        else check(std::fabs(last[kClip]) <= 0.05 * last[kWvF], "|clip_loss| <= 5 % W_v (малая поправка отсечки)");

        // E5: неупругий U_bond (eps_p > 0, q != 0) против независимой формулы; только на шаге строки лога.
        if (broken_pl == 0 && steps % kEvery == 0)
        {
            std::size_t n_eps = 0, n_both = 0;
            for (const auto& b : sc.model.topology.bonds)
            {
                if (b.st.eps_p > 0.0) ++n_eps;
                if (b.st.eps_p > 0.0 && b.st.q != 0.0) ++n_both;
            }
            const double ub_ind = bond_energy_inelastic_indep(sc.model);
            std::printf("  info : E5 U_bond лога %s, независимо %s; связей eps_p>0: %zu, с q!=0: %zu\n",
                        g(last[kUBond]).c_str(), g(ub_ind).c_str(), n_eps, n_both);
            check(n_eps > 0 && (!pl.visco || n_both > 0), "E5: есть связи с eps_p > 0 (и q != 0 для Максвелла)");
            check(close_rel(last[kUBond], ub_ind, 1e-9), "E5: U_bond лога = 1/2 E(eps-eps_p)^2 A l0 + 1/2 q^2/E1 A l0");
        }
    }
    else
    {
        check(last[kWpF] == 0.0 && last[kWvF] == 0.0, "упругая нить: W_p_fiber = W_v_fiber = 0");
        check(last[kClip] == 0.0, "упругая нить: clip_loss = 0");
    }
    return res;
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
    {
        PlasticSetup pl;
        pl.on = true;
        test_balance(Config::T1, time_limit, kResidualTolFabric, 0.0, pl);
        test_balance(Config::T1, time_limit, kResidualTolFabric, 1.0, pl); // l_break = 2 l0: без разрывов, точное совпадение с монитором
        pl.visco = true;
        test_balance(Config::T1, time_limit, kResidualTolFabric, 0.0, pl);
        // Максвелл + разрывы: баланс с D_break и W_v при порванных связях.
        test_balance(Config::T1, time_limit, kResidualTolFabric, 0.005, pl);
    }
    // Законы контакта нитей 1-3: энергия контакта -- в канале D_contact_fiber, баланс сходится.
    for (int ct = 1; ct <= 3; ++ct) test_balance(Config::T1, time_limit, kResidualTolFabric, 0.0, PlasticSetup(), ct);
    test_balance(Config::T1, time_limit, kResidualTolFabric, 0.005, PlasticSetup(), 1);
    test_no_perturbation(Config::A_T1);
    std::printf("\n%s: провалов %d\n", e2e::failures == 0 ? "OK" : "FAIL", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
