// tests/test_fiber_plasticity.cpp
// Unit-тесты чистой функции bond_constitutive_update / bond_stress
// (include/fiber/fiber_plasticity.hpp; notes/fiber-model/fiber-plasticity-design.md,
// раздел 8, тесты Т1-Т7). Решатель не используется.
#include "fiber/fiber_plasticity.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool same(double a, double b) { return !(a < b) && !(a > b); } // точное равенство без -Wfloat-equal

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1e-300, std::fabs(b)); }

constexpr double kE = 2.0e9;
constexpr double kSy = 60.0e6;
constexpr double kH = 0.4e9;

// Владеющая обёртка для тестов; в вызовы функций уходит невладеющий BondLaw.
struct TestLaw
{
    double E = 0.0;
    fiber::PlasticParams plastic;
    fiber::RateParams rate;
    fiber::ViscoParams visco;
    fiber::BondLaw view() const { return fiber::BondLaw(E, plastic, rate, visco); }
};

fiber::BondState upd(const TestLaw& l, const fiber::BondState& s, double e, double dt)
{
    return fiber::bond_constitutive_update(l.view(), s, e, dt);
}

double stress(const TestLaw& l, const fiber::BondState& s, double e, double edot)
{
    return fiber::bond_stress(l.view(), s, e, edot);
}

TestLaw make_law_linear(double H = kH)
{
    TestLaw l;
    l.E = kE;
    l.plastic.type = fiber::HardeningType::Linear;
    l.plastic.sy = kSy;
    l.plastic.H = H;
    return l;
}

// Прогон по списку значений ε (шаг dt); возвращает конечное состояние.
fiber::BondState run_path(const TestLaw& law, const std::vector<double>& eps, double dt,
                          fiber::BondState s = {})
{
    for (double e : eps) s = upd(law, s, e, dt);
    return s;
}

std::vector<double> ramp(double from, double to, int n)
{
    std::vector<double> v;
    for (int k = 1; k <= n; ++k) v.push_back(from + (to - from) * static_cast<double>(k) / static_cast<double>(n));
    return v;
}

// --- Т0: всё выключено => упругий закон E*eps при eps>0, 0 иначе ------------
void test_elastic_equivalence()
{
    std::printf("T0: выключенные параметры == упругий закон\n");
    TestLaw law;
    law.E = kE;
    fiber::BondState s;
    bool ok = true;
    for (double e : {0.001, 0.01, 0.05, 0.2, 0.05, -0.01, 0.0, 0.03})
    {
        s = upd(law, s, e, 1e-6);
        const double sig = stress(law, s, e, 123.0); // eps_dot не влияет при eta = 0
        const double ref = e > 0.0 ? kE * e : 0.0;
        if (sig < ref || sig > ref || s.eps_p > 0.0 || s.q > 0.0 || s.w_p > 0.0 || s.w_v > 0.0) ok = false;
    }
    check(ok, "sigma == E*eps (eps>0), eps_p = q = w = 0, побитово");
}

// --- Т1: монотонное растяжение, линейное упрочнение -------------------------
void test_monotonic_linear()
{
    std::printf("T1: монотонное растяжение, линейное упрочнение\n");
    const TestLaw law = make_law_linear();
    const double dt = 1e-6;
    fiber::BondState s;
    double max_err_ep = 0.0, max_err_sig = 0.0;
    const int n = 500;
    const double emax = 0.05;
    for (int k = 1; k <= n; ++k)
    {
        const double e = emax * static_cast<double>(k) / n;
        s = upd(law, s, e, dt);
        const double ep_ref = std::max(0.0, (kE * e - kSy) / (kE + kH));
        const double sig_ref = kE * (e - ep_ref);
        const double sig = stress(law, s, e, 0.0);
        max_err_ep = std::max(max_err_ep, std::fabs(s.eps_p - ep_ref) / 1e-3);
        max_err_sig = std::max(max_err_sig, std::fabs(sig - sig_ref) / kSy);
    }
    check(max_err_ep < 1e-12, "eps_p == (E*eps - sy)/(E+H), невязка < 1e-12 (отн. 1e-3)");
    check(max_err_sig < 1e-12, "sigma == E(eps - eps_p), невязка < 1e-12 (отн. sy)");
    check(near_rel(stress(law, s, emax, 0.0), kSy + kH * s.eps_p, 1e-12), "sigma = sy + H*eps_p");
    // Работа: int sigma d eps_p = sy*eps_p + H*eps_p^2/2 (трапеция точна)
    check(near_rel(s.w_p, kSy * s.eps_p + 0.5 * kH * s.eps_p * s.eps_p, 1e-12), "w_p = sy*eps_p + H*eps_p^2/2");
}

// --- Т2: разгрузка и повторная нагрузка -------------------------------------
void test_unload_reload()
{
    std::printf("T2: разгрузка / повторная нагрузка\n");
    const TestLaw law = make_law_linear();
    const double dt = 1e-6;
    fiber::BondState s = run_path(law, ramp(0.0, 0.05, 100), dt);
    const double ep = s.eps_p;
    const double sig_max = stress(law, s, 0.05, 0.0);
    check(ep > 0.0, "после нагрузки есть остаточная деформация");

    const double de = 0.004;
    s = upd(law, s, 0.05 - de, dt);
    check(near_rel(s.eps_p, ep, 1e-14), "при разгрузке eps_p не меняется");
    check(near_rel(stress(law, s, 0.05 - de, 0.0), sig_max - kE * de, 1e-12), "наклон разгрузки == E");

    s = upd(law, s, ep * 0.5, dt);
    check(stress(law, s, ep * 0.5, 0.0) < 1e-300, "при eps < eps_p сила == 0");
    s = upd(law, s, ep * 0.5 - 0.01, dt); // даже сжатие
    check(stress(law, s, ep * 0.5 - 0.01, 0.0) < 1e-300 && near_rel(s.eps_p, ep, 1e-14),
          "сжатие не меняет eps_p, сила 0");

    // Повторная нагрузка: упруго до достигнутой eps = 0.05, затем течёт по той же кривой
    s = run_path(law, ramp(ep * 0.5 - 0.01, 0.05 - 1e-9, 200), dt, s);
    check(near_rel(s.eps_p, ep, 1e-12), "до достигнутой eps пластика не возобновляется");
    s = upd(law, s, 0.06, dt);
    const double ep2 = (kE * 0.06 - kSy) / (kE + kH);
    check(near_rel(s.eps_p, ep2, 1e-10), "дальше течёт по той же кривой (упрочнение запомнено)");
}

// --- Т3: идеальная пластичность ---------------------------------------------
void test_perfect_plasticity()
{
    std::printf("T3: идеальная пластичность (H = 0)\n");
    const TestLaw law = make_law_linear(0.0);
    fiber::BondState s;
    double max_sig = 0.0;
    for (double e : ramp(0.0, 0.3, 300))
    {
        s = upd(law, s, e, 1e-6);
        max_sig = std::max(max_sig, stress(law, s, e, 0.0));
    }
    check(max_sig <= kSy * (1.0 + 1e-12), "sigma <= sy при любом eps");
    check(near_rel(s.eps_p, 0.3 - kSy / kE, 1e-12), "eps_p = eps - sy/E");
}

// --- Т4: независимость от dt без скоростной зависимости ---------------------
std::vector<double> t4_path(int scale)
{
    std::vector<double> v = ramp(0.0, 0.08, 40 * scale);
    const auto b = ramp(0.08, 0.03, 20 * scale);
    const auto c = ramp(0.03, 0.15, 60 * scale);
    v.insert(v.end(), b.begin(), b.end());
    v.insert(v.end(), c.begin(), c.end());
    return v;
}

void test_dt_independence()
{
    std::printf("T4: независимость от dt (rate-independent)\n");
    std::vector<TestLaw> laws;
    laws.push_back(make_law_linear());
    {
        TestLaw l = make_law_linear();
        l.plastic.type = fiber::HardeningType::Power;
        l.plastic.K = 300e6;
        l.plastic.n = 0.4;
        laws.push_back(l);
    }
    {
        TestLaw l;
        l.E = kE;
        l.plastic.type = fiber::HardeningType::Table;
        l.plastic.table_eps = {0.0, 0.02, 0.1, 0.3};
        l.plastic.table_sigma = {60e6, 72e6, 95e6, 120e6};
        laws.push_back(l);
    }
    const char* names[] = {"линейное", "степенное", "таблица"};
    for (std::size_t k = 0; k < laws.size(); ++k)
    {
        // Траектория: нагрузка до 0.08, разгрузка до 0.03, нагрузка до 0.15.
        // Возврат на монотонных участках точный => результат не зависит от разбиения и dt.
        const fiber::BondState a = run_path(laws[k], t4_path(1), 1e-6);
        const fiber::BondState b = run_path(laws[k], t4_path(10), 1e-7);
        const double w_tol = k == 0 ? 1e-12 : 2e-3; // линейное: трапеция точна
        check(near_rel(a.eps_p, b.eps_p, 1e-10) && near_rel(a.w_p, b.w_p, w_tol),
              std::string(names[k]) + ": eps_p при dt и dt/10 совпадают (eps_p=" + std::to_string(a.eps_p) + ")");
    }
}

// --- Т5: скоростная зависимость ---------------------------------------------
double flow_stress(const TestLaw& law, double edot)
{
    // Идеальная пластичность по статике: при установившемся течении eps_p' = eps' => sigma = sigma_y(eps').
    const double de = 5e-4; // 2000 шагов (до eps = 1) при любом eps_dot: переходный процесс затухает медленно
    const double dt = de / edot;
    fiber::BondState s;
    double e = 0.0;
    for (int k = 0; k < 2000; ++k)
    {
        e += de;
        s = upd(law, s, e, dt);
    }
    return stress(law, s, e, 0.0);
}

struct RateCase
{
    const char* name;
    fiber::RateParams rate;
    double (*ref)(const fiber::RateParams&, double);
};

double ref_cs(const fiber::RateParams& r, double ed) { return kSy * (1.0 + std::pow(ed / r.a, 1.0 / r.b)); }
double ref_jc(const fiber::RateParams& r, double ed) { return kSy * (1.0 + r.a * std::log(std::max(ed / r.b, 1.0))); }
double ref_pz(const fiber::RateParams& r, double ed) { return kSy + r.a * std::pow(ed, r.b); }

void test_rate_dependence()
{
    std::printf("T5: скоростная зависимость (CS, JC, Perzyna)\n");
    const TestLaw base = make_law_linear(0.0);
    const RateCase cases[] = {
        {"Cowper-Symonds", {fiber::RateModel::CowperSymonds, 40.0, 5.0}, ref_cs},
        {"Johnson-Cook", {fiber::RateModel::JohnsonCook, 0.05, 1e-3}, ref_jc},
        {"Perzyna", {fiber::RateModel::Perzyna, 1e6, 0.5}, ref_pz},
    };
    for (const RateCase& c : cases)
    {
        TestLaw law = base;
        law.rate = c.rate;
        bool formula_ok = true, monotone = true;
        double prev = 0.0;
        for (double ed : {1e-4, 1e-2, 1.0, 10.0, 100.0, 1000.0})
        {
            const double sig = flow_stress(law, ed);
            if (!near_rel(sig, c.ref(c.rate, ed), 1e-6)) formula_ok = false;
            if (sig < prev) monotone = false;
            prev = sig;
        }
        check(formula_ok, std::string(c.name) + ": sigma_flow == sigma_y(eps_dot) по формуле (отн. 1e-6)");
        check(monotone, std::string(c.name) + ": sigma_flow не убывает с eps_dot");
        const double slow = flow_stress(law, 1e-9);
        check(std::fabs(slow - kSy) / kSy < 0.1 && slow >= kSy * (1.0 - 1e-12),
              std::string(c.name) + ": eps_dot->0 стремится к статическому sy (sigma=" + std::to_string(slow) + ")");
    }
}

// --- Т6: точность решения скалярного уравнения -----------------------------
void test_scalar_solver()
{
    std::printf("T6: невязка g(delta) и граница упругости\n");
    TestLaw law = make_law_linear();
    law.plastic.type = fiber::HardeningType::Power;
    law.plastic.K = 300e6;
    law.plastic.n = 0.4;
    law.rate = {fiber::RateModel::CowperSymonds, 40.0, 5.0};

    double worst = 0.0;
    bool unique = true;
    for (double dt : {1e-7, 1e-5, 1e-3})
        for (double e : {0.06, 0.1, 0.4, 0.8})
        {
            fiber::BondState s0;
            s0.eps_p = 0.002;
            s0.eps_prev = e * 0.9;
            const fiber::BondState s1 = upd(law, s0, e, dt);
            const double delta = s1.eps_p - s0.eps_p;
            if (!(delta > 0.0)) unique = false;
            const double sig_tr = kE * (e - s0.eps_p);
            auto g = [&](double d)
            { return sig_tr - kE * d - fiber::plastic_yield_stress(law.plastic, law.rate, s0.eps_p + d, d / dt); };
            worst = std::max(worst, std::fabs(g(delta)) / sig_tr);
            if (!(g(delta * (1.0 + 1e-6) + 1e-15) < 0.0 && g(delta * (1.0 - 1e-6)) > 0.0)) unique = false;
        }
    check(worst < 1e-12, "|g(delta)|/sigma_tr < 1e-12 (CS + степенное; худшая " + std::to_string(worst) + ")");
    check(unique, "g меняет знак строго в корне");

    const fiber::BondState s0;
    const double e_y = kSy / kE;
    const fiber::BondState below = upd(law, s0, e_y * (1.0 - 1e-9), 1e-6);
    check(below.eps_p < 1e-300, "sigma_tr чуть ниже sy => delta == 0");
    const fiber::BondState above = upd(law, s0, e_y * (1.0 + 1e-9), 1e-6);
    check(above.eps_p >= 0.0 && above.eps_p < 1e-9,
          "sigma_tr чуть выше sy => delta малое (" + std::to_string(above.eps_p) + ")");
    const fiber::BondState la = upd(make_law_linear(), s0, e_y * (1.0 + 1e-9), 1e-6);
    check(la.eps_p > 0.0 && la.eps_p < 1e-9, "линейный: непрерывность на границе упругости");
}

// --- Т7: вязкоупругость -------------------------------------------------------
void test_viscoelastic()
{
    std::printf("T7: вязкоупругость (Максвелл + Кельвин)\n");
    TestLaw law;
    law.E = kE;
    law.visco.E1 = 1.0e9;
    law.visco.tau = 1.0e-5;

    // (a) Релаксация: быстрая ступень e0 (за 1e-9*tau), затем выдержка шагами разного dt.
    const double e0 = 0.01;
    for (double dt_hold : {1e-6, 1e-5, 3e-5})
    {
        fiber::BondState s;
        s = upd(law, s, e0, 1e-9 * law.visco.tau);
        const double q0 = s.q;
        check(near_rel(q0, law.visco.E1 * e0, 1e-8), "ступень: q(0) == E1*e0 (dt_hold=" + std::to_string(dt_hold) + ")");
        double t = 0.0, worst = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            s = upd(law, s, e0, dt_hold);
            t += dt_hold;
            worst = std::max(worst, std::fabs(s.q - q0 * std::exp(-t / law.visco.tau)) / q0);
        }
        check(worst < 1e-13, "релаксация q == q0*exp(-t/tau), ошибка не зависит от dt (" + std::to_string(worst) + ")");
        for (int k = 0; k < 400; ++k) s = upd(law, s, e0, dt_hold);
        check(near_rel(s.w_v, 0.5 * law.visco.E1 * e0 * e0, 1e-6), "w_v -> E1*e0^2/2 после полной релаксации");
    }

    // (b) Гармоническая деформация eps = em + ea*sin(wt): E', E'', tan(delta).
    law.visco.eta = 1.0e3; // Па*с
    const double E1 = law.visco.E1, tau = law.visco.tau, eta = law.visco.eta;
    for (double wt : {0.3, 1.0, 3.0}) // omega*tau
    {
        const double w = wt / tau;
        const double T = 2.0 * 3.14159265358979323846 / w;
        const int per = 4000;
        const double dt = T / per;
        const double em = 0.02, ea = 0.002;
        fiber::BondState s;
        s.eps_prev = em;
        double t = 0.0, sig = 0.0;
        auto step = [&]()
        {
            t += dt;
            const double e = em + ea * std::sin(w * t);
            const double ed = ea * w * std::cos(w * t);
            s = upd(law, s, e, dt);
            sig = stress(law, s, e, ed);
        };
        for (int k = 0; k < 20 * per; ++k) step(); // переходный процесс
        double S = 0.0, C = 0.0;
        for (int k = 0; k < per; ++k)
        {
            step();
            S += sig * std::sin(w * t);
            C += sig * std::cos(w * t);
        }
        S *= 2.0 / per / ea;
        C *= 2.0 / per / ea;
        const double w2 = wt * wt;
        const double Ep = kE + E1 * w2 / (1.0 + w2);
        const double Epp = E1 * wt / (1.0 + w2) + w * eta;
        const std::string tag = "w*tau=" + std::to_string(wt);
        check(near_rel(S, Ep, 2e-4), tag + ": E' = " + std::to_string(S) + " (ожид. " + std::to_string(Ep) + ")");
        check(near_rel(C, Epp, 2e-4), tag + ": E'' = " + std::to_string(C) + " (ожид. " + std::to_string(Epp) + ")");
        check(near_rel(C / S, Epp / Ep, 4e-4), tag + ": tan(delta) = E''/E'");
    }
}

// --- Т8: работа w_p при скоростной зависимости (установившееся течение) -------
void test_wp_rate()
{
    std::printf("T8: w_p при скоростной зависимости == int sigma_flow d eps_p\n");
    const TestLaw base = make_law_linear(0.0);
    const RateCase cases[] = {
        {"Cowper-Symonds", {fiber::RateModel::CowperSymonds, 40.0, 5.0}, ref_cs},
        {"Johnson-Cook", {fiber::RateModel::JohnsonCook, 0.05, 1e-3}, ref_jc},
        {"Perzyna", {fiber::RateModel::Perzyna, 1e6, 0.5}, ref_pz},
    };
    for (const RateCase& c : cases)
        for (double ed : {1.0, 100.0})
        {
            TestLaw law = base;
            law.rate = c.rate;
            const double de = 5e-4;
            const double dt = de / ed;
            fiber::BondState s, mid;
            double e = 0.0;
            for (int k = 0; k < 2000; ++k)
            {
                e += de;
                s = upd(law, s, e, dt);
                if (k == 999) mid = s;
            }
            const double ref = c.ref(c.rate, ed) * (s.eps_p - mid.eps_p); // sigma_flow = const
            check(near_rel(s.w_p - mid.w_p, ref, 1e-3),
                  std::string(c.name) + " ed=" + std::to_string(ed) + ": dw_p = " + std::to_string(s.w_p - mid.w_p) +
                      " (ожид. " + std::to_string(ref) + ")");
        }
}

// --- Т9: диссипация Кельвина и Максвелла (энергобаланс) -----------------------
void test_dissipation()
{
    std::printf("T9: w_v Кельвина и Максвелла\n");
    // Кельвин: eps_dot = const => w_v = eta*eps_dot^2*t, для любого dt
    {
        TestLaw law;
        law.E = kE;
        law.visco.eta = 2.0e3;
        for (double dt : {1e-7, 1e-5, 1e-3})
        {
            fiber::BondState s;
            s.eps_prev = 0.01;
            const double ed = 3.0;
            double e = 0.01, t = 0.0;
            for (int k = 0; k < 100; ++k)
            {
                e += ed * dt;
                t += dt;
                s = upd(law, s, e, dt);
            }
            check(near_rel(s.w_v, law.visco.eta * ed * ed * t, 1e-12), "Кельвин: w_v = eta*edot^2*t (dt=" + std::to_string(dt) + ")");
        }
    }
    // Максвелл: работа int q d eps = w_v + q^2/(2 E1) на шаге с линейным eps(t); dt/tau от 1e-3 (ветка малого r) до 20
    {
        TestLaw law;
        law.E = kE;
        law.visco.E1 = 1.0e9;
        law.visco.tau = 1.0e-5;
        const double tau = law.visco.tau, E1 = law.visco.E1;
        for (double r : {1e-3, 1e-2, 0.019, 0.021, 0.3, 1.0, 20.0})
        {
            const double dt = r * tau;
            fiber::BondState s;
            s.eps_prev = 0.02;
            double e = 0.02, work = 0.0;
            for (int k = 0; k < 400; ++k)
            {
                const double de = 1e-4 * std::sin(0.37 * k) + 2e-5; // знакопеременное приращение, eps > 0
                const double a = s.q, b = E1 * (de / dt) * tau, sx = std::exp(-r);
                const double int_q = a * tau * (1.0 - sx) + b * (dt - tau * (1.0 - sx));
                work += (de / dt) * int_q;
                e += de;
                s = upd(law, s, e, dt);
            }
            const double stored = 0.5 * s.q * s.q / E1;
            check(near_rel(work, s.w_v + stored, 1e-9) && s.w_v > 0.0,
                  "Максвелл: работа = w_v + q^2/2E1 (dt/tau=" + std::to_string(r) + ", W=" + std::to_string(work) + ")");
        }
    }
}

// --- Т10: провисшая связь -----------------------------------------------------
void test_slack()
{
    std::printf("T10: провисание (eps <= eps_p)\n");
    TestLaw law = make_law_linear();
    law.visco.E1 = 1.0e9;
    law.visco.tau = 1.0e-5;
    law.visco.eta = 1.0e3;
    const double dt = 1e-6;
    fiber::BondState s = run_path(law, ramp(0.0, 0.05, 100), dt);
    const double ep = s.eps_p;
    check(ep > 0.0 && s.q > 0.0, "после нагрузки eps_p > 0 и q > 0");

    s = upd(law, s, 0.5 * ep, dt); // пересечение границы провисания
    const double q0 = s.q;
    double wv0 = s.w_v;
    double t = 0.0;
    bool ok_q = true;
    for (double e : {0.4 * ep, -0.02, 0.9 * ep, 0.0, 0.7 * ep, ep})
    {
        s = upd(law, s, e, dt);
        t += dt;
        if (!near_rel(s.q, q0 * std::exp(-t / law.visco.tau), 1e-12)) ok_q = false;
        // в провисании энергия берётся только из q: w_v растёт на (q_old^2 - q_new^2)/(2 E1)
        const double q_old = q0 * std::exp(-(t - dt) / law.visco.tau);
        if (!near_rel(s.w_v - wv0, (q_old * q_old - s.q * s.q) / (2.0 * law.visco.E1), 1e-9)) ok_q = false;
        wv0 = s.w_v;
    }
    check(ok_q, "в провисании q только релаксирует экспонентой; w_v растёт лишь на релаксирующую часть q^2/2E1");

    // Без Максвелла (только Кельвин) w_v в провисании постоянно
    TestLaw lk = make_law_linear();
    lk.visco.eta = 1.0e3;
    fiber::BondState k = run_path(lk, ramp(0.0, 0.05, 100), dt);
    k = upd(lk, k, 0.5 * k.eps_p, dt);
    const double wk = k.w_v;
    for (double e : {0.3 * k.eps_p, -0.05, 0.99 * k.eps_p, k.eps_p}) k = upd(lk, k, e, dt);
    check(same(k.w_v, wk), "Кельвин: w_v не растёт при eps <= eps_p");

    // Выход из провисания: учитывается только натянутая часть (eps - eps_p), переход непрерывен
    TestLaw lm;
    lm.E = kE;
    lm.plastic = make_law_linear().plastic;
    lm.visco.E1 = 1.0e9;
    lm.visco.tau = 1.0e3; // почти без релаксации
    fiber::BondState r;
    r.eps_p = 0.01;
    r.eps_prev = 0.0; // провисла
    r = upd(lm, r, 0.012, 1e-6);
    check(near_rel(r.q, lm.visco.E1 * 0.002, 1e-6), "выход из провисания: q = E1*(eps - eps_p)");
    // Вход в провисание: q получает только натянутую часть
    fiber::BondState r2;
    r2.eps_p = 0.01;
    r2.eps_prev = 0.012;
    r2.q = 5.0e6;
    r2 = upd(lm, r2, 0.0, 1e-6);
    check(near_rel(r2.q, 5.0e6 + lm.visco.E1 * (0.01 - 0.012), 1e-6), "вход в провисание: Delta_eps_eff = eps_p - eps_prev");
}

// --- Т11: вырожденные dt и параметры ------------------------------------------
bool finite_state(const fiber::BondState& s)
{
    return std::isfinite(s.eps_p) && std::isfinite(s.q) && std::isfinite(s.eps_prev) && std::isfinite(s.w_p) &&
           std::isfinite(s.w_v);
}

void test_degenerate()
{
    std::printf("T11: вырожденные dt, E1 = 0, tau = 0\n");
    TestLaw law = make_law_linear();
    law.visco = {1.0e3, 0.0, 1.0e9, 1.0e-5};
    fiber::BondState s0 = run_path(law, ramp(0.0, 0.05, 50), 1e-6);
    for (double dt : {0.0, -1e-6})
    {
        const fiber::BondState s = upd(law, s0, 0.07, dt);
        check(same(s.eps_p, s0.eps_p) && same(s.q, s0.q) && same(s.eps_prev, s0.eps_prev) && same(s.w_p, s0.w_p) &&
                  same(s.w_v, s0.w_v),
              "dt <= 0 => состояние без изменений (dt=" + std::to_string(dt) + ")");
    }

    // dt = 1e-300 для разных законов
    TestLaw pz = make_law_linear(0.0);
    pz.rate = {fiber::RateModel::Perzyna, 1.0e6, 0.5};
    pz.visco = {1.0e3, 0.0, 1.0e9, 1.0e-5};
    TestLaw cs = make_law_linear();
    cs.rate = {fiber::RateModel::CowperSymonds, 40.0, 5.0};
    cs.visco = pz.visco;
    TestLaw lin = law;
    const TestLaw* all[] = {&pz, &cs, &lin};
    bool fin = true, sane = true;
    for (const TestLaw* l : all)
        for (double dt : {1e-300, 1e-200, 1e-30})
        {
            fiber::BondState s;
            s.eps_prev = 0.05;
            const fiber::BondState r = upd(*l, s, 0.05 + 1e-9, dt);
            if (!finite_state(r)) fin = false;
            if (!(r.eps_p >= 0.0 && r.eps_p <= 0.05 + 1e-9 && r.w_p >= 0.0 && r.w_v >= 0.0)) sane = false;
            // и «упругая» ступень из нуля
            fiber::BondState s2;
            const fiber::BondState r2 = upd(*l, s2, 1e-9, dt);
            if (!finite_state(r2)) fin = false;
        }
    check(fin, "dt = 1e-300 / 1e-200 / 1e-30: нет NaN/inf в состоянии");
    check(sane, "dt малый: eps_p в [0, eps], w_p, w_v >= 0");
    {
        fiber::BondState s;
        s.eps_prev = 0.05;
        const fiber::BondState r = upd(pz, s, 0.05 + 1e-9, 1e-300);
        // Perzyna m<1, dt -> 0: течение подавлено, работа не раздута (<= sigma_tr*delta)
        const double sig_tr = kE * (0.05 + 1e-9);
        check(r.eps_p < 1e-12 && r.w_p <= sig_tr * r.eps_p * (1.0 + 1e-12), "Perzyna m<1, dt=1e-300: течение подавлено, w_p не раздута");
    }

    // E1 = 0 и tau = 0: ветви Максвелла нет
    for (int variant = 0; variant < 2; ++variant)
    {
        TestLaw l;
        l.E = kE;
        l.visco.E1 = variant == 0 ? 0.0 : 1.0e9;
        l.visco.tau = variant == 0 ? 1.0e-5 : 0.0;
        fiber::BondState s = run_path(l, ramp(0.01, 0.05, 50), 1e-6, [] {
            fiber::BondState b;
            b.eps_prev = 0.01;
            return b;
        }());
        check(same(s.q, 0.0) && same(s.w_v, 0.0) && finite_state(s), variant == 0 ? "E1 = 0: q = 0, w_v = 0" : "tau = 0: ветви Максвелла нет");
    }

    // bond_stress: eta = 0 и бесконечная скорость -- не NaN
    TestLaw el;
    el.E = kE;
    fiber::BondState se;
    se.eps_prev = 0.01;
    check(same(stress(el, se, 0.01, std::numeric_limits<double>::infinity()), kE * 0.01), "eta = 0, eps_dot = inf => sigma = E*eps (не NaN)");
}

// --- Т12: табличное упрочнение: узлы и комбинация с RATE ----------------------
void test_table()
{
    std::printf("T12: табличное упрочнение: узлы, RATE\n");
    TestLaw law;
    law.E = kE;
    law.plastic.type = fiber::HardeningType::Table;
    law.plastic.table_eps = {0.0, 0.02, 0.1, 0.3};
    law.plastic.table_sigma = {60e6, 72e6, 95e6, 120e6};

    bool exact = true;
    for (std::size_t k = 0; k < law.plastic.table_eps.size(); ++k)
        if (!same(fiber::static_yield_stress(law.plastic, law.plastic.table_eps[k]), law.plastic.table_sigma[k])) exact = false;
    exact = exact && same(fiber::static_yield_stress(law.plastic, 5.0), 120e6) &&
            near_rel(fiber::static_yield_stress(law.plastic, 0.06), 83.5e6, 1e-14);
    check(exact, "sigma_s ровно в узлах таблицы, линейно между, константа после последнего");

    // Нагружение до eps, при котором eps_p ровно на узле: eps = eps_p_node + sigma_node/E
    bool on_node = true;
    for (std::size_t k = 1; k < law.plastic.table_eps.size(); ++k)
    {
        const double epn = law.plastic.table_eps[k];
        const double e = epn + law.plastic.table_sigma[k] / kE;
        fiber::BondState s = upd(law, fiber::BondState{}, e, 1e-6); // один шаг через несколько узлов
        fiber::BondState s2 = run_path(law, ramp(0.0, e, 300), 1e-6);
        if (!(std::fabs(s.eps_p - epn) < 1e-11 && std::fabs(s2.eps_p - epn) < 1e-11)) on_node = false;
    }
    check(on_node, "eps_p попадает ровно в узел (шаг через несколько узлов и мелкими шагами)");

    // Табличное + CS: невязка, монотонность по dt, предел dt -> inf
    TestLaw lr = law;
    lr.rate = {fiber::RateModel::CowperSymonds, 40.0, 5.0};
    bool resid = true, mono = true;
    for (double e : {0.03, 0.0805, 0.1203, 0.2, 0.45, 0.7}) // между узлами и за ними
    {
        double prev_ep = -1.0;
        for (double dt : {1e-7, 1e-5, 1e-3, 1e-1})
        {
            const fiber::BondState s1 = upd(lr, fiber::BondState{}, e, dt);
            const double sig_tr = kE * e;
            const double g = sig_tr - kE * s1.eps_p - fiber::plastic_yield_stress(lr.plastic, lr.rate, s1.eps_p, s1.eps_p / dt);
            if (std::fabs(g) / sig_tr > 1e-12) resid = false;
            if (s1.eps_p < prev_ep) mono = false; // чем больше dt, тем меньше скорость и выше eps_p
            prev_ep = s1.eps_p;
        }
        const fiber::BondState slow = upd(lr, fiber::BondState{}, e, 1e12);
        const fiber::BondState stat = upd(law, fiber::BondState{}, e, 1e-6);
        if (!near_rel(slow.eps_p, stat.eps_p, 5e-3)) resid = false;
    }
    check(resid, "table+CS: |g|/sigma_tr < 1e-12; dt -> inf даёт статический результат");
    check(mono, "table+CS: eps_p не убывает с ростом dt (скорость падает)");
}

// --- Т13: kelvin_eta_from_xi --------------------------------------------------
void test_xi()
{
    std::printf("T13: kelvin_eta_from_xi\n");
    const double xi = 0.05, E = 2.0e9, E1 = 1.0e9, A = 3.0e-6, l0 = 2.0e-3, m = 5.0e-7;
    const double E_inst = E + E1;
    const double eta = fiber::kelvin_eta_from_xi(xi, E_inst, A, l0, m);
    const double k = E_inst * A / l0;
    const double c = eta * A / l0;
    check(near_rel(c, 2.0 * xi * std::sqrt(k * m), 1e-14), "c = eta*A/l0 = 2*xi*sqrt(k*m_red)");
    check(near_rel(eta, 2.0 * xi * std::sqrt(E_inst * m * l0 / A), 1e-14), "eta = 2*xi*sqrt(E_inst*m*l0/A)");
    check(near_rel(c, 2.0 * xi * std::sqrt(k * m), 1e-14) && near_rel(c / (2.0 * std::sqrt(k * m)), xi, 1e-14),
          "доля критического демпфирования c/c_crit == xi");
    check(same(fiber::kelvin_eta_from_xi(0.0, E_inst, A, l0, m), 0.0), "xi = 0 => eta = 0");
    // «удвоение числа бусин»: l0/2, m/2 => c тот же, eta меняется как sqrt(m*l0)
    const double eta2 = fiber::kelvin_eta_from_xi(xi, E_inst, A, 0.5 * l0, 0.5 * m);
    check(near_rel(eta2 * A / (0.5 * l0), 2.0 * xi * std::sqrt(2.0 * k * 0.5 * m), 1e-14), "масштабирование c ~ sqrt(k*m) при дроблении");
}

} // namespace

int main()
{
    test_elastic_equivalence();
    test_monotonic_linear();
    test_unload_reload();
    test_perfect_plasticity();
    test_dt_independence();
    test_rate_dependence();
    test_scalar_solver();
    test_viscoelastic();
    test_wp_rate();
    test_dissipation();
    test_slack();
    test_degenerate();
    test_table();
    test_xi();

    if (failures == 0) std::printf("test_fiber_plasticity: OK\n");
    else std::printf("test_fiber_plasticity: %d failure(s)\n", failures);
    return failures;
}
