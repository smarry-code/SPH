// tests/test_fiber_damage.cpp
// Unit-тесты модели повреждённости связи с разупрочнением (DAMAGE; notes/fiber-model/fiber-damage-plan.md):
//   * чистые функции bond_constitutive_update / bond_stress (include/fiber/fiber_plasticity.hpp);
//   * deps_f связи (interactions/fiber_bond.hpp::fiber_bond_damage_deps_f);
//   * критерий разрыва D >= 1 (interactions/postprocessing.hpp::find_fiber_break_events);
//   * разбор директивы DAMAGE в .fib (fiber::load_fiber_file).
// Решатель не используется. Собирается с SPH_DIM=3 (нити).
#include "fiber/fiber_input.hpp"
#include "fiber/fiber_plasticity.hpp"
#include "fiber/fiber_system.hpp"
#include "fiber/model.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/postprocessing.hpp"
#include "param.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_damage требует сборки с SPH_DIM=3");

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

bool same_state(const fiber::BondState& a, const fiber::BondState& b)
{
    return same(a.eps_p, b.eps_p) && same(a.q, b.q) && same(a.eps_prev, b.eps_prev) && same(a.w_p, b.w_p) &&
           same(a.w_v, b.w_v) && same(a.d, b.d) && same(a.kappa, b.kappa) && same(a.eps0, b.eps0) &&
           same(a.w_d, b.w_d);
}

// Поля прежнего закона (без повреждённости): на шаге начала повреждения eps0/kappa уже заданы.
bool same_legacy(const fiber::BondState& a, const fiber::BondState& b)
{
    return same(a.eps_p, b.eps_p) && same(a.q, b.q) && same(a.eps_prev, b.eps_prev) && same(a.w_p, b.w_p) &&
           same(a.w_v, b.w_v) && same(a.d, b.d) && same(a.w_d, b.w_d);
}

// Параметры закона: владеющая обёртка (BondLaw хранит ссылки).
struct P
{
    double E = 80.0e9;
    fiber::PlasticParams plastic;
    fiber::RateParams rate;
    fiber::ViscoParams visco;
    fiber::FailureParams failure;
    fiber::DamageParams damage;
    double deps_f = 0.0;
    fiber::BondLaw law() const { return fiber::BondLaw(E, plastic, rate, visco, failure, damage, deps_f); }
};

constexpr double kSu = 2.5e9;
constexpr double kDepsF = 0.01;
constexpr double kDe = 2.0e-5; // шаг по деформации
constexpr double kDt = 1.0e-8;

P make_plastic_damaged()
{
    P p;
    p.plastic.type = fiber::HardeningType::Linear;
    p.plastic.sy = 2.04e9;
    p.plastic.H = 100.0e9;
    p.failure.sigma_u = kSu;
    p.damage.deps_f = kDepsF;
    p.deps_f = kDepsF;
    return p;
}

P without_damage(P p)
{
    p.damage = fiber::DamageParams{};
    p.deps_f = 0.0;
    return p;
}

fiber::BondState upd(const P& p, const fiber::BondState& s, double e, double dt = kDt)
{
    return fiber::bond_constitutive_update(p.law(), s, e, dt);
}

double stress(const P& p, const fiber::BondState& s, double e, double edot = 0.0)
{
    return fiber::bond_stress(p.law(), s, e, edot);
}

double envelope(double e, double eps0) { return kSu * std::max(0.0, 1.0 - (e - eps0) / kDepsF); }

// ---------------------------------------------------------------------------
// A1: монотонное растяжение с пластикой
// ---------------------------------------------------------------------------
void test_a1_monotone()
{
    std::printf("A1: монотонное растяжение, пластика + разупрочнение\n");
    const P pd = make_plastic_damaged();
    const P pn = without_damage(pd);
    fiber::BondState s, sn;
    int pre_mismatch = 0, post_steps = 0;
    double max_err = 0.0;
    double eps0 = -1.0, ep_onset = 0.0, sig_eff_before = 0.0, sig_eff_at = 0.0;
    bool onset_cap_ok = false;
    bool ep_frozen = true, d_ok_before = true, d_ok_after = true, zero_after = true, d_monotone = true;
    double d_prev = 0.0;
    const int n = 2600;
    for (int k = 1; k <= n; ++k)
    {
        const double e = kDe * static_cast<double>(k);
        const bool started_before = s.eps0 >= 0.0;
        const double sig_eff_prev = pd.E * ((e - kDe) - s.eps_p);
        s = upd(pd, s, e);
        sn = upd(pn, sn, e);
        if (!started_before)
        {
            if (s.eps0 < 0.0 && (!same_legacy(s, sn) || !same(stress(pd, s, e), stress(pn, sn, e)))) ++pre_mismatch;
            if (s.eps0 >= 0.0)
            {
                // На шаге начала eps0 -- точка, где sigma_eff ровно = sigma_u (внутри шага), перелёт
                // уже идёт в разупрочнение: sigma = sigma_u(1 - (e - eps0)/deps_f) <= sigma_u,
                // и eps0 = eps_p + sigma_u/E (q = 0).
                onset_cap_ok = near_rel(stress(pd, s, e), kSu * (1.0 - (e - s.eps0) / kDepsF), 1e-12) &&
                               near_rel(s.eps0, s.eps_p + kSu / pd.E, 1e-12) && stress(pd, s, e) <= kSu;
                eps0 = s.eps0;
                ep_onset = s.eps_p;
                sig_eff_before = sig_eff_prev;
                sig_eff_at = pd.E * (e - s.eps_p);
            }
        }
        else
        {
            ++post_steps;
            const double sg = stress(pd, s, e);
            max_err = std::max(max_err, std::fabs(sg - envelope(e, eps0)) / kSu);
            if (!same(s.eps_p, ep_onset)) ep_frozen = false;
            if (e >= eps0 + kDepsF * (1.0 + 1e-9))
            {
                if (!same(s.d, 1.0)) d_ok_after = false;
                if (!same(sg, 0.0)) zero_after = false;
            }
            if (e < eps0 + kDepsF * (1.0 - 1e-9) && !(s.d < 1.0)) d_ok_before = false;
        }
        if (s.d < d_prev) d_monotone = false;
        d_prev = s.d;
    }
    check(pre_mismatch == 0, "до начала повреждения состояние и напряжение побитово равны закону без DAMAGE");
    check(eps0 >= 0.0 && ep_onset > 0.0, "повреждение началось после пластического течения (eps_p > 0)");
    check(sig_eff_before < kSu && sig_eff_at >= kSu, "начало: E(eps-eps_p) пересекает sigma_u на этом шаге");
    check(sig_eff_at - kSu < 0.1 * kSu, "перескок sigma_eff через sigma_u не превышает шаг дискретизации");
    check(onset_cap_ok, "на шаге начала: eps0 = eps_p + sigma_u/E (точка пересечения), sigma на огибающей и <= sigma_u");
    check(max_err < 1e-9, "после начала sigma = sigma_u(1-(eps-eps0)/deps_f), невязка/sigma_u < 1e-9 (" +
                              std::to_string(max_err) + ")");
    check(ep_frozen, "eps_p не меняется после начала повреждения (побитово)");
    check(d_ok_before, "d < 1 до eps0 + deps_f");
    check(d_ok_after, "d == 1 при eps >= eps0 + deps_f");
    check(zero_after, "sigma == 0 после полного разрушения");
    check(d_monotone, "d не убывает");
    check(post_steps > 500, "пройден весь участок разупрочнения");
}

// ---------------------------------------------------------------------------
// A2: разгрузка / повторная нагрузка
// ---------------------------------------------------------------------------
void test_a2_unload_reload()
{
    std::printf("A2: разгрузка и повторная нагрузка на участке разупрочнения\n");
    const P p = make_plastic_damaged();
    fiber::BondState s;
    int k = 0;
    double e = 0.0;
    while (s.eps0 < 0.0 || e < s.eps0 + 0.8 * kDepsF)
    {
        ++k;
        e = kDe * static_cast<double>(k);
        s = upd(p, s, e);
    }
    const double eps0 = s.eps0;
    const double d_top = s.d, kappa_top = s.kappa, e_top = e;
    check(d_top > 0.5 && d_top < 1.0, "d в разупрочнении: " + std::to_string(d_top));
    check(near_rel(stress(p, s, e), envelope(e, eps0), 1e-9), "в вершине напряжение на огибающей");

    // Разгрузка до половины (kappa - eps_p)
    const double e_low = s.eps_p + 0.5 * (kappa_top - s.eps_p);
    bool d_const = true, kappa_const = true, secant = true;
    const int nu = 200;
    for (int i = 1; i <= nu; ++i)
    {
        const double ei = e_top + (e_low - e_top) * static_cast<double>(i) / nu;
        s = upd(p, s, ei);
        if (!same(s.d, d_top)) d_const = false;
        if (!same(s.kappa, kappa_top)) kappa_const = false;
        if (!near_rel(stress(p, s, ei), (1.0 - s.d) * p.E * (ei - s.eps_p), 1e-12)) secant = false;
    }
    check(d_const, "при разгрузке d не меняется");
    check(kappa_const, "при разгрузке kappa не меняется");
    check(secant, "разгрузка по секущей: sigma = (1-d) E (eps - eps_p)");
    check(stress(p, s, e_low) < envelope(e_top, eps0), "разгруженное напряжение ниже огибающей в вершине");

    // Повторная нагрузка: до kappa по секущей, дальше по огибающей
    bool d_never_down = true, sec_up = true, env_again = true;
    double d_prev = s.d;
    const int nr = 400;
    const double e_end = e_top + 0.1 * kDepsF;
    for (int i = 1; i <= nr; ++i)
    {
        const double ei = e_low + (e_end - e_low) * static_cast<double>(i) / nr;
        s = upd(p, s, ei);
        if (s.d < d_prev) d_never_down = false;
        d_prev = s.d;
        const double sg = stress(p, s, ei);
        if (ei <= kappa_top)
        {
            if (!same(s.d, d_top) || !near_rel(sg, (1.0 - s.d) * p.E * (ei - s.eps_p), 1e-12)) sec_up = false;
        }
        else if (std::fabs(sg - envelope(ei, eps0)) > 1e-9 * kSu)
        {
            env_again = false;
        }
    }
    check(sec_up, "повторная нагрузка до kappa: секущая, d прежнее");
    check(env_again, "за kappa напряжение снова на огибающей");
    check(d_never_down, "d не убывает при повторной нагрузке");
    check(s.d > d_top, "за kappa d растёт");
}

// ---------------------------------------------------------------------------
// A3: энергия
// ---------------------------------------------------------------------------
void test_a3_energy()
{
    std::printf("A3: энергетический баланс при полном разрушении\n");
    const P p = make_plastic_damaged();
    fiber::BondState s;
    double W = 0.0, sig_prev = 0.0;
    double w_p_onset = 0.0, ee_onset = 0.0;
    const int n = 2800;
    for (int k = 1; k <= n; ++k)
    {
        const double e = kDe * static_cast<double>(k);
        const bool started_before = s.eps0 >= 0.0;
        s = upd(p, s, e);
        if (!started_before && s.eps0 >= 0.0)
        {
            w_p_onset = s.w_p;
            ee_onset = s.eps0 - s.eps_p;
        }
        const double sg = stress(p, s, e);
        W += 0.5 * (sig_prev + sg) * kDe;
        sig_prev = sg;
    }
    check(same(s.d, 1.0), "связь полностью разрушена");
    check(s.w_d > 0.0 && s.w_p > 0.0, "w_d > 0 и w_p > 0");
    check(same(s.w_p, w_p_onset), "w_p не растёт после начала повреждения");
    const double diss = s.w_p + s.w_d;
    std::printf("    W = %.6e, w_p + w_d = %.6e (w_p = %.6e, w_d = %.6e)\n", W, diss, s.w_p, s.w_d);
    check(near_rel(diss, W, 5e-3), "внешняя работа == w_p + w_d с точностью 0.5 %");
    const double w_d_ref = 0.5 * kSu * ee_onset + 0.5 * kSu * kDepsF;
    check(near_rel(s.w_d, w_d_ref, 1e-2), "w_d == 0.5 sigma_u (eps0-eps_p) + 0.5 sigma_u deps_f с точностью 1 %");
}

// ---------------------------------------------------------------------------
// A4: энергия разрушения GF и независимость от длины связи
// ---------------------------------------------------------------------------
double softening_energy_per_bond(double gf, double l0, double A, double& deps_f_out, bool& broke)
{
    fiber::FiberMaterial m;
    m.E = 80.0e9;
    m.failure.sigma_u = kSu;
    m.damage.gf = gf;
    deps_f_out = interactions::fiber_bond_damage_deps_f(m, l0);
    P p;
    p.E = m.E;
    p.failure = m.failure;
    p.damage = m.damage;
    p.deps_f = deps_f_out;
    fiber::BondState s;
    double ee_onset = 0.0;
    for (int k = 1; k <= 4000; ++k)
    {
        const double e = kDe * static_cast<double>(k);
        const bool started_before = s.eps0 >= 0.0;
        s = upd(p, s, e);
        if (!started_before && s.eps0 >= 0.0) ee_onset = s.eps0 - s.eps_p;
    }
    broke = same(s.d, 1.0);
    return (s.w_d - 0.5 * p.E * ee_onset * ee_onset) * A * l0;
}

void test_a4_gf()
{
    std::printf("A4: GF, независимость энергии разрушения от длины связи\n");
    const double l0 = 0.01, A = 1.0e-6;
    const double gf = 0.5 * kSu * kDepsF * l0; // даёт deps_f = 0.01 при l0
    fiber::FiberMaterial m;
    m.failure.sigma_u = kSu;
    m.damage.gf = gf;
    check(near_rel(interactions::fiber_bond_damage_deps_f(m, l0), 2.0 * gf / (kSu * l0), 1e-14),
          "deps_f = 2 GF / (sigma_u l0)");
    check(near_rel(interactions::fiber_bond_damage_deps_f(m, 0.5 * l0), 4.0 * gf / (kSu * l0), 1e-14),
          "для l0/2 deps_f вдвое больше");
    fiber::FiberMaterial m2;
    m2.failure.sigma_u = kSu;
    m2.damage.deps_f = 0.03;
    check(same(interactions::fiber_bond_damage_deps_f(m2, l0), 0.03) &&
              same(interactions::fiber_bond_damage_deps_f(m2, 0.1 * l0), 0.03),
          "DEPS_F задан: берётся как есть, не зависит от l0");
    fiber::FiberMaterial m3;
    m3.failure.sigma_u = kSu;
    check(same(interactions::fiber_bond_damage_deps_f(m3, l0), 0.0), "DAMAGE выключен -> deps_f = 0");
    fiber::FiberMaterial m4;
    m4.damage.gf = gf;
    check(same(interactions::fiber_bond_damage_deps_f(m4, l0), 0.0), "без SIGMA_U -> deps_f = 0");

    double df1 = 0.0, df2 = 0.0;
    bool b1 = false, b2 = false;
    const double g1 = softening_energy_per_bond(gf, l0, A, df1, b1);
    const double g2 = softening_energy_per_bond(gf, 0.5 * l0, A, df2, b2);
    std::printf("    deps_f: %.6g / %.6g, энергия: %.6e / %.6e, GF*A = %.6e\n", df1, df2, g1, g2, gf * A);
    check(b1 && b2, "обе связи разрушены");
    check(near_rel(g1, g2, 1e-2), "энергия разупрочнения связи одинакова для l0 и l0/2 (1 %)");
    check(near_rel(g1, gf * A, 1e-2), "энергия разупрочнения == GF * A (1 %)");
}

// ---------------------------------------------------------------------------
// A5: критерий разрыва
// ---------------------------------------------------------------------------
struct OneBond
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

OneBond make_one_bond(bool damage, double length_factor)
{
    OneBond ob;
    const double l0 = 0.01;
    fiber::FiberMaterial m;
    m.E = 80.0e9;
    m.failure.sigma_u = kSu;
    if (damage) m.damage.deps_f = kDepsF;
    ob.materials.push_back(m);
    ob.state.resize(2);
    ob.topology.resize(2);
    for (std::size_t d = 0; d < dim; ++d) ob.state.x[d][0] = ob.state.x[d][1] = 0.0;
    ob.state.x[0][1] = length_factor * l0;
    ob.state.mass[0] = ob.state.mass[1] = 1.0e-6;
    ob.state.object_id[0] = ob.state.object_id[1] = 0;
    ob.state.material_id[0] = ob.state.material_id[1] = 0;
    ob.topology.next[0] = 1;
    ob.topology.prev[1] = 0;
    fiber::FiberBond b;
    b.i = 0;
    b.j = 1;
    b.l0 = l0;
    b.l_break = 1.5 * l0;
    ob.topology.bonds.push_back(b);
    return ob;
}

std::vector<interactions::FiberBreakEvent> find_events(OneBond& ob)
{
    rhs::CpuExecutor exec(1);
    std::vector<std::vector<interactions::FiberBreakEvent>> buf;
    return interactions::find_fiber_break_events(exec, ob.state, ob.topology, ob.materials, buf);
}

void test_a5_break()
{
    std::printf("A5: разрыв по повреждённости\n");
    {
        OneBond ob = make_one_bond(true, 1.04); // eps = 0.04, E*eps = 3.2e9 > sigma_u
        ob.topology.bonds[0].st.d = 1.0 - 1e-6;
        check(find_events(ob).empty(), "d чуть меньше 1 -> события нет");
        ob.topology.bonds[0].st.d = 1.0;
        const auto ev = find_events(ob);
        check(ev.size() == 1 && ev[0].reason == interactions::FiberBreakReason::Damage, "d == 1 -> разрыв, причина Damage");
    }
    {
        OneBond ob = make_one_bond(true, 1.04);
        ob.topology.bonds[0].st.d = 0.0;
        check(find_events(ob).empty(), "DAMAGE включён, sigma > sigma_u, d < 1 -> события нет (SIGMA_U-критерий отключён)");
    }
    {
        OneBond ob = make_one_bond(false, 1.04);
        const auto ev = find_events(ob);
        check(ev.size() == 1 && ev[0].reason == interactions::FiberBreakReason::Stress,
              "без DAMAGE sigma > sigma_u -> разрыв, причина Stress (прежнее поведение)");
    }
    {
        OneBond ob = make_one_bond(true, 1.6); // l > l_break
        ob.topology.bonds[0].st.d = 0.0;
        const auto ev = find_events(ob);
        check(ev.size() == 1 && ev[0].reason == interactions::FiberBreakReason::Strain,
              "l > l_break при включённом DAMAGE -> разрыв, причина Strain");
    }
}

// ---------------------------------------------------------------------------
// A6: ветви Максвелла и Кельвина
// ---------------------------------------------------------------------------
// ОГОВОРКА: точный энергетический баланс здесь не проверяется. Повреждение
// масштабирует Максвелл-ветвь вместе с упругой частью ((1-D) sigma), тогда как
// диссипация w_d считается только по упругой пружине; поэтому проверяются
// качественные свойства (монотонность d, неотрицательность sigma, d -> 1,
// sigma -> 0, w_d > 0), а соотношение работы и диссипации только печатается.
void test_a6_visco()
{
    std::printf("A6: Максвелл + Кельвин\n");
    P p = make_plastic_damaged();
    p.visco.E1 = 0.5e9;
    p.visco.tau = 1.0e-5;
    p.visco.eta = 1.0e3;
    fiber::BondState s;
    double W = 0.0, sig_prev = 0.0, d_prev = 0.0;
    bool d_mono = true, nonneg = true, finite = true;
    const double edot = kDe / kDt;
    for (int k = 1; k <= 2800; ++k)
    {
        const double e = kDe * static_cast<double>(k);
        s = upd(p, s, e);
        const double sg = stress(p, s, e, edot);
        if (s.d < d_prev) d_mono = false;
        d_prev = s.d;
        if (sg < 0.0) nonneg = false;
        if (!std::isfinite(sg) || !std::isfinite(s.w_d) || !std::isfinite(s.w_v)) finite = false;
        W += 0.5 * (sig_prev + sg) * kDe;
        sig_prev = sg;
    }
    std::printf("    W = %.6e, w_p + w_d + w_v = %.6e (w_p %.4e, w_d %.4e, w_v %.4e)\n", W, s.w_p + s.w_d + s.w_v,
                s.w_p, s.w_d, s.w_v);
    check(d_mono, "d не убывает");
    check(nonneg, "sigma никогда не отрицательна");
    check(finite, "все величины конечны");
    check(same(s.d, 1.0), "d достигает 1");
    check(same(stress(p, s, 2800 * kDe, edot), 0.0), "в конце sigma == 0");
    check(s.w_d > 0.0, "w_d > 0");
}

// ---------------------------------------------------------------------------
// A8: разбор .fib
// ---------------------------------------------------------------------------
const char* kPath = "__fiber_damage_fixture.fib";
const char* kBase = "MATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.30\nFIBER 0 3\n0 0 0\n0.01 0 0\n0.02 0 0\n";

fiber::Model load(const std::string& content)
{
    {
        std::ofstream out(kPath);
        out << content;
    }
    fiber::Model m;
    fiber::load_fiber_file(kPath, m);
    return m;
}

void expect_error(const std::string& name, const std::string& content)
{
    try
    {
        load(content);
        check(false, name + ": ожидалось исключение, его нет");
    }
    catch (const std::exception& e)
    {
        check(true, name + ": '" + e.what() + "'");
    }
}

void test_a8_parsing()
{
    std::printf("A8: директива DAMAGE в .fib\n");
    try
    {
        const fiber::Model m = load(std::string("FAILURE 0 SIGMA_U=90e6\nDAMAGE 0 DEPS_F=0.02\n") + kBase);
        const fiber::FiberMaterial& mat = m.materials[0];
        check(near_rel(mat.damage.deps_f, 0.02, 1e-14) && same(mat.damage.gf, 0.0), "DEPS_F принят и сохранён");
        check(mat.damage.enabled() && mat.inelastic(), "damage.enabled() и inelastic()");
        check(near_rel(mat.failure.sigma_u, 90e6, 1e-14), "SIGMA_U не затронут");
    }
    catch (const std::exception& e)
    {
        check(false, std::string("DEPS_F: неожиданное исключение: ") + e.what());
    }
    try
    {
        const fiber::Model m = load(kBase + std::string("DAMAGE 0 gf=150.5\nFAILURE 0 SIGMA_U=90e6\n"));
        check(near_rel(m.materials[0].damage.gf, 150.5, 1e-14) && same(m.materials[0].damage.deps_f, 0.0),
              "GF принят (порядок директив и регистр свободны, DAMAGE до FAILURE)");
    }
    catch (const std::exception& e)
    {
        check(false, std::string("GF: неожиданное исключение: ") + e.what());
    }
    try
    {
        const fiber::Model m = load(kBase);
        check(!m.materials[0].damage.enabled(), "без директивы DAMAGE выключен");
    }
    catch (const std::exception& e)
    {
        check(false, std::string("без DAMAGE: неожиданное исключение: ") + e.what());
    }
    const std::string fail = "FAILURE 0 SIGMA_U=90e6\n";
    expect_error("оба ключа", fail + "DAMAGE 0 DEPS_F=0.02 GF=100\n" + kBase);
    expect_error("нет ключа", fail + "DAMAGE 0\n" + kBase);
    expect_error("DEPS_F = 0", fail + "DAMAGE 0 DEPS_F=0\n" + kBase);
    expect_error("DEPS_F < 0", fail + "DAMAGE 0 DEPS_F=-0.01\n" + kBase);
    expect_error("GF = 0", fail + "DAMAGE 0 GF=0\n" + kBase);
    expect_error("GF < 0", fail + "DAMAGE 0 GF=-5\n" + kBase);
    expect_error("THETA_C", fail + "DAMAGE 0 THETA_C=0.5\n" + kBase);
    expect_error("THETA_C вместе с DEPS_F", fail + "DAMAGE 0 DEPS_F=0.02 THETA_C=0.5\n" + kBase);
    expect_error("без FAILURE", std::string("DAMAGE 0 DEPS_F=0.02\n") + kBase);
    expect_error("FAILURE без SIGMA_U", std::string("FAILURE 0 EPS_P=0.2\nDAMAGE 0 DEPS_F=0.02\n") + kBase);
    expect_error("повтор DAMAGE", fail + "DAMAGE 0 DEPS_F=0.02\nDAMAGE 0 GF=10\n" + kBase);
    expect_error("нечисловое значение", fail + "DAMAGE 0 DEPS_F=abc\n" + kBase);
    std::remove(kPath);
}

// ---------------------------------------------------------------------------
// A9: сжатие / провисание после начала повреждения
// ---------------------------------------------------------------------------
void test_a9_slack()
{
    std::printf("A9: сжатие после начала повреждения\n");
    const P p = make_plastic_damaged();
    fiber::BondState s;
    int k = 0;
    double e = 0.0;
    while (s.eps0 < 0.0 || e < s.eps0 + 0.3 * kDepsF)
    {
        ++k;
        e = kDe * static_cast<double>(k);
        s = upd(p, s, e);
    }
    const double d0 = s.d, kap0 = s.kappa, ep0 = s.eps_p;
    check(d0 > 0.0 && d0 < 1.0 && ep0 > 0.0, "исходно d в (0,1), eps_p > 0");
    bool zero = true, d_const = true;
    const double targets[] = {0.5 * ep0, ep0, 0.0, -0.01};
    for (double et : targets)
    {
        s = upd(p, s, et);
        if (!same(stress(p, s, et, -1.0e3), 0.0)) zero = false;
        if (!same(s.d, d0) || !same(s.kappa, kap0) || !same(s.eps_p, ep0)) d_const = false;
    }
    check(zero, "при eps <= eps_p sigma == 0");
    check(d_const, "d, kappa, eps_p не меняются при сжатии");
}

// ---------------------------------------------------------------------------
// A10: детерминизм числа потоков
// ---------------------------------------------------------------------------
void test_a10_determinism()
{
    std::printf("A10: детерминизм find_fiber_break_events (1 и 4 потока)\n");
    const std::size_t nb = 60, np = nb + 1;
    const double l0 = 0.01;
    fiber::FiberState state;
    fiber::FiberTopology topology;
    state.resize(np);
    topology.resize(np);
    fiber::FiberMaterial m;
    m.E = 80.0e9;
    m.failure.sigma_u = kSu;
    m.damage.deps_f = kDepsF;
    const std::vector<fiber::FiberMaterial> materials{m};
    for (std::size_t k = 0; k < np; ++k)
    {
        for (std::size_t d = 0; d < dim; ++d) state.x[d][k] = 0.0;
        state.x[0][k] = l0 * static_cast<double>(k);
        state.mass[k] = 1.0e-6;
        state.object_id[k] = 0;
        state.material_id[k] = 0;
        if (k + 1 < np) topology.next[k] = static_cast<long long>(k + 1);
        if (k > 0) topology.prev[k] = static_cast<long long>(k - 1);
    }
    for (std::size_t k = 0; k < nb; ++k)
    {
        fiber::FiberBond b;
        b.i = k;
        b.j = k + 1;
        b.l0 = l0;
        b.l_break = (k % 7 == 0) ? 0.5 * l0 : 1.5 * l0; // Strain
        if (k % 3 == 0) b.st.d = 1.0;                    // Damage
        else if (k % 5 == 0) b.st.d = 0.5;               // не рвётся
        topology.bonds.push_back(b);
    }
    rhs::CpuExecutor ex1(1), ex4(4);
    std::vector<std::vector<interactions::FiberBreakEvent>> buf1, buf4;
    const auto e1 = interactions::find_fiber_break_events(ex1, state, topology, materials, buf1);
    const auto e4 = interactions::find_fiber_break_events(ex4, state, topology, materials, buf4);
    bool eq = e1.size() == e4.size();
    for (std::size_t i = 0; eq && i < e1.size(); ++i)
        eq = e1[i].bond_index == e4[i].bond_index && e1[i].reason == e4[i].reason;
    check(!e1.empty(), "события есть (" + std::to_string(e1.size()) + ")");
    check(eq, "списки событий одинаковы для 1 и 4 потоков");
    bool sorted = true;
    for (std::size_t i = 1; i < e4.size(); ++i)
        if (e4[i - 1].bond_index >= e4[i].bond_index) sorted = false;
    check(sorted, "события упорядочены по bond_index");
    // Ожидаемое число: k%3==0 или k%7==0
    std::size_t expected = 0;
    for (std::size_t k = 0; k < nb; ++k)
        if (k % 3 == 0 || k % 7 == 0) ++expected;
    check(e4.size() == expected, "число событий совпадает с ожидаемым (" + std::to_string(expected) + ")");
}

// ---------------------------------------------------------------------------
// R1: регрессия, выключенный DAMAGE побитово не меняет закон
// ---------------------------------------------------------------------------
void test_r1_regression()
{
    std::printf("R1: выключенный DAMAGE побитово равен прежнему закону\n");
    P p;
    p.plastic.type = fiber::HardeningType::Linear;
    p.plastic.sy = 2.04e9;
    p.plastic.H = 100.0e9;
    p.failure.sigma_u = kSu;
    p.visco.E1 = 0.5e9;
    p.visco.tau = 1.0e-5;
    p.visco.eta = 1.0e3;
    const fiber::DamageParams off{};
    // (a) прежний вызов (без damage), (b) damage по умолчанию, (c) damage включён, но SIGMA_U нет
    P pc = p;
    pc.failure = fiber::FailureParams{};
    pc.damage.deps_f = kDepsF;
    pc.deps_f = kDepsF;
    P pc_ref = pc;
    pc_ref.damage = fiber::DamageParams{};
    pc_ref.deps_f = 0.0;

    fiber::BondState sa, sb, sc, sr;
    bool states_eq = true, stress_eq = true, c_eq = true;
    // путь: нагрузка с течением, разгрузка, сжатие, повторная нагрузка выше
    std::vector<double> path;
    for (int i = 1; i <= 1500; ++i) path.push_back(kDe * static_cast<double>(i));
    for (int i = 1; i <= 500; ++i) path.push_back(kDe * 1500.0 - kDe * static_cast<double>(i));
    for (int i = 1; i <= 200; ++i) path.push_back(kDe * 1000.0 - 1.0e-4 * static_cast<double>(i));
    for (int i = 1; i <= 2000; ++i) path.push_back(-0.02 + kDe * static_cast<double>(i));
    double e_prev = 0.0;
    for (double e : path)
    {
        const double edot = (e - e_prev) / kDt;
        e_prev = e;
        sa = fiber::bond_constitutive_update(fiber::BondLaw(p.E, p.plastic, p.rate, p.visco, p.failure), sa, e, kDt);
        sb = fiber::bond_constitutive_update(
            fiber::BondLaw(p.E, p.plastic, p.rate, p.visco, p.failure, off, 0.0), sb, e, kDt);
        const double ta = fiber::bond_stress(fiber::BondLaw(p.E, p.plastic, p.rate, p.visco, p.failure), sa, e, edot);
        const double tb = fiber::bond_stress(
            fiber::BondLaw(p.E, p.plastic, p.rate, p.visco, p.failure, off, 0.0), sb, e, edot);
        if (!same_state(sa, sb)) states_eq = false;
        if (!same(ta, tb)) stress_eq = false;
        sc = upd(pc, sc, e);
        sr = upd(pc_ref, sr, e);
        if (!same_state(sc, sr) || !same(stress(pc, sc, e, edot), stress(pc_ref, sr, e, edot))) c_eq = false;
    }
    check(states_eq, "BondState: старый вызов == damage по умолчанию на всём пути");
    check(stress_eq, "bond_stress: побитово равен");
    check(sa.w_p > 0.0 && sa.w_v > 0.0 && sa.eps0 < 0.0 && same(sa.d, 0.0), "путь содержал пластику и вязкость, d не тронуто");
    check(c_eq, "DAMAGE включён, но SIGMA_U = 0 (softening() == false) -> закон прежний");
}

} // namespace

int main()
{
    test_a1_monotone();
    test_a2_unload_reload();
    test_a3_energy();
    test_a4_gf();
    test_a5_break();
    test_a6_visco();
    test_a8_parsing();
    test_a9_slack();
    test_a10_determinism();
    test_r1_regression();

    if (failures == 0)
        std::printf("test_fiber_damage: OK\n");
    else
        std::printf("test_fiber_damage: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
