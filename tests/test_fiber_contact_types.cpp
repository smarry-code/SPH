// tests/test_fiber_contact_types.cpp
// Законы нормальной силы pinball / Swegle / импульсный для Fiber-Fiber-Contact
// (fiber/plan-fiber-contact-types.md, этап 2, шаг 2.4;
// include/interactions/fiber_contact_law.hpp).
//   1. тип 0 с явными параметрами == вызов без параметров (побитово);
//   2. аналитика: крест двух сегментов, сила ядра == contact::penalty_force;
//   3. знак pdot, симметрия A<->B, ноль в покое и при разлёте;
//   4. пересечение осей D = 0: Swegle даёт fn > 0 (g.r = r_a + r_b);
//   5. баланс импульса, s, t in {0, 0.5, 1}, точечные сегменты;
//   6. согласованность ядра, предела шага и диагностики;
//   7. трение: slot.N == fn, при fn = 0 слот неактивен.
// e2e-метрики (пункт 8) -- в tests/test_e2e_fiber_contact.cpp.
#include "contact_force.hpp"
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_contact_law.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_contact_types требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1e-300, std::fabs(b)); }

const char* type_name(int t)
{
    switch (t)
    {
        case 0: return "Hertz";
        case 1: return "pinball";
        case 2: return "Swegle";
        case 3: return "impulse";
        case 4: return "Hertz+damper";
        default: return "?";
    }
}

constexpr double kE = 2.0e9, kNu = 0.35, kD = 1.0e-3, kRho = 1100.0;
constexpr double kL = 2.0e-3;

// Крест: сегмент A вдоль x на высоте +zA, сегмент B вдоль y на высоте -zA;
// ближайшие точки -- s на A, t на B. A движется по z со скоростью vA, B -- vB.
// point_a / point_b: сторона -- одна изолированная бусина.
struct Cross
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    std::vector<double> bead_dt_ref;
    interactions::FiberContactPrimitive prim;
    std::size_t na = 0; // число бусин стороны A (бусины A: [0, na), B: [na, n))
};

Cross make_cross(double s, double t, double D, double vA, double vB, bool point_a = false, bool point_b = false,
                 double mu_ff = 0.0)
{
    Cross c;
    fiber::FiberMaterial m;
    m.E = kE; m.nu = kNu; m.d = kD; m.rho = kRho; m.mu_ff = mu_ff;
    c.materials.push_back(m);
    const std::size_t na = point_a ? 1 : 2, nb = point_b ? 1 : 2;
    c.na = na;
    const std::size_t n = na + nb;
    c.state.resize(n);
    c.topology.resize(n);
    const double A = interactions::fiber_cross_section_area(m);
    const double m0 = kRho * A * kL;
    const double xc = s * kL; // точка пересечения по x
    const double yc = t * kL; // точка пересечения по y относительно b0
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d) c.state.vx[d][i] = 0.0;
        c.state.material_id[i] = 0;
    }
    // A: a0 = (0, 0, +D/2), a1 = (L, 0, +D/2); при point_a -- бусина в точке пересечения.
    c.state.x[0][0] = point_a ? xc : 0.0; c.state.x[1][0] = 0.0; c.state.x[2][0] = 0.5 * D;
    c.state.mass[0] = m0; c.state.object_id[0] = 0; c.state.vx[2][0] = vA;
    if (!point_a)
    {
        c.state.x[0][1] = kL; c.state.x[1][1] = 0.0; c.state.x[2][1] = 0.5 * D;
        c.state.mass[1] = 1.6 * m0; c.state.object_id[1] = 0; c.state.vx[2][1] = vA;
        c.topology.next[0] = 1; c.topology.prev[1] = 0;
        fiber::FiberBond b; b.i = 0; b.j = 1; b.l0 = kL; b.l_break = 1e9;
        c.topology.bonds.push_back(b);
    }
    // B: b0 = (xc, -yc, -D/2), b1 = (xc, L - yc, -D/2).
    const std::size_t b0 = na;
    c.state.x[0][b0] = xc; c.state.x[1][b0] = point_b ? 0.0 : -yc; c.state.x[2][b0] = -0.5 * D;
    c.state.mass[b0] = 0.8 * m0; c.state.object_id[b0] = 1; c.state.vx[2][b0] = vB;
    if (!point_b)
    {
        const std::size_t b1 = b0 + 1;
        c.state.x[0][b1] = xc; c.state.x[1][b1] = kL - yc; c.state.x[2][b1] = -0.5 * D;
        c.state.mass[b1] = 1.3 * m0; c.state.object_id[b1] = 1; c.state.vx[2][b1] = vB;
        c.topology.next[b0] = static_cast<long long>(b1); c.topology.prev[b1] = static_cast<long long>(b0);
        fiber::FiberBond b; b.i = b0; b.j = b1; b.l0 = kL; b.l_break = 1e9;
        c.topology.bonds.push_back(b);
    }
    interactions::compute_fiber_bead_dt_ref(c.state, c.topology, c.materials, c.bead_dt_ref);
    c.prim.a0 = 0;
    c.prim.a1 = point_a ? -1 : 1;
    c.prim.b0 = static_cast<rhs::index_t>(b0);
    c.prim.b1 = point_b ? -1 : static_cast<long long>(b0 + 1);
    return c;
}

interactions::FiberContactStepParams make_params(const Cross& c, int type)
{
    interactions::FiberContactStepParams p;
    p.law.type = type;
    p.law.kp = 1.0;
    p.law.gamma = 0.8;
    p.dt_prev = 1.0e-8;
    p.bead_dt_ref = &c.bead_dt_ref;
    return p;
}

// Вручную (по §2 плана) собранные g, pr; s, t -- параметры точек.
void manual_inputs(const Cross& c, double s, double t, contact::PairGeom& g, contact::PairProps& pr)
{
    const double r = 0.5 * kD;
    const double D = c.state.x[2][0] - c.state.x[2][c.na];
    g.p = 2.0 * r - D;
    auto side_v = [&](rhs::index_t k0, long long k1, double w) {
        return k1 == -1 ? c.state.vx[2][k0] : (1.0 - w) * c.state.vx[2][k0] + w * c.state.vx[2][static_cast<std::size_t>(k1)];
    };
    auto side_inv_m = [&](rhs::index_t k0, long long k1, double w) {
        return k1 == -1 ? 1.0 / c.state.mass[k0]
                        : (1.0 - w) * (1.0 - w) / c.state.mass[k0] + w * w / c.state.mass[static_cast<std::size_t>(k1)];
    };
    const double vA = side_v(c.prim.a0, c.prim.a1, s), vB = side_v(c.prim.b0, c.prim.b1, t);
    g.pdot = -(vA - vB); // n = +z (от B к A)
    g.v_rel = std::fabs(vA - vB);
    g.r = 2.0 * r;
    const double mA = 1.0 / side_inv_m(c.prim.a0, c.prim.a1, s), mB = 1.0 / side_inv_m(c.prim.b0, c.prim.b1, t);
    pr.m_red = mA * mB / (mA + mB);
    const double G = kE / (2.0 * (1.0 + kNu));
    pr.g_eff = 0.5 * G;
    pr.r_eff = 0.5 * r;
    pr.c_pair = std::sqrt(kE / kRho);
    auto side_dt = [&](rhs::index_t k0, long long k1) {
        return k1 == -1 ? c.bead_dt_ref[k0] : std::min(c.bead_dt_ref[k0], c.bead_dt_ref[static_cast<std::size_t>(k1)]);
    };
    pr.dt_ref = std::min(side_dt(c.prim.a0, c.prim.a1), side_dt(c.prim.b0, c.prim.b1));
}

struct KernelOut
{
    std::vector<fiber::FiberDelta> dv;
    std::vector<interactions::FiberFrictionContact> slots;
    std::vector<double> pair_dt;
    bool any_friction = false;
    std::size_t emitted = 0;
};

KernelOut run_kernel(const Cross& c, const std::vector<interactions::FiberContactPrimitive>& prims,
                     interactions::FiberContactStepParams p)
{
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buf;
    KernelOut out;
    p.pair_dt_limit = &out.pair_dt;
    out.any_friction =
        interactions::compute_fiber_contact_contributions(exec, c.state, c.materials, prims, buf, out.slots, p);
    for (const auto& b : buf) out.emitted += b.size();
    out.dv = rhs::reduce_contributions(exec, buf, c.state.size());
    return out;
}

// Сила на стороне A вдоль +z (= fn при n = +z).
double force_on_a(const Cross& c, const KernelOut& k)
{
    double f = 0.0;
    for (std::size_t i = 0; i < c.na; ++i) f += c.state.mass[i] * k.dv[i].dv[2];
    return f;
}

// ---------------------------------------------------------------------------
void test_type0_unchanged()
{
    std::printf("test_type0_unchanged:\n");
    const Cross c = make_cross(0.3, 0.6, 0.9 * kD, -0.5, 0.4, false, false, 0.2);
    rhs::CpuExecutor exec{1};
    const std::vector<interactions::FiberContactPrimitive> prims{c.prim};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> b0, b1;
    std::vector<interactions::FiberFrictionContact> s0, s1;
    interactions::compute_fiber_contact_contributions(exec, c.state, c.materials, prims, b0, s0);
    interactions::FiberContactStepParams p0;
    interactions::compute_fiber_contact_contributions(exec, c.state, c.materials, prims, b1, s1, p0);
    const auto d0 = rhs::reduce_contributions(exec, b0, c.state.size());
    const auto d1 = rhs::reduce_contributions(exec, b1, c.state.size());
    bool same = s0[0].N == s1[0].N && s0[0].N > 0.0;
    for (std::size_t i = 0; i < c.state.size(); ++i)
        for (std::size_t d = 0; d < dim; ++d) same = same && d0[i].dv[d] == d1[i].dv[d];
    check(same, "тип 0 с явными параметрами побитово == вызов без параметров");
    check(interactions::compute_fiber_contact_dt_limit(exec, c.state, c.materials, prims) ==
              interactions::compute_fiber_contact_dt_limit(exec, c.state, c.materials, prims, p0),
          "тип 0: предел шага побитово прежний");
}

// ---------------------------------------------------------------------------
void test_analytic_sign_symmetry()
{
    std::printf("test_analytic_sign_symmetry:\n");
    const double v = 0.6;
    for (int type = 1; type <= 3; ++type)
    {
        const std::string tag = type_name(type);
        const Cross c = make_cross(0.3, 0.6, 0.97 * kD, -v, 0.5 * v);
        contact::PairGeom g;
        contact::PairProps pr;
        manual_inputs(c, 0.3, 0.6, g, pr);
        const auto p = make_params(c, type);
        const contact::Result ref = contact::penalty_force(type, g, pr, p.law.kp, p.dt_prev, p.law.gamma);
        const KernelOut k = run_kernel(c, {c.prim}, p);
        check(g.pdot > 0.0 && ref.fn > 0.0, tag + ": сближение -> pdot > 0, fn > 0");
        check(near_rel(force_on_a(c, k), ref.fn, 1e-12), tag + ": сила на A == contact::penalty_force, по +n");
        check(k.pair_dt.size() == 1 && near_rel(k.pair_dt[0], ref.dt_limit, 1e-12), tag + ": слот предела == res.dt_limit");

        // Симметрия A <-> B: тот же primitive с переставленными сторонами.
        interactions::FiberContactPrimitive sw = c.prim;
        std::swap(sw.a0, sw.b0); std::swap(sw.a1, sw.b1); std::swap(sw.a_prev, sw.b_prev);
        std::swap(sw.a_has_next_segment, sw.b_has_next_segment);
        const KernelOut ks = run_kernel(c, {sw}, p);
        bool same = true;
        for (std::size_t i = 0; i < c.state.size(); ++i)
            for (std::size_t d = 0; d < dim; ++d) same = same && near_rel(ks.dv[i].dv[d], k.dv[i].dv[d], 1e-12);
        check(same && near_rel(ks.pair_dt[0], k.pair_dt[0], 1e-12),
              tag + ": A<->B -- тот же модуль, нормаль противоположна, вклады те же");

        for (const double vz : {0.0, 1.0})
        {
            // vz = 0 -- покой; vz = 1 -- разлёт (A вверх, B вниз).
            const Cross c0 = make_cross(0.3, 0.6, 0.97 * kD, vz * v, -vz * v);
            const KernelOut k0 = run_kernel(c0, {c0.prim}, make_params(c0, type));
            check(k0.emitted == 0, tag + (vz == 0.0 ? ": покой" : ": разлёт") + " -> fn = 0, вкладов нет");
        }
    }
}

// ---------------------------------------------------------------------------
void test_axes_cross()
{
    std::printf("test_axes_cross:\n");
    // D = 0: оси пересекаются, нормаль из e_A x e_B = +z.
    const Cross c = make_cross(0.5, 0.5, 0.0, -0.3, 0.3);
    rhs::CpuExecutor exec{1};
    const auto geom = interactions::compute_fiber_contact_geometry(c.state, c.materials, c.prim);
    check(geom.contact && geom.D == 0.0 && geom.n[2] > 0.99, "самопроверка: контакт при D = 0, n = +z");
    for (int type = 1; type <= 3; ++type)
    {
        const KernelOut k = run_kernel(c, {c.prim}, make_params(c, type));
        check(force_on_a(c, k) > 0.0, std::string(type_name(type)) + ": D = 0 -> fn > 0");
    }
    contact::PairGeom g;
    contact::PairProps pr;
    interactions::fiber_fiber_pair_inputs(c.state, c.materials, c.prim, geom, make_params(c, 2), g, pr);
    check(g.r == kD, "g.r = r_a + r_b при D = 0");
}

// ---------------------------------------------------------------------------
void test_momentum()
{
    std::printf("test_momentum:\n");
    for (int type = 1; type <= 3; ++type)
    {
        int bad = 0, total = 0;
        for (const double s : {0.0, 0.5, 1.0})
            for (const double t : {0.0, 0.5, 1.0})
                for (int pa = 0; pa <= 1; ++pa)
                    for (int pb = 0; pb <= 1; ++pb)
                    {
                        Cross c = make_cross(s, t, 0.95 * kD, -0.4, 0.7, pa != 0, pb != 0);
                        if (!pa) c.state.vx[2][1] = -0.9; // разные скорости бусин A
                        const auto geom = interactions::compute_fiber_contact_geometry(c.state, c.materials, c.prim);
                        if (!geom.contact) { ++bad; ++total; continue; }
                        const KernelOut k = run_kernel(c, {c.prim}, make_params(c, type));
                        double pz = 0.0, scale = std::fabs(force_on_a(c, k));
                        for (std::size_t i = 0; i < c.state.size(); ++i) pz += c.state.mass[i] * k.dv[i].dv[2];
                        ++total;
                        if (!(scale > 0.0) || std::fabs(pz) > 1e-12 * scale) ++bad;
                    }
        check(bad == 0, std::string(type_name(type)) + ": импульс сохраняется (" + std::to_string(total) +
                            " случаев s,t in {0,0.5,1}, точечные стороны)");
    }
}

// ---------------------------------------------------------------------------
void test_consistency_and_friction()
{
    std::printf("test_consistency_and_friction:\n");
    rhs::CpuExecutor exec{1};
    for (int type = 1; type <= 3; ++type)
    {
        const std::string tag = type_name(type);
        // Два примитива с разным прониканием: две пары крестов в одной сцене
        // не строим -- берём два независимых креста и сравниваем по каждому.
        double kernel_min = std::numeric_limits<double>::max();
        double limit_slots = std::numeric_limits<double>::max();
        double limit_recalc = std::numeric_limits<double>::max();
        double diag_min = std::numeric_limits<double>::max();
        bool fn_match = true;
        for (const double D : {0.97 * kD, 0.9 * kD})
        {
            const Cross c = make_cross(0.4, 0.7, D, -0.5, 0.2, false, false, 0.3);
            auto p = make_params(c, type);
            const KernelOut k = run_kernel(c, {c.prim}, p);
            kernel_min = std::min(kernel_min, k.pair_dt[0]);
            p.pair_dt_limit = const_cast<std::vector<double>*>(&k.pair_dt);
            limit_slots = std::min(limit_slots,
                                   interactions::compute_fiber_contact_dt_limit(exec, c.state, c.materials, {c.prim}, p));
            p.pair_dt_limit = nullptr;
            limit_recalc = std::min(limit_recalc,
                                    interactions::compute_fiber_contact_dt_limit(exec, c.state, c.materials, {c.prim}, p));
            interactions::FiberContactSearchResult search;
            search.primitives = {c.prim};
            diag_min = std::min(diag_min, interactions::compute_fiber_contact_diagnostics(exec, c.state, c.materials,
                                                                                          search, p)
                                              .min_dt);
            fn_match = fn_match && k.any_friction && k.slots[0].active &&
                       near_rel(k.slots[0].N, force_on_a(c, k), 1e-12);

            // Покой -> слот неактивен.
            const Cross c0 = make_cross(0.4, 0.7, D, 0.0, 0.0, false, false, 0.3);
            const KernelOut k0 = run_kernel(c0, {c0.prim}, make_params(c0, type));
            fn_match = fn_match && !k0.any_friction && !k0.slots[0].active;
        }
        check(kernel_min < 1e30 && limit_slots == kernel_min && limit_recalc == kernel_min && diag_min == kernel_min,
              tag + ": ядро, предел шага (слоты и пересчёт) и диагностика дают один dt_limit");
        check(fn_match, tag + ": slot.N == fn; при fn = 0 слот неактивен");
    }
}

// ---------------------------------------------------------------------------
// Тип 4 (Герц + демпфер) на кресте: F = F_H + c*max(pdot, 0); в покое и при
// разлёте F = F_H; D = 0 -- тоже F > 0; импульс; согласованность ядра,
// предела шага и диагностики.
void test_type4_hertz_damped()
{
    std::printf("test_type4_hertz_damped:\n");
    rhs::CpuExecutor exec{1};
    for (const double vz : {-0.6, 0.0, 0.6})
    {
        const Cross c = make_cross(0.3, 0.6, 0.97 * kD, vz, -0.5 * vz, false, false, 0.2);
        auto p = make_params(c, interactions::kFiberLawHertzDamped);
        contact::PairGeom g;
        contact::PairProps pr;
        manual_inputs(c, 0.3, 0.6, g, pr);
        const auto geom = interactions::compute_fiber_contact_geometry(c.state, c.materials, c.prim);
        const double c_d = p.law.kp * pr.m_red / pr.dt_ref;
        const double f_ref = geom.F + (g.pdot > 0.0 ? c_d * g.pdot : 0.0);
        const KernelOut k = run_kernel(c, {c.prim}, p);
        const std::string mode = vz < 0.0 ? "сближение" : (vz == 0.0 ? "покой" : "разлёт");
        check(geom.F > 0.0 && near_rel(force_on_a(c, k), f_ref, 1e-12), "Hertz+damper, " + mode + ": F = F_H + c*max(pdot,0)");
        check(k.slots[0].active && near_rel(k.slots[0].N, f_ref, 1e-12), "Hertz+damper, " + mode + ": слот трения активен");
        double pz = 0.0;
        for (std::size_t i = 0; i < c.state.size(); ++i) pz += c.state.mass[i] * k.dv[i].dv[2];
        check(std::fabs(pz) <= 1e-12 * f_ref, "Hertz+damper, " + mode + ": импульс сохраняется");
        p.pair_dt_limit = nullptr;
        const double lim = interactions::compute_fiber_contact_dt_limit(exec, c.state, c.materials, {c.prim}, p);
        interactions::FiberContactSearchResult search;
        search.primitives = {c.prim};
        const double diag = interactions::compute_fiber_contact_diagnostics(exec, c.state, c.materials, search, p).min_dt;
        check(k.pair_dt[0] == lim && lim == diag, "Hertz+damper, " + mode + ": ядро, предел шага и диагностика согласованы");
    }
    const Cross cx = make_cross(0.5, 0.5, 0.0, 0.0, 0.0);
    const KernelOut kx = run_kernel(cx, {cx.prim}, make_params(cx, interactions::kFiberLawHertzDamped));
    check(force_on_a(cx, kx) > 0.0, "Hertz+damper: D = 0, покой -> F = F_H > 0");
}

// ---------------------------------------------------------------------------
// Тип 3 с выталкиванием (FiberContactLaw::push_tau): глубокий крест p = 0.3 d.
//   - push_tau = 0 -- побитово прежний тип 3;
//   - мелкое проникание (p < push_depth*r) -- выталкивание не включается;
//   - глубокое, медленный разлёт: прежний тип 3 даёт 0, с выталкиванием
//     fn = gamma*m_red*(pdot + excess/tau)/dt_prev > 0, слот трения активен;
//   - быстрый разлёт (pdot + excess/tau <= 0): fn = 0 (контакт не тянет).
void test_type3_push()
{
    std::printf("test_type3_push:\n");
    const double tau = 1.0e-6, depth = 0.1;
    auto with_push = [&](const Cross& c) {
        auto p = make_params(c, contact::Impulse);
        p.law.push_tau = tau;
        p.law.push_depth = depth;
        return p;
    };
    {
        const Cross c = make_cross(0.3, 0.6, 0.7 * kD, -0.5, 0.4, false, false, 0.2);
        auto p0 = make_params(c, contact::Impulse);
        auto p1 = p0;
        p1.law.push_tau = 0.0;
        p1.law.push_depth = 0.37;
        const KernelOut k0 = run_kernel(c, {c.prim}, p0), k1 = run_kernel(c, {c.prim}, p1);
        bool same = k0.slots[0].N == k1.slots[0].N;
        for (std::size_t i = 0; i < c.state.size(); ++i)
            for (std::size_t d = 0; d < dim; ++d) same = same && k0.dv[i].dv[d] == k1.dv[i].dv[d];
        check(same, "push_tau = 0: побитово прежний тип 3");
    }
    {
        const Cross c = make_cross(0.3, 0.6, 0.97 * kD, -0.5, 0.4);
        const double f_old = force_on_a(c, run_kernel(c, {c.prim}, make_params(c, contact::Impulse)));
        const double f_new = force_on_a(c, run_kernel(c, {c.prim}, with_push(c)));
        check(f_old > 0.0 && f_new == f_old, "p < push_depth*r: сила как у прежнего типа 3");
    }
    {
        const Cross c = make_cross(0.3, 0.6, 0.7 * kD, 0.1, -0.05, false, false, 0.2);
        contact::PairGeom g;
        contact::PairProps pr;
        manual_inputs(c, 0.3, 0.6, g, pr);
        auto p = with_push(c);
        const double excess = g.p - depth * g.r;
        const double f_ref = p.law.gamma * pr.m_red * (g.pdot + excess / tau) / p.dt_prev;
        const double f_old = force_on_a(c, run_kernel(c, {c.prim}, make_params(c, contact::Impulse)));
        const KernelOut k = run_kernel(c, {c.prim}, p);
        check(g.pdot < 0.0 && f_old == 0.0, "глубоко, медленный разлёт: прежний тип 3 даёт 0");
        check(f_ref > 0.0 && near_rel(force_on_a(c, k), f_ref, 1e-12),
              "глубоко, медленный разлёт: fn = gamma*m_red*(pdot + excess/tau)/dt_prev");
        check(k.slots[0].active && near_rel(k.slots[0].N, f_ref, 1e-12), "выталкивание: слот трения активен, N = fn");
        double pz = 0.0;
        for (std::size_t i = 0; i < c.state.size(); ++i) pz += c.state.mass[i] * k.dv[i].dv[2];
        check(std::fabs(pz) <= 1e-12 * f_ref, "выталкивание: импульс сохраняется");
    }
    {
        // Предел скорости выталкивания: excess/tau = 200 м/с, предел 20 м/с ->
        // fn = gamma*m_red*(pdot + 20)/dt_prev; без предела -- прежняя формула.
        const Cross c = make_cross(0.3, 0.6, 0.7 * kD, 0.1, -0.05, false, false, 0.2);
        contact::PairGeom g;
        contact::PairProps pr;
        manual_inputs(c, 0.3, 0.6, g, pr);
        auto p = with_push(c);
        p.law.push_vmax = 20.0;
        const double f_ref = p.law.gamma * pr.m_red * (g.pdot + 20.0) / p.dt_prev;
        check(near_rel(force_on_a(c, run_kernel(c, {c.prim}, p)), f_ref, 1e-12),
              "push_vmax = 20 м/с: скорость выталкивания ограничена (200 -> 20 м/с)");
        p.law.push_vmax = 1e9;
        const double f_free = p.law.gamma * pr.m_red * (g.pdot + (g.p - depth * g.r) / tau) / p.dt_prev;
        check(near_rel(force_on_a(c, run_kernel(c, {c.prim}, p)), f_free, 1e-12),
              "push_vmax больше excess/tau: формула без изменений");
    }
    {
        // excess/tau = 0.2e-3/1e-6 = 200 м/с; разлёт 300 м/с перекрывает выталкивание.
        const Cross c = make_cross(0.3, 0.6, 0.7 * kD, 150.0, -150.0);
        check(force_on_a(c, run_kernel(c, {c.prim}, with_push(c))) == 0.0, "быстрый разлёт: fn = 0, контакт не тянет");
    }
}

} // namespace

int main()
{
    test_type0_unchanged();
    test_analytic_sign_symmetry();
    test_axes_cross();
    test_momentum();
    test_consistency_and_friction();
    test_type4_hertz_damped();
    test_type3_push();
    std::printf("\n%s: %d failure(s)\n", failures == 0 ? "PASSED" : "FAILED", failures);
    return failures == 0 ? 0 : 1;
}
