// tests/test_fiber_solid_types.cpp
// Законы нормальной силы pinball / Swegle / импульсный для Fiber-Solid-Contact
// (fiber/plan-fiber-contact-types.md, этап 1, шаг 1.6;
// include/interactions/fiber_contact_law.hpp).
//   1. тип 0 с явными параметрами == вызов без параметров (побитово);
//   2. аналитика: сила ядра == contact::penalty_force на вручную собранных g, pr;
//   3. знак: сближение -> fn > 0; покой и разлёт -> fn = 0 у всех типов 1-3;
//   4. баланс импульса F_a + F_b + F_s = 0, lambda in {0, 0.3, 1} и точка;
//   5. трение: slot.N == fn, при fn = 0 слот неактивен;
//   6. dt_ref нити: связь, изолированная бусина, вязкая ветка pinball;
//      предел шага = минимум res.dt_limit по принятым контактам;
//   7. парсер ключей;
//   8. e2e: нить падает на пластину, типы 0-3, метрики в вывод.
#include "contact_force.hpp"
#include "fiber/fiber_system.hpp"
#include "input.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact_law.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "materials.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_solid_types требует сборки с SPH_DIM=3");

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
        default: return "?";
    }
}

// --- Сцена «одна пара»: сегмент (или точка) нити + одна Solid-частица ------
constexpr double kE = 2.0e9, kNu = 0.35, kD = 1.0e-3, kRho = 1100.0;
constexpr double kL = 2.0e-3;      // длина сегмента
constexpr double kMs = 2.0e-6;     // масса Solid-частицы
constexpr double kRhoS = 2700.0;   // Al
constexpr int kItypeS = 4;         // Al (materials.hpp)
constexpr double kCs = 5000.0;     // sys.c частицы

struct PairScene
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    ParticleSystem sys;
    Materials smats;
    SimulationConfig cfg;
    std::vector<double> bead_dt_ref;
    interactions::FiberSolidContactPrimitive prim;
    double r_s = 0.0;
};

// segment: две бусины (0,0,0)-(L,0,0) со связью; иначе одна изолированная
// бусина в начале координат. Частица над точкой lambda*L на расстоянии
// D = (r_f + r_s)*(1 - pen) по z. Нить движется по z со скоростью vz_f.
PairScene make_pair_scene(bool segment, double lambda, double pen, double vz_f, double mu_fs = 0.0)
{
    PairScene sc;
    fiber::FiberMaterial m;
    m.E = kE; m.nu = kNu; m.d = kD; m.rho = kRho; m.mu_fs = mu_fs;
    sc.materials.push_back(m);
    const std::size_t nb = segment ? 2 : 1;
    sc.state.resize(nb);
    sc.topology.resize(nb);
    const double A = interactions::fiber_cross_section_area(m);
    for (std::size_t k = 0; k < nb; ++k)
    {
        sc.state.x[0][k] = static_cast<double>(k) * kL;
        sc.state.x[1][k] = 0.0;
        sc.state.x[2][k] = 0.0;
        sc.state.vx[0][k] = 0.0;
        sc.state.vx[1][k] = 0.0;
        sc.state.vx[2][k] = vz_f;
        // Разные массы бусин -- чтобы веса lambda проверялись нетривиально.
        sc.state.mass[k] = kRho * A * kL * (k == 0 ? 1.0 : 1.7);
        sc.state.object_id[k] = 0;
        sc.state.material_id[k] = 0;
    }
    if (segment)
    {
        sc.topology.next[0] = 1;
        sc.topology.prev[1] = 0;
        fiber::FiberBond b;
        b.i = 0; b.j = 1; b.l0 = kL; b.l_break = 1e9;
        sc.topology.bonds.push_back(b);
    }

    sc.sys.resize(1);
    sc.sys.inzializ_index(0);
    sc.sys.mass[0] = kMs;
    sc.sys.rho[0] = kRhoS;
    sc.sys.rho_0[0] = kRhoS;
    sc.sys.itype[0] = kItypeS;
    sc.sys.body[0] = kItypeS;
    sc.sys.hsml[0] = 1.0e-3;
    sc.sys.c[0] = kCs;
    sc.r_s = interactions::solid_contact_radius(sc.sys, 0);
    const double D = (0.5 * kD + sc.r_s) * (1.0 - pen);
    sc.sys.x[0][0] = segment ? lambda * kL : 0.0;
    sc.sys.x[1][0] = 0.0;
    sc.sys.x[2][0] = D;
    for (std::size_t d = 0; d < dim; ++d) sc.sys.vx[d][0] = 0.0;

    sc.cfg.apply_defaults();
    interactions::compute_fiber_bead_dt_ref(sc.state, sc.topology, sc.materials, sc.bead_dt_ref);

    if (segment)
    {
        sc.prim.a = 0; sc.prim.b = 1; sc.prim.s = 0;
        sc.prim.has_b = true; sc.prim.a_prev = -1; sc.prim.b_is_last = true;
    }
    else
    {
        sc.prim.a = 0; sc.prim.s = 0; sc.prim.has_b = false; sc.prim.a_prev = -1;
    }
    return sc;
}

interactions::FiberContactStepParams make_params(const PairScene& sc, int type, double dt_prev = 1.0e-8)
{
    interactions::FiberContactStepParams p;
    p.law.type = type;
    p.law.kp = 1.0;
    p.law.gamma = 0.8;
    p.dt_prev = dt_prev;
    p.cfg = &sc.cfg;
    p.bead_dt_ref = &sc.bead_dt_ref;
    return p;
}

// Вручную (по §2 плана, независимо от ядра) собранные g, pr.
void manual_inputs(const PairScene& sc, double lambda, contact::PairGeom& g, contact::PairProps& pr)
{
    const fiber::FiberMaterial& m = sc.materials[0];
    const double r_f = 0.5 * m.d;
    const double D = sc.sys.x[2][0];
    g.p = r_f + sc.r_s - D;
    const bool seg = sc.state.size() == 2;
    const double vq = seg ? (1.0 - lambda) * sc.state.vx[2][0] + lambda * sc.state.vx[2][1] : sc.state.vx[2][0];
    // n = +z (от нити к Solid), v_s = 0.
    g.pdot = -(0.0 - vq);
    g.v_rel = std::fabs(vq);
    g.r = r_f + sc.r_s;
    const double inv_mq = seg ? (1.0 - lambda) * (1.0 - lambda) / sc.state.mass[0] + lambda * lambda / sc.state.mass[1]
                              : 1.0 / sc.state.mass[0];
    const double mq = 1.0 / inv_mq;
    pr.m_red = mq * kMs / (mq + kMs);
    const double Gf = kE / (2.0 * (1.0 + kNu));
    const double Gs = shear_modulus_by_itype(sc.smats, kItypeS);
    pr.g_eff = Gf * Gs / (Gf + Gs);
    pr.r_eff = r_f * sc.r_s / (r_f + sc.r_s);
    pr.c_pair = std::min(std::sqrt(kE / kRho), kCs);
    const double dt_f = seg ? std::min(sc.bead_dt_ref[0], sc.bead_dt_ref[1]) : sc.bead_dt_ref[0];
    const double dt_s = sc.cfg.CFL * sc.sys.hsml[0] / kCs; // v_s = 0
    pr.dt_ref = std::min(dt_f, dt_s);
}

struct KernelOut
{
    std::vector<fiber::FiberDelta> fd;
    std::vector<interactions::SolidDelta> sd;
    std::vector<interactions::FiberSolidFrictionContact> slots;
    std::vector<double> pair_dt;
    bool any_friction = false;
    std::size_t emitted = 0;
};

KernelOut run_kernel(const PairScene& sc, interactions::FiberContactStepParams p)
{
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fb;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sb;
    KernelOut out;
    p.pair_dt_limit = &out.pair_dt;
    const std::vector<interactions::FiberSolidContactPrimitive> prims{sc.prim};
    out.any_friction = interactions::compute_fiber_solid_contact_contributions(
        exec, sc.smats, sc.sys, sc.state, sc.materials, prims, fb, sb, out.slots, p);
    for (const auto& b : fb) out.emitted += b.size();
    for (const auto& b : sb) out.emitted += b.size();
    out.fd = rhs::reduce_contributions(exec, fb, sc.state.size());
    out.sd = rhs::reduce_contributions(exec, sb, sc.sys.size());
    return out;
}

// ---------------------------------------------------------------------------
void test_type0_unchanged()
{
    std::printf("test_type0_unchanged:\n");
    const PairScene sc = make_pair_scene(true, 0.3, 0.05, 0.7, 0.2);
    rhs::CpuExecutor exec{1};
    const std::vector<interactions::FiberSolidContactPrimitive> prims{sc.prim};

    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fb0, fb1;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sb0, sb1;
    std::vector<interactions::FiberSolidFrictionContact> s0, s1;
    interactions::compute_fiber_solid_contact_contributions(exec, sc.smats, sc.sys, sc.state, sc.materials, prims,
                                                             fb0, sb0, s0);
    interactions::FiberContactStepParams p0; // type 0 явно, без bead_dt_ref/cfg
    interactions::compute_fiber_solid_contact_contributions(exec, sc.smats, sc.sys, sc.state, sc.materials, prims,
                                                             fb1, sb1, s1, p0);
    const auto f0 = rhs::reduce_contributions(exec, fb0, 2), f1 = rhs::reduce_contributions(exec, fb1, 2);
    const auto d0 = rhs::reduce_contributions(exec, sb0, 1), d1 = rhs::reduce_contributions(exec, sb1, 1);
    bool same = s0[0].N == s1[0].N;
    for (std::size_t d = 0; d < dim; ++d)
    {
        same = same && f0[0].dv[d] == f1[0].dv[d] && f0[1].dv[d] == f1[1].dv[d] && d0[0].dv[d] == d1[0].dv[d];
    }
    check(same, "тип 0 с явными параметрами побитово == вызов без параметров");
    const double l0 = interactions::compute_fiber_solid_contact_dt_limit(exec, sc.smats, sc.sys, sc.state, sc.materials, prims);
    const double l1 =
        interactions::compute_fiber_solid_contact_dt_limit(exec, sc.smats, sc.sys, sc.state, sc.materials, prims, p0);
    check(l0 == l1, "тип 0: предел шага побитово прежний");
}

// ---------------------------------------------------------------------------
void test_analytic_and_sign()
{
    std::printf("test_analytic_and_sign:\n");
    const double v = 0.9;
    for (int type = 1; type <= 3; ++type)
    {
        for (int seg = 0; seg <= 1; ++seg)
        {
            const double lambda = seg ? 0.3 : 0.0;
            const std::string tag = std::string(type_name(type)) + (seg ? ", сегмент" : ", точка");

            // Сближение.
            const PairScene sc = make_pair_scene(seg != 0, lambda, 0.02, v);
            const auto p = make_params(sc, type);
            contact::PairGeom g;
            contact::PairProps pr;
            manual_inputs(sc, lambda, g, pr);
            const contact::Result ref = contact::penalty_force(type, g, pr, p.law.kp, p.dt_prev, p.law.gamma);
            const KernelOut k = run_kernel(sc, p);
            check(g.pdot > 0.0 && ref.fn > 0.0, tag + ": сближение -> pdot > 0, fn > 0");
            check(near_rel(k.sd[0].dv[2] * kMs, ref.fn, 1e-12), tag + ": сила на Solid == contact::penalty_force");
            check(k.sd[0].dv[0] == 0.0 && k.sd[0].dv[1] == 0.0, tag + ": сила направлена по n (+z)");
            check(k.pair_dt.size() == 1 && near_rel(k.pair_dt[0], ref.dt_limit, 1e-12),
                  tag + ": слот предела шага == res.dt_limit");

            // Баланс импульса.
            double px = 0.0, pz = 0.0, scale = ref.fn;
            for (std::size_t i = 0; i < sc.state.size(); ++i)
            {
                px += sc.state.mass[i] * k.fd[i].dv[0];
                pz += sc.state.mass[i] * k.fd[i].dv[2];
            }
            px += kMs * k.sd[0].dv[0];
            pz += kMs * k.sd[0].dv[2];
            check(std::fabs(px) <= 1e-12 * scale && std::fabs(pz) <= 1e-12 * scale, tag + ": F_a + F_b + F_s = 0");

            // Покой и разлёт.
            for (const double vz : {0.0, -v})
            {
                const PairScene sc0 = make_pair_scene(seg != 0, lambda, 0.02, vz);
                const KernelOut k0 = run_kernel(sc0, make_params(sc0, type));
                check(k0.emitted == 0, tag + (vz == 0.0 ? ": покой" : ": разлёт") + " -> fn = 0, вкладов нет");
            }
        }
    }
}

// ---------------------------------------------------------------------------
void test_momentum_lambdas()
{
    std::printf("test_momentum_lambdas:\n");
    for (int type = 1; type <= 3; ++type)
        for (const double lambda : {0.0, 0.3, 1.0})
        {
            PairScene sc = make_pair_scene(true, lambda, 0.03, 0.5);
            // Разные скорости бусин: v_q зависит от lambda.
            sc.state.vx[2][1] = 1.1;
            const KernelOut k = run_kernel(sc, make_params(sc, type));
            contact::PairGeom g;
            contact::PairProps pr;
            manual_inputs(sc, lambda, g, pr);
            const auto p = make_params(sc, type);
            const contact::Result ref = contact::penalty_force(type, g, pr, p.law.kp, p.dt_prev, p.law.gamma);
            double pz = kMs * k.sd[0].dv[2];
            for (std::size_t i = 0; i < 2; ++i) pz += sc.state.mass[i] * k.fd[i].dv[2];
            char buf[96];
            std::snprintf(buf, sizeof buf, "%s, lambda=%.1f", type_name(type), lambda);
            check(std::fabs(pz) <= 1e-12 * ref.fn, std::string(buf) + ": импульс сохраняется");
            check(near_rel(sc.state.mass[0] * k.fd[0].dv[2], -(1.0 - lambda) * ref.fn, 1e-12) &&
                      near_rel(sc.state.mass[1] * k.fd[1].dv[2], -lambda * ref.fn, 1e-12),
                  std::string(buf) + ": раскладка по весам (1-lambda, lambda)");
        }
}

// ---------------------------------------------------------------------------
void test_friction_slot()
{
    std::printf("test_friction_slot:\n");
    for (int type = 1; type <= 3; ++type)
    {
        const PairScene sc = make_pair_scene(true, 0.4, 0.02, 0.6, 0.3);
        const KernelOut k = run_kernel(sc, make_params(sc, type));
        check(k.any_friction && k.slots[0].active && near_rel(k.slots[0].N, k.sd[0].dv[2] * kMs, 1e-12),
              std::string(type_name(type)) + ": slot.N == fn");
        const PairScene sc0 = make_pair_scene(true, 0.4, 0.02, 0.0, 0.3);
        const KernelOut k0 = run_kernel(sc0, make_params(sc0, type));
        check(!k0.any_friction && !k0.slots[0].active, std::string(type_name(type)) + ": fn = 0 -> слот неактивен");
    }
}

// ---------------------------------------------------------------------------
void test_dt_ref_and_limit()
{
    std::printf("test_dt_ref_and_limit:\n");
    const PairScene seg = make_pair_scene(true, 0.5, 0.02, 0.5);
    const double A = interactions::fiber_cross_section_area(seg.materials[0]);
    const double m0 = seg.state.mass[0], m1 = seg.state.mass[1];
    const double dt_bond = 0.8 / std::sqrt((kE * A / kL) / (m0 * m1 / (m0 + m1)));
    check(near_rel(seg.bead_dt_ref[0], dt_bond, 1e-14) && near_rel(seg.bead_dt_ref[1], dt_bond, 1e-14),
          "бусины со связью: dt_ref_f == dt_bond");
    rhs::CpuExecutor exec{1};
    const double bond_limit = interactions::compute_fiber_bond_dt_limit(exec, seg.state, seg.topology, seg.materials);
    check(seg.bead_dt_ref[0] == bond_limit, "dt_bead побитово == compute_fiber_bond_dt_limit (одна связь)");

    const PairScene pt = make_pair_scene(false, 0.0, 0.02, 0.5);
    const double mp = pt.state.mass[0];
    const double dt_iso = 0.8 / std::sqrt((kE * A / kD) / (0.5 * mp));
    check(near_rel(pt.bead_dt_ref[0], dt_iso, 1e-14), "изолированная бусина: виртуальная связь длины d");

    // Вязкая ветка pinball при Kp = 1 не режет шаг сильнее связей:
    // 0.8*2*dt_ref >= dt_bond при dt_ref = dt_bond (Solid-сторона не мельче).
    contact::PairGeom g;
    g.p = 1e-12; g.pdot = 1e3; g.r = 1e-3; g.v_rel = 1e3; // F1 << F2 невозможно -> вязкая ветка
    contact::PairProps pr;
    pr.m_red = 1e-6; pr.c_pair = 100.0; pr.g_eff = 1e9; pr.r_eff = 1e-4; pr.dt_ref = dt_bond;
    const contact::Result rp = contact::pinball_force(g, pr, 1.0);
    check(rp.dt_limit >= dt_bond, "pinball, Kp = 1: вязкая ветка не ограничивает шаг сильнее связей");

    // Предел шага ядра = минимум res.dt_limit по принятым контактам: две
    // частицы с разным прониканием над одним сегментом.
    for (int type = 1; type <= 3; ++type)
    {
        PairScene sc = make_pair_scene(true, 0.3, 0.02, 0.8);
        sc.sys.resize(2);
        sc.sys.inzializ_index(1);
        sc.sys.mass[1] = 0.5 * kMs;
        sc.sys.rho[1] = kRhoS; sc.sys.rho_0[1] = kRhoS;
        sc.sys.itype[1] = kItypeS; sc.sys.body[1] = kItypeS;
        sc.sys.hsml[1] = 1.0e-3; sc.sys.c[1] = kCs;
        sc.sys.x[0][1] = 0.8 * kL; sc.sys.x[1][1] = 0.0; sc.sys.x[2][1] = (0.5 * kD + sc.r_s) * 0.9;
        sc.sys.x[0][0] = 0.3 * kL; // первая -- как в make_pair_scene
        std::vector<interactions::FiberSolidContactPrimitive> prims{sc.prim, sc.prim};
        prims[1].s = 1;
        auto p = make_params(sc, type);
        std::vector<double> slots;
        p.pair_dt_limit = &slots;
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fb;
        std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sb;
        std::vector<interactions::FiberSolidFrictionContact> fs;
        interactions::compute_fiber_solid_contact_contributions(exec, sc.smats, sc.sys, sc.state, sc.materials, prims,
                                                                 fb, sb, fs, p);
        double expect = std::numeric_limits<double>::max();
        for (const auto& prim : prims)
        {
            const auto geom = interactions::compute_fiber_solid_contact_geometry(sc.smats, sc.sys, sc.state,
                                                                               sc.materials, prim);
            if (!geom.contact) continue;
            expect = std::min(expect, interactions::fiber_solid_contact_penalty(sc.smats, sc.sys, sc.state,
                                                                                 sc.materials, prim, geom, p)
                                          .dt_limit);
        }
        const double with_slots =
            interactions::compute_fiber_solid_contact_dt_limit(exec, sc.smats, sc.sys, sc.state, sc.materials, prims, p);
        auto p_noslot = p;
        p_noslot.pair_dt_limit = nullptr;
        const double recomputed = interactions::compute_fiber_solid_contact_dt_limit(exec, sc.smats, sc.sys, sc.state,
                                                                                     sc.materials, prims, p_noslot);
        check(expect < 1e30 && with_slots == expect && recomputed == expect,
              std::string(type_name(type)) + ": предел = минимум res.dt_limit (по слотам и пересчётом)");
    }
}

// ---------------------------------------------------------------------------
void test_parser()
{
    std::printf("test_parser:\n");
    const std::filesystem::path p = std::filesystem::temp_directory_path() / "sph_test_fiber_solid_types_cfg.txt";
    auto parse = [&](const std::string& body, SimulationConfig& cfg) -> std::string {
        {
            std::ofstream f(p);
            f << body;
        }
        std::ostringstream err;
        std::streambuf* old = std::cerr.rdbuf(err.rdbuf());
        cfg.apply_defaults();
        parse_sph_config(cfg, p.string());
        std::cerr.rdbuf(old);
        return err.str();
    };

    SimulationConfig cfg;
    std::string e = parse("", cfg);
    check(cfg.fiber_solid_contact_type == 0 && cfg.fiber_fiber_contact_type == 0 &&
              cfg.fiber_solid_contact_penalty_scale == 1.0 && cfg.fiber_solid_contact_impulse_gamma == 1.0 &&
              cfg.fiber_fiber_contact_penalty_scale == 1.0 && cfg.fiber_fiber_contact_impulse_gamma == 1.0,
          "умолчания: типы 0 (Герц), Kp = 1, gamma = 1");

    bool all_ok = true;
    for (int t = 0; t <= 3; ++t)
    {
        e = parse("fiber_solid_contact_type = " + std::to_string(t) + "\n", cfg);
        all_ok = all_ok && cfg.fiber_solid_contact_type == t && e.empty();
    }
    check(all_ok, "fiber_solid_contact_type 0..3 принимаются без предупреждений");

    e = parse("fiber_solid_contact_type = 5\nfiber_fiber_contact_type = -1\n", cfg);
    check(cfg.fiber_solid_contact_type == 0 && cfg.fiber_fiber_contact_type == 0 && !e.empty(),
          "вне диапазона -> сброс в 0 с сообщением");

    e = parse("fiber_fiber_contact_type = 1\n", cfg);
    check(cfg.fiber_fiber_contact_type == 0 && e.find("не реализован") != std::string::npos,
          "fiber_fiber_contact_type = 1 (этап 1) -> предупреждение «не реализован», Герц");

    e = parse("fiber_solid_contact_penalty_scale = 2.5\nfiber_fiber_contact_penalty_scale = 0.5\n", cfg);
    check(cfg.fiber_solid_contact_penalty_scale == 2.5 && cfg.fiber_fiber_contact_penalty_scale == 0.5,
          "penalty_scale разбирается");

    e = parse("fiber_solid_contact_impulse_gamma = 1.7\n", cfg);
    check(cfg.fiber_solid_contact_impulse_gamma == 1.7 && e.find("1.6") != std::string::npos,
          "gamma > 1.6 -> предупреждение, значение принято");
    e = parse("fiber_fiber_contact_impulse_gamma = 0\n", cfg);
    check(cfg.fiber_fiber_contact_impulse_gamma == 0.0 && !e.empty(), "gamma <= 0 -> предупреждение");

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

// ---------------------------------------------------------------------------
// 8. e2e: нить падает на неподвижную пластину (сцена test_e2e_fiber_solid.cpp,
// (Б)), локальный PC-драйвер с законом params. Пластина: c = 5000 м/с,
// hsml = 1 м -- CFL-шаг Solid-стороны (2e-5 с) крупнее шага связей, так что
// dt_ref пары задаёт нить, как в связке с реальной Solid-частью, где шаг уже
// ограничен CFL Solid.
// ---------------------------------------------------------------------------
struct E2eMetrics
{
    bool finite = true;
    double max_pen_over_d = 0.0;   // max delta / d_f
    double end_pen_over_d = 0.0;   // delta / d_f в конце (статический контакт)
    double min_dt = std::numeric_limits<double>::max();
    double max_dt_ratio = 0.0;     // max dt_new/dt_prev
    int steps = 0;
    int steps_limited_by_fs = 0;
    int contact_steps = 0;         // шаги с геометрическим контактом
    int contact_steps_zero_fn = 0; // из них шаги, где у всех контактов fn = 0
    double e0 = 0.0, e_end = 0.0, e_max = 0.0;
    double vz_end = 0.0;
    double pen_at[4] = {}; // delta/d_f в моменты t_end/4, t_end/2, 3t_end/4, t_end
    double dt_bond0 = 0.0; // safety * шаг связей в начале
};

void build_plate(ParticleSystem& sys, int ncols_half, double dxp, double rho)
{
    const int side = 2 * ncols_half;
    sys.resize(static_cast<std::size_t>(side * side));
    std::size_t k = 0;
    for (int ix = -ncols_half; ix < ncols_half; ++ix)
        for (int iy = -ncols_half; iy < ncols_half; ++iy)
        {
            sys.inzializ_index(k);
            sys.x[0][k] = (ix + 0.5) * dxp;
            sys.x[1][k] = (iy + 0.5) * dxp;
            sys.x[2][k] = 0.0;
            sys.itype[k] = 4;
            sys.body[k] = 4;
            sys.rho[k] = rho;
            sys.rho_0[k] = rho;
            sys.mass[k] = rho * dxp * dxp * dxp;
            sys.hsml[k] = 1.0;
            sys.c[k] = 5000.0;
            ++k;
        }
}

E2eMetrics run_plate(int type)
{
    const double dxp = 1.5e-3;
    ParticleSystem plate;
    build_plate(plate, 5, dxp, 2700.0);
    const std::size_t nplate = plate.size();
    Materials smats;
    SimulationConfig cfg;
    cfg.apply_defaults();
    const double r_s = 0.5 * dxp;

    const double E_f = 2.0e7, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    const double r_f = 0.5 * d_f, l0 = 2.0e-3;
    const int N = 5;
    fiber::FiberMaterial fm;
    fm.E = E_f; fm.nu = nu_f; fm.d = d_f; fm.rho = rho_f;
    std::vector<fiber::FiberMaterial> fmats{fm};
    const double A = interactions::fiber_cross_section_area(fm);
    const double mass_each = rho_f * A * l0;
    const double v0 = 0.5, gap0 = 2.0e-3;
    const double z0 = r_s + r_f + gap0;

    fiber::FiberState fb;
    fiber::FiberTopology topo;
    fb.resize(N);
    topo.resize(N);
    for (int k = 0; k < N; ++k)
    {
        const std::size_t i = static_cast<std::size_t>(k);
        fb.x[0][i] = (k - N / 2) * l0; fb.x[1][i] = 0.0; fb.x[2][i] = z0;
        fb.vx[0][i] = 0.0; fb.vx[1][i] = 0.0; fb.vx[2][i] = -v0;
        fb.mass[i] = mass_each; fb.object_id[i] = 0; fb.material_id[i] = 0;
        if (k > 0) topo.prev[i] = k - 1;
        if (k + 1 < N)
        {
            topo.next[i] = k + 1;
            fiber::FiberBond b;
            b.i = i; b.j = i + 1; b.l0 = l0; b.l_break = 1e9;
            topo.bonds.push_back(b);
        }
    }

    auto energy = [&]() {
        double e = 0.0;
        for (std::size_t i = 0; i < fb.size(); ++i)
            for (std::size_t d = 0; d < dim; ++d) e += 0.5 * fb.mass[i] * fb.vx[d][i] * fb.vx[d][i];
        for (const auto& b : topo.bonds)
        {
            double l2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double dx = fb.x[d][b.j] - fb.x[d][b.i];
                l2 += dx * dx;
            }
            const double s = std::sqrt(l2) - b.l0;
            if (s > 0.0) e += 0.5 * (E_f * A / b.l0) * s * s;
        }
        return e;
    };

    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> bond_buffers, fs_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fs_solid_buffers;
    std::vector<interactions::FiberSolidFrictionContact> slots;
    interactions::FiberSolidContactSearchResult search;
    std::vector<double> bead_dt_ref;
    std::vector<double> dv_prev[dim], v_min[dim];
    for (std::size_t d = 0; d < dim; ++d) { dv_prev[d].assign(N, 0.0); v_min[d].assign(N, 0.0); }
    bool first = true;
    const double safety = 0.3;
    double dt_prev = 0.0;

    E2eMetrics m;
    m.e0 = energy();
    m.e_max = m.e0;
    const double t_end = 2.4e-2;
    int next_sample = 0;
    double time = 0.0;
    const std::size_t center = static_cast<std::size_t>(N / 2);
    while (time < t_end && m.steps < 40000)
    {
        interactions::FiberContactStepParams p;
        p.law.type = type;
        p.dt_prev = dt_prev > 0.0 ? dt_prev : 1e-7;
        p.cfg = &cfg;
        p.bead_dt_ref = &bead_dt_ref;
        p.pair_dt_limit = &search.pair_dt_limit;
        if (type != 0) interactions::compute_fiber_bead_dt_ref(fb, topo, fmats, bead_dt_ref);

        const auto bond_delta = interactions::compute_fiber_bond(exec, fb, topo, fmats, bond_buffers);
        const double bond_dt = interactions::compute_fiber_bond_dt_limit(exec, fb, topo, fmats);
        interactions::build_fiber_solid_contact_search(exec, smats, plate, nplate, fb, topo, fmats, search);
        interactions::compute_fiber_solid_contact_contributions(exec, smats, plate, fb, fmats, search.primitives,
                                                                 fs_fiber_buffers, fs_solid_buffers, slots, p);
        const auto fs_delta = rhs::reduce_contributions(exec, fs_fiber_buffers, fb.size());
        const double fs_dt =
            interactions::compute_fiber_solid_contact_dt_limit(exec, smats, plate, fb, fmats, search.primitives, p);

        double pen = 0.0;
        bool any_contact = false, any_force = false;
        for (const auto& prim : search.primitives)
        {
            const auto g = interactions::compute_fiber_solid_contact_geometry(smats, plate, fb, fmats, prim);
            if (!g.contact) continue;
            any_contact = true;
            pen = std::max(pen, g.delta);
            const double fn = (type == 0) ? g.F
                                          : interactions::fiber_solid_contact_penalty(smats, plate, fb, fmats, prim,
                                                                                       g, p).fn;
            if (fn > 0.0) any_force = true;
        }
        m.max_pen_over_d = std::max(m.max_pen_over_d, pen / d_f);
        m.end_pen_over_d = pen / d_f;
        if (any_contact)
        {
            ++m.contact_steps;
            if (!any_force) ++m.contact_steps_zero_fn;
        }

        // Запас 0.3 -- только на пределы с постоянной жёсткостью (связи, Герц).
        // Предел типа 3 равен 1.6*dt_prev/gamma и сам зависит от шага: запас
        // поверх него дал бы петлю dt -> 0.48 dt за шаг (§6 плана, «жёсткость
        // -> шаг»). В продовом цикле запаса поверх пределов нет.
        const double dt = std::min(safety * bond_dt, type == 0 ? safety * fs_dt : fs_dt);
        if (m.steps == 0) m.dt_bond0 = safety * bond_dt;
        if (fs_dt < bond_dt) ++m.steps_limited_by_fs;
        m.min_dt = std::min(m.min_dt, dt);
        if (dt_prev > 0.0) m.max_dt_ratio = std::max(m.max_dt_ratio, dt / dt_prev);

        // PC-драйвер, как FiberOnPlateIntegrator в test_e2e_fiber_solid.cpp.
        std::vector<double> dv_now[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_now[d].assign(N, 0.0);
            for (std::size_t i = 0; i < fb.size(); ++i) dv_now[d][i] = bond_delta[i].dv[d] + fs_delta[i].dv[d];
        }
        if (!first)
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < fb.size(); ++i)
                {
                    v_min[d][i] = fb.vx[d][i];
                    fb.vx[d][i] += 0.5 * dt * dv_prev[d][i];
                }
        for (std::size_t d = 0; d < dim; ++d) dv_prev[d] = dv_now[d];
        const double coeff = first ? 0.5 * dt : dt;
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < fb.size(); ++i)
            {
                const double base = first ? fb.vx[d][i] : v_min[d][i];
                fb.vx[d][i] = base + coeff * dv_prev[d][i];
                fb.x[d][i] += dt * fb.vx[d][i];
                if (!std::isfinite(fb.x[d][i]) || !std::isfinite(fb.vx[d][i])) m.finite = false;
            }
        first = false;
        dt_prev = dt;
        time += dt;
        ++m.steps;
        const double e = energy();
        m.e_max = std::max(m.e_max, e);
        while (next_sample < 4 && time >= 0.25 * (next_sample + 1) * t_end) m.pen_at[next_sample++] = pen / d_f;
        if (!m.finite) break;
    }
    m.e_end = energy();
    m.vz_end = fb.vx[2][center];
    return m;
}

void test_e2e_plate()
{
    std::printf("test_e2e_plate:\n");
    std::printf("    %-8s %6s %9s %9s %10s %7s %8s %9s %9s %9s %9s\n", "type", "steps", "max d/d", "end d/d",
                "min dt", "fs-lim", "fn=0", "E_end/E0", "E_max/E0", "vz_end", "dtN/dtP");
    for (int type = 0; type <= 3; ++type)
    {
        const E2eMetrics m = run_plate(type);
        std::printf("    %-8s %6d %9.3e %9.3e %10.3e %6.1f%% %7.1f%% %9.4f %9.4f %9.4f %9.4f\n", type_name(type),
                    m.steps, m.max_pen_over_d, m.end_pen_over_d, m.min_dt,
                    100.0 * m.steps_limited_by_fs / std::max(1, m.steps),
                    100.0 * m.contact_steps_zero_fn / std::max(1, m.contact_steps), m.e_end / m.e0, m.e_max / m.e0,
                    m.vz_end, m.max_dt_ratio);
        const std::string tag = type_name(type);
        check(m.finite, tag + ": нет NaN/Inf");
        check(m.max_pen_over_d > 0.0, tag + ": контакт произошёл");
        // Порог как у Hertz-сцены test_e2e_fiber_solid: энергия не растёт
        // больше чем на 1 % (у типов 1-3 ожидается потеря).
        check(m.e_max <= 1.01 * m.e0, tag + ": энергия не растёт больше 1 %");
        check(m.min_dt > 0.1 * m.dt_bond0, tag + ": шаг не схлопывается (min dt > 0.1 шага связей)");
        std::printf("      delta/d по времени (t/4, t/2, 3t/4, t): %.3e %.3e %.3e %.3e\n", m.pen_at[0], m.pen_at[1],
                    m.pen_at[2], m.pen_at[3]);
    }
}

} // namespace

int main()
{
    test_type0_unchanged();
    test_analytic_and_sign();
    test_momentum_lambdas();
    test_friction_slot();
    test_dt_ref_and_limit();
    test_parser();
    test_e2e_plate();
    std::printf("\n%s: %d failure(s)\n", failures == 0 ? "PASSED" : "FAILED", failures);
    return failures == 0 ? 0 : 1;
}
