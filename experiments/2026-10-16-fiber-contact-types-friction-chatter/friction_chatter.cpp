// friction_chatter.cpp -- диагностика дребезга силы трения нить-нить при законах
// нормальной силы 0-3 (Герц, pinball, Swegle, импульсный). Отчёт -- README.md
// рядом. Продовый код не меняется: сцена и PC-драйвер повторяют
// tests/test_e2e_fiber_contact.cpp (build_cross_scene, FiberBondContactIntegrator
// с полем law), сюда добавлены только запись по шагам и режим hold_n
// (контрфакт для гипотезы H1: на шагах с геометрическим контактом и fn = 0 слот
// трения заполняется последним ненулевым N этого primitive).
//
// Сборка из корня репозитория (см. README.md):
//   clang++ -std=c++20 -O2 -fopenmp -DSPH_DIM=3 -Iinclude <этот файл> -o friction_chatter.exe
// Запуск: friction_chatter.exe <каталог CSV>
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_contact_law.hpp"
#include "interactions/fiber_friction.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>
#include <vector>

static_assert(dim == 3, "нужна сборка с SPH_DIM=3");

namespace
{

// --- Сцена: копия build_cross_scene (tests/test_e2e_fiber_contact.cpp) -----
struct CrossScene
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    std::size_t n_p = 0, center_p = 0, center_q = 0;
    double d = 0.0;
};

void append_chain(fiber::FiberState& state, fiber::FiberTopology& topology, std::size_t base, int object_id,
                  const std::vector<std::array<double, 3>>& beads, double mass_each)
{
    const std::size_t n = beads.size();
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t gi = base + k;
        for (int d3 = 0; d3 < 3; ++d3) state.x[d3][gi] = beads[k][d3];
        for (std::size_t dd = 0; dd < dim; ++dd) state.vx[dd][gi] = 0.0;
        state.mass[gi] = mass_each;
        state.object_id[gi] = object_id;
        state.material_id[gi] = 0;
        if (k > 0) topology.prev[gi] = static_cast<long long>(gi - 1);
        if (k + 1 < n)
        {
            topology.next[gi] = static_cast<long long>(gi + 1);
            fiber::FiberBond b;
            b.i = gi; b.j = gi + 1;
            double l2 = 0.0;
            for (int d3 = 0; d3 < 3; ++d3)
            {
                const double dx = beads[k + 1][d3] - beads[k][d3];
                l2 += dx * dx;
            }
            b.l0 = std::sqrt(l2);
            b.l_break = b.l0 * 1e9;
            topology.bonds.push_back(b);
        }
    }
}

CrossScene build_cross_scene(double mu_ff)
{
    CrossScene sc;
    const double E = 2.0e7, nu = 0.3, d = 1.0e-3, rho = 1200.0, l0 = 2.0e-3;
    const double mass = rho * (pi * d * d / 4.0) * l0;
    const double z0 = 1.0e-3, v0 = 0.5;
    const int N = 7;
    std::vector<std::array<double, 3>> P, Q;
    for (int k = -(N / 2); k <= N / 2; ++k)
    {
        P.push_back({static_cast<double>(k) * l0, 0.0, z0});
        Q.push_back({0.0, static_cast<double>(k) * l0, -z0});
    }
    fiber::FiberMaterial m;
    m.E = E; m.nu = nu; m.d = d; m.rho = rho; m.eps_break = 1.0; m.mu_ff = mu_ff;
    sc.materials.push_back(m);
    sc.state.resize(P.size() + Q.size());
    sc.topology.resize(sc.state.size());
    append_chain(sc.state, sc.topology, 0, 0, P, mass);
    append_chain(sc.state, sc.topology, P.size(), 1, Q, mass);
    sc.n_p = P.size();
    sc.center_p = static_cast<std::size_t>(N / 2);
    sc.center_q = P.size() + static_cast<std::size_t>(N / 2);
    sc.d = d;
    for (std::size_t i = 0; i < P.size(); ++i) sc.state.vx[2][i] = -v0;
    for (std::size_t i = 0; i < Q.size(); ++i) sc.state.vx[2][P.size() + i] = v0;
    return sc;
}

// --- Запись шага ------------------------------------------------------------
struct StepRec
{
    int step = 0;
    double t = 0.0, dt = 0.0;
    std::size_t n_geom = 0;   // принятые геометрические контакты
    std::size_t n_force = 0;  // из них с fn > 0 (до hold)
    std::size_t n_held = 0;   // слотов, заполненных удержанием N (hold_n)
    double sumN = 0.0;        // сумма N по активным слотам (после hold)
    double imp_q = 0.0;       // импульс трения на Q по y за шаг (сила трения * h)
    double cap = 0.0;         // mu * sumN * h -- предел импульса Кулона за шаг (масштаб)
    double ut = 0.0;          // max ||u_t|| точки контакта после шага
    std::size_t stick = 0, slip = 0;
    double vy_rel = 0.0;      // vy центра P - vy центра Q
    int branch = -1;          // pinball, контакт с max N: 0 -- вязкая F1, 1 -- герцевская F2; -1 -- нет
};

// --- PC-драйвер: копия FiberBondContactIntegrator::step (тот же порядок) ----
struct Driver
{
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> bond_buffers, contact_buffers, friction_buffers;
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> count_buffers;
    std::vector<double> bead_count;
    std::vector<interactions::FiberFrictionContact> slots;
    std::vector<double> dv_prev[dim], v_min[dim];
    std::vector<double> bead_dt_ref;
    bool first = true;
    double dt_prev = 0.0;
    double safety = 0.5;
    interactions::FiberContactLaw law;
    bool hold_n = false;
    std::map<std::tuple<long long, long long, long long, long long>, double> last_n; // hold_n

    StepRec step(fiber::FiberState& fb, const fiber::FiberTopology& topo, const std::vector<fiber::FiberMaterial>& mats,
                 std::size_t q_first, std::size_t q_count)
    {
        StepRec r;
        const std::size_t n = fb.size();
        if (dv_prev[0].size() != n)
            for (std::size_t d = 0; d < dim; ++d) { dv_prev[d].assign(n, 0.0); v_min[d].assign(n, 0.0); }

        if (!first)
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i)
                {
                    v_min[d][i] = fb.vx[d][i];
                    fb.vx[d][i] += 0.5 * dt_prev * dv_prev[d][i];
                }

        const auto bond_delta = interactions::compute_fiber_bond(exec, fb, topo, mats, bond_buffers);
        const double bond_dt = interactions::compute_fiber_bond_dt_limit(exec, fb, topo, mats);

        interactions::FiberContactSearchResult search;
        interactions::build_fiber_contact_search(exec, fb, topo, mats, search);
        interactions::FiberContactStepParams params;
        if (law.type != 0)
        {
            params.law = law;
            params.dt_prev = dt_prev > 0.0 ? dt_prev : 1e-7;
            interactions::compute_fiber_bead_dt_ref(fb, topo, mats, bead_dt_ref);
            params.bead_dt_ref = &bead_dt_ref;
            params.pair_dt_limit = &search.pair_dt_limit;
        }
        bool has_friction = interactions::compute_fiber_contact_contributions(exec, fb, mats, search.primitives,
                                                                             contact_buffers, slots, params);
        const auto contact_delta = rhs::reduce_contributions(exec, contact_buffers, n);
        const double contact_dt = interactions::compute_fiber_contact_dt_limit(exec, fb, mats, search.primitives, params);

        // Геометрия, fn > 0, ветка pinball, hold_n.
        double best_n = -1.0;
        for (std::size_t k = 0; k < search.primitives.size(); ++k)
        {
            const auto& prim = search.primitives[k];
            const auto geom = interactions::compute_fiber_contact_geometry(fb, mats, prim);
            if (!geom.contact) continue;
            ++r.n_geom;
            const auto key = std::make_tuple(static_cast<long long>(prim.a0), prim.a1,
                                             static_cast<long long>(prim.b0), prim.b1);
            if (slots[k].active)
            {
                ++r.n_force;
                last_n[key] = slots[k].N;
                if (law.type == 1 && slots[k].N > best_n)
                {
                    best_n = slots[k].N;
                    contact::PairGeom g;
                    contact::PairProps pr;
                    interactions::fiber_fiber_pair_inputs(fb, mats, prim, geom, params, g, pr);
                    const double f1 = pr.m_red * g.pdot / pr.dt_ref;
                    const double f2 = pr.g_eff * std::sqrt(pr.r_eff) * g.p * std::sqrt(g.p);
                    r.branch = (f2 < f1) ? 1 : 0;
                }
            }
            else if (hold_n)
            {
                const auto it = last_n.find(key);
                if (it == last_n.end()) continue;
                auto& s = slots[k];
                s.active = true;
                s.a0 = prim.a0; s.a1 = prim.a1; s.b0 = prim.b0; s.b1 = prim.b1;
                s.s = geom.s; s.t = geom.t;
                for (std::size_t d = 0; d < dim; ++d) s.n[d] = geom.n[d];
                s.N = it->second;
                s.mu = interactions::fiber_friction_coefficient(mats[0], mats[0]);
                ++r.n_held;
                has_friction = true;
            }
        }

        const double dt = (law.type == 0) ? safety * std::min(bond_dt, contact_dt)
                                          : std::min(safety * bond_dt, contact_dt);

        std::vector<double> dv_now[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_now[d].assign(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) dv_now[d][i] = bond_delta[i].dv[d] + contact_delta[i].dv[d];
        }
        const double coeff = first ? dt * 0.5 : dt;

        if (has_friction)
        {
            interactions::compute_fiber_friction_bead_counts(exec, n, slots, count_buffers, bead_count);
            ode::FiberVelocityUpdate upd{};
            for (std::size_t d = 0; d < dim; ++d) upd.base[d] = first ? &fb.vx[d] : &v_min[d];
            upd.h = coeff;
            const auto fd = interactions::compute_fiber_friction_diagnostics(exec, fb, upd, dv_now, slots, bead_count);
            r.stick = fd.stick_count;
            r.slip = fd.slip_count;
            for (const auto& c : slots)
                if (c.active) r.sumN += c.N;
            for (const auto& c : slots)
                if (c.active) { r.cap = c.mu * r.sumN * coeff; break; }
            interactions::compute_fiber_friction_contributions(exec, fb, upd, dv_now, slots, bead_count,
                                                               friction_buffers);
            const auto fdel = rhs::reduce_contributions(exec, friction_buffers, n);
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i) dv_now[d][i] += fdel[i].dv[d];
            for (std::size_t i = q_first; i < q_first + q_count; ++i) r.imp_q += fb.mass[i] * coeff * fdel[i].dv[1];
        }

        for (std::size_t d = 0; d < dim; ++d) dv_prev[d] = dv_now[d];
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < n; ++i)
            {
                const double base = first ? fb.vx[d][i] : v_min[d][i];
                fb.vx[d][i] = base + coeff * dv_prev[d][i];
                fb.x[d][i] += dt * fb.vx[d][i];
            }

        // ||u_t|| по итоговым скоростям, max по активным слотам.
        for (const auto& c : slots)
        {
            if (!c.active) continue;
            double u[dim], un = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                double vA = fb.vx[d][c.a0];
                if (c.a1 != -1) vA = (1.0 - c.s) * vA + c.s * fb.vx[d][static_cast<std::size_t>(c.a1)];
                double vB = fb.vx[d][c.b0];
                if (c.b1 != -1) vB = (1.0 - c.t) * vB + c.t * fb.vx[d][static_cast<std::size_t>(c.b1)];
                u[d] = vA - vB;
                un += u[d] * c.n[d];
            }
            double ut2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d) ut2 += (u[d] - un * c.n[d]) * (u[d] - un * c.n[d]);
            r.ut = std::max(r.ut, std::sqrt(ut2));
        }

        first = false;
        dt_prev = dt;
        r.dt = dt;
        return r;
    }
};

// --- Прогон и сводка --------------------------------------------------------
struct Summary
{
    int steps = 0, contact_steps = 0, zero_fn_steps = 0;
    double contact_time = 0.0;
    int sgn_v = 0;            // смены знака vy_rel на шагах контакта (как в тесте)
    double amp_v = 0.0;       // max |vy_rel| после первой смены
    int sgn_f = 0;            // смены знака импульса трения на Q
    int sgn_f_after_zero = 0; // из них -- на шаге сразу после шага с fn = 0 (или на нём)
    double amp_f = 0.0;       // max |imp_q|/cap после первой смены знака силы
    double ut_max_after = 0.0;
    double n_cv = 0.0;        // коэффициент вариации sumN по шагам с силой
    int branch_switches = 0;
    int slip_steps = 0, stick_steps = 0;
    double imp_total = 0.0;
};

Summary run(int type, double mu, double v_t0, bool hold, int max_contact_steps, const std::string& csv)
{
    CrossScene sc = build_cross_scene(mu);
    for (std::size_t i = 0; i < sc.n_p; ++i) sc.state.vx[1][i] = v_t0;
    Driver drv;
    drv.law.type = type;
    drv.hold_n = hold;
    FILE* f = csv.empty() ? nullptr : std::fopen(csv.c_str(), "w");
    if (f) std::fprintf(f, "step,t,dt,n_geom,n_force,n_held,sumN,imp_q,cap,ut,stick,slip,vy_rel,branch\n");

    Summary s;
    double t = 0.0;
    int last_sv = 0, last_sf = 0, last_branch = -1;
    bool first_v = false, first_f = false, prev_zero = false;
    double n_sum = 0.0, n_sum2 = 0.0;
    int n_cnt = 0;
    while (t < 8.0e-3 && s.steps < 40000)
    {
        StepRec r = drv.step(sc.state, sc.topology, sc.materials, sc.n_p, sc.n_p);
        t += r.dt;
        r.t = t;
        r.step = s.steps++;
        r.vy_rel = sc.state.vx[1][sc.center_p] - sc.state.vx[1][sc.center_q];
        if (f)
            std::fprintf(f, "%d,%.9e,%.6e,%zu,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%zu,%zu,%.6e,%d\n", r.step, r.t, r.dt,
                         r.n_geom, r.n_force, r.n_held, r.sumN, r.imp_q, r.cap, r.ut, r.stick, r.slip, r.vy_rel,
                         r.branch);
        if (r.n_geom == 0) continue;
        if (max_contact_steps > 0 && s.contact_steps >= max_contact_steps) continue;
        ++s.contact_steps;
        s.contact_time += r.dt;
        const bool zero = (r.n_force == 0);
        if (zero) ++s.zero_fn_steps;
        s.imp_total += r.imp_q;
        if (r.slip > 0) ++s.slip_steps;
        if (r.stick > 0 && r.slip == 0) ++s.stick_steps;

        const int sv = (r.vy_rel > 0.0) - (r.vy_rel < 0.0);
        if (sv != 0 && last_sv != 0 && sv != last_sv) { ++s.sgn_v; first_v = true; }
        if (sv != 0) last_sv = sv;
        if (first_v) s.amp_v = std::max(s.amp_v, std::fabs(r.vy_rel));

        const int sf = (r.imp_q > 0.0) - (r.imp_q < 0.0);
        if (sf != 0 && last_sf != 0 && sf != last_sf)
        {
            ++s.sgn_f;
            first_f = true;
            if (prev_zero) ++s.sgn_f_after_zero;
        }
        if (sf != 0) last_sf = sf;
        if (first_f && r.cap > 0.0) s.amp_f = std::max(s.amp_f, std::fabs(r.imp_q) / r.cap);
        if (first_f) s.ut_max_after = std::max(s.ut_max_after, r.ut);
        // «после fn = 0»: этот шаг без силы или любой из двух предыдущих.
        prev_zero = zero || (sf == 0 && prev_zero);

        if (!zero) { n_sum += r.sumN; n_sum2 += r.sumN * r.sumN; ++n_cnt; }
        if (r.branch >= 0 && last_branch >= 0 && r.branch != last_branch) ++s.branch_switches;
        if (r.branch >= 0) last_branch = r.branch;
    }
    if (n_cnt > 1)
    {
        const double mean = n_sum / n_cnt;
        const double var = std::max(0.0, n_sum2 / n_cnt - mean * mean);
        s.n_cv = mean > 0.0 ? std::sqrt(var) / mean : 0.0;
    }
    if (f) std::fclose(f);
    return s;
}

void print_header()
{
    std::printf("%-8s %-5s %6s %5s %8s %7s %5s %9s %5s %7s %9s %9s %6s %6s %6s %6s %10s\n", "type", "hold", "v_t0",
                "Kstp", "cont_us", "fn0%", "sgnV", "ampV", "sgnF", "sgnF0", "ampF/cap", "ut_aft", "N_cv", "brSw",
                "slip", "stick", "imp_tot");
}

void print_row(const char* name, bool hold, double v_t0, int K, const Summary& s)
{
    std::printf("%-8s %-5s %6.2f %5d %8.1f %6.1f%% %5d %9.2e %5d %7d %9.3f %9.2e %6.3f %6d %6d %6d %10.3e\n", name,
                hold ? "yes" : "no", v_t0, K, s.contact_time * 1e6,
                100.0 * s.zero_fn_steps / std::max(1, s.contact_steps), s.sgn_v, s.amp_v, s.sgn_f, s.sgn_f_after_zero,
                s.amp_f, s.ut_max_after, s.n_cv, s.branch_switches, s.slip_steps, s.stick_steps, s.imp_total);
}

} // namespace

int main(int argc, char** argv)
{
    const std::string out = argc > 1 ? argv[1] : "";
    static const char* names[] = {"Hertz", "pinball", "Swegle", "impulse"};
    const double mu = 0.2;

    // Шагов контакта у Герца (окно для сравнения при равной экспозиции, H4).
    int hertz_contact_steps = 0;
    {
        const Summary h = run(0, mu, 0.05, false, 0, "");
        hertz_contact_steps = h.contact_steps;
    }
    std::printf("Окно равной экспозиции: K = %d шагов контакта (Герц, v_t0 = 0.05)\n\n", hertz_contact_steps);

    print_header();
    for (const double v_t0 : {0.05, 2.0})
    {
        for (int type = 0; type <= 3; ++type)
        {
            for (int hold = 0; hold <= (type == 0 ? 0 : 1); ++hold)
            {
                char csv[512] = "";
                if (!out.empty())
                    std::snprintf(csv, sizeof csv, "%s/%s_v%.2f%s.csv", out.c_str(), names[type], v_t0,
                                  hold ? "_hold" : "");
                const Summary full = run(type, mu, v_t0, hold != 0, 0, csv);
                const Summary win = run(type, mu, v_t0, hold != 0, hertz_contact_steps, "");
                print_row(names[type], hold != 0, v_t0, 0, full);
                print_row(names[type], hold != 0, v_t0, hertz_contact_steps, win);
            }
        }
        std::printf("\n");
    }
    std::printf("Столбцы: Kstp -- окно шагов контакта (0 -- весь прогон); cont_us -- время в контакте, мкс;\n"
                "fn0%% -- доля шагов контакта, где ни у одного контакта нет fn > 0; sgnV/ampV -- смены знака\n"
                "vy_rel центров и max |vy_rel| после первой смены, м/с; sgnF -- смены знака импульса трения на Q;\n"
                "sgnF0 -- из них на шаге с fn = 0 или сразу после; ampF/cap -- max |импульс|/(mu*sumN*h) после\n"
                "первой смены; ut_aft -- max ||u_t|| после первой смены, м/с; N_cv -- коэф. вариации sumN;\n"
                "brSw -- переключения ветви pinball F1<->F2; slip/stick -- шаги со скольжением / только прилипание;\n"
                "imp_tot -- суммарный импульс трения на Q, кг*м/с.\n");
    return 0;
}
