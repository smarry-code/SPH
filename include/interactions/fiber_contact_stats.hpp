#pragma once
// interactions/fiber_contact_stats.hpp
// Диагностическая статистика нормальной силы и трения контактов нитей
// (Fiber-Fiber и Fiber-Solid) для проверки гипотез о законах pinball / Swegle
// (contact/12-fiber-contact-types.md; experiments/*-guo2019-1t-contact-types-friction).
// Только в сборке с -DSPH_CONTACT_STATS=ON (src/time_integration.cpp).
//
// Без hold_n модуль только ЧИТАЕТ состояние: траектория побитово та же, что
// без статистики. hold_n = true -- контрфакт «удержание N»: если у контакта есть
// геометрический контакт, а закон дал fn = 0, слот трения заполняется последним
// ненулевым N этого контакта и текущей геометрией (меняет физику, только для
// диагностики).
//
// Что считается (по шагам и итог в конце):
//   G1  контакты с fn = 0 среди геометрических: по шагам и с весом силы Герца
//       той же геометрии (geom.F);
//   G2  прилипание/скольжение активных слотов (||u_t|| по итоговым скоростям
//       после шага), «потерянные» слоты (на прошлом шаге активен, сейчас fn = 0);
//   G3  N против силы Герца и скорости проникания (корреляции), доля вязкой
//       ветки pinball (F1 < F2) и потолка Swegle;
//   G4  смены направления u_t одного контакта между шагами (u_t . u_t_prev < 0),
//       длины серий одного направления, смены на шаге с fn = 0 или после него;
//   G5  max delta/d.
// Ключ контакта: Fiber-Fiber (a0, a1, b0, b1), Fiber-Solid (a, b, s).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "contact_force.hpp"
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_contact_law.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "materials.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"

namespace interactions
{
namespace contact_stats
{

using Key = std::array<long long, 4>;
struct KeyHash
{
    std::size_t operator()(const Key& k) const noexcept
    {
        std::uint64_t h = 1469598103934665603ull;
        for (long long v : k)
        {
            h ^= static_cast<std::uint64_t>(v) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        }
        return static_cast<std::size_t>(h);
    }
};

struct KeyState
{
    long long last_step = -10;       // последний шаг с геометрическим контактом
    long long last_zero_step = -10;  // последний шаг с fn = 0 при геометрическом контакте
    long long last_active_step = -10;
    long long run_start = -1;        // начало текущей серии одного направления u_t
    double last_n = 0.0;             // последнее ненулевое N (hold_n)
    double ut[3] = {0.0, 0.0, 0.0};  // u_t прошлого активного шага
    bool has_ut = false;
    long long contact_steps = 0;
};

// Агрегаты одного вида контакта.
struct Side
{
    // по шагу
    std::size_t n_geom = 0, n_force = 0, n_zero = 0, n_held = 0, n_lost = 0;
    std::size_t n_active = 0, n_stick = 0, n_near = 0, n_slip = 0, n_flip = 0, n_flip_after_zero = 0;
    std::size_t n_branch_visc = 0, n_capped = 0;
    double sum_n = 0.0, sum_n_hertz = 0.0, sum_n_hertz_zero = 0.0, max_pen_d = 0.0;
    void reset_step()
    {
        *this = Side{};
    }
};

struct Totals
{
    double contact_steps = 0, force_steps = 0, zero_steps = 0, held = 0, lost = 0;
    double active = 0, stick = 0, near = 0, slip = 0, flips = 0, flips_after_zero = 0;
    double branch_visc = 0, capped = 0;
    double sum_n = 0, sum_n_hertz = 0, sum_n_hertz_zero = 0;
    double max_pen_d = 0;
    // корреляции по контакт-шагам с fn > 0: (N, N_hertz) и (N, pdot)
    double sx = 0, sy = 0, sz = 0, sxx = 0, syy = 0, szz = 0, sxy = 0, sxz = 0, cnt = 0;
    // N_cv по шагам (сумма N шага)
    double step_n = 0, step_n2 = 0, step_cnt = 0;
    std::array<double, 6> run_hist{}; // длины серий: 1, 2, 3-5, 6-10, 11-30, >30
    double keys = 0, key_steps = 0;

    void add(const Side& s)
    {
        contact_steps += static_cast<double>(s.n_geom);
        force_steps += static_cast<double>(s.n_force);
        zero_steps += static_cast<double>(s.n_zero);
        held += static_cast<double>(s.n_held);
        lost += static_cast<double>(s.n_lost);
        branch_visc += static_cast<double>(s.n_branch_visc);
        capped += static_cast<double>(s.n_capped);
        sum_n += s.sum_n;
        sum_n_hertz += s.sum_n_hertz;
        sum_n_hertz_zero += s.sum_n_hertz_zero;
        max_pen_d = std::max(max_pen_d, s.max_pen_d);
        if (s.n_force > 0)
        {
            step_n += s.sum_n;
            step_n2 += s.sum_n * s.sum_n;
            step_cnt += 1.0;
        }
    }
    void add_after(const Side& s)
    {
        active += static_cast<double>(s.n_active);
        stick += static_cast<double>(s.n_stick);
        near += static_cast<double>(s.n_near);
        slip += static_cast<double>(s.n_slip);
        flips += static_cast<double>(s.n_flip);
        flips_after_zero += static_cast<double>(s.n_flip_after_zero);
    }
    void corr(double n, double nh, double pdot)
    {
        sx += n; sy += nh; sz += pdot;
        sxx += n * n; syy += nh * nh; szz += pdot * pdot;
        sxy += n * nh; sxz += n * pdot;
        cnt += 1.0;
    }
    void run(long long len)
    {
        const std::size_t b = len <= 1 ? 0 : len == 2 ? 1 : len <= 5 ? 2 : len <= 10 ? 3 : len <= 30 ? 4 : 5;
        run_hist[b] += 1.0;
    }
};

inline double pearson(double n, double sx, double sy, double sxx, double syy, double sxy)
{
    if (n < 2.0) return 0.0;
    const double vx = sxx - sx * sx / n, vy = syy - sy * sy / n;
    if (!(vx > 0.0) || !(vy > 0.0)) return 0.0;
    return (sxy - sx * sy / n) / std::sqrt(vx * vy);
}

class Collector
{
public:
    Collector(const SimulationConfig& cfg, std::size_t every, std::string csv_path, std::string summary_path, bool hold_n)
        : cfg_(cfg), every_(every), summary_path_(std::move(summary_path)), hold_n_(hold_n)
    {
        csv_ = std::fopen(csv_path.c_str(), "w");
        if (csv_)
        {
            std::fprintf(csv_, "step,t,dt");
            for (const char* p : {"ff", "fs"})
                std::fprintf(csv_,
                             ",%s_geom,%s_force,%s_zero,%s_held,%s_lost,%s_active,%s_stick,%s_slip,%s_flip,%s_flip0,"
                             "%s_visc,%s_capped,%s_sumN,%s_sumNH,%s_sumNH_zero,%s_pen_d",
                             p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p);
            std::fprintf(csv_, "\n");
        }
    }
    Collector(const Collector&) = delete;
    Collector& operator=(const Collector&) = delete;
    ~Collector()
    {
        if (csv_) std::fclose(csv_);
        write_summary();
    }

    // После RHS обоих контактов, ДО стадии трения. При hold_n дополняет слоты;
    // возвращает true, если какой-либо слот стал активным (флаги трения шага
    // нужно поднять).
    bool after_rhs(long long step, const Materials& mats, const ParticleSystem& sys, const fiber::FiberState& fb,
                   const fiber::FiberTopology& topo, const std::vector<fiber::FiberMaterial>& materials,
                   const FiberContactSearchResult& ff_search, std::vector<FiberFrictionContact>& ff_slots,
                   const FiberContactLaw& ff_law, const FiberSolidContactSearchResult& fs_search,
                   std::vector<FiberSolidFrictionContact>& fs_slots, const FiberContactLaw& fs_law, double dt_prev,
                   bool& ff_any, bool& fs_any)
    {
        step_ = step;
        ff_.reset_step();
        fs_.reset_step();
        compute_fiber_bead_dt_ref(fb, topo, materials, bead_dt_ref_);
        FiberContactStepParams pff, pfs;
        pff.law = ff_law; pff.dt_prev = dt_prev > 0.0 ? dt_prev : 1e-12; pff.cfg = &cfg_; pff.bead_dt_ref = &bead_dt_ref_;
        pfs.law = fs_law; pfs.dt_prev = pff.dt_prev; pfs.cfg = &cfg_; pfs.bead_dt_ref = &bead_dt_ref_;

        // --- Fiber-Fiber ---------------------------------------------------
        for (std::size_t k = 0; k < ff_search.primitives.size(); ++k)
        {
            const FiberContactPrimitive& prim = ff_search.primitives[k];
            const FiberContactGeometry geom = compute_fiber_contact_geometry(fb, materials, prim);
            if (!geom.contact) continue;
            const Key key{static_cast<long long>(prim.a0), prim.a1, static_cast<long long>(prim.b0), prim.b1};
            KeyState& ks = ff_keys_[key];
            const double d = materials[static_cast<std::size_t>(fb.material_id[prim.a0])].d;
            contact::PairGeom g;
            contact::PairProps pr;
            fiber_fiber_pair_inputs(fb, materials, prim, geom, pff, g, pr);
            FiberFrictionContact& slot = ff_slots[k];
            // fn > 0 определяется законом, а не слотом (слот есть только при mu > 0).
            const double fn = (ff_law.type == 0) ? geom.F : fiber_fiber_contact_penalty(fb, materials, prim, geom, pff).fn;
            account(ff_, ff_tot_, ks, geom.F, fn, g, pr, ff_law.type, geom.delta / d);
            if (fn > 0.0)
            {
                ks.last_n = fn;
            }
            else if (hold_n_ && ks.last_n > 0.0 && !slot.active)
            {
                const double mu = fiber_friction_coefficient(materials[static_cast<std::size_t>(fb.material_id[prim.a0])],
                                                             materials[static_cast<std::size_t>(fb.material_id[prim.b0])]);
                if (mu > 0.0)
                {
                    slot.active = true;
                    slot.a0 = prim.a0; slot.a1 = prim.a1; slot.b0 = prim.b0; slot.b1 = prim.b1;
                    slot.s = geom.s; slot.t = geom.t;
                    for (std::size_t dd = 0; dd < dim; ++dd) slot.n[dd] = geom.n[dd];
                    slot.N = ks.last_n;
                    slot.mu = mu;
                    ++ff_.n_held;
                    ff_any = true;
                }
            }
        }

        // --- Fiber-Solid ---------------------------------------------------
        for (std::size_t k = 0; k < fs_search.primitives.size(); ++k)
        {
            const FiberSolidContactPrimitive& prim = fs_search.primitives[k];
            const FiberSolidContactGeometry geom = compute_fiber_solid_contact_geometry(mats, sys, fb, materials, prim);
            if (!geom.contact) continue;
            const Key key{static_cast<long long>(prim.a), prim.has_b ? static_cast<long long>(prim.b) : -1,
                          static_cast<long long>(prim.s), -2};
            KeyState& ks = fs_keys_[key];
            const fiber::FiberMaterial& mat_a = materials[static_cast<std::size_t>(fb.material_id[prim.a])];
            contact::PairGeom g;
            contact::PairProps pr;
            fiber_solid_pair_inputs(mats, sys, fb, materials, prim, geom, pfs, g, pr);
            const double fn =
                (fs_law.type == 0) ? geom.F : fiber_solid_contact_penalty(mats, sys, fb, materials, prim, geom, pfs).fn;
            account(fs_, fs_tot_, ks, geom.F, fn, g, pr, fs_law.type, geom.delta / mat_a.d);
            FiberSolidFrictionContact& slot = fs_slots[k];
            if (fn > 0.0)
            {
                ks.last_n = fn;
            }
            else if (hold_n_ && ks.last_n > 0.0 && !slot.active)
            {
                const double mu = fiber_solid_friction_coefficient(mat_a, sys.itype[prim.s]);
                if (mu > 0.0)
                {
                    slot.active = true;
                    slot.a = prim.a;
                    slot.b = prim.has_b ? static_cast<long long>(prim.b) : -1;
                    slot.s = prim.s;
                    slot.lambda = prim.has_b ? geom.lambda : 0.0;
                    for (std::size_t dd = 0; dd < dim; ++dd) slot.n[dd] = geom.n[dd];
                    slot.N = ks.last_n;
                    slot.mu = mu;
                    ++fs_.n_held;
                    fs_any = true;
                }
            }
        }
        ff_tot_.add(ff_);
        fs_tot_.add(fs_);
        return ff_.n_held > 0 || fs_.n_held > 0;
    }

    // После solver.advance: u_t активных слотов по итоговым скоростям.
    void after_advance(double time, double dt, const fiber::FiberState& fb, const ParticleSystem& sys,
                       const std::vector<FiberFrictionContact>& ff_slots,
                       const std::vector<FiberSolidFrictionContact>& fs_slots)
    {
        for (const auto& c : ff_slots)
        {
            if (!c.active) continue;
            double u[3] = {0.0, 0.0, 0.0};
            for (std::size_t d = 0; d < dim; ++d)
            {
                double vA = fb.vx[d][c.a0];
                if (c.a1 != -1) vA = (1.0 - c.s) * vA + c.s * fb.vx[d][static_cast<std::size_t>(c.a1)];
                double vB = fb.vx[d][c.b0];
                if (c.b1 != -1) vB = (1.0 - c.t) * vB + c.t * fb.vx[d][static_cast<std::size_t>(c.b1)];
                u[d] = vA - vB;
            }
            const Key key{static_cast<long long>(c.a0), c.a1, static_cast<long long>(c.b0), c.b1};
            tangential(ff_, ff_tot_, ff_keys_[key], u, c.n);
        }
        for (const auto& c : fs_slots)
        {
            if (!c.active) continue;
            double u[3] = {0.0, 0.0, 0.0};
            for (std::size_t d = 0; d < dim; ++d)
            {
                double vq = fb.vx[d][c.a];
                if (c.b != -1) vq = (1.0 - c.lambda) * vq + c.lambda * fb.vx[d][static_cast<std::size_t>(c.b)];
                u[d] = vq - sys.vx[d][c.s];
            }
            const Key key{static_cast<long long>(c.a), c.b, static_cast<long long>(c.s), -2};
            tangential(fs_, fs_tot_, fs_keys_[key], u, c.n);
        }
        ff_tot_.add_after(ff_);
        fs_tot_.add_after(fs_);
        ++steps_;
        if (csv_ && every_ > 0 && (step_ % static_cast<long long>(every_)) == 0)
        {
            std::fprintf(csv_, "%lld,%.9e,%.6e", step_, time, dt);
            for (const Side* s : {&ff_, &fs_})
                std::fprintf(csv_, ",%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%.6e,%.6e,%.6e,%.6e", s->n_geom,
                             s->n_force, s->n_zero, s->n_held, s->n_lost, s->n_active, s->n_stick, s->n_slip,
                             s->n_flip, s->n_flip_after_zero, s->n_branch_visc, s->n_capped, s->sum_n, s->sum_n_hertz,
                             s->sum_n_hertz_zero, s->max_pen_d);
            std::fprintf(csv_, "\n");
        }
    }

private:
    // Прилипание: ||u_t|| ниже порога (проекция гасит u_t до ~1e-16 относительно).
    static constexpr double kStickUt = 1e-9; // м/с
    // «Почти прилипание»: связанные контакты один проход Якоби гасит не до нуля.
    static constexpr double kNearUt = 1e-3; // м/с

    void account(Side& s, Totals& tot, KeyState& ks, double n_hertz, double fn, const contact::PairGeom& g,
                 const contact::PairProps& pr, int type, double pen_d)
    {
        ++s.n_geom;
        s.sum_n_hertz += n_hertz;
        s.max_pen_d = std::max(s.max_pen_d, pen_d);
        ++ks.contact_steps;
        if (fn > 0.0)
        {
            ++s.n_force;
            s.sum_n += fn;
            tot.corr(fn, n_hertz, g.pdot);
            if (type == 1 && pr.dt_ref > 0.0)
            {
                const double f1 = pr.m_red * g.pdot / pr.dt_ref;
                const double f2 = (pr.g_eff > 0.0 && pr.r_eff > 0.0) ? pr.g_eff * std::sqrt(pr.r_eff) * g.p * std::sqrt(g.p)
                                                                      : 1e300;
                if (f1 < f2) ++s.n_branch_visc;
            }
            if (type == 2 && pr.dt_ref > 0.0 && g.r > 0.0)
            {
                const double f1 = pr.m_red * g.p / (pr.dt_ref * pr.dt_ref);
                const double f2 = 2.0 * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r);
                const double cap = pr.m_red * g.pdot / pr.dt_ref;
                if (std::max(f1, f2) > cap) ++s.n_capped;
            }
        }
        else
        {
            ++s.n_zero;
            s.sum_n_hertz_zero += n_hertz;
            ks.last_zero_step = step_;
            if (ks.last_active_step == step_ - 1) ++s.n_lost;
        }
        ks.last_step = step_;
    }

    void tangential(Side& s, Totals& tot, KeyState& ks, const double u[3], const double n[dim])
    {
        ++s.n_active;
        double un = 0.0;
        for (std::size_t d = 0; d < dim; ++d) un += u[d] * n[d];
        double ut[3] = {0.0, 0.0, 0.0}, ut2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            ut[d] = u[d] - un * n[d];
            ut2 += ut[d] * ut[d];
        }
        const double nut = std::sqrt(ut2);
        if (nut < kNearUt) ++s.n_near;
        if (nut < kStickUt)
        {
            ++s.n_stick;
        }
        else
        {
            ++s.n_slip;
            if (ks.has_ut && ks.last_active_step == step_ - 1)
            {
                double dot = 0.0;
                for (std::size_t d = 0; d < dim; ++d) dot += ut[d] * ks.ut[d];
                if (dot < 0.0)
                {
                    ++s.n_flip;
                    if (ks.last_zero_step >= step_ - 2) ++s.n_flip_after_zero;
                    if (ks.run_start >= 0) tot.run(step_ - ks.run_start);
                    ks.run_start = step_;
                }
            }
            if (ks.run_start < 0 || ks.last_active_step != step_ - 1) ks.run_start = step_;
            for (std::size_t d = 0; d < dim; ++d) ks.ut[d] = ut[d];
            ks.has_ut = true;
        }
        ks.last_active_step = step_;
    }

    void write_side(FILE* f, const char* name, const Totals& t,
                    const std::unordered_map<Key, KeyState, KeyHash>& keys) const
    {
        double key_steps = 0.0;
        for (const auto& kv : keys) key_steps += static_cast<double>(kv.second.contact_steps);
        const double nk = static_cast<double>(keys.size());
        auto frac = [](double a, double b) { return b > 0.0 ? a / b : 0.0; };
        const double mean = frac(t.step_n, t.step_cnt);
        const double var = std::max(0.0, frac(t.step_n2, t.step_cnt) - mean * mean);
        double runs = 0.0;
        for (double v : t.run_hist) runs += v;
        std::fprintf(f, "[%s]\n", name);
        std::fprintf(f, "contact_steps = %.0f   (контакт-шагов: геометрический контакт)\n", t.contact_steps);
        std::fprintf(f, "G1 zero_fn_frac = %.6f   (доля контакт-шагов с fn = 0)\n", frac(t.zero_steps, t.contact_steps));
        std::fprintf(f, "G1 zero_fn_frac_hertz_weighted = %.6f   (доля силы Герца на контактах с fn = 0)\n",
                     frac(t.sum_n_hertz_zero, t.sum_n_hertz));
        std::fprintf(f, "G1 held_slots = %.0f   (заполнено удержанием N)\n", t.held);
        std::fprintf(f,
                     "G2 active_slots = %.0f, stick_frac = %.6f (|u_t| < 1e-9), near_stick_frac = %.6f (|u_t| < 1e-3), "
                     "slip_frac = %.6f\n",
                     t.active, frac(t.stick, t.active), frac(t.near, t.active), frac(t.slip, t.active));
        std::fprintf(f, "G2 lost_slots = %.0f   (активен на прошлом шаге, сейчас fn = 0), lost_per_active = %.6f\n",
                     t.lost, frac(t.lost, t.active));
        std::fprintf(f, "G2 sumN / sumN_hertz = %.6f\n", frac(t.sum_n, t.sum_n_hertz));
        std::fprintf(f, "G3 corr(N, N_hertz) = %.4f, corr(N, pdot) = %.4f   (по контакт-шагам с fn > 0, n = %.0f)\n",
                     pearson(t.cnt, t.sx, t.sy, t.sxx, t.syy, t.sxy), pearson(t.cnt, t.sx, t.sz, t.sxx, t.szz, t.sxz),
                     t.cnt);
        std::fprintf(f, "G3 pinball_visc_frac = %.6f, swegle_capped_frac = %.6f   (среди контакт-шагов с fn > 0)\n",
                     frac(t.branch_visc, t.force_steps), frac(t.capped, t.force_steps));
        std::fprintf(f, "G4 flips = %.0f, flips_per_active = %.6f, flips_after_zero_frac = %.6f\n", t.flips,
                     frac(t.flips, t.active), frac(t.flips_after_zero, t.flips));
        std::fprintf(f, "G4 run_hist (1, 2, 3-5, 6-10, 11-30, >30) = %.0f %.0f %.0f %.0f %.0f %.0f   (всего %.0f)\n",
                     t.run_hist[0], t.run_hist[1], t.run_hist[2], t.run_hist[3], t.run_hist[4], t.run_hist[5], runs);
        std::fprintf(f, "G4 N_cv = %.6f   (коэф. вариации суммарного N по шагам с силой)\n", mean > 0.0 ? std::sqrt(var) / mean : 0.0);
        std::fprintf(f, "G4 keys = %.0f, contact_steps_per_key = %.3f   (экспозиция)\n", nk, frac(key_steps, nk));
        std::fprintf(f, "G5 max_pen_d = %.6f\n\n", t.max_pen_d);
    }

    void write_summary() const
    {
        FILE* f = std::fopen(summary_path_.c_str(), "w");
        if (!f) return;
        std::fprintf(f, "steps = %lld, hold_n = %d, fs_type = %d, ff_type = %d\n\n", steps_, hold_n_ ? 1 : 0,
                     cfg_.fiber_solid_contact_type, cfg_.fiber_fiber_contact_type);
        write_side(f, "fiber-fiber", ff_tot_, ff_keys_);
        write_side(f, "fiber-solid", fs_tot_, fs_keys_);
        std::fclose(f);
    }

    const SimulationConfig& cfg_;
    std::size_t every_ = 0;
    std::string summary_path_;
    bool hold_n_ = false;
    FILE* csv_ = nullptr;
    long long step_ = 0, steps_ = 0;
    std::vector<double> bead_dt_ref_;
    Side ff_, fs_;
    Totals ff_tot_, fs_tot_;
    std::unordered_map<Key, KeyState, KeyHash> ff_keys_, fs_keys_;
};

} // namespace contact_stats
} // namespace interactions
