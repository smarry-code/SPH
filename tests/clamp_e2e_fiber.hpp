#pragma once
// tests/clamp_e2e_fiber.hpp
// e2e-тесты закрепления бусин Fiber (Clamp/CLAMP_TESTS.md §3, Clamp/
// FIBER_CLAMP_SPEC.md §8): EF1 (рывок нити от зажима) и EF1-b (разрыв у кромки
// зажима), EF2 (натянутая струна), EF4 (трение о полностью зажатую нить),
// EF5 (частичная маска yz), EF6 (разбор CLAMP_BOX). EF3 -- в Solid-части.
// Числа кейсов -- tests/clamp_case.hpp, namespace clamp_case::fib.
//
// Знаковая конвенция реакции (clamp_report.fiber_force): R = -sum m*a_raw по
// зажатым компонентам, т.е. СИЛА ОПОРЫ на зажатые бусины. Баланс импульса
// монитора: P^{n+1} - P^n = h*R (опора -- единственная внешняя сила сцены),
// поэтому в EF1 R_x = -sqrt(km)*v0 (опора тянет нить назад), а ∫R dt равен
// приращению импульса нити.
//
// Сцена без Solid: фиктивная Solid-частица, ntotal = 0 (цикл пишет log.txt по
// sys.vx[moni_particle]); все прогоны -- одним непрерывным time_integration
// (RunGuard::set_step_observer), дополнительные величины снимает наблюдатель
// каждого кейса поверх ClampMonitor::on_step.
#include "clamp_case.hpp"
#include "clamp_e2e_common.hpp"

#include "fiber/fiber_input.hpp"

#include <cstring>

namespace clamp_e2e
{
namespace fibtest
{

namespace fc = clamp_case::fib;

// ---------------------------------------------------------------------------
// Материал и сцена (по образцу tests/clamp_integration_fiber.hpp: make_mat,
// add_chain, init_model -- заголовок не включается, он зависит от анонимного
// namespace test_clamp.cpp)
// ---------------------------------------------------------------------------
inline fiber::FiberMaterial make_mat(double mu_ff = 0.0, double eps_break = fc::eps_break)
{
    fiber::FiberMaterial m;
    m.E = fc::E;
    m.nu = fc::nu;
    m.d = fc::d;
    m.rho = fc::rho;
    m.eps_break = eps_break;
    m.mu_ff = mu_ff;
    return m;
}

inline double mat_area(const fiber::FiberMaterial& m) { return interactions::fiber_cross_section_area(m); }
inline double wave_speed(const fiber::FiberMaterial& m) { return std::sqrt(m.E / m.rho); }

inline void init_model(fiber::Model& model)
{
    model.solid.resize(1); // фиктивная частица, ntotal = 0
    model.solid.inzializ_index(0);
    model.fiber.resize(0);
    model.topology.resize(0);
}

// Прямая нить из n бусин: x_k = origin + k*spacing*dir, связи (k,k+1) с l0,
// l_break = (1+eps_break)*l0; масса m = rho*A*l0 (у концов половина), как в
// src/fiber_input.cpp. Возвращает индекс первой бусины.
inline std::size_t add_chain(fiber::Model& model, std::size_t n, const double (&origin)[3], const double (&dir)[3],
                             double spacing, double l0, int material_id, int object_id)
{
    const std::size_t first = model.fiber.size();
    model.fiber.resize(first + n);
    model.topology.prev.resize(first + n, -1);
    model.topology.next.resize(first + n, -1);
    const fiber::FiberMaterial& mat = model.materials[static_cast<std::size_t>(material_id)];
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t g = first + k;
        for (std::size_t dd = 0; dd < 3; ++dd)
        {
            model.fiber.x[dd][g] = origin[dd] + static_cast<double>(k) * spacing * dir[dd];
            model.fiber.vx[dd][g] = 0.0;
        }
        const double seg = (k == 0 || k + 1 == n) ? 0.5 : 1.0;
        model.fiber.mass[g] = mat.rho * mat_area(mat) * l0 * seg;
        model.fiber.object_id[g] = object_id;
        model.fiber.material_id[g] = material_id;
        if (k > 0) model.topology.prev[g] = static_cast<long long>(g - 1);
        if (k + 1 < n) model.topology.next[g] = static_cast<long long>(g + 1);
    }
    for (std::size_t k = 0; k + 1 < n; ++k)
    {
        fiber::FiberBond b;
        b.i = first + k;
        b.j = first + k + 1;
        b.l0 = l0;
        b.l_break = (1.0 + mat.eps_break) * l0;
        model.topology.bonds.push_back(b);
    }
    return first;
}

// «Пейсмейкер шага» для явного Эйлера. Явный Эйлер неустойчив на связях нитей
// при dt_limit (|λ| ~ 1.5 за шаг на старших осевых модах, цепочка разлетается
// за ~20 шагов; дефект схемы, не механизма закрепления). Шаг цикла --
// min(dt_limit по ВСЕМ связям), поэтому сцену Эйлера дополняет изолированная
// двухбусинная нить с коротким l0/factor: её dt_limit в factor раз меньше, и
// весь расчёт идёт с шагом dt_limit/factor. Нить лежит в 4d от оси сцены (контакта нет, а поисковая сетка контактов остаётся компактной),
// покоится и недеформирована (сила ровно 0.0, импульс и энергия не вносит).
// Индексы бусин сцены остаются прежними: пейсмейкер добавляется ПОСЛЕДНИМ.
inline void add_pacemaker(fiber::Model& model, double factor, int material_id, int object_id)
{
    if (!(factor > 1.0)) return;
    const double o[3] = {0.0, 0.0, 4.0 * fc::d}, dr[3] = {1, 0, 0}; // 4d над осью сцены: зазор > d, контакта нет
    const double la = fc::l0 / factor;
    add_chain(model, 2, o, dr, la, la, material_id, object_id);
}

// Маска бусин [first, first+count): размер fiber_mask всегда == nfiber.
inline void set_mask(fiber::Model& model, std::size_t first, std::size_t count,
                     std::uint8_t mk = interactions::kClampAll)
{
    if (model.clamps.fiber_mask.size() != model.fiber.size()) model.clamps.fiber_mask.resize(model.fiber.size(), 0);
    for (std::size_t k = 0; k < count; ++k) model.clamps.fiber_mask[first + k] = mk;
}

// Сила связей на каждую бусину по заданным положениям x (тот же закон, что
// interactions/fiber_bond.hpp: F = EA*max(l/l0 - 1, 0) вдоль связи, только
// растяжение). Независимая от механизма закрепления проверка реакции.
inline void bond_forces(const fiber::Model& m, const std::vector<double> (&x)[3], std::vector<double> (&f)[3])
{
    for (std::size_t dd = 0; dd < 3; ++dd) f[dd].assign(m.fiber.size(), 0.0);
    for (const auto& b : m.topology.bonds)
    {
        const fiber::FiberMaterial& mat = m.materials[static_cast<std::size_t>(m.fiber.material_id[b.i])];
        double r[3], l2 = 0.0;
        for (std::size_t dd = 0; dd < 3; ++dd)
        {
            r[dd] = x[dd][b.j] - x[dd][b.i];
            l2 += r[dd] * r[dd];
        }
        const double l = std::sqrt(l2);
        const double strain = l / b.l0 - 1.0;
        if (!(strain > 0.0) || !(l > 0.0)) continue;
        const double F = mat.E * mat_area(mat) * strain;
        for (std::size_t dd = 0; dd < 3; ++dd)
        {
            f[dd][b.i] += F * r[dd] / l;
            f[dd][b.j] -= F * r[dd] / l;
        }
    }
}

inline double max_bond_strain(const fiber::Model& m)
{
    double s = -1.0;
    for (const auto& b : m.topology.bonds)
    {
        double l2 = 0.0;
        for (std::size_t dd = 0; dd < 3; ++dd)
        {
            const double r = m.fiber.x[dd][b.j] - m.fiber.x[dd][b.i];
            l2 += r * r;
        }
        s = std::max(s, std::sqrt(l2) / b.l0 - 1.0);
    }
    return s;
}

inline double mass_sum(const fiber::Model& m, std::size_t first, std::size_t last_excl)
{
    double s = 0.0;
    for (std::size_t k = first; k < last_excl; ++k) s += m.fiber.mass[k];
    return s;
}

// ---------------------------------------------------------------------------
// Обработка ряда записей монитора
// ---------------------------------------------------------------------------
// Реакция шага n действует на интервале (t_n - dt_n, t_n]; вес h -- тот же
// множитель, что в балансе импульса (первый шаг PC -- полушаг).
inline double step_mid(const StepRecord& r) { return r.t - 0.5 * r.dt; }

// ∫ R_d dt по записям с серединой шага в [a, b) (сумма R*h -- точный интеграл
// шагового баланса).
inline double impulse_window(const std::vector<StepRecord>& rec, std::size_t d, double a, double b)
{
    double s = 0.0;
    for (const auto& r : rec)
        if (step_mid(r) >= a && step_mid(r) < b) s += r.Rf[d] * r.h;
    return s;
}

// Среднее R_d по времени на [a, b): ∫R dt / сумма h.
inline double mean_window(const std::vector<StepRecord>& rec, std::size_t d, double a, double b, std::size_t* nsteps = nullptr)
{
    double s = 0.0, w = 0.0;
    std::size_t n = 0;
    for (const auto& r : rec)
        if (step_mid(r) >= a && step_mid(r) < b)
        {
            s += r.Rf[d] * r.h;
            w += r.h;
            ++n;
        }
    if (nsteps) *nsteps = n;
    return w > 0.0 ? s / w : 0.0;
}

// ---------------------------------------------------------------------------
// Прогон одной сцены: монитор + наблюдатель + боевой цикл
// ---------------------------------------------------------------------------
inline RunOutcome run_scene(const Options& opt, fiber::Model& m, double t_end, const std::string& tag,
                            ClampMonitor& mon, std::function<void(double, double)> extra)
{
    double dt = 1.0e-9, time = 0.0; // dt перезаписывается циклом; PC на первом шаге -- полушаг
    auto obs = [&mon, extra](double d, double t)
    {
        mon.on_step(d, t);
        if (extra) extra(d, t);
    };
    return run_case(load_config(opt.solver, 1), m, 0, dt, time, t_end, 0, opt, tag, obs);
}

// ===========================================================================
// EF1. Рывок нити от зажима (+ EF1-b: разрыв у кромки зажима)
// ===========================================================================
struct Ef1Scene
{
    fiber::Model model;
    std::size_t nseg = 0, nclamp = 0;
    double v0 = 0.0;
    double c = 0.0, k = 0.0, m = 0.0, A = 0.0;
    double L = 0.0;       // свободная длина (последняя зажатая бусина -> свободный конец)
    double Mfree = 0.0;   // масса свободных бусин
    double plateau = 0.0; // -sqrt(km)*v0
};

// Нить из nseg сегментов (nseg+1 бусин) вдоль x, зажаты бусины 0..nclamp-1
// (xyz), остальным v_x = v0. code_scene == false -- только параметры (для EF6).
inline void ef1_fill_params(Ef1Scene& s, std::size_t nseg, double v0_factor, double eps_break = fc::eps_break)
{
    const fiber::FiberMaterial mat = make_mat(0.0, eps_break);
    s.nseg = nseg;
    s.nclamp = clamp_case::fib::ef1_nclamp;
    s.c = wave_speed(mat);
    s.A = mat_area(mat);
    s.k = mat.E * s.A / fc::l0;
    s.m = mat.rho * s.A * fc::l0;
    s.v0 = v0_factor * fc::eps_break * s.c; // v0 -- от штатного eps_break, даже если у Euler он поднят
    s.L = static_cast<double>(nseg - (s.nclamp - 1)) * fc::l0;
    s.plateau = -std::sqrt(s.k * s.m) * s.v0;
}

inline void ef1_set_velocities(Ef1Scene& s)
{
    for (std::size_t k = s.nclamp; k <= s.nseg; ++k) s.model.fiber.vx[0][k] = s.v0;
    s.Mfree = mass_sum(s.model, s.nclamp, s.nseg + 1);
}

inline Ef1Scene build_ef1(std::size_t nseg, double v0_factor, double eps_break = fc::eps_break, double pacemaker = 1.0)
{
    Ef1Scene s;
    ef1_fill_params(s, nseg, v0_factor, eps_break);
    init_model(s.model);
    s.model.materials = {make_mat(0.0, eps_break)};
    const double o[3] = {0, 0, 0}, dr[3] = {1, 0, 0};
    add_chain(s.model, nseg + 1, o, dr, fc::l0, fc::l0, 0, 0);
    set_mask(s.model, 0, s.nclamp);
    ef1_set_velocities(s);
    add_pacemaker(s.model, pacemaker, 0, 1);
    s.model.clamps.fiber_mask.resize(s.model.fiber.size(), 0);
    return s;
}

struct Ef1Extra
{
    std::vector<double> maxstrain; // max деформация связей на шаге
    std::vector<std::size_t> nbonds;
};

inline void run_EF1_impl(const Options& opt, bool variant_b)
{
    const bool pc = opt.solver == 0;
    const std::size_t nseg = opt.full ? clamp_case::fib::ef1_nseg_full
                                       : (pc ? clamp_case::fib::ef1_nseg_smoke : clamp_case::fib::ef1_nseg_euler);
    const double v0_factor = variant_b ? clamp_case::fib::ef1b_v0_factor : clamp_case::fib::ef1_v0_factor;
    const double pace = pc ? 1.0 : clamp_case::fib::euler_pacemaker;
    Ef1Scene s = build_ef1(nseg, v0_factor, fc::eps_break, pace);
    const double t_end = clamp_case::fib::ef1_tend_factor * s.L / s.c;
    std::printf("EF1%s: нить %zu сегментов, зажато %zu бусин, L = %.4g м, c = %.4g м/с, v0 = %.4g м/с, "
                "-sqrt(km) v0 = %s Н, 2L/c = %s с, t_end = %s с\n",
                variant_b ? "-b" : "", nseg, s.nclamp, s.L, s.c, s.v0, fmt(s.plateau).c_str(),
                fmt(2 * s.L / s.c).c_str(), fmt(t_end).c_str());

    const std::vector<std::uint8_t> mask0 = s.model.clamps.fiber_mask;
    const std::size_t nb0 = s.model.topology.bonds.size();
    const std::vector<fiber::FiberBond> bonds0 = s.model.topology.bonds;
    const int id_free0 = s.model.fiber.object_id[s.nclamp];
    ClampMonitor mon(s.model, 0, opt.solver);
    Ef1Extra ex;
    std::size_t prev_nb = nb0;
    std::vector<std::size_t> break_steps; // номера шагов (с 0), на которых число связей уменьшилось
    std::vector<std::pair<std::size_t, std::size_t>> lost_first;
    auto extra = [&](double, double)
    {
        ex.maxstrain.push_back(max_bond_strain(s.model));
        const std::size_t nb = s.model.topology.bonds.size();
        ex.nbonds.push_back(nb);
        if (nb != prev_nb)
        {
            break_steps.push_back(ex.nbonds.size() - 1);
            if (lost_first.empty())
                for (const auto& b0 : bonds0)
                {
                    bool found = false;
                    for (const auto& b : s.model.topology.bonds)
                        if (b.i == b0.i && b.j == b0.j) { found = true; break; }
                    if (!found) lost_first.emplace_back(b0.i, b0.j);
                }
            prev_nb = nb;
        }
    };
    const RunOutcome run = run_scene(opt, s.model, t_end, variant_b ? "ef1b" : "ef1", mon, extra);
    maybe_write_monitor(opt, mon, variant_b ? ".ef1b" : ".ef1");
    const auto& rec = mon.records();

    // §1.3 п. 3: только нити -- 1e-12 * max|P| (округление накопленных сумм
    // обсуждается у clamp_case::fib::balance_rel).
    const double bal_tol = clamp_case::fib::balance_rel * mon.max_abs_momentum();
    common_checks(variant_b ? " EF1-b" : " EF1", s.model, 0, mon, run, bal_tol);
    if (run.threw || rec.empty()) return;

    const double Lc = s.L / s.c;
    std::printf("\n[EF1%s] Критерии теста\n", variant_b ? "-b" : "");
    note("шагов " + std::to_string(rec.size()) + ", dt = " + fmt(run.dt_first, 3) + ", wall " + fmt(run.wall_s, 2) + " с");

    if (!variant_b)
    {
        // п. 6: нет разрывов
        double smax = -1.0;
        for (double v : ex.maxstrain) smax = std::max(smax, v);
        check(s.model.topology.bonds.size() == nb0 && break_steps.empty(),
              "нет разрывов: связей " + std::to_string(s.model.topology.bonds.size()) + " из " + std::to_string(nb0) +
                  ", max деформация " + fmt(smax, 3) + " < eps_break = " + fmt(fc::eps_break, 3));
        check(s.model.clamps.fiber_mask == mask0, "маска не изменилась");

        const double tw_a = 0.2 * Lc, tw_b = 1.8 * Lc;
        // п. 2: плато
        std::size_t nw = 0;
        const double Rm = mean_window(rec, 0, tw_a, tw_b, &nw);
        double ring = 0.0;
        for (const auto& r : rec)
            if (step_mid(r) >= tw_a && step_mid(r) < tw_b) ring = std::max(ring, std::fabs(r.Rf[0] - s.plateau));
        const double perr = (Rm - s.plateau) / std::fabs(s.plateau);
        const bool plateau_ok = std::fabs(perr) <= clamp_case::fib::ef1_plateau_tol;
        if (pc || opt.full)
            check(plateau_ok, "плато R_x на [0.2, 1.8] L/c: " + fmt(Rm) + " против -sqrt(km) v0 = " + fmt(s.plateau) +
                                  ", отклонение " + fmt(100 * perr, 3) + " % (допуск " +
                                  fmt(100 * clamp_case::fib::ef1_plateau_tol, 2) + " %, " + std::to_string(nw) + " шагов)");
        else
            check(plateau_ok, "плато R_x: " + fmt(Rm) + " против " + fmt(s.plateau) + ", отклонение " +
                                  fmt(100 * perr, 3) + " %");
        note("звон: max |R - плато| в окне = " + fmt(100 * ring / std::fabs(s.plateau), 3) + " % плато");

        // п. 3: длительность
        double t_drop = -1.0;
        for (const auto& r : rec)
            if (step_mid(r) >= 1.5 * Lc && std::fabs(r.Rf[0]) < 0.1 * std::fabs(s.plateau))
            {
                t_drop = step_mid(r);
                break;
            }
        const double t2 = 2.0 * Lc;
        check(t_drop > 0.0 && std::fabs(t_drop - t2) <= 0.03 * t2,
              "падение R_x ниже 10 % плато в t = " + fmt(t_drop) + " с против 2L/c = " + fmt(t2) + " с (" +
                  fmt(100 * (t_drop - t2) / t2, 3) + " %, допуск 3 %)");
        double tail = 0.0;
        std::size_t ntail = 0;
        for (const auto& r : rec)
            if (step_mid(r) >= 2.2 * Lc)
            {
                tail = std::max(tail, std::fabs(r.Rf[0]));
                ++ntail;
            }
        check(ntail > 0 && tail <= 1.0e-3 * std::fabs(s.plateau),
              "после 2.2 L/c: max |R_x| = " + fmt(tail) + " <= 1e-3 плато (" + fmt(1e-3 * std::fabs(s.plateau)) +
                  ", " + std::to_string(ntail) + " шагов)");

        // п. 4: импульс
        double J = 0.0;
        for (const auto& r : rec) J += r.Rf[0] * r.h;
        const double Jexp = -2.0 * s.Mfree * s.v0;
        check(std::fabs(J - Jexp) <= 0.02 * std::fabs(Jexp),
              "импульс ∫R_x dt = " + fmt(J) + " против -2 M v0 = " + fmt(Jexp) + " (" +
                  fmt(100 * (J - Jexp) / std::fabs(Jexp), 3) + " %, допуск 2 %)");

        // п. 5: энергия
        const double E0 = mon.initial_energy();
        double dEmax = 0.0, Emax = 0.0, Emin = 1e300;
        for (const auto& r : rec)
        {
            const double E = r.Ekin + r.Eint + r.Ebond;
            dEmax = std::max(dEmax, std::fabs(E - E0));
            Emax = std::max(Emax, E);
            Emin = std::min(Emin, E);
        }
        const std::string en = "E_kin + E_bond: max |E - E0|/E0 = " + fmt(100 * dEmax / E0, 3) + " % (E0 = " + fmt(E0) + ")";
        if (pc)
            check(dEmax <= clamp_case::fib::ef1_energy_tol * E0, en + ", допуск " + fmt(100 * clamp_case::fib::ef1_energy_tol, 2) + " %");
        else
            note("Euler: " + en);
    }
    else
    {
        // EF1-b
        check(!break_steps.empty(), "связь порвалась (шаг " + (break_steps.empty() ? std::string("-") : std::to_string(break_steps[0] + 1)) + ")");
        if (!break_steps.empty())
        {
            const std::size_t kb = break_steps[0];
            const double t_break = rec[kb].t;
            check(t_break <= clamp_case::fib::ef1b_break_time_frac * Lc,
                  "разрыв в первые шаги: шаг " + std::to_string(kb + 1) + ", t = " + fmt(t_break) + " с <= " +
                      fmt(clamp_case::fib::ef1b_break_time_frac, 2) + " L/c = " + fmt(clamp_case::fib::ef1b_break_time_frac * Lc) + " с");
            check(lost_first.size() == 1 && lost_first[0].first == s.nclamp - 1 && lost_first[0].second == s.nclamp,
                  "первой рвётся связь (" + std::to_string(s.nclamp - 1) + "," + std::to_string(s.nclamp) +
                      ") «последняя зажатая - первая свободная»; на этом шаге пропали: " +
                      (lost_first.empty() ? std::string("-") : std::to_string(lost_first[0].first) + "-" + std::to_string(lost_first[0].second)) +
                      ", всего " + std::to_string(lost_first.size()));
            // R == 0.0 после разрыва до конца прогона
            bool zero_after = true;
            double rmax_after = 0.0;
            for (std::size_t n = kb + 1; n < rec.size(); ++n)
                for (std::size_t dd = 0; dd < 3; ++dd)
                {
                    if (rec[n].Rf[dd] != 0.0) zero_after = false;
                    rmax_after = std::max(rmax_after, std::fabs(rec[n].Rf[dd]));
                }
            check(zero_after && rec.size() > kb + 1, "R == 0.0 на всех " + std::to_string(rec.size() - kb - 1) +
                                                       " шагах после разрыва (max |R| = " + fmt(rmax_after) + ")");
            bool before_nonzero = false;
            for (std::size_t n = 0; n <= kb; ++n)
                if (rec[n].Rf[0] != 0.0) before_nonzero = true;
            check(before_nonzero, "до разрыва реакция ненулевая (проверка не вырождена)");
            check(s.model.topology.bonds.size() == nb0 - 1, "число связей уменьшилось ровно на одну (" +
                                                              std::to_string(s.model.topology.bonds.size()) + " из " + std::to_string(nb0) + ")");
            check(s.model.clamps.fiber_mask == mask0, "маска не изменилась после разрыва");
            check(s.model.fiber.object_id[s.nclamp] != id_free0 &&
                      s.model.fiber.object_id[s.nclamp] != s.model.fiber.object_id[s.nclamp - 1],
                  "у оторвавшейся части новый object_id; зажатая часть сохранила прежний");
            // импульс оторвавшейся нити = M v0 + ∫R dt до разрыва (точный баланс)
            double Pfree = 0.0, Pclamped = 0.0;
            for (std::size_t k = 0; k < s.model.fiber.size(); ++k)
                (k >= s.nclamp ? Pfree : Pclamped) += s.model.fiber.mass[k] * s.model.fiber.vx[0][k];
            double Jb = 0.0;
            for (std::size_t n = 0; n <= kb; ++n) Jb += rec[n].Rf[0] * rec[n].h;
            const double Pexp = s.Mfree * s.v0 + Jb;
            check(std::fabs(Pfree - Pexp) <= 1e-12 * s.Mfree * s.v0,
                  "импульс оторвавшейся нити " + fmt(Pfree) + " == M_free v0 + ∫R dt до разрыва = " + fmt(Pexp) +
                      " (невязка " + fmt(std::fabs(Pfree - Pexp)) + ", допуск " + fmt(1e-12 * s.Mfree * s.v0) + ")");
            check(Pclamped == 0.0, "импульс зажатой части ровно 0.0");
            note("∫R dt до разрыва = " + fmt(Jb) + ", M_free v0 = " + fmt(s.Mfree * s.v0) + " (опора успела отнять " +
                 fmt(100 * std::fabs(Jb) / (s.Mfree * s.v0), 3) + " % импульса)");
        }
    }
}

inline void run_EF1(const Options& opt)
{
    if (opt.variant == "b") run_EF1_impl(opt, true);
    else if (opt.variant.empty()) run_EF1_impl(opt, false);
    else { std::printf("  FAIL : EF1: неизвестный вариант --variant=%s (есть: b)\n", opt.variant.c_str()); ++e2e::failures; }
}

// ===========================================================================
// EF2. Натянутая струна с зажатыми концами
// ===========================================================================
struct Ef2Scene
{
    fiber::Model model;
    std::size_t n = 0;              // сегментов между последними зажатыми бусинами
    std::size_t nbeads_chain = 0;   // n + 3 (по 2 зажатые с каждого конца)
    std::size_t first_free = 2, mid = 0;
    double eps0 = 0.0, l = 0.0, N = 0.0, m = 0.0, A = 0.0;
    double omega1 = 0.0, T = 0.0, amp = 0.0, V = 0.0;
};

inline Ef2Scene build_ef2(std::size_t n, double pacemaker)
{
    Ef2Scene s;
    const fiber::FiberMaterial mat = make_mat();
    s.n = n;
    s.nbeads_chain = n + 3;
    s.mid = 1 + n / 2;
    s.eps0 = clamp_case::fib::ef2_eps0_factor * mat.eps_break;
    s.l = fc::l0 * (1.0 + s.eps0);
    s.A = mat_area(mat);
    s.N = mat.E * s.A * s.eps0;
    s.m = mat.rho * s.A * fc::l0;
    s.omega1 = 2.0 * std::sqrt(s.N / (s.m * s.l)) * std::sin(pi / (2.0 * static_cast<double>(n)));
    s.T = 2.0 * pi / s.omega1;
    s.amp = clamp_case::fib::ef2_amp_over_l * s.l;
    s.V = s.amp * s.omega1;
    init_model(s.model);
    s.model.materials = {mat};
    const double o[3] = {0, 0, 0}, dr[3] = {1, 0, 0};
    add_chain(s.model, s.nbeads_chain, o, dr, s.l, fc::l0, 0, 0);
    set_mask(s.model, 0, 2);
    set_mask(s.model, n + 1, 2);
    for (std::size_t k = 2; k <= n; ++k)
        s.model.fiber.vx[1][k] = s.V * std::sin(pi * static_cast<double>(k - 1) / static_cast<double>(n));
    add_pacemaker(s.model, pacemaker, 0, 1);
    s.model.clamps.fiber_mask.resize(s.model.fiber.size(), 0);
    return s;
}

// Нули ряда (линейная интерполяция по соседним записям; начальный нуль не входит).
inline std::vector<double> zero_crossings(const std::vector<double>& t, const std::vector<double>& y)
{
    std::vector<double> z;
    for (std::size_t i = 1; i < y.size(); ++i)
        if ((y[i - 1] < 0.0) != (y[i] < 0.0) && y[i - 1] != y[i])
            z.push_back(t[i - 1] + (t[i] - t[i - 1]) * (0.0 - y[i - 1]) / (y[i] - y[i - 1]));
    return z;
}

inline void run_EF2(const Options& opt)
{
    const bool pc = opt.solver == 0;
    const std::size_t n = opt.full ? clamp_case::fib::ef2_n_full : clamp_case::fib::ef2_n_smoke;
    Ef2Scene s = build_ef2(n, pc ? 1.0 : clamp_case::fib::euler_pacemaker);
    const double periods = opt.full ? clamp_case::fib::ef2_periods_full : clamp_case::fib::ef2_periods_smoke;
    const double t_end = (periods + clamp_case::fib::ef2_tail) * s.T;
    std::printf("EF2: струна %zu сегментов, натяжение eps0 = %.4g (N = %.5g Н), l = %.6g м, "
                "T_дискр = %.6g с, амплитуда %.4g м, t_end = %.5g с (%.3g T)\n",
                n, s.eps0, s.N, s.l, s.T, s.amp, t_end, t_end / s.T);

    const std::vector<std::uint8_t> mask0 = s.model.clamps.fiber_mask;
    ClampMonitor mon(s.model, 0, opt.solver);
    std::vector<double> ts, ymid, Rxl, Rxr, Rexp_dev;
    std::vector<double> xprev[3];
    for (std::size_t dd = 0; dd < 3; ++dd) xprev[dd] = s.model.fiber.x[dd];
    double sym_max = 0.0, y_abs_max = 0.0, Rdev_max = 0.0, strain_max = -1.0;
    double Rl_sum = 0.0, Rr_sum = 0.0, w_sum = 0.0;
    std::vector<double> Ry_series;
    auto extra = [&](double, double t)
    {
        const auto& r = mon.records().back();
        // независимая реакция: -сила связей на зажатые бусины по положениям НАЧАЛА шага
        std::vector<double> f[3];
        bond_forces(s.model, xprev, f);
        double Rl[3] = {}, Rr[3] = {};
        for (std::size_t dd = 0; dd < 3; ++dd)
        {
            Rl[dd] = -(f[dd][0] + f[dd][1]);
            Rr[dd] = -(f[dd][s.n + 1] + f[dd][s.n + 2]);
            Rdev_max = std::max(Rdev_max, std::fabs(r.Rf[dd] - (Rl[dd] + Rr[dd])));
        }
        Rl_sum += Rl[0] * r.h;
        Rr_sum += Rr[0] * r.h;
        w_sum += r.h;
        ts.push_back(t);
        ymid.push_back(s.model.fiber.x[1][s.mid]);
        Ry_series.push_back(r.Rf[1]);
        y_abs_max = std::max(y_abs_max, std::fabs(ymid.back()));
        for (std::size_t j = 1; j < s.n / 2; ++j)
            sym_max = std::max(sym_max, std::fabs(s.model.fiber.x[1][1 + j] - s.model.fiber.x[1][1 + s.n - j]));
        strain_max = std::max(strain_max, max_bond_strain(s.model));
        for (std::size_t dd = 0; dd < 3; ++dd) xprev[dd] = s.model.fiber.x[dd];
    };
    const RunOutcome run = run_scene(opt, s.model, t_end, "ef2", mon, extra);
    maybe_write_monitor(opt, mon, ".ef2");
    const auto& rec = mon.records();
    const double bal_tol = clamp_case::fib::balance_rel * mon.max_abs_momentum();
    common_checks(" EF2", s.model, 0, mon, run, bal_tol);
    if (run.threw || rec.empty()) return;

    std::printf("\n[EF2] Критерии теста\n");
    note("шагов " + std::to_string(rec.size()) + ", dt = " + fmt(run.dt_first, 3) + ", wall " + fmt(run.wall_s, 2) + " с");
    check(s.model.clamps.fiber_mask == mask0, "маска не изменилась");
    check(strain_max < fc::eps_break, "нет разрывов, max деформация " + fmt(strain_max, 3) + " < eps_break = " + fmt(fc::eps_break, 3));

    // п. 4: реакция
    check(Rdev_max <= clamp_case::fib::ef2_react_dev * s.N,
          "реакция опоры из цикла (R_fiber) совпадает с независимой суммой -F_связей по зажатым бусинам: max |dR| = " +
              fmt(Rdev_max) + " Н <= " + fmt(clamp_case::fib::ef2_react_dev * s.N) + " (допуск " + fmt(clamp_case::fib::ef2_react_dev, 1) + " N)");
    const double Rl_mean = Rl_sum / w_sum, Rr_mean = Rr_sum / w_sum;
    check(std::fabs(Rl_mean + s.N) <= clamp_case::fib::ef2_react_tol * s.N &&
              std::fabs(Rr_mean - s.N) <= clamp_case::fib::ef2_react_tol * s.N,
          "среднее R_x: левый зажим " + fmt(Rl_mean) + ", правый " + fmt(Rr_mean) + " при N = " + fmt(s.N) + " (отклонения " +
              fmt(100 * (Rl_mean + s.N) / s.N, 2) + " %, " + fmt(100 * (Rr_mean - s.N) / s.N, 2) + " %; допуск " +
              fmt(100 * clamp_case::fib::ef2_react_tol, 2) + " %)");
    double Rx_tot_max = 0.0;
    for (const auto& r : rec) Rx_tot_max = std::max(Rx_tot_max, std::fabs(r.Rf[0]));
    note("clamp_report суммирует оба зажима: R_x(левый) + R_x(правый) = 0 до " + fmt(Rx_tot_max) + " Н; по знаку R = сила ОПОРЫ, "
         "левая опора тянет струну к -x (R_x = -N), правая к +x (R_x = +N); в постановке стоит обратный знак");

    // R_y: амплитуда 2 N A sin(pi/n)/l, частота w1
    double Ry_max = 0.0;
    for (double v : Ry_series) Ry_max = std::max(Ry_max, std::fabs(v));
    const double Ry_exp = 2.0 * s.N * s.amp * std::sin(pi / static_cast<double>(s.n)) / s.l;
    check(std::fabs(Ry_max - Ry_exp) <= clamp_case::fib::ef2_Ry_tol * Ry_exp,
          "амплитуда R_y = " + fmt(Ry_max) + " против 2 N A sin(pi/n)/l = " + fmt(Ry_exp) + " (" + fmt(100 * (Ry_max - Ry_exp) / Ry_exp, 3) +
              " %, допуск " + fmt(100 * clamp_case::fib::ef2_Ry_tol, 2) + " %)");

    // п. 5: симметрия. Явный Эйлер в полном прогоне (~1e6 шагов) усиливает округление
    // старших осевых мод и рано или поздно ломает симметрию -- там note.
    {
        const std::string ss = "симметрия y_j = y_{n-j}: max |dy| = " + fmt(sym_max) + " (" + fmt(sym_max / y_abs_max, 2) + " от max|y| = " +
                               fmt(y_abs_max) + ", допуск " + fmt(clamp_case::fib::ef2_sym_tol, 1) + ")";
        if (pc || !opt.full) check(sym_max <= clamp_case::fib::ef2_sym_tol * y_abs_max, ss);
        else note("Euler full: " + ss + "; рост старших мод Эйлера");
    }

    // п. 2, 3: период и энергия
    const double E0 = mon.initial_energy();
    double dEmax = 0.0, Ekin_max = 0.0, Ebond_min = 1e300, Ebond_max = 0.0;
    for (const auto& r : rec)
    {
        dEmax = std::max(dEmax, std::fabs(r.Ekin + r.Eint + r.Ebond - E0));
        Ekin_max = std::max(Ekin_max, r.Ekin);
        Ebond_min = std::min(Ebond_min, r.Ebond);
        Ebond_max = std::max(Ebond_max, r.Ebond);
    }
    const std::string en = "E_kin + E_bond: max |E - E0|/E0 = " + fmt(100 * dEmax / E0, 3) + " % (E0 = " + fmt(E0) +
                           "), относительно динамической части (max E_kin = " + fmt(Ekin_max) + "): " + fmt(100 * dEmax / Ekin_max, 3) + " %";
    if (pc)
        check(dEmax <= clamp_case::fib::ef2_energy_tol * E0, en + "; допуск " + fmt(100 * clamp_case::fib::ef2_energy_tol, 2) + " % E0");
    else
        note("Euler: " + en);

    const auto z = zero_crossings(ts, ymid);
    {
        if (z.size() >= 3)
        {
            const double Tm = 2.0 * (z.back() - z.front()) / static_cast<double>(z.size() - 1);
            const double err = (Tm - s.T) / s.T;
            const std::string ps = "период по нулям y средней бусины (" + std::to_string(z.size()) + " нулей): T = " + fmt(Tm, 6) +
                                   " с против 2π/ω1 = " + fmt(s.T, 6) + " с, отклонение " + fmt(100 * err, 3) + " %";
            if (pc)
                check(std::fabs(err) <= clamp_case::fib::ef2_period_tol, ps + " (допуск " + fmt(100 * clamp_case::fib::ef2_period_tol, 2) + " %)");
            else
                note("Euler: " + ps);
            // пределы струны: T = 2L sqrt(rho A / N), L = n*l (FIBER_CLAMP_SPEC §8 п. 3)
            const double Tcont = 2.0 * static_cast<double>(s.n) * s.l * std::sqrt(s.m / s.l / s.N);
            note("непрерывная струна T = 2L sqrt(mu/N) = " + fmt(Tcont, 6) + " с (отличие дискретной цепочки " + fmt(100 * (s.T - Tcont) / Tcont, 3) + " %)");
        }
        else
            warn("EF2: нулей y средней бусины меньше трёх (" + std::to_string(z.size()) + "): прогон короче периода, проверка периода пропущена");
    }
}

// ===========================================================================
// EF4. Трение о полностью зажатую нить
// ===========================================================================
// Нить B вдоль y зажата целиком, свободная нить A вдоль x над ней падает на B со
// скоростью v_n и скользит ВДОЛЬ B со скоростью v_t (A лежит над серединой
// сегмента A, B -- вдоль y). Скольжение вдоль зажатой нити выбрано потому, что
// нормаль контакта перпендикулярна оси B, т.е. не имеет y-компоненты, а сцена
// зеркально симметрична относительно x = 0: боковая сила в y -- ТОЛЬКО трение,
// нормальный импульс -- ТОЛЬКО по z (и x-составляющие нормалей гасятся
// симметрией). Тогда ∫R_z dt = ∫N dt (N >= 0), ΔP_Ay = -J_f.
struct Ef4Scene
{
    fiber::Model model;
    std::size_t nA = 0, nB = 0, firstB = 0;
    double v_n = 0.0, v_t = 0.0, mu = 0.0, MA = 0.0;
};

inline Ef4Scene build_ef4(double mu)
{
    Ef4Scene s;
    namespace cf = clamp_case::fib;
    s.nA = cf::ef4_nA;
    s.nB = cf::ef4_nB;
    s.v_n = cf::ef4_vn;
    s.v_t = cf::ef4_vt;
    s.mu = mu;
    init_model(s.model);
    s.model.materials = {make_mat(mu)};
    const double zA = fc::d * (1.0 + cf::ef4_gap_over_d);
    const double oA[3] = {-(static_cast<double>(s.nA) / 2.0 - 0.5) * fc::l0, cf::ef4_yA_over_l0 * fc::l0, zA}, dA[3] = {1, 0, 0};
    const std::size_t a0 = add_chain(s.model, s.nA, oA, dA, fc::l0, fc::l0, 0, 0);
    const double oB[3] = {0.0, -(static_cast<double>(s.nB - 1) / 2.0) * fc::l0, 0.0}, dB[3] = {0, 1, 0};
    s.firstB = add_chain(s.model, s.nB, oB, dB, fc::l0, fc::l0, 0, 1);
    s.model.clamps.fiber_mask.assign(s.model.fiber.size(), 0);
    set_mask(s.model, s.firstB, s.nB);
    for (std::size_t k = 0; k < s.nA; ++k)
    {
        s.model.fiber.vx[1][a0 + k] = s.v_t;
        s.model.fiber.vx[2][a0 + k] = -s.v_n;
    }
    s.MA = mass_sum(s.model, a0, a0 + s.nA);
    return s;
}

struct Ef4Result
{
    double P0[3] = {}, P1[3] = {};   // импульс нити A в начале и в конце
    double Jn = 0.0;                 // ∫ R_z dt
    double Rz_min = 0.0;
    double Rsum[3] = {};             // ∫ R dt по осям
    double t_first_contact = -1.0, t_last_contact = -1.0;
    std::size_t contact_steps = 0;
    bool ok = false;
};

inline Ef4Result run_ef4_case(const Options& opt, double mu, const std::string& tag)
{
    Ef4Result res;
    Ef4Scene s = build_ef4(mu);
    const double t_end = clamp_case::fib::ef4_t_end;
    ClampMonitor mon(s.model, 0, opt.solver);
    auto momA = [&](double (&P)[3])
    {
        P[0] = P[1] = P[2] = 0.0;
        for (std::size_t k = 0; k < s.nA; ++k)
            for (std::size_t dd = 0; dd < 3; ++dd) P[dd] += s.model.fiber.mass[k] * s.model.fiber.vx[dd][k];
    };
    momA(res.P0);
    std::vector<double> PAy, PAz, tt;
    auto extra = [&](double, double t)
    {
        double P[3];
        momA(P);
        PAy.push_back(P[1]);
        PAz.push_back(P[2]);
        tt.push_back(t);
    };
    const RunOutcome run = run_scene(opt, s.model, t_end, tag, mon, extra);
    maybe_write_monitor(opt, mon, "." + tag);
    const auto& rec = mon.records();
    common_checks(" EF4 (" + tag + ")", s.model, 0, mon, run, clamp_case::fib::balance_rel * mon.max_abs_momentum());
    if (run.threw || rec.empty()) return res;
    momA(res.P1);
    res.Rz_min = 1e300;
    for (const auto& r : rec)
    {
        for (std::size_t dd = 0; dd < 3; ++dd) res.Rsum[dd] += r.Rf[dd] * r.h;
        res.Jn += r.Rf[2] * r.h;
        res.Rz_min = std::min(res.Rz_min, r.Rf[2]);
        if (r.Rf[2] != 0.0)
        {
            ++res.contact_steps;
            if (res.t_first_contact < 0) res.t_first_contact = step_mid(r);
            res.t_last_contact = step_mid(r);
        }
    }
    res.ok = true;
    return res;
}

inline void run_EF4(const Options& opt)
{
    namespace cf = clamp_case::fib;
    const double mu = cf::ef4_mu;
    std::printf("EF4: нить A (%zu бусин) падает на полностью зажатую нить B (%zu бусин) со скоростью %.3g м/с и скользит вдоль B со %.3g м/с, mu_ff = %.3g\n",
                cf::ef4_nA, cf::ef4_nB, cf::ef4_vn, cf::ef4_vt, mu);
    const Ef4Result r0 = run_ef4_case(opt, 0.0, "ef4_mu0");
    const Ef4Result r1 = run_ef4_case(opt, mu, "ef4_mu");
    if (!r0.ok || !r1.ok) { check(false, "EF4: прогон не состоялся"); return; }
    // Критерий (4) постановки о паре «зажатая-зажатая»: внутри B все бусины одного object_id, а
    // контакт ищется только между разными object_id, так что таких пар в сцене нет и проверить
    // их здесь нельзя. Пара разных полностью зажатых нитей проверяется в F2 (tests/test_clamp.cpp).
    warn("EF4: пар трения «зажатая-зажатая» в сцене нет (одна нить B, контакт только между разными "
         "object_id) -- критерий (4) в этой части не проверяется, см. F2 в test_clamp");
    std::printf("\n[EF4] Критерии теста\n");
    note("mu = 0: контакт " + std::to_string(r0.contact_steps) + " шагов в t = [" + fmt(r0.t_first_contact) + ", " + fmt(r0.t_last_contact) +
         "], J_n = " + fmt(r0.Jn) + " Н*с, P_Az: " + fmt(r0.P0[2]) + " -> " + fmt(r0.P1[2]));
    note("mu > 0: контакт " + std::to_string(r1.contact_steps) + " шагов в t = [" + fmt(r1.t_first_contact) + ", " + fmt(r1.t_last_contact) +
         "], J_n = " + fmt(r1.Jn) + " Н*с, P_Ay: " + fmt(r1.P0[1]) + " -> " + fmt(r1.P1[1]));
    // п. 2
    const double dPy0 = std::fabs(r0.P1[1] - r0.P0[1]) / std::fabs(r0.P0[1]);
    check(r0.contact_steps > 0 && r0.Jn > 0.0, "mu = 0: контакт был (J_n = " + fmt(r0.Jn) + " Н*с > 0), проверка не вырождена");
    check(dPy0 <= cf::ef4_free_tol, "mu = 0: касательный импульс A не меняется: |dP_y|/P_y0 = " + fmt(dPy0) + " <= " + fmt(cf::ef4_free_tol, 1));
    // п. 3
    const double Jf = r1.P0[1] - r1.P1[1];
    check(r1.contact_steps > 0 && r1.Jn > 0.0, "mu > 0: контакт был (J_n = " + fmt(r1.Jn) + " Н*с > 0)");
    check(Jf > 0.0, "mu > 0: трение отнимает касательный импульс: J_f = " + fmt(Jf) + " Н*с > 0");
    check(Jf <= mu * r1.Jn * (1.0 + cf::ef4_coulomb_slack), "закон Кулона: J_f = " + fmt(Jf) + " <= mu*∫N dt = " + fmt(mu * r1.Jn) + " (запас " + fmt(cf::ef4_coulomb_slack, 1) + ")");
    check(Jf >= cf::ef4_sliding_frac * mu * r1.Jn, "скольжение, а не пропуск трения: J_f/(mu J_n) = " + fmt(Jf / (mu * r1.Jn), 3) + " >= " + fmt(cf::ef4_sliding_frac, 2));
    // п. 5
    check(r1.Rz_min >= 0.0 && r0.Rz_min >= 0.0, "R_z >= 0 на каждом шаге (нормальный контакт только отталкивает): min = " + fmt(r1.Rz_min));
}

// ===========================================================================
// EF5. Частичная маска нити: концы на «рельсах» вдоль оси
// ===========================================================================
inline void run_EF5(const Options& opt)
{
    namespace cf = clamp_case::fib;
    const bool pc = opt.solver == 0;
    const std::size_t nseg = opt.full ? cf::ef1_nseg_full : cf::ef1_nseg_smoke;
    Ef1Scene s = build_ef1(nseg, cf::ef1_v0_factor, fc::eps_break, pc ? 1.0 : cf::euler_pacemaker);
    // Сцена EF1 без маски xyz: рельсы -- маска yz на двух концевых бусинах. Всем бусинам нити v_x = v0.
    const std::uint8_t kYZ = 0x6;
    s.model.clamps.fiber_mask.assign(s.model.fiber.size(), 0);
    s.model.clamps.fiber_mask[0] = kYZ;
    s.model.clamps.fiber_mask[nseg] = kYZ;
    for (std::size_t k = 0; k <= nseg; ++k) s.model.fiber.vx[0][k] = s.v0;
    const double v1 = cf::ef5_v1_over_v0 * s.v0;
    s.model.fiber.vx[1][nseg / 2] = v1;       // поперечный толчок средней бусине (постановка)
    s.model.fiber.vx[1][1] = v1;              // и бусинам, соседним с концами: иначе рельсы не нагружены
    s.model.fiber.vx[1][nseg - 1] = v1;       // (R_y, R_z были бы тождественно нулём)
    const double t_end = cf::ef1_tend_factor * s.L / s.c;
    std::printf("EF5: нить %zu сегментов, концы (бусины 0 и %zu) на рельсах (маска yz), v_x = %.4g м/с у всех, v_y = %.3g м/с у средней и соседних с концами; t_end = %.4g с\n",
                nseg, nseg, s.v0, v1, t_end);

    const std::vector<double> x0[3] = {s.model.fiber.x[0], s.model.fiber.x[1], s.model.fiber.x[2]};
    double Mtot = 0.0, xcm0 = 0.0;
    for (std::size_t k = 0; k <= nseg; ++k)
    {
        Mtot += s.model.fiber.mass[k];
        xcm0 += s.model.fiber.mass[k] * s.model.fiber.x[0][k];
    }
    xcm0 /= Mtot;
    ClampMonitor mon(s.model, 0, opt.solver);
    double dx_end_max = 0.0, dyz_end_max = 0.0;
    auto extra = [&](double, double t)
    {
        for (std::size_t e : {std::size_t(0), nseg})
        {
            dx_end_max = std::max(dx_end_max, std::fabs(s.model.fiber.x[0][e] - (x0[0][e] + s.v0 * t)));
            dyz_end_max = std::max(dyz_end_max, std::max(std::fabs(s.model.fiber.x[1][e] - x0[1][e]), std::fabs(s.model.fiber.x[2][e] - x0[2][e])));
        }
    };
    const RunOutcome run = run_scene(opt, s.model, t_end, "ef5", mon, extra);
    maybe_write_monitor(opt, mon, ".ef5");
    const auto& rec = mon.records();
    common_checks(" EF5", s.model, 0, mon, run, cf::balance_rel * mon.max_abs_momentum());
    if (run.threw || rec.empty()) return;

    std::printf("\n[EF5] Критерии теста\n");
    note("шагов " + std::to_string(rec.size()) + ", dt = " + fmt(run.dt_first, 3) + ", wall " + fmt(run.wall_s, 2) + " с, закреплённых компонент " +
         std::to_string(mon.clamped_components()));
    check(mon.clamped_components() == 4, "маска yz у двух концов: 4 зажатые компоненты (x свободна)");
    // п. 1: y, z концов неизменны побитово (то же утверждает монитор), x растёт как x0 + v0 t
    check(dyz_end_max == 0.0, "y, z концевых бусин неизменны побитово (max отклонение " + fmt(dyz_end_max) + ")");
    const double travel = s.v0 * run.time;
    check(dx_end_max <= cf::ef5_x_tol * travel, "x концов = x0 + v0 t: max |dx| = " + fmt(dx_end_max) + " м <= " + fmt(cf::ef5_x_tol, 1) + " * v0 t (" + fmt(travel) + " м)");
    // п. 2: R_x == 0.0; импульс по x
    std::size_t nz = 0;
    double Ryz_max = 0.0;
    for (const auto& r : rec)
    {
        if (r.Rf[0] != 0.0) ++nz;
        Ryz_max = std::max(Ryz_max, std::max(std::fabs(r.Rf[1]), std::fabs(r.Rf[2])));
    }
    check(nz == 0, "R_x == 0.0 на каждом шаге (ненулевых: " + std::to_string(nz) + ")");
    check(Ryz_max > 0.0, "рельсы нагружены: max |R_y,z| = " + fmt(Ryz_max) + " Н > 0 (проверка не вырождена)");
    double Px = 0.0;
    for (std::size_t k = 0; k <= nseg; ++k) Px += s.model.fiber.mass[k] * s.model.fiber.vx[0][k];
    const double eps_acc = std::max(1.0e-12, 4.0 * 2.2e-16 * static_cast<double>(rec.size()));
    check(std::fabs(Px - Mtot * s.v0) <= eps_acc * Mtot * s.v0, "импульс по x сохраняется: P_x = " + fmt(Px) + " против M v0 = " + fmt(Mtot * s.v0) +
                                                                  " (отн. " + fmt(std::fabs(Px - Mtot * s.v0) / (Mtot * s.v0)) + ", допуск " + fmt(eps_acc, 1) + ")");
    // п. 3: центр масс
    double xcm = 0.0;
    for (std::size_t k = 0; k <= nseg; ++k) xcm += s.model.fiber.mass[k] * s.model.fiber.x[0][k];
    xcm /= Mtot;
    const double cm_err = std::fabs((xcm - xcm0) - s.v0 * run.time) / (s.v0 * run.time);
    check(cm_err <= eps_acc, "центр масс по x движется с v0: смещение " + fmt(xcm - xcm0) + " против v0 t = " + fmt(s.v0 * run.time) + " (отн. " + fmt(cm_err) + ", допуск " + fmt(eps_acc, 1) + ")");
}

// ===========================================================================
// EF6. Разбор CLAMP_BOX: сцена EF1 из временного .fib
// ===========================================================================
inline bool same_bits_d(double a, double b)
{
    std::uint64_t ua, ub;
    std::memcpy(&ua, &a, sizeof a);
    std::memcpy(&ub, &b, sizeof b);
    return ua == ub;
}

inline std::size_t count_scene_diffs(const fiber::Model& a, const fiber::Model& b)
{
    std::size_t diff = 0;
    if (a.fiber.size() != b.fiber.size() || a.topology.bonds.size() != b.topology.bonds.size()) return 1;
    for (std::size_t k = 0; k < a.fiber.size(); ++k)
    {
        for (std::size_t dd = 0; dd < 3; ++dd)
            if (!same_bits_d(a.fiber.x[dd][k], b.fiber.x[dd][k]) || !same_bits_d(a.fiber.vx[dd][k], b.fiber.vx[dd][k])) ++diff;
        if (!same_bits_d(a.fiber.mass[k], b.fiber.mass[k])) ++diff;
        if (a.fiber.object_id[k] != b.fiber.object_id[k] || a.fiber.material_id[k] != b.fiber.material_id[k]) ++diff;
        if (a.topology.prev[k] != b.topology.prev[k] || a.topology.next[k] != b.topology.next[k]) ++diff;
    }
    for (std::size_t i = 0; i < a.topology.bonds.size(); ++i)
    {
        const auto& p = a.topology.bonds[i];
        const auto& q = b.topology.bonds[i];
        if (p.i != q.i || p.j != q.j || !same_bits_d(p.l0, q.l0) || !same_bits_d(p.l_break, q.l_break)) ++diff;
    }
    const auto& ma = a.materials[0];
    const auto& mb = b.materials[0];
    if (!same_bits_d(ma.E, mb.E) || !same_bits_d(ma.nu, mb.nu) || !same_bits_d(ma.d, mb.d) || !same_bits_d(ma.rho, mb.rho) ||
        !same_bits_d(ma.eps_break, mb.eps_break) || !same_bits_d(ma.mu_ff, mb.mu_ff))
        ++diff;
    return diff;
}

inline void run_EF6(const Options& opt)
{
    namespace cf = clamp_case::fib;
    const bool pc = opt.solver == 0;
    const std::size_t nseg = opt.full ? cf::ef1_nseg_full : cf::ef1_nseg_smoke;
    const double pace = pc ? 1.0 : cf::euler_pacemaker;
    std::printf("EF6: сцена EF1 (%zu сегментов) из кода и из временного .fib с CLAMP_BOX\n", nseg);

    // --- сцена из кода
    Ef1Scene code = build_ef1(nseg, cf::ef1_v0_factor, fc::eps_break, pace);

    // --- файл: те же координаты (%.17g), тот же материал, CLAMP_BOX на первые nclamp бусин
    namespace fs = std::filesystem;
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const fs::path fib_path = fs::temp_directory_path() / ("sph_e2e_clamp_ef6_" + std::to_string(static_cast<long long>(stamp)) + ".fib");
    {
        std::ofstream f(fib_path);
        char buf[512];
        const fiber::FiberMaterial mat = make_mat();
        std::snprintf(buf, sizeof buf, "# EF6: сцена EF1 для разбора CLAMP_BOX\nMATERIAL 0 %.17g %.17g %.17g %.17g %.17g\n", mat.E, mat.nu, mat.d, mat.rho, mat.eps_break);
        f << buf;
        f << "FIBER 0 " << (nseg + 1) << "\n";
        for (std::size_t k = 0; k <= nseg; ++k)
        {
            std::snprintf(buf, sizeof buf, "%.17g %.17g %.17g\n", code.model.fiber.x[0][k], code.model.fiber.x[1][k], code.model.fiber.x[2][k]);
            f << buf;
        }
        // бокс накрывает бусины 0..nclamp-1 (середина между последней зажатой и первой свободной -- граница)
        const double xhi = (static_cast<double>(code.nclamp) - 0.5) * fc::l0;
        std::snprintf(buf, sizeof buf, "CLAMP_BOX %.17g %.17g %.17g %.17g %.17g %.17g\n", -0.5 * fc::l0, -1.0e-3, -1.0e-3, xhi, 1.0e-3, 1.0e-3);
        f << buf;
    }
    fiber::Model fm;
    bool loaded = true;
    try
    {
        init_model(fm);
        fiber::load_fiber_file(fib_path.string(), fm);
    }
    catch (const std::exception& ex)
    {
        loaded = false;
        check(false, std::string("load_fiber_file: исключение: ") + ex.what());
    }
    std::error_code ec;
    fs::remove(fib_path, ec);
    if (!loaded) return;

    // маска из разбора файла == маска из кода == mark_clamp_box по тем же координатам
    std::vector<std::uint8_t> box_mask(code.model.fiber.size(), 0);
    {
        rhs::CpuExecutor exec(1);
        const double lo[3] = {-0.5 * fc::l0, -1.0e-3, -1.0e-3};
        const double hi[3] = {(static_cast<double>(code.nclamp) - 0.5) * fc::l0, 1.0e-3, 1.0e-3};
        interactions::mark_clamp_box(exec, code.model.fiber.x, code.model.fiber.size(), lo, hi, interactions::kClampAll, box_mask);
    }
    box_mask.resize(fm.fiber.size(), 0);
    std::size_t nclamped_file = 0;
    for (auto b : fm.clamps.fiber_mask) nclamped_file += (b != 0);
    check(fm.clamps.fiber_mask == box_mask, "маска из .fib побитово совпадает с mark_clamp_box по тем же координатам");
    check(fm.clamps.fiber_mask.size() == fm.fiber.size() && nclamped_file == code.nclamp,
          "в файле закреплено " + std::to_string(nclamped_file) + " бусин (ожидалось " + std::to_string(code.nclamp) + "), размер маски == nfiber");
    // маска кода = первые nclamp бусин xyz
    {
        bool same = true;
        for (std::size_t k = 0; k < fm.fiber.size(); ++k)
            if (fm.clamps.fiber_mask[k] != (k < code.nclamp ? interactions::kClampAll : 0)) same = false;
        check(same, "маска из .fib == маска сцены EF1 из кода (бусины 0.." + std::to_string(code.nclamp - 1) + ", xyz)");
    }

    // файл -> сцена: скорости и пейсмейкер, как в build_ef1
    Ef1Scene file_scene;
    ef1_fill_params(file_scene, nseg, cf::ef1_v0_factor, fc::eps_break);
    file_scene.model = std::move(fm);
    ef1_set_velocities(file_scene);
    add_pacemaker(file_scene.model, pace, 0, 1);
    file_scene.model.clamps.fiber_mask.resize(file_scene.model.fiber.size(), 0);

    const std::size_t d0 = count_scene_diffs(code.model, file_scene.model);
    check(d0 == 0, "сцена из .fib побитово совпадает с кодовой: x, v, m, связи (l0, l_break), топология, материал (расхождений " + std::to_string(d0) + ")");
    check(code.model.clamps.fiber_mask == file_scene.model.clamps.fiber_mask, "маски кодовой и файловой сцен (с пейсмейкером) совпадают");

    // --- два прогона
    const double t_end = cf::ef1_tend_factor * code.L / code.c;
    ClampMonitor mon_c(code.model, 0, opt.solver), mon_f(file_scene.model, 0, opt.solver);
    const RunOutcome run_c = run_scene(opt, code.model, t_end, "ef6_code", mon_c, nullptr);
    const RunOutcome run_f = run_scene(opt, file_scene.model, t_end, "ef6_file", mon_f, nullptr);
    common_checks(" EF6 (из кода)", code.model, 0, mon_c, run_c, cf::balance_rel * mon_c.max_abs_momentum());
    common_checks(" EF6 (из .fib)", file_scene.model, 0, mon_f, run_f, cf::balance_rel * mon_f.max_abs_momentum());
    std::printf("\n[EF6] Критерии теста\n");
    const auto& rc = mon_c.records();
    const auto& rf = mon_f.records();
    bool recs_equal = rc.size() == rf.size();
    std::size_t first_bad = 0;
    for (std::size_t n = 0; recs_equal && n < rc.size(); ++n)
    {
        bool eq = same_bits_d(rc[n].t, rf[n].t) && same_bits_d(rc[n].dt, rf[n].dt) && same_bits_d(rc[n].Ekin, rf[n].Ekin) && same_bits_d(rc[n].Ebond, rf[n].Ebond);
        for (std::size_t dd = 0; dd < 3; ++dd) eq = eq && same_bits_d(rc[n].Rf[dd], rf[n].Rf[dd]) && same_bits_d(rc[n].Pf[dd], rf[n].Pf[dd]);
        if (!eq) { recs_equal = false; first_bad = n; }
    }
    check(recs_equal, "монитор побитово: t, dt, R_fiber, P, E_kin, E_bond на всех " + std::to_string(rc.size()) + " шагах" +
                          (recs_equal ? std::string() : " (расхождение на шаге " + std::to_string(first_bad + 1) + ")"));
    const std::size_t d1 = count_scene_diffs(code.model, file_scene.model);
    check(d1 == 0, "конечное состояние побитово: x, v всех бусин, топология (расхождений " + std::to_string(d1) + ")");
    check(rc.size() > 0 && !rc.empty() && rc.back().t > 0.0, "прогон не вырожден: " + std::to_string(rc.size()) + " шагов, max |R_x| = " + fmt([&] {
              double m = 0.0;
              for (const auto& r : rc) m = std::max(m, std::fabs(r.Rf[0]));
              return m;
          }()));
}

} // namespace fibtest

using fibtest::run_EF1;
using fibtest::run_EF2;
using fibtest::run_EF4;
using fibtest::run_EF5;
using fibtest::run_EF6;

} // namespace clamp_e2e
