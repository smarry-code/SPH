#pragma once
// tests/clamp_integration_solid.hpp
// Интеграционные тесты S1-S4 кинематического закрепления Solid через боевой
// time_integration(cfg, fiber::Model&, ...) (fiber/TZ_CLAMP_SOLID_FIBER.md, §10).
// Включается из tests/test_clamp.cpp ПОСЛЕ определения check()/failures/same_bits
// (анонимное пространство имён того файла), поэтому использует их напрямую.
//
//   S1  плита в несколько слоёв, полоса у кромки зажата целиком; PC и Эйлер:
//       x зажатых побитово, v == 0, «живой» зажим (напряжения/плотность
//       продолжают интегрироваться), свободная часть движется.
//   S2  скольжение нити по плите (сцена solid_kinematic из test_e2e_fiber_solid):
//       плита зажата маской и считается боевым циклом; сравнение замедления
//       нити с кинематической плитой.
//   S3  частичные маски (одна компонента).
//   S4  баланс импульса P^{n+1} - P^n = h*R_solid при cfg.epsilon = 0.
//
// Всё, что не относится к API теста, лежит в namespace clamp_solid_it, чтобы
// не пересекаться с хелперами соседних заголовков (clamp_integration_fiber.hpp).

#include "fiber/fiber_system.hpp"
#include "fiber/model.hpp"
#include "int_force.hpp"
#include "interactions/clamp.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_friction.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "materials.hpp"
#include "ode/ode_solver.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
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
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace clamp_solid_it
{

using interactions::kClampAll;

// ---------------------------------------------------------------------------
// Плита и конфиг
// ---------------------------------------------------------------------------
// Параллелепипед nx*ny*nz частиц алюминия (itype = 4, как build_plate), шаг dx,
// hsml = 2*dx (как в bc2_mini3d_geometry.hpp: носитель кубического сплайна 2h).
// Индекс частицы: (ix*ny + iy)*nz + iz. Центр по y, z в нуле, x от 0.
struct Plate
{
    std::size_t nx = 0, ny = 0, nz = 0;
    double dx = 0.0;
    std::size_t n() const { return nx * ny * nz; }
    std::size_t index(std::size_t ix, std::size_t iy, std::size_t iz) const { return (ix * ny + iy) * nz + iz; }
    std::size_t column(std::size_t i) const { return i / (ny * nz); } // ix частицы
};

constexpr double kRho = 2700.0;

inline Plate build_plate_layers(ParticleSystem& sys, std::size_t nx, std::size_t ny, std::size_t nz, double dx,
                                std::size_t maxn)
{
    Plate pl;
    pl.nx = nx; pl.ny = ny; pl.nz = nz; pl.dx = dx;
    sys.resize(maxn);
    for (std::size_t i = 0; i < maxn; ++i) sys.inzializ_index(i);
    for (std::size_t ix = 0; ix < nx; ++ix)
        for (std::size_t iy = 0; iy < ny; ++iy)
            for (std::size_t iz = 0; iz < nz; ++iz)
            {
                const std::size_t k = pl.index(ix, iy, iz);
                sys.x[0][k] = (static_cast<double>(ix) + 0.5) * dx;
                sys.x[1][k] = (static_cast<double>(iy) - 0.5 * static_cast<double>(ny) + 0.5) * dx;
                sys.x[2][k] = (static_cast<double>(iz) - 0.5 * static_cast<double>(nz) + 0.5) * dx;
                sys.itype[k] = 4; // алюминий
                sys.body[k] = 4;
                sys.rho[k] = kRho;
                sys.rho_0[k] = kRho;
                sys.mass[k] = kRho * dx * dx * dx;
                sys.hsml[k] = 2.0 * dx;
            }
    return pl;
}

// Конфиг: значения по умолчанию, Solid-SPH без контакта Solid-Solid, без
// внешних сил (в 3D нет ни виртуальных частиц -- build_ghost_set работает
// только при axis_sim && dim == 2, -- ни гравитации: в single_step/
// time_integration внешних сил нет). save_t/print_step велики: кадр выводится
// только на первом шаге прогона (в рабочий каталог ScopedWorkDir).
inline SimulationConfig make_cfg(int ode_solver, std::size_t maxn)
{
    SimulationConfig cfg;
    cfg.apply_defaults();
    cfg.ode_solver = ode_solver;
    cfg.contact_force_type = 0;
    cfg.maxn = maxn;
    cfg.max_interaction = 100 * maxn;
    cfg.save_t = 1.0e30;
    cfg.print_step = 1000000000;
    cfg.t_total = 1.0e30;
    return cfg;
}

// Подготовка EOS-полей (p, c) как в test_e2e_fiber_solid / smoke.
inline void init_eos(const SimulationConfig& cfg, ParticleSystem& sys, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i) urs(cfg.mats, sys, i);
}

inline bool same_bits_u64(double a, double b)
{
    std::uint64_t ua, ub;
    std::memcpy(&ua, &a, sizeof a);
    std::memcpy(&ub, &b, sizeof b);
    return ua == ub;
}

// Запуск ровно steps шагов боевого цикла (новый решатель на каждый вызов).
// Возвращает true, если цикл завершён без исключения; dt -- как записал цикл.
inline bool run_steps(const SimulationConfig& cfg, fiber::Model& model, std::size_t n, std::size_t steps, double& dt,
                      double& time, std::string& err)
{
    RunGuard guard(0.0, 0.0, 0);
    guard.set_step_limit(steps);
    e2e::ScopedWorkDir wd("sph_clamp_solid", "it");
    e2e::ScopedSilence silence(true);
    try
    {
        time_integration(cfg, model, n, dt, time, /*delta_x=*/0.0, &guard);
    }
    catch (const std::exception& ex)
    {
        err = ex.what();
        return false;
    }
    return guard.reason() == RunGuard::Reason::StepLimit && guard.steps() == steps;
}

inline bool solid_finite(const ParticleSystem& s, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d)
            if (!std::isfinite(s.x[d][i]) || !std::isfinite(s.vx[d][i])) return false;
        if (!std::isfinite(s.rho[i]) || !std::isfinite(s.u[i]) || !std::isfinite(s.sxx[i]) ||
            !std::isfinite(s.syy[i]) || !std::isfinite(s.szz[i]) || !std::isfinite(s.sxy[i]) ||
            !std::isfinite(s.syz[i]) || !std::isfinite(s.sxz[i]))
            return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Сцена S1/S4: консоль -- плита 20x8x4 (640 частиц, dx = 1 мм, 4 слоя по
// толщине), колонны ix < kClampCols зажаты целиком по ВСЕЙ толщине и ширине.
// Ширина полосы 6 dx = 6 мм >= носитель ядра 2h = 4 dx = 4 мм.
// Свободной части (ix >= kClampCols) задана линейно нарастающая от линии зажима
// скорость: поперечная vz = -8*s м/с и продольная vx = 2*s м/с, s = 0 у линии
// зажима и 1 на свободном торце (непрерывное поле: сильного разрыва скорости у
// зажима нет, но изгиб и растяжение консоли есть). c_Al ~ 5.4 км/с, dt ~ 4e-8 с.
// ---------------------------------------------------------------------------
constexpr std::size_t kNx = 20, kNy = 8, kNz = 4;
constexpr std::size_t kClampCols = 6;
constexpr double kDx = 1.0e-3;

struct Scene
{
    fiber::Model model;
    Plate pl;
    std::size_t n = 0, maxn = 0;
    std::vector<double> x0[dim];
    SimulationConfig cfg;
};

inline void build_cantilever(Scene& sc, int ode_solver, double epsilon, bool with_clamps)
{
    sc.n = kNx * kNy * kNz;
    sc.maxn = static_cast<std::size_t>(static_cast<double>(sc.n) * 1.4);
    sc.cfg = make_cfg(ode_solver, sc.maxn);
    sc.cfg.epsilon = epsilon;
    sc.pl = build_plate_layers(sc.model.solid, kNx, kNy, kNz, kDx, sc.maxn);
    sc.model.fiber.resize(0);
    sc.model.topology.resize(0);
    for (std::size_t ix = kClampCols; ix < kNx; ++ix)
    {
        const double s = static_cast<double>(ix - kClampCols + 1) / static_cast<double>(kNx - kClampCols);
        for (std::size_t iy = 0; iy < kNy; ++iy)
            for (std::size_t iz = 0; iz < kNz; ++iz)
            {
                const std::size_t k = sc.pl.index(ix, iy, iz);
                sc.model.solid.vx[2][k] = -8.0 * s;
                sc.model.solid.vx[0][k] = 2.0 * s;
            }
    }
    init_eos(sc.cfg, sc.model.solid, sc.n);
    if (with_clamps)
    {
        sc.model.clamps.solid_mask.assign(sc.n, 0);
        for (std::size_t i = 0; i < sc.n; ++i)
            if (sc.pl.column(i) < kClampCols) sc.model.clamps.solid_mask[i] = kClampAll;
    }
    for (std::size_t d = 0; d < dim; ++d) sc.x0[d].assign(sc.model.solid.x[d].begin(), sc.model.solid.x[d].begin() + sc.n);
}

inline double initial_dt(const Scene& sc)
{
    return compute_cfl_timestep(sc.cfg, const_cast<ParticleSystem&>(sc.model.solid), sc.n);
}

// ---------------------------------------------------------------------------
// S1
// ---------------------------------------------------------------------------
inline void s1_variant(int ode_solver, const char* name)
{
    std::printf("S1 (%s). консоль 20x8x4, зажато %zu колонн, 150 шагов\n", name, kClampCols);
    Scene sc;
    build_cantilever(sc, ode_solver, /*epsilon=*/0.3, /*with_clamps=*/true);
    double dt = initial_dt(sc), time = 0.0;
    std::string err;
    const auto t0 = std::chrono::steady_clock::now();
    const bool ok = run_steps(sc.cfg, sc.model, sc.n, 150, dt, time, err);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    check(ok, std::string("150 шагов без исключений") + (ok ? "" : ": " + err));
    if (!ok) return;
    std::printf("  info : время прогона %.2f с, dt = %.3e, t = %.3e\n", sec, dt, time);

    const ParticleSystem& s = sc.model.solid;
    bool x_ok = true, v_ok = true;
    std::size_t nclamped = 0;
    for (std::size_t i = 0; i < sc.n; ++i)
    {
        if (sc.pl.column(i) >= kClampCols) continue;
        ++nclamped;
        for (std::size_t d = 0; d < dim; ++d)
        {
            x_ok = x_ok && same_bits_u64(s.x[d][i], sc.x0[d][i]);
            v_ok = v_ok && (s.vx[d][i] == 0.0);
        }
    }
    check(nclamped == kClampCols * kNy * kNz, "зажато столько частиц, сколько задано");
    check(x_ok, "x зажатых побитово равны начальным");
    check(v_ok, "v зажатых == 0.0");
    check(solid_finite(s, sc.n), "нет NaN/Inf в x, v, rho, u, напряжениях");

    // Живой зажим: у линии зажима (последняя зажатая колонна) напряжения и
    // плотность интегрируются (ТЗ §3.4).
    double smax = 0.0, drho_max = 0.0;
    const std::size_t line = kClampCols - 1;
    for (std::size_t iy = 0; iy < kNy; ++iy)
        for (std::size_t iz = 0; iz < kNz; ++iz)
        {
            const std::size_t k = sc.pl.index(line, iy, iz);
            smax = std::max({smax, std::fabs(s.sxx[k]), std::fabs(s.syy[k]), std::fabs(s.szz[k]), std::fabs(s.sxy[k]),
                             std::fabs(s.syz[k]), std::fabs(s.sxz[k])});
            drho_max = std::max(drho_max, std::fabs(s.rho[k] - kRho));
        }
    std::printf("  info : у линии зажима max|s| = %.3e Па, max|rho-rho0| = %.3e кг/м3\n", smax, drho_max);
    check(smax > 0.0, "напряжения в зажатых частицах у линии зажима ненулевые");
    check(drho_max > 0.0, "плотность зажатых у линии зажима изменилась (живой зажим)");

    // Свободная часть сдвинулась; торец прогнулся вниз (vz < 0 в начале).
    double zshift = 0.0;
    bool moved = true;
    for (std::size_t iy = 0; iy < kNy; ++iy)
        for (std::size_t iz = 0; iz < kNz; ++iz)
        {
            const std::size_t k = sc.pl.index(kNx - 1, iy, iz);
            zshift = std::min(zshift, s.x[2][k] - sc.x0[2][k]);
        }
    for (std::size_t i = 0; i < sc.n; ++i)
        if (sc.pl.column(i) >= kClampCols && s.x[2][i] == sc.x0[2][i] && s.x[0][i] == sc.x0[0][i]) moved = false;
    std::printf("  info : смещение торца по z = %.3e м\n", zshift);
    check(moved, "все свободные частицы сдвинулись");
    check(zshift < -1.0e-6, "свободный торец прогнулся вниз (dz < -1 мкм)");

    const auto& rep = sc.model.clamp_report;
    const double rn = std::sqrt(rep.solid_force[0] * rep.solid_force[0] + rep.solid_force[1] * rep.solid_force[1] +
                                rep.solid_force[2] * rep.solid_force[2]);
    std::printf("  info : реакция последнего шага R_solid = (%.4e, %.4e, %.4e) Н\n", rep.solid_force[0],
                rep.solid_force[1], rep.solid_force[2]);
    check(std::isfinite(rn) && rn > 0.0, "реакция последнего шага конечна и ненулевая");
    check(rep.solid_clamped == nclamped, "clamp_report.solid_clamped == числу зажатых");
}

inline void test_s1()
{
    s1_variant(0, "PC");
    s1_variant(1, "Эйлер");
}

// ---------------------------------------------------------------------------
// S3: частичные маски. Плита 20x8x4; область A (ix < 6) закреплена по z (бит 2),
// область B (ix >= 14) -- по x (бит 0). В A начальная vx = +3 м/с, vy = +1 м/с
// (свободные компоненты, validate_clamps это допускает), в B vy = -2 м/с;
// середина (6 <= ix < 14) получает vz = -5 м/с и vx = 2 м/с. Проверяется: закреплённая
// координата побитово неизменна, скорость по ней == 0, по остальным частицы движутся.
// ---------------------------------------------------------------------------
inline void s3_variant(int ode_solver, const char* name)
{
    std::printf("S3 (%s). частичные маски: A (ix<6) -- только z, B (ix>=14) -- только x, 100 шагов\n", name);
    Scene sc;
    build_cantilever(sc, ode_solver, 0.3, /*with_clamps=*/false);
    ParticleSystem& s = sc.model.solid;
    const std::uint8_t bit_x = 1u << 0, bit_z = 1u << 2;
    sc.model.clamps.solid_mask.assign(sc.n, 0);
    for (std::size_t i = 0; i < sc.n; ++i)
        for (std::size_t d = 0; d < dim; ++d) s.vx[d][i] = 0.0;
    for (std::size_t i = 0; i < sc.n; ++i)
    {
        const std::size_t c = sc.pl.column(i);
        if (c < 6)
        {
            sc.model.clamps.solid_mask[i] = bit_z;
            s.vx[0][i] = 3.0;
            s.vx[1][i] = 1.0;
        }
        else if (c >= 14)
        {
            sc.model.clamps.solid_mask[i] = bit_x;
            s.vx[1][i] = -2.0;
        }
        else
        {
            s.vx[2][i] = -5.0;
            s.vx[0][i] = 2.0;
        }
    }
    double dt = initial_dt(sc), time = 0.0;
    std::string err;
    const auto t0 = std::chrono::steady_clock::now();
    const bool ok = run_steps(sc.cfg, sc.model, sc.n, 100, dt, time, err);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    check(ok, std::string("100 шагов без исключений") + (ok ? "" : ": " + err));
    if (!ok) return;
    std::printf("  info : время прогона %.2f с, t = %.3e\n", sec, time);

    bool a_z_ok = true, a_vz_ok = true, b_x_ok = true, b_vx_ok = true;
    bool a_x_moved = true, a_y_moved = true, b_y_moved = true, b_z_moved = true;
    for (std::size_t i = 0; i < sc.n; ++i)
    {
        const std::size_t c = sc.pl.column(i);
        if (c < 6)
        {
            a_z_ok = a_z_ok && same_bits_u64(s.x[2][i], sc.x0[2][i]);
            a_vz_ok = a_vz_ok && (s.vx[2][i] == 0.0);
            a_x_moved = a_x_moved && (s.x[0][i] != sc.x0[0][i]);
            a_y_moved = a_y_moved && (s.x[1][i] != sc.x0[1][i]);
        }
        else if (c >= 14)
        {
            b_x_ok = b_x_ok && same_bits_u64(s.x[0][i], sc.x0[0][i]);
            b_vx_ok = b_vx_ok && (s.vx[0][i] == 0.0);
            b_y_moved = b_y_moved && (s.x[1][i] != sc.x0[1][i]);
            b_z_moved = b_z_moved && (s.x[2][i] != sc.x0[2][i]);
        }
    }
    check(a_z_ok && a_vz_ok, "A: z побитово неизменна, vz == 0.0");
    check(a_x_moved && a_y_moved, "A: x и y изменились у всех частиц");
    check(b_x_ok && b_vx_ok, "B: x побитово неизменна, vx == 0.0");
    check(b_y_moved && b_z_moved, "B: y и z изменились у всех частиц");
    check(solid_finite(s, sc.n), "нет NaN/Inf");
}

inline void test_s3()
{
    s3_variant(0, "PC");
    s3_variant(1, "Эйлер");
}

// ---------------------------------------------------------------------------
// S4: баланс импульса. Сцена S1 с cfg.epsilon = 0 (XSPH не входит в v, а av
// зажатых в [1] обнуляется, так что баланс без XSPH точен). Шаги -- по одному
// вызовом time_integration с set_step_limit(1) на той же модели: решатель
// создаётся заново, поэтому у PC КАЖДЫЙ вызов -- «первый шаг», h = dt/2, а у Эйлера
// h = dt; dt берётся из записанного циклом значения. Для шага k:
//     dP_k = P^{k+1} - P^k  (P = sum m v по ВСЕМ частицам, зажатые дают 0),
//     dP_k =? h * R_solid,  R_solid -- model.clamp_report.solid_force.
// Внешних сил нет: contact_force_type = 0, виртуальных частиц в 3D нет, гравитации
// в цикле нет (см. make_cfg).
//
// ДОПУСК. Точный баланс требует sum_i m_i a_i = 0 по ВСЕМ частицам до проекции,
// т.е. попарной антисимметрии SPH-сил. Измерение на том же прогоне БЕЗ зажимов
// (dP_base = h sum m a) показывает, что SPH-силы проекта антисимметричны до
// округления: |dP_base| / S ~ 2e-15 (измерено, 8 шагов PC и Эйлер) при S = sum m|v^k| + sum m|v^{k+1}|
// (масштаб импульсов, по которым берётся разность). Значит остаток -- это
// округление: P -- сумма 640 слагаемых порядка m*v, а dP мал по сравнению с P
// (|hR|/S ~ 1e-5), поэтому остаток относительно |hR| выглядит как 1e-10.
// Допуск поэтому абсолютный от масштаба S: 1e-13*S (~450 машинных эпсилон; запас
// на суммирование 640 слагаемых и на накопление SPH-сумм). Он на ~8 порядков
// меньше |hR|/S ~ 1e-5, то есть пропущенная или неверно посчитанная реакция
// (ошибка порядка 1e-8 от R) тест бы поймала. Запас допуска над измеренным ~50.
struct MomentumStats
{
    double max_err_S = 0.0;      // max_k |dP_k - h R_k| / S_k  (с зажимом)
    double max_err_R = 0.0;      // max_k |dP_k - h R_k| / |h R_k|
    double min_R_over_S = 1e300; // min_k |h R_k| / S_k  (чувствительность проверки)
    double max_base_S = 0.0;     // max_k |dP_k| / S_k  (без зажимов, базовая неконсервативность)
};

inline void momentum(const ParticleSystem& s, std::size_t n, double P[dim], double* abs_sum = nullptr)
{
    for (std::size_t d = 0; d < dim; ++d) P[d] = 0.0;
    double a = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        double v2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            P[d] += s.mass[i] * s.vx[d][i];
            v2 += s.vx[d][i] * s.vx[d][i];
        }
        a += s.mass[i] * std::sqrt(v2);
    }
    if (abs_sum) *abs_sum = a;
}

inline double norm3(const double v[dim])
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

inline MomentumStats s4_run(int ode_solver, int steps, bool with_clamps)
{
    Scene sc;
    build_cantilever(sc, ode_solver, /*epsilon=*/0.0, with_clamps);
    MomentumStats st;
    double dt = initial_dt(sc), time = 0.0;
    std::string err;
    for (int k = 0; k < steps; ++k)
    {
        double P0[dim], P1[dim], S0 = 0.0, S1 = 0.0;
        momentum(sc.model.solid, sc.n, P0, &S0);
        const bool ok = run_steps(sc.cfg, sc.model, sc.n, 1, dt, time, err);
        check(ok, "шаг " + std::to_string(k + 1) + " выполнен" + (ok ? "" : ": " + err));
        if (!ok) return st;
        momentum(sc.model.solid, sc.n, P1, &S1);
        const double h = (ode_solver == 0) ? 0.5 * dt : dt;
        const double S = S0 + S1;
        double dP[dim], diff[dim];
        for (std::size_t d = 0; d < dim; ++d) dP[d] = P1[d] - P0[d];
        if (with_clamps)
        {
            double hR[dim];
            for (std::size_t d = 0; d < dim; ++d)
            {
                hR[d] = h * sc.model.clamp_report.solid_force[d];
                diff[d] = dP[d] - hR[d];
            }
            const double nR = norm3(hR), nD = norm3(diff);
            st.max_err_S = std::max(st.max_err_S, nD / S);
            st.max_err_R = std::max(st.max_err_R, nD / nR);
            st.min_R_over_S = std::min(st.min_R_over_S, nR / S);
            std::printf("  info : шаг %d: h = %.3e, |hR| = %.6e, |dP - hR| = %.3e (к |hR|: %.2e, к S: %.2e)\n", k + 1, h,
                        nR, nD, nD / nR, nD / S);
        }
        else
        {
            st.max_base_S = std::max(st.max_base_S, norm3(dP) / S);
            std::printf("  info : шаг %d (без зажимов): |dP| = %.3e (к S: %.2e)\n", k + 1, norm3(dP), norm3(dP) / S);
        }
    }
    return st;
}

inline void s4_variant(int ode_solver, const char* name)
{
    std::printf("S4 (%s). баланс импульса, cfg.epsilon = 0, 8 отдельных шагов\n", name);
    const MomentumStats c = s4_run(ode_solver, 8, true);
    const MomentumStats b = s4_run(ode_solver, 8, false);
    std::printf("  info : с зажимом: max |dP - hR|/S = %.3e, max |dP - hR|/|hR| = %.3e, min |hR|/S = %.3e; "
                "базовая (без зажимов) max |dP|/S = %.3e\n",
                c.max_err_S, c.max_err_R, c.min_R_over_S, b.max_base_S);
    const double tol = 1.0e-13;
    check(c.min_R_over_S > 1.0e-7, "реакция значима на каждом шаге (|hR|/S > 1e-7)");
    check(c.max_err_S <= tol, "|dP - h R| / S <= 1e-13 на каждом из 8 шагов");
    check(b.max_base_S <= tol, "базовая: без зажимов |dP| / S <= 1e-13 (SPH-силы антисимметричны до округления)");
}

inline void test_s4()
{
    s4_variant(0, "PC");
    s4_variant(1, "Эйлер");
}

// ---------------------------------------------------------------------------
// S2: нить скользит по плите (сцена test_e2e_fiber_solid.cpp, solid_kinematic =
// true, FrictionPlateIntegrator), плита зажата целиком маской и считается
// боевым time_integration.
//
// Сцена (копия build_fs_scene_plate): плита 10x10x1, dx = 1.5 мм, алюминий;
// нить из 5 бусин (d = 1 мм, E = 2e7, l0 = 2 мм), над плитой на зазоре 0.5 мм,
// vz = -0.5 м/с, vx = v_t0 = 2 м/с (нить стартует со сдвигом по x так, что к
// контакту оказывается над центром плиты), mu_fs = 0.2 (быстрое скольжение).
//
// Эталон -- kin_reference(): копия FrictionPlateIntegrator при free_plate =
// false (кинематическая плита: SolidVelocityUpdate{h = 0}, solid_kinematic =
// true), реальный PredictorCorrectorSolver, плита вне решателя.
//
// КАК СОГЛАСОВАНЫ УСЛОВИЯ. (1) Шаг dt. Эталон берёт dt = safety * min(лимит
// связи, лимит контакта Fiber-Solid); боевой цикл -- min(CFL по плите, лимит
// нити, лимит контакта, контактный лимит Solid) БЕЗ safety (встроенные 0.8 уже в
// лимитах). Поэтому (а) в боевом прогоне cfg.CFL = 1e6, чтобы Solid-CFL
// (0.1*h/c ~ 5e-8 с, плита неподвижна и к динамике нити отношения не имеет) не
// определял шаг и решал лимит нити, и (б) основное сравнение -- с эталоном при
// safety = 1, то есть с тем же dt; второе сравнение -- с safety = 0.3 (как в
// самом test_e2e_fiber_solid) показывает чувствительность к шагу.
// (2) SPH-сил в плите нет: состояние покоя однородно (напряжения 0, p(rho0, u=0) = 0
// для EOS алюминия, скорости 0), плита зажата, XSPH и dvx обнуляются.
// (3) Трение: c_s = 0 у зажатых (точка [2]) == solid_kinematic в эталоне.
// (4) Решатель -- PC в обоих прогонах; нить пересчитывается одной схемой.
// Измерение: vx0 - vx1 -- потеря средней (по массе) скорости нити вдоль x за
// T_END; к T_END контакт закончился (vz1 > 0, нить улетает), дальше vx не меняется.
// ---------------------------------------------------------------------------
constexpr double kS2Dxp = 1.5e-3;
constexpr double kS2TEnd = 5.0e-3;

inline fiber::FiberMaterial s2_fiber_material(double mu_fs)
{
    fiber::FiberMaterial m;
    m.E = 2.0e7; m.nu = 0.3; m.d = 1.0e-3; m.rho = 1200.0; m.eps_break = 1.0;
    m.mu_fs = mu_fs;
    return m;
}

struct S2Scene
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    double mass_bead = 0.0;
};

inline void s2_build_fiber(S2Scene& sc, double mu_fs, double v_t0)
{
    const double E_f = 2.0e7, d_f = 1.0e-3, rho_f = 1200.0;
    (void)E_f;
    const double r_s = 0.5 * kS2Dxp, r_f = 0.5 * d_f;
    const double l0 = 2.0e-3;
    const int N = 5;
    const double A = pi * d_f * d_f / 4.0;
    sc.mass_bead = rho_f * A * l0;
    const double gap0 = 0.5e-3, v0 = 0.5;
    const double z0 = r_s + r_f + gap0;
    const double z_first_contact = 1.0e-3;
    const double x_shift = -v_t0 * (z0 - z_first_contact) / v0;
    sc.materials.push_back(s2_fiber_material(mu_fs));
    sc.state.resize(N);
    sc.topology.resize(N);
    for (int k = 0; k < N; ++k)
    {
        const std::size_t gi = static_cast<std::size_t>(k);
        sc.state.x[0][gi] = static_cast<double>(k - N / 2) * l0 + x_shift;
        sc.state.x[1][gi] = 0.0;
        sc.state.x[2][gi] = z0;
        for (std::size_t dd = 0; dd < dim; ++dd) sc.state.vx[dd][gi] = 0.0;
        sc.state.vx[2][gi] = -v0;
        sc.state.vx[0][gi] = v_t0;
        sc.state.mass[gi] = sc.mass_bead;
        sc.state.object_id[gi] = 0;
        sc.state.material_id[gi] = 0;
        if (k > 0) sc.topology.prev[gi] = static_cast<long long>(gi - 1);
        if (k + 1 < N)
        {
            sc.topology.next[gi] = static_cast<long long>(gi + 1);
            fiber::FiberBond b;
            b.i = gi; b.j = gi + 1;
            b.l0 = l0;
            b.l_break = b.l0 * 1e9; // разрыв -- не предмет теста
            sc.topology.bonds.push_back(b);
        }
    }
}

// Плита S2: 10x10x1, центр плиты по x в нуле (как build_plate в эталонной сцене).
inline void s2_build_plate(ParticleSystem& plate, std::size_t maxn)
{
    build_plate_layers(plate, 10, 10, 1, kS2Dxp, maxn);
    for (std::size_t i = 0; i < 100; ++i) plate.x[0][i] -= 0.5 * 10.0 * kS2Dxp;
}

inline double s2_mean_v(const fiber::FiberState& fb, std::size_t d)
{
    double sp = 0.0, sm = 0.0;
    for (std::size_t i = 0; i < fb.size(); ++i) { sp += fb.mass[i] * fb.vx[d][i]; sm += fb.mass[i]; }
    return sp / sm;
}

struct S2Result
{
    double vx0 = 0.0, vx1 = 0.0, vz0 = 0.0, vz1 = 0.0;
    long steps = 0;
    double time = 0.0, sec = 0.0, dt_min = 1e9, dt_max = 0.0;
    double min_z = 1e9, max_pen = 0.0;
    bool ok = true;
    std::string err;
};

// Эталон: кинематическая плита, вне решателя. Копия FrictionPlateIntegrator::step
// без диагностик (free_plate = false).
inline S2Result s2_kin_reference(double mu_fs, double v_t0, double safety)
{
    S2Result res;
    S2Scene sc;
    s2_build_fiber(sc, mu_fs, v_t0);
    ParticleSystem plate;
    s2_build_plate(plate, 100);
    const std::size_t nplate = 100;
    fiber::FiberState& fb = sc.state;
    res.vx0 = s2_mean_v(fb, 0);
    res.vz0 = s2_mean_v(fb, 2);

    Materials smats;
    ode::PredictorCorrectorSolver solver;
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> bond_buffers, fs_fiber_buffers, fr_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fs_solid_buffers, fr_solid_buffers;
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> count_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidCountDelta>> count_solid_buffers;
    std::vector<double> bead_count, solid_count;
    std::vector<interactions::FiberSolidFrictionContact> slots;
    std::vector<interactions::FiberFrictionContact> no_ff_slots;
    std::vector<std::vector<double>> zero_dvx(dim, std::vector<double>(nplate, 0.0));
    const std::size_t n = fb.size();
    solver.reset(0, 0, n);
    double dt_prev = 0.0, time = 0.0;
    while (time < kS2TEnd)
    {
        solver.predict(plate, 0, fb, n, dt_prev);
        std::vector<fiber::FiberDelta> bond_delta = interactions::compute_fiber_bond(exec, fb, sc.topology, sc.materials, bond_buffers);
        const double bond_dt_limit = interactions::compute_fiber_bond_dt_limit(exec, fb, sc.topology, sc.materials);
        interactions::FiberSolidContactSearchResult fs_search;
        interactions::build_fiber_solid_contact_search(exec, smats, plate, nplate, fb, sc.topology, sc.materials, fs_search);
        const bool has_active = interactions::compute_fiber_solid_contact_contributions(
            exec, smats, plate, fb, sc.materials, fs_search.primitives, fs_fiber_buffers, fs_solid_buffers, slots);
        std::vector<fiber::FiberDelta> fs_delta = rhs::reduce_contributions(exec, fs_fiber_buffers, n);
        const double fs_dt_limit = interactions::compute_fiber_solid_contact_dt_limit(exec, smats, plate, fb, sc.materials,
                                                                                        fs_search.primitives);
        for (const auto& p : fs_search.primitives)
        {
            const interactions::FiberSolidContactGeometry g =
                interactions::compute_fiber_solid_contact_geometry(smats, plate, fb, sc.materials, p);
            if (g.contact) res.max_pen = std::max(res.max_pen, g.delta);
        }
        const double dt = safety * std::min(bond_dt_limit, fs_dt_limit);
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < n; ++i) solver.fiber_rhs.dv[d][i] = bond_delta[i].dv[d] + fs_delta[i].dv[d];
        if (has_active)
        {
            const auto fupd = solver.fiber_velocity_update(fb, dt);
            ode::SolidVelocityUpdate supd{};
            for (std::size_t d = 0; d < dim; ++d) supd.base[d] = &plate.vx[d];
            supd.h = 0.0;
            supd.av = nullptr;
            interactions::compute_friction_counts(exec, n, nplate, no_ff_slots, slots, count_fiber_buffers,
                                                   count_solid_buffers, bead_count, solid_count);
            interactions::compute_friction_contributions(exec, fb, fupd, solver.fiber_rhs.dv, no_ff_slots, plate, supd,
                                                          zero_dvx, slots, bead_count, solid_count, /*kin=*/true,
                                                          fr_fiber_buffers, fr_solid_buffers);
            const std::vector<fiber::FiberDelta> fr_f = rhs::reduce_contributions(exec, fr_fiber_buffers, n);
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i) solver.fiber_rhs.dv[d][i] += fr_f[i].dv[d];
        }
        solver.advance(plate, 0, fb, n, dt);
        dt_prev = dt;
        time += dt;
        ++res.steps;
        res.dt_min = std::min(res.dt_min, dt);
        res.dt_max = std::max(res.dt_max, dt);
        res.min_z = std::min(res.min_z, fb.x[2][n / 2]);
    }
    res.time = time;
    res.vx1 = s2_mean_v(fb, 0);
    res.vz1 = s2_mean_v(fb, 2);
    return res;
}

// Боевой прогон: плита в model.solid зажата kClampAll, нить в model.fiber.
inline S2Result s2_production(double mu_fs, double v_t0, int ode_solver, bool with_clamps,
                              std::vector<double> (*x0_out)[dim] = nullptr, fiber::Model* keep = nullptr)
{
    S2Result res;
    S2Scene sc;
    s2_build_fiber(sc, mu_fs, v_t0);
    fiber::Model local;
    fiber::Model& model = keep ? *keep : local;
    const std::size_t n = 100, maxn = 140;
    SimulationConfig cfg = make_cfg(ode_solver, maxn);
    cfg.CFL = 1.0e6; // шаг определяют лимиты нити, см. шапку секции S2
    cfg.t_total = kS2TEnd;
    s2_build_plate(model.solid, maxn);
    init_eos(cfg, model.solid, n);
    model.fiber = sc.state;
    model.topology = sc.topology;
    model.materials = sc.materials;
    if (with_clamps) model.clamps.solid_mask.assign(n, kClampAll);
    if (x0_out)
        for (std::size_t d = 0; d < dim; ++d) (*x0_out)[d].assign(model.solid.x[d].begin(), model.solid.x[d].begin() + n);
    res.vx0 = s2_mean_v(model.fiber, 0);
    res.vz0 = s2_mean_v(model.fiber, 2);

    double dt = compute_cfl_timestep(cfg, model.solid, n), time = 0.0;
    RunGuard guard(0.0, 0.0, 0);
    guard.set_step_limit(2000000);
    const auto t0 = std::chrono::steady_clock::now();
    {
        e2e::ScopedWorkDir wd("sph_clamp_solid", "s2");
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(cfg, model, n, dt, time, 0.0, &guard);
        }
        catch (const std::exception& ex)
        {
            res.ok = false;
            res.err = ex.what();
        }
    }
    res.sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    res.steps = static_cast<long>(guard.steps());
    res.time = time;
    res.vx1 = s2_mean_v(model.fiber, 0);
    res.vz1 = s2_mean_v(model.fiber, 2);
    return res;
}

// Малый боевой прогон сцены S2 с CFL по умолчанию (cfg.apply_defaults): проверяет,
// что зажатые частицы Solid остаются в выборе dt (TZ §3.4) -- в основном прогоне
// CFL = 1e6 отключает ограничение Solid ради совпадения dt с эталоном.
inline void s2_default_cfl_check()
{
    const double mu = 0.2, v_t0 = 2.0;
    const std::size_t n = 100, maxn = 140, nsteps = 20;
    S2Scene sc;
    s2_build_fiber(sc, mu, v_t0);
    fiber::Model model;
    SimulationConfig cfg = make_cfg(0, maxn); // CFL -- значение apply_defaults
    cfg.t_total = 1.0e30;
    s2_build_plate(model.solid, maxn);
    init_eos(cfg, model.solid, n);
    model.fiber = sc.state;
    model.topology = sc.topology;
    model.materials = sc.materials;
    model.clamps.solid_mask.assign(n, kClampAll);
    std::vector<double> x0[dim];
    for (std::size_t d = 0; d < dim; ++d) x0[d].assign(model.solid.x[d].begin(), model.solid.x[d].begin() + n);

    double dt = compute_cfl_timestep(cfg, model.solid, n), time = 0.0;
    const double cfl0 = dt;
    bool ok = true;
    std::string err;
    double dt_last = 0.0, cfl_before_last = 0.0;
    {
        e2e::ScopedWorkDir wd("sph_clamp_solid", "s2cfl");
        e2e::ScopedSilence silence(true);
        try
        {
            // nsteps-1 шагов, затем замер CFL плиты непосредственно перед последним шагом
            // (плита неподвижна, но c, rho зажатых живут -- CFL может дрейфовать).
            RunGuard g1(0.0, 0.0, 0);
            g1.set_step_limit(nsteps - 1);
            time_integration(cfg, model, n, dt, time, 0.0, &g1);
            cfl_before_last = compute_cfl_timestep(cfg, model.solid, n);
            RunGuard g2(0.0, 0.0, 0);
            g2.set_step_limit(1);
            time_integration(cfg, model, n, dt, time, 0.0, &g2);
            dt_last = dt;
        }
        catch (const std::exception& ex)
        {
            ok = false;
            err = ex.what();
        }
    }
    check(ok, std::string("S2 (CFL по умолчанию): прогон без исключений") + (ok ? "" : ": " + err));
    if (!ok) return;
    const double cfl_after = compute_cfl_timestep(cfg, model.solid, n);
    bool x_ok = true, v_ok = true;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < dim; ++d)
        {
            x_ok = x_ok && same_bits_u64(model.solid.x[d][i], x0[d][i]);
            v_ok = v_ok && (model.solid.vx[d][i] == 0.0);
        }
    bool fin = solid_finite(model.solid, n) && std::isfinite(dt_last);
    for (std::size_t i = 0; i < model.fiber.size(); ++i)
        for (std::size_t d = 0; d < dim; ++d)
            fin = fin && std::isfinite(model.fiber.x[d][i]) && std::isfinite(model.fiber.vx[d][i]);
    std::printf("  info : S2, CFL по умолчанию (cfg.CFL = %g): %zu шагов, compute_cfl_timestep(начало) = %.6e, "
                "(перед последним шагом) = %.6e, (после) = %.6e, dt последнего шага = %.6e (dt/CFL = %.6f)\n",
                cfg.CFL, nsteps, cfl0, cfl_before_last, cfl_after, dt_last, dt_last / cfl_before_last);
    check(x_ok && v_ok, "S2 (CFL по умолчанию): плита побитово неподвижна (x == x0, v == 0.0)");
    check(fin, "S2 (CFL по умолчанию): нет NaN/Inf");
    check(dt_last > 0.0 && dt_last <= cfl_before_last * (1.0 + 1e-12),
          "S2 (CFL по умолчанию): dt <= compute_cfl_timestep(плита): зажатые частицы Solid остаются в выборе dt");
    check(std::fabs(dt_last / cfl_before_last - 1.0) <= 1e-9,
          "S2 (CFL по умолчанию): dt == compute_cfl_timestep(плита) -- CFL зажатой плиты здесь и есть ограничивающий (лимиты нити выше)");
}

inline void test_s2()
{
    std::printf("S2. нить скользит по зажатой плите (PC): mu_fs = 0.2, v_t0 = 2 м/с, T = %.1e с\n", kS2TEnd);
    const double mu = 0.2, v_t0 = 2.0;

    const auto t0 = std::chrono::steady_clock::now();
    const S2Result ref1 = s2_kin_reference(mu, v_t0, 1.0);
    const S2Result ref03 = s2_kin_reference(mu, v_t0, 0.3);
    const double sec_ref = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::vector<double> x0[dim];
    fiber::Model model;
    const S2Result prod = s2_production(mu, v_t0, 0, true, &x0, &model);

    const double dec_ref1 = ref1.vx0 - ref1.vx1, dec_ref03 = ref03.vx0 - ref03.vx1, dec_prod = prod.vx0 - prod.vx1;
    std::printf("  info : эталон (кинематическая плита, safety=1.0): %ld шагов, dt=[%.3e,%.3e], vx %.6f -> %.6f (торможение %.6f), "
                "vz %.4f -> %.4f\n", ref1.steps, ref1.dt_min, ref1.dt_max, ref1.vx0, ref1.vx1, dec_ref1, ref1.vz0, ref1.vz1);
    std::printf("  info : эталон (safety=0.3): %ld шагов, vx %.6f -> %.6f (торможение %.6f), vz -> %.4f; %.2f с на оба эталона\n",
                ref03.steps, ref03.vx0, ref03.vx1, dec_ref03, ref03.vz1, sec_ref);
    std::printf("  info : боевой time_integration + зажим: %ld шагов, vx %.6f -> %.6f (торможение %.6f), vz %.4f -> %.4f, %.2f с\n",
                prod.steps, prod.vx0, prod.vx1, dec_prod, prod.vz0, prod.vz1, prod.sec);
    check(prod.ok, std::string("боевой прогон без исключений") + (prod.ok ? "" : ": " + prod.err));
    if (!prod.ok) return;
    check(ref1.vz1 > 0.0 && ref03.vz1 > 0.0 && prod.vz1 > 0.0, "контакт закончен, нить отскочила (vz > 0) во всех прогонах");
    check(dec_prod > 0.1, "трение тормозит нить (торможение > 0.1 м/с)");
    const double rel1 = std::fabs(dec_prod - dec_ref1) / dec_ref1;
    const double rel03 = std::fabs(dec_prod - dec_ref03) / dec_ref03;
    std::printf("  info : |торможение боевое - эталон| / эталон: %.3e (тот же dt), %.3e (safety 0.3)\n", rel1, rel03);
    check(rel1 < 0.05, "замедление нити совпадает с кинематической плитой (тот же dt) в пределах 5 %");
    check(rel03 < 0.05, "замедление нити совпадает с эталоном при safety = 0.3 в пределах 5 %");

    const ParticleSystem& s = model.solid;
    bool x_ok = true, v_ok = true;
    for (std::size_t i = 0; i < 100; ++i)
        for (std::size_t d = 0; d < dim; ++d)
        {
            x_ok = x_ok && same_bits_u64(s.x[d][i], x0[d][i]);
            v_ok = v_ok && (s.vx[d][i] == 0.0);
        }
    check(x_ok, "плита: x побитово равны начальным");
    check(v_ok, "плита: v == 0.0");
    bool fin = solid_finite(s, 100);
    for (std::size_t i = 0; i < model.fiber.size(); ++i)
        for (std::size_t d = 0; d < dim; ++d) fin = fin && std::isfinite(model.fiber.x[d][i]) && std::isfinite(model.fiber.vx[d][i]);
    check(fin, "нет NaN/Inf (плита и нить)");
    const auto& rep = model.clamp_report;
    check(rep.solid_clamped == 100, "clamp_report.solid_clamped == 100");

    s2_default_cfl_check();
}

} // namespace clamp_solid_it
inline void run_clamp_solid_integration_tests()
{
    std::printf("\n=== Solid: интеграционные тесты S1-S4 ===\n");
    clamp_solid_it::test_s1();
    clamp_solid_it::test_s2();
    clamp_solid_it::test_s3();
    clamp_solid_it::test_s4();
}
