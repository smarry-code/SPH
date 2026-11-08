// tests/test_ode_solver.cpp
// Unit-тесты OdeSolver (fiber/PLAN.md, этап Ф8; §3.9 «ODESolver»):
//   1) Euler -- первый порядок сходимости на задаче с аналитическим решением
//      (гармонический осциллятор x'' = -k*x, x(t) = x0*cos(omega*t));
//   2) PredictorCorrector -- совпадение с прежним циклом на BC2/BC6. Это уже
//      доказывается bit-exact-регрессией (tests/test_regression.cpp,
//      `ctest -L regression`) и отдельно -- прямым сравнением старой и новой
//      сигнатуры time_integration на непустом Solid
//      (tests/test_fiber_time_integration_equivalence.cpp), которое ПОСЛЕ
//      унификации Ф8 сравнивает две обёртки над ОДНОЙ и той же реализацией
//      (detail::run_time_integration<ode::PredictorCorrectorSolver>,
//      src/time_integration.cpp), то есть является даже более сильной
//      проверкой, чем до Ф8. Отдельный дублирующий прогон здесь не добавлен
//      осознанно (fiber/PLAN.md §Ф8, раздел G: «если считаешь регрессию
//      достаточным покрытием этого пункта, можешь не дублировать, но явно
//      поясни»).
//   3) нить без контактов: полная механическая энергия (кинетическая +
//      упругая потенциальная в связях) ограничена на PredictorCorrector на
//      dt, при котором Эйлер по формуле (1+omega^2*dt^2) гарантированно
//      расходится (план: «энергия растёт... при любом dt»).
//
// Собирается с SPH_DIM=3 (CMakeLists.txt), по образцу test_fiber_bond.cpp --
// тест 3) использует реальную Fiber-Bond физику (interactions::compute_fiber_bond).
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp"
#include "ode/ode_solver.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

static_assert(dim == 3, "test_ode_solver требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b)); }

// ---------------------------------------------------------------------------
// 1) Euler: первый порядок сходимости.
//
// Гармонический осциллятор x'' = -k*x (единичная "масса" -- dv/dt = -k*x
// прописывается вручную в fiber_rhs.dv, минуя Fiber-Bond физику: цель теста --
// проверить ПОРЯДОК СХОДИМОСТИ САМОГО РЕШАТЕЛЯ, а не Fiber-Bond). Аналитика
// при x(0) = x0, v(0) = 0: x(t) = x0*cos(omega*t), omega = sqrt(k).
//
// ode::EulerSolver::step_fiber -- буквальная реализация формулы плана
// (fiber/PLAN.md §3.9): x^{n+1} = x^n + dt*v^n (СТАРАЯ скорость), затем
// v^{n+1} = v^n + dt*f(x^n). T и dt0 подобраны так, чтобы T/dt0,
// T/(dt0/2), T/(dt0/4) были целыми -- одна и та же конечная точка t = T при
// всех трёх шагах, без ошибки округления числа шагов.
// ---------------------------------------------------------------------------
void test_euler_convergence()
{
    std::printf("test_euler_convergence:\n");
    const double k = 4.0;
    const double omega = std::sqrt(k);
    const double x0 = 1.0;
    const double T = 1.0;
    const double dt0 = 0.01; // T/dt0 = 100 -- целое; /2 -> 200; /4 -> 400.

    auto run = [&](double dt) -> double {
        fiber::FiberState fb;
        fb.resize(1);
        for (std::size_t d = 0; d < dim; ++d)
        {
            fb.x[d][0] = 0.0;
            fb.vx[d][0] = 0.0;
        }
        fb.x[0][0] = x0;
        fb.mass[0] = 1.0;
        fb.object_id[0] = 0;
        fb.material_id[0] = 0;

        ode::EulerSolver solver;
        solver.reset(/*nbuf=*/0, /*ntotal=*/0, /*nfiber=*/1);

        const std::size_t nsteps = static_cast<std::size_t>(std::llround(T / dt));
        for (std::size_t s = 0; s < nsteps; ++s)
        {
            solver.fiber_rhs.dv[0][0] = -k * fb.x[0][0]; // f(y^n) = -k*x^n, компоненты 1,2 остаются 0
            solver.step_fiber(fb, 1, dt);
        }
        return fb.x[0][0];
    };

    const double x_exact = x0 * std::cos(omega * T);

    const double e1 = std::fabs(run(dt0) - x_exact);
    const double e2 = std::fabs(run(dt0 / 2.0) - x_exact);
    const double e3 = std::fabs(run(dt0 / 4.0) - x_exact);

    std::printf("    e(dt0)=%.6e  e(dt0/2)=%.6e  e(dt0/4)=%.6e\n", e1, e2, e3);
    check(e1 > 0.0 && e2 > 0.0 && e3 > 0.0, "самопроверка: ошибки ненулевые (тест не тривиален)");
    check(e2 < e1 && e3 < e2, "ошибка монотонно убывает при уменьшении dt");

    const double ratio1 = e1 / e2;
    const double ratio2 = e2 / e3;
    std::printf("    ratio(e1/e2)=%.3f  ratio(e2/e3)=%.3f  (ожидается ~2.0 -- первый порядок)\n", ratio1, ratio2);
    // Первый порядок (не второй): допуск широкий (1.7..2.3), т.к. T не
    // бесконечно мало и в ошибку вносит вклад не только главный член O(dt).
    check(ratio1 > 1.7 && ratio1 < 2.3, "первый порядок сходимости: e(dt0)/e(dt0/2) ~ 2");
    check(ratio2 > 1.7 && ratio2 < 2.3, "первый порядок сходимости: e(dt0/2)/e(dt0/4) ~ 2");
}

// ---------------------------------------------------------------------------
// 3) Нить без контактов: полная механическая энергия ограничена на
// PredictorCorrector.
//
// Короткая прямая цепочка бусин, связанных Fiber-Bond, БЕЗ Fiber-Fiber/
// Fiber-Solid контакта (топология -- единственный источник сил). Начальное
// состояние -- небольшое неравномерное растяжение каждого сегмента (все
// связи в тензии), нулевая начальная скорость -- вся энергия исходно
// потенциальная (fiber/SPH.md, Fiber-Bond: F = EA*max(l/l0-1,0)*x_ij/l,
// поэтому потенциальная энергия связи в тензии PE = 0.5*E*A*(l-l0)^2/l0,
// получена интегрированием F по l от l0; в сжатии PE = 0, сила отсутствует).
//
// dt выбран как малая доля dt_limit связи (fiber/PLAN.md §3.7:
// 0.8/sqrt(k_eff/m_red), k_eff = EA/l0) -- при ЛЮБОМ dt > 0 явный Эйлер для
// этой линейной-по-участкам системы растёт как (1+omega^2*dt^2) за шаг
// (fiber/PLAN.md §3.9); за достаточное число шагов рост становится огромным
// даже при dt много меньше dt_limit. PredictorCorrector должен остаться
// ограниченным на том же dt и том же числе шагов -- в этом и есть проверяемое
// отличие (fiber/PLAN.md §Ф8, раздел G).
// ---------------------------------------------------------------------------
double cross_section_area(double d) { return pi * d * d / 4.0; }

struct ChainFixture
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    double E = 0.0, A = 0.0, l0 = 0.0;
};

ChainFixture make_chain(std::size_t N, double l0, double E, double d, double rho, double mass)
{
    ChainFixture f;
    f.E = E;
    f.A = cross_section_area(d);
    f.l0 = l0;
    fiber::FiberMaterial m;
    m.E = E;
    m.nu = 0.3;
    m.d = d;
    m.rho = rho;
    m.eps_break = 1.0e9; // разрыв не тестируется здесь (Ф7)
    f.materials.push_back(m);

    f.state.resize(N);
    f.topology.resize(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t dcomp = 0; dcomp < dim; ++dcomp)
        {
            f.state.x[dcomp][i] = 0.0;
            f.state.vx[dcomp][i] = 0.0;
        }
        // Неравномерное растяжение каждого сегмента (1%..3%), все связи в
        // тензии -- по образцу tests/test_fiber_bond.cpp::test_momentum_conservation.
        const double stretch = 1.0 + 0.01 * static_cast<double>(i % 3 + 1);
        f.state.x[0][i] = (i == 0) ? 0.0 : (f.state.x[0][i - 1] + l0 * stretch);
        f.state.mass[i] = mass;
        f.state.object_id[i] = 0;
        f.state.material_id[i] = 0;
        if (i + 1 < N)
        {
            f.topology.next[i] = static_cast<long long>(i + 1);
            f.topology.prev[i + 1] = static_cast<long long>(i);
            fiber::FiberBond b;
            b.i = i;
            b.j = i + 1;
            b.l0 = l0;
            b.l_break = l0 * 1.0e9;
            f.topology.bonds.push_back(b);
        }
    }
    return f;
}

double total_kinetic_energy(const fiber::FiberState& state)
{
    double ke = 0.0;
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        double v2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d) v2 += state.vx[d][i] * state.vx[d][i];
        ke += 0.5 * state.mass[i] * v2;
    }
    return ke;
}

double total_bond_potential_energy(const fiber::FiberState& state, const fiber::FiberTopology& topology, double E,
                                    double A, double l0)
{
    double pe = 0.0;
    for (const fiber::FiberBond& b : topology.bonds)
    {
        double l2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double dx = state.x[d][b.j] - state.x[d][b.i];
            l2 += dx * dx;
        }
        const double l = std::sqrt(l2);
        if (l > l0)
        {
            const double stretch = l - l0;
            pe += 0.5 * E * A * stretch * stretch / l0;
        }
    }
    return pe;
}

void test_predictor_corrector_energy_bound()
{
    std::printf("test_predictor_corrector_energy_bound:\n");
    const std::size_t N = 6;
    const double l0 = 0.01, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double A = cross_section_area(d);
    const double m = rho * A * l0;

    ChainFixture fixture = make_chain(N, l0, E, d, rho, m);

    rhs::CpuExecutor exec(2);
    const double dt_limit =
        interactions::compute_fiber_bond_dt_limit(exec, fixture.state, fixture.topology, fixture.materials);
    const double dt = 0.1 * dt_limit; // много меньше dt_limit -- Эйлер тем не менее расходится (см. шапку файла)
    const std::size_t nsteps = 2000;

    const double e0 = total_kinetic_energy(fixture.state) +
                       total_bond_potential_energy(fixture.state, fixture.topology, E, A, l0);
    check(e0 > 0.0, "самопроверка: начальная энергия положительна (связи в тензии)");

    ode::PredictorCorrectorSolver solver;
    solver.reset(/*nbuf=*/0, /*ntotal=*/0, N);
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;

    double e_max = e0;
    bool finite = true;
    for (std::size_t s = 0; s < nsteps; ++s)
    {
        solver.predict_fiber(fixture.state, N, dt);

        std::vector<fiber::FiberDelta> delta =
            interactions::compute_fiber_bond(exec, fixture.state, fixture.topology, fixture.materials, buffers);
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < N; ++i) solver.fiber_rhs.dv[d][i] = delta[i].dv[d];

        solver.correct_fiber(fixture.state, N, dt);
        solver.mark_step_done();

        const double e = total_kinetic_energy(fixture.state) +
                          total_bond_potential_energy(fixture.state, fixture.topology, E, A, l0);
        if (!std::isfinite(e))
        {
            finite = false;
            break;
        }
        if (e > e_max) e_max = e;
    }

    std::printf("    e0 = %.6e, e_max за %zu шагов PredictorCorrector = %.6e (dt = %.3e = 0.1*dt_limit)\n", e0,
                nsteps, e_max, dt);
    check(finite, "энергия остаётся конечной (нет NaN/Inf) на всём прогоне PredictorCorrector");
    // Явный Эйлер на этом dt и числе шагов расходится в >1e100 раз (см.
    // комментарий у fixture выше: growth ~ exp(nsteps*omega^2*dt^2/2), при
    // dt = 0.1*dt_limit и nsteps = 2000 показатель порядка 10). Порог 5.0 --
    // с большим запасом относительно ожидаемого поведения PredictorCorrector
    // (энергия колеблется у e0, не растёт секулярно), но НА ПОРЯДКИ строже,
    // чем неизбежный взрыв Эйлера -- разделяет два решателя надёжно.
    check(e_max < 5.0 * e0, "энергия PredictorCorrector ограничена (не расходится, в отличие от Эйлера)");
}

// ---------------------------------------------------------------------------
// 4) Интерфейс решателя fiber_velocity_update (fiber/fiber_friction_projection.md
//    §4, тест 8.9): advance() даёт fb.vx ПОБИТОВО равную base + h*dv, где
//    base/h взяты из fiber_velocity_update() ДО advance().
//    Сравнение -- ровно тем же выражением `base + h * dv`, что в
//    pc_finalize_velocity_position (`base + coeff * dv[d][i]`) и в
//    euler_step_position_velocity (`vx += dt * dv`, то есть vx + dt*dv); один
//    компилятор и один режим контракции FMA дают одинаковый результат в обоих
//    местах, поэтому проверка именно на равенство (==) битов.
// ---------------------------------------------------------------------------
struct Lcg
{
    unsigned long long s = 0x9E3779B97F4A7C15ULL;
    double next() // равномерно в [-1, 1)
    {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) * (2.0 / 9007199254740992.0) - 1.0;
    }
};

template <class Solver>
void check_velocity_update_step(Solver& solver, fiber::FiberState& fb, std::size_t n, double dt, Lcg& rng,
                                const char* label, double expect_h, bool expect_base_is_fb)
{
    ParticleSystem sys;
    solver.predict(sys, 0, fb, n, dt);
    // RHS шага: новая dv (в боевом цикле -- между predict и advance).
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t i = 0; i < n; ++i) solver.fiber_rhs.dv[d][i] = 10.0 * rng.next();

    const auto upd = solver.fiber_velocity_update(fb, dt);
    std::vector<double> base[dim];
    std::vector<double> dv[dim];
    bool base_ptr_ok = true;
    for (std::size_t d = 0; d < dim; ++d)
    {
        base[d] = *upd.base[d];
        dv[d] = solver.fiber_rhs.dv[d];
        if (expect_base_is_fb && upd.base[d] != &fb.vx[d]) base_ptr_ok = false;
        if (!expect_base_is_fb && upd.base[d] == &fb.vx[d]) base_ptr_ok = false;
    }

    solver.advance(sys, 0, fb, n, dt);

    bool bitwise = true;
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t i = 0; i < n; ++i)
        {
            const double expected = base[d][i] + upd.h * dv[d][i];
            if (!(fb.vx[d][i] == expected)) bitwise = false;
        }
    check(upd.h == expect_h, std::string(label) + ": h == " + (expect_h == dt ? "dt" : "dt/2"));
    check(base_ptr_ok, std::string(label) + ": база = " + (expect_base_is_fb ? "fb.vx" : "v_min_fiber_ (не fb.vx)"));
    check(bitwise, std::string(label) + ": fb.vx побитово == base + h*dv");
}

void test_fiber_velocity_update()
{
    std::printf("test_fiber_velocity_update:\n");
    const std::size_t n = 37;
    const double dt = 0.0137;
    Lcg rng;

    auto make_fb = [&](fiber::FiberState& fb) {
        fb.resize(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            fb.mass[i] = 1.0 + 0.1 * static_cast<double>(i);
            fb.object_id[i] = 0;
            fb.material_id[i] = 0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                fb.x[d][i] = rng.next();
                fb.vx[d][i] = 3.0 * rng.next();
            }
        }
    };

    {
        fiber::FiberState fb;
        make_fb(fb);
        ode::PredictorCorrectorSolver pc;
        pc.reset(0, 0, n);
        check_velocity_update_step(pc, fb, n, dt, rng, "PC шаг 1", dt / 2.0, /*base_is_fb=*/true);
        // шаг n>=2: predict сдвигает fb.vx (u*), база -- сохранённая v^n
        check_velocity_update_step(pc, fb, n, dt * 1.3, rng, "PC шаг 2", dt * 1.3, false);
        check_velocity_update_step(pc, fb, n, dt * 0.7, rng, "PC шаг 3", dt * 0.7, false);
    }
    {
        fiber::FiberState fb;
        make_fb(fb);
        ode::EulerSolver eu;
        eu.reset(0, 0, n);
        check_velocity_update_step(eu, fb, n, dt, rng, "Euler шаг 1", dt, true);
        check_velocity_update_step(eu, fb, n, dt * 1.3, rng, "Euler шаг 2", dt * 1.3, true);
    }
}

// ---------------------------------------------------------------------------
// 5) Интерфейс решателя solid_velocity_update (fiber/fiber_solid_friction_
//    implementation.md §4, аналог 8.9 для Solid): advance() даёт sys.vx
//    ПОБИТОВО равную (base + h*dvx) + av (PC; av -- поправка XSPH,
//    pc_finalize_velocity_position_with_av: `base + coeff*dv + av`) либо
//    base + h*dvx (Эйлер, av == nullptr), где base/h/av взяты из
//    solid_velocity_update() ДО advance(). Порядок сложения именно такой --
//    (base + h*dvx) + av; проверяется равенство (==) битов.
// ---------------------------------------------------------------------------
template <class Solver>
void check_solid_velocity_update_step(Solver& solver, ParticleSystem& sys, fiber::FiberState& fb, std::size_t ns,
                                      std::size_t nf, double dt, double dt_prev, Lcg& rng, const char* label,
                                      double expect_h, bool expect_base_is_sysvx, bool expect_av)
{
    solver.predict(sys, ns, fb, nf, dt_prev);
    // RHS шага: новые dvx и XSPH-поправка av (в боевом цикле -- между predict и advance).
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t i = 0; i < ns; ++i)
        {
            solver.solid_rhs.dvx[d][i] = 10.0 * rng.next();
            solver.solid_rhs.av[d][i] = 0.05 * rng.next();
        }

    const auto upd = solver.solid_velocity_update(sys, dt);
    std::vector<double> base[dim];
    std::vector<double> dvx[dim];
    std::vector<double> avv[dim];
    bool base_ptr_ok = true;
    for (std::size_t d = 0; d < dim; ++d)
    {
        base[d] = *upd.base[d];
        dvx[d] = solver.solid_rhs.dvx[d];
        avv[d] = solver.solid_rhs.av[d];
        const bool is_sysvx = (upd.base[d] == &sys.vx[d]);
        if (expect_base_is_sysvx != is_sysvx) base_ptr_ok = false;
    }
    const bool av_ok = expect_av ? (upd.av == &solver.solid_rhs.av) : (upd.av == nullptr);

    solver.advance(sys, ns, fb, nf, dt);

    bool bitwise = true;
    bool differs_from_no_av = false; // самопроверка: av на самом деле влияет на результат (PC)
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t i = 0; i < ns; ++i)
        {
            double expected = base[d][i] + upd.h * dvx[d][i];
            const double without_av = expected;
            if (upd.av) expected += avv[d][i];
            if (!(sys.vx[d][i] == expected)) bitwise = false;
            if (sys.vx[d][i] != without_av) differs_from_no_av = true;
        }
    check(upd.h == expect_h, std::string(label) + ": h == " + (expect_h == dt ? "dt" : "dt/2"));
    check(base_ptr_ok, std::string(label) + ": база = " + (expect_base_is_sysvx ? "sys.vx" : "v_min_ (не sys.vx)"));
    check(av_ok, std::string(label) + (expect_av ? ": av == &solid_rhs.av (XSPH)" : ": av == nullptr (у Эйлера XSPH нет)"));
    check(bitwise, std::string(label) + ": sys.vx побитово == (base + h*dvx) + av");
    if (expect_av) check(differs_from_no_av, std::string(label) + ": самопроверка: av влияет на результат (тест не тривиален)");
}

void test_solid_velocity_update()
{
    std::printf("test_solid_velocity_update:\n");
    const std::size_t ns = 41, nf = 3;
    const double dt = 0.0137;
    Lcg rng;

    auto make_sys = [&](ParticleSystem& sys) {
        sys.resize(ns);
        for (std::size_t i = 0; i < ns; ++i)
        {
            sys.inzializ_index(i);
            sys.mass[i] = 1.0 + 0.1 * static_cast<double>(i);
            sys.rho[i] = 1000.0;
            sys.rho_0[i] = 1000.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                sys.x[d][i] = rng.next();
                sys.vx[d][i] = 3.0 * rng.next();
            }
        }
    };
    auto make_fb = [&](fiber::FiberState& fb) {
        fb.resize(nf);
        for (std::size_t i = 0; i < nf; ++i)
        {
            fb.mass[i] = 1.0;
            fb.object_id[i] = 0;
            fb.material_id[i] = 0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                fb.x[d][i] = rng.next();
                fb.vx[d][i] = rng.next();
            }
        }
    };

    {
        ParticleSystem sys;
        fiber::FiberState fb;
        make_sys(sys);
        make_fb(fb);
        ode::PredictorCorrectorSolver pc;
        pc.reset(ns, ns, nf);
        check_solid_velocity_update_step(pc, sys, fb, ns, nf, dt, dt, rng, "PC шаг 1", dt / 2.0, /*base=sys.vx*/ true, true);
        // шаг n>=2: predict сдвигает sys.vx (u*), база -- сохранённая v^n
        check_solid_velocity_update_step(pc, sys, fb, ns, nf, dt * 1.3, dt, rng, "PC шаг 2", dt * 1.3, false, true);
        check_solid_velocity_update_step(pc, sys, fb, ns, nf, dt * 0.7, dt * 1.3, rng, "PC шаг 3", dt * 0.7, false, true);
    }
    {
        ParticleSystem sys;
        fiber::FiberState fb;
        make_sys(sys);
        make_fb(fb);
        ode::EulerSolver eu;
        eu.reset(ns, ns, nf);
        check_solid_velocity_update_step(eu, sys, fb, ns, nf, dt, dt, rng, "Euler шаг 1", dt, true, false);
        check_solid_velocity_update_step(eu, sys, fb, ns, nf, dt * 1.3, dt, rng, "Euler шаг 2", dt * 1.3, true, false);
    }
}

// ---------------------------------------------------------------------------
// 6) xsph_mode = position: скорость = base + h*dvx БЕЗ av (побитово), позиция
//    сдвинута на dt*(v_new + av), solid_velocity_update().av == nullptr.
// ---------------------------------------------------------------------------
void test_xsph_position()
{
    std::printf("test_xsph_position:\n");
    const std::size_t ns = 41, nf = 3;
    const double dt = 0.0137;
    Lcg rng;

    ParticleSystem sys;
    sys.resize(ns);
    for (std::size_t i = 0; i < ns; ++i)
    {
        sys.inzializ_index(i);
        sys.mass[i] = 1.0 + 0.1 * static_cast<double>(i);
        sys.rho[i] = 1000.0;
        sys.rho_0[i] = 1000.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            sys.x[d][i] = rng.next();
            sys.vx[d][i] = 3.0 * rng.next();
        }
    }
    fiber::FiberState fb;
    fb.resize(nf);

    ode::PredictorCorrectorSolver pc;
    pc.xsph_position = true;
    pc.reset(ns, ns, nf);

    const double dts[3] = {dt, dt * 1.3, dt * 0.7};
    double dt_prev = dt;
    for (int step = 0; step < 3; ++step)
    {
        const double h = dts[step];
        pc.predict(sys, ns, fb, nf, dt_prev);
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < ns; ++i)
            {
                pc.solid_rhs.dvx[d][i] = 10.0 * rng.next();
                pc.solid_rhs.av[d][i] = 0.05 * rng.next();
            }
        const auto upd = pc.solid_velocity_update(sys, h);
        std::vector<double> base[dim], x0[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            base[d] = *upd.base[d];
            x0[d] = sys.x[d];
        }
        const double hh = upd.h;
        const bool av_null = (upd.av == nullptr);
        pc.advance(sys, ns, fb, nf, h);

        bool v_ok = true, x_ok = true;
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < ns; ++i)
            {
                const double v_new = base[d][i] + hh * pc.solid_rhs.dvx[d][i];
                if (!(sys.vx[d][i] == v_new)) v_ok = false;
                if (!(sys.x[d][i] == x0[d][i] + h * (v_new + pc.solid_rhs.av[d][i]))) x_ok = false;
            }
        const std::string lbl = "position шаг " + std::to_string(step + 1);
        check(av_null, lbl + ": upd.av == nullptr");
        check(v_ok, lbl + ": sys.vx побитово == base + h*dvx (без av)");
        check(x_ok, lbl + ": x побитово == x0 + dt*(v_new + av)");
        dt_prev = h;
    }
}

} // namespace

int main()
{
    test_xsph_position();
    test_euler_convergence();
    test_predictor_corrector_energy_bound();
    test_fiber_velocity_update();
    test_solid_velocity_update();

    if (failures == 0)
        std::printf("test_ode_solver: OK\n");
    else
        std::printf("test_ode_solver: %d failure(s)\n", failures);

    return failures;
}
