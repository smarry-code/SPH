// tests/test_fiber_friction.cpp
// Unit-тесты стадии трения нить-нить -- проекции на уровне скоростей
// (interactions/fiber_friction.hpp; fiber/fiber_friction_projection.md, §8).
//
// Тесты работают со стадией НАПРЯМУЮ: задаются FiberState (массы), слоты
// FiberFrictionContact, dv_no_friction и ode::FiberVelocityUpdate{base, h};
// вызываются compute_fiber_friction_contributions + rhs::reduce_contributions;
// итоговая скорость v^{n+1} = base + h*(a + a^f) считается той же формулой,
// что в решателе (ode/ode_solver.hpp, §4 спецификации).
//
// Разделы соответствуют §8 спецификации:
//   8.1 прилипание; 8.2 скольжение; 8.3 статическое трение под другими силами;
//   8.4 касательность F_t и чисто нормальная скорость; 8.5 импульс и
//   распределение; 8.6 диссипация на 1000 случайных контактов; 8.7 независимость
//   прилипания от dt; 8.8 mu=0; плюс диагностика FiberFrictionDiagnostics.
//   (8.9 -- интерфейс решателя -- в tests/test_ode_solver.cpp.)
//
// Вторая часть файла (fiber/fiber_solid_friction_implementation.md, задача 4,
// §7.1) -- трение НИТЬ-SOLID: единая стадия compute_friction_contributions по
// слотам FiberSolidFrictionContact и ode::SolidVelocityUpdate. Сцена SolidScene
// (нить + частицы Solid), применение стадии apply_solid_stage; разделы:
// test_fs_stick, test_fs_slip, test_fs_static_under_force, test_fs_xsph,
// test_fs_signs (таблица §3), test_fs_tangential,
// test_fs_momentum_distribution, test_fs_dissipation, test_fs_zero_mu,
// test_fs_kinematic, test_fs_counts_and_diagnostics, test_fs_mass_splitting
// (три частицы Solid под сегментом; реальная плита из kernel'а; бусина в двух
// видах контактов; одиночный контакт == формула задачи 3),
// test_fs_solver_interface (реальные PredictorCorrector/Euler + стадия).
// Собирается с SPH_DIM=3 (CMakeLists.txt).
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_friction.hpp"
#include "interactions/fiber_normal_projection.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "materials.hpp"
#include "ode/ode_solver.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_friction требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

double dot3(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm3(const double* a) { return std::sqrt(dot3(a, a)); }

// Детерминированный генератор (splitmix64) -- не зависит от стандартной
// библиотеки, результат одинаков на всех платформах.
struct Rng
{
    std::uint64_t s;
    explicit Rng(std::uint64_t seed) : s(seed) {}
    std::uint64_t next()
    {
        std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uni() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); } // [0,1)
    double range(double a, double b) { return a + (b - a) * uni(); }
};

// -----------------------------------------------------------------------
// Сцена: n бусин, base, dv_no_friction, слоты; upd указывает на base сцены
// (сцену НЕ копировать после init_upd).
// -----------------------------------------------------------------------
struct Scene
{
    fiber::FiberState state;
    std::vector<double> base[dim];
    std::vector<double> dv[dim];
    std::vector<interactions::FiberFrictionContact> slots;
    ode::FiberVelocityUpdate upd{};

    void init(std::size_t n, double h)
    {
        state.resize(n);
        for (std::size_t d = 0; d < dim; ++d)
        {
            base[d].assign(n, 0.0);
            dv[d].assign(n, 0.0);
            upd.base[d] = &base[d];
        }
        upd.h = h;
    }
};

// Результат применения стадии.
struct StageResult
{
    std::vector<fiber::FiberDelta> delta; // a^f по бусинам
    std::vector<double> vtrial[dim];      // base + h*dv (без трения)
    std::vector<double> vnew[dim];        // base + h*(dv + a^f)
};

// c_k (mass splitting) сцены -- той же функцией, что в стадии.
std::vector<double> bead_counts(const rhs::CpuExecutor& exec, const Scene& sc)
{
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> cbuf;
    std::vector<double> counts;
    interactions::compute_fiber_friction_bead_counts(exec, sc.state.size(), sc.slots, cbuf, counts);
    return counts;
}

StageResult apply_stage(const rhs::CpuExecutor& exec, const Scene& sc)
{
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
    const std::vector<double> counts = bead_counts(exec, sc);
    interactions::compute_fiber_friction_contributions(exec, sc.state, sc.upd, sc.dv, sc.slots, counts, buffers);
    StageResult r;
    r.delta = rhs::reduce_contributions(exec, buffers, sc.state.size());
    const std::size_t n = sc.state.size();
    for (std::size_t d = 0; d < dim; ++d)
    {
        r.vtrial[d].resize(n);
        r.vnew[d].resize(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            r.vtrial[d][i] = (*sc.upd.base[d])[i] + sc.upd.h * sc.dv[d][i];
            r.vnew[d][i] = (*sc.upd.base[d])[i] + sc.upd.h * (sc.dv[d][i] + r.delta[i].dv[d]);
        }
    }
    return r;
}

// Касательная относительная скорость точки контакта по массивам скоростей v.
// Возвращает |u_t|; ut -- вектор.
double contact_ut(const std::vector<double> (&v)[dim], const interactions::FiberFrictionContact& c, double ut[3])
{
    double u[3];
    for (std::size_t d = 0; d < dim; ++d)
    {
        double vA = v[d][c.a0];
        if (c.a1 != -1) vA = (1.0 - c.s) * vA + c.s * v[d][static_cast<std::size_t>(c.a1)];
        double vB = v[d][c.b0];
        if (c.b1 != -1) vB = (1.0 - c.t) * vB + c.t * v[d][static_cast<std::size_t>(c.b1)];
        u[d] = vA - vB;
    }
    const double un = dot3(u, c.n);
    for (std::size_t d = 0; d < dim; ++d) ut[d] = u[d] - un * c.n[d];
    return norm3(ut);
}

// Кинетическая энергия первых 4 бусин.
double kinetic(const Scene& sc, const std::vector<double> (&v)[dim])
{
    double ke = 0.0;
    for (std::size_t i = 0; i < sc.state.size(); ++i)
        ke += 0.5 * sc.state.mass[i] * (v[0][i] * v[0][i] + v[1][i] * v[1][i] + v[2][i] * v[2][i]);
    return ke;
}

void unit(double v[3])
{
    const double L = norm3(v);
    for (int k = 0; k < 3; ++k) v[k] /= L;
}

// Одиночный контакт: бусины 0,1 -- сторона A, 2,3 -- сторона B.
void set_single_slot(Scene& sc, double s, double t, const double n[3], double N, double mu)
{
    sc.slots.assign(1, interactions::FiberFrictionContact{});
    auto& c = sc.slots[0];
    c.active = true;
    c.a0 = 0;
    c.a1 = 1;
    c.b0 = 2;
    c.b1 = 3;
    c.s = s;
    c.t = t;
    for (int k = 0; k < 3; ++k) c.n[k] = n[k];
    c.N = N;
    c.mu = mu;
}

// Сцена 8.1: s=0.3, t=0.7, разные массы, ненулевые базовые скорости и dv.
void build_scene_8_1(Scene& sc, double h, double N, double mu)
{
    sc.init(4, h);
    const double m[4] = {1.0, 3.0, 0.5, 2.0};
    for (int i = 0; i < 4; ++i) sc.state.mass[i] = m[i];
    const double v[4][3] = {{1.0, 0.2, -0.4}, {0.5, -0.3, 0.7}, {-0.2, 0.9, 0.1}, {0.4, 0.1, -0.6}};
    const double a[4][3] = {{5.0, -2.0, 1.0}, {-1.0, 3.0, 0.5}, {2.0, 2.0, -4.0}, {0.3, -0.7, 1.1}};
    for (int i = 0; i < 4; ++i)
        for (int d = 0; d < 3; ++d)
        {
            sc.base[d][i] = v[i][d];
            sc.dv[d][i] = a[i][d];
        }
    double n[3] = {0.3, -0.5, 0.8};
    unit(n);
    set_single_slot(sc, 0.3, 0.7, n, N, mu);
}

double slot_W(const Scene& sc, const interactions::FiberFrictionContact& c)
{
    return interactions::fiber_contact_side_inv_mass(sc.state, c.a0, c.a1, c.s) +
           interactions::fiber_contact_side_inv_mass(sc.state, c.b0, c.b1, c.t);
}

// ----------------------------------------------------------------------
void test_8_1_stick(double& err_out)
{
    std::printf("test_8_1_stick:\n");
    rhs::CpuExecutor exec(2);
    Scene sc;
    build_scene_8_1(sc, 1.0e-3, 2.0, 1000.0); // muN = 2000 -- заведомо достаточно
    const StageResult r = apply_stage(exec, sc);

    double ut0[3], ut1[3];
    const double n0 = contact_ut(r.vtrial, sc.slots[0], ut0);
    const double n1 = contact_ut(r.vnew, sc.slots[0], ut1);
    err_out = n1 / n0;
    std::printf("    |u_t| до трения %.6e, после %.3e, отношение %.3e\n", n0, n1, n1 / n0);
    check(n0 > 1e-3, "исходная касательная скорость ненулевая");
    check(n1 <= 1e-12 * n0, "после обновления |u_t| <= 1e-12*|u_t| (точное прилипание)");

    interactions::FiberFrictionImpulse imp;
    interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, sc.slots[0], imp);
    const double Ft = norm3(imp.Ft);
    const double muN = sc.slots[0].mu * sc.slots[0].N;
    check(imp.stick, "выбрана ветка прилипания");
    check(Ft < muN, "|F_t| < mu*N");
    // Ft*h*W == |u_t|
    check(std::fabs(Ft * sc.upd.h * imp.W - n0) <= 1e-12 * n0, "|F_t|*h*W == |u_t| (остановка ровно)");
}

// ----------------------------------------------------------------------
void test_8_2_slip(double& ratio_err_out)
{
    std::printf("test_8_2_slip:\n");
    rhs::CpuExecutor exec(2);
    Scene sc;
    build_scene_8_1(sc, 1.0e-3, 2.0, 0.25); // muN = 0.5, muN*h*W << |u_t|
    const StageResult r = apply_stage(exec, sc);

    double ut0[3], ut1[3];
    const double n0 = contact_ut(r.vtrial, sc.slots[0], ut0);
    const double n1 = contact_ut(r.vnew, sc.slots[0], ut1);
    const double muN = sc.slots[0].mu * sc.slots[0].N;
    const double W = slot_W(sc, sc.slots[0]);
    const double h = sc.upd.h;
    std::printf("    |u_t| до %.6e, после %.6e; muN*h*W = %.6e\n", n0, n1, muN * h * W);

    interactions::FiberFrictionImpulse imp;
    interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, sc.slots[0], imp);
    check(!imp.stick, "выбрана ветка скольжения");
    const double Ft = norm3(imp.Ft);
    ratio_err_out = std::fabs(Ft / muN - 1.0);
    std::printf("    |F_t|/muN - 1 = %.3e\n", Ft / muN - 1.0);
    check(ratio_err_out <= 1e-14, "|F_t| == mu*N (отн. 1e-14)");

    // |F_t| по реально применённому ускорению бусины a0: a^f = (1-s)*F_t/m
    const auto& c = sc.slots[0];
    double fa0[3];
    for (int d = 0; d < 3; ++d) fa0[d] = r.delta[c.a0].dv[d] * sc.state.mass[c.a0] / (1.0 - c.s);
    check(std::fabs(norm3(fa0) / muN - 1.0) <= 1e-14, "|F_t| == mu*N по вкладу бусины a0");

    // направление -u_t/|u_t|
    double uhat[3] = {ut0[0] / n0, ut0[1] / n0, ut0[2] / n0};
    check(std::fabs(dot3(fa0, uhat) / (-norm3(fa0)) - 1.0) <= 1e-13, "направление F_t == -u_t/|u_t|");

    // Δu_t == -muN*h*W*û_t
    double dut[3], expect[3];
    for (int d = 0; d < 3; ++d)
    {
        dut[d] = ut1[d] - ut0[d];
        expect[d] = -muN * h * W * uhat[d];
    }
    double diff[3] = {dut[0] - expect[0], dut[1] - expect[1], dut[2] - expect[2]};
    check(norm3(diff) <= 1e-10 * muN * h * W, "касательная скорость уменьшилась ровно на muN*h*W");
    check(dot3(ut1, uhat) > 0.0, "касательная относительная скорость не сменила знак");
    check(std::fabs(n0 - n1 - muN * h * W) <= 1e-10 * muN * h * W, "|u_t| уменьшился ровно на muN*h*W");
}

// ----------------------------------------------------------------------
void test_8_3_static_under_force(double& stick_err_out)
{
    std::printf("test_8_3_static_under_force:\n");
    rhs::CpuExecutor exec(2);
    const double h = 1.0e-3, N = 2.5, mu = 0.4;
    const double muN = mu * N; // 1.0

    for (int variant = 0; variant < 2; ++variant)
    {
        Scene sc;
        build_scene_8_1(sc, h, N, mu);
        // Нулевая базовая скорость и нулевые dv, кроме заданного относительного ускорения g.
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < 4; ++i)
            {
                sc.base[d][i] = 0.0;
                sc.dv[d][i] = 0.0;
            }
        const auto& c = sc.slots[0];
        const double W = slot_W(sc, c);
        const double m_eff = 1.0 / W;
        // касательный единичный вектор t1 ⟂ n
        double tt[3] = {c.n[1] * 1.0 - c.n[2] * 0.0, c.n[2] * 1.0 - c.n[0] * 0.0 - 0.0, c.n[0] * 0.0 - c.n[1] * 1.0};
        // cross(n, x̂) = (0, n_z, -n_y)
        tt[0] = 0.0;
        tt[1] = c.n[2];
        tt[2] = -c.n[1];
        unit(tt);
        const double frac = (variant == 0) ? 0.5 : 2.0; // m_eff|g_t| = frac*muN
        const double gt = frac * muN / m_eff;
        // g = gt*t1 + 3*n; ускорение прикладываем к обеим бусинам стороны A (точка контакта
        // получает ровно g, независимо от s), сторона B неподвижна.
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double g = gt * tt[d] + 3.0 * c.n[d];
            sc.dv[d][0] = g;
            sc.dv[d][1] = g;
        }
        const StageResult r = apply_stage(exec, sc);

        double ut0[3], ut1[3];
        const double n0 = contact_ut(r.vtrial, sc.slots[0], ut0);
        const double n1 = contact_ut(r.vnew, sc.slots[0], ut1);
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, sc.slots[0], imp);
        const double Ft = norm3(imp.Ft);
        std::printf("    m_eff|g_t|/muN = %.2f: |u_t| до %.6e, после %.6e, |F_t|/muN = %.15f\n", frac, n0, n1, Ft / muN);
        check(std::fabs(n0 - h * gt) <= 1e-12 * h * gt, "пробная касательная скорость == h*g_t");
        if (variant == 0)
        {
            stick_err_out = n1 / n0;
            check(n1 <= 1e-12 * n0, "m_eff|g| < muN: относительная касательная скорость после шага 0");
            check(imp.stick && Ft < muN, "прилипание, |F_t| < muN");
            check(std::fabs(Ft - m_eff * gt) <= 1e-12 * muN, "F_t уравновешивает m_eff*g_t");
        }
        else
        {
            check(!imp.stick, "m_eff|g| > muN: скольжение");
            check(std::fabs(Ft / muN - 1.0) <= 1e-14, "|F_t| == mu*N");
            const double expected = h * (gt - W * muN);
            check(std::fabs(n1 - expected) <= 1e-10 * expected, "остаточная скорость h*(g_t - W*muN)");
            double gh[3] = {tt[0], tt[1], tt[2]};
            check(dot3(ut1, gh) > 0.0, "остаточная скорость направлена вдоль g_t");
        }
    }
}

// ----------------------------------------------------------------------
void test_8_4_tangential(void)
{
    std::printf("test_8_4_tangential:\n");
    rhs::CpuExecutor exec(2);

    // F_t . n == 0 для наклонной нормали, оба режима.
    for (int mode = 0; mode < 2; ++mode)
    {
        Scene sc;
        build_scene_8_1(sc, 1.0e-3, 2.0, mode == 0 ? 1000.0 : 0.25);
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, sc.slots[0], imp);
        const double Ft = norm3(imp.Ft);
        check(Ft > 0.0 && std::fabs(dot3(imp.Ft, sc.slots[0].n)) <= 1e-14 * Ft,
              mode == 0 ? "F_t . n == 0 (прилипание)" : "F_t . n == 0 (скольжение)");
        // F_t . u_t <= 0
        const StageResult r = apply_stage(exec, sc);
        double ut0[3];
        contact_ut(r.vtrial, sc.slots[0], ut0);
        check(dot3(imp.Ft, ut0) <= 0.0, "F_t . u_t <= 0");
    }

    // u чисто нормальная (нормаль вдоль оси -> u_t точно 0): F_t == 0, dv не меняется.
    {
        Scene sc;
        build_scene_8_1(sc, 1.0e-3, 2.0, 1000.0);
        const double n[3] = {0.0, 0.0, 1.0};
        set_single_slot(sc, 0.3, 0.7, n, 2.0, 1000.0);
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < 4; ++i)
            {
                sc.base[d][i] = 0.0;
                sc.dv[d][i] = 0.0;
            }
        for (std::size_t i = 0; i < 2; ++i) sc.base[2][i] = 1.7; // A летит вдоль n, B неподвижна
        const StageResult r = apply_stage(exec, sc);
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, sc.slots[0], imp);
        bool all_zero = true;
        for (const auto& dl : r.delta)
            for (int d = 0; d < 3; ++d)
                if (dl.dv[d] != 0.0) all_zero = false;
        check(imp.ut_norm == 0.0 && !imp.has_force, "u чисто нормальная: u_t == 0, силы нет");
        check(norm3(imp.Ft) == 0.0 && all_zero, "u чисто нормальная: F_t == 0, вклады нулевые");
    }

    // Наклонная нормаль: u = c*n с точностью до округления -> F_t пренебрежимо мала.
    {
        Scene sc;
        build_scene_8_1(sc, 1.0e-3, 2.0, 1000.0);
        const auto& c = sc.slots[0];
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < 4; ++i)
            {
                sc.base[d][i] = 0.0;
                sc.dv[d][i] = 0.0;
            }
        for (std::size_t d = 0; d < dim; ++d)
        {
            sc.base[d][0] = 1.7 * c.n[d];
            sc.base[d][1] = 1.7 * c.n[d];
        }
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, sc.slots[0], imp);
        const double h = sc.upd.h;
        check(imp.ut_norm <= 1e-15 * 1.7 * 2.0, "наклонная нормаль, u ~ нормальная: |u_t| на уровне округления");
        check(norm3(imp.Ft) * h * imp.W <= 1e-15 * 1.7 * 2.0, "F_t*h*W на уровне округления (нет заметной силы)");
    }
}

// ----------------------------------------------------------------------
// Общая проверка 8.5 для одной конфигурации. Возвращает false при нарушении.
struct Check85
{
    double mom_err = 0.0;   // max |sum m a^f| / sum |m a^f|
    double dist_err = 0.0;  // max относительное отклонение от ±вес*F_t/m
    bool zero_weight_ok = true;
};

void check_85_config(Check85& acc, const rhs::CpuExecutor& exec, Rng& rng, double s, double t, bool a_point, bool b_point)
{
    Scene sc;
    sc.init(4, rng.range(1e-4, 1e-2));
    for (int i = 0; i < 4; ++i) sc.state.mass[i] = rng.range(0.1, 10.0);
    for (std::size_t d = 0; d < dim; ++d)
        for (int i = 0; i < 4; ++i)
        {
            sc.base[d][i] = rng.range(-1.0, 1.0);
            sc.dv[d][i] = rng.range(-5.0, 5.0);
        }
    double n[3] = {rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)};
    unit(n);
    set_single_slot(sc, a_point ? 0.0 : s, b_point ? 0.0 : t, n, rng.range(0.5, 5.0), rng.range(0.05, 1.0));
    if (a_point) sc.slots[0].a1 = -1;
    if (b_point) sc.slots[0].b1 = -1;
    const auto& c = sc.slots[0];

    const StageResult r = apply_stage(exec, sc);
    interactions::FiberFrictionImpulse imp;
    interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, c, imp);

    // импульс
    double sum[3] = {0, 0, 0};
    double sabs = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int d = 0; d < 3; ++d)
        {
            const double f = sc.state.mass[i] * r.delta[i].dv[d];
            sum[d] += f;
            sabs += std::fabs(f);
        }
    if (sabs > 0.0) acc.mom_err = std::max(acc.mom_err, norm3(sum) / sabs);

    // ожидаемое распределение
    double exp_a[4][3] = {};
    const double wa0 = (c.a1 == -1) ? 1.0 : 1.0 - c.s;
    const double wa1 = (c.a1 == -1) ? 0.0 : c.s;
    const double wb0 = (c.b1 == -1) ? 1.0 : 1.0 - c.t;
    const double wb1 = (c.b1 == -1) ? 0.0 : c.t;
    for (int d = 0; d < 3; ++d)
    {
        exp_a[0][d] = wa0 * imp.Ft[d] / sc.state.mass[0];
        exp_a[1][d] = wa1 * imp.Ft[d] / sc.state.mass[1];
        exp_a[2][d] = -wb0 * imp.Ft[d] / sc.state.mass[2];
        exp_a[3][d] = -wb1 * imp.Ft[d] / sc.state.mass[3];
    }
    const double Ftn = norm3(imp.Ft);
    for (int i = 0; i < 4; ++i)
        for (int d = 0; d < 3; ++d)
        {
            const double scale = Ftn / sc.state.mass[i];
            if (scale > 0.0) acc.dist_err = std::max(acc.dist_err, std::fabs(r.delta[i].dv[d] - exp_a[i][d]) / scale);
        }
    // бусины с нулевым весом (s==0 -> a1, s==1 -> a0, ...) или отсутствующие -- ровно 0
    auto is_zero = [&](int i) { return r.delta[i].dv[0] == 0.0 && r.delta[i].dv[1] == 0.0 && r.delta[i].dv[2] == 0.0; };
    if ((c.a1 == -1 || c.s == 0.0) && !is_zero(1)) acc.zero_weight_ok = false;
    if (c.s == 1.0 && c.a1 != -1 && !is_zero(0)) acc.zero_weight_ok = false;
    if ((c.b1 == -1 || c.t == 0.0) && !is_zero(3)) acc.zero_weight_ok = false;
    if (c.t == 1.0 && c.b1 != -1 && !is_zero(2)) acc.zero_weight_ok = false;
}

void test_8_5_momentum_distribution()
{
    std::printf("test_8_5_momentum_distribution:\n");
    rhs::CpuExecutor exec(2);
    Rng rng(20260930);

    Check85 gen;
    for (int k = 0; k < 500; ++k) check_85_config(gen, exec, rng, rng.uni(), rng.uni(), false, false);
    std::printf("    произвольные s,t (500): max|sum m a|/sum|m a| = %.3e, max откл. распределения = %.3e\n",
                gen.mom_err, gen.dist_err);
    check(gen.mom_err <= 1e-14, "sum m_k a^f_k == 0 при произвольных s,t и массах");
    check(gen.dist_err <= 1e-14, "вклад бусины == ±вес*F_t/m_k (веса 1-s, s, 1-t, t; B со знаком -)");

    // Точечные сегменты: вся сила на a0 (и/или b0).
    Check85 pt;
    for (int k = 0; k < 100; ++k)
    {
        check_85_config(pt, exec, rng, 0.0, rng.uni(), true, false);
        check_85_config(pt, exec, rng, rng.uni(), 0.0, false, true);
        check_85_config(pt, exec, rng, 0.0, 0.0, true, true);
    }
    check(pt.mom_err <= 1e-14, "точечный сегмент: сумма m a^f == 0");
    check(pt.dist_err <= 1e-14, "точечный сегмент: вся сила на a0 (b0) с весом 1");
    check(pt.zero_weight_ok, "точечный сегмент: несуществующая бусина a1 (b1) не получает вклада");

    // Нулевые веса: s=0, s=1, t=0, t=1 (в разных комбинациях).
    Check85 zw;
    const double edge[2] = {0.0, 1.0};
    for (double s : edge)
        for (double t : edge)
            for (int k = 0; k < 20; ++k) check_85_config(zw, exec, rng, s, t, false, false);
    for (int k = 0; k < 20; ++k)
    {
        check_85_config(zw, exec, rng, 0.0, rng.uni(), false, false);
        check_85_config(zw, exec, rng, 1.0, rng.uni(), false, false);
        check_85_config(zw, exec, rng, rng.uni(), 0.0, false, false);
        check_85_config(zw, exec, rng, rng.uni(), 1.0, false, false);
    }
    check(zw.mom_err <= 1e-14, "s,t на границах: сумма m a^f == 0");
    check(zw.dist_err <= 1e-14, "s,t на границах: распределение верное");
    check(zw.zero_weight_ok, "бусины с нулевым весом (s=0/1, t=0/1) не получают вклада (ровно 0)");
}

// ----------------------------------------------------------------------
void test_8_6_dissipation(double& max_growth_out)
{
    std::printf("test_8_6_dissipation:\n");
    rhs::CpuExecutor exec(2);
    Rng rng(8600001);

    int bad = 0, n_stick = 0, n_slip = 0, n_nofrc = 0;
    double max_growth = 0.0; // max (KE_after - KE_trial)/KE_trial (отрицательное -- убыль)
    double max_drop = 0.0;
    for (int k = 0; k < 1000; ++k)
    {
        Scene sc;
        sc.init(4, rng.range(1e-3, 1e-1));
        for (int i = 0; i < 4; ++i) sc.state.mass[i] = std::pow(10.0, rng.range(-1.0, 1.0));
        for (std::size_t d = 0; d < dim; ++d)
            for (int i = 0; i < 4; ++i)
            {
                sc.base[d][i] = rng.range(-1.0, 1.0);
                sc.dv[d][i] = rng.range(-5.0, 5.0);
            }
        double n[3] = {rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)};
        unit(n);
        const double muN = std::pow(10.0, rng.range(-3.0, 1.5));
        set_single_slot(sc, rng.uni(), rng.uni(), n, 1.0, muN);

        const StageResult r = apply_stage(exec, sc);
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, sc.slots[0], imp);
        if (!imp.has_force) ++n_nofrc;
        else if (imp.stick) ++n_stick;
        else ++n_slip;

        const double ke0 = kinetic(sc, r.vtrial);
        const double ke1 = kinetic(sc, r.vnew);
        const double growth = (ke1 - ke0) / ke0;
        max_growth = (k == 0) ? growth : std::max(max_growth, growth);
        max_drop = std::min(max_drop, growth);
        if (!(ke1 <= ke0 * (1.0 + 1e-12))) ++bad;
    }
    max_growth_out = max_growth;
    std::printf("    stick %d, slip %d, без силы %d; max (KE_после-KE_проб)/KE_проб = %.3e (в ед. eps=%.3e), min = %.3e\n",
                n_stick, n_slip, n_nofrc, max_growth, max_growth / 2.220446049250313e-16, max_drop);
    check(n_stick > 100 && n_slip > 100, "выборка содержит и прилипание, и скольжение");
    check(bad == 0, "KE с трением <= KE без трения (отн. допуск 1e-12) во всех 1000 контактах");
}

// ----------------------------------------------------------------------
void test_8_7_dt_independence(double& err_h_out, double& err_h10_out)
{
    std::printf("test_8_7_dt_independence:\n");
    rhs::CpuExecutor exec(2);
    const double hs[2] = {1.0e-3, 1.0e-4};
    double errs[2];
    for (int k = 0; k < 2; ++k)
    {
        Scene sc;
        build_scene_8_1(sc, hs[k], 2.0, 1.0e6); // muN*h >> |u_t|/W в обоих случаях
        const StageResult r = apply_stage(exec, sc);
        double ut0[3], ut1[3];
        const double n0 = contact_ut(r.vtrial, sc.slots[0], ut0);
        const double n1 = contact_ut(r.vnew, sc.slots[0], ut1);
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, sc.slots[0], imp);
        errs[k] = n1 / n0;
        std::printf("    h = %.1e: |u_t| до %.6e, после %.3e (отн. %.3e)\n", hs[k], n0, n1, n1 / n0);
        check(imp.stick, "прилипание");
        check(n1 <= 1e-12 * n0, "итоговая касательная относительная скорость == 0");
    }
    err_h_out = errs[0];
    err_h10_out = errs[1];
}

// ----------------------------------------------------------------------
// Сцена с контактом двух крестом лежащих нитей (2 бусины на нить, проникание
// D=0.9d): через compute_fiber_contact_contributions со слотами.
void test_8_8_zero_mu()
{
    std::printf("test_8_8_zero_mu:\n");
    rhs::CpuExecutor exec(2);

    const double d = 1.0e-3;
    fiber::FiberState state;
    fiber::FiberTopology topology;
    state.resize(4);
    topology.resize(4);
    const double X[4][3] = {{-0.005, 0.0, 0.0}, {0.005, 0.0, 0.0}, {0.0, -0.005, 0.9 * d}, {0.0, 0.005, 0.9 * d}};
    for (int i = 0; i < 4; ++i)
    {
        for (int k = 0; k < 3; ++k)
        {
            state.x[k][i] = X[i][k];
            state.vx[k][i] = 0.0;
        }
        state.mass[i] = 1.0e-3 * (i + 1);
        state.object_id[i] = (i < 2) ? 0 : 1;
        state.material_id[i] = 0;
    }
    topology.next[0] = 1;
    topology.prev[1] = 0;
    topology.next[2] = 3;
    topology.prev[3] = 2;
    for (int b = 0; b < 2; ++b)
    {
        fiber::FiberBond bond;
        bond.i = 2 * b;
        bond.j = 2 * b + 1;
        bond.l0 = 0.01;
        bond.l_break = 1e9;
        topology.bonds.push_back(bond);
    }

    for (int variant = 0; variant < 2; ++variant)
    {
        std::vector<fiber::FiberMaterial> materials(1);
        materials[0].E = 2.0e9;
        materials[0].nu = 0.3;
        materials[0].d = d;
        materials[0].rho = 1200.0;
        materials[0].eps_break = 1.0;
        materials[0].mu_ff = (variant == 0) ? 0.0 : 0.5;

        interactions::FiberContactSearchResult search;
        interactions::build_fiber_contact_search(exec, state, topology, materials, search);
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> nb;
        std::vector<interactions::FiberFrictionContact> slots;
        const bool has_active =
            interactions::compute_fiber_contact_contributions(exec, state, materials, search.primitives, nb, slots);

        if (variant == 1)
        {
            check(has_active, "контрольный прогон mu=0.5: has_active == true");
            bool any = false;
            for (const auto& sl : slots) any = any || sl.active;
            check(any, "контрольный прогон mu=0.5: слот активен");
            continue;
        }

        check(!search.primitives.empty(), "в сцене есть примитив контакта");
        check(!has_active, "mu_ff=0: compute_fiber_contact_contributions возвращает has_active == false");
        bool any_active = false;
        for (const auto& sl : slots) any_active = any_active || sl.active;
        check(!any_active, "mu_ff=0: все слоты неактивны");
        check(slots.size() == search.primitives.size(), "слотов по одному на primitive");

        // Стадия на неактивных слотах: ничего не эмитирует.
        Scene sc;
        sc.init(4, 1.0e-4);
        for (int i = 0; i < 4; ++i) sc.state.mass[i] = state.mass[i];
        Rng rng(88);
        for (std::size_t dd = 0; dd < dim; ++dd)
            for (int i = 0; i < 4; ++i)
            {
                sc.base[dd][i] = rng.range(-1.0, 1.0);
                sc.dv[dd][i] = rng.range(-5.0, 5.0);
            }
        sc.slots = slots;
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
        const std::vector<double> counts0 = bead_counts(exec, sc);
        interactions::compute_fiber_friction_contributions(exec, sc.state, sc.upd, sc.dv, sc.slots, counts0, buffers);
        std::size_t emitted = 0;
        for (const auto& b : buffers) emitted += b.size();
        check(emitted == 0, "стадия на неактивных слотах ничего не эмитирует");

        std::vector<double> dv_before[dim];
        for (std::size_t dd = 0; dd < dim; ++dd) dv_before[dd] = sc.dv[dd];
        const auto delta = rhs::reduce_contributions(exec, buffers, sc.state.size());
        for (std::size_t dd = 0; dd < dim; ++dd)
            for (std::size_t i = 0; i < 4; ++i) sc.dv[dd][i] += delta[i].dv[dd];
        bool bitwise = true;
        for (std::size_t dd = 0; dd < dim; ++dd)
            bitwise = bitwise && std::memcmp(dv_before[dd].data(), sc.dv[dd].data(), 4 * sizeof(double)) == 0;
        check(bitwise, "dv после стадии побитово равен dv до неё");

        // Диагностика на неактивных слотах: пусто.
        const auto dg = interactions::compute_fiber_friction_diagnostics(exec, sc.state, sc.upd, sc.dv, sc.slots,
                                                                        bead_counts(exec, sc));
        check(dg.active_count == 0 && dg.stick_count == 0 && dg.slip_count == 0 && dg.shared_bead_count == 0,
              "диагностика на неактивных слотах: нули");
    }
}

// ----------------------------------------------------------------------
// Mass splitting (перенесён из задачи 4): два активных слота в ОДНОЙ точке у
// общей бусины 3 -- сегменты (2,3) при s=1 и (3,4) при s=0 против одного
// сегмента B (10,11) -- сцена test_e2e_cross_friction (c)/(d). Без mass
// splitting (c_k=1) каждый контакт берёт полный импульс прилипания и
// суммарная поправка вдвое больше нужной (перескок); с ним каждый видит
// половину массы общих бусин (c=2), суммарный импульс ровно останавливает
// касательное скольжение.
// ----------------------------------------------------------------------
void test_mass_splitting()
{
    std::printf("test_mass_splitting:\n");
    rhs::CpuExecutor exec(2);
    Scene sc;
    sc.init(12, 1.0e-3);
    const double m[12] = {1.0, 1.0, 2.0, 1.5, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.7, 1.9};
    for (int i = 0; i < 12; ++i) sc.state.mass[i] = m[i];
    Rng rng(555);
    for (std::size_t d = 0; d < dim; ++d)
        for (int i = 0; i < 12; ++i)
        {
            sc.base[d][i] = rng.range(-1.0, 1.0);
            sc.dv[d][i] = rng.range(-3.0, 3.0);
        }
    double n[3] = {0.3, -0.5, 0.8};
    unit(n);
    auto make = [&](std::size_t a0, long long a1, double s) {
        interactions::FiberFrictionContact c;
        c.active = true;
        c.a0 = a0;
        c.a1 = a1;
        c.b0 = 10;
        c.b1 = 11;
        c.s = s;
        c.t = 0.4;
        for (int k = 0; k < 3; ++k) c.n[k] = n[k];
        c.N = 2.0;
        c.mu = 1.0e6; // muN*h заведомо больше J_stick: режим прилипания
        return c;
    };
    sc.slots.clear();
    sc.slots.push_back(make(2, 3, 1.0)); // точка контакта в бусине 3 (вес бусины 2 равен нулю)
    sc.slots.push_back(make(3, 4, 0.0)); // та же точка (вес бусины 4 равен нулю)

    // c_k: бусины 3, 10, 11 -- по 2 слота; бусины 2 и 4 (нулевой вес) -- 0.
    const std::vector<double> cnt = bead_counts(exec, sc);
    check(cnt[3] == 2.0 && cnt[10] == 2.0 && cnt[11] == 2.0, "c_3 = c_10 = c_11 = 2");
    check(cnt[2] == 0.0 && cnt[4] == 0.0, "бусины с нулевым весом (2 и 4) не считаются: c = 0");
    check(cnt[0] == 0.0 && cnt[5] == 0.0, "бусины вне слотов: c = 0");

    double ut0[3];
    std::vector<double> vtrial[dim];
    for (std::size_t d = 0; d < dim; ++d)
    {
        vtrial[d].resize(12);
        for (int i = 0; i < 12; ++i) vtrial[d][i] = sc.base[d][i] + sc.upd.h * sc.dv[d][i];
    }
    const double ut0n = contact_ut(vtrial, sc.slots[0], ut0);
    double ut0b[3];
    const double ut0bn = contact_ut(vtrial, sc.slots[1], ut0b);
    check(std::fabs(ut0n - ut0bn) <= 1e-14 * ut0n, "оба слота -- одна и та же точка контакта (равные u_t)");

    // (1) с mass splitting -- штатная стадия
    const StageResult r = apply_stage(exec, sc);
    double ut1[3];
    const double ut1n = contact_ut(r.vnew, sc.slots[0], ut1);
    const double along = dot3(ut1, ut0) / ut0n; // компонента u_t^{n+1} вдоль исходной u_t
    std::printf("    |u_t| до %.6e, после (c=2) %.3e, вдоль исходной %.3e\n", ut0n, ut1n, along);
    check(along >= -1e-12 * ut0n, "с mass splitting: касательная скорость не сменила знак (нет перескока; допуск округления 1e-12)");
    check(ut1n <= ut0n, "с mass splitting: |u_t| не больше исходной");
    check(ut1n <= 1e-12 * ut0n, "с mass splitting: точное прилипание (|u_t| после <= 1e-12 исходной)");

    // (2) самопроверка сцены: c_k=1 -- прямой вызов формулы, вклады суммируются вручную
    std::vector<double> dv1[dim];
    for (std::size_t d = 0; d < dim; ++d) dv1[d] = sc.dv[d];
    for (const auto& c : sc.slots)
    {
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, c, imp); // bead_count == nullptr: c_k = 1
        interactions::fiber_friction_for_each_bead(
            c, [&](rhs::index_t bead, double weight, double sign, std::uint16_t) {
                for (std::size_t d = 0; d < dim; ++d) dv1[d][bead] += sign * weight * imp.Ft[d] / sc.state.mass[bead];
            });
    }
    std::vector<double> v1[dim];
    for (std::size_t d = 0; d < dim; ++d)
    {
        v1[d].resize(12);
        for (int i = 0; i < 12; ++i) v1[d][i] = sc.base[d][i] + sc.upd.h * dv1[d][i];
    }
    double ut2[3];
    const double ut2n = contact_ut(v1, sc.slots[0], ut2);
    const double along1 = dot3(ut2, ut0) / ut0n;
    std::printf("    без mass splitting (c=1): |u_t| после %.6e, вдоль исходной %.6e (отношение %.3f)\n", ut2n, along1,
                along1 / ut0n);
    check(along1 < 0.0, "самопроверка: при c_k=1 в этой сцене происходит перескок (знак u_t сменился)");
    check(std::fabs(ut2n - ut0n) <= 1e-9 * ut0n, "самопроверка: при c_k=1 перекомпенсация вдвое (|u_t| после == |u_t| до)");

    // (3) одиночный контакт: c_k = 1 из bead_counts даёт побитово то же, что nullptr
    Scene one;
    build_scene_8_1(one, 1.0e-3, 2.0, 1000.0);
    const std::vector<double> c1 = bead_counts(exec, one);
    interactions::FiberFrictionImpulse ia, ib;
    interactions::fiber_friction_impulse(one.state, one.upd, one.dv, one.slots[0], ia);
    interactions::fiber_friction_impulse(one.state, one.upd, one.dv, one.slots[0], ib, c1.data());
    check(c1[0] == 1.0 && c1[1] == 1.0 && c1[2] == 1.0 && c1[3] == 1.0, "одиночный контакт: c_k = 1 на всех четырёх бусинах");
    check((ia.W == ib.W && ia.J == ib.J && ia.Ft[0] == ib.Ft[0] &&
                                                     ia.Ft[1] == ib.Ft[1] && ia.Ft[2] == ib.Ft[2]),
          "одиночный контакт: результат с c_k=1 побитово совпадает с формулой без mass splitting");
}

// ----------------------------------------------------------------------
void test_diagnostics()
{
    std::printf("test_diagnostics:\n");
    rhs::CpuExecutor exec(2);
    Rng rng(4242);

    Scene sc;
    sc.init(15, 1.0e-3);
    for (int i = 0; i < 15; ++i) sc.state.mass[i] = rng.range(0.5, 4.0);
    for (std::size_t d = 0; d < dim; ++d)
        for (int i = 0; i < 15; ++i)
        {
            sc.base[d][i] = rng.range(-1.0, 1.0);
            sc.dv[d][i] = rng.range(-3.0, 3.0);
        }

    auto make = [&](std::size_t a0, long long a1, std::size_t b0, long long b1, double s, double t, double N,
                    double mu, bool active) {
        interactions::FiberFrictionContact c;
        c.active = active;
        c.a0 = a0;
        c.a1 = a1;
        c.b0 = b0;
        c.b1 = b1;
        c.s = s;
        c.t = t;
        double n[3] = {rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)};
        unit(n);
        for (int k = 0; k < 3; ++k) c.n[k] = n[k];
        c.N = N;
        c.mu = mu;
        return c;
    };
    sc.slots.clear();
    sc.slots.push_back(make(0, 1, 2, 3, 0.3, 0.6, 2.0, 1.0e6, true));  // c0: прилипание
    sc.slots.push_back(make(1, 4, 5, 6, 0.5, 0.5, 2.0, 0.05, true));    // c1: скольжение; общая бусина 1 с c0
    sc.slots.push_back(make(7, 8, 5, 9, 0.4, 0.2, 2.0, 1.0e6, true));  // c2: прилипание; общая бусина 5 с c1
    sc.slots.push_back(make(10, 11, 12, 13, 0.5, 0.5, 0.0, 0.5, true)); // c3: N=0 -- активен, F_t=0
    sc.slots.push_back(make(0, 1, 2, 3, 0.5, 0.5, 2.0, 0.5, false));   // c4: неактивен (бусины c0 не считать)

    const auto dg = interactions::compute_fiber_friction_diagnostics(exec, sc.state, sc.upd, sc.dv, sc.slots,
                                                                        bead_counts(exec, sc));
    std::printf("    active %zu, stick %zu, slip %zu, shared %zu, max_force_ratio %.15f\n", dg.active_count,
                dg.stick_count, dg.slip_count, dg.shared_bead_count, dg.max_force_ratio);
    check(dg.active_count == 4, "active_count == 4 (неактивный слот не считается)");

    // ожидаемое stick/slip -- независимо, через fiber_friction_impulse по слотам
    // (с теми же c_k, что и диагностика)
    const std::vector<double> dg_counts = bead_counts(exec, sc);
    std::size_t st = 0, sl = 0;
    for (const auto& c : sc.slots)
    {
        if (!c.active) continue;
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.upd, sc.dv, c, imp, dg_counts.data());
        if (imp.has_force) (imp.stick ? st : sl) += 1;
    }
    check(st == 2 && sl == 1, "по слотам: 2 прилипания, 1 скольжение (контакт с N=0 без силы)");
    check(dg.stick_count == 2, "stick_count == 2");
    check(dg.slip_count == 1, "slip_count == 1");
    check(dg.shared_bead_count == 2, "shared_bead_count == 2 (бусины 1 и 5)");
    check(dg.max_force_ratio <= 1.0 + 1e-14, "max_force_ratio <= 1");
    // скользящий контакт даёт ровно 1
    check(std::fabs(dg.max_force_ratio - 1.0) <= 1e-14, "max_force_ratio == 1 (есть скользящий контакт)");

    // Без общих бусин: shared == 0
    sc.slots.erase(sc.slots.begin() + 2);
    const auto dg2 = interactions::compute_fiber_friction_diagnostics(exec, sc.state, sc.upd, sc.dv, sc.slots,
                                                                        bead_counts(exec, sc));
    check(dg2.shared_bead_count == 1, "после удаления c2 остаётся одна общая бусина (1)");
    sc.slots.erase(sc.slots.begin() + 1);
    const auto dg3 = interactions::compute_fiber_friction_diagnostics(exec, sc.state, sc.upd, sc.dv, sc.slots,
                                                                        bead_counts(exec, sc));
    check(dg3.shared_bead_count == 0, "без пересекающихся по бусинам контактов shared_bead_count == 0");
}

// ============================================================================
// ТРЕНИЕ НИТЬ-SOLID (fiber/fiber_solid_friction_implementation.md, задача 4,
// §7.1). Стадия -- единый проход compute_friction_contributions (оба вида
// контактов), c_k/c_s -- compute_friction_counts (mass splitting), пробная
// скорость Solid -- ode::SolidVelocityUpdate: (base + h*dvx) + av. Скорость
// после обновления считается ТОЙ ЖЕ арифметикой, что в решателе:
// нить -- base + h*(dv + a^f), Solid -- (base + h*(dvx + a^f)) + av.
// ============================================================================

// Сцена: nf бусин, ns частиц Solid. fupd/supd указывают на массивы сцены
// (сцену НЕ копировать после init).
struct SolidScene
{
    fiber::FiberState state;
    ParticleSystem sys;
    std::vector<double> base[dim];  // база скорости нити
    std::vector<double> dv[dim];    // ускорение нити без трения
    std::vector<double> sbase[dim]; // база скорости Solid
    std::vector<std::vector<double>> dvx; // ускорение Solid без трения
    std::vector<std::vector<double>> av;  // XSPH-поправка
    std::vector<interactions::FiberFrictionContact> ff_slots;
    std::vector<interactions::FiberSolidFrictionContact> fs_slots;
    ode::FiberVelocityUpdate fupd{};
    ode::SolidVelocityUpdate supd{};

    void init(std::size_t nf, std::size_t ns, double hf, double hs, bool with_av = true)
    {
        state.resize(nf);
        sys.resize(ns);
        for (std::size_t d = 0; d < dim; ++d)
        {
            base[d].assign(nf, 0.0);
            dv[d].assign(nf, 0.0);
            sbase[d].assign(ns, 0.0);
            fupd.base[d] = &base[d];
            supd.base[d] = &sbase[d];
        }
        dvx.assign(dim, std::vector<double>(ns, 0.0));
        av.assign(dim, std::vector<double>(ns, 0.0));
        for (std::size_t i = 0; i < ns; ++i) sys.mass[i] = 1.0;
        fupd.h = hf;
        supd.h = hs;
        supd.av = with_av ? &av : nullptr;
    }
};

struct SolidStage
{
    std::vector<fiber::FiberDelta> fdelta;      // a^f нити
    std::vector<interactions::SolidDelta> sdelta; // a^f Solid
    std::size_t solid_emitted = 0;              // число вкладов в Solid (type 13)
    std::vector<double> bcnt, scnt;             // c_k, c_s
    std::vector<double> fvtrial[dim], fvnew[dim]; // нить: без трения / после
    std::vector<double> svtrial[dim], svnew[dim]; // Solid: без трения / после
};

// no_split=true -- самопроверка: все c=1 (как если бы mass splitting не было).
SolidStage apply_solid_stage(const rhs::CpuExecutor& exec, const SolidScene& sc, bool kinematic = false,
                             bool no_split = false)
{
    SolidStage r;
    const std::size_t nf = sc.state.size();
    const std::size_t ns = sc.dvx[0].size();
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> fcb;
    std::vector<rhs::ContributionBuffer<interactions::SolidCountDelta>> scb;
    interactions::compute_friction_counts(exec, nf, ns, sc.ff_slots, sc.fs_slots, fcb, scb, r.bcnt, r.scnt);
    if (no_split)
    {
        r.bcnt.assign(nf, 1.0);
        r.scnt.assign(ns, 1.0);
    }
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fb;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sb;
    interactions::compute_friction_contributions(exec, sc.state, sc.fupd, sc.dv, sc.ff_slots, sc.sys, sc.supd, sc.dvx,
                                                 sc.fs_slots, r.bcnt, r.scnt, kinematic, fb, sb);
    r.fdelta = rhs::reduce_contributions(exec, fb, nf);
    r.sdelta = rhs::reduce_contributions(exec, sb, ns);
    for (const auto& b : sb) r.solid_emitted += b.size();
    for (std::size_t d = 0; d < dim; ++d)
    {
        r.fvtrial[d].resize(nf);
        r.fvnew[d].resize(nf);
        r.svtrial[d].resize(ns);
        r.svnew[d].resize(ns);
        for (std::size_t i = 0; i < nf; ++i)
        {
            r.fvtrial[d][i] = (*sc.fupd.base[d])[i] + sc.fupd.h * sc.dv[d][i];
            r.fvnew[d][i] = (*sc.fupd.base[d])[i] + sc.fupd.h * (sc.dv[d][i] + r.fdelta[i].dv[d]);
        }
        for (std::size_t i = 0; i < ns; ++i)
        {
            double t = (*sc.supd.base[d])[i] + sc.supd.h * sc.dvx[d][i];
            double n = (*sc.supd.base[d])[i] + sc.supd.h * (sc.dvx[d][i] + r.sdelta[i].dv[d]);
            if (sc.supd.av)
            {
                t += (*sc.supd.av)[d][i];
                n += (*sc.supd.av)[d][i];
            }
            r.svtrial[d][i] = t;
            r.svnew[d][i] = n;
        }
    }
    return r;
}

// Касательная относительная скорость точки контакта Fiber-Solid (нить
// относительно Solid) по массивам скоростей; возвращает |u_t|.
double fs_ut(const std::vector<double> (&fv)[dim], const std::vector<double> (&sv)[dim],
             const interactions::FiberSolidFrictionContact& c, double ut[3])
{
    double u[3];
    for (std::size_t d = 0; d < dim; ++d)
    {
        double vF = fv[d][c.a];
        if (c.b != -1) vF = (1.0 - c.lambda) * vF + c.lambda * fv[d][static_cast<std::size_t>(c.b)];
        u[d] = vF - sv[d][c.s];
    }
    const double un = dot3(u, c.n);
    for (std::size_t d = 0; d < dim; ++d) ut[d] = u[d] - un * c.n[d];
    return norm3(ut);
}

interactions::FiberSolidFrictionContact make_fs_slot(std::size_t a, long long b, std::size_t s, double lambda,
                                                     const double n[3], double N, double mu)
{
    interactions::FiberSolidFrictionContact c;
    c.active = true;
    c.a = a;
    c.b = b;
    c.s = s;
    c.lambda = lambda;
    for (int k = 0; k < 3; ++k) c.n[k] = n[k];
    c.N = N;
    c.mu = mu;
    return c;
}

// Одиночный контакт Fiber-Solid: бусины 0 (a), 1 (b), частица Solid 0.
// m3 = {m_a, m_b, m_s}; ненулевые базовые скорости, dv, dvx и av.
void build_fs_scene(SolidScene& sc, double h, double N, double mu, double lambda = 0.3, const double* m3 = nullptr,
                    bool with_av = true)
{
    const double dm[3] = {1.0, 3.0, 2.0};
    const double* m = m3 ? m3 : dm;
    sc.init(2, 1, h, h, with_av);
    sc.state.mass[0] = m[0];
    sc.state.mass[1] = m[1];
    sc.sys.mass[0] = m[2];
    const double va[3] = {1.0, 0.2, -0.4}, vb[3] = {0.5, -0.3, 0.7}, vs[3] = {-0.2, 0.9, 0.1};
    const double aa[3] = {5.0, -2.0, 1.0}, ab[3] = {-1.0, 3.0, 0.5}, as_[3] = {2.0, 2.0, -4.0};
    const double avs[3] = {0.03, -0.02, 0.01};
    for (std::size_t d = 0; d < dim; ++d)
    {
        sc.base[d][0] = va[d];
        sc.base[d][1] = vb[d];
        sc.dv[d][0] = aa[d];
        sc.dv[d][1] = ab[d];
        sc.sbase[d][0] = vs[d];
        sc.dvx[d][0] = as_[d];
        sc.av[d][0] = avs[d];
    }
    double n[3] = {0.3, -0.5, 0.8};
    unit(n);
    sc.fs_slots.assign(1, make_fs_slot(0, 1, 0, lambda, n, N, mu));
}

interactions::FiberFrictionImpulse fs_impulse(const SolidScene& sc, const SolidStage& r, std::size_t k,
                                              bool kinematic = false)
{
    interactions::FiberFrictionImpulse imp;
    interactions::fiber_solid_friction_impulse(sc.state, sc.fupd, sc.dv, sc.sys, sc.supd, sc.dvx, sc.fs_slots[k], imp,
                                               r.bcnt.data(), r.scnt.data(), kinematic);
    return imp;
}

// ----------------------------------------------------------------------
// FS.1 Прилипание: после обновления касательная относительная скорость точки
// контакта равна нулю (отн. 1e-12), |F_t| < mu*N. Разные массы (в том числе
// крайние: лёгкий/тяжёлый Solid) и lambda.
// ----------------------------------------------------------------------
void test_fs_stick(double& worst_out)
{
    std::printf("test_fs_stick:\n");
    rhs::CpuExecutor exec(2);
    const double masses[4][3] = {{1.0, 3.0, 2.0}, {1.0e-3, 1.0e-3, 10.0}, {5.0, 0.2, 1.0e-2}, {1.0, 1.0, 1.0e6}};
    const double lambdas[3] = {0.3, 0.9, 0.5};
    double worst = 0.0;
    int n_stick = 0, n_cases = 0, n_exact = 0, n_lt = 0, n_stop = 0, n_big = 0;
    for (const auto& m : masses)
        for (double lam : lambdas)
        {
            SolidScene sc;
            build_fs_scene(sc, 1.0e-3, 2.0, 1.0e6, lam, m);
            const SolidStage r = apply_solid_stage(exec, sc);
            double ut0[3], ut1[3];
            const double n0 = fs_ut(r.fvtrial, r.svtrial, sc.fs_slots[0], ut0);
            const double n1 = fs_ut(r.fvnew, r.svnew, sc.fs_slots[0], ut1);
            const auto imp = fs_impulse(sc, r, 0);
            const double Ft = norm3(imp.Ft);
            const double muN = sc.fs_slots[0].mu * sc.fs_slots[0].N;
            ++n_cases;
            if (n0 > 1e-3) ++n_big;
            if (imp.stick) ++n_stick;
            if (n1 <= 1e-12 * n0) ++n_exact;
            if (Ft < muN) ++n_lt;
            if (std::fabs(Ft * sc.fupd.h * imp.W - n0) <= 1e-12 * n0) ++n_stop;
            worst = std::max(worst, n1 / n0);
        }
    worst_out = worst;
    std::printf("    %d случаев (4 набора масс x 3 lambda), max |u_t|после/до = %.3e\n", n_cases, worst);
    check(n_big == n_cases, "исходная касательная скорость ненулевая во всех случаях");
    check(n_stick == n_cases, "выбрана ветка прилипания во всех случаях");
    check(n_exact == n_cases, "после обновления |u_t| <= 1e-12*|u_t до| (точное прилипание) во всех случаях");
    check(n_lt == n_cases, "|F_t| < mu*N во всех случаях");
    check(n_stop == n_cases, "|F_t|*h*W == |u_t| (остановка ровно) во всех случаях");
}

// ----------------------------------------------------------------------
// FS.2 Скольжение: |F_t| = mu*N, направление -u_t/|u_t|, u_t не сменила знак;
// сила на Solid -F_t.
// ----------------------------------------------------------------------
void test_fs_slip(double& ratio_err_out)
{
    std::printf("test_fs_slip:\n");
    rhs::CpuExecutor exec(2);
    SolidScene sc;
    build_fs_scene(sc, 1.0e-3, 2.0, 0.25); // muN = 0.5, muN*h*W << |u_t|
    const SolidStage r = apply_solid_stage(exec, sc);
    const auto& c = sc.fs_slots[0];

    double ut0[3], ut1[3];
    const double n0 = fs_ut(r.fvtrial, r.svtrial, c, ut0);
    const double n1 = fs_ut(r.fvnew, r.svnew, c, ut1);
    const double muN = c.mu * c.N;
    const double h = sc.fupd.h;
    const auto imp = fs_impulse(sc, r, 0);
    const double W = imp.W;
    std::printf("    |u_t| до %.6e, после %.6e; muN*h*W = %.6e\n", n0, n1, muN * h * W);
    check(!imp.stick, "выбрана ветка скольжения");
    const double Ft = norm3(imp.Ft);
    ratio_err_out = std::fabs(Ft / muN - 1.0);
    std::printf("    |F_t|/muN - 1 = %.3e\n", Ft / muN - 1.0);
    check(ratio_err_out <= 1e-14, "|F_t| == mu*N (отн. 1e-14)");

    // сила по фактическому вкладу бусины a: m_a a^f = (1-lambda) F_t
    double fa[3], fs_[3];
    for (int d = 0; d < 3; ++d)
    {
        fa[d] = r.fdelta[c.a].dv[d] * sc.state.mass[c.a] / (1.0 - c.lambda);
        fs_[d] = -r.sdelta[c.s].dv[d] * sc.sys.mass[c.s]; // на Solid действует -F_t
    }
    check(std::fabs(norm3(fa) / muN - 1.0) <= 1e-14, "|F_t| == mu*N по вкладу бусины a");
    double uhat[3] = {ut0[0] / n0, ut0[1] / n0, ut0[2] / n0};
    check(std::fabs(dot3(fa, uhat) / (-norm3(fa)) - 1.0) <= 1e-13, "направление F_t == -u_t/|u_t|");
    double dfs[3] = {fa[0] - fs_[0], fa[1] - fs_[1], fa[2] - fs_[2]};
    check(norm3(dfs) <= 1e-13 * muN, "сила на Solid == -F_t (по вкладу частицы Solid)");

    double expect[3], diff[3];
    for (int d = 0; d < 3; ++d)
    {
        expect[d] = ut0[d] - muN * h * W * uhat[d];
        diff[d] = ut1[d] - expect[d];
    }
    check(norm3(diff) <= 1e-10 * muN * h * W, "касательная скорость уменьшилась ровно на muN*h*W");
    check(dot3(ut1, uhat) > 0.0, "касательная относительная скорость не сменила знак");
    check(std::fabs(n0 - n1 - muN * h * W) <= 1e-10 * muN * h * W, "|u_t| уменьшился ровно на muN*h*W");
}

// ----------------------------------------------------------------------
// FS.3 Статическое трение: базовая относительная скорость ноль, тянущее
// касательное ускорение g_t задано ускорением бусин (вариант A) либо частицы
// Solid в противоположную сторону (вариант B). Ниже порога mu*N -- прилипание
// (F_t уравновешивает m_eff*g_t), выше -- скольжение с mu*N.
// ----------------------------------------------------------------------
void test_fs_static_under_force(double& stick_err_out)
{
    std::printf("test_fs_static_under_force:\n");
    rhs::CpuExecutor exec(2);
    const double h = 1.0e-3, N = 2.5, mu = 0.4;
    const double muN = mu * N;
    for (int who = 0; who < 2; ++who)
        for (int variant = 0; variant < 2; ++variant)
        {
            SolidScene sc;
            build_fs_scene(sc, h, N, mu, 0.35, nullptr, /*with_av=*/false);
            for (std::size_t d = 0; d < dim; ++d)
            {
                for (std::size_t i = 0; i < 2; ++i)
                {
                    sc.base[d][i] = 0.0;
                    sc.dv[d][i] = 0.0;
                }
                sc.sbase[d][0] = 0.0;
                sc.dvx[d][0] = 0.0;
            }
            const auto& c = sc.fs_slots[0];
            const auto imp0 = [&] {
                const SolidStage r0 = apply_solid_stage(exec, sc);
                return fs_impulse(sc, r0, 0);
            }();
            const double W = imp0.W;
            const double m_eff = 1.0 / W;
            double tt[3] = {0.0, c.n[2], -c.n[1]}; // n x x_hat, касательный
            unit(tt);
            const double frac = (variant == 0) ? 0.5 : 2.0; // m_eff|g_t| = frac*muN
            const double gt = frac * muN / m_eff;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double g = gt * tt[d] + 3.0 * c.n[d];
                if (who == 0)
                {
                    sc.dv[d][0] = g; // обе бусины ускоряются одинаково: точка контакта получает g при любой lambda
                    sc.dv[d][1] = g;
                }
                else
                    sc.dvx[d][0] = -g; // Solid ускоряется в противоположную сторону: относительное ускорение g
            }
            const SolidStage r = apply_solid_stage(exec, sc);
            double ut0[3], ut1[3];
            const double n0 = fs_ut(r.fvtrial, r.svtrial, c, ut0);
            const double n1 = fs_ut(r.fvnew, r.svnew, c, ut1);
            const auto imp = fs_impulse(sc, r, 0);
            const double Ft = norm3(imp.Ft);
            std::printf("    %s, m_eff|g_t|/muN = %.2f: |u_t| до %.6e, после %.6e, |F_t|/muN = %.15f\n",
                        who == 0 ? "тянет нить" : "тянет Solid", frac, n0, n1, Ft / muN);
            check(std::fabs(n0 - h * gt) <= 1e-12 * h * gt, "пробная касательная скорость == h*g_t");
            if (variant == 0)
            {
                if (who == 0) stick_err_out = n1 / n0;
                check(n1 <= 1e-12 * n0, "m_eff|g| < muN: относительная касательная скорость после шага 0");
                check(imp.stick && Ft < muN, "прилипание, |F_t| < muN");
                check(std::fabs(Ft - m_eff * gt) <= 1e-12 * muN, "F_t уравновешивает m_eff*g_t");
            }
            else
            {
                check(!imp.stick, "m_eff|g| > muN: скольжение");
                check(std::fabs(Ft / muN - 1.0) <= 1e-14, "|F_t| == mu*N");
                const double expected = h * (gt - W * muN);
                check(std::fabs(n1 - expected) <= 1e-10 * expected, "остаточная скорость h*(g_t - W*muN)");
                check(dot3(ut1, tt) > 0.0, "остаточная скорость направлена вдоль g_t");
            }
        }
}

// ----------------------------------------------------------------------
// FS.4 XSPH: ненулевой av частицы Solid входит в пробную скорость; прилипание
// точное С УЧЁТОМ av. Самопроверка: если пробная скорость не знает об av
// (supd.av = nullptr), а решатель av прибавит, касательная остаётся ненулевой.
// ----------------------------------------------------------------------
void test_fs_xsph(double& err_out)
{
    std::printf("test_fs_xsph:\n");
    rhs::CpuExecutor exec(2);
    const double avbig[3] = {0.4, -0.7, 0.2};
    SolidScene sc;
    build_fs_scene(sc, 1.0e-3, 2.0, 1.0e6);
    for (std::size_t d = 0; d < dim; ++d) sc.av[d][0] = avbig[d];
    const auto& c = sc.fs_slots[0];

    // Пробная скорость Solid ровно (base + h*dvx) + av; u_t считается независимо.
    const SolidStage r = apply_solid_stage(exec, sc);
    double vs_expect[3];
    for (std::size_t d = 0; d < dim; ++d) vs_expect[d] = (sc.sbase[d][0] + sc.supd.h * sc.dvx[d][0]) + sc.av[d][0];
    bool trial_ok = true;
    for (std::size_t d = 0; d < dim; ++d) trial_ok = trial_ok && (r.svtrial[d][0] == vs_expect[d]);
    check(trial_ok, "пробная скорость Solid в тесте == (base + h*dvx) + av");

    double ut0[3], ut1[3];
    const double n0 = fs_ut(r.fvtrial, r.svtrial, c, ut0);
    const double n1 = fs_ut(r.fvnew, r.svnew, c, ut1);
    const auto imp = fs_impulse(sc, r, 0);
    check(std::fabs(imp.ut_norm - n0) <= 1e-13 * n0, "|u_t| стадии == |u_t| по (base + h*dvx) + av (av учтён)");
    check(imp.stick, "прилипание");
    err_out = n1 / n0;
    std::printf("    |u_t| до %.6e, после %.3e (отношение %.3e); |av_t| = %.3e\n", n0, n1, n1 / n0,
                [&] {
                    double t[3];
                    const double un = dot3(avbig, c.n);
                    for (int d = 0; d < 3; ++d) t[d] = avbig[d] - un * c.n[d];
                    return norm3(t);
                }());
    check(n1 <= 1e-12 * n0, "прилипание точное с учётом av: |u_t после| <= 1e-12*|u_t до|");

    // Влияние av: без av пробная u_t другая.
    SolidScene sc2;
    build_fs_scene(sc2, 1.0e-3, 2.0, 1.0e6);
    for (std::size_t d = 0; d < dim; ++d) sc2.av[d][0] = avbig[d];
    sc2.supd.av = nullptr;
    const SolidStage r2 = apply_solid_stage(exec, sc2);
    const auto imp2 = fs_impulse(sc2, r2, 0);
    check(std::fabs(imp2.ut_norm - imp.ut_norm) > 1e-2 * imp.ut_norm, "av существенно меняет пробную u_t (тест не тривиален)");
    // Решатель av прибавит: остаточная касательная == касательная часть av (пробная про av не знала).
    std::vector<double> svfinal[dim];
    for (std::size_t d = 0; d < dim; ++d) svfinal[d] = {r2.svnew[d][0] + avbig[d]};
    double ut2[3];
    const double n2 = fs_ut(r2.fvnew, svfinal, sc2.fs_slots[0], ut2);
    check(n2 > 1e-2, "самопроверка: без учёта av в пробной скорости прилипание НЕ точное (остаток порядка |av_t|)");
}

// ----------------------------------------------------------------------
// Сцена на реальной геометрии: сегмент нити (2 бусины вдоль x) и одна частица
// Solid над ним (Al, itype 4); слоты заполняет kernel нормального контакта.
// ----------------------------------------------------------------------
struct KernelScene
{
    fiber::FiberState state;
    fiber::FiberTopology topo;
    std::vector<fiber::FiberMaterial> materials;
    ParticleSystem sys;
    Materials mats;
    std::vector<interactions::FiberSolidContactPrimitive> prims;
    std::vector<interactions::FiberSolidFrictionContact> slots;
    std::vector<fiber::FiberDelta> fd;         // нормальная сила: нить
    std::vector<interactions::SolidDelta> sd;  // нормальная сила: Solid
    bool has_active = false;
};

void build_kernel_scene(KernelScene& ks, rhs::CpuExecutor& exec, double mu_fs, bool with_contact = true)
{
    const double d_f = 1.0e-3, l0 = 2.0e-3, rho_f = 1200.0;
    const double mass_f = rho_f * pi * d_f * d_f / 4.0 * l0;
    ks.state.resize(2);
    ks.topo.resize(2);
    for (int i = 0; i < 2; ++i)
    {
        ks.state.x[0][i] = i * l0;
        ks.state.x[1][i] = 0.0;
        ks.state.x[2][i] = 0.0;
        for (std::size_t d = 0; d < dim; ++d) ks.state.vx[d][i] = 0.0;
        ks.state.vx[0][i] = 0.3;
        ks.state.mass[i] = mass_f;
        ks.state.object_id[i] = 0;
        ks.state.material_id[i] = 0;
    }
    ks.topo.next[0] = 1;
    ks.topo.prev[1] = 0;
    fiber::FiberMaterial m;
    m.E = 2.0e7;
    m.nu = 0.3;
    m.d = d_f;
    m.rho = rho_f;
    m.eps_break = 1.0;
    m.mu_fs = mu_fs;
    ks.materials.assign(1, m);

    const double dxp = 1.5e-3, rho_s = 2700.0;
    ks.sys.resize(1);
    ks.sys.inzializ_index(0);
    ks.sys.x[0][0] = 0.7 * l0;                     // lambda = 0.7
    ks.sys.x[1][0] = 0.0;
    ks.sys.x[2][0] = with_contact ? 1.0e-3 : 5.0e-3; // над нитью: D = 1e-3 < r_f + r_s = 1.25e-3 (или далеко)
    ks.sys.mass[0] = rho_s * dxp * dxp * dxp;
    ks.sys.rho[0] = rho_s;
    ks.sys.rho_0[0] = rho_s;
    ks.sys.itype[0] = 4;
    ks.sys.hsml[0] = 1.0;

    interactions::FiberSolidContactPrimitive p;
    p.a = 0;
    p.b = 1;
    p.s = 0;
    p.has_b = true;
    p.a_prev = -1;
    p.b_is_last = true;
    ks.prims.assign(1, p);

    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
    ks.has_active = interactions::compute_fiber_solid_contact_contributions(exec, ks.mats, ks.sys, ks.state,
                                                                             ks.materials, ks.prims, fbuf, sbuf, ks.slots);
    ks.fd = rhs::reduce_contributions(exec, fbuf, ks.state.size());
    ks.sd = rhs::reduce_contributions(exec, sbuf, ks.sys.size());
}

// SolidScene из KernelScene: dv/dvx -- нормальные ускорения kernel'а, базы --
// скорости нити (Solid покоится), av = 0.
void solid_scene_from_kernel(SolidScene& sc, const KernelScene& ks, double h)
{
    sc.init(2, 1, h, h, /*with_av=*/false);
    for (std::size_t i = 0; i < 2; ++i) sc.state.mass[i] = ks.state.mass[i];
    sc.sys.mass[0] = ks.sys.mass[0];
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (std::size_t i = 0; i < 2; ++i)
        {
            sc.base[d][i] = ks.state.vx[d][i];
            sc.dv[d][i] = ks.fd[i].dv[d];
        }
        sc.dvx[d][0] = ks.sd[0].dv[d];
    }
    sc.fs_slots = ks.slots;
}

// ----------------------------------------------------------------------
// FS.5 Знаки (таблица §3): нормальная и тангенциальная части на a, b, s.
//   бусина a: -(1-lambda) N n   и  +(1-lambda) F_t
//   бусина b: -lambda N n       и  +lambda F_t
//   Solid s:  +N n              и  -F_t
// n -- от нити к Solid. Оба режима: скольжение (mu мало) и прилипание.
// ----------------------------------------------------------------------
void test_fs_signs()
{
    std::printf("test_fs_signs:\n");
    rhs::CpuExecutor exec(2);
    for (int mode = 0; mode < 2; ++mode)
    {
        KernelScene ks;
        build_kernel_scene(ks, exec, mode == 0 ? 0.1 : 100.0);
        const char* label = mode == 0 ? "скольжение" : "прилипание";
        check(ks.has_active && ks.slots.size() == 1 && ks.slots[0].active,
              std::string(label) + ": kernel заполнил активный слот");
        const auto& c = ks.slots[0];
        check(c.a == 0 && c.b == 1 && c.s == 0, std::string(label) + ": слот: a=0, b=1, s=0");
        check(std::fabs(c.lambda - 0.7) <= 1e-9, std::string(label) + ": lambda == 0.7");
        check(std::fabs(c.n[2] - 1.0) <= 1e-12 && std::fabs(c.n[0]) <= 1e-12 && std::fabs(c.n[1]) <= 1e-12,
              std::string(label) + ": n = +z (от нити к Solid, Solid над нитью)");
        check(c.N > 0.0, std::string(label) + ": N > 0 (N = " + std::to_string(c.N) + ")");

        // --- нормальная часть по таблице -------------------------------
        const double lam = c.lambda;
        double ea[3], eb[3], es[3], ra[3], rb[3], rs[3];
        for (std::size_t d = 0; d < dim; ++d)
        {
            ea[d] = -(1.0 - lam) * c.N * c.n[d];
            eb[d] = -lam * c.N * c.n[d];
            es[d] = +c.N * c.n[d];
            ra[d] = ks.state.mass[0] * ks.fd[0].dv[d];
            rb[d] = ks.state.mass[1] * ks.fd[1].dv[d];
            rs[d] = ks.sys.mass[0] * ks.sd[0].dv[d];
        }
        auto dist = [](const double* x, const double* y) {
            double dd[3] = {x[0] - y[0], x[1] - y[1], x[2] - y[2]};
            return norm3(dd);
        };
        check(dist(ra, ea) <= 1e-9 * c.N, std::string(label) + ": нормальная часть на a: -(1-lambda) N n");
        check(dist(rb, eb) <= 1e-9 * c.N, std::string(label) + ": нормальная часть на b: -lambda N n");
        check(dist(rs, es) <= 1e-9 * c.N, std::string(label) + ": нормальная часть на Solid: +N n");

        // --- тангенциальная часть --------------------------------------
        SolidScene sc;
        solid_scene_from_kernel(sc, ks, 1.0e-7);
        const SolidStage r = apply_solid_stage(exec, sc);
        const auto imp = fs_impulse(sc, r, 0);
        const double Ft = norm3(imp.Ft);
        check(imp.has_force && (mode == 0 ? !imp.stick : imp.stick),
              std::string(label) + ": ожидаемая ветка (" + (mode == 0 ? "slip" : "stick") + ")");
        if (mode == 0)
            check(std::fabs(Ft / (c.mu * c.N) - 1.0) <= 1e-14, "скольжение: |F_t| == mu*N");
        check(imp.Ft[0] < 0.0 && std::fabs(imp.Ft[1]) <= 1e-9 * Ft && std::fabs(imp.Ft[2]) <= 1e-9 * Ft,
              std::string(label) + ": нить скользит по +x относительно Solid -> F_t на нити вдоль -x");
        double fa[3], fb_[3], fsol[3], sum[3];
        for (std::size_t d = 0; d < dim; ++d)
        {
            fa[d] = ks.state.mass[0] * r.fdelta[0].dv[d];
            fb_[d] = ks.state.mass[1] * r.fdelta[1].dv[d];
            fsol[d] = sc.sys.mass[0] * r.sdelta[0].dv[d];
            sum[d] = fa[d] + fb_[d] + fsol[d];
        }
        double ta[3] = {(1.0 - lam) * imp.Ft[0], (1.0 - lam) * imp.Ft[1], (1.0 - lam) * imp.Ft[2]};
        double tb[3] = {lam * imp.Ft[0], lam * imp.Ft[1], lam * imp.Ft[2]};
        double ts[3] = {-imp.Ft[0], -imp.Ft[1], -imp.Ft[2]};
        check(dist(fa, ta) <= 1e-12 * Ft, std::string(label) + ": трение на a: +(1-lambda) F_t");
        check(dist(fb_, tb) <= 1e-12 * Ft, std::string(label) + ": трение на b: +lambda F_t");
        check(dist(fsol, ts) <= 1e-12 * Ft, std::string(label) + ": трение на Solid: -F_t");
        check(fsol[0] > 0.0, std::string(label) + ": Solid увлекается вдоль +x (движения нити)");
        check(norm3(sum) <= 1e-12 * Ft, std::string(label) + ": сумма m a^f по a, b, s == 0");
        check(std::fabs(dot3(imp.Ft, c.n)) <= 1e-12 * Ft, std::string(label) + ": F_t . n == 0");

        // нормаль + трение вместе: импульс сохраняется
        double tot[3], sabs = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            tot[d] = ra[d] + rb[d] + rs[d] + sum[d];
            sabs += std::fabs(ra[d]) + std::fabs(rb[d]) + std::fabs(rs[d]);
        }
        check(norm3(tot) <= 1e-12 * sabs, std::string(label) + ": нормаль + трение: полный импульс контакта сохраняется");
        std::printf("    %s: N = %.4f, mu*N = %.4f, |F_t| = %.4f, m_a a_x = %.4e, m_b a_x = %.4e, m_s a_x = %.4e\n",
                    label, c.N, c.mu * c.N, Ft, fa[0], fb_[0], fsol[0]);
    }
}

// ----------------------------------------------------------------------
// FS.6 Касательность F_t . n == 0 и F_t . u_t <= 0 (оба режима, наклонная
// нормаль, точечный контакт и lambda=0/1).
// ----------------------------------------------------------------------
void test_fs_tangential()
{
    std::printf("test_fs_tangential:\n");
    rhs::CpuExecutor exec(2);
    for (int mode = 0; mode < 2; ++mode)
        for (int shape = 0; shape < 3; ++shape)
        {
            SolidScene sc;
            build_fs_scene(sc, 1.0e-3, 2.0, mode == 0 ? 1.0e6 : 0.25, shape == 1 ? 0.0 : (shape == 2 ? 1.0 : 0.4));
            if (shape == 1) sc.fs_slots[0].b = -1; // точечный контакт (b = -1, lambda = 0)
            const SolidStage r = apply_solid_stage(exec, sc);
            const auto imp = fs_impulse(sc, r, 0);
            double ut0[3];
            fs_ut(r.fvtrial, r.svtrial, sc.fs_slots[0], ut0);
            const double Ft = norm3(imp.Ft);
            const std::string tag = std::string(mode == 0 ? "прилипание" : "скольжение") +
                                    (shape == 0 ? ", lambda=0.4" : (shape == 1 ? ", точечный контакт" : ", lambda=1"));
            check(Ft > 0.0 && std::fabs(dot3(imp.Ft, sc.fs_slots[0].n)) <= 1e-14 * Ft, tag + ": F_t . n == 0");
            check(dot3(imp.Ft, ut0) <= 0.0, tag + ": F_t . u_t <= 0");
        }
}

// ----------------------------------------------------------------------
// FS.7 Импульс: m_a a_a + m_b a_b + m_s a_s == 0; распределение по бусинам
// (веса 1-lambda, lambda; точечный контакт -- всё на a; нулевой вес -- ровно 0).
// 500 случайных контактов.
// ----------------------------------------------------------------------
void test_fs_momentum_distribution(double& mom_err_out)
{
    std::printf("test_fs_momentum_distribution:\n");
    rhs::CpuExecutor exec(2);
    Rng rng(20261004);
    double mom_err = 0.0, dist_err = 0.0;
    bool zero_ok = true;
    int n_stick = 0, n_slip = 0;
    for (int k = 0; k < 500; ++k)
    {
        SolidScene sc;
        const double h = rng.range(1e-3, 1e-1);
        sc.init(2, 1, h, h, true);
        sc.state.mass[0] = std::pow(10.0, rng.range(-1.0, 1.0));
        sc.state.mass[1] = std::pow(10.0, rng.range(-1.0, 1.0));
        sc.sys.mass[0] = std::pow(10.0, rng.range(-1.0, 1.0));
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < 2; ++i)
            {
                sc.base[d][i] = rng.range(-1.0, 1.0);
                sc.dv[d][i] = rng.range(-5.0, 5.0);
            }
            sc.sbase[d][0] = rng.range(-1.0, 1.0);
            sc.dvx[d][0] = rng.range(-5.0, 5.0);
            sc.av[d][0] = rng.range(-0.05, 0.05);
        }
        double n[3] = {rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)};
        unit(n);
        const int shape = k % 4; // 0 -- произвольная lambda, 1 -- lambda=0, 2 -- lambda=1, 3 -- точечный контакт
        const double lam = shape == 0 ? rng.uni() : (shape == 1 ? 0.0 : (shape == 2 ? 1.0 : 0.0));
        sc.fs_slots.assign(1, make_fs_slot(0, shape == 3 ? -1 : 1, 0, lam, n, rng.range(0.5, 5.0),
                                           std::pow(10.0, rng.range(-1.0, 1.5))));
        const SolidStage r = apply_solid_stage(exec, sc);
        const auto imp = fs_impulse(sc, r, 0);
        if (!imp.has_force) continue;
        (imp.stick ? n_stick : n_slip) += 1;

        double sum[3] = {0, 0, 0}, sabs = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            const double fa = sc.state.mass[0] * r.fdelta[0].dv[d];
            const double fb = sc.state.mass[1] * r.fdelta[1].dv[d];
            const double fs = sc.sys.mass[0] * r.sdelta[0].dv[d];
            sum[d] = fa + fb + fs;
            sabs += std::fabs(fa) + std::fabs(fb) + std::fabs(fs);
            const double wa = shape == 3 ? 1.0 : 1.0 - lam;
            const double wb = shape == 3 ? 0.0 : lam;
            const double Ftn = norm3(imp.Ft);
            dist_err = std::max(dist_err, std::fabs(fa - wa * imp.Ft[d]) / Ftn);
            dist_err = std::max(dist_err, std::fabs(fb - wb * imp.Ft[d]) / Ftn);
            dist_err = std::max(dist_err, std::fabs(fs + imp.Ft[d]) / Ftn);
        }
        mom_err = std::max(mom_err, norm3(sum) / sabs);
        if (shape == 1 || shape == 3) // вес бусины b нулевой: ровно 0
            zero_ok = zero_ok && r.fdelta[1].dv[0] == 0.0 && r.fdelta[1].dv[1] == 0.0 && r.fdelta[1].dv[2] == 0.0;
        if (shape == 2)
            zero_ok = zero_ok && r.fdelta[0].dv[0] == 0.0 && r.fdelta[0].dv[1] == 0.0 && r.fdelta[0].dv[2] == 0.0;
    }
    mom_err_out = mom_err;
    std::printf("    stick %d, slip %d; max|sum m a|/sum|m a| = %.3e, max откл. распределения = %.3e\n", n_stick, n_slip,
                mom_err, dist_err);
    check(n_stick > 50 && n_slip > 50, "выборка содержит и прилипание, и скольжение");
    check(mom_err <= 1e-14, "m_a a_a + m_b a_b + m_s a_s == 0 (500 случайных контактов, все формы)");
    check(dist_err <= 1e-14, "вклады: a -- (1-lambda)F_t, b -- lambda F_t, Solid -- -F_t (точечный: всё на a)");
    check(zero_ok, "бусина с нулевым весом (lambda=0 -> b, lambda=1 -> a, точечный -> b) не получает вклада (ровно 0)");
}

// ----------------------------------------------------------------------
// FS.8 Диссипация: 1000 случайных контактов, KE (2 бусины + частица Solid)
// после трения не больше, чем без него.
// ----------------------------------------------------------------------
void test_fs_dissipation(double& max_growth_out)
{
    std::printf("test_fs_dissipation:\n");
    rhs::CpuExecutor exec(2);
    Rng rng(8600002);
    int bad = 0, n_stick = 0, n_slip = 0;
    double max_growth = -1.0e300;
    for (int k = 0; k < 1000; ++k)
    {
        SolidScene sc;
        const double h = rng.range(1e-3, 1e-1);
        sc.init(2, 1, h, h, true);
        sc.state.mass[0] = std::pow(10.0, rng.range(-1.0, 1.0));
        sc.state.mass[1] = std::pow(10.0, rng.range(-1.0, 1.0));
        sc.sys.mass[0] = std::pow(10.0, rng.range(-1.0, 1.0));
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < 2; ++i)
            {
                sc.base[d][i] = rng.range(-1.0, 1.0);
                sc.dv[d][i] = rng.range(-5.0, 5.0);
            }
            sc.sbase[d][0] = rng.range(-1.0, 1.0);
            sc.dvx[d][0] = rng.range(-5.0, 5.0);
            sc.av[d][0] = rng.range(-0.1, 0.1);
        }
        double n[3] = {rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)};
        unit(n);
        const double muN = std::pow(10.0, rng.range(-3.0, 1.5));
        sc.fs_slots.assign(1, make_fs_slot(0, 1, 0, rng.uni(), n, 1.0, muN));
        const SolidStage r = apply_solid_stage(exec, sc);
        const auto imp = fs_impulse(sc, r, 0);
        (imp.stick ? n_stick : n_slip) += 1;

        auto ke = [&](const std::vector<double> (&fv)[dim], const std::vector<double> (&sv)[dim]) {
            double e = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
                e += 0.5 * (sc.state.mass[0] * fv[d][0] * fv[d][0] + sc.state.mass[1] * fv[d][1] * fv[d][1] +
                            sc.sys.mass[0] * sv[d][0] * sv[d][0]);
            return e;
        };
        const double e0 = ke(r.fvtrial, r.svtrial);
        const double e1 = ke(r.fvnew, r.svnew);
        max_growth = std::max(max_growth, (e1 - e0) / e0);
        if (!(e1 <= e0 * (1.0 + 1e-12))) ++bad;
    }
    max_growth_out = max_growth;
    std::printf("    stick %d, slip %d; max (KE_после-KE_проб)/KE_проб = %.3e\n", n_stick, n_slip, max_growth);
    check(n_stick > 100 && n_slip > 100, "выборка содержит и прилипание, и скольжение");
    check(bad == 0, "KE с трением <= KE без трения (отн. допуск 1e-12) во всех 1000 контактах");
}

// ----------------------------------------------------------------------
// FS.9 mu_fs = 0: слот kernel'а неактивен, стадия не эмитирует вкладов, dv и
// dvx побитово не изменились (стадия в боевом цикле при этом вообще не
// вызывается -- has_active == false); слот с mu=0 не даёт силы и при принудительной
// активации.
// ----------------------------------------------------------------------
void test_fs_zero_mu()
{
    std::printf("test_fs_zero_mu:\n");
    rhs::CpuExecutor exec(2);
    KernelScene ks;
    build_kernel_scene(ks, exec, 0.0);
    check(!ks.has_active, "mu_fs=0: compute_fiber_solid_contact_contributions возвращает has_active == false");
    check(ks.slots.size() == ks.prims.size() && !ks.slots[0].active, "mu_fs=0: слот на каждый primitive, все неактивны");
    check(ks.fd[0].dv[2] != 0.0 || ks.sd[0].dv[2] != 0.0, "самопроверка: нормальный контакт при этом есть");

    SolidScene sc;
    solid_scene_from_kernel(sc, ks, 1.0e-7);
    std::vector<double> dv_before[dim];
    std::vector<double> dvx_before[dim];
    for (std::size_t d = 0; d < dim; ++d)
    {
        dv_before[d] = sc.dv[d];
        dvx_before[d] = sc.dvx[d];
    }
    const SolidStage r = apply_solid_stage(exec, sc);
    std::size_t femit = 0;
    {
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fb;
        std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sb;
        interactions::compute_friction_contributions(exec, sc.state, sc.fupd, sc.dv, sc.ff_slots, sc.sys, sc.supd,
                                                     sc.dvx, sc.fs_slots, r.bcnt, r.scnt, false, fb, sb);
        for (const auto& b : fb) femit += b.size();
    }
    check(femit == 0 && r.solid_emitted == 0, "неактивные слоты: стадия ничего не эмитирует (ни нить, ни Solid)");
    bool bitwise = true;
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (std::size_t i = 0; i < 2; ++i) sc.dv[d][i] += r.fdelta[i].dv[d];
        sc.dvx[d][0] += r.sdelta[0].dv[d];
        bitwise = bitwise && std::memcmp(dv_before[d].data(), sc.dv[d].data(), 2 * sizeof(double)) == 0 &&
                  std::memcmp(dvx_before[d].data(), sc.dvx[d].data(), 1 * sizeof(double)) == 0;
    }
    check(bitwise, "dv и dvx после стадии побитово равны значениям до неё");
    check(r.bcnt[0] == 0.0 && r.bcnt[1] == 0.0 && r.scnt[0] == 0.0, "c_k и c_s == 0 (неактивные слоты не считаются)");

    const auto dg = interactions::compute_friction_diagnostics(exec, sc.state, sc.fupd, sc.dv, sc.ff_slots, sc.sys,
                                                               sc.supd, sc.dvx, sc.fs_slots, r.bcnt, r.scnt, false);
    check(dg.active_count == 0 && dg.stick_count == 0 && dg.slip_count == 0 && dg.shared_bead_count == 0 &&
              dg.shared_solid_count == 0,
          "диагностика на неактивных слотах: нули");

    // Принудительно активный слот с mu = 0: mu*N == 0 -> силы нет.
    SolidScene s2;
    build_fs_scene(s2, 1.0e-3, 2.0, 0.0);
    const SolidStage r2 = apply_solid_stage(exec, s2);
    const auto imp2 = fs_impulse(s2, r2, 0);
    bool zero = true;
    for (std::size_t d = 0; d < dim; ++d) zero = zero && r2.fdelta[0].dv[d] == 0.0 && r2.fdelta[1].dv[d] == 0.0 &&
                                                  r2.sdelta[0].dv[d] == 0.0;
    check(!imp2.has_force && zero && r2.solid_emitted == 0, "слот с mu=0 (принудительно активный): силы и вкладов нет");
}

// ----------------------------------------------------------------------
// FS.10 Кинематическое тело: solid_kinematic=true -- вклад в Solid не
// эмитируется, c_s/m_s = 0; SolidVelocityUpdate{base = скорость пластины,
// h = 0, av = nullptr} -- пробная скорость Solid == base ТОЧНО (dvx игнорируется);
// прилипание точное относительно заданной скорости Solid; результат совпадает с
// пределом бесконечной массы Solid.
// ----------------------------------------------------------------------
void test_fs_kinematic()
{
    std::printf("test_fs_kinematic:\n");
    rhs::CpuExecutor exec(2);
    const double vplate[3] = {0.7, -0.2, 0.1};

    SolidScene sc;
    build_fs_scene(sc, 1.0e-3, 2.0, 1.0e6, 0.4, nullptr, /*with_av=*/false);
    sc.supd.h = 0.0; // пластина не интегрируется
    for (std::size_t d = 0; d < dim; ++d) sc.sbase[d][0] = vplate[d];
    // dvx (ускорение Solid без трения) остаётся большим: при h=0 оно не должно влиять.
    const SolidStage r = apply_solid_stage(exec, sc, /*kinematic=*/true);
    bool trial_exact = true;
    for (std::size_t d = 0; d < dim; ++d) trial_exact = trial_exact && r.svtrial[d][0] == vplate[d];
    check(trial_exact, "пробная скорость Solid == base ровно (h=0, av=nullptr; dvx игнорируется)");
    check(r.solid_emitted == 0, "solid_kinematic: вклады в Solid (type 13) не эмитируются");
    bool sd_zero = true;
    for (std::size_t d = 0; d < dim; ++d) sd_zero = sd_zero && r.sdelta[0].dv[d] == 0.0;
    check(sd_zero, "solid_kinematic: суммарный вклад в Solid == 0");

    const auto imp = fs_impulse(sc, r, 0, /*kinematic=*/true);
    const double Wfib = interactions::fiber_contact_side_inv_mass(sc.state, 0, 1, sc.fs_slots[0].lambda);
    check(imp.W == Wfib || std::fabs(imp.W - Wfib) <= 1e-15 * Wfib,
          "W = 1/m_eff только нити (c_s/m_s = 0, бесконечная масса) -- c_k=1 у одиночного контакта");
    double ut0[3], ut1[3];
    const double n0 = fs_ut(r.fvtrial, r.svtrial, sc.fs_slots[0], ut0);
    const double n1 = fs_ut(r.fvnew, r.svnew, sc.fs_slots[0], ut1);
    std::printf("    |u_t| до %.6e, после %.3e\n", n0, n1);
    check(n0 > 1e-3 && imp.stick, "исходная касательная скорость ненулевая, ветка прилипания");
    check(n1 <= 1e-12 * n0, "прилипание точное относительно заданной скорости Solid");

    // Предел бесконечной массы: Solid с m_s = 1e15 и той же скоростью даёт ту же F_t.
    SolidScene heavy;
    const double mh[3] = {sc.state.mass[0], sc.state.mass[1], 1.0e15};
    build_fs_scene(heavy, 1.0e-3, 2.0, 1.0e6, 0.4, mh, false);
    heavy.supd.h = 0.0;
    for (std::size_t d = 0; d < dim; ++d) heavy.sbase[d][0] = vplate[d];
    const SolidStage rh = apply_solid_stage(exec, heavy, /*kinematic=*/false);
    const auto imph = fs_impulse(heavy, rh, 0, false);
    double dFt[3] = {imp.Ft[0] - imph.Ft[0], imp.Ft[1] - imph.Ft[1], imp.Ft[2] - imph.Ft[2]};
    check(norm3(dFt) <= 1e-9 * norm3(imp.Ft), "solid_kinematic совпадает с пределом m_s -> inf (F_t, отн. 1e-9)");

    // Конечная масса даёт меньший импульс прилипания (W больше) -- флаг реально меняет W.
    SolidScene fin;
    build_fs_scene(fin, 1.0e-3, 2.0, 1.0e6, 0.4, nullptr, false);
    fin.supd.h = 0.0;
    for (std::size_t d = 0; d < dim; ++d) fin.sbase[d][0] = vplate[d];
    const SolidStage rf = apply_solid_stage(exec, fin, false);
    const auto impf = fs_impulse(fin, rf, 0, false);
    check(impf.W > imp.W && norm3(impf.Ft) < norm3(imp.Ft),
          "с конечной массой Solid: W больше, F_t меньше, чем у кинематического");

    // Три контакта на одной неподвижной пластине (mass splitting по бусинам, c_s/m_s = 0).
    SolidScene m3;
    m3.init(2, 3, 1.0e-3, 0.0, false);
    m3.state.mass[0] = 1.0;
    m3.state.mass[1] = 1.4;
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t i = 0; i < 2; ++i) m3.base[d][i] = (d == 0) ? 1.0 : (d == 1 ? 0.1 : 0.0);
    const double nz[3] = {0.0, 0.0, 1.0};
    const double lams[3] = {0.2, 0.5, 0.8};
    m3.fs_slots.clear();
    for (std::size_t k = 0; k < 3; ++k) m3.fs_slots.push_back(make_fs_slot(0, 1, k, lams[k], nz, 2.0, 1.0e6));
    const SolidStage r3 = apply_solid_stage(exec, m3, true);
    bool no_flip = true, dec = true;
    double worst_along = 1.0;
    for (const auto& c : m3.fs_slots)
    {
        double a0[3], a1[3];
        const double p0 = fs_ut(r3.fvtrial, r3.svtrial, c, a0);
        const double p1 = fs_ut(r3.fvnew, r3.svnew, c, a1);
        const double along = dot3(a1, a0) / p0;
        worst_along = std::min(worst_along, along / p0);
        if (along < -1e-12 * p0) no_flip = false;
        if (!(p1 < p0)) dec = false;
    }
    std::printf("    три контакта на неподвижной пластине: c_k = %.0f, %.0f; min(along/|u_t до|) = %.4f\n", r3.bcnt[0],
                r3.bcnt[1], worst_along);
    check(r3.bcnt[0] == 3.0 && r3.bcnt[1] == 3.0, "c_a = c_b = 3 (три контакта на бусину)");
    check(no_flip && dec, "три контакта на неподвижной пластине: нет перескока, |u_t| уменьшилась");
}

// ----------------------------------------------------------------------
// Подсчёт c_k / c_s (mass splitting, §5.2): оба вида контактов, нулевые веса,
// неактивные слоты; и диагностика единой стадии.
// ----------------------------------------------------------------------
void test_fs_counts_and_diagnostics()
{
    std::printf("test_fs_counts_and_diagnostics:\n");
    rhs::CpuExecutor exec(2);
    Rng rng(4243);
    SolidScene sc;
    sc.init(6, 3, 1.0e-3, 1.0e-3, true);
    for (std::size_t i = 0; i < 6; ++i) sc.state.mass[i] = rng.range(0.5, 4.0);
    for (std::size_t i = 0; i < 3; ++i) sc.sys.mass[i] = rng.range(0.5, 4.0);
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (std::size_t i = 0; i < 6; ++i)
        {
            sc.base[d][i] = rng.range(-1.0, 1.0);
            sc.dv[d][i] = rng.range(-3.0, 3.0);
        }
        for (std::size_t i = 0; i < 3; ++i)
        {
            sc.sbase[d][i] = rng.range(-1.0, 1.0);
            sc.dvx[d][i] = rng.range(-3.0, 3.0);
            sc.av[d][i] = rng.range(-0.05, 0.05);
        }
    }
    auto rnd_n = [&](double n[3]) {
        n[0] = rng.range(-1, 1);
        n[1] = rng.range(-1, 1);
        n[2] = rng.range(-1, 1);
        unit(n);
    };
    auto make_ff = [&](std::size_t a0, long long a1, std::size_t b0, long long b1, double s, double t, double N, double mu,
                       bool active) {
        interactions::FiberFrictionContact c;
        c.active = active;
        c.a0 = a0;
        c.a1 = a1;
        c.b0 = b0;
        c.b1 = b1;
        c.s = s;
        c.t = t;
        rnd_n(c.n);
        c.N = N;
        c.mu = mu;
        return c;
    };
    // нить-нить: f0 -- обычный, f1 -- неактивный, f2 -- t=1 (вес b0 нулевой)
    sc.ff_slots.push_back(make_ff(0, 1, 2, 3, 0.5, 0.5, 2.0, 1.0e6, true));
    sc.ff_slots.push_back(make_ff(0, 1, 2, 3, 0.5, 0.5, 2.0, 0.5, false));
    sc.ff_slots.push_back(make_ff(1, -1, 4, 5, 0.0, 1.0, 2.0, 0.05, true));
    // Fiber-Solid: s0 (lambda=0: вес b нулевой), s1 (lambda=1: вес a нулевой), s2 (точечный, N=0), s3 неактивный
    double n[3];
    rnd_n(n);
    sc.fs_slots.push_back(make_fs_slot(0, 1, 0, 0.0, n, 2.0, 1.0e6));
    rnd_n(n);
    sc.fs_slots.push_back(make_fs_slot(1, 2, 0, 1.0, n, 2.0, 0.05));
    rnd_n(n);
    sc.fs_slots.push_back(make_fs_slot(3, -1, 2, 0.0, n, 0.0, 0.5));
    rnd_n(n);
    sc.fs_slots.push_back(make_fs_slot(3, -1, 1, 0.0, n, 2.0, 0.5));
    sc.fs_slots.back().active = false;

    const SolidStage r = apply_solid_stage(exec, sc);
    const double ec[6] = {2, 2, 2, 2, 0, 1};
    const double es[3] = {2, 0, 1};
    std::printf("    c_k = ");
    for (double v : r.bcnt) std::printf("%.0f ", v);
    std::printf("; c_s = ");
    for (double v : r.scnt) std::printf("%.0f ", v);
    std::printf("\n");
    bool ck_ok = r.bcnt.size() == 6, cs_ok = r.scnt.size() == 3;
    for (int i = 0; i < 6 && ck_ok; ++i) ck_ok = r.bcnt[i] == ec[i];
    for (int i = 0; i < 3 && cs_ok; ++i) cs_ok = r.scnt[i] == es[i];
    check(ck_ok, "c_k по обоим видам: нулевые веса и неактивные слоты не считаются (ожидается 2 2 2 2 0 1)");
    check(cs_ok, "c_s: частица 0 -- два слота, частица 1 -- ни одного (неактивный не считается), частица 2 -- один");

    const auto dg = interactions::compute_friction_diagnostics(exec, sc.state, sc.fupd, sc.dv, sc.ff_slots, sc.sys,
                                                               sc.supd, sc.dvx, sc.fs_slots, r.bcnt, r.scnt, false);
    // независимая пересборка stick/slip по слотам (те же c_k/c_s)
    std::size_t st = 0, sl = 0, act = 0;
    double max_ratio = 0.0;
    for (const auto& c : sc.ff_slots)
    {
        if (!c.active) continue;
        ++act;
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(sc.state, sc.fupd, sc.dv, c, imp, r.bcnt.data());
        if (imp.has_force) (imp.stick ? st : sl) += 1;
        if (imp.has_force) max_ratio = std::max(max_ratio, norm3(imp.Ft) / (c.mu * c.N));
    }
    for (std::size_t k = 0; k < sc.fs_slots.size(); ++k)
    {
        if (!sc.fs_slots[k].active) continue;
        ++act;
        const auto imp = fs_impulse(sc, r, k);
        if (imp.has_force) (imp.stick ? st : sl) += 1;
        if (imp.has_force) max_ratio = std::max(max_ratio, norm3(imp.Ft) / (sc.fs_slots[k].mu * sc.fs_slots[k].N));
    }
    std::printf("    active %zu, stick %zu, slip %zu, shared_bead %zu, shared_solid %zu, max_force_ratio %.15f\n",
                dg.active_count, dg.stick_count, dg.slip_count, dg.shared_bead_count, dg.shared_solid_count,
                dg.max_force_ratio);
    check(dg.active_count == 5 && act == 5, "active_count == 5 (2 нить-нить + 3 Fiber-Solid; неактивные не считаются)");
    check(dg.stick_count == st && dg.slip_count == sl, "stick_count/slip_count совпадают с пересборкой по слотам");
    check(st >= 1 && sl >= 1 && st + sl == 4, "в сцене есть и прилипание, и скольжение; контакт с N=0 без силы");
    check(dg.shared_bead_count == 4, "shared_bead_count == 4 (бусины 0..3 с c_k=2)");
    check(dg.shared_solid_count == 1, "shared_solid_count == 1 (частица Solid 0)");
    check(dg.max_force_ratio <= 1.0 + 1e-14 && std::fabs(dg.max_force_ratio - max_ratio) <= 1e-14,
          "max_force_ratio <= 1 и равен независимой пересборке");
}

// ----------------------------------------------------------------------
// Mass splitting, Fiber-Solid (§5.2, §7.1).
// (1) Сегмент из двух бусин на трёх частицах Solid (lambda = 0.2, 0.5, 0.8),
//     общее скольжение, все контакты в режиме прилипания: с c_k после шага
//     касательная скорость КАЖДОЙ точки контакта не сменила знак и уменьшилась;
//     без mass splitting (c = 1) -- перескок (самопроверка сцены).
// (2) Бусина одновременно в контакте нить-нить и нить-Solid: перескока нет.
// (3) Одиночные контакты (c = 1): побитово формула задачи 3.
// ----------------------------------------------------------------------
void test_fs_mass_splitting()
{
    std::printf("test_fs_mass_splitting:\n");
    rhs::CpuExecutor exec(2);

    // ---- (1) три частицы Solid под одним сегментом ----
    SolidScene sc;
    sc.init(2, 3, 1.0e-3, 1.0e-3, false);
    sc.state.mass[0] = 1.0;
    sc.state.mass[1] = 1.4;
    const double ms[3] = {50.0, 80.0, 120.0};
    for (std::size_t k = 0; k < 3; ++k) sc.sys.mass[k] = ms[k];
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t i = 0; i < 2; ++i) sc.base[d][i] = (d == 0) ? 1.0 : (d == 1 ? 0.1 : 0.0);
    const double nz[3] = {0.0, 0.0, 1.0};
    const double lams[3] = {0.2, 0.5, 0.8};
    for (std::size_t k = 0; k < 3; ++k) sc.fs_slots.push_back(make_fs_slot(0, 1, k, lams[k], nz, 2.0, 1.0e6));

    const SolidStage r = apply_solid_stage(exec, sc);
    check(r.bcnt[0] == 3.0 && r.bcnt[1] == 3.0 && r.scnt[0] == 1.0 && r.scnt[1] == 1.0 && r.scnt[2] == 1.0,
          "c_a = c_b = 3 (три контакта на бусину), c_s = 1 у каждой частицы Solid");
    bool no_flip = true, decreased = true;
    double min_along = 1e300, max_along = -1e300;
    for (const auto& c : sc.fs_slots)
    {
        double a0[3], a1[3];
        const double p0 = fs_ut(r.fvtrial, r.svtrial, c, a0);
        const double p1 = fs_ut(r.fvnew, r.svnew, c, a1);
        const double along = dot3(a1, a0) / p0; // компонента u_t после вдоль исходной u_t
        min_along = std::min(min_along, along / p0);
        max_along = std::max(max_along, along / p0);
        if (along < -1e-12 * p0) no_flip = false;
        if (!(p1 < p0 && along < p0)) decreased = false;
    }
    std::printf("    с mass splitting: u_t после / u_t до вдоль исходной: [%.4f, %.4f]\n", min_along, max_along);
    check(no_flip, "с mass splitting: касательная скорость каждой точки контакта не сменила знак (нет перескока)");
    check(decreased, "с mass splitting: касательная скорость каждой точки контакта уменьшилась");

    // самопроверка: c = 1 (без mass splitting) -- перескок
    const SolidStage r1 = apply_solid_stage(exec, sc, false, /*no_split=*/true);
    int flipped = 0;
    double min_along1 = 1e300;
    for (const auto& c : sc.fs_slots)
    {
        double a0[3], a1[3];
        const double p0 = fs_ut(r1.fvtrial, r1.svtrial, c, a0);
        fs_ut(r1.fvnew, r1.svnew, c, a1);
        const double along = dot3(a1, a0) / p0;
        min_along1 = std::min(min_along1, along / p0);
        if (along < 0.0) ++flipped;
    }
    std::printf("    без mass splitting (c=1): min u_t после / u_t до вдоль исходной = %.4f, перескок в %d из 3 контактов\n",
                min_along1, flipped);
    check(flipped == 3, "самопроверка: при c = 1 в этой сцене перескок (знак u_t сменился) во всех трёх контактах");

    // ---- (1b) сегмент на реальной пластине: слоты kernel'а ----
    {
        rhs::CpuExecutor exec1(1);
        Materials mats;
        const double dxp = 1.5e-3, rho_s = 2700.0, d_f = 1.0e-3, rho_f = 1200.0, l0 = 2.0e-3;
        const int nbeads = 4;
        fiber::FiberState state;
        fiber::FiberTopology topo;
        state.resize(nbeads);
        topo.resize(nbeads);
        for (int i = 0; i < nbeads; ++i)
        {
            state.x[0][i] = i * l0;
            state.x[1][i] = 0.0;
            state.x[2][i] = 0.9e-3;
            for (std::size_t d = 0; d < dim; ++d) state.vx[d][i] = 0.0;
            state.vx[0][i] = 0.3;
            state.mass[i] = rho_f * pi * d_f * d_f / 4.0 * l0;
            state.object_id[i] = 0;
            state.material_id[i] = 0;
            if (i > 0) topo.prev[i] = i - 1;
            if (i + 1 < nbeads)
            {
                topo.next[i] = i + 1;
                // связь нужна для радиуса поиска (L_max берётся по topology.bonds)
                fiber::FiberBond bond;
                bond.i = i;
                bond.j = i + 1;
                bond.l0 = l0;
                bond.l_break = l0 * 1.0e9;
                topo.bonds.push_back(bond);
            }
        }
        std::vector<fiber::FiberMaterial> materials(1);
        materials[0].E = 2.0e7;
        materials[0].nu = 0.3;
        materials[0].d = d_f;
        materials[0].rho = rho_f;
        materials[0].eps_break = 1.0;
        materials[0].mu_fs = 1.0e6;
        // пластина: колонны x = (ix+0.5)*dxp, два ряда y = +-0.5*dxp (как в e2e), z = 0
        ParticleSystem plate;
        std::size_t np = 0;
        const int nx = 8;
        plate.resize(static_cast<std::size_t>(nx) * 2);
        for (int ix = -2; ix < nx - 2; ++ix)
            for (int iy = -1; iy <= 0; ++iy)
            {
                plate.inzializ_index(np);
                plate.x[0][np] = (ix + 0.5) * dxp;
                plate.x[1][np] = (iy + 0.5) * dxp;
                plate.x[2][np] = 0.0;
                plate.itype[np] = 4;
                plate.body[np] = 4;
                plate.rho[np] = rho_s;
                plate.rho_0[np] = rho_s;
                plate.mass[np] = rho_s * dxp * dxp * dxp;
                plate.hsml[np] = 1.0;
                ++np;
            }
        interactions::FiberSolidContactSearchResult search;
        interactions::build_fiber_solid_contact_search(exec1, mats, plate, np, state, topo, materials, search);
        std::vector<interactions::FiberSolidFrictionContact> slots;
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
        std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
        const bool has = interactions::compute_fiber_solid_contact_contributions(
            exec1, mats, plate, state, materials, search.primitives, fbuf, sbuf, slots);
        const auto fd = rhs::reduce_contributions(exec1, fbuf, state.size());
        const auto sd = rhs::reduce_contributions(exec1, sbuf, np);

        SolidScene ps;
        ps.init(nbeads, np, 1.0e-7, 1.0e-7, false);
        for (int i = 0; i < nbeads; ++i) ps.state.mass[i] = state.mass[i];
        for (std::size_t i = 0; i < np; ++i) ps.sys.mass[i] = plate.mass[i];
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (int i = 0; i < nbeads; ++i)
            {
                ps.base[d][i] = state.vx[d][i];
                ps.dv[d][i] = fd[i].dv[d];
            }
            for (std::size_t i = 0; i < np; ++i) ps.dvx[d][i] = sd[i].dv[d];
        }
        for (const auto& s : slots)
            if (s.active) ps.fs_slots.push_back(s);
        check(has && ps.fs_slots.size() >= 6, "плита: kernel дал не менее 6 активных контактов нить-Solid (получено " +
                                                   std::to_string(ps.fs_slots.size()) + ")");

        const SolidStage rs = apply_solid_stage(exec1, ps);
        double cmin = 1e300, cmax = 0.0, csmax = 0.0;
        for (int i = 0; i < nbeads; ++i)
        {
            cmin = std::min(cmin, rs.bcnt[i]);
            cmax = std::max(cmax, rs.bcnt[i]);
        }
        for (std::size_t i = 0; i < np; ++i) csmax = std::max(csmax, rs.scnt[i]);
        bool nf = true, dc = true, small_ut = true;
        double mina = 1e300, maxa = -1e300;
        for (const auto& c : ps.fs_slots)
        {
            double a0[3], a1[3];
            const double p0 = fs_ut(rs.fvtrial, rs.svtrial, c, a0);
            const double p1 = fs_ut(rs.fvnew, rs.svnew, c, a1);
            const double along = dot3(a1, a0) / p0;
            mina = std::min(mina, along / p0);
            maxa = std::max(maxa, along / p0);
            if (p0 < 1e-3) small_ut = false;
            if (along < -1e-12 * p0) nf = false;
            if (!(along < p0)) dc = false;
            (void)p1;
        }
        std::printf("    плита: %zu контактов, c_k в [%.0f, %.0f], max c_s = %.0f; с mass splitting: u_t после/до "
                    "вдоль исходной в [%.4f, %.4f]\n",
                    ps.fs_slots.size(), cmin, cmax, csmax, mina, maxa);
        check(small_ut, "плита: исходная касательная скорость во всех контактах ненулевая");
        check(cmax >= 4.0, "плита: бусина участвует не менее чем в 4 контактах (c_max >= 4)");
        check(nf, "плита с mass splitting: нет перескока в каждой точке контакта");
        check(dc, "плита с mass splitting: касательная скорость вдоль исходной уменьшилась в каждой точке контакта");
        const SolidStage rs1 = apply_solid_stage(exec1, ps, false, true);
        int fl = 0;
        double mina1 = 1e300;
        for (const auto& c : ps.fs_slots)
        {
            double a0[3], a1[3];
            const double p0 = fs_ut(rs1.fvtrial, rs1.svtrial, c, a0);
            fs_ut(rs1.fvnew, rs1.svnew, c, a1);
            const double along = dot3(a1, a0) / p0;
            mina1 = std::min(mina1, along / p0);
            if (along < 0.0) ++fl;
        }
        std::printf("    плита без mass splitting (c=1): min u_t после/до вдоль исходной = %.4f, перескок в %d из %zu\n",
                    mina1, fl, ps.fs_slots.size());
        check(fl >= 1, "самопроверка плиты: без mass splitting (c = 1) перескок есть хотя бы в одном контакте");
    }

    // ---- (2) бусина одновременно в контакте нить-нить и нить-Solid ----
    // Бусина 1 (масса 1) входит и в контакт с бусиной 2 другой нити, и в контакт
    // с частицей Solid; обе стороны-партнёры покоятся. Массы партнёров m_2 и m_s:
    // при m_2 == m_s оба контакта с c_1 = 2 останавливают скольжение ТОЧНО (сумма
    // двух долей ровно 1); при m_2 != m_s один из контактов перекомпенсирует на
    // доли процента (контакт не видит движения партнёра другого контакта) --
    // перескока "вдвое", как при c = 1, нет.
    for (int variant = 0; variant < 2; ++variant)
    {
        const double m2 = 50.0, msol = (variant == 0) ? 50.0 : 80.0;
        SolidScene s2;
        s2.init(4, 1, 1.0e-3, 1.0e-3, false);
        s2.state.mass[0] = 1.0;
        s2.state.mass[1] = 1.0;
        s2.state.mass[2] = m2;
        s2.state.mass[3] = 1.0;
        s2.sys.mass[0] = msol;
        s2.base[0][1] = 1.0; // бусина 1 скользит вдоль x и y относительно бусины 2 нити B и частицы Solid
        s2.base[1][1] = 0.3;
        interactions::FiberFrictionContact ff;
        ff.active = true;
        ff.a0 = 0;
        ff.a1 = 1;
        ff.s = 1.0; // вес бусины 0 нулевой, вся сторона A -- бусина 1
        ff.b0 = 2;
        ff.b1 = -1;
        ff.t = 0.0;
        ff.n[2] = 1.0;
        ff.N = 2.0;
        ff.mu = 1.0e6;
        s2.ff_slots.assign(1, ff);
        s2.fs_slots.assign(1, make_fs_slot(1, -1, 0, 0.0, nz, 2.0, 1.0e6));
        const std::string tag = variant == 0 ? "m_2 == m_s" : "m_2 != m_s";
        const SolidStage rr = apply_solid_stage(exec, s2);
        check(rr.bcnt[1] == 2.0 && rr.bcnt[0] == 0.0 && rr.bcnt[2] == 1.0 && rr.scnt[0] == 1.0,
              tag + ": c_1 = 2 (нить-нить + нить-Solid), c_0 = 0 (нулевой вес), c_2 = 1, c_s = 1");
        double u_ff0[3], u_fs0[3], u_ff1[3], u_fs1[3];
        // касательная скорость нить-нить: бусина 1 относительно бусины 2 (вес A -- бусина 1)
        auto ff_ut = [&](const std::vector<double> (&v)[dim], double ut[3]) {
            double u[3];
            for (std::size_t d = 0; d < dim; ++d) u[d] = v[d][1] - v[d][2];
            const double un = dot3(u, ff.n);
            for (std::size_t d = 0; d < dim; ++d) ut[d] = u[d] - un * ff.n[d];
            return norm3(ut);
        };
        const double f0 = ff_ut(rr.fvtrial, u_ff0), f1 = ff_ut(rr.fvnew, u_ff1);
        const double g0 = fs_ut(rr.fvtrial, rr.svtrial, s2.fs_slots[0], u_fs0);
        const double g1 = fs_ut(rr.fvnew, rr.svnew, s2.fs_slots[0], u_fs1);
        const double aff = dot3(u_ff1, u_ff0) / f0, afs = dot3(u_fs1, u_fs0) / g0;
        std::printf("    %s: нить-нить |u_t| %.4f -> вдоль исходной %.4e (%.2e от исходной); нить-Solid %.4f -> %.4e (%.2e)\n",
                    tag.c_str(), f0, aff, aff / f0, g0, afs, afs / g0);
        if (variant == 0)
        {
            check(aff >= -1e-12 * f0 && afs >= -1e-12 * g0, tag + ": с mass splitting: перескока нет ни в одном из двух контактов");
            check(f1 <= 1e-12 * f0 && g1 <= 1e-12 * g0, tag + ": с mass splitting: прилипание точное в обоих контактах");
        }
        else
        {
            // Mass splitting по счёту контактов точен только при равных массах партнёров
            // (известное ограничение): при m_2 != m_s остаётся малая перекомпенсация (факт ~0.18 %).
            check(aff >= -0.01 * f0 && afs >= -0.01 * g0,
                  tag + ": с mass splitting: перекомпенсация не более 1% (перескока вдвое нет)");
        }
        const SolidStage r0 = apply_solid_stage(exec, s2, false, true);
        double v_ff0[3], v_fs0[3], v_ff1[3], v_fs1[3];
        ff_ut(r0.fvtrial, v_ff0);
        ff_ut(r0.fvnew, v_ff1);
        fs_ut(r0.fvtrial, r0.svtrial, s2.fs_slots[0], v_fs0);
        fs_ut(r0.fvnew, r0.svnew, s2.fs_slots[0], v_fs1);
        const double bff = dot3(v_ff1, v_ff0) / f0, bfs = dot3(v_fs1, v_fs0) / g0;
        std::printf("    %s, без mass splitting (c=1): нить-нить %.4f, нить-Solid %.4f (вдоль исходной, |u_t до| = %.4f / %.4f)\n",
                    tag.c_str(), bff, bfs, f0, g0);
        check(bff < 0.0 && bfs < 0.0, tag + ": самопроверка: при c = 1 обе относительные касательные скорости сменили знак");
    }

    // ---- (3) одиночные контакты: c = 1 -- побитово формула задачи 3 ----
    {
        SolidScene s3;
        build_fs_scene(s3, 1.0e-3, 2.0, 1.0e6);
        const SolidStage r3 = apply_solid_stage(exec, s3);
        check(r3.bcnt[0] == 1.0 && r3.bcnt[1] == 1.0 && r3.scnt[0] == 1.0, "одиночный контакт: c_a = c_b = c_s = 1");
        const auto with_counts = fs_impulse(s3, r3, 0);
        interactions::FiberFrictionImpulse without;
        interactions::fiber_solid_friction_impulse(s3.state, s3.fupd, s3.dv, s3.sys, s3.supd, s3.dvx, s3.fs_slots[0],
                                                   without);
        // формула задачи 3: W = fiber_contact_side_inv_mass + 1/m_s, проекция без счётчиков
        const auto& c = s3.fs_slots[0];
        const double W3 = interactions::fiber_contact_side_inv_mass(s3.state, 0, 1, c.lambda) + 1.0 / s3.sys.mass[0];
        double u[3];
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double va = s3.base[d][0] + s3.fupd.h * s3.dv[d][0];
            const double vb = s3.base[d][1] + s3.fupd.h * s3.dv[d][1];
            const double vs = (s3.sbase[d][0] + s3.supd.h * s3.dvx[d][0]) + s3.av[d][0];
            u[d] = ((1.0 - c.lambda) * va + c.lambda * vb) - vs;
        }
        interactions::FiberFrictionImpulse ref;
        interactions::friction_projection_impulse(u, c.n, W3, c.mu, c.N, s3.fupd.h, ref);
        auto same = [](const interactions::FiberFrictionImpulse& x, const interactions::FiberFrictionImpulse& y) {
            return x.W == y.W && x.J == y.J && x.ut_norm == y.ut_norm && x.Ft[0] == y.Ft[0] && x.Ft[1] == y.Ft[1] &&
                   x.Ft[2] == y.Ft[2];
        };
        check(same(with_counts, without), "одиночный контакт: результат с c = 1 побитово совпадает с вызовом без счётчиков");
        check(same(with_counts, ref), "одиночный контакт: результат побитово совпадает с формулой задачи 3 (W3, проекция)");
    }
}

// ----------------------------------------------------------------------
// Интерфейс решателя в связке со стадией: реальные PredictorCorrectorSolver
// (первый шаг и шаги n>=2) и EulerSolver, XSPH-поправка av ненулевая. Стадия
// берёт пробные скорости из solver.fiber_velocity_update /
// solver.solid_velocity_update, добавляет a^f в solver.*_rhs, advance() даёт
// скорости, после которых относительная касательная скорость точки контакта
// равна нулю (прилипание точное, одиночный контакт).
// (Побитовость пробной скорости Solid против advance() без стадии --
// tests/test_ode_solver.cpp.)
// ----------------------------------------------------------------------
template <class Solver>
double run_solver_stick(const char* label, int nsteps)
{
    rhs::CpuExecutor exec(2);
    Rng rng(777);
    const std::size_t nf = 2, ns = 1;
    fiber::FiberState fb;
    fb.resize(nf);
    fb.mass[0] = 1.0;
    fb.mass[1] = 3.0;
    for (std::size_t i = 0; i < nf; ++i)
    {
        fb.object_id[i] = 0;
        fb.material_id[i] = 0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            fb.x[d][i] = rng.range(-1, 1);
            fb.vx[d][i] = rng.range(-1, 1);
        }
    }
    ParticleSystem sys;
    sys.resize(ns);
    sys.inzializ_index(0);
    sys.mass[0] = 2.0;
    sys.rho[0] = 1000.0;
    sys.rho_0[0] = 1000.0;
    for (std::size_t d = 0; d < dim; ++d) sys.vx[d][0] = rng.range(-1, 1);

    Solver solver;
    solver.reset(ns, ns, nf);
    double n[3] = {0.3, -0.5, 0.8};
    unit(n);
    std::vector<interactions::FiberSolidFrictionContact> slots(1, make_fs_slot(0, 1, 0, 0.3, n, 2.0, 1.0e6));
    std::vector<interactions::FiberFrictionContact> none;

    double worst = 0.0;
    bool bit_ok = true;
    for (int k = 0; k < nsteps; ++k)
    {
        const double dt = 1.0e-3 * (1.0 + 0.3 * k);
        solver.predict(sys, ns, fb, nf, k == 0 ? dt : 1.0e-3 * (1.0 + 0.3 * (k - 1)));
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < nf; ++i) solver.fiber_rhs.dv[d][i] = rng.range(-5.0, 5.0);
            solver.solid_rhs.dvx[d][0] = rng.range(-3.0, 3.0);
            solver.solid_rhs.av[d][0] = rng.range(-0.02, 0.02); // XSPH (Эйлер его игнорирует)
        }
        const auto fupd = solver.fiber_velocity_update(fb, dt);
        const auto supd = solver.solid_velocity_update(sys, dt);
        std::vector<double> ftrial[dim], strial[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < nf; ++i) ftrial[d].push_back((*fupd.base[d])[i] + fupd.h * solver.fiber_rhs.dv[d][i]);
            double t = (*supd.base[d])[0] + supd.h * solver.solid_rhs.dvx[d][0];
            if (supd.av) t += (*supd.av)[d][0];
            strial[d].push_back(t);
        }
        std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> fcb;
        std::vector<rhs::ContributionBuffer<interactions::SolidCountDelta>> scb;
        std::vector<double> bc, sc_;
        interactions::compute_friction_counts(exec, nf, ns, none, slots, fcb, scb, bc, sc_);
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fb_buf;
        std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sb_buf;
        interactions::compute_friction_contributions(exec, fb, fupd, solver.fiber_rhs.dv, none, sys, supd,
                                                     solver.solid_rhs.dvx, slots, bc, sc_, false, fb_buf, sb_buf);
        const auto fd = rhs::reduce_contributions(exec, fb_buf, nf);
        const auto sd = rhs::reduce_contributions(exec, sb_buf, ns);
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < nf; ++i) solver.fiber_rhs.dv[d][i] += fd[i].dv[d];
            solver.solid_rhs.dvx[d][0] += sd[0].dv[d];
        }
        double ut0[3];
        const double n0 = fs_ut(ftrial, strial, slots[0], ut0);
        solver.advance(sys, ns, fb, nf, dt);
        std::vector<double> fvn[dim], svn[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            fvn[d] = fb.vx[d];
            svn[d] = {sys.vx[d][0]};
        }
        double ut1[3];
        const double n1 = fs_ut(fvn, svn, slots[0], ut1);
        worst = std::max(worst, n1 / n0);
        bit_ok = bit_ok && n0 > 1e-3;
    }
    std::printf("    %s: %d шагов, max |u_t после advance| / |u_t пробная| = %.3e\n", label, nsteps, worst);
    check(bit_ok, std::string(label) + ": пробная касательная скорость ненулевая на всех шагах");
    check(worst <= 1e-12, std::string(label) + ": после advance() прилипание точное (пробная скорость Solid = формула advance)");
    return worst;
}

void test_fs_solver_interface()
{
    std::printf("test_fs_solver_interface:\n");
    run_solver_stick<ode::PredictorCorrectorSolver>("PredictorCorrector (шаг 1, затем n>=2, XSPH)", 4);
    run_solver_stick<ode::EulerSolver>("Euler (без XSPH)", 4);
}

// ============================================================================
// Режим списка активных слотов трения (FrictionActiveRef) против полного прохода: два параллельных
// конвейера (ядра контактов -> [fill из нормальных слотов] -> счётчики c_k/c_s -> вклады трения -> reduce)
// на сцене из сетки нитей и частиц Solid, 5 «шагов» с возмущением положений (множество активных слотов
// меняется), при 1, 2, 7, 8 потоках. Сравнение побитовое: слоты, флаги, c_k/c_s, итоговые a^f.
// ============================================================================

bool bits_eq(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

bool slot_eq(const interactions::FiberFrictionContact& a, const interactions::FiberFrictionContact& b)
{
    bool eq = a.active == b.active && a.a0 == b.a0 && a.a1 == b.a1 && a.b0 == b.b0 && a.b1 == b.b1 &&
              bits_eq(a.s, b.s) && bits_eq(a.t, b.t) && bits_eq(a.N, b.N) && bits_eq(a.mu, b.mu);
    for (std::size_t d = 0; d < dim; ++d) eq = eq && bits_eq(a.n[d], b.n[d]);
    return eq;
}

bool slot_eq(const interactions::FiberSolidFrictionContact& a, const interactions::FiberSolidFrictionContact& b)
{
    bool eq = a.active == b.active && a.a == b.a && a.b == b.b && a.s == b.s && bits_eq(a.lambda, b.lambda) &&
              bits_eq(a.N, b.N) && bits_eq(a.mu, b.mu);
    for (std::size_t d = 0; d < dim; ++d) eq = eq && bits_eq(a.n[d], b.n[d]);
    return eq;
}

template <class Slot>
bool slots_eq(const std::vector<Slot>& a, const std::vector<Slot>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!slot_eq(a[i], b[i])) return false;
    return true;
}

bool vec_bits_eq(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!bits_eq(a[i], b[i])) return false;
    return true;
}

// Сцена: nf x-нитей и nf y-нитей по 12 бусин (пересечения на разной глубине) и частицы Solid рядом.
struct ListScene
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    ParticleSystem sys;
    Materials mats;
    std::vector<double> base[dim], dv[dim], sbase[dim];
    std::vector<std::vector<double>> dvx;
    ode::FiberVelocityUpdate fupd{};
    ode::SolidVelocityUpdate supd{};
    int nf = 6;
    int nb = 12;
    std::size_t np = 14;

    static constexpr double l0 = 0.01;

    void add_chain(std::size_t base_idx, int oid, const std::vector<std::array<double, 3>>& beads)
    {
        const std::size_t n = beads.size();
        for (std::size_t k = 0; k < n; ++k)
        {
            const std::size_t gi = base_idx + k;
            for (int d3 = 0; d3 < 3; ++d3) state.x[d3][gi] = beads[k][d3];
            state.mass[gi] = 1.0e-6 * (1.0 + 0.05 * static_cast<double>(gi % 5));
            state.object_id[gi] = oid;
            state.material_id[gi] = 0;
            if (k > 0) topology.prev[gi] = static_cast<long long>(gi - 1);
            if (k + 1 < n)
            {
                topology.next[gi] = static_cast<long long>(gi + 1);
                fiber::FiberBond b;
                b.i = gi;
                b.j = gi + 1;
                b.l0 = l0;
                b.l_break = l0 * 1e9;
                topology.bonds.push_back(b);
            }
        }
    }

    void init()
    {
        fiber::FiberMaterial m;
        m.E = 2.0e9;
        m.nu = 0.3;
        m.d = 1.0e-3;
        m.rho = 1200.0;
        m.eps_break = 1.0;
        m.mu_ff = 0.3;
        m.mu_fs = 0.2;
        materials.assign(1, m);
        const std::size_t nfib = static_cast<std::size_t>(2 * nf * nb);
        state.resize(nfib);
        topology.resize(nfib);
        std::size_t idx = 0;
        int oid = 0;
        for (int j = 0; j < nf; ++j, idx += static_cast<std::size_t>(nb))
        {
            std::vector<std::array<double, 3>> beads;
            for (int i = 0; i < nb; ++i) beads.push_back({(i - 5.5) * l0, 0.5 * l0 + j * 0.9 * l0, 0.0});
            add_chain(idx, oid++, beads);
        }
        for (int k = 0; k < nf; ++k, idx += static_cast<std::size_t>(nb))
        {
            std::vector<std::array<double, 3>> beads;
            for (int i = 0; i < nb; ++i) beads.push_back({0.5 * l0 + k * 0.9 * l0 - 2.5 * l0, (i - 2.0) * l0, 0.0});
            add_chain(idx, oid++, beads);
        }
        const double dxp = 1.5e-3, rho_s = 2700.0;
        sys.resize(np);
        for (std::size_t i = 0; i < np; ++i)
        {
            sys.inzializ_index(i);
            sys.mass[i] = rho_s * dxp * dxp * dxp;
            sys.rho[i] = rho_s;
            sys.rho_0[i] = rho_s;
            sys.itype[i] = 4;
            sys.hsml[i] = 1.0;
        }
        for (std::size_t d = 0; d < dim; ++d)
        {
            base[d].resize(nfib);
            dv[d].resize(nfib);
            sbase[d].resize(np);
            for (std::size_t i = 0; i < nfib; ++i)
            {
                base[d][i] = 0.3 * std::sin(0.7 * static_cast<double>(i + 1) + static_cast<double>(d));
                dv[d][i] = 50.0 * std::cos(0.37 * static_cast<double>(i) + 1.3 * static_cast<double>(d));
                state.vx[d][i] = base[d][i];
            }
            for (std::size_t i = 0; i < np; ++i)
                sbase[d][i] = 0.2 * std::cos(0.9 * static_cast<double>(i) + static_cast<double>(d));
            fupd.base[d] = &base[d];
            supd.base[d] = &sbase[d];
        }
        dvx.assign(dim, std::vector<double>(np, 0.0));
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < np; ++i) dvx[d][i] = 30.0 * std::sin(0.51 * static_cast<double>(i) + static_cast<double>(d));
        fupd.h = 1.0e-6;
        supd.h = 1.0e-6;
        supd.av = nullptr;
        set_step(0);
    }

    // Положения на «шаге» st: y-нити (object_id >= nf) смещены по z, частицы -- по высоте над x-нитями.
    void set_step(int st)
    {
        static const double shift[6] = {0.0, 0.25e-3, -0.1e-3, 0.5e-3, -0.3e-3, 0.8e-3};
        const double sh = shift[st % 6];
        for (std::size_t i = 0; i < state.size(); ++i)
        {
            const int oid = state.object_id[i];
            if (oid < nf)
                state.x[2][i] = 1.0e-4 * oid;
            else
            {
                const int k = oid - nf;
                state.x[2][i] = 0.95e-3 - 0.7e-4 * k + sh * (1.0 + 0.5 * (k % 3));
            }
        }
        for (std::size_t i = 0; i < np; ++i)
        {
            const int j = static_cast<int>(i % static_cast<std::size_t>(nf));
            sys.x[0][i] = (-0.045 + 0.0071 * static_cast<double>(i));
            sys.x[1][i] = 0.5 * l0 + j * 0.9 * l0 + 0.2e-3 * static_cast<double>((i + static_cast<std::size_t>(st)) % 3);
            sys.x[2][i] = 1.0e-4 * j + (0.4e-3 + 0.2e-3 * static_cast<double>((i + 2 * static_cast<std::size_t>(st)) % 7));
        }
    }
};

struct PipeOut
{
    bool has_ff = false, has_fs = false, stage = false;
    std::vector<interactions::FiberFrictionContact> ff;
    std::vector<interactions::FiberSolidFrictionContact> fs;
    std::vector<double> bcnt, scnt;
    std::vector<fiber::FiberDelta> fdelta;
    std::vector<interactions::SolidDelta> sdelta;
    std::size_t n_active_ff = 0, n_active_fs = 0;
};

struct Pipe
{
    bool list = false;
    bool proj = false;
    std::vector<interactions::FiberFrictionContact> ff_slots;
    std::vector<interactions::FiberSolidFrictionContact> fs_slots;
    std::vector<interactions::FrictionActiveRef> ff_act, fs_act;
    std::vector<rhs::CpuExecutor::SegPart<interactions::FrictionActiveRef>> ff_parts, fs_parts;
    std::vector<interactions::FiberNormalContact> normal_ff;
    std::vector<rhs::CpuExecutor::SegPart<interactions::FiberNormalContact>> normal_parts;
    std::vector<double> bead_dt;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> k_fbuf, k_fbuf2, st_fbuf;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> k_sbuf, st_sbuf;
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> fcb;
    std::vector<rhs::ContributionBuffer<interactions::SolidCountDelta>> scb;
    std::vector<interactions::FiberCountDelta> fcd;
    std::vector<interactions::SolidCountDelta> scd;
    rhs::ReduceScratch scratch;

    PipeOut step(const rhs::CpuExecutor& exec, ListScene& sc, double h)
    {
        PipeOut out;
        interactions::FiberContactSearchResult ffs;
        interactions::build_fiber_contact_search(exec, sc.state, sc.topology, sc.materials, ffs);
        interactions::FiberSolidContactSearchResult fss;
        interactions::build_fiber_solid_contact_search(exec, sc.mats, sc.sys, sc.sys.size(), sc.state, sc.topology,
                                                        sc.materials, fss);
        interactions::FiberContactStepParams pf, ps;
        if (proj)
        {
            // Проекция: ядро нить-нить даёт нормальные слоты вместо силы (закон 1, только геометрия).
            interactions::compute_fiber_bead_dt_ref(exec, sc.state, sc.topology, sc.materials, bead_dt);
            pf.law.type = 1;
            pf.bead_dt_ref = &bead_dt;
            pf.normal_ff = &normal_ff;
            pf.normal_ff_parts = &normal_parts;
        }
        if (list)
        {
            pf.fric_active = &ff_act;
            pf.fric_active_parts = &ff_parts;
            ps.fric_active = &fs_act;
            ps.fric_active_parts = &fs_parts;
        }
        out.has_ff = interactions::compute_fiber_contact_contributions(exec, sc.state, sc.materials, ffs.primitives,
                                                                       k_fbuf, ff_slots, pf);
        out.has_fs = interactions::compute_fiber_solid_contact_contributions(
            exec, sc.mats, sc.sys, sc.state, sc.materials, fss.primitives, k_fbuf2, k_sbuf, fs_slots, ps);
        if (proj)
        {
            // Нормальные импульсы P -- синтетические (часть нулевых: слот трения не создаётся).
            for (std::size_t j = 0; j < normal_ff.size(); ++j)
                normal_ff[j].P = (j % 4 == 0) ? 0.0 : 2.0e-4 * static_cast<double>(1 + j % 3);
            out.has_ff = interactions::fill_friction_slots_from_normal_ff(exec, h, normal_ff, ff_slots,
                                                                         list ? &ff_act : nullptr);
        }
        const std::size_t nf = sc.state.size(), ns = sc.sys.size();
        const auto* ffa = list ? &ff_act : nullptr;
        const auto* fsa = list ? &fs_act : nullptr;
        if (nf > 0 && (out.has_ff || out.has_fs))
        {
            out.stage = true;
            interactions::compute_friction_counts(exec, nf, ns, ff_slots, fs_slots, fcb, scb, out.bcnt, out.scnt,
                                                  fcd, scd, scratch, ffa, fsa);
            interactions::compute_friction_contributions(exec, sc.state, sc.fupd, sc.dv, ff_slots, sc.sys, sc.supd,
                                                         sc.dvx, fs_slots, out.bcnt, out.scnt, false, st_fbuf,
                                                         st_sbuf, ffa, fsa);
            out.fdelta = rhs::reduce_contributions(exec, st_fbuf, nf);
            out.sdelta = rhs::reduce_contributions(exec, st_sbuf, ns);
        }
        out.ff = ff_slots;
        out.fs = fs_slots;
        for (const auto& c : ff_slots) out.n_active_ff += c.active ? 1u : 0u;
        for (const auto& c : fs_slots) out.n_active_fs += c.active ? 1u : 0u;
        if (list)
        {
            // Инварианты списка: возрастание k, число записей == число активных слотов, только активные.
            bool ok = ff_act.size() == out.n_active_ff && fs_act.size() == out.n_active_fs;
            for (std::size_t j = 0; j < ff_act.size(); ++j)
                ok = ok && ff_act[j].k < ff_slots.size() && ff_slots[ff_act[j].k].active && (j == 0 || ff_act[j - 1].k < ff_act[j].k);
            for (std::size_t j = 0; j < fs_act.size(); ++j)
                ok = ok && fs_act[j].k < fs_slots.size() && fs_slots[fs_act[j].k].active && (j == 0 || fs_act[j - 1].k < fs_act[j].k);
            if (!ok) ++failures;
        }
        return out;
    }
};

bool pipe_out_eq(const PipeOut& a, const PipeOut& b)
{
    bool eq = a.has_ff == b.has_ff && a.has_fs == b.has_fs && a.stage == b.stage && slots_eq(a.ff, b.ff) &&
              slots_eq(a.fs, b.fs);
    if (a.stage && b.stage)
    {
        eq = eq && vec_bits_eq(a.bcnt, b.bcnt) && vec_bits_eq(a.scnt, b.scnt) && a.fdelta.size() == b.fdelta.size() &&
             a.sdelta.size() == b.sdelta.size();
        if (eq)
        {
            for (std::size_t i = 0; i < a.fdelta.size(); ++i)
                for (std::size_t d = 0; d < dim; ++d) eq = eq && bits_eq(a.fdelta[i].dv[d], b.fdelta[i].dv[d]);
            for (std::size_t i = 0; i < a.sdelta.size(); ++i)
                for (std::size_t d = 0; d < dim; ++d) eq = eq && bits_eq(a.sdelta[i].dv[d], b.sdelta[i].dv[d]);
        }
    }
    return eq;
}

void test_active_list_vs_full(bool projection)
{
    std::printf("test_active_list_vs_full (%s):\n", projection ? "проекция нить-нить" : "силовой путь");
    for (int nt : {1, 2, 7, 8})
    {
        rhs::CpuExecutor exec(nt);
        ListScene sc;
        sc.init();
        Pipe full, lst;
        full.proj = lst.proj = projection;
        lst.list = true;
        bool all_eq = true;
        std::size_t min_active = static_cast<std::size_t>(-1), max_active = 0;
        bool changed = false;
        std::size_t prev_active = 0;
        for (int st = 0; st < 6; ++st)
        {
            sc.set_step(st);
            const PipeOut a = full.step(exec, sc, 1.0e-6);
            const PipeOut b = lst.step(exec, sc, 1.0e-6);
            if (!pipe_out_eq(a, b)) all_eq = false;
            const std::size_t act = a.n_active_ff + a.n_active_fs;
            min_active = std::min(min_active, act);
            max_active = std::max(max_active, act);
            if (st > 0 && act != prev_active) changed = true;
            prev_active = act;
        }
        const std::string lab = "потоков " + std::to_string(nt) + ": ";
        std::printf("    %sактивных слотов за 6 шагов: min=%zu, max=%zu\n", lab.c_str(), min_active, max_active);
        check(all_eq, lab + "слоты, флаги, c_k/c_s и a^f в режиме списка побитово равны полному проходу");
        check(max_active > 0, lab + "в сценах есть активные слоты (предпосылка)");
        check(changed, lab + "множество активных слотов меняется между шагами (предпосылка)");
    }
}

// compute_fiber_bead_dt_ref с executor побитово равен прежней последовательной версии.
void test_bead_dt_ref_exec()
{
    std::printf("test_bead_dt_ref_exec:\n");
    ListScene sc;
    sc.init();
    sc.materials.push_back(sc.materials[0]);
    sc.materials[1].E = 7.0e8;
    sc.materials[1].d = 0.8e-3;
    for (std::size_t i = 0; i < sc.state.size(); ++i) sc.state.material_id[i] = static_cast<int>((i / 12) % 2);
    std::vector<double> ref;
    interactions::compute_fiber_bead_dt_ref(sc.state, sc.topology, sc.materials, ref);
    bool eq = true;
    for (int nt : {1, 2, 7, 8})
    {
        rhs::CpuExecutor exec(nt);
        std::vector<double> got;
        interactions::compute_fiber_bead_dt_ref(exec, sc.state, sc.topology, sc.materials, got);
        eq = eq && vec_bits_eq(ref, got);
    }
    check(eq && !ref.empty(), "compute_fiber_bead_dt_ref(exec, ...) == прежняя версия (побитово, 1/2/7/8 потоков)");
}

} // namespace

int main()
{
    double e81 = 0, e82 = 0, e83 = 0, growth86 = 0, e87a = 0, e87b = 0;
    test_8_1_stick(e81);
    test_8_2_slip(e82);
    test_8_3_static_under_force(e83);
    test_8_4_tangential();
    test_8_5_momentum_distribution();
    test_8_6_dissipation(growth86);
    test_8_7_dt_independence(e87a, e87b);
    test_8_8_zero_mu();
    test_diagnostics();
    test_mass_splitting();

    double fs_stick = 0, fs_slip = 0, fs_static = 0, fs_xsph = 0, fs_mom = 0, fs_growth = 0;
    test_fs_stick(fs_stick);
    test_fs_slip(fs_slip);
    test_fs_static_under_force(fs_static);
    test_fs_xsph(fs_xsph);
    test_fs_signs();
    test_fs_tangential();
    test_fs_momentum_distribution(fs_mom);
    test_fs_dissipation(fs_growth);
    test_fs_zero_mu();
    test_fs_kinematic();
    test_fs_counts_and_diagnostics();
    test_fs_mass_splitting();
    test_fs_solver_interface();
    test_active_list_vs_full(false);
    test_active_list_vs_full(true);
    test_bead_dt_ref_exec();

    std::printf("ключевые числа: 8.1 |u_t|после/до = %.3e; 8.2 |F_t|/muN-1 = %.3e; 8.3 |u_t|после/до = %.3e; "
                "8.6 max прирост KE = %.3e; 8.7 = %.3e / %.3e\n",
                e81, e82, e83, growth86, e87a, e87b);

    std::printf("ключевые числа Fiber-Solid: FS.1 |u_t|после/до = %.3e; FS.2 |F_t|/muN-1 = %.3e; FS.3 = %.3e; "
                "FS.4 (XSPH) = %.3e; FS.7 max|sum m a|/sum|m a| = %.3e; FS.8 max прирост KE = %.3e\n",
                fs_stick, fs_slip, fs_static, fs_xsph, fs_mom, fs_growth);

    if (failures == 0)
        std::printf("test_fiber_friction: OK\n");
    else
        std::printf("test_fiber_friction: %d failure(s)\n", failures);
    return failures;
}
