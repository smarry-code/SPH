// tests/test_fiber_plasticity_integration.cpp
// Этап 2 дизайна пластичности нити (notes/fiber-model/fiber-plasticity-design.md, раздел
// 10, п. 2): поля состояния в FiberBond, RHS-ядро с bond_stress и проход
// Postprocessing apply_fiber_bond_constitutive. Тесты раздела 8:
//   Т0   упругий материал даёт ровно прежние силы и dt-лимит (побитово), проход
//        Postprocessing упругую связь не трогает (межверсионное побитовое
//        сравнение целых прогонов -- cases/guo2019/experiment/fiber_t0_dump.cpp, см. отчёт);
//   Т8   затухающие колебания пары бусин с Кельвином: ω_d = ω·sqrt(1-ξ²),
//        огибающая e^{-ξωt}, логарифмический декремент;
//   Т9   предел шага при демпфере (пара и цепочка, точная формула
//        критического шага схемы проекта, см. комментарий к тесту);
//   Т10  энергетический баланс цикла нагрузка-разгрузка по отдельным ветвям
//        (пластика, Максвелл, Кельвин, всё вместе), W_p >= 0, W_v >= 0; остаток
//        баланса объяснён и разложен (O(dt) + потеря на односторонней отсечке);
//   Т11  сохранение импульса и момента импульса с пластикой и демпфером;
//   начальное eps_prev (init_fiber_bond_state);
//   БОЕВОЙ ЦИКЛ: time_integration(cfg, fiber::Model&) на малой сцене с
//   неупругой нитью: проход apply_fiber_bond_constitutive вызывается (внутри
//   run_fiber_topology_postprocessing), eps_p/w_p/w_v ненулевые и монотонные,
//   состояние выживших связей переживает разрыв;
//   Т14  критерии разрыва (eps_break, FAILURE.EPS_P, FAILURE.SIGMA_U, в т.ч. без
//        пластики), topology.version, сохранение состояния уцелевших связей;
//   плюс полноцикловая проверка пластики (остаточное удлинение, сила 0 при
//   провисании).
//
// Динамика -- локальный драйвер PredictorCorrector (как FiberPcIntegrator в
// test_fiber_bond.cpp, формулы src/time_integration.cpp), к которому добавлены
// тот же порядок шага, что в боевом цикле, и вызов прохода Postprocessing после
// advance. Кинематически заданное движение конца (driven) реализовано как в
// clamp: ускорение бусины обнуляется, скорость задаётся; реакция
// m·a_связи записывается, работа внешней силы = -реакция·v·dt.
// Собирается с SPH_DIM=3 (CMakeLists.txt), по образцу test_fiber_bond.
#include "fiber/fiber_plasticity.hpp"
#include "fiber/fiber_system.hpp"
#include "fiber/model.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/postprocessing.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "run_guard.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"

#include "e2e_harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_plasticity_integration требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

void info(const std::string& what) { std::printf("  info : %s\n", what.c_str()); }

std::string fmt(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.6g", v);
    return buf;
}

double area_of(double d) { return pi * d * d / 4.0; }

fiber::FiberMaterial make_material(double E, double d, double rho)
{
    fiber::FiberMaterial m;
    m.E = E;
    m.nu = 0.3;
    m.d = d;
    m.rho = rho;
    m.eps_break = 1.0e9; // разрыв в этих тестах не нужен
    return m;
}

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// --- Цепочка бусин ----------------------------------------------------------
struct Chain
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

// N бусин; координаты задаёт вызывающий (x[0..2][i]), связи i--i+1 с длиной l0.
Chain make_chain(std::size_t N, double l0, const fiber::FiberMaterial& mat, double mass)
{
    Chain c;
    c.materials.push_back(mat);
    c.state.resize(N);
    c.topology.resize(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d)
        {
            c.state.x[d][i] = 0.0;
            c.state.vx[d][i] = 0.0;
        }
        c.state.x[0][i] = static_cast<double>(i) * l0;
        c.state.mass[i] = mass;
        c.state.object_id[i] = 0;
        c.state.material_id[i] = 0;
        if (i + 1 < N)
        {
            c.topology.next[i] = static_cast<long long>(i + 1);
            c.topology.prev[i + 1] = static_cast<long long>(i);
            c.topology.bonds.push_back(fiber::FiberBond{i, i + 1, l0, l0 * 1.0e9});
        }
    }
    return c;
}

// --- Драйвер шага -----------------------------------------------------------
// Порядок шага как в боевом цикле: предиктор (полушаг по прошлому ускорению) ->
// RHS (Fiber-Bond) -> advance (корректор, x += dt·v) -> проход Postprocessing.
struct Driver
{
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
    std::vector<double> dv_prev[dim];
    std::vector<double> v_min[dim];
    std::vector<double> ext[dim];  // постоянное внешнее ускорение (предварительное натяжение)
    std::vector<char> driven;      // бусина задана кинематически вдоль x
    std::vector<double> vdrive;    // её скорость на текущем шаге
    std::vector<double> react_x;   // сила связей на заданную бусину (x-компонента), последний шаг
    bool first = true;
    bool constitutive = true;
    double work_ext = 0.0;         // работа внешних сил над заданными бусинами [Дж]
    // Диагностика односторонности связи: энергия, "потерянная" отсечкой
    // σ = max(σ_модели, 0) (ветви связи сжимают σ_модели < 0 при ε > ε_p --
    // быстрая разгрузка Максвелла/Кельвина -- а сила вдоль связи не бывает
    // отрицательной). Внутренние ветви уже отдали (запасённую) энергию, а бусинам
    // она не досталась: Σ σ_модели·Δε·V по этим шагам [Дж]. Не входит в W_p, W_v.
    double clip_loss = 0.0;

    void ensure(std::size_t n)
    {
        if (dv_prev[0].size() == n) return;
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_prev[d].assign(n, 0.0);
            v_min[d].assign(n, 0.0);
            if (ext[d].size() != n) ext[d].assign(n, 0.0);
        }
        if (driven.size() != n) driven.assign(n, 0);
        if (vdrive.size() != n) vdrive.assign(n, 0.0);
        react_x.assign(n, 0.0);
    }

    void step(fiber::FiberState& fb, fiber::FiberTopology& topo, const std::vector<fiber::FiberMaterial>& mats,
              double dt)
    {
        const std::size_t n = fb.size();
        ensure(n);
        // Как в начале run_time_integration: eps_prev = текущая ε (связь могла
        // быть создана при l != l0).
        if (first && constitutive) interactions::init_fiber_bond_state(fb, topo, mats);
        for (std::size_t i = 0; i < n; ++i)
            if (driven[i])
            {
                fb.vx[0][i] = vdrive[i];
                fb.vx[1][i] = 0.0;
                fb.vx[2][i] = 0.0;
            }
        if (!first)
        {
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i)
                {
                    v_min[d][i] = fb.vx[d][i];
                    fb.vx[d][i] += 0.5 * dt * dv_prev[d][i];
                }
        }
        std::vector<double> sig_neg(topo.bonds.size(), 0.0), eps_before(topo.bonds.size(), 0.0);
        for (std::size_t kb = 0; kb < topo.bonds.size(); ++kb)
        {
            const fiber::FiberBond& b = topo.bonds[kb];
            const fiber::FiberMaterial& mat = mats[static_cast<std::size_t>(fb.material_id[b.i])];
            if (!mat.inelastic()) continue;
            const double A = interactions::fiber_cross_section_area(mat);
            const double l = interactions::fiber_bond_length(fb, b);
            const double eps = l / b.l0 - 1.0;
            eps_before[kb] = eps;
            const fiber::ViscoParams vp = interactions::fiber_bond_visco(mat, A, b.l0, fb.mass[b.i], fb.mass[b.j]);
            double edot = 0.0;
            if (vp.eta != 0.0)
            {
                double vrel = 0.0;
                for (std::size_t d = 0; d < dim; ++d) vrel += (fb.vx[d][b.j] - fb.vx[d][b.i]) * (fb.x[d][b.j] - fb.x[d][b.i]);
                edot = (vrel / l) / b.l0;
            }
            const double sig = mat.E * (eps - b.st.eps_p) + b.st.q + vp.eta * edot;
            if (eps > b.st.eps_p && sig < 0.0) sig_neg[kb] = sig;
        }
        const std::vector<fiber::FiberDelta> delta = interactions::compute_fiber_bond(exec, fb, topo, mats, buffers);
        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t d = 0; d < dim; ++d) dv_prev[d][i] = delta[i].dv[d] + ext[d][i];
            if (driven[i])
            {
                react_x[i] = fb.mass[i] * delta[i].dv[0];
                work_ext += -react_x[i] * vdrive[i] * dt;
                for (std::size_t d = 0; d < dim; ++d) dv_prev[d][i] = 0.0;
            }
        }
        const double coeff = first ? dt * 0.5 : dt;
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < n; ++i)
            {
                const double base = first ? fb.vx[d][i] : v_min[d][i];
                fb.vx[d][i] = base + coeff * dv_prev[d][i];
                fb.x[d][i] += dt * fb.vx[d][i];
            }
        first = false;
        if (constitutive) interactions::apply_fiber_bond_constitutive(exec, fb, topo, mats, dt);
        for (std::size_t kb = 0; kb < topo.bonds.size(); ++kb)
            if (sig_neg[kb] != 0.0)
            {
                const fiber::FiberBond& b = topo.bonds[kb];
                const fiber::FiberMaterial& mat = mats[static_cast<std::size_t>(fb.material_id[b.i])];
                const double V = interactions::fiber_cross_section_area(mat) * b.l0;
                clip_loss += sig_neg[kb] * (interactions::fiber_bond_length(fb, b) / b.l0 - 1.0 - eps_before[kb]) * V;
            }
    }
};

struct Energies
{
    double kin = 0.0, el = 0.0, wp = 0.0, wv = 0.0;
    double total() const { return kin + el + wp + wv; }
};

// E_kin + E_упр (пружина E в натянутом состоянии + пружина Максвелла) + накопленные W_p, W_v [Дж].
Energies energies(const fiber::FiberState& fb, const fiber::FiberTopology& topo,
                  const std::vector<fiber::FiberMaterial>& mats)
{
    Energies e;
    for (std::size_t i = 0; i < fb.size(); ++i)
    {
        double v2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d) v2 += fb.vx[d][i] * fb.vx[d][i];
        e.kin += 0.5 * fb.mass[i] * v2;
    }
    for (const fiber::FiberBond& b : topo.bonds)
    {
        const fiber::FiberMaterial& mat = mats[static_cast<std::size_t>(fb.material_id[b.i])];
        const double V = interactions::fiber_cross_section_area(mat) * b.l0;
        const double l = interactions::fiber_bond_length(fb, b);
        const double eps = l / b.l0 - 1.0;
        double w = 0.0;
        if (eps > b.st.eps_p) w += 0.5 * mat.E * (eps - b.st.eps_p) * (eps - b.st.eps_p);
        if (mat.visco.E1 > 0.0 && mat.visco.tau > 0.0) w += 0.5 * b.st.q * b.st.q / mat.visco.E1;
        e.el += w * V;
        e.wp += b.st.w_p * V;
        e.wv += b.st.w_v * V;
    }
    return e;
}

// ----------------------------------------------------------------------
// Т0. Упругий материал: ровно прежняя формула и прежний dt-лимит; проход
// Postprocessing упругую связь не меняет.
// ----------------------------------------------------------------------
void test_t0_elastic_bitwise()
{
    std::printf("T0 (упругий материал побитово прежний):\n");
    const double l0 = 1.0e-3, E = 2.0e9, d = 1.0e-4, rho = 1140.0;
    const double A = area_of(d);
    const std::size_t N = 9;
    fiber::FiberMaterial mat{E, 0.35, d, rho, 0.05};
    Chain c = make_chain(N, l0, mat, rho * A * l0);
    // Неравномерные сжатие/растяжение и 3D-излом, чтобы покрыть обе ветви ядра.
    for (std::size_t i = 0; i < N; ++i)
    {
        const double s = 1.0 + 0.013 * static_cast<double>(static_cast<int>(i % 4) - 1);
        c.state.x[0][i] = (i == 0 ? 0.0 : c.state.x[0][i - 1] + l0 * s);
        c.state.x[1][i] = 2.0e-5 * static_cast<double>(i % 3);
        c.state.x[2][i] = -1.0e-5 * static_cast<double>(i % 2);
        c.state.mass[i] = rho * A * l0 * (1.0 + 0.2 * static_cast<double>(i % 3));
        c.state.vx[0][i] = 3.0 * static_cast<double>(i) - 7.0; // скорость не должна влиять на упругую силу
        c.state.vx[1][i] = (i % 2) ? 1.5 : -2.0;
    }
    rhs::CpuExecutor exec(2);
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
    const std::vector<fiber::FiberDelta> got = interactions::compute_fiber_bond(exec, c.state, c.topology, c.materials, buffers);

    // Эталон: дословно прежняя формула ядра.
    std::vector<fiber::FiberDelta> ref(N);
    for (const fiber::FiberBond& b : c.topology.bonds)
    {
        double xij[dim];
        double l2 = 0.0;
        for (std::size_t dd = 0; dd < dim; ++dd)
        {
            xij[dd] = c.state.x[dd][b.j] - c.state.x[dd][b.i];
            l2 += xij[dd] * xij[dd];
        }
        const double l = std::sqrt(l2);
        if (l > b.l0 && l > 0.0)
        {
            const double strain = l / b.l0 - 1.0;
            const double coeff = mat.E * A * strain / l;
            for (std::size_t dd = 0; dd < dim; ++dd)
            {
                const double f = coeff * xij[dd];
                ref[b.i].dv[dd] += f / c.state.mass[b.i];
                ref[b.j].dv[dd] += -f / c.state.mass[b.j];
            }
        }
    }
    bool same = true;
    bool nonzero = false;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t dd = 0; dd < dim; ++dd)
        {
            if (!bits_equal(got[i].dv[dd], ref[i].dv[dd])) same = false;
            if (ref[i].dv[dd] != 0.0) nonzero = true;
        }
    check(nonzero, "самопроверка: в сцене есть ненулевые силы");
    check(same, "ускорения бусин совпадают с прежней формулой побитово");

    // dt-лимит: дословно прежняя формула.
    const double dt_got = interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
    double dt_ref = std::numeric_limits<double>::max();
    for (const fiber::FiberBond& b : c.topology.bonds)
    {
        const double k_eff = mat.E * A / b.l0;
        const double mi = c.state.mass[b.i], mj = c.state.mass[b.j];
        const double m_red = (mi * mj) / (mi + mj);
        dt_ref = std::min(dt_ref, 0.8 / std::sqrt(k_eff / m_red));
    }
    check(bits_equal(dt_got, dt_ref), "dt-лимит связей совпадает с прежней формулой побитово");

    // Проход Postprocessing: состояние упругих связей остаётся нулевым.
    for (int s = 0; s < 3; ++s) interactions::apply_fiber_bond_constitutive(exec, c.state, c.topology, c.materials, 1.0e-8);
    bool zero_state = true;
    for (const fiber::FiberBond& b : c.topology.bonds)
        if (b.st.eps_p != 0.0 || b.st.q != 0.0 || b.st.eps_prev != 0.0 || b.st.w_p != 0.0 || b.st.w_v != 0.0)
            zero_state = false;
    check(zero_state, "проход apply_fiber_bond_constitutive не меняет состояние упругих связей");
    check(!mat.inelastic(), "FiberMaterial по умолчанию упругий (inelastic() == false)");

    // Агрегатная инициализация старого вида компилируется и даёт выключенные блоки.
    fiber::FiberMaterial m7{E, 0.35, d, rho, 0.05, 0.1, 0.2};
    check(!m7.plastic.enabled() && m7.visco.eta == 0.0 && m7.visco.xi == 0.0 && m7.visco.E1 == 0.0 && m7.rate.model == fiber::RateModel::None,
          "FiberMaterial{E,nu,d,rho,eps_break,mu_ff,mu_fs}: неупругие блоки выключены");

    // Пластика включена, но связь ещё в упругой зоне (eps_p = 0, q = 0, eta = 0):
    // сила совпадает с упругой с точностью до округления (другой порядок операций).
    fiber::FiberMaterial pm = mat;
    pm.plastic.sy = 1.0e12; // недостижимый предел
    pm.plastic.H = 0.0;
    Chain cp = c;
    cp.materials[0] = pm;
    const std::vector<fiber::FiberDelta> gp = interactions::compute_fiber_bond(exec, cp.state, cp.topology, cp.materials, buffers);
    double worst = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t dd = 0; dd < dim; ++dd)
        {
            worst = std::max(worst, std::fabs(gp[i].dv[dd] - ref[i].dv[dd]));
            scale = std::max(scale, std::fabs(ref[i].dv[dd]));
        }
    check(worst <= 1.0e-13 * scale, "пластический путь в упругой зоне совпадает с упругим силой до округления (" + fmt(worst / scale) + ")");
}

// ----------------------------------------------------------------------
// Полноцикловая проверка пластики: пара бусин, заданное движение конца.
// ----------------------------------------------------------------------
void test_plastic_full_cycle()
{
    std::printf("Полный цикл пластики (пара, заданное движение конца):\n");
    const double l0 = 1.0e-3, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double sy = 60.0e6, H = 0.4e9;
    const double A = area_of(d);
    fiber::FiberMaterial mat = make_material(E, d, rho);
    mat.plastic.sy = sy;
    mat.plastic.H = H;
    Chain c = make_chain(2, l0, mat, rho * A * l0);

    Driver drv;
    drv.ensure(2);
    drv.driven[0] = 1; // закреплена (v = 0)
    drv.driven[1] = 1; // заданное движение
    const double eps_max = 0.08;
    const double V = 0.5; // м/с
    const double dt = 1.0e-8;
    const double L_max = l0 * (1.0 + eps_max);

    // Нагрузка до eps_max.
    drv.vdrive[1] = V;
    double F_last = 0.0;
    while (c.state.x[0][1] - c.state.x[0][0] < L_max)
    {
        drv.step(c.state, c.topology, c.materials, dt);
        F_last = -drv.react_x[1]; // тянущая сила, приложенная снаружи
    }
    const double eps_p_exp = (eps_max - sy / E) * E / (E + H);
    const fiber::FiberBond& b = c.topology.bonds[0];
    const double eps_end = (c.state.x[0][1] - c.state.x[0][0]) / l0 - 1.0;
    const double eps_p_exp_end = (eps_end - sy / E) * E / (E + H);
    check(std::fabs(b.st.eps_p - eps_p_exp_end) < 1.0e-12, "ε_p после нагрузки = (ε - σ_y0/E)·E/(E+H) (" + fmt(b.st.eps_p) + " против " + fmt(eps_p_exp_end) + "), ε_max = " + fmt(eps_p_exp));
    // Сила последнего шага считалась по x^n и ε_p^n: A·σ_s(ε_p^n) (σ_ep лежит на поверхности текучести).
    check(std::fabs(F_last - A * (sy + H * (b.st.eps_p - 0.0))) < 3.0 * A * E * V * dt / l0 + 1e-12 * A * sy,
          "сила в пластике = A·(σ_y0 + H·ε_p) с точностью до одного шага (" + fmt(F_last) + " Н)");
    check(b.st.w_p > 0.0, "пластическая работа w_p > 0");

    // Разгрузка: двигаем конец назад, ищем момент обнуления силы.
    drv.vdrive[1] = -V;
    const double l_p = l0 * (1.0 + b.st.eps_p);
    double l_zero = -1.0;
    double F_prev = -1.0; // первый шаг разгрузки считает силу в точке поворота (наклон там иной)
    bool slope_ok = true;
    int guard = 0;
    while (c.state.x[0][1] - c.state.x[0][0] > l0 * 0.9 && guard++ < 1000000)
    {
        drv.step(c.state, c.topology, c.materials, dt);
        const double F = -drv.react_x[1];
        const double l_now = c.state.x[0][1] - c.state.x[0][0];
        if (l_zero < 0.0 && F == 0.0) l_zero = l_now + V * dt; // длина на момент расчёта силы (x^n)
        if (l_zero < 0.0 && F_prev > 0.0 && F > 0.0)
        {
            // наклон разгрузки = E (с учётом сдвига шага ε_p неизменна)
            const double slope = (F_prev - F) / (A * E * V * dt / l0);
            if (std::fabs(slope - 1.0) > 1.0e-6) slope_ok = false;
        }
        if (l_zero > 0.0 && F != 0.0) slope_ok = false; // после провисания сила обязана быть нулевой
        F_prev = F;
    }
    check(slope_ok, "наклон разгрузки равен E, после провисания сила строго 0 (до сжатия l = 0.9·l0)");
    check(l_zero > 0.0 && std::fabs(l_zero - l_p) <= 1.5 * V * dt,
          "сила обнуляется при l = l0(1+ε_p) (" + fmt(l_zero) + " против " + fmt(l_p) + ", шаг " + fmt(V * dt) + ")");
    check(l_zero > l0 * 1.03, "остаточное удлинение заметно: l_p/l0 - 1 = " + fmt(b.st.eps_p));

    // Повторная нагрузка: сила появляется только за l_p, текучесть на достигнутом σ_y.
    drv.vdrive[1] = V;
    const double eps_p_before = b.st.eps_p;
    double first_F_l = -1.0;
    while (c.state.x[0][1] - c.state.x[0][0] < L_max)
    {
        drv.step(c.state, c.topology, c.materials, dt);
        if (first_F_l < 0.0 && drv.react_x[1] != 0.0) first_F_l = c.state.x[0][1] - c.state.x[0][0];
    }
    check(first_F_l > 0.0 && std::fabs(first_F_l - l_p) <= 2.5 * V * dt, "при повторной нагрузке сила появляется за l0(1+ε_p)");
    check(std::fabs(b.st.eps_p - eps_p_before) < 1.0e-5, "повторная нагрузка до прежнего максимума ε_p не растит (запомнено упрочнение): Δε_p = " + fmt(b.st.eps_p - eps_p_before));
}

// ----------------------------------------------------------------------
// Т8. Свободные затухающие колебания пары бусин с демпфером Кельвина.
// Постоянная внешняя сила сдвигает равновесие в зону натяжения (s_eq > 0), и
// колебания s(t) = s_eq + a·e^{-ξωt}(cos ω_d t + ξω/ω_d sin ω_d t) идут целиком
// в натянутой фазе, где закон ЛИНЕЕН (см. вывод в test_fiber_bond.cpp).
// ----------------------------------------------------------------------
struct PairOsc
{
    double omega = 0.0, omega_d_meas = 0.0, delta_meas = 0.0, slope_meas = 0.0;
    int peaks = 0;
};

PairOsc run_damped_pair(fiber::FiberMaterial mat, double l0, double m, double dt_over_T)
{
    PairOsc out;
    const double A = area_of(mat.d);
    const double E_inst = interactions::fiber_instant_modulus(mat);
    const double k = E_inst * A / l0;
    const double m_red = m / 2.0;
    const double omega = std::sqrt(k / m_red);
    out.omega = omega;
    const double T = 2.0 * pi / omega;
    const double s_eq = 0.5 * l0;
    const double a0 = 0.3 * s_eq;
    Chain c = make_chain(2, l0, mat, m);
    c.state.x[0][1] = l0 + s_eq + a0;
    Driver drv;
    drv.ensure(2);
    const double g = k * s_eq; // сила тяги, уравновешивающая s_eq
    drv.ext[0][1] = g / m;
    drv.ext[0][0] = -g / m;
    const double dt = T * dt_over_T;
    const std::size_t nsteps = static_cast<std::size_t>(std::round(8.0 * T / dt));
    std::vector<double> t(nsteps + 1), y(nsteps + 1);
    t[0] = 0.0;
    y[0] = a0;
    for (std::size_t n = 1; n <= nsteps; ++n)
    {
        drv.step(c.state, c.topology, c.materials, dt);
        t[n] = static_cast<double>(n) * dt;
        y[n] = (c.state.x[0][1] - c.state.x[0][0]) - l0 - s_eq;
    }
    // Максимумы y: t = 0 (старт из покоя -- максимум) и локальные с параболической интерполяцией.
    std::vector<double> tp{0.0}, ap{a0};
    for (std::size_t n = 1; n + 1 <= nsteps; ++n)
    {
        if (y[n] > y[n - 1] && y[n] >= y[n + 1] && y[n] > 0.0)
        {
            const double den = y[n - 1] - 2.0 * y[n] + y[n + 1];
            const double off = den != 0.0 ? 0.5 * (y[n - 1] - y[n + 1]) / den : 0.0;
            tp.push_back(t[n] + off * dt);
            ap.push_back(y[n] - 0.25 * (y[n - 1] - y[n + 1]) * off);
        }
    }
    out.peaks = static_cast<int>(tp.size());
    if (tp.size() >= 4)
    {
        const std::size_t K = 3;
        const double Td = (tp[K] - tp[0]) / static_cast<double>(K);
        out.omega_d_meas = 2.0 * pi / Td;
        out.delta_meas = std::log(ap[0] / ap[K]) / static_cast<double>(K);
        // наклон ln(a_k) по t_k (линейная регрессия) = -ξω
        double st = 0, sy = 0, stt = 0, sty = 0;
        for (std::size_t kk = 0; kk <= K; ++kk)
        {
            st += tp[kk];
            sy += std::log(ap[kk]);
            stt += tp[kk] * tp[kk];
            sty += tp[kk] * std::log(ap[kk]);
        }
        const double nn = static_cast<double>(K + 1);
        out.slope_meas = (nn * sty - st * sy) / (nn * stt - st * st);
    }
    return out;
}

void test_t8_damped_oscillation()
{
    std::printf("Т8 (затухающие колебания пары с Кельвином):\n");
    const double l0 = 0.01, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double m = 1.0e-6;
    const double A = area_of(d);
    for (const double xi : {0.1, 0.2, 0.3})
    {
        // (а) демпфирование задано через xi
        fiber::FiberMaterial mat = make_material(E, d, rho);
        mat.visco.xi = xi;
        const PairOsc r = run_damped_pair(mat, l0, m, 1.0 / 4000.0);
        const double wd_exact = r.omega * std::sqrt(1.0 - xi * xi);
        const double delta_exact = 2.0 * pi * xi / std::sqrt(1.0 - xi * xi);
        info("ξ = " + fmt(xi) + ": ω = " + fmt(r.omega) + ", ω_d = " + fmt(r.omega_d_meas) + " (точное " + fmt(wd_exact) +
             "), δ = " + fmt(r.delta_meas) + " (точное " + fmt(delta_exact) + "), наклон огибающей " + fmt(r.slope_meas) +
             " (точное " + fmt(-xi * r.omega) + "), максимумов " + std::to_string(r.peaks));
        check(r.peaks >= 4, "ξ = " + fmt(xi) + ": найдено не менее 4 максимумов");
        check(std::fabs(r.omega_d_meas / wd_exact - 1.0) < 2.0e-3, "ξ = " + fmt(xi) + ": ω_d = ω·sqrt(1-ξ²) (отн. ошибка " + fmt(r.omega_d_meas / wd_exact - 1.0) + ")");
        check(std::fabs(r.delta_meas / delta_exact - 1.0) < 1.0e-2, "ξ = " + fmt(xi) + ": логарифмический декремент 2πξ/sqrt(1-ξ²) (отн. ошибка " + fmt(r.delta_meas / delta_exact - 1.0) + ")");
        check(std::fabs(r.slope_meas / (-xi * r.omega) - 1.0) < 1.0e-2, "ξ = " + fmt(xi) + ": огибающая ∝ e^{-ξωt} (отн. ошибка наклона " + fmt(r.slope_meas / (-xi * r.omega) - 1.0) + ")");

        // (б) то же демпфирование, заданное вязкостью eta [Па·с]: c = η A/l0 = 2ξ sqrt(k m_red)
        fiber::FiberMaterial mat_eta = make_material(E, d, rho);
        const double k = E * A / l0;
        mat_eta.visco.eta = 2.0 * xi * std::sqrt(k * m / 2.0) * l0 / A;
        const PairOsc re = run_damped_pair(mat_eta, l0, m, 1.0 / 4000.0);
        check(std::fabs(re.omega_d_meas / r.omega_d_meas - 1.0) < 1.0e-9 && std::fabs(re.delta_meas / r.delta_meas - 1.0) < 1.0e-9,
              "ξ = " + fmt(xi) + ": задание через eta и через xi дают один и тот же результат");
    }
}

// ----------------------------------------------------------------------
// Т9. Предел шага при демпфере. Линейная пара/цепочка под предварительным
// натяжением, малое возмущение δ0 = 1e-8·l0; шаг k·dt_bond,
// dt_bond = compute_fiber_bond_dt_limit. Рост = max|s - s_eq| за прогон / δ0.
//
// Критический шаг схемы проекта (линейная часть, спектральный радиус 3x3-матрицы
// шага «предиктор скорости -> RHS -> corrector»: v_pred = v + h/2·a_prev,
// a = -ω²s - 2ξω·v_pred, v' = v + h·a, s' = s + h·v'; h = ω·dt) выходит на 1
// через -1 и в замкнутой форме равен
//     h_кр(ξ) = 2(sqrt(1+4ξ²) - 2ξ)
// (формула Белычко с УДВОЕННЫМ ξ: демпфер берёт запаздывающую скорость
// предиктора): 2.00 (ξ=0), 0.829 (0.5), 0.472 (1), 0.246 (2), 0.0998 (5),
// h_кр·ξ -> 1/2. Формулу проверяют прогоны ниже, а числа сверены с прямым
// расчётом спектрального радиуса матрицы.
//
// Старшая мода ЦЕПОЧКИ в sqrt2 выше по частоте, чем у пары, и (демпфер
// пропорционален жёсткости) по ξ: ω_цепи = sqrt2·ω_пары, ξ_цепи = sqrt2·ξ_пары.
// Прежний лимит (Белычко для пары, ξ-множитель sqrt(1+ξ²)-ξ) при ξ ≳ 1.2 давал
// неустойчивую цепочку даже при 0.9·dt_lim. Принят
//     dt_bond = min( 0.8/ω,  0.8·h_кр(sqrt2·ξ)/(sqrt2·ω) ),
// множитель 0.8 -- тот же запас, что в упругом лимите.
// Результаты (growth цепочки из 12 бусин, коэффициенты k от dt_bond):
//  * цепочка с половинными концевыми массами (старшая мода строго чередующаяся,
//    лимит задаёт концевая связь с m_red = m/3, поэтому запас лишний в sqrt(3/2)
//    раз): устойчиво до 1.5·dt_bond, неустойчиво с 1.7·dt_bond (рост ~1e7);
//  * цепочка с одинаковыми массами: устойчиво до 1.25·dt_bond (1.25 = 1/0.8 --
//    предел формулы), неустойчиво при 1.5·dt_bond (рост ~1e7..1e8). Здесь лимит
//    почти точный (не чрезмерно консервативен);
//  * пара: устойчиво при 0.9 и 2.0·dt_bond, неустойчиво при 3.0·dt_bond.
// 0.9·dt_bond устойчиво везде, ξ = 0.5 ... 5.
// ----------------------------------------------------------------------
double pair_growth(double xi, double dt_factor, std::size_t nsteps, double* dt_bond_out = nullptr)
{
    const double l0 = 0.01, E = 2.0e9, d = 1.0e-3, rho = 1200.0, m = 1.0e-6;
    const double A = area_of(d);
    fiber::FiberMaterial mat = make_material(E, d, rho);
    mat.visco.xi = xi;
    const double k = E * A / l0;
    const double s_eq = 0.5 * l0;
    const double delta0 = 1.0e-8 * l0;
    Chain c = make_chain(2, l0, mat, m);
    c.state.x[0][1] = l0 + s_eq + delta0;
    Driver drv;
    drv.ensure(2);
    drv.ext[0][1] = k * s_eq / m;
    drv.ext[0][0] = -k * s_eq / m;
    rhs::CpuExecutor exec(1);
    const double dt_bond = interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
    if (dt_bond_out) *dt_bond_out = dt_bond;
    const double dt = dt_factor * dt_bond;
    double amp = delta0;
    for (std::size_t n = 0; n < nsteps; ++n)
    {
        drv.step(c.state, c.topology, c.materials, dt);
        amp = std::max(amp, std::fabs((c.state.x[0][1] - c.state.x[0][0]) - l0 - s_eq));
    }
    return amp / delta0;
}

// То же для цепочки из N бусин (внутренние массы m, концевые m/2: старшая мода
// строго чередующаяся, omega = 2 sqrt(k/m) = sqrt2 * omega_пары). Предварительно
// растянута равномерно (s_eq на связь, концевые силы k·s_eq), возмущение --
// чередующееся, амплитуда delta0 = 1e-8 l0 (почти чистая старшая мода).
double chain_growth(std::size_t N, double xi, double dt_factor, std::size_t nsteps, bool half_end_masses = true)
{
    const double l0 = 0.01, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double A = area_of(d);
    const double m = rho * A * l0;
    fiber::FiberMaterial mat = make_material(E, d, rho);
    mat.visco.xi = xi;
    Chain c = make_chain(N, l0, mat, m);
    if (half_end_masses) c.state.mass[0] = c.state.mass[N - 1] = m / 2.0;
    const double s_eq = 0.5 * l0;
    const double k = E * A / l0;
    for (std::size_t i = 0; i < N; ++i) c.state.x[0][i] = static_cast<double>(i) * (l0 + s_eq);
    Driver drv;
    drv.ensure(N);
    drv.ext[0][0] = -k * s_eq / c.state.mass[0];
    drv.ext[0][N - 1] = k * s_eq / c.state.mass[N - 1];
    rhs::CpuExecutor exec(1);
    const double dt = dt_factor * interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
    const double delta0 = 1.0e-8 * l0;
    for (std::size_t i = 0; i < N; ++i) c.state.x[0][i] += ((i % 2) ? delta0 : -delta0) * 0.5;
    double amp = delta0;
    for (std::size_t n = 0; n < nsteps; ++n)
    {
        drv.step(c.state, c.topology, c.materials, dt);
        for (std::size_t b = 0; b + 1 < N; ++b)
            amp = std::max(amp, std::fabs(c.state.x[0][b + 1] - c.state.x[0][b] - l0 - s_eq));
        if (!std::isfinite(amp) || amp > 1.0) break;
    }
    return amp / delta0;
}

void test_t9_dt_limit()
{
    std::printf("Т9 (предел шага при демпфере):\n");
    // Формула лимита против ручного расчёта.
    {
        const double l0 = 0.01, E = 2.0e9, d = 1.0e-3, rho = 1200.0, m = 1.0e-6, xi = 2.0;
        const double A = area_of(d);
        fiber::FiberMaterial mat = make_material(E, d, rho);
        mat.visco.xi = xi;
        Chain c = make_chain(2, l0, mat, m);
        rhs::CpuExecutor exec(1);
        const double got = interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
        const double omega = std::sqrt(E * A / l0 / (m / 2.0));
        const double xc = std::sqrt(2.0) * xi;
        const double want = std::min(0.8 / omega, 0.8 * 2.0 * (std::sqrt(1.0 + 4.0 * xc * xc) - 2.0 * xc) / (std::sqrt(2.0) * omega));
        check(std::fabs(got / want - 1.0) < 1.0e-12, "dt_bond = min(0.8/ω, 0.8·h_кр(√2ξ)/(√2ω)), h_кр = 2(sqrt(1+4x²)-2x) (ξ = 2)");
        check(got < 0.8 / omega, "демпфер ужесточает лимит (ξ = 2)");
        // Малое ξ: демпфер лимит не ослабляет (min с упругим 0.8/ω) и не меняет.
        mat.visco.xi = 1.0e-3;
        c.materials[0] = mat;
        const double got_small = interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
        check(std::fabs(got_small / (0.8 / omega) - 1.0) < 1.0e-14, "малое ξ: лимит равен упругому (min с упругим)");
        // eta, пересчитанная из xi, даёт тот же лимит; E1 увеличивает жёсткость.
        fiber::FiberMaterial me = make_material(E, d, rho);
        me.visco.eta = fiber::kelvin_eta_from_xi(xi, E, A, l0, m / 2.0);
        c.materials[0] = me;
        const double got_eta = interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
        check(std::fabs(got_eta / want - 1.0) < 1.0e-12, "лимит при заданной eta совпадает с лимитом при xi");
        fiber::FiberMaterial m1 = make_material(E, d, rho);
        m1.visco.E1 = 3.0 * E;
        m1.visco.tau = 1.0e-5;
        c.materials[0] = m1;
        const double got_e1 = interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
        const double want_e1 = 0.8 / std::sqrt((4.0 * E) * A / l0 / (m / 2.0));
        check(std::fabs(got_e1 / want_e1 - 1.0) < 1.0e-14, "E_inst = E + E1 в лимите шага");
        { // материал twaron_plastic_visco_rate.fib: E1 = 0.1 E, tau, xi = 0.05 => лимит в sqrt(E/(E+E1)) раз меньше упругого
            fiber::FiberMaterial mv = make_material(E, d, rho);
            mv.visco.E1 = 0.1 * E;
            mv.visco.tau = 20.0e-6;
            mv.visco.xi = 0.05;
            c.materials[0] = mv;
            const double got_v = interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
            const double want_v = 0.8 / omega * std::sqrt(E / (1.1 * E));
            check(std::fabs(got_v / want_v - 1.0) < 1.0e-12, "E1 = 0.1E, xi = 0.05: лимит = упругий * sqrt(E/(E+E1)) = 0.9535 (Кельвин не ослабляет)");
        }
        fiber::FiberMaterial mp = make_material(E, d, rho);
        mp.plastic.sy = 1.0e7;
        c.materials[0] = mp;
        const double got_p = interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
        const double want_p = 0.8 / std::sqrt(E * A / l0 / (m / 2.0));
        check(bits_equal(got_p, want_p), "пластика лимит не ужесточает (значение прежнее побитово)");
    }
    for (const double xi : {0.5, 1.0, 1.2, 2.0, 5.0})
    {
        const std::size_t nsteps = xi > 3.0 ? 30000 : 12000;
        double dt_bond = 0.0;
        const double p09 = pair_growth(xi, 0.9, nsteps, &dt_bond);
        const double p30 = pair_growth(xi, 3.0, nsteps);
        const double c09 = chain_growth(12, xi, 0.9, nsteps, true);
        const double c20 = chain_growth(12, xi, 2.0, nsteps, true);
        const double u09 = chain_growth(12, xi, 0.9, nsteps, false);
        const double u125 = chain_growth(12, xi, 1.25, nsteps, false);
        const double u15 = chain_growth(12, xi, 1.5, nsteps, false);
        info("ξ = " + fmt(xi) + ": рост возмущения: пара 0.9/3.0·dt_bond -> " + fmt(p09) + "/" + fmt(p30) + "; цепочка (концы m/2) 0.9/2.0 -> " + fmt(c09) + "/" + fmt(c20) +
             "; цепочка (равные массы) 0.9/1.25/1.5 -> " + fmt(u09) + "/" + fmt(u125) + "/" + fmt(u15));
        const std::string tag = "ξ = " + fmt(xi) + ": ";
        check(p09 < 2.0, tag + "пара, 0.9·dt_bond устойчиво");
        check(p30 > 1.0e4, tag + "пара, 3.0·dt_bond неустойчиво (лимит не чрезмерно консервативен)");
        check(c09 < 2.0, tag + "цепочка из 12 бусин (концы m/2), 0.9·dt_bond устойчиво");
        check(c20 > 1.0e4, tag + "цепочка (концы m/2), 2.0·dt_bond неустойчиво");
        check(u09 < 2.0 && u125 < 2.0, tag + "цепочка (равные массы), 0.9 и 1.25·dt_bond устойчиво");
        check(u15 > 1.0e4, tag + "цепочка (равные массы), 1.5·dt_bond неустойчиво: лимит почти точный");
    }
}

// ----------------------------------------------------------------------
// Т10. Энергобаланс цикла нагрузка-разгрузка. Цепочка из 5 бусин, бусина 0
// движется по заданной скорости (полусинус вперёд, потом назад -- возврат в
// исходное положение), бусина 4 закреплена; материал с пластикой, ветвью
// Максвелла и демпфером Кельвина. Сравнивается E_кин + E_упр + W_p + W_v с
// работой внешних сил над заданной бусиной.
// ----------------------------------------------------------------------
fiber::FiberMaterial full_inelastic_material(double d, double rho)
{
    fiber::FiberMaterial mat = make_material(2.0e9, d, rho);
    mat.plastic.sy = 60.0e6;
    mat.plastic.H = 0.4e9;
    mat.visco.E1 = 0.5e9;
    mat.visco.tau = 30.0e-6;
    mat.visco.xi = 0.05;
    return mat;
}

struct CycleResult
{
    double work = 0.0, total = 0.0, resid_rel = 0.0, resid_max_rel = 0.0, clip = 0.0, resid_clip_rel = 0.0;
    double wp = 0.0, wv = 0.0, kin = 0.0, el = 0.0;
    bool monotone = true;
    bool all_plastic = false;
};

CycleResult run_cycle(const fiber::FiberMaterial& mat, double dt_factor)
{
    const double l0 = 1.0e-3, d = 1.0e-3, rho = 1200.0;
    const double A = area_of(d);
    const std::size_t N = 5;
    Chain c = make_chain(N, l0, mat, rho * A * l0);
    Driver drv;
    drv.ensure(N);
    drv.driven[0] = 1;
    drv.driven[N - 1] = 1;
    // бусина 0 уходит в -x (растягивает цепочку), бусина N-1 закреплена
    const double T_up = 60.0e-6, T_rest = 40.0e-6;
    const double eps_peak = 0.10;
    const double ext_total = eps_peak * l0 * static_cast<double>(N - 1);
    const double V = ext_total * pi / (2.0 * T_up);
    rhs::CpuExecutor exec(1);
    const double dt = dt_factor * interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
    const std::size_t nsteps = static_cast<std::size_t>((2.0 * T_up + T_rest) / dt);
    CycleResult res;
    std::vector<double> wp_prev(N - 1, 0.0), wv_prev(N - 1, 0.0), ep_prev(N - 1, 0.0);
    double peak_work = 0.0;
    double t = 0.0;
    for (std::size_t n = 0; n < nsteps; ++n)
    {
        double v = 0.0;
        if (t < T_up) v = -V * std::sin(pi * t / T_up);
        else if (t < 2.0 * T_up) v = V * std::sin(pi * (t - T_up) / T_up);
        drv.vdrive[0] = v;
        drv.step(c.state, c.topology, c.materials, dt);
        t += dt;
        for (std::size_t kb = 0; kb < N - 1; ++kb)
        {
            const fiber::BondState& s = c.topology.bonds[kb].st;
            if (!(s.w_p >= wp_prev[kb]) || !(s.w_v >= wv_prev[kb]) || !(s.w_p >= 0.0) || !(s.w_v >= 0.0) || !(s.eps_p >= ep_prev[kb])) res.monotone = false;
            wp_prev[kb] = s.w_p;
            wv_prev[kb] = s.w_v;
            ep_prev[kb] = s.eps_p;
        }
        peak_work = std::max(peak_work, std::fabs(drv.work_ext));
        const Energies e = energies(c.state, c.topology, c.materials);
        res.resid_max_rel = std::max(res.resid_max_rel, std::fabs(drv.work_ext - e.total())); // абсолютная, нормируется ниже
    }
    const Energies e = energies(c.state, c.topology, c.materials);
    res.work = drv.work_ext;
    res.total = e.total();
    res.wp = e.wp;
    res.wv = e.wv;
    res.kin = e.kin;
    res.el = e.el;
    res.resid_rel = std::fabs(drv.work_ext - e.total()) / peak_work;
    res.clip = drv.clip_loss;
    res.resid_clip_rel = std::fabs(drv.work_ext - (e.total() - drv.clip_loss)) / peak_work;
    res.resid_max_rel /= peak_work;
    res.all_plastic = true;
    for (const fiber::FiberBond& b : c.topology.bonds) res.all_plastic = res.all_plastic && b.st.eps_p > 0.0;
    return res;
}

// Невязка баланса -- свойство схемы, а не ошибка учёта энергии. Два источника
// (разобраны численно, прогоны одной и той же сцены с шагом 0.4 ... 0.05 от
// лимита, по отдельным материалам):
//  1) O(dt): сила в шаге берётся по состоянию на x^n (явная схема), Кельвин --
//     по скорости предиктора (запаздывает на dt/2), а энергии связи считаются
//     по x^{n+1} трапецией; невязка убывает ПРОПОРЦИОНАЛЬНО dt (проверяется:
//     при вдвое меньшем шаге вдвое меньше, у пластики и Кельвина это вся невязка);
//  2) НЕ зависит от dt: односторонность связи. Сила вдоль связи не бывает
//     отрицательной, поэтому σ = max(σ_модели, 0). При быстрой разгрузке ветвь
//     Максвелла (q < 0) и Кельвин (η·ε̇ < 0) делают σ_модели < 0 при ε > ε_p: связь
//     отдаёт внутреннюю энергию, а бусинам сила не передаётся -- энергия
//     "теряется" на отсечке. Драйвер теста считает её отдельно (clip_loss =
//     Σ σ_модели·Δε·V по шагам с отсечкой): после её учёта невязка снова убывает
//     пропорционально dt. Это физика односторонней модели (дизайн, разд. 1.3 и
//     9, риск 4), а не неточность интегратора.
struct EnergyCase
{
    const char* name;
    fiber::FiberMaterial mat;
    bool plastic, maxwell, kelvin;
    // Допуски (доли пика внешней работы) для шага 0.1·лимит, ~1.5-2 фактического значения.
    double tol_end, tol_max, tol_end_clip;
};

void test_t10_energy_balance()
{
    std::printf("Т10 (энергобаланс цикла нагрузка-разгрузка, по отдельным ветвям):\n");
    const double d = 1.0e-3, rho = 1200.0;
    std::vector<EnergyCase> cases;
    {
        fiber::FiberMaterial m = make_material(2.0e9, d, rho);
        m.plastic.sy = 60.0e6;
        m.plastic.H = 0.4e9;
        cases.push_back({"только пластика", m, true, false, false, 2.0e-3, 6.0e-3, 2.0e-3});
        m = make_material(2.0e9, d, rho);
        m.visco.E1 = 0.5e9;
        m.visco.tau = 30.0e-6;
        cases.push_back({"только Максвелл", m, false, true, false, 1.0e-2, 1.3e-2, 2.5e-3});
        m = make_material(2.0e9, d, rho);
        m.visco.xi = 0.05;
        cases.push_back({"только Кельвин", m, false, false, true, 3.0e-3, 5.0e-3, 3.0e-3});
        cases.push_back({"всё вместе", full_inelastic_material(d, rho), true, true, true, 4.0e-3, 7.0e-3, 1.5e-3});
    }
    for (const EnergyCase& ec : cases)
    {
        const CycleResult r1 = run_cycle(ec.mat, 0.1);
        const CycleResult r2 = run_cycle(ec.mat, 0.05);
        const std::string tag = std::string(ec.name) + ": ";
        info(tag + "A_внеш = " + fmt(r1.work) + " Дж, W_p = " + fmt(r1.wp) + ", W_v = " + fmt(r1.wv) + ", E_упр = " + fmt(r1.el) + ", E_кин = " + fmt(r1.kin) +
             "; невязка (конец/макс) " + fmt(r1.resid_rel) + "/" + fmt(r1.resid_max_rel) + "; потери на отсечке " + fmt(r1.clip) + " Дж, невязка за вычетом отсечки " + fmt(r1.resid_clip_rel));
        info(tag + "шаг вдвое меньше: невязка " + fmt(r2.resid_rel) + ", за вычетом отсечки " + fmt(r2.resid_clip_rel));
        check(r1.work > 0.0, tag + "внешняя работа за цикл положительна");
        check(r1.monotone, tag + "W_p >= 0, W_v >= 0, ε_p не убывает: для каждой связи и каждого шага");
        check(ec.plastic == (r1.wp > 0.0) && ec.plastic == r1.all_plastic, tag + "пластика есть ровно там, где включена (W_p, ε_p)");
        check((ec.maxwell || ec.kelvin) == (r1.wv > 0.0), tag + "вязкая диссипация есть ровно там, где включены Максвелл/Кельвин");
        check(r1.resid_rel < ec.tol_end, tag + "невязка баланса в конце цикла < " + fmt(ec.tol_end) + " пика внешней работы");
        check(r1.resid_max_rel < ec.tol_max, tag + "невязка по всему прогону < " + fmt(ec.tol_max) + " пика внешней работы");
        check(r1.resid_clip_rel < ec.tol_end_clip, tag + "невязка за вычетом потерь на отсечке < " + fmt(ec.tol_end_clip));
        check(r2.resid_clip_rel < 0.65 * r1.resid_clip_rel, tag + "за вычетом отсечки невязка убывает ~ пропорционально dt (вдвое меньший шаг -- < 0.65)");
        if (!ec.maxwell && !ec.kelvin) check(r1.clip == 0.0, tag + "без вязких ветвей отсечка не срабатывает");
        else if (ec.maxwell) check(r1.clip > 1.0e-3 * r1.work, tag + "у Максвелла потери на отсечке заметны (источник невязки, не зависящей от dt)");
    }
}

// Случайная кинематика: монотонность W_p, W_v, ε_p и неотрицательность на нерегулярных траекториях.
void test_t10_random()
{
    std::printf("Т10 (случайные траектории: W_p, W_v >= 0):\n");
    const double l0 = 1.0e-3, d = 1.0e-3, rho = 1200.0;
    const double A = area_of(d);
    const fiber::FiberMaterial mat = full_inelastic_material(d, rho);
    bool ok = true;
    bool any_plastic = false;
    for (std::uint64_t seed = 1; seed <= 5; ++seed)
    {
        Chain c = make_chain(6, l0, mat, rho * A * l0);
        Driver drv;
        drv.ensure(6);
        drv.driven[0] = 1;
        drv.driven[5] = 1; // второй конец закреплён
        std::uint64_t st = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        auto rnd = [&st]() {
            st = st * 6364136223846793005ULL + 1442695040888963407ULL;
            return static_cast<double>(st >> 11) / 9007199254740992.0;
        };
        rhs::CpuExecutor exec(1);
        const double dt = 0.3 * interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
        std::vector<double> wp_prev(5, 0.0), wv_prev(5, 0.0), ep_prev(5, 0.0);
        for (int n = 0; n < 4000; ++n)
        {
            if (n % 17 == 0) drv.vdrive[0] = 5.0 - 16.0 * rnd(); // м/с, в обе стороны (в среднем растяжение)
            drv.step(c.state, c.topology, c.materials, dt);
            for (std::size_t kb = 0; kb < 5; ++kb)
            {
                const fiber::BondState& s = c.topology.bonds[kb].st;
                if (!std::isfinite(s.w_p) || !std::isfinite(s.w_v) || !std::isfinite(s.q)) ok = false;
                if (!(s.w_p >= wp_prev[kb]) || !(s.w_v >= wv_prev[kb]) || !(s.eps_p >= ep_prev[kb])) ok = false;
                wp_prev[kb] = s.w_p;
                wv_prev[kb] = s.w_v;
                ep_prev[kb] = s.eps_p;
                if (s.eps_p > 0.0) any_plastic = true;
            }
        }
    }
    check(any_plastic, "самопроверка: на случайных траекториях возникает пластика");
    check(ok, "W_p, W_v, ε_p монотонны и неотрицательны на 5 случайных траекториях по 4000 шагов");
}

// ----------------------------------------------------------------------
// Т11. Импульс и момент импульса с пластикой и демпфером (сила вдоль связи,
// третий закон Ньютона: ни пластика, ни демпфер центральность не нарушают).
// Схема: v^{n+1} = v^n + dt·a(x^n, ·), x^{n+1} = x^n + dt·v^{n+1}, поэтому
// L = Σ m x×v сохраняется до округления.
// ----------------------------------------------------------------------
void test_t11_momentum()
{
    std::printf("Т11 (импульс и момент импульса с пластикой и демпфером):\n");
    const double l0 = 0.02, d = 8.0e-4, rho = 1000.0;
    fiber::FiberMaterial mat = make_material(1.5e9, d, rho);
    mat.plastic.sy = 4.0e6;
    mat.plastic.H = 0.1e9;
    mat.visco.E1 = 0.3e9;
    mat.visco.tau = 5.0e-5;
    mat.visco.xi = 0.1;
    const std::size_t N = 5;
    const double masses[N] = {1.0e-3, 1.6e-3, 0.9e-3, 2.2e-3, 1.3e-3};
    Chain c = make_chain(N, l0, mat, 1.0e-3);
    for (std::size_t i = 0; i < N; ++i)
    {
        c.state.mass[i] = masses[i];
        const double stretch = 1.0 + 0.04 * static_cast<double>(i + 1) / static_cast<double>(N);
        c.state.x[0][i] = (i == 0) ? 0.0 : (c.state.x[0][i - 1] + l0 * stretch);
        c.state.x[1][i] = 0.003 * static_cast<double>(i % 2 == 0 ? 1 : -1) * static_cast<double>(i);
        c.state.x[2][i] = 0.002 * static_cast<double>(i * i % 3);
        c.state.vx[0][i] = 0.1 * static_cast<double>(i) - 0.2;
        c.state.vx[1][i] = (i % 2 == 0) ? 0.05 : -0.03;
        c.state.vx[2][i] = 0.02 * static_cast<double>(i);
    }
    auto totals = [&](double p[dim], double L[dim]) {
        for (std::size_t dd = 0; dd < dim; ++dd) p[dd] = L[dd] = 0.0;
        for (std::size_t i = 0; i < N; ++i)
        {
            const double m = c.state.mass[i];
            for (std::size_t dd = 0; dd < dim; ++dd) p[dd] += m * c.state.vx[dd][i];
            L[0] += m * (c.state.x[1][i] * c.state.vx[2][i] - c.state.x[2][i] * c.state.vx[1][i]);
            L[1] += m * (c.state.x[2][i] * c.state.vx[0][i] - c.state.x[0][i] * c.state.vx[2][i]);
            L[2] += m * (c.state.x[0][i] * c.state.vx[1][i] - c.state.x[1][i] * c.state.vx[0][i]);
        }
    };
    double p0[dim], L0[dim], p1[dim], L1[dim];
    totals(p0, L0);
    rhs::CpuExecutor exec(1);
    const double dt = 0.3 * interactions::compute_fiber_bond_dt_limit(exec, c.state, c.topology, c.materials);
    Driver drv;
    drv.ensure(N);
    // Предиктор шага добавляет к скорости 0.5·dt·a_prev, а корректор возвращает
    // v_min + dt·a: импульс по v_min сохраняется точно. Сравниваем после шага,
    // добавив вклад полушага: P = Σ m (v - 0.5·dt·a_prev) -- проще: измерять по v_min.
    const int nsteps = 300;
    for (int s = 0; s < nsteps; ++s) drv.step(c.state, c.topology, c.materials, dt);
    totals(p1, L1);
    // Скорости в конце шага = v_min + dt·a (v^{n+1}) -- именно они входят в закон сохранения.
    double pscale = 0.0, Lscale = 0.0;
    for (std::size_t dd = 0; dd < dim; ++dd)
    {
        pscale = std::max(pscale, std::fabs(p0[dd]));
        Lscale = std::max(Lscale, std::fabs(L0[dd]));
    }
    double mom_scale_abs = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t dd = 0; dd < dim; ++dd) mom_scale_abs += c.state.mass[i] * std::fabs(c.state.vx[dd][i]);
    bool plastic = false, visc = false;
    for (const fiber::FiberBond& b : c.topology.bonds)
    {
        if (b.st.eps_p > 0.0) plastic = true;
        if (b.st.w_v > 0.0) visc = true;
    }
    check(plastic && visc, "самопроверка: пластика и вязкая диссипация действительно включались");
    {
        double dp = 0.0, dL = 0.0;
        for (std::size_t dd = 0; dd < dim; ++dd)
        {
            dp = std::max(dp, std::fabs(p1[dd] - p0[dd]));
            dL = std::max(dL, std::fabs(L1[dd] - L0[dd]));
        }
        info("max|Δp|/Σm|v| = " + fmt(dp / mom_scale_abs) + ", max|ΔL|/max|L| = " + fmt(dL / std::max(Lscale, 1e-300)));
    }
    // Сохранение импульса в схеме: Σ m·v^{n+1} = Σ m·v^n (+ полушаговый сдвиг предиктора у обоих полушагов
    // пропорционален Σ m·a = 0): допуск относительно масштаба Σ m|v|.
    for (std::size_t dd = 0; dd < dim; ++dd)
        check(std::fabs(p1[dd] - p0[dd]) < 1.0e-12 * mom_scale_abs,
              "импульс сохраняется (компонента " + std::to_string(dd) + "): " + fmt(p0[dd]) + " -> " + fmt(p1[dd]));
    double Lm = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t dd = 0; dd < dim; ++dd)
            Lm += c.state.mass[i] * std::fabs(c.state.x[(dd + 1) % dim][i] * c.state.vx[(dd + 2) % dim][i]);
    for (std::size_t dd = 0; dd < dim; ++dd)
        check(std::fabs(L1[dd] - L0[dd]) < 1.0e-10 * std::max(Lm, Lscale),
              "момент импульса сохраняется (компонента " + std::to_string(dd) + "): " + fmt(L0[dd]) + " -> " + fmt(L1[dd]));
}

// ----------------------------------------------------------------------
// Начальное eps_prev (init_fiber_bond_state). Связь, созданная при l != l0, без
// инициализации давала бы на первом шаге фиктивное Δε = ε - 0.
// ----------------------------------------------------------------------
void test_init_eps_prev()
{
    std::printf("Начальное eps_prev неупругих связей:\n");
    const double l0 = 1.0e-3, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double A = area_of(d);
    fiber::FiberMaterial elastic = make_material(E, d, rho);
    fiber::FiberMaterial visc = make_material(E, d, rho);
    visc.visco.xi = 0.1;
    visc.visco.E1 = 0.5e9;
    visc.visco.tau = 20.0e-6;
    // бусины 0-2 -- вязкоупругий материал 1, бусины 3-4 -- упругий материал 0
    Chain c = make_chain(5, l0, elastic, rho * A * l0);
    c.materials.push_back(visc);
    for (std::size_t i = 0; i < 3; ++i) c.state.material_id[i] = 1;
    c.topology.bonds.erase(c.topology.bonds.begin() + 2); // связь 2-3 убираем (разные материалы)
    const double eps0 = 0.05;
    for (std::size_t i = 0; i < 5; ++i) c.state.x[0][i] = static_cast<double>(i) * l0 * (1.0 + eps0);
    Chain without = c; // состояние нулевое, как до инициализации
    interactions::init_fiber_bond_state(c.state, c.topology, c.materials);
    bool ok_visc = true, ok_el = true;
    for (const fiber::FiberBond& b : c.topology.bonds)
    {
        const bool inel = c.materials[static_cast<std::size_t>(c.state.material_id[b.i])].inelastic();
        const double eps = interactions::fiber_bond_length(c.state, b) / b.l0 - 1.0;
        if (inel)
            ok_visc = ok_visc && b.st.eps_prev == eps && b.st.eps_p == 0.0 && b.st.q == 0.0 && b.st.w_p == 0.0 && b.st.w_v == 0.0;
        else
            ok_el = ok_el && bits_equal(b.st.eps_prev, 0.0) && b.st.eps_p == 0.0 && b.st.q == 0.0 && b.st.w_p == 0.0 && b.st.w_v == 0.0;
    }
    check(ok_visc, "неупругие связи: eps_prev = текущая ε, остальное состояние нулевое");
    check(ok_el, "упругие связи не тронуты (st нулевое)");

    // Первый шаг связи в покое при ε0 != 0: с инициализацией приращение Δε = 0,
    // без неё -- фиктивный скачок Δε = ε0 и огромные q и w_v.
    const double dt = 1.0e-8;
    rhs::CpuExecutor exec(1);
    interactions::apply_fiber_bond_constitutive(exec, c.state, c.topology, c.materials, dt);
    interactions::apply_fiber_bond_constitutive(exec, without.state, without.topology, without.materials, dt);
    const fiber::BondState& sw = c.topology.bonds[0].st;
    const fiber::BondState& so = without.topology.bonds[0].st;
    info("первый шаг при ε0 = " + fmt(eps0) + ": с init w_v = " + fmt(sw.w_v) + ", q = " + fmt(sw.q) + "; без init w_v = " + fmt(so.w_v) + ", q = " + fmt(so.q));
    check(sw.w_v == 0.0 && sw.q == 0.0, "с init: связь в покое не диссипирует и не копит q на первом шаге (Δε = 0)");
    check(so.w_v > 1.0e3 && so.q > 1.0e6, "без init: фиктивный скачок даёт большие w_v и q (самопроверка теста)");
}

// ----------------------------------------------------------------------
// Полный боевой цикл: time_integration(cfg, fiber::Model&, ...) на малой
// сцене из свободной нити с неупругим материалом (без Solid, ntotal = 0, как
// tests/clamp_integration_fiber.hpp). Свободная нить растягивается линейным
// профилем скорости (концы разлетаются, центр в покое), одна связь ослаблена и
// рвётся (l_break = 1.02·l0), остальные до разрыва проходят предел текучести.
// Проверяется то, что драйверы тестов с ручным шагом не покрывают: проход
// apply_fiber_bond_constitutive реально вызывается в цикле (внутри
// run_fiber_topology_postprocessing), eps_p, w_p, w_v ненулевые, они
// монотонны, eps_prev на каждом шаге = текущей ε, состояние выживших связей
// сохраняется при разрыве.
// ----------------------------------------------------------------------
fiber::FiberMaterial loop_material(double d, double rho)
{
    fiber::FiberMaterial m = make_material(1.0e7, d, rho);
    m.plastic.sy = 1.0e5; // ε_y = 1%
    m.plastic.H = 2.0e6;
    m.visco.xi = 0.02;
    m.visco.E1 = 2.0e6;
    m.visco.tau = 2.0e-3;
    return m;
}

struct LoopScene
{
    fiber::Model model;
    SimulationConfig cfg;
};

LoopScene make_loop_scene(std::size_t n, double l0, double stretch0, double v_max)
{
    const double d = 1.0e-3, rho = 1000.0;
    LoopScene sc;
    sc.model.solid.resize(1); // фиктивная частица, ntotal = 0
    sc.model.solid.inzializ_index(0);
    sc.model.materials.push_back(loop_material(d, rho));
    const fiber::FiberMaterial& mat = sc.model.materials[0];
    sc.model.fiber.resize(n);
    sc.model.topology.resize(n);
    const double A = area_of(d);
    for (std::size_t k = 0; k < n; ++k)
    {
        for (std::size_t dd = 0; dd < dim; ++dd)
        {
            sc.model.fiber.x[dd][k] = 0.0;
            sc.model.fiber.vx[dd][k] = 0.0;
        }
        sc.model.fiber.x[0][k] = static_cast<double>(k) * l0 * (1.0 + stretch0);
        // линейный профиль скорости, центр в покое: концы разлетаются, постоянная скорость деформации
        sc.model.fiber.vx[0][k] = v_max * (2.0 * static_cast<double>(k) / static_cast<double>(n - 1) - 1.0);
        sc.model.fiber.mass[k] = mat.rho * A * l0 * ((k == 0 || k + 1 == n) ? 0.5 : 1.0);
        sc.model.fiber.object_id[k] = 0;
        sc.model.fiber.material_id[k] = 0;
        if (k > 0) sc.model.topology.prev[k] = static_cast<long long>(k - 1);
        if (k + 1 < n) sc.model.topology.next[k] = static_cast<long long>(k + 1);
    }
    for (std::size_t k = 0; k + 1 < n; ++k)
    {
        fiber::FiberBond b;
        b.i = k;
        b.j = k + 1;
        b.l0 = l0;
        b.l_break = l0 * 1.0e9;
        sc.model.topology.bonds.push_back(b);
    }
    sc.cfg.apply_defaults();
    sc.cfg.maxn = 1;
    sc.cfg.max_interaction = 100;
    sc.cfg.moni_particle = 0;
    sc.cfg.t_total = 1.0e30;
    sc.cfg.save_t = 1.0e30;
    sc.cfg.print_step = 1000000000;
    return sc;
}

void test_full_loop()
{
    std::printf("Боевой цикл time_integration с неупругой нитью (разрыв одной связи):\n");
    const double l0 = 0.01;
    const std::size_t n = 8;
    {
        // --- нить без предварительного натяжения: ненулевые eps_p, w_p, w_v; монотонность; разрыв ---
        LoopScene sc = make_loop_scene(n, l0, 0.0, 5.0);
        sc.model.topology.bonds[3].l_break = l0 * 1.02; // слабая связь: рвётся при ε = 2% (предел текучести 1%)
        struct Rec
        {
            std::size_t i, j;
            fiber::BondState st;
        };
        std::vector<Rec> prev;
        bool monotone = true, prev_ok = true, survivors_kept = true, broke = false;
        std::size_t nb_prev = sc.model.topology.bonds.size();
        std::size_t nsteps_seen = 0;
        e2e::ScopedWorkDir wd("sph_plast_loop", "a");
        RunGuard guard(0.0, 0.0, 0);
        guard.set_step_limit(600);
        guard.set_step_observer([&](double, double) {
            ++nsteps_seen;
            const fiber::Model& m = sc.model;
            const std::size_t nb = m.topology.bonds.size();
            if (nb < nb_prev) broke = true;
            for (const fiber::FiberBond& b : m.topology.bonds)
            {
                const double eps = interactions::fiber_bond_length(m.fiber, b) / b.l0 - 1.0;
                if (!bits_equal(b.st.eps_prev, eps)) prev_ok = false;
                for (const Rec& r : prev)
                    if (r.i == b.i && r.j == b.j)
                    {
                        if (b.st.eps_p < r.st.eps_p || b.st.w_p < r.st.w_p || b.st.w_v < r.st.w_v) monotone = false;
                        if (nb < nb_prev && b.st.eps_p < r.st.eps_p) survivors_kept = false; // шаг с разрывом
                    }
            }
            prev.clear();
            for (const fiber::FiberBond& b : m.topology.bonds) prev.push_back({b.i, b.j, b.st});
            nb_prev = nb;
        });
        double dt = 1.0e-9, time = 0.0;
        {
            e2e::ScopedSilence silence(true);
            time_integration(sc.cfg, sc.model, 0, dt, time, 0.0, &guard);
        }
        info("шагов " + std::to_string(nsteps_seen) + ", t = " + fmt(time) + ", связей осталось " + std::to_string(sc.model.topology.bonds.size()));
        check(broke && sc.model.topology.bonds.size() == n - 2, "слабая связь порвалась (осталось n-2 связи)");
        double ep_max = 0.0, wp = 0.0, wv = 0.0;
        std::size_t n_plastic = 0;
        for (const fiber::FiberBond& b : sc.model.topology.bonds)
        {
            ep_max = std::max(ep_max, b.st.eps_p);
            wp += b.st.w_p;
            wv += b.st.w_v;
            if (b.st.eps_p > 0.0) ++n_plastic;
        }
        info("в конце: max eps_p = " + fmt(ep_max) + ", ΣW_p = " + fmt(wp) + ", ΣW_v = " + fmt(wv) + " Дж/м^3");
        check(ep_max > 0.01 && wp > 0.0 && wv > 0.0, "в боевом цикле накоплены eps_p, w_p, w_v (проход Postprocessing вызывается)");
        check(n_plastic >= 4, "не менее четырёх выживших связей получили eps_p > 0 (концевые связи разгружаются свободным концом)");
        check(prev_ok, "eps_prev == текущая ε на каждом шаге для каждой связи (проход вызван каждый шаг)");
        check(monotone, "eps_p, w_p, w_v не убывают на каждом шаге у каждой выжившей связи");
        check(survivors_kept, "при разрыве состояние оставшихся связей сохраняется (не обнуляется)");
        check(sc.model.topology.version >= 1, "topology.version увеличена разрывом");
    }
    {
        // --- нить, созданная при l != l0 (ε0 = 0.5%, упругая зона): нет фиктивного скачка на первом шаге ---
        const double eps0 = 0.005;
        LoopScene sc = make_loop_scene(n, l0, eps0, 0.0);
        e2e::ScopedWorkDir wd("sph_plast_loop", "b");
        RunGuard guard(0.0, 0.0, 0);
        guard.set_step_limit(1);
        double dt = 1.0e-9, time = 0.0;
        {
            e2e::ScopedSilence silence(true);
            time_integration(sc.cfg, sc.model, 0, dt, time, 0.0, &guard);
        }
        const fiber::FiberMaterial& mat = sc.model.materials[0];
        const double A = area_of(mat.d);
        const double m_red = sc.model.fiber.mass[1] / 2.0;
        const double eta = fiber::kelvin_eta_from_xi(mat.visco.xi, interactions::fiber_instant_modulus(mat), A, l0, m_red);
        const double spurious = eta * eps0 * eps0 / dt; // что дал бы скачок Δε = ε0 за один шаг, Дж/м^3
        double wv_max = 0.0, q_max = 0.0;
        for (const fiber::FiberBond& b : sc.model.topology.bonds)
        {
            wv_max = std::max(wv_max, b.st.w_v);
            q_max = std::max(q_max, std::fabs(b.st.q));
        }
        info("первый шаг при ε0 = " + fmt(eps0) + ": max w_v = " + fmt(wv_max) + " (скачок Δε = ε0 дал бы ~" + fmt(spurious) + "), max |q| = " + fmt(q_max));
        check(guard.steps() == 1, "сделан ровно один шаг");
        // Концевые связи на первом шаге реально меняют длину (свободный конец релаксирует: Δε ~ 0.2·ε0),
        // внутренние (1..5) не меняют её вовсе -- у них w_v и q должны остаться нулевыми.
        double wv_in = 0.0, q_in = 0.0;
        for (std::size_t k = 1; k + 2 < n; ++k)
        {
            wv_in = std::max(wv_in, sc.model.topology.bonds[k].st.w_v);
            q_in = std::max(q_in, std::fabs(sc.model.topology.bonds[k].st.q));
        }
        check(wv_in < 1.0e-9 * spurious && q_in < 1.0e-9 * mat.visco.E1 * eps0, "первый шаг: внутренние связи (Δε = 0) не получили ни w_v, ни q (без init дали бы ~" + fmt(spurious) + ")");
        check(wv_max < 0.1 * spurious && q_max < 0.3 * mat.visco.E1 * eps0, "первый шаг: даже концевые связи (реальное Δε) далеки от фиктивного скачка Δε = ε0");
    }
}

} // namespace

// ----------------------------------------------------------------------
// Т14. Критерии разрыва (раздел 7 дизайна): ε > eps_break (прежний), ε_p >=
// FAILURE.EPS_P, σ_full = E(ε - ε_p) + q >= FAILURE.SIGMA_U -- по КОНЕЧНОМУ состоянию
// шага (после apply_fiber_bond_constitutive). Цепочка из 4 бусин (3 связи):
// бусина 0 закреплена, бусина 3 тянется с постоянной скоростью; материал связи
// k -- материал бусины k (связь 0 -- "слабая", связи 1, 2 -- уцелевшие).
// На КАЖДОМ шаге независимо вычисленное предсказание (критерий по состоянию
// связи) сравнивается с найденными событиями, поэтому разрыв срабатывает ровно
// на том шаге, где конечное состояние перешло порог (ни раньше, ни позже).
// ----------------------------------------------------------------------
struct T14Result
{
    bool broke = false;
    std::size_t steps = 0;
};

template <class Pred>
T14Result run_t14(const std::string& name, std::vector<fiber::FiberMaterial> mats, Pred crossed,
                  interactions::FiberBreakReason expect_reason, bool expect_survivor_plastic)
{
    std::printf("%s:\n", name.c_str());
    const double l0 = 1.0e-3, d = 1.0e-3, rho = 1200.0;
    const double A = area_of(d);
    Chain c = make_chain(4, l0, mats[0], rho * A * l0);
    c.materials = mats;
    c.state.material_id[0] = 0;
    for (std::size_t i = 1; i < 4; ++i) c.state.material_id[i] = 1;
    for (fiber::FiberBond& b : c.topology.bonds)
        b.l_break = (1.0 + mats[static_cast<std::size_t>(c.state.material_id[b.i])].eps_break) * b.l0;

    Driver drv;
    drv.ensure(4);
    drv.driven[0] = 1;
    drv.driven[3] = 1;
    drv.vdrive[3] = 1.0;
    const double dt = 2.0e-8;
    std::vector<std::vector<interactions::FiberBreakEvent>> buf;
    T14Result res;
    bool equiv_ok = true, version_ok = true;
    for (std::size_t step = 0; step < 400000 && !res.broke; ++step)
    {
        const unsigned long long ver0 = c.topology.version;
        drv.step(c.state, c.topology, c.materials, dt);
        ++res.steps;
        const auto events = interactions::find_fiber_break_events(drv.exec, c.state, c.topology, c.materials, buf);
        std::vector<char> expect(c.topology.bonds.size(), 0);
        for (std::size_t k = 0; k < c.topology.bonds.size(); ++k)
        {
            const fiber::FiberBond& b = c.topology.bonds[k];
            const fiber::FiberMaterial& m = c.materials[static_cast<std::size_t>(c.state.material_id[b.i])];
            const double l = interactions::fiber_bond_length(c.state, b);
            expect[k] = (l > b.l_break || crossed(b, m, l / b.l0 - 1.0)) ? 1 : 0;
        }
        std::vector<char> got(c.topology.bonds.size(), 0);
        for (const auto& e : events) got[e.bond_index] = 1;
        if (got != expect) equiv_ok = false;
        if (events.empty())
        {
            if (c.topology.version != ver0) version_ok = false;
            continue;
        }
        res.broke = true;
        check(events.size() == 1 && events[0].bond_index == 0, "рвётся ровно связь 0 (слабая)");
        check(events[0].reason == expect_reason, "причина разрыва корректна");
        // Прежняя перегрузка (без материалов) FAILURE не видит.
        if (expect_reason != interactions::FiberBreakReason::Strain)
            check(interactions::find_fiber_break_events(drv.exec, c.state, c.topology, buf).empty(),
                  "перегрузка без материалов не учитывает FAILURE");
        const fiber::BondState s1 = c.topology.bonds[1].st, s2 = c.topology.bonds[2].st;
        const std::size_t nb = c.topology.bonds.size();
        interactions::apply_fiber_break_events(c.state, c.topology, events);
        check(c.topology.version == ver0 + 1, "topology.version +1 за проход");
        check(c.topology.bonds.size() == nb - 1, "связь удалена из topology.bonds");
        check(c.state.object_id[0] != c.state.object_id[1] && c.state.object_id[1] == c.state.object_id[3],
              "хвост получил новый object_id");
        auto same = [](const fiber::BondState& a, const fiber::BondState& b) {
            return bits_equal(a.eps_p, b.eps_p) && bits_equal(a.q, b.q) && bits_equal(a.eps_prev, b.eps_prev) &&
                   bits_equal(a.w_p, b.w_p) && bits_equal(a.w_v, b.w_v);
        };
        check(c.topology.bonds[0].i == 1 && same(c.topology.bonds[0].st, s1) && same(c.topology.bonds[1].st, s2),
              "eps_p/q/eps_prev/w_p/w_v уцелевших связей побитово сохранены");
        if (expect_survivor_plastic)
            check(s1.eps_p > 0.0 && s1.q != 0.0 && s1.eps_prev != 0.0 && s1.w_p > 0.0 && s1.w_v > 0.0,
                  "состояние уцелевшей связи нетривиально (eps_p " + fmt(s1.eps_p) + ", q " + fmt(s1.q) + ")");
        else
            check(s1.eps_p == 0.0 && s1.q == 0.0 && s1.w_p == 0.0 && s1.w_v == 0.0,
                  "упругая уцелевшая связь: st нулевой");
    }
    check(res.broke, "разрыв произошёл (" + std::to_string(res.steps) + " шагов)");
    check(equiv_ok, "на каждом шаге события == предсказание по конечному состоянию (разрыв на шаге перехода порога)");
    check(version_ok, "без разрыва topology.version не меняется");
    return res;
}

fiber::FiberMaterial t14_plastic(double sy)
{
    fiber::FiberMaterial m = make_material(2.0e9, 1.0e-3, 1200.0);
    m.plastic.sy = sy;
    m.plastic.H = 0.4e9;
    m.visco.E1 = 0.5e9; // ветвь Максвелла: q != 0 у уцелевших связей
    m.visco.tau = 2.0e-5;
    return m;
}

void test_t14_failure_criteria()
{
    using Reason = interactions::FiberBreakReason;
    { // EPS_P при пластике
        fiber::FiberMaterial weak = t14_plastic(60.0e6), strong = t14_plastic(62.0e6);
        weak.failure.eps_p = 0.03;
        run_t14("Т14a: разрыв по FAILURE.EPS_P", {weak, strong},
                [](const fiber::FiberBond& b, const fiber::FiberMaterial& m, double) {
                    return m.failure.eps_p > 0.0 && b.st.eps_p >= m.failure.eps_p;
                },
                Reason::PlasticStrain, true);
    }
    { // SIGMA_U при пластике
        fiber::FiberMaterial weak = t14_plastic(60.0e6), strong = t14_plastic(62.0e6);
        weak.failure.sigma_u = 70.0e6;
        run_t14("Т14b: разрыв по FAILURE.SIGMA_U (с пластикой)", {weak, strong},
                [](const fiber::FiberBond& b, const fiber::FiberMaterial& m, double eps) {
                    return m.failure.sigma_u > 0.0 && m.E * (eps - b.st.eps_p) + b.st.q >= m.failure.sigma_u;
                },
                Reason::Stress, true);
    }
    { // SIGMA_U без PLASTIC/VISCO: хрупкий упругий материал
        fiber::FiberMaterial weak = make_material(2.0e9, 1.0e-3, 1200.0), strong = weak;
        weak.failure.sigma_u = 40.0e6;
        check(!weak.inelastic() && weak.has_failure() && !strong.has_failure(),
              "FAILURE без PLASTIC/VISCO: inelastic() == false, has_failure() == true");
        run_t14("Т14c: разрыв по FAILURE.SIGMA_U без пластики (хрупкий)", {weak, strong},
                [](const fiber::FiberBond&, const fiber::FiberMaterial& m, double eps) {
                    return m.failure.sigma_u > 0.0 && m.E * eps >= m.failure.sigma_u;
                },
                Reason::Stress, false);
        // Ядро остаётся упругим: силы побитово такие же, как без FAILURE.
        fiber::FiberMaterial plain = weak;
        plain.failure = fiber::FailureParams{};
        Chain ca = make_chain(3, 1.0e-3, weak, 1.0e-6), cb = make_chain(3, 1.0e-3, plain, 1.0e-6);
        for (std::size_t i = 0; i < 3; ++i)
        {
            ca.state.x[0][i] = cb.state.x[0][i] =
                1.0e-3 * static_cast<double>(i) * (1.0 + 0.01 * static_cast<double>(i));
            ca.state.vx[0][i] = cb.state.vx[0][i] = 0.3 * static_cast<double>(i);
        }
        rhs::CpuExecutor ex(1);
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> b1, b2;
        const auto da = interactions::compute_fiber_bond(ex, ca.state, ca.topology, ca.materials, b1);
        const auto db = interactions::compute_fiber_bond(ex, cb.state, cb.topology, cb.materials, b2);
        bool same = da.size() == db.size();
        for (std::size_t i = 0; same && i < da.size(); ++i)
            for (std::size_t k = 0; k < dim; ++k) same = same && bits_equal(da[i].dv[k], db[i].dv[k]);
        check(same, "силы RHS для материала только с FAILURE побитово равны упругим");
    }
    { // Старый eps_break при включённой пластике (и далёком FAILURE.EPS_P)
        fiber::FiberMaterial weak = t14_plastic(60.0e6), strong = t14_plastic(62.0e6);
        weak.eps_break = 0.05;
        weak.failure.eps_p = 0.5;
        run_t14("Т14d: разрыв по eps_break при включённой пластике", {weak, strong},
                [](const fiber::FiberBond&, const fiber::FiberMaterial&, double) { return false; }, Reason::Strain,
                true);
    }
}

// ----------------------------------------------------------------------
// Юнит-кейсы find_fiber_break_events без динамики: координаты бусин и st связей
// заданы вручную, связь 0 -- между бусинами 0 и 1, если не сказано иначе.
// ----------------------------------------------------------------------
namespace
{

const char* reason_name(interactions::FiberBreakReason r)
{
    switch (r)
    {
    case interactions::FiberBreakReason::Strain: return "Strain";
    case interactions::FiberBreakReason::PlasticStrain: return "PlasticStrain";
    case interactions::FiberBreakReason::Stress: return "Stress";
    }
    return "?";
}

// Цепочка бусин с заданными x (ось x), связи i--i+1 с номинальной длиной l0 и
// порогом l_break = lb_factor * l0 (по умолчанию -- недостижимый).
Chain manual_chain(const std::vector<double>& xs, double l0, const std::vector<fiber::FiberMaterial>& mats,
                   double lb_factor = 1.0e9)
{
    Chain c = make_chain(xs.size(), l0, mats[0], 1.0e-6);
    c.materials = mats;
    for (std::size_t i = 0; i < xs.size(); ++i) c.state.x[0][i] = xs[i];
    for (fiber::FiberBond& b : c.topology.bonds) b.l_break = lb_factor * l0;
    return c;
}

std::vector<interactions::FiberBreakEvent> find_with_threads(const Chain& c, int nthreads)
{
    rhs::CpuExecutor ex(nthreads);
    std::vector<std::vector<interactions::FiberBreakEvent>> buf;
    return interactions::find_fiber_break_events(ex, c.state, c.topology, c.materials, buf);
}

} // namespace

void test_break_unit_cases()
{
    using R = interactions::FiberBreakReason;
    const double l0 = 1.0e-3;

    std::printf("Unit (a): SIGMA_U при eps < eps_p -- связь не рвётся:\n");
    {
        // eps = 0.01 < eps_p = 0.02 => sigma_ep = E(eps - eps_p) = -2e7 Па < SIGMA_U = 5e7.
        // EPS_P = 0.05 задан, но не достигнут (0.02 < 0.05).
        fiber::FiberMaterial m = t14_plastic(60.0e6);
        m.failure.eps_p = 0.05;
        m.failure.sigma_u = 5.0e7;
        Chain c = manual_chain({0.0, 1.01 * l0}, l0, {m});
        c.topology.bonds[0].st.eps_p = 0.02;
        check(find_with_threads(c, 1).empty(), "eps=0.01 < eps_p=0.02: SIGMA_U и EPS_P не срабатывают");
    }

    std::printf("Unit (b): EPS_P без PLASTIC (st нулевой) -- не рвёт:\n");
    {
        fiber::FiberMaterial m = make_material(2.0e9, 1.0e-3, 1200.0); // упругий материал
        m.failure.eps_p = 0.01;
        check(!m.inelastic() && m.has_failure(), "материал упругий, FAILURE.EPS_P задан");
        Chain c = manual_chain({0.0, 1.5 * l0}, l0, {m}); // eps = 0.5, st нулевой
        check(find_with_threads(c, 1).empty(), "eps=0.5, st.eps_p=0: EPS_P не срабатывает (eps_p == 0)");
    }

    std::printf("Unit (c): одна связь с двумя причинами (l > l_break и eps_p >= EPS_P):\n");
    {
        fiber::FiberMaterial m = t14_plastic(60.0e6);
        m.failure.eps_p = 0.03;
        m.failure.sigma_u = 1.0e7;
        Chain c = manual_chain({0.0, 1.02 * l0}, l0, {m}, 1.01); // l = 1.02 l0 > l_break = 1.01 l0
        c.topology.bonds[0].st.eps_p = 0.05;                     // eps_p >= EPS_P
        const auto ev = find_with_threads(c, 1);
        check(ev.size() == 1, "ровно одно событие");
        check(!ev.empty() && ev[0].bond_index == 0 && ev[0].reason == R::Strain,
              std::string("reason = Strain (приоритет l > l_break), получено ") +
                  (ev.empty() ? "нет события" : reason_name(ev[0].reason)));
    }

    std::printf("Unit (d): несколько одновременных разрывов -- все найдены, отсортированы:\n");
    {
        // Связи: 0: l0 (цела), 1: 1.1 l0 (l > l_break), 2: l0 (цела),
        //        3: 1.1 l0 (l > l_break), 4: l0 с eps_p = 0.05 >= EPS_P = 0.03.
        fiber::FiberMaterial m = t14_plastic(60.0e6);
        m.failure.eps_p = 0.03;
        const double seg[5] = {1.0, 1.1, 1.0, 1.1, 1.0};
        std::vector<double> xs{0.0};
        for (double f : seg) xs.push_back(xs.back() + f * l0);
        Chain c = manual_chain(xs, l0, {m}, 1.05);
        c.topology.bonds[4].st.eps_p = 0.05;

        const auto ev1 = find_with_threads(c, 1);
        const auto ev4 = find_with_threads(c, 4);
        const std::vector<std::size_t> want{1, 3, 4};
        bool idx_ok = ev1.size() == want.size();
        for (std::size_t k = 0; idx_ok && k < want.size(); ++k) idx_ok = ev1[k].bond_index == want[k];
        check(idx_ok, "найдены все 3 разрыва, индексы 1, 3, 4 по возрастанию");
        check(ev1.size() == 3 && ev1[0].reason == R::Strain && ev1[1].reason == R::Strain &&
                  ev1[2].reason == R::PlasticStrain,
              "причины: Strain, Strain, PlasticStrain");
        bool same = ev1.size() == ev4.size();
        for (std::size_t k = 0; same && k < ev1.size(); ++k)
            same = ev1[k].bond_index == ev4[k].bond_index && ev1[k].reason == ev4[k].reason;
        check(same, "1 поток и 4 потока дают одинаковый список событий");
    }

    std::printf("Unit (e): SIGMA_U численно (E, l0, l, eps_p заданы вручную):\n");
    {
        // E = 2e9 Па, l0 = 1e-3 м, l = 1.03e-3 м => eps = 0.03; eps_p = 0.01.
        // sigma_ep = E(eps - eps_p) = 2e9 * 0.02 = 4.0e7 Па (расчёт на бумаге).
        const double kSigmaEpPaper = 4.0e7;
        const double eps = 0.03, eps_p = 0.01;
        fiber::FiberMaterial base = t14_plastic(60.0e6); // E = 2e9
        base.failure.eps_p = 0.5;                        // EPS_P не достигается
        const double sigma_calc = base.E * (eps - eps_p);
        check(std::fabs(sigma_calc - kSigmaEpPaper) < 1.0e-6 * kSigmaEpPaper,
              "E(eps - eps_p) = " + fmt(sigma_calc) + " Па совпадает с расчётом на бумаге 4.0e7");

        auto run = [&](double sigma_u) {
            fiber::FiberMaterial m = base;
            m.failure.sigma_u = sigma_u;
            Chain c = manual_chain({0.0, 1.03 * l0}, l0, {m}, 1.5);
            c.topology.bonds[0].st.eps_p = eps_p;
            return find_with_threads(c, 1);
        };
        const auto below = run(3.99e7); // порог чуть ниже 4.0e7 -- рвётся
        const auto above = run(4.01e7); // порог чуть выше 4.0e7 -- не рвётся
        check(below.size() == 1 && below[0].reason == R::Stress, "SIGMA_U = 3.99e7 < 4.0e7: рвётся, reason = Stress");
        check(above.empty(), "SIGMA_U = 4.01e7 > 4.0e7: не рвётся");
    }

    std::printf("Unit (f): SIGMA_U по полному напряжению E(eps - eps_p) + q, q != 0:\n");
    {
        // E = 2e9, eps = 0.03, eps_p = 0.01 => E(eps - eps_p) = 4.0e7; q = 1.0e7 (ветвь Максвелла)
        // => sigma_full = 5.0e7 Па (на бумаге). Кельвин (eta*eps_dot) в критерий не входит.
        const double eps_p = 0.01, q = 1.0e7;
        auto run = [&](double sigma_u, double E1, double qv, double eta) {
            fiber::FiberMaterial m = t14_plastic(60.0e6);
            m.failure.eps_p = 0.5;
            m.visco.E1 = E1;
            m.visco.eta = eta;
            m.failure.sigma_u = sigma_u;
            Chain c = manual_chain({0.0, 1.03 * l0}, l0, {m}, 1.5);
            c.topology.bonds[0].st.eps_p = eps_p;
            c.topology.bonds[0].st.q = qv;
            c.topology.bonds[0].st.eps_prev = 0.0; // большой "eps_dot" -- не должен влиять
            return find_with_threads(c, 1);
        };
        const auto below = run(4.99e7, 0.5e9, q, 0.0); // чуть ниже 4.0e7 + 1.0e7 -- рвётся
        const auto above = run(5.01e7, 0.5e9, q, 0.0); // чуть выше -- не рвётся
        check(below.size() == 1 && below[0].reason == R::Stress, "SIGMA_U = 4.99e7 < E(eps-eps_p)+q = 5.0e7: рвётся, reason = Stress");
        check(above.empty(), "SIGMA_U = 5.01e7 > 5.0e7: не рвётся");
        check(run(4.5e7, 0.5e9, q, 0.0).size() == 1, "SIGMA_U = 4.5e7: с q рвётся (без q: 4.0e7 < 4.5e7 не рвалась бы)");
        check(run(4.5e7, 0.5e9, 0.0, 0.0).empty(), "q = 0 (Максвелл не нагружен): прежнее поведение, 4.0e7 < 4.5e7 не рвётся");
        check(run(4.5e7, 0.0, 0.0, 0.0).empty() && run(3.99e7, 0.0, 0.0, 0.0).size() == 1,
              "без Максвелла (E1 = 0, q = 0): прежнее поведение E(eps-eps_p) vs SIGMA_U");
        check(run(5.01e7, 0.5e9, q, 1.0e3).empty(), "член Кельвина (eta > 0) в критерий не входит");
    }
}

int main()
{
    test_t0_elastic_bitwise();
    test_plastic_full_cycle();
    test_t8_damped_oscillation();
    test_t9_dt_limit();
    test_t10_energy_balance();
    test_t10_random();
    test_t11_momentum();
    test_init_eps_prev();
    test_full_loop();
    test_t14_failure_criteria();
    test_break_unit_cases();

    if (failures == 0)
        std::printf("test_fiber_plasticity_integration: OK\n");
    else
        std::printf("test_fiber_plasticity_integration: %d failure(s)\n", failures);
    return failures;
}
