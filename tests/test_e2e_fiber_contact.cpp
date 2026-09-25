// tests/test_e2e_fiber_contact.cpp
// E2E-критерий готовности этапа Ф5 (fiber/PLAN.md, «Fiber-Fiber-Contact»,
// «Готово, когда»): «две нити, летящие навстречу крест-накрест, отскакивают
// без взаимного проникания».
//
// Сцена: нить P вдоль оси x на высоте z=+z0, летит вниз (vz=-v0); нить Q
// вдоль оси y на высоте z=-z0, летит вверх (vz=+v0) -- проекции пересекаются
// в окрестности (x,y)=(0,0). Обе нити прямые, несколько бусин, связаны
// Fiber-Bond (interactions/fiber_bond.hpp) и взаимодействуют через
// Fiber-Fiber-Contact (interactions/fiber_contact.hpp). Интегрирование --
// тот же PredictorCorrector, что в src/time_integration.cpp (перегрузка
// fiber::Model&, этап Ф3/Ф5), локальным драйвером -- по образцу
// FiberPcIntegrator из tests/test_fiber_bond.cpp (там же обоснование, почему
// не используется полноценный time_integration(fiber::Model&,...): требует
// Solid-подсистему, не относящуюся к предмету теста), расширенным суммой
// Fiber-Bond + Fiber-Fiber-Contact в один dv (как в src/time_integration.cpp,
// перегрузка fiber::Model&, блок «RHS: Fiber-Fiber-Contact»).
//
// Критерии готовности:
//   (а) нити НЕ проходят друг сквозь друга -- минимальное наблюдённое
//       расстояние между осями (в контакте) не опускается намного ниже d;
//   (б) после контакта нити расходятся -- представительные скорости вдоль
//       оси сближения (z) меняют знак (были навстречу -- становятся врозь).
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

static_assert(dim == 3, "test_e2e_fiber_contact требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

// -----------------------------------------------------------------------
// Локальный PredictorCorrector-драйвер: Fiber-Bond + Fiber-Fiber-Contact
// суммируются в один dv_fiber, шаг dt -- ДИНАМИЧЕСКИЙ (пересчитывается
// каждый шаг как в продовом цикле, а не фиксированный, как в более простых
// тестах test_fiber_bond.cpp -- здесь нужен, т.к. в момент контакта
// жёсткость Герца резко растёт и фиксированный dt из состояния "до
// контакта" перестал бы быть устойчивым).
// -----------------------------------------------------------------------
struct FiberBondContactIntegrator
{
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> bond_buffers;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> contact_buffers;
    std::vector<double> dv_prev[dim];
    std::vector<double> v_min[dim];
    bool first = true;

    // safety -- дополнительный запас поверх встроенных 0.8 в формулах
    // dt_limit (fiber/PLAN.md §3.7): здесь ДВЕ независимые силы (bond и
    // contact) суммируются в одну схему интегрирования одновременно,
    // формулы dt_limit каждой физики выведены без учёта другой -- запас
    // компенсирует это (аналогично 0.5*dt_limit в
    // test_fiber_bond.cpp::test_momentum_conservation).
    double safety = 0.5;
    double last_max_penetration = 0.0; // диагностика последнего шага (a)

    double step(fiber::FiberState& fb, const fiber::FiberTopology& topo, const std::vector<fiber::FiberMaterial>& mats)
    {
        const std::size_t n = fb.size();
        if (dv_prev[0].size() != n)
        {
            for (std::size_t d = 0; d < dim; ++d)
            {
                dv_prev[d].assign(n, 0.0);
                v_min[d].assign(n, 0.0);
            }
        }

        // --- RHS в ТЕКУЩЕЙ конфигурации: dt для ЭТОГО шага берём из dt_limit,
        // вычисленного на текущем состоянии (как в src/time_integration.cpp:
        // "Общий шаг" вычисляется ПОСЛЕ RHS, из уже посчитанных dt_limit'ов). ---
        std::vector<fiber::FiberDelta> bond_delta = interactions::compute_fiber_bond(exec, fb, topo, mats, bond_buffers);
        const double bond_dt_limit = interactions::compute_fiber_bond_dt_limit(exec, fb, topo, mats);

        interactions::FiberContactSearchResult search;
        interactions::build_fiber_contact_search(exec, fb, topo, mats, search);
        interactions::compute_fiber_contact_contributions(exec, fb, mats, search.primitives, contact_buffers);
        std::vector<fiber::FiberDelta> contact_delta = rhs::reduce_contributions(exec, contact_buffers, n);
        const double contact_dt_limit = interactions::compute_fiber_contact_dt_limit(exec, fb, mats, search.primitives);

        interactions::FiberContactDiagnostics diag =
            interactions::compute_fiber_contact_diagnostics(exec, fb, mats, search.primitives);
        last_max_penetration = diag.max_penetration;

        double dt = safety * std::min(bond_dt_limit, contact_dt_limit);

        std::vector<double> dv_now[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_now[d].assign(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) dv_now[d][i] = bond_delta[i].dv[d] + contact_delta[i].dv[d];
        }

        if (!first)
        {
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i)
                {
                    v_min[d][i] = fb.vx[d][i];
                    fb.vx[d][i] += 0.5 * dt * dv_prev[d][i];
                }
            // Пересчитываем RHS ПОСЛЕ полушага предиктора скорости -- в
            // продовом цикле предиктор изменяет скорость, но НЕ координату
            // (см. pc_predict_velocity_half_step), поэтому RHS (функция
            // координат, не скоростей, для Fiber-Bond/Fiber-Fiber-Contact --
            // обе силы -- чисто позиционные) пересчитывать не нужно: dv_now
            // уже посчитан выше ДО предиктора, что соответствует
            // src/time_integration.cpp ("single_step" вызывается ДО
            // предиктора следующего шага, RHS этого шага использует
            // координаты, оставшиеся от корректора предыдущего шага).
        }

        for (std::size_t d = 0; d < dim; ++d) dv_prev[d] = dv_now[d];

        const double coeff = first ? dt * 0.5 : dt;
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                const double base = first ? fb.vx[d][i] : v_min[d][i];
                fb.vx[d][i] = base + coeff * dv_prev[d][i];
                fb.x[d][i] += dt * fb.vx[d][i];
            }
        }
        first = false;
        return dt;
    }
};

fiber::FiberMaterial make_material(double E, double nu, double d, double rho, double eps_break = 1.0)
{
    fiber::FiberMaterial m;
    m.E = E; m.nu = nu; m.d = d; m.rho = rho; m.eps_break = eps_break;
    return m;
}

void append_chain(fiber::FiberState& state, fiber::FiberTopology& topology, std::size_t base, int object_id,
                   int material_id, const std::vector<std::array<double, 3>>& beads, double mass_each)
{
    const std::size_t n = beads.size();
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t gi = base + k;
        for (int d3 = 0; d3 < 3; ++d3) state.x[d3][gi] = beads[k][d3];
        for (std::size_t dd = 0; dd < dim; ++dd) state.vx[dd][gi] = 0.0;
        state.mass[gi] = mass_each;
        state.object_id[gi] = object_id;
        state.material_id[gi] = material_id;
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

void test_e2e_cross_bounce()
{
    std::printf("test_e2e_cross_bounce:\n");

    // Материал -- намеренно "мягкий" (E ниже типичного для реальных волокон)
    // -- уменьшает k_eff Fiber-Bond и, соответственно, требуемое число шагов
    // dt_limit-ограниченного интегрирования при разумной длительности теста;
    // на физику критерия готовности (нет проникания / есть отскок) это не
    // влияет -- Hertz-контакт по-прежнему растёт как delta^1.5 и останавливает
    // сближение, просто с более мягкой, но качественно той же реакцией.
    const double E = 2.0e7, nu = 0.3, d = 1.0e-3, rho = 1200.0;
    const double l0 = 2.0e-3; // шаг бусин вдоль своей нити
    const double A = pi * d * d / 4.0;
    const double mass = rho * A * l0;

    const double z0 = 1.0e-3;  // половина начального зазора по оси сближения (z)
    const double v0 = 0.5;     // м/с, скорость каждой нити навстречу (по z)
    const int N = 7;           // бусин на нить (нечётное -- есть центральная)

    std::vector<std::array<double, 3>> beadsP, beadsQ;
    for (int k = -(N / 2); k <= N / 2; ++k)
    {
        beadsP.push_back({static_cast<double>(k) * l0, 0.0, z0});
        beadsQ.push_back({0.0, static_cast<double>(k) * l0, -z0});
    }

    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    materials.push_back(make_material(E, nu, d, rho));

    state.resize(beadsP.size() + beadsQ.size());
    topology.resize(state.size());
    append_chain(state, topology, 0, 0, 0, beadsP, mass);
    append_chain(state, topology, beadsP.size(), 1, 0, beadsQ, mass);

    const std::size_t centerP = static_cast<std::size_t>(N / 2);         // центральная бусина P
    const std::size_t centerQ = beadsP.size() + static_cast<std::size_t>(N / 2); // центральная бусина Q

    for (std::size_t i = 0; i < beadsP.size(); ++i) state.vx[2][i] = -v0;
    for (std::size_t i = 0; i < beadsQ.size(); ++i) state.vx[2][beadsP.size() + i] = v0;

    const double vz_center_initial_P = state.vx[2][centerP];
    const double vz_center_initial_Q = state.vx[2][centerQ];
    check(vz_center_initial_P < 0.0 && vz_center_initial_Q > 0.0, "самопроверка: нити изначально летят навстречу друг другу");

    FiberBondContactIntegrator integ;
    double time = 0.0;
    double max_penetration_observed = 0.0;
    double running_min_gap = std::numeric_limits<double>::max();

    // Критерий (б) читается СРАЗУ после подтверждённого отскока (зазор
    // восстановился на recovery_threshold сверх наблюдённого минимума), а НЕ
    // в конце длинного прогона: нити -- конечные упругие тела (Fiber-Bond),
    // после отскока в точке контакта волна растяжения многократно отражается
    // от свободных концов (длина нити ~ N*l0 ~ 1.2e-2 м, скорость волны
    // sqrt(E/rho) ~ 129 м/с -- время прохода нити ~ 1e-4 с, за длительный
    // прогон -- десятки отражений), и скорость ЦЕНТРАЛЬНОЙ бусины на большом
    // времени отражает СЛОЖНУЮ вторичную динамику всей упругой нити, а не
    // простой сигнал "оттолкнулись/нет" сразу у места контакта.
    const double recovery_threshold = 2.0e-4; // 0.2 мм -- восстановление зазора после минимума
    bool bounce_confirmed = false;
    double vz_bounce_P = 0.0, vz_bounce_Q = 0.0;
    double gap_at_bounce = 0.0;

    const int max_steps = 4000;
    const double max_time = 8.0e-3; // запас на сближение + отскок, если он запаздывает
    int steps_done = 0;
    for (; steps_done < max_steps && time < max_time; ++steps_done)
    {
        const double dt = integ.step(state, topology, materials);
        time += dt;
        max_penetration_observed = std::max(max_penetration_observed, integ.last_max_penetration);
        const double gap = state.x[2][centerP] - state.x[2][centerQ];
        running_min_gap = std::min(running_min_gap, gap);

        if (!bounce_confirmed && gap > running_min_gap + recovery_threshold)
        {
            bounce_confirmed = true;
            vz_bounce_P = state.vx[2][centerP];
            vz_bounce_Q = state.vx[2][centerQ];
            gap_at_bounce = gap;
            break; // читаем состояние сразу после отскока, не продолжаем дальше
        }
    }

    std::printf("    steps=%d time=%.4e s bounce_confirmed=%d\n", steps_done, time, bounce_confirmed ? 1 : 0);
    std::printf("    max_penetration_observed=%.4e m (d_eff=%.4e m)\n", max_penetration_observed, d);
    std::printf("    running_min_gap(z_centerP - z_centerQ)=%.4e m\n", running_min_gap);
    std::printf("    сразу после отскока: gap=%.4e m, vz_centerP: %.4f -> %.4f м/с;  vz_centerQ: %.4f -> %.4f м/с\n",
                gap_at_bounce, vz_center_initial_P, vz_bounce_P, vz_center_initial_Q, vz_bounce_Q);

    // --- (а) нет взаимного проникания: наблюдённое проникновение << d ------
    check(max_penetration_observed > 0.0, "контакт действительно произошёл (max_penetration_observed > 0)");
    check(max_penetration_observed < 0.5 * d,
          "проникновение остаётся малым относительно диаметра нити (< 0.5*d) -- нет 'прохождения насквозь'");

    // --- (б) отскок: скорости вдоль оси сближения меняют знак сразу после
    //     минимума зазора (до вторичной динамики отражений волны) ----------
    check(bounce_confirmed, "отскок (восстановление зазора) зафиксирован в пределах отведённого времени");
    check(vz_bounce_P > 0.0, "центральная бусина P отскочила -- vz сменил знак с - на + сразу после сближения");
    check(vz_bounce_Q < 0.0, "центральная бусина Q отскочила -- vz сменил знак с + на - сразу после сближения");
}

} // namespace

int main()
{
    test_e2e_cross_bounce();

    if (failures == 0)
        std::printf("test_e2e_fiber_contact: OK\n");
    else
        std::printf("test_e2e_fiber_contact: %d failure(s)\n", failures);

    return failures;
}
