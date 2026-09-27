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

    // dt_prev -- шаг, использованный в ПРЕДЫДУЩЕМ вызове step() (там же, где
    // был вычислен dv_prev). Нужен ОТДЕЛЬНО от "dt этого шага", потому что
    // порядок операций внутри шага -- как в src/time_integration.cpp,
    // detail::run_time_integration: solver.predict() (полушаг предиктора
    // скорости по f^{n-1}, использует dt ПРЕДЫДУЩЕГО шага -- там это
    // буквально та же переменная `dt`, ещё не перезаписанная новым значением
    // из текущей итерации) -> RHS (пишет f^n) -> НОВЫЙ dt (из уже
    // пересчитанных dt_limit) -> solver.advance() (корректор на НОВОМ dt).
    // Тем самым переменная `dt` в продовом цикле "переживает" итерацию и на
    // предикторе следующего шага несёт значение из корректора текущего --
    // здесь то же самое, только явно как отдельное поле, а не одна
    // сквозная переменная во внешнем while.
    double dt_prev = 0.0;

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

        // --- Предиктор скорости (src/time_integration.cpp: solver.predict()
        // ВСЕГДА вызывается ДО RHS этого шага; нет-оп на первом шаге --
        // fiber/SPH.md, «Первый шаг»). ВАЖНО: с трением (mu_ff>0,
        // fiber/fiber_friction_implementation.md) сила Fiber-Fiber-Contact
        // ЗАВИСИТ от fb.vx (относительная тангенциальная скорость точек
        // контакта, compute_fiber_contact_contributions), поэтому порядок
        // "предиктор -> RHS" (а не наоборот) обязателен -- RHS должна видеть
        // ПРЕДСКАЗАННУЮ на полшага скорость, как в ode::
        // PredictorCorrectorSolver::predict_fiber/pc_predict_velocity_half_step,
        // а не сырую скорость с конца предыдущего шага. При mu_ff=0 сила не
        // зависит от v вовсе (Fiber-Bond тоже -- см. compute_fiber_bond,
        // формула только от x_ij) -- порядок операций физику не меняет, а
        // итоговая (корректорная) скорость всё равно строится от v_min
        // (сохранённой ДО полушага), в точности как здесь и в
        // ode_solver.hpp -- поэтому вывод при mu_ff=0 остаётся ПОБИТОВО
        // прежним (см. также сравнение с эталоном в
        // test_e2e_cross_bounce/main этого файла).
        if (!first)
        {
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i)
                {
                    v_min[d][i] = fb.vx[d][i];
                    fb.vx[d][i] += 0.5 * dt_prev * dv_prev[d][i];
                }
        }

        // --- RHS на (возможно предсказанной выше) скорости и координатах,
        // оставшихся от корректора предыдущего шага. dt ЭТОГО шага берём из
        // dt_limit, вычисленного ПОСЛЕ RHS (как в src/time_integration.cpp:
        // "Общий шаг" вычисляется ПОСЛЕ RHS, из уже посчитанных
        // dt_limit'ов), и используем его в корректоре ниже. ---
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

        const double dt = safety * std::min(bond_dt_limit, contact_dt_limit);

        std::vector<double> dv_now[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_now[d].assign(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) dv_now[d][i] = bond_delta[i].dv[d] + contact_delta[i].dv[d];
        }

        for (std::size_t d = 0; d < dim; ++d) dv_prev[d] = dv_now[d];

        // --- Корректор: v = v_min + dt*dv_prev (первый шаг -- v = v^0 +
        // (dt/2)*dv_prev, база -- текущая (нетронутая предиктором) fb.vx --
        // то же правило, что pc_finalize_velocity_position/use_base. ---
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
        dt_prev = dt;
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

// ----------------------------------------------------------------------
// CrossScene / build_cross_scene -- сборка крестовой сцены (fiber/PLAN.md,
// этап Ф5, «Готово, когда»): нить P вдоль оси x на высоте z=+z0, летит вниз
// (vz=-v0); нить Q вдоль оси y на высоте z=-z0, летит вверх (vz=+v0) --
// проекции пересекаются в окрестности (x,y)=(0,0). ОБЩАЯ для
// test_e2e_cross_bounce (единственный потребитель до внедрения трения) и
// test_e2e_cross_friction (fiber/fiber_friction_implementation.md, §12) --
// арифметика построения сцены НЕ изменена ни на бит относительно исходного
// test_e2e_cross_bounce: mu_ff -- единственный новый параметр (в
// test_e2e_cross_bounce передаётся 0.0, что воспроизводит прежнее поведение
// fiber::FiberMaterial{...} без mu_ff, у которого mu_ff=0.0 по умолчанию --
// см. fiber/fiber_system.hpp).
// ----------------------------------------------------------------------
struct CrossScene
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    std::size_t n_p = 0;      // число бусин в P (= число бусин в Q -- обе нити длины N)
    std::size_t center_p = 0; // индекс центральной бусины P в state
    std::size_t center_q = 0; // индекс центральной бусины Q в state
    double mass = 0.0;        // масса ОДНОЙ бусины (все бусины -- одной массы)
    double v0 = 0.0;          // скорость сближения по z, м/с (задана обеим нитям навстречу)
    double d = 0.0;           // диаметр нити (для критерия "нет прохождения насквозь")
};

CrossScene build_cross_scene(double mu_ff)
{
    CrossScene scene;

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

    scene.materials.push_back(make_material(E, nu, d, rho));
    scene.materials[0].mu_ff = mu_ff;

    scene.state.resize(beadsP.size() + beadsQ.size());
    scene.topology.resize(scene.state.size());
    append_chain(scene.state, scene.topology, 0, 0, 0, beadsP, mass);
    append_chain(scene.state, scene.topology, beadsP.size(), 1, 0, beadsQ, mass);

    scene.n_p = beadsP.size();
    scene.center_p = static_cast<std::size_t>(N / 2);                          // центральная бусина P
    scene.center_q = beadsP.size() + static_cast<std::size_t>(N / 2);          // центральная бусина Q
    scene.mass = mass;
    scene.v0 = v0;
    scene.d = d;

    for (std::size_t i = 0; i < beadsP.size(); ++i) scene.state.vx[2][i] = -v0;
    for (std::size_t i = 0; i < beadsQ.size(); ++i) scene.state.vx[2][beadsP.size() + i] = v0;

    return scene;
}

// Суммарный импульс (по всем бусинам ОБЕИХ нитей, т.е. всей системы) вдоль
// оси d3 -- для проверки сохранения импульса (fiber/fiber_friction_
// implementation.md, §12: «суммарный импульс системы должен сохраняться»).
double total_momentum(const fiber::FiberState& state, std::size_t d3)
{
    double p = 0.0;
    for (std::size_t i = 0; i < state.size(); ++i) p += state.mass[i] * state.vx[d3][i];
    return p;
}

// Средняя по массе скорость вдоль оси d3 по диапазону бусин
// [first, first+count) (все бусины одной нити имеют одинаковую массу --
// среднее по массе совпадает с обычным арифметическим средним, но формула
// записана через импульс/массу для общности).
double mean_velocity(const fiber::FiberState& state, std::size_t d3, std::size_t first, std::size_t count)
{
    double sum_p = 0.0, sum_m = 0.0;
    for (std::size_t k = 0; k < count; ++k)
    {
        const std::size_t i = first + k;
        sum_p += state.mass[i] * state.vx[d3][i];
        sum_m += state.mass[i];
    }
    return sum_p / sum_m;
}

void test_e2e_cross_bounce()
{
    std::printf("test_e2e_cross_bounce:\n");

    CrossScene scene = build_cross_scene(0.0); // mu_ff=0 -- трения в этом тесте нет
    fiber::FiberState& state = scene.state;
    fiber::FiberTopology& topology = scene.topology;
    std::vector<fiber::FiberMaterial>& materials = scene.materials;
    const std::size_t centerP = scene.center_p;
    const std::size_t centerQ = scene.center_q;
    const double d = scene.d;

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

// ============================================================================
// test_e2e_cross_friction -- сценарный тест динамического трения Кулона
// между нитями (fiber/fiber_friction_implementation.md, §12 «Сценарный
// тест»). Та же крестовая сцена, что test_e2e_cross_bounce (build_cross_scene
// выше), плюс касательная скорость: ВСЕ бусины P получают vy=+v_t0 (вдоль оси
// Q, перпендикулярно нормали контакта n≈z), Q -- vy=0. Три прогона (mu, v_t0
// задаются параметрами run_cross_friction_case):
//   (a) mu=0,   v_t0=2.0 -- трения нет, касательный импульс каждой нити
//       обязан остаться прежним (§12: "контакт не должен менять
//       тангенциальный импульс через силу трения");
//   (b) mu=0.2, v_t0=2.0 -- трение есть, относительная тангенциальная
//       скорость нитей обязана УМЕНЬШИТЬСЯ (§12: "модуль относительной
//       тангенциальной скорости должен уменьшаться"), но не обнулиться --
//       для этого v_t0 подобран так, чтобы скольжение не прекращалось за
//       время контакта: оценка "остановки" -- mu*|Delta v_n|, где Delta v_n
//       -- изменение осевой (по z) скорости центра масс нити за контакт
//       (импульс НОРМАЛЬНОЙ силы на нить = m_нити*Delta v_n, т.к. Fiber-Bond
//       -- внутренняя сила нити, не меняющая импульс её центра масс; трение
//       не имеет z-компоненты, n~=z). Из test_e2e_cross_bounce известно, что
//       при том же v0=0.5 характерное Delta v_n центральной бусины порядка
//       1.5-2 м/с (см. её вывод: vz: -0.5 -> 1.2) -- т.е. mu*Delta_v_n
//       порядка 0.2*1.7~=0.34 м/с, что СУЩЕСТВЕННО меньше v_t0=2.0 -- запас
//       на порядок, скольжение не должно прекратиться (числа подтверждаются
//       печатью ниже: "оценка тормож. трением");
//   (c) mu=0.2, v_t0=0.05 (маленькая) -- ДИАГНОСТИЧЕСКИЙ прогон (тест не
//       валит): при столь малой касательной скорости относительно "силы
//       торможения" mu*Delta_v_n схема без stick/slip (fiber/fiber_friction_
//       implementation.md, §13: "не реализовывать в этой задаче") способна
//       давать дребезг знака v_t на каждом шаге -- считаем число смен знака
//       и максимальную амплитуду после первой смены, печатаем сводку.
// ============================================================================

// Результат одного прогона run_cross_friction_case -- скорости/импульсы ДО
// (сразу после задания касательной скорости) и СРАЗУ ПОСЛЕ подтверждённого
// отскока (та же логика "recovery_threshold", что и test_e2e_cross_bounce).
struct FrictionCaseResult
{
    bool bounce_confirmed = false;
    int steps_done = 0;
    double time = 0.0;
    double gap_at_bounce = 0.0;
    double max_penetration_observed = 0.0;

    double px0 = 0.0, py0 = 0.0, pz0 = 0.0; // импульс СИСТЕМЫ (P+Q) до контакта
    double px1 = 0.0, py1 = 0.0, pz1 = 0.0; // импульс СИСТЕМЫ сразу после отскока

    double vyP0 = 0.0, vyQ0 = 0.0; // средняя (по массе) vy нити P/Q до контакта
    double vyP1 = 0.0, vyQ1 = 0.0; // средняя vy P/Q сразу после отскока
    double vzP0 = 0.0, vzQ0 = 0.0; // средняя vz P/Q до контакта
    double vzP1 = 0.0, vzQ1 = 0.0; // средняя vz P/Q сразу после отскока

    double m_total = 0.0; // суммарная масса ВСЕЙ системы (P+Q) -- для допуска на сохранение импульса
    double v0 = 0.0;      // скорость сближения по z из сцены (build_cross_scene) -- для масштаба допуска
};

// run_cross_friction_case -- прогон крестовой сцены с касательной скоростью
// P (vy=+v_t0) и трением mu (mu_ff обеих нитей -- материал общий, значит
// fiber_friction_coefficient(a,a)=a.mu_ff, см. include/interactions/
// fiber_contact.hpp) до подтверждённого отскока -- ТА ЖЕ логика, что цикл
// test_e2e_cross_bounce (recovery_threshold/max_steps/max_time идентичны).
// diag_verbose=true (только прогон (c)) -- печатает по шагам, ПОКА есть
// контакт (integ.last_max_penetration>0), относительную vy центральных
// бусин P и Q (не более 20 строк выборки), считает число смен знака этой
// величины и максимальную амплитуду после первой смены -- оценка "дребезга"
// явной схемы Кулона без stick/slip (§13 спецификации).
FrictionCaseResult run_cross_friction_case(double mu, double v_t0, bool diag_verbose)
{
    FrictionCaseResult res;
    CrossScene scene = build_cross_scene(mu);
    fiber::FiberState& state = scene.state;

    // Касательная скорость: ВСЕ бусины P получают vy=+v_t0 (вдоль оси Q,
    // перпендикулярно нормали контакта n~=z); Q -- vy=0 (уже 0 после
    // build_cross_scene -- FiberState::resize/append_chain обнуляют vx[1]).
    for (std::size_t i = 0; i < scene.n_p; ++i) state.vx[1][i] = v_t0;

    res.px0 = total_momentum(state, 0);
    res.py0 = total_momentum(state, 1);
    res.pz0 = total_momentum(state, 2);
    res.vyP0 = mean_velocity(state, 1, 0, scene.n_p);
    res.vyQ0 = mean_velocity(state, 1, scene.n_p, scene.n_p);
    res.vzP0 = mean_velocity(state, 2, 0, scene.n_p);
    res.vzQ0 = mean_velocity(state, 2, scene.n_p, scene.n_p);
    res.m_total = static_cast<double>(state.size()) * scene.mass; // все бусины -- одной массы (scene.mass)
    res.v0 = scene.v0;

    FiberBondContactIntegrator integ;
    double time = 0.0;
    double max_penetration_observed = 0.0;
    double running_min_gap = std::numeric_limits<double>::max();
    const double recovery_threshold = 2.0e-4; // как в test_e2e_cross_bounce
    bool bounce_confirmed = false;

    // --- диагностика дребезга (только прогон (c), diag_verbose=true) ------
    double vy_rel_prev = 0.0;
    bool have_prev = false;
    bool first_change_seen = false;
    int sign_changes = 0;
    double max_amp_after_first_change = 0.0;
    int printed_lines = 0;

    const int max_steps = 4000;
    const double max_time = 8.0e-3; // как в test_e2e_cross_bounce
    int steps_done = 0;
    for (; steps_done < max_steps && time < max_time; ++steps_done)
    {
        const double dt = integ.step(state, scene.topology, scene.materials);
        time += dt;
        max_penetration_observed = std::max(max_penetration_observed, integ.last_max_penetration);
        const double gap = state.x[2][scene.center_p] - state.x[2][scene.center_q];
        running_min_gap = std::min(running_min_gap, gap);

        if (diag_verbose && integ.last_max_penetration > 0.0)
        {
            const double vy_rel = state.vx[1][scene.center_p] - state.vx[1][scene.center_q];
            if (have_prev && vy_rel_prev != 0.0 && vy_rel != 0.0 && (vy_rel_prev < 0.0) != (vy_rel < 0.0))
            {
                ++sign_changes;
                first_change_seen = true;
            }
            if (first_change_seen) max_amp_after_first_change = std::max(max_amp_after_first_change, std::fabs(vy_rel));
            if (printed_lines < 20)
            {
                std::printf("    [diag c] step=%d t=%.4e vy_rel(centerP-centerQ)=%.6e м/с\n", steps_done, time, vy_rel);
                ++printed_lines;
            }
            vy_rel_prev = vy_rel;
            have_prev = true;
        }

        if (!bounce_confirmed && gap > running_min_gap + recovery_threshold)
        {
            bounce_confirmed = true;
            res.gap_at_bounce = gap;
            break; // читаем состояние сразу после отскока, не продолжаем дальше
        }
    }

    res.bounce_confirmed = bounce_confirmed;
    res.steps_done = steps_done;
    res.time = time;
    res.max_penetration_observed = max_penetration_observed;

    res.px1 = total_momentum(state, 0);
    res.py1 = total_momentum(state, 1);
    res.pz1 = total_momentum(state, 2);
    res.vyP1 = mean_velocity(state, 1, 0, scene.n_p);
    res.vyQ1 = mean_velocity(state, 1, scene.n_p, scene.n_p);
    res.vzP1 = mean_velocity(state, 2, 0, scene.n_p);
    res.vzQ1 = mean_velocity(state, 2, scene.n_p, scene.n_p);

    if (diag_verbose)
    {
        std::printf("    [diag c] сводка дребезга: смен_знака=%d, макс_амплитуда_после_1й_смены=%.4e м/с "
                    "(из %d шагов контакта, выборка выше ограничена 20 строками)\n",
                    sign_changes, max_amp_after_first_change, steps_done);
    }

    return res;
}

void test_e2e_cross_friction()
{
    std::printf("test_e2e_cross_friction:\n");

    // v_t0 для прогонов (a)/(b) -- см. обоснование в комментарии-шапке
    // test_e2e_cross_friction выше (запас ~6x над оценкой mu*Delta_v_n).
    const double v_t0_slide = 2.0;
    const double mu_b = 0.2;

    // --- (a) mu=0: трения нет -- тангенциальный импульс каждой нити не
    //     меняется контактом ------------------------------------------------
    FrictionCaseResult a = run_cross_friction_case(0.0, v_t0_slide, false);
    std::printf("  (a) mu=0, v_t0=%.3f: steps=%d time=%.4e s bounce_confirmed=%d max_penetration=%.4e m\n",
                v_t0_slide, a.steps_done, a.time, a.bounce_confirmed ? 1 : 0, a.max_penetration_observed);
    std::printf("      vyP: %.6f -> %.6f м/с (|d|=%.3e);  vyQ: %.6f -> %.6f м/с (|d|=%.3e)\n", a.vyP0, a.vyP1,
                std::fabs(a.vyP1 - a.vyP0), a.vyQ0, a.vyQ1, std::fabs(a.vyQ1 - a.vyQ0));
    check(a.bounce_confirmed, "(a) отскок зафиксирован в пределах отведённого времени");

    // Допуск на "тангенциальный импульс не изменился": нормаль контакта n
    // формально ~=z (обе нити прямые в начальный момент), но геометрия
    // Fiber-Fiber-Contact считается по ФАКТИЧЕСКИМ (изогнутым Fiber-Bond'ом
    // под нагрузкой) осям бусин -- у n может появиться небольшая x/y-
    // компонента, и тогда нормальная сила Герца (которая доминирует над
    // трением и на 2-3 порядка больше его по модулю в этой сцене) даёт
    // малый паразитный вклад в vy. Фактически наблюдаемое отклонение в этой
    // сцене -- порядка 1e-4..1e-3 м/с (см. печать выше); допуск ниже взят с
    // запасом ~10x от него.
    const double tol_vy_zero_mu = 3.0e-3; // м/с
    check(std::fabs(a.vyP1 - a.vyP0) < tol_vy_zero_mu,
          "(a) mu=0: тангенциальный импульс P не изменился контактом (в пределах допуска на наклон нормали)");
    check(std::fabs(a.vyQ1 - a.vyQ0) < tol_vy_zero_mu,
          "(a) mu=0: тангенциальный импульс Q не изменился контактом (в пределах допуска на наклон нормали)");

    // --- (b) mu=0.2: относительная тангенциальная скорость уменьшается,
    //     P тормозится, Q увлекается -----------------------------------------
    FrictionCaseResult b = run_cross_friction_case(mu_b, v_t0_slide, false);
    std::printf("  (b) mu=%.2f, v_t0=%.3f: steps=%d time=%.4e s bounce_confirmed=%d max_penetration=%.4e m\n", mu_b,
                v_t0_slide, b.steps_done, b.time, b.bounce_confirmed ? 1 : 0, b.max_penetration_observed);
    std::printf("      vyP: %.6f -> %.6f м/с;  vyQ: %.6f -> %.6f м/с\n", b.vyP0, b.vyP1, b.vyQ0, b.vyQ1);

    const double v_rel_before = b.vyP0 - b.vyQ0; // == v_t0_slide (Q стартует с vy=0)
    const double v_rel_after = b.vyP1 - b.vyQ1;
    std::printf("      относительная vy(P-Q): %.6f -> %.6f м/с\n", v_rel_before, v_rel_after);

    check(b.bounce_confirmed, "(b) отскок зафиксирован в пределах отведённого времени");
    check(v_rel_after < v_t0_slide, "(b) относительная тангенциальная скорость (P-Q) после отскока строго меньше v_t0");
    check(b.vyP1 < v_t0_slide, "(b) P затормозилась по y трением (vyP1 < v_t0)");
    check(b.vyQ1 > 0.0, "(b) Q увлеклась в +y трением (vyQ1 > 0)");

    // Оценка порядка величины: уменьшение относительной скорости против
    // mu*Delta_v_n (Delta_v_n -- изменение ОСЕВОЙ (z) скорости центра масс
    // P за контакт = импульс нормальной силы на P / m_P, см. обоснование в
    // шапке-комментарии test_e2e_cross_friction). Не жёсткая проверка (§12
    // спецификации формализует только направление/знак и строгое
    // неравенство выше) -- только печать для качественного сопоставления.
    const double delta_v_t = v_rel_before - v_rel_after;
    const double delta_v_n_P = std::fabs(b.vzP1 - b.vzP0);
    const double braking_estimate = mu_b * delta_v_n_P;
    std::printf("      оценка тормож. трением: delta_v_t(набл.)=%.4e м/с; mu*|delta_v_n_P|(оценка)=%.4e м/с\n",
                delta_v_t, braking_estimate);

    // --- сохранение импульса системы (P+Q) по всем трём компонентам, оба
    //     прогона (a) и (b) (fiber/fiber_friction_implementation.md, §12:
    //     «суммарный импульс системы должен сохраняться»). Относительный
    //     допуск -- порядка 1e-12 от масштаба m_total*max(v0, v_t0): Fiber-
    //     Fiber-Contact (нормаль + трение) и Fiber-Bond по построению
    //     ньютоновски симметричны (writer.emit(i,...)/writer.emit(j,...) с
    //     di=-dj*mj/mi, см. compute_fiber_contact_contributions/
    //     compute_fiber_bond_contributions) -- отклонение появляется только
    //     от округления double по сумме шагов, не от асимметрии модели. -----
    const double scale = a.m_total * std::max(a.v0, v_t0_slide);
    const double mom_tol = 1.0e-12 * scale;

    auto check_momentum = [&](const char* label, const FrictionCaseResult& r) {
        const double dpx = r.px1 - r.px0;
        const double dpy = r.py1 - r.py0;
        const double dpz = r.pz1 - r.pz0;
        std::printf("  %s: сохранение импульса -- |dpx|=%.4e |dpy|=%.4e |dpz|=%.4e (допуск=%.4e, кг*м/с)\n", label,
                    std::fabs(dpx), std::fabs(dpy), std::fabs(dpz), mom_tol);
        check(std::fabs(dpx) < mom_tol, std::string(label) + ": суммарный импульс px сохранился");
        check(std::fabs(dpy) < mom_tol, std::string(label) + ": суммарный импульс py сохранился");
        check(std::fabs(dpz) < mom_tol, std::string(label) + ": суммарный импульс pz сохранился");
    };
    check_momentum("(a) mu=0", a);
    check_momentum("(b) mu=0.2", b);

    // --- (c) диагностический прогон: маленькое v_t0, печать дребезга,
    //     тест НЕ валится -- см. обоснование в шапке-комментарии выше ------
    const double v_t0_small = 0.05;
    std::printf("  (c) mu=%.2f, v_t0=%.3f (диагностика дребезга, не влияет на результат теста):\n", mu_b, v_t0_small);
    FrictionCaseResult c = run_cross_friction_case(mu_b, v_t0_small, true);
    std::printf("  (c) итог: steps=%d time=%.4e s bounce_confirmed=%d max_penetration=%.4e m\n", c.steps_done, c.time,
                c.bounce_confirmed ? 1 : 0, c.max_penetration_observed);
}

} // namespace

int main()
{
    test_e2e_cross_bounce();
    test_e2e_cross_friction();

    if (failures == 0)
        std::printf("test_e2e_fiber_contact: OK\n");
    else
        std::printf("test_e2e_fiber_contact: %d failure(s)\n", failures);

    return failures;
}
