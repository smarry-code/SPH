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
#include "interactions/fiber_contact_law.hpp"
#include "interactions/fiber_friction.hpp"
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
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> friction_buffers;
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> friction_count_buffers;
    std::vector<double> friction_bead_count; // c_k (mass splitting)
    std::vector<interactions::FiberFrictionContact> friction_contacts; // слоты трения последнего шага
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

    // Закон нормальной силы (fiber/plan-fiber-contact-types.md, этап 2).
    // По умолчанию тип 0 (Герц) -- прежний путь, все вызовы с params = {}.
    interactions::FiberContactLaw law;
    std::vector<double> bead_dt_ref;
    double last_contact_dt_limit = 0.0;  // предел шага контакта последнего шага
    double last_bond_dt_limit = 0.0;     // предел шага связей последнего шага
    std::size_t last_contact_count = 0;  // принятые контакты (геометрия) последнего шага

    // --- Режимы и диагностика для сценарных прогонов (c)/(d) (§8.10). По
    // умолчанию выключены -- прогоны bounce/(a)/(b) идут по прежнему пути. ---
    // euler=true -- Эйлер-подобная схема, согласованная с ode::EulerSolver
    // (include/ode/ode_solver.hpp: step_fiber -> euler_step_position_velocity):
    // без предиктора; x^{n+1} = x^n + dt*v^n (СТАРАЯ скорость, именно так
    // делает euler_step_position_velocity), v^{n+1} = v^n + dt*a; стадия
    // трения: base = fb.vx (v^n), h = dt (EulerSolver::fiber_velocity_update).
    bool euler = false;
    // collect_friction_info=true -- после каждого шага заполняется last_info.
    bool collect_friction_info = false;

    struct StepFrictionInfo
    {
        bool has_slots = false;            // есть активные слоты с mu*N>0 на этом шаге
        std::size_t active_count = 0, stick_count = 0, slip_count = 0, shared_bead_count = 0;
        double ut_max = 0.0;               // max по слотам ||u_t|| точки контакта по ИТОГОВЫМ скоростям шага
        double normal_term = 0.0;          // max_слоты|u_n| * max_{i,j}|n_i-n_j| (по слотам с mu*N>0; 0 при одном слоте)
        std::size_t force_slots = 0;       // число слотов с mu*N>0
        double capacity_sum = 0.0;         // sum mu*N*h*W по активным слотам с mu*N>0 (м/с)
        double impulse_y_on_Q = 0.0;       // импульс трения на нить Q (y) за шаг: sum m*h*a_f_y
        double h = 0.0;
    };
    StepFrictionInfo last_info;
    std::size_t q_first = 0, q_count = 0; // диапазон бусин нити Q (для impulse_y_on_Q); q_count=0 -- не считать

    // ||u_t|| точки контакта по ТЕКУЩИМ (итоговым) скоростям fb.vx, с весами
    // s,t и нормалью каждого активного слота с mu*N>0 ЭТОГО шага; max по слотам.
    double measure_ut_max(const fiber::FiberState& fb)
    {
        double best = 0.0;
        bool any = false;
        double un_max = 0.0;
        std::size_t nslots = 0;
        const interactions::FiberFrictionContact* first_slot = nullptr;
        double dn_max = 0.0;
        for (const auto& c : friction_contacts)
        {
            if (!c.active || !(c.mu * c.N > 0.0)) continue;
            any = true;
            ++nslots;
            if (first_slot)
            {
                // max |n_i - n_j| по парам слотов
                for (const auto& c2 : friction_contacts)
                {
                    if (&c2 == &c) break;
                    if (!c2.active || !(c2.mu * c2.N > 0.0)) continue;
                    double d2 = 0.0;
                    for (std::size_t d = 0; d < dim; ++d) d2 += (c.n[d] - c2.n[d]) * (c.n[d] - c2.n[d]);
                    dn_max = std::max(dn_max, std::sqrt(d2));
                }
            }
            else
                first_slot = &c;
            double u[dim];
            for (std::size_t d = 0; d < dim; ++d)
            {
                double vA = fb.vx[d][c.a0];
                if (c.a1 != -1) vA = (1.0 - c.s) * vA + c.s * fb.vx[d][static_cast<std::size_t>(c.a1)];
                double vB = fb.vx[d][c.b0];
                if (c.b1 != -1) vB = (1.0 - c.t) * vB + c.t * fb.vx[d][static_cast<std::size_t>(c.b1)];
                u[d] = vA - vB;
            }
            double un = 0.0;
            for (std::size_t d = 0; d < dim; ++d) un += u[d] * c.n[d];
            un_max = std::max(un_max, std::fabs(un));
            double ut2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d) ut2 += (u[d] - un * c.n[d]) * (u[d] - un * c.n[d]);
            best = std::max(best, std::sqrt(ut2));
        }
        last_info.has_slots = any;
        last_info.force_slots = nslots;
        last_info.normal_term = un_max * dn_max;
        return best;
    }

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
        // fiber/SPH.md, «Первый шаг»). RHS нитей (Fiber-Bond, нормальный
        // Fiber-Fiber-Contact) от скорости НЕ зависит (только от координат;
        // трение -- отдельная стадия ПОСЛЕ вычисления dt, см. ниже,
        // fiber/fiber_friction_projection.md, §1), поэтому порядок
        // "предиктор -> RHS" физику не определяет: он сохранён только чтобы
        // повторять порядок продового цикла. Итоговая (корректорная) скорость
        // строится от v_min (сохранённой ДО полушага), в точности как здесь и
        // в ode_solver.hpp.
        if (!first && !euler)
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
        interactions::FiberContactStepParams params;
        if (law.type != 0)
        {
            params.law = law;
            params.dt_prev = dt_prev > 0.0 ? dt_prev : 1e-7; // шаг прошлого вызова, как в продовом цикле
            interactions::compute_fiber_bead_dt_ref(fb, topo, mats, bead_dt_ref);
            params.bead_dt_ref = &bead_dt_ref;
            params.pair_dt_limit = &search.pair_dt_limit;
        }
        const bool has_friction = interactions::compute_fiber_contact_contributions(
            exec, fb, mats, search.primitives, contact_buffers, friction_contacts, params);
        std::vector<fiber::FiberDelta> contact_delta = rhs::reduce_contributions(exec, contact_buffers, n);
        const double contact_dt_limit =
            interactions::compute_fiber_contact_dt_limit(exec, fb, mats, search.primitives, params);

        interactions::FiberContactDiagnostics diag =
            interactions::compute_fiber_contact_diagnostics(exec, fb, mats, search, params);
        last_max_penetration = diag.max_penetration;
        last_contact_dt_limit = contact_dt_limit;
        last_bond_dt_limit = bond_dt_limit;
        last_contact_count = diag.contact_count;

        // Тип 0: прежнее выражение. Типы 1-3: запас safety -- только на связи;
        // предел типа 3 = 1.6*dt_prev/gamma зависит от шага, и запас поверх него
        // дал бы петлю dt -> 0.8 dt за шаг (в продовом цикле запаса нет).
        const double dt = (law.type == 0) ? safety * std::min(bond_dt_limit, contact_dt_limit)
                                          : std::min(safety * bond_dt_limit, contact_dt_limit);

        std::vector<double> dv_now[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_now[d].assign(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) dv_now[d][i] = bond_delta[i].dv[d] + contact_delta[i].dv[d];
        }

        // База и коэффициент корректора: v = v_min + dt*dv_prev (первый шаг --
        // v = v^0 + (dt/2)*dv_prev, база -- текущая (нетронутая предиктором)
        // fb.vx) -- то же правило, что pc_finalize_velocity_position/use_base.
        const double coeff = (first && !euler) ? dt * 0.5 : dt;

        // --- Стадия трения нить-нить (src/time_integration.cpp: после dt, до
        // advance; fiber/fiber_friction_projection.md, §5) -- проекция по
        // ровно этой базе и коэффициенту. Ускорение трения добавляется в
        // dv_now, поэтому в dv_prev (f^{n-1} следующего предиктора) попадает
        // полное применённое ускорение. При mu_ff=0 стадия не выполняется. ---
        last_info = StepFrictionInfo{};
        last_info.h = coeff;
        if (has_friction)
        {
            interactions::compute_fiber_friction_bead_counts(exec, n, friction_contacts, friction_count_buffers,
                                                              friction_bead_count);
            ode::FiberVelocityUpdate upd{};
            for (std::size_t d = 0; d < dim; ++d)
                upd.base[d] = euler ? &fb.vx[d] : (first ? &fb.vx[d] : &v_min[d]);
            upd.h = coeff;
            if (collect_friction_info)
            {
                // Диагностика -- по dv ДО трения (как и сама стадия), только чтение.
                const interactions::FiberFrictionDiagnostics fd =
                    interactions::compute_fiber_friction_diagnostics(exec, fb, upd, dv_now, friction_contacts,
                                                                        friction_bead_count);
                last_info.active_count = fd.active_count;
                last_info.stick_count = fd.stick_count;
                last_info.slip_count = fd.slip_count;
                last_info.shared_bead_count = fd.shared_bead_count;
                for (const auto& c : friction_contacts)
                {
                    if (!c.active || !(c.mu * c.N > 0.0)) continue;
                    const double W = interactions::fiber_contact_side_inv_mass(fb, c.a0, c.a1, c.s) +
                                     interactions::fiber_contact_side_inv_mass(fb, c.b0, c.b1, c.t);
                    last_info.capacity_sum += c.mu * c.N * coeff * W;
                }
            }
            interactions::compute_fiber_friction_contributions(exec, fb, upd, dv_now, friction_contacts,
                                                                friction_bead_count, friction_buffers);
            const std::vector<fiber::FiberDelta> friction_delta = rhs::reduce_contributions(exec, friction_buffers, n);
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i) dv_now[d][i] += friction_delta[i].dv[d];
            for (std::size_t i = q_first; i < q_first + q_count; ++i)
                last_info.impulse_y_on_Q += fb.mass[i] * coeff * friction_delta[i].dv[1];
        }

        for (std::size_t d = 0; d < dim; ++d) dv_prev[d] = dv_now[d];

        if (euler)
        {
            // ode::euler_step_position_velocity: x += dt*v^n, затем v += dt*a.
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i) fb.x[d][i] += dt * fb.vx[d][i];
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i) fb.vx[d][i] += coeff * dv_prev[d][i];
        }
        else
        {
            for (std::size_t d = 0; d < dim; ++d)
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    const double base = first ? fb.vx[d][i] : v_min[d][i];
                    fb.vx[d][i] = base + coeff * dv_prev[d][i];
                    fb.x[d][i] += dt * fb.vx[d][i];
                }
            }
        }

        if (collect_friction_info) last_info.ut_max = measure_ut_max(fb);
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
//   (a) mu=0,   v_t0=2.0 -- трения нет (§12: "контакт не должен менять
//       тангенциальный импульс через силу трения"); y-импульс передаётся
//       только нормальной силой через наклон нормали на изогнутой под
//       нагрузкой нити -- меньше, чем в (b), и ограничен грубой границей
//       (см. проверки после прогона (b));
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

// ----------------------------------------------------------------------
// run_slow_contact -- ПРОВЕРЯЕМЫЙ прогон медленного скольжения (§8.10 (c)/(d)
// fiber/fiber_friction_projection.md). Та же сцена и тот же критерий останова
// (подтверждённый отскок), что у run_cross_friction_case, но интегратор
// собирает диагностику стадии трения (collect_friction_info) и записывает по
// шагам: ||u_t|| точки контакта по ИТОГОВЫМ скоростям шага (после корректора,
// веса s,t и нормаль активных слотов ЭТОГО шага, max по слотам), stick/slip,
// ёмкость прилипания sum mu*N*h*W и импульс трения на нить Q.
// euler=true -- Эйлер-подобный интегратор (см. FiberBondContactIntegrator).
//
// Импульс трения на Q (dpy_Q) -- сумма по шагам вкладов САМОЙ стадии трения:
// sum_{i in Q} m_i*h*a^f_{y,i}, где h -- коэффициент обновления скорости
// (для PC: dt, на первом шаге dt/2; для Эйлера: dt). Это ровно приращение
// y-импульса нити Q, вызванное трением: вклад нормальной силы не смешивается
// (разность Delta(sum m vy) - нормаль потребовала бы отдельно считать
// y-компоненту нормального импульса из-за наклона нормали).
// ----------------------------------------------------------------------
struct SlowContactResult
{
    bool bounce_confirmed = false;
    int steps_done = 0;
    int contact_steps = 0;     // шагов с активными слотами (mu*N>0)
    int stick_steps = 0;       // шагов со stick_count>0
    int first_stick_step = -1; // первый шаг, где ВСЕ силовые слоты прилипли (stick>0, slip==0)
    int stick_phase_steps = 0; // подряд идущие такие шаги с первого прилипания (проверяются)
    int slip_steps = 0;        // шагов со slip_count>0
    double worst_ut_in_stick = 0.0; // max ||u_t|| по шагам stick-фазы (абсолютная, м/с)
    double worst_ratio_in_stick = 0.0; // max ||u_t||/tol_step по шагам stick-фазы
    double max_normal_term = 0.0;      // max второго члена допуска по stick-фазе
    double max_growth_excess = -1.0e300; // max(рост - tol_step) до конца stick-фазы
    int steps_one_slot = 0, steps_two_slots = 0; // шагов stick-фазы с 1 / >=2 силовыми слотами
    double ut_max_contact = 0.0;    // max ||u_t|| за весь контакт
    double max_growth = 0.0;        // max (||u_t||^k - ||u_t||^{k-1}) от начала контакта до конца stick-фазы
    double max_growth_unload = 0.0; // то же после конца stick-фазы (разгрузка: только печать)
    double cap_peak = 0.0;          // max по шагам sum(mu*N*h*W) по слотам
    double cap_total = 0.0;         // sum по шагам контакта
    std::size_t shared_max = 0;     // max shared_bead_count за контакт
    std::size_t active_max = 0;     // max active_count
    double dpy_Q = 0.0;             // импульс трения на нить Q (y), за весь прогон
    double dt_min = 1.0e9, dt_max = 0.0;
    double px0 = 0, py0 = 0, pz0 = 0, px1 = 0, py1 = 0, pz1 = 0;
    double m_total = 0.0;
    std::vector<double> ut_series; // ||u_t|| по шагам контакта (для печати)
    std::vector<int> ut_steps;
};

SlowContactResult run_slow_contact(double mu, double v_t0, double safety, bool euler)
{
    SlowContactResult res;
    CrossScene scene = build_cross_scene(mu);
    fiber::FiberState& state = scene.state;
    for (std::size_t i = 0; i < scene.n_p; ++i) state.vx[1][i] = v_t0;
    res.px0 = total_momentum(state, 0);
    res.py0 = total_momentum(state, 1);
    res.pz0 = total_momentum(state, 2);
    res.m_total = static_cast<double>(state.size()) * scene.mass;

    FiberBondContactIntegrator integ;
    integ.safety = safety;
    integ.euler = euler;
    integ.collect_friction_info = true;
    integ.q_first = scene.n_p;
    integ.q_count = scene.n_p;

    double time = 0.0;
    double running_min_gap = std::numeric_limits<double>::max();
    const double recovery_threshold = 2.0e-4;
    bool stick_phase_done = false;
    bool prev_valid = false;
    double ut_prev = 0.0;

    const int max_steps = 40000; // при safety/5 шагов в ~5 раз больше
    const double max_time = 8.0e-3;
    int k = 0;
    for (; k < max_steps && time < max_time; ++k)
    {
        const double dt = integ.step(state, scene.topology, scene.materials);
        time += dt;
        res.dt_min = std::min(res.dt_min, dt);
        res.dt_max = std::max(res.dt_max, dt);
        const auto& info = integ.last_info;
        res.dpy_Q += info.impulse_y_on_Q;

        if (info.has_slots)
        {
            ++res.contact_steps;
            res.ut_series.push_back(info.ut_max);
            res.ut_steps.push_back(k);
            res.ut_max_contact = std::max(res.ut_max_contact, info.ut_max);
            const double growth = prev_valid ? info.ut_max - ut_prev : 0.0;
            ut_prev = info.ut_max;
            prev_valid = true;
            res.cap_peak = std::max(res.cap_peak, info.capacity_sum);
            res.cap_total += info.capacity_sum;
            res.shared_max = std::max(res.shared_max, info.shared_bead_count);
            res.active_max = std::max(res.active_max, info.active_count);
            if (info.stick_count > 0) ++res.stick_steps;
            if (info.slip_count > 0) ++res.slip_steps;

            // Фаза прилипания: с первого шага, где все силовые слоты прилипли,
            // и ПОКА так (stick>0, slip==0); первый шаг иного вида -- конец.
            const bool stage_stick = info.stick_count > 0 && info.slip_count == 0;
            if (!stick_phase_done)
            {
                if (stage_stick)
                {
                    if (res.first_stick_step < 0) res.first_stick_step = k;
                    ++res.stick_phase_steps;
                    res.worst_ut_in_stick = std::max(res.worst_ut_in_stick, info.ut_max);
                    const double tol_step = 1.0e-9 * v_t0 + info.normal_term;
                    res.worst_ratio_in_stick = std::max(res.worst_ratio_in_stick, info.ut_max / tol_step);
                    res.max_normal_term = std::max(res.max_normal_term, info.normal_term);
                    if (info.force_slots >= 2) ++res.steps_two_slots; else ++res.steps_one_slot;
                }
                else if (res.first_stick_step >= 0)
                    stick_phase_done = true;
            }
            // Рост ||u_t|| за шаг: до конца stick-фазы (шаг, закончивший фазу,
            // уже относится к разгрузке) и после неё -- раздельно.
            if (prev_valid)
            {
                if (!stick_phase_done)
                {
                    res.max_growth = std::max(res.max_growth, growth);
                    res.max_growth_excess =
                        std::max(res.max_growth_excess, growth - (1.0e-9 * v_t0 + info.normal_term));
                }
                else
                    res.max_growth_unload = std::max(res.max_growth_unload, growth);
            }
        }
        else
            prev_valid = false;

        const double gap = state.x[2][scene.center_p] - state.x[2][scene.center_q];
        running_min_gap = std::min(running_min_gap, gap);
        if (gap > running_min_gap + recovery_threshold)
        {
            res.bounce_confirmed = true;
            break;
        }
    }
    res.steps_done = k;
    res.px1 = total_momentum(state, 0);
    res.py1 = total_momentum(state, 1);
    res.pz1 = total_momentum(state, 2);
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
    // Проверки передачи y-импульса при mu=0 -- ниже, после прогона (b): одна
    // из них сравнивает (a) с (b).

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

    // --- (a) mu=0: передача y-импульса нормальной силой ---------------------
    // Без трения y-импульс между нитями передаётся только нормальной силой
    // Герца через наклон нормали n (n_y != 0), когда P скользит по изогнутой
    // под нагрузкой Q. В схеме «пара сегментов -- одна сила»
    // (fiber/fiber_contact_segment_pairs.md) сила пересечения сосредоточена
    // на 1-2 бусинах, и Q продавливается острым изломом в узле (центральная
    // бусина Q прогибается ~8e-4 м при ~1-2e-4 м у соседей; d=1e-3, l0=2e-3
    // -- редкая дискретизация). P, скользя по y через этот излом, видит
    // max|n_y|~0.3, и |d vy| каждой нити ~4.7e-2 м/с. До перехода на пары
    // сегментов ~7 primitives на пересечение размазывали силу, max|n_y|~0.017,
    // |d vy|~4.8e-4 м/с, и проверка была "|d vy| < 3e-3 м/с". Это не
    // нарушение сохранения импульса (проверяется ниже отдельно, ~1e-20) и не
    // трение, а геометрия контакта дискретной нити. Поэтому проверяем:
    //   (1) при mu=0 относительная y-скорость меняется МЕНЬШЕ, чем при
    //       mu=0.2 (сейчас 0.093 против 0.241 м/с, запас ~2.6x);
    //   (2) грубая граница: |d vy| каждой нити < 0.5*|d vz_P| (средний наклон
    //       нормали меньше ~27 градусов; сейчас отношение ~0.157, запас ~3x).
    const double dvrel_a = std::fabs((a.vyP1 - a.vyQ1) - (a.vyP0 - a.vyQ0));
    const double dvrel_b = std::fabs((b.vyP1 - b.vyQ1) - (b.vyP0 - b.vyQ0));
    const double dvz_P_a = std::fabs(a.vzP1 - a.vzP0);
    std::printf("  (a) mu=0: |d vrel_y|=%.4e м/с (при mu=%.2f: %.4e);  |d vyP|=%.4e |d vyQ|=%.4e;  |d vzP|=%.4e м/с"
                " (отношение |d vyP|/|d vzP|=%.4f)\n",
                dvrel_a, mu_b, dvrel_b, std::fabs(a.vyP1 - a.vyP0), std::fabs(a.vyQ1 - a.vyQ0), dvz_P_a,
                dvz_P_a > 0.0 ? std::fabs(a.vyP1 - a.vyP0) / dvz_P_a : 0.0);
    check(dvrel_a < dvrel_b,
          "(a) mu=0: относительная y-скорость меняется меньше, чем при mu>0 (передача только через наклон нормали)");
    const double max_tilt_ratio = 0.5;
    check(std::fabs(a.vyP1 - a.vyP0) < max_tilt_ratio * dvz_P_a,
          "(a) mu=0: |d vyP| < 0.5*|d vzP| (средний наклон нормали контакта мал)");
    check(std::fabs(a.vyQ1 - a.vyQ0) < max_tilt_ratio * dvz_P_a,
          "(a) mu=0: |d vyQ| < 0.5*|d vzP| (средний наклон нормали контакта мал)");

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

    // --- (c)/(d) ПРОВЕРЯЕМОЕ медленное скольжение (§8.10) -------------------
    // Ёмкость прилипания за шаг: sum mu*N*h*W (м/с) -- максимум, на который
    // стадия трения способна изменить относительную касательную скорость
    // точки контакта за один шаг. Прилипание возможно, когда ёмкость >= |u_t|.
    // Зондирование при mu=0.2 (как в (b)) печатает ёмкость; выбор mu_slow для
    // (c)/(d) и обоснование -- в комментарии у mu_slow ниже. (a)/(b) не тронуты.
    const double v_t0_slow = v_t0_small;
    {
        SlowContactResult probe = run_slow_contact(mu_b, v_t0_slow, 0.5, false);
        std::printf("  (c0) зондирование mu=%.2f, v_t0=%.3f: ёмкость прилипания sum(mu*N*h*W): пик за шаг=%.4e, "
                    "суммарно за контакт=%.4e м/с (v_t0=%.3f); stick-шагов=%d\n",
                    mu_b, v_t0_slow, probe.cap_peak, probe.cap_total, v_t0_slow, probe.stick_steps);
    }

    auto report_slow = [&](const char* label, double mu, const SlowContactResult& r) {
        std::printf("  %s mu=%.2f: steps=%d bounce=%d dt=[%.3e,%.3e] контакт-шагов=%d (активных слотов max=%zu)\n", label,
                    mu, r.steps_done, r.bounce_confirmed ? 1 : 0, r.dt_min, r.dt_max, r.contact_steps, r.active_max);
        std::printf("      ёмкость прилипания: пик за шаг=%.4e, суммарно за контакт=%.4e м/с (v_t0=%.3f)\n", r.cap_peak,
                    r.cap_total, v_t0_slow);
        std::printf("      первый шаг прилипания=%d, шагов со stick=%d (stick-фаза подряд=%d), шагов со slip=%d\n",
                    r.first_stick_step, r.stick_steps, r.stick_phase_steps, r.slip_steps);
        std::printf("      max||u_t||=%.4e (v_t0=%.3f), max||u_t|| в stick-фазе=%.4e (порог %.1e), max рост за шаг до конца "
                    "stick-фазы=%.4e, max рост за шаг в разгрузке=%.4e, shared_bead_count max=%zu\n",
                    r.ut_max_contact, v_t0_slow, r.worst_ut_in_stick, 1.0e-9 * v_t0_slow, r.max_growth,
                    r.max_growth_unload, r.shared_max);
        std::printf("      допуск с нормалями: max ||u_t||/tol_step в stick-фазе=%.4e, max второго члена |u_n|*|dn|=%.4e, "
                    "шагов stick-фазы с 1 слотом=%d, с >=2 слотами=%d, max(рост-tol_step)=%.4e\n",
                    r.worst_ratio_in_stick, r.max_normal_term, r.steps_one_slot, r.steps_two_slots,
                    r.max_growth_excess);
        std::printf("      импульс трения на Q (y) за контакт: %.6e кг*м/с\n", r.dpy_Q);
        std::printf("      ||u_t|| по шагам контакта:");
        for (std::size_t i = 0; i < r.ut_series.size(); ++i)
        {
            if (i % 8 == 0) std::printf("\n        ");
            std::printf(" %d:%.2e", r.ut_steps[i], r.ut_series[i]);
        }
        std::printf("\n");
    };

    auto check_slow = [&](const char* label, const SlowContactResult& r) {
        const std::string L(label);
        check(r.bounce_confirmed, L + ": отскок зафиксирован");
        check(r.contact_steps > 0, L + ": контакт с трением был (есть слоты с mu*N>0)");
        check(r.first_stick_step >= 0, L + ": прилипание наступило (все силовые слоты stick хотя бы на одном шаге)");
        // Критерий с учётом нормалей: tol_step = 1e-9*v_t0 + |u_n|*max|n_i-n_j| по активным
        // слотам шага. Второй член: два слота «седла» (сегменты (2,3) s~1 и (3,4) s~0 против
        // одного сегмента B) лежат почти в одной точке, но имеют РАЗНЫЕ нормали (излом у общей
        // бусины). Трение каждого слота перпендикулярно его нормали и не может обнулить
        // касательные проекции нормального расхождения u_n в обеих плоскостях сразу; остаток
        // ~ |u_n|*|dn| (при n_x=+-eps: |u_n|*eps в x). Итерации Якоби/Гаусса-Зейделя это не
        // исправят: ограничения несовместны. Открытый вопрос задачи 1: два контакта седла почти
        // в одной точке (вероятно, удвоение нормальной силы на пологом изломе). Для одного
        // слота второй член = 0, критерий прежний 1e-9*v_t0.
        check(r.worst_ratio_in_stick <= 1.0,
              L + ": ||u_t|| <= 1e-9*v_t0 + |u_n|*max|dn| на всех шагах stick-фазы (stick>0, slip==0)");
        check(r.ut_max_contact <= v_t0_slow, L + ": max||u_t|| за контакт <= v_t0 (роста нет)");
        // Монотонность ||u_t|| -- только от начала контакта до конца stick-фазы
        // (конец -- первый шаг, где не stick>0 && slip==0). В разгрузке после
        // прилипания N падает при расхождении нитей, остальные бусины P тянут
        // точку контакта, ёмкость muN*h*W мала, и скольжение с |F_t|=muN --
        // допустимое поведение (§8.10 не требует нулевой скорости после
        // отскока). Там проверяется только max||u_t|| <= v_t0 (выше), рост за
        // шаг печатается (report_slow), но не проверяется.
        check(r.max_growth_excess <= 0.0,
              L + ": ||u_t||_n <= ||u_t||_{n-1} + tol_step_n от начала контакта до конца stick-фазы");
        const double ptol = 1.0e-12 * r.m_total * std::max(a.v0, v_t0_slow);
        check(std::fabs(r.px1 - r.px0) < ptol && std::fabs(r.py1 - r.py0) < ptol && std::fabs(r.pz1 - r.pz0) < ptol,
              L + ": суммарный импульс системы сохранился");
        std::printf("      сохранение импульса: |dpx|=%.3e |dpy|=%.3e |dpz|=%.3e (допуск %.3e)\n",
                    std::fabs(r.px1 - r.px0), std::fabs(r.py1 - r.py0), std::fabs(r.pz1 - r.pz0), ptol);
    };

    // mu в (c)/(d) НЕ увеличено: по печати зондирования (c0) ёмкость
    // прилипания в этой сцене уже достаточна -- за контакт суммарно ~0.4 м/с
    // (в 8 раз больше v_t0=0.05), а на пике нагрузки ёмкость одного шага
    // ~4.75e-2 м/с ~ v_t0, так что прилипание наступает (8 шагов со stick при
    // safety=0.5). Физически ожидаемое прилипание получается при mu=0.2, как в (b).
    const double mu_slow = mu_b;
    SlowContactResult c1 = run_slow_contact(mu_slow, v_t0_slow, 0.5, false);
    SlowContactResult c5 = run_slow_contact(mu_slow, v_t0_slow, 0.1, false);
    std::printf("  (c) PredictorCorrector, mu=%.2f, v_t0=%.3f:\n", mu_slow, v_t0_slow);
    report_slow("(c) safety   ", mu_slow, c1);
    report_slow("(c) safety/5 ", mu_slow, c5);
    check_slow("(c) safety", c1);
    check_slow("(c) safety/5", c5);
    {
        const double rel = std::fabs(c1.dpy_Q - c5.dpy_Q) / std::max(std::fabs(c1.dpy_Q), std::fabs(c5.dpy_Q));
        std::printf("      сходимость по шагу: dp_y(Q) safety=%.6e, safety/5=%.6e, расхождение=%.2f %%\n", c1.dpy_Q,
                    c5.dpy_Q, 100.0 * rel);
        check(rel <= 0.10, "(c) импульс трения на Q при safety и safety/5 отличается не более чем на 10 %");
    }

    SlowContactResult d1 = run_slow_contact(mu_slow, v_t0_slow, 0.5, true);
    SlowContactResult d5 = run_slow_contact(mu_slow, v_t0_slow, 0.1, true);
    std::printf("  (d) Эйлер-подобный интегратор, mu=%.2f, v_t0=%.3f:\n", mu_slow, v_t0_slow);
    report_slow("(d) safety   ", mu_slow, d1);
    report_slow("(d) safety/5 ", mu_slow, d5);
    check_slow("(d) safety", d1);
    check_slow("(d) safety/5", d5);
    {
        const double rel = std::fabs(d1.dpy_Q - d5.dpy_Q) / std::max(std::fabs(d1.dpy_Q), std::fabs(d5.dpy_Q));
        std::printf("      сходимость по шагу (d): dp_y(Q) safety=%.6e, safety/5=%.6e, расхождение=%.2f %%\n", d1.dpy_Q,
                    d5.dpy_Q, 100.0 * rel);
    }
}

} // namespace

// ============================================================================
// test_e2e_cross_types -- крестовая сцена для законов 0-3
// (fiber/plan-fiber-contact-types.md, шаг 2.4, п. 8). Метрики печатаются,
// критерии: нет NaN, контакт случился, энергия не растёт больше 1 % (типы
// 1-3), шаг не схлопывается. Столбец sgn -- смены знака относительной
// касательной скорости центров в прогоне с трением (mu = 0.2, v_t0 = 0.05). Для типов 1-3 ожидается потеря энергии (отскока нет, §0 п.3);
// ориентир Герца +0.6 % -- fiber/FRICTION_REPORT.md.
// ============================================================================
double cross_energy(const fiber::FiberState& fb, const fiber::FiberTopology& topo,
                    const std::vector<fiber::FiberMaterial>& mats)
{
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
        const double strain = std::sqrt(l2) - b.l0;
        const fiber::FiberMaterial& m = mats[static_cast<std::size_t>(fb.material_id[b.i])];
        if (strain > 0.0) e += 0.5 * (m.E * pi * m.d * m.d / 4.0 / b.l0) * strain * strain;
    }
    return e;
}

struct CrossTypeMetrics
{
    bool finite = true;
    int steps = 0;
    double max_pen_over_d = 0.0;
    double min_dt = std::numeric_limits<double>::max();
    double dt0 = 0.0;
    int limited_by_contact = 0;
    double e0 = 0.0, e_max = 0.0, e_end = 0.0;
    double vz_p_end = 0.0, vz_q_end = 0.0;
    double max_dt_ratio = 0.0;
    int sign_changes = 0; // смены знака относительной касательной скорости центров (прогон с трением)
};

CrossTypeMetrics run_cross_type(int type, double mu_ff, double v_t0)
{
    CrossScene scene = build_cross_scene(mu_ff);
    for (std::size_t i = 0; i < scene.n_p; ++i) scene.state.vx[1][i] = v_t0;
    FiberBondContactIntegrator integ;
    integ.law.type = type;
    CrossTypeMetrics m;
    m.e0 = cross_energy(scene.state, scene.topology, scene.materials);
    m.e_max = m.e0;
    double time = 0.0, dt_prev = 0.0;
    int last_sign = 0;
    const double t_end = 8.0e-3;
    while (time < t_end && m.steps < 40000)
    {
        const double dt = integ.step(scene.state, scene.topology, scene.materials);
        if (m.steps == 0) m.dt0 = dt;
        time += dt;
        ++m.steps;
        m.min_dt = std::min(m.min_dt, dt);
        if (dt_prev > 0.0) m.max_dt_ratio = std::max(m.max_dt_ratio, dt / dt_prev);
        dt_prev = dt;
        if (integ.last_contact_count > 0 && integ.last_contact_dt_limit < integ.safety * integ.last_bond_dt_limit)
            ++m.limited_by_contact;
        m.max_pen_over_d = std::max(m.max_pen_over_d, integ.last_max_penetration / scene.d);
        const double e = cross_energy(scene.state, scene.topology, scene.materials);
        m.e_max = std::max(m.e_max, e);
        if (!std::isfinite(e)) { m.finite = false; break; }
        if (v_t0 > 0.0 && integ.last_contact_count > 0)
        {
            const double u = scene.state.vx[1][scene.center_p] - scene.state.vx[1][scene.center_q];
            const int sg = (u > 0.0) - (u < 0.0);
            if (sg != 0 && last_sign != 0 && sg != last_sign) ++m.sign_changes;
            if (sg != 0) last_sign = sg;
        }
    }
    m.e_end = cross_energy(scene.state, scene.topology, scene.materials);
    m.vz_p_end = scene.state.vx[2][scene.center_p];
    m.vz_q_end = scene.state.vx[2][scene.center_q];
    return m;
}

void test_e2e_cross_types()
{
    std::printf("test_e2e_cross_types:\n");
    static const char* names[] = {"Hertz", "pinball", "Swegle", "impulse", "Hertz+damp"};
    std::printf("    %-8s %6s %9s %10s %7s %9s %9s %9s %9s %8s %6s\n", "type", "steps", "max d/d", "min dt",
                "ct-lim", "E_end/E0", "E_max/E0", "vzP_end", "vzQ_end", "dtN/dtP", "sgn");
    for (int type = 0; type <= 4; ++type)
    {
        const CrossTypeMetrics m = run_cross_type(type, 0.0, 0.0);
        const CrossTypeMetrics mf = run_cross_type(type, 0.2, 0.05); // трение, малая касательная скорость
        std::printf("    %-8s %6d %9.3e %10.3e %6.1f%% %9.4f %9.4f %9.4f %9.4f %8.4f %6d\n", names[type], m.steps,
                    m.max_pen_over_d, m.min_dt, 100.0 * m.limited_by_contact / std::max(1, m.steps), m.e_end / m.e0,
                    m.e_max / m.e0, m.vz_p_end, m.vz_q_end, m.max_dt_ratio, mf.sign_changes);
        const std::string tag = names[type];
        check(m.finite && mf.finite, tag + ": нет NaN/Inf (без трения и с трением)");
        check(m.max_pen_over_d > 0.0, tag + ": контакт произошёл");
        // Герц -- только справочно (прежнее поведение, рост энергии ~1 % в
        // этой сцене, fiber/FRICTION_REPORT.md); типы 1-3 -- критерий.
        if (type != 0) check(m.e_max <= 1.01 * m.e0, tag + ": энергия не растёт больше 1 %");
        check(m.min_dt > 0.1 * m.dt0 && mf.min_dt > 0.1 * mf.dt0, tag + ": шаг не схлопывается");
    }
}

int main()
{
    test_e2e_cross_bounce();
    test_e2e_cross_friction();
    test_e2e_cross_types();

    if (failures == 0)
        std::printf("test_e2e_fiber_contact: OK\n");
    else
        std::printf("test_e2e_fiber_contact: %d failure(s)\n", failures);

    return failures;
}
