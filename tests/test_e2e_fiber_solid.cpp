// tests/test_e2e_fiber_solid.cpp
// E2E-критерии готовности этапа Ф6 (notes/fiber-model/fiber-plan.md, «Fiber-Solid-Contact»,
// «Готово, когда»): «тесты проходят, e2e без проникания нити в пластину».
//
// Два сценария:
//
//  (А) smoke одного Solid в 3D БЕЗ нитей -- регрессионная защита: подключение
//      Ф6 (interactions/fiber_solid_contact.hpp, вызов из
//      src/time_integration.cpp, перегрузка fiber::Model&) не должно ломать
//      обычный Solid-путь, когда нитей нет. Кейс -- уменьшенная 3D-геометрия
//      BC2 (tests/bc2_mini3d_geometry.hpp), та же, что использует
//      test_fiber_time_integration_equivalence.cpp. Модель нитей -- ПУСТАЯ
//      (fb.size()==0): по построению src/time_integration.cpp весь блок
//      Ф6 (как и Ф3/Ф5) обёрнут в `if (nfiber > 0)` и при пустых нитях не
//      исполняется вовсе -- ноль накладных расходов и byte-for-byte то же
//      поведение, что до Ф6. Тест -- смоук (несколько шагов без исключений,
//      состояние конечно), НЕ повторяет количественное сравнение старой и
//      новой сигнатур -- это уже делает test_fiber_time_integration_equivalence.cpp.
//
//  (Б) нить, падающая на неподвижную Solid-пластину -- прямая физика Ф6.
//      Пластина -- маленький плоский слой частиц Solid (алюминий, itype=4),
//      ЗАФИКСИРОВАННЫЙ (позиции и скорости не меняются в ходе прогона -- так
//      прямо разрешено постановкой: "нить, падающая... на НЕПОДВИЖНУЮ
//      Solid-пластину"). Нить -- прямая горизонтальная цепочка бусин над
//      пластиной с начальной скоростью вниз. Интегрируется ЛОКАЛЬНЫМ
//      PredictorCorrector-драйвером (тот же приём, что
//      tests/test_e2e_fiber_contact.cpp::FiberBondContactIntegrator) -- полный
//      time_integration(fiber::Model&,...) не нужен: пластина не должна
//      двигаться, а он обязательно прогоняет Solid-SPH и по ней тоже.
//      Драйвер суммирует Fiber-Bond (interactions/fiber_bond.hpp, держит нить
//      прямой/упругой) и Fiber-Solid-Contact (interactions/
//      fiber_solid_contact.hpp) в dv_fiber; вклад в Solid (SolidDelta)
//      вычисляется (тем же kernel'ом, что и в продовом пути), но НЕ
//      применяется к пластине -- она зафиксирована по определению сценария.
//
//  (В) трение нить-Solid (notes/fiber-contact/fiber-solid-friction-implementation.md, §7.2):
//      та же сцена (нить над пластиной) с касательной скоростью v_t0 вдоль x,
//      интегрирование РЕАЛЬНЫМ решателем (ode::PredictorCorrectorSolver) в порядке
//      шага src/time_integration.cpp, включая единую стадию трения:
//        (a) mu_fs = 0 -- касательная скорость нити контактом не меняется;
//        (b) быстрое скольжение -- торможение порядка mu*|d v_n|;
//        (c) медленное скольжение -- прилипание, сходимость импульса по dt;
//        (d) свободная (интегрируемая) пластина -- сохранение суммарного импульса.
//      Кинематическая пластина -- SolidVelocityUpdate{h=0} + solid_kinematic=true.
//
// Критерии готовности:
//   (а) нить НЕ проходит сквозь пластину -- наблюдённое проникновение
//       (delta_fs = r_f+r_s^contact-D_fs) остаётся малым относительно
//       r_f+r_s^contact, не растёт неограниченно;
//   (б) есть отклик -- скорость центральной бусины вдоль нормали к пластине
//       (z) меняет знак (падала -- отскочила) в пределах отведённого времени.
#include "fiber/fiber_system.hpp"
#include "fiber/model.hpp"
#include "input.hpp"
#include "int_force.hpp"
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

#include "bc2_mini3d_geometry.hpp"
#include "e2e_harness.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

static_assert(dim == 3, "test_e2e_fiber_solid требует сборки с SPH_DIM=3");

#ifndef SPH_TEST_DATA_DIR
#define SPH_TEST_DATA_DIR "tests/data"
#endif

namespace fs = std::filesystem;
using e2e::check;
using e2e::note;

namespace
{

// =========================================================================
// (А) smoke: один Solid в 3D, БЕЗ нитей.
// =========================================================================
void test_smoke_solid_only_3d()
{
    std::printf("test_smoke_solid_only_3d:\n");

    fs::path config_path = fs::path(SPH_TEST_DATA_DIR) / "bc2_initial_conditions.txt";
    if (!fs::exists(config_path))
    {
        check(false, "не найден конфиг теста: " + config_path.string());
        return;
    }

    SimulationConfig cfg;
    cfg.apply_defaults();
    parse_sph_config(cfg, config_path.string());
    cfg.contact_force_type = 1; // pinball -- та же ветка, что штатно включена в BC2

    fiber::Model model;
    bc2m3::Geometry geom = bc2m3::build_geometry(model.solid, bc2m3::kDxSmoke);
    const std::size_t n = geom.ntotal;
    cfg.moni_particle = geom.moni;
    cfg.maxn = static_cast<size_t>(static_cast<double>(n) * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    const std::size_t prev = model.solid.size();
    model.solid.resize(cfg.maxn);
    for (std::size_t i = prev; i < model.solid.size(); ++i) model.solid.inzializ_index(i);
    for (std::size_t i = 0; i < n; ++i) urs(cfg.mats, model.solid, i);

    // Модель нитей -- ПУСТАЯ (fiber_file не задан, fb.size()==0): весь блок
    // Ф6 в src/time_integration.cpp обёрнут в `if (nfiber > 0)` -- код Ф6
    // не исполняется вовсе (см. шапку файла).
    model.fiber.resize(0);
    model.topology.resize(0);

    double dt = compute_cfl_timestep(cfg, model.solid, n);
    double time = 0.0;
    RunGuard guard(0.0, 0.0, 0);
    guard.set_step_limit(8);

    bool threw = false;
    std::string what;
    {
        e2e::ScopedWorkDir wd("sph_e2e_fiber_solid", "smoke");
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(cfg, model, n, dt, time, bc2m3::kDxSmoke, &guard);
        }
        catch (const std::exception& ex)
        {
            threw = true;
            what = ex.what();
        }
    }

    check(!threw, "time_integration(fiber::Model&, ...) без нитей не бросает исключений" +
                       (threw ? (": " + what) : std::string()));
    check(guard.reason() == RunGuard::Reason::StepLimit, "остановлена штатно по StepLimit (не авария/застой)");
    check(guard.steps() == 8, "сделано ровно 8 шагов");

    bool all_finite = true;
    for (std::size_t i = 0; i < n && all_finite; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d)
            if (!std::isfinite(model.solid.x[d][i]) || !std::isfinite(model.solid.vx[d][i])) all_finite = false;
        if (!std::isfinite(model.solid.rho[i]) || !std::isfinite(model.solid.u[i])) all_finite = false;
    }
    check(all_finite, "состояние Solid конечно после 8 шагов (x, vx, rho, u -- нет NaN/Inf)");
    check(std::isfinite(dt) && dt > 0.0, "итоговый dt конечен и положителен");
}

// =========================================================================
// (Б) нить падает на неподвижную Solid-пластину.
// =========================================================================
fiber::FiberMaterial make_fiber_material(double E, double nu, double d, double rho, double eps_break = 1.0)
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
            b.l_break = b.l0 * 1e9; // разрыв -- Ф7, не предмет этого теста
            topology.bonds.push_back(b);
        }
    }
}

// Локальный PredictorCorrector-драйвер: Fiber-Bond + Fiber-Solid-Contact
// суммируются в dv_fiber; пластина (plate) НЕ меняется -- см. шапку файла.
// По образцу tests/test_e2e_fiber_contact.cpp::FiberBondContactIntegrator.
struct FiberOnPlateIntegrator
{
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> bond_buffers;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fs_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fs_solid_buffers; // вычислен, но отброшен
    std::vector<double> dv_prev[dim];
    std::vector<double> v_min[dim];
    bool first = true;
    double safety = 0.3; // запас поверх встроенных 0.8 в dt_limit -- две независимые физики (bond+contact) в одной схеме
    double last_max_penetration = 0.0;

    double step(fiber::FiberState& fb, const fiber::FiberTopology& topo, const std::vector<fiber::FiberMaterial>& fmats,
                const Materials& smats, const ParticleSystem& plate, std::size_t nplate)
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

        std::vector<fiber::FiberDelta> bond_delta = interactions::compute_fiber_bond(exec, fb, topo, fmats, bond_buffers);
        const double bond_dt_limit = interactions::compute_fiber_bond_dt_limit(exec, fb, topo, fmats);

        interactions::FiberSolidContactSearchResult fs_search;
        interactions::build_fiber_solid_contact_search(exec, smats, plate, nplate, fb, topo, fmats, fs_search);
        interactions::compute_fiber_solid_contact_contributions(exec, smats, plate, fb, fmats, fs_search.primitives,
                                                                  fs_fiber_buffers, fs_solid_buffers);
        std::vector<fiber::FiberDelta> fs_delta = rhs::reduce_contributions(exec, fs_fiber_buffers, n);
        const double fs_dt_limit = interactions::compute_fiber_solid_contact_dt_limit(exec, smats, plate, fb, fmats,
                                                                                        fs_search.primitives);

        double max_pen = 0.0;
        for (const auto& p : fs_search.primitives)
        {
            const interactions::FiberSolidContactGeometry g =
                interactions::compute_fiber_solid_contact_geometry(smats, plate, fb, fmats, p);
            if (g.contact) max_pen = std::max(max_pen, g.delta);
        }
        last_max_penetration = max_pen;

        double dt = safety * std::min(bond_dt_limit, fs_dt_limit);

        std::vector<double> dv_now[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_now[d].assign(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) dv_now[d][i] = bond_delta[i].dv[d] + fs_delta[i].dv[d];
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

// Плоская решётка Solid-частиц (алюминий, itype=4) -- "пластина".
// nz=1 (один слой), центр по z в 0.
void build_plate(ParticleSystem& sys, int ncols_half, double dxp, double rho)
{
    const int side = 2 * ncols_half;
    const std::size_t ntot = static_cast<std::size_t>(side) * static_cast<std::size_t>(side);
    sys.resize(ntot);
    std::size_t k = 0;
    for (int ix = -ncols_half; ix < ncols_half; ++ix)
        for (int iy = -ncols_half; iy < ncols_half; ++iy)
        {
            sys.inzializ_index(k);
            sys.x[0][k] = (ix + 0.5) * dxp;
            sys.x[1][k] = (iy + 0.5) * dxp;
            sys.x[2][k] = 0.0;
            sys.itype[k] = 4; // алюминий (materials.hpp)
            sys.body[k] = 4;
            sys.rho[k] = rho;
            sys.rho_0[k] = rho;
            sys.mass[k] = rho * dxp * dxp * dxp; // куб dxp x dxp x dxp -- толщина слоя == шагу решётки
            sys.hsml[k] = 1.0; // не используется Fiber-Solid-Contact (нет SPH-ядра)
            ++k;
        }
}

void test_fiber_falls_on_plate()
{
    std::printf("test_fiber_falls_on_plate:\n");

    // --- Пластина: неподвижная, itype=4 (Al), dxp -- шаг решётки. ----------
    const double dxp = 1.5e-3;
    const double rho_s = 2700.0;
    ParticleSystem plate;
    build_plate(plate, 5, dxp, rho_s); // 10x10 колонн, half-width 7.5e-3 м
    const std::size_t nplate = plate.size();
    Materials smats;

    const double r_s = 0.5 * dxp; // solid_contact_radius при rho==rho_0: 0.5*(dxp^3)^(1/3)=0.5*dxp
    const double plate_top = 0.0; // одна ячейка толщиной dxp, центр в 0 -- по построению contribution идёт от x[2]=0

    // --- Нить: прямая, N бусин вдоль x, y=0, летит вниз (vz=-v0). ----------
    const double E_f = 2.0e7, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0; // "мягкий" материал -- см. test_e2e_fiber_contact.cpp
    const double r_f = 0.5 * d_f;
    const double l0 = 2.0e-3;
    const int N = 5;
    const double A = pi * d_f * d_f / 4.0;
    const double mass_each = rho_f * A * l0;

    const double gap0 = 2.0e-3;     // начальный зазор до контактного порога r_f+r_s
    const double v0 = 0.5;          // м/с, вниз
    const double z0 = plate_top + r_s + r_f + gap0;

    std::vector<std::array<double, 3>> beads;
    for (int k = -(N / 2); k <= N / 2; ++k) beads.push_back({static_cast<double>(k) * l0, 0.0, z0});

    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    materials.push_back(make_fiber_material(E_f, nu_f, d_f, rho_f));
    state.resize(beads.size());
    topology.resize(beads.size());
    append_chain(state, topology, 0, 0, 0, beads, mass_each);
    for (std::size_t i = 0; i < beads.size(); ++i) state.vx[2][i] = -v0;

    const std::size_t center = static_cast<std::size_t>(N / 2);
    check(state.vx[2][center] < 0.0, "самопроверка: нить изначально падает вниз (vz<0)");
    check(z0 > r_f + r_s, "самопроверка: начальный зазор положителен (нет проникновения на старте)");

    FiberOnPlateIntegrator integ;
    double time = 0.0;
    double max_penetration_observed = 0.0;
    double running_min_z = std::numeric_limits<double>::max();

    const double recovery_threshold = 2.0e-4; // 0.2 мм -- восстановление после минимума z
    bool bounce_confirmed = false;
    double vz_before = state.vx[2][center], vz_after = 0.0;

    const int max_steps = 6000;
    const double max_time = 2.0e-2;
    int steps_done = 0;
    for (; steps_done < max_steps && time < max_time; ++steps_done)
    {
        const double dt = integ.step(state, topology, materials, smats, plate, nplate);
        time += dt;
        max_penetration_observed = std::max(max_penetration_observed, integ.last_max_penetration);
        const double z = state.x[2][center];
        running_min_z = std::min(running_min_z, z);

        if (!bounce_confirmed && z > running_min_z + recovery_threshold)
        {
            bounce_confirmed = true;
            vz_after = state.vx[2][center];
            break;
        }
    }

    std::printf("    steps=%d time=%.4e s bounce_confirmed=%d\n", steps_done, time, bounce_confirmed ? 1 : 0);
    std::printf("    max_penetration_observed=%.4e m (r_f+r_s=%.4e m)\n", max_penetration_observed, r_f + r_s);
    std::printf("    running_min_z(center)=%.4e m (contact threshold z=%.4e m)\n", running_min_z, r_f + r_s);
    std::printf("    vz(center): до контакта %.4f -> после отскока %.4f м/с\n", vz_before, vz_after);

    check(max_penetration_observed > 0.0, "контакт с пластиной действительно произошёл (max_penetration_observed > 0)");
    check(max_penetration_observed < 0.5 * (r_f + r_s),
          "проникновение остаётся малым относительно r_f+r_s (< 50%) -- нить не проходит сквозь пластину");
    check(bounce_confirmed, "отскок (восстановление высоты) зафиксирован в пределах отведённого времени/шагов");
    check(vz_after > 0.0, "скорость центральной бусины сменила знак (падала -- отскочила от пластины)");
}

// =========================================================================
// (В) трение нить-Solid на пластине (notes/fiber-contact/fiber-solid-friction-implementation.md,
// §7.2). Та же сцена, что test_fiber_falls_on_plate (нить над пластиной падает
// вниз), плюс касательная скорость v_t0 вдоль x (вдоль нити). Интегратор
// повторяет порядок шага src/time_integration.cpp на РЕАЛЬНОМ решателе
// (ode::PredictorCorrectorSolver / EulerSolver):
//   predict -> RHS (Fiber-Bond + Fiber-Solid нормаль + слоты трения) -> dt ->
//   единая стадия трения (пробные скорости из solver.fiber_velocity_update /
//   solver.solid_velocity_update; c_k, c_s -- compute_friction_counts) ->
//   a^f в solver.*_rhs -> advance.
//   - кинематическая пластина (a, b, c): скорость Solid не интегрируется
//     (solver.reset(0, 0, n): решатель не трогает Solid), в стадию передаётся
//     SolidVelocityUpdate{base = скорость пластины (нули), h = 0, av = nullptr} и
//     solid_kinematic = true (c_s/m_s = 0, вкладов в Solid нет);
//   - свободная пластина (d): каждая частица пластины интегрируется тем же
//     решателем как свободная масса (solver.reset(nplate, nplate, n); Solid-SPH
//     нет, XSPH-поправка нулевая), solid_kinematic = false.
// =========================================================================
struct FsStepInfo
{
    bool has_slots = false;           // были активные слоты трения на шаге
    std::size_t active = 0, stick = 0, slip = 0, shared_bead = 0, shared_solid = 0;
    double ut_max = 0.0;              // max |u_t| точек контакта по итоговым скоростям шага (слоты ЭТОГО шага)
    double ut_x_min = 1.0e300;        // min x-компоненты u_t по точкам контакта (проекция на исходное направление скольжения +x)
    double capacity_fiber = 0.0;      // sum mu*N*h / M_fiber -- макс. изменение средней касательной скорости нити за шаг
    double capacity_contact = 0.0;    // max mu*N*h*W (с c_k, c_s) по контактам, м/с
    double imp_f[3] = {0, 0, 0};      // импульс трения на нить за шаг: sum m*h*a^f
    double imp_s[3] = {0, 0, 0};      // импульс трения на пластину за шаг
    double N_h = 0.0;                 // sum N*h по активным слотам (импульс нормальной силы, модуль)
    double h = 0.0;
};

template <class Solver>
struct FrictionPlateIntegrator
{
    Solver solver;
    rhs::CpuExecutor exec{1};
    bool free_plate = false;
    double safety = 0.3;
    double dt_prev = 0.0;
    bool inited = false;

    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> bond_buffers;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fs_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fs_solid_buffers;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fr_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fr_solid_buffers;
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> count_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidCountDelta>> count_solid_buffers;
    std::vector<double> bead_count, solid_count;
    std::vector<interactions::FiberSolidFrictionContact> slots;
    std::vector<interactions::FiberFrictionContact> no_ff_slots; // нить-нить контактов нет
    std::vector<std::vector<double>> zero_dvx;                    // ускорение Solid без трения у кинематической пластины

    double last_max_penetration = 0.0;
    FsStepInfo info;

    double step(fiber::FiberState& fb, const fiber::FiberTopology& topo, const std::vector<fiber::FiberMaterial>& fmats,
                const Materials& smats, ParticleSystem& plate, std::size_t nplate)
    {
        const std::size_t n = fb.size();
        const std::size_t nsolid = free_plate ? nplate : 0; // сколько частиц интегрирует решатель
        if (!inited)
        {
            solver.reset(nsolid, nsolid, n);
            zero_dvx.assign(dim, std::vector<double>(nplate, 0.0));
            inited = true;
        }

        solver.predict(plate, nsolid, fb, n, dt_prev);

        std::vector<fiber::FiberDelta> bond_delta = interactions::compute_fiber_bond(exec, fb, topo, fmats, bond_buffers);
        const double bond_dt_limit = interactions::compute_fiber_bond_dt_limit(exec, fb, topo, fmats);

        interactions::FiberSolidContactSearchResult fs_search;
        interactions::build_fiber_solid_contact_search(exec, smats, plate, nplate, fb, topo, fmats, fs_search);
        const bool has_active = interactions::compute_fiber_solid_contact_contributions(
            exec, smats, plate, fb, fmats, fs_search.primitives, fs_fiber_buffers, fs_solid_buffers, slots);
        std::vector<fiber::FiberDelta> fs_delta = rhs::reduce_contributions(exec, fs_fiber_buffers, n);
        std::vector<interactions::SolidDelta> fs_solid_delta = rhs::reduce_contributions(exec, fs_solid_buffers, nplate);
        const double fs_dt_limit = interactions::compute_fiber_solid_contact_dt_limit(exec, smats, plate, fb, fmats,
                                                                                        fs_search.primitives);
        double max_pen = 0.0;
        for (const auto& p : fs_search.primitives)
        {
            const interactions::FiberSolidContactGeometry g =
                interactions::compute_fiber_solid_contact_geometry(smats, plate, fb, fmats, p);
            if (g.contact) max_pen = std::max(max_pen, g.delta);
        }
        last_max_penetration = max_pen;

        const double dt = safety * std::min(bond_dt_limit, fs_dt_limit);

        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < n; ++i) solver.fiber_rhs.dv[d][i] = bond_delta[i].dv[d] + fs_delta[i].dv[d];
        if (free_plate)
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < nplate; ++i) solver.solid_rhs.dvx[d][i] = fs_solid_delta[i].dv[d];

        // --- единая стадия трения (после dt, до advance) ------------------------
        info = FsStepInfo{};
        info.h = 0.0;
        if (has_active)
        {
            const auto fupd = solver.fiber_velocity_update(fb, dt);
            ode::SolidVelocityUpdate supd{};
            if (free_plate)
                supd = solver.solid_velocity_update(plate, dt);
            else
            {
                for (std::size_t d = 0; d < dim; ++d) supd.base[d] = &plate.vx[d]; // скорость пластины (нули)
                supd.h = 0.0;
                supd.av = nullptr;
            }
            const std::vector<std::vector<double>>& dvx_ref = free_plate ? solver.solid_rhs.dvx : zero_dvx;
            const bool kin = !free_plate;

            interactions::compute_friction_counts(exec, n, nplate, no_ff_slots, slots, count_fiber_buffers,
                                                   count_solid_buffers, bead_count, solid_count);
            const interactions::FiberFrictionDiagnostics dg = interactions::compute_friction_diagnostics(
                exec, fb, fupd, solver.fiber_rhs.dv, no_ff_slots, plate, supd, dvx_ref, slots, bead_count, solid_count, kin);
            info.has_slots = true;
            info.active = dg.active_count;
            info.stick = dg.stick_count;
            info.slip = dg.slip_count;
            info.shared_bead = dg.shared_bead_count;
            info.shared_solid = dg.shared_solid_count;
            info.h = fupd.h;

            double M = 0.0;
            for (std::size_t i = 0; i < n; ++i) M += fb.mass[i];
            for (const auto& c : slots)
            {
                if (!c.active || !(c.mu * c.N > 0.0)) continue;
                interactions::FiberFrictionImpulse imp;
                interactions::fiber_solid_friction_impulse(fb, fupd, solver.fiber_rhs.dv, plate, supd, dvx_ref, c, imp,
                                                            bead_count.data(), solid_count.data(), kin);
                info.capacity_fiber += c.mu * c.N * fupd.h / M;
                info.capacity_contact = std::max(info.capacity_contact, imp.mu_N_h * imp.W);
                info.N_h += c.N * fupd.h;
            }

            interactions::compute_friction_contributions(exec, fb, fupd, solver.fiber_rhs.dv, no_ff_slots, plate, supd,
                                                          dvx_ref, slots, bead_count, solid_count, kin, fr_fiber_buffers,
                                                          fr_solid_buffers);
            const std::vector<fiber::FiberDelta> fr_f = rhs::reduce_contributions(exec, fr_fiber_buffers, n);
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i)
                {
                    solver.fiber_rhs.dv[d][i] += fr_f[i].dv[d];
                    info.imp_f[d] += fb.mass[i] * fupd.h * fr_f[i].dv[d];
                }
            if (free_plate)
            {
                const std::vector<interactions::SolidDelta> fr_s = rhs::reduce_contributions(exec, fr_solid_buffers, nplate);
                for (std::size_t d = 0; d < dim; ++d)
                    for (std::size_t i = 0; i < nplate; ++i)
                    {
                        solver.solid_rhs.dvx[d][i] += fr_s[i].dv[d];
                        info.imp_s[d] += plate.mass[i] * supd.h * fr_s[i].dv[d];
                    }
            }
        }

        solver.advance(plate, nsolid, fb, n, dt);

        // |u_t| точек контакта по итоговым скоростям шага (слоты и нормали этого шага)
        if (info.has_slots)
        {
            for (const auto& c : slots)
            {
                if (!c.active || !(c.mu * c.N > 0.0)) continue;
                double u[3];
                for (std::size_t d = 0; d < dim; ++d)
                {
                    double vF = fb.vx[d][c.a];
                    if (c.b != -1) vF = (1.0 - c.lambda) * vF + c.lambda * fb.vx[d][static_cast<std::size_t>(c.b)];
                    u[d] = vF - plate.vx[d][c.s];
                }
                double un = 0.0;
                for (std::size_t d = 0; d < dim; ++d) un += u[d] * c.n[d];
                double ut2 = 0.0;
                for (std::size_t d = 0; d < dim; ++d) ut2 += (u[d] - un * c.n[d]) * (u[d] - un * c.n[d]);
                info.ut_max = std::max(info.ut_max, std::sqrt(ut2));
                info.ut_x_min = std::min(info.ut_x_min, u[0] - un * c.n[0]);
            }
        }
        dt_prev = dt;
        return dt;
    }
};

// Пластина (как в test_fiber_falls_on_plate); density -- плотность частиц (для
// свободной пластины «малая масса»: меньшая плотность при том же объёме, чтобы не
// менять радиус контакта r_s = 0.5*(m/rho)^(1/3) = 0.5*dxp).
struct FsScene
{
    ParticleSystem plate;
    std::size_t nplate = 0;
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    std::size_t center = 0;
    double r_f = 0.0, r_s = 0.0, mass_bead = 0.0;
};

void build_fs_scene_plate(FsScene& sc, double mu_fs, double v_t0, double plate_density, double y_off)
{
    const double dxp = 1.5e-3;
    build_plate(sc.plate, 5, dxp, plate_density);
    sc.nplate = sc.plate.size();
    sc.r_s = 0.5 * dxp;

    const double E_f = 2.0e7, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    sc.r_f = 0.5 * d_f;
    const double l0 = 2.0e-3;
    const int N = 5;
    const double A = pi * d_f * d_f / 4.0;
    sc.mass_bead = rho_f * A * l0;
    // Зазор меньше, чем в test_fiber_falls_on_plate (0.5 мм против 2 мм), а нить стартует
    // со сдвигом по x на -v_t0*t_подхода: за время падения до контакта (t_подхода ~ 1.5 мс,
    // первый контакт при z ~ 1.0 мм) касательная скорость сместила бы нить на v_t0*t_подхода
    // -- при v_t0 = 2 м/с на 3 мм -- и она приходила бы к контакту не над серединой пластины
    // (полуширина пластины 7.5 мм). Сдвиг -- геометрия сцены, а не подбор результата: к
    // моменту контакта нить оказывается над центром пластины при любом v_t0.
    const double gap0 = 0.5e-3;
    const double v0 = 0.5;
    const double z0 = sc.r_s + sc.r_f + gap0;
    const double z_first_contact = 1.0e-3; // оценка высоты центра нити при первом касании (seam-геометрия)
    const double t_approach = (z0 - z_first_contact) / v0;
    const double x_shift = -v_t0 * t_approach;

    std::vector<std::array<double, 3>> beads;
    for (int k = -(N / 2); k <= N / 2; ++k) beads.push_back({static_cast<double>(k) * l0 + x_shift, y_off, z0});
    sc.materials.push_back(make_fiber_material(E_f, nu_f, d_f, rho_f));
    sc.materials[0].mu_fs = mu_fs;
    sc.state.resize(beads.size());
    sc.topology.resize(beads.size());
    append_chain(sc.state, sc.topology, 0, 0, 0, beads, sc.mass_bead);
    for (std::size_t i = 0; i < beads.size(); ++i)
    {
        sc.state.vx[2][i] = -v0;
        sc.state.vx[0][i] = v_t0;
    }
    sc.center = static_cast<std::size_t>(N / 2);
}

struct FsResult
{
    bool bounce = false;   // отскок подтверждён (кинематическая пластина)
    bool contact_seen = false;
    int steps = 0;
    double time = 0.0;
    double max_pen = 0.0;
    double vx0 = 0.0, vx1 = 0.0, vz0 = 0.0, vz1 = 0.0; // средняя скорость нити (по массе): до / в конце
    double fiber_p0[3] = {}, fiber_p1[3] = {}, plate_p0[3] = {}, plate_p1[3] = {};
    double imp_f[3] = {}, imp_s[3] = {}; // импульс трения на нить / пластину за прогон
    double int_N = 0.0;                  // sum N*h по шагам (модуль импульса нормальной силы)
    int contact_steps = 0, stick_steps = 0, slip_steps = 0, first_stick_step = -1;
    int stick_only_steps = 0;            // шаги stick>0 && slip==0
    double worst_ut_stick_any = 0.0;     // max |u_t| по шагам с stick>0, начиная с первого прилипания
    double worst_ut_stick_only = 0.0;    // то же по шагам stick>0 && slip==0
    int steps_to_converge = -1;          // шагов после первого прилипания до |u_t| <= 1e-9 v_t0 (первый такой шаг)
    int conv_step = -1;                  // номер шага сходимости
    double t_conv = 0.0;                 // время от первого прилипания до сходимости, с
    double worst_ut_after_conv = 0.0;    // max |u_t| на шагах stick>0 && slip==0 ПОСЛЕ сходимости
    int after_conv_stick_steps = 0;      // число таких шагов
    double min_utx_stick = 1.0e300;      // min проекции u_t на +x (исходное направление) на шагах stick>0 && slip==0
    double max_growth_stick = 0.0;       // max (|u_t|_k - |u_t|_{k-1}) на подряд идущих шагах stick>0 && slip==0
    double ut_max_contact = 0.0;
    double cap_fiber_peak = 0.0, cap_fiber_total = 0.0, cap_contact_peak = 0.0;
    std::size_t active_max = 0, shared_bead_max = 0, shared_solid_max = 0;
    double m_total = 0.0;
    double dt_min = 1e9, dt_max = 0.0;
    std::vector<double> ut_series;
    std::vector<int> ut_steps;
};

static void fs_momentum(const fiber::FiberState& fb, const ParticleSystem& plate, std::size_t np, double pf[3], double pp[3])
{
    for (int d = 0; d < 3; ++d)
    {
        pf[d] = 0.0;
        pp[d] = 0.0;
        for (std::size_t i = 0; i < fb.size(); ++i) pf[d] += fb.mass[i] * fb.vx[d][i];
        for (std::size_t i = 0; i < np; ++i) pp[d] += plate.mass[i] * plate.vx[d][i];
    }
}

static double fs_mean_v(const fiber::FiberState& fb, std::size_t d)
{
    double sp = 0.0, sm = 0.0;
    for (std::size_t i = 0; i < fb.size(); ++i)
    {
        sp += fb.mass[i] * fb.vx[d][i];
        sm += fb.mass[i];
    }
    return sp / sm;
}

// Прогон. free_plate=false: до подтверждённого отскока (как test_fiber_falls_on_plate);
// free_plate=true: до конца контакта (нет проникания после того, как контакт был) либо
// лимита шагов/времени.
template <class Solver>
FsResult run_fs_case(double mu_fs, double v_t0, double safety, bool free_plate, double plate_density = 2700.0,
                     double y_off = 0.0)
{
    FsResult res;
    FsScene sc;
    build_fs_scene_plate(sc, mu_fs, v_t0, plate_density, y_off);
    fiber::FiberState& fb = sc.state;
    res.vx0 = fs_mean_v(fb, 0);
    res.vz0 = fs_mean_v(fb, 2);
    fs_momentum(fb, sc.plate, sc.nplate, res.fiber_p0, res.plate_p0);
    for (std::size_t i = 0; i < fb.size(); ++i) res.m_total += fb.mass[i];
    if (free_plate)
        for (std::size_t i = 0; i < sc.nplate; ++i) res.m_total += sc.plate.mass[i];

    Materials smats;
    FrictionPlateIntegrator<Solver> integ;
    integ.free_plate = free_plate;
    integ.safety = safety;

    double time = 0.0;
    double running_min_z = std::numeric_limits<double>::max();
    const double recovery_threshold = 2.0e-4;
    const int max_steps = free_plate ? 60000 : 80000;
    const double max_time = 2.0e-2;
    bool prev_contact = false;
    double t_first_stick = 0.0, ut_prev_stick = -1.0;
    int k = 0;
    for (; k < max_steps && time < max_time; ++k)
    {
        const double dt = integ.step(fb, sc.topology, sc.materials, smats, sc.plate, sc.nplate);
        time += dt;
        res.dt_min = std::min(res.dt_min, dt);
        res.dt_max = std::max(res.dt_max, dt);
        res.max_pen = std::max(res.max_pen, integ.last_max_penetration);
        const auto& info = integ.info;
        for (int d = 0; d < 3; ++d)
        {
            res.imp_f[d] += info.imp_f[d];
            res.imp_s[d] += info.imp_s[d];
        }
        res.int_N += info.N_h;

        if (info.has_slots)
        {
            ++res.contact_steps;
            res.ut_series.push_back(info.ut_max);
            res.ut_steps.push_back(k);
            res.ut_max_contact = std::max(res.ut_max_contact, info.ut_max);
            res.cap_fiber_peak = std::max(res.cap_fiber_peak, info.capacity_fiber);
            res.cap_fiber_total += info.capacity_fiber;
            res.cap_contact_peak = std::max(res.cap_contact_peak, info.capacity_contact);
            res.active_max = std::max(res.active_max, info.active);
            res.shared_bead_max = std::max(res.shared_bead_max, info.shared_bead);
            res.shared_solid_max = std::max(res.shared_solid_max, info.shared_solid);
            if (info.stick > 0)
            {
                ++res.stick_steps;
                if (res.first_stick_step < 0)
                {
                    res.first_stick_step = k;
                    t_first_stick = time;
                }
                res.worst_ut_stick_any = std::max(res.worst_ut_stick_any, info.ut_max);
                if (res.steps_to_converge < 0 && info.ut_max <= 1.0e-9 * v_t0)
                {
                    res.steps_to_converge = k - res.first_stick_step;
                    res.conv_step = k;
                    res.t_conv = time - t_first_stick;
                }
                if (info.slip == 0)
                {
                    ++res.stick_only_steps;
                    res.worst_ut_stick_only = std::max(res.worst_ut_stick_only, info.ut_max);
                    if (res.conv_step >= 0 && k > res.conv_step)
                    {
                        ++res.after_conv_stick_steps;
                        res.worst_ut_after_conv = std::max(res.worst_ut_after_conv, info.ut_max);
                    }
                    res.min_utx_stick = std::min(res.min_utx_stick, info.ut_x_min);
                    if (ut_prev_stick >= 0.0) res.max_growth_stick = std::max(res.max_growth_stick, info.ut_max - ut_prev_stick);
                    ut_prev_stick = info.ut_max;
                }
                else
                    ut_prev_stick = -1.0;
            }
            else
                ut_prev_stick = -1.0;
            if (info.slip > 0) ++res.slip_steps;
        }

        const bool in_contact = integ.last_max_penetration > 0.0;
        if (in_contact) res.contact_seen = true;
        const double z = fb.x[2][sc.center];
        running_min_z = std::min(running_min_z, z);
        if (!free_plate)
        {
            if (!res.bounce && z > running_min_z + recovery_threshold)
            {
                res.bounce = true;
                break;
            }
        }
        else
        {
            if (res.contact_seen && !in_contact && prev_contact) break; // контакт закончился
        }
        prev_contact = in_contact;
    }
    res.steps = k;
    res.time = time;
    res.vx1 = fs_mean_v(fb, 0);
    res.vz1 = fs_mean_v(fb, 2);
    fs_momentum(fb, sc.plate, sc.nplate, res.fiber_p1, res.plate_p1);
    return res;
}

// Печать сводки прогона.
static void report_fs(const char* label, const FsResult& r, double v_t0)
{
    std::printf("  %s: steps=%d time=%.4e с dt=[%.3e,%.3e] bounce=%d max_pen=%.4e м контакт-шагов=%d "
                "(активных слотов max=%zu)\n",
                label, r.steps, r.time, r.dt_min, r.dt_max, r.bounce ? 1 : 0, r.max_pen, r.contact_steps, r.active_max);
    std::printf("      vx нити: %.6f -> %.6f м/с; vz нити: %.6f -> %.6f м/с\n", r.vx0, r.vx1, r.vz0, r.vz1);
    std::printf("      импульс трения на нить (x): %.6e; на пластину (x): %.6e; int N dt = %.6e кг*м/с\n", r.imp_f[0],
                r.imp_s[0], r.int_N);
    std::printf("      stick-шагов=%d (из них без slip=%d), slip-шагов=%d, первый stick=%d, shared_bead max=%zu, "
                "shared_solid max=%zu\n",
                r.stick_steps, r.stick_only_steps, r.slip_steps, r.first_stick_step, r.shared_bead_max,
                r.shared_solid_max);
    std::printf("      ёмкость прилипания: нить (sum mu*N*h/M) пик за шаг=%.4e, суммарно за контакт=%.4e м/с (v_t0=%.3f); "
                "контакт mu*N*h*W пик=%.4e м/с\n",
                r.cap_fiber_peak, r.cap_fiber_total, v_t0, r.cap_contact_peak);
    std::printf("      max|u_t| за контакт=%.4e, max|u_t| на шагах со stick>0=%.4e (без slip: %.4e), шагов до |u_t|<=1e-9*v_t0 "
                "после первого stick=%d (%.3e с)\n",
                r.ut_max_contact, r.worst_ut_stick_any, r.worst_ut_stick_only, r.steps_to_converge, r.t_conv);
    std::printf("      после сходимости: шагов stick без slip=%d, max|u_t| на них=%.3e; max рост |u_t| за шаг на подряд идущих "
                "stick-шагах=%.3e\n",
                r.after_conv_stick_steps, r.worst_ut_after_conv, r.max_growth_stick);
}

void test_fiber_plate_friction()
{
    std::printf("test_fiber_plate_friction:\n");
    const double v_t0_slide = 2.0;
    const double mu_b = 0.2;

    // --- (a) mu_fs = 0: касательная скорость нити контактом не меняется ------
    FsResult a = run_fs_case<ode::PredictorCorrectorSolver>(0.0, v_t0_slide, 0.3, false);
    report_fs("(a) mu_fs=0, v_t0=2.0", a, v_t0_slide);
    check(a.bounce, "(а) отскок зафиксирован");
    check(a.contact_steps == 0, "(а) mu_fs=0: активных слотов трения нет ни на одном шаге (стадия не выполняется)");
    const double dvx_a = std::fabs(a.vx1 - a.vx0), dvz_a = std::fabs(a.vz1 - a.vz0);
    std::printf("      |d vx| = %.4e, |d vz| = %.4e (отношение %.4f)\n", dvx_a, dvz_a, dvx_a / dvz_a);
    // Допуск на наклон нормали (как в test_e2e_cross_friction (a)): без трения x-импульс
    // передаётся только x-компонентой нормальной силы (n_x != 0 у крайних контактов).
    check(dvx_a < 0.5 * dvz_a, "(а) mu_fs=0: |d vx| < 0.5*|d vz| (касательный импульс передаёт только наклон нормали)");

    // --- (b) быстрое скольжение -------------------------------------------
    FsResult b = run_fs_case<ode::PredictorCorrectorSolver>(mu_b, v_t0_slide, 0.3, false);
    report_fs("(b) mu_fs=0.2, v_t0=2.0", b, v_t0_slide);
    check(b.bounce, "(b) отскок зафиксирован");
    const double dvx_b = b.vx0 - b.vx1, dvn_b = std::fabs(b.vz1 - b.vz0);
    const double est = mu_b * dvn_b;
    std::printf("      торможение d vx = %.4e м/с; оценка mu*|d v_n| = %.4e м/с (отношение %.3f)\n", dvx_b, est, dvx_b / est);
    check(dvx_b > 0.0, "(b) касательная скорость нити уменьшилась (vx1 < vx0)");
    check(b.vx1 > 0.0, "(b) скольжение не прекратилось (vx1 > 0): v_t0 много больше ёмкости прилипания");
    check(b.cap_fiber_total < 0.5 * v_t0_slide, "(b) v_t0 больше суммарной ёмкости прилипания за контакт в 2 и более раза");
    check(dvx_b > 0.5 * est && dvx_b < 2.0 * est, "(b) торможение по порядку величины равно mu*|d v_n| (0.5..2)");
    const double fr = std::fabs(b.imp_f[0]);
    check(fr <= mu_b * b.int_N * (1.0 + 1e-12), "(b) импульс трения не превышает mu * (импульс нормальной силы)");
    check(fr >= 0.8 * mu_b * b.int_N, "(b) при быстром скольжении трение почти всё время на пределе: >= 0.8 mu*int N dt");
    check(std::fabs(b.imp_f[1]) < 0.5 * fr, "(b) импульс трения направлен вдоль скольжения (y-компонента мала)");
}

// --- (c) медленное скольжение --------------------------------------------
// v_t0 = 0.05 м/с много меньше ёмкости прилипания за контакт (печатается:
// sum mu*N*h/M по шагам контакта, ~0.2 м/с). mu НЕ поднимался: при mu_fs = 0.2
// (значение пункта (b)) ёмкость уже достаточна, а причина расхождения с
// критерием спецификации не в ёмкости (см. ниже).
//
// Две геометрии одной и той же сцены (нить над плитой dxp = 1.5 мм, первый контакт
// при z ~ 1 мм):
//   y_off = 0      -- как в test_fiber_falls_on_plate: нить над ШВОМ между двумя рядами
//                     частиц пластины (y = +-0.75 мм); контакты идут с частицами обоих
//                     рядов, нормали наклонены на ~40 градусов от вертикали. Скорость
//                     падения (0.5 м/с) даёт на наклонных нормалях касательную
//                     проекцию u_t ~ 0.3 м/с -- на порядок больше v_t0, поэтому ||u_t||
//                     точек контакта в этой геометрии определяется падением, а не v_t0,
//                     и критерии на ||u_t||/v_t0 из §7.2 (c) неприменимы. Здесь
//                     проверяются только отскок, ёмкость и сходимость импульса по dt.
//   y_off = 0.75 мм -- нить над РЯДОМ частиц: нормаль вертикальна, ||u_t|| == касательная
//                     скорость вдоль x (и y); критерии §7.2 (c) применимы буквально.
void check_slow_case(const char* label, const FsResult& r, double v_t0, bool vertical_normals)
{
    const std::string L(label);
    check(r.bounce, L + ": отскок зафиксирован");
    check(r.contact_steps > 0, L + ": контакт с трением был (активные слоты)");
    check(r.cap_fiber_total > v_t0, L + ": суммарная ёмкость прилипания за контакт больше v_t0 (" +
                                        std::to_string(r.cap_fiber_total) + " м/с)");
    if (!vertical_normals) return;
    check(r.stick_steps > 0, L + ": прилипание наступило (stick_count > 0 хотя бы на одном шаге)");
    check(r.ut_max_contact <= v_t0, L + ": max||u_t|| за контакт <= v_t0 (max = " + std::to_string(r.ut_max_contact) + ")");
    // Критерий (i) (вместо «||u_t|| <= 1e-9*v_t0 сразу после первого шага прилипания», решение
    // пользователя): после времени релаксации T_relax (от первого шага с stick_count > 0 до
    // первого шага с ||u_t|| <= 1e-9*v_t0) и до конца фазы прилипания ||u_t|| <= 1e-9*v_t0.
    // Причина релаксации: mass splitting (§5.2, §6.4 спецификации) -- бусина входит в c = 5-6
    // контактов, каждый берёт долю 1/c массы, поэтому точное прилипание достигается не за один
    // шаг, а за несколько (недокомпенсирование); итерации Якоби/Гаусса-Зейделя запрещены §9.
    check(r.steps_to_converge >= 0, L + ": ||u_t|| достигает 1e-9*v_t0 за конечное время T_relax = " +
                                        std::to_string(r.t_conv) + " с (" + std::to_string(r.steps_to_converge) + " шагов)");
    check(r.worst_ut_after_conv <= 1.0e-9 * v_t0,
          L + ": (i) после T_relax ||u_t|| <= 1e-9*v_t0 на шагах stick без slip (шагов: " +
              std::to_string(r.after_conv_stick_steps) + ", max = " + std::to_string(r.worst_ut_after_conv) + ")");
    // (iii) в фазе прилипания (подряд идущие шаги stick>0 && slip==0): рост ||u_t|| за шаг <= 0.01*v_t0
    // и перескок проекции u_t на исходное направление (+x) не глубже -0.05*v_t0.
    // Допуск ОТ v_t0, а не от округления, потому что один проход Якоби (с mass splitting) на
    // связанных контактах на общих бусинах при неравномерной u_t вдоль нити даёт малые остатки:
    // контакт с малой u_t получает поправки соседей (на шаге ||u_t|| этого контакта меняет
    // знак на малую величину), плюс меняется нагрузка между шагами (изгибная волна, N падает).
    // Итерации запрещены §9, и даже 30 развёрток Гаусса-Зейделя не дают 1e-9*v_t0 (диагностика:
    // рост 5.5e-5, перескок -7.4e-5 м/с при dt). Измеренные остатки: рост ~0.2 % v_t0, перескок
    // ~2.3 % v_t0. Проверка при этом ловит возврат дребезга: без mass splitting (c = 1) перескок
    // проекции 160-220 % (см. test_fs_mass_splitting), то есть на два порядка больше допуска.
    const double grow_tol = 0.01 * v_t0, flip_tol = 0.05 * v_t0;
    std::printf("      (iii) %s: рост ||u_t|| за шаг max = %.3e (допуск %.3e), перескок проекции min = %.3e (допуск %.3e)\n",
                label, r.max_growth_stick, grow_tol, r.min_utx_stick, -flip_tol);
    check(r.max_growth_stick <= grow_tol, L + ": (iii) рост ||u_t|| за шаг в фазе прилипания <= 0.01*v_t0 (max = " +
                                              std::to_string(r.max_growth_stick) + ")");
    check(r.min_utx_stick >= -flip_tol, L + ": (iii) перескок проекции u_t на исходное направление >= -0.05*v_t0 (min = " +
                                              std::to_string(r.min_utx_stick) + ")");
}

void test_fiber_plate_friction_slow(double y_off, bool vertical_normals, const char* tag)
{
    std::printf("test_fiber_plate_friction_slow (%s):\n", tag);
    const double v_t0 = 0.05;
    const double mu = 0.2;

    FsResult c1 = run_fs_case<ode::PredictorCorrectorSolver>(mu, v_t0, 0.3, false, 2700.0, y_off);
    FsResult c5 = run_fs_case<ode::PredictorCorrectorSolver>(mu, v_t0, 0.06, false, 2700.0, y_off);
    report_fs("(c) dt (safety 0.3)", c1, v_t0);
    report_fs("(c) dt/5 (safety 0.06)", c5, v_t0);

    check_slow_case("(c) dt", c1, v_t0, vertical_normals);
    check_slow_case("(c) dt/5", c5, v_t0, vertical_normals);
    const double a1 = std::fabs(c1.imp_f[0]), a5 = std::fabs(c5.imp_f[0]);
    const double rel = std::fabs(a1 - a5) / std::max(a1, a5);
    std::printf("      сходимость по шагу: импульс трения (x) на нить: dt %.6e, dt/5 %.6e, расхождение %.2f %%\n", a1, a5,
                100.0 * rel);
    check(rel <= 0.10, "(c) импульс, переданный трением, при dt и dt/5 отличается не более чем на 10 %");

    if (!vertical_normals) return;
    // (ii) T_relax при dt и dt/5 совпадает в пределах 20 %.
    if (c1.steps_to_converge >= 0 && c5.steps_to_converge >= 0)
    {
        const double ratio = c1.t_conv / c5.t_conv;
        std::printf("      T_relax (первый stick -> первый шаг с ||u_t|| <= 1e-9*v_t0): dt %.4e с, dt/5 %.4e с (отношение %.3f)\n",
                    c1.t_conv, c5.t_conv, ratio);
        check(std::fabs(ratio - 1.0) <= 0.20, "(c) (ii) T_relax при dt и dt/5 совпадает в пределах 20 %");
    }
    else
        check(false, "(c) (ii) T_relax определён при dt и dt/5");
}

// --- (d) свободная пластина ------------------------------------------------
void test_fiber_plate_friction_free()
{
    std::printf("test_fiber_plate_friction_free:\n");
    const double v_t0 = 2.0;
    const double mu = 0.2;
    // Плотность частиц пластины в 10 раз меньше (объём тот же, радиус контакта тот же):
    // масса частицы 9.1e-7 кг, чуть меньше массы бусины нити (1.9e-6 кг).
    FsResult d = run_fs_case<ode::PredictorCorrectorSolver>(mu, v_t0, 0.3, true, 270.0);
    report_fs("(d) свободная пластина, mu_fs=0.2, v_t0=2.0", d, v_t0);
    check(d.contact_seen && d.contact_steps > 0, "(d) контакт с трением был");
    check(d.slip_steps > 0, "(d) трение работало (есть скользящие контакты)");

    const double scale = d.m_total * std::max(0.5, v_t0);
    double dp[3];
    for (int k = 0; k < 3; ++k) dp[k] = (d.fiber_p1[k] + d.plate_p1[k]) - (d.fiber_p0[k] + d.plate_p0[k]);
    std::printf("      суммарный импульс нити+пластины: |dp| = (%.3e, %.3e, %.3e), допуск %.3e\n", std::fabs(dp[0]),
                std::fabs(dp[1]), std::fabs(dp[2]), 1.0e-12 * scale);
    check(std::fabs(dp[0]) < 1.0e-12 * scale && std::fabs(dp[1]) < 1.0e-12 * scale && std::fabs(dp[2]) < 1.0e-12 * scale,
          "(d) суммарный линейный импульс нити и пластины сохраняется (все три компоненты)");

    // Импульс трения: на пластину равен и противоположен импульсу на нить.
    std::printf("      импульс трения: на нить (%.6e, %.3e, %.3e), на пластину (%.6e, %.3e, %.3e)\n", d.imp_f[0], d.imp_f[1],
                d.imp_f[2], d.imp_s[0], d.imp_s[1], d.imp_s[2]);
    double sum_fr[3];
    for (int k = 0; k < 3; ++k) sum_fr[k] = d.imp_f[k] + d.imp_s[k];
    const double fr = std::fabs(d.imp_f[0]);
    check(fr > 0.0, "(d) импульс трения на нить ненулевой");
    check(std::fabs(sum_fr[0]) <= 1e-12 * fr && std::fabs(sum_fr[1]) <= 1e-12 * fr && std::fabs(sum_fr[2]) <= 1e-12 * fr,
          "(d) импульс трения на пластину равен по модулю и противоположен по знаку импульсу на нить (1e-12)");
    check(d.imp_s[0] > 0.0, "(d) пластина увлекается трением вдоль +x (движения нити)");
    // Касательный (x) импульс пластины в целом: получен трением (доминирует) плюс x-компонента нормальной силы.
    const double dpl_x = d.plate_p1[0] - d.plate_p0[0];
    const double dfb_x = d.fiber_p1[0] - d.fiber_p0[0];
    std::printf("      x-импульс: пластина %.6e, нить %.6e (сумма %.3e)\n", dpl_x, dfb_x, dpl_x + dfb_x);
    check(dpl_x > 0.0 && std::fabs(dpl_x + dfb_x) < 1e-12 * scale,
          "(d) пластина получила касательный x-импульс, равный и противоположный изменению x-импульса нити");
    check(std::fabs(dpl_x - d.imp_s[0]) <= 0.2 * dpl_x, "(d) x-импульс пластины в основном (>= 80 %) создан трением");
}

} // namespace

int main()
{
    test_smoke_solid_only_3d();
    test_fiber_falls_on_plate();
    test_fiber_plate_friction();
    test_fiber_plate_friction_slow(0.0, false, "нить над швом между рядами, как в test_fiber_falls_on_plate");
    test_fiber_plate_friction_slow(0.75e-3, true, "нить над рядом частиц, нормаль вертикальна");
    test_fiber_plate_friction_free();

    std::printf("\n=== %s (провалов: %d) ===\n", e2e::failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
