// tests/test_e2e_fiber_solid.cpp
// E2E-критерии готовности этапа Ф6 (fiber/PLAN.md, «Fiber-Solid-Contact»,
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
#include "interactions/fiber_solid_contact.hpp"
#include "materials.hpp"
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
            if (g.contact) max_pen = std::max(max_pen, static_cast<double>(g.delta));
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

} // namespace

int main()
{
    test_smoke_solid_only_3d();
    test_fiber_falls_on_plate();

    std::printf("\n=== %s (провалов: %d) ===\n", e2e::failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
