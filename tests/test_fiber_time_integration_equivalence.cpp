// test_fiber_time_integration_equivalence.cpp
// Проверка того, что скопированный Solid-цикл внутри
// time_integration(cfg, fiber::Model&, ...) (src/time_integration.cpp, этап
// Ф3, fiber/PLAN.md, раздел «Ф3. Fiber-Bond») физически ЭКВИВАЛЕНТЕН старому
// циклу time_integration(cfg, ParticleSystem&, ...) при пустых нитях.
//
// ЗАЧЕМ. Новая перегрузка не переиспользует старый цикл -- она содержит его
// построчную копию (u/rho/s -- буквально; v/x -- через общие шаблоны
// pc_predict_velocity_half_step/pc_finalize_velocity_position*, см. комментарий
// TODO(Ф8) в src/time_integration.cpp). Ни один существующий тест её не
// вызывает с НЕПУСТЫМ Solid: test_e2e_bc2/test_e2e_bc6 используют старую
// сигнатуру напрямую, а test_fiber_bond.cpp упражняет только Fiber-физику при
// ntotal == 0 (Solid-путь при этом не проходит НИ ОДНОЙ итерации). Скопированный
// Solid-цикл внутри новой перегрузки, таким образом, до этого теста не был
// проверен вообще ничем -- риск скрытой регрессии, которая проявится только
// на Ф4+, когда что-то реально начнёт вызывать новую перегрузку с непустым
// Solid.
//
// КАК. Берём одну и ту же начальную конфигурацию Solid -- уменьшенную 3D-
// геометрию BC2 (tests/bc2_mini3d_geometry.hpp, тот же кейс, что и
// test_e2e_3d.cpp), клонируем её в две независимые копии и прогоняем
// одинаковое число шагов (RunGuard::set_step_limit) через обе перегрузки:
//   А) time_integration(cfg, sysA, ...)               -- старая сигнатура;
//   Б) time_integration(cfg, modelB, ...), modelB.fiber -- ПУСТОЙ (resize(0)),
//      modelB.topology -- ПУСТАЯ (resize(0), bonds.empty()).
// При пустых нитях новая перегрузка не должна вносить в Solid никакой
// дополнительной физики (Fiber-Bond kernel по пустой topology.bonds не
// производит ни одного contribution -- fiber_dt_limit == +inf, на итоговый
// dt не влияет), поэтому modelB.solid после N шагов должен совпасть с sysA
// после тех же N шагов с точностью до порядка округления (не обязательно
// побитно: порядок суммирования в OpenMP-редукциях парных сил может отличаться
// между двумя процессами/прогонами).
//
// Прогон повторяется для нескольких contact_force_type, чтобы проверить
// также ветку dt = min(CFL, wb.min_dt_contact) -- она продублирована в новой
// перегрузке отдельно от Solid-SPH части.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "input.hpp"
#include "int_force.hpp"
#include "time_integration.hpp"
#include "run_guard.hpp"
#include "fiber/model.hpp"

#include "e2e_harness.hpp"
#include "bc2_mini3d_geometry.hpp"

static_assert(dim == 3, "test_fiber_time_integration_equivalence требует сборки с -DSPH_DIM=3");

#ifndef SPH_TEST_DATA_DIR
#define SPH_TEST_DATA_DIR "tests/data"
#endif

namespace fs = std::filesystem;
using e2e::check;
using e2e::note;

namespace
{

constexpr size_t kSteps = 6; // достаточно, чтобы пройти и первый шаг (особая ветка), и n>=2 (предиктор-корректор)

fs::path config_path() { return fs::path(SPH_TEST_DATA_DIR) / "bc2_initial_conditions.txt"; }

// Готовит ParticleSystem и SimulationConfig ровно как test_e2e_3d.cpp
// (tests/test_e2e_3d.cpp, run_case): один и тот же детерминированный
// генератор геометрии, поэтому два независимых вызова с тем же delta_x дают
// побитно одинаковое начальное состояние.
void build_case(SimulationConfig& cfg, ParticleSystem& sys, size_t& n, int contact_force_type)
{
    cfg.apply_defaults();
    if (!fs::exists(config_path()))
        throw std::runtime_error("не найден конфиг теста: " + config_path().string());
    parse_sph_config(cfg, config_path().string());
    cfg.contact_force_type = contact_force_type;

    bc2m3::Geometry geom = bc2m3::build_geometry(sys, bc2m3::kDxSmoke);
    n = geom.ntotal;
    cfg.moni_particle = geom.moni;
    cfg.maxn = static_cast<size_t>(static_cast<double>(n) * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    const size_t prev = sys.size();
    sys.resize(cfg.maxn);
    for (size_t i = prev; i < sys.size(); ++i) sys.inzializ_index(i);
    for (size_t i = 0; i < n; ++i) urs(cfg.mats, sys, i);
}

// Максимум |a|, |b|, пол scale -- чтобы сравнение не требовало абсурдной
// точности у полей, которые в этой постановке остаются нулевыми (например,
// компоненты скорости, перпендикулярные направлению удара).
bool close(double a, double b, double floor_scale, double tol)
{
    const double scale = std::max({std::fabs(a), std::fabs(b), floor_scale});
    return std::fabs(a - b) <= tol * scale;
}

struct FieldDiff
{
    double worst = 0.0;
    std::string field;
    size_t particle = 0;
};

// Сравнивает все поля Solid, изменяемые time_integration, на первых n
// частицах. tol -- относительный допуск (см. шапку файла: coordinator
// просил диапазон 1e-9..1e-12; берём 1e-9 -- прогон использует OpenMP-
// редукции парных сил, порядок суммирования которых между двумя процессами
// не гарантирован, поэтому побитное совпадение не предполагается).
FieldDiff compare_solid(const ParticleSystem& a, const ParticleSystem& b, size_t n, double tol)
{
    FieldDiff worst;
    auto probe = [&](const char* name, double va, double vb, double floor_scale) {
        if (!close(va, vb, floor_scale, tol))
        {
            const double rel = std::fabs(va - vb) / std::max({std::fabs(va), std::fabs(vb), floor_scale});
            if (rel > worst.worst)
            {
                worst.worst = rel;
                worst.field = name;
            }
        }
    };
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            probe(d == 0 ? "x0" : (d == 1 ? "x1" : "x2"), a.x[d][i], b.x[d][i], 1e-3);
            probe(d == 0 ? "v0" : (d == 1 ? "v1" : "v2"), a.vx[d][i], b.vx[d][i], 1.0);
        }
        probe("u", a.u[i], b.u[i], 1.0);
        probe("rho", a.rho[i], b.rho[i], 1.0);
        probe("sxx", a.sxx[i], b.sxx[i], 1e6);
        probe("syy", a.syy[i], b.syy[i], 1e6);
        probe("szz", a.szz[i], b.szz[i], 1e6);
        probe("sxy", a.sxy[i], b.sxy[i], 1e6);
        probe("syz", a.syz[i], b.syz[i], 1e6);
        probe("sxz", a.sxz[i], b.sxz[i], 1e6);
        if (worst.worst > 0.0 && worst.particle == 0) worst.particle = i;
    }
    return worst;
}

void run_case(int contact_force_type)
{
    std::printf("\n[contact_force_type = %d]\n", contact_force_type);

    SimulationConfig cfg;
    ParticleSystem sysA;
    size_t n = 0;
    build_case(cfg, sysA, n, contact_force_type);

    double dtA = compute_cfl_timestep(cfg, sysA, n);
    ParticleSystem sysB = sysA; // побитная копия начального состояния
    double dtB = dtA;

    fiber::Model modelB;
    modelB.solid = sysB;
    modelB.fiber.resize(0);
    modelB.topology.resize(0);
    // modelB.materials остаётся пустым -- Fiber-Bond kernel не обращается к
    // нему при topology.bonds.empty() (compute_fiber_bond_contributions/
    // compute_fiber_bond_dt_limit возвращают рано при nbonds == 0, см.
    // include/interactions/fiber_bond.hpp).

    double timeA = 0.0, timeB = 0.0;
    RunGuard guardA(0.0, 0.0, 0);
    guardA.set_step_limit(kSteps);
    RunGuard guardB(0.0, 0.0, 0);
    guardB.set_step_limit(kSteps);

    {
        e2e::ScopedWorkDir wd("sph_fiber_ti_equiv", "old_c" + std::to_string(contact_force_type));
        e2e::ScopedSilence silence(true);
        time_integration(cfg, sysA, n, dtA, timeA, bc2m3::kDxSmoke, &guardA);
    }
    {
        e2e::ScopedWorkDir wd("sph_fiber_ti_equiv", "new_c" + std::to_string(contact_force_type));
        e2e::ScopedSilence silence(true);
        time_integration(cfg, modelB, n, dtB, timeB, bc2m3::kDxSmoke, &guardB);
    }

    check(guardA.reason() == RunGuard::Reason::StepLimit,
          "старая сигнатура: остановлена штатно по StepLimit (не авария)");
    check(guardB.reason() == RunGuard::Reason::StepLimit,
          "новая сигнатура (fiber::Model, пустые нити): остановлена штатно по StepLimit");
    check(guardA.steps() == kSteps && guardB.steps() == kSteps,
          "обе сигнатуры сделали одинаковое число шагов (" + std::to_string(kSteps) + ")");

    check(timeA == timeB, "время расчёта совпадает после N шагов: " + std::to_string(timeA) +
                               " == " + std::to_string(timeB));
    check(dtA == dtB, "итоговый dt совпадает после N шагов (CFL/контакт -- общая формула в обеих сигнатурах)");

    const double tol = 1e-9;
    const FieldDiff diff = compare_solid(sysA, modelB.solid, n, tol);
    std::ostringstream os;
    os << "Solid-состояние совпадает со старой сигнатурой (<= " << tol << " относительно): худшее поле "
       << (diff.field.empty() ? "-" : diff.field) << ", частица " << diff.particle
       << ", относительное расхождение " << diff.worst;
    check(diff.field.empty(), os.str());
}

} // namespace

int main()
{
    std::printf("=== Эквивалентность Solid-цикла: time_integration(ParticleSystem&) vs "
                "time_integration(fiber::Model&) при пустых нитях ===\n");
    std::printf("(fiber/PLAN.md, этап Ф3 -- закрывает риск непроверенной копии Solid-цикла,\n"
                " отмеченный автором этапа как компромисс)\n");

    for (int contact_type : {0, 1, 3})
    {
        run_case(contact_type);
    }

    std::printf("\n=== %s (провалов: %d) ===\n", e2e::failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН",
                e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
