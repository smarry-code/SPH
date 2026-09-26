// time_integration.cpp
// Определения функций временной интеграции.
// Объявления — в time_integration.hpp.
#include "time_integration.hpp"
#include "output.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "interactions/postprocessing.hpp"
#include "ode/ode_solver.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include <fstream>
#include <iostream>
#include <limits>

#ifdef _WIN32
#define NOMINMAX // std::min/std::max (rhs/contributions.hpp, interactions/fiber_bond.hpp) конфликтуют с макросами min/max windows.h
#include <windows.h>
#endif

// Вычисление в real_t (precision.hpp), результат -- double: dt живёт в часах
// расчёта (double), float -> double без потерь.
double compute_cfl_timestep(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal)
{
    const real_t cfl = static_cast<real_t>(cfg.CFL);
    real_t min_dt = real_t(1e6);              // шаг по куранту
    for (size_t i = 0; i < ntotal; ++i)
    {
        real_t v = real_t(0.0);
        for (size_t d = 0; d < dim; ++d)
        {
             v +=  sys.vx[d][i] * sys.vx[d][i];
        }
        v = std::sqrt(v);  // Скорость частицы
        real_t dt = cfl * sys.hsml[i] / (sys.c[i] + v);
        if (dt < min_dt) min_dt = dt;
    }
    return static_cast<double>(min_dt);
}

namespace detail
{
// ----------------------------------------------------------------------
// run_time_integration<Solver> -- ЕДИНСТВЕННАЯ реализация главного цикла по
// времени (fiber/PLAN.md, этап Ф8, раздел E «Унификация time_integration»).
// Обе публичные перегрузки time_integration() ниже -- тонкие обёртки над этим
// шаблоном: разница между "Solid-only" (старая сигнатура, ParticleSystem&) и
// "Solid+Fiber" (fiber::Model&) -- только в том, какие fb/topology/materials
// передаются (пустые по умолчанию -> nfiber == 0 -> вся Fiber-физика ниже
// no-op по существующим `if (nfiber > 0)`-веткам, как и раньше в перегрузке
// для fiber::Model). Никакой копии ParticleSystem не возникает — оба
// вызывающих места передают его по ссылке.
//
// Solver -- ode::PredictorCorrectorSolver (по умолчанию, fiber/PLAN.md §Ф8)
// или ode::EulerSolver (include/ode/ode_solver.hpp); выбор — по
// cfg.ode_solver, ОДИН раз перед циклом while, на стороне вызывающих
// перегрузок. Интерфейс обоих решателей единообразен (predict()/advance()),
// поэтому тело цикла не содержит if/switch по типу решателя (fiber/PLAN.md
// §3.9: «Дизайн — выбор через if/switch по enum на этапе вызова» — здесь этот
// выбор целиком сведён к выбору ШАБЛОННОГО ПАРАМЕТРА при инстанцировании).
//
// Порядок одного шага (fiber/SPH.md, «Описание алгоритмов и структур
// данных» / PLAN.md §3.9): solver.predict() -> RHS (single_step + Fiber-Bond/
// Fiber-Fiber-Contact/Fiber-Solid-Contact, ОДНО вычисление на шаг) -> новый
// dt -> solver.advance() -> Postprocessing (apply_material_strength, разрыв
// связей Fiber-Bond) -> вывод. update_smoothing_length -- ИСКЛЮЧЕНИЕ, см.
// большой комментарий в include/interactions/postprocessing.hpp (шапка
// файла) -- почему он остаётся МЕЖДУ RHS и вычислением нового dt, а не после
// advance().
// ----------------------------------------------------------------------
template <class Solver>
void run_time_integration(const SimulationConfig& cfg, ParticleSystem& sys, fiber::FiberState& fb,
                           fiber::FiberTopology& topology, std::vector<fiber::FiberMaterial>& materials,
                           size_t ntotal, double& dt, double& time, double /*delta_x*/, RunGuard* guard)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    const size_t nfiber = fb.size();

    const double t_total = cfg.t_total;
    const size_t print_step = cfg.print_step;
    const double save_t = cfg.save_t;
    const size_t moni_particle = cfg.moni_particle;

    // Solid RHS-буферы охватывают реальные И виртуальные частицы: single_step
    // вызывает узлы конвейера с ntotal + nvirt, и int_force/art_visc/av_vel
    // индексируют эти массивы всем диапазоном. Число зеркальных частиц
    // определяется каждый шаг заново, поэтому берётся запас maxn (main.cpp:
    // maxn = ntotal * 1.4 «под возможные виртуальные частицы»). Интегрируются
    // по-прежнему только реальные частицы [0, ntotal).
    const size_t nbuf = (cfg.maxn > ntotal) ? cfg.maxn : ntotal;

    Solver solver;
    solver.reset(nbuf, ntotal, nfiber);

    // переиспользуемые рабочие буферы (живут всю симуляцию, передаются в single_step)
    WorkBuffers wb;

    // --- Fiber RHS-инфраструктура (fiber/PLAN.md, этапы Ф3/Ф5/Ф6/Ф7):
    // per-thread буферы contributions/событий, переиспользуются между шагами
    // ("память только растёт", как и остальные буферы RhsWorkspace). Не
    // относятся к состоянию решателя (не f^{n-1}) -- это инфраструктура
    // ВЫЧИСЛЕНИЯ ОДНОГО RHS, а не история интегрирования по времени.
    rhs::CpuExecutor fiber_exec;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_thread_buffers;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_contact_thread_buffers;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_solid_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fiber_solid_solid_buffers;
    std::vector<std::vector<interactions::FiberBreakEvent>> fiber_break_thread_buffers;

    size_t current_ts = 0;
    double next_output_time = 0.0;

    // Основной цикл по времени
    while (time < t_total)
    {
        current_ts += 1;
        // Вывод прогресса
        if (current_ts % print_step == 0)
        {
            std::wcout << "______________________________________" << std::endl;
            std::wcout << "  Текущий щаг по времени = " << dt
                      << "    текущее время расчета = " << time + dt
                      << "    текущий цикл = " << current_ts << std::endl;
            if (cfg.contact_force_type != 0)
            {
                const bool limited_by_contact = (wb.min_dt_contact < 1e30) &&
                    (wb.min_dt_contact <= compute_cfl_timestep(cfg, sys, ntotal));
                std::wcout << "  Шаг ограничен " << (limited_by_contact ? "КОНТАКТОМ" : "CFL")
                          << "    contact_pair_count = " << wb.contact_pair_count
                          << "    contact_max_delta = " << wb.contact_max_delta << std::endl;
            }
            std::wcout << "______________________________________" << std::endl;
        }

        // OdeSolver: предиктор (State::save + полушаг по f^{n-1}; нет-оп у
        // Эйлера и на первом шаге PredictorCorrector).
        solver.predict(sys, ntotal, fb, nfiber, static_cast<real_t>(dt));

        // RHS: Solid (SPH+контакт, как раньше) -- пишет f^n в solver.solid_rhs.
        single_step(cfg, sys, current_ts, static_cast<real_t>(dt), ntotal, solver.solid_rhs.t, solver.solid_rhs.tdsdt, solver.solid_rhs.dx,
                    solver.solid_rhs.dvx, solver.solid_rhs.du, solver.solid_rhs.ds, solver.solid_rhs.drho,
                    solver.solid_rhs.av, solver.solid_rhs.dsxx, solver.solid_rhs.dsyy, solver.solid_rhs.dszz,
                    solver.solid_rhs.dsxy, solver.solid_rhs.dsyz, solver.solid_rhs.dsxz, wb);

        // Postprocessing, этап Ф7 (fiber/PLAN.md §3.8): sle=1 -- пересчёт h.
        // МЕСТО ВЫЗОВА -- сразу после single_step(), ДО compute_cfl_timestep():
        // см. большой комментарий у interactions::update_smoothing_length
        // (include/interactions/postprocessing.hpp) -- почему это обязательно
        // для бит-в-бит совпадения с прежним поведением (h_upgrade раньше
        // вызывался отсюда же, изнутри single_step.cpp).
        interactions::update_smoothing_length(cfg, sys, static_cast<real_t>(dt), ntotal);

        // RHS: Fiber-Bond (fiber/PLAN.md, этап Ф3) — readonly search по
        // topology.bonds, contributions -> reduce (fiber/SPH.md, «Основные
        // требования к interaction kernel»).
        // RHS: Fiber-Fiber-Contact (fiber/PLAN.md, этап Ф5) — readonly
        // spatial search по бусинам РАЗНЫХ нитей (interactions/fiber_contact.hpp),
        // тот же FiberDelta reduce, СКЛАДЫВАЕТСЯ с вкладом Fiber-Bond (обе
        // физики независимо производят ускорение бусин — fiber/PLAN.md §3.3:
        // «FiberDelta: dv[dim]», сумма нескольких interaction types в общий dv).
        double fiber_dt_limit = std::numeric_limits<double>::max();
        if (nfiber > 0)
        {
            std::vector<fiber::FiberDelta> fiber_bond_delta =
                interactions::compute_fiber_bond(fiber_exec, fb, topology, materials, fiber_thread_buffers);

            interactions::FiberContactSearchResult fiber_contact_search;
            interactions::build_fiber_contact_search(fiber_exec, fb, topology, materials, fiber_contact_search);
            interactions::compute_fiber_contact_contributions(fiber_exec, fb, materials,
                                                               fiber_contact_search.primitives, fiber_contact_thread_buffers);
            std::vector<fiber::FiberDelta> fiber_contact_delta =
                rhs::reduce_contributions(fiber_exec, fiber_contact_thread_buffers, fb.size());

            for (size_t d = 0; d < dim; ++d)
                for (size_t i = 0; i < nfiber; ++i)
                    solver.fiber_rhs.dv[d][i] = fiber_bond_delta[i].dv[d] + fiber_contact_delta[i].dv[d];

            // dt_limit взаимодействий -- real_t, итоговый dt -- double (precision.hpp).
            fiber_dt_limit = static_cast<double>(interactions::compute_fiber_bond_dt_limit(fiber_exec, fb, topology, materials));
            const double fiber_contact_dt_limit = static_cast<double>(interactions::compute_fiber_contact_dt_limit(
                fiber_exec, fb, materials, fiber_contact_search.primitives));
            if (fiber_contact_dt_limit < fiber_dt_limit) fiber_dt_limit = fiber_contact_dt_limit;

            // RHS: Fiber-Solid-Contact (fiber/PLAN.md, этап Ф6) -- механический
            // контакт бусины нити с реальной Solid-частицей. Вклад в Solid идёт
            // ЧЕРЕЗ ТОТ ЖЕ МЕХАНИЗМ, что и Solid-SPH/Solid-Solid-Contact внутри
            // single_step (SolidDelta + rhs::reduce_contributions) — здесь он
            // лишь ДОБАВЛЯЕТСЯ к уже посчитанному single_step() dvx
            // (interactions/fiber_solid_contact.hpp, шапка файла, решения 1-6).
            // Вклад в Fiber складывается с Fiber-Bond/Fiber-Fiber-Contact в тот
            // же dv, как и они между собой.
            const interactions::FiberSolidContactResult fiber_solid_result = interactions::compute_fiber_solid_contact(
                fiber_exec, cfg.mats, sys, ntotal, fb, topology, materials, fiber_solid_fiber_buffers,
                fiber_solid_solid_buffers);

            for (size_t d = 0; d < dim; ++d)
                for (size_t i = 0; i < nfiber; ++i)
                    solver.fiber_rhs.dv[d][i] += fiber_solid_result.fiber_delta[i].dv[d];
            for (size_t d = 0; d < dim; ++d)
                for (size_t i = 0; i < ntotal; ++i)
                    solver.solid_rhs.dvx[d][i] += fiber_solid_result.solid_delta[i].dv[d];

            const double fiber_solid_dt_limit = static_cast<double>(interactions::compute_fiber_solid_contact_dt_limit(
                fiber_exec, cfg.mats, sys, fb, materials, fiber_solid_result.primitives));
            if (fiber_solid_dt_limit < fiber_dt_limit) fiber_dt_limit = fiber_solid_dt_limit;
        }

        // Общий шаг (fiber/SPH.md, «Шаг по времени»): min(CFL, контакт Solid, dt_limit связей).
        dt = compute_cfl_timestep(cfg, sys, ntotal);
        if (wb.min_dt_contact < dt) dt = wb.min_dt_contact;
        if (fiber_dt_limit < dt) dt = fiber_dt_limit;

        // OdeSolver: State::advance -- обновление на полный шаг (Solid: u/rho/s/
        // v/x с правилами, специфичными для Solid (u>=0, XSPH), — внутри
        // решателя, fiber/PLAN.md §3.9; Fiber: только x, v).
        solver.advance(sys, ntotal, fb, nfiber, static_cast<real_t>(dt));

        // Postprocessing, этап Ф7 (fiber/PLAN.md §3.8) -- ПОСЛЕ шага
        // ODESolver (predictor-corrector/Euler Solid И Fiber уже завершены
        // выше):
        //   1) Мизес/JH-2 -- interactions::apply_material_strength, буквально
        //      прежний код, читает финальные sys.sxx.. (место вызова не
        //      изменилось);
        //   2) разрыв связей Fiber-Bond -- interactions::
        //      run_fiber_topology_postprocessing: поиск событий (readonly
        //      parallel_for по topology.bonds, ТЕКУЩИЕ fb.x после advance()
        //      выше) -> применение строго последовательно (prev/next,
        //      новый object_id j-стороне, topology.bonds/version). Fiber-Bond
        //      kernel и Fiber-Fiber-Contact search строят primitives заново
        //      из topology/object_id на КАЖДОМ вызове RHS -- отдельной
        //      инвалидации кэша не требуется (fiber/SPH.md, «Postprocessing и
        //      изменение topology»).
        interactions::apply_material_strength(cfg, sys, ntotal, static_cast<real_t>(dt), solver.solid_rhs.drho);
        if (nfiber > 0)
        {
            interactions::run_fiber_topology_postprocessing(fiber_exec, fb, topology, fiber_break_thread_buffers);
        }
        if (time >= next_output_time)
        {
            output(sys, ntotal, time, current_ts, dt);
            next_output_time += save_t;
        }
        time += dt;

        // Вывод информации о мониторинге
        std::ofstream file("log.txt", std::ios::app);
        // Пишется ОСЕВАЯ компонента. В осесимметричной постановке ось удара —
        // x[1] (z), и именно осевую скорость поверхности измеряет интерферометр,
        // с которым сверяется эталон. Осевая компонента — последняя координата:
        // z в 3D, z (ось симметрии) в осесимметричном 2D, x в 1D.
        const size_t moni_axis = dim - 1;
        file << time << " " << sys.vx[moni_axis][moni_particle] << " " << sys.itype[moni_particle] << std::endl;
        file.close();

        // Сторож времени/застоя (run_guard.hpp). Без него цикл при коллапсе
        // шага крутится бесконечно: время задачи практически не растёт, а
        // признаков зависания для ОС нет. Причина обрыва остаётся в guard.
        if (guard != nullptr && guard->should_stop(dt, time, t_total))
        {
            // Штатный останов по числу шагов (StepLimit) аварией не является.
            if (guard->failed()) std::cerr << guard->message() << std::endl;
            break;
        }
    }
}

} // namespace detail

// Основная функция временной интеграции (fiber/PLAN.md, этап Ф8: тонкая
// обёртка над detail::run_time_integration с пустыми Fiber-структурами по
// умолчанию -- FiberState::n == 0, FiberTopology::bonds.empty(), поэтому вся
// Fiber-физика внутри общего цикла no-op, как и раньше в перегрузке для
// fiber::Model). Выбор решателя -- по cfg.ode_solver (0 = PredictorCorrector,
// по умолчанию; 1 = Euler), инстанцирует run_time_integration нужным типом.
void time_integration(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, double& dt, double& time, double delta_x,
                      RunGuard* guard)
{
    fiber::FiberState empty_fb;
    fiber::FiberTopology empty_topology;
    std::vector<fiber::FiberMaterial> empty_materials;
    if (cfg.ode_solver == 1)
        detail::run_time_integration<ode::EulerSolver>(cfg, sys, empty_fb, empty_topology, empty_materials, ntotal, dt,
                                                        time, delta_x, guard);
    else
        detail::run_time_integration<ode::PredictorCorrectorSolver>(cfg, sys, empty_fb, empty_topology,
                                                                     empty_materials, ntotal, dt, time, delta_x, guard);
}

// Перегрузка для fiber::Model (этап Ф2: fiber/PLAN.md пункт 7 — сигнатура;
// этап Ф3+: физика Fiber-Bond/Fiber-Fiber-Contact/Fiber-Solid-Contact/
// Postprocessing). Тонкая обёртка над detail::run_time_integration — та же
// функция, что и у старой сигнатуры выше, различаются только передаваемые
// model.fiber/model.topology/model.materials вместо пустых структур.
void time_integration(const SimulationConfig& cfg, fiber::Model& model, size_t ntotal, double& dt, double& time,
                      double delta_x, RunGuard* guard)
{
    if (cfg.ode_solver == 1)
        detail::run_time_integration<ode::EulerSolver>(cfg, model.solid, model.fiber, model.topology, model.materials,
                                                        ntotal, dt, time, delta_x, guard);
    else
        detail::run_time_integration<ode::PredictorCorrectorSolver>(cfg, model.solid, model.fiber, model.topology,
                                                                     model.materials, ntotal, dt, time, delta_x, guard);
}
