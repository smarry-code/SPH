// time_integration.cpp
// Определения функций временной интеграции.
// Объявления — в time_integration.hpp.
#include "time_integration.hpp"
#include "output.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "interactions/postprocessing.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include <fstream>
#include <iostream>
#include <limits>

#ifdef _WIN32
#define NOMINMAX // std::min/std::max (rhs/contributions.hpp, interactions/fiber_bond.hpp) конфликтуют с макросами min/max windows.h
#include <windows.h>
#endif

namespace
{
// ----------------------------------------------------------------------
// Общая часть PredictorCorrector-схемы (fiber/SPH.md, «ODESolver»;
// fiber/PLAN.md §3.9) для произвольного набора x[dim]/vx[dim]: сохранение
// v^n и полушаг предиктора, затем финальное обновление скорости/координаты.
// Используется НИЖЕ, в перегрузке time_integration(fiber::Model&, ...), и
// для model.solid, и для model.fiber (fiber/PLAN.md, этап Ф3, пункт 3:
// «нити интегрируются... той же схемой, что Solid»).
//
// Шаблон по типу контейнера: у Solid внутри time_integration производные
// (dvx, av) — std::vector<std::vector<double>>, у ParticleSystem/FiberState
// x/vx — std::vector<double>[dim]; обоим достаточно operator[][], поэтому
// один шаблон обслуживает оба случая без адаптеров.
//
// TODO(Ф8): старый цикл выше (перегрузка time_integration(ParticleSystem&,
// ...), регрессионно протестированная BC2/BC6) НЕ использует эти функции —
// трогать его формулы ради унификации сейчас рискованно (малейшая
// перестановка операций способна изменить последний бит результата и
// сломать `ctest -L regression`). Поэтому Solid-цикл внутри
// time_integration(fiber::Model&, ...) ниже — временная копия структуры
// старого цикла (без u/rho/s — они специфичны для Solid и не переносятся в
// общий helper). Настоящая унификация обоих циклов через единую структуру
// OdeSolver::step — задача этапа Ф8 (fiber/PLAN.md §3.9, §4).
// ----------------------------------------------------------------------

// v_min = vx; vx += (dt_prev/2) * dv. Координаты не трогаются.
template <class N, class VxArr, class DvArr>
void pc_predict_velocity_half_step(N n, double dt_prev, VxArr& vx, const DvArr& dv, VxArr& v_min)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            v_min[d][i] = vx[d][i];
            vx[d][i] += (dt_prev / 2.0) * dv[d][i];
        }
    }
}

// vx = base + coeff*dv + av; x += dt*vx (Solid — с поправкой XSPH; av и dv
// могут иметь свой тип контейнера, отдельный от VxArr самих x/vx, — у Solid
// внутри time_integration это std::vector<std::vector<double>>, привязанный
// сигнатурой single_step).
// use_base == false -> base = текущая vx (первый шаг, current_ts == 1;
//   fiber/SPH.md «Первый шаг»: v += dt/2*v̇^0 + ṽ, coeff = dt/2);
// use_base == true  -> base = v_base (обычно v_min, сохранённая v^n;
//   fiber/SPH.md «Корректор»: v = v^n + dt*v̇^n + ṽ, coeff = dt).
template <class N, class XArr, class VxArr, class DvArr, class AvArr>
void pc_finalize_velocity_position_with_av(N n, double dt, double coeff, XArr& x, VxArr& vx, const DvArr& dv,
                                            const VxArr& v_base, bool use_base, const AvArr& av)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            const double base = use_base ? v_base[d][i] : vx[d][i];
            vx[d][i] = base + coeff * dv[d][i] + av[d][i];
            x[d][i] += dt * vx[d][i];
        }
    }
}

// То же самое без поправки XSPH (Fiber — fiber/SPH.md, «finalize»: у Fiber
// XSPH не определён, только ẋ = v).
template <class N, class XArr, class VxArr, class DvArr>
void pc_finalize_velocity_position(N n, double dt, double coeff, XArr& x, VxArr& vx, const DvArr& dv,
                                    const VxArr& v_base, bool use_base)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            const double base = use_base ? v_base[d][i] : vx[d][i];
            vx[d][i] = base + coeff * dv[d][i];
            x[d][i] += dt * vx[d][i];
        }
    }
}
} // namespace

double compute_cfl_timestep(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal)
{
    double min_dt = 1e6;                      // шаг по куранту
    for (size_t i = 0; i < ntotal; ++i)
    {
        double v = 0.0;
        for (size_t d = 0; d < dim; ++d)
        {
             v +=  sys.vx[d][i] * sys.vx[d][i];
        }
        v = std::sqrt(v);  // Скорость частицы
        double dt = cfg.CFL * sys.hsml[i] / (sys.c[i] + v);
        if (dt < min_dt) min_dt = dt;
    }
    return min_dt;
}
// Основная функция временной интеграции
void time_integration(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, double& dt, double& time, double delta_x,
                      RunGuard* guard)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    // Локальные псевдонимы параметров (тело цикла не изменяется):
    // t_total, print_step, save_t, moni_particle ранее — глобалы param.hpp.
    const double t_total = cfg.t_total;
    const size_t print_step = cfg.print_step;
    const double save_t = cfg.save_t;
    const size_t moni_particle = cfg.moni_particle;
    // Инициализация переменных.
    // Массивы производных охватывают реальные И виртуальные частицы: single_step
    // вызывает узлы конвейера с ntotal + nvirt, и int_force/art_visc/av_vel
    // индексируют эти массивы всем диапазоном. Число зеркальных частиц
    // определяется каждый шаг заново, поэтому берётся запас maxn — он для того
    // и заводится (main.cpp: maxn = ntotal * 1.4 «под возможные виртуальные
    // частицы»). Интегрируются по-прежнему только реальные частицы [0, ntotal).
    const size_t nbuf = (cfg.maxn > ntotal) ? cfg.maxn : ntotal;
    std::vector<std::vector<double>> av(dim, std::vector<double>(nbuf, 0.0));
    std::vector<std::vector<double>> dx(dim, std::vector<double>(nbuf, 0.0));
    std::vector<std::vector<double>> dvx(dim, std::vector<double>(nbuf, 0.0)); // Полное ускорение
    std::vector<double> du(nbuf, 0.0);  // Скорость изменения удельной внутренней
    std::vector<double> drho(nbuf, 0.0);
    std::vector<double> ds(nbuf, 0.0);
    std::vector<double> t(nbuf, 0.0);
    std::vector<double> tdsdt(nbuf, 0.0);

    std::vector<double> dsxx(nbuf, 0.0);
    std::vector<double> dsyy(nbuf, 0.0);
    std::vector<double> dszz(nbuf, 0.0);
    std::vector<double> dsxy(nbuf, 0.0);
    std::vector<double> dsyz(nbuf, 0.0);
    std::vector<double> dsxz(nbuf, 0.0);

    std::vector<double> sxx_min(ntotal, 0.0);
    std::vector<double> syy_min(ntotal, 0.0);
    std::vector<double> szz_min(ntotal, 0.0);
    std::vector<double> sxy_min(ntotal, 0.0);
    std::vector<double> syz_min(ntotal, 0.0);
    std::vector<double> sxz_min(ntotal, 0.0);

    std::vector<std::vector<double>> x_min(dim, std::vector<double>(ntotal, 0.0));
    std::vector<std::vector<double>> v_min(dim, std::vector<double>(ntotal, 0.0));
    std::vector<double> u_min(ntotal, 0.0);
    std::vector<double> rho_min(ntotal, 0.0);

    // переиспользуемые рабочие буферы (живут всю симуляцию, передаются в single_step)
    WorkBuffers wb;
    // sys (ParticleSystem) — первичное хранилище частиц, передаётся из main.

    size_t current_ts = 0;

    // Основной цикл по времени
    double next_output_time = 0.0;
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

        // Обновление термальной энергии, плотности и скорости (половина шага)
        if (current_ts != 1)
        {
            for (size_t i = 0; i < ntotal; ++i)
            {
                u_min[i] = sys.u[i];
                double temp_u = 0.0;

                if (dim == 1)
                {
                    temp_u = -nsym * sys.p[i] * sys.vx[0][i] / sys.x[0][i] / sys.rho[i];
                }

                sys.u[i] += (dt/2.0) * (du[i] + temp_u);
                if (sys.u[i] < 0) sys.u[i] = 0.0;

                if (!summation_density)
                {
                    rho_min[i] = sys.rho[i];
                    double temp_rho = 0.0;

                    if (dim == 1)
                    {
                        temp_rho = -nsym * sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
                    }

                    sys.rho[i] += (dt/2.0) * (drho[i] + temp_rho);
                }

                for (size_t d = 0; d < dim; ++d)
                {
                    v_min[d][i] = sys.vx[d][i];
                    sys.vx[d][i] += (dt/2.0) * dvx[d][i];
                }
                sxx_min[i] = sys.sxx[i];
                sxy_min[i] = sys.sxy[i];
                syy_min[i] = sys.syy[i];
                syz_min[i] = sys.syz[i];
                sxz_min[i] = sys.sxz[i];
                szz_min[i] = sys.szz[i];
                sys.sxx[i] += (dt/2.0) * dsxx[i];
                sys.syy[i] += (dt/2.0) * dsyy[i];
                sys.szz[i] += (dt/2.0) * dszz[i];
                sys.sxy[i] += (dt/2.0) * dsxy[i];
                sys.syz[i] += (dt/2.0) * dsyz[i];
                sys.sxz[i] += (dt/2.0) * dsxz[i];
            }
        }

        // Вычисление сил и производных
        single_step(cfg, sys, current_ts, dt, ntotal, t, tdsdt, dx, dvx, du, ds, drho, av, dsxx, dsyy, dszz, dsxy, dsyz, dsxz, wb);
        // Postprocessing, этап Ф7 (fiber/PLAN.md §3.8): sle=1 -- пересчёт h.
        // МЕСТО ВЫЗОВА -- сразу после single_step(), ДО compute_cfl_timestep():
        // см. большой комментарий у interactions::update_smoothing_length
        // (include/interactions/postprocessing.hpp) -- почему это обязательно
        // для бит-в-бит совпадения с прежним поведением (h_upgrade раньше
        // вызывался отсюда же, изнутри single_step.cpp).
        interactions::update_smoothing_length(cfg, sys, dt, ntotal);
        dt = compute_cfl_timestep(cfg, sys, ntotal);
        if (wb.min_dt_contact < dt) dt = wb.min_dt_contact;   // ADR-007: контактный предел
        // Обновление свойств частиц (полный временной шаг)
        if (current_ts == 1)
        {
            for (size_t i = 0; i < ntotal; ++i)
            {
                double temp_u = 0.0;

                if (dim == 1)
                {
                    temp_u = -nsym * sys.p[i] * sys.vx[0][i] / sys.x[0][i] / sys.rho[i];
                }

                sys.u[i] += (dt/2.0) * (du[i] + temp_u);
                if (sys.u[i] < 0) sys.u[i] = 0.0;

                if (!summation_density)
                {
                    double temp_rho = 0.0;

                    if (dim == 1)
                    {
                        temp_rho = -nsym * sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
                    }

                    sys.rho[i] += (dt/2.0) * (drho[i] + temp_rho);
                }

                for (size_t d = 0; d < dim; ++d)
                {                                       
                    if (d == 0)
                    {
                        if (sys.axis[i] == true)
                        {
                            sys.vx[d][i] += (dt/2.0) * dvx[d][i] + av[d][i];
                            sys.x[d][i] += dt * sys.vx[d][i];
                            /*if (sys.vx[d][i] < 0)
                            {
                                sys.vx[d][i] = 0.0;
                                sys.x[d][i] = 0.5 * delta_x;
                            }*/
                        }
                        else
                        {
                            sys.vx[d][i] += (dt/2.0) * dvx[d][i] + av[d][i];
                            sys.x[d][i] += dt * sys.vx[d][i];
                        }
                    }
                    else
                    {
                        sys.vx[d][i] += (dt/2.0) * dvx[d][i] + av[d][i];
                        sys.x[d][i] += dt * sys.vx[d][i];
                    }
                }
                sys.sxx[i] += (dt/2.0) * dsxx[i];
                sys.syy[i] += (dt/2.0) * dsyy[i];
                sys.szz[i] += (dt/2.0) * dszz[i];
                sys.sxy[i] += (dt/2.0) * dsxy[i];
                sys.syz[i] += (dt/2.0) * dsyz[i];
                sys.sxz[i] += (dt/2.0) * dsxz[i];
            }
        }
        else
        {
            for (size_t i = 0; i < ntotal; ++i)
            {
                double temp_u = 0.0;

                if (dim == 1)
                {
                    temp_u = -nsym * sys.p[i] * sys.vx[0][i] / sys.x[0][i] / sys.rho[i];
                }

                sys.u[i] = u_min[i] + dt * (du[i] + temp_u);
                if (sys.u[i] < 0) sys.u[i] = 0.0;

                if (!summation_density)
                {
                    double temp_rho = 0.0;

                    if (dim == 1)
                    {
                        temp_rho = -nsym * sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
                    }

                    sys.rho[i] = rho_min[i] + dt * (drho[i] + temp_rho);
                }

                for (size_t d = 0; d < dim; ++d)
                {                   
                    if (d == 0)
                    {
                        if (sys.axis[i] == true)
                        {
                            sys.vx[d][i] = v_min[d][i] + dt * dvx[d][i] + av[d][i];
                            sys.x[d][i] += dt * sys.vx[d][i];
                            /*if (sys.vx[d][i] < 0)
                            {
                                sys.vx[d][i] = 0.0;
                                sys.x[d][i] = 0.5 * delta_x;
                            }*/
                        }
                        else
                        {
                            sys.vx[d][i] = v_min[d][i] + dt * dvx[d][i] + av[d][i];
                            sys.x[d][i] += dt * sys.vx[d][i];
                        }
                    }
                    else
                    {
                        sys.vx[d][i] = v_min[d][i] + dt * dvx[d][i] + av[d][i];
                        sys.x[d][i] += dt * sys.vx[d][i];
                    }
                }
                sys.sxx[i] = sxx_min[i] + (dt) * dsxx[i];
                sys.syy[i] = syy_min[i] + (dt) * dsyy[i];
                sys.szz[i] = szz_min[i] + (dt) * dszz[i];
                sys.sxy[i] = sxy_min[i] + (dt) * dsxy[i];
                sys.syz[i] = syz_min[i] + (dt) * dsyz[i];
                sys.sxz[i] = sxz_min[i] + (dt) * dsxz[i];
            }
        }
        // Postprocessing, этап Ф7: Мизес/JH-2 (interactions::apply_material_strength,
        // include/interactions/postprocessing.hpp) -- буквально прежний код,
        // читает ФИНАЛЬНЫЕ (скорректированные) sys.sxx.. и т.д., место вызова
        // не изменилось (сразу после corrector'а, как и раньше).
        interactions::apply_material_strength(cfg, sys, ntotal, dt, drho);
        if (time >= next_output_time)
        {
            output(sys, ntotal, time, current_ts, dt);
            next_output_time += save_t;
        }
        time += dt;

        // Сохранение данных
        /*if (itimestep % save_step == 0)
        {
            output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal, time, current_ts);
        }*/

        // Вывод информации о мониторинге
        std::ofstream file("log.txt", std::ios::app);
        // Пишется ОСЕВАЯ компонента. В осесимметричной постановке ось удара —
        // x[1] (z), и именно осевую скорость поверхности измеряет интерферометр,
        // с которым сверяется эталон. Раньше писалась vx[0] (радиальная), из-за
        // чего расчёт и эксперимент были несопоставимы в принципе.
        // Прежние Керамика/**/log_*.txt содержат vx[0] и с новыми несравнимы.
        // Осевая компонента — последняя координата: z в 3D, z (ось симметрии)
        // в осесимметричном 2D, x в 1D. Прежняя запись (dim > 1 ? 1 : 0)
        // совпадает с этой при dim <= 2, а в 3D снимала бы поперечную v_y.
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

// Перегрузка для fiber::Model (этап Ф2: fiber/PLAN.md пункт 7 — сигнатура и
// тонкая обёртка; этап Ф3: fiber/PLAN.md, раздел «Ф3. Fiber-Bond», пункт 3 —
// «нити интегрируются в текущем цикле time_integration той же схемой, что
// Solid; шаг — общий минимум»).
//
// Структура цикла — копия цикла выше (Solid, старая сигнатура), с двумя
// добавлениями: (1) после single_step считается RHS связей Fiber-Bond
// (interactions/fiber_bond.hpp) и её dt_limit; (2) итоговый dt — минимум
// CFL/контакта Solid И dt_limit связей (fiber/SPH.md, «Шаг по времени»).
// Предиктор/корректор скорости и координаты (Solid и Fiber) используют общие
// шаблоны pc_predict_velocity_half_step/pc_finalize_velocity_position выше
// (TODO(Ф8) там же — почему старый Solid-цикл их не использует).
void time_integration(const SimulationConfig& cfg, fiber::Model& model, size_t ntotal, double& dt, double& time,
                      double delta_x, RunGuard* guard)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    ParticleSystem& sys = model.solid;
    fiber::FiberState& fb = model.fiber;
    const size_t nfiber = fb.size();

    const double t_total = cfg.t_total;
    const size_t print_step = cfg.print_step;
    const double save_t = cfg.save_t;
    const size_t moni_particle = cfg.moni_particle;

    // --- Solid: те же буферы, что в старом цикле (см. комментарии там) ---
    const size_t nbuf = (cfg.maxn > ntotal) ? cfg.maxn : ntotal;
    std::vector<std::vector<double>> av(dim, std::vector<double>(nbuf, 0.0));
    std::vector<std::vector<double>> dx(dim, std::vector<double>(nbuf, 0.0));
    std::vector<std::vector<double>> dvx(dim, std::vector<double>(nbuf, 0.0));
    std::vector<double> du(nbuf, 0.0);
    std::vector<double> drho(nbuf, 0.0);
    std::vector<double> ds(nbuf, 0.0);
    std::vector<double> t(nbuf, 0.0);
    std::vector<double> tdsdt(nbuf, 0.0);

    std::vector<double> dsxx(nbuf, 0.0);
    std::vector<double> dsyy(nbuf, 0.0);
    std::vector<double> dszz(nbuf, 0.0);
    std::vector<double> dsxy(nbuf, 0.0);
    std::vector<double> dsyz(nbuf, 0.0);
    std::vector<double> dsxz(nbuf, 0.0);

    std::vector<double> sxx_min(ntotal, 0.0);
    std::vector<double> syy_min(ntotal, 0.0);
    std::vector<double> szz_min(ntotal, 0.0);
    std::vector<double> sxy_min(ntotal, 0.0);
    std::vector<double> syz_min(ntotal, 0.0);
    std::vector<double> sxz_min(ntotal, 0.0);

    // v_min — того же типа, что sys.vx (std::vector<double>[dim]), а не
    // std::vector<std::vector<double>>, как в старом цикле выше: нужен общий
    // тип VxArr для pc_predict_velocity_half_step/pc_finalize_velocity_position*,
    // которыми пользуются и Solid (здесь), и Fiber (v_min_fiber ниже).
    std::vector<double> v_min[dim];
    for (size_t d = 0; d < dim; ++d) v_min[d].assign(ntotal, 0.0);
    std::vector<double> u_min(ntotal, 0.0);
    std::vector<double> rho_min(ntotal, 0.0);

    WorkBuffers wb;

    // --- Fiber: ускорение связей на шаге, v^n сохранённая для корректора ---
    std::vector<double> dv_fiber[dim];
    std::vector<double> v_min_fiber[dim];
    for (size_t d = 0; d < dim; ++d)
    {
        dv_fiber[d].assign(nfiber, 0.0);
        v_min_fiber[d].assign(nfiber, 0.0);
    }
    rhs::CpuExecutor fiber_exec;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_thread_buffers;
    // Fiber-Fiber-Contact (fiber/PLAN.md, этап Ф5) -- отдельный буфер
    // (interaction_type отличается, kFiberFiberContactInteractionType !=
    // kFiberBondInteractionType), переиспользуется между шагами так же, как
    // fiber_thread_buffers выше ("память только растёт").
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_contact_thread_buffers;
    // Fiber-Solid-Contact (fiber/PLAN.md, этап Ф6) -- два буфера (сторона
    // Fiber и сторона Solid, разные Delta), переиспользуются между шагами по
    // тому же принципу, что буферы выше.
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_solid_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fiber_solid_solid_buffers;
    // Postprocessing, этап Ф7 (fiber/PLAN.md §3.8) -- разрыв связей Fiber-Bond:
    // per-thread буфер событий поиска (interactions/postprocessing.hpp,
    // find_fiber_break_events), переиспользуется между шагами по тому же
    // принципу, что буферы выше.
    std::vector<std::vector<interactions::FiberBreakEvent>> fiber_break_thread_buffers;

    size_t current_ts = 0;
    double next_output_time = 0.0;

    while (time < t_total)
    {
        current_ts += 1;
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

        // Предиктор на полшага (по производным ПРЕДЫДУЩЕГО шага) — Solid и Fiber.
        if (current_ts != 1)
        {
            for (size_t i = 0; i < ntotal; ++i)
            {
                u_min[i] = sys.u[i];
                double temp_u = 0.0;
                if (dim == 1) temp_u = -nsym * sys.p[i] * sys.vx[0][i] / sys.x[0][i] / sys.rho[i];
                sys.u[i] += (dt / 2.0) * (du[i] + temp_u);
                if (sys.u[i] < 0) sys.u[i] = 0.0;

                if (!summation_density)
                {
                    rho_min[i] = sys.rho[i];
                    double temp_rho = 0.0;
                    if (dim == 1) temp_rho = -nsym * sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
                    sys.rho[i] += (dt / 2.0) * (drho[i] + temp_rho);
                }

                sxx_min[i] = sys.sxx[i];
                sxy_min[i] = sys.sxy[i];
                syy_min[i] = sys.syy[i];
                syz_min[i] = sys.syz[i];
                sxz_min[i] = sys.sxz[i];
                szz_min[i] = sys.szz[i];
                sys.sxx[i] += (dt / 2.0) * dsxx[i];
                sys.syy[i] += (dt / 2.0) * dsyy[i];
                sys.szz[i] += (dt / 2.0) * dszz[i];
                sys.sxy[i] += (dt / 2.0) * dsxy[i];
                sys.syz[i] += (dt / 2.0) * dsyz[i];
                sys.sxz[i] += (dt / 2.0) * dsxz[i];
            }
            pc_predict_velocity_half_step(ntotal, dt, sys.vx, dvx, v_min);
            pc_predict_velocity_half_step(nfiber, dt, fb.vx, dv_fiber, v_min_fiber);
        }

        // RHS: Solid (SPH+контакт, как раньше).
        single_step(cfg, sys, current_ts, dt, ntotal, t, tdsdt, dx, dvx, du, ds, drho, av, dsxx, dsyy, dszz, dsxy, dsyz, dsxz, wb);

        // Postprocessing, этап Ф7 (fiber/PLAN.md §3.8): sle=1 -- пересчёт h.
        // МЕСТО ВЫЗОВА -- сразу после single_step(), ДО compute_cfl_timestep():
        // см. большой комментарий у interactions::update_smoothing_length
        // (include/interactions/postprocessing.hpp) -- почему это обязательно
        // для бит-в-бит совпадения с прежним поведением.
        interactions::update_smoothing_length(cfg, sys, dt, ntotal);

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
                interactions::compute_fiber_bond(fiber_exec, fb, model.topology, model.materials, fiber_thread_buffers);

            interactions::FiberContactSearchResult fiber_contact_search;
            interactions::build_fiber_contact_search(fiber_exec, fb, model.topology, model.materials, fiber_contact_search);
            interactions::compute_fiber_contact_contributions(fiber_exec, fb, model.materials,
                                                               fiber_contact_search.primitives, fiber_contact_thread_buffers);
            std::vector<fiber::FiberDelta> fiber_contact_delta =
                rhs::reduce_contributions(fiber_exec, fiber_contact_thread_buffers, fb.size());

            for (size_t d = 0; d < dim; ++d)
                for (size_t i = 0; i < nfiber; ++i)
                    dv_fiber[d][i] = fiber_bond_delta[i].dv[d] + fiber_contact_delta[i].dv[d];

            fiber_dt_limit = interactions::compute_fiber_bond_dt_limit(fiber_exec, fb, model.topology, model.materials);
            const double fiber_contact_dt_limit = interactions::compute_fiber_contact_dt_limit(
                fiber_exec, fb, model.materials, fiber_contact_search.primitives);
            if (fiber_contact_dt_limit < fiber_dt_limit) fiber_dt_limit = fiber_contact_dt_limit;

            // RHS: Fiber-Solid-Contact (fiber/PLAN.md, этап Ф6) -- механический
            // контакт бусины нити с реальной Solid-частицей. Ф4 уже завершён,
            // поэтому вклад в Solid идёт ЧЕРЕЗ ТОТ ЖЕ МЕХАНИЗМ, что и Solid-SPH/
            // Solid-Solid-Contact внутри single_step (SolidDelta + rhs::
            // reduce_contributions) — здесь он лишь ДОБАВЛЯЕТСЯ к уже
            // посчитанному single_step() dvx (interactions/fiber_solid_contact.hpp,
            // шапка файла, решения 1-6). Вклад в Fiber складывается с Fiber-Bond/
            // Fiber-Fiber-Contact в тот же dv_fiber, как и они между собой.
            const interactions::FiberSolidContactResult fiber_solid_result = interactions::compute_fiber_solid_contact(
                fiber_exec, cfg.mats, sys, ntotal, fb, model.topology, model.materials, fiber_solid_fiber_buffers,
                fiber_solid_solid_buffers);

            for (size_t d = 0; d < dim; ++d)
                for (size_t i = 0; i < nfiber; ++i)
                    dv_fiber[d][i] += fiber_solid_result.fiber_delta[i].dv[d];
            for (size_t d = 0; d < dim; ++d)
                for (size_t i = 0; i < ntotal; ++i)
                    dvx[d][i] += fiber_solid_result.solid_delta[i].dv[d];

            const double fiber_solid_dt_limit = interactions::compute_fiber_solid_contact_dt_limit(
                fiber_exec, cfg.mats, sys, fb, model.materials, fiber_solid_result.primitives);
            if (fiber_solid_dt_limit < fiber_dt_limit) fiber_dt_limit = fiber_solid_dt_limit;
        }

        // Общий шаг (fiber/SPH.md, «Шаг по времени»): min(CFL, контакт Solid, dt_limit связей).
        dt = compute_cfl_timestep(cfg, sys, ntotal);
        if (wb.min_dt_contact < dt) dt = wb.min_dt_contact;
        if (fiber_dt_limit < dt) dt = fiber_dt_limit;

        // Обновление на полный шаг — Solid (u/rho/s по-прежнему построчно, как
        // в старом цикле) и Fiber (только x, v — fiber/SPH.md, «ODESolver»,
        // «Нити интегрируются той же схемой (только x, v)»).
        const bool first_step = (current_ts == 1);
        for (size_t i = 0; i < ntotal; ++i)
        {
            double temp_u = 0.0;
            if (dim == 1) temp_u = -nsym * sys.p[i] * sys.vx[0][i] / sys.x[0][i] / sys.rho[i];

            if (first_step)
            {
                sys.u[i] += (dt / 2.0) * (du[i] + temp_u);
            }
            else
            {
                sys.u[i] = u_min[i] + dt * (du[i] + temp_u);
            }
            if (sys.u[i] < 0) sys.u[i] = 0.0;

            if (!summation_density)
            {
                double temp_rho = 0.0;
                if (dim == 1) temp_rho = -nsym * sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
                if (first_step)
                    sys.rho[i] += (dt / 2.0) * (drho[i] + temp_rho);
                else
                    sys.rho[i] = rho_min[i] + dt * (drho[i] + temp_rho);
            }

            if (first_step)
            {
                sys.sxx[i] += (dt / 2.0) * dsxx[i];
                sys.syy[i] += (dt / 2.0) * dsyy[i];
                sys.szz[i] += (dt / 2.0) * dszz[i];
                sys.sxy[i] += (dt / 2.0) * dsxy[i];
                sys.syz[i] += (dt / 2.0) * dsyz[i];
                sys.sxz[i] += (dt / 2.0) * dsxz[i];
            }
            else
            {
                sys.sxx[i] = sxx_min[i] + dt * dsxx[i];
                sys.syy[i] = syy_min[i] + dt * dsyy[i];
                sys.szz[i] = szz_min[i] + dt * dszz[i];
                sys.sxy[i] = sxy_min[i] + dt * dsxy[i];
                sys.syz[i] = syz_min[i] + dt * dsyz[i];
                sys.sxz[i] = sxz_min[i] + dt * dsxz[i];
            }
        }
        pc_finalize_velocity_position_with_av(ntotal, dt, first_step ? dt / 2.0 : dt, sys.x, sys.vx, dvx, v_min,
                                               /*use_base=*/!first_step, av);
        pc_finalize_velocity_position(nfiber, dt, first_step ? dt / 2.0 : dt, fb.x, fb.vx, dv_fiber, v_min_fiber,
                                       /*use_base=*/!first_step);

        // Postprocessing, этап Ф7 (fiber/PLAN.md §3.8) -- ПОСЛЕ шага
        // ODESolver (predictor-corrector Solid И Fiber уже завершены выше):
        //   1) Мизес/JH-2 -- interactions::apply_material_strength, буквально
        //      прежний код, читает финальные sys.sxx.. (место вызова не
        //      изменилось);
        //   2) разрыв связей Fiber-Bond -- interactions::
        //      run_fiber_topology_postprocessing: поиск событий (readonly
        //      parallel_for по topology.bonds, ТЕКУЩИЕ fb.x после corrector'а
        //      выше) -> применение строго последовательно (prev/next,
        //      новый object_id j-стороне, topology.bonds/version). Fiber-Bond
        //      kernel и Fiber-Fiber-Contact search строят primitives заново
        //      из topology/object_id на КАЖДОМ вызове RHS -- отдельной
        //      инвалидации кэша не требуется (fiber/SPH.md, «Postprocessing и
        //      изменение topology»).
        interactions::apply_material_strength(cfg, sys, ntotal, dt, drho);
        if (nfiber > 0)
        {
            interactions::run_fiber_topology_postprocessing(fiber_exec, fb, model.topology, fiber_break_thread_buffers);
        }
        if (time >= next_output_time)
        {
            output(sys, ntotal, time, current_ts, dt);
            next_output_time += save_t;
        }
        time += dt;

        std::ofstream file("log.txt", std::ios::app);
        const size_t moni_axis = dim - 1;
        file << time << " " << sys.vx[moni_axis][moni_particle] << " " << sys.itype[moni_particle] << std::endl;
        file.close();

        if (guard != nullptr && guard->should_stop(dt, time, t_total))
        {
            if (guard->failed()) std::cerr << guard->message() << std::endl;
            break;
        }
    }
}
