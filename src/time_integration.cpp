// time_integration.cpp
// Определения функций временной интеграции.
// Объявления — в time_integration.hpp.
#include "time_integration.hpp"
#include "output.hpp"
#include <fstream>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

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

    // переиспользуемые рабочие буферы (живут всю симуляцию, передаются в int_force)
    WorkBuffers wb;
    // граф пар + окраска (живёт всю симуляцию, перестраивается каждый шаг в single_step)
    PairGraph graph;
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
        single_step(cfg, sys, current_ts, dt, ntotal, t, tdsdt, dx, dvx, du, ds, drho, av, dsxx, dsyy, dszz, dsxy, dsyz, dsxz, wb, graph);
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
        for (size_t i = 0; i < ntotal; ++i)
        {
            if (sys.itype[i] >= 3 && sys.itype[i] <= 9)
            {
                prochnost_mises(cfg.mats, static_cast<size_t>(sys.itype[i]), sys.sxx[i], sys.syy[i], sys.szz[i], sys.sxy[i], sys.syz[i], sys.sxz[i]);
            }
            else if (sys.itype[i] == 10)
            {
                prochnost_cem_JH2(cfg.mats, sys, i, dt, drho[i]);
            }
            else if (sys.itype[i] == 11)
            {
                prochnost_cem_JH2(cfg.mats, sys, i, dt, drho[i]);
            }
        }
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
        const size_t moni_axis = (dim > 1) ? 1 : 0;
        file << time << " " << sys.vx[moni_axis][moni_particle] << " " << sys.itype[moni_particle] << std::endl;
        file.close();

        // Сторож времени/застоя (run_guard.hpp). Без него цикл при коллапсе
        // шага крутится бесконечно: время задачи практически не растёт, а
        // признаков зависания для ОС нет. Причина обрыва остаётся в guard.
        if (guard != nullptr && guard->should_stop(dt, time, t_total))
        {
            std::cerr << guard->message() << std::endl;
            break;
        }
    }
}
