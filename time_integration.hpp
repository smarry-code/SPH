#include "single_step.hpp"
#include "ppochnost.hpp"
double compute_cfl_timestep(std::vector<Particle>& Elem, size_t ntotal)
{
    double min_dt = 1e6;
    for (size_t i = 0; i < ntotal; ++i)
    {
        double v = 0.0;
        for (size_t d = 0; d < dim; ++d)
        {
             v +=  Elem[i].vx[d] * Elem[i].vx[d];
        }
        v = std::sqrt(v);  // Скорость частицы
        double dt = CFL * Elem[i].hsml / (Elem[i].c + v);
        if (dt < min_dt) min_dt = dt;
    }
    return min_dt;
}
// Основная функция временной интеграции
void time_integration(std::vector<Particle>& Elem, size_t ntotal, double& dt, double& time, double delta_x)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    // Инициализация переменных
    std::vector<std::vector<double>> av(dim, std::vector<double>(maxn, 0.0));
    std::vector<std::vector<double>> dx(dim, std::vector<double>(maxn, 0.0));
    std::vector<std::vector<double>> dvx(dim, std::vector<double>(maxn, 0.0));
    std::vector<double> du(maxn, 0.0);
    std::vector<double> drho(maxn, 0.0);
    std::vector<double> ds(maxn, 0.0);
    std::vector<double> t(maxn, 0.0);
    std::vector<double> tdsdt(maxn, 0.0);

    std::vector<double> dsxx(maxn, 0.0);
    std::vector<double> dsyy(maxn, 0.0);
    std::vector<double> dszz(maxn, 0.0);
    std::vector<double> dsxy(maxn, 0.0);
    std::vector<double> dsyz(maxn, 0.0);
    std::vector<double> dsxz(maxn, 0.0);

    std::vector<double> sxx_min(maxn, 0.0);
    std::vector<double> syy_min(maxn, 0.0);
    std::vector<double> szz_min(maxn, 0.0);
    std::vector<double> sxy_min(maxn, 0.0);
    std::vector<double> syz_min(maxn, 0.0);
    std::vector<double> sxz_min(maxn, 0.0);

    std::vector<std::vector<double>> x_min(dim, std::vector<double>(maxn, 0.0));
    std::vector<std::vector<double>> v_min(dim, std::vector<double>(maxn, 0.0));
    std::vector<double> u_min(maxn, 0.0);
    std::vector<double> rho_min(maxn, 0.0);

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
            std::wcout << "______________________________________" << std::endl;
        }

        // Обновление термальной энергии, плотности и скорости (половина шага)
        if (current_ts != 1)
        {
            for (size_t i = 0; i < ntotal; ++i)
            {
                u_min[i] = Elem[i].u;
                double temp_u = 0.0;

                if (dim == 1)
                {
                    temp_u = -nsym * Elem[i].p * Elem[i].vx[0] / Elem[i].x[0] / Elem[i].rho;
                }

                Elem[i].u += (dt/2.0) * (du[i] + temp_u);
                if (Elem[i].u < 0) Elem[i].u = 0.0;

                if (!summation_density)
                {
                    rho_min[i] = Elem[i].rho;
                    double temp_rho = 0.0;

                    if (dim == 1)
                    {
                        temp_rho = -nsym * Elem[i].rho * Elem[i].vx[0] / Elem[i].x[0];
                    }

                    Elem[i].rho += (dt/2.0) * (drho[i] + temp_rho);
                }

                for (size_t d = 0; d < dim; ++d)
                {
                    v_min[d][i] = Elem[i].vx[d];
                    Elem[i].vx[d] += (dt/2.0) * dvx[d][i];
                }
                sxx_min[i] = Elem[i].sxx;
                sxy_min[i] = Elem[i].sxy;
                syy_min[i] = Elem[i].syy;
                syz_min[i] = Elem[i].syz;
                sxz_min[i] = Elem[i].sxz;
                szz_min[i] = Elem[i].szz;
                Elem[i].sxx += (dt/2.0) * dsxx[i];
                Elem[i].syy += (dt/2.0) * dsyy[i];
                Elem[i].szz += (dt/2.0) * dszz[i];
                Elem[i].sxy += (dt/2.0) * dsxy[i];
                Elem[i].syz += (dt/2.0) * dsyz[i];
                Elem[i].sxz += (dt/2.0) * dsxz[i];
            }
        }

        // Вычисление сил и производных
        single_step(Elem, current_ts, dt, ntotal, t, tdsdt, dx, dvx, du, ds, drho, av, dsxx, dsyy, dszz, dsxy, dsyz, dsxz);
        dt = compute_cfl_timestep(Elem, ntotal);
        // Обновление свойств частиц (полный временной шаг)
        if (current_ts == 1)
        {
            for (size_t i = 0; i < ntotal; ++i)
            {
                double temp_u = 0.0;

                if (dim == 1)
                {
                    temp_u = -nsym * Elem[i].p * Elem[i].vx[0] / Elem[i].x[0] / Elem[i].rho;
                }

                Elem[i].u += (dt/2.0) * (du[i] + temp_u);
                if (Elem[i].u < 0) Elem[i].u = 0.0;

                if (!summation_density)
                {
                    double temp_rho = 0.0;

                    if (dim == 1)
                    {
                        temp_rho = -nsym * Elem[i].rho * Elem[i].vx[0] / Elem[i].x[0];
                    }

                    Elem[i].rho += (dt/2.0) * (drho[i] + temp_rho);
                }

                for (size_t d = 0; d < dim; ++d)
                {                                       
                    if (d == 0)
                    {
                        if (Elem[i].axis == true)
                        {
                            Elem[i].vx[d] += (dt/2.0) * dvx[d][i] + av[d][i];
                            Elem[i].x[d] += dt * Elem[i].vx[d];
                            /*if (Elem[i].vx[d] < 0)
                            {
                                Elem[i].vx[d] = 0.0;
                                Elem[i].x[d] = 0.5 * delta_x;
                            }*/
                        }
                        else
                        {
                            Elem[i].vx[d] += (dt/2.0) * dvx[d][i] + av[d][i];
                            Elem[i].x[d] += dt * Elem[i].vx[d];
                        }
                    }
                    else
                    {
                        Elem[i].vx[d] += (dt/2.0) * dvx[d][i] + av[d][i];
                        Elem[i].x[d] += dt * Elem[i].vx[d];
                    }
                }
                Elem[i].sxx += (dt/2.0) * dsxx[i];
                Elem[i].syy += (dt/2.0) * dsyy[i];
                Elem[i].szz += (dt/2.0) * dszz[i];
                Elem[i].sxy += (dt/2.0) * dsxy[i];
                Elem[i].syz += (dt/2.0) * dsyz[i];
                Elem[i].sxz += (dt/2.0) * dsxz[i];
            }
        }
        else
        {
            for (size_t i = 0; i < ntotal; ++i)
            {
                double temp_u = 0.0;

                if (dim == 1)
                {
                    temp_u = -nsym * Elem[i].p * Elem[i].vx[0] / Elem[i].x[0] / Elem[i].rho;
                }

                Elem[i].u = u_min[i] + dt * (du[i] + temp_u);
                if (Elem[i].u < 0) Elem[i].u = 0.0;

                if (!summation_density)
                {
                    double temp_rho = 0.0;

                    if (dim == 1)
                    {
                        temp_rho = -nsym * Elem[i].rho * Elem[i].vx[0] / Elem[i].x[0];
                    }

                    Elem[i].rho = rho_min[i] + dt * (drho[i] + temp_rho);
                }

                for (size_t d = 0; d < dim; ++d)
                {                   
                    if (d == 0)
                    {
                        if (Elem[i].axis == true)
                        {
                            Elem[i].vx[d] = v_min[d][i] + dt * dvx[d][i] + av[d][i];
                            Elem[i].x[d] += dt * Elem[i].vx[d];
                            /*if (Elem[i].vx[d] < 0)
                            {
                                Elem[i].vx[d] = 0.0;
                                Elem[i].x[d] = 0.5 * delta_x;
                            }*/
                        }
                        else
                        {
                            Elem[i].vx[d] = v_min[d][i] + dt * dvx[d][i] + av[d][i];
                            Elem[i].x[d] += dt * Elem[i].vx[d];
                        }
                    }
                    else
                    {
                        Elem[i].vx[d] = v_min[d][i] + dt * dvx[d][i] + av[d][i];
                        Elem[i].x[d] += dt * Elem[i].vx[d];
                    }
                }
                Elem[i].sxx = sxx_min[i] + (dt) * dsxx[i];
                Elem[i].syy = syy_min[i] + (dt) * dsyy[i];
                Elem[i].szz = szz_min[i] + (dt) * dszz[i];
                Elem[i].sxy = sxy_min[i] + (dt) * dsxy[i];
                Elem[i].syz = syz_min[i] + (dt) * dsyz[i];
                Elem[i].sxz = sxz_min[i] + (dt) * dsxz[i];
            }
        }
        for (size_t i = 0; i < ntotal; ++i)
        {
            if (Elem[i].itype >= 3 && Elem[i].itype <= 9)
            {
                prochnost_mises(Elem[i].itype, Elem[i].sxx, Elem[i].syy, Elem[i].szz, Elem[i].sxy, Elem[i].syz, Elem[i].sxz);
            }
            else if (Elem[i].itype == 10)
            {
                prochnost_cem_JH2(Elem[i], dt, drho[i]);
            }
            else if (Elem[i].itype == 11)
            {
                prochnost_cem_JH2(Elem[i], dt, drho[i]);
            }
        }
        if (time >= next_output_time)
        {
            output(Elem, ntotal, time, current_ts, dt);
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
        file << time << " " << Elem[moni_particle].vx[0] << " " << Elem[moni_particle].itype << std::endl;
        file.close();
    }
}
