// int_force.cpp
void urs(std::vector<Particle> &Elem, size_t i)
{
    if (Elem[i].itype == 1) // идеальный газ
    {
        p_gas(Elem[i].rho, Elem[i].u, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 2) // вода
    {
        p_art_water(Elem[i].rho, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 3) //
    {
        eos_mi_grunaizen(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].u, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 4) // алюминий
    {
        eos_mi_grunaizen(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].u, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 5) // LiF
    {
        eos_ua_teta(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 6) // полиуретан 160
    {
        eos_ua_teta(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 7) // полиуретан 320
    {
        eos_ua_teta(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 8) // полиуретан 626
    {
        eos_ua_teta(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 9) // PMMA
    {
        eos_ua_teta(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].p, Elem[i].c);
    }
    else if (Elem[i].itype == 10) // B4C
    {
        urs_cem_JH2(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].dp, Elem[i].c, Elem[i].p, Elem[i].D);
    }
    else if (Elem[i].itype == 11) // Al2O3
    {
        urs_cem_JH2(Elem[i].itype, Elem[i].rho, Elem[i].rho_0, Elem[i].dp, Elem[i].c, Elem[i].p, Elem[i].D);
    }
}
// Основная функция int_force
void int_force(std::vector<Particle>& Elem, size_t ntotal, size_t niac, std::vector<double>& eta, std::vector<size_t>& pair_i,
               std::vector<size_t>& pair_j, std::vector<std::vector<double>>& dwdx, const std::vector<double>& w, std::vector<double>& t,
               std::vector<std::vector<double>>& dvxdt, std::vector<double>& tdsdt, std::vector<double>& dedt, std::vector<double>& dsxx, std::vector<double>& dsyy,
               std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, std::vector<double> drho)

{
    // Инициализация массивов
    std::vector<double> txx(maxn, 0.0);
    std::vector<double> tyy(maxn, 0.0);
    std::vector<double> tzz(maxn, 0.0);
    std::vector<double> txy(maxn, 0.0);
    std::vector<double> txz(maxn, 0.0);
    std::vector<double> tyz(maxn, 0.0);
    std::vector<double> vcc(maxn, 0.0);
    // Обнуление выходных массивов
    for (size_t i = 0; i < ntotal; ++i)
    {
        tdsdt[i] = 0.0;
        dedt[i] = 0.0;
        //Скорость деформаций рассчитывается из кинематических соотношений. Не эволюционный параметр
        Elem[i].exx = 0.0;
        Elem[i].exy = 0.0;
        Elem[i].exz = 0.0;
        Elem[i].eyy = 0.0;
        Elem[i].eyz = 0.0;
        Elem[i].ezz = 0.0;
        Elem[i].rxy = 0.0;
        Elem[i].rxz = 0.0;
        Elem[i].ryz = 0.0;
        for (size_t d = 0; d < dim; ++d)
        {
            dvxdt[d][i] = 0.0;
        }
    }
    // Вычисление тензора сдвиговых напряжений
    if (visc)
    {
        for (size_t k = 0; k < niac; ++k)
        {
            size_t i = pair_i[k] - 1; // Переводим в C++ индексацию
            size_t j = pair_j[k] - 1;

            std::vector<double> dvx(dim);
            for (size_t d = 0; d < dim; ++d)
            {
                dvx[d] = Elem[j].vx[d] - Elem[i].vx[d];
            }

            double hxx = 0.0, hyy = 0.0, hzz = 0.0;
            double hxy = 0.0, hxz = 0.0, hyz = 0.0;

            if (dim == 1)
            {
                hxx = 2.0 * dvx[0] * dwdx[0][k];
            }
            else if (dim == 2)
            {
                hxx = 2.0 * dvx[0] * dwdx[0][k] - dvx[1] * dwdx[1][k];
                hyy = 2.0 * dvx[1] * dwdx[1][k] - dvx[0] * dwdx[0][k];
                hxy = dvx[0] * dwdx[1][k] + dvx[1] * dwdx[0][k];
            }
            else if (dim == 3)
            {
                hxx = 2.0 * dvx[0] * dwdx[0][k] - dvx[1] * dwdx[1][k] - dvx[2] * dwdx[2][k];
                hyy = 2.0 * dvx[1] * dwdx[1][k] - dvx[0] * dwdx[0][k] - dvx[2] * dwdx[2][k];
                hzz = 2.0 * dvx[2] * dwdx[2][k] - dvx[0] * dwdx[0][k] - dvx[1] * dwdx[1][k];
                hxy = dvx[0] * dwdx[1][k] + dvx[1] * dwdx[0][k];
                hxz = dvx[0] * dwdx[2][k] + dvx[2] * dwdx[0][k];
                hyz = dvx[1] * dwdx[2][k] + dvx[2] * dwdx[1][k];
            }

            double factor = 2.0 / 3.0;
            hxx *= factor;
            hyy *= factor;
            hzz *= factor;

            if (dim == 1)
            {
                txx[i] += Elem[j].mass * hxx / Elem[j].rho;
                txx[j] += Elem[i].mass * hxx / Elem[i].rho;
            }
            else if (dim == 2)
            {
                txx[i] += Elem[j].mass * hxx / Elem[j].rho;
                txx[j] += Elem[i].mass * hxx / Elem[i].rho;
                txy[i] += Elem[j].mass * hxy / Elem[j].rho;
                txy[j] += Elem[i].mass * hxy / Elem[i].rho;
                tyy[i] += Elem[j].mass * hyy / Elem[j].rho;
                tyy[j] += Elem[i].mass * hyy / Elem[i].rho;
            }
            else if (dim == 3)
            {
                txx[i] += Elem[j].mass * hxx / Elem[j].rho;
                txx[j] += Elem[i].mass * hxx / Elem[i].rho;
                txy[i] += Elem[j].mass * hxy / Elem[j].rho;
                txy[j] += Elem[i].mass * hxy / Elem[i].rho;
                txz[i] += Elem[j].mass * hxz / Elem[j].rho;
                txz[j] += Elem[i].mass * hxz / Elem[i].rho;
                tyy[i] += Elem[j].mass * hyy / Elem[j].rho;
                tyy[j] += Elem[i].mass * hyy / Elem[i].rho;
                tyz[i] += Elem[j].mass * hyz / Elem[j].rho;
                tyz[j] += Elem[i].mass * hyz / Elem[i].rho;
                tzz[i] += Elem[j].mass * hzz / Elem[j].rho;
                tzz[j] += Elem[i].mass * hzz / Elem[i].rho;
            }

            // Расчёт дивергенции скорости
            double hvcc = 0.0;
            for (size_t d = 0; d < dim; ++d)
            {
                hvcc += dvx[d] * dwdx[d][k];
            }

            vcc[i] += Elem[j].mass * hvcc / Elem[j].rho;
            vcc[j] += Elem[i].mass * hvcc / Elem[i].rho;
        }
    }

    // Вычисление вязкой энтропии и давления
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (visc)
        {
            if (dim == 1)
            {
                tdsdt[i] = txx[i] * txx[i];
            }
            else if (dim == 2)
            {
                tdsdt[i] = txx[i] * txx[i] + 2.0 * txy[i] * txy[i] + tyy[i] * tyy[i];
            }
            else if (dim == 3)
            {
                tdsdt[i] = txx[i] * txx[i] + 2.0 * txy[i] * txy[i] + 2.0 * txz[i] * txz[i] + tyy[i] * tyy[i] + 2.0 * tyz[i] * tyz[i] + tzz[i] * tzz[i];
            }
            tdsdt[i] = 0.5 * eta[i] / Elem[i].rho * tdsdt[i];
        }

        // Вычисление давления через уравнение состояния
        //std::wcout << "Вычисление УРС частицы № " << i << std::endl;
        urs(Elem, i);
    }
    // Скорость деформаций и тензор поворота
#pragma omp parallel
    {
#pragma omp for schedule(auto)        
        for (size_t k = 0; k < niac; ++k)
        {
            size_t i = pair_i[k] - 1;
            size_t j = pair_j[k] - 1;
            std::vector<double> dvx(dim);
            for (size_t d = 0; d < dim; ++d)
            {
                dvx[d] = Elem[j].vx[d] - Elem[i].vx[d];
            }
            if (dim == 1)
            {
                double exx_i = (Elem[j].mass * dvx[0] * dwdx[0][k]) / Elem[j].rho;
                double exx_j = (Elem[i].mass * dvx[0] * dwdx[0][k]) / Elem[i].rho;
#pragma omp atomic update
                Elem[i].exx += exx_i;
#pragma omp atomic update
                Elem[j].exx += exx_j;
            }
            else if (dim == 2)
            {
                double exx_i = (Elem[j].mass * dvx[0] * dwdx[0][k]) / Elem[j].rho;
                double eyy_i = (Elem[j].mass * dvx[1] * dwdx[1][k]) / Elem[j].rho;
                double exy_i = 0.5 * (((Elem[j].mass * dvx[0] * dwdx[1][k]) / Elem[j].rho) + ((Elem[j].mass * dvx[1] * dwdx[0][k]) / Elem[j].rho));
                double exx_j = (Elem[i].mass * dvx[0] * dwdx[0][k]) / Elem[i].rho;
                double eyy_j = (Elem[i].mass * dvx[1] * dwdx[1][k]) / Elem[i].rho;
                double exy_j = 0.5 * (((Elem[i].mass * dvx[0] * dwdx[1][k]) / Elem[i].rho) + ((Elem[i].mass * dvx[1] * dwdx[0][k]) / Elem[i].rho));
                double rxy_i = 0.5 * (((Elem[j].mass * dvx[0] * dwdx[1][k]) / Elem[j].rho) - ((Elem[j].mass * dvx[1] * dwdx[0][k]) / Elem[j].rho));
                double rxy_j = 0.5 * (((Elem[i].mass * dvx[0] * dwdx[1][k]) / Elem[i].rho) - ((Elem[i].mass * dvx[1] * dwdx[0][k]) / Elem[i].rho));
#pragma omp atomic update
                Elem[i].exx += exx_i;
#pragma omp atomic update
                Elem[i].eyy += eyy_i;
#pragma omp atomic update
                Elem[i].exy += exy_i;
#pragma omp atomic update
                Elem[j].exx += exx_j;
#pragma omp atomic update
                Elem[j].eyy += eyy_j;
#pragma omp atomic update
                Elem[j].exy += exy_j;

#pragma omp atomic update
                Elem[i].rxy += rxy_i;
#pragma omp atomic update
                Elem[j].rxy += rxy_j;
            }
            else if (dim == 3)
            {
                double exx_i = (Elem[j].mass * dvx[0] * dwdx[0][k]) / Elem[j].rho;
                double eyy_i = (Elem[j].mass * dvx[1] * dwdx[1][k]) / Elem[j].rho;
                double ezz_i = (Elem[j].mass * dvx[2] * dwdx[2][k]) / Elem[j].rho;
                double exy_i = 0.5 * (((Elem[j].mass * dvx[0] * dwdx[1][k]) / Elem[j].rho) + ((Elem[j].mass * dvx[1] * dwdx[0][k]) / Elem[j].rho));
                double eyz_i = 0.5 * (((Elem[j].mass * dvx[1] * dwdx[2][k]) / Elem[j].rho) + ((Elem[j].mass * dvx[2] * dwdx[1][k]) / Elem[j].rho));
                double exz_i = 0.5 * (((Elem[j].mass * dvx[0] * dwdx[2][k]) / Elem[j].rho) + ((Elem[j].mass * dvx[2] * dwdx[0][k]) / Elem[j].rho));
                double exx_j = (Elem[i].mass * dvx[0] * dwdx[0][k]) / Elem[i].rho;
                double eyy_j = (Elem[i].mass * dvx[1] * dwdx[1][k]) / Elem[i].rho;
                double ezz_j = (Elem[i].mass * dvx[2] * dwdx[2][k]) / Elem[i].rho;
                double exy_j = 0.5 * (((Elem[i].mass * dvx[0] * dwdx[1][k]) / Elem[i].rho) + ((Elem[i].mass * dvx[1] * dwdx[0][k]) / Elem[i].rho));
                double eyz_j = 0.5 * (((Elem[i].mass * dvx[1] * dwdx[2][k]) / Elem[i].rho) + ((Elem[i].mass * dvx[2] * dwdx[1][k]) / Elem[i].rho));
                double exz_j = 0.5 * (((Elem[i].mass * dvx[0] * dwdx[2][k]) / Elem[i].rho) + ((Elem[i].mass * dvx[2] * dwdx[0][k]) / Elem[i].rho));
                double rxy_i = 0.5 * (((Elem[j].mass * dvx[0] * dwdx[1][k]) / Elem[j].rho) - ((Elem[j].mass * dvx[1] * dwdx[0][k]) / Elem[j].rho));
                double ryz_i = 0.5 * (((Elem[j].mass * dvx[1] * dwdx[2][k]) / Elem[j].rho) - ((Elem[j].mass * dvx[2] * dwdx[1][k]) / Elem[j].rho));
                double rxz_i = 0.5 * (((Elem[j].mass * dvx[0] * dwdx[2][k]) / Elem[j].rho) - ((Elem[j].mass * dvx[2] * dwdx[0][k]) / Elem[j].rho));
                double rxy_j = 0.5 * (((Elem[i].mass * dvx[0] * dwdx[1][k]) / Elem[i].rho) - ((Elem[i].mass * dvx[1] * dwdx[0][k]) / Elem[i].rho));
                double ryz_j = 0.5 * (((Elem[i].mass * dvx[1] * dwdx[2][k]) / Elem[i].rho) - ((Elem[i].mass * dvx[2] * dwdx[1][k]) / Elem[i].rho));
                double rxz_j = 0.5 * (((Elem[i].mass * dvx[0] * dwdx[2][k]) / Elem[i].rho) - ((Elem[i].mass * dvx[2] * dwdx[0][k]) / Elem[i].rho));
#pragma omp atomic update
                Elem[i].exx += exx_i;
#pragma omp atomic update
                Elem[i].eyy += eyy_i;
#pragma omp atomic update
                Elem[i].ezz += ezz_i;
#pragma omp atomic update
                Elem[i].exy += exy_i;
#pragma omp atomic update
                Elem[i].eyz += eyz_i;
#pragma omp atomic update
                Elem[i].exz += exz_i;
#pragma omp atomic update
                Elem[j].exx += exx_j;
#pragma omp atomic update
                Elem[j].eyy += eyy_j;
#pragma omp atomic update
                Elem[j].ezz += ezz_j;
#pragma omp atomic update
                Elem[j].exy += exy_j;
#pragma omp atomic update
                Elem[j].eyz += eyz_j;
#pragma omp atomic update
                Elem[j].exz += exz_j;

#pragma omp atomic update
                Elem[i].rxy += rxy_i;
#pragma omp atomic update
                Elem[i].ryz += ryz_i;
#pragma omp atomic update
                Elem[i].rxz += rxz_i;
#pragma omp atomic update
                Elem[j].rxy += rxy_j;
#pragma omp atomic update
                Elem[j].ryz += ryz_j;
#pragma omp atomic update
                Elem[j].rxz += rxz_j;
            }
        }
        // Вычисление сил давления и вязкости
#pragma omp for schedule(auto)
        for (size_t k = 0; k < niac; ++k)
        {
            size_t i = pair_i[k] - 1;
            size_t j = pair_j[k] - 1;
            std::vector<double> dvx(dim);
            for (size_t d = 0; d < dim; ++d)
            {
                dvx[d] = Elem[j].vx[d] - Elem[i].vx[d];
            }
            //---Расчет направляющих косинусов
            std::vector<double> dxiac(dim);
            std::vector<double> l(dim);
            double driac = 0.0;
            for (size_t d = 0; d < dim; ++d)
            {
                dxiac[d] = Elem[i].x[d] - Elem[j].x[d];
                driac += dxiac[d] * dxiac[d];
            }
            double r = std::sqrt(driac);
            for (size_t d = 0; d < dim; ++d)
            {
                l[d] = dxiac[d] / r;
            }
            double U_i, U_j = 0.0;
            if (dim == 1)
            {
                U_i = Elem[i].vx[0] * l[0];
                U_j = Elem[j].vx[0] * l[0];
            }
            else if (dim == 2)
            {
                U_i = Elem[i].vx[0] * l[0] + Elem[i].vx[1] * l[1];
                U_j = Elem[j].vx[0] * l[0] + Elem[j].vx[1] * l[1];
            }
            //--------------------------
            double he = 0.0;
            double rhoij = 1.0 / (Elem[i].rho * Elem[j].rho);
            // Алгоритм SPH для силы давления и вязкости
            if (pa_sph == 1)
            {
                for (size_t d = 0; d < dim; ++d)
                {
                    if (Elem[i].itype != Elem[j].itype)
                    {
                        /*double p_= (Elem[i].rho * Elem[i].c * Elem[j].p + Elem[j].rho * Elem[j].c * Elem[i].p + Elem[i].rho * Elem[i].c * Elem[j].rho * Elem[j].c * (U_i - U_j)) /
                                   (Elem[i].rho * Elem[i].c + Elem[j].rho * Elem[j].c);
                        Elem[i].p = 1.0 * p_;
                        Elem[j].p = 1.0 * p_;  */
                        double p_i = 0.5 * (Elem[i].p + Elem[j].p);
                        double p_j = 0.5 * (Elem[i].p + Elem[j].p);
#pragma omp atomic write
                        Elem[i].p = p_i;
#pragma omp atomic write
                        Elem[j].p = p_j;
                    }
                    double h = -(Elem[i].p + Elem[j].p) * dwdx[d][k];

                    he += (Elem[j].vx[d] - Elem[i].vx[d]) * h;
                    if (d == 0)
                    {
                        double h_sxx = (Elem[i].sxx + Elem[j].sxx) * dwdx[0][k];
                        h += h_sxx;
                        if (dim >= 2)
                        {
                            double h_sxy = (Elem[i].sxy + Elem[j].sxy) * dwdx[1][k];
                            h += h_sxy;
                            if (dim == 3)
                            {
                                double h_sxz = (Elem[i].sxz + Elem[j].sxz) * dwdx[2][k];
                                h += h_sxz;
                            }
                        }
                    }
                    else if (d == 1)
                    {
                        double h_sxy = (Elem[i].sxy + Elem[j].sxy) * dwdx[0][k];
                        double h_syy = (Elem[i].syy + Elem[j].syy) * dwdx[1][k];
                        h += h_sxy + h_syy;
                        if (dim == 3)
                        {
                            double h_syz = (Elem[i].syz + Elem[j].syz) * dwdx[2][k];
                            h += h_syz;
                        }
                    }
                    else if (d == 2)
                    {
                        double h_szz = (Elem[i].szz + Elem[j].szz) * dwdx[2][k];
                        double h_syz = (Elem[i].syz + Elem[j].syz) * dwdx[1][k];
                        double h_sxz = (Elem[i].sxz + Elem[j].sxz) * dwdx[0][k];
                        h += h_sxz + h_syz + h_szz;
                    }
                    // Добавление вязких сил
                    if (visc)
                    {
                        if (d == 0) // x-компонента
                        {
                            h += (eta[i] * txx[i] + eta[j] * txx[j]) * dwdx[0][k];
                            if (dim >= 2)
                            {
                                h += (eta[i] * txy[i] + eta[j] * txy[j]) * dwdx[1][k];
                                if (dim == 3)
                                {
                                    h += (eta[i] * txz[i] + eta[j] * txz[j]) * dwdx[2][k];
                                }
                            }
                        }
                        else if (d == 1) // y-компонента
                        {
                            h += (eta[i] * txy[i] + eta[j] * txy[j]) * dwdx[0][k] +
                                 (eta[i] * tyy[i] + eta[j] * tyy[j]) * dwdx[1][k];
                            if (dim == 3)
                            {
                                h += (eta[i] * tyz[i] + eta[j] * tyz[j]) * dwdx[2][k];
                            }
                        }
                        else if (d == 2) // z-компонента
                        {
                            h += (eta[i] * txz[i] + eta[j] * txz[j]) * dwdx[0][k] +
                                 (eta[i] * tyz[i] + eta[j] * tyz[j]) * dwdx[1][k] +
                                 (eta[i] * tzz[i] + eta[j] * tzz[j]) * dwdx[2][k];
                        }
                    }

                    h *= rhoij;
                    double dvxdt_i = Elem[j].mass * h;
                    double dvxdt_j = Elem[i].mass * h;
#pragma omp atomic update
                    dvxdt[d][i] += dvxdt_i;
#pragma omp atomic update
                    dvxdt[d][j] -= dvxdt_j;
                    if (axis_sim == true && dim == 2)
                    {
                        if (d == 0)
                        {
                            double sym_0_i = rhoij * Elem[j].mass * ((2 * Elem[j].sxx + Elem[j].szz) / Elem[j].x[0]) * w[k];
                            double sym_0_j = rhoij * Elem[i].mass * ((2 * Elem[i].sxx + Elem[i].szz) / Elem[i].x[0]) * w[k];
#pragma omp atomic update
                            dvxdt[d][i] += sym_0_i;
#pragma omp atomic update
                            dvxdt[d][j] += sym_0_j;
                        }
                        else if (d == 1)
                        {
                            double sym_1_i = rhoij * Elem[j].mass * (Elem[j].sxz / Elem[j].x[0]) * w[k];
                            double sym_1_j = rhoij * Elem[i].mass * (Elem[i].sxz / Elem[i].x[0]) * w[k];
#pragma omp atomic update
                            dvxdt[d][i] += sym_1_i;
#pragma omp atomic update
                            dvxdt[d][j] += sym_1_j;
                        }
                    }
                }
                he *= rhoij;
                double dedt_i = Elem[j].mass * he;
                double dedt_j = Elem[i].mass * he;
#pragma omp atomic update
                dedt[i] += dedt_i;
#pragma omp atomic update
                dedt[j] += dedt_j;
            }
            else if (pa_sph == 2)
            {
                for (size_t d = 0; d < dim; ++d)
                {
                    double h = -(Elem[i].p / std::pow(Elem[i].rho, 2) + Elem[j].p / std::pow(Elem[j].rho, 2)) * dwdx[d][k];
                    he += (Elem[j].vx[d] - Elem[i].vx[d]) * h;
                    if (d == 0)
                    {
                        double h_sxx = (Elem[i].sxx / std::pow(Elem[i].rho, 2) + Elem[j].sxx / std::pow(Elem[j].rho, 2)) * dwdx[0][k];
                        h += h_sxx;
                        if (dim >= 2)
                        {
                            double h_sxy = (Elem[i].sxy / std::pow(Elem[i].rho, 2) + Elem[j].sxy / std::pow(Elem[j].rho, 2)) * dwdx[1][k];
                            h += h_sxy;
                            if (dim == 3)
                            {
                                double h_sxz = (Elem[i].sxz / std::pow(Elem[i].rho, 2) + Elem[j].sxz / std::pow(Elem[j].rho, 2)) * dwdx[2][k];
                                h += h_sxz;
                            }
                        }
                    }
                    else if (d == 1)
                    {
                        double h_sxy = (Elem[i].sxy / std::pow(Elem[i].rho, 2) + Elem[j].sxy / std::pow(Elem[j].rho, 2)) * dwdx[0][k];
                        double h_syy = (Elem[i].syy / std::pow(Elem[i].rho, 2) + Elem[j].syy / std::pow(Elem[j].rho, 2)) * dwdx[1][k];
                        h += h_sxy + h_syy;
                        if (dim == 3)
                        {
                            double h_syz = (Elem[i].syz / std::pow(Elem[i].rho, 2) + Elem[j].syz / std::pow(Elem[j].rho, 2)) * dwdx[2][k];
                            h += h_syz;
                        }
                    }
                    else if (d == 2)
                    {
                        double h_szz = (Elem[i].szz / std::pow(Elem[i].rho, 2) + Elem[j].szz / std::pow(Elem[j].rho, 2)) * dwdx[2][k];
                        double h_syz = (Elem[i].syz / std::pow(Elem[i].rho, 2) + Elem[j].syz / std::pow(Elem[j].rho, 2)) * dwdx[1][k];
                        double h_sxz = (Elem[i].sxz / std::pow(Elem[i].rho, 2) + Elem[j].sxz / std::pow(Elem[j].rho, 2)) * dwdx[0][k];
                        h += h_sxz + h_syz + h_szz;
                    }
                    if (visc)
                    {
                        if (d == 0) // x-компонента
                        {
                            h += (eta[i] * txx[i] / std::pow(Elem[i].rho, 2) + eta[j] * txx[j] / std::pow(Elem[j].rho, 2)) * dwdx[0][k];
                            if (dim >= 2)
                            {
                                h += (eta[i] * txy[i] / std::pow(Elem[i].rho, 2) + eta[j] * txy[j] / std::pow(Elem[j].rho, 2)) * dwdx[1][k];
                                if (dim == 3)
                                {
                                    h += (eta[i] * txz[i] / std::pow(Elem[i].rho, 2) + eta[j] * txz[j] / std::pow(Elem[j].rho, 2)) * dwdx[2][k];
                                }
                            }
                        }
                        else if (d == 1) // y-компонента
                        {
                            h += (eta[i] * txy[i] / std::pow(Elem[i].rho, 2) + eta[j] * txy[j] / std::pow(Elem[j].rho, 2)) * dwdx[0][k] +
                                 (eta[i] * tyy[i] / std::pow(Elem[i].rho, 2) + eta[j] * tyy[j] / std::pow(Elem[j].rho, 2)) * dwdx[1][k];
                            if (dim == 3)
                            {
                                h += (eta[i] * tyz[i] / std::pow(Elem[i].rho, 2) + eta[j] * tyz[j] / std::pow(Elem[j].rho, 2)) * dwdx[2][k];
                            }
                        }
                        else if (d == 2) // z-компонента
                        {
                            h += (eta[i] * txz[i] / std::pow(Elem[i].rho, 2) + eta[j] * txz[j] / std::pow(Elem[j].rho, 2)) * dwdx[0][k] +
                                 (eta[i] * tyz[i] / std::pow(Elem[i].rho, 2) + eta[j] * tyz[j] / std::pow(Elem[j].rho, 2)) * dwdx[1][k] +
                                 (eta[i] * tzz[i] / std::pow(Elem[i].rho, 2) + eta[j] * tzz[j] / std::pow(Elem[j].rho, 2)) * dwdx[2][k];
                        }
                    }
                    double dvxdt_i = Elem[j].mass * h;
                    double dvxdt_j = Elem[i].mass * h;
#pragma omp atomic update
                    dvxdt[d][i] += dvxdt_i;
#pragma omp atomic update
                    dvxdt[d][j] -= dvxdt_j;
                    if (axis_sim == true && dim == 2)
                    {
                        if (d == 0)
                        {
                            double sym_0_i = rhoij * Elem[j].mass * ((2 * Elem[j].sxx + Elem[j].szz) / Elem[j].x[0]) * w[k];
                            double sym_0_j = rhoij * Elem[i].mass * ((2 * Elem[i].sxx + Elem[i].szz) / Elem[i].x[0]) * w[k];
#pragma omp atomic update
                            dvxdt[d][i] += sym_0_i;
#pragma omp atomic update
                            dvxdt[d][j] += sym_0_j;
                            // dvxdt[d][i] += ((2 * Elem[i].sxx + Elem[i].szz) /  Elem[i].x[0]);
                            // dvxdt[d][j] += ((2 * Elem[j].sxx + Elem[j].szz) /  Elem[j].x[0]);
                        }
                        else if (d == 1)
                        {
                            double sym_1_i = rhoij * Elem[j].mass * (Elem[j].sxz / Elem[j].x[0]) * w[k];
                            double sym_1_j = rhoij * Elem[i].mass * (Elem[i].sxz / Elem[i].x[0]) * w[k];
#pragma omp atomic update
                            dvxdt[d][i] += sym_1_i;
#pragma omp atomic update
                            dvxdt[d][j] += sym_1_j;
                            // dvxdt[d][i] += (Elem[i].sxz /  Elem[i].x[0]);
                            // dvxdt[d][j] += (Elem[j].sxz /  Elem[j].x[0]);
                        }
                    }
                }
                double dedt_i = Elem[j].mass * he;
                double dedt_j = Elem[i].mass * he;
#pragma omp atomic update
                dedt[i] += dedt_i;
#pragma omp atomic update
                dedt[j] += dedt_j;
            }
        }
    }

    // Обновление внутренней энергии
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (axis_sim == true && dim == 2)
        {
            dedt[i] = tdsdt[i] + 0.5 * dedt[i] + (Elem[i].sxx * (Elem[i].exx - (Elem[i].vx[0]) / (Elem[i].x[0])) + Elem[i].szz * (Elem[i].ezz - (Elem[i].vx[0]) / (Elem[i].x[0])) + 2 * Elem[i].sxy * Elem[i].exy) / Elem[i].rho;
        }
        else
        {
            dedt[i] = tdsdt[i] + 0.5 * dedt[i] + (Elem[i].sxx * Elem[i].exx + Elem[i].syy * Elem[i].eyy + Elem[i].szz * Elem[i].ezz + 2 * Elem[i].sxy * Elem[i].exy + 2 * Elem[i].syz * Elem[i].eyz + 2 * Elem[i].sxz * Elem[i].exz) / Elem[i].rho;
        }

    }
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (Elem[i].itype > 2)
        {
            size_t nmat = Elem[i].itype;
            double G_ = 0.0;
            if (nmat >= 3 && nmat <= 9)
            {
                G_ = G[nmat - n_Me];
            }            
            else if (Elem[i].itype == 10)
            {
                G_ = const_JH_2[Elem[i].itype - 10][1];
            }
            else if (Elem[i].itype == 11)
            {
                G_ = const_JH_2[Elem[i].itype - 10][1];
            }
            if (dim > 1)
            {
                //dsxx[i] = G_ * (Elem[i].exx - ((Elem[i].exx + Elem[i].eyy + Elem[i].ezz) / 3.0)) + 2 * (Elem[i].sxy * Elem[i].rxy + Elem[i].sxz * Elem[i].rxz);
                //dsyy[i] = G_ * (Elem[i].eyy - ((Elem[i].exx + Elem[i].eyy + Elem[i].ezz) / 3.0)) + 2 * (Elem[i].syz * Elem[i].ryz - Elem[i].sxy * Elem[i].rxy);
                //dszz[i] = G_ * (Elem[i].ezz - ((Elem[i].exx + Elem[i].eyy + Elem[i].ezz) / 3.0)) - 2 * (Elem[i].sxz * Elem[i].rxz - Elem[i].syz * Elem[i].ryz);
                dsxx[i] = G_ * (Elem[i].exx + (drho[i] / (3.0 * Elem[i].rho))) + 2 * (Elem[i].sxy * Elem[i].rxy + Elem[i].sxz * Elem[i].rxz);
                dsyy[i] = G_ * (Elem[i].eyy + (drho[i] / (3.0 * Elem[i].rho))) + 2 * (Elem[i].syz * Elem[i].ryz - Elem[i].sxy * Elem[i].rxy);
                dszz[i] = G_ * (Elem[i].ezz + (drho[i] / (3.0 * Elem[i].rho))) - 2 * (Elem[i].sxz * Elem[i].rxz - Elem[i].syz * Elem[i].ryz);
                dsxy[i] = G_ * Elem[i].exy + Elem[i].rxy * (Elem[i].syy - Elem[i].sxx) + Elem[i].sxz * Elem[i].ryz + Elem[i].syz * Elem[i].rxz;
                dsyz[i] = G_ * Elem[i].eyz - Elem[i].sxy * Elem[i].rxz - Elem[i].sxz * Elem[i].rxy + Elem[i].ryz * (Elem[i].szz - Elem[i].syy);
                dsxz[i] = G_ * Elem[i].exz + Elem[i].rxz * (Elem[i].szz - Elem[i].sxx) - Elem[i].sxy * Elem[i].ryz + Elem[i].syz * Elem[i].rxy;
            }
            else if (dim == 1)
            {
                //dsxx[i] = G_ * (Elem[i].exx - ((Elem[i].exx + Elem[i].eyy + Elem[i].ezz) / 3.0));
                //dsyy[i] = G_ * ( - ((Elem[i].exx + Elem[i].eyy + Elem[i].ezz) / 3.0));
                //dszz[i] = G_ * ( - ((Elem[i].exx + Elem[i].eyy + Elem[i].ezz) / 3.0));
                dsxx[i] = G_ * (Elem[i].exx - (drho[i] / (3.0 * Elem[i].rho)));
                dsyy[i] = G_ * ( - (drho[i] / (3.0 * Elem[i].rho)));
                dszz[i] = G_ * ( - (drho[i] / (3.0 * Elem[i].rho)));
            }
        }
    }
    // Проверка массива dvxdt
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            if (std::isnan(dvxdt[d][i]) || std::isinf(dvxdt[d][i]))
            {
                std::cerr << "ERROR: dvxdt[" << d << "][" << i << "] = "
                          << dvxdt[d][i] << " (NaN or Inf detected)" << std::endl;
            }
        }
    }

    // Проверка массива dedt
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (std::isnan(dedt[i]) || std::isinf(dedt[i]))
        {
                    std::cerr << "ERROR: dedt[" << i << "] = "
                              << dedt[i] << " (NaN or Inf detected)" << std::endl;
                }
            }

    // Проверка массива давления p
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (std::isnan(Elem[i].p) || std::isinf(Elem[i].p))
        {
            std::cerr << "ERROR: p[" << i << "] = "
                      << Elem[i].p << " (NaN or Inf detected)" << std::endl;
        }
    }

    // Проверка массива скорости звука c
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (std::isnan(Elem[i].c) || std::isinf(Elem[i].c))
        {
            std::cerr << "ERROR: c[" << i << "] = "
                      << Elem[i].c << " (NaN or Inf detected)" << std::endl;
        }
    }

    // Проверка массива tdsdt
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (std::isnan(tdsdt[i]) || std::isinf(tdsdt[i]))
        {
            std::cerr << "ERROR: tdsdt[" << i << "] = "
                      << tdsdt[i] << " (NaN or Inf detected)" << std::endl;
        }
    }
}
