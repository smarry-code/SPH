// int_force.cpp
// Определения функций внутренних сил SPH (urs — диспетчер УРС, int_force —
// основной силовой расчёт). Объявления — в int_force.hpp.
#include "int_force.hpp"
#include <omp.h>
#include <fstream>
#include <cmath>

// Диспетчеризация уравнения состояния по типу материала (itype).
// цепочка if/else без else заменена на switch с default —
// управляемый отказ. Ранее при неизвестном itype УРС молчаливо не пересчитывался
// (p/c оставались прежними => некорректные силы, возможен NaN). Теперь — явное
// исключение с диагностикой. Материалы (mats) передаются по ссылке (ADR-4).
void urs(const Materials& mats, ParticleSystem& sys, size_t i)
{
    switch (sys.itype[i])
    {
        case 1: // идеальный газ
            p_gas(sys.rho[i], sys.u[i], sys.p[i], sys.c[i]);
            break;
        case 2: // вода
            p_art_water(sys.rho[i], sys.p[i], sys.c[i]);
            break;
        case 3: case 4: // тантал, алюминий
            eos_mi_grunaizen(mats, static_cast<size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.u[i], sys.p[i], sys.c[i]);
            break;
        case 5: case 6: case 7: case 8: case 9: // LiF, полиуретаны, PMMA
            eos_ua_teta(mats, static_cast<size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.p[i], sys.c[i]);
            break;
        case 10: case 11: // B4C, Al2O3
            urs_cem_JH2(mats, static_cast<size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.dp[i], sys.c[i], sys.p[i], sys.D[i]);
            break;
        default:
            throw std::runtime_error("urs(): неизвестный itype=" + std::to_string(sys.itype[i]));
    }
}
// Основная функция int_force
void int_force(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t niac,
               const std::vector<double>& eta, const std::vector<size_t>& pair_i,
               const std::vector<size_t>& pair_j, const std::vector<std::vector<double>>& dwdx,
               const std::vector<double>& w, std::vector<double>& t,
               std::vector<std::vector<double>>& dvxdt, std::vector<double>& tdsdt, std::vector<double>& dedt, std::vector<double>& dsxx, std::vector<double>& dsyy,
               std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, const std::vector<double>& drho,
               WorkBuffers& wb, const PairGraph& graph)
 
{
    // Локальные псевдонимы параметров (горячие лямбды body1/body2 и пост-циклы
    // НЕ изменяются). pa_sph, print_step ранее — глобалы param.hpp; материалы —
    // из cfg.mats (ADR-3 этапа 6.3).
    const int pa_sph = cfg.pa_sph;
    const size_t print_step = cfg.print_step;
    const Materials& mats = cfg.mats;
    // Граничные инварианты: размерности пар и буферов согласованы
    // с конфигом. Активны без -DNDEBUG (сборка -O2 без NDEBUG, как в baseline).
    assert(ntotal <= cfg.maxn);
    assert(niac <= cfg.max_interaction);
    assert(pair_i.size() >= niac && pair_j.size() >= niac);
    assert(dwdx[0].size() >= niac);
    // Инициализация массивов
    // переиспользуемые тензорные буферы из WorkBuffers (без аллокации
    // на каждом шаге). ensure() — no-op после первого шага; clear_tensor()
    // обнуляет диапазон [0, ntotal). Локальные ссылки-псевдонимы сохраняют
    // читаемость тела функции (дальнейший код не изменяется).
    wb.ensure(ntotal);
    wb.clear_tensor(ntotal);

    // Защитная проверка ключевых величин ПЕРЕД использованием: плотность и масса
    // каждой частицы должны быть конечными и положительными. Нулевая/отрицательная/
    // NaN плотность здесь немедленно даёт деление на ноль (p/rho^2, 1/rho) и NaN
    // во всех ускорениях dvxdt и в dedt — лучше упасть с понятной диагностикой.
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (!std::isfinite(sys.rho[i]) || sys.rho[i] <= 0.0 ||
            !std::isfinite(sys.mass[i]) || sys.mass[i] <= 0.0)
        {
            throw std::runtime_error(
                "int_force: некорректная частица " + std::to_string(i) +
                " (rho=" + std::to_string(sys.rho[i]) +
                ", mass=" + std::to_string(sys.mass[i]) +
                ", itype=" + std::to_string(sys.itype[i]) + ")");
        }
    }

    std::vector<double>& txx = wb.txx;
    std::vector<double>& tyy = wb.tyy;
    std::vector<double>& tzz = wb.tzz;
    std::vector<double>& txy = wb.txy;
    std::vector<double>& txz = wb.txz;
    std::vector<double>& tyz = wb.tyz;
    // Обнуление выходных массивов
    for (size_t i = 0; i < ntotal; ++i)
    {
        tdsdt[i] = 0.0;
        dedt[i] = 0.0;
        //Скорость деформаций рассчитывается из кинематических соотношений. Не эволюционный параметр
        sys.exx[i] = 0.0;
        sys.exy[i] = 0.0;
        sys.exz[i] = 0.0;
        sys.eyy[i] = 0.0;
        sys.eyz[i] = 0.0;
        sys.ezz[i] = 0.0;
        sys.rxy[i] = 0.0;
        sys.rxz[i] = 0.0;
        sys.ryz[i] = 0.0;
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
            size_t i = pair_i[k]; // 0-базовая индексация
            size_t j = pair_j[k];

            std::array<double, dim> dvx{};
            for (size_t d = 0; d < dim; ++d)
            {
                dvx[d] = sys.vx[d][j] - sys.vx[d][i];
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
                txx[i] += sys.mass[j] * hxx / sys.rho[j];
                txx[j] += sys.mass[i] * hxx / sys.rho[i];
            }
            else if (dim == 2)
            {
                txx[i] += sys.mass[j] * hxx / sys.rho[j];
                txx[j] += sys.mass[i] * hxx / sys.rho[i];
                txy[i] += sys.mass[j] * hxy / sys.rho[j];
                txy[j] += sys.mass[i] * hxy / sys.rho[i];
                tyy[i] += sys.mass[j] * hyy / sys.rho[j];
                tyy[j] += sys.mass[i] * hyy / sys.rho[i];
            }
            else if (dim == 3)
            {
                txx[i] += sys.mass[j] * hxx / sys.rho[j];
                txx[j] += sys.mass[i] * hxx / sys.rho[i];
                txy[i] += sys.mass[j] * hxy / sys.rho[j];
                txy[j] += sys.mass[i] * hxy / sys.rho[i];
                txz[i] += sys.mass[j] * hxz / sys.rho[j];
                txz[j] += sys.mass[i] * hxz / sys.rho[i];
                tyy[i] += sys.mass[j] * hyy / sys.rho[j];
                tyy[j] += sys.mass[i] * hyy / sys.rho[i];
                tyz[i] += sys.mass[j] * hyz / sys.rho[j];
                tyz[j] += sys.mass[i] * hyz / sys.rho[i];
                tzz[i] += sys.mass[j] * hzz / sys.rho[j];
                tzz[j] += sys.mass[i] * hzz / sys.rho[i];
            }

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
            tdsdt[i] = 0.5 * eta[i] / sys.rho[i] * tdsdt[i];
        }

        // Вычисление давления через уравнение состояния
        //std::wcout << "Вычисление УРС частицы № " << i << std::endl;
        urs(mats, sys, i);
    }

    // === Мост AoS -> SoA для горячих циклов (Этап 4.3.2). ====================
    // gather копирует поля частиц в плоские непрерывные массивы sys. Горячие
    // циклы body1/body2 читают/пишут через __restrict__-указатели на эти массивы:
    // одно разыменование, потоковый доступ к кэшу, доказательство неалиасинга.
    // После циклов scatter_tensor возвращает тензор скоростей деформаций/поворота
    // (единственные поля, изменяемые в цикле body1) обратно в Elem; остальные
    // поля в sys тождественны Elem и не требуют scatter.
    // Указатели на читаемые поля (const __restrict__ — доказывают неалиасинг).
    const double* __restrict__ pmass = sys.mass.data();
    const double* __restrict__ prho  = sys.rho.data();
    const double* __restrict__ pp    = sys.p.data();
    const double* __restrict__ psxx  = sys.sxx.data();
    const double* __restrict__ psyy  = sys.syy.data();
    const double* __restrict__ pszz  = sys.szz.data();
    const double* __restrict__ psxy  = sys.sxy.data();
    const double* __restrict__ psyz  = sys.syz.data();
    const double* __restrict__ psxz  = sys.sxz.data();
    const double* __restrict__ pvx[dim];
    const double* __restrict__ px[dim];
    for (size_t d = 0; d < dim; ++d)
    {
        pvx[d] = sys.vx[d].data();
        px[d] = sys.x[d].data();
    }
    // Указатели на записываемые поля тензора деформаций/поворота (body1).
    double* __restrict__ pexx = sys.exx.data();
    double* __restrict__ peyy = sys.eyy.data();
    double* __restrict__ pezz = sys.ezz.data();
    double* __restrict__ pexy = sys.exy.data();
    double* __restrict__ pexz = sys.exz.data();
    double* __restrict__ peyz = sys.eyz.data();
    double* __restrict__ prxy = sys.rxy.data();
    double* __restrict__ prxz = sys.rxz.data();
    double* __restrict__ pryz = sys.ryz.data();

    // Скорость деформаций и тензор поворота
#ifdef SPH_PROFILE
    // --- Замер времени модернизируемых параллельных циклов (baseline: atomic) ---
    static double g_force_par_time = 0.0;
    static size_t g_force_par_calls = 0;
#endif
    // Указатели SoA захватываются по значению (pvx/px — массивы, остальные —
    // скаляры), остальные переменные (dvxdt, dedt, eta, txx/wb, dwdx, w) — по
    // ссылке. Захват по значению сохраняет __restrict__-доказательство.
    auto body1 = [&]<class Acc>(size_t i, size_t j, size_t k)
    {
                    std::array<double, dim> dvx{};
                    for (size_t d = 0; d < dim; ++d)
                    {
                        dvx[d] = pvx[d][j] - pvx[d][i];
                    }
                    if (dim == 1)
                    {
                        double exx_i = (pmass[j] * dvx[0] * dwdx[0][k]) / prho[j];
                        double exx_j = (pmass[i] * dvx[0] * dwdx[0][k]) / prho[i];
                        Acc::add(pexx[i], exx_i);
                        Acc::add(pexx[j], exx_j);
                    }
                    else if (dim == 2)
                    {
                        double exx_i = (pmass[j] * dvx[0] * dwdx[0][k]) / prho[j];
                        double eyy_i = (pmass[j] * dvx[1] * dwdx[1][k]) / prho[j];
                        double exy_i = 0.5 * (((pmass[j] * dvx[0] * dwdx[1][k]) / prho[j]) + ((pmass[j] * dvx[1] * dwdx[0][k]) / prho[j]));
                        double exx_j = (pmass[i] * dvx[0] * dwdx[0][k]) / prho[i];
                        double eyy_j = (pmass[i] * dvx[1] * dwdx[1][k]) / prho[i];
                        double exy_j = 0.5 * (((pmass[i] * dvx[0] * dwdx[1][k]) / prho[i]) + ((pmass[i] * dvx[1] * dwdx[0][k]) / prho[i]));
                        double rxy_i = 0.5 * (((pmass[j] * dvx[0] * dwdx[1][k]) / prho[j]) - ((pmass[j] * dvx[1] * dwdx[0][k]) / prho[j]));
                        double rxy_j = 0.5 * (((pmass[i] * dvx[0] * dwdx[1][k]) / prho[i]) - ((pmass[i] * dvx[1] * dwdx[0][k]) / prho[i]));
                        Acc::add(pexx[i], exx_i);
                        Acc::add(peyy[i], eyy_i);
                        Acc::add(pexy[i], exy_i);
                        Acc::add(pexx[j], exx_j);
                        Acc::add(peyy[j], eyy_j);
                        Acc::add(pexy[j], exy_j);

                        Acc::add(prxy[i], rxy_i);
                        Acc::add(prxy[j], rxy_j);
                    }
                    else if (dim == 3)
                    {
                        double exx_i = (pmass[j] * dvx[0] * dwdx[0][k]) / prho[j];
                        double eyy_i = (pmass[j] * dvx[1] * dwdx[1][k]) / prho[j];
                        double ezz_i = (pmass[j] * dvx[2] * dwdx[2][k]) / prho[j];
                        double exy_i = 0.5 * (((pmass[j] * dvx[0] * dwdx[1][k]) / prho[j]) + ((pmass[j] * dvx[1] * dwdx[0][k]) / prho[j]));
                        double eyz_i = 0.5 * (((pmass[j] * dvx[1] * dwdx[2][k]) / prho[j]) + ((pmass[j] * dvx[2] * dwdx[1][k]) / prho[j]));
                        double exz_i = 0.5 * (((pmass[j] * dvx[0] * dwdx[2][k]) / prho[j]) + ((pmass[j] * dvx[2] * dwdx[0][k]) / prho[j]));
                        double exx_j = (pmass[i] * dvx[0] * dwdx[0][k]) / prho[i];
                        double eyy_j = (pmass[i] * dvx[1] * dwdx[1][k]) / prho[i];
                        double ezz_j = (pmass[i] * dvx[2] * dwdx[2][k]) / prho[i];
                        double exy_j = 0.5 * (((pmass[i] * dvx[0] * dwdx[1][k]) / prho[i]) + ((pmass[i] * dvx[1] * dwdx[0][k]) / prho[i]));
                        double eyz_j = 0.5 * (((pmass[i] * dvx[1] * dwdx[2][k]) / prho[i]) + ((pmass[i] * dvx[2] * dwdx[1][k]) / prho[i]));
                        double exz_j = 0.5 * (((pmass[i] * dvx[0] * dwdx[2][k]) / prho[i]) + ((pmass[i] * dvx[2] * dwdx[0][k]) / prho[i]));
                        double rxy_i = 0.5 * (((pmass[j] * dvx[0] * dwdx[1][k]) / prho[j]) - ((pmass[j] * dvx[1] * dwdx[0][k]) / prho[j]));
                        double ryz_i = 0.5 * (((pmass[j] * dvx[1] * dwdx[2][k]) / prho[j]) - ((pmass[j] * dvx[2] * dwdx[1][k]) / prho[j]));
                        double rxz_i = 0.5 * (((pmass[j] * dvx[0] * dwdx[2][k]) / prho[j]) - ((pmass[j] * dvx[2] * dwdx[0][k]) / prho[j]));
                        double rxy_j = 0.5 * (((pmass[i] * dvx[0] * dwdx[1][k]) / prho[i]) - ((pmass[i] * dvx[1] * dwdx[0][k]) / prho[i]));
                        double ryz_j = 0.5 * (((pmass[i] * dvx[1] * dwdx[2][k]) / prho[i]) - ((pmass[i] * dvx[2] * dwdx[1][k]) / prho[i]));
                        double rxz_j = 0.5 * (((pmass[i] * dvx[0] * dwdx[2][k]) / prho[i]) - ((pmass[i] * dvx[2] * dwdx[0][k]) / prho[i]));
                        Acc::add(pexx[i], exx_i);
                        Acc::add(peyy[i], eyy_i);
                        Acc::add(pezz[i], ezz_i);
                        Acc::add(pexy[i], exy_i);
                        Acc::add(peyz[i], eyz_i);
                        Acc::add(pexz[i], exz_i);
                        Acc::add(pexx[j], exx_j);
                        Acc::add(peyy[j], eyy_j);
                        Acc::add(pezz[j], ezz_j);
                        Acc::add(pexy[j], exy_j);
                        Acc::add(peyz[j], eyz_j);
                        Acc::add(pexz[j], exz_j);

                        Acc::add(prxy[i], rxy_i);
                        Acc::add(pryz[i], ryz_i);
                        Acc::add(prxz[i], rxz_i);
                        Acc::add(prxy[j], rxy_j);
                        Acc::add(pryz[j], ryz_j);
                        Acc::add(prxz[j], rxz_j);
                    }
        
    };
    auto body2 = [&]<class Acc>(size_t i, size_t j, size_t k)
    {
                    std::array<double, dim> dvx{};
                    for (size_t d = 0; d < dim; ++d)
                    {
                        dvx[d] = pvx[d][j] - pvx[d][i];
                    }
                    double he = 0.0;
                    double rhoij = 1.0 / (prho[i] * prho[j]);
                    // Алгоритм SPH для силы давления и вязкости
                    if (pa_sph == 1)
                    {
                        for (size_t d = 0; d < dim; ++d)
                        {
                            double h = -(pp[i] + pp[j]) * dwdx[d][k];

                            he += (pvx[d][j] - pvx[d][i]) * h;
                            if (d == 0)
                            {
                                double h_sxx = (psxx[i] + psxx[j]) * dwdx[0][k];
                                h += h_sxx;
                                if (dim >= 2)
                                {
                                    double h_sxy = (psxy[i] + psxy[j]) * dwdx[1][k];
                                    h += h_sxy;
                                    if (dim == 3)
                                    {
                                        double h_sxz = (psxz[i] + psxz[j]) * dwdx[2][k];
                                        h += h_sxz;
                                    }
                                }
                            }
                            else if (d == 1)
                            {
                                double h_sxy = (psxy[i] + psxy[j]) * dwdx[0][k];
                                double h_syy = (psyy[i] + psyy[j]) * dwdx[1][k];
                                h += h_sxy + h_syy;
                                if (dim == 3)
                                {
                                    double h_syz = (psyz[i] + psyz[j]) * dwdx[2][k];
                                    h += h_syz;
                                }
                            }
                            else if (d == 2)
                            {
                                double h_szz = (pszz[i] + pszz[j]) * dwdx[2][k];
                                double h_syz = (psyz[i] + psyz[j]) * dwdx[1][k];
                                double h_sxz = (psxz[i] + psxz[j]) * dwdx[0][k];
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
                            double dvxdt_i = pmass[j] * h;
                            double dvxdt_j = pmass[i] * h;
                            Acc::add(dvxdt[d][i], dvxdt_i);
                            Acc::add(dvxdt[d][j], (-(dvxdt_j)));
                            if (axis_sim == true && dim == 2)
                            {
                                if (d == 0)
                                {
                                    double sym_0_i = rhoij * pmass[j] * ((psxx[j] - pszz[j]) / px[0][j]) * w[k];
                                    double sym_0_j = rhoij * pmass[i] * ((psxx[i] - pszz[i]) / px[0][i]) * w[k];
                                    Acc::add(dvxdt[d][i], sym_0_i);
                                    Acc::add(dvxdt[d][j], sym_0_j);
                                }
                                else if (d == 1)
                                {
                                    double sym_1_i = rhoij * pmass[j] * (psxy[j] / px[0][j]) * w[k];
                                    double sym_1_j = rhoij * pmass[i] * (psxy[i] / px[0][i]) * w[k];
                                    Acc::add(dvxdt[d][i], sym_1_i);
                                    Acc::add(dvxdt[d][j], sym_1_j);
                                }
                            }
                        }
                        he *= rhoij;
                        double dedt_i = pmass[j] * he;
                        double dedt_j = pmass[i] * he;
                        Acc::add(dedt[i], dedt_i);
                        Acc::add(dedt[j], dedt_j);
                    }
                    else if (pa_sph == 2)
                    {
                        for (size_t d = 0; d < dim; ++d)
                        {
                            double h = -(pp[i] / (prho[i] * prho[i]) + pp[j] / (prho[j] * prho[j])) * dwdx[d][k];
                            he += (pvx[d][j] - pvx[d][i]) * h;
                            if (d == 0)
                            {
                                double h_sxx = (psxx[i] / (prho[i] * prho[i]) + psxx[j] / (prho[j] * prho[j])) * dwdx[0][k];
                                h += h_sxx;
                                if (dim >= 2)
                                {
                                    double h_sxy = (psxy[i] / (prho[i] * prho[i]) + psxy[j] / (prho[j] * prho[j])) * dwdx[1][k];
                                    h += h_sxy;
                                    if (dim == 3)
                                    {
                                        double h_sxz = (psxz[i] / (prho[i] * prho[i]) + psxz[j] / (prho[j] * prho[j])) * dwdx[2][k];
                                        h += h_sxz;
                                    }
                                }
                            }
                            else if (d == 1)
                            {
                                double h_sxy = (psxy[i] / (prho[i] * prho[i]) + psxy[j] / (prho[j] * prho[j])) * dwdx[0][k];
                                double h_syy = (psyy[i] / (prho[i] * prho[i]) + psyy[j] / (prho[j] * prho[j])) * dwdx[1][k];
                                h += h_sxy + h_syy;
                                if (dim == 3)
                                {
                                    double h_syz = (psyz[i] / (prho[i] * prho[i]) + psyz[j] / (prho[j] * prho[j])) * dwdx[2][k];
                                    h += h_syz;
                                }
                            }
                            else if (d == 2)
                            {
                                double h_szz = (pszz[i] / (prho[i] * prho[i]) + pszz[j] / (prho[j] * prho[j])) * dwdx[2][k];
                                double h_syz = (psyz[i] / (prho[i] * prho[i]) + psyz[j] / (prho[j] * prho[j])) * dwdx[1][k];
                                double h_sxz = (psxz[i] / (prho[i] * prho[i]) + psxz[j] / (prho[j] * prho[j])) * dwdx[0][k];
                                h += h_sxz + h_syz + h_szz;
                            }
                            if (visc)
                            {
                                if (d == 0) // x-компонента
                                {
                                    h += (eta[i] * txx[i] / (prho[i] * prho[i]) + eta[j] * txx[j] / (prho[j] * prho[j])) * dwdx[0][k];
                                    if (dim >= 2)
                                    {
                                        h += (eta[i] * txy[i] / (prho[i] * prho[i]) + eta[j] * txy[j] / (prho[j] * prho[j])) * dwdx[1][k];
                                        if (dim == 3)
                                        {
                                            h += (eta[i] * txz[i] / (prho[i] * prho[i]) + eta[j] * txz[j] / (prho[j] * prho[j])) * dwdx[2][k];
                                        }
                                    }
                                }
                                else if (d == 1) // y-компонента
                                {
                                    h += (eta[i] * txy[i] / (prho[i] * prho[i]) + eta[j] * txy[j] / (prho[j] * prho[j])) * dwdx[0][k] +
                                         (eta[i] * tyy[i] / (prho[i] * prho[i]) + eta[j] * tyy[j] / (prho[j] * prho[j])) * dwdx[1][k];
                                    if (dim == 3)
                                    {
                                        h += (eta[i] * tyz[i] / (prho[i] * prho[i]) + eta[j] * tyz[j] / (prho[j] * prho[j])) * dwdx[2][k];
                                    }
                                }
                                else if (d == 2) // z-компонента
                                {
                                    h += (eta[i] * txz[i] / (prho[i] * prho[i]) + eta[j] * txz[j] / (prho[j] * prho[j])) * dwdx[0][k] +
                                         (eta[i] * tyz[i] / (prho[i] * prho[i]) + eta[j] * tyz[j] / (prho[j] * prho[j])) * dwdx[1][k] +
                                         (eta[i] * tzz[i] / (prho[i] * prho[i]) + eta[j] * tzz[j] / (prho[j] * prho[j])) * dwdx[2][k];
                                }
                            }
                            double dvxdt_i = pmass[j] * h;
                            double dvxdt_j = pmass[i] * h;
                            Acc::add(dvxdt[d][i], dvxdt_i);
                            Acc::add(dvxdt[d][j], (-(dvxdt_j)));
                            if (axis_sim == true && dim == 2)
                            {
                                if (d == 0)
                                {
                                    double sym_0_i = rhoij * pmass[j] * ((psxx[j] - pszz[j]) / px[0][j]) * w[k];
                                    double sym_0_j = rhoij * pmass[i] * ((psxx[i] - pszz[i]) / px[0][i]) * w[k];
                                    Acc::add(dvxdt[d][i], sym_0_i);
                                    Acc::add(dvxdt[d][j], sym_0_j);
                                    // dvxdt[d][i] += ((2 * sys.sxx[i] + sys.szz[i]) /  sys.x[0][i]);
                                    // dvxdt[d][j] += ((2 * sys.sxx[j] + sys.szz[j]) /  sys.x[0][j]);
                                }
                                else if (d == 1)
                                {
                                    double sym_1_i = rhoij * pmass[j] * (psxy[j] / px[0][j]) * w[k];
                                    double sym_1_j = rhoij * pmass[i] * (psxy[i] / px[0][i]) * w[k];
                                    Acc::add(dvxdt[d][i], sym_1_i);
                                    Acc::add(dvxdt[d][j], sym_1_j);
                                    // dvxdt[d][i] += (sys.sxz[i] /  sys.x[0][i]);
                                    // dvxdt[d][j] += (sys.sxz[j] /  sys.x[0][j]);
                                }
                            }
                        }
                        double dedt_i = pmass[j] * he;
                        double dedt_j = pmass[i] * he;
                        Acc::add(dedt[i], dedt_i);
                        Acc::add(dedt[j], dedt_j);
                    }
        
    };
#ifdef SPH_PROFILE
    double _fpar_t0 = omp_get_wtime();
#endif
#pragma omp parallel
    {
        // ---- Цикл 1: тензор скоростей деформаций и поворота (окраска графа пар) ----
        // Пары одного цвета не разделяют частиц => конкурентной записи нет => atomic не нужен.
        // Все потоки итерируют цвета синхронно; неявный барьер #pragma omp for
        // в конце каждого цвета разделяет цвета (c+1 может делить частицу с c).
        // ---- Цикл 1 (условный путь: coloring при num_colors>0, иначе atomic) ----
        if (graph.num_colors > 0)
        {
            for (int _col = 0; _col < graph.num_colors; ++_col)
            {
#pragma omp for schedule(static)
                for (size_t _p = graph.color_offset[static_cast<size_t>(_col)]; _p < graph.color_offset[static_cast<size_t>(_col) + 1]; ++_p)
                {
                    size_t k = graph.pairs_by_color[_p];
                    body1.template operator()<AccPlain>(pair_i[k], pair_j[k], k);
                }
            }
        }
        else
        {
#pragma omp for schedule(auto)
            for (size_t k = 0; k < niac; ++k)
            {
                body1.template operator()<AccAtomic>(pair_i[k], pair_j[k], k);
            }
        }
                // ---- Цикл 2 (условный путь: coloring при num_colors>0, иначе atomic) ----
        if (graph.num_colors > 0)
        {
            for (int _col = 0; _col < graph.num_colors; ++_col)
            {
#pragma omp for schedule(static)
                for (size_t _p = graph.color_offset[static_cast<size_t>(_col)]; _p < graph.color_offset[static_cast<size_t>(_col) + 1]; ++_p)
                {
                    size_t k = graph.pairs_by_color[_p];
                    body2.template operator()<AccPlain>(pair_i[k], pair_j[k], k);
                }
            }
        }
        else
        {
#pragma omp for schedule(auto)
            for (size_t k = 0; k < niac; ++k)
            {
                body2.template operator()<AccAtomic>(pair_i[k], pair_j[k], k);
            }
        }
    }

    // === Мост SoA -> AoS: вернуть тензор деформаций/поворота в Elem. ==========
    // body1 — единственный параллельный цикл, изменяющий поля частиц (exx..eyz,
    // rxy..ryz). Остальные поля в sys тождественны Elem (после gather), поэтому
    // scatter_tensor копирует только 9 массивов (вместо ~30 у scatter).
#ifdef SPH_PROFILE
    double _fpar_t1 = omp_get_wtime();
    g_force_par_time += (_fpar_t1 - _fpar_t0);
    ++g_force_par_calls;
    // Вывод тайминга в отдельный CSV-файл (надёжнее, чем std::wcout при редиректе).
    // Столбцы: версия, № вызова, среднее время параллельных циклов (мкс/вызов).
    if (g_force_par_calls % print_step == 0)
    {
        std::ofstream _tf("force_timing.csv", std::ios::app);
        if (_tf)
        {
            _tf << (graph.num_colors > 0 ? "cycles_coloring," : "cycles_atomic,")
                << g_force_par_calls << "," << graph.num_colors << ","
                << (g_force_par_time / static_cast<double>(g_force_par_calls)) * 1e6
                << "\n";
        }
    }
#endif

    // Окружная скорость деформации в осесимметричном случае.
    // Парный цикл выше заполняет только exx = dv_r/dr и eyy = dv_z/dz: он
    // ходит по индексам координат 0 и 1, а третьей координаты в плоскости
    // меридиана нет. Компонента ε_θθ = v_r/r относится к ЧАСТИЦЕ, а не к
    // паре, поэтому задаётся здесь — так же, как осевой член неразрывности
    // в density.cpp и осевой член энергии ниже.
    //
    // Без неё ezz остаётся нулём, и это расходится по трём местам:
    //   - dszz = 2G(ezz + drho/(3rho)) недосчитывает 2G*v_r/r, из-за чего
    //     девиатор перестаёт быть бесследовым, хотя sigma_i() для
    //     осесимметрии считает третью компоненту зависимой;
    //   - JH-2 берёт ezz в eps_p_z и копит по нему повреждение D
    //     (ppochnost.cpp, prochnost_cem_JH2);
    //   - девиаторная работа в энергии теряет окружной вклад.
    if (axis_sim == true && dim == 2)
    {
        for (size_t i = 0; i < ntotal; ++i)
        {
            sys.ezz[i] = sys.vx[0][i] / sys.x[0][i];
        }
    }

    // Обновление внутренней энергии
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (axis_sim == true && dim == 2)
        {
            // Девиаторная работа S:eps. Компоненты в осесимметричном случае:
            // sxx = S_rr, syy = S_zz, szz = S_theta, sxy = S_rz; окружная
            // скорость деформации ezz = v_r/r задана выше. Гидростатическую
            // часть работы (-p*div v) добавляет time_integration, включая её
            // осевой член -p*v_r/r/rho.
            dedt[i] = tdsdt[i] + 0.5 * dedt[i] + (sys.sxx[i] * sys.exx[i] + sys.syy[i] * sys.eyy[i] + sys.szz[i] * sys.ezz[i] + 2 * sys.sxy[i] * sys.exy[i]) / sys.rho[i];
        }
        else
        {
            dedt[i] = tdsdt[i] + 0.5 * dedt[i] + (sys.sxx[i] * sys.exx[i] + sys.syy[i] * sys.eyy[i] + sys.szz[i] * sys.ezz[i] + 2 * sys.sxy[i] * sys.exy[i] + 2 * sys.syz[i] * sys.eyz[i] + 2 * sys.sxz[i] * sys.exz[i]) / sys.rho[i];
        }

    }
    // Материальные константы (ранее — глобалы param.hpp): G, const_JH_2, n_Me.
    const auto& G = mats.G;
    const auto& const_JH_2 = mats.const_JH_2;
    const size_t n_Me = mats.n_Me;
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (sys.itype[i] > 2)
        {
            size_t nmat = static_cast<size_t>(sys.itype[i]);
            double G_ = 0.0;
            if (nmat >= 3 && nmat <= 9)
            {
                G_ = G[nmat - n_Me];
            }
            else if (sys.itype[i] == 10)
            {
                G_ = const_JH_2[static_cast<size_t>(sys.itype[i] - 10)][1];
            }
            else if (sys.itype[i] == 11)
            {
                G_ = const_JH_2[static_cast<size_t>(sys.itype[i] - 10)][1];
            }
            if (dim > 1)
            {
                //dsxx[i] = G_ * (sys.exx[i] - ((sys.exx[i] + sys.eyy[i] + sys.ezz[i]) / 3.0)) + 2 * (sys.sxy[i] * sys.rxy[i] + sys.sxz[i] * sys.rxz[i]);
                //dsyy[i] = G_ * (sys.eyy[i] - ((sys.exx[i] + sys.eyy[i] + sys.ezz[i]) / 3.0)) + 2 * (sys.syz[i] * sys.ryz[i] - sys.sxy[i] * sys.rxy[i]);
                //dszz[i] = G_ * (sys.ezz[i] - ((sys.exx[i] + sys.eyy[i] + sys.ezz[i]) / 3.0)) - 2 * (sys.sxz[i] * sys.rxz[i] - sys.syz[i] * sys.ryz[i]);
                dsxx[i] = 2 * G_ * (sys.exx[i] + (drho[i] / (3.0 * sys.rho[i]))) + 2 * (sys.sxy[i] * sys.rxy[i] + sys.sxz[i] * sys.rxz[i]);
                dsyy[i] = 2 * G_ * (sys.eyy[i] + (drho[i] / (3.0 * sys.rho[i]))) + 2 * (sys.syz[i] * sys.ryz[i] - sys.sxy[i] * sys.rxy[i]);
                // Поворотные члены производной Яуманна: dS/dt = ... + (R·S − S·R),
                // R_ab = ½(∂v_a/∂x_b − ∂v_b/∂x_a). Для zz это −2(sxz·rxz + syz·ryz);
                // прежний знак «− syz·ryz» делал след поворотной части равным
                // 4·syz·ryz, и при вращении вокруг x девиатор переставал быть
                // бесследовым (tests/test_kinematics_3d.cpp, проверки C2–C4).
                // При dim = 2 ryz = 0, и правка осесимметричную ветку не меняет.
                dszz[i] = 2 * G_ * (sys.ezz[i] + (drho[i] / (3.0 * sys.rho[i]))) - 2 * (sys.sxz[i] * sys.rxz[i] + sys.syz[i] * sys.ryz[i]);
                dsxy[i] = 2 * G_ * sys.exy[i] + sys.rxy[i] * (sys.syy[i] - sys.sxx[i]) + sys.sxz[i] * sys.ryz[i] + sys.syz[i] * sys.rxz[i];
                dsyz[i] = 2 * G_ * sys.eyz[i] - sys.sxy[i] * sys.rxz[i] - sys.sxz[i] * sys.rxy[i] + sys.ryz[i] * (sys.szz[i] - sys.syy[i]);
                dsxz[i] = 2 * G_ * sys.exz[i] + sys.rxz[i] * (sys.szz[i] - sys.sxx[i]) - sys.sxy[i] * sys.ryz[i] + sys.syz[i] * sys.rxy[i];
            }
            else if (dim == 1)
            {
                //dsxx[i] = G_ * (sys.exx[i] - ((sys.exx[i] + sys.eyy[i] + sys.ezz[i]) / 3.0));
                //dsyy[i] = G_ * ( - ((sys.exx[i] + sys.eyy[i] + sys.ezz[i]) / 3.0));
                //dszz[i] = G_ * ( - ((sys.exx[i] + sys.eyy[i] + sys.ezz[i]) / 3.0));
                dsxx[i] = 2 * G_ * (sys.exx[i] + (drho[i] / (3.0 * sys.rho[i])));
                dsyy[i] = 2 * G_ * (drho[i] / (3.0 * sys.rho[i]));
                dszz[i] = 2 * G_ * (drho[i] / (3.0 * sys.rho[i]));
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
        if (std::isnan(sys.p[i]) || std::isinf(sys.p[i]))
        {
            std::cerr << "ERROR: p[" << i << "] = "
                      << sys.p[i] << " (NaN or Inf detected)" << std::endl;
        }
    }

    // Проверка массива скорости звука c
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (std::isnan(sys.c[i]) || std::isinf(sys.c[i]))
        {
            std::cerr << "ERROR: c[" << i << "] = "
                      << sys.c[i] << " (NaN or Inf detected)" << std::endl;
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
