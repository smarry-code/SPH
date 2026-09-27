#pragma once
// physics/solid_sph_pair.hpp
// Вклад одного SolidSphPrimitive (i, j) в Solid -- fiber/CUDA_PLAN.md §3.2,
// этап Г2: перевод математики interaction type Solid-SPH
// (interactions/solid_sph.hpp::compute_solid_sph_contributions) в
// SPH_HD-функцию, общую для CPU- и (в будущих частях Г2/Г4) GPU-путей.
//
// Построчная копия арифметики compute_solid_sph_contributions (порядок
// операций и типы промежуточных величин не изменены -- требование побитной
// регрессии, план §3.2), которая сама была построчной копией:
//   - неразрывность: src/density.cpp, con_density (все три density_method);
//   - скорости деформаций/поворота: src/int_force.cpp, body1 (dim 1..3);
//   - импульс+энергия (pa_sph 1/2) и осесимметричные члены:
//     src/int_force.cpp, body2;
//   - искусственная вязкость: src/art_visc.cpp;
//   - XSPH: src/av_vel.cpp.
// Физическая вязкость (visc) по-прежнему не перенесена (см. комментарий в
// начале interactions/solid_sph.hpp).
//
// SolidSphPairInput/SolidSphPairOutput -- НЕ SolidView/SolidDelta из плана
// §3.2: это чистый вход/выход одной пары (без std::vector/индексов),
// собираемый вызывающей стороной (CPU: interactions/solid_sph.hpp из
// SolidExtended; GPU в будущих частях Г2/Г4: DeviceSolidState). Поля
// SolidSphPairOutput дублируют поля interactions::SolidDelta (кроме dh,
// который Solid-SPH не заполняет, см. комментарий у SolidDelta) отдельно
// для стороны i и стороны j, поскольку арифметика оригинала пишет в di и dj
// как в независимые накопители.
#include <cstddef>
#include "param.hpp"       // dim, axis_sim
#include "physics/hd.hpp"  // SPH_HD

namespace physics
{

struct SolidSphPairInput
{
    real_t x_i[dim] = {}, x_j[dim] = {};
    real_t vx_i[dim] = {}, vx_j[dim] = {};
    real_t mass_i = 0, mass_j = 0;
    real_t rho_i = 0, rho_j = 0;
    real_t p_i = 0, p_j = 0;
    real_t c_i = 0, c_j = 0;
    real_t hsml_i = 0, hsml_j = 0;
    real_t sxx_i = 0, syy_i = 0, szz_i = 0, sxy_i = 0, syz_i = 0, sxz_i = 0;
    real_t sxx_j = 0, syy_j = 0, szz_j = 0, sxy_j = 0, syz_j = 0, sxz_j = 0;

    // Metadata primitive (build_solid_sph_search, kernel()/sph_kernel_w_dwdx).
    real_t w = 0;
    real_t dwdx[dim] = {};

    // Параметры SimulationConfig, нужные этой математике (в real_t, как в
    // compute_solid_sph_contributions -- параметры конфига там же приведены
    // из double к real_t один раз до цикла по primitives).
    int density_method = 1;
    int pa_sph = 1;
    real_t xi = 0, epsilon_h = 0;
    real_t alpha = 0, beta = 0, etq = 0;
};

struct SolidSphPairOutput
{
    real_t dv_i[dim] = {}, dv_j[dim] = {};
    real_t drho_i = 0, drho_j = 0;
    real_t du_i = 0, du_j = 0;
    real_t dexx_i = 0, deyy_i = 0, dezz_i = 0, dexy_i = 0, dexz_i = 0, deyz_i = 0;
    real_t dexx_j = 0, deyy_j = 0, dezz_j = 0, dexy_j = 0, dexz_j = 0, deyz_j = 0;
    real_t drxy_i = 0, drxz_i = 0, dryz_i = 0;
    real_t drxy_j = 0, drxz_j = 0, dryz_j = 0;
    real_t xsph_i[dim] = {}, xsph_j[dim] = {};
};

SPH_HD inline void compute_solid_sph_pair(const SolidSphPairInput& in, SolidSphPairOutput& out)
{
    const real_t* dwdx = in.dwdx;
    const real_t w = in.w;

    // --- неразрывность (density.cpp::con_density, тело body) -----------
    {
        real_t dvx_ij[dim]; // v_i - v_j, как в density.cpp
        for (std::size_t d = 0; d < dim; ++d) dvx_ij[d] = in.vx_i[d] - in.vx_j[d];
        real_t vcc = 0;
        for (std::size_t d = 0; d < dim; ++d) vcc += dvx_ij[d] * dwdx[d];

        if (in.density_method == 1)
        {
            out.drho_i = in.mass_j * vcc;
            out.drho_j = in.mass_i * vcc;
        }
        else if (in.density_method == 2)
        {
            out.drho_i = (in.rho_i / in.rho_j) * in.mass_j * vcc;
            out.drho_j = (in.rho_j / in.rho_i) * in.mass_i * vcc;
        }
        else // density_method == 3 (искусственная диффузия, Bonet-Lok/Molteni-Colagrossi)
        {
            const real_t q_i = in.mass_j * vcc;
            const real_t q_j = in.mass_i * vcc;
            const real_t h_ij = real_t(0.5) * (in.hsml_i + in.hsml_j);
            real_t rr = 0;
            real_t dx[dim];
            for (std::size_t d = 0; d < dim; ++d) { dx[d] = in.x_i[d] - in.x_j[d]; rr += dx[d] * dx[d]; }
            const real_t vi = real_t(1.0) / in.rho_i, vj = real_t(1.0) / in.rho_j;
            const real_t denom = rr + in.epsilon_h * h_ij * h_ij;
            const real_t coef_i = real_t(2.0) * (vi / vj - real_t(1.0)) / denom;
            const real_t coef_j = real_t(2.0) * (vj / vi - real_t(1.0)) / denom;
            real_t psi_i = 0, psi_j = 0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                psi_i += (coef_i * dx[d]) * dwdx[d];
                psi_j += (coef_j * (-dx[d])) * dwdx[d];
            }
            const real_t dVi = in.mass_i * vi, dVj = in.mass_j * vj;
            const real_t c0 = real_t(0.5) * (in.c_i + in.c_j);
            const real_t D_i = in.xi * h_ij * c0 * psi_i * dVj;
            const real_t D_j = in.xi * h_ij * c0 * psi_j * dVi;
            out.drho_i = q_i - in.rho_i * D_i;
            out.drho_j = q_j - in.rho_j * D_j;
        }
    }

    // --- скорости деформаций/поворота (int_force.cpp body1) ------------
    {
        real_t dvx_ji[dim]; // v_j - v_i, как в int_force.cpp body1
        for (std::size_t d = 0; d < dim; ++d) dvx_ji[d] = in.vx_j[d] - in.vx_i[d];

        if constexpr (dim == 1)
        {
            out.dexx_i = (in.mass_j * dvx_ji[0] * dwdx[0]) / in.rho_j;
            out.dexx_j = (in.mass_i * dvx_ji[0] * dwdx[0]) / in.rho_i;
        }
        else if constexpr (dim == 2)
        {
            out.dexx_i = (in.mass_j * dvx_ji[0] * dwdx[0]) / in.rho_j;
            out.deyy_i = (in.mass_j * dvx_ji[1] * dwdx[1]) / in.rho_j;
            out.dexy_i = real_t(0.5) * (((in.mass_j * dvx_ji[0] * dwdx[1]) / in.rho_j)
                            + ((in.mass_j * dvx_ji[1] * dwdx[0]) / in.rho_j));
            out.dexx_j = (in.mass_i * dvx_ji[0] * dwdx[0]) / in.rho_i;
            out.deyy_j = (in.mass_i * dvx_ji[1] * dwdx[1]) / in.rho_i;
            out.dexy_j = real_t(0.5) * (((in.mass_i * dvx_ji[0] * dwdx[1]) / in.rho_i)
                            + ((in.mass_i * dvx_ji[1] * dwdx[0]) / in.rho_i));
            out.drxy_i = real_t(0.5) * (((in.mass_j * dvx_ji[0] * dwdx[1]) / in.rho_j)
                            - ((in.mass_j * dvx_ji[1] * dwdx[0]) / in.rho_j));
            out.drxy_j = real_t(0.5) * (((in.mass_i * dvx_ji[0] * dwdx[1]) / in.rho_i)
                            - ((in.mass_i * dvx_ji[1] * dwdx[0]) / in.rho_i));
        }
        else if constexpr (dim == 3)
        {
            out.dexx_i = (in.mass_j * dvx_ji[0] * dwdx[0]) / in.rho_j;
            out.deyy_i = (in.mass_j * dvx_ji[1] * dwdx[1]) / in.rho_j;
            out.dezz_i = (in.mass_j * dvx_ji[2] * dwdx[2]) / in.rho_j;
            out.dexy_i = real_t(0.5) * (((in.mass_j * dvx_ji[0] * dwdx[1]) / in.rho_j) + ((in.mass_j * dvx_ji[1] * dwdx[0]) / in.rho_j));
            out.deyz_i = real_t(0.5) * (((in.mass_j * dvx_ji[1] * dwdx[2]) / in.rho_j) + ((in.mass_j * dvx_ji[2] * dwdx[1]) / in.rho_j));
            out.dexz_i = real_t(0.5) * (((in.mass_j * dvx_ji[0] * dwdx[2]) / in.rho_j) + ((in.mass_j * dvx_ji[2] * dwdx[0]) / in.rho_j));
            out.dexx_j = (in.mass_i * dvx_ji[0] * dwdx[0]) / in.rho_i;
            out.deyy_j = (in.mass_i * dvx_ji[1] * dwdx[1]) / in.rho_i;
            out.dezz_j = (in.mass_i * dvx_ji[2] * dwdx[2]) / in.rho_i;
            out.dexy_j = real_t(0.5) * (((in.mass_i * dvx_ji[0] * dwdx[1]) / in.rho_i) + ((in.mass_i * dvx_ji[1] * dwdx[0]) / in.rho_i));
            out.deyz_j = real_t(0.5) * (((in.mass_i * dvx_ji[1] * dwdx[2]) / in.rho_i) + ((in.mass_i * dvx_ji[2] * dwdx[1]) / in.rho_i));
            out.dexz_j = real_t(0.5) * (((in.mass_i * dvx_ji[0] * dwdx[2]) / in.rho_i) + ((in.mass_i * dvx_ji[2] * dwdx[0]) / in.rho_i));
            out.drxy_i = real_t(0.5) * (((in.mass_j * dvx_ji[0] * dwdx[1]) / in.rho_j) - ((in.mass_j * dvx_ji[1] * dwdx[0]) / in.rho_j));
            out.dryz_i = real_t(0.5) * (((in.mass_j * dvx_ji[1] * dwdx[2]) / in.rho_j) - ((in.mass_j * dvx_ji[2] * dwdx[1]) / in.rho_j));
            out.drxz_i = real_t(0.5) * (((in.mass_j * dvx_ji[0] * dwdx[2]) / in.rho_j) - ((in.mass_j * dvx_ji[2] * dwdx[0]) / in.rho_j));
            out.drxy_j = real_t(0.5) * (((in.mass_i * dvx_ji[0] * dwdx[1]) / in.rho_i) - ((in.mass_i * dvx_ji[1] * dwdx[0]) / in.rho_i));
            out.dryz_j = real_t(0.5) * (((in.mass_i * dvx_ji[1] * dwdx[2]) / in.rho_i) - ((in.mass_i * dvx_ji[2] * dwdx[1]) / in.rho_i));
            out.drxz_j = real_t(0.5) * (((in.mass_i * dvx_ji[0] * dwdx[2]) / in.rho_i) - ((in.mass_i * dvx_ji[2] * dwdx[0]) / in.rho_i));
        }
    }

    // --- импульс + энергия (int_force.cpp body2, pa_sph 1/2) -----------
    {
        real_t dvx_ji[dim];
        for (std::size_t d = 0; d < dim; ++d) dvx_ji[d] = in.vx_j[d] - in.vx_i[d];
        const real_t rhoij = real_t(1.0) / (in.rho_i * in.rho_j);
        real_t he = 0;

        if (in.pa_sph == 1)
        {
            for (std::size_t d = 0; d < dim; ++d)
            {
                real_t h = -(in.p_i + in.p_j) * dwdx[d];
                he += dvx_ji[d] * h; // he использует h ДО добавления члена напряжений (int_force.cpp)

                if (d == 0)
                {
                    h += (in.sxx_i + in.sxx_j) * dwdx[0];
                    if constexpr (dim >= 2) h += (in.sxy_i + in.sxy_j) * dwdx[1];
                    if constexpr (dim == 3) h += (in.sxz_i + in.sxz_j) * dwdx[2];
                }
                if constexpr (dim >= 2)
                {
                    if (d == 1)
                    {
                        h += (in.sxy_i + in.sxy_j) * dwdx[0] + (in.syy_i + in.syy_j) * dwdx[1];
                        if constexpr (dim == 3) h += (in.syz_i + in.syz_j) * dwdx[2];
                    }
                }
                if constexpr (dim == 3)
                {
                    if (d == 2)
                    {
                        h += (in.sxz_i + in.sxz_j) * dwdx[0] + (in.syz_i + in.syz_j) * dwdx[1]
                           + (in.szz_i + in.szz_j) * dwdx[2];
                    }
                }

                // int_force.cpp (pa_sph==1): "h *= rhoij;" применяется к
                // ПОЛНОМУ h (давление+напряжения) ПЕРЕД домножением на
                // массу -- в отличие от he, который использует h ДО этого
                // домножения (см. "he += ... * h;" выше, до if(d==0)/...).
                h *= rhoij;
                const real_t dvxdt_i = in.mass_j * h;
                const real_t dvxdt_j = in.mass_i * h;
                out.dv_i[d] += dvxdt_i;
                out.dv_j[d] += -dvxdt_j;

                if constexpr (axis_sim && dim == 2)
                {
                    if (d == 0)
                    {
                        out.dv_i[0] += rhoij * in.mass_j * ((in.sxx_j - in.szz_j) / in.x_j[0]) * w;
                        out.dv_j[0] += rhoij * in.mass_i * ((in.sxx_i - in.szz_i) / in.x_i[0]) * w;
                    }
                    else if (d == 1)
                    {
                        out.dv_i[1] += rhoij * in.mass_j * (in.sxy_j / in.x_j[0]) * w;
                        out.dv_j[1] += rhoij * in.mass_i * (in.sxy_i / in.x_i[0]) * w;
                    }
                }
            }
            he *= rhoij;
            out.du_i += in.mass_j * he;
            out.du_j += in.mass_i * he;
        }
        else // pa_sph == 2
        {
            for (std::size_t d = 0; d < dim; ++d)
            {
                real_t h = -(in.p_i / (in.rho_i * in.rho_i) + in.p_j / (in.rho_j * in.rho_j)) * dwdx[d];
                he += dvx_ji[d] * h;

                if (d == 0)
                {
                    h += (in.sxx_i / (in.rho_i * in.rho_i) + in.sxx_j / (in.rho_j * in.rho_j)) * dwdx[0];
                    if constexpr (dim >= 2)
                        h += (in.sxy_i / (in.rho_i * in.rho_i) + in.sxy_j / (in.rho_j * in.rho_j)) * dwdx[1];
                    if constexpr (dim == 3)
                        h += (in.sxz_i / (in.rho_i * in.rho_i) + in.sxz_j / (in.rho_j * in.rho_j)) * dwdx[2];
                }
                if constexpr (dim >= 2)
                {
                    if (d == 1)
                    {
                        h += (in.sxy_i / (in.rho_i * in.rho_i) + in.sxy_j / (in.rho_j * in.rho_j)) * dwdx[0]
                           + (in.syy_i / (in.rho_i * in.rho_i) + in.syy_j / (in.rho_j * in.rho_j)) * dwdx[1];
                        if constexpr (dim == 3)
                            h += (in.syz_i / (in.rho_i * in.rho_i) + in.syz_j / (in.rho_j * in.rho_j)) * dwdx[2];
                    }
                }
                if constexpr (dim == 3)
                {
                    if (d == 2)
                    {
                        h += (in.sxz_i / (in.rho_i * in.rho_i) + in.sxz_j / (in.rho_j * in.rho_j)) * dwdx[0]
                           + (in.syz_i / (in.rho_i * in.rho_i) + in.syz_j / (in.rho_j * in.rho_j)) * dwdx[1]
                           + (in.szz_i / (in.rho_i * in.rho_i) + in.szz_j / (in.rho_j * in.rho_j)) * dwdx[2];
                    }
                }

                const real_t dvxdt_i = in.mass_j * h;
                const real_t dvxdt_j = in.mass_i * h;
                out.dv_i[d] += dvxdt_i;
                out.dv_j[d] += -dvxdt_j;

                if constexpr (axis_sim && dim == 2)
                {
                    if (d == 0)
                    {
                        out.dv_i[0] += rhoij * in.mass_j * ((in.sxx_j - in.szz_j) / in.x_j[0]) * w;
                        out.dv_j[0] += rhoij * in.mass_i * ((in.sxx_i - in.szz_i) / in.x_i[0]) * w;
                    }
                    else if (d == 1)
                    {
                        out.dv_i[1] += rhoij * in.mass_j * (in.sxy_j / in.x_j[0]) * w;
                        out.dv_j[1] += rhoij * in.mass_i * (in.sxy_i / in.x_i[0]) * w;
                    }
                }
            }
            out.du_i += in.mass_j * he;
            out.du_j += in.mass_i * he;
        }
    }

    // --- искусственная вязкость (art_visc.cpp) --------------------------
    {
        const real_t mhsml = real_t(0.5) * (in.hsml_i + in.hsml_j);
        real_t dvx_iminusj[dim], dx[dim];
        real_t vr = 0, rr = 0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            dvx_iminusj[d] = in.vx_i[d] - in.vx_j[d];
            dx[d] = in.x_i[d] - in.x_j[d];
            vr += dvx_iminusj[d] * dx[d];
            rr += dx[d] * dx[d];
        }
        if (vr < 0)
        {
            const real_t muv = mhsml * vr / (rr + mhsml * mhsml * in.etq * in.etq);
            const real_t mc = real_t(0.5) * (in.c_i + in.c_j);
            const real_t mrho = real_t(0.5) * (in.rho_i + in.rho_j);
            const real_t piv = (in.beta * muv - in.alpha * mc) * muv / mrho;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const real_t h = -piv * dwdx[d];
                const real_t dvxdt_i = in.mass_j * h;
                const real_t dvxdt_j = in.mass_i * h;
                out.dv_i[d] += dvxdt_i;
                out.dv_j[d] += -dvxdt_j;
                const real_t dedt_i = in.mass_j * dvx_iminusj[d] * h;
                const real_t dedt_j = in.mass_i * dvx_iminusj[d] * h;
                out.du_i += -dedt_i; // Acc::add(dedt[i], -(dedt_i)) в art_visc.cpp
                out.du_j += -dedt_j;
            }
        }
    }

    // --- XSPH (av_vel.cpp) ----------------------------------------------
    {
        const real_t rho_sum = in.rho_i + in.rho_j;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const real_t dvx_ij = in.vx_i[d] - in.vx_j[d];
            const real_t av_i = real_t(2.0) * in.mass_j * dvx_ij / rho_sum * w;
            const real_t av_j = real_t(2.0) * in.mass_i * dvx_ij / rho_sum * w;
            out.xsph_i[d] += -av_i;
            out.xsph_j[d] += av_j;
        }
    }
}

} // namespace physics
