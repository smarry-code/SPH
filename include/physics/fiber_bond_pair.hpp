#pragma once
// physics/fiber_bond_pair.hpp
// Вклад одной связи Fiber-Bond -- fiber/CUDA_PLAN.md §3.2, этап Г2 (часть
// "Fiber-Bond"): перевод математики
// interactions/fiber_bond.hpp::compute_fiber_bond_contributions /
// compute_fiber_bond_dt_limit в SPH_HD-функции, общие для CPU- и (в будущих
// частях Г2/Г6) GPU-путей.
//
// Построчная копия арифметики (порядок операций и типы промежуточных
// величин не изменены -- требование побитной регрессии, план §3.2): сила
// связи (только растяжение, fiber/SPH.md, «Постановка проекта») и dt_limit
// связи (fiber/SPH.md, «Шаг по времени»).
//
// fiber::FiberMaterial (fiber/fiber_system.hpp) уже POD (E, nu, d, rho,
// eps_break -- пять real_t, без std::vector) -- отдельная "FiberMaterialTable"
// (упомянутая в fiber/CUDA_PLAN.md §3.2 как файл общего назначения) для ЭТИХ
// двух функций не нужна: вызывающая сторона (interactions/fiber_bond.hpp)
// просто читает нужные материалу поля (E, d) и передаёт их сюда как
// скаляры -- ровно как SolidSphPairInput получает уже готовые p_i/c_i, а не
// таблицу материалов Solid.
#include <cstddef>
#include <cmath>

#include "param.hpp"       // dim, pi
#include "physics/hd.hpp"  // SPH_HD

namespace physics
{

// Побитная копия interactions::fiber_cross_section_area: A = pi*d^2/4.
SPH_HD inline real_t fiber_cross_section_area(real_t d)
{
    return pi * d * d / real_t(4.0);
}

struct FiberBondPairInput
{
    real_t x_i[dim] = {}, x_j[dim] = {};
    real_t mass_i = 0, mass_j = 0;
    real_t l0 = 0;
    real_t E = 0; // fiber::FiberMaterial::E частицы i (материал связи -- по i, как в исходнике)
    real_t A = 0; // fiber_cross_section_area(mat.d)
};

struct FiberBondPairOutput
{
    real_t dv_i[dim] = {}, dv_j[dim] = {};
};

// Сила связи (i, j): только растяжение (l > l0). x_ij = x_j - x_i, l = |x_ij|,
// strain = l/l0 - 1, F_i<-j = E*A*strain/l * x_ij, F_j<-i = -F_i<-j,
// dv_i = F_i<-j/m_i, dv_j = F_j<-i/m_j. l_break здесь не используется --
// критерий разрыва проверяется в Postprocessing (не в этом kernel).
SPH_HD inline void compute_fiber_bond_pair(const FiberBondPairInput& in, FiberBondPairOutput& out)
{
    out = FiberBondPairOutput{};

    real_t xij[dim];
    real_t l2 = 0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        xij[d] = in.x_j[d] - in.x_i[d];
        l2 += xij[d] * xij[d];
    }
    const real_t l = std::sqrt(l2);

    if (l > in.l0 && l > 0) // только растяжение (fiber/SPH.md, "Постановка проекта")
    {
        const real_t strain = l / in.l0 - real_t(1.0);
        const real_t coeff = in.E * in.A * strain / l; // F_i<-j = coeff * xij
        for (std::size_t d = 0; d < dim; ++d)
        {
            const real_t f = coeff * xij[d];
            out.dv_i[d] = f / in.mass_i;  // +F_i<-j / m_i
            out.dv_j[d] = -f / in.mass_j; // -F_i<-j / m_j == F_j<-i / m_j
        }
    }
}

// dt_limit одной связи (fiber/SPH.md, "Шаг по времени"; fiber/PLAN.md §3.7):
//     dt_lim = 0.8 / sqrt(k_eff / m_red), k_eff = E*A/l0,
//     m_red = m_i*m_j / (m_i+m_j).
// Считается для ВСЕХ связей независимо от знака деформации (см.
// interactions/fiber_bond.hpp::compute_fiber_bond_dt_limit).
SPH_HD inline real_t fiber_bond_pair_dt_limit(real_t E, real_t A, real_t l0, real_t mass_i, real_t mass_j)
{
    const real_t k_eff = E * A / l0;
    const real_t m_red = (mass_i * mass_j) / (mass_i + mass_j);
    return real_t(0.8) / std::sqrt(k_eff / m_red);
}

} // namespace physics
