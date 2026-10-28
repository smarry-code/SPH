#pragma once
// physics/smoothing_length.hpp
// update_smoothing_length (sle=1) -- fiber/CUDA_PLAN.md §3.2/§4 Г2д: этот
// файл не входит в перечень §3.2 (он перечисляет hd.hpp/material_table.hpp/
// sph_kernel.hpp/eos.hpp/solid_sph_pair.hpp/solid_finalize.hpp/
// solid_contact_pair.hpp/fiber_bond_pair.hpp/fiber_contact_geom.hpp/
// fiber_solid_geom.hpp/strength.hpp/ode_update.hpp) -- добавлен отдельным
// файлом, т.к. ни один из перечисленных не подходит по смыслу (это не EOS,
// не прочность, не решатель ODE, а геометрическая величина h). Задание Г2д
// прямо требует перенести "update_smoothing_length" в physics/ (пункт 4).
//
// Перенесена ТОЛЬКО формула метода sle=1 (h = fac*(mass/rho)^(1/dim),
// src/h_upgrade.cpp) -- она единственная используется текущими конфигами
// (BC2/BC6 -> sle=1) и единственная, которую GPU обязан поддержать
// (fiber/CUDA_PLAN.md §1.3: "sle=2 ни один текущий конфиг не использует").
// sle=0 (no-op) и sle=2 (сумма по парам соседей, div v -- не per-частичная
// формула, не переносится в SPH_HD) остаются в src/h_upgrade.cpp как есть.
#include <cmath>

#include "param.hpp"      // dim
#include "physics/hd.hpp" // SPH_HD

namespace physics
{

// src/h_upgrade.cpp, метод sle=1: h = fac*(mass/rho)^(1/dim), fac = 1.3.
// rho_floor из оригинала НЕ перенесён -- он объявлен, но нигде не
// используется в формуле (мёртвая константа, комментарий "защита плотности"
// не подключён ни к одному вычислению) -- перенос дал бы функцию, которая
// делает больше, чем оригинал; поведение (в т.ч. деление на rho<=0 без
// клампа) сохранено побитно.
SPH_HD inline real_t smoothing_length_sle1(real_t mass, real_t rho)
{
    constexpr real_t fac = real_t(1.3);
    return fac * std::pow(mass / rho, real_t(1.0) / dim);
}

} // namespace physics
