#pragma once
// physics/material_table.hpp
// POD-таблица материалов для Solid-Contact -- fiber/CUDA_PLAN.md §3.2, этап
// Г2 (часть "Solid-Contact"): materials.hpp::Materials хранит константы в
// std::vector/std::vector<std::vector<...>>, что запрещено в SPH_HD-коде
// (план §3.2, "Правила кода"). MaterialTable -- фиксированного размера копия
// РОВНО того, что нужно compute_solid_contact_pair (physics/
// solid_contact_pair.hpp): модуль сдвига G по itype (используется в pr.g_eff
// -- ветка Герца типа 1). Упругие константы (E, nu; materials.hpp::
// elastic_by_itype) сюда не входят -- они нужны Fiber-Solid-Contact, который
// в объём этой части Г2 не входит (fiber/CUDA_PLAN.md, задание координатора:
// "Не трогай Fiber-Fiber, Fiber-Solid").
//
// build_material_table() -- НЕ SPH_HD: собирает таблицу из Materials
// (std::vector) один раз при подготовке параметров RHS (как MaterialTable в
// целом по плану §3.1: "передаётся в ядра по значению или через
// __constant__"), не вызывается из device-кода.
//
// material_shear_modulus() -- построчная копия диспетчера
// materials.hpp::shear_modulus_by_itype над этой POD-таблицей вместо
// Materials; та же логика (int_force.cpp:636-652), тот же результат бит в
// бит для тех же itype.
#include <cstddef>

#include "materials.hpp"  // Materials (сборка таблицы, не SPH_HD)
#include "physics/hd.hpp" // SPH_HD

namespace physics
{

struct MaterialTable
{
    // itype 3..9 (Materials::G) -- металлы/полимеры, модуль сдвига Мизеса.
    static constexpr std::size_t kMetalCount = 7;
    // itype 10..11 (Materials::const_JH_2[.][1]) -- керамика JH-2, столбец G.
    static constexpr std::size_t kJh2Count = 2;

    real_t G[kMetalCount] = {};
    real_t jh2_g[kJh2Count] = {};
};

inline MaterialTable build_material_table(const Materials& m)
{
    MaterialTable t;
    for (std::size_t k = 0; k < MaterialTable::kMetalCount; ++k)
        t.G[k] = (k < m.G.size()) ? m.G[k] : real_t(0.0);
    for (std::size_t k = 0; k < MaterialTable::kJh2Count; ++k)
        t.jh2_g[k] = (k < m.const_JH_2.size() && m.const_JH_2[k].size() > 1) ? m.const_JH_2[k][1] : real_t(0.0);
    return t;
}

// Побитная копия materials.hpp::shear_modulus_by_itype над MaterialTable.
SPH_HD inline real_t material_shear_modulus(const MaterialTable& t, int itype)
{
    if (itype >= 3 && itype <= 9)
    {
        const int k = itype - static_cast<int>(Materials::n_Me);
        return (k >= 0 && static_cast<std::size_t>(k) < MaterialTable::kMetalCount) ? t.G[static_cast<std::size_t>(k)]
                                                                                     : real_t(0.0);
    }
    if (itype == 10 || itype == 11)
    {
        const int k = itype - 10;
        return (k >= 0 && static_cast<std::size_t>(k) < MaterialTable::kJh2Count) ? t.jh2_g[static_cast<std::size_t>(k)]
                                                                                    : real_t(0.0);
    }
    return real_t(0.0);
}

} // namespace physics
