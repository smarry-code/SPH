#pragma once
// jc_state.hpp
// Начальное состояние частиц модели Джонсона-Кука (itype 12-15).
#include <cstddef>
#include "materials.hpp"
#include "particle_system.hpp"

// T = T_r для JC-частиц, T = 0 и edot = 0 для остальных. Вызывается ОДИН раз
// до цикла по времени и НЕ вызывается при рестарте из чекпоинта (там T, edot
// восстанавливаются из файла).
inline void init_jc_state(const Materials& mats, ParticleSystem& sys, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        sys.T[i] = is_jc_itype(sys.itype[i]) ? jc_material(mats, sys.itype[i]).s.Tr : 0.0;
        sys.edot[i] = 0.0;
    }
}
