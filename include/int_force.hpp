#pragma once
// int_force.hpp
// Объявление urs() — диспетчера уравнения состояния (УРС) по типу материала
// (itype). Определение — в int_force.cpp.
//
// fiber/PLAN.md, этап Ф4-B2: основной силовой расчёт (int_force(), горячие
// циклы body1/body2 — старый путь SPH) удалён вместе с direct_find/
// linked_list/PairGraph/density/art_visc/av_vel/virt_part; урс()
// используется НЕ ТОЛЬКО через int_force() (interactions::solid_eos,
// include/interactions/solid_sph.hpp, — точная копия того же диспетчера для
// нового пути), но и напрямую в нескольких тестах для инициализации
// начального состояния (p, c) ПЕРЕД первым шагом расчёта:
// tests/regression_cases.hpp, tests/test_e2e_bc2.cpp, tests/test_e2e_bc6.cpp,
// tests/test_e2e_3d.cpp, tests/test_fiber_time_integration_equivalence.cpp —
// поэтому файл не осиротел целиком и не удалён, удалена только мёртвая
// часть (int_force()).
#include <cstddef>
#include "particle_system.hpp"
#include "eos.hpp"

// Диспетчеризация уравнения состояния по типу материала (itype) через switch
// с default — управляемый отказ.
void urs(const Materials& mats, ParticleSystem& sys, size_t i);
