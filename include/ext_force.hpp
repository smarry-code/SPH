#pragma once
// ext_force.hpp
// Объявления функций внешних сил (отталкивание от границы, самогравитация,
// контакт «частица-частица» — part_0.md, ADR-007). Определения — в ext_force.cpp.
//
// Контракт contact_force_type (SimulationConfig): 0 - выключен (поведение
// бит-в-бит как до внедрения контакта), 1 - пинбол Белицко, 2 - Свигл,
// 3 - импульсный (см. include/contact_force.hpp).
#include <cstddef>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "buffers.hpp"
#include "contact_force.hpp"

// dvxdt -- double при любом SPH_PRECISION (include/precision.hpp): ext_force
// вне продового пути (его вызывают только tests/test_contact.cpp и
// tests/test_e2e_bc2.cpp), и сохранение типа избавляет тесты от правок.
void virtual_symm_2d(ParticleSystem& sys, size_t i, size_t j, std::vector<std::vector<double>>& dvxdt);

// ntotal — полный диапазон частиц в парах (реальные + зеркальные), по нему
//          размерены dvxdt и ведётся обнуление/проверка на NaN.
// nreal  — из них реальных; зеркальные занимают [nreal, ntotal). Контактные
//          пары строятся только среди реальных (см. is_contact_material_pair
//          в ext_force.cpp). Когда зеркальных нет, nreal == ntotal.
void ext_force(const SimulationConfig& cfg, ParticleSystem& sys, real_t dt,
               size_t ntotal, size_t nreal, size_t niac,
               const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
               std::vector<std::vector<double>>& dvxdt, WorkBuffers& wb);
