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

void virtual_symm_2d(ParticleSystem& sys, size_t i, size_t j, std::vector<std::vector<double>>& dvxdt);

void ext_force(const SimulationConfig& cfg, ParticleSystem& sys, double dt,
               size_t ntotal, size_t niac,
               const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
               std::vector<std::vector<double>>& dvxdt, WorkBuffers& wb);
