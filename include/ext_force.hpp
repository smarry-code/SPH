#pragma once
// ext_force.hpp
// Объявления функций внешних сил (отталкивание от границы, самогравитация).
// Определения — в ext_force.cpp.
#include <cstddef>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "particle_system.hpp"

void virtual_symm_2d(ParticleSystem& sys, size_t i, size_t j, std::vector<std::vector<double>>& dvxdt);

void ext_force(ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j, std::vector<std::vector<double>>& dvxdt);
