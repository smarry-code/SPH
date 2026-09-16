#pragma once
// direct_find.hpp
// Объявление функции прямого поиска соседних частиц.
// Определение — в direct_find.cpp.
#include <cstddef>
#include <vector>
#include <array>
#include <stdexcept>
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "kernel.hpp"

void direct_find(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t& niac, std::vector<size_t>& pair_i,
                 std::vector<size_t>& pair_j, std::vector<double>& w, std::vector<std::vector<double>>& dwdx, std::vector<int>& countiac);
