#pragma once
// h_upgrade.hpp
// Объявление функции обновления длины сглаживания.
// Определение — в h_upgrade.cpp.
#include <cstddef>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "simulation_config.hpp"
#include "particle_system.hpp"

// Обновление длины сглаживания (метод задаётся cfg.sle).
void h_upgrade(const SimulationConfig& cfg, ParticleSystem& sys, double dt, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
               const std::vector<std::vector<double>>& dwdx);
