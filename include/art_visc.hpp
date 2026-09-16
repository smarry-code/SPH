#pragma once
// art_visc.hpp
// Объявление функции искусственной вязкости.
// Определение — в art_visc.cpp.
#include <cstddef>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "pair_graph.hpp"
#include <omp.h>

void art_visc(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
              const std::vector<std::vector<double>>& dwdx, std::vector<std::vector<double>>& dvxdt, std::vector<double>& dedt, const PairGraph& graph);
