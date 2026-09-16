#pragma once
// av_vel.hpp
// Объявление функции усреднения скорости для предотвращения проникания.
// Определение — в av_vel.cpp.
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

void av_vel(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i,
            const std::vector<size_t>& pair_j, const std::vector<double>& w, std::vector<std::vector<double>>& av, const PairGraph& graph);
