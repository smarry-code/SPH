#pragma once
// density.hpp
// Объявления функций расчёта плотности (суммирование, уравнение непрерывности).
// Определения — в density.cpp.
// УРС (p_gas/p_art_water/eos_mi_grunaizen/eos_ua_teta/urs_cem_JH2) перенесены в eos.hpp.
#include <cstddef>
#include <vector>
#include <array>
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "eos.hpp"
#include "kernel.hpp"
#include "pair_graph.hpp"
#include <omp.h>

// Вычисление плотности через SPH-суммирование.
void sum_density(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j, const std::vector<double>& w);

// Вычисление изменения плотности через уравнение непрерывности.
void con_density(const SimulationConfig& cfg, const ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
                 const std::vector<std::vector<double>>& dwdx, std::vector<double>& drhodt, const PairGraph& graph);
