#pragma once
// single_step.hpp
// Объявления функции печати текущего времени (time_print) и оркестратора
// одного шага симуляции (single_step). Определения — в single_step.cpp.
#include <cstddef>
#include <vector>
#include "pair_graph.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "diagnostics.hpp"
#include "buffers.hpp"

// Печать текущего системного времени (используется в main.cpp).
void time_print();

// Оркестратор одного шага: поиск соседей, плотность, внутренние/внешние силы,
// обновление длины сглаживания, усреднение скорости, объединение приращений.
void single_step(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, double dt, size_t ntotal, std::vector<double>& t, std::vector<double>& tdsdt,
                 std::vector<std::vector<double>>& dx, std::vector<std::vector<double>>& dvx, std::vector<double>& du, std::vector<double>& ds,
                 std::vector<double>& drho, std::vector<std::vector<double>>& av, std::vector<double>& dsxx, std::vector<double>& dsyy,
                 std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, WorkBuffers& wb,
                 PairGraph& graph);
