#pragma once
// linked_list.hpp
// Объявления функций поиска соседей на сетке (link_list) и замера времени
// выполнения программы (time_elapsed). Определения — в linked_list.cpp.
//
// time_elapsed отсчитывается от момента первого включения linked_list.cpp в
// сборку (единая точка отсчёта program_start, определена в .cpp).
#include <cstddef>
#include <vector>
#include <array>
#include <chrono>
#include "flat_grid.hpp"
#include "simulation_config.hpp"
#include "particle_system.hpp"
#include "kernel.hpp"

// Замер секунд, прошедших с момента старта программы (program_start в .cpp).
void time_elapsed(double& elapsed);

// Поиск соседних частиц через связный список / плотную сетку (FlatGrid).
void link_list(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, size_t ntotal, size_t& niac, std::vector<size_t>& pair_i,
               std::vector<size_t>& pair_j, std::vector<double>& w, std::vector<std::vector<double>>& dwdx);
