#pragma once
// time_integration.hpp
// Объявления функций временной интеграции (compute_cfl_timestep — шаг по Куранту,
// time_integration — главный цикл по времени). Определения — в time_integration.cpp.
#include <cstddef>
#include <vector>
#include "single_step.hpp"
#include "ppochnost.hpp"
#include "buffers.hpp"
#include "simulation_config.hpp"

// Вычисление шага по времени из условия Куранта.
double compute_cfl_timestep(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal);

// Основной цикл временной интеграции.
void time_integration(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, double& dt, double& time, double delta_x);
