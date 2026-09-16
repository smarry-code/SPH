#pragma once
// input.hpp
// Объявления функций ввода/инициализации данных (parse_sph_config — разбор
// конфигурации, input — точка входа инициализации частиц).
// Определения (включая shock_tube/shear_cavity/parse_k_file и вспомогательные
// функции разбора) — в input.cpp.
#include <cstddef>
#include <string>
#include "particle_system.hpp"
#include "simulation_config.hpp"

// Разбор конфигурации симуляции из Initial_Conditions.txt в cfg.
void parse_sph_config(SimulationConfig& cfg, const std::string& filename);

// Точка входа инициализации: загрузка/генерация начальной конфигурации частиц.
void input(SimulationConfig& cfg, ParticleSystem& sys, size_t& ntotal, double& dx);
