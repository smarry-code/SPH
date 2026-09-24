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

// Внутренние парсеры/генераторы input.cpp (внешняя линковка; объявлены здесь,
// чтобы не дублировать сигнатуры в tests/test_fiber_model.cpp — тест body[]
// из fiber/PLAN.md, этап Ф2, проверяет их результат напрямую).
void parse_k_file(SimulationConfig& cfg, const std::string& filename, ParticleSystem& sys);
void shock_tube(SimulationConfig& cfg, ParticleSystem& sys, size_t& ntotal, double& delta_x);
