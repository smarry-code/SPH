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
#include "run_guard.hpp"
#include "fiber/model.hpp"

// Вычисление шага по времени из условия Куранта.
double compute_cfl_timestep(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal);

// Основной цикл временной интеграции.
// guard — необязательный сторож времени/застоя (см. run_guard.hpp). nullptr
// (по умолчанию) — цикл идёт до t_total, как и прежде. Если сторож сработал,
// цикл прерывается досрочно, а причина остаётся в самом guard: вызывающая
// сторона обязана проверить guard->triggered() — сам по себе выход из
// time_integration об успехе не свидетельствует.
void time_integration(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, double& dt, double& time, double delta_x,
                      RunGuard* guard = nullptr);

// Перегрузка для fiber::Model (fiber/PLAN.md, этап Ф2, пункт 7). На этом
// этапе — тонкая обёртка: вызывает прежнюю логику для Solid-части
// (model.solid) и не делает ничего физического с нитями (model.fiber/
// model.topology) — их интеграция в цикл начинается с этапа Ф3. Прежняя
// сигнатура выше не изменена и продолжает работать без правок для
// test_e2e_bc2/test_e2e_bc6.
void time_integration(const SimulationConfig& cfg, fiber::Model& model, size_t ntotal, double& dt, double& time,
                      double delta_x, RunGuard* guard = nullptr);
