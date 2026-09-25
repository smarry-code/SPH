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

// Перегрузка для fiber::Model (fiber/PLAN.md, этап Ф2, пункт 7 — сигнатура;
// этап Ф3 — физика). Интегрирует Solid (model.solid) тем же
// PredictorCorrector, что и старая сигнатура, и одновременно нити
// (model.fiber/model.topology): на каждом шаге считается RHS связей
// Fiber-Bond (include/interactions/fiber_bond.hpp) и её dt_limit, итоговый
// шаг — общий минимум CFL/контакта Solid и dt_limit связей (fiber/SPH.md,
// «Шаг по времени»). Solid-цикл внутри этой перегрузки — копия цикла старой
// сигнатуры (см. TODO(Ф8) в src/time_integration.cpp); эквивалентность
// проверяет tests/test_fiber_time_integration_equivalence.cpp. Прежняя
// сигнатура выше не изменена и продолжает работать без правок для
// test_e2e_bc2/test_e2e_bc6.
void time_integration(const SimulationConfig& cfg, fiber::Model& model, size_t ntotal, double& dt, double& time,
                      double delta_x, RunGuard* guard = nullptr);
