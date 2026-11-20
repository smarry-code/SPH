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

namespace checkpoint
{
struct UserState;
}

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

// Перегрузка для fiber::Model (notes/fiber-model/fiber-plan.md, этап Ф2, пункт 7 — сигнатура;
// этап Ф3 — физика; этап Ф8 — унификация). Интегрирует Solid (model.solid) и
// нити (model.fiber/model.topology) одновременно: на каждом шаге считается
// RHS Fiber-Bond/Fiber-Fiber-Contact/Fiber-Solid-Contact
// (include/interactions/) и их dt_limit, итоговый шаг — общий минимум
// CFL/контакта Solid и dt_limit нитей (notes/fiber-model/fiber-problem-statement.md, «Шаг по времени»).
// Решатель — ode::PredictorCorrectorSolver (по умолчанию) или ode::EulerSolver
// (include/ode/ode_solver.hpp), выбор — по cfg.ode_solver.
//
// Эта перегрузка и перегрузка выше (ParticleSystem&) — тонкие обёртки над
// ОДНОЙ общей реализацией (detail::run_time_integration<Solver>,
// src/time_integration.cpp): различаются только передаваемые
// fb/topology/materials (пустые по умолчанию у старой сигнатуры -> вся
// Fiber-физика no-op). Эквивалентность прежней сигнатуры при пустых нитях
// проверяет tests/test_fiber_time_integration_equivalence.cpp.
// user_state (необязательно): состояние вызывающего в блоке USER чекпоинта (include/checkpoint.hpp,
// checkpoint::UserState): save -- при каждой записи чекпоинта, restore -- после чтения cfg.restart_file.
// nullptr -- блок пуст, поведение прежнее.
void time_integration(const SimulationConfig& cfg, fiber::Model& model, size_t ntotal, double& dt, double& time,
                      double delta_x, RunGuard* guard = nullptr, const checkpoint::UserState* user_state = nullptr);
