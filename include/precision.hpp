#pragma once
// precision.hpp
// Единственное место выбора вещественного типа расчёта (fiber/CUDA_PLAN.md,
// этап Г0-FP32).
//
// real_t — тип всех переменных состояния частиц (Solid и Fiber), производных
// RHS, primitives, материальных констант и арифметики interaction kernels.
// Выбирается при сборке макросом SPH_PRECISION, по образцу SPH_DIM (ADR-002:
// каждая цель компилирует свои объектники конвейера со своими макросами):
//
//     SPH_PRECISION=64 (по умолчанию) — double, эталонный расчёт проекта;
//     SPH_PRECISION=32                — float.
//
// Сменить точность всего расчёта = сменить значение макроса; править
// остальной код не нужно.
//
// Что НАМЕРЕННО остаётся double при любом SPH_PRECISION:
//   - часы расчёта: time, накопление time += dt, t_total, save_t, RunGuard.
//     Это скаляры (стоимость нулевая), а накопление ~10^3 шагов во float
//     дало бы относительную ошибку времени ~1e-4 — время кадров и строка
//     log.txt разъехались бы с fp64 без всякой физической причины;
//   - dt хранится в double, но вычисляется как минимум величин real_t
//     (float -> double без потерь) и в решатель передаётся приведённым к
//     real_t (см. include/ode/ode_solver.hpp);
//   - параметры SimulationConfig (CFL, alpha, beta, ...) — разбираются из
//     файла как double и приводятся к real_t в месте использования;
//   - ввод/вывод и регрессионный снимок (tests/regression_snapshot.hpp):
//     real_t -> double без потерь.
//
// ПРАВИЛО ДЛЯ ЛИТЕРАЛОВ. Выражение float * 2.0 вычисляется в double (2.0 —
// литерал double), поэтому в арифметике над real_t вещественные литералы
// пишутся как real_t(2.0). При SPH_PRECISION=64 это тождественно 2.0, и
// fp64-расчёт остаётся побитно прежним. Проверка: fp32-цели собираются с
// -Wdouble-promotion (CMakeLists.txt, SPH_FP32_WARNINGS) — любое пропущенное
// повышение до double в конвейере видно предупреждением.

#ifndef SPH_PRECISION
#define SPH_PRECISION 64
#endif

#if SPH_PRECISION == 64
using real_t = double;
#elif SPH_PRECISION == 32
using real_t = float;
#else
#error "SPH_PRECISION должен быть 32 или 64"
#endif

inline constexpr int precision_bits = SPH_PRECISION;
