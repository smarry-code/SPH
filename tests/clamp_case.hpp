#pragma once
// tests/clamp_case.hpp
// ВСЕ числа e2e-кейсов закрепления (Clamp/CLAMP_TESTS.md): размеры в долях dx,
// материалы, скорости, длительности, допуски приёмки и их обоснования.
// Построение сцен и проверки -- tests/clamp_e2e_solid.hpp, tests/clamp_e2e_fiber.hpp;
// параметры решателя -- tests/data/clamp_initial_conditions.txt.
//
// ЕДИНИЦЫ: метры, секунды, м/с, кг/м^3, Па.

#include <cstddef>

namespace clamp_case
{

// ===========================================================================
// Solid (ES1-ES4, EK1)
// ===========================================================================
namespace solid
{
} // namespace solid

// ===========================================================================
// Fiber (EF1-EF6)
// ===========================================================================
namespace fib
{
} // namespace fib

} // namespace clamp_case
