#pragma once
// tests/clamp_e2e_fiber.hpp
// ЗАГОТОВКА: сцены заполняются отдельно.
#include "clamp_case.hpp"
#include "clamp_e2e_common.hpp"

namespace clamp_e2e
{
inline void run_EF1(const Options&) { warn("EF1: не реализован"); ++e2e::failures; }
inline void run_EF2(const Options&) { warn("EF2: не реализован"); ++e2e::failures; }
inline void run_EF4(const Options&) { warn("EF4: не реализован"); ++e2e::failures; }
inline void run_EF5(const Options&) { warn("EF5: не реализован"); ++e2e::failures; }
inline void run_EF6(const Options&) { warn("EF6: не реализован"); ++e2e::failures; }
} // namespace clamp_e2e
