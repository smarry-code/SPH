#pragma once
// tests/clamp_e2e_solid.hpp
// ЗАГОТОВКА: сцены заполняются отдельно.
#include "clamp_case.hpp"
#include "clamp_e2e_common.hpp"

namespace clamp_e2e
{
inline void run_ES1(const Options&) { warn("ES1: не реализован"); ++e2e::failures; }
inline void run_ES2(const Options&) { warn("ES2: не реализован"); ++e2e::failures; }
inline void run_ES3(const Options&) { warn("ES3: не реализован"); ++e2e::failures; }
inline void run_ES4(const Options&) { warn("ES4: не реализован"); ++e2e::failures; }
inline void run_EF3(const Options&) { warn("EF3: не реализован"); ++e2e::failures; }
inline void run_EK1(const Options&) { warn("EK1: не реализован"); ++e2e::failures; }
} // namespace clamp_e2e
