#pragma once
// output.hpp
// Объявление функции вывода состояния частиц в файл.
// Определение — в output.cpp.
#include <cstddef>
#include "particle_system.hpp"

void output(ParticleSystem& sys, size_t ntotal, double time, size_t itimestep, double dt);
