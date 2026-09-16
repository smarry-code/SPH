#pragma once
// virt_part.hpp
// Объявление функции генерации виртуальных частиц.
// Определение — в virt_part.cpp.
#include <cstddef>
#include "particle_system.hpp"

void virt_part(ParticleSystem& sys, size_t ntotal, size_t& nvirt);
