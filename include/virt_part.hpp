#pragma once
// virt_part.hpp
// Объявление функции генерации виртуальных частиц.
// Определение — в virt_part.cpp.
#include <cstddef>
#include "particle_system.hpp"

// maxn — верхняя граница запаса под зеркальные частицы (cfg.maxn).
// 0 — не задана, тогда ограничением служит только размер ParticleSystem.
void virt_part(ParticleSystem& sys, size_t ntotal, size_t& nvirt, size_t maxn = 0);
