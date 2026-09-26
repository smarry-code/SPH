#pragma once
// ppochnost.hpp
// Объявления моделей прочности (Мизес, JH-2) и интенсивности напряжений.
// Определения — в ppochnost.cpp.
#include <cstddef>
#include <cmath>
#include "particle_system.hpp"
#include "materials.hpp"

real_t sigma_i(real_t &sxx, real_t &syy, real_t &sxy, real_t &syz, real_t &sxz);

void prochnost_mises(const Materials& mats, size_t nmat, real_t &sxx, real_t &syy, real_t &szz, real_t &sxy, real_t &syz, real_t &sxz);

void prochnost_cem_JH2(const Materials& mats, ParticleSystem& sys, size_t i, real_t dt, real_t drho);
