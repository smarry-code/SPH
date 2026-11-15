#pragma once
// ppochnost.hpp
// Объявления моделей прочности (Мизес, JH-2) и интенсивности напряжений.
// Определения — в ppochnost.cpp.
#include <cstddef>
#include <cmath>
#include "particle_system.hpp"
#include "materials.hpp"

// J2 девиатора напряжений с явными dim_/axis_ (sigma_i -- обёртка над ней).
double sigma_i_general(int dim_, bool axis_, double sxx, double syy, double sxy, double syz, double sxz);

double sigma_i(double &sxx, double &syy, double &sxy, double &syz, double &sxz);

void prochnost_mises(const Materials& mats, size_t nmat, double &sxx, double &syy, double &szz, double &sxy, double &syz, double &sxz);

void prochnost_cem_JH2(const Materials& mats, ParticleSystem& sys, size_t i, double dt, double drho);
