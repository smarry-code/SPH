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

// Предел текучести JH-2 sigma_eq для P*, D, eps (c -- строка const_JH_2[nmat]).
double jh2_strength(const std::vector<double>& c, double P_n, double D, double eps_rate, bool clamp_low_rate = false);

double sigma_i(double &sxx, double &syy, double &sxy, double &syz, double &sxz);

void prochnost_mises(const Materials& mats, size_t nmat, double &sxx, double &syy, double &szz, double &sxy, double &syz, double &sxz);

void prochnost_cem_JH2(const Materials& mats, ParticleSystem& sys, size_t i, double dt, double drho);
