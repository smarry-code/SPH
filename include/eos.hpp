#pragma once
// eos.hpp
// Уравнения состояния (УРС) и функции давления.
// Выделены из density.hpp,
// чтобы density.hpp содержал только расчёт плотности (sum_density/con_density).
//
// УРС принимают материальные константы по ссылке (const Materials& mats) —
// чистые функции, тестируемые изолированно, без глобального
// состояния. Ранее читали векторы aa/zz/gamm/const_JH_2 и константу
// n_Me из неконстантных глобалов param.hpp.
//
//объявления здесь, определения — в eos.cpp.
#include <cmath>
#include <cstddef>
#include "materials.hpp"

void p_gas(real_t rho, real_t u, real_t& p, real_t& c);

void p_art_water(real_t rho, real_t& p, real_t& c);

void eos_mi_grunaizen(const Materials& mats, size_t nmat, real_t rho, real_t rho_0, real_t u, real_t& p, real_t& c);

void eos_ua_teta(const Materials& mats, size_t nmat, real_t rho, real_t rho_0, real_t &p, real_t &c);

void urs_cem_JH2(const Materials& mats, size_t nmat, real_t qo, real_t qo0, real_t dp, real_t &as, real_t &pv, real_t &D);
