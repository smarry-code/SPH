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

void p_gas(double rho, double u, double& p, double& c);

void p_art_water(double rho, double& p, double& c);

void eos_mi_grunaizen(const Materials& mats, size_t nmat, double rho, double rho_0, double u, double& p, double& c);

void eos_ua_teta(const Materials& mats, size_t nmat, double rho, double rho_0, double &p, double &c);

void urs_cem_JH2(const Materials& mats, size_t nmat, double qo, double qo0, double dp, double &as, double &pv, double &D);
