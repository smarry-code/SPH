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
#include "eos_types.hpp"

void p_gas(double rho, double u, double& p, double& c);

void p_art_water(double rho, double& p, double& c);

void eos_mi_grunaizen(const Materials& mats, size_t nmat, double rho, double rho_0, double u, double& p, double& c);

void eos_ua_teta(const Materials& mats, size_t nmat, double rho, double rho_0, double &p, double &c);

void urs_cem_JH2(const Materials& mats, size_t nmat, double qo, double qo0, double dp, double &as, double &pv, double &D);

// Autodyn Linear: p = K·(ρ/ρ0 − 1), c = sqrt(K/ρ0).
void eos_linear(const LinearEos& eos, double rho, double rho0, double& p, double& c);

// Autodyn Shock (Ми-Грюнайзен по Us = c0 + S1·up); u — удельная внутренняя энергия.
void eos_shock(const ShockEos& eos, double rho, double rho0, double u, double& p, double& c);

// Единый УРС материалов Джонсона-Кука (itype 12-15) для ОБОИХ путей (urs() в
// int_force.cpp и solid_eos() в solid_sph.hpp), чтобы они не могли разойтись.
// Диспетчер по m.eos (Linear/Shock). rho0 -- sys.rho_0[i] частицы (как в
// остальных УРС), а не m.rho0 из таблицы. D -- повреждение частицы: при
// D >= 1 (разрушена) p = max(p, 0) -- растяжение не воспринимается (§2.6).
// Лежит в eos.cpp: eos.hpp уже зависит от materials.hpp (ADR-005: зависимость
// УРС -> Materials допустима, обратной нет), JcMaterial живёт в materials.hpp.
void eos_jc(const JcMaterial& m, double rho, double rho0, double u, double D, double& p, double& c);
