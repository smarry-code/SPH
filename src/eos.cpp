// eos.cpp
// Определения уравнений состояния (УРС) и функций давления.
// Объявления — в eos.hpp.
//
// fiber/CUDA_PLAN.md §3.2, этап Г2д: сама арифметика перенесена в
// physics/eos.hpp (SPH_HD-функции над physics::MaterialTable) -- эти функции
// остались тонкими обёртками: строят MaterialTable из Materials (дёшево --
// фиксированные размеры kMetalCount=7/kJh2Count=2, см. materials.hpp) и
// вызывают physics:: версию. Сигнатуры (принимают Materials&, не
// MaterialTable) не изменены -- их используют напрямую несколько тестов
// (int_force.hpp, шапка файла) и src/int_force.cpp::urs(), менять сигнатуру
// вне объёма Г2д.
#include "eos.hpp"
#include "physics/eos.hpp"
#include "physics/material_table.hpp"

void p_gas(real_t rho, real_t u, real_t& p, real_t& c)
{
    physics::p_gas(rho, u, p, c);
}

void p_art_water(real_t rho, real_t& p, real_t& c)
{
    physics::p_art_water(rho, p, c);
}

void eos_mi_grunaizen(const Materials& mats, size_t nmat, real_t rho, real_t rho_0, real_t u, real_t& p, real_t& c)
{
    const physics::MaterialTable t = physics::build_material_table(mats);
    physics::eos_mi_grunaizen(t, static_cast<int>(nmat), rho, rho_0, u, p, c);
}

void eos_ua_teta(const Materials& mats, size_t nmat, real_t rho, real_t rho_0, real_t &p, real_t &c)
{
    const physics::MaterialTable t = physics::build_material_table(mats);
    physics::eos_ua_teta(t, static_cast<int>(nmat), rho, rho_0, p, c);
}

void urs_cem_JH2(const Materials& mats, size_t nmat, real_t qo, real_t qo0, real_t dp, real_t &as, real_t &pv, real_t &D)
{
    const physics::MaterialTable t = physics::build_material_table(mats);
    physics::urs_cem_JH2(t, static_cast<int>(nmat), qo, qo0, dp, as, pv, D);
}
