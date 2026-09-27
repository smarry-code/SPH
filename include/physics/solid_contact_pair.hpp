#pragma once
// physics/solid_contact_pair.hpp
// Вклад одной контактной пары Solid-Solid -- fiber/CUDA_PLAN.md §3.2, этап
// Г2 (часть "Solid-Contact"): перевод математики
// interactions/solid_contact.hpp::compute_solid_contact_contributions /
// compute_solid_contact_diagnostics в SPH_HD-функции, общие для CPU- и
// (в будущих частях Г2/Г4) GPU-путей.
//
// Построчная копия арифметики compute_solid_contact_contributions (порядок
// операций и типы промежуточных величин не изменены -- требование побитной
// регрессии, план §3.2): геометрия пары (проникание, скорость сближения),
// contact::penalty_force (contact_force.hpp -- получил SPH_HD этой же
// частью Г2, см. комментарий там), раскладка силы по осям.
//
// Отличие compute_solid_contact_contributions / compute_solid_contact_diagnostics
// от общей математики Solid-SPH: ДВА разных вызывающих кода читают ЧАСТИЧНО
// разный набор результатов одной и той же геометрии (contributions -- dv,
// diagnostics -- проникание p и dt_limit, БЕЗ фильтра "res.fn <= 0", который
// diagnostics в старом коде никогда не применял -- см. interactions/
// solid_contact.hpp, compute_solid_contact_diagnostics). Чтобы не
// дублировать геометрию и не разойтись между двумя копиями, ОДНА функция
// compute_solid_contact_pair считает всё (геометрию, PairGeom/PairProps,
// penalty_force, силу по осям) и заполняет SolidContactPairOutput целиком;
// CPU-обёртки (interactions/solid_contact.hpp) сами решают, какие поля им
// нужны и какие свои проверки (fn <= 0 для contributions) применить поверх --
// ровно то же поведение, что было при раздельных копиях геометрии, потому
// что PairOutput.fn/dt_limit заполняются ДО проверки "fn <= 0" внутри этой
// функции (как и dv, зануляемое этой же веткой).
#include <cstddef>
#include <cmath>

#include "contact_force.hpp" // contact::PairGeom/PairProps/Result/penalty_force (SPH_HD)
#include "param.hpp"         // dim
#include "physics/hd.hpp"    // SPH_HD, sph_min/sph_max

namespace physics
{

// ----------------------------------------------------------------------
// Геометрические помощники ОДНОЙ частицы (interactions/solid_contact.hpp::
// solid_contact_radius / solid_contact_pair_dt_ref, ставшие тонкими
// обёртками над этими функциями -- см. комментарий там же о причине
// сохранения прежних имён и сигнатур: solid_contact_radius напрямую
// переиспользуется include/interactions/fiber_solid_contact.hpp, часть Г2
// его не трогает).
// ----------------------------------------------------------------------

// Побитная копия particle_radius (исходно src/ext_force.cpp, анонимный
// namespace; интерфейс interactions/solid_contact.hpp::solid_contact_radius):
// R = 0.5*(m/max(rho,rho_0))^(1/dim).
SPH_HD inline real_t solid_contact_particle_radius(real_t mass, real_t rho, real_t rho_0)
{
    if (!(rho > 0)) return real_t(0.0);
    const real_t rho_eff = sph_max(rho, rho_0);
    const real_t volume = mass / rho_eff;
    return real_t(0.5) * std::pow(volume, real_t(1.0) / static_cast<real_t>(dim));
}

// Побитная копия pair_dt_ref (исходно src/ext_force.cpp, анонимный
// namespace; интерфейс interactions/solid_contact.hpp::
// solid_contact_pair_dt_ref) для ОДНОЙ частицы k: CFL -- уже приведённый к
// real_t cfg.CFL (как в исходнике, где static_cast<real_t>(cfg.CFL)
// вычисляется в месте использования).
SPH_HD inline real_t solid_contact_particle_dt_ref(real_t CFL, real_t hsml, const real_t vx[dim], real_t c)
{
    real_t v2 = 0;
    for (std::size_t d = 0; d < dim; ++d) v2 += vx[d] * vx[d];
    const real_t denom = c + std::sqrt(v2);
    return denom > 0 ? CFL * hsml / denom : real_t(1e30);
}

// ----------------------------------------------------------------------
// Пара.
// ----------------------------------------------------------------------
struct SolidContactPairInput
{
    real_t x_i[dim] = {}, x_j[dim] = {};
    real_t vx_i[dim] = {}, vx_j[dim] = {};
    real_t mass_i = 0, mass_j = 0;

    // Радиусы частиц (solid_contact_particle_radius) и опорные dt пары
    // (solid_contact_particle_dt_ref) -- считаются ВЫЗЫВАЮЩЕЙ стороной по
    // частице ДО построения пары (как Ri/Rj/dt_ref в исходном kernel).
    real_t Ri = 0, Rj = 0;
    real_t c_i = 0, c_j = 0;
    real_t G_i = 0, G_j = 0; // модуль сдвига (materials.hpp::shear_modulus_by_itype)
    real_t dt_ref_i = 0, dt_ref_j = 0;

    int contact_type = 0; // contact::ForceType
    real_t penalty_scale = 0;
    real_t dt_applied = 0;
    real_t impulse_gamma = real_t(1.0);
};

struct SolidContactPairOutput
{
    bool contact = false;              // r > 0 && p > 0 (проникание есть)
    real_t p = 0;                      // проникание (0, если !contact)
    real_t fn = 0;                     // penalty_force::Result::fn (0, если !contact)
    real_t dt_limit = contact::kNoLimit;
    real_t dv_i[dim] = {}, dv_j[dim] = {}; // ускорение (0, если fn <= 0)
};

SPH_HD inline void compute_solid_contact_pair(const SolidContactPairInput& in, SolidContactPairOutput& out)
{
    out = SolidContactPairOutput{};

    real_t dx[dim];
    real_t r2 = 0;
    for (std::size_t d = 0; d < dim; ++d) { dx[d] = in.x_i[d] - in.x_j[d]; r2 += dx[d] * dx[d]; }
    const real_t r = std::sqrt(r2);
    const real_t p = (in.Ri + in.Rj) - r;
    if (!(r > 0) || !(p > 0)) return;
    out.contact = true;
    out.p = p;

    real_t dot = 0, v2 = 0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        const real_t dv = in.vx_i[d] - in.vx_j[d];
        dot += dx[d] * dv;
        v2 += dv * dv;
    }
    contact::PairGeom g;
    g.r = r; g.p = p; g.pdot = -dot / r; g.v_rel = std::sqrt(v2);

    contact::PairProps pr;
    pr.m_red = in.mass_i * in.mass_j / (in.mass_i + in.mass_j);
    pr.c_pair = sph_min(in.c_i, in.c_j);
    pr.g_eff = (in.G_i > 0 && in.G_j > 0) ? in.G_i * in.G_j / (in.G_i + in.G_j) : real_t(0.0);
    pr.r_eff = in.Ri * in.Rj / (in.Ri + in.Rj);
    pr.dt_ref = sph_min(in.dt_ref_i, in.dt_ref_j);

    const contact::Result res = contact::penalty_force(in.contact_type, g, pr, in.penalty_scale, in.dt_applied, in.impulse_gamma);
    out.fn = res.fn;
    out.dt_limit = res.dt_limit;
    if (res.fn <= 0) return;

    const real_t inv_r = real_t(1.0) / r;
    for (std::size_t d = 0; d < dim; ++d)
    {
        const real_t f = res.fn * dx[d] * inv_r;
        out.dv_i[d] = f / in.mass_i;
        out.dv_j[d] = -f / in.mass_j;
    }
}

} // namespace physics
