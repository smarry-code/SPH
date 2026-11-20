#pragma once
// energy.hpp
// Интегральные величины состояния сцены Solid + Fiber: импульс, кинетическая
// энергия, внутренняя энергия Solid, упругая энергия связей нитей.
// Общая библиотека постановок (cases/README.md): ей пользуются и тесты
// (tests/clamp_e2e_common.hpp -- через using), и постановки кейсов
// (cases/guo2019/*/guo2019_monitor.hpp). Код перенесён из
// tests/clamp_e2e_common.hpp без изменений.

#include "fiber/fiber_system.hpp"
#include "fiber/model.hpp"
#include "interactions/fiber_bond.hpp"
#include "param.hpp"
#include "particle_system.hpp"

#include <cmath>
#include <cstddef>

namespace sph_energy
{

inline void momentum_solid(const ParticleSystem& s, std::size_t n, double (&P)[3])
{
    P[0] = P[1] = P[2] = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < dim; ++d) P[d] += s.mass[i] * s.vx[d][i];
}

inline void momentum_fiber(const fiber::FiberState& f, double (&P)[3])
{
    P[0] = P[1] = P[2] = 0.0;
    for (std::size_t k = 0; k < f.size(); ++k)
        for (std::size_t d = 0; d < dim; ++d) P[d] += f.mass[k] * f.vx[d][k];
}

inline double kinetic_solid(const ParticleSystem& s, std::size_t n)
{
    double e = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < dim; ++d) e += 0.5 * s.mass[i] * s.vx[d][i] * s.vx[d][i];
    return e;
}

inline double kinetic_fiber(const fiber::FiberState& f)
{
    double e = 0.0;
    for (std::size_t k = 0; k < f.size(); ++k)
        for (std::size_t d = 0; d < dim; ++d) e += 0.5 * f.mass[k] * f.vx[d][k] * f.vx[d][k];
    return e;
}

// Внутренняя энергия Solid sum m*u (§1.2).
inline double internal_solid(const ParticleSystem& s, std::size_t n)
{
    double e = 0.0;
    for (std::size_t i = 0; i < n; ++i) e += s.mass[i] * s.u[i];
    return e;
}

// Упругая энергия связей sum 1/2 (EA/l0)(l - l0)^2 [l > l0] (связь работает только на растяжение);
// для неупругих материалов -- с вычетом eps_p и пружиной Максвелла (см. ниже).
inline double bond_energy(const fiber::Model& m)
{
    double e = 0.0;
    for (const auto& b : m.topology.bonds)
    {
        const fiber::FiberMaterial& mat = m.materials[static_cast<std::size_t>(fiber::bond_material_id(b, m.fiber))];
        double l2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double r = m.fiber.x[d][b.j] - m.fiber.x[d][b.i];
            l2 += r * r;
        }
        const double l = std::sqrt(l2);
        if (mat.inelastic())
        {
            // Упруго-пластичная связь: запасённая энергия 0.5·E·(ε-ε_p)²·A·l0 (ε > ε_p) плюс энергия пружины
            // ветви Максвелла 0.5·q²/E1·A·l0. Пластическая работа W_p и вязкая W_v в U_b НЕ входят (см. fiber_system).
            const double A = interactions::fiber_cross_section_area(mat);
            const double eps = (l - b.l0) / b.l0;
            if (eps > b.st.eps_p) e += 0.5 * mat.E * (eps - b.st.eps_p) * (eps - b.st.eps_p) * A * b.l0;
            if (mat.visco.E1 > 0.0) e += 0.5 * b.st.q * b.st.q / mat.visco.E1 * A * b.l0;
        }
        else if (l > b.l0)
        {
            const double k = mat.E * interactions::fiber_cross_section_area(mat) / b.l0;
            e += 0.5 * k * (l - b.l0) * (l - b.l0);
        }
    }
    return e;
}

} // namespace sph_energy
