#pragma once
// interactions/jc_strength.hpp
// Прочность и разрушение Джонсона-Кука для одной частицы (itype 12-15):
// fiber/jc-implementation-plan.md, разделы 2.1-2.3, 2.5, 2.6, этапы J4/J5.
// Вызывается из interactions::apply_material_strength (postprocessing.hpp)
// ПОСЛЕ corrector'а -- читает финальные sys.sxx.. и sys.exx.. шага. Вынесено
// в отдельный лёгкий заголовок, чтобы тесты звали её без конвейера.
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <stdexcept>

#include "materials.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "strength/johnson_cook.hpp"

// Эффективный модуль сдвига для закона Гука (§2.6): разрушенная JC-частица
// (D >= 1) не набирает девиатор заново, G = 0. Для всех остальных частиц --
// ровно shear_modulus_by_itype (то же значение, то же вычисление в Гуке).
// shear_modulus_by_itype НЕ меняется: контакт использует табличный G.
inline double effective_shear_modulus(const Materials& mats, const ParticleSystem& sys, std::size_t i)
{
    if (is_jc_itype(sys.itype[i]) && sys.D[i] >= 1.0) return 0.0;
    return shear_modulus_by_itype(mats, sys.itype[i]);
}

namespace interactions
{

// Геометрия сборки (dim/axis_sim -- constexpr из param.hpp).
constexpr jc::Geom jc_build_geom()
{
    if (dim == 1) return jc::Geom::D1;
    if (dim == 2) return axis_sim ? jc::Geom::D2Axi : jc::Geom::D2Plane;
    return jc::Geom::D3;
}

// Один шаг прочности JC для частицы i. Порядок:
//   1) T = Tr + u/cp (heating = Energy, §2.3), sys.T[i];
//   2) edot: Total -- эквивалентная скорость полной девиаторной деформации из
//      sys.exx..eyz (в осесимметрии sys.ezz = vr/r: finalize_solid пишет
//      окружную компоненту в f.ezz, single_step переносит её в sys.ezz);
//      Plastic -- Δε_p/dt прошлого шага из sys.edot;
//   3) D >= 1: девиатор = 0, ε_p не растёт, выход;
//   4) radial return по σ_eq_trial = sqrt(3 J2); масштаб на все 6 компонент;
//      sys.e_pl += Δε_p;
//   5) повреждение (§2.5): σ* = -p/Y, D += Δε_p/ε_f, D <= 1; при D = 1 --
//      девиатор сразу обнуляется;
//   6) sys.edot: Plastic -- Δε_p/dt этого шага, Total -- уже записана в п. 2.
// plastic_acc (может быть nullptr): += m*Y*Δε_p/rho -- пластическая работа
// для лога энергии (SPH_ENERGY_LOG).
inline void apply_jc_strength(const Materials& mats, ParticleSystem& sys, std::size_t i, double dt,
                              double* plastic_acc = nullptr)
{
    constexpr jc::Geom geom = jc_build_geom();
    const JcMaterial& m = jc_material(mats, sys.itype[i]);
    if (m.th.heating != jc::Heating::Energy)
        throw std::runtime_error("JC: heating = TaylorQuinney not implemented");

    // 1) Температура из внутренней энергии.
    const double T = m.s.Tr + sys.u[i] / m.th.cp;
    sys.T[i] = T;

    // 2) Скорость деформации.
    double edot;
    if (m.s.rate_source == jc::RateSource::Plastic)
        edot = sys.edot[i];
    else
    {
        edot = jc::equivalent_strain_rate(geom, sys.exx[i], sys.eyy[i], sys.ezz[i], sys.exy[i], sys.eyz[i], sys.exz[i]);
        sys.edot[i] = edot;
    }

    // 3) Разрушенная частица: прочности нет.
    auto zero_deviator = [&]()
    {
        sys.sxx[i] = 0.0; sys.syy[i] = 0.0; sys.szz[i] = 0.0;
        sys.sxy[i] = 0.0; sys.syz[i] = 0.0; sys.sxz[i] = 0.0;
    };
    if (sys.D[i] >= 1.0)
    {
        zero_deviator();
        if (m.s.rate_source == jc::RateSource::Plastic) sys.edot[i] = 0.0;
        return;
    }

    // 4) Возврат на поверхность текучести.
    const double sig_tr = std::sqrt(3.0 * jc::j2_of_deviator(geom, sys.sxx[i], sys.syy[i], sys.sxy[i], sys.syz[i], sys.sxz[i]));
    const jc::JcReturnResult r = jc::jc_radial_return(m.s, sig_tr, sys.e_pl[i], edot, T);
    if (r.scale != 1.0)
    {
        sys.sxx[i] *= r.scale; sys.syy[i] *= r.scale; sys.szz[i] *= r.scale;
        sys.sxy[i] *= r.scale; sys.syz[i] *= r.scale; sys.sxz[i] *= r.scale;
    }
    sys.e_pl[i] += r.d_eps_p;
    if (!r.converged)
    {
        static bool warned = false; // одно предупреждение на запуск
        if (!warned)
        {
            warned = true;
            std::fprintf(stderr, "JC: radial return did not converge (particle %zu), continuing\n", i);
        }
    }
    if (m.s.rate_source == jc::RateSource::Plastic) sys.edot[i] = dt > 0.0 ? r.d_eps_p / dt : 0.0;

    // 5) Накопление повреждения.
    if (jc::jc_damage_active(m.d) && r.d_eps_p > 0.0 && r.Y > 1.0e-12 * m.s.A)
    {
        const double triax = -sys.p[i] / r.Y; // p > 0 при сжатии
        const double ef = jc::jc_fracture_strain(m.d, m.s, triax, edot, T);
        sys.D[i] = std::min(sys.D[i] + r.d_eps_p / ef, 1.0);
        if (sys.D[i] >= 1.0) zero_deviator();
    }

    // 6) Лог энергии. Для JC логируется именно пластическая работа Y*Δε_p/rho
    // (а не упругая энергия, снятая возвратом, как у Мизеса): при упрочнении и
    // чувствительности к скорости они различаются, а тест I5 сверяет Y*Δε_p.
    if (plastic_acc && sys.rho[i] > 0.0) *plastic_acc += sys.mass[i] * r.Y * r.d_eps_p / sys.rho[i];
}

} // namespace interactions
