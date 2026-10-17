#pragma once
// interactions/fiber_contact_law.hpp
// Законы нормальной силы pinball / Swegle / импульсный для контактов нитей
// (Fiber-Solid и Fiber-Fiber; fiber/plan-fiber-contact-types.md, §2;
// fiber/SPH.md, «Законы нормальной силы»).
//
// Физика силы НЕ дублируется: fiber_penalty_force -- тонкая обёртка над
// contact::penalty_force (include/contact_force.hpp), той же функцией, что у
// Solid-Solid-Contact. Здесь только то, что специфично для нитей: параметры
// шага и опорный шаг dt_ref стороны нити.
//
// ТИП 0 = ГЕРЦ (не «выключено», в отличие от contact_force_type). Ветки
// «тип 0» здесь нет: тип 0 остаётся прежним кодом в ядрах
// (fiber_solid_contact.hpp, fiber_contact.hpp), так результат при всех ключах
// = 0 побитово прежний.
//
// СВОЙСТВО ВСЕХ ТИПОВ 1-3: при скорости проникания pdot <= 0 (покой или
// разлёт) сила равна нулю (у pinball min(F1, F2) = 0 при F1 = 0). Контакт
// неупругий: отскока нет, у неподвижного контакта нет нормальной силы, а
// значит и трения на этом шаге.
//
// dt_ref НИТИ (§2 плана). У связей критерий шага осцилляторный, без CFL
// (compute_fiber_bond_dt_limit), и нити не входят в compute_cfl_timestep.
// Поэтому опорный шаг стороны нити -- тот же предел шага связи:
//   dt_bead[k] = min по связям, касающимся k, от fiber_bond_pair_dt;
//   у бусины без связей -- виртуальная связь длины d (l0 = d, m_red = m_k/2).
// Глобальный шаг уже не больше min dt_bond, поэтому dt_applied <= dt_ref, как
// у Solid-Solid. Исключение -- изолированные бусины: их виртуальная связь в
// глобальный шаг не входит.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "contact_force.hpp"
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp" // fiber_bond_pair_dt
#include "rhs/cpu_executor.hpp"

struct SimulationConfig;

namespace interactions
{

// Закон нормальной силы. type: 0 Герц, 1 pinball, 2 Swegle, 3 импульсный,
// 4 Герц + демпфер (fiber_hertz_damped_force).
struct FiberContactLaw
{
    int type = 0;       // 0 = Герц (не «выключено»)
    double kp = 1.0;    // Kp (типы 1, 2; у типа 4 -- масштаб демпфера)
    double gamma = 1.0; // доля гашения (тип 3)
};

struct FiberContactStepParams
{
    FiberContactLaw law;
    // Шаг ПРЕДЫДУЩЕГО цикла (тип 3). Новый шаг в момент RHS ещё не известен:
    // dt пересчитывается после RHS (src/time_integration.cpp), так же как у
    // Solid-Solid (single_step.cpp передаёт в RHS шаг прошлого цикла). Если
    // шаг уменьшился, гасится меньше gamma; если вырос -- больше.
    double dt_prev = 0.0;
    // Для solid_contact_pair_dt_ref (CFL стороны Solid, только Fiber-Solid).
    const SimulationConfig* cfg = nullptr;
    // dt_bead (compute_fiber_bead_dt_ref), обязателен при type != 0.
    const std::vector<double>* bead_dt_ref = nullptr;
    // Слот предела шага на primitive (выход ядра при type != 0). Если задан,
    // ядро пишет туда res.dt_limit, а *_dt_limit берёт минимум по слотам --
    // сила и предел из одного расчёта. Если не задан, *_dt_limit пересчитывает
    // ту же функцию по primitive.
    std::vector<double>* pair_dt_limit = nullptr;
    // Только для типа 4 и лога энергии: ядро считает лишь вязкую часть силы.
    bool damper_only = false;
};

// Есть ли у закона новые (не Герц) типы -- нужен ли dt_bead.
inline bool fiber_contact_law_is_penalty(const FiberContactLaw& law) { return law.type != 0; }

// ----------------------------------------------------------------------
// compute_fiber_bead_dt_ref -- dt_bead (см. шапку файла). Один раз на шаг,
// O(nbonds), последовательно: минимум не зависит от порядка, результат
// детерминирован при любом числе потоков.
// ----------------------------------------------------------------------
inline std::vector<double>& compute_fiber_bead_dt_ref(const fiber::FiberState& state,
                                                       const fiber::FiberTopology& topology,
                                                       const std::vector<fiber::FiberMaterial>& materials,
                                                       std::vector<double>& out)
{
    const double kInf = std::numeric_limits<double>::max();
    const std::size_t n = state.size();
    out.assign(n, kInf);
    for (const fiber::FiberBond& b : topology.bonds)
    {
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(state.material_id[b.i])];
        const double dt = fiber_bond_pair_dt(mat, b.l0, state.mass[b.i], state.mass[b.j]);
        out[b.i] = std::min(out[b.i], dt);
        out[b.j] = std::min(out[b.j], dt);
    }
    for (std::size_t k = 0; k < n; ++k)
    {
        if (out[k] < kInf) continue;
        // Виртуальная связь длины d: E*A/d, m_red = m_k*m_k/(2 m_k) = m_k/2.
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(state.material_id[k])];
        out[k] = fiber_bond_pair_dt(mat, mat.d, state.mass[k], state.mass[k]);
    }
    return out;
}

// dt_ref стороны нити: точка на сегменте (k0, k1) или бусина k0 (k1 = -1).
inline double fiber_side_dt_ref(const std::vector<double>& bead_dt_ref, std::size_t k0, long long k1)
{
    const double d0 = bead_dt_ref[k0];
    if (k1 < 0) return d0;
    return std::min(d0, bead_dt_ref[static_cast<std::size_t>(k1)]);
}

// Модуль сдвига материала нити из E, nu (изотропный).
inline double fiber_shear_modulus(const fiber::FiberMaterial& m)
{
    return m.E / (2.0 * (1.0 + m.nu));
}

// Стержневая скорость звука материала нити sqrt(E/rho).
inline double fiber_sound_speed(const fiber::FiberMaterial& m)
{
    return (m.E > 0.0 && m.rho > 0.0) ? std::sqrt(m.E / m.rho) : 0.0;
}

// G_1 G_2/(G_1+G_2); 0, если хотя бы один модуль не положителен (ветка F2
// pinball тогда отключается, как у Solid-Solid).
inline double fiber_contact_g_eff(double g1, double g2)
{
    if (!(g1 > 0.0) || !(g2 > 0.0)) return 0.0;
    return g1 * g2 / (g1 + g2);
}

// = contact::penalty_force(type, g, pr, kp, dt_prev, gamma). Только для type 1-3.
inline contact::Result fiber_penalty_force(const FiberContactStepParams& params, const contact::PairGeom& g,
                                           const contact::PairProps& pr)
{
    return contact::penalty_force(params.law.type, g, pr, params.law.kp, params.dt_prev, params.law.gamma);
}

// ----------------------------------------------------------------------
// Тип 4 -- «Герц + демпфер» (модель Кельвина -- Фойгта / Ханта -- Кроссли):
//     F = F_H(delta) + c * max(pdot, 0),   c = Kp * m_red / dt_ref,
// F_H -- сила Герца той же геометрии (geom.F, как у типа 0). В отличие от
// типов 1-3 упругая часть есть и в покое, и при разлёте: контакт держит
// нагрузку, слот трения не пропадает. Демпфер действует только на сближении
// и при pdot -> 0 обращается в ноль непрерывно (без релейного включения).
// Kp -- ключ fiber_*_contact_penalty_scale (масштаб демпфера).
//
// Предел шага (явная схема с демпфером; множитель Белычко):
//     omega = sqrt(k_eff/m_red), k_eff = 1.5*F_H/delta (как у Герца),
//     xi    = c / (2 sqrt(k_eff m_red)),
//     dt    = min(0.8/omega, 1.6/omega * (sqrt(1+xi^2) - xi)).
// При xi = 0 это прежний предел Герца 0.8/omega; при xi -> inf --
// предел демпфера 1.6*m_red/c (с тем же запасом 0.8 от 2 m/c).
//
// damper_only = true -- только вязкая часть (для канала энергии
// D_contact_fiber: упругая часть потенциальна и учитывается в U_hertz_*).
// ----------------------------------------------------------------------
constexpr int kFiberLawHertzDamped = 4;

inline contact::Result fiber_hertz_damped_force(double f_hertz, double delta, const contact::PairGeom& g,
                                                const contact::PairProps& pr, double kp, bool damper_only = false)
{
    contact::Result out;
    if (!(delta > 0.0) || !(pr.m_red > 0.0)) return out;
    const double c = (pr.dt_ref > 0.0 && kp > 0.0) ? kp * pr.m_red / pr.dt_ref : 0.0;
    const double f_d = (g.pdot > 0.0) ? c * g.pdot : 0.0;
    out.fn = damper_only ? f_d : f_hertz + f_d;
    if (!(out.fn > 0.0)) return out;

    const double k_eff = (f_hertz > 0.0) ? 1.5 * f_hertz / delta : 0.0;
    if (k_eff > 0.0)
    {
        const double omega = std::sqrt(k_eff / pr.m_red);
        const double xi = c / (2.0 * std::sqrt(k_eff * pr.m_red));
        out.dt_limit = std::min(contact::kOscSafety / omega,
                                2.0 * contact::kOscSafety / omega * (std::sqrt(1.0 + xi * xi) - xi));
    }
    else if (c > 0.0)
    {
        out.dt_limit = contact::kOscSafety * 2.0 * pr.m_red / c;
    }
    return out;
}

// Закон типа t имеет упругий потенциал Герца (энергия контакта в U_hertz_*).
inline bool fiber_contact_law_has_hertz_potential(int type) { return type == 0 || type == kFiberLawHertzDamped; }

} // namespace interactions
