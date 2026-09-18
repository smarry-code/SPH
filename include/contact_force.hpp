#pragma once
// contact_force.hpp
// Штрафные контактные силы «частица-частица» (part_0.md, разд. 2.3; ADR-007).
// Чистые функции: вход — скалярные характеристики пары, выход — модуль
// нормальной силы и ограничитель шага. Направление и раскладка по осям —
// в оркестраторе (ext_force.cpp), поэтому файл не зависит от dim.
#include <algorithm>
#include <cmath>

namespace contact
{

// Значение «ограничитель не установлен» (совместимо с min() в time_integration).
inline constexpr double kNoLimit = 1e30;
// Запас устойчивости явной схемы: dt*omega < 2, берём 0.8.
inline constexpr double kOscSafety = 0.8;

enum ForceType : int { Off = 0, Pinball = 1, Swegle = 2, Impulse = 3 };

// Геометрия пары по (2.1). Заполняется оркестратором.
struct PairGeom
{
    double r     = 0.0;  // |x_i - x_j|, > 0
    double p     = 0.0;  // проникание 0.5*(h_i+h_j) - r; контакт при p > 0
    double pdot  = 0.0;  // скорость проникания, > 0 при сближении
    double v_rel = 0.0;  // ||v_i - v_j||, полная относительная скорость (2.10)
};

// Свойства пары, не зависящие от геометрии.
struct PairProps
{
    double m_red  = 0.0; // приведённая масса m_i m_j/(m_i+m_j)
    double c_pair = 0.0; // min(c_i, c_j)
    double g_eff  = 0.0; // G_i G_j/(G_i+G_j); 0 => G неизвестен (ветка F2 отключается)
    double r_eff  = 0.0; // R_i R_j/(R_i+R_j), R = h/2
    double dt_ref = 0.0; // локальный CFL-шаг пары; НЕ применяемый шаг
};

struct Result
{
    double fn       = 0.0;      // модуль нормальной силы, Н (>= 0)
    double dt_limit = kNoLimit; // ограничитель шага; kNoLimit => не ограничивает
};

// Тип 1, (2.3)-(2.5). Множитель p из буквальной (2.3) опущен: F1 и F2 —
// силы, и p делал бы размерность Н*м (решение D1 в contact/02-decisions.md).
inline Result pinball_force(const PairGeom& g, const PairProps& pr, double kp)
{
    Result out;
    if (g.p <= 0.0 || pr.m_red <= 0.0 || pr.dt_ref <= 0.0 || kp <= 0.0) return out;

    const double f1 = (g.pdot > 0.0) ? pr.m_red * g.pdot / pr.dt_ref : 0.0;
    const bool   hertz_known = (pr.g_eff > 0.0 && pr.r_eff > 0.0);
    const double f2 = hertz_known
                    ? pr.g_eff * std::sqrt(pr.r_eff) * g.p * std::sqrt(g.p)
                    : 0.0;

    // Откат при неизвестном G: F2 исключается из min, а не обнуляет силу.
    const bool hertz_active = hertz_known && (f2 < f1);
    const double f = hertz_known ? std::min(f1, f2) : f1;
    out.fn = kp * f;
    if (out.fn <= 0.0) return out;

    if (hertz_active)
    {
        // Герц: F ~ p^{3/2} => k_eff = dF/dp = 1.5 * F / p.
        const double k_eff = 1.5 * out.fn / g.p;
        out.dt_limit = kOscSafety / std::sqrt(k_eff / pr.m_red);
    }
    else
    {
        // Вязкая ветка: F = c_d * pdot, c_d = kp*m_red/dt_ref.
        // Предел явной схемы для демпфера: dt < 2*m_red/c_d = 2*dt_ref/kp.
        out.dt_limit = kOscSafety * 2.0 * pr.dt_ref / kp;
    }
    return out;
}

// Тип 2, (2.6)-(2.8). Во второй формуле первоисточника опечатка: она помечена
// F1, но по (2.6) это F2.
inline Result swegle_force(const PairGeom& g, const PairProps& pr, double kp)
{
    Result out;
    if (g.p <= 0.0 || g.pdot <= 0.0 || pr.m_red <= 0.0 ||
        pr.dt_ref <= 0.0 || g.r <= 0.0 || kp <= 0.0) return out;

    const double f1 = pr.m_red * g.p / (pr.dt_ref * pr.dt_ref);
    const double f2 = 2.0 * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r);
    out.fn = kp * std::max(f1, f2);
    if (out.fn <= 0.0) return out;

    const double k_eff = out.fn / g.p;   // F линейна по p
    out.dt_limit = kOscSafety / std::sqrt(k_eff / pr.m_red);
    return out;
}

// Тип 3, (2.10). От ПРИМЕНЯЕМОГО шага: импульс F*dt = m_red*||v_ij|| от шага
// не зависит, поэтому ограничитель не ставится и обратной связи нет.
inline Result impulse_force(const PairGeom& g, const PairProps& pr, double dt_applied)
{
    Result out;
    if (g.p <= 0.0 || g.pdot <= 0.0 || pr.m_red <= 0.0 || dt_applied <= 0.0) return out;
    out.fn = pr.m_red * g.v_rel / dt_applied;   // dt_limit остаётся kNoLimit
    return out;
}

inline Result penalty_force(int type, const PairGeom& g, const PairProps& pr,
                            double kp, double dt_applied)
{
    switch (type)
    {
        case Pinball: return pinball_force(g, pr, kp);
        case Swegle:  return swegle_force(g, pr, kp);
        case Impulse: return impulse_force(g, pr, dt_applied);
        default:      return Result{};
    }
}

} // namespace contact
