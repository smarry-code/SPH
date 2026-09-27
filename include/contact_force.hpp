#pragma once
// contact_force.hpp
// Штрафные контактные силы «частица-частица» (part_0.md, разд. 2.3; ADR-007).
// Чистые функции: вход — скалярные характеристики пары, выход — модуль
// нормальной силы и ограничитель шага. Направление и раскладка по осям —
// в оркестраторе (ext_force.cpp), поэтому файл не зависит от dim.
//
// SPH_HD (fiber/CUDA_PLAN.md §3.2, этап Г2, часть Solid-Contact): эти
// функции вызываются из physics/solid_contact_pair.hpp::
// compute_solid_contact_pair, поэтому обязаны собираться и для host, и для
// device (иначе __host__ __device__ функция не может вызвать чисто хостовую).
// std::min/std::max заменены на sph_min/sph_max (правило SPH_HD-кода, план
// §3.2) — та же операция сравнения, тот же результат.
#include <cmath>

#include "physics/hd.hpp" // SPH_HD, sph_min/sph_max
#include "precision.hpp"

namespace contact
{

// Значение «ограничитель не установлен» (совместимо с min() в time_integration).
inline constexpr real_t kNoLimit = real_t(1e30);
// Запас устойчивости явной схемы: dt*omega < 2, берём 0.8.
inline constexpr real_t kOscSafety = real_t(0.8);

enum ForceType : int { Off = 0, Pinball = 1, Swegle = 2, Impulse = 3 };

// Геометрия пары по (2.1). Заполняется оркестратором.
struct PairGeom
{
    real_t r     = 0;  // |x_i - x_j|, > 0
    real_t p     = 0;  // проникание 0.5*(h_i+h_j) - r; контакт при p > 0
    real_t pdot  = 0;  // скорость проникания, > 0 при сближении
    // Полная относительная скорость ||v_i - v_j||. Справочная: буквальная
    // (2.10) строила на ней силу типа 3, реализация перешла на нормальную
    // составляющую pdot (см. impulse_force). Оркестратор поле по-прежнему
    // заполняет — им меряется доля тангенциального скольжения в паре.
    real_t v_rel = 0;
};

// Свойства пары, не зависящие от геометрии.
struct PairProps
{
    real_t m_red  = 0; // приведённая масса m_i m_j/(m_i+m_j)
    real_t c_pair = 0; // min(c_i, c_j)
    real_t g_eff  = 0; // G_i G_j/(G_i+G_j); 0 => G неизвестен (ветка F2 отключается)
    real_t r_eff  = 0; // R_i R_j/(R_i+R_j), R = h/2
    real_t dt_ref = 0; // локальный CFL-шаг пары; НЕ применяемый шаг
};

struct Result
{
    real_t fn       = 0;      // модуль нормальной силы, Н (>= 0)
    real_t dt_limit = kNoLimit; // ограничитель шага; kNoLimit => не ограничивает
};

// Тип 1, (2.3)-(2.5). Множитель p из буквальной (2.3) опущен: F1 и F2 —
// силы, и p делал бы размерность Н*м (решение D1 в contact/02-decisions.md).
SPH_HD inline Result pinball_force(const PairGeom& g, const PairProps& pr, real_t kp)
{
    Result out;
    if (g.p <= 0 || pr.m_red <= 0 || pr.dt_ref <= 0 || kp <= 0) return out;

    const real_t f1 = (g.pdot > 0) ? pr.m_red * g.pdot / pr.dt_ref : real_t(0.0);
    const bool   hertz_known = (pr.g_eff > 0 && pr.r_eff > 0);
    const real_t f2 = hertz_known
                    ? pr.g_eff * std::sqrt(pr.r_eff) * g.p * std::sqrt(g.p)
                    : real_t(0.0);

    // Откат при неизвестном G: F2 исключается из min, а не обнуляет силу.
    const bool hertz_active = hertz_known && (f2 < f1);
    const real_t f = hertz_known ? sph_min(f1, f2) : f1;
    out.fn = kp * f;
    if (out.fn <= 0) return out;

    if (hertz_active)
    {
        // Герц: F ~ p^{3/2} => k_eff = dF/dp = 1.5 * F / p.
        const real_t k_eff = real_t(1.5) * out.fn / g.p;
        out.dt_limit = kOscSafety / std::sqrt(k_eff / pr.m_red);
    }
    else
    {
        // Вязкая ветка: F = c_d * pdot, c_d = kp*m_red/dt_ref.
        // Предел явной схемы для демпфера: dt < 2*m_red/c_d = 2*dt_ref/kp.
        out.dt_limit = kOscSafety * real_t(2.0) * pr.dt_ref / kp;
    }
    return out;
}

// Тип 2, (2.6)-(2.8). Во второй формуле первоисточника опечатка: она помечена
// F1, но по (2.6) это F2.
SPH_HD inline Result swegle_force(const PairGeom& g, const PairProps& pr, real_t kp)
{
    Result out;
    if (g.p <= 0 || g.pdot <= 0 || pr.m_red <= 0 ||
        pr.dt_ref <= 0 || g.r <= 0 || kp <= 0) return out;

    const real_t f1 = pr.m_red * g.p / (pr.dt_ref * pr.dt_ref);
    const real_t f2 = real_t(2.0) * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r);

    // Потолок: сила не больше той, что гасит сближение пары за один опорный
    // шаг. (2.6) выведена для p, накопленного ЗА ШАГ (p = pdot*dt_ref); в этом
    // режиме f1 = m_red*pdot/dt_ref ровно совпадает с потолком, поэтому в
    // области применимости формулы он неактивен. Превышение равно в точности
    // p/(pdot*dt_ref) — во сколько раз проникание больше накопленного
    // сближением, — и потолок вступает ровно тогда, когда допущение нарушено.
    //
    // Без потолка это разрушало кейс BC2 при contact_force_type = 2: приосевая
    // частица полиуретана разгружалась, particle_radius по mass/rho раздувал
    // порог R_i + R_j, и p росло у РАСХОДЯЩИХСЯ частиц. К шагу 1665 p = 1.7e-4
    // при pdot*dt_ref ~ 2.7e-8, то есть сила была завышена в ~6e3 раз: 4.5e8 Н
    // и 3.5e4 м/с приращения за шаг, плотность уходила в минус.
    const real_t f_cap  = pr.m_red * g.pdot / pr.dt_ref;
    const real_t f_raw  = sph_max(f1, f2);
    const bool   capped = (f_raw > f_cap);

    out.fn = kp * (capped ? f_cap : f_raw);
    if (out.fn <= 0) return out;

    if (capped)
    {
        // На потолке F = c_d*pdot, c_d = kp*m_red/dt_ref: та же вязкая ветка,
        // что у типа 1, и тот же предел явной схемы dt < 2*m_red/c_d.
        out.dt_limit = kOscSafety * real_t(2.0) * pr.dt_ref / kp;
    }
    else
    {
        const real_t k_eff = out.fn / g.p;   // F линейна по p
        out.dt_limit = kOscSafety / std::sqrt(k_eff / pr.m_red);
    }
    return out;
}

// Тип 3, (2.10). Импульсный контакт: за шаг снимается доля gamma нормальной
// скорости сближения пары.
//
// ДВА ОТСТУПЛЕНИЯ ОТ БУКВАЛЬНОЙ (2.10). Оба — следствие измерений на BC6
// (contact/06-e2e-bc6.md, разбор «зашумлённости» типа 3).
//
// 1. Вместо полной относительной скорости ||v_i - v_j|| берётся её НОРМАЛЬНАЯ
//    составляющая pdot. Буквальная (2.10) прикладывает вдоль нормали модуль,
//    посчитанный по полному вектору, поэтому:
//      - тангенциальное скольжение пары превращается в нормальное отталкивание
//        (силы трения в модели нет, гасить сдвиг нечем, а импульс уходит
//        в отрыв поверхностей);
//      - сила РАЗРЫВНА в точке включения: условие срабатывания pdot > 0, а
//        величина от pdot не зависит, поэтому при сколь угодно малом сближении
//        сила скачком принимает конечное значение m_red*||v_ij||/dt. Это
//        релейное управление без зоны нечувствительности; на BC6 оно давало
//        широкополосный дребезг на границе B4C/PMMA после прихода разгрузки:
//        СКО ряби кривой монитора 5.1 м/с против 0.12 м/с у типов 0-1,
//        средний скачок скорости 8.6 м/с за шаг против 0.04 м/с.
//    С pdot сила непрерывна: pdot -> +0 => fn -> 0.
//
// 2. Доля гашения gamma (отложенное решение ADR-007, «доля гашения gamma < 1
//    для типа 3»). gamma = 1 — буквальная (2.10): относительное сближение
//    гасится за шаг полностью, схема работает на границе апериодичности.
//    gamma < 1 оставляет паре (1-gamma)*pdot, то есть превращает релейный удар
//    в пропорциональный регулятор с запасом.
//
// ЧЕГО ЭТО НЕ ЧИНИТ. Приращение скорости за шаг
//   dv_i = fn*dt/m_i = gamma*(m_red/m_i)*pdot
// от шага НЕ зависит: fn ~ 1/dt, и dt сокращается. Это свойство самой
// постановки «погасить за один шаг», а не дефект реализации: измельчение CFL
// амплитуду контактного воздействия не уменьшает, и по шагу тип 3 не сходится.
// У типа 1 на герцевской ветви fn от dt не зависит и dv ~ dt — поэтому он
// гладкий. Считать тип 3 сеточно-сходящимся нельзя ни при каком gamma.
SPH_HD inline Result impulse_force(const PairGeom& g, const PairProps& pr, real_t dt_applied,
                            real_t gamma = real_t(1.0))
{
    Result out;
    if (g.p <= 0 || g.pdot <= 0 || pr.m_red <= 0 || dt_applied <= 0 ||
        gamma <= 0) return out;
    out.fn = gamma * pr.m_red * g.pdot / dt_applied;
    if (out.fn <= 0) return out;

    // Вязкий демпфер с c_d = gamma*m_red/dt_applied. Предел явной схемы
    // dt < 2*m_red/c_d = 2*dt_applied/gamma, с тем же запасом kOscSafety, что
    // у типов 1-2. При gamma <= 1.6 это >= dt_applied, то есть ограничитель не
    // срабатывает никогда и остаётся страховкой; при gamma > 1.6 он режет шаг,
    // а от урезанного шага зависит и сама жёсткость — петля «жёсткость -> шаг».
    // Поэтому gamma > 1.6 вне рабочего диапазона, и разбор конфигурации
    // предупреждает об этом (input.cpp).
    out.dt_limit = kOscSafety * real_t(2.0) * dt_applied / gamma;
    return out;
}

SPH_HD inline Result penalty_force(int type, const PairGeom& g, const PairProps& pr,
                            real_t kp, real_t dt_applied, real_t gamma = real_t(1.0))
{
    switch (type)
    {
        case Pinball: return pinball_force(g, pr, kp);
        case Swegle:  return swegle_force(g, pr, kp);
        case Impulse: return impulse_force(g, pr, dt_applied, gamma);
        default:      return Result{};
    }
}

} // namespace contact
