#pragma once
// interactions/fiber_contact_law.hpp
// Законы нормальной силы pinball / Swegle / импульсный для контактов нитей
// (Fiber-Solid и Fiber-Fiber; notes/fiber-contact/plan-fiber-contact-types.md, §2;
// notes/fiber-model/fiber-problem-statement.md, «Законы нормальной силы»).
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
#include <cstdint>
#include <limits>
#include <vector>

#include "contact_force.hpp"
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp" // fiber_bond_pair_dt
#include "interactions/fiber_normal_projection_types.hpp" // FiberNormalContact, FiberSolidNormalContact
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp" // dt_scale (impulse_unscaled_limit)

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
    // НЕ ИСПОЛЬЗОВАТЬ ДЛЯ РАСЧЁТОВ (решение владельца 2026-10-21; ключ оставлен для
    // воспроизведения серий 2026-10-20/21). Выталкивание на уровне скоростей -- источник
    // энергии: закон 3 гасит скорость пары за шаг при любой силе, и с выталкиванием пара,
    // которую натянутая нить прижимает обратно, разводится с заданной скоростью сотни шагов
    // подряд, совершая работу против натяжения. 1T 71 м/с: -19,8 % E0 за 5 мкс без предела,
    // -37,5 % E0 с пределом push_vmax = 20 м/с (предел ограничивает скорость, а не работу),
    // пачки по 12-16 ложных разрывов. Вернуть проникание без источника энергии можно только
    // потенциальной (упругой) силой от глубины -- это закон 4 (Герц + демпфер).
    //
    // Выталкивание типа 3 (стабилизация Баумгарте на уровне скоростей): при
    // p > push_depth*r закон гасит не pdot, а pdot + (p - push_depth*r)/push_tau,
    // то есть выводит лишнее проникание с постоянной времени push_tau, с.
    // Без него тип 3 держит нагрузку (F ~ gamma*P), но не возвращает
    // накопленное проникание: при pdot <= 0 сила ноль, и нити копят проникание
    // до сквозного прохода. Время, а не доля за шаг: доля beta за шаг дала бы
    // скорость выталкивания beta*excess/dt ~ сотни м/с при dt ~ 2e-8 с и
    // зависела бы от шага. Требование устойчивости: push_tau >> dt.
    // push_tau = 0 -- выключено (побитово прежний тип 3).
    double push_tau = 0.0;
    double push_depth = 0.1; // доля g.r (суммы радиусов пары)
    // Предел скорости выталкивания, м/с (0 -- без предела). Без него при внезапно глубоком
    // проникании (порванный конец влетает в перекрестье) (p - depth*r)/tau достигает сотен м/с,
    // и контакт выбрасывает энергию: 1T 71 м/с, 2026-10-21, -19,8 % E0 за 5 мкс и 12 ложных
    // разрывов. Под нагрузкой P пара за шаг набирает P*dt/m_red (4,5-9 м/с при 25-50 Н и
    // dt = 9e-9 с), выталкивание должно её перекрывать -- поэтому предел порядка 20 м/с, а не 1 м/с.
    // За одно срабатывание контакт может отдать не больше m_red*vmax^2/2 (~1e-5 Дж при 20 м/с).
    double push_vmax = 0.0;
};

// Ссылка на активный слот трения: индекс primitive k и номер отрезка seg (kDetSegments), в котором k стоял
// при обходе parallel_for_segments по ВСЕМ primitive. Списки таких ссылок (по возрастанию k) позволяют
// стадии трения обходить только активные слоты с прежним разбиением на отрезки -- порядок суммирования
// на частицу (буфер seg -> позиция в буфере) тот же, что при обходе плотного массива слотов.
struct FrictionActiveRef
{
    std::size_t k = 0;
    std::uint16_t seg = 0;
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
    // РЕЖИМ ПРОЕКЦИИ (законы 1-3, notes/fiber-contact/law3-projection-tz.md): если задан, ядро НЕ считает
    // нормальную силу, а заполняет нормальные слоты для стадии run_normal_projection -- сжатым
    // массивом только активных контактов (delta > 0) по возрастанию primitive, с индексом
    // primitive (prim) и номером отрезка (seg); вкладов в ускорение, слота трения и предела шага
    // (в pair_dt_limit -- kNoLimit) нет. normal_ff -- ядро нить-нить, normal_fs -- нить-Solid.
    // nullptr -- прежний путь, побитово.
    std::vector<FiberNormalContact>* normal_ff = nullptr;
    std::vector<FiberSolidNormalContact>* normal_fs = nullptr;
    // Рабочая память сбора нормальных слотов по отрезкам (живёт между шагами, память только
    // растёт; заголовки векторов -- в своих кэш-линиях). nullptr -- ядро заводит временную.
    std::vector<rhs::CpuExecutor::SegPart<FiberNormalContact>>* normal_ff_parts = nullptr;
    std::vector<rhs::CpuExecutor::SegPart<FiberSolidNormalContact>>* normal_fs_parts = nullptr;
    // Только тип 0: если задан, ядро отдаёт предел шага контакта (минимум fiber_contact_pair_dt_limit по
    // принятым контактам; max(), если контактов нет) -- тот же результат, что compute_*_dt_limit для
    // типа 0 (та же свёртка по отрезкам kDetSegments, побитово), без повторного обхода primitive.
    // nullptr -- не считается.
    double* law0_dt_min_out = nullptr;
    // РЕЖИМ СПИСКА слотов трения: если fric_active задан, плотный массив слотов friction_contacts не
    // обнуляется целиком на каждом шаге -- ядро сбрасывает только слоты из списка ПРОШЛОГО вызова, а
    // после прохода заполняет список активными слотами (по возрастанию k, с номером отрезка).
    // Инвариант: каждый слот вне списка равен FiberFrictionContact{} / FiberSolidFrictionContact{}.
    // Список ведёт тот же вызывающий код, что и массив слотов (в режиме проекции его перезаполняет
    // fill_friction_slots_from_normal_*). nullptr -- прежний путь с полным сбросом, побитово.
    std::vector<FrictionActiveRef>* fric_active = nullptr;
    // Рабочая память сбора списка по отрезкам (память только растёт); nullptr -- временная.
    std::vector<rhs::CpuExecutor::SegPart<FrictionActiveRef>>* fric_active_parts = nullptr;
};

// Есть ли у закона новые (не Герц) типы -- нужен ли dt_bead.
inline bool fiber_contact_law_is_penalty(const FiberContactLaw& law) { return law.type != 0; }

// ----------------------------------------------------------------------
// compute_fiber_bead_dt_ref -- dt_bead (см. шапку файла). O(nbonds),
// последовательно: минимум не зависит от порядка, результат детерминирован
// при любом числе потоков. В основном шаге не вызывается каждый шаг: значения
// хранит FiberDtRefCache (вариант с executor ниже), пересчёт -- только при смене
// топологии; этот последовательный вариант зовёт диагностика fiber_contact_stats.
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
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(fiber::bond_material_id(b, state))];
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

// Вариант с executor: dt по связям -- параллельно (та же fiber_bond_pair_dt, те же аргументы), затем
// последовательная раскладка минимума по бусинам в том же порядке связей -- результат побитово прежний.
inline std::vector<double>& compute_fiber_bead_dt_ref(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                                       const fiber::FiberTopology& topology,
                                                       const std::vector<fiber::FiberMaterial>& materials,
                                                       std::vector<double>& out)
{
    const double kInf = std::numeric_limits<double>::max();
    const std::size_t n = state.size();
    out.assign(n, kInf);
    const std::size_t nb = topology.bonds.size();
    std::vector<double> bond_dt(nb);
    exec.parallel_for(nb, [&](std::size_t q) {
        const fiber::FiberBond& b = topology.bonds[q];
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(fiber::bond_material_id(b, state))];
        bond_dt[q] = fiber_bond_pair_dt(mat, b.l0, state.mass[b.i], state.mass[b.j]);
    });
    for (std::size_t q = 0; q < nb; ++q)
    {
        const fiber::FiberBond& b = topology.bonds[q];
        const double dt = bond_dt[q];
        out[b.i] = std::min(out[b.i], dt);
        out[b.j] = std::min(out[b.j], dt);
    }
    exec.parallel_for(n, [&](std::size_t k) {
        if (out[k] < kInf) return;
        // Виртуальная связь длины d: E*A/d, m_red = m_k*m_k/(2 m_k) = m_k/2.
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(state.material_id[k])];
        out[k] = fiber_bond_pair_dt(mat, mat.d, state.mass[k], state.mass[k]);
    });
    return out;
}

// ----------------------------------------------------------------------
// FiberDtRefCache -- кэш опорных пределов шага нитей: bond_dt_limit
// (compute_fiber_bond_dt_limit) и bead_dt_ref (compute_fiber_bead_dt_ref).
// Обе функции зависят ТОЛЬКО от topology.bonds[*].{i,j,l0}, state.material_id,
// state.mass и materials (не от x, v, состояния связи, повреждения, dt), поэтому
// пересчёт нужен лишь при смене топологии (topology.version растёт в
// apply_fiber_break_events), размеров или числа материалов. Значения считаются
// теми же функциями -- побитово равны прямому вызову. Масса/материалы бусин и l0
// связей в ходе расчёта не меняются; при их ручной правке вызывать invalidate().
// ----------------------------------------------------------------------
struct FiberDtRefCache
{
    bool valid = false;
    unsigned long long version = 0;
    std::size_t nbonds = 0, n = 0, nmat = 0;
    bool have_bead = false;
    std::vector<double> bead_dt_ref;
    double bond_dt_limit = std::numeric_limits<double>::max();

    void ensure(const rhs::CpuExecutor& exec, const fiber::FiberState& state, const fiber::FiberTopology& topology,
                const std::vector<fiber::FiberMaterial>& materials, bool need_bead)
    {
        const bool same_key = valid && version == topology.version && nbonds == topology.bonds.size() &&
                              n == state.size() && nmat == materials.size();
        if (same_key && (!need_bead || have_bead)) return;
        bond_dt_limit = compute_fiber_bond_dt_limit(exec, state, topology, materials);
        have_bead = false;
        if (need_bead)
        {
            compute_fiber_bead_dt_ref(exec, state, topology, materials, bead_dt_ref);
            have_bead = true;
        }
        version = topology.version;
        nbonds = topology.bonds.size();
        n = state.size();
        nmat = materials.size();
        valid = true;
    }

    void invalidate() { valid = false; }
};

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

// Предел шага типа 3 (2*kOscSafety*dt_prev/gamma) выражен через сам прошлый шаг. При
// cfg.dt_scale = s < 1 общий шаг умножается на s ПОСЛЕ минимума по пределам, и этот
// предел замыкал бы петлю dt_new = s*1.6*dt_prev: шаг падал бы в 0.8 раза за шаг до
// нуля (B-этап 2026-10-20: dt_min = 3e-29 за 2000 шагов). Поэтому предел типа 3
// делится на s: после умножения он снова 1.6*dt_prev/gamma. При s = 1 -- без изменений.
inline contact::Result impulse_unscaled_limit(const FiberContactStepParams& params, contact::Result out)
{
    if (params.cfg != nullptr && params.cfg->dt_scale < 1.0 && params.cfg->dt_scale > 0.0 &&
        out.dt_limit < contact::kNoLimit)
        out.dt_limit /= params.cfg->dt_scale;
    return out;
}

// = contact::penalty_force(type, g, pr, kp, dt_prev, gamma). Только для type 1-3.
// У типа 3 с push_tau > 0 к pdot добавляется скорость выталкивания (см.
// FiberContactLaw::push_tau); сама формула impulse_force не меняется.
inline contact::Result fiber_penalty_force(const FiberContactStepParams& params, const contact::PairGeom& g,
                                           const contact::PairProps& pr)
{
    if (params.law.type == contact::Impulse && params.law.push_tau > 0.0)
    {
        const double excess = g.p - params.law.push_depth * g.r;
        if (excess > 0.0)
        {
            contact::PairGeom gp = g;
            double v_push = excess / params.law.push_tau;
            if (params.law.push_vmax > 0.0 && v_push > params.law.push_vmax) v_push = params.law.push_vmax;
            gp.pdot += v_push;
            return impulse_unscaled_limit(params, contact::impulse_force(gp, pr, params.dt_prev, params.law.gamma));
        }
    }
    contact::Result out = contact::penalty_force(params.law.type, g, pr, params.law.kp, params.dt_prev, params.law.gamma);
    return params.law.type == contact::Impulse ? impulse_unscaled_limit(params, out) : out;
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
