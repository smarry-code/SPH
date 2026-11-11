#pragma once
// interactions/fiber_normal_projection.hpp
// Нормальный контакт нитей (законы 1-3) как ПРОЕКЦИЯ НА УРОВНЕ СКОРОСТЕЙ
// (fiber/LAW3_VARIANT1_VELOCITY_PROJECTION.md, fiber/TZ_LAW3_PROJECTION.md,
// этап A3). Стадия стоит между выбором шага Δt и стадией трения
// (src/time_integration.cpp), по образцу interactions/fiber_friction.hpp:
// слоты, счётчики c_k/c_s, пробные скорости от БАЗЫ решателя
// (ode::FiberVelocityUpdate / ode::SolidVelocityUpdate), вклады через
// ContributionWriter + reduce_contributions_into.
//
// ----------------------------------------------------------------------------
// ПОСТАНОВКА (док. §2-§4). Контакт c: сторона A (нить: бусины a0,a1 с весами
// (1-s, s); Solid: частица s), сторона B (нить: b0,b1 с весами (1-t, t)),
// нормаль n_c от B к A (нить-нить) / от нити к Solid (нить-Solid, тогда A -- Solid).
// Скорость сближения
//     pdot_c = -n_c . (v_A - v_B).
// Неизвестный импульс P_c >= 0 за шаг (N_c = P_c/h). Импульс +P n действует на
// A, -P n на B: бусина k получает Δv_k = sigma w P n / m_k (sigma = +1 на A,
// -1 на B); для нить-Solid частица Solid получает +P n/m_s, бусины нити -w P n/m_k.
// Условие Синьорини: 0 <= P_c  ⊥  -pdot_c^{n+1} >= 0.
//
// ПРОБНАЯ СКОРОСТЬ -- ровно та же формула, что у advance() и у трения:
//     нить:  v~_k = base_k + h * dv_k;
//     Solid: v~_s = base_s + h * dvx_s + av_s   (av -- XSPH; nullptr -- нет).
// dv_fiber / dvx_solid -- ускорения шага БЕЗ нормального контакта законов 1-3
// (после проекции зажимов, до трения). h для нити и Solid ОБЯЗАН совпадать
// (fupd.h == supd.h: импульс P -- один на контакт); в стадии это проверяется assert.
//
// ИТЕРАЦИИ ЯКОБИ С РАЗДЕЛЕНИЕМ МАССЫ (док. §4.2). Для каждого активного слота
//     W~_c = sum_{бусины} c_k w^2/m_k  (+ c_s/m_s для нить-Solid),
//     P_new = max(0, P_old + pdot_c/W~_c),     dP = P_new - P_old.
// Накопленный импульс зажимается, приращение dP -- нет (Catto): импульс, лишний
// после соседних контактов, можно снять. Бусина получает ПОЛНОЕ dv от всех
// контактов, а контакт при расчёте dP видит только долю 1/c_k массы -- поэтому
// число контактов на бусину не вызывает перекомпенсации. Порядок стадии:
//   0. пробные скорости; 1. начальные импульсы P^(0) (тёплый старт, уже лежат в
//   slot.P), один проход «применить»; 1'. энергетическая безопасность тёплого
//   старта: если его работа W0 > 0, все P^(0) и Δv умножаются на lambda < 1 так,
//   что W0 = 0 (stats.warm_scale; вывод -- в run_normal_projection); 2. K итераций
//   (parallel_for_segments по слотам -> вклады -> reduce -> Δv += ), каждая не
//   увеличивает работу стадии (D >= A); 3. работа стадии W <= 0; 4. диагностика;
//   5. Δv/h добавляется в dv_fiber / dvx_solid (ускорение), slot.P -- итоговый
//   накопленный импульс. Остаток после стадии трения -- normal_projection_post_residual.
//
// ЗАЖИМЫ. Счётчик c = 0 у бусины/частицы, входящей в слот, означает бесконечную
// массу: вклад в W~ нулевой, Δv не получает (пробная скорость остаётся скоростью
// опоры). Вызывающая сторона получает c из compute_normal_counts и затем
// вызывает interactions::apply_clamp_friction_counts (interactions/clamp.hpp),
// которая зануляет c у полностью закреплённых бусин и частиц -- ровно как у трения.
//
// ДЕТЕРМИНИЗМ. Только parallel_for_segments (фиксированные отрезки,
// det_segments()), ContributionWriter + reduce_contributions_into и
// parallel_reduce (фиксированные отрезки, затем последовательное объединение).
// Ни atomics, ни mutex, ни неупорядоченных контейнеров: результат побитово
// одинаков при любом числе потоков. Слот меняется только своим собственным
// вызовом parallel_for_segments (без гонок).
//
// ТЁПЛЫЙ СТАРТ -- NormalWarmStart: отсортированный по ключу пары массив
// (бинарный поиск), P^(0) = P_prev * h / h_prev (импульс под нагрузкой
// пропорционален шагу). Таблица экспортируется/импортируется для чекпоинта.
//
// ТАБЛИЦА interaction_type: домен FiberDelta -- 1 Bond, 2 нить-нить контакт,
// 3 нить-Solid контакт, 4/5 трение, 6 -- нормаль нить-нить, 7 -- нормаль
// нить-Solid (сторона нити); домен Solid -- 10/11/12 заняты, 13 -- трение
// Fiber-Solid, 14 -- нормаль нить-Solid (сторона Solid).
// ============================================================================
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "interactions/fiber_contact.hpp"       // FiberFrictionContact
#include "interactions/fiber_friction.hpp"      // FiberCountDelta, SolidCountDelta, detail_friction::prepare
#include "interactions/fiber_normal_projection_types.hpp"
#include "interactions/fiber_solid_contact.hpp" // FiberSolidFrictionContact
#include "ode/ode_solver.hpp"                   // ode::FiberVelocityUpdate, ode::SolidVelocityUpdate
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"

namespace interactions
{

constexpr std::uint32_t kFiberNormalFiberFiberType = 6;
constexpr std::uint32_t kFiberNormalFiberSolidFiberType = 7;
constexpr std::uint32_t kFiberNormalFiberSolidSolidType = 14;

// Вклад итерации в приращение скорости бусины / частицы Solid (одна структура
// на оба домена: нужны только три компоненты Δv).
struct NormalDelta
{
    double dv[dim] = {};
    NormalDelta& operator+=(const NormalDelta& o)
    {
        for (std::size_t d = 0; d < dim; ++d) dv[d] += o.dv[d];
        return *this;
    }
};

// ----------------------------------------------------------------------
// normal_projection_impulse -- формулы ОДНОГО контакта (без исключений и
// аллокаций, кандидат на __host__ __device__). pdot -- текущая скорость
// сближения, Wtilde -- W~_c с разделением массы, P_old -- накопленный импульс.
//   P_new = max(0, P_old + pdot/W~),  dP = P_new - P_old.
// Wtilde <= 0 (обе стороны с бесконечной массой): P_new = P_old, dP = 0.
// ----------------------------------------------------------------------
inline void normal_projection_impulse(double pdot, double Wtilde, double P_old, double& P_new, double& dP)
{
    if (!(Wtilde > 0.0))
    {
        P_new = P_old;
        dP = 0.0;
        return;
    }
    P_new = std::max(0.0, P_old + pdot / Wtilde);
    dP = P_new - P_old;
}

// ----------------------------------------------------------------------
// *_for_each_bead -- бусины слота с НЕНУЛЕВЫМ весом: f(bead, weight, sigma,
// slot). sigma -- знак относительно стороны A (его стороны A: +1).
// Нить-нить: slot 0..3 = a0,a1,b0,b1; sigma +1 на A, -1 на B (как у трения).
// Нить-Solid: slot 0..1 = a,b; sigma = -1 (нить -- сторона B, Solid -- A;
// Solid обрабатывается отдельно с sigma = +1). Единственное место, где задано
// «нулевой вес»: им пользуются и счётчики, и скорость, и эмиссия вкладов.
// ----------------------------------------------------------------------
template <class F>
inline void fiber_normal_for_each_bead(const FiberNormalContact& c, F&& f)
{
    if (c.s < 1.0) f(c.a0, 1.0 - c.s, 1.0, static_cast<std::uint16_t>(0));
    if (c.a1 != -1 && c.s > 0.0) f(static_cast<rhs::index_t>(c.a1), c.s, 1.0, static_cast<std::uint16_t>(1));
    if (c.t < 1.0) f(c.b0, 1.0 - c.t, -1.0, static_cast<std::uint16_t>(2));
    if (c.b1 != -1 && c.t > 0.0) f(static_cast<rhs::index_t>(c.b1), c.t, -1.0, static_cast<std::uint16_t>(3));
}

template <class F>
inline void fiber_solid_normal_for_each_bead(const FiberSolidNormalContact& c, F&& f)
{
    if (c.lambda < 1.0) f(c.a, 1.0 - c.lambda, -1.0, static_cast<std::uint16_t>(0));
    if (c.b != -1 && c.lambda > 0.0) f(static_cast<rhs::index_t>(c.b), c.lambda, -1.0, static_cast<std::uint16_t>(1));
}

// ----------------------------------------------------------------------
// Рабочая память стадии (владеет вызывающая сторона, переиспользуется между
// шагами, память только растёт).
//   bead_count / solid_count -- c_k, c_s: заполняет compute_normal_counts, затем
//       зажимы (apply_clamp_friction_counts); run_normal_projection их ЧИТАЕТ;
//   dvel_fiber / dvel_solid  -- накопленные Δv стадии (размер nfiber / nsolid);
//   остальное -- буферы вкладов и reduce.
// ----------------------------------------------------------------------
struct NormalProjectionWorkspace
{
    std::vector<double> bead_count, solid_count;

    std::vector<double> dvel_fiber[dim];
    std::vector<double> dvel_solid[dim];

    std::vector<rhs::ContributionBuffer<FiberCountDelta>> fiber_count_buffers;
    std::vector<rhs::ContributionBuffer<SolidCountDelta>> solid_count_buffers;
    std::vector<FiberCountDelta> fiber_counts;
    std::vector<SolidCountDelta> solid_counts;

    std::vector<rhs::ContributionBuffer<NormalDelta>> fiber_buffers, solid_buffers;
    std::vector<NormalDelta> fiber_sum, solid_sum;
    rhs::ReduceScratch scratch;

    // Нулевые Δv для normal_projection_post_residual (скорость = база + h*dv итоговое).
    std::vector<double> zero_fiber[dim];
    std::vector<double> zero_solid[dim];
};

// ----------------------------------------------------------------------
// Счётчики c_k, c_s по активным нормальным слотам ОБОИХ видов: c_k -- число
// слотов нить-нить и нить-Solid, в которых бусина k имеет ненулевой вес; c_s --
// число слотов нить-Solid с частицей s. Результат -- ws.bead_count (nfiber),
// ws.solid_count (nsolid). После вызова нужно применить зажимы
// (apply_clamp_friction_counts(exec, clamp_ws, ws.bead_count, ws.solid_count)).
// ----------------------------------------------------------------------
namespace detail_normal
{

inline void emit_counts_ff(const rhs::CpuExecutor& exec, const std::vector<FiberNormalContact>& slots,
                            std::vector<rhs::ContributionBuffer<FiberCountDelta>>& buffers)
{
    const std::size_t n = slots.size();
    if (n == 0) return;
    exec.parallel_for_segments(n, [&](std::size_t seg, std::size_t k) {
        const FiberNormalContact& c = slots[k];
        if (!c.active) return;
        rhs::ContributionWriter<FiberCountDelta> writer(buffers[seg], kFiberNormalFiberFiberType);
        FiberCountDelta one;
        one.c = 1.0;
        fiber_normal_for_each_bead(c, [&](rhs::index_t bead, double, double, std::uint16_t slot) {
            writer.emit(bead, one, static_cast<std::uint64_t>(k), slot);
        });
    });
}

inline void emit_counts_fs(const rhs::CpuExecutor& exec, const std::vector<FiberSolidNormalContact>& slots,
                            std::vector<rhs::ContributionBuffer<FiberCountDelta>>& fiber_buffers,
                            std::vector<rhs::ContributionBuffer<SolidCountDelta>>& solid_buffers)
{
    const std::size_t n = slots.size();
    if (n == 0) return;
    exec.parallel_for_segments(n, [&](std::size_t seg, std::size_t k) {
        const FiberSolidNormalContact& c = slots[k];
        if (!c.active) return;
        rhs::ContributionWriter<FiberCountDelta> fw(fiber_buffers[seg], kFiberNormalFiberSolidFiberType);
        rhs::ContributionWriter<SolidCountDelta> sw(solid_buffers[seg], kFiberNormalFiberSolidSolidType);
        FiberCountDelta one;
        one.c = 1.0;
        fiber_solid_normal_for_each_bead(c, [&](rhs::index_t bead, double, double, std::uint16_t slot) {
            fw.emit(bead, one, static_cast<std::uint64_t>(k), slot);
        });
        SolidCountDelta ones;
        ones.c = 1.0;
        sw.emit(c.s, ones, static_cast<std::uint64_t>(k), 0);
    });
}

} // namespace detail_normal

inline void compute_normal_counts(const rhs::CpuExecutor& exec, std::size_t nfiber, std::size_t nsolid,
                                   const std::vector<FiberNormalContact>& ff_slots,
                                   const std::vector<FiberSolidNormalContact>& fs_slots,
                                   NormalProjectionWorkspace& ws)
{
    detail_friction::prepare(exec, ws.fiber_count_buffers);
    detail_friction::prepare(exec, ws.solid_count_buffers);
    ws.bead_count.assign(nfiber, 0.0);
    ws.solid_count.assign(nsolid, 0.0);

    detail_normal::emit_counts_ff(exec, ff_slots, ws.fiber_count_buffers);
    detail_normal::emit_counts_fs(exec, fs_slots, ws.fiber_count_buffers, ws.solid_count_buffers);

    if (!ff_slots.empty() || !fs_slots.empty())
    {
        rhs::reduce_contributions_into(exec, ws.fiber_count_buffers, nfiber, ws.fiber_counts, ws.scratch);
        for (std::size_t i = 0; i < nfiber; ++i) ws.bead_count[i] = ws.fiber_counts[i].c;
    }
    if (!fs_slots.empty())
    {
        rhs::reduce_contributions_into(exec, ws.solid_count_buffers, nsolid, ws.solid_counts, ws.scratch);
        for (std::size_t i = 0; i < nsolid; ++i) ws.solid_count[i] = ws.solid_counts[i].c;
    }
}

// ----------------------------------------------------------------------
// Текущие скорости стадии и вычисление pdot / W~ одного слота. Скорость =
// пробная (база решателя + h*dv [+ av]) + накопленное Δv стадии.
// ----------------------------------------------------------------------
struct NormalStageView
{
    const fiber::FiberState& state;
    const ode::FiberVelocityUpdate& fupd;
    const std::vector<double> (&dv_fiber)[dim];
    const std::vector<double> (&acc_fiber)[dim];
    const ParticleSystem& sys;
    const ode::SolidVelocityUpdate& supd;
    const std::vector<std::vector<double>>& dvx_solid;
    const std::vector<double> (&acc_solid)[dim];
    const double* bead_count;
    const double* solid_count;

    double trial_fiber(std::size_t d, rhs::index_t k) const { return (*fupd.base[d])[k] + fupd.h * dv_fiber[d][k]; }
    double trial_solid(std::size_t d, rhs::index_t s) const
    {
        double v = (*supd.base[d])[s] + supd.h * dvx_solid[d][s];
        if (supd.av) v += (*supd.av)[d][s];
        return v;
    }
    double fiber_v(std::size_t d, rhs::index_t k) const { return trial_fiber(d, k) + acc_fiber[d][k]; }
    double solid_v(std::size_t d, rhs::index_t s) const { return trial_solid(d, s) + acc_solid[d][s]; }

    // pdot_c = -n.(v_A - v_B) и W~_c нить-нить. Необязательный выход Wplain --
    // W_c = sum w^2/m_k БЕЗ счётчиков c_k (= 1/m_red контакта, только диагностика
    // остатка): бусина с c_k = 0 (зажим, бесконечная масса) в сумму не входит.
    void eval_ff(const FiberNormalContact& c, double& pdot, double& Wt, double* Wplain = nullptr) const
    {
        double rel[dim] = {};
        double W = 0.0;
        double Wp = 0.0;
        fiber_normal_for_each_bead(c, [&](rhs::index_t k, double w, double sg, std::uint16_t) {
            for (std::size_t d = 0; d < dim; ++d) rel[d] += sg * w * fiber_v(d, k);
            W += bead_count[k] * (w * w) / state.mass[k];
            if (bead_count[k] > 0.0) Wp += (w * w) / state.mass[k];
        });
        double rn = 0.0;
        for (std::size_t d = 0; d < dim; ++d) rn += c.n[d] * rel[d];
        pdot = -rn;
        Wt = W;
        if (Wplain != nullptr) *Wplain = Wp;
    }

    // Нить-Solid: A = Solid, B = точка нити; n от нити к Solid. Wplain -- как у eval_ff
    // (частица Solid с c_s = 0 тоже не входит).
    void eval_fs(const FiberSolidNormalContact& c, double& pdot, double& Wt, double* Wplain = nullptr) const
    {
        double rel[dim];
        for (std::size_t d = 0; d < dim; ++d) rel[d] = solid_v(d, c.s);
        double W = solid_count[c.s] / sys.mass[c.s];
        double Wp = solid_count[c.s] > 0.0 ? 1.0 / sys.mass[c.s] : 0.0;
        fiber_solid_normal_for_each_bead(c, [&](rhs::index_t k, double w, double sg, std::uint16_t) {
            for (std::size_t d = 0; d < dim; ++d) rel[d] += sg * w * fiber_v(d, k);
            W += bead_count[k] * (w * w) / state.mass[k];
            if (bead_count[k] > 0.0) Wp += (w * w) / state.mass[k];
        });
        double rn = 0.0;
        for (std::size_t d = 0; d < dim; ++d) rn += c.n[d] * rel[d];
        pdot = -rn;
        Wt = W;
        if (Wplain != nullptr) *Wplain = Wp;
    }
};

namespace detail_normal
{

// Агрегат диагностики по слотам (после последней итерации).
struct SlotAgg
{
    std::size_t active = 0;
    std::size_t warm = 0;
    double resid = 0.0;       // max pdot^+
    double resid_rel = 0.0;   // max pdot^+ / (W_c P_c) по слотам с P_c > 0 и W_c > 0 (W_c без c_k)
    double sum_missing = 0.0; // sum_c pdot_c^+ / W~_c по активным слотам с W~_c > 0 (недостающий импульс)
    double sum_P = 0.0;       // sum_c P_c по тем же слотам
};

inline SlotAgg agg_combine(SlotAgg a, const SlotAgg& b)
{
    a.active += b.active;
    a.warm += b.warm;
    a.resid = std::max(a.resid, b.resid);
    a.resid_rel = std::max(a.resid_rel, b.resid_rel);
    a.sum_missing += b.sum_missing;
    a.sum_P += b.sum_P;
    return a;
}

// Wt -- W~_c (с разделением массы), Wc -- W_c = sum w^2/m_k без c_k (1/m_red).
inline SlotAgg agg_one(bool warm, double pdot, double Wt, double Wc, double P)
{
    SlotAgg a;
    a.active = 1;
    a.warm = warm ? 1u : 0u;
    const double pp = std::max(0.0, pdot);
    a.resid = pp;
    if (P > 0.0 && Wc > 0.0) a.resid_rel = pp / (Wc * P);
    if (Wt > 0.0)
    {
        a.sum_missing = pp / Wt;
        a.sum_P = P;
    }
    return a;
}

} // namespace detail_normal

// ----------------------------------------------------------------------
// run_normal_projection -- вся стадия (см. шапку).
//   fupd/supd              -- ode::FiberVelocityUpdate / SolidVelocityUpdate (как у трения);
//   dv_fiber / dvx_solid   -- ускорения шага без нормали; на выходе += Δv/h стадии;
//                             dvx_solid[d][s], размер >= nsolid = ws.solid_count.size();
//   ff_slots / fs_slots    -- нормальные слоты; на входе slot.P = P^(0) (тёплый старт,
//                             NormalWarmStart::apply или 0), на выходе -- накопленный
//                             импульс шага. Неактивные слоты не меняются;
//   iters                  -- K (0 -- только начальные импульсы);
//   ws                     -- рабочая память; ws.bead_count/solid_count УЖЕ посчитаны
//                             compute_normal_counts и скорректированы зажимами;
//   stats                  -- диагностика шага (заполняется целиком, enabled = true).
// Диагностика:
//   resid     = max_c (pdot_c^{n+1})^+ после K итераций (один дополнительный проход);
//   resid_rel = max_c pdot_c^+/(W_c P_c) по слотам с P_c > 0, W_c > 0, где
//               W_c = sum w^2/m_k (+ 1/m_s) БЕЗ счётчиков c_k (= 1/m_red контакта; бусины
//               и частицы с c = 0, то есть зажимы, в сумму не входят). Это остаток в долях
//               h N_c/m_red (ТЗ A4): импульс pdot^+/W_c, которого не хватает «честному»
//               контакту двух масс, делённый на накопленный P_c. Метрика НЕ МЕНЬШЕ
//               прежней (с W~_c): при c_k >= 1 W~_c >= W_c, то есть строже;
//   resid_load = (sum_c pdot_c^+ / W~_c) / (sum_c P_c) по ВСЕМ активным слотам с W~_c > 0:
//               недостающий импульс (в терминах итерации Якоби, с W~_c) в долях полной
//               нагрузки; в отличие от resid_rel учитывает и слоты с P_c = 0 (контакт,
//               который итерации не успели «включить»). sum P = 0 -- 0. Суммы идут по
//               фиксированным отрезкам parallel_reduce (детерминированно);
//   warm_frac = (число активных слотов с slot.warm) / (число активных слотов);
//   W_stage   = sum m (|v|^2 - |v~|^2)/2 = sum m (v~.Δv + |Δv|^2/2) по бусинам и частицам,
//               фиксированные отрезки parallel_reduce, затем последовательное сложение;
//   max_ck/max_cs -- максимумы ws.bead_count / ws.solid_count (после зажимов).
// ----------------------------------------------------------------------
inline void run_normal_projection(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                   const ode::FiberVelocityUpdate& fupd, std::vector<double> (&dv_fiber)[dim],
                                   std::vector<FiberNormalContact>& ff_slots, const ParticleSystem& sys,
                                   const ode::SolidVelocityUpdate& supd,
                                   std::vector<std::vector<double>>& dvx_solid,
                                   std::vector<FiberSolidNormalContact>& fs_slots, int iters,
                                   NormalProjectionWorkspace& ws, NormalProjectionStats& stats)
{
    stats = NormalProjectionStats{};
    stats.enabled = true;
    stats.iters = iters;

    const std::size_t nff = ff_slots.size();
    const std::size_t nfs = fs_slots.size();
    const std::size_t nfiber = state.size();
    const std::size_t nsolid = ws.solid_count.size();
    assert(ws.bead_count.size() == nfiber);
    assert(std::fabs(fupd.h - supd.h) <= 1e-12 * std::fabs(fupd.h)); // один h на нить и Solid
    assert(fupd.h > 0.0);

    for (std::size_t d = 0; d < dim; ++d)
    {
        ws.dvel_fiber[d].assign(nfiber, 0.0);
        ws.dvel_solid[d].assign(nsolid, 0.0);
    }
    if (nff == 0 && nfs == 0) return;

    const NormalStageView view{state,    fupd, dv_fiber, ws.dvel_fiber,   sys,
                               supd,     dvx_solid, ws.dvel_solid, ws.bead_count.data(), ws.solid_count.data()};

    // --- один проход: эмиссия вкладов Δv. apply_initial -- применить slot.P как есть
    //     (тёплый старт), иначе -- обновить P по текущим скоростям.
    auto emit_pass = [&](bool apply_initial) {
        detail_friction::prepare(exec, ws.fiber_buffers);
        if (nfs > 0) detail_friction::prepare(exec, ws.solid_buffers);

        exec.parallel_for_segments(nff, [&](std::size_t seg, std::size_t k) {
            FiberNormalContact& c = ff_slots[k];
            if (!c.active) return;
            double dP;
            if (apply_initial)
            {
                if (!(c.P > 0.0)) c.P = 0.0; // P >= 0, NaN -> 0
                dP = c.P;
            }
            else
            {
                double pdot, Wt, Pn;
                view.eval_ff(c, pdot, Wt);
                normal_projection_impulse(pdot, Wt, c.P, Pn, dP);
                c.P = Pn;
            }
            if (dP == 0.0) return;
            rhs::ContributionWriter<NormalDelta> writer(ws.fiber_buffers[seg], kFiberNormalFiberFiberType);
            fiber_normal_for_each_bead(c, [&](rhs::index_t bead, double w, double sg, std::uint16_t slot) {
                if (!(ws.bead_count[bead] > 0.0)) return; // бесконечная масса (зажим)
                const double m = state.mass[bead];
                NormalDelta dl;
                for (std::size_t d = 0; d < dim; ++d) dl.dv[d] = (sg * w * dP) * c.n[d] / m;
                writer.emit(bead, dl, static_cast<std::uint64_t>(k), slot);
            });
        });

        exec.parallel_for_segments(nfs, [&](std::size_t seg, std::size_t k) {
            FiberSolidNormalContact& c = fs_slots[k];
            if (!c.active) return;
            double dP;
            if (apply_initial)
            {
                if (!(c.P > 0.0)) c.P = 0.0;
                dP = c.P;
            }
            else
            {
                double pdot, Wt, Pn;
                view.eval_fs(c, pdot, Wt);
                normal_projection_impulse(pdot, Wt, c.P, Pn, dP);
                c.P = Pn;
            }
            if (dP == 0.0) return;
            rhs::ContributionWriter<NormalDelta> fw(ws.fiber_buffers[seg], kFiberNormalFiberSolidFiberType);
            fiber_solid_normal_for_each_bead(c, [&](rhs::index_t bead, double w, double sg, std::uint16_t slot) {
                if (!(ws.bead_count[bead] > 0.0)) return;
                const double m = state.mass[bead];
                NormalDelta dl;
                for (std::size_t d = 0; d < dim; ++d) dl.dv[d] = (sg * w * dP) * c.n[d] / m;
                fw.emit(bead, dl, static_cast<std::uint64_t>(k), slot);
            });
            if (ws.solid_count[c.s] > 0.0)
            {
                rhs::ContributionWriter<NormalDelta> sw(ws.solid_buffers[seg], kFiberNormalFiberSolidSolidType);
                const double ms = sys.mass[c.s];
                NormalDelta dl;
                for (std::size_t d = 0; d < dim; ++d) dl.dv[d] = dP * c.n[d] / ms;
                sw.emit(c.s, dl, static_cast<std::uint64_t>(k), 0);
            }
        });

        // reduce вкладов и накопление Δv
        rhs::reduce_contributions_into(exec, ws.fiber_buffers, nfiber, ws.fiber_sum, ws.scratch);
        exec.parallel_for(nfiber, [&](std::size_t i) {
            for (std::size_t d = 0; d < dim; ++d) ws.dvel_fiber[d][i] += ws.fiber_sum[i].dv[d];
        });
        if (nfs > 0)
        {
            rhs::reduce_contributions_into(exec, ws.solid_buffers, nsolid, ws.solid_sum, ws.scratch);
            exec.parallel_for(nsolid, [&](std::size_t i) {
                for (std::size_t d = 0; d < dim; ++d) ws.dvel_solid[d][i] += ws.solid_sum[i].dv[d];
            });
        }
    };

    // --- активные слоты: если их нет, ничего не делаем (стадия не меняет dv)
    const std::size_t n_active_ff = exec.parallel_reduce(
        nff, std::size_t(0), [&](std::size_t k) -> std::size_t { return ff_slots[k].active ? 1u : 0u; },
        [](std::size_t a, std::size_t b) { return a + b; });
    const std::size_t n_active_fs = exec.parallel_reduce(
        nfs, std::size_t(0), [&](std::size_t k) -> std::size_t { return fs_slots[k].active ? 1u : 0u; },
        [](std::size_t a, std::size_t b) { return a + b; });
    stats.n_slots_ff = n_active_ff;
    stats.n_slots_fs = n_active_fs;
    if (n_active_ff + n_active_fs == 0) return;

    // Работа стадии W = sum m (v~.Δv + |Δv|^2/2) по текущим Δv (без вычитания близких величин).
    //
    // ПОЧЕМУ W <= 0 (ревью Codex, W_proj_pos_max > 0 в серии C). С v = v~ + M^-1 J^T P стадия минимизирует
    // f(P) = b^T P + P^T A P / 2, A = J M^-1 J^T, b = -pdot по пробной скорости; f(0) = 0. Без движущихся
    // опор f(P) = W (работа над свободными частицами). Шаг Якоби P <- max(0, P - D^-1 grad f), D = diag(W~_c), --
    // минимизация мажоранты f(P) + grad f^T d + d^T D d / 2, и он не увеличивает f, если D >= A. С разделением
    // массы это так: x^T A x = sum_k |sum_{c э k} sigma w x_c n_c|^2 / m_k <= sum_k c_k sum_{c э k} w^2 x_c^2 / m_k
    // = x^T D x (Коши-Буняковский по c_k слотам бусины; частицы Solid -- так же, закреплённые не входят ни туда,
    // ни туда). Значит, итерации f не наращивают, и f > 0 может прийти только из НАЧАЛЬНОЙ точки -- тёплого
    // старта P^(0) с f(P^(0)) > 0 (импульс прошлого шага у пары, которая уже расходится, или разгруженной
    // перестройкой соседей). Поэтому f(P^(0)) = A0 + B0/2 считается до итераций и при f > 0 весь тёплый старт
    // умножается на lambda = -2 A0/B0 (наибольший множитель с f(lambda P^(0)) <= 0; A0 >= 0 -- холодный старт).
    // Импульс сохраняется (Δv линейны по P), P >= 0 остаётся; при f(P^(0)) <= 0 (равновесие, цепочки) ничего
    // не меняется. ДВИЖУЩАЯСЯ опора (закреплённая с ненулевой заданной скоростью) может законно совершать
    // положительную работу W над свободными частицами: W = f - работа опоры, гарантия -- f <= 0.
    // Суммы -- по фиксированным отрезкам parallel_reduce: lambda от числа потоков не зависит.
    const auto stage_work = [&]() {
        const double W_fiber = exec.parallel_reduce(
            nfiber, 0.0,
            [&](std::size_t i) -> double {
                double s = 0.0;
                for (std::size_t d = 0; d < dim; ++d)
                {
                    const double a = ws.dvel_fiber[d][i];
                    if (a != 0.0) s += view.trial_fiber(d, static_cast<rhs::index_t>(i)) * a + 0.5 * a * a;
                }
                return state.mass[i] * s;
            },
            [](double a, double b) { return a + b; });
        const double W_solid = exec.parallel_reduce(
            nsolid, 0.0,
            [&](std::size_t i) -> double {
                double s = 0.0;
                for (std::size_t d = 0; d < dim; ++d)
                {
                    const double a = ws.dvel_solid[d][i];
                    if (a != 0.0) s += view.trial_solid(d, static_cast<rhs::index_t>(i)) * a + 0.5 * a * a;
                }
                return sys.mass[i] * s;
            },
            [](double a, double b) { return a + b; });
        return W_fiber + W_solid;
    };

    // 1. начальные импульсы (тёплый старт) и их энергетическая безопасность.
    //    g0 = f(P^(0)) = A0 + B0/2: A0 = sum_c P_c^(0) * (-pdot_c) по пробной скорости (Δv ещё нулевые),
    //    B0 = P^T A P = sum m |Δv0|^2 после прохода «применить». A0 считается по ОТНОСИТЕЛЬНЫМ скоростям
    //    контактов, а не по скоростям частиц: у движущейся опоры (закреплённая бусина с заданной скоростью)
    //    работа стадии над свободными частицами законно может быть > 0 (опора толкает), а убывает при
    //    итерациях именно f. Без движущихся опор f = W_stage.
    const double A0_ff = exec.parallel_reduce(
        nff, 0.0,
        [&](std::size_t k) -> double {
            const FiberNormalContact& c = ff_slots[k];
            if (!c.active || !(c.P > 0.0)) return 0.0;
            double pdot, Wt;
            view.eval_ff(c, pdot, Wt);
            return -pdot * c.P;
        },
        [](double x, double y) { return x + y; });
    const double A0_fs = exec.parallel_reduce(
        nfs, 0.0,
        [&](std::size_t k) -> double {
            const FiberSolidNormalContact& c = fs_slots[k];
            if (!c.active || !(c.P > 0.0)) return 0.0;
            double pdot, Wt;
            view.eval_fs(c, pdot, Wt);
            return -pdot * c.P;
        },
        [](double x, double y) { return x + y; });
    const double A0 = A0_ff + A0_fs;
    emit_pass(true);
    const auto dv_square = [&]() {
        const double Bf = exec.parallel_reduce(
            nfiber, 0.0,
            [&](std::size_t i) -> double {
                double q = 0.0;
                for (std::size_t d = 0; d < dim; ++d) q += ws.dvel_fiber[d][i] * ws.dvel_fiber[d][i];
                return state.mass[i] * q;
            },
            [](double x, double y) { return x + y; });
        const double Bs = exec.parallel_reduce(
            nsolid, 0.0,
            [&](std::size_t i) -> double {
                double q = 0.0;
                for (std::size_t d = 0; d < dim; ++d) q += ws.dvel_solid[d][i] * ws.dvel_solid[d][i];
                return sys.mass[i] * q;
            },
            [](double x, double y) { return x + y; });
        return Bf + Bs;
    };
    const double B0 = dv_square();
    stats.warm_scale = 1.0;
    stats.W_warm = A0 + 0.5 * B0; // без тёплого старта A0 = B0 = 0
    if (stats.W_warm > 0.0)
    {
        // Наибольший lambda с f(lambda P^(0)) = lambda A0 + lambda^2 B0/2 <= 0; A0 >= 0 -- холодный старт.
        const double lambda = (A0 < 0.0 && B0 > 0.0) ? std::min(1.0, -2.0 * A0 / B0) : 0.0;
        stats.warm_scale = lambda;
        exec.parallel_for_segments(nff, [&](std::size_t, std::size_t k) {
            if (ff_slots[k].active) ff_slots[k].P *= lambda;
        });
        exec.parallel_for_segments(nfs, [&](std::size_t, std::size_t k) {
            if (fs_slots[k].active) fs_slots[k].P *= lambda;
        });
        // Проход «применить» заново с lambda P^(0): та же формула вкладов, что и у итераций
        // (умножение уже накопленных Δv на lambda нарушало бы точную антисимметрию вкладов пары).
        for (std::size_t d = 0; d < dim; ++d)
        {
            ws.dvel_fiber[d].assign(nfiber, 0.0);
            ws.dvel_solid[d].assign(nsolid, 0.0);
        }
        emit_pass(true);
    }

    // 2. K итераций Якоби (каждая не увеличивает W, см. выше)
    for (int it = 0; it < iters; ++it) emit_pass(false);

    // 3. работа стадии за шаг
    stats.W_stage = stage_work();

    // 4. диагностика после последней итерации
    const detail_normal::SlotAgg agg_ff = exec.parallel_reduce(
        nff, detail_normal::SlotAgg{},
        [&](std::size_t k) -> detail_normal::SlotAgg {
            const FiberNormalContact& c = ff_slots[k];
            if (!c.active) return detail_normal::SlotAgg{};
            double pdot, Wt, Wc;
            view.eval_ff(c, pdot, Wt, &Wc);
            return detail_normal::agg_one(c.warm, pdot, Wt, Wc, c.P);
        },
        detail_normal::agg_combine);
    const detail_normal::SlotAgg agg_fs = exec.parallel_reduce(
        nfs, detail_normal::SlotAgg{},
        [&](std::size_t k) -> detail_normal::SlotAgg {
            const FiberSolidNormalContact& c = fs_slots[k];
            if (!c.active) return detail_normal::SlotAgg{};
            double pdot, Wt, Wc;
            view.eval_fs(c, pdot, Wt, &Wc);
            return detail_normal::agg_one(c.warm, pdot, Wt, Wc, c.P);
        },
        detail_normal::agg_combine);
    const detail_normal::SlotAgg agg = detail_normal::agg_combine(agg_ff, agg_fs);
    stats.resid = agg.resid;
    stats.resid_rel = agg.resid_rel;
    stats.resid_load = agg.sum_P > 0.0 ? agg.sum_missing / agg.sum_P : 0.0;
    stats.warm_frac = agg.active > 0 ? static_cast<double>(agg.warm) / static_cast<double>(agg.active) : 0.0;

    stats.max_ck = exec.parallel_reduce(
        nfiber, 0.0, [&](std::size_t i) { return ws.bead_count[i]; },
        [](double a, double b) { return std::max(a, b); });
    stats.max_cs = exec.parallel_reduce(
        nsolid, 0.0, [&](std::size_t i) { return ws.solid_count[i]; },
        [](double a, double b) { return std::max(a, b); });

    // 5. выход: ускорение Δv/h
    const double hf = fupd.h;
    const double hs = supd.h;
    exec.parallel_for(nfiber, [&](std::size_t i) {
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double a = ws.dvel_fiber[d][i];
            if (a != 0.0) dv_fiber[d][i] += a / hf;
        }
    });
    exec.parallel_for(nsolid, [&](std::size_t i) {
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double a = ws.dvel_solid[d][i];
            if (a != 0.0) dvx_solid[d][i] += a / hs;
        }
    });
}

// ----------------------------------------------------------------------
// normal_projection_post_residual -- остаток условия Синьорини по ИТОГОВЫМ
// скоростям шага (после стадии трения и повторной проекции зажимов, до
// advance): v = база + h*dv (+ av у Solid), dv_fiber / dvx_solid -- полные
// ускорения шага. Трение идёт после нормальной стадии и может дать нормальную
// составляющую в соседних контактах (ревью Codex), поэтому stats.resid
// (после K итераций) его не видит. Только диагностика: слоты и скорости не
// меняются. ws.bead_count / ws.solid_count -- те же, что у стадии этого шага.
// Заполняет stats.resid_post, stats.resid_load_post, stats.resid_post_valid.
// ----------------------------------------------------------------------
inline void normal_projection_post_residual(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                            const ode::FiberVelocityUpdate& fupd,
                                            const std::vector<double> (&dv_fiber)[dim],
                                            const std::vector<FiberNormalContact>& ff_slots, const ParticleSystem& sys,
                                            const ode::SolidVelocityUpdate& supd,
                                            const std::vector<std::vector<double>>& dvx_solid,
                                            const std::vector<FiberSolidNormalContact>& fs_slots,
                                            NormalProjectionWorkspace& ws, NormalProjectionStats& stats)
{
    const std::size_t nfiber = state.size();
    const std::size_t nsolid = ws.solid_count.size();
    for (std::size_t d = 0; d < dim; ++d)
    {
        ws.zero_fiber[d].assign(nfiber, 0.0);
        ws.zero_solid[d].assign(nsolid, 0.0);
    }
    const NormalStageView view{state, fupd,      dv_fiber,      ws.zero_fiber,        sys,
                               supd,  dvx_solid, ws.zero_solid, ws.bead_count.data(), ws.solid_count.data()};
    const detail_normal::SlotAgg agg_ff = exec.parallel_reduce(
        ff_slots.size(), detail_normal::SlotAgg{},
        [&](std::size_t k) -> detail_normal::SlotAgg {
            const FiberNormalContact& c = ff_slots[k];
            if (!c.active) return detail_normal::SlotAgg{};
            double pdot, Wt, Wc;
            view.eval_ff(c, pdot, Wt, &Wc);
            return detail_normal::agg_one(c.warm, pdot, Wt, Wc, c.P);
        },
        detail_normal::agg_combine);
    const detail_normal::SlotAgg agg_fs = exec.parallel_reduce(
        fs_slots.size(), detail_normal::SlotAgg{},
        [&](std::size_t k) -> detail_normal::SlotAgg {
            const FiberSolidNormalContact& c = fs_slots[k];
            if (!c.active) return detail_normal::SlotAgg{};
            double pdot, Wt, Wc;
            view.eval_fs(c, pdot, Wt, &Wc);
            return detail_normal::agg_one(c.warm, pdot, Wt, Wc, c.P);
        },
        detail_normal::agg_combine);
    const detail_normal::SlotAgg agg = detail_normal::agg_combine(agg_ff, agg_fs);
    stats.resid_post = agg.resid;
    stats.resid_load_post = agg.sum_P > 0.0 ? agg.sum_missing / agg.sum_P : 0.0;
    stats.resid_post_valid = true;
}

// ----------------------------------------------------------------------
// NormalWarmStart -- таблица тёплого старта: P прошлого шага по ключу пары.
// Хранится отсортированный (stable_sort по ключу) массив без повторов ключей
// (при повторе остаётся первая запись в порядке слотов) и h_prev.
//   apply   -- до стадии: P^(0) = P_prev*h/h_prev, warm = true, если пара найдена
//              (и enabled, h_prev > 0); иначе P = 0, warm = false. Неактивные
//              слоты обнуляются. Чтение таблицы параллельно безопасно.
//   rebuild -- после стадии: записи активных слотов с P > 0 (обход в порядке
//              слотов, нить-нить затем нить-Solid), h_prev = h.
//   entries()/h_prev()/set() -- экспорт/импорт для чекпоинта.
// ----------------------------------------------------------------------
class NormalWarmStart
{
public:
    void clear()
    {
        table_.clear();
        h_prev_ = 0.0;
    }

    std::size_t size() const { return table_.size(); }
    const std::vector<NormalWarmEntry>& entries() const { return table_; }
    double h_prev() const { return h_prev_; }

    // Импорт (чекпоинт): сортирует и убирает повторы, как rebuild.
    void set(std::vector<NormalWarmEntry> entries, double h_prev)
    {
        table_ = std::move(entries);
        normalize();
        h_prev_ = h_prev;
    }

    void apply(const rhs::CpuExecutor& exec, std::vector<FiberNormalContact>& ff_slots,
               std::vector<FiberSolidNormalContact>& fs_slots, double h, bool enabled) const
    {
        const bool use = enabled && h_prev_ > 0.0 && !table_.empty();
        const double scale = use ? h / h_prev_ : 0.0;
        exec.parallel_for_segments(ff_slots.size(), [&](std::size_t, std::size_t k) {
            FiberNormalContact& c = ff_slots[k];
            lookup(use, scale, c.active, c.key, c.P, c.warm);
        });
        exec.parallel_for_segments(fs_slots.size(), [&](std::size_t, std::size_t k) {
            FiberSolidNormalContact& c = fs_slots[k];
            lookup(use, scale, c.active, c.key, c.P, c.warm);
        });
    }

    void rebuild(const std::vector<FiberNormalContact>& ff_slots, const std::vector<FiberSolidNormalContact>& fs_slots,
                 double h)
    {
        table_.clear();
        for (const FiberNormalContact& c : ff_slots)
            if (c.active && c.P > 0.0) table_.push_back(NormalWarmEntry{c.key, c.P});
        for (const FiberSolidNormalContact& c : fs_slots)
            if (c.active && c.P > 0.0) table_.push_back(NormalWarmEntry{c.key, c.P});
        normalize();
        h_prev_ = h;
    }

private:
    void normalize()
    {
        std::stable_sort(table_.begin(), table_.end(),
                         [](const NormalWarmEntry& a, const NormalWarmEntry& b) { return a.key < b.key; });
        table_.erase(std::unique(table_.begin(), table_.end(),
                                 [](const NormalWarmEntry& a, const NormalWarmEntry& b) { return a.key == b.key; }),
                     table_.end());
    }

    void lookup(bool use, double scale, bool active, const NormalPairKey& key, double& P, bool& warm) const
    {
        P = 0.0;
        warm = false;
        if (!use || !active) return;
        const auto it = std::lower_bound(table_.begin(), table_.end(), key,
                                         [](const NormalWarmEntry& e, const NormalPairKey& k) { return e.key < k; });
        if (it == table_.end() || !(it->key == key)) return;
        P = it->P * scale;
        warm = true;
    }

    std::vector<NormalWarmEntry> table_;
    double h_prev_ = 0.0;
};

// ----------------------------------------------------------------------
// Слоты трения из нормальных (после стадии): слот k активен, если нормальный
// слот активен, P > 0 и mu > 0; геометрия копируется, N = P/h. Неактивные слоты
// сбрасываются в значение по умолчанию. Размер выходного вектора приводится к
// размеру нормального (индекс k сохраняется). Возвращает true, если активен хотя
// бы один слот (детерминированно: фиксированные отрезки parallel_reduce).
// Два вида заполняются НЕЗАВИСИМО (fill_friction_slots_from_normal_ff / _fs):
// если проекция включена только для одного вида, слоты трения другого вида
// (их заполняет ядро нормальной силы) не трогаются.
// ----------------------------------------------------------------------
inline bool fill_friction_slots_from_normal_ff(const rhs::CpuExecutor& exec, double h,
                                                const std::vector<FiberNormalContact>& ff_normal,
                                                std::vector<FiberFrictionContact>& ff_friction)
{
    ff_friction.resize(ff_normal.size());
    exec.parallel_for_segments(ff_normal.size(), [&](std::size_t, std::size_t k) {
        const FiberNormalContact& c = ff_normal[k];
        FiberFrictionContact f;
        if (c.active && c.P > 0.0 && c.mu > 0.0)
        {
            f.active = true;
            f.a0 = c.a0;
            f.a1 = c.a1;
            f.b0 = c.b0;
            f.b1 = c.b1;
            f.s = c.s;
            f.t = c.t;
            for (std::size_t d = 0; d < dim; ++d) f.n[d] = c.n[d];
            f.N = c.P / h;
            f.mu = c.mu;
        }
        ff_friction[k] = f;
    });
    const std::size_t nff = exec.parallel_reduce(
        ff_friction.size(), std::size_t(0),
        [&](std::size_t k) -> std::size_t { return ff_friction[k].active ? 1u : 0u; },
        [](std::size_t a, std::size_t b) { return a + b; });
    return nff > 0;
}

inline bool fill_friction_slots_from_normal_fs(const rhs::CpuExecutor& exec, double h,
                                                const std::vector<FiberSolidNormalContact>& fs_normal,
                                                std::vector<FiberSolidFrictionContact>& fs_friction)
{
    fs_friction.resize(fs_normal.size());
    exec.parallel_for_segments(fs_normal.size(), [&](std::size_t, std::size_t k) {
        const FiberSolidNormalContact& c = fs_normal[k];
        FiberSolidFrictionContact f;
        if (c.active && c.P > 0.0 && c.mu > 0.0)
        {
            f.active = true;
            f.a = c.a;
            f.b = c.b;
            f.s = c.s;
            f.lambda = c.lambda;
            for (std::size_t d = 0; d < dim; ++d) f.n[d] = c.n[d];
            f.N = c.P / h;
            f.mu = c.mu;
        }
        fs_friction[k] = f;
    });
    const std::size_t nfs = exec.parallel_reduce(
        fs_friction.size(), std::size_t(0),
        [&](std::size_t k) -> std::size_t { return fs_friction[k].active ? 1u : 0u; },
        [](std::size_t a, std::size_t b) { return a + b; });
    return nfs > 0;
}

// Оба вида сразу (оба выходных вектора перезаписываются).
inline bool fill_friction_slots_from_normal(const rhs::CpuExecutor& exec, double h,
                                             const std::vector<FiberNormalContact>& ff_normal,
                                             const std::vector<FiberSolidNormalContact>& fs_normal,
                                             std::vector<FiberFrictionContact>& ff_friction,
                                             std::vector<FiberSolidFrictionContact>& fs_friction)
{
    const bool a = fill_friction_slots_from_normal_ff(exec, h, ff_normal, ff_friction);
    const bool b = fill_friction_slots_from_normal_fs(exec, h, fs_normal, fs_friction);
    return a || b;
}

} // namespace interactions
