#pragma once
// interactions/fiber_friction.hpp
// Единая стадия трения нити как ПРОЕКЦИЯ НА УРОВНЕ СКОРОСТЕЙ (ограничение
// импульса трения) для ДВУХ видов контактов: нить-нить (fiber/
// fiber_friction_projection.md, задача 3 из 4) и нить-Solid (fiber/
// fiber_solid_friction_implementation.md, задача 4 из 4). Заменяет
// разрывный закон Кулона F_t = -mu*N*v_t/|v_t| внутри RHS, который в явной
// схеме давал дребезг (fiber/friction_chattering_analysis.md).
//
// ----------------------------------------------------------------------------
// МОДЕЛЬ НИТЬ-НИТЬ (§2 спецификации задачи 3). Для принятого контакта известны
// бусины a0,a1 и b0,b1 (a1/b1 == -1 для точечного сегмента), параметры s,t,
// нормаль n (от B к A), нормальная сила N>=0 и mu (слот FiberFrictionContact,
// его заполняет kernel нормального контакта, interactions/fiber_contact.hpp).
// Веса точки контакта: (1-s, s) на A и (1-t, t) на B.
//
//   Решатель обновит скорость бусины как v^{n+1} = base + h*a  (база и h --
//   ode::FiberVelocityUpdate, зависят от решателя). Пробная скорость -- БЕЗ
//   трения: v~_k = base_k + h*a_k^{без трения}. Относительная скорость точек
//   контакта u = sum_A alpha_k v~_k - sum_B beta_k v~_k, касательная часть
//   u_t = u - (u.n) n. W = 1/m_eff = sum_A alpha_k^2/m_k + sum_B beta_k^2/m_k
//   (fiber_contact_side_inv_mass, обе стороны).
//
//   При mu*N>0 и |u_t|>0:
//       J_stick = |u_t|/W,   J = min(mu*N*h, J_stick),   F_t = -(J/h) u_t/|u_t|.
//   Иначе F_t = 0. F_t действует на A (a_k = alpha_k F_t/m_k), -F_t на B
//   (a_k = -beta_k F_t/m_k). Бусина с нулевым весом не эмитируется.
//   После обновления Δu_t = h*W*F_t; при J=J_stick это ровно -u_t
//   (точное прилипание). Свойства: F_t.n=0, |F_t|<=mu*N, F_t.u_t<=0,
//   sum m_k a_k^f = 0.
//
// МОДЕЛЬ НИТЬ-SOLID (§3-§4 спецификации задачи 4). Слот FiberSolidFrictionContact
// (interactions/fiber_solid_contact.hpp): бусины a (вес 1-lambda), b (вес lambda,
// b == -1 для точечного контакта), частица Solid s, нормаль n ОТ НИТИ К SOLID,
// N, mu. Пробная скорость нити -- как выше; пробная скорость частицы Solid --
// ode::SolidVelocityUpdate: v~_s = (base + h*dvx) + av (av -- поправка XSPH,
// nullptr у Эйлера), ровно арифметика advance(). Относительная скорость (нить
// относительно Solid) u = [(1-lambda) v~_a + lambda v~_b] - v~_s, W =
// sum_{a,b} alpha_k^2/m_k + 1/m_s; проекция та же формула. F_t действует на
// НИТЬ (a: +(1-lambda)F_t, b: +lambda F_t), на Solid -F_t. Линейный импульс
// сохраняется точно; момент импульса -- нет (§6.1 спецификации: точки
// приложения разнесены по нормали, у частиц нет вращения).
//
// ЕДИНЫЙ ПРОХОД. Формулы проекции вынесены в friction_projection_impulse
// (принимает относительную скорость, нормаль, W, mu, N, h -- не знает вида
// контакта); fiber_friction_impulse и fiber_solid_friction_impulse лишь
// собирают u и W. Оба вида контактов читают ОДНО И ТО ЖЕ ускорение без трения
// (Якоби) и пишут вклады в общие буферы; два последовательных прохода не
// допускаются (бусина, касающаяся и нити, и пластины, получила бы двойную
// поправку).
//
// ПОЧЕМУ ПРОЕКЦИЯ. Явный Кулон перебрасывает скорость через ноль, если за шаг
// сила mu*N способна изменить относительную скорость больше её величины
// (двухшаговый цикл, трение не рассеивает, а иногда закачивает энергию), а
// PredictorCorrector добавляет смещение: знак берётся от предсказанной u*, а
// не от скорости, от которой считает корректор. Регуляризации (tanh, вязкий
// участок) не подходят: прилипания нет, результат зависит от dt. Проекция
// берёт РОВНО такой импульс, чтобы остановить касательное скольжение, но не
// больше mu*N*h: точное прилипание (в том числе под действием других сил --
// проекция видит полное ускорение шага без трения), скольжение с |F_t|=mu*N
// иначе, без перескока, память контакта не нужна. Поэтому ОДИН коэффициент
// (mu_ff, mu_fs) работает и как статический, и как кинематический.
//
// ПОРЯДОК СТАДИИ в шаге (src/time_integration.cpp):
//   predict -> RHS (single_step, Fiber-Bond, Fiber-Fiber нормаль + слоты,
//   Fiber-Solid нормаль + слоты) -> solver.fiber_rhs.dv, solver.solid_rhs.dvx
//   (всё, кроме трения) -> dt=min(CFL, контакт Solid, fiber_dt_limit) -> ЭТА
//   СТАДИЯ (fupd=solver.fiber_velocity_update(fb,dt), supd=solver.
//   solid_velocity_update(sys,dt); вклады -> reduce -> += в solver.fiber_rhs.dv
//   и solver.solid_rhs.dvx) -> advance. Скорость берётся от БАЗЫ решателя
//   (upd.base), а не из fb.vx/sys.vx: после предиктора PC там u*, и именно это
//   смещало цикл дребезга. Вклады пишутся через ContributionWriter -- порядок
//   обхода на результат не влияет, детерминизм и совместимость с CUDA как у
//   остальных interaction kernels. При mu_ff=mu_fs=0 у всех материалов
//   активных слотов нет, стадия НЕ выполняется, результат побитово прежний.
//   Трение не имеет собственной жёсткости и на dt не влияет.
//
// ----------------------------------------------------------------------------
// MASS SPLITTING (несколько контактов на бусину / частицу; §7 спецификации
// задачи 3, §5.2 задачи 4). Проекция по Якоби считает импульс каждого контакта
// так, будто остальные контакты на этом шаге не действуют; если бусина входит в
// k>1 активных контактов (узел сетки, плетение, пучок, два соседних сегмента
// одной нити, сегмент на нескольких частицах Solid), их суммарная поправка
// превышала бы нужную примерно в k раз -- перекомпенсация и перескок
// относительной скорости. Поэтому перед проходом считаются c_k -- число
// активных слотов ОБОИХ видов, в которых бусина k имеет НЕНУЛЕВОЙ вес, и c_s --
// число активных слотов Fiber-Solid с частицей s: compute_friction_counts,
// contributions + reduce с типами FiberCountDelta/SolidCountDelta, без atomics.
// В эффективной массе 1/m заменяется на c/m:
//   W = sum_A c_k alpha_k^2/m_k + sum_B c_k beta_k^2/m_k        (нить-нить),
//   W = sum_{a,b} c_k alpha_k^2/m_k + c_s/m_s                   (Fiber-Solid).
// Каждый контакт видит только свою долю массы. При c=1 (одиночные контакты)
// умножение на 1.0 точное -- результат побитово прежний. Векторы c живут у
// вызывающей стороны (переиспользуются между шагами); на CUDA-стороне это
// обычные массивы double[nbeads], double[nsolid].
//
// ИЗВЕСТНОЕ ОГРАНИЧЕНИЕ: недокомпенсирование. Контакт считает, что ему доступна
// лишь 1/c массы бусины (частицы), поэтому в узле с многими контактами
// прилипание может достигаться не за один шаг, а за несколько (медленная
// остаточная ползучесть); перескок при нескольких контактах на бусину
// устранён. Точное решение требует итераций Якоби/Гаусса-Зейделя (дальнейшая
// работа, не реализовано). Диагностика: shared_bead_count -- число бусин с
// c_k>1, shared_solid_count -- число частиц Solid с c_s>1.
// Сцена «нить на пластине» (связанные контакты, неравномерная u_t): рост |u_t|
// за шаг <= 0.21 % v_t0, перескок проекции >= -2.3 % v_t0; прилипание --
// релаксацией за конечное время. Второе ограничение: бусина и в контакте
// нить-нить, и в нить-Solid при РАЗНЫХ массах партнёров -- c по счёту
// контактов точен лишь при равных массах (перекомпенсация 0.18 %, без
// splitting 103 %). Третье: момент импульса не сохраняется (§6.1 спецификации
// задачи 4).
//
// КИНЕМАТИЧЕСКОЕ SOLID-ТЕЛО (только тесты, §5.4): флаг solid_kinematic=true --
// c_s/m_s=0 (бесконечная масса), вклады в Solid не эмитируются; скорость Solid
// -- та, что даёт переданный SolidVelocityUpdate (для неподвижной пластины --
// h=0 и base = её скорость). В боевом цикле флаг всегда false.
// Зажатым частицам Solid тот же эффект по частице даёт c_s = 0 (interactions/clamp.hpp).
//
// ПАРАЛЛЕЛЬНОСТЬ / CUDA (§10): kernel -- parallel_for по слотам фиксированного
// размера + reduce_contributions, без mutex/atomics и аллокаций внутри
// kernel. Формулы вынесены в friction_projection_impulse и
// *_friction_impulse (без CpuExecutor и исключений) -- кандидаты на
// __host__ __device__.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "interactions/fiber_contact.hpp"       // FiberFrictionContact, fiber_contact_side_inv_mass
#include "interactions/fiber_solid_contact.hpp" // FiberSolidFrictionContact
#include "interactions/solid_sph.hpp"           // SolidDelta
#include "ode/ode_solver.hpp"                   // ode::FiberVelocityUpdate, ode::SolidVelocityUpdate
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"

namespace interactions
{

// Метки interaction_type. Домен FiberDelta: 1 -- Fiber-Bond, 2 --
// Fiber-Fiber-Contact, 3 -- Fiber-Solid-Contact, 4 -- трение нить-нить,
// 5 -- трение Fiber-Solid (сторона нити). Домен SolidDelta: 10, 11, 12 заняты
// (Solid-SPH, Solid-Solid, Fiber-Solid нормаль), 13 -- трение Fiber-Solid
// (сторона Solid).
constexpr std::uint32_t kFiberFrictionInteractionType = 4;
constexpr std::uint32_t kFiberSolidFrictionFiberType = 5;
constexpr std::uint32_t kFiberSolidFrictionSolidType = 13;

// Результат формул одного контакта.
struct FiberFrictionImpulse
{
    double Ft[dim] = {};      // сила трения на сторону A (нить); на B (Solid) -- минус
    double ut_norm = 0.0;     // |u_t| пробной относительной касательной скорости
    double W = 0.0;           // 1/m_eff точки контакта
    double mu_N_h = 0.0;      // mu*N*h -- предел импульса
    double J_stick = 0.0;     // |u_t|/W -- импульс, останавливающий скольжение
    double J = 0.0;           // min(mu_N_h, J_stick); 0, если силы нет
    bool has_force = false;   // mu*N>0 && |u_t|>0
    bool stick = false;       // J_stick <= mu_N_h (выбрана ветка прилипания)
};

// ----------------------------------------------------------------------
// friction_projection_impulse -- формулы проекции для ОДНОГО контакта, НЕ
// привязанные к виду контакта: u -- пробная относительная скорость точек
// контакта (сторона A относительно B), n -- нормаль (знак не важен), W -- 1/m_eff
// (с учётом mass splitting), mu, N -- коэффициент трения и нормальная сила,
// h -- коэффициент при ускорении в обновлении скорости. Без исключений и
// аллокаций.
// ----------------------------------------------------------------------
inline void friction_projection_impulse(const double (&u)[dim], const double (&n)[dim], double W, double mu, double N,
                                         double h, FiberFrictionImpulse& out)
{
    out = FiberFrictionImpulse{};

    double un = 0.0;
    for (std::size_t d = 0; d < dim; ++d) un += u[d] * n[d];
    double ut[dim];
    double ut2 = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        ut[d] = u[d] - un * n[d];
        ut2 += ut[d] * ut[d];
    }
    out.ut_norm = std::sqrt(ut2);
    out.W = W;
    out.mu_N_h = mu * N * h;

    if (mu * N > 0.0 && out.ut_norm > 0.0)
    {
        out.has_force = true;
        out.J_stick = out.ut_norm / out.W;
        out.stick = (out.J_stick <= out.mu_N_h);
        out.J = out.stick ? out.J_stick : out.mu_N_h;
        const double scale = -(out.J / h) / out.ut_norm;
        for (std::size_t d = 0; d < dim; ++d) out.Ft[d] = scale * ut[d];
    }
}

// 1/m_eff точки контакта на одной стороне с mass splitting: fiber_contact_side_inv_mass
// с множителем c_k при каждой бусине; bead_count == nullptr -- все c_k=1. При
// c_k=1 умножение точное -- результат побитово как у fiber_contact_side_inv_mass.
// k1 == -1 -- точечный сегмент; s -- вес бусины k1 (вес k0 -- 1-s).
inline double friction_side_inv_mass(const fiber::FiberState& state, rhs::index_t k0, long long k1, double s,
                                      const double* bead_count)
{
    auto cnt = [&](rhs::index_t k) -> double { return bead_count ? bead_count[k] : 1.0; };
    if (k1 == -1) return cnt(k0) / state.mass[k0];
    const double w0 = 1.0 - s;
    const double w1 = s;
    const rhs::index_t kk1 = static_cast<rhs::index_t>(k1);
    return (cnt(k0) * (w0 * w0)) / state.mass[k0] + (cnt(kk1) * (w1 * w1)) / state.mass[kk1];
}

// ----------------------------------------------------------------------
// fiber_friction_impulse -- нить-нить: собирает u и W и вызывает
// friction_projection_impulse.
// bead_count -- c_k (число активных слотов с ненулевым весом бусины k, mass
// splitting, см. шапку); nullptr -- все c_k=1 (одиночный контакт).
// Скорость бусины k: (*upd.base[d])[k] + upd.h * dv_no_friction[d][k]
// (решатель обновит её именно так; dv_no_friction -- solver.fiber_rhs.dv ДО
// добавления трения). Скорость точки контакта стороны: (1-s)v_k0 + s*v_k1
// (k1==-1 -- v_k0).
// ----------------------------------------------------------------------
inline void fiber_friction_impulse(const fiber::FiberState& state, const ode::FiberVelocityUpdate& upd,
                                    const std::vector<double> (&dv_no_friction)[dim],
                                    const FiberFrictionContact& c, FiberFrictionImpulse& out,
                                    const double* bead_count = nullptr)
{
    auto trial_v = [&](std::size_t d, rhs::index_t k) -> double {
        return (*upd.base[d])[k] + upd.h * dv_no_friction[d][k];
    };

    double u[dim];
    for (std::size_t d = 0; d < dim; ++d)
    {
        double vA = trial_v(d, c.a0);
        if (c.a1 != -1) vA = (1.0 - c.s) * vA + c.s * trial_v(d, static_cast<rhs::index_t>(c.a1));
        double vB = trial_v(d, c.b0);
        if (c.b1 != -1) vB = (1.0 - c.t) * vB + c.t * trial_v(d, static_cast<rhs::index_t>(c.b1));
        u[d] = vA - vB;
    }

    const double W = friction_side_inv_mass(state, c.a0, c.a1, c.s, bead_count) +
                     friction_side_inv_mass(state, c.b0, c.b1, c.t, bead_count);
    friction_projection_impulse(u, c.n, W, c.mu, c.N, upd.h, out);
}

// ----------------------------------------------------------------------
// fiber_solid_friction_impulse -- нить-Solid: u = [(1-lambda) v~_a + lambda v~_b]
// - v~_s, W = side(a,b) + c_s/m_s. Пробная скорость частицы Solid:
// (base + h*dvx) + av (av -- XSPH, nullptr у Эйлера), ровно как advance().
// fiber_count/solid_count -- c_k, c_s (nullptr -- все 1). solid_kinematic:
// c_s/m_s = 0 (бесконечная масса Solid, только тесты).
// ----------------------------------------------------------------------
inline void fiber_solid_friction_impulse(const fiber::FiberState& state, const ode::FiberVelocityUpdate& fupd,
                                          const std::vector<double> (&dv_fiber_no_friction)[dim],
                                          const ParticleSystem& sys, const ode::SolidVelocityUpdate& supd,
                                          const std::vector<std::vector<double>>& dvx_solid_no_friction,
                                          const FiberSolidFrictionContact& c, FiberFrictionImpulse& out,
                                          const double* fiber_count = nullptr, const double* solid_count = nullptr,
                                          bool solid_kinematic = false)
{
    auto trial_vf = [&](std::size_t d, rhs::index_t k) -> double {
        return (*fupd.base[d])[k] + fupd.h * dv_fiber_no_friction[d][k];
    };
    auto trial_vs = [&](std::size_t d) -> double {
        double v = (*supd.base[d])[c.s] + supd.h * dvx_solid_no_friction[d][c.s];
        if (supd.av) v += (*supd.av)[d][c.s];
        return v;
    };

    double u[dim];
    for (std::size_t d = 0; d < dim; ++d)
    {
        double vF = trial_vf(d, c.a);
        if (c.b != -1) vF = (1.0 - c.lambda) * vF + c.lambda * trial_vf(d, static_cast<rhs::index_t>(c.b));
        u[d] = vF - trial_vs(d);
    }

    const double cs = solid_count ? solid_count[c.s] : 1.0;
    const double solid_inv_mass = solid_kinematic ? 0.0 : cs / sys.mass[c.s];
    const double W = friction_side_inv_mass(state, c.a, c.b, c.lambda, fiber_count) + solid_inv_mass;
    friction_projection_impulse(u, c.n, W, c.mu, c.N, fupd.h, out);
}

// ----------------------------------------------------------------------
// *_for_each_bead -- бусины слота с НЕНУЛЕВЫМ весом: f(bead, weight,
// side_sign, slot). Нить-нить: slot 0..3 = a0,a1,b0,b1 (ContributionMeta).
// Fiber-Solid: slot 0..1 = a,b (сторона нити, sign +1). Единственное место,
// где определено «ненулевой вес»: им пользуются и подсчёт c_k, и эмиссия
// вкладов.
// ----------------------------------------------------------------------
template <class F>
inline void fiber_friction_for_each_bead(const FiberFrictionContact& c, F&& f)
{
    if (c.s < 1.0) f(c.a0, 1.0 - c.s, 1.0, static_cast<std::uint16_t>(0));
    if (c.a1 != -1 && c.s > 0.0) f(static_cast<rhs::index_t>(c.a1), c.s, 1.0, static_cast<std::uint16_t>(1));
    if (c.t < 1.0) f(c.b0, 1.0 - c.t, -1.0, static_cast<std::uint16_t>(2));
    if (c.b1 != -1 && c.t > 0.0) f(static_cast<rhs::index_t>(c.b1), c.t, -1.0, static_cast<std::uint16_t>(3));
}

template <class F>
inline void fiber_solid_friction_for_each_bead(const FiberSolidFrictionContact& c, F&& f)
{
    if (c.lambda < 1.0) f(c.a, 1.0 - c.lambda, 1.0, static_cast<std::uint16_t>(0));
    if (c.b != -1 && c.lambda > 0.0) f(static_cast<rhs::index_t>(c.b), c.lambda, 1.0, static_cast<std::uint16_t>(1));
}

// Маленькие типы Delta для подсчёта контактов на бусину / частицу через
// contributions + reduce (без atomics).
struct FiberCountDelta
{
    double c = 0.0;
    FiberCountDelta& operator+=(const FiberCountDelta& o)
    {
        c += o.c;
        return *this;
    }
};

struct SolidCountDelta
{
    double c = 0.0;
    SolidCountDelta& operator+=(const SolidCountDelta& o)
    {
        c += o.c;
        return *this;
    }
};

// ----------------------------------------------------------------------
// Подсчёт c_k / c_s. Внутренние эмиттеры не чистят буферы (см. публичные
// функции ниже).
// ----------------------------------------------------------------------
namespace detail_friction
{

template <class Delta>
inline void prepare(const rhs::CpuExecutor& exec, std::vector<rhs::ContributionBuffer<Delta>>& buffers)
{
    const std::size_t nworkers = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    if (buffers.size() < nworkers) buffers.resize(nworkers);
    for (auto& buf : buffers) buf.clear();
}

template <class Delta>
inline std::size_t worker_slot(const std::vector<rhs::ContributionBuffer<Delta>>& buffers)
{
    int w = rhs::CpuExecutor::worker_id();
    if (w < 0 || static_cast<std::size_t>(w) >= buffers.size()) w = 0;
    return static_cast<std::size_t>(w);
}

inline void emit_counts_fiber_fiber(const rhs::CpuExecutor& exec, const std::vector<FiberFrictionContact>& contacts,
                                     std::vector<rhs::ContributionBuffer<FiberCountDelta>>& count_buffers)
{
    const std::size_t n = contacts.size();
    if (n == 0) return;
    exec.parallel_for(n, [&](std::size_t k) {
        const FiberFrictionContact& c = contacts[k];
        if (!c.active) return;
        rhs::ContributionWriter<FiberCountDelta> writer(count_buffers[worker_slot(count_buffers)],
                                                          kFiberFrictionInteractionType);
        FiberCountDelta one;
        one.c = 1.0;
        fiber_friction_for_each_bead(c, [&](rhs::index_t bead, double, double, std::uint16_t slot) {
            writer.emit(bead, one, static_cast<std::uint64_t>(k), slot);
        });
    });
}

inline void emit_counts_fiber_solid(const rhs::CpuExecutor& exec, const std::vector<FiberSolidFrictionContact>& contacts,
                                     std::vector<rhs::ContributionBuffer<FiberCountDelta>>& fiber_count_buffers,
                                     std::vector<rhs::ContributionBuffer<SolidCountDelta>>& solid_count_buffers)
{
    const std::size_t n = contacts.size();
    if (n == 0) return;
    exec.parallel_for(n, [&](std::size_t k) {
        const FiberSolidFrictionContact& c = contacts[k];
        if (!c.active) return;
        rhs::ContributionWriter<FiberCountDelta> fw(fiber_count_buffers[worker_slot(fiber_count_buffers)],
                                                      kFiberSolidFrictionFiberType);
        rhs::ContributionWriter<SolidCountDelta> sw(solid_count_buffers[worker_slot(solid_count_buffers)],
                                                      kFiberSolidFrictionSolidType);
        FiberCountDelta one;
        one.c = 1.0;
        fiber_solid_friction_for_each_bead(c, [&](rhs::index_t bead, double, double, std::uint16_t slot) {
            fw.emit(bead, one, static_cast<std::uint64_t>(k), slot);
        });
        SolidCountDelta ones;
        ones.c = 1.0;
        sw.emit(c.s, ones, static_cast<std::uint64_t>(k), 0);
    });
}

} // namespace detail_friction

// ----------------------------------------------------------------------
// compute_fiber_friction_bead_counts -- c_k для mass splitting ТОЛЬКО по
// контактам нить-нить: число активных слотов, в которых бусина k имеет
// ненулевой вес. bead_count переразмеряется до nbeads (ёмкость только растёт),
// count_buffers -- переиспользуемые per-thread буферы вызывающей стороны.
// Вызывается ПЕРЕД проходом стадии трения (после заполнения слотов).
// ----------------------------------------------------------------------
inline void compute_fiber_friction_bead_counts(const rhs::CpuExecutor& exec, std::size_t nbeads,
                                                const std::vector<FiberFrictionContact>& contacts,
                                                std::vector<rhs::ContributionBuffer<FiberCountDelta>>& count_buffers,
                                                std::vector<double>& bead_count)
{
    detail_friction::prepare(exec, count_buffers);
    bead_count.assign(nbeads, 0.0);
    if (contacts.empty()) return;

    detail_friction::emit_counts_fiber_fiber(exec, contacts, count_buffers);
    const std::vector<FiberCountDelta> counts = rhs::reduce_contributions(exec, count_buffers, nbeads);
    for (std::size_t i = 0; i < nbeads; ++i) bead_count[i] = counts[i].c;
}

// ----------------------------------------------------------------------
// compute_friction_counts -- c_k и c_s по контактам ОБОИХ видов (§5.2
// спецификации задачи 4): c_k -- число активных слотов нить-нить и Fiber-Solid,
// в которых бусина k имеет ненулевой вес; c_s -- число активных слотов
// Fiber-Solid с частицей s. bead_count размера nfiber, solid_count размера
// nsolid (реальные частицы [0, ntotal)); буферы переиспользуются. Результаты
// reduce (fiber_counts, solid_counts) и рабочая память reduce (scratch) --
// тоже переменные вызывающей стороны (память только растёт).
// ----------------------------------------------------------------------
inline void compute_friction_counts(const rhs::CpuExecutor& exec, std::size_t nfiber, std::size_t nsolid,
                                     const std::vector<FiberFrictionContact>& fiber_contacts,
                                     const std::vector<FiberSolidFrictionContact>& solid_contacts,
                                     std::vector<rhs::ContributionBuffer<FiberCountDelta>>& fiber_count_buffers,
                                     std::vector<rhs::ContributionBuffer<SolidCountDelta>>& solid_count_buffers,
                                     std::vector<double>& bead_count, std::vector<double>& solid_count,
                                     std::vector<FiberCountDelta>& fiber_counts,
                                     std::vector<SolidCountDelta>& solid_counts, rhs::ReduceScratch& scratch)
{
    detail_friction::prepare(exec, fiber_count_buffers);
    detail_friction::prepare(exec, solid_count_buffers);
    bead_count.assign(nfiber, 0.0);
    solid_count.assign(nsolid, 0.0);

    detail_friction::emit_counts_fiber_fiber(exec, fiber_contacts, fiber_count_buffers);
    detail_friction::emit_counts_fiber_solid(exec, solid_contacts, fiber_count_buffers, solid_count_buffers);

    if (!fiber_contacts.empty() || !solid_contacts.empty())
    {
        rhs::reduce_contributions_into(exec, fiber_count_buffers, nfiber, fiber_counts, scratch);
        for (std::size_t i = 0; i < nfiber; ++i) bead_count[i] = fiber_counts[i].c;
    }
    if (!solid_contacts.empty())
    {
        rhs::reduce_contributions_into(exec, solid_count_buffers, nsolid, solid_counts, scratch);
        for (std::size_t i = 0; i < nsolid; ++i) solid_count[i] = solid_counts[i].c;
    }
}

// Обёртка со старой сигнатурой: локальные результаты reduce и рабочая память.
inline void compute_friction_counts(const rhs::CpuExecutor& exec, std::size_t nfiber, std::size_t nsolid,
                                     const std::vector<FiberFrictionContact>& fiber_contacts,
                                     const std::vector<FiberSolidFrictionContact>& solid_contacts,
                                     std::vector<rhs::ContributionBuffer<FiberCountDelta>>& fiber_count_buffers,
                                     std::vector<rhs::ContributionBuffer<SolidCountDelta>>& solid_count_buffers,
                                     std::vector<double>& bead_count, std::vector<double>& solid_count)
{
    std::vector<FiberCountDelta> fiber_counts;
    std::vector<SolidCountDelta> solid_counts;
    rhs::ReduceScratch scratch;
    compute_friction_counts(exec, nfiber, nsolid, fiber_contacts, solid_contacts, fiber_count_buffers,
                            solid_count_buffers, bead_count, solid_count, fiber_counts, solid_counts, scratch);
}

// ----------------------------------------------------------------------
// Эмиссия вкладов a^f (не чистит буферы -- см. публичные функции ниже).
// ----------------------------------------------------------------------
namespace detail_friction
{

inline void emit_fiber_fiber(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                              const ode::FiberVelocityUpdate& upd, const std::vector<double> (&dv_no_friction)[dim],
                              const std::vector<FiberFrictionContact>& contacts, const std::vector<double>& bead_count,
                              std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers)
{
    const std::size_t n = contacts.size();
    if (n == 0) return;

    exec.parallel_for(n, [&](std::size_t k) {
        const FiberFrictionContact& c = contacts[k];
        if (!c.active) return;

        FiberFrictionImpulse imp;
        fiber_friction_impulse(state, upd, dv_no_friction, c, imp, bead_count.data());
        if (!imp.has_force) return;

        rhs::ContributionWriter<fiber::FiberDelta> writer(thread_buffers[worker_slot(thread_buffers)],
                                                            kFiberFrictionInteractionType);

        // Сторона A: +F_t, веса (1-s, s); сторона B: -F_t, веса (1-t, t).
        fiber_friction_for_each_bead(c, [&](rhs::index_t bead, double weight, double sign, std::uint16_t slot) {
            const double m = state.mass[bead];
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = sign * weight * imp.Ft[dd] / m;
            writer.emit(bead, d, static_cast<std::uint64_t>(k), slot);
        });
    });
}

inline void emit_fiber_solid(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                              const ode::FiberVelocityUpdate& fupd, const std::vector<double> (&dv_fiber)[dim],
                              const ParticleSystem& sys, const ode::SolidVelocityUpdate& supd,
                              const std::vector<std::vector<double>>& dvx_solid,
                              const std::vector<FiberSolidFrictionContact>& contacts,
                              const std::vector<double>& bead_count, const std::vector<double>& solid_count,
                              bool solid_kinematic,
                              std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& fiber_buffers,
                              std::vector<rhs::ContributionBuffer<SolidDelta>>& solid_buffers)
{
    const std::size_t n = contacts.size();
    if (n == 0) return;

    exec.parallel_for(n, [&](std::size_t k) {
        const FiberSolidFrictionContact& c = contacts[k];
        if (!c.active) return;

        FiberFrictionImpulse imp;
        fiber_solid_friction_impulse(state, fupd, dv_fiber, sys, supd, dvx_solid, c, imp, bead_count.data(),
                                      solid_count.data(), solid_kinematic);
        if (!imp.has_force) return;

        rhs::ContributionWriter<fiber::FiberDelta> fw(fiber_buffers[worker_slot(fiber_buffers)],
                                                        kFiberSolidFrictionFiberType);

        // Нить: +F_t, веса (1-lambda, lambda); Solid: -F_t.
        fiber_solid_friction_for_each_bead(c, [&](rhs::index_t bead, double weight, double sign, std::uint16_t slot) {
            const double m = state.mass[bead];
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = sign * weight * imp.Ft[dd] / m;
            fw.emit(bead, d, static_cast<std::uint64_t>(k), slot);
        });

        if (!solid_kinematic)
        {
            rhs::ContributionWriter<SolidDelta> sw(solid_buffers[worker_slot(solid_buffers)],
                                                     kFiberSolidFrictionSolidType);
            const double ms = sys.mass[c.s];
            SolidDelta ds{};
            for (std::size_t dd = 0; dd < dim; ++dd) ds.dv[dd] = -imp.Ft[dd] / ms;
            sw.emit(c.s, ds, static_cast<std::uint64_t>(k), 0);
        }
    });
}

} // namespace detail_friction

// ----------------------------------------------------------------------
// compute_fiber_friction_contributions -- kernel стадии трения ТОЛЬКО для
// нить-нить: parallel_for по слотам, неактивный пропускается, вклады a_k^f
// эмитируются по слотам ContributionMeta 0..3 (a0,a1,b0,b1); бусины с нулевым
// весом (s==1 -> a0, s==0 -> a1, аналогично для b) не эмитируются.
// Суммирование -- rhs::reduce_contributions вызывающей стороной. dv_no_friction
// только читается. bead_count -- c_k (размер >= state.size()).
// ----------------------------------------------------------------------
inline void compute_fiber_friction_contributions(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                                  const ode::FiberVelocityUpdate& upd,
                                                  const std::vector<double> (&dv_no_friction)[dim],
                                                  const std::vector<FiberFrictionContact>& contacts,
                                                  const std::vector<double>& bead_count,
                                                  std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers)
{
    detail_friction::prepare(exec, thread_buffers);
    detail_friction::emit_fiber_fiber(exec, state, upd, dv_no_friction, contacts, bead_count, thread_buffers);
}

// ----------------------------------------------------------------------
// compute_friction_contributions -- ЕДИНЫЙ проход стадии трения по контактам
// обоих видов (§5 спецификации задачи 4). Оба вида читают одни и те же
// dv_fiber/dvx_solid (без трения, Якоби). Вклады нити (нить-нить: type 4,
// Fiber-Solid: type 5) -- в fiber_buffers; вклады Solid (type 13) -- в
// solid_buffers (при solid_kinematic=true не эмитируются). Суммирование --
// rhs::reduce_contributions вызывающей стороной.
// ----------------------------------------------------------------------
inline void compute_friction_contributions(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const ode::FiberVelocityUpdate& fupd,
    const std::vector<double> (&dv_fiber)[dim], const std::vector<FiberFrictionContact>& fiber_contacts,
    const ParticleSystem& sys, const ode::SolidVelocityUpdate& supd,
    const std::vector<std::vector<double>>& dvx_solid, const std::vector<FiberSolidFrictionContact>& solid_contacts,
    const std::vector<double>& bead_count, const std::vector<double>& solid_count, bool solid_kinematic,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& fiber_buffers,
    std::vector<rhs::ContributionBuffer<SolidDelta>>& solid_buffers)
{
    detail_friction::prepare(exec, fiber_buffers);
    detail_friction::prepare(exec, solid_buffers);
    detail_friction::emit_fiber_fiber(exec, state, fupd, dv_fiber, fiber_contacts, bead_count, fiber_buffers);
    detail_friction::emit_fiber_solid(exec, state, fupd, dv_fiber, sys, supd, dvx_solid, solid_contacts, bead_count,
                                       solid_count, solid_kinematic, fiber_buffers, solid_buffers);
}

// ----------------------------------------------------------------------
// Диагностика (§6). Не вызывается в горячем цикле; для тестов и демо, по
// аналогии с compute_fiber_contact_diagnostics.
//   active_count      -- активные слоты (contact && mu>0);
//   stick_count       -- J == J_stick среди контактов с ненулевым F_t;
//   slip_count        -- J == mu*N*h среди контактов с ненулевым F_t
//                        (active_count-stick-slip -- контакты с F_t=0:
//                        N=0 либо |u_t|=0);
//   shared_bead_count -- бусины (с ненулевым весом), входящие в >1 активный
//                        контакт, c_k>1 (см. «MASS SPLITTING» в шапке);
//                        диагностика использует те же c_k (bead_count), что и стадия;
//   shared_solid_count -- частицы Solid с c_s>1 (только compute_friction_
//                        diagnostics);
//   max_force_ratio   -- max |F_t|/(mu*N), должно быть <= 1.
// ----------------------------------------------------------------------
struct FiberFrictionDiagnostics
{
    std::size_t active_count = 0;
    std::size_t stick_count = 0;
    std::size_t slip_count = 0;
    std::size_t shared_bead_count = 0;
    std::size_t shared_solid_count = 0;
    double max_force_ratio = 0.0;
};

namespace detail_friction
{

inline FiberFrictionDiagnostics diag_combine(FiberFrictionDiagnostics a, const FiberFrictionDiagnostics& b)
{
    a.active_count += b.active_count;
    a.stick_count += b.stick_count;
    a.slip_count += b.slip_count;
    a.max_force_ratio = std::max(a.max_force_ratio, b.max_force_ratio);
    return a;
}

inline FiberFrictionDiagnostics diag_one(const FiberFrictionImpulse& imp, double mu, double N)
{
    FiberFrictionDiagnostics d1;
    d1.active_count = 1;
    if (!imp.has_force) return d1;
    if (imp.stick) d1.stick_count = 1;
    else d1.slip_count = 1;
    double f2 = 0.0;
    for (std::size_t d = 0; d < dim; ++d) f2 += imp.Ft[d] * imp.Ft[d];
    d1.max_force_ratio = std::sqrt(f2) / (mu * N);
    return d1;
}

inline FiberFrictionDiagnostics diag_fiber_fiber(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                                  const ode::FiberVelocityUpdate& upd,
                                                  const std::vector<double> (&dv_no_friction)[dim],
                                                  const std::vector<FiberFrictionContact>& contacts,
                                                  const std::vector<double>& bead_count)
{
    const std::size_t n = contacts.size();
    if (n == 0) return FiberFrictionDiagnostics{};
    return exec.parallel_reduce(
        n, FiberFrictionDiagnostics{},
        [&](std::size_t k) -> FiberFrictionDiagnostics {
            const FiberFrictionContact& c = contacts[k];
            if (!c.active) return FiberFrictionDiagnostics{};
            FiberFrictionImpulse imp;
            fiber_friction_impulse(state, upd, dv_no_friction, c, imp, bead_count.data());
            return diag_one(imp, c.mu, c.N);
        },
        diag_combine);
}

} // namespace detail_friction

inline FiberFrictionDiagnostics compute_fiber_friction_diagnostics(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const ode::FiberVelocityUpdate& upd,
    const std::vector<double> (&dv_no_friction)[dim], const std::vector<FiberFrictionContact>& contacts,
    const std::vector<double>& bead_count)
{
    FiberFrictionDiagnostics result =
        detail_friction::diag_fiber_fiber(exec, state, upd, dv_no_friction, contacts, bead_count);

    // shared_bead_count: бусины с c_k>1 (>1.5 с запасом на округление).
    for (std::size_t i = 0; i < state.size() && i < bead_count.size(); ++i)
        if (bead_count[i] > 1.5) ++result.shared_bead_count;

    return result;
}

// Диагностика единой стадии (оба вида контактов); c_k/c_s -- из
// compute_friction_counts.
inline FiberFrictionDiagnostics compute_friction_diagnostics(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const ode::FiberVelocityUpdate& fupd,
    const std::vector<double> (&dv_fiber)[dim], const std::vector<FiberFrictionContact>& fiber_contacts,
    const ParticleSystem& sys, const ode::SolidVelocityUpdate& supd,
    const std::vector<std::vector<double>>& dvx_solid, const std::vector<FiberSolidFrictionContact>& solid_contacts,
    const std::vector<double>& bead_count, const std::vector<double>& solid_count, bool solid_kinematic)
{
    FiberFrictionDiagnostics result =
        detail_friction::diag_fiber_fiber(exec, state, fupd, dv_fiber, fiber_contacts, bead_count);

    const std::size_t n = solid_contacts.size();
    if (n > 0)
    {
        const FiberFrictionDiagnostics fs = exec.parallel_reduce(
            n, FiberFrictionDiagnostics{},
            [&](std::size_t k) -> FiberFrictionDiagnostics {
                const FiberSolidFrictionContact& c = solid_contacts[k];
                if (!c.active) return FiberFrictionDiagnostics{};
                FiberFrictionImpulse imp;
                fiber_solid_friction_impulse(state, fupd, dv_fiber, sys, supd, dvx_solid, c, imp, bead_count.data(),
                                              solid_count.data(), solid_kinematic);
                return detail_friction::diag_one(imp, c.mu, c.N);
            },
            detail_friction::diag_combine);
        result = detail_friction::diag_combine(result, fs);
    }

    for (std::size_t i = 0; i < state.size() && i < bead_count.size(); ++i)
        if (bead_count[i] > 1.5) ++result.shared_bead_count;
    for (std::size_t i = 0; i < solid_count.size(); ++i)
        if (solid_count[i] > 1.5) ++result.shared_solid_count;

    return result;
}

} // namespace interactions
