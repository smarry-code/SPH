#pragma once
// interactions/fiber_contact.hpp
// Fiber-Fiber-Contact: пространственный механический контакт бусин РАЗНЫХ
// нитей (fiber/SPH.md, раздел «Fiber-Fiber-Contact»; fiber/PLAN.md, этап Ф5,
// «3.2. Типы взаимодействий и primitives», «3.3», «3.7»).
//
// Как и Fiber-Bond (interactions/fiber_bond.hpp) и Solid-Solid-Contact
// (interactions/solid_contact.hpp), kernel только ЧИТАЕТ FiberState/
// FiberTopology/FiberMaterial и производит typed contributions (FiberDelta)
// через rhs::ContributionWriter — суммирование отдельно, в
// rhs::reduce_contributions. Никаких atomics/mutex.
//
// ============================================================================
// РАЗРЕШЁННЫЕ НЕОДНОЗНАЧНОСТИ (SPH.md формулирует физику, но оставляет знаки
// и обобщение на разные материалы без явного вывода — см. также LaTeX-дефект
// «\$$6pt]» вместо «\\[6pt]» в формулах F^c/w_ij, чисто вёрсточный, физику не
// затрагивает):
//
// 1. Знак P_ij и ориентация p_i/p_j. SPH.md задаёт матричную систему
//    [p_i, n, -p_j]*[P_ij,D_ij,P_ji]^T = x_ij ОДНОВременно с проекционной
//    формулой Surface-End P_ij = -p_i . x_ij; при x_ij = x_i - x_j (принятая
//    здесь конвенция, стандартная для SPH: r_ij = x_i - x_j) эти две формулы
//    дают P_ij с ПРОТИВОПОЛОЖНЫМИ знаками друг относительно друга (проверено
//    прямым выводом: см. также fiber/PLAN.md, задание Ф5 — «сделай
//    самостоятельный вывод знака из требования непрерывности веса»). Здесь
//    используется САМОСОГЛАСОВАННОЕ определение, выведенное из требований:
//      - P_ij>0 <=> точка контакта на оси i смещена в сторону i_prev
//        (непрерывность w_ij на стыке P_ij=0: обе ветви дают w=1);
//      - n_ij направлена "от нити j к нити i" (проверено на End-End: там
//        единственная однозначная формула n=x_ij/|x_ij| = (x_i-x_j)/|..| —
//        сила на i вдоль +n отталкивает i от j, физически корректно).
//    Отсюда: p_i определяется как ЕДИНИЧНЫЙ вектор normalize(x[i_prev] -
//    x[i_next]) ("назад" по нити, в сторону prev); P_ij — параметр ближайшей
//    точки на оси i, Q_i = x_i + P_ij*p_i, найденный проекцией/решением
//    системы «ближайшая точка между двумя прямыми»; n_ij = (Q_i-подобная
//    точка на i) - (Q_j-подобная точка на j), нормированная. Все 4 случая
//    (Surface-Surface/Surface-End/End-Surface/End-End) реализуют это ЕДИНОЕ
//    правило — см. compute_fiber_contact_geometry ниже, где вывод расписан
//    по шагам для каждого случая.
//
// 2. d/E*/R* при разных материалах нитей i и j. SPH.md даёт формулу Meyer
//    только для ОДИНАКОВЫХ нитей (E* = E/(2(1-nu^2)), R* = sqrt(d/4)), но в
//    разделе «Fiber-Solid-Contact», «Постановка проекта: коэффициент K_fs»
//    для ПАРЫ РАЗНЫХ материалов (Fiber-Solid) уже выведена корректная
//    Hertz-обобщённая формула:
//      1/E* = (1-nu_a^2)/E_a + (1-nu_b^2)/E_b,   1/R* = 1/r_a + 1/r_b.
//    Та же логика (эффективные модуль/кривизна пары контактирующих тел из
//    контактной механики Герца) применена здесь к паре нитей i,j: E_i,nu_i
//    из FiberMaterial нити i, E_j,nu_j — нити j; r_i=d_i/2, r_j=d_j/2.
//    Порог проникновения обобщён СИММЕТРИЧНО тем же способом, каким в
//    Fiber-Solid-Contact используется delta_fs = r_f+r_s-D_fs (а не
//    сравнение с «диаметром одной стороны»): d_eff = r_i+r_j,
//    delta = d_eff - D_ij, контакт при D_ij < d_eff.
//    САМОПРОВЕРКА (одинаковые нити, E_i=E_j=E, nu_i=nu_j=nu, d_i=d_j=d):
//      1/E* = 2(1-nu^2)/E => E* = E/(2(1-nu^2))                    — совпадает с Meyer;
//      r_i=r_j=d/2 => d_eff=d (порог D<d, как у Meyer);
//      1/R* = 1/(d/2)+1/(d/2) = 4/d => R*=d/4, sqrt(R*)=sqrt(d/4)  — совпадает с Meyer.
//    (Проверено также юнит-тестом test_fiber_contact.cpp: perpendicular
//    одинаковых нитей.)
//
// 3. Радиус поиска d+l_max. l_max и d берутся как ГЛОБАЛЬНЫЕ максимумы по
//    ВСЕМ материалам/связям расчёта (не по каждой паре отдельно) — SpatialIndex
//    строится с ОДНИМ фиксированным радиусом на весь набор Fiber-частиц
//    (rhs/spatial_index.hpp: build(coords, n, cell_size) принимает единственное
//    число), поэтому радиус обязан быть консервативной ВЕРХНЕЙ оценкой для
//    любой пары, а не точным значением для каждой — недостаточный радиус для
//    одной пары означал бы пропущенный контакт, избыточный для другой — лишь
//    более широкую (но корректную) фазу поиска, точный критерий контакта
//    (D_ij<d_eff) всё равно проверяется в kernel. l_max — максимум l0 по
//    topology.bonds (номинальная, недеформированная длина связи; PLAN.md Ф5
//    явно допускает этот вариант).
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "param.hpp"
#include "physics/fiber_contact_geom.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/spatial_index.hpp"
#include "rhs/types.hpp"

namespace interactions
{

// Метка interaction_type для ContributionMeta/детерминированного reduce
// (rhs/contributions.hpp, «Детерминизм»); отличается от kFiberBondInteractionType.
constexpr std::uint32_t kFiberFiberContactInteractionType = 2;

// Порог |p_i x p_j| < eps ("почти параллельные нити") и защита от деления на
// ~0 при нормировке n -- fiber/CUDA_PLAN.md §3.2, этап Г2 (часть "Fiber-
// Fiber-Contact"): значения и обоснование (диапазон PLAN.md 1e-6..1e-3,
// устойчивость решения, test_fiber_contact.cpp "почти параллельные нити")
// перенесены БЕЗ ИЗМЕНЕНИЙ в physics::kFiberContactParallelEps/
// kFiberContactMinD (physics/fiber_contact_geom.hpp) вместе с самой
// геометрией (compute_fiber_contact_geometry ниже) -- используются только
// там, здесь оставлен только комментарий-указатель.

// ----------------------------------------------------------------------
// FiberContactCase, FiberContactPrimitive -- primitive Fiber-Fiber-Contact
// (fiber/PLAN.md §3.2, таблица; решение Р8: вид контакта определяется ПРИ
// ПОСТРОЕНИИ primitive по prev/next, kernel далее использует только
// сохранённых здесь соседей, topology повторно не читает).
// ----------------------------------------------------------------------
enum class FiberContactCase : std::uint8_t
{
    SurfaceSurface, // обе частицы внутренние (prev!=-1 && next!=-1)
    SurfaceEnd,     // i внутренняя, j -- конец нити
    EndSurface,     // i -- конец нити, j внутренняя
    EndEnd          // обе частицы -- концы своих нитей
};

struct FiberContactPrimitive
{
    rhs::index_t i = 0, j = 0;
    long long i_prev = -1, i_next = -1;
    long long j_prev = -1, j_next = -1;
    FiberContactCase case_ = FiberContactCase::EndEnd;
};

// ----------------------------------------------------------------------
// fiber_contact_search_radius -- d_max + l_max (см. пункт 3 в шапке файла).
// ----------------------------------------------------------------------
inline real_t fiber_contact_search_radius(const fiber::FiberTopology& topology,
                                           const std::vector<fiber::FiberMaterial>& materials)
{
    real_t d_max = 0;
    for (const auto& m : materials) d_max = std::max(d_max, m.d);
    real_t l_max = 0;
    for (const auto& b : topology.bonds) l_max = std::max(l_max, b.l0);
    return d_max + l_max;
}

// ----------------------------------------------------------------------
// build_fiber_contact_search -- SpatialIndex<dim> (dim==3, Р9) по бусинам
// Fiber, фильтр i<j (канонический, из SpatialIndex) И object_id[i]!=
// object_id[j] (разные нити -- самоконтакт одной нити не формулируется,
// SPH.md, «Fiber-Fiber-Contact», последний абзац). Для каждой найденной пары
// читает topology ОДИН раз и классифицирует case (Р8).
// ----------------------------------------------------------------------
struct FiberContactSearchResult
{
    std::vector<FiberContactPrimitive> primitives;
};

inline void build_fiber_contact_search(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                        const fiber::FiberTopology& topology,
                                        const std::vector<fiber::FiberMaterial>& materials,
                                        FiberContactSearchResult& out)
{
    out.primitives.clear();
    const std::size_t n = state.size();
    if (n == 0) return;

    if constexpr (dim != 3)
    {
        // Нити физически существуют только при dim==3 (fiber/SPH.md; Р9);
        // граница проверяется fiber::check_fiber_requires_dim3() при разборе
        // конфигурации -- сюда с непустым state в сборке dim!=3 попасть
        // невозможно. Заглушка нужна только для компиляции заголовка во ВСЕХ
        // SPH_DIM-вариантах пайплайна (см. include/interactions/fiber_bond.hpp,
        // include/fiber/fiber_system.hpp -- тот же паттерн).
        return;
    }
    else
    {
        const real_t radius = fiber_contact_search_radius(topology, materials);

        rhs::SpatialIndex<dim> index;
        std::array<const real_t*, dim> coords{};
        for (std::size_t d = 0; d < dim; ++d) coords[d] = state.x[d].data();
        index.build(coords, n, radius > 0 ? radius : real_t(1.0));

        auto accept = [&](rhs::index_t a, rhs::index_t b) -> bool {
            return state.object_id[a] != state.object_id[b];
        };

        const auto pairs = index.find_pairs_parallel(exec, radius, accept);

        out.primitives.reserve(pairs.size());
        for (const auto& pr : pairs)
        {
            FiberContactPrimitive prim;
            prim.i = pr.first;
            prim.j = pr.second;
            prim.i_prev = topology.prev[prim.i];
            prim.i_next = topology.next[prim.i];
            prim.j_prev = topology.prev[prim.j];
            prim.j_next = topology.next[prim.j];

            const bool i_end = (prim.i_prev == -1 || prim.i_next == -1);
            const bool j_end = (prim.j_prev == -1 || prim.j_next == -1);
            if (!i_end && !j_end) prim.case_ = FiberContactCase::SurfaceSurface;
            else if (!i_end && j_end) prim.case_ = FiberContactCase::SurfaceEnd;
            else if (i_end && !j_end) prim.case_ = FiberContactCase::EndSurface;
            else prim.case_ = FiberContactCase::EndEnd;

            out.primitives.push_back(prim);
        }
    }
}

// ----------------------------------------------------------------------
// FiberContactGeometry -- итог геометрического решения ОДНОГО primitive:
// нормаль/расстояние/сила Герца (если есть контакт) и веса/соседи для
// распределения силы (SPH.md, "w_ij"). Общая точка входа для
// compute_fiber_contact_contributions, compute_fiber_contact_dt_limit и
// диагностики -- геометрия вычисляется РОВНО ОДИН РАЗ на каждый вызывающий
// проход (аналог PairGeom/PairProps в contact_force.hpp), но не кэшируется
// между тремя проходами (SPH.md о кэшировании geometрии между проходами не
// требует; см. также interactions/solid_contact.hpp -- там геометрия тоже
// пересчитывается отдельно в diagnostics).
// ----------------------------------------------------------------------
struct FiberContactGeometry
{
    bool contact = false;
    real_t n[dim] = {};
    real_t D = 0;
    real_t F = 0;     // сила Герца, >=0
    real_t delta = 0; // проникновение d_eff - D_ij, >0 при контакте

    // P_ij/P_ji -- параметры ближайших точек вдоль p_i/p_j (см. шапку файла,
    // пункт 1); заполняются ВСЕГДА, когда геометрия найдена (даже без
    // контакта, D>=d_eff) -- удобно для тестов и диагностики. Остаются 0.0
    // для End-роли соответствующей стороны (см. side_weight ниже -- там P не
    // применяется).
    real_t P_ij = 0;
    real_t P_ji = 0;

    real_t w_i = 0;
    bool i_has_neighbor = false;
    rhs::index_t i_neighbor = 0;

    real_t w_j = 0;
    bool j_has_neighbor = false;
    rhs::index_t j_neighbor = 0;
};

// fiber/CUDA_PLAN.md §3.2, этап Г2 (часть "Fiber-Fiber-Contact"): вся
// геометрия (4 случая, вырожденная ветка, сила Герца, веса "w_ij") вынесена
// в physics::compute_fiber_contact_geometry (physics/fiber_contact_geom.hpp)
// -- SPH_HD-функция, общая для CPU и (в будущей части Г7) GPU. Здесь --
// только сборка входа из FiberState/FiberTopology/FiberMaterial (points,
// материалы) и раскладка выхода в FiberContactGeometry (в т.ч. перевод
// i_neighbor_is_prev/j_neighbor_is_prev в rhs::index_t по prim.i_prev/
// i_next/j_prev/j_next -- индексы знает только эта, topology-осведомлённая,
// сторона).
inline FiberContactGeometry compute_fiber_contact_geometry(const fiber::FiberState& state,
                                                             const std::vector<fiber::FiberMaterial>& materials,
                                                             const FiberContactPrimitive& prim)
{
    FiberContactGeometry geom{};
    if constexpr (dim != 3)
    {
        return geom; // см. build_fiber_contact_search -- недостижимо при реальном использовании
    }
    else
    {
        const rhs::index_t i = prim.i, j = prim.j;
        const fiber::FiberMaterial& mat_i = materials[static_cast<std::size_t>(state.material_id[i])];
        const fiber::FiberMaterial& mat_j = materials[static_cast<std::size_t>(state.material_id[j])];

        physics::FiberContactGeomInput in{};
        in.case_ = static_cast<int>(prim.case_); // та же нумерация 0..3, см. physics::FiberContactCase
        for (std::size_t d = 0; d < dim; ++d)
        {
            in.x_i[d] = state.x[d][i];
            in.x_j[d] = state.x[d][j];
        }
        // x_*_prev/x_*_next читаются physics-функцией ТОЛЬКО когда case_
        // отводит этой стороне Surface-роль (prev!=-1 && next!=-1
        // гарантировано классификацией case_, см. build_fiber_contact_search) --
        // заполняем безусловно там, где индекс валиден.
        if (prim.i_prev != -1)
            for (std::size_t d = 0; d < dim; ++d) in.x_i_prev[d] = state.x[d][static_cast<rhs::index_t>(prim.i_prev)];
        if (prim.i_next != -1)
            for (std::size_t d = 0; d < dim; ++d) in.x_i_next[d] = state.x[d][static_cast<rhs::index_t>(prim.i_next)];
        if (prim.j_prev != -1)
            for (std::size_t d = 0; d < dim; ++d) in.x_j_prev[d] = state.x[d][static_cast<rhs::index_t>(prim.j_prev)];
        if (prim.j_next != -1)
            for (std::size_t d = 0; d < dim; ++d) in.x_j_next[d] = state.x[d][static_cast<rhs::index_t>(prim.j_next)];

        in.r_i = real_t(0.5) * mat_i.d;
        in.r_j = real_t(0.5) * mat_j.d;
        in.E_i = mat_i.E; in.nu_i = mat_i.nu;
        in.E_j = mat_j.E; in.nu_j = mat_j.nu;

        physics::FiberContactGeomOutput out{};
        physics::compute_fiber_contact_geometry(in, out);

        geom.contact = out.contact;
        for (std::size_t d = 0; d < dim; ++d) geom.n[d] = out.n[d];
        geom.D = out.D;
        geom.F = out.F;
        geom.delta = out.delta;
        geom.P_ij = out.P_ij;
        geom.P_ji = out.P_ji;

        geom.w_i = out.w_i;
        geom.i_has_neighbor = out.i_has_neighbor;
        if (out.i_has_neighbor)
            geom.i_neighbor = static_cast<rhs::index_t>(out.i_neighbor_is_prev ? prim.i_prev : prim.i_next);

        geom.w_j = out.w_j;
        geom.j_has_neighbor = out.j_has_neighbor;
        if (out.j_has_neighbor)
            geom.j_neighbor = static_cast<rhs::index_t>(out.j_neighbor_is_prev ? prim.j_prev : prim.j_next);

        return geom;
    }
}

// ----------------------------------------------------------------------
// compute_fiber_contact_contributions -- kernel: для каждого primitive с
// контактом (geom.contact) эмитирует вклады FiberDelta на i, j и (если есть)
// их соседей i_prev/i_next, j_prev/j_next (до 4 частиц на primitive, SPH.md
// §3.2). Момент сохраняется точно: сумма вкладов*mass по всем ДО ЧЕТЫРЁХ
// частиц одного primitive равна нулю (F*n на стороне i, -F*n на стороне j,
// см. test_fiber_contact.cpp, "сумма импульсов").
// ----------------------------------------------------------------------
inline void compute_fiber_contact_contributions(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const std::vector<fiber::FiberMaterial>& materials,
    const std::vector<FiberContactPrimitive>& primitives,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers)
{
    const std::size_t nworkers = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const std::size_t n = primitives.size();
    if (n == 0) return;

    exec.parallel_for(n, [&](std::size_t k) {
        const FiberContactPrimitive& prim = primitives[k];
        const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
        if (!geom.contact) return;

        int w = rhs::CpuExecutor::worker_id();
        if (w < 0 || static_cast<std::size_t>(w) >= thread_buffers.size()) w = 0;
        rhs::ContributionWriter<fiber::FiberDelta> writer(thread_buffers[static_cast<std::size_t>(w)],
                                                            kFiberFiberContactInteractionType);

        const real_t mi = state.mass[prim.i];
        const real_t mj = state.mass[prim.j];

        fiber::FiberDelta di{}, dj{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            di.dv[d] = geom.w_i * geom.F * geom.n[d] / mi;
            dj.dv[d] = -geom.w_j * geom.F * geom.n[d] / mj;
        }
        writer.emit(prim.i, di, static_cast<std::uint64_t>(k), 0);
        writer.emit(prim.j, dj, static_cast<std::uint64_t>(k), 2);

        if (geom.i_has_neighbor)
        {
            const real_t m_nb = state.mass[geom.i_neighbor];
            const real_t rem = real_t(1.0) - geom.w_i;
            fiber::FiberDelta dnb{};
            for (std::size_t d = 0; d < dim; ++d) dnb.dv[d] = rem * geom.F * geom.n[d] / m_nb;
            writer.emit(geom.i_neighbor, dnb, static_cast<std::uint64_t>(k), 1);
        }
        if (geom.j_has_neighbor)
        {
            const real_t m_nb = state.mass[geom.j_neighbor];
            const real_t rem = real_t(1.0) - geom.w_j;
            fiber::FiberDelta dnb{};
            for (std::size_t d = 0; d < dim; ++d) dnb.dv[d] = -rem * geom.F * geom.n[d] / m_nb;
            writer.emit(geom.j_neighbor, dnb, static_cast<std::uint64_t>(k), 3);
        }
    });
}

// compute_fiber_contact -- удобная обёртка: search + contributions + reduce.
inline std::vector<fiber::FiberDelta> compute_fiber_contact(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const fiber::FiberTopology& topology,
    const std::vector<fiber::FiberMaterial>& materials,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers)
{
    FiberContactSearchResult search;
    build_fiber_contact_search(exec, state, topology, materials, search);
    compute_fiber_contact_contributions(exec, state, materials, search.primitives, thread_buffers);
    return rhs::reduce_contributions(exec, thread_buffers, state.size());
}

// ----------------------------------------------------------------------
// compute_fiber_contact_dt_limit -- fiber/PLAN.md §3.7:
//     dt_lim = 0.8/sqrt(k_eff/m_red),  k_eff = 1.5*F/delta,
//     m_red = m_i*m_j/(m_i+m_j)  (бусины i,j -- "якоря" primitive, по
//     аналогии с Fiber-Bond/Solid-Solid-Contact; соседи, получающие часть
//     силы через w_ij, в m_red не участвуют -- как и там).
// Считается ТОЛЬКО при фактическом контакте (geom.contact, т.е. D<d_eff) --
// В ОТЛИЧИЕ от Fiber-Bond (где dt_limit считается для ВСЕХ связей независимо
// от знака деформации, т.к. связь существует всегда, а её жёсткость k=EA/l0
// не зависит от текущего удлинения). Здесь контакта в принципе МОЖЕТ не быть
// (D_ij>=d_eff -- k_eff не определена, силы нет), поэтому естественный
// аналог -- Solid-Solid-Contact (dt_limit только по РЕАЛЬНО контактирующим
// парам, contact_force.hpp/ADR-007), а не Fiber-Bond. +inf, если контактов
// нет -- вызывающая сторона берёт min(...) безусловно.
// ----------------------------------------------------------------------
inline real_t compute_fiber_contact_dt_limit(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                              const std::vector<fiber::FiberMaterial>& materials,
                                              const std::vector<FiberContactPrimitive>& primitives)
{
    const real_t kInf = std::numeric_limits<real_t>::max();
    const std::size_t n = primitives.size();
    if (n == 0) return kInf;

    return exec.parallel_reduce(
        n, kInf,
        [&](std::size_t k) -> real_t {
            const FiberContactPrimitive& prim = primitives[k];
            const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
            if (!geom.contact || !(geom.delta > 0)) return kInf;

            const real_t mi = state.mass[prim.i];
            const real_t mj = state.mass[prim.j];
            const real_t m_red = (mi * mj) / (mi + mj);
            const real_t k_eff = real_t(1.5) * geom.F / geom.delta;
            if (!(k_eff > 0)) return kInf;
            return real_t(0.8) / std::sqrt(k_eff / m_red);
        },
        [](real_t a, real_t b) { return std::min(a, b); });
}

// ----------------------------------------------------------------------
// Диагностика (fiber/PLAN.md, этап Ф5: «диагностика: число пар на одно
// пересечение»). Физическое пересечение (или сближение) двух нитей, как
// правило, порождает НЕСКОЛЬКО spatial-пар бусин одновременно (соседние
// бусины обеих нитей одновременно попадают в радиус поиска d+l_max друг у
// друга) -- это распознанное свойство модели (fiber-design-decisions),
// НЕ баг: contact_count в норме кратно превышает число геометрических
// "пересечений" нитей, что и проверяет test_fiber_contact.cpp
// ("диагностика: число пар/пересечение", разумный небольшой множитель, а не
// 1:1 и не абсурдно много).
//   pair_count    -- все spatial-пары после фильтра object_id[i]!=object_id[j]
//                     (до проверки D_ij<d_eff -- широкая фаза);
//   contact_count -- подмножество с фактическим проникновением (D_ij<d_eff);
//   max_penetration, min_dt -- как в SolidContactDiagnostics
//                     (interactions/solid_contact.hpp), для единообразия.
// ----------------------------------------------------------------------
struct FiberContactDiagnostics
{
    std::size_t pair_count = 0;
    std::size_t contact_count = 0;
    real_t max_penetration = 0;
    real_t min_dt = std::numeric_limits<real_t>::max();
};

inline FiberContactDiagnostics compute_fiber_contact_diagnostics(const rhs::CpuExecutor& exec,
                                                                   const fiber::FiberState& state,
                                                                   const std::vector<fiber::FiberMaterial>& materials,
                                                                   const std::vector<FiberContactPrimitive>& primitives)
{
    FiberContactDiagnostics zero;
    zero.pair_count = primitives.size();
    const std::size_t n = primitives.size();
    if (n == 0) return zero;

    auto combine = [](FiberContactDiagnostics a, FiberContactDiagnostics b) {
        FiberContactDiagnostics r;
        r.pair_count = a.pair_count + b.pair_count;
        r.contact_count = a.contact_count + b.contact_count;
        r.max_penetration = std::max(a.max_penetration, b.max_penetration);
        r.min_dt = std::min(a.min_dt, b.min_dt);
        return r;
    };

    FiberContactDiagnostics init;
    init.pair_count = 0; // складывается через combine ниже по одному на primitive
    FiberContactDiagnostics acc = exec.parallel_reduce(
        n, init,
        [&](std::size_t k) -> FiberContactDiagnostics {
            FiberContactDiagnostics d1;
            d1.pair_count = 1;
            const FiberContactPrimitive& prim = primitives[k];
            const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
            if (geom.contact)
            {
                d1.contact_count = 1;
                d1.max_penetration = geom.delta;
                const real_t mi = state.mass[prim.i];
                const real_t mj = state.mass[prim.j];
                const real_t m_red = (mi * mj) / (mi + mj);
                const real_t k_eff = (geom.delta > 0) ? real_t(1.5) * geom.F / geom.delta : real_t(0.0);
                d1.min_dt = (k_eff > 0) ? real_t(0.8) / std::sqrt(k_eff / m_red) : std::numeric_limits<real_t>::max();
            }
            return d1;
        },
        combine);
    return acc;
}

} // namespace interactions
