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
//
// 4. Трение. Модель — динамическое сухое трение Кулона БЕЗ памяти контакта
//    (fiber/fiber_friction_implementation.md, целиком). Скорость точки
//    контакта на каждой стороне строится ИЗ ТЕХ ЖЕ весов/соседей, что уже
//    использует распределение нормальной силы (w_i/i_has_neighbor/
//    i_neighbor и симметрично для j, см. FiberContactGeometry и
//    fiber_contact_point_velocity ниже): v_c = w*v_bead + (1-w)*v_neighbor;
//    для End-стороны (has_neighbor=false, w=1) это вырождается в v_end —
//    отдельной формулы не вводится. В вырожденной почти-параллельной ветке
//    Surface-Surface (см. пункт 1) j трактуется как точка (w_j=1, соседа
//    нет) — та же формула fiber_contact_point_velocity даёт v_c=v_j без
//    дополнительного разбора случаев. Итоговая сила F = N*n + F_t
//    распределяется по бусинам ОДНИМИ И ТЕМИ ЖЕ весами, что и N*n (§5-6
//    спецификации) — никакого отдельного распределения для трения.
//    Правило смешения mu для разных материалов вынесено в отдельную функцию
//    fiber_friction_coefficient (§4/§8 спецификации) — менять только её.
//    dt_limit (compute_fiber_contact_dt_limit) трением не затрагивается: у
//    закона Кулона нет собственной жёсткости (в отличие от нормального
//    закона Герца), характерное время трения не вводит более жёсткого
//    ограничения, чем уже посчитанное по нормальному контакту.
//    Одно геометрическое пересечение нитей обычно даёт НЕСКОЛЬКО contact
//    primitives (см. диагностику ниже) — суммарная сила трения по всем
//    primitives одного пересечения равна mu*(сумма их N), что согласовано с
//    суммарной нормальной силой тем же способом (никакой отдельной
//    нормировки трения на число primitives не вводится).
//    Известный риск явной схемы интегрирования: при v_t, близкой к нулю,
//    сила трения на каждом шаге пересчитывается заново по текущему знаку
//    v_t (модель без памяти, без stick/slip) — это может давать дребезг
//    (осцилляцию) v_t около нуля с амплитудой порядка mu*N*dt/m за шаг. В
//    этой итерации это НЕ регуляризуется (нет сглаживания или регуляризации
//    около v_t=0 сверх порога eps_v, нет тангенциальной пружины) —
//    статическое трение и stick/slip state machine — отдельная задача
//    (см. §13 спецификации). mu_ff=0 (значение по умолчанию в
//    FiberMaterial) обязан давать ПОБИТОВО прежнее поведение — это
//    единственный жёсткий инвариант регрессии, проверяемый
//    test_friction_zero_mu_regression.
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
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/spatial_index.hpp"
#include "rhs/types.hpp"

namespace interactions
{

// Метка interaction_type для ContributionMeta/детерминированного reduce
// (rhs/contributions.hpp, «Детерминизм»); отличается от kFiberBondInteractionType.
constexpr std::uint32_t kFiberFiberContactInteractionType = 2;

// Порог |p_i x p_j| < eps -- нити считаются "почти параллельными" (Surface-
// Surface система [p_i n -p_j] вырождена, знаменатель 2x2-решения ~ eps^2).
// eps -- безразмерный (|p_i x p_j| = sin(угол между p_i,p_j) для единичных
// векторов). Выбрано ближе к консервативному концу диапазона, предложенного
// PLAN.md (1e-6..1e-3): при eps=1e-3 знаменатель решения >= 1e-6, что на
// много порядков выше плавающего эпсилон (~1e-16) -- решение остаётся
// численно устойчивым (см. разбор в верхнем комментарии файла, пункт 1),
// P_ij/P_ji при cr_norm вблизи eps могут получаться большими, но это не
// приводит к NaN/Inf: итоговая проверка D_ij>=d_eff (нет контакта) либо
// weight-формула (branch "иначе -> 0") гасят большие |P| корректно
// (test_fiber_contact.cpp, "почти параллельные нити").
constexpr double kFiberContactParallelEps = 1e-3;

// Защита от деления на ~0 при нормировке n (совпадающие точки/оси).
constexpr double kFiberContactMinD = 1e-12;

// Порог |v_t| в м/с для закона трения Кулона (fiber/fiber_friction_implementation.md,
// §4/§10): нужен ТОЛЬКО как защита от деления на ноль при нормировке v_t и
// от численного шума около v_t=0 -- не физический параметр модели. Выбран
// малым по сравнению с характерными скоростями задач проекта (десятые доли
// -- сотни м/с).
constexpr double kFiberContactFrictionMinVt = 1e-8;

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
inline double fiber_contact_search_radius(const fiber::FiberTopology& topology,
                                           const std::vector<fiber::FiberMaterial>& materials)
{
    double d_max = 0.0;
    for (const auto& m : materials) d_max = std::max(d_max, m.d);
    double l_max = 0.0;
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
        const double radius = fiber_contact_search_radius(topology, materials);

        rhs::SpatialIndex<dim> index;
        std::array<const double*, dim> coords{};
        for (std::size_t d = 0; d < dim; ++d) coords[d] = state.x[d].data();
        index.build(coords, n, radius > 0.0 ? radius : 1.0);

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
    double n[dim] = {};
    double D = 0.0;
    double F = 0.0;     // сила Герца, >=0
    double delta = 0.0; // проникновение d_eff - D_ij, >0 при контакте

    // P_ij/P_ji -- параметры ближайших точек вдоль p_i/p_j (см. шапку файла,
    // пункт 1); заполняются ВСЕГДА, когда геометрия найдена (даже без
    // контакта, D>=d_eff) -- удобно для тестов и диагностики. Остаются 0.0
    // для End-роли соответствующей стороны (см. side_weight ниже -- там P не
    // применяется).
    double P_ij = 0.0;
    double P_ji = 0.0;

    double w_i = 0.0;
    bool i_has_neighbor = false;
    rhs::index_t i_neighbor = 0;

    double w_j = 0.0;
    bool j_has_neighbor = false;
    rhs::index_t j_neighbor = 0;
};

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
        using Vec3 = std::array<double, 3>;
        auto sub = [](const Vec3& a, const Vec3& b) -> Vec3 { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
        auto dot = [](const Vec3& a, const Vec3& b) -> double { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
        auto norm = [&](const Vec3& a) -> double { return std::sqrt(dot(a, a)); };
        auto cross = [](const Vec3& a, const Vec3& b) -> Vec3 {
            return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
        };
        auto point = [&](rhs::index_t k) -> Vec3 { return {state.x[0][k], state.x[1][k], state.x[2][k]}; };
        // "Назад" по нити -- normalize(x[prev]-x[next]); знак/ориентация
        // выведены и задокументированы в шапке файла, пункт 1.
        auto direction_backward = [&](long long prev, long long next) -> Vec3 {
            const Vec3 xp = point(static_cast<rhs::index_t>(prev));
            const Vec3 xn = point(static_cast<rhs::index_t>(next));
            Vec3 v = sub(xp, xn);
            const double len = norm(v);
            if (len > 0.0) { v[0] /= len; v[1] /= len; v[2] /= len; }
            return v;
        };

        const rhs::index_t i = prim.i, j = prim.j;
        const Vec3 xi = point(i);
        const Vec3 xj = point(j);
        const Vec3 r = sub(xi, xj); // x_ij := x_i - x_j (см. шапку файла)

        double P_ij = 0.0, P_ji = 0.0;
        bool i_surface = false; // используется ли Surface-формула (P_ij) для i
        bool j_surface = false;
        Vec3 n{};
        double D = 0.0;

        if (prim.case_ == FiberContactCase::EndEnd)
        {
            const Vec3 diff = r;
            D = norm(diff);
            if (D < kFiberContactMinD) return geom; // совпадающие бусины -- пропуск, защита от NaN
            n = {diff[0] / D, diff[1] / D, diff[2] / D};
        }
        else if (prim.case_ == FiberContactCase::SurfaceEnd)
        {
            const Vec3 p_i = direction_backward(prim.i_prev, prim.i_next);
            P_ij = dot(p_i, sub(xj, xi));
            const Vec3 Q_i = {xi[0] + P_ij * p_i[0], xi[1] + P_ij * p_i[1], xi[2] + P_ij * p_i[2]};
            const Vec3 diff = sub(Q_i, xj);
            D = norm(diff);
            if (D < kFiberContactMinD) return geom;
            n = {diff[0] / D, diff[1] / D, diff[2] / D};
            i_surface = true;
        }
        else if (prim.case_ == FiberContactCase::EndSurface)
        {
            const Vec3 p_j = direction_backward(prim.j_prev, prim.j_next);
            P_ji = dot(p_j, sub(xi, xj));
            const Vec3 Q_j = {xj[0] + P_ji * p_j[0], xj[1] + P_ji * p_j[1], xj[2] + P_ji * p_j[2]};
            const Vec3 diff = sub(xi, Q_j);
            D = norm(diff);
            if (D < kFiberContactMinD) return geom;
            n = {diff[0] / D, diff[1] / D, diff[2] / D};
            j_surface = true;
        }
        else // SurfaceSurface
        {
            const Vec3 p_i = direction_backward(prim.i_prev, prim.i_next);
            const Vec3 p_j = direction_backward(prim.j_prev, prim.j_next);
            const Vec3 cr = cross(p_i, p_j);
            const double cr_norm = norm(cr);

            if (cr_norm < kFiberContactParallelEps)
            {
                // Нити почти параллельны -- система [p_i n -p_j] вырождена
                // (fiber/PLAN.md Ф5, fiber/SPH.md «Fiber-Fiber-Contact»):
                // переход на формулу Surface-End. Детерминированный выбор --
                // i остаётся "поверхностной" стороной (i<j по построению
                // поиска), j трактуется как точка: w_j=1, вклад j НЕ
                // распределяется на его соседей в этой редкой вырожденной
                // ветке (сознательное упрощение, см. шапку файла).
                P_ij = dot(p_i, sub(xj, xi));
                const Vec3 Q_i = {xi[0] + P_ij * p_i[0], xi[1] + P_ij * p_i[1], xi[2] + P_ij * p_i[2]};
                const Vec3 diff = sub(Q_i, xj);
                D = norm(diff);
                if (D < kFiberContactMinD) return geom;
                n = {diff[0] / D, diff[1] / D, diff[2] / D};
                i_surface = true;
            }
            else
            {
                // Ближайшие точки двух прямых (2x2 система из условия
                // перпендикулярности Q_i-Q_j обеим p_i,p_j -- эквивалентна
                // решению [p_i n -p_j]*[P_ij,D_ij,P_ji]^T=x_ij относительно
                // P_ij,P_ji при уже известной n=+-(p_i x p_j)/|..|):
                const double c = dot(p_i, p_j);
                const double denom = 1.0 - c * c; // == cr_norm^2 (>= eps^2 здесь)
                const double rpi = dot(r, p_i);
                const double rpj = dot(r, p_j);
                P_ji = (rpj - c * rpi) / denom;
                P_ij = c * P_ji - rpi;
                const Vec3 Q_i = {xi[0] + P_ij * p_i[0], xi[1] + P_ij * p_i[1], xi[2] + P_ij * p_i[2]};
                const Vec3 Q_j = {xj[0] + P_ji * p_j[0], xj[1] + P_ji * p_j[1], xj[2] + P_ji * p_j[2]};
                const Vec3 diff = sub(Q_i, Q_j);
                D = norm(diff);
                if (D < kFiberContactMinD)
                {
                    // Оси фактически пересекаются (D==0) -- нормаль берём из
                    // p_i x p_j, ориентированную "от j к i" (dot(n,r)>=0).
                    n = {cr[0] / cr_norm, cr[1] / cr_norm, cr[2] / cr_norm};
                    if (dot(n, r) < 0.0) { n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2]; }
                    D = 0.0;
                }
                else
                {
                    n = {diff[0] / D, diff[1] / D, diff[2] / D};
                }
                i_surface = true;
                j_surface = true;
            }
        }

        // P_ij/P_ji заполняются ВСЕГДА, когда геометрия дошла до этой точки
        // (даже если ниже окажется, что контакта нет, D>=d_eff) -- см.
        // комментарий у полей FiberContactGeometry.
        geom.P_ij = P_ij;
        geom.P_ji = P_ji;

        // --- Hertz-сила (см. пункт 2 в шапке файла) -------------------------
        const fiber::FiberMaterial& mat_i = materials[static_cast<std::size_t>(state.material_id[i])];
        const fiber::FiberMaterial& mat_j = materials[static_cast<std::size_t>(state.material_id[j])];
        const double ri = 0.5 * mat_i.d;
        const double rj = 0.5 * mat_j.d;
        const double d_eff = ri + rj;
        if (!(D < d_eff)) return geom; // нет проникновения

        const double delta = d_eff - D;
        const double inv_Estar = (mat_i.E > 0.0 ? (1.0 - mat_i.nu * mat_i.nu) / mat_i.E : 0.0) +
                                  (mat_j.E > 0.0 ? (1.0 - mat_j.nu * mat_j.nu) / mat_j.E : 0.0);
        if (!(inv_Estar > 0.0)) return geom; // некорректный материал (E<=0 с обеих сторон) -- защита
        const double E_star = 1.0 / inv_Estar;
        const double R_star = (ri > 0.0 && rj > 0.0) ? (ri * rj) / (ri + rj) : 0.0;
        const double F = (4.0 / 3.0) * E_star * std::sqrt(R_star) * std::pow(delta, 1.5);

        geom.contact = true;
        geom.D = D;
        geom.delta = delta;
        geom.F = F;
        geom.n[0] = n[0];
        geom.n[1] = n[1];
        geom.n[2] = n[2];

        // --- Веса и соседи (SPH.md, "w_ij") ----------------------------------
        // surface_role==false (End) -> w=1, соседа нет (SPH.md: контакт с
        // концом нити прикладывается целиком к концевой бусине).
        auto side_weight = [&](bool surface_role, double P, long long prev, long long next, const Vec3& xs,
                                double& w_out, bool& has_nb, rhs::index_t& nb_out) {
            if (!surface_role)
            {
                w_out = 1.0;
                has_nb = false;
                return;
            }
            // surface_role==true подразумевает Surface-роль по классификации
            // primitive (Р8) => частица внутренняя => prev!=-1 && next!=-1
            // ГАРАНТИРОВАННО (build_fiber_contact_search формирует case
            // именно из этого условия) -- static_cast ниже безопасен.
            const double d_p = norm(sub(xs, point(static_cast<rhs::index_t>(prev))));
            const double d_n = norm(sub(xs, point(static_cast<rhs::index_t>(next))));
            double w;
            if (P > 0.0 && P < d_p) w = (d_p - P) / d_p;
            else if (P <= 0.0 && P > -d_n) w = (d_n + P) / d_n;
            else w = 0.0; // контакт проецируется дальше ближайшего соседа -- см. шапку файла, диагностика "пар на пересечение"
            w_out = w;
            has_nb = true;
            nb_out = static_cast<rhs::index_t>(P > 0.0 ? prev : next);
        };

        side_weight(i_surface, P_ij, prim.i_prev, prim.i_next, xi, geom.w_i, geom.i_has_neighbor, geom.i_neighbor);
        side_weight(j_surface, P_ji, prim.j_prev, prim.j_next, xj, geom.w_j, geom.j_has_neighbor, geom.j_neighbor);

        return geom;
    }
}

// ----------------------------------------------------------------------
// fiber_friction_coefficient -- правило смешения mu_ff двух (возможно
// разных) материалов нитей (fiber/fiber_friction_implementation.md, §4/§8).
// Вынесено в ОТДЕЛЬНУЮ функцию сознательно -- если проект примет другое
// правило смешения, менять нужно ТОЛЬКО её, а не compute_fiber_contact_
// contributions. Для одинаковых mu_ff возвращается точное (без std::sqrt)
// значение -- избегаем лишней ошибки округления в самом частом случае (обе
// нити одного материала).
// ----------------------------------------------------------------------
inline double fiber_friction_coefficient(const fiber::FiberMaterial& a, const fiber::FiberMaterial& b)
{
    if (a.mu_ff == b.mu_ff) return a.mu_ff;
    return std::sqrt(a.mu_ff * b.mu_ff);
}

// ----------------------------------------------------------------------
// coulomb_friction_force -- чистый (без обращения к FiberState/geometry)
// закон динамического трения Кулона без памяти контакта
// (fiber/fiber_friction_implementation.md, §4): v_t -- тангенциальная
// (перпендикулярная n) составляющая относительной скорости точек контакта;
// при mu*N>0 и |v_t|>eps_v сила Ft = -mu*N*v_t/|v_t| (направлена против
// скольжения, |Ft|=mu*N, Ft.n=0 по построению -- v_t уже ортогональна n);
// иначе Ft=0 (нет скольжения либо нет нормальной силы/трения). dim --
// compile-time константа проекта (param.hpp); функция компилируется при
// любом SPH_DIM -- циклы по d<dim корректны для любого dim (тот же паттерн,
// что и в остальном файле).
// ----------------------------------------------------------------------
inline void coulomb_friction_force(const double n[dim], double N, double mu, const double v_rel[dim], double eps_v,
                                    double Ft[dim])
{
    double vn = 0.0;
    for (std::size_t d = 0; d < dim; ++d) vn += v_rel[d] * n[d];

    double vt[dim];
    double vt_norm2 = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        vt[d] = v_rel[d] - vn * n[d];
        vt_norm2 += vt[d] * vt[d];
    }
    const double vt_norm = std::sqrt(vt_norm2);

    if (mu * N > 0.0 && vt_norm > eps_v)
    {
        const double scale = -mu * N / vt_norm;
        for (std::size_t d = 0; d < dim; ++d) Ft[d] = scale * vt[d];
    }
    else
    {
        for (std::size_t d = 0; d < dim; ++d) Ft[d] = 0.0;
    }
}

// ----------------------------------------------------------------------
// fiber_contact_point_velocity -- скорость точки контакта на ОДНОЙ стороне
// (fiber/fiber_friction_implementation.md, §2): v = w*v_bead + (1-w)*v_nb
// при has_nb, иначе v = v_bead (End-случай -- вся скорость приходится на
// концевую бусину, w=1 там уже гарантирован compute_fiber_contact_geometry,
// но эта функция не полагается на конкретное значение w при has_nb=false --
// читает только v_bead). Использует ТЕ ЖЕ w/сосед, что FiberContactGeometry
// -- никакой независимой геометрии для трения не вводится (§2 спецификации,
// шапка файла, пункт 4).
// ----------------------------------------------------------------------
inline void fiber_contact_point_velocity(const fiber::FiberState& state, rhs::index_t bead, double w, bool has_nb,
                                          rhs::index_t nb, double v[dim])
{
    if (has_nb)
    {
        for (std::size_t d = 0; d < dim; ++d) v[d] = w * state.vx[d][bead] + (1.0 - w) * state.vx[d][nb];
    }
    else
    {
        for (std::size_t d = 0; d < dim; ++d) v[d] = state.vx[d][bead];
    }
}

// ----------------------------------------------------------------------
// compute_fiber_contact_contributions -- kernel: для каждого primitive с
// контактом (geom.contact) эмитирует вклады FiberDelta на i, j и (если есть)
// их соседей i_prev/i_next, j_prev/j_next (до 4 частиц на primitive, SPH.md
// §3.2). Момент сохраняется точно: сумма вкладов*mass по всем ДО ЧЕТЫРЁХ
// частиц одного primitive равна нулю -- итоговая сила F = F_n + F_t
// (F_n = geom.F*n -- нормальный закон Герца, не изменён; F_t -- тангенциальное
// трение Кулона, fiber/fiber_friction_implementation.md) на стороне i,
// -F на стороне j, распределены ОДНИМИ И ТЕМИ ЖЕ весами w_i/w_j (см.
// test_fiber_contact.cpp, "сумма импульсов"). При mu_ff=0 (или mu=0 после
// смешения) блок трения не выполняется вообще -- выражения для нормальной
// части остаются буквально теми же, что до внедрения трения (побитовая
// регрессия, test_friction_zero_mu_regression).
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

        const double mi = state.mass[prim.i];
        const double mj = state.mass[prim.j];

        // --- Трение (fiber/fiber_friction_implementation.md, §2-§5) ---------
        // Вычисляется ПОСЛЕ геометрии/нормальной силы Герца (geom.F, geom.n),
        // но ДО распределения итоговой силы по бусинам -- как требует §9
        // спецификации. При mu==0 (mu_ff=0 у обеих нитей -- значение по
        // умолчанию FiberMaterial) has_ft остаётся false и блок ниже не
        // добавляет НИЧЕГО к di/dj/dnb -- нормальная часть остаётся
        // побитово такой же, как до внедрения трения.
        const double mu = fiber_friction_coefficient(materials[static_cast<std::size_t>(state.material_id[prim.i])],
                                                       materials[static_cast<std::size_t>(state.material_id[prim.j])]);
        double Ft[dim] = {};
        bool has_ft = false;
        if (mu > 0.0)
        {
            double v_ci[dim], v_cj[dim];
            fiber_contact_point_velocity(state, prim.i, geom.w_i, geom.i_has_neighbor, geom.i_neighbor, v_ci);
            fiber_contact_point_velocity(state, prim.j, geom.w_j, geom.j_has_neighbor, geom.j_neighbor, v_cj);
            double v_rel[dim];
            for (std::size_t d = 0; d < dim; ++d) v_rel[d] = v_ci[d] - v_cj[d];
            coulomb_friction_force(geom.n, geom.F, mu, v_rel, kFiberContactFrictionMinVt, Ft);
            has_ft = true;
        }

        fiber::FiberDelta di{}, dj{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            di.dv[d] = geom.w_i * geom.F * geom.n[d] / mi;
            dj.dv[d] = -geom.w_j * geom.F * geom.n[d] / mj;
        }
        // Ft -- сила на стороне i (n направлена "от j к i", v_rel=v_ci-v_cj,
        // §5 спецификации: F=F_n+F_t на i, -F на j), поэтому на j действует
        // -Ft -- тот же знак, что уже применён выше к нормальной части.
        if (has_ft)
        {
            for (std::size_t d = 0; d < dim; ++d)
            {
                di.dv[d] += geom.w_i * Ft[d] / mi;
                dj.dv[d] -= geom.w_j * Ft[d] / mj;
            }
        }
        writer.emit(prim.i, di, static_cast<std::uint64_t>(k), 0);
        writer.emit(prim.j, dj, static_cast<std::uint64_t>(k), 2);

        if (geom.i_has_neighbor)
        {
            const double m_nb = state.mass[geom.i_neighbor];
            const double rem = 1.0 - geom.w_i;
            fiber::FiberDelta dnb{};
            for (std::size_t d = 0; d < dim; ++d) dnb.dv[d] = rem * geom.F * geom.n[d] / m_nb;
            if (has_ft)
                for (std::size_t d = 0; d < dim; ++d) dnb.dv[d] += rem * Ft[d] / m_nb;
            writer.emit(geom.i_neighbor, dnb, static_cast<std::uint64_t>(k), 1);
        }
        if (geom.j_has_neighbor)
        {
            const double m_nb = state.mass[geom.j_neighbor];
            const double rem = 1.0 - geom.w_j;
            fiber::FiberDelta dnb{};
            for (std::size_t d = 0; d < dim; ++d) dnb.dv[d] = -rem * geom.F * geom.n[d] / m_nb;
            if (has_ft)
                for (std::size_t d = 0; d < dim; ++d) dnb.dv[d] -= rem * Ft[d] / m_nb;
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
inline double compute_fiber_contact_dt_limit(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                              const std::vector<fiber::FiberMaterial>& materials,
                                              const std::vector<FiberContactPrimitive>& primitives)
{
    const double kInf = std::numeric_limits<double>::max();
    const std::size_t n = primitives.size();
    if (n == 0) return kInf;

    return exec.parallel_reduce(
        n, kInf,
        [&](std::size_t k) -> double {
            const FiberContactPrimitive& prim = primitives[k];
            const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
            if (!geom.contact || !(geom.delta > 0.0)) return kInf;

            const double mi = state.mass[prim.i];
            const double mj = state.mass[prim.j];
            const double m_red = (mi * mj) / (mi + mj);
            const double k_eff = 1.5 * geom.F / geom.delta;
            if (!(k_eff > 0.0)) return kInf;
            return 0.8 / std::sqrt(k_eff / m_red);
        },
        [](double a, double b) { return std::min(a, b); });
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
    double max_penetration = 0.0;
    double min_dt = std::numeric_limits<double>::max();
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
                const double mi = state.mass[prim.i];
                const double mj = state.mass[prim.j];
                const double m_red = (mi * mj) / (mi + mj);
                const double k_eff = (geom.delta > 0.0) ? 1.5 * geom.F / geom.delta : 0.0;
                d1.min_dt = (k_eff > 0.0) ? 0.8 / std::sqrt(k_eff / m_red) : std::numeric_limits<double>::max();
            }
            return d1;
        },
        combine);
    return acc;
}

} // namespace interactions
