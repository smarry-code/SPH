#pragma once
// interactions/fiber_contact.hpp
// Fiber-Fiber-Contact: пространственный механический контакт бусин РАЗНЫХ
// нитей (fiber/SPH.md, раздел «Fiber-Fiber-Contact»; fiber/PLAN.md, этап Ф5,
// «3.2. Типы взаимодействий и primitives», «3.3», «3.7»;
// fiber/fiber_contact_segment_pairs.md -- переход на схему «пара сегментов —
// одна сила», описанную ниже).
//
// Как и Fiber-Bond (interactions/fiber_bond.hpp) и Solid-Solid-Contact
// (interactions/solid_contact.hpp), kernel только ЧИТАЕТ FiberState/
// FiberTopology/FiberMaterial и производит typed contributions (FiberDelta)
// через rhs::ContributionWriter — суммирование отдельно, в
// rhs::reduce_contributions. Никаких atomics/mutex.
//
// ============================================================================
// РЕШЕНИЯ (fiber/fiber_contact_segment_pairs.md формулирует физику и
// алгоритм, но оставляет несколько инженерных деталей на усмотрение
// реализации):
//
// 1. Primitive -- ПАРА СЕГМЕНТОВ, а не пара бусин. Сегмент нити A задаётся
//    первой бусиной a0 и её соседом a1=next(a0), q_A(s)=x_{a0}+s*(x_{a1}-x_{a0}),
//    s in [0,1] (симметрично для B, параметр t). Четыре случая Surface-
//    Surface/Surface-End/End-Surface/End-End старой схемы ИСЧЕЗЛИ вместе с
//    enum FiberContactCase: контакт с концом нити — это просто ближайшая
//    точка при s=0 или s=1 на крайнем сегменте, отдельной формулы не нужно.
//
//    Построение primitive из spatial-пары бусин (i,j) (find_pairs_parallel +
//    фильтр object_id, как раньше) отображает КАЖДУЮ сторону пары в сегмент
//    по правилу «сегмент принадлежит первой бусине» (тот же приём, что
//    fiber_solid_contact.hpp, решение 2): next(k)!=-1 -> сегмент (k,next(k));
//    next(k)==-1 && prev(k)!=-1 (k — последняя бусина нити) -> ПАРА
//    ОТБРАСЫВАЕТСЯ целиком (эту зону покроет пара (prev(k),*), которая
//    гарантированно найдётся — см. вывод неравенства треугольника у
//    fiber_contact_search_radius ниже); next(k)==-1 && prev(k)==-1
//    (изолированная бусина) -> точечный сегмент (k,-1). Пара бусин {i,j}
//    приходит из поиска РОВНО ОДИН РАЗ (канонический i<j в SpatialIndex),
//    поэтому пара сегментов тоже строится ровно один раз — дублей на уровне
//    построения нет, сортировка/удаление повторов не нужны.
//
//    Ближайшие точки двух отрезков ищет ЕДИНСТВЕННАЯ функция
//    fiber_segment_closest_params (Ericson, «Real-Time Collision Detection»,
//    §5.1.9, ClosestPtSegmentSegment) со SNAPPING результата к 0/1 в
//    η=1e-9-окрестности границы (без него пересечение ровно в узле сети то
//    классифицируется как внутренняя точка сегмента, то как его конец —
//    прототип на 20000 случайных скрещиваний без snapping дал ~6% неверной
//    классификации, со snapping — 0, fiber/fiber_contact_segment_pairs.md,
//    §2).
//
//    ПРАВИЛО ВЛАДЕНИЯ (§3 спецификации) гарантирует ровно один принятый
//    контакт на физическое пересечение. Идея: две соседние вершины сегментов
//    нити A делят узел x_{a0} — он одновременно конец S_{a_prev} (при s=1) и
//    начало S_{a0} (при s=0). Каждый сегмент, у которого ближайшая точка
//    попала на ЕГО СОБСТВЕННЫЙ конец, а этот конец продолжается дальше по
//    нити, ОТКЛОНЯЕТ примитив (условия A1/B1: s=1 и a_has_next_segment, или
//    t=1 и b_has_next_segment) — эту точку "видит" сосед. Обратная сторона
//    (условия A2/B2: s=0 и a_prev!=-1, либо t=0 и b_prev!=-1) должна
//    подтвердить, что она — законный владелец узла: та же функция
//    fiber_segment_closest_params пересчитывается для пары (S_a_prev,S_B)
//    (сегмент S_a_prev=(a_prev,a0)) или (S_A,S_b_prev), и точка принимается,
//    только если у ЭТОГО соседнего сегмента ближайшая точка оказалась на ЕГО
//    дальнем конце (s'=1 или t'=1 соответственно) — иначе минимум на самом
//    деле лежит ВНУТРИ соседнего сегмента, и владеет им он. Первый/последний
//    сегмент нити (a_prev==-1 либо !a_has_next_segment) принимает свою
//    концевую точку БЕЗ условий — там соревноваться не с кем.
//    Почему этого достаточно (полный разбор — fiber/fiber_contact_segment_
//    pairs.md, §3): если минимум лежит строго внутри сегмента, его "видит"
//    только этот сегмент (сосед видит ту же зону как свой конец и
//    отклоняется по A1/B1). Если минимум ровно в узле, "верхний" по ходу
//    нити сегмент (для которого это s=1) отклоняется по A1, а "нижний" (для
//    которого это s=0) принимает его по A2 — после проверки, что более
//    ранний сосед не претендует на ту же точку изнутри. Если сегмент
//    отклонён по A1, то по неравенству треугольника dist(S_{a0},S_B) <=
//    dist(x_{a0},S_B) = D — контакт не теряется, а переходит к следующему
//    сегменту, который гарантированно есть среди primitives (тот же радиус
//    поиска). Изогнутые нити с НЕСКОЛЬКИМИ локальными минимумами (седло) —
//    единственный случай, когда несколько primitives одного и того же
//    физического пересечения остаются принятыми: это разные физические
//    точки контакта, не дубликат одной и той же.
//    Правило живёт ЦЕЛИКОМ внутри compute_fiber_contact_geometry — kernel,
//    dt_limit и диагностика видят ОДНО И ТО ЖЕ множество принятых контактов.
//    Отклонённый primitive возвращает geom.contact=false,
//    geom.rejected_by_ownership=true (для диагностики).
//
// 2. d/E*/R* при разных материалах нитей a,b. SPH.md даёт формулу Meyer
//    только для ОДИНАКОВЫХ нитей (E* = E/(2(1-nu^2)), R* = sqrt(d/4)), но в
//    разделе «Fiber-Solid-Contact», «Постановка проекта: коэффициент K_fs»
//    для ПАРЫ РАЗНЫХ материалов (Fiber-Solid) уже выведена корректная
//    Hertz-обобщённая формула:
//      1/E* = (1-nu_a^2)/E_a + (1-nu_b^2)/E_b,   1/R* = 1/r_a + 1/r_b.
//    Та же логика (эффективные модуль/кривизна пары контактирующих тел из
//    контактной механики Герца) применена здесь к паре нитей a,b: E_a,nu_a
//    из FiberMaterial стороны A (материал первой бусины сегмента, a0),
//    E_b,nu_b — стороны B (по бусине b0), r_a=d_a/2, r_b=d_b/2 (fiber_solid_
//    contact.hpp, решение 4 — материал стороны берётся по её "якорной"
//    бусине). Порог проникновения: d_eff = r_a+r_b, delta = d_eff - D,
//    контакт при D < d_eff.
//    САМОПРОВЕРКА (одинаковые нити, E_a=E_b=E, nu_a=nu_b=nu, d_a=d_b=d):
//      1/E* = 2(1-nu^2)/E => E* = E/(2(1-nu^2))                    — совпадает с Meyer;
//      r_a=r_b=d/2 => d_eff=d (порог D<d, как у Meyer);
//      1/R* = 1/(d/2)+1/(d/2) = 4/d => R*=d/4, sqrt(R*)=sqrt(d/4)  — совпадает с Meyer.
//
// 3. Радиус поиска d_max + 2*L_max (fiber_contact_search_radius). Пара
//    сегментов в контакте должна гарантированно породить spatial-пару своих
//    первых бусин a0,b0. Если D_AB < d_eff, то по неравенству треугольника
//        |x_{a0}-x_{b0}| <= |x_{a0}-q_A| + D_AB + |q_B-x_{b0}|
//                         <= L_A + d_eff + L_B <= 2*L_max + d_max,
//    где L_A,L_B — ТЕКУЩИЕ (не номинальные l0!) длины сегментов A,B: связи
//    Fiber-Bond работают только на растяжение и рвутся лишь в
//    Postprocessing ПОСЛЕ шага (fiber/SPH.md), поэтому внутри шага сегмент
//    может быть длиннее l0 и даже длиннее l_break — номинальное l0 как
//    оценка не годится. L_max — максимум ТЕКУЩЕЙ длины |x_i-x_j| по всем
//    topology.bonds расчёта (fiber_contact_search_radius получает state и
//    exec именно для этого — считает максимум через exec.parallel_reduce).
//
// 4. Трение. Kernel этого файла считает ТОЛЬКО нормальную силу Герца.
//    Трение нить-нить -- отдельная стадия шага между вычислением dt и
//    ODESolver::advance (проекция на уровне скоростей, ограничение импульса
//    трения; interactions/fiber_friction.hpp, fiber/fiber_friction_projection.md).
//    Kernel лишь заполняет по одному слоту FiberFrictionContact на primitive
//    (геометрия принятого контакта, N=geom.F, mu из fiber_friction_coefficient
//    -- правило смешения mu для разных материалов живёт здесь, менять только
//    его). Слот active только при contact && mu>0, поэтому при mu_ff=0 у всех
//    материалов стадия трения не выполняется и результат побитово прежний.
//    dt_limit трением не затрагивается.
//
// 5. Закон нормальной силы (fiber/plan-fiber-contact-types.md, этап 2;
//    interactions/fiber_contact_law.hpp): ключ fiber_fiber_contact_type,
//    0 = Герц (прежний код, побитово), 1-3 -- contact::penalty_force на g, pr
//    из fiber_fiber_pair_inputs (одна функция для ядра, предела шага и
//    диагностики). Геометрия не менялась. Предел шага типов 1-3 -- минимум
//    res.dt_limit по слотам FiberContactSearchResult::pair_dt_limit.
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "contact_force.hpp"
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_contact_common.hpp" // kFiberContactSnapEta, fiber_max_current_bond_length
#include "interactions/fiber_contact_law.hpp"    // FiberContactStepParams, fiber_penalty_force
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

// Порог |ê_A x ê_B| (единичные направления сегментов), ниже которого пара
// осей в вырожденной ветке D<kFiberContactMinD (см. compute_fiber_contact_
// geometry, §4 спецификации) считается "совпадающей параллельной" —
// направление нормали через e_A x e_B неопределимо (числитель ~0 при
// ненулевом знаменателе, либо оба направления вырождены сами по себе —
// точечный сегмент). eps -- безразмерный (|ê_A x ê_B| = sin(угол между
// направлениями) для единичных векторов). Контакт в этой ветке НЕ
// формируется (geom.degenerate=true, fiber/fiber_contact_segment_pairs.md,
// §4) -- в отличие от старой схемы, здесь это не "почти параллельная"
// аппроксимация нормали, а честный отказ: у двух совпадающих/параллельных
// осей нормаль контакта в принципе не определена этой геометрией.
constexpr double kFiberContactParallelEps = 1e-3;

// Защита от деления на ~0 при нормировке n (совпадающие точки/оси), а также
// (в compute_fiber_contact_geometry) индикатор того, что хотя бы один из
// сегментов вырожден в точку (|e|~0) -- в этом случае направление e_A x e_B
// не определяет ось однозначно, см. kFiberContactParallelEps выше.
constexpr double kFiberContactMinD = 1e-12;

// Относительный порог вырожденности отрезка в fiber_segment_closest_params
// (§2 спецификации): |d|^2 <= eps*L_ref^2, L_ref^2 = max(|d1|^2,|d2|^2) --
// больший из квадратов длин ДВУХ сравниваемых отрезков (в функции нет
// отдельного параметра "d" для L_ref, поэтому в его роли выступает масштаб
// самой пары отрезков; если оба почти вырождены, L_ref~0 и оба сразу
// попадают в общую ветку "обе стороны -- точки", деления на 0 не возникает).
constexpr double kFiberSegmentDegenerateEpsRel = 1e-12;

// Snapping параметров ближайших точек s,t к 0/1 (§2 спецификации, η) --
// kFiberContactSnapEta теперь в interactions/fiber_contact_common.hpp (общая
// с Fiber-Solid-Contact, fiber/fiber_solid_contact_fix.md, §1 п.2).

// ----------------------------------------------------------------------
// FiberContactPrimitive -- пара сегментов осей нитей A,B (fiber/
// fiber_contact_segment_pairs.md, §1). Сегмент A = (a0, a1=next(a0));
// a1==-1 -- точечный сегмент (изолированная бусина без соседей). a_prev,
// a_has_next_segment -- нужны правилу владения (см. шапку файла, пункт 1) и
// сохранены при построении primitive, чтобы kernel/геометрия не перечитывали
// topology повторно (тот же приём, что FiberSolidContactPrimitive::
// a_prev/b_is_last в fiber_solid_contact.hpp). Симметрично для B.
// ----------------------------------------------------------------------
struct FiberContactPrimitive
{
    // Сторона A
    rhs::index_t a0 = 0;
    long long    a1 = -1;      // next(a0); -1 -> точечный сегмент (изолированная бусина)
    long long    a_prev = -1;  // prev(a0); нужен правилу владения
    bool         a_has_next_segment = false; // a1 != -1 && next(a1) != -1
    // Сторона B -- симметрично
    rhs::index_t b0 = 0;
    long long    b1 = -1;
    long long    b_prev = -1;
    bool         b_has_next_segment = false;
};

// ----------------------------------------------------------------------
// fiber_contact_search_radius -- d_max + 2*L_max (см. пункт 3 в шапке
// файла: вывод неравенства треугольника там же). L_max -- ТЕКУЩАЯ (не
// номинальная l0!) максимальная длина связи, interactions/
// fiber_contact_common.hpp::fiber_max_current_bond_length (общая с Fiber-
// Solid-Contact, fiber/fiber_solid_contact_fix.md, §2 -- вычисление ровно
// одно, побитово то же, что было здесь до выноса).
// ----------------------------------------------------------------------
inline double fiber_contact_search_radius(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                           const fiber::FiberTopology& topology,
                                           const std::vector<fiber::FiberMaterial>& materials)
{
    double d_max = 0.0;
    for (const auto& m : materials) d_max = std::max(d_max, m.d);

    const double L_max = fiber_max_current_bond_length(exec, state, topology);

    return d_max + 2.0 * L_max;
}

// ----------------------------------------------------------------------
// build_fiber_contact_search -- SpatialIndex<dim> (dim==3, Р9) по бусинам
// Fiber, фильтр i<j (канонический, из SpatialIndex) И object_id[i]!=
// object_id[j] (разные нити -- самоконтакт одной нити не формулируется,
// SPH.md, «Fiber-Fiber-Contact», последний абзац). Каждая сторона найденной
// пары отображается в сегмент по правилу «сегмент принадлежит первой
// бусине» (шапка файла, пункт 1); если ХОТЯ БЫ ОДНА сторона -- последняя
// бусина нити (next==-1, prev!=-1), пара отбрасывается ЦЕЛИКОМ (эту зону
// покроет пара (prev(k),*)).
// ----------------------------------------------------------------------
struct FiberContactSearchResult
{
    std::vector<FiberContactPrimitive> primitives;
    std::size_t pair_count = 0; // spatial-пары после фильтра object_id (широкая фаза, до отбрасывания по таблице сторон)
    // Постоянная рабочая память поиска: результат живёт между шагами, поэтому
    // индекс и вектор пар не создаются заново на каждом вызове (память только
    // растёт).
    rhs::SpatialIndex<dim> index;
    std::vector<std::pair<rhs::index_t, rhs::index_t>> pairs;
    // Предел шага на primitive (законы 1-3, FiberContactStepParams::pair_dt_limit):
    // заполняет ядро, минимум берёт compute_fiber_contact_dt_limit.
    std::vector<double> pair_dt_limit;
};

inline void build_fiber_contact_search(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                        const fiber::FiberTopology& topology,
                                        const std::vector<fiber::FiberMaterial>& materials,
                                        FiberContactSearchResult& out)
{
    out.primitives.clear();
    out.pair_count = 0;
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
        const double radius = fiber_contact_search_radius(exec, state, topology, materials);

        std::array<const double*, dim> coords{};
        for (std::size_t d = 0; d < dim; ++d) coords[d] = state.x[d].data();
        out.index.build(coords, n, radius > 0.0 ? radius : 1.0);

        auto accept = [&](rhs::index_t a, rhs::index_t b) -> bool {
            return state.object_id[a] != state.object_id[b];
        };

        out.index.find_pairs_parallel(exec, radius, accept, out.pairs);
        out.pair_count = out.pairs.size();

        out.primitives.reserve(out.pairs.size());
        for (const auto& pr : out.pairs)
        {
            const rhs::index_t k_a = pr.first;
            const rhs::index_t k_b = pr.second;

            // Сторона A -- по таблице §1 спецификации.
            const long long a_prev = topology.prev[k_a];
            const long long a_next = topology.next[k_a];
            long long a1 = -1;
            if (a_next != -1)
            {
                a1 = a_next;
            }
            else if (a_prev != -1)
            {
                continue; // k_a -- последняя бусина своей нити: пара отбрасывается целиком
            }
            // иначе: k_a -- изолированная бусина, точечный сегмент (a1 остаётся -1)

            // Сторона B -- симметрично.
            const long long b_prev = topology.prev[k_b];
            const long long b_next = topology.next[k_b];
            long long b1 = -1;
            if (b_next != -1)
            {
                b1 = b_next;
            }
            else if (b_prev != -1)
            {
                continue;
            }

            FiberContactPrimitive prim;
            prim.a0 = k_a;
            prim.a1 = a1;
            prim.a_prev = a_prev;
            prim.a_has_next_segment = (a1 != -1) && (topology.next[static_cast<std::size_t>(a1)] != -1);
            prim.b0 = k_b;
            prim.b1 = b1;
            prim.b_prev = b_prev;
            prim.b_has_next_segment = (b1 != -1) && (topology.next[static_cast<std::size_t>(b1)] != -1);

            out.primitives.push_back(prim);
        }
    }
}

// ----------------------------------------------------------------------
// fiber_segment_closest_params -- ближайшие точки двух отрезков [p0,p1] и
// [r0,r1] (C. Ericson, Real-Time Collision Detection, §5.1.9,
// ClosestPtSegmentSegment). ЕДИНСТВЕННАЯ функция на ВСЕ вызовы (§2
// спецификации): и для "собственного" primitive (S_A,S_B), и для проверок
// правила владения (S_a_prev,S_B)/(S_A,S_b_prev) в compute_fiber_contact_
// geometry ниже -- значит результаты должны быть побитово воспроизводимы
// для одинаковых входов (правило владения опирается на точное сравнение
// s/t с 0.0/1.0 после snapping). Без исключений, аллокаций и std::function,
// без захваченных лямбд и ручного переупорядочивания выражений -- годится
// для переноса в CUDA kernel добавлением __host__ __device__. Проект
// собирается Clang с -ffp-contract=on по умолчанию (FMA-сжатие внутри
// одного выражения детерминировано и не зависит от места inlining) -- флаги
// не меняются.
//
// Вырожденные отрезки (длина~0) обрабатываются как точки: относительный
// порог kFiberSegmentDegenerateEpsRel*L_ref^2 (комментарий у константы
// выше). Snapping (η=kFiberContactSnapEta): без него пересечение ровно в
// узле сети то классифицируется как внутренняя точка сегмента (s~1-1e-16),
// то как его конец -- см. шапку файла, пункт 1.
// ----------------------------------------------------------------------
inline void fiber_segment_closest_params(const double p0[3], const double p1[3], const double r0[3],
                                          const double r1[3], double& s, double& t)
{
    const double d1[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    const double d2[3] = {r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2]};
    const double rr[3] = {p0[0] - r0[0], p0[1] - r0[1], p0[2] - r0[2]};

    const double a = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]; // |d1|^2
    const double e = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]; // |d2|^2
    const double f = d2[0] * rr[0] + d2[1] * rr[1] + d2[2] * rr[2];

    const double Lref2 = (a > e) ? a : e;
    const bool a_degenerate = (a <= kFiberSegmentDegenerateEpsRel * Lref2);
    const bool e_degenerate = (e <= kFiberSegmentDegenerateEpsRel * Lref2);

    if (a_degenerate && e_degenerate)
    {
        // Оба отрезка -- фактически точки.
        s = 0.0;
        t = 0.0;
    }
    else if (a_degenerate)
    {
        // Первый отрезок вырожден в точку -- проекция на второй.
        s = 0.0;
        t = f / e;
        t = (t < 0.0) ? 0.0 : ((t > 1.0) ? 1.0 : t);
    }
    else
    {
        const double c = d1[0] * rr[0] + d1[1] * rr[1] + d1[2] * rr[2];
        if (e_degenerate)
        {
            // Второй отрезок вырожден в точку -- проекция на первый.
            t = 0.0;
            s = -c / a;
            s = (s < 0.0) ? 0.0 : ((s > 1.0) ? 1.0 : s);
        }
        else
        {
            // Общий невырожденный случай.
            const double b = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2];
            const double denom = a * e - b * b; // всегда >= 0 (неравенство Коши-Буняковского)

            if (denom > 0.0)
            {
                s = (b * f - c * e) / denom;
                s = (s < 0.0) ? 0.0 : ((s > 1.0) ? 1.0 : s);
            }
            else
            {
                // Отрезки параллельны -- s не определяется однозначно, берём s=0
                // (Ericson, ClosestPtSegmentSegment) и находим t по нему ниже.
                s = 0.0;
            }

            t = (b * s + f) / e;

            if (t < 0.0)
            {
                t = 0.0;
                s = -c / a;
                s = (s < 0.0) ? 0.0 : ((s > 1.0) ? 1.0 : s);
            }
            else if (t > 1.0)
            {
                t = 1.0;
                s = (b - c) / a;
                s = (s < 0.0) ? 0.0 : ((s > 1.0) ? 1.0 : s);
            }
        }
    }

    // --- Snapping к 0/1 (см. комментарий у функции выше) --------------------
    if (s < kFiberContactSnapEta) s = 0.0;
    else if (s > 1.0 - kFiberContactSnapEta) s = 1.0;
    if (t < kFiberContactSnapEta) t = 0.0;
    else if (t > 1.0 - kFiberContactSnapEta) t = 1.0;
}

// ----------------------------------------------------------------------
// FiberContactGeometry -- итог геометрического решения ОДНОГО primitive:
// параметры ближайших точек s,t (после snapping и правила владения),
// нормаль/расстояние/сила Герца (если контакт принят). Общая точка входа
// для compute_fiber_contact_contributions, compute_fiber_contact_dt_limit и
// диагностики -- геометрия (включая правило владения) вычисляется РОВНО
// ОДИН РАЗ на каждый вызывающий проход (аналог PairGeom/PairProps в
// contact_force.hpp), но не кэшируется между проходами (см. также
// interactions/solid_contact.hpp -- там геометрия тоже пересчитывается
// отдельно в diagnostics).
// ----------------------------------------------------------------------
struct FiberContactGeometry
{
    bool contact = false;
    double n[dim] = {};
    double D = 0.0;
    double F = 0.0;     // сила Герца, >=0
    double delta = 0.0; // проникновение d_eff - D, >0 при контакте

    double s = 0.0; // параметр ближайшей точки на сегменте A (после snapping)
    double t = 0.0; // параметр ближайшей точки на сегменте B (после snapping)

    bool rejected_by_ownership = false; // отклонён правилом владения (§3) -- точку "видит" сосед
    bool degenerate = false;            // D~0 при совпадающих/параллельных (или неопределённых) осях -- контакт не сформирован (§4)
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

        const Vec3 xa0 = point(prim.a0);
        const Vec3 xa1 = (prim.a1 == -1) ? xa0 : point(static_cast<rhs::index_t>(prim.a1));
        const Vec3 xb0 = point(prim.b0);
        const Vec3 xb1 = (prim.b1 == -1) ? xb0 : point(static_cast<rhs::index_t>(prim.b1));

        double s = 0.0, t = 0.0;
        fiber_segment_closest_params(xa0.data(), xa1.data(), xb0.data(), xb1.data(), s, t);

        // --- Правило владения (§3 спецификации; вывод -- шапка файла, пункт 1) --
        bool accept = true;
        if (s >= 1.0 && prim.a_has_next_segment)
        {
            // A1: узел x_{a1} принадлежит следующему сегменту A -- отклонить.
            accept = false;
        }
        else if (s <= 0.0 && prim.a_prev != -1)
        {
            // A2: узел x_{a0} -- общий минимум, владелец решается пересчётом
            // для (S_a_prev, S_B) = ((a_prev,a0), (b0,b1)); принимаем, только
            // если у СОСЕДА минимум оказался на ЕГО дальнем конце (s'=1) --
            // иначе минимум внутри S_a_prev, и владеет им он.
            const Vec3 xap = point(static_cast<rhs::index_t>(prim.a_prev));
            double sp = 0.0, tp = 0.0;
            fiber_segment_closest_params(xap.data(), xa0.data(), xb0.data(), xb1.data(), sp, tp);
            if (!(sp >= 1.0)) accept = false;
        }
        if (accept)
        {
            if (t >= 1.0 && prim.b_has_next_segment)
            {
                // B1 -- симметрично A1.
                accept = false;
            }
            else if (t <= 0.0 && prim.b_prev != -1)
            {
                // B2 -- симметрично A2: (S_A,S_b_prev) = ((a0,a1), (b_prev,b0)),
                // владелец подтверждён, если t'=1.
                const Vec3 xbp = point(static_cast<rhs::index_t>(prim.b_prev));
                double sp = 0.0, tp = 0.0;
                fiber_segment_closest_params(xa0.data(), xa1.data(), xbp.data(), xb0.data(), sp, tp);
                if (!(tp >= 1.0)) accept = false;
            }
        }

        geom.s = s;
        geom.t = t;
        if (!accept)
        {
            geom.rejected_by_ownership = true;
            return geom; // geom.contact остаётся false
        }

        // --- Геометрия принятого контакта (§4 спецификации) ---------------------
        const Vec3 qA = {xa0[0] + s * (xa1[0] - xa0[0]), xa0[1] + s * (xa1[1] - xa0[1]),
                          xa0[2] + s * (xa1[2] - xa0[2])};
        const Vec3 qB = {xb0[0] + t * (xb1[0] - xb0[0]), xb0[1] + t * (xb1[1] - xb0[1]),
                          xb0[2] + t * (xb1[2] - xb0[2])};
        const Vec3 diff = sub(qA, qB);
        double D = norm(diff);
        Vec3 n{};

        if (D < kFiberContactMinD)
        {
            // Оси (почти) пересекаются -- нормаль из e_A x e_B, ориентация
            // "от B к A" по серединам сегментов (§4: не по q_A/q_B -- они
            // здесь совпадают и ориентацию не несут).
            const Vec3 eA = sub(xa1, xa0);
            const Vec3 eB = sub(xb1, xb0);
            const double eAn = norm(eA);
            const double eBn = norm(eB);
            const Vec3 cr = cross(eA, eB);
            const double crn = norm(cr);

            // |ê_A x ê_B| = crn/(eAn*eBn) -- сравниваем без деления (крайние
            // случаи eAn~0/eBn~0 -- точечный сегмент, направление одной из
            // осей не определено вовсе, тоже считаем вырожденным случаем).
            const bool degenerate_axes =
                (eAn < kFiberContactMinD) || (eBn < kFiberContactMinD) || (crn < kFiberContactParallelEps * eAn * eBn);
            if (degenerate_axes)
            {
                geom.degenerate = true;
                return geom; // geom.contact остаётся false
            }

            n = {cr[0] / crn, cr[1] / crn, cr[2] / crn};
            const Vec3 mA = {0.5 * (xa0[0] + xa1[0]), 0.5 * (xa0[1] + xa1[1]), 0.5 * (xa0[2] + xa1[2])};
            const Vec3 mB = {0.5 * (xb0[0] + xb1[0]), 0.5 * (xb0[1] + xb1[1]), 0.5 * (xb0[2] + xb1[2])};
            if (dot(n, sub(mA, mB)) < 0.0) { n = {-n[0], -n[1], -n[2]}; }
            D = 0.0;
        }
        else
        {
            n = {diff[0] / D, diff[1] / D, diff[2] / D};
        }

        // --- Hertz-сила (пункт 2 в шапке файла); материал стороны -- по её
        //     "якорной" бусине a0/b0 -------------------------------------------
        const fiber::FiberMaterial& mat_a = materials[static_cast<std::size_t>(state.material_id[prim.a0])];
        const fiber::FiberMaterial& mat_b = materials[static_cast<std::size_t>(state.material_id[prim.b0])];
        const double ra = 0.5 * mat_a.d;
        const double rb = 0.5 * mat_b.d;
        const double d_eff = ra + rb;
        if (!(D < d_eff)) return geom; // нет проникновения

        const double delta = d_eff - D;
        const double inv_Estar = (mat_a.E > 0.0 ? (1.0 - mat_a.nu * mat_a.nu) / mat_a.E : 0.0) +
                                  (mat_b.E > 0.0 ? (1.0 - mat_b.nu * mat_b.nu) / mat_b.E : 0.0);
        if (!(inv_Estar > 0.0)) return geom; // некорректный материал (E<=0 с обеих сторон) -- защита
        const double E_star = 1.0 / inv_Estar;
        const double R_star = (ra > 0.0 && rb > 0.0) ? (ra * rb) / (ra + rb) : 0.0;
        const double F = (4.0 / 3.0) * E_star * std::sqrt(R_star) * std::pow(delta, 1.5);

        geom.contact = true;
        geom.D = D;
        geom.delta = delta;
        geom.F = F;
        geom.n[0] = n[0];
        geom.n[1] = n[1];
        geom.n[2] = n[2];
        return geom;
    }
}

// ----------------------------------------------------------------------
// fiber_contact_side_inv_mass -- 1/m эффективной массы точки контакта на
// ОДНОЙ стороне (§5 спецификации): (1-s)^2/m0 + s^2/m1 при наличии соседа,
// точечный сегмент (k1==-1) -> 1/m0. Общая функция для dt_limit (пункт 6 в
// шапке файла) и стадии трения (fiber_friction.hpp).
// ----------------------------------------------------------------------
inline double fiber_contact_side_inv_mass(const fiber::FiberState& st, rhs::index_t k0, long long k1, double s)
{
    if (k1 == -1) return 1.0 / st.mass[k0];
    const double w0 = 1.0 - s;
    const double w1 = s;
    return (w0 * w0) / st.mass[k0] + (w1 * w1) / st.mass[static_cast<rhs::index_t>(k1)];
}

// ----------------------------------------------------------------------
// fiber_friction_coefficient -- правило смешения mu_ff двух (возможно
// разных) материалов нитей (fiber/fiber_friction_implementation.md, §4/§8).
// Вынесено в ОТДЕЛЬНУЮ функцию сознательно -- если проект примет другое
// правило смешения, менять нужно ТОЛЬКО её. Для одинаковых mu_ff
// возвращается точное (без std::sqrt) значение -- избегаем лишней ошибки
// округления в самом частом случае (обе нити одного материала).
// ----------------------------------------------------------------------
inline double fiber_friction_coefficient(const fiber::FiberMaterial& a, const fiber::FiberMaterial& b)
{
    if (a.mu_ff == b.mu_ff) return a.mu_ff;
    return std::sqrt(a.mu_ff * b.mu_ff);
}

// ----------------------------------------------------------------------
// FiberFrictionContact -- слот принятого контакта для стадии трения
// (interactions/fiber_friction.hpp; fiber/fiber_friction_projection.md, §3).
// ПО ОДНОМУ слоту на primitive, запись по индексу k -- без гонок и без
// аллокаций внутри kernel. POD без указателей и std::vector: переносится в
// device-память как есть. Объявлен здесь (а не в fiber_friction.hpp), потому
// что его заполняет kernel нормального контакта, а fiber_friction.hpp
// включает этот файл (fiber_contact_side_inv_mass) -- иначе цикл include.
// active = контакт принят && mu>0; неактивный слот остаётся value-initialized.
// ----------------------------------------------------------------------
struct FiberFrictionContact
{
    bool active = false;
    rhs::index_t a0 = 0;
    long long a1 = -1; // -1 -- точечный сегмент
    rhs::index_t b0 = 0;
    long long b1 = -1;
    double s = 0.0;
    double t = 0.0;
    double n[dim] = {}; // нормаль от B к A
    double N = 0.0;     // нормальная сила (Герц geom.F или закон 1-3), >=0
    double mu = 0.0;
};

// ----------------------------------------------------------------------
// fiber_fiber_pair_inputs -- g, pr пары для законов 1-3
// (fiber/plan-fiber-contact-types.md, §2; fiber/SPH.md, Fiber-Fiber-Contact,
// «Законы нормальной силы»). ЕДИНСТВЕННАЯ сборка входов: ею пользуются ядро,
// предел шага и диагностика (одна формула на одном множестве контактов).
//   n от B к A (geom.n; при D = 0 -- из e_A x e_B, уже в геометрии);
//   g.p = delta = r_a + r_b - D; v_A = (1-s) v_a0 + s v_a1 (a1 = -1: v_a0),
//   v_B аналогично; g.pdot = -n.(v_A - v_B); g.v_rel = |v_A - v_B|;
//   g.r = r_a + r_b (> 0 и при пересечении осей D = 0, где Swegle иначе
//   отключился бы); m_red по эффективным массам точек
//   (fiber_contact_side_inv_mass); g_eff = G_a G_b/(G_a+G_b) по a0, b0;
//   r_eff = R_star; c_pair = min(c_a, c_b); dt_ref = min(dt_ref_A, dt_ref_B).
// Скорости -- фаза предиктора (solver.predict до RHS).
// ----------------------------------------------------------------------
inline void fiber_fiber_pair_inputs(const fiber::FiberState& state, const std::vector<fiber::FiberMaterial>& materials,
                                    const FiberContactPrimitive& prim, const FiberContactGeometry& geom,
                                    const FiberContactStepParams& params, contact::PairGeom& g,
                                    contact::PairProps& pr)
{
    const fiber::FiberMaterial& mat_a = materials[static_cast<std::size_t>(state.material_id[prim.a0])];
    const fiber::FiberMaterial& mat_b = materials[static_cast<std::size_t>(state.material_id[prim.b0])];
    const double ra = 0.5 * mat_a.d;
    const double rb = 0.5 * mat_b.d;

    double pdot = 0.0, v2 = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        const double vA = (prim.a1 == -1)
                              ? state.vx[d][prim.a0]
                              : (1.0 - geom.s) * state.vx[d][prim.a0] +
                                    geom.s * state.vx[d][static_cast<rhs::index_t>(prim.a1)];
        const double vB = (prim.b1 == -1)
                              ? state.vx[d][prim.b0]
                              : (1.0 - geom.t) * state.vx[d][prim.b0] +
                                    geom.t * state.vx[d][static_cast<rhs::index_t>(prim.b1)];
        const double vrel = vA - vB;
        pdot -= geom.n[d] * vrel;
        v2 += vrel * vrel;
    }
    g.p = geom.delta;
    g.pdot = pdot;
    g.v_rel = std::sqrt(v2);
    g.r = ra + rb;

    const double m_A = 1.0 / fiber_contact_side_inv_mass(state, prim.a0, prim.a1, geom.s);
    const double m_B = 1.0 / fiber_contact_side_inv_mass(state, prim.b0, prim.b1, geom.t);
    pr.m_red = (m_A * m_B) / (m_A + m_B);
    pr.g_eff = fiber_contact_g_eff(fiber_shear_modulus(mat_a), fiber_shear_modulus(mat_b));
    pr.r_eff = (ra > 0.0 && rb > 0.0) ? (ra * rb) / (ra + rb) : 0.0;
    pr.c_pair = std::min(fiber_sound_speed(mat_a), fiber_sound_speed(mat_b));
    pr.dt_ref = std::min(fiber_side_dt_ref(*params.bead_dt_ref, prim.a0, prim.a1),
                         fiber_side_dt_ref(*params.bead_dt_ref, prim.b0, prim.b1));
}

// Сила и предел шага закона 1-3 для принятого контакта (geom.contact).
inline contact::Result fiber_fiber_contact_penalty(const fiber::FiberState& state,
                                                   const std::vector<fiber::FiberMaterial>& materials,
                                                   const FiberContactPrimitive& prim, const FiberContactGeometry& geom,
                                                   const FiberContactStepParams& params)
{
    contact::PairGeom g;
    contact::PairProps pr;
    fiber_fiber_pair_inputs(state, materials, prim, geom, params, g, pr);
    if (params.law.type == kFiberLawHertzDamped)
        return fiber_hertz_damped_force(geom.F, geom.delta, g, pr, params.law.kp, params.damper_only);
    return fiber_penalty_force(params, g, pr);
}

// Проверка входов законов 1-3 до parallel_for.
inline void fiber_fiber_check_params(const FiberContactStepParams& params, const fiber::FiberState& state)
{
    if (params.law.type == 0) return;
    if (params.law.type < 0 || params.law.type > kFiberLawHertzDamped)
        throw std::invalid_argument("Fiber-Fiber: тип закона вне 0..4");
    if (params.bead_dt_ref == nullptr || params.bead_dt_ref->size() != state.size())
        throw std::invalid_argument("Fiber-Fiber: тип 1-3 требует bead_dt_ref размера state.size()");
}

// ----------------------------------------------------------------------
// compute_fiber_contact_contributions -- kernel: для каждого primitive с
// принятым контактом (geom.contact) эмитирует вклады FiberDelta НОРМАЛЬНОЙ
// силы Герца на a0, (a1, если существует и вес не нулевой), b0, (b1 --
// симметрично). Слоты ContributionMeta: 0=a0, 1=a1, 2=b0, 3=b1 (§5
// спецификации). Сила F_n = geom.F*n прикладывается на стороне A весами
// (1-s,s), -F_n на стороне B весами (1-t,t) -- момент сохраняется точно.
// Вклад с нулевым весом (s==0 -> a1, s==1 -> a0, симметрично для b0/b1) НЕ
// эмитируется -- сравнение с границей делается через строгое неравенство
// (geom.s/geom.t гарантированно в [0,1] после clamp в
// fiber_segment_closest_params).
//
// Трения здесь нет (fiber/fiber_friction_projection.md, §1, §3): kernel
// ДОПОЛНИТЕЛЬНО заполняет friction_contacts[k] для каждого primitive.
// friction_contacts.resize(primitives.size()) выполняется в начале, ДО
// parallel_for (внутри kernel аллокаций нет; память только растёт, если
// вызывающая сторона переиспользует вектор). Возвращает true, если есть хотя
// бы один активный слот (contact && mu>0); флаг считается отдельным
// parallel_reduce по слотам (целочисленная сумма -- детерминированно).
// ----------------------------------------------------------------------
inline bool compute_fiber_contact_contributions(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const std::vector<fiber::FiberMaterial>& materials,
    const std::vector<FiberContactPrimitive>& primitives,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers,
    std::vector<FiberFrictionContact>& friction_contacts, const FiberContactStepParams& params = {})
{
    fiber_fiber_check_params(params, state);
    const std::size_t nworkers = exec.det_segments(); // детерминированно: не зависит от числа потоков
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const std::size_t n = primitives.size();
    friction_contacts.resize(n);
    const bool penalty = params.law.type != 0;
    std::vector<double>* pair_dt = penalty ? params.pair_dt_limit : nullptr;
    if (pair_dt) pair_dt->assign(n, std::numeric_limits<double>::max());
    if (n == 0) return false;

    exec.parallel_for_segments(n, [&](std::size_t seg, std::size_t k) {
        friction_contacts[k] = FiberFrictionContact{};
        const FiberContactPrimitive& prim = primitives[k];
        const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
        if (!geom.contact) return;

        // --- Нормальная сила: тип 0 -- Герц из геометрии, 1-3 -- закон ---------
        double F = geom.F;
        if (penalty)
        {
            const contact::Result res = fiber_fiber_contact_penalty(state, materials, prim, geom, params);
            if (pair_dt) (*pair_dt)[k] = res.dt_limit;
            F = res.fn;
            if (!(F > 0.0)) return; // покой/разлёт: ни вкладов, ни слота трения
        }

        // --- Слот для стадии трения --------------------------------------------
        const double mu = fiber_friction_coefficient(materials[static_cast<std::size_t>(state.material_id[prim.a0])],
                                                       materials[static_cast<std::size_t>(state.material_id[prim.b0])]);
        if (mu > 0.0)
        {
            FiberFrictionContact& slot = friction_contacts[k];
            slot.active = true;
            slot.a0 = prim.a0;
            slot.a1 = prim.a1;
            slot.b0 = prim.b0;
            slot.b1 = prim.b1;
            slot.s = geom.s;
            slot.t = geom.t;
            for (std::size_t d = 0; d < dim; ++d) slot.n[d] = geom.n[d];
            slot.N = F;
            slot.mu = mu;
        }

        rhs::ContributionWriter<fiber::FiberDelta> writer(thread_buffers[seg],
                                                            kFiberFiberContactInteractionType);

        // --- Сторона A: +F, веса (1-s, s) -------------------------------------
        if (geom.s < 1.0)
        {
            const double m = state.mass[prim.a0];
            const double wA = 1.0 - geom.s;
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = wA * F * geom.n[dd] / m;
            writer.emit(prim.a0, d, static_cast<std::uint64_t>(k), 0);
        }
        if (prim.a1 != -1 && geom.s > 0.0)
        {
            const rhs::index_t a1 = static_cast<rhs::index_t>(prim.a1);
            const double m = state.mass[a1];
            const double wA = geom.s;
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = wA * F * geom.n[dd] / m;
            writer.emit(a1, d, static_cast<std::uint64_t>(k), 1);
        }

        // --- Сторона B: -F, веса (1-t, t) -------------------------------------
        if (geom.t < 1.0)
        {
            const double m = state.mass[prim.b0];
            const double wB = 1.0 - geom.t;
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = -wB * F * geom.n[dd] / m;
            writer.emit(prim.b0, d, static_cast<std::uint64_t>(k), 2);
        }
        if (prim.b1 != -1 && geom.t > 0.0)
        {
            const rhs::index_t b1 = static_cast<rhs::index_t>(prim.b1);
            const double m = state.mass[b1];
            const double wB = geom.t;
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = -wB * F * geom.n[dd] / m;
            writer.emit(b1, d, static_cast<std::uint64_t>(k), 3);
        }
    });

    const std::size_t active_count = exec.parallel_reduce(
        n, std::size_t{0}, [&](std::size_t k) -> std::size_t { return friction_contacts[k].active ? 1u : 0u; },
        [](std::size_t a, std::size_t b) { return a + b; });
    return active_count > 0;
}

// Перегрузка без слотов трения (тесты нормального контакта): слоты во
// временном векторе, флаг активности отбрасывается.
inline void compute_fiber_contact_contributions(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const std::vector<fiber::FiberMaterial>& materials,
    const std::vector<FiberContactPrimitive>& primitives,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers, const FiberContactStepParams& params = {})
{
    std::vector<FiberFrictionContact> unused_slots;
    (void)compute_fiber_contact_contributions(exec, state, materials, primitives, thread_buffers, unused_slots,
                                              params);
}

// compute_fiber_contact -- удобная обёртка: search + contributions + reduce
// (только нормальная сила, без трения).
inline std::vector<fiber::FiberDelta> compute_fiber_contact(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const fiber::FiberTopology& topology,
    const std::vector<fiber::FiberMaterial>& materials,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers, const FiberContactStepParams& params = {})
{
    FiberContactSearchResult search;
    build_fiber_contact_search(exec, state, topology, materials, search);
    compute_fiber_contact_contributions(exec, state, materials, search.primitives, thread_buffers, params);
    return rhs::reduce_contributions(exec, thread_buffers, state.size());
}

// ----------------------------------------------------------------------
// fiber_contact_pair_dt_limit -- Δt_lim ОДНОГО primitive (пункт 6 в шапке
// файла), общая точка входа для compute_fiber_contact_dt_limit и
// диагностики (min_dt считается ТОЙ ЖЕ формулой, что и настоящий dt_limit,
// fiber/fiber_contact_segment_pairs.md, §7). m_red строится по ЭФФЕКТИВНЫМ
// массам точек контакта (fiber_contact_side_inv_mass), а не по "якорным"
// бусинам a0/b0, как в старой схеме -- теперь на пересечение приходится
// ровно одна сила, и реальная жёсткость узла совпадает с той, по которой
// считается Δt. +inf, если контакта нет или проникновение не положительно.
//
// ИЗВЕСТНОЕ ОГРАНИЧЕНИЕ (то же, что у Fiber-Solid-Contact и Solid-Solid-
// Contact): если одна бусина одновременно участвует в НЕСКОЛЬКИХ контактах
// (узлы плотной сетки, плетение), их жёсткости складываются, а эта оценка
// по-прежнему считается по ОДНОМУ контакту -- в редких сценариях с высокой
// кратностью контактов на бусину Δt может оказаться недостаточно
// консервативным.
// ----------------------------------------------------------------------
inline double fiber_contact_pair_dt_limit(const fiber::FiberState& state, const FiberContactPrimitive& prim,
                                           const FiberContactGeometry& geom)
{
    const double kInf = std::numeric_limits<double>::max();
    if (!geom.contact || !(geom.delta > 0.0)) return kInf;

    const double m_A = 1.0 / fiber_contact_side_inv_mass(state, prim.a0, prim.a1, geom.s);
    const double m_B = 1.0 / fiber_contact_side_inv_mass(state, prim.b0, prim.b1, geom.t);
    const double m_red = (m_A * m_B) / (m_A + m_B);
    const double k_eff = 1.5 * geom.F / geom.delta;
    if (!(k_eff > 0.0)) return kInf;
    return 0.8 / std::sqrt(k_eff / m_red);
}

// Законы 1-3: минимум res.dt_limit. Если params.pair_dt_limit заполнен ядром
// на тех же primitives (размер совпадает) -- редукция по слотам; иначе тот же
// расчёт (fiber_fiber_contact_penalty) повторяется по primitive. Тип 0 --
// прежний fiber_contact_pair_dt_limit.
inline double compute_fiber_contact_dt_limit(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                              const std::vector<fiber::FiberMaterial>& materials,
                                              const std::vector<FiberContactPrimitive>& primitives,
                                              const FiberContactStepParams& params = {})
{
    const double kInf = std::numeric_limits<double>::max();
    const std::size_t n = primitives.size();
    if (n == 0) return kInf;

    if (params.law.type != 0)
    {
        const std::vector<double>* slots = params.pair_dt_limit;
        if (slots != nullptr && slots->size() == n)
            return exec.parallel_reduce(
                n, kInf, [&](std::size_t k) -> double { return (*slots)[k]; },
                [](double a, double b) { return std::min(a, b); });
        fiber_fiber_check_params(params, state);
        return exec.parallel_reduce(
            n, kInf,
            [&](std::size_t k) -> double {
                const FiberContactPrimitive& prim = primitives[k];
                const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
                if (!geom.contact) return kInf;
                return fiber_fiber_contact_penalty(state, materials, prim, geom, params).dt_limit;
            },
            [](double a, double b) { return std::min(a, b); });
    }

    return exec.parallel_reduce(
        n, kInf,
        [&](std::size_t k) -> double {
            const FiberContactPrimitive& prim = primitives[k];
            const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
            return fiber_contact_pair_dt_limit(state, prim, geom);
        },
        [](double a, double b) { return std::min(a, b); });
}

// ----------------------------------------------------------------------
// Диагностика (fiber/fiber_contact_segment_pairs.md, §7). В НОВОЙ схеме
// «пара сегментов — одна сила» физическое пересечение двух нитей в общем
// случае даёт РОВНО ОДИН принятый контакт (contact_count), а не несколько,
// как в старой Meyer-подобной схеме, — это и есть цель перехода (см. шапку
// файла, пункт 1). pair_count/segment_pair_count позволяют отличить широкую
// фазу поиска (после фильтра object_id) от построенных primitives (после
// отбрасывания пар с последней бусиной стороны), rejected_by_ownership и
// degenerate_count -- сколько из построенных primitives отклонены правилом
// владения (§3) и сколько признаны геометрически вырожденными (§4)
// соответственно.
//   pair_count            -- spatial-пары после фильтра object_id[i]!=
//                             object_id[j] (широкая фаза, из search.pair_count);
//   segment_pair_count    -- построенные primitives (search.primitives.size());
//   rejected_by_ownership -- primitives, отклонённые правилом владения;
//   degenerate_count      -- primitives с совпадающими/неопределёнными осями
//                             (D~0, контакт не сформирован);
//   contact_count         -- принятые контакты с фактическим проникновением
//                             (D < d_eff);
//   max_penetration, min_dt -- как в SolidContactDiagnostics
//                             (interactions/solid_contact.hpp), min_dt --
//                             формулой fiber_contact_pair_dt_limit выше.
// ----------------------------------------------------------------------
struct FiberContactDiagnostics
{
    std::size_t pair_count = 0;
    std::size_t segment_pair_count = 0;
    std::size_t rejected_by_ownership = 0;
    std::size_t degenerate_count = 0;
    std::size_t contact_count = 0;
    double max_penetration = 0.0;
    double min_dt = std::numeric_limits<double>::max();
};

inline FiberContactDiagnostics compute_fiber_contact_diagnostics(const rhs::CpuExecutor& exec,
                                                                   const fiber::FiberState& state,
                                                                   const std::vector<fiber::FiberMaterial>& materials,
                                                                   const FiberContactSearchResult& search,
                                                                   const FiberContactStepParams& params = {})
{
    fiber_fiber_check_params(params, state);
    FiberContactDiagnostics result;
    result.pair_count = search.pair_count;
    result.segment_pair_count = search.primitives.size();

    const std::size_t n = search.primitives.size();
    if (n == 0) return result;

    auto combine = [](FiberContactDiagnostics a, FiberContactDiagnostics b) {
        FiberContactDiagnostics r;
        r.pair_count = a.pair_count + b.pair_count; // остаются 0 в свёртке -- проставляются после reduce
        r.segment_pair_count = a.segment_pair_count + b.segment_pair_count;
        r.rejected_by_ownership = a.rejected_by_ownership + b.rejected_by_ownership;
        r.degenerate_count = a.degenerate_count + b.degenerate_count;
        r.contact_count = a.contact_count + b.contact_count;
        r.max_penetration = std::max(a.max_penetration, b.max_penetration);
        r.min_dt = std::min(a.min_dt, b.min_dt);
        return r;
    };

    FiberContactDiagnostics init;
    FiberContactDiagnostics acc = exec.parallel_reduce(
        n, init,
        [&](std::size_t k) -> FiberContactDiagnostics {
            FiberContactDiagnostics d1;
            const FiberContactPrimitive& prim = search.primitives[k];
            const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
            if (geom.rejected_by_ownership)
            {
                d1.rejected_by_ownership = 1;
            }
            else if (geom.degenerate)
            {
                d1.degenerate_count = 1;
            }
            else if (geom.contact)
            {
                d1.contact_count = 1;
                d1.max_penetration = geom.delta;
                // min_dt -- той же формулой, что настоящий предел шага.
                d1.min_dt = (params.law.type == 0)
                                ? fiber_contact_pair_dt_limit(state, prim, geom)
                                : fiber_fiber_contact_penalty(state, materials, prim, geom, params).dt_limit;
            }
            return d1;
        },
        combine);

    acc.pair_count = search.pair_count;
    acc.segment_pair_count = search.primitives.size();
    return acc;
}

} // namespace interactions
