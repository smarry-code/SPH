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
                                           const std::vector<fiber::FiberMaterial>& materials,
                                           const double* l_max = nullptr)
{
    double d_max = 0.0;
    for (const auto& m : materials) d_max = std::max(d_max, m.d);

    // l_max != nullptr -- L_max уже посчитан вызывающей стороной на этом шаге (то же значение).
    const double L_max = l_max ? *l_max : fiber_max_current_bond_length(exec, state, topology);

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
// Определена ниже (консервативный отсев primitive); нужна build в режиме отсечки.
inline bool fiber_contact_surely_apart(const fiber::FiberState& state,
                                        const std::vector<fiber::FiberMaterial>& materials,
                                        const FiberContactPrimitive& prim);

// Выжившая пара в сборке: primitive и её номер kl внутри отрезка ячеек среди ВСЕХ построенных primitive
// (до отсечки). Глобальный индекс k = смещение отрезка + kl.
struct FiberContactSurvivor
{
    FiberContactPrimitive prim;
    std::size_t kl = 0;
};

// Счётчики отрезка ячеек: alignas(64) -- каждый отрезок пишет свои, без ложного разделения кэш-линий.
struct alignas(64) FiberContactSegCounts
{
    std::size_t all = 0;   // построенные primitive отрезка (до отсечки)
    std::size_t pairs = 0; // spatial-пары отрезка после фильтра object_id
};

struct FiberContactBuildOptions
{
    // true -- в primitives остаются только пары, не прошедшие fiber_contact_surely_apart (остальные не дают ни
    // контакта, ни вкладов, ни предела шага); индекс каждой выжившей в полном списке -- prim_k.
    bool cull_apart = false;
    const double* l_max = nullptr; // L_max уже посчитан вызывающей стороной (то же значение), иначе считается здесь
};

struct FiberContactSearchResult
{
    // Полный режим: все построенные primitive, индекс j == k. Режим отсечки (culled): только выжившие,
    // j-й выживший имеет полный индекс k = prim_k[j] (строго возрастает).
    std::vector<FiberContactPrimitive> primitives;
    std::size_t pair_count = 0; // spatial-пары после фильтра object_id (широкая фаза, до отбрасывания по таблице сторон)
    // Постоянная рабочая память поиска: результат живёт между шагами, поэтому индекс не создаётся заново на
    // каждом вызове (память только растёт). Вектор пар не строится: пары сразу превращаются в primitive
    // (SpatialIndex::for_each_pair_segmented).
    rhs::SpatialIndex<dim> index;
    // Предел шага на primitive (законы 1-3, FiberContactStepParams::pair_dt_limit):
    // заполняет ядро, минимум берёт compute_fiber_contact_dt_limit. Размер -- primitives.size() (плотный по j).
    std::vector<double> pair_dt_limit;

    // Режим отсечки. k = индекс primitive в полном (неотсечённом) списке -- он задаёт разбиение на отрезки
    // kDetSegments и порядок суммирования, поэтому нужен ядру, даже когда сами отсечённые пары не хранятся.
    std::vector<std::size_t> prim_k; // пуст в полном режиме
    std::size_t n_full = 0;          // число primitive до отсечки
    // seg_bound[s]..seg_bound[s+1] -- j, чей k лежит в отрезке s разбиения parallel_for_segments(n_full):
    // в полном режиме n_full*s/64, в режиме отсечки lower_bound(prim_k, n_full*s/64); seg_bound[64] = primitives.size().
    std::size_t seg_bound[rhs::CpuExecutor::kDetSegments + 1] = {};
    bool culled = false;

    // Рабочая память сборки: выжившие по отрезкам ячеек и счётчики отрезков.
    std::vector<rhs::CpuExecutor::SegPart<FiberContactSurvivor>> surv_parts;
    FiberContactSegCounts seg_cnt[rhs::CpuExecutor::kDetSegments];
};

// Пара бусин (ka, kb) -> primitive по правилу «сегмент принадлежит первой бусине» (шапка файла, пункт 1).
// false -- пара отбрасывается целиком (одна из сторон -- последняя бусина нити).
inline bool make_fiber_contact_prim(const fiber::FiberTopology& topology, rhs::index_t k_a, rhs::index_t k_b,
                                     FiberContactPrimitive& out)
{
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
        return false; // k_a -- последняя бусина своей нити: пара отбрасывается целиком
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
        return false;
    }

    out.a0 = k_a;
    out.a1 = a1;
    out.a_prev = a_prev;
    out.a_has_next_segment = (a1 != -1) && (topology.next[static_cast<std::size_t>(a1)] != -1);
    out.b0 = k_b;
    out.b1 = b1;
    out.b_prev = b_prev;
    out.b_has_next_segment = (b1 != -1) && (topology.next[static_cast<std::size_t>(b1)] != -1);
    return true;
}

// Сборка без вектора пар. Пары обходятся SpatialIndex::for_each_pair_segmented в том же порядке, в каком
// их выдавал find_pairs_parallel: отрезок ячеек seg по возрастанию, внутри -- порядок обхода. Каждый отрезок
// строит свои primitive (kl -- номер внутри отрезка среди ВСЕХ построенных, даже если пара затем отсечена),
// смещения -- префиксная сумма по 64 отрезкам. Поэтому полный индекс k = off_all[seg] + kl совпадает с
// индексом после прежней склейки concat_segments, а при cull_apart == false primitives и порядок -- те же.
// Отсечённые пары в primitives не попадают, но их k учтён в n_full и prim_k выживших.
inline void build_fiber_contact_search(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                        const fiber::FiberTopology& topology,
                                        const std::vector<fiber::FiberMaterial>& materials,
                                        FiberContactSearchResult& out, const FiberContactBuildOptions& opt)
{
    constexpr std::size_t S = rhs::CpuExecutor::kDetSegments;
    const bool cull = opt.cull_apart;
    // primitives очищается только на ранних выходах: при построении он перезаписывается целиком
    // (resize инициализирует лишь прирост).
    out.pair_count = 0;
    out.n_full = 0;
    out.culled = cull;
    out.prim_k.clear();
    for (std::size_t s = 0; s <= S; ++s) out.seg_bound[s] = 0;
    const std::size_t n = state.size();
    if (n == 0)
    {
        out.primitives.clear();
        return;
    }

    if constexpr (dim != 3)
    {
        // Нити физически существуют только при dim==3 (fiber/SPH.md; Р9);
        // граница проверяется fiber::check_fiber_requires_dim3() при разборе
        // конфигурации -- сюда с непустым state в сборке dim!=3 попасть
        // невозможно. Заглушка нужна только для компиляции заголовка во ВСЕХ
        // SPH_DIM-вариантах пайплайна (см. include/interactions/fiber_bond.hpp,
        // include/fiber/fiber_system.hpp -- тот же паттерн).
        out.primitives.clear();
        return;
    }
    else
    {
        const double radius = fiber_contact_search_radius(exec, state, topology, materials, opt.l_max);

        std::array<const double*, dim> coords{};
        for (std::size_t d = 0; d < dim; ++d) coords[d] = state.x[d].data();
        out.index.build(exec, coords, n, radius > 0.0 ? radius : 1.0);

        auto accept = [&](rhs::index_t a, rhs::index_t b) -> bool {
            return state.object_id[a] != state.object_id[b];
        };

        exec.prepare_segment_parts(out.surv_parts);
        for (std::size_t s = 0; s < S; ++s) out.seg_cnt[s] = FiberContactSegCounts{};
        out.index.for_each_pair_segmented(exec, radius, accept, [&](std::size_t seg, std::size_t ka, std::size_t kb) {
            FiberContactSegCounts& cnt = out.seg_cnt[seg];
            ++cnt.pairs;
            FiberContactSurvivor sv;
            if (!make_fiber_contact_prim(topology, static_cast<rhs::index_t>(ka), static_cast<rhs::index_t>(kb),
                                         sv.prim))
                return;
            sv.kl = cnt.all++;
            if (cull && fiber_contact_surely_apart(state, materials, sv.prim)) return;
            out.surv_parts[seg].v.push_back(sv);
        });

        // Префиксные суммы по 64 отрезкам (последовательно): смещения в полном списке и в списке выживших.
        std::size_t off_all[S + 1], off_surv[S + 1];
        off_all[0] = 0;
        off_surv[0] = 0;
        for (std::size_t s = 0; s < S; ++s)
        {
            off_all[s + 1] = off_all[s] + out.seg_cnt[s].all;
            off_surv[s + 1] = off_surv[s] + out.surv_parts[s].v.size();
            out.pair_count += out.seg_cnt[s].pairs;
        }
        out.n_full = off_all[S];
        const std::size_t m = off_surv[S];
        out.primitives.resize(m);
        if (cull) out.prim_k.resize(m);
        exec.parallel_for(S, std::size_t(1), [&](std::size_t s) {
            const auto& v = out.surv_parts[s].v;
            for (std::size_t i = 0; i < v.size(); ++i)
            {
                out.primitives[off_surv[s] + i] = v[i].prim;
                if (cull) out.prim_k[off_surv[s] + i] = off_all[s] + v[i].kl;
            }
        });

        // Границы отрезков разбиения parallel_for_segments(n_full) по j.
        for (std::size_t s = 0; s < S; ++s)
        {
            const std::size_t b = rhs::CpuExecutor::seg_begin(out.n_full, s, S);
            out.seg_bound[s] =
                cull ? static_cast<std::size_t>(std::lower_bound(out.prim_k.begin(), out.prim_k.end(), b) -
                                                 out.prim_k.begin())
                     : b;
        }
        out.seg_bound[S] = m;
    }
}

inline void build_fiber_contact_search(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                        const fiber::FiberTopology& topology,
                                        const std::vector<fiber::FiberMaterial>& materials,
                                        FiberContactSearchResult& out, const double* l_max = nullptr)
{
    FiberContactBuildOptions opt;
    opt.cull_apart = false;
    opt.l_max = l_max;
    build_fiber_contact_search(exec, state, topology, materials, out, opt);
}

// Вид на набор primitive для ядер: m primitive, j-й имеет полный индекс k_of(j) (k == nullptr -- j == k),
// bound -- границы отрезков по j (см. FiberContactSearchResult::seg_bound); в полном виде (culled == false)
// ядра считают границы сами (m*s/64), указатель bound не используется и может быть nullptr -- так вид
// остаётся верным, даже если вызывающий изменил primitives после сборки.
struct FiberContactPrimView
{
    const FiberContactPrimitive* prim = nullptr;
    std::size_t m = 0;
    const std::size_t* k = nullptr;
    std::size_t n_full = 0;
    const std::size_t* bound = nullptr;
    bool culled = false;

    // Пользовательский конструктор: вид не агрегат, поэтому вызов вида f(..., {prim}) с braced-списком
    // однозначно выбирает перегрузку по std::vector<FiberContactPrimitive>.
    FiberContactPrimView() {}

    std::size_t k_of(std::size_t j) const { return k ? k[j] : j; }

    // Границы отрезков по j: в режиме отсечки bound, иначе m*s/64 в local.
    const std::size_t* seg_bounds(std::size_t* local) const
    {
        if (culled) return bound;
        for (std::size_t s = 0; s <= rhs::CpuExecutor::kDetSegments; ++s)
            local[s] = rhs::CpuExecutor::seg_begin(m, s);
        return local;
    }
};

inline FiberContactPrimView fiber_contact_view(const FiberContactSearchResult& r)
{
    FiberContactPrimView v;
    v.prim = r.primitives.data();
    v.m = r.primitives.size();
    v.culled = r.culled;
    v.k = r.culled ? r.prim_k.data() : nullptr;
    v.n_full = r.culled ? r.n_full : r.primitives.size();
    v.bound = r.seg_bound;
    return v;
}

// Полный вид над вектором (старые вызовы): j == k.
inline FiberContactPrimView fiber_contact_view(const std::vector<FiberContactPrimitive>& primitives)
{
    FiberContactPrimView v;
    v.prim = primitives.data();
    v.m = primitives.size();
    v.n_full = primitives.size();
    return v;
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

// ----------------------------------------------------------------------
// fiber_contact_surely_apart -- дешёвый КОНСЕРВАТИВНЫЙ отсев primitive: true, только если контакт
// заведомо невозможен (тогда compute_fiber_contact_geometry вернула бы contact == false).
// Расстояние между отрезками >= |m_a - m_b| - (L_a + L_b)/2 (m -- середины, L -- длины; точечный
// сегмент: a1 == -1, L = 0), а контакт требует D < d_eff = 0.5*(d_a + d_b) -- те же материалы по
// якорным бусинам a0/b0, что в геометрии. Радиус R = d_eff + (L_a + L_b)/2 берётся с запасом
// 1e-8*(d_eff + L_a + L_b) на snapping s,t (kFiberContactSnapEta) и округление. Некорректный индекс
// материала -> false (решает полная геометрия). NaN в координатах -> сравнение ложно -> false.
// ----------------------------------------------------------------------
inline bool fiber_contact_surely_apart(const fiber::FiberState& state,
                                        const std::vector<fiber::FiberMaterial>& materials,
                                        const FiberContactPrimitive& prim)
{
    if constexpr (dim != 3)
    {
        (void)state;
        (void)materials;
        (void)prim;
        return false;
    }
    else
    {
        const long long ma_id = static_cast<long long>(state.material_id[prim.a0]);
        const long long mb_id = static_cast<long long>(state.material_id[prim.b0]);
        const long long nmat = static_cast<long long>(materials.size());
        if (ma_id < 0 || ma_id >= nmat || mb_id < 0 || mb_id >= nmat) return false;
        const double d_eff = 0.5 * materials[static_cast<std::size_t>(ma_id)].d +
                             0.5 * materials[static_cast<std::size_t>(mb_id)].d;

        double xa0[3], xa1[3], xb0[3], xb1[3];
        for (std::size_t d = 0; d < 3; ++d)
        {
            xa0[d] = state.x[d][prim.a0];
            xa1[d] = (prim.a1 == -1) ? xa0[d] : state.x[d][static_cast<rhs::index_t>(prim.a1)];
            xb0[d] = state.x[d][prim.b0];
            xb1[d] = (prim.b1 == -1) ? xb0[d] : state.x[d][static_cast<rhs::index_t>(prim.b1)];
        }
        double la2 = 0.0, lb2 = 0.0, dm2 = 0.0;
        double xabs = 0.0; // max |x| по четырём точкам: масштаб абсолютной ошибки округления
        for (std::size_t d = 0; d < 3; ++d)
        {
            xabs = std::max({xabs, std::fabs(xa0[d]), std::fabs(xa1[d]), std::fabs(xb0[d]), std::fabs(xb1[d])});
            const double ea = xa1[d] - xa0[d];
            const double eb = xb1[d] - xb0[d];
            la2 += ea * ea;
            lb2 += eb * eb;
            const double dm = 0.5 * (xa0[d] + xa1[d]) - 0.5 * (xb0[d] + xb1[d]);
            dm2 += dm * dm;
        }
        const double la = std::sqrt(la2);
        const double lb = std::sqrt(lb2);
        double R = d_eff + 0.5 * (la + lb);
        // Запас: относительный -- на snapping s,t и округление в масштабе отрезков; абсолютный -- на округление
        // середин, разностей и расстояния в геометрии при больших координатах (ошибка ~ ulp(max|x|), не зависит
        // от длин: при |x| ~ 2^52 и длинах ~ 1 без него отсев ложен). 64*eps*max|x| покрывает обе стороны оценки.
        R += 1.0e-8 * (d_eff + la + lb) + 64.0 * std::numeric_limits<double>::epsilon() * xabs;
        return std::isfinite(R) && dm2 > R * R;
    }
}

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
// friction_contacts.resize(m) (m -- число primitive вида) выполняется в начале, ДО
// parallel_for (внутри kernel аллокаций нет; память только растёт, если
// вызывающая сторона переиспользует вектор).
//
// Вид FiberContactPrimView (режим отсечки поиска): primitive -- только выжившие после
// fiber_contact_surely_apart, j-й имеет полный индекс k = view.k_of(j). Плотные массивы (слоты трения,
// pair_dt_limit, нормальные слоты c.prim) индексируются по j. Побитовость результата держится на том, что
// номер отрезка seg вкладов -- отрезок ПОЛНОГО индекса k в разбиении n_full (его даёт view.bound), а не
// отрезок j: буферы вкладов, dt_seg_min и нормальные слоты видят тот же seg и тот же порядок (по
// возрастанию k), что и при обходе полного списка с пропуском отсечённых, которые ничего не выпускают.
// Выжившие уже прошли fiber_contact_surely_apart при сборке, повторно её вызывают только в полном виде. Возвращает true, если есть хотя
// бы один активный слот (contact && mu>0); флаг считается отдельным
// parallel_reduce по слотам (целочисленная сумма -- детерминированно).
// ----------------------------------------------------------------------
// Определена ниже (предел шага одного primitive); ядро вызывает её для закона 0 (law0_dt_min_out).
inline double fiber_contact_pair_dt_limit(const fiber::FiberState& state, const FiberContactPrimitive& prim,
                                          const FiberContactGeometry& geom);

inline bool compute_fiber_contact_contributions(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const std::vector<fiber::FiberMaterial>& materials,
    const FiberContactPrimView& view,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers,
    std::vector<FiberFrictionContact>& friction_contacts, const FiberContactStepParams& params = {})
{
    fiber_fiber_check_params(params, state);
    const std::size_t nworkers = exec.det_segments(); // детерминированно: не зависит от числа потоков
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const std::size_t m = view.m;
    friction_contacts.resize(m);
    const bool penalty = params.law.type != 0;
    std::vector<double>* pair_dt = penalty ? params.pair_dt_limit : nullptr;
    if (pair_dt) pair_dt->resize(m); // слот j инициализируется в начале тела цикла (покрывает все j)
    // Режим проекции: нормальные слоты только активных контактов (delta > 0), сжатым массивом по возрастанию
    // primitive: каждый отрезок пишет свои слоты в local[seg] (по возрастанию k), затем отрезки склеиваются
    // по порядку. prim слота -- индекс j в плотных массивах вида, seg -- номер отрезка k (см. FiberNormalContact).
    std::vector<FiberNormalContact>* const normal = penalty ? params.normal_ff : nullptr;
    if (m == 0)
    {
        if (normal) normal->clear();
        if (params.fric_active) params.fric_active->clear();
        if (!penalty && params.law0_dt_min_out) *params.law0_dt_min_out = std::numeric_limits<double>::max();
        return false;
    }
    std::size_t bound_local[rhs::CpuExecutor::kDetSegments + 1];
    const std::size_t* const bound = view.seg_bounds(bound_local);
    std::vector<rhs::CpuExecutor::SegPart<FiberNormalContact>> parts_tmp;
    auto& normal_local = params.normal_ff_parts ? *params.normal_ff_parts : parts_tmp;
    if (normal) exec.prepare_segment_parts(normal_local);

    // Режим списка слотов трения: сброс только слотов прошлого вызова (вместо сброса каждого j в цикле).
    // Инвариант: после resize(m) все слоты, кроме перечисленных в прошлом списке, пусты ({}).
    std::vector<FrictionActiveRef>* const fric_list = params.fric_active;
    std::vector<rhs::CpuExecutor::SegPart<FrictionActiveRef>> fric_parts_tmp;
    auto& fric_local = params.fric_active_parts ? *params.fric_active_parts : fric_parts_tmp;
    if (fric_list)
    {
        exec.parallel_for(fric_list->size(), [&](std::size_t j) {
            const std::size_t kk = (*fric_list)[j].k;
            if (kk < m) friction_contacts[kk] = FiberFrictionContact{};
        });
        exec.prepare_segment_parts(fric_local);
    }

    // Предел шага закона 0: свёртка min по отрезкам, как в parallel_reduce compute_fiber_contact_dt_limit
    // (неконтактные элементы дают max() -- нейтральный элемент, поэтому достаточно принятых контактов).
    const double kInfDt = std::numeric_limits<double>::max();
    double dt_seg_min[rhs::CpuExecutor::kDetSegments];
    double* const law0_dt_out = (!penalty) ? params.law0_dt_min_out : nullptr;
    if (law0_dt_out)
        for (std::size_t s = 0; s < rhs::CpuExecutor::kDetSegments; ++s) dt_seg_min[s] = kInfDt;

    exec.parallel_for_seg_bounds(bound, [&](std::size_t seg, std::size_t j) {
        const std::size_t k = view.k_of(j); // полный индекс: ключ meta вкладов
        if (pair_dt) (*pair_dt)[j] = std::numeric_limits<double>::max();
        if (!fric_list) friction_contacts[j] = FiberFrictionContact{};
        const FiberContactPrimitive& prim = view.prim[j];
        if (!view.culled && fiber_contact_surely_apart(state, materials, prim)) return;
        const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
        if (!geom.contact) return;

        // --- Режим проекции: только нормальный слот, силы и слота трения нет ---
        if (normal)
        {
            if (!(geom.delta > 0.0)) return;
            FiberNormalContact slot;
            slot.active = true;
            slot.prim = j;
            slot.seg = static_cast<std::uint16_t>(seg);
            slot.a0 = prim.a0;
            slot.a1 = prim.a1;
            slot.b0 = prim.b0;
            slot.b1 = prim.b1;
            slot.s = geom.s;
            slot.t = geom.t;
            for (std::size_t d = 0; d < dim; ++d) slot.n[d] = geom.n[d];
            slot.delta = geom.delta;
            slot.mu = fiber_friction_coefficient(materials[static_cast<std::size_t>(state.material_id[prim.a0])],
                                                 materials[static_cast<std::size_t>(state.material_id[prim.b0])]);
            slot.key = make_ff_pair_key(prim.a0, prim.a1, prim.b0, prim.b1);
            normal_local[seg].v.push_back(slot);
            if (pair_dt) (*pair_dt)[j] = contact::kNoLimit;
            return;
        }

        // --- Нормальная сила: тип 0 -- Герц из геометрии, 1-3 -- закон ---------
        double F = geom.F;
        if (law0_dt_out)
            dt_seg_min[seg] = std::min(dt_seg_min[seg], fiber_contact_pair_dt_limit(state, prim, geom));
        if (penalty)
        {
            const contact::Result res = fiber_fiber_contact_penalty(state, materials, prim, geom, params);
            if (pair_dt) (*pair_dt)[j] = res.dt_limit;
            F = res.fn;
            if (!(F > 0.0)) return; // покой/разлёт: ни вкладов, ни слота трения
        }

        // --- Слот для стадии трения --------------------------------------------
        const double mu = fiber_friction_coefficient(materials[static_cast<std::size_t>(state.material_id[prim.a0])],
                                                       materials[static_cast<std::size_t>(state.material_id[prim.b0])]);
        if (mu > 0.0)
        {
            FiberFrictionContact& slot = friction_contacts[j];
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
            if (fric_list) fric_local[seg].v.push_back(FrictionActiveRef{j, static_cast<std::uint16_t>(seg)});
        }

        rhs::ContributionWriter<fiber::FiberDelta> writer(thread_buffers[seg],
                                                            kFiberFiberContactInteractionType);

        // --- Сторона A: +F, веса (1-s, s) -------------------------------------
        if (geom.s < 1.0)
        {
            const double mass = state.mass[prim.a0];
            const double wA = 1.0 - geom.s;
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = wA * F * geom.n[dd] / mass;
            writer.emit(prim.a0, d, static_cast<std::uint64_t>(k), 0);
        }
        if (prim.a1 != -1 && geom.s > 0.0)
        {
            const rhs::index_t a1 = static_cast<rhs::index_t>(prim.a1);
            const double mass = state.mass[a1];
            const double wA = geom.s;
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = wA * F * geom.n[dd] / mass;
            writer.emit(a1, d, static_cast<std::uint64_t>(k), 1);
        }

        // --- Сторона B: -F, веса (1-t, t) -------------------------------------
        if (geom.t < 1.0)
        {
            const double mass = state.mass[prim.b0];
            const double wB = 1.0 - geom.t;
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = -wB * F * geom.n[dd] / mass;
            writer.emit(prim.b0, d, static_cast<std::uint64_t>(k), 2);
        }
        if (prim.b1 != -1 && geom.t > 0.0)
        {
            const rhs::index_t b1 = static_cast<rhs::index_t>(prim.b1);
            const double mass = state.mass[b1];
            const double wB = geom.t;
            fiber::FiberDelta d{};
            for (std::size_t dd = 0; dd < dim; ++dd) d.dv[dd] = -wB * F * geom.n[dd] / mass;
            writer.emit(b1, d, static_cast<std::uint64_t>(k), 3);
        }
    });

    if (normal) exec.concat_segments(normal_local, *normal);

    if (law0_dt_out)
    {
        double res = kInfDt; // init свёртки parallel_reduce, затем отрезки по порядку
        for (std::size_t s = 0; s < rhs::CpuExecutor::kDetSegments; ++s) res = std::min(res, dt_seg_min[s]);
        *law0_dt_out = res;
    }

    if (fric_list)
    {
        exec.concat_segments(fric_local, *fric_list);
        return !fric_list->empty();
    }

    const std::size_t active_count = exec.parallel_reduce(
        m, std::size_t{0}, [&](std::size_t j) -> std::size_t { return friction_contacts[j].active ? 1u : 0u; },
        [](std::size_t a, std::size_t b) { return a + b; });
    return active_count > 0;
}

// Перегрузки по вектору primitive (полный вид, j == k): тесты и прежние вызовы.
inline bool compute_fiber_contact_contributions(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const std::vector<fiber::FiberMaterial>& materials,
    const std::vector<FiberContactPrimitive>& primitives,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers,
    std::vector<FiberFrictionContact>& friction_contacts, const FiberContactStepParams& params = {})
{
    return compute_fiber_contact_contributions(exec, state, materials, fiber_contact_view(primitives), thread_buffers,
                                               friction_contacts, params);
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
// на тех же primitives (размер совпадает) -- редукция по слотам; иначе тот
// же расчёт (fiber_fiber_contact_penalty) повторяется по primitive. Тип 0 --
// прежний fiber_contact_pair_dt_limit.
// Вид (FiberContactPrimView): свёртка идёт по границам отрезков ПОЛНОГО индекса k (view.bound), тот же
// порядок min по отрезкам, что при обходе полного списка; отсечённые пары дают +inf (нейтральны для min).
// В полном виде закон 0 сперва вызывает fiber_contact_surely_apart, в режиме отсечки -- нет (выжившие
// её уже прошли; законы 1-4 её не вызывали и раньше).
inline double compute_fiber_contact_dt_limit(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                              const std::vector<fiber::FiberMaterial>& materials,
                                              const FiberContactPrimView& view,
                                              const FiberContactStepParams& params = {})
{
    const double kInf = std::numeric_limits<double>::max();
    const std::size_t m = view.m;
    if (m == 0) return kInf;
    std::size_t bound_local[rhs::CpuExecutor::kDetSegments + 1];
    const std::size_t* const bound = view.seg_bounds(bound_local);
    auto fmin = [](double a, double b) { return std::min(a, b); };

    if (params.law.type != 0)
    {
        const std::vector<double>* slots = params.pair_dt_limit;
        if (slots != nullptr && slots->size() == m)
            return exec.parallel_reduce_seg_bounds(
                bound, kInf, [&](std::size_t j) -> double { return (*slots)[j]; }, fmin);
        fiber_fiber_check_params(params, state);
        return exec.parallel_reduce_seg_bounds(
            bound, kInf,
            [&](std::size_t j) -> double {
                const FiberContactPrimitive& prim = view.prim[j];
                const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
                if (!geom.contact) return kInf;
                return fiber_fiber_contact_penalty(state, materials, prim, geom, params).dt_limit;
            },
            fmin);
    }

    return exec.parallel_reduce_seg_bounds(
        bound, kInf,
        [&](std::size_t j) -> double {
            const FiberContactPrimitive& prim = view.prim[j];
            if (!view.culled && fiber_contact_surely_apart(state, materials, prim)) return kInf;
            const FiberContactGeometry geom = compute_fiber_contact_geometry(state, materials, prim);
            return fiber_contact_pair_dt_limit(state, prim, geom);
        },
        fmin);
}

inline double compute_fiber_contact_dt_limit(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                              const std::vector<fiber::FiberMaterial>& materials,
                                              const std::vector<FiberContactPrimitive>& primitives,
                                              const FiberContactStepParams& params = {})
{
    return compute_fiber_contact_dt_limit(exec, state, materials, fiber_contact_view(primitives), params);
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
    // Диагностика считает по ВСЕМ построенным primitive (segment_pair_count, rejected_by_ownership, degenerate_count);
    // в режиме отсечки их часть не хранится -- результат был бы молча неверным.
    if (search.culled)
        throw std::logic_error("compute_fiber_contact_diagnostics: поиск построен в режиме отсечки (cull_apart)");
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
