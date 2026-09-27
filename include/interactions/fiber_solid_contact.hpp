#pragma once
// interactions/fiber_solid_contact.hpp
// Fiber-Solid-Contact: механический (Hertz penalty) контакт бусины нити с
// РЕАЛЬНОЙ Solid-частицей (fiber/SPH.md, раздел «Fiber-Solid-Contact»;
// fiber/PLAN.md, этап Ф6, «3.2. Типы взаимодействий и primitives», «3.3»,
// «3.7»).
//
// Ф4 (перевод Solid на новую RHS, SolidDelta + rhs::reduce_contributions) на
// момент реализации Ф6 УЖЕ ЗАВЕРШЁН, поэтому вклад в Solid здесь формируется
// строго тем же способом, что Solid-SPH/Solid-Solid-Contact — как
// interactions::SolidDelta через rhs::ContributionWriter<SolidDelta>, а НЕ
// временным прибавлением к ускорению внутри single_step (развилка
// PLAN.md была на случай Ф6-раньше-Ф4, здесь не актуальна).
//
// ============================================================================
// РЕШЕНИЯ (задокументированы явно, т.к. постановка оставляла их на усмотрение
// реализации — fiber/PLAN.md, этап Ф6, последний абзац):
//
// 1. Поиск Fiber -> Solid использует rhs::SpatialIndex::find_pairs_cross
//    (rhs/spatial_index.hpp, ИМЕННО для этого случая: "Fiber -> Solid само
//    обеспечивает однократную генерацию пары"). Индексируется («this»)
//    подмножество РЕАЛЬНЫХ Solid-частиц, для которых elastic_by_itype()
//    определён (materials.hpp) — не индексировать частицы, для которых
//    Fiber-Solid-Contact заведомо не формируется (газ/вода, itype 1-2),
//    дешевле, чем отбрасывать их в kernel (fiber/PLAN.md, этап Ф6, пункт 4:
//    "предпочтительно фильтровать раньше"). Запрос («other») — все бусины
//    Fiber. Результат — пары (индекс бусины Fiber, локальный индекс в
//    подмножестве Solid), однократно, без дублей и без сортировки.
//
// 2. Радиус поиска -- r_f_max + r_s_max^contact + L_max, где L_max -- ТЕКУЩАЯ
//    (не номинальная l0!) максимальная длина связи topology.bonds
//    (interactions/fiber_contact_common.hpp::fiber_max_current_bond_length,
//    общая с Fiber-Fiber-Contact, fiber/fiber_solid_contact_fix.md, §2).
//    Связи Fiber-Bond работают только на растяжение и рвутся лишь в
//    Postprocessing ПОСЛЕ шага (fiber/SPH.md), поэтому внутри шага сегмент
//    может быть длиннее l0 -- номинальное l0 как оценка не годится (иначе
//    контакт у дальнего конца растянутого сегмента может выпасть из поиска,
//    fiber_solid_contact_fix.md, §2). "Правило «сегмент принадлежит первой
//    бусине»" (PLAN.md §3.2) по-прежнему реализовано на этапе ПОСТРОЕНИЯ
//    primitive, а не в kernel:
//      - next(f) существует  -> primitive-сегмент (f, next(f), s);
//      - next(f) НЕ существует (f -- последняя бусина нити):
//          - prev(f) существует  -> primitive НЕ строится вовсе: сегмент
//            (prev(f), f) уже покрывает эту зону контакта, ОН СТРОИТСЯ, когда
//            поиск встречает пару (prev(f), s), — при радиусе поиска
//            r_f_max + r_s_max^contact + L_max это ГАРАНТИРОВАНО происходит
//            всегда, когда контакт вообще возможен рядом с f (неравенство
//            треугольника: |x_s-x_{prev(f)}| <= |x_s-x_f| + |x_f-x_{prev(f)}|
//            <= (r_f+r_s) + |x_f-x_{prev(f)}| <= r_f_max+r_s_max^contact+L_max,
//            т.к. |x_f-x_{prev(f)}| -- ТЕКУЩАЯ длина связи (prev(f),f) и,
//            значит, <= L_max = максимуму по ВСЕМ связям расчёта -- см. вывод
//            ниже, у fiber_solid_contact_search_radius);
//          - prev(f) тоже НЕ существует (f -- изолированная бусина без
//            соседей, PLAN.md: "одиночная бусина без соседей даёт точечный
//            контакт") -> primitive-точка (f, s), сила целиком на f.
//    Это исключает ситуацию, разобранную в задании (двойной учёт через
//    ОДНОВРЕМЕННОЕ появление сегмента (prev(f),f) с λ=1 у последнего сегмента
//    И отдельного точечного контакта в f): точечный примитив в f строится
//    ТОЛЬКО когда сегмент (prev(f),f) в принципе не существует.
//
// 3. Правило владения узлом (fiber/fiber_solid_contact_fix.md, §1) --
//    заменяет старое "полуинтервал λ∈[0,1) против [0,1] по СЫРОМУ параметру
//    t", которое теряло контакт целиком, если нить изломана в бусине a, а
//    Solid-частица лежит во ВНЕШНЕМ секторе излома (оба соседних сегмента
//    отбрасывали контакт, т.к. |t|>=1 и |t|<0 одновременно у "своего" конца).
//    Параметр t вычисляется ЕДИНОЙ функцией fiber_solid_segment_param(xa,xb,xs)
//    -- используется и для "своего" сегмента (a,b), и для соседнего
//    (prev(a),a) -- с прижатием к границе (η=kFiberContactSnapEta,
//    interactions/fiber_contact_common.hpp, общая с Fiber-Fiber-Contact):
//      if (t > 1.0 - eta) t = std::max(t, 1.0);
//      if (t < eta)        t = std::min(t, 0.0);
//    Значения СТРОГО за пределами [0,1] сохраняются -- прижимаются только
//    значения у самих границ (защита от рассогласованного округления у двух
//    соседних сегментов, которые иначе оба могли бы принять один и тот же
//    узел; для прямой нити это на практике не воспроизводится, но прижатие
//    делает правило независимым от округления). После прижатия -- правило по
//    таблице (primitive (a,b,s), b=next(a)):
//      0 < t < 1                                 -> контакт, λ=t;
//      t >= 1, b_is_last (next(b)==-1)            -> контакт, λ=1;
//      t >= 1, !b_is_last                         -> пропуск (узел b -- у
//                                                     следующего сегмента);
//      t <= 0, первый сегмент (a_prev==-1)        -> контакт, λ=0;
//      t <= 0, НЕ первый (a_prev!=-1)              -> вычислить t_prev той же
//                                                     функцией для сегмента
//                                                     (a_prev,a) с тем же
//                                                     прижатием; t_prev>=1 ->
//                                                     контакт, λ=0; иначе ->
//                                                     пропуск.
//    Почему это даёт ровно один контакт: прямая нить, частица над узлом --
//    сегмент (prev(a),a) получает t>=1 и отказывается (не последний), а
//    сегмент (a,b) получает t<=0, проверяет соседа, видит t_prev>=1 и
//    принимает; ВНЕШНИЙ сектор излома -- то же самое (оба соседних сегмента
//    участвуют в том же t<=0/t>=1 разборе, контакт с λ=0 у вершины вместо
//    потери); частица проецируется внутрь ПРЕДЫДУЩЕГО сегмента (0<t_prev<1)
//    -- его принимает предыдущий сегмент, текущий отказывается по проверке
//    t_prev; вогнутый излом, частица проецируется внутрь ОБОИХ сегментов --
//    два контакта в разных точках (это физически две точки касания, так и
//    должно быть). Тест test_fiber_solid.cpp, "узел не получает двойную
//    силу" проверяет частный случай (прямая нить) численно: Solid-частица в
//    узле между двумя сегментами, сумма вкладов от ОБОИХ primitives равна
//    вкладу ОДНОГО контакта.
//    FiberSolidContactPrimitive хранит a_prev=prev(a) (long long, -1 если нет)
//    вместо старого bool a_is_first -- нужен САМ индекс бусины prev(a) для
//    вызова fiber_solid_segment_param на соседнем сегменте. b_is_last не
//    изменился. Точечный primitive (изолированная бусина) не меняется по
//    смыслу.
//
// 4. K_fs берёт FiberMaterial ИМЕННО бусины a (PLAN.md, этап Ф6, пункт 4,
//    дословно: "использует elastic_by_itype(itype_s) и FiberMaterial
//    материала нити a") — не b. Для однородной нити (все бусины одного
//    материала, как во всех текущих сценариях) выбор a/b не имеет значения.
//
// 5. dt_limit: m_red = m_a*m_s/(m_a+m_s) -- по аналогии с Fiber-Fiber-Contact
//    (interactions/fiber_contact.hpp: там, после перехода на схему «пара
//    сегментов -- одна сила» (fiber/fiber_contact_segment_pairs.md), m_red
//    строится по ЭФФЕКТИВНЫМ массам точек контакта на каждой стороне пары
//    сегментов, а не по одной "якорной" бусине) и с Solid-Solid-Contact.
//    Здесь, в Fiber-Solid, схема проще: бусина b в m_red не участвует
//    (получает лишь долю λ силы, ровно как соседняя бусина сегмента в
//    Fiber-Fiber-Contact).
//
// 6. Предупреждение "Fiber-Solid-Contact не формируется для itype=..(газ/
//    вода)" -- напечатано ОДИН РАЗ ЗА ПРОЦЕСС ВЫПОЛНЕНИЯ на itype (static
//    std::set в build_fiber_solid_contact_search, см. ниже), при первой
//    встрече такой частицы среди реально используемых Solid при построении
//    поиска (вызывается только когда нити есть, fiber_file непустой -- иначе
//    fb.size()==0 и время_integration вообще не входит в этот блок,
//    src/time_integration.cpp). Минимально достаточный вариант, явно
//    допущенный постановкой (fiber/PLAN.md, этап Ф6, пункт 1): глубже
//    встраивать проверку в parse_sph_config избыточно, т.к. на этапе разбора
//    конфигурации ещё не построен список "реально используемых" Solid-частиц
//    в удобной для этой проверки форме.
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <vector>
#include <iostream>

#include "fiber/fiber_system.hpp"
#include "interactions/fiber_contact_common.hpp" // kFiberContactSnapEta, fiber_max_current_bond_length
#include "interactions/solid_contact.hpp"        // solid_contact_radius
#include "interactions/solid_sph.hpp"            // SolidDelta
#include "materials.hpp"                         // elastic_by_itype
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/spatial_index.hpp"
#include "rhs/types.hpp"

namespace interactions
{

// Метки interaction_type для ContributionMeta -- разные Delta (FiberDelta,
// SolidDelta) суммируются в НЕЗАВИСИМЫХ reduce_contributions, поэтому
// коллизия численного значения между "доменами" разных Delta не критична, но
// внутри своего домена значение обязано быть уникальным:
//   FiberDelta:  1 = Fiber-Bond, 2 = Fiber-Fiber-Contact, 3 = ЗДЕСЬ (сторона Fiber);
//   SolidDelta: 10 = Solid-SPH, 11 = Solid-Solid-Contact, 12 = ЗДЕСЬ (сторона Solid).
constexpr std::uint32_t kFiberSolidContactFiberType = 3;
constexpr std::uint32_t kFiberSolidContactSolidType = 12;

constexpr double kFiberSolidContactMinD = 1e-12;

// ----------------------------------------------------------------------
// FiberSolidContactPrimitive = (a, b, s) -- fiber/SPH.md, «Fiber-Solid-
// Contact»: a, b -- соседние бусины сегмента нити (b может отсутствовать --
// вырожденный точечный контакт, has_b==false), s -- РЕАЛЬНАЯ Solid-частица.
// a_prev/b_is_last -- топологические данные, нужные kernel'у для правила
// владения узлом (см. пункт 3 в шапке файла): b_is_last решает "клэмпить t в
// 1 или передать ответственность следующему сегменту"; a_prev=prev(a) (не
// просто флаг "первый ли сегмент", как раньше a_is_first) -- нужен САМ
// индекс бусины prev(a), т.к. при t<=0 у не первого сегмента правило
// пересчитывает параметр для соседнего сегмента (prev(a),a) функцией
// fiber_solid_segment_param (fiber/fiber_solid_contact_fix.md, §1 п.3-4).
// Сохранены в primitive при построении (тот же приём, что
// FiberContactPrimitive::a_prev/a_has_next_segment (и симметрично для
// стороны B) в interactions/fiber_contact.hpp — kernel не перечитывает
// topology).
// ----------------------------------------------------------------------
struct FiberSolidContactPrimitive
{
    rhs::index_t a = 0;
    rhs::index_t s = 0;
    bool has_b = false;
    rhs::index_t b = 0;
    long long a_prev = -1;  // prev[a]; -1, если нет (значимо всегда, даже при !has_b)
    bool b_is_last = false; // has_b && next[b] == -1
};

// ----------------------------------------------------------------------
// fiber_solid_contact_search_radius -- r_f_max + r_s_max^contact + L_max
// (fiber/PLAN.md §3.2, таблица "Радиусы поиска"; fiber/
// fiber_solid_contact_fix.md, §2).
//
// L_max -- ТЕКУЩАЯ (не номинальная l0!) максимальная длина связи
// (interactions/fiber_contact_common.hpp::fiber_max_current_bond_length,
// общая с Fiber-Fiber-Contact): связи Fiber-Bond работают только на
// растяжение и рвутся лишь в Postprocessing ПОСЛЕ шага (fiber/SPH.md),
// поэтому внутри шага сегмент может быть длиннее номинального l0, и l0 как
// оценка не годится.
//
// Эта же оценка обосновывает решение (в шапке файла, пункт 2) о том, что
// primitive НЕ строится для пары (f,s), если next(f)==-1 и prev(f)!=-1: если
// контакт существует у Solid-частицы s рядом с концом f сегмента (prev(f),f),
// т.е. |x_s - x_f| <= r_f+r_s, то по неравенству треугольника
//   |x_s - x_{prev(f)}| <= |x_s-x_f| + |x_f-x_{prev(f)}|
//                        <= (r_f+r_s) + |x_f-x_{prev(f)}|
//                        <= r_f_max + r_s_max^contact + L_max = radius,
// то есть пара (prev(f), s) гарантированно попадает в результат поиска
// (|x_f-x_{prev(f)}| -- ТЕКУЩАЯ длина связи (prev(f),f), значит она <= L_max
// = максимуму ТЕКУЩЕЙ длины ПО ВСЕМ связям расчёта).
// ----------------------------------------------------------------------
inline double fiber_solid_contact_search_radius(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                                  const fiber::FiberTopology& topology,
                                                  const std::vector<fiber::FiberMaterial>& materials,
                                                  double r_s_max_contact)
{
    double r_f_max = 0.0;
    for (const auto& m : materials) r_f_max = std::max(r_f_max, 0.5 * m.d);
    const double l_max = fiber_max_current_bond_length(exec, state, topology);
    return r_f_max + r_s_max_contact + l_max;
}

// Предупреждение "Fiber-Solid-Contact не формируется для газа/воды" -- один
// раз за процесс на itype (решение 6 в шапке файла).
inline void warn_fiber_solid_contact_undefined_itype_once(int itype)
{
    static std::set<int> warned;
    if (warned.insert(itype).second)
    {
        std::wcout << "ПРЕДУПРЕЖДЕНИЕ: Fiber-Solid-Contact не формируется для "
                      "Solid-частиц itype=" << itype << " -- модуль Юнга/Пуассона "
                      "не определены (нет модуля объёмного сжатия и/или сдвига, "
                      "напр. газ/вода itype 1-2; fiber/SPH.md, \"Fiber-Solid-"
                      "Contact\", \"Постановка проекта: коэффициент K_fs\")."
                   << std::endl;
    }
}

struct FiberSolidContactSearchResult
{
    std::vector<FiberSolidContactPrimitive> primitives;
};

// ----------------------------------------------------------------------
// build_fiber_solid_contact_search -- собирает подмножество contact-eligible
// РЕАЛЬНЫХ Solid-частиц ([0,ntotal), elastic_by_itype определён), строит по
// нему SpatialIndex<dim> и запрашивает "чужой тип" (Fiber -> Solid,
// find_pairs_cross_parallel), затем классифицирует каждую найденную пару в
// primitive по топологии (решение 2 в шапке файла).
// ----------------------------------------------------------------------
inline void build_fiber_solid_contact_search(const rhs::CpuExecutor& exec, const Materials& mats,
                                              const ParticleSystem& sys, std::size_t ntotal,
                                              const fiber::FiberState& state, const fiber::FiberTopology& topology,
                                              const std::vector<fiber::FiberMaterial>& materials,
                                              FiberSolidContactSearchResult& out)
{
    out.primitives.clear();
    const std::size_t nfiber = state.size();
    if (nfiber == 0 || ntotal == 0) return;

    if constexpr (dim != 3)
    {
        // Нити физически существуют только при dim==3 (см. fiber_contact.hpp,
        // тот же паттерн заглушки для компиляции во всех SPH_DIM-вариантах).
        return;
    }
    else
    {
        std::vector<rhs::index_t> eligible;
        eligible.reserve(ntotal);
        double r_s_max = 0.0;
        for (std::size_t i = 0; i < ntotal; ++i)
        {
            const auto es = elastic_by_itype(mats, sys.itype[i]);
            if (!es.has_value())
            {
                warn_fiber_solid_contact_undefined_itype_once(sys.itype[i]);
                continue;
            }
            eligible.push_back(static_cast<rhs::index_t>(i));
            r_s_max = std::max(r_s_max, solid_contact_radius(sys, i));
        }
        if (eligible.empty()) return;

        const double radius = fiber_solid_contact_search_radius(exec, state, topology, materials, r_s_max);

        std::vector<double> sx[dim];
        for (std::size_t d = 0; d < dim; ++d) sx[d].resize(eligible.size());
        for (std::size_t k = 0; k < eligible.size(); ++k)
            for (std::size_t d = 0; d < dim; ++d) sx[d][k] = sys.x[d][eligible[k]];

        rhs::SpatialIndex<dim> index;
        std::array<const double*, dim> scoords{};
        for (std::size_t d = 0; d < dim; ++d) scoords[d] = sx[d].data();
        index.build(scoords, eligible.size(), radius > 0.0 ? radius : 1.0);

        std::array<const double*, dim> fcoords{};
        for (std::size_t d = 0; d < dim; ++d) fcoords[d] = state.x[d].data();

        const auto pairs = index.find_pairs_cross_parallel(exec, fcoords, nfiber, radius);

        out.primitives.reserve(pairs.size());
        for (const auto& pr : pairs)
        {
            const rhs::index_t f = pr.first;        // бусина Fiber ("other")
            const rhs::index_t local_s = pr.second;  // локальный индекс в eligible
            const rhs::index_t s = eligible[local_s];

            FiberSolidContactPrimitive prim;
            prim.a = f;
            prim.s = s;
            prim.a_prev = topology.prev[f]; // -1, если бусина первая/изолированная

            const long long nxt = topology.next[f];
            if (nxt == -1)
            {
                if (prim.a_prev == -1)
                {
                    // Изолированная бусина без соседей -- точечный контакт
                    // (решение 2 в шапке файла).
                    prim.has_b = false;
                    out.primitives.push_back(prim);
                }
                // else: сегмент (prev(f),f) покрывает эту зону -- пропуск
                // (доказано в комментарии у fiber_solid_contact_search_radius).
                continue;
            }

            prim.has_b = true;
            prim.b = static_cast<rhs::index_t>(nxt);
            prim.b_is_last = (topology.next[prim.b] == -1);
            out.primitives.push_back(prim);
        }
    }
}

// ----------------------------------------------------------------------
// FiberSolidSegmentParam -- результат fiber_solid_segment_param ниже:
// valid=false <=> сегмент (xa,xb) вырожден (e2<=0, совпадающие бусины) --
// параметр t не определён (защита от деления на 0, то же поведение, что
// раньше -- прямой отказ от контакта на этом сегменте).
// ----------------------------------------------------------------------
struct FiberSolidSegmentParam
{
    bool valid = false;
    double t = 0.0;
};

// ----------------------------------------------------------------------
// fiber_solid_segment_param -- ЕДИНСТВЕННАЯ функция вычисления параметра t
// проекции Solid-частицы xs на сегмент нити (xa,xb) (fiber/
// fiber_solid_contact_fix.md, §1 п.1): используется и для "своего" сегмента
// primitive (a,b), и для соседнего (prev(a),a) при проверке правила
// владения узлом в compute_fiber_solid_contact_geometry -- поэтому
// результат должен быть побитово воспроизводим для одинаковых входов (то же
// выражение, тот же порядок операций, что раньше был инлайнен в
// compute_fiber_solid_contact_geometry -- для прямой нити БЕЗ контактов у
// узлов результат остаётся побитово прежним):
//   t = dot(xs-xa, e)/dot(e,e),  e = xb-xa.
// Прижатие к границам (η=kFiberContactSnapEta, interactions/
// fiber_contact_common.hpp, общая с Fiber-Fiber-Contact; fiber/
// fiber_solid_contact_fix.md, §1 п.2): значения СТРОГО за пределами [0,1]
// сохраняются -- прижимаются только значения у самих границ (см. шапку
// файла, решение 3, для полного разбора и обоснования):
//   if (t > 1.0 - eta) t = std::max(t, 1.0);
//   if (t < eta)       t = std::min(t, 0.0);
// ----------------------------------------------------------------------
inline FiberSolidSegmentParam fiber_solid_segment_param(const double xa[3], const double xb[3], const double xs[3])
{
    FiberSolidSegmentParam out;
    const double e[3] = {xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]};
    const double e2 = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    if (!(e2 > 0.0)) return out; // вырожденный сегмент (совпадающие бусины) -- valid остаётся false

    const double d[3] = {xs[0] - xa[0], xs[1] - xa[1], xs[2] - xa[2]};
    double t = (d[0] * e[0] + d[1] * e[1] + d[2] * e[2]) / e2;

    if (t > 1.0 - kFiberContactSnapEta) t = std::max(t, 1.0);
    if (t < kFiberContactSnapEta) t = std::min(t, 0.0);

    out.valid = true;
    out.t = t;
    return out;
}

// ----------------------------------------------------------------------
// FiberSolidContactGeometry -- итог геометрического решения ОДНОГО primitive
// (fiber/SPH.md, «Fiber-Solid-Contact»): λ, ближайшая точка q, D_fs, n_fs,
// δ_fs, сила Герца F^c_fs. Общая точка входа для contributions и dt_limit
// (аналог FiberContactGeometry в interactions/fiber_contact.hpp).
// ----------------------------------------------------------------------
struct FiberSolidContactGeometry
{
    bool contact = false;
    double n[dim] = {};
    double D = 0.0;
    double delta = 0.0;
    double F = 0.0;      // >= 0
    double lambda = 0.0; // значим только при has_b
};

inline FiberSolidContactGeometry compute_fiber_solid_contact_geometry(
    const Materials& mats, const ParticleSystem& sys, const fiber::FiberState& state,
    const std::vector<fiber::FiberMaterial>& materials, const FiberSolidContactPrimitive& prim)
{
    FiberSolidContactGeometry geom{};
    if constexpr (dim != 3)
    {
        return geom; // недостижимо при реальном использовании, см. build_fiber_solid_contact_search
    }
    else
    {
        using Vec3 = std::array<double, 3>;
        auto sub = [](const Vec3& x, const Vec3& y) -> Vec3 { return {x[0] - y[0], x[1] - y[1], x[2] - y[2]}; };
        auto dot = [](const Vec3& x, const Vec3& y) -> double { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; };
        auto norm = [&](const Vec3& x) -> double { return std::sqrt(dot(x, x)); };
        auto point_f = [&](rhs::index_t k) -> Vec3 { return {state.x[0][k], state.x[1][k], state.x[2][k]}; };
        auto point_s = [&](rhs::index_t k) -> Vec3 { return {sys.x[0][k], sys.x[1][k], sys.x[2][k]}; };

        const Vec3 xs = point_s(prim.s);
        const fiber::FiberMaterial& mat_a = materials[static_cast<std::size_t>(state.material_id[prim.a])];
        const double r_f = 0.5 * mat_a.d;

        Vec3 q{};
        double lambda = 0.0;
        if (!prim.has_b)
        {
            q = point_f(prim.a);
        }
        else
        {
            const Vec3 xa = point_f(prim.a);
            const Vec3 xb = point_f(prim.b);
            const Vec3 e = sub(xb, xa);

            // Правило владения узлом (fiber/fiber_solid_contact_fix.md, §1;
            // шапка файла, решение 3) -- параметр t (после прижатия к
            // границе) вычисляется общей функцией fiber_solid_segment_param.
            const FiberSolidSegmentParam pt = fiber_solid_segment_param(xa.data(), xb.data(), xs.data());
            if (!pt.valid) return geom; // вырожденный сегмент (совпадающие бусины) -- защита от деления на 0
            const double t = pt.t;

            if (t >= 1.0)
            {
                if (!prim.b_is_last) return geom; // узел b -- у следующего сегмента
                lambda = 1.0;
            }
            else if (t <= 0.0)
            {
                if (prim.a_prev == -1)
                {
                    lambda = 0.0; // первый сегмент нити -- клэмпить некому
                }
                else
                {
                    // Узел a -- общая точка с соседним сегментом (a_prev,a);
                    // владелец решается пересчётом t_prev ТОЙ ЖЕ функцией с
                    // тем же прижатием -- принимаем, только если у соседа
                    // минимум оказался на его дальнем конце (t_prev>=1),
                    // иначе минимум внутри соседнего сегмента, и владеет им
                    // он (см. таблицу в шапке файла, решение 3).
                    const Vec3 xap = point_f(static_cast<rhs::index_t>(prim.a_prev));
                    const FiberSolidSegmentParam pt_prev = fiber_solid_segment_param(xap.data(), xa.data(), xs.data());
                    if (!(pt_prev.valid && pt_prev.t >= 1.0)) return geom;
                    lambda = 0.0;
                }
            }
            else
            {
                lambda = t;
            }
            q = {xa[0] + lambda * e[0], xa[1] + lambda * e[1], xa[2] + lambda * e[2]};
        }

        const Vec3 diff = sub(xs, q);
        const double D = norm(diff);
        if (D < kFiberSolidContactMinD) return geom; // совпадающие точки -- пропуск, защита от NaN
        const Vec3 n = {diff[0] / D, diff[1] / D, diff[2] / D};

        const double r_s = solid_contact_radius(sys, prim.s);
        const double delta = r_f + r_s - D;
        if (!(delta > 0.0)) return geom; // нет проникновения

        const auto es = elastic_by_itype(mats, sys.itype[prim.s]);
        if (!es.has_value()) return geom; // защита (search уже отфильтровал такие s)

        const double inv_Estar = (mat_a.E > 0.0 ? (1.0 - mat_a.nu * mat_a.nu) / mat_a.E : 0.0) +
                                  (es->E > 0.0 ? (1.0 - es->nu * es->nu) / es->E : 0.0);
        if (!(inv_Estar > 0.0)) return geom;
        const double E_star = 1.0 / inv_Estar;
        const double R_star = (r_f > 0.0 && r_s > 0.0) ? (r_f * r_s) / (r_f + r_s) : 0.0;
        const double K_fs = (4.0 / 3.0) * E_star * std::sqrt(R_star);
        const double F = K_fs * std::pow(delta, 1.5);

        geom.contact = true;
        geom.n[0] = n[0]; geom.n[1] = n[1]; geom.n[2] = n[2];
        geom.D = D;
        geom.delta = delta;
        geom.F = F;
        geom.lambda = lambda;
        return geom;
    }
}

// ----------------------------------------------------------------------
// compute_fiber_solid_contact_contributions -- kernel: эмитирует вклад в
// Solid s (SolidDelta, только dv[dim] -- нет неразрывности/энергии/
// деформации от контакта, как и Solid-Solid-Contact, interactions/
// solid_contact.hpp) и в Fiber a (+b, если has_b) (FiberDelta). Сохранение
// импульса на primitive: F_a+F_b+F_s=0 (F_s=F*n, F_a=-(1-λ)F*n,
// F_b=-λF*n) -- см. test_fiber_solid.cpp.
// ----------------------------------------------------------------------
inline void compute_fiber_solid_contact_contributions(
    const rhs::CpuExecutor& exec, const Materials& mats, const ParticleSystem& sys, const fiber::FiberState& state,
    const std::vector<fiber::FiberMaterial>& materials, const std::vector<FiberSolidContactPrimitive>& primitives,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& fiber_buffers,
    std::vector<rhs::ContributionBuffer<SolidDelta>>& solid_buffers)
{
    const std::size_t nworkers = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    if (fiber_buffers.size() < nworkers) fiber_buffers.resize(nworkers);
    if (solid_buffers.size() < nworkers) solid_buffers.resize(nworkers);
    for (auto& buf : fiber_buffers) buf.clear();
    for (auto& buf : solid_buffers) buf.clear();

    const std::size_t n = primitives.size();
    if (n == 0) return;

    exec.parallel_for(n, [&](std::size_t k) {
        const FiberSolidContactPrimitive& prim = primitives[k];
        const FiberSolidContactGeometry geom = compute_fiber_solid_contact_geometry(mats, sys, state, materials, prim);
        if (!geom.contact) return;

        int w = rhs::CpuExecutor::worker_id();
        if (w < 0 || static_cast<std::size_t>(w) >= fiber_buffers.size()) w = 0;
        rhs::ContributionWriter<fiber::FiberDelta> fw(fiber_buffers[static_cast<std::size_t>(w)],
                                                        kFiberSolidContactFiberType);
        rhs::ContributionWriter<SolidDelta> sw(solid_buffers[static_cast<std::size_t>(w)],
                                                kFiberSolidContactSolidType);

        const double ms = sys.mass[prim.s];
        SolidDelta ds{};
        for (std::size_t d = 0; d < dim; ++d) ds.dv[d] = geom.F * geom.n[d] / ms;
        sw.emit(prim.s, ds, static_cast<std::uint64_t>(k), 0);

        if (!prim.has_b)
        {
            const double ma = state.mass[prim.a];
            fiber::FiberDelta da{};
            for (std::size_t d = 0; d < dim; ++d) da.dv[d] = -geom.F * geom.n[d] / ma;
            fw.emit(prim.a, da, static_cast<std::uint64_t>(k), 1);
        }
        else
        {
            const double ma = state.mass[prim.a];
            const double mb = state.mass[prim.b];
            fiber::FiberDelta da{}, db{};
            for (std::size_t d = 0; d < dim; ++d)
            {
                da.dv[d] = -(1.0 - geom.lambda) * geom.F * geom.n[d] / ma;
                db.dv[d] = -geom.lambda * geom.F * geom.n[d] / mb;
            }
            fw.emit(prim.a, da, static_cast<std::uint64_t>(k), 1);
            fw.emit(prim.b, db, static_cast<std::uint64_t>(k), 2);
        }
    });
}

// compute_fiber_solid_contact -- удобная обёртка: search + contributions +
// reduce обеих сторон (Fiber и Solid).
struct FiberSolidContactResult
{
    std::vector<fiber::FiberDelta> fiber_delta;
    std::vector<SolidDelta> solid_delta;
    std::vector<FiberSolidContactPrimitive> primitives; // для dt_limit/диагностики без повторного поиска
};

inline FiberSolidContactResult compute_fiber_solid_contact(
    const rhs::CpuExecutor& exec, const Materials& mats, const ParticleSystem& sys, std::size_t ntotal,
    const fiber::FiberState& state, const fiber::FiberTopology& topology,
    const std::vector<fiber::FiberMaterial>& materials,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& fiber_buffers,
    std::vector<rhs::ContributionBuffer<SolidDelta>>& solid_buffers)
{
    FiberSolidContactSearchResult search;
    build_fiber_solid_contact_search(exec, mats, sys, ntotal, state, topology, materials, search);
    compute_fiber_solid_contact_contributions(exec, mats, sys, state, materials, search.primitives, fiber_buffers,
                                               solid_buffers);
    FiberSolidContactResult out;
    out.fiber_delta = rhs::reduce_contributions(exec, fiber_buffers, state.size());
    out.solid_delta = rhs::reduce_contributions(exec, solid_buffers, ntotal);
    out.primitives = search.primitives;
    return out;
}

// ----------------------------------------------------------------------
// compute_fiber_solid_contact_dt_limit -- fiber/PLAN.md §3.7:
//     dt_lim = 0.8/sqrt(k_eff/m_red),  k_eff = 1.5*F/delta,
//     m_red = m_a*m_s/(m_a+m_s) (решение 5 в шапке файла). Считается ТОЛЬКО
// при фактическом контакте (geom.contact, delta>0) -- как Fiber-Fiber-
// Contact/Solid-Solid-Contact. +inf, если контактов нет.
// ----------------------------------------------------------------------
inline double compute_fiber_solid_contact_dt_limit(const rhs::CpuExecutor& exec, const Materials& mats,
                                                     const ParticleSystem& sys, const fiber::FiberState& state,
                                                     const std::vector<fiber::FiberMaterial>& materials,
                                                     const std::vector<FiberSolidContactPrimitive>& primitives)
{
    const double kInf = std::numeric_limits<double>::max();
    const std::size_t n = primitives.size();
    if (n == 0) return kInf;

    return exec.parallel_reduce(
        n, kInf,
        [&](std::size_t k) -> double {
            const FiberSolidContactPrimitive& prim = primitives[k];
            const FiberSolidContactGeometry geom = compute_fiber_solid_contact_geometry(mats, sys, state, materials, prim);
            if (!geom.contact || !(geom.delta > 0.0)) return kInf;

            const double ma = state.mass[prim.a];
            const double ms = sys.mass[prim.s];
            const double m_red = (ma * ms) / (ma + ms);
            const double k_eff = 1.5 * geom.F / geom.delta;
            if (!(k_eff > 0.0)) return kInf;
            return 0.8 / std::sqrt(k_eff / m_red);
        },
        [](double x, double y) { return std::min(x, y); });
}

} // namespace interactions
