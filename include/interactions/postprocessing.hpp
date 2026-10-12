#pragma once
// interactions/postprocessing.hpp
// Postprocessing (fiber/PLAN.md, этап Ф7; fiber/SPH.md, разделы
// «Postprocessing» и «Postprocessing и изменение topology»): дискретные
// пороговые правила, выполняемые ПОСЛЕ шага ODESolver -- в две фазы, где это
// применимо (поиск событий readonly параллельно -> применение
// последовательно):
//
//   - Мизес (itype 3-9) / JH-2 (itype 10-11)  -- apply_material_strength;
//   - sle = 1 (пересчёт h)                     -- update_smoothing_length;
//   - разрыв связей Fiber-Bond по l > l_break  -- find_fiber_break_events /
//                                                  apply_fiber_break_events /
//                                                  run_fiber_topology_postprocessing;
//   - обновление внутренних переменных связи   -- apply_fiber_bond_constitutive
//     (пластика, Максвелл, диссипация; ПЕРЕД поиском разрыва).
//
// Единая функция для ОБЕИХ перегрузок time_integration
// (src/time_integration.cpp) -- было два независимых инлайновых блока кода
// (Мизес/JH-2) и один блок внутри single_step.cpp (h_upgrade), теперь общий
// код.
//
// ============================================================================
// МЕСТО ВЫЗОВА update_smoothing_length В ЦИКЛЕ -- ВАЖНО, ЧИТАТЬ ПЕРЕД
// ПЕРЕСТАНОВКОЙ:
//
// В отличие от apply_material_strength/фазы разрыва связей (которые
// действительно можно и нужно вызывать ПОСЛЕ corrector'а -- они используют
// финальное, скорректированное состояние шага n), update_smoothing_length
// (sle=1: h = fac*(mass/rho)^(1/dim)) ОБЯЗАН вызываться на прежнем месте --
// сразу после single_step(), ДО compute_cfl_timestep()/corrector'а, как это
// было раньше ВНУТРИ single_step.cpp. Причина -- численная, не
// косметическая:
//
//   sys.rho в момент вызова single_step(n) -- это ПРЕДСКАЗАННАЯ (predictor,
//   половинный шаг от rho^(n-1) по производной f^(n-1)) плотность, а НЕ
//   финальная (corrector) плотность шага n. Раньше h_upgrade(sle=1) внутри
//   single_step(n) считал h именно от этой предсказанной rho -- ДО того, как
//   compute_cfl_timestep() использует получившийся h для CFL-шага ЭТОГО ЖЕ
//   шага n (и раньше, чем corrector вообще посчитает финальную rho^(n)).
//
//   Если бы update_smoothing_length вызывался буквально "после шага
//   ODESolver" (после corrector'а, от финальной rho^(n)), это было бы ДРУГОЕ
//   число (rho^(n) != rho*_n в общем случае, они отличаются на слагаемое
//   порядка (dt/2)*(f^(n)-f^(n-1))) -- и, что важнее, это НЕ УСПЕВАЛО бы к
//   compute_cfl_timestep() ЭТОГО ЖЕ шага (та вызывается ДО corrector'а),
//   так что CFL-шаг n использовал бы УСТАРЕВШИЙ h (с предыдущего шага) --
//   двойное расхождение с прежним поведением. Допуски tests/test_regression.cpp
//   откалиброваны на уровне ~1e-11..1e-10 (шум порядка потоков/округления
//   представления), а не на физическую грубую точность -- ЛЮБОЕ из этих
//   расхождений на порядки больше допуска и валит регрессию (проверено;
//   см. отчёт по этапу Ф7 в истории коммитов).
//
// Поэтому update_smoothing_length вызывается в time_integration.cpp СРАЗУ
// после single_step(), ДО compute_cfl_timestep() -- то есть на том самом
// месте в графе вычислений, где раньше физически проявлялся эффект вызова
// h_upgrade изнутри single_step.cpp. Формально это "внутри шага ODESolver",
// а не строго "после" него в смысле fiber/SPH.md -- сознательное отступление
// от буквальной формулировки ради нулевого изменения численного результата
// Solid-физики (жёсткое требование этапа Ф7 fiber/PLAN.md: «BC2/BC6 не
// изменились»). Мизес/JH-2 и разрыв связей такому ограничению не подчиняются
// (они и раньше считались после corrector'а) и вызываются в буквальном месте
// -- после corrector'а, вместе.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "h_upgrade.hpp"
#include "interactions/fiber_bond.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "ppochnost.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"

namespace interactions
{

// ----------------------------------------------------------------------
// apply_material_strength -- Мизес (itype 3-9) и JH-2 (itype 10-11),
// fiber/PLAN.md §3.8. Код -- буквальная копия цикла, ранее продублированного
// в обеих перегрузках time_integration (src/time_integration.cpp), вызывается
// ПОСЛЕ corrector'а -- читает ФИНАЛЬНЫЕ (скорректированные) sys.sxx.. и т.д.
// Шим совместимости (interactions/solid_sph.hpp, комментарий у single_step):
// prochnost_cem_JH2 по-прежнему читает 9 величин деформации/поворота
// (exx..eyz, rxy..ryz) из sys, а не явными параметрами -- решение
// координатора этапа Ф7: сигнатуру prochnost_cem_JH2 (src/ppochnost.cpp) не
// меняем, риск регрессии от смены источника данных не оправдан архитектурной
// чистотой на этом этапе (задача допускает оба варианта, см. fiber/PLAN.md
// §Ф7, пункт 5 обоснования координатора).
// ----------------------------------------------------------------------
inline void apply_material_strength(const SimulationConfig& cfg, ParticleSystem& sys, std::size_t ntotal, double dt,
                                     const std::vector<double>& drho)
{
    for (std::size_t i = 0; i < ntotal; ++i)
    {
        if (sys.itype[i] >= 3 && sys.itype[i] <= 9)
        {
            prochnost_mises(cfg.mats, static_cast<size_t>(sys.itype[i]), sys.sxx[i], sys.syy[i], sys.szz[i],
                             sys.sxy[i], sys.syz[i], sys.sxz[i]);
        }
        else if (sys.itype[i] == 10 || sys.itype[i] == 11)
        {
            prochnost_cem_JH2(cfg.mats, sys, i, dt, drho[i]);
        }
    }
}

// ----------------------------------------------------------------------
// update_smoothing_length -- sle=1: h = fac*(mass/rho)^(1/dim)
// (src/h_upgrade.cpp), fiber/PLAN.md §3.8 "sle=1 -- Postprocessing". sle=0 --
// no-op (h_upgrade не вызывается вовсе, как и раньше). sle=2 -- явный отказ:
// путь RHS (fiber/PLAN.md, этап Ф4) не строит сумму по парам соседей (div v)
// в формате старого поиска, которую требует эта ветка h_upgrade; ни один
// текущий конфиг (BC2/BC6 -> sle=1) sle=2 не использует. ПРО МЕСТО ВЫЗОВА --
// см. большой комментарий в шапке файла.
// ----------------------------------------------------------------------
inline void update_smoothing_length(const SimulationConfig& cfg, ParticleSystem& sys, double dt, std::size_t ntotal)
{
    if (cfg.sle == 2)
    {
        throw std::runtime_error(
            "update_smoothing_length: sle=2 требует сумму по парам соседей (div v), которую "
            "путь RHS (fiber/PLAN.md, этап Ф4) не строит в формате старого "
            "поиска для h_upgrade; ни один текущий конфиг (BC2/BC6 -> sle=1) "
            "sle=2 не использует.");
    }
    if (cfg.sle == 1)
    {
        static const std::vector<size_t> kEmptyPairIdx;
        static const std::vector<std::vector<double>> kEmptyDwdx(dim);
        h_upgrade(cfg, sys, dt, ntotal, 0, kEmptyPairIdx, kEmptyPairIdx, kEmptyDwdx);
    }
    // sle == 0 (и любое иное непредусмотренное значение) -- h_upgrade не
    // вызывается, как и раньше (single_step.cpp вызывал h_upgrade ТОЛЬКО при
    // sle == 1).
}

// Текущая длина связи |x_j - x_i| (общая для поиска разрыва и прохода
// apply_fiber_bond_constitutive; порядок операций прежний).
inline double fiber_bond_length(const fiber::FiberState& state, const fiber::FiberBond& b)
{
    double l2 = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        const double dx = state.x[d][b.j] - state.x[d][b.i];
        l2 += dx * dx;
    }
    return std::sqrt(l2);
}

// ----------------------------------------------------------------------
// apply_fiber_bond_constitutive -- обновление внутренних переменных связей
// (eps_p, q, eps_prev, w_p, w_v) по определяющему закону
// fiber::bond_constitutive_update (fiber/fiber-plasticity-design.md, раздел 4).
// Вызывается в time_integration ПОСЛЕ advance (x^{n+1} финальны) и ПЕРЕД
// поиском разрыва, с dt текущего шага. Параллельно по связям через exec: связь k
// читает только state.x своих двух бусин и пишет только topology.bonds[k].st --
// гонок нет, результат не зависит от числа потоков.
//
// Упругий материал (FiberMaterial::inelastic() == false) НЕ трогается: его st
// остаётся нулевым, проход ничего не меняет. Если неупругого материала нет
// вовсе, проход выходит сразу после одного обхода списка материалов.
// ----------------------------------------------------------------------
inline void apply_fiber_bond_constitutive(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                          fiber::FiberTopology& topology,
                                          const std::vector<fiber::FiberMaterial>& materials, double dt)
{
    if (!(dt > 0.0)) return;
    const std::size_t nbonds = topology.bonds.size();
    if (nbonds == 0) return;
    const bool any_inelastic = std::any_of(materials.begin(), materials.end(),
                                           [](const fiber::FiberMaterial& m) { return m.inelastic(); });
    if (!any_inelastic) return;

    exec.parallel_for(nbonds, [&](std::size_t k) {
        fiber::FiberBond& b = topology.bonds[k];
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(state.material_id[b.i])];
        if (!mat.inelastic()) return;
        const double eps = fiber_bond_length(state, b) / b.l0 - 1.0;
        const double A = fiber_cross_section_area(mat);
        const fiber::ViscoParams vp = fiber_bond_visco(mat, A, b.l0, state.mass[b.i], state.mass[b.j]);
        const fiber::BondLaw law(mat.E, mat.plastic, mat.rate, vp, mat.failure);
        b.st = fiber::bond_constitutive_update(law, b.st, eps, dt);
    });
}

// ----------------------------------------------------------------------
// init_fiber_bond_state -- начальное eps_prev неупругих связей. Связь, созданная
// при l != l0 (предварительное натяжение, сцена со сдвигом бусин), иначе дала бы
// на первом шаге фиктивное приращение Δε = ε - 0 и огромные w_v / вязкую силу.
// Выставляет eps_prev = текущей ε = l/l0 - 1 для связей неупругих материалов
// (остальные поля состояния не трогает: eps_p = q = 0 -- мгновенная упругая
// нагрузка). Упругие связи НЕ меняются (их st остаётся нулевым). Вызывается
// один раз перед первым шагом (run_time_integration в начале; драйверы тестов --
// сами); повторный вызов безвреден, если x не менялись после последнего шага.
// ----------------------------------------------------------------------
inline void init_fiber_bond_state(const fiber::FiberState& state, fiber::FiberTopology& topology,
                                  const std::vector<fiber::FiberMaterial>& materials)
{
    for (fiber::FiberBond& b : topology.bonds)
    {
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(state.material_id[b.i])];
        if (!mat.inelastic()) continue;
        b.st.eps_prev = fiber_bond_length(state, b) / b.l0 - 1.0;
    }
}

// ----------------------------------------------------------------------
// Разрыв связей Fiber-Bond (fiber/SPH.md, «Fiber-Bond», «Разрыв связи»;
// fiber/PLAN.md §3.8): критерий по расстоянию l_ij > l_break = (1+eps_b)*l0,
// проверяется ПОСЛЕ шага ODESolver (координаты бусин -- уже финальные,
// скорректированные fb.x этого шага).
//
// Две фазы (fiber/SPH.md, «Postprocessing и изменение topology»):
//   1. find_fiber_break_events  -- parallel_for ТОЛЬКО НА ЧТЕНИЕ по
//      topology.bonds, потокобезопасные локальные буферы событий (тот же
//      паттерн "per-thread buffer -> merge", что в
//      interactions/fiber_bond.hpp::compute_fiber_bond_contributions и
//      rhs/contributions.hpp);
//   2. apply_fiber_break_events -- строго последовательно.
// ----------------------------------------------------------------------
// Причина разрыва (только диагностика; на применение разрыва не влияет).
// Если сработало несколько критериев, берётся первый по порядку перечисления.
enum class FiberBreakReason : std::int64_t
{
    Strain = 0,  // l > l_break (eps_break, прежний критерий)
    PlasticStrain, // eps_p >= FAILURE.EPS_P
    Stress         // sigma_ep = E(eps - eps_p) >= FAILURE.SIGMA_U
};

struct FiberBreakEvent
{
    std::size_t bond_index = 0;
    FiberBreakReason reason = FiberBreakReason::Strain;
};

// find_fiber_break_events -- READONLY: читает ТЕКУЩИЕ (после шага ODESolver)
// координаты бусин state.x[][] и l_break каждой связи, ничего не
// мутирует. thread_buffers переиспользуется между вызовами ("память только
// растёт", как и другие per-thread буферы Ф1-Ф6).
namespace detail
{
// Общая реализация поиска. materials == nullptr => только прежний критерий
// l > l_break (побитово прежнее поведение).
inline std::vector<FiberBreakEvent> find_fiber_break_events_impl(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const fiber::FiberTopology& topology,
    const std::vector<fiber::FiberMaterial>* materials, std::vector<std::vector<FiberBreakEvent>>& thread_buffers)
{
    const std::size_t nworkers = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const std::size_t nbonds = topology.bonds.size();
    if (nbonds == 0) return {};

    const bool any_failure = materials && std::any_of(materials->begin(), materials->end(),
                                                      [](const fiber::FiberMaterial& m) { return m.has_failure(); });

    exec.parallel_for(nbonds, [&](std::size_t k) {
        const fiber::FiberBond& b = topology.bonds[k];
        const double l = fiber_bond_length(state, b);
        bool broken = false;
        FiberBreakReason reason = FiberBreakReason::Strain;
        if (l > b.l_break)
        {
            broken = true;
        }
        else if (any_failure)
        {
            const fiber::FiberMaterial& mat = (*materials)[static_cast<std::size_t>(state.material_id[b.i])];
            if (mat.failure.eps_p > 0.0 && b.st.eps_p >= mat.failure.eps_p)
            {
                broken = true;
                reason = FiberBreakReason::PlasticStrain;
            }
            else if (mat.failure.sigma_u > 0.0)
            {
                // Конечное состояние шага (после возврата). У материала без
                // неупругости st нулевой => sigma_ep = E*eps.
                const double eps = l / b.l0 - 1.0;
                if (mat.E * (eps - b.st.eps_p) >= mat.failure.sigma_u)
                {
                    broken = true;
                    reason = FiberBreakReason::Stress;
                }
            }
        }
        if (broken)
        {
            int w = rhs::CpuExecutor::worker_id();
            if (w < 0 || static_cast<std::size_t>(w) >= thread_buffers.size()) w = 0;
            thread_buffers[static_cast<std::size_t>(w)].push_back(FiberBreakEvent{k, reason});
        }
    });

    std::size_t total = 0;
    for (const auto& buf : thread_buffers) total += buf.size();
    std::vector<FiberBreakEvent> merged;
    merged.reserve(total);
    for (const auto& buf : thread_buffers) merged.insert(merged.end(), buf.begin(), buf.end());
    // Порядок не должен зависеть от числа потоков и распределения связей по буферам.
    std::sort(merged.begin(), merged.end(),
              [](const FiberBreakEvent& a, const FiberBreakEvent& b) { return a.bond_index < b.bond_index; });
    return merged;
}
} // namespace detail

// Прежний критерий (l > l_break) без материалов.
// ТОЛЬКО критерий l > l_break; FAILURE (EPS_P, SIGMA_U) здесь игнорируется.
// НЕ для боевого пути -- используйте перегрузку с materials / run_fiber_topology_postprocessing.
// Оставлена для драйверов, которым нужен только геометрический критерий
// (tests/test_fiber_break.cpp, tests/fiber_t0_dump.cpp).
inline std::vector<FiberBreakEvent> find_fiber_break_events(const rhs::CpuExecutor& exec,
                                                              const fiber::FiberState& state,
                                                              const fiber::FiberTopology& topology,
                                                              std::vector<std::vector<FiberBreakEvent>>& thread_buffers)
{
    return detail::find_fiber_break_events_impl(exec, state, topology, nullptr, thread_buffers);
}

// С критериями FAILURE (раздел 7 дизайна): l > l_break ИЛИ eps_p >= EPS_P ИЛИ
// E(eps - eps_p) >= SIGMA_U по состоянию b.st ПОСЛЕ apply_fiber_bond_constitutive.
// Для материалов без FAILURE результат идентичен перегрузке без материалов.
inline std::vector<FiberBreakEvent> find_fiber_break_events(const rhs::CpuExecutor& exec,
                                                              const fiber::FiberState& state,
                                                              const fiber::FiberTopology& topology,
                                                              const std::vector<fiber::FiberMaterial>& materials,
                                                              std::vector<std::vector<FiberBreakEvent>>& thread_buffers)
{
    return detail::find_fiber_break_events_impl(exec, state, topology, &materials, thread_buffers);
}

// apply_fiber_break_events -- строго последовательное применение,
// fiber/SPH.md: «ни один RHS kernel никогда не работает одновременно с
// изменением FiberTopology».
//
// Решения (задокументированы по месту, см. также сопроводительный отчёт
// этапа Ф7):
//
//   - Порядок обработки нескольких одновременных событий -- по возрастанию
//     bond_index (детерминизм; повторные индексы, если вдруг встретятся,
//     схлопываются в один).
//
//   - Какая сторона получает новый object_id -- ВСЕГДА j-сторона (частица j
//     связи (i,j) и вся цепочка next[] от неё, т.е. "хвост" нити после
//     разрыва). Простой, детерминированный, всегда корректный вариант
//     (fiber/SPH.md прямо допускает: «простой вариант "всегда j-сторона
//     получает новый id" тоже корректен, если он последователен»).
//     Корректно работает и при КАСКАДНОМ разрыве -- если после первого
//     разрыва в обработке (по возрастанию bond_index) встречается второй
//     разрыв ДАЛЬШЕ по цепочке j-стороны, его "хвост" (начиная со своего j)
//     получает СВОЙ, ЕЩЁ более новый id, "перекрывая" уже присвоенный первым
//     разрывом -- итоговое разбиение на отрезки не зависит от порядка
//     обработки (независимая проверка рассуждением для обоих порядков
//     возрастания/убывания bond_index -- см. отчёт этапа).
//
//   - Новый id -- (текущий максимум object_id по всем бусинам) + 1,
//     пересчитывается ОДИН раз в начале вызова и далее просто
//     инкрементируется на каждое обрабатываемое событие (без выделенного
//     счётчика в Model/FiberState -- разрывы редки, а число бусин конечно,
//     полный проход дёшев и не требует персистентного состояния между
//     шагами/сериализации).
//
//   - Связи удаляются ИЗ topology.bonds (не помечаются флагом active) --
//     дешёвый одноразовый rebuild-фильтр по всему вектору вместо
//     erase()/swap-pop в цикле (который потребовал бы аккуратной обработки
//     инвалидации индексов при НЕСКОЛЬКИХ одновременных событиях). Fiber-Bond
//     kernel (interactions/fiber_bond.hpp) не хранит собственного кэша
//     primitives между шагами -- он строится заново из topology.bonds на
//     КАЖДОМ вызове RHS, поэтому удаление элемента из bonds здесь уже
//     достаточно ("interaction primitives... построены заново перед
//     следующим RHS", fiber/PLAN.md §3.8) -- fiber_bond.hpp трогать не
//     потребовалось.
//
//   - topology.version увеличивается ОДИН раз на весь проход Postprocessing
//     (не на каждое отдельное событие) -- fiber/PLAN.md §3.8 говорит об
//     инвалидации топологии как о ЕДИНОМ факте "топология этого прохода
//     изменилась", а не о подсчёте числа разрывов; читателей topology.version
//     в проекте пока нет (задел на будущее кеширование), поэтому выбор не
//     наблюдаем ни одним тестом -- документируется как осознанный.
inline void apply_fiber_break_events(fiber::FiberState& state, fiber::FiberTopology& topology,
                                      const std::vector<FiberBreakEvent>& events)
{
    if (events.empty()) return;

    std::vector<std::size_t> break_indices;
    break_indices.reserve(events.size());
    for (const auto& e : events) break_indices.push_back(e.bond_index);
    std::sort(break_indices.begin(), break_indices.end());
    break_indices.erase(std::unique(break_indices.begin(), break_indices.end()), break_indices.end());

    int next_object_id = 0;
    for (int id : state.object_id) next_object_id = std::max(next_object_id, id);
    ++next_object_id;

    // 1) Обновить prev/next и присвоить новый object_id j-стороне каждой
    //    рвущейся связи (в порядке возрастания bond_index -- см. комментарий
    //    выше).
    for (std::size_t bidx : break_indices)
    {
        const fiber::FiberBond& b = topology.bonds[bidx];
        const std::size_t i = b.i;
        const std::size_t j = b.j;

        topology.next[i] = -1;
        topology.prev[j] = -1;

        const int new_id = next_object_id++;
        long long cur = static_cast<long long>(j);
        while (cur != -1)
        {
            state.object_id[static_cast<std::size_t>(cur)] = new_id;
            cur = topology.next[static_cast<std::size_t>(cur)];
        }
    }

    // 2) Удалить порвавшиеся связи из topology.bonds (rebuild-фильтр,
    //    корректно для любого числа одновременных событий).
    std::vector<char> removed(topology.bonds.size(), 0);
    for (std::size_t bidx : break_indices) removed[bidx] = 1;
    std::vector<fiber::FiberBond> kept;
    kept.reserve(topology.bonds.size() - break_indices.size());
    for (std::size_t k = 0; k < topology.bonds.size(); ++k)
    {
        if (!removed[k]) kept.push_back(topology.bonds[k]);
    }
    topology.bonds.swap(kept);

    ++topology.version;
}

// run_fiber_topology_postprocessing -- полный проход Postprocessing нитей за
// один вызов: (1) apply_fiber_bond_constitutive (внутренние переменные связей,
// с dt шага, по x^{n+1}), (2) поиск и применение разрывов (по образцу
// compute_fiber_bond в interactions/fiber_bond.hpp). Конститутивный проход
// ВСТРОЕН сюда, чтобы любой драйвер, вызывающий обёртку (боевой цикл,
// tests/fiber_t0_dump), автоматически обновлял eps_p/q/w_p/w_v и не мог его
// забыть. Драйверы, вызывающие find_fiber_break_events/apply_fiber_break_events
// напрямую (tests/test_fiber_break.cpp), конститутивный проход пропускают --
// им нужен явный вызов apply_fiber_bond_constitutive (для упругих материалов он
// и так no-op). Для упругих материалов поведение обёртки прежнее побитово.
inline void run_fiber_topology_postprocessing(const rhs::CpuExecutor& exec, fiber::FiberState& state,
                                               fiber::FiberTopology& topology,
                                               const std::vector<fiber::FiberMaterial>& materials, double dt,
                                               std::vector<std::vector<FiberBreakEvent>>& thread_buffers)
{
    apply_fiber_bond_constitutive(exec, state, topology, materials, dt);
    const std::vector<FiberBreakEvent> events =
        find_fiber_break_events(exec, state, topology, materials, thread_buffers);
    apply_fiber_break_events(state, topology, events);
}

} // namespace interactions
