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
//                                                  run_fiber_topology_postprocessing.
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
                                     const std::vector<double>& drho
#ifdef SPH_ENERGY_LOG
                                     ,
                                     double* plastic_acc = nullptr
#endif
)
{
    for (std::size_t i = 0; i < ntotal; ++i)
    {
        if (sys.itype[i] >= 3 && sys.itype[i] <= 9)
        {
#ifdef SPH_ENERGY_LOG
            double ss_before = 0.0;
            if (plastic_acc)
                ss_before = sys.sxx[i] * sys.sxx[i] + sys.syy[i] * sys.syy[i] + sys.szz[i] * sys.szz[i] +
                            2.0 * (sys.sxy[i] * sys.sxy[i] + sys.syz[i] * sys.syz[i] + sys.sxz[i] * sys.sxz[i]);
#endif
            prochnost_mises(cfg.mats, static_cast<size_t>(sys.itype[i]), sys.sxx[i], sys.syy[i], sys.szz[i],
                             sys.sxy[i], sys.syz[i], sys.sxz[i]);
#ifdef SPH_ENERGY_LOG
            // Радиальный возврат Мизеса: упругая энергия сдвига (s:s)/(4G rho)
            // переходит в тепло. Она уже входит в Σm·u (справочная величина,
            // в баланс не добавляется). JH-2 (itype 10-11) не учитывается.
            if (plastic_acc)
            {
                const double ss_after =
                    sys.sxx[i] * sys.sxx[i] + sys.syy[i] * sys.syy[i] + sys.szz[i] * sys.szz[i] +
                    2.0 * (sys.sxy[i] * sys.sxy[i] + sys.syz[i] * sys.syz[i] + sys.sxz[i] * sys.sxz[i]);
                const double G = shear_modulus_by_itype(cfg.mats, sys.itype[i]);
                if (G > 0.0 && sys.rho[i] > 0.0)
                    *plastic_acc += sys.mass[i] * (ss_before - ss_after) / (4.0 * G * sys.rho[i]);
            }
#endif
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
struct FiberBreakEvent
{
    std::size_t bond_index = 0;
};

// find_fiber_break_events -- READONLY: читает ТЕКУЩИЕ (после шага ODESolver)
// координаты бусин state.x[][] и l_break каждой связи, ничего не
// мутирует. thread_buffers переиспользуется между вызовами ("память только
// растёт", как и другие per-thread буферы Ф1-Ф6).
inline std::vector<FiberBreakEvent> find_fiber_break_events(const rhs::CpuExecutor& exec,
                                                              const fiber::FiberState& state,
                                                              const fiber::FiberTopology& topology,
                                                              std::vector<std::vector<FiberBreakEvent>>& thread_buffers)
{
    const std::size_t nworkers = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const std::size_t nbonds = topology.bonds.size();
    if (nbonds == 0) return {};

    exec.parallel_for(nbonds, [&](std::size_t k) {
        const fiber::FiberBond& b = topology.bonds[k];
        double l2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double dx = state.x[d][b.j] - state.x[d][b.i];
            l2 += dx * dx;
        }
        const double l = std::sqrt(l2);
        if (l > b.l_break)
        {
            int w = rhs::CpuExecutor::worker_id();
            if (w < 0 || static_cast<std::size_t>(w) >= thread_buffers.size()) w = 0;
            thread_buffers[static_cast<std::size_t>(w)].push_back(FiberBreakEvent{k});
        }
    });

    std::size_t total = 0;
    for (const auto& buf : thread_buffers) total += buf.size();
    std::vector<FiberBreakEvent> merged;
    merged.reserve(total);
    for (const auto& buf : thread_buffers) merged.insert(merged.end(), buf.begin(), buf.end());
    return merged;
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
                                      const std::vector<FiberBreakEvent>& events
#ifdef SPH_ENERGY_LOG
                                      ,
                                      const std::vector<fiber::FiberMaterial>* materials = nullptr,
                                      double* break_acc = nullptr
#endif
)
{
    if (events.empty()) return;

    std::vector<std::size_t> break_indices;
    break_indices.reserve(events.size());
    for (const auto& e : events) break_indices.push_back(e.bond_index);
    std::sort(break_indices.begin(), break_indices.end());
    break_indices.erase(std::unique(break_indices.begin(), break_indices.end()), break_indices.end());

#ifdef SPH_ENERGY_LOG
    // Упругая энергия порвавшихся связей, 0.5*k*(l-l0)^2, k = E*A/l0: она исчезает
    // из модели при удалении связи. Считается по дедуплицированным индексам ДО
    // любой мутации topology.bonds (индексы ещё валидны).
    if (materials && break_acc)
    {
        for (std::size_t bidx : break_indices)
        {
            const fiber::FiberBond& b = topology.bonds[bidx];
            double l2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double dx = state.x[d][b.j] - state.x[d][b.i];
                l2 += dx * dx;
            }
            const double l = std::sqrt(l2);
            if (l > b.l0)
            {
                const fiber::FiberMaterial& mat = (*materials)[static_cast<std::size_t>(state.material_id[b.i])];
                const double k = mat.E * fiber_cross_section_area(mat) / b.l0;
                *break_acc += 0.5 * k * (l - b.l0) * (l - b.l0);
            }
        }
    }
#endif

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

// run_fiber_topology_postprocessing -- удобная обёртка: поиск + применение
// за один вызов (по образцу compute_fiber_bond в interactions/fiber_bond.hpp).
inline void run_fiber_topology_postprocessing(const rhs::CpuExecutor& exec, fiber::FiberState& state,
                                               fiber::FiberTopology& topology,
                                               std::vector<std::vector<FiberBreakEvent>>& thread_buffers
#ifdef SPH_ENERGY_LOG
                                               ,
                                               const std::vector<fiber::FiberMaterial>* materials = nullptr,
                                               double* break_acc = nullptr
#endif
)
{
    const std::vector<FiberBreakEvent> events = find_fiber_break_events(exec, state, topology, thread_buffers);
#ifdef SPH_ENERGY_LOG
    apply_fiber_break_events(state, topology, events, materials, break_acc);
#else
    apply_fiber_break_events(state, topology, events);
#endif
}

} // namespace interactions
