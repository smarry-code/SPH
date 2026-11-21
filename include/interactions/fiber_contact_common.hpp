#pragma once
// interactions/fiber_contact_common.hpp
// Общие константы/утилиты для Fiber-Fiber-Contact (interactions/
// fiber_contact.hpp) и Fiber-Solid-Contact (interactions/
// fiber_solid_contact.hpp): обе задачи (notes/fiber-contact/fiber-contact-segment-pairs.md
// и notes/fiber-contact/fiber-solid-contact-fix.md) используют одно и то же прижатие
// параметра к границе сегмента (η) и один и тот же вывод неравенства
// треугольника для радиуса поиска (текущая, не номинальная l0!, максимальная
// длина связи L_max) -- вынесены сюда, чтобы не заводить вторую константу и
// не дублировать вычисление (notes/fiber-contact/fiber-solid-contact-fix.md, §1 п.2 и §2).
#include <algorithm>
#include <cstddef>

#include "fiber/fiber_system.hpp"
#include "param.hpp"
#include "rhs/cpu_executor.hpp"

namespace interactions
{

// Snapping параметров ближайших точек к 0/1 у границы сегмента нити
// (перенесено из fiber_contact.hpp без изменений; fiber/
// fiber-contact-segment-pairs.md, §2 спецификации, η): без него пересечение
// ровно в узле сети то классифицируется как внутренняя точка сегмента
// (s~1-1e-16), то как его конец -- см. шапки fiber_contact.hpp и
// fiber_solid_contact.hpp.
constexpr double kFiberContactSnapEta = 1e-9;

// ----------------------------------------------------------------------
// fiber_max_current_bond_length -- L_max, максимум ТЕКУЩЕЙ (не номинальной
// l0!) длины |x_i-x_j| по всем связям topology.bonds расчёта. Связи Fiber-
// Bond работают только на растяжение и рвутся лишь в Postprocessing ПОСЛЕ
// шага (notes/fiber-model/fiber-problem-statement.md), поэтому внутри шага сегмент может быть длиннее l0 (и
// даже длиннее l_break) -- номинальное l0 как оценка радиуса поиска не
// годится (notes/fiber-contact/fiber-contact-segment-pairs.md, пункт 3 в шапке файла;
// notes/fiber-contact/fiber-solid-contact-fix.md, §2). Используется ОБОИМИ контактами
// (Fiber-Fiber, interactions/fiber_contact.hpp; Fiber-Solid, interactions/
// fiber_solid_contact.hpp) -- вычисление ровно одно, порядок операций и
// reduce тот же, что был в fiber_contact_search_radius ДО выноса сюда
// (побитовая воспроизводимость -- Fiber-Fiber-Contact не должен измениться
// ни на бит).
// ----------------------------------------------------------------------
inline double fiber_max_current_bond_length(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                             const fiber::FiberTopology& topology)
{
    const std::size_t nb = topology.bonds.size();
    return exec.parallel_reduce(
        nb, 0.0,
        [&](std::size_t k) -> double {
            const fiber::FiberBond& b = topology.bonds[k];
            if (b.mat >= 0) return 0.0; // связь-стежок -- не сегмент нити: радиус поиска от неё не зависит
            double d2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double dx = state.x[d][b.i] - state.x[d][b.j];
                d2 += dx * dx;
            }
            return std::sqrt(d2);
        },
        [](double a, double b) { return std::max(a, b); });
}

} // namespace interactions
