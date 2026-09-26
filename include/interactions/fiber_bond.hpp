#pragma once
// interactions/fiber_bond.hpp
// Fiber-Bond: сила упругой связи между топологическими соседями одной нити
// (fiber/SPH.md, раздел «Fiber-Bond», «Постановка проекта: только
// растяжение»; fiber/PLAN.md, этап Ф3, §3.2/§3.3/§3.7).
//
// Primitive — уже существующая fiber::FiberBond (i, j, l0, l_break),
// fiber/fiber_system.hpp, этап Ф2: список связей строится из topology-графа
// нити (FiberTopology::bonds), а не пространственным поиском (PLAN.md §3.2,
// таблица "Типы взаимодействий и primitives").
//
// Kernel только ЧИТАЕТ FiberState/FiberMaterial (readonly) и производит
// typed contributions (FiberDelta) через rhs::ContributionWriter —
// суммирование происходит отдельно, в rhs::reduce_contributions. Никаких
// atomics/mutex здесь нет (fiber/SPH.md, «Основные требования к interaction
// kernel», пп. 1-5).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"

namespace interactions
{

// Метка interaction_type для ContributionMeta/детерминированного reduce
// (rhs/contributions.hpp, «Детерминизм») — произвольная константа, уникальная
// среди типов взаимодействий, вносящих вклад в Fiber.
constexpr std::uint32_t kFiberBondInteractionType = 1;

// Площадь поперечного сечения нити (fiber/SPH.md, Fiber-Bond): A = pi*d^2/4.
// ВНИМАНИЕ: модуль Юнга нити — FiberMaterial.E, не FiberMaterial.nu
// (коэффициент Пуассона, в Fiber-Bond не используется — постановка проекта
// работает только с осевой жёсткостью EA, изгиб не моделируется).
inline real_t fiber_cross_section_area(const fiber::FiberMaterial& mat)
{
    return pi * mat.d * mat.d / real_t(4.0);
}

// ----------------------------------------------------------------------
// compute_fiber_bond_contributions — заполняет thread_buffers (по одному
// на CPU worker, как в rhs/contributions.hpp) contributions Fiber-Bond для
// ВСЕХ связей topology.bonds. Буферы переиспользуются между вызовами
// (изменяют размер под exec.num_workers() при необходимости, иначе только
// clear() — «память только растёт», fiber/SPH.md «Временная память»).
//
// Для каждой связи (i, j, l0, l_break):
//   x_ij = x_j - x_i,  l_ij = |x_ij|
//   F_i<-j = E*A * max(l_ij/l0 - 1, 0) * x_ij/l_ij      (только растяжение)
//   F_j<-i = -F_i<-j
//   dv_i = F_i<-j / m_i,  dv_j = F_j<-i / m_j = -F_i<-j / m_j
// l_break здесь не используется — критерий разрыва проверяется в
// Postprocessing (этап Ф7, fiber/SPH.md «Разрыв связи»), НЕ в этом kernel.
// ----------------------------------------------------------------------
inline void compute_fiber_bond_contributions(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const fiber::FiberTopology& topology,
    const std::vector<fiber::FiberMaterial>& materials,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers)
{
    const std::size_t nworkers = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const std::size_t nbonds = topology.bonds.size();
    if (nbonds == 0) return;

    exec.parallel_for(nbonds, [&](std::size_t k) {
        const fiber::FiberBond& b = topology.bonds[k];
        const std::size_t i = b.i;
        const std::size_t j = b.j;
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(state.material_id[i])];
        const real_t A = fiber_cross_section_area(mat);

        real_t xij[dim];
        real_t l2 = 0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            xij[d] = state.x[d][j] - state.x[d][i];
            l2 += xij[d] * xij[d];
        }
        const real_t l = std::sqrt(l2);

        fiber::FiberDelta delta_i{};
        fiber::FiberDelta delta_j{};

        if (l > b.l0 && l > 0) // только растяжение (fiber/SPH.md, "Постановка проекта")
        {
            const real_t strain = l / b.l0 - real_t(1.0);
            const real_t coeff = mat.E * A * strain / l; // F_i<-j = coeff * xij
            const real_t mi = state.mass[i];
            const real_t mj = state.mass[j];
            for (std::size_t d = 0; d < dim; ++d)
            {
                const real_t f = coeff * xij[d];
                delta_i.dv[d] = f / mi;  // +F_i<-j / m_i
                delta_j.dv[d] = -f / mj; // -F_i<-j / m_j == F_j<-i / m_j
            }
        }

        int w = rhs::CpuExecutor::worker_id();
        if (w < 0 || static_cast<std::size_t>(w) >= thread_buffers.size()) w = 0;
        rhs::ContributionWriter<fiber::FiberDelta> writer(thread_buffers[static_cast<std::size_t>(w)],
                                                            kFiberBondInteractionType);
        writer.emit(i, delta_i, /*primitive_id=*/static_cast<std::uint64_t>(k), /*slot=*/0);
        writer.emit(j, delta_j, /*primitive_id=*/static_cast<std::uint64_t>(k), /*slot=*/1);
    });
}

// compute_fiber_bond — удобная обёртка: contributions + reduce за один
// вызов. Возвращает вектор FiberDelta длины state.size() (индекс — частица
// Fiber). thread_buffers передаётся вызывающей стороной (обычно —
// переиспользуемый буфер из RhsWorkspace/локальной переменной цикла) для
// соблюдения "память только растёт".
inline std::vector<fiber::FiberDelta> compute_fiber_bond(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const fiber::FiberTopology& topology,
    const std::vector<fiber::FiberMaterial>& materials,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers)
{
    compute_fiber_bond_contributions(exec, state, topology, materials, thread_buffers);
    return rhs::reduce_contributions(exec, thread_buffers, state.size());
}

// ----------------------------------------------------------------------
// compute_fiber_bond_dt_limit — dt_limit связей (fiber/SPH.md, "Шаг по
// времени"; fiber/PLAN.md §3.7):
//     dt_lim = 0.8 / sqrt(k_eff / m_red),  k_eff = E*A/l0,
//     m_red = m_i*m_j / (m_i+m_j).
// Считается для ВСЕХ связей независимо от знака деформации ("чтобы шаг не
// прыгал при переходе через l0") — минимум по всем связям через
// exec.parallel_reduce. Возвращает +inf (numeric_limits::max), если связей
// нет — так вызывающая сторона может безусловно брать min(..., это
// значение) без отдельной проверки на пустоту.
// ----------------------------------------------------------------------
inline real_t compute_fiber_bond_dt_limit(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                           const fiber::FiberTopology& topology,
                                           const std::vector<fiber::FiberMaterial>& materials)
{
    const std::size_t nbonds = topology.bonds.size();
    const real_t kInf = std::numeric_limits<real_t>::max();
    if (nbonds == 0) return kInf;

    return exec.parallel_reduce(
        nbonds, kInf,
        [&](std::size_t k) -> real_t {
            const fiber::FiberBond& b = topology.bonds[k];
            const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(state.material_id[b.i])];
            const real_t A = fiber_cross_section_area(mat);
            const real_t k_eff = mat.E * A / b.l0;
            const real_t mi = state.mass[b.i];
            const real_t mj = state.mass[b.j];
            const real_t m_red = (mi * mj) / (mi + mj);
            return real_t(0.8) / std::sqrt(k_eff / m_red);
        },
        [](real_t a, real_t b2) { return std::min(a, b2); });
}

} // namespace interactions
