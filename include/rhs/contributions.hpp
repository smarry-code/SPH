#ifndef RHS_CONTRIBUTIONS_HPP
#define RHS_CONTRIBUTIONS_HPP

#include <algorithm>
#include <cstdint>
#include <vector>

#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"

// ============================================================================
// rhs/contributions.hpp — Contribution/ContributionBuffer/ContributionWriter
// и reduce по шардам owner-computes (fiber/SPH.md, разделы "Compute
// interaction kernel", "Reduce", "Детерминизм"; fiber/PLAN.md §3.5).
//
// Interaction kernel не пишет Delta напрямую в общий RHS частицы (гонка
// данных при нескольких primitives, влияющих на одну частицу). Вместо этого
// kernel формирует Contribution через ContributionWriter, принадлежащий
// своему worker'у; contributions складываются на отдельной стадии Reduce по
// схеме owner-computes: диапазон particle id режется на shard'ы, каждый
// shard — единственный владелец своего диапазона и суммирует его без
// mutex/atomic.
// ============================================================================

namespace rhs
{

// Метаданные вклада, нужные ТОЛЬКО для канонического порядка суммирования в
// детерминированном режиме (SPH.md, "Детерминизм": interaction_type,
// primitive_id, contribution_slot).
struct ContributionMeta
{
    std::uint32_t interaction_type = 0;
    std::uint64_t primitive_id = 0;
    std::uint16_t slot = 0;
};

template <class Delta>
struct Contribution
{
    index_t target;
    Delta delta;
    ContributionMeta meta;
};

// ----------------------------------------------------------------------
// ContributionBuffer<Delta> — буфер вкладов ОДНОГО worker'а для ОДНОЙ
// compute-структуры. Часть RhsWorkspace: переиспользуется между вызовами
// RHS — clear() сбрасывает размер, но не освобождает память (SPH.md,
// "Временная память": "не следует выполнять большое количество небольших
// heap allocations внутри... interaction kernel").
// ----------------------------------------------------------------------
template <class Delta>
class ContributionBuffer
{
public:
    void clear() { data_.clear(); }
    void push(index_t target, const Delta& delta, ContributionMeta meta)
    {
        data_.push_back(Contribution<Delta>{target, delta, meta});
    }
    std::size_t size() const { return data_.size(); }
    bool empty() const { return data_.empty(); }
    std::size_t capacity() const { return data_.capacity(); }

    const Contribution<Delta>& operator[](std::size_t i) const { return data_[i]; }
    Contribution<Delta>& operator[](std::size_t i) { return data_[i]; }

private:
    std::vector<Contribution<Delta>> data_;
};

// ----------------------------------------------------------------------
// ContributionWriter<Delta> — интерфейс kernel'а на одну compute-структуру:
//     output.emit(target, delta, primitive_id[, slot]);
// соответствует SPH.md: "output.emit<Fiber>(fiber_id, fiber_delta)".
//
// Инженерное отклонение от буквального текста SPH.md: там один writer умеет
// emit<Fiber>/emit<Solid> сразу для нескольких compute-структур. На Ф1 нет ни
// Fiber, ни Solid (задача явно запрещает физику) — вводить полиморфный по
// compute-структуре writer нечем наполнить и нечем протестировать.
// ContributionWriter<Delta> — writer ОДНОЙ структуры; kernel, вносящий вклад
// сразу в несколько структур (например, Fiber-Solid-Contact — в Fiber и в
// Solid), получает несколько таких writer'ов, по одному на структуру. Эта
// комбинация тривиальна и будет собрана в Ф2+/Ф6, когда появятся реальные
// FiberDelta/SolidDelta.
// ----------------------------------------------------------------------
template <class Delta>
class ContributionWriter
{
public:
    ContributionWriter(ContributionBuffer<Delta>& buf, std::uint32_t interaction_type)
        : buf_(buf), interaction_type_(interaction_type)
    {
    }

    void emit(index_t target, const Delta& delta, std::uint64_t primitive_id, std::uint16_t slot = 0)
    {
        buf_.push(target, delta, ContributionMeta{interaction_type_, primitive_id, slot});
    }

private:
    ContributionBuffer<Delta>& buf_;
    std::uint32_t interaction_type_;
};

// ----------------------------------------------------------------------
// reduce_contributions — CPU reduce по схеме owner-computes (SPH.md,
// "Reduce"): диапазон [0, particle_count) режется на shard'ы (заметно
// мельче числа worker'ов), вклады раскладываются по shard'ам двухпроходной
// схемой count -> prefix sum -> copy (та же схема, что у Search), затем
// каждый shard суммирует свой диапазон независимо и параллельно другим
// shard'ам.
//
// В детерминированном режиме (options.deterministic) вклады внутри shard'а
// предварительно сортируются по (target, interaction_type, primitive_id,
// slot) — результат не зависит от того, сколько потоков и в каком порядке
// произвели contributions (SPH.md, "Детерминизм").
//
// Delta должен быть default-конструируем в "нулевой" элемент суммы и
// поддерживать operator+=.
// ----------------------------------------------------------------------
struct ReduceOptions
{
    std::size_t shard_count = 0; // 0 => подобрать автоматически
    bool deterministic = false;
};

template <class Delta>
std::vector<Delta> reduce_contributions(const CpuExecutor& exec,
                                         const std::vector<ContributionBuffer<Delta>>& thread_buffers,
                                         std::size_t particle_count,
                                         const ReduceOptions& options = ReduceOptions{})
{
    std::vector<Delta> output(particle_count, Delta{});
    if (particle_count == 0) return output;

    std::size_t total = 0;
    for (const auto& buf : thread_buffers) total += buf.size();
    if (total == 0) return output;

    std::size_t shard_count = options.shard_count;
    if (shard_count == 0)
    {
        shard_count = static_cast<std::size_t>(std::max(1, exec.num_workers())) * 4;
    }
    shard_count = std::max<std::size_t>(1, std::min(shard_count, particle_count));

    // Непрерывные границы shard'ов: [0,b1), [b1,b2), ... (SPH.md, "Reduce").
    std::vector<std::size_t> shard_bound(shard_count + 1);
    for (std::size_t s = 0; s <= shard_count; ++s)
    {
        shard_bound[s] = (particle_count * s) / shard_count;
    }

    auto shard_of = [&shard_bound, shard_count](index_t target) -> std::size_t {
        auto it = std::upper_bound(shard_bound.begin(), shard_bound.end(), target);
        std::size_t s = static_cast<std::size_t>(it - shard_bound.begin());
        return (s == 0) ? 0 : std::min(s - 1, shard_count - 1);
    };

    const std::size_t nthreads = thread_buffers.size();

    // Проход 1: подсчёт (thread, shard) -> количество.
    std::vector<std::size_t> counts(nthreads * shard_count, 0);
    for (std::size_t t = 0; t < nthreads; ++t)
    {
        const auto& buf = thread_buffers[t];
        for (std::size_t k = 0; k < buf.size(); ++k)
        {
            const std::size_t s = shard_of(buf[k].target);
            counts[t * shard_count + s]++;
        }
    }

    // Префиксная сумма: сначала непрерывный блок shard'а 0 (по всем thread'ам
    // по возрастанию thread id), затем shard'а 1 и т.д.
    std::vector<std::size_t> shard_total(shard_count, 0);
    for (std::size_t s = 0; s < shard_count; ++s)
    {
        for (std::size_t t = 0; t < nthreads; ++t) shard_total[s] += counts[t * shard_count + s];
    }

    std::vector<std::size_t> shard_start(shard_count + 1, 0);
    for (std::size_t s = 0; s < shard_count; ++s) shard_start[s + 1] = shard_start[s] + shard_total[s];

    std::vector<std::size_t> dest_offset(nthreads * shard_count, 0);
    for (std::size_t s = 0; s < shard_count; ++s)
    {
        std::size_t running = shard_start[s];
        for (std::size_t t = 0; t < nthreads; ++t)
        {
            dest_offset[t * shard_count + s] = running;
            running += counts[t * shard_count + s];
        }
    }

    // Проход 2: scatter в непрерывный массив merged, сгруппированный по shard'ам.
    std::vector<Contribution<Delta>> merged(total);
    {
        std::vector<std::size_t> cursor = dest_offset;
        for (std::size_t t = 0; t < nthreads; ++t)
        {
            const auto& buf = thread_buffers[t];
            for (std::size_t k = 0; k < buf.size(); ++k)
            {
                const std::size_t s = shard_of(buf[k].target);
                merged[cursor[t * shard_count + s]++] = buf[k];
            }
        }
    }

    // Проход 3: суммирование по shard'ам, параллельно между shard'ами.
    exec.parallel_for(shard_count, std::size_t(1), [&](std::size_t s) {
        const std::size_t begin = shard_start[s];
        const std::size_t end = shard_start[s + 1];
        if (options.deterministic)
        {
            std::sort(merged.begin() + static_cast<std::ptrdiff_t>(begin),
                      merged.begin() + static_cast<std::ptrdiff_t>(end),
                      [](const Contribution<Delta>& a, const Contribution<Delta>& b) {
                          if (a.target != b.target) return a.target < b.target;
                          if (a.meta.interaction_type != b.meta.interaction_type)
                              return a.meta.interaction_type < b.meta.interaction_type;
                          if (a.meta.primitive_id != b.meta.primitive_id)
                              return a.meta.primitive_id < b.meta.primitive_id;
                          return a.meta.slot < b.meta.slot;
                      });
        }
        for (std::size_t k = begin; k < end; ++k)
        {
            output[merged[k].target] += merged[k].delta;
        }
    });

    return output;
}

} // namespace rhs

#endif // RHS_CONTRIBUTIONS_HPP
