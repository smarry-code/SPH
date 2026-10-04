#ifndef RHS_CONTRIBUTIONS_HPP
#define RHS_CONTRIBUTIONS_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
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
    const Contribution<Delta>* data() const { return data_.data(); }

private:
    std::vector<Contribution<Delta>> data_;
};

// ----------------------------------------------------------------------
// IndexedContributions<Delta> — постоянный массив вкладов с записью по
// индексу, для kernel'ов, у которых число вкладов известно заранее (например,
// Solid-SPH: примитив даёт вклад для i и, если j реальна, для j; позиции
// считаются префиксной суммой). Нет общего заголовка, в который пишут все
// потоки (у ContributionBuffer заголовки vector'ов соседних потоков делят
// кэш-линию — false sharing), а порядок вкладов канонический и не зависит от
// числа потоков. prepare(n) только растит ёмкость и не обнуляет данные:
// память переиспользуется между шагами.
//
// Рядом с вкладами хранится компактный массив целей targets_ (4 Б на вклад):
// проходы 1-2 свёртки читают из вклада только target, а вклад занимает 176 Б,
// т.е. целую кэш-линию на 4 полезных байта; чтение компактного массива снимает
// этот упор в память. Писать вклады можно только через writer().set(...), чтобы
// цель и вклад не разошлись; неконстантных data()/operator[] нет.
// ----------------------------------------------------------------------
template <class Delta>
class IndexedContributions
{
public:
    void prepare(std::size_t n)
    {
        if (n > data_.size()) data_.resize(n);
        if (n > targets_.size()) targets_.resize(n);
        size_ = n;
    }
    std::size_t size() const { return size_; }
    std::size_t capacity() const { return data_.size(); }

    const Contribution<Delta>* data() const { return data_.data(); }
    const Contribution<Delta>& operator[](std::size_t i) const { return data_[i]; }
    const index_t* targets() const { return targets_.data(); }

    // Единственный способ записи: вклад и его цель пишутся вместе. Сырые
    // указатели -- чтобы цикл kernel'а не перечитывал заголовки vector.
    struct Writer
    {
        Contribution<Delta>* d;
        index_t* t;
        void set(std::size_t pos, const Contribution<Delta>& c) const
        {
            d[pos] = c;
            t[pos] = c.target;
        }
    };
    Writer writer() { return Writer{data_.data(), targets_.data()}; }

    // Постоянные вспомогательные массивы kernel'а (число вкладов на примитив
    // и смещения); тоже только растут.
    std::vector<std::size_t> counts, offsets;

private:
    std::vector<Contribution<Delta>> data_;
    std::vector<index_t> targets_; // targets_[k] == data_[k].target
    std::size_t size_ = 0;
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
// мельче числа worker'ов), ссылки на вклады раскладываются по shard'ам
// двухпроходной схемой count -> prefix sum -> scatter (та же схема, что у
// Search; сами вклады не копируются), затем
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

// Постоянная рабочая память reduce: ReduceScratch живёт между вызовами (в
// RhsWorkspace/решателе), все массивы только растут — повторные вызовы с тем
// же или меньшим числом вкладов не обращаются к аллокатору (память только
// растёт). Не зависит от Delta: буферы и вклады хранятся как const void*.
struct ReduceScratch
{
    // Отрезки вкладов в порядке суммирования: начало (Contribution<Delta>*) и размер.
    std::vector<const void*> seg_ptr;
    std::vector<std::size_t> seg_size;
    // Компактные цели отрезка (const index_t*); nullptr -- цель читается из самого вклада.
    std::vector<const index_t*> seg_tgt;
    std::vector<std::size_t> shard_bound; // границы shard'ов, shard_count + 1
    std::vector<std::size_t> shard_start; // начало блока shard'а в refs, shard_count + 1
    std::vector<std::size_t> counts;      // (буфер, shard) -> количество вкладов
    std::vector<std::size_t> dest_offset; // (буфер, shard) -> позиция записи в refs
    // Ссылки на вклады (указатели на Contribution<Delta> внутри буферов),
    // сгруппированные по shard'ам; вместо прежнего массива полных копий вкладов.
    std::vector<const void*> refs;
};

// Номер shard'а по target: наибольшее s с shard_bound[s] <= target, но не
// больше shard_count - 1 (target >= particle_count -> последний shard) —
// ровно семантика прежнего std::upper_bound. Начальная оценка s через double
// уточняется по shard_bound шагами +-1; произведения target * shard_count
// (возможное переполнение) нет.
inline std::size_t reduce_shard_of(const std::size_t* shard_bound, std::size_t shard_count,
                                    std::size_t particle_count, double inv_width, std::size_t target)
{
    if (target >= particle_count) return shard_count - 1;
    std::size_t s = std::min(static_cast<std::size_t>(static_cast<double>(target) * inv_width), shard_count - 1);
    while (s > 0 && shard_bound[s] > target) --s;
    while (s + 1 < shard_count && shard_bound[s + 1] <= target) ++s;
    return s;
}

// Ядро свёртки: отрезки вкладов seg_ptr[0..seg_count) / seg_size[...] в
// scratch (заполняет вызывающий) логически конкатенированы в порядке номеров.
// output не должен совпадать с отрезками; ссылки в scratch живут только
// внутри вызова. Порядок суммирования на каждую частицу: shard -> номер
// отрезка -> позиция в отрезке, поэтому результат не зависит от числа потоков.
template <class Delta>
void reduce_segments_into(const CpuExecutor& exec, std::size_t nbuf, std::size_t particle_count,
                          std::vector<Delta>& output, ReduceScratch& scratch, const ReduceOptions& options)
{
    output.assign(particle_count, Delta{});
    if (particle_count == 0) return;

    std::size_t total = 0;
    for (std::size_t g = 0; g < nbuf; ++g) total += scratch.seg_size[g];
    if (total == 0) return;

    std::size_t shard_count = options.shard_count;
    if (shard_count == 0)
    {
        shard_count = static_cast<std::size_t>(std::max(1, exec.num_workers())) * 4;
    }
    shard_count = std::max<std::size_t>(1, std::min(shard_count, particle_count));

    // Непрерывные границы shard'ов: [0,b1), [b1,b2), ... (SPH.md, "Reduce").
    scratch.shard_bound.resize(shard_count + 1);
    for (std::size_t s = 0; s <= shard_count; ++s)
    {
        scratch.shard_bound[s] = (particle_count * s) / shard_count;
    }
    const std::size_t* const bound = scratch.shard_bound.data();
    const double inv_width = static_cast<double>(shard_count) / static_cast<double>(particle_count);

    const auto* const seg_ptr = scratch.seg_ptr.data();
    const std::size_t* const seg_size = scratch.seg_size.data();
    const index_t* const* const seg_tgt = scratch.seg_tgt.data();

    // Проход 1: подсчёт (буфер, shard) -> количество; параллельно по буферам,
    // каждый буфер пишет только свою строку counts.
    scratch.counts.assign(nbuf * shard_count, 0);
    std::size_t* const counts = scratch.counts.data();
    exec.parallel_for(nbuf, std::size_t(1), [&](std::size_t g) {
        const auto* const seg = static_cast<const Contribution<Delta>*>(seg_ptr[g]);
        const std::size_t n = seg_size[g];
        std::size_t* const row = counts + g * shard_count;
        const auto count_all = [&](auto target_of) {
            for (std::size_t k = 0; k < n; ++k)
                row[reduce_shard_of(bound, shard_count, particle_count, inv_width, target_of(k))]++;
        };
        // Выбор источника целей один раз на отрезок, не на элемент.
        if (const index_t* const tg = seg_tgt[g])
            count_all([tg](std::size_t k) { return tg[k]; });
        else
            count_all([seg](std::size_t k) { return seg[k].target; });
    });

    // Префиксная сумма: сначала непрерывный блок shard'а 0 (по всем буферам по
    // возрастанию номера), затем shard'а 1 и т.д.
    scratch.shard_start.resize(shard_count + 1);
    scratch.dest_offset.resize(nbuf * shard_count);
    {
        std::size_t running = 0;
        for (std::size_t s = 0; s < shard_count; ++s)
        {
            scratch.shard_start[s] = running;
            for (std::size_t g = 0; g < nbuf; ++g)
            {
                scratch.dest_offset[g * shard_count + s] = running;
                running += counts[g * shard_count + s];
            }
        }
        scratch.shard_start[shard_count] = running;
    }

    // Проход 2: scatter ссылок в refs, сгруппированный по shard'ам; пары
    // (буфер, shard) пишут непересекающиеся диапазоны, каждый буфер двигает
    // только свою строку dest_offset.
    if (scratch.refs.size() < total) scratch.refs.resize(total);
    const void** const refs = scratch.refs.data();
    std::size_t* const dest = scratch.dest_offset.data();
    exec.parallel_for(nbuf, std::size_t(1), [&](std::size_t g) {
        const auto* const seg = static_cast<const Contribution<Delta>*>(seg_ptr[g]);
        const std::size_t n = seg_size[g];
        std::size_t* const cursor = dest + g * shard_count;
        const auto scatter_all = [&](auto target_of) {
            for (std::size_t k = 0; k < n; ++k)
            {
                const std::size_t s = reduce_shard_of(bound, shard_count, particle_count, inv_width, target_of(k));
                refs[cursor[s]++] = &seg[k];
            }
        };
        if (const index_t* const tg = seg_tgt[g])
            scatter_all([tg](std::size_t k) { return tg[k]; });
        else
            scatter_all([seg](std::size_t k) { return seg[k].target; });
    });

    // Проход 3: суммирование по shard'ам, параллельно между shard'ами.
    using Contrib = Contribution<Delta>;
    const std::size_t* const start = scratch.shard_start.data();
    Delta* const out = output.data();
    exec.parallel_for(shard_count, std::size_t(1), [&](std::size_t s) {
        const void** const begin = refs + start[s];
        const void** const end = refs + start[s + 1];
        if (options.deterministic)
        {
            std::sort(begin, end, [](const void* pa, const void* pb) {
                const Contrib& a = *static_cast<const Contrib*>(pa);
                const Contrib& b = *static_cast<const Contrib*>(pb);
                if (a.target != b.target) return a.target < b.target;
                if (a.meta.interaction_type != b.meta.interaction_type)
                    return a.meta.interaction_type < b.meta.interaction_type;
                if (a.meta.primitive_id != b.meta.primitive_id)
                    return a.meta.primitive_id < b.meta.primitive_id;
                return a.meta.slot < b.meta.slot;
            });
        }
        for (const void** p = begin; p != end; ++p)
        {
            const Contrib& c = *static_cast<const Contrib*>(*p);
            out[c.target] += c.delta;
        }
    });
}

// Наборы буферов sets[0..set_count), логически конкатенированные в порядке
// перечисления (как если бы все буферы были скопированы в один vector);
// каждый буфер — один отрезок (пустой допустим).
template <class Delta>
void reduce_contributions_into(const CpuExecutor& exec,
                                const std::vector<ContributionBuffer<Delta>>* const* sets, std::size_t set_count,
                                std::size_t particle_count, std::vector<Delta>& output, ReduceScratch& scratch,
                                const ReduceOptions& options)
{
    std::size_t nbuf = 0;
    for (std::size_t k = 0; k < set_count; ++k) nbuf += sets[k]->size();
    scratch.seg_ptr.resize(nbuf);
    scratch.seg_size.resize(nbuf);
    scratch.seg_tgt.assign(nbuf, nullptr);
    std::size_t g = 0;
    for (std::size_t k = 0; k < set_count; ++k)
    {
        for (const auto& buf : *sets[k])
        {
            scratch.seg_ptr[g] = buf.data();
            scratch.seg_size[g] = buf.size();
            ++g;
        }
    }
    reduce_segments_into<Delta>(exec, nbuf, particle_count, output, scratch, options);
}

// Индексированный массив first (режется на max(1, num_workers()) смежных
// отрезков почти равного размера в порядке индексов) и затем наборы буферов
// rest в порядке перечисления. Порядок суммирования на частицу — порядок
// вкладов в first, затем rest; не зависит от числа потоков.
template <class Delta>
void reduce_contributions_into(const CpuExecutor& exec, const IndexedContributions<Delta>& first,
                                std::initializer_list<const std::vector<ContributionBuffer<Delta>>*> rest,
                                std::size_t particle_count, std::vector<Delta>& output, ReduceScratch& scratch,
                                const ReduceOptions& options = ReduceOptions{})
{
    const std::size_t m = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    std::size_t nbuf = m;
    for (const auto* set : rest) nbuf += set->size();
    scratch.seg_ptr.resize(nbuf);
    scratch.seg_size.resize(nbuf);
    scratch.seg_tgt.assign(nbuf, nullptr);
    const std::size_t n = first.size();
    for (std::size_t q = 0; q < m; ++q)
    {
        const std::size_t b = (n * q) / m, e = (n * (q + 1)) / m;
        scratch.seg_ptr[q] = first.data() + b;
        scratch.seg_size[q] = e - b;
        scratch.seg_tgt[q] = first.targets() + b;
    }
    std::size_t g = m;
    for (const auto* set : rest)
    {
        for (const auto& buf : *set)
        {
            scratch.seg_ptr[g] = buf.data();
            scratch.seg_size[g] = buf.size();
            ++g;
        }
    }
    reduce_segments_into<Delta>(exec, nbuf, particle_count, output, scratch, options);
}

// Один набор буферов.
template <class Delta>
void reduce_contributions_into(const CpuExecutor& exec, const std::vector<ContributionBuffer<Delta>>& thread_buffers,
                                std::size_t particle_count, std::vector<Delta>& output, ReduceScratch& scratch,
                                const ReduceOptions& options = ReduceOptions{})
{
    const std::vector<ContributionBuffer<Delta>>* const set = &thread_buffers;
    reduce_contributions_into(exec, &set, std::size_t(1), particle_count, output, scratch, options);
}

// Несколько наборов буферов (конкатенация в порядке перечисления).
template <class Delta>
void reduce_contributions_into(
    const CpuExecutor& exec, std::initializer_list<const std::vector<ContributionBuffer<Delta>>*> buffer_sets,
    std::size_t particle_count, std::vector<Delta>& output, ReduceScratch& scratch,
    const ReduceOptions& options = ReduceOptions{})
{
    reduce_contributions_into(exec, buffer_sets.begin(), buffer_sets.size(), particle_count, output, scratch,
                              options);
}

// Обёртка со старой сигнатурой: локальная рабочая память, результат по значению.
template <class Delta>
std::vector<Delta> reduce_contributions(const CpuExecutor& exec,
                                         const std::vector<ContributionBuffer<Delta>>& thread_buffers,
                                         std::size_t particle_count,
                                         const ReduceOptions& options = ReduceOptions{})
{
    std::vector<Delta> output;
    ReduceScratch scratch;
    reduce_contributions_into(exec, thread_buffers, particle_count, output, scratch, options);
    return output;
}

// Обёртка для одного индексированного массива: локальная рабочая память.
template <class Delta>
std::vector<Delta> reduce_contributions(const CpuExecutor& exec, const IndexedContributions<Delta>& contribs,
                                         std::size_t particle_count,
                                         const ReduceOptions& options = ReduceOptions{})
{
    std::vector<Delta> output;
    ReduceScratch scratch;
    reduce_contributions_into(exec, contribs, {}, particle_count, output, scratch, options);
    return output;
}

} // namespace rhs

#endif // RHS_CONTRIBUTIONS_HPP
