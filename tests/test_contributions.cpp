// tests/test_contributions.cpp
// Этап Ф1 (fiber/PLAN.md): unit-тест ContributionWriter/ContributionBuffer и
// reduce_contributions (include/rhs/contributions.hpp), а также
// RhsWorkspace (include/rhs/rhs_workspace.hpp) как контейнера буферов.
//
// Проверяет:
//  1) сумма вкладов, записанных через ContributionWriter и сведённых
//     reduce_contributions, совпадает с последовательным суммированием;
//  2) в детерминированном режиме результат ПОБИТНО одинаков при 1, 2 и 7
//     потоках (omp-число потоков задаётся через CpuExecutor(N), без
//     глобального omp_set_num_threads);
//  3) RhsWorkspace отдаёт один и тот же переиспользуемый буфер между
//     "вызовами RHS", capacity не уменьшается после clear();
//  4) reduce_contributions_into (без копий вкладов, с постоянной рабочей
//     памятью ReduceScratch) побитно совпадает с прежней реализацией
//     (reduce_contributions_reference, копия старого кода) при любом числе
//     потоков, числе буферов/наборов буферов, при повторных вызовах на одной
//     ReduceScratch с растущими и убывающими размерами; reduce_shard_of
//     совпадает с std::upper_bound.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>
#ifdef _WIN32
#define NOMINMAX // std::max/std::min из <rhs/...> конфликтуют с макросами min/max windows.h
#include <windows.h>
#endif

#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/rhs_workspace.hpp"
#include "rhs/types.hpp"

using rhs::Contribution;
using rhs::ContributionBuffer;
using rhs::ContributionWriter;
using rhs::CpuExecutor;
using rhs::index_t;
using rhs::reduce_contributions;
using rhs::reduce_contributions_into;
using rhs::reduce_shard_of;
using rhs::ReduceOptions;
using rhs::ReduceScratch;
using rhs::RhsWorkspace;

static int failures = 0;

static void check(bool cond, const char* msg)
{
    if (!cond)
    {
        std::printf("FAIL: %s\n", msg);
        ++failures;
    }
}

// Тестовая Delta — НЕ физическая величина (Ф1 не содержит Fiber/Solid),
// просто вектор из 3 double с определённым сложением, как того требует
// SPH.md ("операция сложения Delta определена на уровне compute-структуры").
struct TestDelta
{
    double v[3] = {0.0, 0.0, 0.0};
    TestDelta& operator+=(const TestDelta& o)
    {
        v[0] += o.v[0];
        v[1] += o.v[1];
        v[2] += o.v[2];
        return *this;
    }
};

static bool exact_equal(const std::vector<TestDelta>& a, const std::vector<TestDelta>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (int k = 0; k < 3; ++k)
            if (a[i].v[k] != b[i].v[k]) return false;
    return true;
}

static double max_abs_diff(const std::vector<TestDelta>& a, const std::vector<TestDelta>& b)
{
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (int k = 0; k < 3; ++k) m = std::max(m, std::fabs(a[i].v[k] - b[i].v[k]));
    return m;
}

// Один "primitive" (синтетическое взаимодействие) с фиксированным по номеру
// primitive содержимым — не зависит от того, какой поток его обработал, так
// что итоговое МНОЖЕСТВО contributions одно и то же при любом числе потоков.
struct SyntheticPrimitive
{
    index_t targets[2];
    TestDelta delta_a, delta_b; // на targets[0] и targets[1] соответственно
};

static std::vector<SyntheticPrimitive> make_primitives(std::size_t count, std::size_t particle_count)
{
    std::vector<SyntheticPrimitive> prims(count);
    std::mt19937 rng(777);
    std::uniform_int_distribution<std::size_t> pid(0, particle_count - 1);
    std::uniform_real_distribution<double> val(-10.0, 10.0);
    for (std::size_t k = 0; k < count; ++k)
    {
        prims[k].targets[0] = pid(rng);
        prims[k].targets[1] = pid(rng);
        for (int c = 0; c < 3; ++c)
        {
            prims[k].delta_a.v[c] = val(rng);
            prims[k].delta_b.v[c] = -prims[k].delta_a.v[c]; // как в Fiber-Bond: F_i = -F_j
        }
    }
    return prims;
}

// Прогоняет primitives через CpuExecutor::parallel_for с nthreads worker'ов,
// каждый worker пишет в собственный ContributionBuffer через ContributionWriter,
// затем reduce_contributions сводит их в итоговый вектор.
static std::vector<TestDelta> run_reduce(const std::vector<SyntheticPrimitive>& prims, std::size_t particle_count,
                                          int nthreads, bool deterministic)
{
    CpuExecutor exec(nthreads);
    std::vector<ContributionBuffer<TestDelta>> buffers(static_cast<std::size_t>(exec.num_workers()));

    const std::uint32_t interaction_type = 42; // произвольная фиксированная метка типа взаимодействия

    exec.parallel_for(prims.size(), [&](std::size_t k) {
        int w = CpuExecutor::worker_id();
        if (w < 0 || w >= exec.num_workers()) w = 0;
        ContributionWriter<TestDelta> writer(buffers[static_cast<std::size_t>(w)], interaction_type);
        const auto& p = prims[k];
        writer.emit(p.targets[0], p.delta_a, /*primitive_id=*/k, /*slot=*/0);
        writer.emit(p.targets[1], p.delta_b, /*primitive_id=*/k, /*slot=*/1);
    });

    ReduceOptions opts;
    opts.deterministic = deterministic;
    return reduce_contributions(exec, buffers, particle_count, opts);
}

static std::vector<TestDelta> sequential_reference(const std::vector<SyntheticPrimitive>& prims,
                                                     std::size_t particle_count)
{
    std::vector<TestDelta> out(particle_count, TestDelta{});
    for (const auto& p : prims)
    {
        out[p.targets[0]] += p.delta_a;
        out[p.targets[1]] += p.delta_b;
    }
    return out;
}

// -----------------------------------------------------------------------
// 1) Сумма вкладов (обычный режим) совпадает с последовательным суммированием
//    с точностью до порядка сложения (плавающая точка) — проверяем допуском.
// -----------------------------------------------------------------------
static void test_sum_matches_sequential()
{
    const std::size_t particle_count = 500;
    auto prims = make_primitives(3000, particle_count);
    auto reference = sequential_reference(prims, particle_count);

    for (int nt : {1, 2, 7})
    {
        auto result = run_reduce(prims, particle_count, nt, /*deterministic=*/false);
        double diff = max_abs_diff(result, reference);
        check(diff < 1e-9, "reduce_contributions (обычный режим): сумма расходится с последовательным эталоном");
    }

    // Пустой набор primitives.
    {
        std::vector<SyntheticPrimitive> empty;
        auto result = run_reduce(empty, particle_count, 4, false);
        auto reference_empty = sequential_reference(empty, particle_count);
        check(exact_equal(result, reference_empty), "reduce_contributions: пустой набор primitives должен давать нулевой RHS");
    }

    // particle_count == 0.
    {
        std::vector<SyntheticPrimitive> empty;
        auto result = run_reduce(empty, 0, 4, false);
        check(result.empty(), "reduce_contributions: particle_count == 0 должен давать пустой результат");
    }
}

// -----------------------------------------------------------------------
// 2) Детерминированный режим: результат побитно одинаков при 1, 2 и 7
//    потоках.
// -----------------------------------------------------------------------
static void test_deterministic_bitexact()
{
    const std::size_t particle_count = 733; // не кратно ни 1*4, ни 2*4, ни 7*4 -- разные shard_count
    auto prims = make_primitives(5000, particle_count);

    auto result_1 = run_reduce(prims, particle_count, 1, /*deterministic=*/true);
    auto result_2 = run_reduce(prims, particle_count, 2, /*deterministic=*/true);
    auto result_7 = run_reduce(prims, particle_count, 7, /*deterministic=*/true);

    check(exact_equal(result_1, result_2),
          "детерминированный режим: результат при 1 и 2 потоках не совпадает побитно");
    check(exact_equal(result_1, result_7),
          "детерминированный режим: результат при 1 и 7 потоках не совпадает побитно");

    // Дополнительно: сумма (по модулю порядка) совпадает с последовательным эталоном.
    auto reference = sequential_reference(prims, particle_count);
    check(max_abs_diff(result_1, reference) < 1e-9,
          "детерминированный режим: сумма расходится с последовательным эталоном");

    // Повторный прогон с теми же данными и тем же числом потоков даёт тот же
    // результат (воспроизводимость одного и того же запуска, не только
    // между разными N потоков).
    auto result_2_again = run_reduce(prims, particle_count, 2, true);
    check(exact_equal(result_2, result_2_again), "детерминированный режим: повторный прогон при 2 потоках дал другой результат");
}

// Копия прежней реализации reduce_contributions (до перехода на ссылки и
// ReduceScratch) — эталон для побитового сравнения.
template <class Delta>
static std::vector<Delta> reduce_contributions_reference(const CpuExecutor& exec,
                                                          const std::vector<ContributionBuffer<Delta>>& thread_buffers,
                                                          std::size_t particle_count, const ReduceOptions& options)
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

using BufferSet = std::vector<ContributionBuffer<TestDelta>>;

// Побитовое сравнение (memcmp: различает и +0/-0), в отличие от exact_equal.
static bool bitwise_equal(const std::vector<TestDelta>& a, const std::vector<TestDelta>& b)
{
    if (a.size() != b.size()) return false;
    return a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(TestDelta)) == 0;
}

// Набор из nbuf буферов со случайным (в том числе нулевым) числом вкладов в
// каждом. primitive_id уникален в пределах next_id, так что ключи
// детерминированной сортировки (target, type, primitive_id, slot) уникальны.
static BufferSet make_buffers(std::size_t nbuf, std::size_t max_per_buffer, std::size_t particle_count,
                               std::mt19937& rng, std::uint64_t& next_id)
{
    BufferSet set(nbuf);
    std::uniform_int_distribution<std::size_t> cnt(0, max_per_buffer);
    std::uniform_int_distribution<std::size_t> pid(0, particle_count - 1);
    std::uniform_real_distribution<double> val(-10.0, 10.0);
    for (auto& buf : set)
    {
        const std::size_t n = cnt(rng);
        for (std::size_t k = 0; k < n; ++k)
        {
            TestDelta d;
            for (int c = 0; c < 3; ++c) d.v[c] = val(rng);
            buf.push(pid(rng), d, rhs::ContributionMeta{static_cast<std::uint32_t>(next_id % 3), next_id, 0});
            ++next_id;
        }
    }
    return set;
}

static BufferSet concat_sets(std::initializer_list<const BufferSet*> sets)
{
    BufferSet all;
    for (const BufferSet* set : sets) all.insert(all.end(), set->begin(), set->end());
    return all;
}

// -----------------------------------------------------------------------
// 4a) reduce_shard_of == std::upper_bound со старой семантикой, включая
//     target >= particle_count (последний shard).
// -----------------------------------------------------------------------
static void test_shard_of_matches_upper_bound()
{
    bool ok = true;
    for (std::size_t n : {1u, 2u, 3u, 7u, 100u, 733u, 4096u, 100003u})
    {
        for (std::size_t S : {1u, 2u, 3u, 4u, 28u, 64u, 1000u, 200000u})
        {
            const std::size_t shards = std::max<std::size_t>(1, std::min(S, n));
            std::vector<std::size_t> bound(shards + 1);
            for (std::size_t s = 0; s <= shards; ++s) bound[s] = (n * s) / shards;
            const double inv = static_cast<double>(shards) / static_cast<double>(n);
            for (std::size_t t = 0; t < n + 5; ++t)
            {
                auto it = std::upper_bound(bound.begin(), bound.end(), t);
                std::size_t s = static_cast<std::size_t>(it - bound.begin());
                s = (s == 0) ? 0 : std::min(s - 1, shards - 1);
                if (reduce_shard_of(bound.data(), shards, n, inv, t) != s) ok = false;
            }
        }
    }
    check(ok, "reduce_shard_of расходится с std::upper_bound");
}

// -----------------------------------------------------------------------
// 4b) reduce_contributions_into побитно == прежняя реализация: 1/2/7 потоков,
//     разное число буферов, shard_count, обычный и детерминированный режимы.
// -----------------------------------------------------------------------
static void test_reduce_into_matches_reference()
{
    std::mt19937 rng(12345);
    std::uint64_t next_id = 0;
    bool ok = true, ok_wrapper = true;
    for (int nt : {1, 2, 7})
    {
        CpuExecutor exec(nt);
        for (std::size_t particle_count : {1u, 2u, 5u, 733u, 5000u})
        {
            for (std::size_t nbuf : {0u, 1u, 3u, 16u, 40u})
            {
                BufferSet buffers = make_buffers(nbuf, 400, particle_count, rng, next_id);
                for (std::size_t shard_count : {0u, 1u, 3u, 100000u})
                {
                    for (bool det : {false, true})
                    {
                        ReduceOptions opts;
                        opts.shard_count = shard_count;
                        opts.deterministic = det;
                        auto expected = reduce_contributions_reference(exec, buffers, particle_count, opts);
                        ReduceScratch scratch;
                        std::vector<TestDelta> actual;
                        reduce_contributions_into(exec, buffers, particle_count, actual, scratch, opts);
                        if (!bitwise_equal(actual, expected)) ok = false;
                        if (!bitwise_equal(reduce_contributions(exec, buffers, particle_count, opts), expected))
                            ok_wrapper = false;
                    }
                }
            }
        }
    }
    check(ok, "reduce_contributions_into не совпадает побитно с прежней реализацией");
    check(ok_wrapper, "reduce_contributions (обёртка) не совпадает побитно с прежней реализацией");

    // Граничные случаи: total == 0, particle_count == 0, нет буферов.
    {
        CpuExecutor exec(4);
        ReduceScratch scratch;
        std::vector<TestDelta> out(10, TestDelta{{1.0, 2.0, 3.0}});
        BufferSet none;
        reduce_contributions_into(exec, none, 8, out, scratch);
        check(out.size() == 8 && bitwise_equal(out, std::vector<TestDelta>(8)),
              "reduce_contributions_into: без буферов должен быть нулевой вектор нужного размера");
        BufferSet empty_bufs(5);
        reduce_contributions_into(exec, empty_bufs, 8, out, scratch);
        check(out.size() == 8 && bitwise_equal(out, std::vector<TestDelta>(8)),
              "reduce_contributions_into: total == 0 должен давать нулевой вектор");
        BufferSet some = make_buffers(3, 50, 8, rng, next_id);
        reduce_contributions_into(exec, some, 8, out, scratch);
        check(out.size() == 8, "reduce_contributions_into: размер результата");
        reduce_contributions_into(exec, some, 0, out, scratch);
        check(out.empty(), "reduce_contributions_into: particle_count == 0 должен давать пустой результат");
    }
}

// -----------------------------------------------------------------------
// 4c) Несколько наборов буферов == прежняя реализация на их конкатенации.
// -----------------------------------------------------------------------
static void test_reduce_into_multiple_sets()
{
    std::mt19937 rng(54321);
    std::uint64_t next_id = 0;
    bool ok = true, ok_single = true;
    for (int nt : {1, 2, 7})
    {
        CpuExecutor exec(nt);
        for (std::size_t particle_count : {1u, 37u, 1500u})
        {
            BufferSet a = make_buffers(static_cast<std::size_t>(nt), 300, particle_count, rng, next_id);
            BufferSet b = make_buffers(static_cast<std::size_t>(nt) + 2, 500, particle_count, rng, next_id);
            BufferSet empty_set;
            BufferSet c = make_buffers(3, 200, particle_count, rng, next_id);
            const BufferSet all = concat_sets({&a, &b, &empty_set, &c});
            for (bool det : {false, true})
            {
                ReduceOptions opts;
                opts.deterministic = det;
                auto expected = reduce_contributions_reference(exec, all, particle_count, opts);
                ReduceScratch scratch;
                std::vector<TestDelta> actual;
                reduce_contributions_into(exec, {&a, &b, &empty_set, &c}, particle_count, actual, scratch, opts);
                if (!bitwise_equal(actual, expected)) ok = false;

                // Порядок наборов важен: {c, a} совпадает с эталоном на concat(c, a).
                auto expected_ca = reduce_contributions_reference(exec, concat_sets({&c, &a}), particle_count, opts);
                reduce_contributions_into(exec, {&c, &a}, particle_count, actual, scratch, opts);
                if (!bitwise_equal(actual, expected_ca)) ok = false;

                // Один набор через initializer_list.
                auto expected_a = reduce_contributions_reference(exec, a, particle_count, opts);
                reduce_contributions_into(exec, {&a}, particle_count, actual, scratch, opts);
                if (!bitwise_equal(actual, expected_a)) ok_single = false;
            }
        }
    }
    check(ok, "reduce_contributions_into (несколько наборов) не совпадает побитно с эталоном на конкатенации");
    check(ok_single, "reduce_contributions_into (initializer_list из одного набора) не совпадает с эталоном");
}

// -----------------------------------------------------------------------
// 4d) Повторные вызовы на одной ReduceScratch и одном output: размеры
//     растут и убывают, результат каждый раз побитно равен эталону; память
//     scratch не уменьшается.
// -----------------------------------------------------------------------
static void test_reduce_into_scratch_reuse()
{
    std::mt19937 rng(99);
    std::uint64_t next_id = 0;
    CpuExecutor exec(3);
    ReduceScratch scratch;
    std::vector<TestDelta> output;
    bool ok = true, mem_ok = true;
    std::size_t prev_refs_cap = 0, prev_counts_cap = 0, prev_out_cap = 0;
    // Размеры: рост, спад, снова рост, смена particle_count и числа буферов.
    const std::size_t particles[] = {50, 2000, 300, 7, 2000, 1, 900, 2000};
    const std::size_t nbufs[] = {3, 12, 5, 1, 12, 2, 20, 12};
    const std::size_t maxper[] = {20, 600, 100, 5, 700, 3, 250, 600};
    for (std::size_t it = 0; it < sizeof(particles) / sizeof(particles[0]); ++it)
    {
        BufferSet x = make_buffers(nbufs[it], maxper[it], particles[it], rng, next_id);
        BufferSet y = make_buffers(2, maxper[it] / 2, particles[it], rng, next_id);
        for (bool det : {false, true})
        {
            ReduceOptions opts;
            opts.deterministic = det;
            auto expected = reduce_contributions_reference(exec, concat_sets({&x, &y}), particles[it], opts);
            reduce_contributions_into(exec, {&x, &y}, particles[it], output, scratch, opts);
            if (!bitwise_equal(output, expected)) ok = false;
        }
        if (scratch.refs.capacity() < prev_refs_cap || scratch.counts.capacity() < prev_counts_cap ||
            output.capacity() < prev_out_cap)
            mem_ok = false;
        prev_refs_cap = scratch.refs.capacity();
        prev_counts_cap = scratch.counts.capacity();
        prev_out_cap = output.capacity();
    }
    check(ok, "reduce_contributions_into: повторные вызовы на одной ReduceScratch дали результат, отличный от эталона");
    check(mem_ok, "ReduceScratch/output: capacity не должна уменьшаться (память только растёт)");

    // После самого большого вызова меньший не перевыделяет память.
    const void* const* refs_data = scratch.refs.data();
    const TestDelta* out_data = output.data();
    BufferSet shrunk = make_buffers(2, 10, 20, rng, next_id);
    reduce_contributions_into(exec, shrunk, 20, output, scratch);
    check(scratch.refs.data() == refs_data && output.data() == out_data,
          "reduce_contributions_into: меньший вызов после большего не должен перевыделять память");
}

// -----------------------------------------------------------------------
// 4e) IndexedContributions + наборы буферов == прежняя реализация на одном
//     буфере (1 поток), где лежит вся последовательность "индексированные,
//     затем буферы контакта"; при 1/2/7 потоках, оба режима. Плюс prepare:
//     ёмкость только растёт, уменьшение не перевыделяет память.
// -----------------------------------------------------------------------
static void fill_indexed(rhs::IndexedContributions<TestDelta>& idx, std::size_t n, std::size_t particle_count,
                          std::mt19937& rng, std::uint64_t& next_id)
{
    idx.prepare(n);
    std::uniform_int_distribution<std::size_t> pid(0, particle_count - 1);
    std::uniform_real_distribution<double> val(-10.0, 10.0);
    for (std::size_t k = 0; k < n; ++k)
    {
        TestDelta d;
        for (int c = 0; c < 3; ++c) d.v[c] = val(rng);
        idx[k] = rhs::Contribution<TestDelta>{pid(rng), d,
                                              rhs::ContributionMeta{static_cast<std::uint32_t>(next_id % 3), next_id, 0}};
        ++next_id;
    }
}

static void test_indexed_contributions()
{
    std::mt19937 rng(777);
    std::uint64_t next_id = 0;
    bool ok = true;
    CpuExecutor exec1(1);
    for (int nt : {1, 2, 7})
    {
        CpuExecutor exec(nt);
        rhs::IndexedContributions<TestDelta> idx;
        ReduceScratch scratch;
        std::vector<TestDelta> actual;
        // Рост, затем спад и снова рост размера на одном массиве и одной scratch.
        for (std::size_t n : {0u, 1u, 5u, 1000u, 37u, 0u, 2500u})
        {
            for (std::size_t particle_count : {1u, 40u, 1200u})
            {
                fill_indexed(idx, n, particle_count, rng, next_id);
                BufferSet contact = make_buffers(static_cast<std::size_t>(nt), 200, particle_count, rng, next_id);
                BufferSet one(1);
                for (std::size_t k = 0; k < idx.size(); ++k) one[0].push(idx[k].target, idx[k].delta, idx[k].meta);
                const BufferSet all = concat_sets({&one, &contact});
                for (bool det : {false, true})
                {
                    ReduceOptions opts;
                    opts.deterministic = det;
                    // Эталон: прежняя схема, один поток, вся последовательность в одном буфере.
                    BufferSet merged(1);
                    for (const auto& buf : all)
                        for (std::size_t k = 0; k < buf.size(); ++k)
                            merged[0].push(buf[k].target, buf[k].delta, buf[k].meta);
                    auto expected = reduce_contributions_reference(exec1, merged, particle_count, opts);
                    reduce_contributions_into(exec, idx, {&contact}, particle_count, actual, scratch, opts);
                    // Побитово при любом числе потоков: порядок слагаемых на частицу
                    // совпадает с эталоном (first по индексам, затем буферы контакта по порядку).
                    if (!bitwise_equal(actual, expected)) ok = false;
                    if (!bitwise_equal(reduce_contributions(exec, idx, particle_count, opts),
                                       reduce_contributions_reference(exec1, one, particle_count, opts)))
                        ok = false;
                }
            }
        }
    }
    check(ok, "IndexedContributions + буферы: результат не совпадает побитно с прежней схемой на одном буфере");

    // prepare: ёмкость только растёт, уменьшение не перевыделяет и не затирает.
    rhs::IndexedContributions<TestDelta> idx;
    idx.prepare(100);
    const auto* p = idx.data();
    const std::size_t cap = idx.capacity();
    idx.prepare(10);
    check(idx.size() == 10 && idx.data() == p && idx.capacity() == cap,
          "IndexedContributions::prepare: уменьшение размера не должно перевыделять память");
    idx.prepare(100);
    check(idx.size() == 100 && idx.data() == p && idx.capacity() == cap,
          "IndexedContributions::prepare: возврат к прежнему размеру не должен перевыделять память");
    idx.prepare(500);
    check(idx.size() == 500 && idx.capacity() >= 500, "IndexedContributions::prepare: рост размера");
}

// -----------------------------------------------------------------------
// 3) RhsWorkspace: переиспользуемый буфер между "вызовами RHS".
// -----------------------------------------------------------------------
static void test_rhs_workspace_reuse()
{
    RhsWorkspace ws;
    using BufVec = std::vector<ContributionBuffer<TestDelta>>;

    BufVec& buf1 = ws.get<BufVec>("fiber.contrib");
    buf1.resize(4);
    for (auto& b : buf1)
        for (int i = 0; i < 100; ++i) b.push(static_cast<index_t>(i), TestDelta{}, rhs::ContributionMeta{});
    std::size_t cap_after_fill = buf1[0].capacity();
    check(cap_after_fill >= 100, "RhsWorkspace: буфер должен вместить записанные вклады");

    // Имитация конца одного RHS: kernel очищает буферы, не освобождая память.
    for (auto& b : buf1) b.clear();

    BufVec& buf2 = ws.get<BufVec>("fiber.contrib");
    check(&buf1 == &buf2, "RhsWorkspace: повторный get<T> с тем же ключом должен вернуть тот же объект");
    check(buf2[0].size() == 0, "RhsWorkspace: clear() должен сбрасывать размер буфера");
    check(buf2[0].capacity() == cap_after_fill,
          "RhsWorkspace: clear() не должен уменьшать capacity (память только растёт, SPH.md)");

    // Другой ключ / другой тип — независимый буфер.
    int& counter = ws.get<int>("some.counter");
    check(counter == 0, "RhsWorkspace: новый буфер должен быть проинициализирован значением по умолчанию");
    counter = 5;
    check(ws.get<int>("some.counter") == 5, "RhsWorkspace: значение должно сохраняться между обращениями");
}

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    std::printf("=== test_contributions: ContributionWriter/reduce_contributions (fiber/PLAN.md, этап Ф1) ===\n");

    test_sum_matches_sequential();
    test_deterministic_bitexact();
    test_rhs_workspace_reuse();
    test_shard_of_matches_upper_bound();
    test_reduce_into_matches_reference();
    test_reduce_into_multiple_sets();
    test_reduce_into_scratch_reuse();
    test_indexed_contributions();

    if (failures == 0)
        std::printf("test_contributions: ПРОЙДЕН\n");
    else
        std::printf("test_contributions: ПРОВАЛЕН (%d ошибок)\n", failures);

    return failures == 0 ? 0 : 1;
}
