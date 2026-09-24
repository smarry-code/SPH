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
//     "вызовами RHS", capacity не уменьшается после clear().
#include <cmath>
#include <cstdio>
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
using rhs::ReduceOptions;
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

    if (failures == 0)
        std::printf("test_contributions: ПРОЙДЕН\n");
    else
        std::printf("test_contributions: ПРОВАЛЕН (%d ошибок)\n", failures);

    return failures == 0 ? 0 : 1;
}
