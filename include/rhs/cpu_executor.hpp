#ifndef RHS_CPU_EXECUTOR_HPP
#define RHS_CPU_EXECUTOR_HPP

#include <algorithm>
#include <cstddef>
#include <future>
#include <iterator>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// rhs/cpu_executor.hpp — CpuExecutor: единая абстракция параллельности CPU
// backend (fiber/SPH.md, раздел "CPU executor"; fiber/PLAN.md §3.5).
//
// Один объект на весь расчёт. Interaction kernels никогда не видят
// #pragma omp / mutex / atomic напрямую (SPH.md, "Основные требования к
// interaction kernel", п.3) — только эти примитивы:
//
//     parallel_for(n, grain, f)
//     parallel_reduce(n, init, map, combine)
//     parallel_scan(in, out, n)
//     parallel_sort(first, last, less)
//     task_group: run(f), wait()
//     num_workers(), worker_id()
//
// Реализация — на OpenMP (проект уже линкуется с -fopenmp; второй пул
// потоков конкурировал бы с ним за ядра, PLAN.md §3.5). Без -fopenmp на
// этапе компиляции (нет макроса _OPENMP) executor вырождается в
// последовательное выполнение на вызывающем потоке — так же, как
// include/pair_graph.hpp обходится без OpenMP (см. комментарий в
// CMakeLists.txt у цели test_pair_graph).
// ============================================================================

namespace rhs
{

class CpuExecutor
{
public:
    // num_threads == 0 -> взять текущее значение omp_get_max_threads()
    // (или 1, если OpenMP не подключён на этапе компиляции). Явное значение
    // задаётся через num_threads(N) в omp-директивах ниже, а не через
    // глобальный omp_set_num_threads() — так несколько CpuExecutor с разным
    // числом потоков не мешают друг другу (используется в test_contributions
    // для детерминизма при 1/2/7 потоках).
    explicit CpuExecutor(int num_threads = 0)
        : num_threads_(resolve_thread_count(num_threads))
    {
    }

    int num_workers() const { return num_threads_; }

    // Индекс текущего потока внутри активного parallel-региона; вне региона
    // (последовательный контекст) — 0.
    static int worker_id()
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

    // ------------------------------------------------------------------
    // parallel_for(n, grain, f) — f(std::size_t i) для каждого i из [0, n)
    // независимо; порядок и разбиение по потокам не определены (SPH.md:
    // "Kernel не содержит #pragma omp"; диапазон режется на blocks с
    // dynamic scheduling).
    // ------------------------------------------------------------------
    template <class Func>
    void parallel_for(std::size_t n, std::size_t grain, Func&& f) const
    {
        if (n == 0) return;
        const std::size_t chunk = std::max<std::size_t>(grain, 1);
#ifdef _OPENMP
        if (num_threads_ > 1)
        {
            const long long nn = static_cast<long long>(n);
            const long long lchunk = static_cast<long long>(chunk);
#pragma omp parallel for schedule(dynamic, lchunk) num_threads(num_threads_)
            for (long long i = 0; i < nn; ++i)
            {
                f(static_cast<std::size_t>(i));
            }
            return;
        }
#endif
        for (std::size_t i = 0; i < n; ++i) f(i);
    }

    template <class Func>
    void parallel_for(std::size_t n, Func&& f) const
    {
        parallel_for(n, default_grain(n), std::forward<Func>(f));
    }

    // ------------------------------------------------------------------
    // parallel_reduce(n, init, map, combine): T map(i), T combine(a, b).
    // Контракт (как в TBB/std parallel reduce): `init` — нейтральный элемент
    // combine (combine(init, x) == x), combine ассоциативна и, для
    // воспроизводимости порядка объединения частичных сумм по потокам,
    // ожидается коммутативной. Каждый worker сворачивает свою часть
    // диапазона в локальный аккумулятор, финальное объединение — по числу
    // потоков (SPH.md, "ODE solver"/"Reduce" — общий паттерн CPU backend).
    // Раздача итераций — schedule(dynamic, lchunk), lchunk = n/(8*nthreads)
    // (около 8 порций на поток): при порции 1 каждая итерация — атомарная
    // операция над общим счётчиком libgomp, и на лёгком map (единицы нс)
    // потоки упираются в этот счётчик (speedup/contact.md).
    // ------------------------------------------------------------------
    template <class T, class MapFunc, class CombineFunc>
    T parallel_reduce(std::size_t n, T init, MapFunc&& map, CombineFunc&& combine) const
    {
        if (n == 0) return init;
#ifdef _OPENMP
        if (num_threads_ > 1)
        {
            std::vector<T> partial(static_cast<std::size_t>(num_threads_), init);
            const long long nn = static_cast<long long>(n);
            const long long lchunk = static_cast<long long>(
                std::max<std::size_t>(1, n / (8 * static_cast<std::size_t>(num_threads_))));
#pragma omp parallel num_threads(num_threads_)
            {
                const int tid = omp_get_thread_num();
                T local = init;
#pragma omp for schedule(dynamic, lchunk)
                for (long long i = 0; i < nn; ++i)
                {
                    local = combine(local, map(static_cast<std::size_t>(i)));
                }
                partial[static_cast<std::size_t>(tid)] = local;
            }
            T result = init;
            for (auto& p : partial) result = combine(result, p);
            return result;
        }
#endif
        T result = init;
        for (std::size_t i = 0; i < n; ++i) result = combine(result, map(i));
        return result;
    }

    // ------------------------------------------------------------------
    // parallel_scan(in, out, n) — исключающая префиксная сумма:
    // out[0] = 0 (T{}), out[i+1] = out[i] + in[i], out[n] = полная сумма.
    // `out` должен иметь размер n+1. Двухпроходная блочная схема (локальные
    // суммы блоков -> последовательная свёртка по блокам -> параллельная
    // запись со сдвигом) — тот же паттерн "count -> prefix sum -> copy",
    // что описан в SPH.md для Search и построения primitives.
    // ------------------------------------------------------------------
    template <class T>
    void parallel_scan(const T* in, T* out, std::size_t n) const
    {
        if (n == 0)
        {
            out[0] = T{};
            return;
        }
#ifdef _OPENMP
        if (num_threads_ > 1 && n >= static_cast<std::size_t>(num_threads_) * 2)
        {
            const int nt = num_threads_;
            std::vector<std::size_t> bounds(static_cast<std::size_t>(nt) + 1);
            for (int t = 0; t <= nt; ++t)
            {
                bounds[static_cast<std::size_t>(t)] =
                    (n * static_cast<std::size_t>(t)) / static_cast<std::size_t>(nt);
            }
            std::vector<T> block_sum(static_cast<std::size_t>(nt), T{});
#pragma omp parallel for num_threads(nt)
            for (int t = 0; t < nt; ++t)
            {
                T s{};
                const std::size_t b = bounds[static_cast<std::size_t>(t)];
                const std::size_t e = bounds[static_cast<std::size_t>(t) + 1];
                for (std::size_t i = b; i < e; ++i) s = s + in[i];
                block_sum[static_cast<std::size_t>(t)] = s;
            }
            std::vector<T> block_offset(static_cast<std::size_t>(nt) + 1, T{});
            for (int t = 0; t < nt; ++t)
            {
                block_offset[static_cast<std::size_t>(t) + 1] =
                    block_offset[static_cast<std::size_t>(t)] + block_sum[static_cast<std::size_t>(t)];
            }
#pragma omp parallel for num_threads(nt)
            for (int t = 0; t < nt; ++t)
            {
                T running = block_offset[static_cast<std::size_t>(t)];
                const std::size_t b = bounds[static_cast<std::size_t>(t)];
                const std::size_t e = bounds[static_cast<std::size_t>(t) + 1];
                for (std::size_t i = b; i < e; ++i)
                {
                    out[i] = running;
                    running = running + in[i];
                }
            }
            out[n] = block_offset[static_cast<std::size_t>(nt)];
            return;
        }
#endif
        T running{};
        for (std::size_t i = 0; i < n; ++i)
        {
            out[i] = running;
            running = running + in[i];
        }
        out[n] = running;
    }

    // ------------------------------------------------------------------
    // parallel_sort(first, last, comp) — диапазон делится на num_workers()
    // смежных блоков, каждый сортируется параллельно (std::sort), затем
    // блоки последовательно сливаются std::inplace_merge. Простая и
    // заведомо корректная схема; не самая быстрая асимптотика слияния —
    // для объёма данных этого этапа (buffers пар/primitives) это не
    // критично, а корректность важнее.
    // ------------------------------------------------------------------
    template <class RandomIt, class Compare>
    void parallel_sort(RandomIt first, RandomIt last, Compare comp) const
    {
        const auto n = std::distance(first, last);
        if (n < 2) return;
#ifdef _OPENMP
        if (num_threads_ > 1 && n >= static_cast<decltype(n)>(num_threads_) * 64)
        {
            const int nt = num_threads_;
            std::vector<RandomIt> bounds(static_cast<std::size_t>(nt) + 1);
            for (int t = 0; t <= nt; ++t)
            {
                bounds[static_cast<std::size_t>(t)] = first + (n * t) / nt;
            }
#pragma omp parallel for num_threads(nt)
            for (int t = 0; t < nt; ++t)
            {
                std::sort(bounds[static_cast<std::size_t>(t)], bounds[static_cast<std::size_t>(t) + 1], comp);
            }
            for (int t = 1; t < nt; ++t)
            {
                std::inplace_merge(bounds[0], bounds[static_cast<std::size_t>(t)],
                                    bounds[static_cast<std::size_t>(t) + 1], comp);
            }
            return;
        }
#endif
        std::sort(first, last, comp);
    }

    // ------------------------------------------------------------------
    // task_group — запуск независимых ветвей графа RHS (SPH.md: "После
    // построения необходимых индексов branches различных interaction
    // types могут выполняться независимо"; PLAN.md §3.4 — Fiber/Solid
    // ветви и их reduce).
    //
    // Инженерное отклонение: реализовано через std::async, а не через
    // собственный persistent thread pool executor'а. Причина — граф RHS
    // (Fiber/Solid ветви, их reduce) появляется только в Ф2+ вместе с
    // самими compute-структурами; на Ф1 task_group не участвует ни в одном
    // из требуемых тестов (test_executor проверяет только
    // for/reduce/scan/sort). Интерфейс run()/wait() соответствует
    // спецификации; механизм будет пересмотрен вместе с подключением
    // реального RHS pipeline, если понадобится.
    // ------------------------------------------------------------------
    class TaskGroup
    {
    public:
        template <class Func>
        void run(Func&& f)
        {
            futures_.push_back(std::async(std::launch::async, std::forward<Func>(f)));
        }

        void wait()
        {
            for (auto& fut : futures_) fut.get();
            futures_.clear();
        }

        ~TaskGroup() { wait(); }

    private:
        std::vector<std::future<void>> futures_;
    };

    TaskGroup make_task_group() const { return TaskGroup{}; }

private:
    static std::size_t default_grain(std::size_t n)
    {
        return std::max<std::size_t>(std::size_t(1), n / 64);
    }

    static int resolve_thread_count(int requested)
    {
        if (requested > 0) return requested;
#ifdef _OPENMP
        return std::max(1, omp_get_max_threads());
#else
        return 1;
#endif
    }

    int num_threads_;
};

} // namespace rhs

#endif // RHS_CPU_EXECUTOR_HPP
