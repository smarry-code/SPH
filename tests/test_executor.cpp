// tests/test_executor.cpp
// Этап Ф1 (fiber/PLAN.md): unit-тест CpuExecutor (include/rhs/cpu_executor.hpp).
// Проверяет parallel_for/parallel_reduce/parallel_scan/parallel_sort против
// последовательных эталонов на разных N, включая пустой диапазон и один
// элемент, при 1, 2 и 7 потоках.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <random>
#include <vector>
#ifdef _WIN32
#define NOMINMAX // std::max/std::min из <rhs/...> конфликтуют с макросами min/max windows.h
#include <windows.h>
#endif

#include "rhs/cpu_executor.hpp"

using rhs::CpuExecutor;

static int failures = 0;

static void check(bool cond, const char* msg)
{
    if (!cond)
    {
        std::printf("FAIL: %s\n", msg);
        ++failures;
    }
}

// -----------------------------------------------------------------------
// parallel_for: заполнение массива f(i) = i*i сравнивается с последовательным
// эталоном для нескольких N (включая 0 и 1) и чисел потоков.
// -----------------------------------------------------------------------
static void test_parallel_for()
{
    const std::vector<std::size_t> sizes = {0, 1, 2, 7, 1000, 100000};
    const std::vector<int> thread_counts = {1, 2, 7};

    for (int nt : thread_counts)
    {
        CpuExecutor exec(nt);
        check(exec.num_workers() == nt, "num_workers() возвращает запрошенное число потоков");

        for (std::size_t n : sizes)
        {
            std::vector<long long> actual(n, -1);
            exec.parallel_for(n, [&](std::size_t i) { actual[i] = static_cast<long long>(i) * static_cast<long long>(i); });

            std::vector<long long> expected(n);
            for (std::size_t i = 0; i < n; ++i) expected[i] = static_cast<long long>(i) * static_cast<long long>(i);

            check(actual == expected, "parallel_for: результат не совпадает с последовательным эталоном");
        }
    }

    // Явный grain, отличный от default.
    {
        CpuExecutor exec(4);
        const std::size_t n = 5000;
        std::vector<int> actual(n, 0);
        exec.parallel_for(n, std::size_t(3), [&](std::size_t i) { actual[i] = 1; });
        check(std::all_of(actual.begin(), actual.end(), [](int v) { return v == 1; }),
              "parallel_for с явным grain: каждый элемент обработан ровно один раз");
    }
}

// -----------------------------------------------------------------------
// parallel_reduce: сумма квадратов как int64 (точная, ассоциативная и
// коммутативная операция — сравнение побитно корректно независимо от
// порядка объединения по потокам).
// -----------------------------------------------------------------------
static void test_parallel_reduce()
{
    const std::vector<std::size_t> sizes = {0, 1, 2, 7, 1000, 50000};
    const std::vector<int> thread_counts = {1, 2, 7};

    for (int nt : thread_counts)
    {
        CpuExecutor exec(nt);
        for (std::size_t n : sizes)
        {
            long long expected = 0;
            for (std::size_t i = 0; i < n; ++i)
                expected += static_cast<long long>(i) * static_cast<long long>(i);

            long long actual = exec.parallel_reduce(
                n, static_cast<long long>(0),
                [](std::size_t i) { return static_cast<long long>(i) * static_cast<long long>(i); },
                [](long long a, long long b) { return a + b; });

            check(actual == expected, "parallel_reduce: сумма квадратов не совпадает с последовательным эталоном");
        }
    }

    // reduce с init != 0 (нейтральный элемент combine): max с init = минимально возможным.
    {
        CpuExecutor exec(7);
        std::vector<int> data = {3, -1, 42, 17, 0, -100, 8};
        int expected = *std::max_element(data.begin(), data.end());
        int actual = exec.parallel_reduce(
            data.size(), std::numeric_limits<int>::min(), [&](std::size_t i) { return data[i]; },
            [](int a, int b) { return std::max(a, b); });
        check(actual == expected, "parallel_reduce: max по данным не совпадает с эталоном");
    }
}

// -----------------------------------------------------------------------
// parallel_scan: исключающая префиксная сумма против последовательного
// накопления. out размера n+1, out[n] == полная сумма.
// -----------------------------------------------------------------------
static void test_parallel_scan()
{
    const std::vector<std::size_t> sizes = {0, 1, 2, 7, 1000, 33333};
    const std::vector<int> thread_counts = {1, 2, 7};

    for (int nt : thread_counts)
    {
        CpuExecutor exec(nt);
        for (std::size_t n : sizes)
        {
            std::vector<long long> in(n);
            for (std::size_t i = 0; i < n; ++i) in[i] = static_cast<long long>(i % 13) + 1;

            std::vector<long long> expected(n + 1, 0);
            for (std::size_t i = 0; i < n; ++i) expected[i + 1] = expected[i] + in[i];

            std::vector<long long> actual(n + 1, -1);
            exec.parallel_scan(in.data(), actual.data(), n);

            check(actual == expected, "parallel_scan: результат не совпадает с последовательным эталоном");
        }
    }
}

// -----------------------------------------------------------------------
// parallel_sort: сортировка случайного вектора int (с повторами) против
// std::sort. Проверяется как порядок, так и сохранность мультимножества
// элементов.
// -----------------------------------------------------------------------
static void test_parallel_sort()
{
    const std::vector<std::size_t> sizes = {0, 1, 2, 7, 1000, 20000};
    const std::vector<int> thread_counts = {1, 2, 7};
    std::mt19937 rng(12345);

    for (int nt : thread_counts)
    {
        CpuExecutor exec(nt);
        for (std::size_t n : sizes)
        {
            std::uniform_int_distribution<int> dist(-1000, 1000);
            std::vector<int> data(n);
            for (auto& v : data) v = dist(rng);

            std::vector<int> expected = data;
            std::sort(expected.begin(), expected.end());

            std::vector<int> actual = data;
            exec.parallel_sort(actual.begin(), actual.end(), std::less<int>());

            check(actual == expected, "parallel_sort: результат не совпадает с std::sort");
        }
    }

    // Сортировка по убыванию (проверка произвольного компаратора).
    {
        CpuExecutor exec(4);
        std::vector<int> data = {5, 3, 8, 1, 9, 2, 7, 4, 6, 0};
        std::vector<int> expected = data;
        std::sort(expected.begin(), expected.end(), std::greater<int>());
        exec.parallel_sort(data.begin(), data.end(), std::greater<int>());
        check(data == expected, "parallel_sort: компаратор std::greater не соблюдён");
    }
}

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    std::printf("=== test_executor: CpuExecutor (fiber/PLAN.md, этап Ф1) ===\n");

    test_parallel_for();
    test_parallel_reduce();
    test_parallel_scan();
    test_parallel_sort();

    if (failures == 0)
        std::printf("test_executor: ПРОЙДЕН\n");
    else
        std::printf("test_executor: ПРОВАЛЕН (%d ошибок)\n", failures);

    return failures == 0 ? 0 : 1;
}
