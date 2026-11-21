// tests/test_executor.cpp
// Этап Ф1 (notes/fiber-model/fiber-plan.md): unit-тест CpuExecutor (include/rhs/cpu_executor.hpp).
// Проверяет parallel_for/parallel_reduce/parallel_scan/parallel_sort против
// последовательных эталонов на разных N, включая пустой диапазон и один
// элемент, при 1, 2 и 7 потоках.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <random>
#include <cstring>
#include <cmath>
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


// ---------------------------------------------------------------------------
// Подмножество с прежним разбиением на отрезки (parallel_for_seg_groups / parallel_reduce_seg_groups):
// свёртка по выбранным индексам побитово равна parallel_reduce по полному диапазону с нейтральным элементом
// у невыбранных; обход -- каждый выбранный ровно один раз, в своём отрезке, по возрастанию. concat_segments --
// порядок как у последовательного цикла.
// ---------------------------------------------------------------------------
static void test_seg_groups()
{
    std::mt19937 rng(77u);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    std::uniform_real_distribution<double> mag(-30.0, 30.0);
    bool ok_sum = true, ok_max = true, ok_cnt = true, ok_for = true, ok_cat = true;
    for (int nt : {1, 3, 8})
    {
        rhs::CpuExecutor exec(nt);
        const std::size_t S = exec.det_segments();
        for (std::size_t n : {std::size_t(0), std::size_t(1), std::size_t(5), std::size_t(63), std::size_t(64),
                              std::size_t(65), std::size_t(1000), std::size_t(12345)})
        {
            for (double frac : {0.0, 0.01, 0.3, 1.0})
            {
                std::vector<double> x(n);
                for (auto& v : x) v = val(rng) * std::pow(10.0, mag(rng) / 3.0); // разные порядки величин и знаки
                std::vector<std::size_t> sel;
                std::vector<std::uint16_t> seg;
                for (std::size_t s = 0; s < S; ++s)
                    for (std::size_t i = n * s / S; i < n * (s + 1) / S; ++i)
                        if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < frac)
                        {
                            sel.push_back(i);
                            seg.push_back(static_cast<std::uint16_t>(s));
                        }
                std::vector<char> is_sel(n, 0);
                for (std::size_t i : sel) is_sel[i] = 1;
                const auto plus = [](double a, double b) { return a + b; };
                const double full = exec.parallel_reduce(
                    n, 0.0, [&](std::size_t i) { return is_sel[i] ? x[i] : 0.0; }, plus);
                const double sub = exec.parallel_reduce_seg_groups(
                    sel.size(), [&](std::size_t j) { return seg[j]; }, 0.0, [&](std::size_t j) { return x[sel[j]]; },
                    plus);
                ok_sum = ok_sum && (std::memcmp(&full, &sub, sizeof(double)) == 0);
                const auto mx = [](double a, double b) { return std::max(a, b); };
                const double fm = exec.parallel_reduce(
                    n, 0.0, [&](std::size_t i) { return is_sel[i] ? std::fabs(x[i]) : 0.0; }, mx);
                const double sm = exec.parallel_reduce_seg_groups(
                    sel.size(), [&](std::size_t j) { return seg[j]; }, 0.0,
                    [&](std::size_t j) { return std::fabs(x[sel[j]]); }, mx);
                ok_max = ok_max && (fm == sm);
                const std::size_t c = exec.parallel_reduce_seg_groups(
                    sel.size(), [&](std::size_t j) { return seg[j]; }, std::size_t(0),
                    [&](std::size_t) { return std::size_t(1); }, [](std::size_t a, std::size_t b) { return a + b; });
                ok_cnt = ok_cnt && (c == sel.size());
                // обход: отрезок верный, каждый j один раз, внутри отрезка -- по возрастанию
                std::vector<std::vector<std::size_t>> seen(S);
                exec.parallel_for_seg_groups(sel.size(), [&](std::size_t j) { return seg[j]; },
                                             [&](std::size_t s, std::size_t j) { seen[s].push_back(j); });
                std::size_t tot = 0;
                for (std::size_t s = 0; s < S; ++s)
                {
                    tot += seen[s].size();
                    for (std::size_t q = 0; q < seen[s].size(); ++q)
                        ok_for = ok_for && seg[seen[s][q]] == s && (q == 0 || seen[s][q] > seen[s][q - 1]);
                }
                ok_for = ok_for && tot == sel.size();
                // concat_segments: порядок последовательного цикла
                std::vector<rhs::CpuExecutor::SegPart<std::size_t>> parts;
                exec.prepare_segment_parts(parts);
                exec.parallel_for_segments(n, [&](std::size_t s, std::size_t i) {
                    if (is_sel[i]) parts[s].v.push_back(i);
                });
                std::vector<std::size_t> cat(7, 999); // непустой out: прежнее содержимое перезаписывается
                exec.concat_segments(parts, cat);
                ok_cat = ok_cat && (cat == sel);
            }
        }
    }
    check(ok_sum, "parallel_reduce_seg_groups: сумма по подмножеству побитово = parallel_reduce с нулями");
    check(ok_max, "parallel_reduce_seg_groups: max по подмножеству = parallel_reduce с нулями");
    check(ok_cnt, "parallel_reduce_seg_groups: счётчик");
    check(ok_for, "parallel_for_seg_groups: каждый элемент один раз, в своём отрезке, по возрастанию");
    check(ok_cat, "concat_segments: порядок последовательного цикла");
}

// ---------------------------------------------------------------------------
// seg_of_index(k, n) -- отрезок индекса k при parallel_for_segments(n): совпадает с отрезком, который
// реально передаёт parallel_for_segments (все k для n = 1..3000, выборка k при n ~ 2^40 по границам
// отрезков). parallel_for_seg_bounds / parallel_reduce_seg_bounds при границах из seg_group_bounds --
// то же, что parallel_for_seg_groups / parallel_reduce_seg_groups (побитово для суммы double).
// ---------------------------------------------------------------------------
static void test_seg_bounds()
{
    constexpr std::size_t S = CpuExecutor::kDetSegments;
    bool ok_all = true, ok_big = true, ok_cover = true;
    {
        CpuExecutor exec(1);
        for (std::size_t n = 1; n <= 3000; ++n)
        {
            std::vector<std::size_t> seen(n, S + 1);
            exec.parallel_for_segments(n, [&](std::size_t s, std::size_t i) { seen[i] = s; });
            for (std::size_t k = 0; k < n; ++k)
                if (CpuExecutor::seg_of_index(k, n) != seen[k]) ok_all = false;
        }
    }
    // n ~ 2^40: проверяем границы каждого отрезка (b-1, b, b+1) и крайние индексы по определению отрезка.
    for (std::size_t n : {(std::size_t(1) << 40), (std::size_t(1) << 40) + 12345, (std::size_t(1) << 40) - 7})
    {
        for (std::size_t s = 0; s < S; ++s)
        {
            const std::size_t b = n * s / S, e = n * (s + 1) / S; // [b, e) -- отрезок s
            for (std::size_t k : {b, b + 1, e - 1, e, b == 0 ? std::size_t(0) : b - 1})
            {
                if (k >= n) continue;
                const std::size_t sk = CpuExecutor::seg_of_index(k, n);
                const bool in_seg = n * sk / S <= k && k < n * (sk + 1) / S;
                if (!in_seg) ok_big = false;
            }
            if (b < e && CpuExecutor::seg_of_index(b, n) != s) ok_big = false;
            if (b < e && CpuExecutor::seg_of_index(e - 1, n) != s) ok_big = false;
        }
        if (CpuExecutor::seg_of_index(n - 1, n) != S - 1) ok_big = false;
        if (CpuExecutor::seg_of_index(0, n) != 0) ok_big = false;
    }
    check(ok_all, "seg_of_index == seg из parallel_for_segments (n = 1..3000, все k)");
    check(ok_big, "seg_of_index при n ~ 2^40: k лежит в своём отрезке");

    // n у верхней границы size_t: n*s и 64*k+63 переполнились бы. seg_begin == floor(n*s/S) (эталон 128-битный),
    // seg_of_index кладёт границы отрезков в свой отрезок.
    bool ok_huge = true;
    for (std::size_t n : {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max() - 12345,
                          (std::size_t(1) << 63) + 1, (std::size_t(1) << 58) + 3})
    {
        for (std::size_t s = 0; s <= S; ++s)
        {
#if defined(__SIZEOF_INT128__)
            const auto ref = static_cast<std::size_t>(static_cast<unsigned __int128>(n) * s / S);
            if (CpuExecutor::seg_begin(n, s) != ref) ok_huge = false;
#endif
            if (s == S && CpuExecutor::seg_begin(n, s) != n) ok_huge = false;
        }
        for (std::size_t s = 0; s < S; ++s)
        {
            const std::size_t b = CpuExecutor::seg_begin(n, s), e = CpuExecutor::seg_begin(n, s + 1);
            if (b < e && CpuExecutor::seg_of_index(b, n) != s) ok_huge = false;
            if (b < e && CpuExecutor::seg_of_index(e - 1, n) != s) ok_huge = false;
        }
        if (CpuExecutor::seg_of_index(n - 1, n) != S - 1) ok_huge = false;
    }
    check(ok_huge, "seg_begin/seg_of_index при n ~ SIZE_MAX: без переполнения, границы в своих отрезках");

    std::mt19937 rng(91u);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    std::uniform_real_distribution<double> mag(-30.0, 30.0);
    bool ok_sum = true, ok_for = true, ok_cnt = true;
    for (int nt : {1, 3, 7, 8})
    {
        CpuExecutor exec(nt);
        for (std::size_t n : {std::size_t(0), std::size_t(1), std::size_t(5), std::size_t(64), std::size_t(1000),
                              std::size_t(12345)})
        {
            std::vector<double> x(n);
            for (auto& v : x) v = val(rng) * std::pow(10.0, mag(rng) / 3.0);
            for (double frac : {0.0, 0.05, 0.5, 1.0})
            {
                std::vector<std::size_t> sel;
                std::vector<std::uint16_t> seg;
                for (std::size_t s = 0; s < S; ++s)
                    for (std::size_t i = n * s / S; i < n * (s + 1) / S; ++i)
                        if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < frac)
                        {
                            sel.push_back(i);
                            seg.push_back(static_cast<std::uint16_t>(s));
                        }
                const std::size_t m = sel.size();
                std::size_t bound[S + 1];
                exec.seg_group_bounds(m, [&](std::size_t j) { return seg[j]; }, bound);
                // границы выбранных через seg_of_index (как у поиска с отсечкой): тот же результат
                std::size_t bound2[S + 1];
                for (std::size_t s = 0; s <= S; ++s)
                    bound2[s] = static_cast<std::size_t>(std::lower_bound(sel.begin(), sel.end(), n * s / S) - sel.begin());
                for (std::size_t s = 0; s <= S; ++s) ok_cover = ok_cover && bound[s] == bound2[s];
                for (std::size_t j = 0; j < m; ++j) ok_cover = ok_cover && CpuExecutor::seg_of_index(sel[j], n) == seg[j];

                const auto plus = [](double a, double b) { return a + b; };
                const double g = m ? exec.parallel_reduce_seg_groups(
                                         m, [&](std::size_t j) { return seg[j]; }, 0.0,
                                         [&](std::size_t j) { return x[sel[j]]; }, plus)
                                   : 0.0;
                const double b = exec.parallel_reduce_seg_bounds(bound, 0.0, [&](std::size_t j) { return x[sel[j]]; }, plus);
                ok_sum = ok_sum && std::memcmp(&g, &b, sizeof(double)) == 0;
                const std::size_t c = exec.parallel_reduce_seg_bounds(
                    bound, std::size_t(0), [&](std::size_t) { return std::size_t(1); },
                    [](std::size_t p, std::size_t q) { return p + q; });
                ok_cnt = ok_cnt && c == m;

                std::vector<std::vector<std::size_t>> seen_g(S), seen_b(S);
                exec.parallel_for_seg_groups(m, [&](std::size_t j) { return seg[j]; },
                                             [&](std::size_t s, std::size_t j) { seen_g[s].push_back(j); });
                exec.parallel_for_seg_bounds(bound, [&](std::size_t s, std::size_t j) { seen_b[s].push_back(j); });
                ok_for = ok_for && seen_g == seen_b;
            }
        }
    }
    check(ok_cover, "границы seg_group_bounds == lower_bound по границам отрезков; seg_of_index == seg выбранных");
    check(ok_sum, "parallel_reduce_seg_bounds: сумма double побитово = parallel_reduce_seg_groups (1/3/7/8 потоков)");
    check(ok_cnt, "parallel_reduce_seg_bounds: счётчик");
    check(ok_for, "parallel_for_seg_bounds: те же (s, j), что parallel_for_seg_groups");
}

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    std::printf("=== test_executor: CpuExecutor (notes/fiber-model/fiber-plan.md, этап Ф1) ===\n");

    test_parallel_for();
    test_parallel_reduce();
    test_parallel_scan();
    test_parallel_sort();
    test_seg_groups();
    test_seg_bounds();

    if (failures == 0)
        std::printf("test_executor: ПРОЙДЕН\n");
    else
        std::printf("test_executor: ПРОВАЛЕН (%d ошибок)\n", failures);

    return failures == 0 ? 0 : 1;
}
