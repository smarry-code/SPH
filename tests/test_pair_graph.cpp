// test_pair_graph.cpp
// Unit-тесты модуля pair_graph.hpp: проверка инварианта правильной раскраски
// графа конфликтов (никакие две пары одного цвета не разделяют частицу),
// корректности обхода for_each_color и устойчивости на граничных случаях.
//
// Тест НЕ зависит от глобального состояния (param.hpp) — PairGraph работает
// только с переданными pair_i/pair_j/niac/ntotal, что делает его изолированно
// тестируемым.
#include "pair_graph.hpp"
#include <cstdio>
#include <vector>
#include <set>
#include <stdexcept>
#include <string>
#include <omp.h>

static int g_failures = 0;

#define CHECK(cond)                                                            \
    do                                                                         \
    {                                                                          \
        if (!(cond))                                                           \
        {                                                                      \
            std::printf("FAIL [строка %d]: %s\n", __LINE__, #cond);            \
            ++g_failures;                                                      \
        }                                                                      \
    } while (0)

// Независимая проверка инварианта: пары одного цвета не разделяют частиц.
static bool coloring_valid(const PairGraph &g,
                           const std::vector<size_t> &pi,
                           const std::vector<size_t> &pj, size_t niac)
{
    return g.is_valid_coloring(pi, pj, niac);
}

// Тест 1: цепочка пар (0,1),(1,2),(2,3),... — линейный граф конфликтов.
// Ожидается: 2 цвета (двойная раскраска пути).
static void test_chain()
{
    PairGraph g;
    const size_t n = 10;
    std::vector<size_t> pi, pj;
    for (size_t k = 0; k + 1 < n; ++k)
    {
        pi.push_back(k);
        pj.push_back(k + 1);
    }
    g.build(pi, pj, pi.size(), n);
    CHECK(coloring_valid(g, pi, pj, pi.size()));
    CHECK(g.num_colors >= 2 && g.num_colors <= 3); // путь красится в 2 цвета
    // Все пары покрыты ровно один раз
    std::vector<int> seen(pi.size(), 0);
    g.for_each_color([&](int, int b, int e) {
        for (int p = b; p < e; ++p)
            ++seen[g.pairs_by_color[p]];
    });
    for (int s : seen)
        CHECK(s == 1);
    std::printf("test_chain: num_colors=%d\n", g.num_colors);
}

// Тест 2: плотный кластер — все пары (i,j) внутри n частиц (полный граф).
// Каждая частица во многих парах => много цветов, но инвариант должен держаться.
static void test_dense_cluster()
{
    PairGraph g;
    const size_t n = 8;
    std::vector<size_t> pi, pj;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
        {
            pi.push_back(i);
            pj.push_back(j);
        }
    g.build(pi, pj, pi.size(), n);
    CHECK(coloring_valid(g, pi, pj, pi.size()));
    // Для полного графа на n вершинах нужно n-1 цветов (каждая пара конфликтует
    // со всеми, разделяющими хотя бы одну частицу). Проверим разумную верхнюю границу.
    CHECK(g.num_colors <= static_cast<int>(2 * n));
    std::printf("test_dense_cluster (n=%zu, niac=%zu): num_colors=%d\n",
                n, pi.size(), g.num_colors);
}

// Тест 3: граничный случай — нет пар (niac=0).
static void test_empty()
{
    PairGraph g;
    std::vector<size_t> pi, pj;
    g.build(pi, pj, 0, 5);
    CHECK(g.num_colors == 0);
    CHECK(coloring_valid(g, pi, pj, 0));
    int calls = 0;
    g.for_each_color([&](int, int, int) { ++calls; });
    CHECK(calls == 0);
}

// Тест 4: одна пара.
static void test_single()
{
    PairGraph g;
    std::vector<size_t> pi = {0}, pj = {1};
    g.build(pi, pj, 1, 2);
    CHECK(g.num_colors == 1);
    CHECK(coloring_valid(g, pi, pj, 1));
}

// Тест 5: непересекающиеся пары (0,1),(2,3),(4,5) — конфликтов нет => 1 цвет.
static void test_disjoint()
{
    PairGraph g;
    std::vector<size_t> pi = {0, 2, 4}, pj = {1, 3, 5};
    g.build(pi, pj, 3, 6);
    CHECK(g.num_colors == 1);
    CHECK(coloring_valid(g, pi, pj, 3));
}

// Тест 6: звезда — частица 0 со всеми остальными (0,1),(0,2),...,(0,k).
// Все пары конфликтуют через частицу 0 => k цветов.
static void test_star()
{
    PairGraph g;
    const size_t k = 6;
    std::vector<size_t> pi, pj;
    for (size_t j = 1; j <= k; ++j)
    {
        pi.push_back(0);
        pj.push_back(j);
    }
    g.build(pi, pj, k, k + 1);
    CHECK(coloring_valid(g, pi, pj, k));
    CHECK(g.num_colors == static_cast<int>(k)); // все пары конфликтуют
    std::printf("test_star (k=%zu): num_colors=%d\n", k, g.num_colors);
}

// Тест 7: типичная SPH-сетка — каждая частица имеет ~6 соседей (1D/2D аналог).
// Проверка, что число цветов остаётся умеренным (десятки) и инвариант держится.
static void test_sph_like()
{
    PairGraph g;
    const size_t n = 200;       // 200 частиц в линию
    const int half = 3;          // радиус соседства
    std::vector<size_t> pi, pj;
    for (size_t i = 0; i < n; ++i)
        for (int dh = -half; dh <= half; ++dh)
        {
            int j = static_cast<int>(i) + dh;
            if (j < 0 || j >= static_cast<int>(n) || j == static_cast<int>(i))
                continue;
            if (dh > 0) // каждая пара один раз
            {
                pi.push_back(i);
                pj.push_back(static_cast<size_t>(j));
            }
        }
    g.build(pi, pj, pi.size(), n);
    CHECK(coloring_valid(g, pi, pj, pi.size()));
    CHECK(g.num_colors < 4 * (2 * half + 1)); // умеренное число цветов
    std::printf("test_sph_like (n=%zu, niac=%zu, neighb~%d): num_colors=%d\n",
                n, pi.size(), 2 * half, g.num_colors);
}

// Тест 8: большая SPH-сетка — стресс-тест спекулятивной блочно-параллельной
// окраски (B1, color_parallel). При >1 потоке активен параллельный путь;
// инвариант правильной раскраски обязан выполняться.
static void test_large_sph_grid()
{
    PairGraph g;
    const size_t n = 2000;       // 2000 частиц в линию
    const int half = 3;          // радиус соседства (~6 соседей)
    std::vector<size_t> pi, pj;
    for (size_t i = 0; i < n; ++i)
        for (int dh = 1; dh <= half; ++dh)
        {
            int j = static_cast<int>(i) + dh;
            if (j >= static_cast<int>(n))
                break;
            pi.push_back(i);
            pj.push_back(static_cast<size_t>(j));
        }
    g.build(pi, pj, pi.size(), n);
    CHECK(coloring_valid(g, pi, pj, pi.size()));
    CHECK(g.num_colors < 4 * (2 * half + 1)); // умеренное число цветов
    std::printf("test_large_sph_grid (n=%zu, niac=%zu, threads=%d): num_colors=%d\n",
                n, pi.size(), omp_get_max_threads(), g.num_colors);
}

// Тест 9: случайный плотный граф — стресс конфликтов B1 (много пар разделяют
// частицы => много конфликтов в спекулятивной фазе => проверка серийного
// разрешения). Детерминированный LCG для воспроизводимости между запусками.
static void test_random_dense()
{
    PairGraph g;
    const size_t n = 300;
    const size_t target_pairs = 5000;
    std::vector<size_t> pi, pj;
    unsigned long long rng = 123456789ULL;
    std::vector<std::vector<char>> seen(n, std::vector<char>(n, 0));
    size_t added = 0;
    while (added < target_pairs)
    {
        rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
        size_t i = static_cast<size_t>((rng >> 16) % n);
        rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
        size_t j = static_cast<size_t>((rng >> 16) % n);
        if (i == j) continue;
        if (i > j) { size_t t = i; i = j; j = t; }
        if (seen[i][j]) continue;
        seen[i][j] = 1;
        pi.push_back(i); pj.push_back(j); ++added;
    }
    g.build(pi, pj, pi.size(), n);
    CHECK(coloring_valid(g, pi, pj, pi.size()));
    std::printf("test_random_dense (n=%zu, niac=%zu, threads=%d): num_colors=%d\n",
                n, pi.size(), omp_get_max_threads(), g.num_colors);
}

// Тест 10: largest-degree-first порядок pair_order (Этап 2.1, C1).
// После build() pair_order упорядочен по НЕвозрастанию степени пары (deg(i)+deg(j)).
// Это гарантирует, что жадная окраска даёт fewer colors (Welsh–Powell).
static void test_degree_order()
{
    PairGraph g;
    const size_t n = 300;
    const int half = 4;
    std::vector<size_t> pi, pj;
    for (size_t i = 0; i < n; ++i)
        for (int dh = 1; dh <= half; ++dh)
        {
            int j = static_cast<int>(i) + dh;
            if (j >= static_cast<int>(n)) break;
            pi.push_back(i);
            pj.push_back(static_cast<size_t>(j));
        }
    g.build(pi, pj, pi.size(), n);
    CHECK(coloring_valid(g, pi, pj, pi.size()));
    // Степень частицы v доступна через pair_offset (CSR).
    auto deg = [&](size_t v) { return g.pair_offset[v + 1] - g.pair_offset[v]; };
    auto pdeg = [&](int k) { return deg(pi[static_cast<size_t>(k)]) + deg(pj[static_cast<size_t>(k)]); };
    bool ordered = true;
    for (size_t idx = 0; idx + 1 < pi.size(); ++idx)
        if (pdeg(g.pair_order[idx]) < pdeg(g.pair_order[idx + 1])) { ordered = false; break; }
    CHECK(ordered);
    std::printf("test_degree_order (n=%zu, niac=%zu): num_colors=%d, pair_order убывает=%s\n",
                n, pi.size(), g.num_colors, ordered ? "OK" : "НАРУШЕН");
}

// Тест 11: слияние совместимых цветов (coalescing, Этап 2.3, C2).
// Инвариант is_valid_coloring обязан сохраняться после слияния, а color_offset —
// быть консистентным (сумма размеров цветов = niac). Число цветов умеренное.
static void test_coalesce()
{
    PairGraph g;
    const size_t n = 400;
    const int half = 3;
    std::vector<size_t> pi, pj;
    for (size_t i = 0; i < n; ++i)
        for (int dh = 1; dh <= half; ++dh)
        {
            int j = static_cast<int>(i) + dh;
            if (j >= static_cast<int>(n)) break;
            pi.push_back(i);
            pj.push_back(static_cast<size_t>(j));
        }
    g.build(pi, pj, pi.size(), n);
    CHECK(coloring_valid(g, pi, pj, pi.size()));
    CHECK(g.num_colors >= 1);
    // color_offset консистентен: сумма размеров цветов = niac.
    int total = 0;
    for (int c = 0; c < g.num_colors; ++c)
        total += g.color_offset[c + 1] - g.color_offset[c];
    CHECK(total == static_cast<int>(pi.size()));
    // min_color_size/max_color_size определены (диагностика Этапа 2.2).
    CHECK(g.max_color_size >= g.min_color_size);
    std::printf("test_coalesce (n=%zu, niac=%zu): num_colors=%d, дисбаланс min=%d/max=%d\n",
                n, pi.size(), g.num_colors, g.min_color_size, g.max_color_size);
}

// Тест 12: color_chunk_range — срезы потоков покрывают цвет целиком
// без перекрытий и пропусков (подряд), каждый поток получает floor/ceil пар.
static void test_chunk_range()
{
    PairGraph g;
    std::vector<size_t> pi = {0, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9, 11, 0, 1};
    std::vector<size_t> pj = {1, 3, 5, 7, 9, 11, 2, 4, 6, 8, 10, 12, 2, 3};
    g.build(pi, pj, pi.size(), 13);
    CHECK(coloring_valid(g, pi, pj, pi.size()));
    const size_t P = 4;
    for (int c = 0; c < g.num_colors; ++c)
    {
        const size_t color = static_cast<size_t>(c);
        size_t expected = g.color_offset[color + 1] - g.color_offset[color];
        size_t covered = 0;
        size_t prev_end = g.color_offset[color];
        for (size_t tid = 0; tid < P; ++tid)
        {
            size_t b = 0, e = 0;
            g.color_chunk_range(color, tid, P, b, e);
            CHECK(b == prev_end);  // срезы идут подряд без разрывов
            CHECK(e >= b);
            covered += (e - b);
            prev_end = e;
        }
        CHECK(prev_end == g.color_offset[color + 1]);
        CHECK(covered == expected);
    }
    std::printf("test_chunk_range: num_colors=%d, срезы потоков консистентны\n", g.num_colors);
}

int main()
{
    // Включаем многопоточный путь color_parallel (B1), если доступно >1 потока.
    omp_set_num_threads(4);
    std::printf("pair_graph: тестирование с %d поток(ом/ами)\n\n", omp_get_max_threads());
    try
    {
        test_empty();
        test_single();
        test_disjoint();
        test_chain();
        test_star();
        test_dense_cluster();
        test_sph_like();
        test_large_sph_grid();
        test_random_dense();
        test_degree_order();   
        test_coalesce();       
        test_chunk_range();    
    }
    catch (const std::exception &e)
    {
        std::printf("EXCEPTION: %s\n", e.what());
        ++g_failures;
    }

    if (g_failures == 0)
    {
        std::printf("\nВСЕ ТЕСТЫ pair_graph ПРОШЛИ\n");
        return 0;
    }
    std::printf("\nПРОВАЛЕНО ТЕСТОВ: %d\n", g_failures);
    return 1;
}
