// tests/test_spatial_index.cpp
// Этап Ф1 (fiber/PLAN.md): unit-тест SpatialIndex<Dim> (include/rhs/spatial_index.hpp).
// Проверяет, что пары, найденные SpatialIndex, совпадают с результатом
// полного перебора (brute force) — для "своего типа" (внутри одного набора)
// и для "чужого типа" (между двумя разными наборами), на нескольких
// случайных конфигурациях разной плотности, включая последовательные и
// параллельные (через CpuExecutor) запросы.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>
#ifdef _WIN32
#define NOMINMAX // std::max/std::min из <rhs/...> конфликтуют с макросами min/max windows.h
#include <windows.h>
#endif

#include "rhs/cpu_executor.hpp"
#include "rhs/spatial_index.hpp"

using rhs::CpuExecutor;
using rhs::index_t;
using rhs::SpatialIndex;

static int failures = 0;

static void check(bool cond, const char* msg)
{
    if (!cond)
    {
        std::printf("FAIL: %s\n", msg);
        ++failures;
    }
}

using PairVec = std::vector<std::pair<index_t, index_t>>;

static PairVec sorted_copy(PairVec v)
{
    std::sort(v.begin(), v.end());
    return v;
}

template <std::size_t Dim>
struct Cloud
{
    std::array<std::vector<double>, Dim> coord;
    std::size_t n = 0;

    std::array<const double*, Dim> ptrs() const
    {
        std::array<const double*, Dim> p{};
        for (std::size_t d = 0; d < Dim; ++d) p[d] = coord[d].data();
        return p;
    }
};

template <std::size_t Dim>
static Cloud<Dim> random_cloud(std::size_t n, double box, std::mt19937& rng)
{
    Cloud<Dim> c;
    c.n = n;
    std::uniform_real_distribution<double> dist(-box, box);
    for (std::size_t d = 0; d < Dim; ++d)
    {
        c.coord[d].resize(n);
        for (std::size_t i = 0; i < n; ++i) c.coord[d][i] = dist(rng);
    }
    return c;
}

template <std::size_t Dim>
static double distance(const Cloud<Dim>& c, std::size_t i, std::size_t j)
{
    double d2 = 0.0;
    for (std::size_t d = 0; d < Dim; ++d)
    {
        double dd = c.coord[d][i] - c.coord[d][j];
        d2 += dd * dd;
    }
    return std::sqrt(d2);
}

template <std::size_t Dim, class Accept>
static PairVec brute_force_same(const Cloud<Dim>& c, double radius, Accept accept)
{
    PairVec out;
    for (std::size_t i = 0; i < c.n; ++i)
        for (std::size_t j = i + 1; j < c.n; ++j)
            if (distance(c, i, j) <= radius && accept(i, j)) out.emplace_back(i, j);
    return out;
}

template <std::size_t Dim, class Accept>
static PairVec brute_force_cross(const Cloud<Dim>& a, const Cloud<Dim>& b, double radius, Accept accept)
{
    PairVec out;
    for (std::size_t ia = 0; ia < a.n; ++ia)
    {
        for (std::size_t ib = 0; ib < b.n; ++ib)
        {
            double d2 = 0.0;
            for (std::size_t d = 0; d < Dim; ++d)
            {
                double dd = a.coord[d][ia] - b.coord[d][ib];
                d2 += dd * dd;
            }
            if (std::sqrt(d2) <= radius && accept(ia, ib)) out.emplace_back(ia, ib);
        }
    }
    return out;
}

// -----------------------------------------------------------------------
// "Свой тип": несколько случайных облаков разной плотности (варьируем число
// частиц в фиксированном объёме и размер ячейки относительно radius).
// -----------------------------------------------------------------------
template <std::size_t Dim>
static void test_same_type()
{
    std::mt19937 rng(1000 + static_cast<unsigned>(Dim));
    const std::vector<std::size_t> sizes = {0, 1, 2, 50, 400};
    const std::vector<double> boxes = {1.0, 5.0, 20.0}; // разная плотность при фиксированном N
    const std::vector<double> radii = {0.3, 1.0, 3.0};
    const std::vector<double> cell_scales = {0.5, 1.0, 2.5}; // размер ячейки относительно radius

    for (std::size_t n : sizes)
    {
        for (double box : boxes)
        {
            Cloud<Dim> cloud = random_cloud<Dim>(n, box, rng);
            auto coords = cloud.ptrs();
            for (double radius : radii)
            {
                for (double cs : cell_scales)
                {
                    SpatialIndex<Dim> index;
                    index.build(coords, n, radius * cs);

                    // Фильтр "object_id": чётность индекса — имитация Fiber-Fiber-Contact
                    // (object_id[i] != object_id[j]).
                    auto accept = [](index_t i, index_t j) { return (i % 2) != (j % 2); };

                    PairVec expected = sorted_copy(brute_force_same<Dim>(cloud, radius, accept));
                    PairVec actual = sorted_copy(index.find_pairs(radius, accept));
                    check(actual == expected, "SpatialIndex::find_pairs (свой тип, с фильтром) расходится с brute force");

                    // Без фильтра (все пары).
                    auto accept_all = [](index_t, index_t) { return true; };
                    PairVec expected_all = sorted_copy(brute_force_same<Dim>(cloud, radius, accept_all));
                    PairVec actual_all = sorted_copy(index.find_pairs(radius, accept_all));
                    check(actual_all == expected_all,
                          "SpatialIndex::find_pairs (свой тип, без фильтра) расходится с brute force");

                    // Параллельный вариант должен находить то же множество пар.
                    CpuExecutor exec(4);
                    PairVec actual_par = sorted_copy(index.find_pairs_parallel(exec, radius, accept_all));
                    check(actual_par == expected_all, "SpatialIndex::find_pairs_parallel расходится с sequential find_pairs");
                }
            }
        }
    }
}

// -----------------------------------------------------------------------
// "Чужой тип": два независимых облака (например, Fiber vs Solid), с разной
// плотностью и разным числом частиц в каждом из наборов.
// -----------------------------------------------------------------------
template <std::size_t Dim>
static void test_cross_type()
{
    std::mt19937 rng(2000 + static_cast<unsigned>(Dim));
    const std::vector<std::pair<std::size_t, std::size_t>> size_pairs = {
        {0, 0}, {0, 10}, {10, 0}, {1, 1}, {5, 60}, {200, 30}};
    const std::vector<double> boxes = {1.0, 8.0};
    const std::vector<double> radii = {0.4, 2.0};

    for (auto [na, nb] : size_pairs)
    {
        for (double box : boxes)
        {
            Cloud<Dim> a = random_cloud<Dim>(na, box, rng);
            Cloud<Dim> b = random_cloud<Dim>(nb, box, rng);
            auto coords_b = b.ptrs();
            auto coords_a = a.ptrs();

            for (double radius : radii)
            {
                SpatialIndex<Dim> index_b;
                index_b.build(coords_b, nb, radius);

                auto accept = [](index_t ia, index_t ib) { return (ia + ib) % 3 != 0; };

                PairVec expected = sorted_copy(brute_force_cross<Dim>(a, b, radius, accept));
                PairVec actual = sorted_copy(index_b.find_pairs_cross(coords_a, na, radius, accept));
                check(actual == expected, "SpatialIndex::find_pairs_cross расходится с brute force");

                CpuExecutor exec(3);
                PairVec actual_par = sorted_copy(index_b.find_pairs_cross_parallel(exec, coords_a, na, radius, accept));
                check(actual_par == expected, "SpatialIndex::find_pairs_cross_parallel расходится с sequential");
            }
        }
    }
}

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    std::printf("=== test_spatial_index: SpatialIndex<Dim> (fiber/PLAN.md, этап Ф1) ===\n");

    test_same_type<2>();
    test_same_type<3>();
    test_cross_type<2>();
    test_cross_type<3>();

    if (failures == 0)
        std::printf("test_spatial_index: ПРОЙДЕН\n");
    else
        std::printf("test_spatial_index: ПРОВАЛЕН (%d ошибок)\n", failures);

    return failures == 0 ? 0 : 1;
}
