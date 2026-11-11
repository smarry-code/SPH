// tests/test_spatial_index.cpp
// Этап Ф1 (fiber/PLAN.md): unit-тест SpatialIndex<Dim> (include/rhs/spatial_index.hpp).
// Проверяет, что пары, найденные SpatialIndex, совпадают с результатом
// полного перебора (brute force) — для "своего типа" (внутри одного набора)
// и для "чужого типа" (между двумя разными наборами), на нескольких
// случайных конфигурациях разной плотности, включая последовательные и
// параллельные (через CpuExecutor) запросы, а также повторные запросы с
// выходным параметром (find_pairs_parallel(..., out&)) на одном и том же
// индексе, перестраиваемом с разным N.
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

// -----------------------------------------------------------------------
// Переиспользование: один SpatialIndex и один out, build() с разным N
// (рост/спад), find_pairs_parallel/find_pairs_cross_parallel с out& дают те
// же пары, что варианты по значению. При одном потоке порядок пар должен
// совпадать точно (worker-major, порядок обхода), при нескольких — как
// множество.
// -----------------------------------------------------------------------
template <std::size_t Dim>
static void test_reuse_out()
{
    std::mt19937 rng(3000 + static_cast<unsigned>(Dim));
    const std::vector<std::size_t> sizes = {300, 20, 0, 600, 5, 150, 600};
    const double radius = 1.5;
    auto accept = [](index_t i, index_t j) { return (i + j) % 5 != 0; };

    for (int nthreads : {1, 4})
    {
        CpuExecutor exec(nthreads);
        SpatialIndex<Dim> index; // один и тот же индекс на все перестроения
        PairVec out_same, out_cross;
        bool ok_same = true, ok_cross = true, ok_ref = true;
        for (std::size_t n : sizes)
        {
            Cloud<Dim> cloud = random_cloud<Dim>(n, 6.0, rng);
            Cloud<Dim> other = random_cloud<Dim>(n / 2 + 1, 6.0, rng);
            auto coords = cloud.ptrs();
            auto other_coords = other.ptrs();
            index.build(coords, n, radius);

            // Свежий индекс — эталон "по значению".
            SpatialIndex<Dim> fresh;
            fresh.build(coords, n, radius);
            PairVec by_value = fresh.find_pairs_parallel(exec, radius, accept);
            PairVec by_value_cross = fresh.find_pairs_cross_parallel(exec, other_coords, other.n, radius, accept);

            // Повторять вызов на одном индексе и одном out: результат не зависит от прошлых вызовов.
            for (int rep = 0; rep < 2; ++rep)
            {
                index.find_pairs_parallel(exec, radius, accept, out_same);
                index.find_pairs_cross_parallel(exec, other_coords, other.n, radius, accept, out_cross);
                if (nthreads == 1)
                {
                    if (out_same != by_value) ok_same = false;
                    if (out_cross != by_value_cross) ok_cross = false;
                }
                else
                {
                    if (sorted_copy(out_same) != sorted_copy(by_value)) ok_same = false;
                    if (sorted_copy(out_cross) != sorted_copy(by_value_cross)) ok_cross = false;
                }
            }
            if (sorted_copy(out_same) != sorted_copy(brute_force_same<Dim>(cloud, radius, accept))) ok_ref = false;
        }
        check(ok_same, "SpatialIndex::find_pairs_parallel(out&) при повторных вызовах расходится с вариантом по значению");
        check(ok_cross,
              "SpatialIndex::find_pairs_cross_parallel(out&) при повторных вызовах расходится с вариантом по значению");
        check(ok_ref, "SpatialIndex::find_pairs_parallel(out&) после перестроений расходится с brute force");
    }
}


// -----------------------------------------------------------------------
// Порядок пар после отсечения ячеек (for_each_near_cell) -- побитово прежний: эталон -- копия прежнего обхода
// (все (2*reach+1)^Dim ячеек вокруг ячейки точки, ось 0 быстрее всех, без отсечения). Вектор пар сравнивается
// целиком, с порядком. Облака: случайные, решётка с шагом, кратным ячейке (расстояния ровно = radius), точки
// вне сетки для поиска чужого типа; radius = h, меньше и больше h (reach 1..3).
// -----------------------------------------------------------------------
template <std::size_t Dim>
struct RefGrid
{
    double h = 1.0;
    std::array<double, Dim> origin{};
    std::array<int, Dim> nc{};
    std::vector<std::size_t> offset, sorted;
    const Cloud<Dim>* c = nullptr;

    std::array<int, Dim> cell_of(const std::array<double, Dim>& p) const
    {
        std::array<int, Dim> r{};
        for (std::size_t d = 0; d < Dim; ++d)
            r[d] = std::clamp(static_cast<int>(std::floor((p[d] - origin[d]) / h)), 0, nc[d] - 1);
        return r;
    }
    std::size_t lin(const std::array<int, Dim>& q) const
    {
        std::size_t l = 0;
        for (std::size_t k = Dim; k-- > 0;) l = l * static_cast<std::size_t>(nc[k]) + static_cast<std::size_t>(q[k]);
        return l;
    }
    void build(const Cloud<Dim>& cloud, double cell)
    {
        c = &cloud;
        h = cell;
        for (std::size_t d = 0; d < Dim; ++d)
        {
            double lo = cloud.coord[d][0], hi = lo;
            for (std::size_t i = 1; i < cloud.n; ++i)
            {
                lo = std::min(lo, cloud.coord[d][i]);
                hi = std::max(hi, cloud.coord[d][i]);
            }
            origin[d] = lo;
            nc[d] = std::max(static_cast<int>(std::ceil((hi - lo) / h)) + 1, 1);
        }
        std::size_t total = 1;
        for (std::size_t d = 0; d < Dim; ++d) total *= static_cast<std::size_t>(nc[d]);
        std::vector<std::size_t> cnt(total, 0), cell_i(cloud.n);
        for (std::size_t i = 0; i < cloud.n; ++i)
        {
            std::array<double, Dim> p{};
            for (std::size_t d = 0; d < Dim; ++d) p[d] = cloud.coord[d][i];
            cell_i[i] = lin(cell_of(p));
            ++cnt[cell_i[i]];
        }
        offset.assign(total + 1, 0);
        for (std::size_t q = 0; q < total; ++q) offset[q + 1] = offset[q] + cnt[q];
        std::vector<std::size_t> cur(offset.begin(), offset.end() - 1);
        sorted.resize(cloud.n);
        for (std::size_t i = 0; i < cloud.n; ++i) sorted[cur[cell_i[i]]++] = i;
    }
    template <class F>
    void neighbors(const std::array<int, Dim>& center, int reach, F&& f) const
    {
        const std::size_t span = static_cast<std::size_t>(2 * reach + 1);
        std::size_t combos = 1;
        for (std::size_t d = 0; d < Dim; ++d) combos *= span;
        for (std::size_t k = 0; k < combos; ++k)
        {
            std::size_t kk = k;
            std::array<int, Dim> q{};
            bool ok = true;
            for (std::size_t d = 0; d < Dim; ++d)
            {
                const int v = center[d] + static_cast<int>(kk % span) - reach;
                kk /= span;
                if (v < 0 || v >= nc[d])
                {
                    ok = false;
                    break;
                }
                q[d] = v;
            }
            if (ok) f(lin(q));
        }
    }
    PairVec same(double radius) const
    {
        PairVec out;
        const double r2 = radius * radius;
        const int reach = static_cast<int>(std::floor(radius / h)) + 1;
        std::size_t total = offset.size() - 1;
        for (std::size_t q = 0; q < total; ++q)
        {
            if (offset[q] == offset[q + 1]) continue;
            std::array<int, Dim> center{};
            std::size_t l = q;
            for (std::size_t d = 0; d < Dim; ++d)
            {
                center[d] = static_cast<int>(l % static_cast<std::size_t>(nc[d]));
                l /= static_cast<std::size_t>(nc[d]);
            }
            for (std::size_t pi = offset[q]; pi < offset[q + 1]; ++pi)
            {
                const std::size_t i = sorted[pi];
                neighbors(center, reach, [&](std::size_t nb) {
                    for (std::size_t pj = offset[nb]; pj < offset[nb + 1]; ++pj)
                    {
                        const std::size_t j = sorted[pj];
                        if (j <= i) continue;
                        double d2 = 0.0;
                        for (std::size_t d = 0; d < Dim; ++d)
                        {
                            const double dd = c->coord[d][i] - c->coord[d][j];
                            d2 += dd * dd;
                        }
                        if (d2 <= r2) out.emplace_back(i, j);
                    }
                });
            }
        }
        return out;
    }
    PairVec cross(const Cloud<Dim>& other, double radius) const
    {
        PairVec out;
        const double r2 = radius * radius;
        const int reach = static_cast<int>(std::floor(radius / h)) + 1;
        for (std::size_t a = 0; a < other.n; ++a)
        {
            std::array<double, Dim> p{};
            for (std::size_t d = 0; d < Dim; ++d) p[d] = other.coord[d][a];
            neighbors(cell_of(p), reach, [&](std::size_t nb) {
                for (std::size_t pb = offset[nb]; pb < offset[nb + 1]; ++pb)
                {
                    const std::size_t b = sorted[pb];
                    double d2 = 0.0;
                    for (std::size_t d = 0; d < Dim; ++d)
                    {
                        const double dd = p[d] - c->coord[d][b];
                        d2 += dd * dd;
                    }
                    if (d2 <= r2) out.emplace_back(a, b);
                }
            });
        }
        return out;
    }
};

template <std::size_t Dim>
static Cloud<Dim> lattice_cloud(int m, double step, double shift)
{
    Cloud<Dim> c;
    std::size_t n = 1;
    for (std::size_t d = 0; d < Dim; ++d) n *= static_cast<std::size_t>(m);
    c.n = n;
    for (std::size_t d = 0; d < Dim; ++d) c.coord[d].resize(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::size_t l = i;
        for (std::size_t d = 0; d < Dim; ++d)
        {
            c.coord[d][i] = shift + step * static_cast<double>(l % static_cast<std::size_t>(m));
            l /= static_cast<std::size_t>(m);
        }
    }
    return c;
}

template <std::size_t Dim>
static void test_pruned_order()
{
    std::mt19937 rng(20261112u);
    rhs::CpuExecutor exec1(1), exec8(8);
    const auto all = [](index_t, index_t) { return true; };
    std::vector<Cloud<Dim>> clouds;
    clouds.push_back(random_cloud<Dim>(Dim == 3 ? 600 : 400, 1.0, rng));
    clouds.push_back(lattice_cloud<Dim>(Dim == 3 ? 7 : 15, 0.25, 0.1));     // расстояния ровно 0.25, 0.5 ...
    clouds.push_back(lattice_cloud<Dim>(Dim == 3 ? 6 : 12, 0.3, -0.731));   // шаг не кратен ячейке
    const double cells[] = {0.25, 0.5, 0.3};
    const double radii_k[] = {1.0, 0.7, 1.6, 2.5};
    bool ok_same = true, ok_cross = true, ok_par = true;
    for (const auto& cl : clouds)
    {
        for (double cell : cells)
        {
            for (double rk : radii_k)
            {
                const double radius = rk * cell;
                rhs::SpatialIndex<Dim> idx;
                idx.build(cl.ptrs(), cl.n, cell);
                RefGrid<Dim> ref;
                ref.build(cl, cell);
                const PairVec expect = ref.same(radius);
                ok_same = ok_same && (idx.find_pairs(radius, all) == expect);
                PairVec par;
                idx.find_pairs_parallel(exec8, radius, all, par);
                PairVec par1;
                idx.find_pairs_parallel(exec1, radius, all, par1);
                // parallel: отрезки ячеек склеиваются по порядку -- тот же вектор, что у последовательного.
                ok_par = ok_par && (par == expect) && (par1 == expect);

                // чужой тип: точки облака b внутри и далеко вне сетки, на границах ячеек
                Cloud<Dim> other = random_cloud<Dim>(300, 3.0, rng);
                for (std::size_t i = 0; i < 40 && i < other.n; ++i)
                    for (std::size_t d = 0; d < Dim; ++d) other.coord[d][i] = cl.coord[d][i % cl.n] + (i % 3) * cell;
                const PairVec expect_x = ref.cross(other, radius);
                ok_cross = ok_cross && (idx.find_pairs_cross(other.ptrs(), other.n, radius, all) == expect_x);
                PairVec px, px1;
                idx.find_pairs_cross_parallel(exec8, other.ptrs(), other.n, radius, all, px);
                idx.find_pairs_cross_parallel(exec1, other.ptrs(), other.n, radius, all, px1);
                ok_par = ok_par && (px == expect_x) && (px1 == expect_x);
            }
        }
    }
    check(ok_same, "find_pairs: пары и их ПОРЯДОК совпадают с прежним обходом всех соседних ячеек");
    check(ok_cross, "find_pairs_cross: пары и порядок совпадают с прежним обходом (точки вне сетки)");
    check(ok_par, "parallel-варианты на 1 и 8 потоках: тот же вектор пар, что у прежнего обхода");
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
    test_reuse_out<2>();
    test_reuse_out<3>();
    test_pruned_order<2>();
    test_pruned_order<3>();

    if (failures == 0)
        std::printf("test_spatial_index: ПРОЙДЕН\n");
    else
        std::printf("test_spatial_index: ПРОВАЛЕН (%d ошибок)\n", failures);

    return failures == 0 ? 0 : 1;
}
