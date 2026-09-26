#ifndef RHS_SPATIAL_INDEX_HPP
#define RHS_SPATIAL_INDEX_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"
#include "precision.hpp"

// ============================================================================
// rhs/spatial_index.hpp — SpatialIndex<Dim>: обобщение прежней фиксированной
// сетки поиска соседей (FlatGrid, удалена вместе со старым путём на
// fiber/PLAN.md, этапе Ф4-B2) под требования fiber/SPH.md
// ("SpatialIndex<Fiber>", "SpatialIndex<Solid>", раздел "Search") и
// fiber/PLAN.md §3.2/§3.5.
//
// В отличие от прежней сетки:
//   - строится по SoA-массивам координат (не по ParticleSystem);
//   - размер ячейки — явный параметр build(), не выводится внутри из
//     hsml/scale_k;
//   - размерность — параметр шаблона (не завязан на глобальный constexpr
//     `dim` из param.hpp), поэтому один и тот же код обслуживает и Solid
//     (dim=2/3 проекта), и Fiber (dim=3, PLAN.md Р9), без привязки к
//     конкретной compute-структуре;
//   - умеет запросы "свой тип" (пары внутри одного набора, канонический
//     порядок i<j) и "чужой тип" (пары между двумя разными наборами,
//     например Fiber -> Solid).
//
// Алгоритм ячеек (counting sort: count -> prefix sum -> scatter) — тот же,
// что был у прежней сетки; сохранена конвенция происхождения (origin = min
// координата, число ячеек по оси = ceil(range/cell_size)+1).
// ============================================================================

namespace rhs
{

template <std::size_t Dim>
class SpatialIndex
{
public:
    using Coords = std::array<const real_t*, Dim>;

    SpatialIndex() = default;

    // Строит сетку по SoA-массивам координат coords[d][0..n) и явному размеру
    // ячейки cell_size. Координаты копируются во внутреннее хранилище:
    // SpatialIndex — часть RhsWorkspace и должен пережить возможную
    // перестройку readonly-состояния между вызовами RHS без dangling-
    // указателей (владение, а не заимствование).
    void build(Coords coords, std::size_t n, real_t cell_size)
    {
        n_ = n;
        cell_size_ = (cell_size > 0) ? cell_size : real_t(1.0);

        for (std::size_t d = 0; d < Dim; ++d)
        {
            coord_[d].assign(coords[d], coords[d] + n);
        }

        if (n_ == 0)
        {
            ncell_.fill(1);
            origin_.fill(real_t(0.0));
            offset_.assign(2, 0);
            sorted_idx_.clear();
            return;
        }

        std::array<real_t, Dim> lo{}, hi{};
        for (std::size_t d = 0; d < Dim; ++d)
        {
            lo[d] = hi[d] = coord_[d][0];
            for (std::size_t i = 1; i < n_; ++i)
            {
                lo[d] = std::min(lo[d], coord_[d][i]);
                hi[d] = std::max(hi[d], coord_[d][i]);
            }
            origin_[d] = lo[d];
            int nc_d = static_cast<int>(std::ceil((hi[d] - lo[d]) / cell_size_)) + 1;
            ncell_[d] = std::max(nc_d, 1);
        }

        const std::size_t nc = total_cells();
        std::vector<std::size_t> count(nc, 0);
        offset_.assign(nc + 1, 0);
        sorted_idx_.resize(n_);

        for (std::size_t i = 0; i < n_; ++i) count[cell_of(i)]++;
        for (std::size_t c = 0; c < nc; ++c) offset_[c + 1] = offset_[c] + count[c];

        std::vector<std::size_t> cursor(offset_.begin(), offset_.begin() + static_cast<std::ptrdiff_t>(nc));
        for (std::size_t i = 0; i < n_; ++i)
        {
            const std::size_t c = cell_of(i);
            sorted_idx_[cursor[c]++] = i;
        }
    }

    std::size_t size() const { return n_; }
    real_t cell_size() const { return cell_size_; }

    // "Свой тип": пары (i, j), i < j, |x_i - x_j| <= radius, accept(i, j).
    // accept — дополнительный фильтр (например, object_id[i] != object_id[j]
    // для Fiber-Fiber-Contact); по умолчанию пропускает все пары.
    template <class Accept = bool (*)(index_t, index_t)>
    std::vector<std::pair<index_t, index_t>> find_pairs(
        real_t radius, Accept accept = [](index_t, index_t) { return true; }) const
    {
        std::vector<std::pair<index_t, index_t>> out;
        if (n_ == 0) return out;
        const std::size_t nc = total_cells();
        for (std::size_t c = 0; c < nc; ++c) collect_same_type_cell(c, radius, accept, out);
        return out;
    }

    // "Чужой тип": пары (a, b), a — индекс во внешнем множестве other_coords
    // (n_other частиц), b — индекс в this (проиндексированном множестве).
    // |x_a - x_b| <= radius, accept(a, b). Направление a -> b не создаёт
    // дублей естественно, так как a и b принадлежат разным наборам (SPH.md:
    // "Fiber -> Solid само обеспечивает однократную генерацию пары").
    template <class Accept = bool (*)(index_t, index_t)>
    std::vector<std::pair<index_t, index_t>> find_pairs_cross(
        Coords other_coords, std::size_t n_other, real_t radius,
        Accept accept = [](index_t, index_t) { return true; }) const
    {
        std::vector<std::pair<index_t, index_t>> out;
        if (n_ == 0 || n_other == 0) return out;
        for (std::size_t a = 0; a < n_other; ++a) collect_cross_one(other_coords, a, radius, accept, out);
        return out;
    }

    // Параллельные варианты (SPH.md, "Search": каждый worker пишет найденные
    // пары в собственный буфер, затем буферы объединяются в один непрерывный
    // массив). Итоговое МНОЖЕСТВО пар совпадает с последовательными
    // find_pairs*/find_pairs_cross с точностью до порядка следования.
    template <class Accept = bool (*)(index_t, index_t)>
    std::vector<std::pair<index_t, index_t>> find_pairs_parallel(
        const CpuExecutor& exec, real_t radius,
        Accept accept = [](index_t, index_t) { return true; }) const
    {
        const std::size_t nc = total_cells();
        return parallel_merge(exec, nc, [&](std::size_t c, std::vector<std::pair<index_t, index_t>>& local) {
            collect_same_type_cell(c, radius, accept, local);
        });
    }

    template <class Accept = bool (*)(index_t, index_t)>
    std::vector<std::pair<index_t, index_t>> find_pairs_cross_parallel(
        const CpuExecutor& exec, Coords other_coords, std::size_t n_other, real_t radius,
        Accept accept = [](index_t, index_t) { return true; }) const
    {
        return parallel_merge(exec, n_other, [&](std::size_t a, std::vector<std::pair<index_t, index_t>>& local) {
            collect_cross_one(other_coords, a, radius, accept, local);
        });
    }

private:
    std::array<std::vector<real_t>, Dim> coord_;
    std::size_t n_ = 0;
    real_t cell_size_ = real_t(1.0);
    std::array<real_t, Dim> origin_{};
    std::array<int, Dim> ncell_{};
    std::vector<std::size_t> offset_;
    std::vector<std::size_t> sorted_idx_;

    std::size_t total_cells() const
    {
        std::size_t nc = 1;
        for (std::size_t d = 0; d < Dim; ++d) nc *= static_cast<std::size_t>(ncell_[d]);
        return nc;
    }

    std::array<int, Dim> cell_of_point(const std::array<real_t, Dim>& p) const
    {
        std::array<int, Dim> c{};
        for (std::size_t d = 0; d < Dim; ++d)
        {
            const int idx = static_cast<int>(std::floor((p[d] - origin_[d]) / cell_size_));
            c[d] = std::clamp(idx, 0, ncell_[d] - 1);
        }
        return c;
    }

    std::array<int, Dim> cell_of_index(std::size_t i) const
    {
        std::array<real_t, Dim> p{};
        for (std::size_t d = 0; d < Dim; ++d) p[d] = coord_[d][i];
        return cell_of_point(p);
    }

    std::size_t cell_of(std::size_t i) const { return linear_index(cell_of_index(i)); }

    std::size_t linear_index(const std::array<int, Dim>& c) const
    {
        std::size_t lin = 0;
        for (std::size_t k = Dim; k-- > 0;)
        {
            lin = lin * static_cast<std::size_t>(ncell_[k]) + static_cast<std::size_t>(c[k]);
        }
        return lin;
    }

    std::array<int, Dim> decode(std::size_t lin) const
    {
        std::array<int, Dim> c{};
        for (std::size_t d = 0; d < Dim; ++d)
        {
            c[d] = static_cast<int>(lin % static_cast<std::size_t>(ncell_[d]));
            lin /= static_cast<std::size_t>(ncell_[d]);
        }
        return c;
    }

    std::pair<std::size_t, std::size_t> cell_range(std::size_t lin) const
    {
        return {offset_[lin], offset_[lin + 1]};
    }

    real_t dist2(std::size_t i, std::size_t j) const
    {
        real_t d2 = 0;
        for (std::size_t d = 0; d < Dim; ++d)
        {
            const real_t dd = coord_[d][i] - coord_[d][j];
            d2 += dd * dd;
        }
        return d2;
    }

    // Сколько слоёв соседних ячеек в каждую сторону нужно обойти, чтобы
    // гарантированно найти все точки в пределах radius. Если cell_size >=
    // radius, достаточно одного слоя (offset -1..0..1, как раньше). Но
    // SpatialIndex строится с ПРОИЗВОЛЬНЫМ cell_size (build() принимает его
    // явным параметром — в отличие от прежней сетки, где ячейка выводилась
    // из hsml), и ничто не гарантирует cell_size >= radius: например,
    // радиус поиска может отличаться от размера ячейки, используемого при
    // построении индекса. Тогда пара точек на расстоянии <= radius может
    // лежать в ячейках, отстоящих друг от друга более чем на один шаг
    // сетки — фиксированный 3^Dim стенсил такие пары молча теряет. Formula
    // floor(radius/cell_size)+1 — стандартная безопасная (не обязательно
    // самая узкая) оценка: |floor(a/h)-floor(b/h)| <= floor(|a-b|/h)+1 для
    // любых a,b,h>0; использована вместо ceil(...), чтобы не зависеть от
    // округления чисел с плавающей точкой around точных кратных cell_size.
    int neighbor_reach(real_t radius) const
    {
        return static_cast<int>(std::floor(radius / cell_size_)) + 1;
    }

    // Обходит все существующие ячейки в кубе (2*reach+1)^Dim вокруг center
    // (clamp по границам сетки — несуществующие соседи просто пропускаются,
    // дублей не возникает: разным смещениям всегда соответствуют разные
    // многоиндексы).
    template <class F>
    void for_each_neighbor_cell(const std::array<int, Dim>& center, int reach, F&& f) const
    {
        const std::size_t span = static_cast<std::size_t>(2 * reach + 1);
        std::size_t combos = 1;
        for (std::size_t d = 0; d < Dim; ++d) combos *= span;
        for (std::size_t k = 0; k < combos; ++k)
        {
            std::size_t kk = k;
            std::array<int, Dim> c{};
            bool valid = true;
            for (std::size_t d = 0; d < Dim; ++d)
            {
                const int off = static_cast<int>(kk % span) - reach;
                kk /= span;
                const int v = center[d] + off;
                if (v < 0 || v >= ncell_[d])
                {
                    valid = false;
                    break;
                }
                c[d] = v;
            }
            if (valid) f(linear_index(c));
        }
    }

    template <class Accept>
    void collect_same_type_cell(std::size_t cell_lin, real_t radius, const Accept& accept,
                                 std::vector<std::pair<index_t, index_t>>& out) const
    {
        const real_t r2 = radius * radius;
        auto [b, e] = cell_range(cell_lin);
        if (b == e) return;
        const std::array<int, Dim> center = decode(cell_lin);
        const int reach = neighbor_reach(radius);
        for (std::size_t pi = b; pi < e; ++pi)
        {
            const std::size_t i = sorted_idx_[pi];
            for_each_neighbor_cell(center, reach, [&](std::size_t nb_lin) {
                auto [nb, ne] = cell_range(nb_lin);
                for (std::size_t pj = nb; pj < ne; ++pj)
                {
                    const std::size_t j = sorted_idx_[pj];
                    if (j <= i) continue; // канонический порядок i<j, без дублей
                    if (dist2(i, j) <= r2 && accept(i, j)) out.emplace_back(i, j);
                }
            });
        }
    }

    template <class Accept>
    void collect_cross_one(Coords other_coords, std::size_t a, real_t radius, const Accept& accept,
                            std::vector<std::pair<index_t, index_t>>& out) const
    {
        const real_t r2 = radius * radius;
        std::array<real_t, Dim> p{};
        for (std::size_t d = 0; d < Dim; ++d) p[d] = other_coords[d][a];
        const std::array<int, Dim> cell = cell_of_point(p);
        const int reach = neighbor_reach(radius);

        for_each_neighbor_cell(cell, reach, [&](std::size_t nb_lin) {
            auto [nb, ne] = cell_range(nb_lin);
            for (std::size_t pb = nb; pb < ne; ++pb)
            {
                const std::size_t b = sorted_idx_[pb];
                real_t d2 = 0;
                for (std::size_t d = 0; d < Dim; ++d)
                {
                    const real_t dd = p[d] - coord_[d][b];
                    d2 += dd * dd;
                }
                if (d2 <= r2 && accept(a, b)) out.emplace_back(a, b);
            }
        });
    }

    template <class Body>
    std::vector<std::pair<index_t, index_t>> parallel_merge(const CpuExecutor& exec, std::size_t range_n,
                                                              Body&& body) const
    {
        const int nw = std::max(1, exec.num_workers());
        std::vector<std::vector<std::pair<index_t, index_t>>> local(static_cast<std::size_t>(nw));
        exec.parallel_for(range_n, [&](std::size_t idx) {
            int w = CpuExecutor::worker_id();
            if (w < 0 || w >= nw) w = 0;
            body(idx, local[static_cast<std::size_t>(w)]);
        });
        std::size_t total = 0;
        for (auto& v : local) total += v.size();
        std::vector<std::pair<index_t, index_t>> out;
        out.reserve(total);
        for (auto& v : local) out.insert(out.end(), v.begin(), v.end());
        return out;
    }
};

} // namespace rhs

#endif // RHS_SPATIAL_INDEX_HPP
