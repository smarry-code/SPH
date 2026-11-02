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
    using Coords = std::array<const double*, Dim>;

    SpatialIndex() = default;

    // Строит сетку по SoA-массивам координат coords[d][0..n) и явному размеру
    // ячейки cell_size. Координаты копируются во внутреннее хранилище:
    // SpatialIndex — часть RhsWorkspace и должен пережить возможную
    // перестройку readonly-состояния между вызовами RHS без dangling-
    // указателей (владение, а не заимствование).
    void build(Coords coords, std::size_t n, double cell_size)
    {
        n_ = n;
        cell_size_ = (cell_size > 0.0) ? cell_size : 1.0;

        for (std::size_t d = 0; d < Dim; ++d)
        {
            coord_[d].assign(coords[d], coords[d] + n);
        }

        if (n_ == 0)
        {
            ncell_.fill(1);
            origin_.fill(0.0);
            offset_.assign(2, 0);
            sorted_idx_.clear();
            return;
        }

        std::array<double, Dim> lo{}, hi{};
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
        count_.assign(nc, 0);
        offset_.assign(nc + 1, 0);
        sorted_idx_.resize(n_);

        for (std::size_t i = 0; i < n_; ++i) count_[cell_of(i)]++;
        for (std::size_t c = 0; c < nc; ++c) offset_[c + 1] = offset_[c] + count_[c];

        cursor_.assign(offset_.begin(), offset_.begin() + static_cast<std::ptrdiff_t>(nc));
        for (std::size_t i = 0; i < n_; ++i)
        {
            const std::size_t c = cell_of(i);
            sorted_idx_[cursor_[c]++] = i;
        }
    }

    std::size_t size() const { return n_; }
    double cell_size() const { return cell_size_; }

    // "Свой тип": пары (i, j), i < j, |x_i - x_j| <= radius, accept(i, j).
    // accept — дополнительный фильтр (например, object_id[i] != object_id[j]
    // для Fiber-Fiber-Contact); по умолчанию пропускает все пары.
    template <class Accept = bool (*)(index_t, index_t)>
    std::vector<std::pair<index_t, index_t>> find_pairs(
        double radius, Accept accept = [](index_t, index_t) { return true; }) const
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
        Coords other_coords, std::size_t n_other, double radius,
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
    //
    // Варианты с выходным параметром out переиспользуют буферы worker'ов
    // (члены индекса) и память out — вызов не обращается к аллокатору, если
    // емкости хватает; порядок пар тот же, что у вариантов по значению.
    // Они не потокобезопасны относительно друг друга: один индекс — один
    // параллельный запрос одновременно.
    template <class Accept = bool (*)(index_t, index_t)>
    void find_pairs_parallel(const CpuExecutor& exec, double radius, Accept accept,
                              std::vector<std::pair<index_t, index_t>>& out) const
    {
        const std::size_t nc = total_cells();
        parallel_merge(exec, nc, out, [&](std::size_t c, std::vector<std::pair<index_t, index_t>>& local) {
            collect_same_type_cell(c, radius, accept, local);
        });
    }

    template <class Accept = bool (*)(index_t, index_t)>
    std::vector<std::pair<index_t, index_t>> find_pairs_parallel(
        const CpuExecutor& exec, double radius,
        Accept accept = [](index_t, index_t) { return true; }) const
    {
        std::vector<std::pair<index_t, index_t>> out;
        find_pairs_parallel(exec, radius, accept, out);
        return out;
    }

    template <class Accept = bool (*)(index_t, index_t)>
    void find_pairs_cross_parallel(const CpuExecutor& exec, Coords other_coords, std::size_t n_other, double radius,
                                    Accept accept, std::vector<std::pair<index_t, index_t>>& out) const
    {
        parallel_merge(exec, n_other, out, [&](std::size_t a, std::vector<std::pair<index_t, index_t>>& local) {
            collect_cross_one(other_coords, a, radius, accept, local);
        });
    }

    template <class Accept = bool (*)(index_t, index_t)>
    std::vector<std::pair<index_t, index_t>> find_pairs_cross_parallel(
        const CpuExecutor& exec, Coords other_coords, std::size_t n_other, double radius,
        Accept accept = [](index_t, index_t) { return true; }) const
    {
        std::vector<std::pair<index_t, index_t>> out;
        find_pairs_cross_parallel(exec, other_coords, n_other, radius, accept, out);
        return out;
    }

private:
    std::array<std::vector<double>, Dim> coord_;
    std::size_t n_ = 0;
    double cell_size_ = 1.0;
    std::array<double, Dim> origin_{};
    std::array<int, Dim> ncell_{};
    std::vector<std::size_t> offset_;
    std::vector<std::size_t> sorted_idx_;
    // Рабочая память build() и параллельных запросов: переиспользуется между
    // вызовами, только растёт (индекс живёт между шагами).
    std::vector<std::size_t> count_;
    std::vector<std::size_t> cursor_;
    // alignas(64): заголовок вектора каждого потока -- в своей кэш-линии (false sharing).
    struct alignas(64) LocalPairs { std::vector<std::pair<index_t, index_t>> v; };
    mutable std::vector<LocalPairs> local_pairs_;

    std::size_t total_cells() const
    {
        std::size_t nc = 1;
        for (std::size_t d = 0; d < Dim; ++d) nc *= static_cast<std::size_t>(ncell_[d]);
        return nc;
    }

    std::array<int, Dim> cell_of_point(const std::array<double, Dim>& p) const
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
        std::array<double, Dim> p{};
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

    double dist2(std::size_t i, std::size_t j) const
    {
        double d2 = 0.0;
        for (std::size_t d = 0; d < Dim; ++d)
        {
            const double dd = coord_[d][i] - coord_[d][j];
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
    int neighbor_reach(double radius) const
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
    void collect_same_type_cell(std::size_t cell_lin, double radius, const Accept& accept,
                                 std::vector<std::pair<index_t, index_t>>& out) const
    {
        const double r2 = radius * radius;
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
    void collect_cross_one(Coords other_coords, std::size_t a, double radius, const Accept& accept,
                            std::vector<std::pair<index_t, index_t>>& out) const
    {
        const double r2 = radius * radius;
        std::array<double, Dim> p{};
        for (std::size_t d = 0; d < Dim; ++d) p[d] = other_coords[d][a];
        const std::array<int, Dim> cell = cell_of_point(p);
        const int reach = neighbor_reach(radius);

        for_each_neighbor_cell(cell, reach, [&](std::size_t nb_lin) {
            auto [nb, ne] = cell_range(nb_lin);
            for (std::size_t pb = nb; pb < ne; ++pb)
            {
                const std::size_t b = sorted_idx_[pb];
                double d2 = 0.0;
                for (std::size_t d = 0; d < Dim; ++d)
                {
                    const double dd = p[d] - coord_[d][b];
                    d2 += dd * dd;
                }
                if (d2 <= r2 && accept(a, b)) out.emplace_back(a, b);
            }
        });
    }

    // Буферы worker'ов — mutable-член local_pairs_ (clear() сохраняет память);
    // результат в порядке отрезков дописывается в out, предварительно очищенный.
    template <class Body>
    void parallel_merge(const CpuExecutor& exec, std::size_t range_n,
                         std::vector<std::pair<index_t, index_t>>& out, Body&& body) const
    {
        // Детерминированно: отрезки диапазона фиксированы (не зависят от числа потоков), слияние -- в порядке отрезков.
        const std::size_t nw = exec.det_segments();
        if (local_pairs_.size() < nw) local_pairs_.resize(nw);
        for (auto& lp : local_pairs_) lp.v.clear();
        exec.parallel_for_segments(range_n, [&](std::size_t seg, std::size_t idx) { body(idx, local_pairs_[seg].v); });
        std::size_t total = 0;
        for (auto& lp : local_pairs_) total += lp.v.size();
        out.clear();
        out.reserve(total);
        for (auto& lp : local_pairs_) out.insert(out.end(), lp.v.begin(), lp.v.end());
    }
};

} // namespace rhs

#endif // RHS_SPATIAL_INDEX_HPP
