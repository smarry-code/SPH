#ifndef FLAT_GRID_HPP
#define FLAT_GRID_HPP

#include <vector>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cassert>
#include "param.hpp"
#include "particle_system.hpp"

// ============================================================================
// FlatGrid — плотная сетка с прямым индексированием ячеек через counting sort
// Замена std::unordered_map<CellKey, std::vector<size_t>> для ускорения
// построения ячеек и парного поиска соседей.
// ============================================================================

struct FlatGrid
{
    int nx, ny, nz;   // размеры сетки по каждой оси
    double dx;        // размер ячейки (= scale_k * max_hmsl)
    double origin[3]; // сдвиг (xmin, ymin, zmin)

    std::vector<size_t> count;      // [n_cells] — число частиц в ячейке
    std::vector<size_t> offset;     // [n_cells+1] — префиксная сумма
    std::vector<size_t> sorted_idx; // [N] — индексы частиц, сгруппированные по ячейкам

    // Полное число ячеек
    int n_cells() const { return nx * ny * nz; }

    // -----------------------------------------------------------------------
    // cell_index — вычисление линейного индекса ячейки для точки (x, y, z)
    // Со сдвигом на origin и clamp в границы сетки.
    // Для dim=1: y и z игнорируются (iy=0, iz=0).
    // Для dim=2: z игнорируется (iz=0).
    // -----------------------------------------------------------------------
    int cell_index(double x, double y, double z) const
    {
        int ix = static_cast<int>(std::floor((x - origin[0]) / dx));
        int iy = (dim >= 2) ? static_cast<int>(std::floor((y - origin[1]) / dx)) : 0;
        int iz = (dim >= 3) ? static_cast<int>(std::floor((z - origin[2]) / dx)) : 0;

        // Clamp в границы сетки (защита от вылетевших/граничных частиц)
        ix = std::clamp(ix, 0, nx - 1);
        iy = std::clamp(iy, 0, ny - 1);
        iz = std::clamp(iz, 0, nz - 1);

        return ix + nx * (iy + ny * iz);
    }

    // -----------------------------------------------------------------------
    // cell_range — диапазон индексов в sorted_idx для ячейки c
    // Возвращает [begin, end): частицы ячейки c лежат в sorted_idx[begin..end)
    // -----------------------------------------------------------------------
    std::pair<size_t, size_t> cell_range(int c) const
    {
        const size_t cell = static_cast<size_t>(c);
        return {offset[cell], offset[cell + 1]};
    }

    // -----------------------------------------------------------------------
    // build — построение сетки за 3 прохода (counting sort)
    //
    // Pass 1 (Bounds):   поиск min/max координат, вычисление origin, dx, nx/ny/nz
    // Pass 2 (Count):    подсчёт числа частиц в каждой ячейке
    // Pass 3 (Prefix+Scatter): префиксная сумма + раскладка индексов по ячейкам
    //
    // Параметры:
    //   Elem    — массив частиц (вектор Particle)
    //   ntotal  — число частиц
    //   scale_k — масштабный коэффициент ядра (2 для куб. сплайна, 3 для Гаусса)
    // -----------------------------------------------------------------------
    void build(const ParticleSystem& sys, size_t ntotal, int scale_k)
    {
        if (ntotal == 0)
        {
            nx = ny = nz = 1;
            dx = 1.0;
            origin[0] = origin[1] = origin[2] = 0.0;
            count.assign(1, 0);
            offset.assign(2, 0);
            sorted_idx.clear();
            return;
        }

        // ==== Pass 1: Bounds ====
        // Поиск min/max координат и max_hsml
        double xmin = sys.x[0][0], xmax = sys.x[0][0];
        double ymin = 0.0, ymax = 0.0;
        double zmin = 0.0, zmax = 0.0;
        double max_hmsl = sys.hsml[0];

        for (size_t i = 1; i < ntotal; ++i)
        {
            if (sys.x[0][i] < xmin)
                xmin = sys.x[0][i];
            if (sys.x[0][i] > xmax)
                xmax = sys.x[0][i];
            if (sys.hsml[i] > max_hmsl)
                max_hmsl = sys.hsml[i];

            if (dim >= 2)
            {
                if (i == 1)
                {
                    ymin = sys.x[1][i];
                    ymax = sys.x[1][i];
                }
                if (sys.x[1][i] < ymin)
                    ymin = sys.x[1][i];
                if (sys.x[1][i] > ymax)
                    ymax = sys.x[1][i];
            }
            if (dim >= 3)
            {
                if (i == 1)
                {
                    zmin = sys.x[2][i];
                    zmax = sys.x[2][i];
                }
                if (sys.x[2][i] < zmin)
                    zmin = sys.x[2][i];
                if (sys.x[2][i] > zmax)
                    zmax = sys.x[2][i];
            }
        }

        // Также учтём первую частицу для y/z
        if (dim >= 2)
        {
            if (sys.x[1][0] < ymin)
                ymin = sys.x[1][0];
            if (sys.x[1][0] > ymax)
                ymax = sys.x[1][0];
        }
        if (dim >= 3)
        {
            if (sys.x[2][0] < zmin)
                zmin = sys.x[2][0];
            if (sys.x[2][0] > zmax)
                zmax = sys.x[2][0];
        }

        // Размер ячейки (сохраняем формулу текущего кода)
        dx = scale_k * max_hmsl;
        if (dx <= 0.0)
            dx = 1.0; // защита от нулевого dx

        // Начало сетки = минимальные координаты
        origin[0] = xmin;
        origin[1] = (dim >= 2) ? ymin : 0.0;
        origin[2] = (dim >= 3) ? zmin : 0.0;

        // Число ячеек по каждой оси (+1 для гарантии попадания крайней частицы)
        nx = static_cast<int>(std::ceil((xmax - xmin) / dx)) + 1;
        ny = (dim >= 2) ? static_cast<int>(std::ceil((ymax - ymin) / dx)) + 1 : 1;
        nz = (dim >= 3) ? static_cast<int>(std::ceil((zmax - zmin) / dx)) + 1 : 1;

        const size_t nc = static_cast<size_t>(n_cells());

        // Выделение памяти
        count.assign(nc, 0);
        offset.resize(nc + 1, 0);
        sorted_idx.resize(ntotal);

        // ==== Pass 2: Count ====
        // Подсчёт числа частиц в каждой ячейке
        for (size_t i = 0; i < ntotal; ++i)
        {
            int c = cell_index(sys.x[0][i],
                               (dim >= 2) ? sys.x[1][i] : 0.0,
                               (dim >= 3) ? sys.x[2][i] : 0.0);
            count[static_cast<size_t>(c)]++;
        }

        // Контрольный инвариант: сумма count == N
        size_t total_counted = 0;
        for (size_t c = 0; c < nc; ++c)
        {
            total_counted += count[c];
        }
        assert(total_counted == ntotal &&
               "FlatGrid: сумма count[] не равна ntotal!");

        // ==== Pass 3: Prefix sum + Scatter ====
        // Префиксная сумма: offset[c] = начало ячейки c в sorted_idx
        offset[0] = 0;
        for (size_t c = 0; c < nc; ++c)
        {
            offset[c + 1] = offset[c] + count[c];
        }
        // offset[nc] должен равняться ntotal
        assert(offset[nc] == ntotal &&
               "FlatGrid: offset[n_cells] != ntotal!");

        // Копия offset для scatter (pos[c] — текущая позиция записи в ячейке c)
        std::vector<size_t> pos(offset.begin(), offset.begin() + static_cast<std::ptrdiff_t>(nc));

        // Раскладка индексов частиц по ячейкам
        for (size_t i = 0; i < ntotal; ++i)
        {
            int c = cell_index(sys.x[0][i],
                               (dim >= 2) ? sys.x[1][i] : 0.0,
                               (dim >= 3) ? sys.x[2][i] : 0.0);
            sorted_idx[pos[static_cast<size_t>(c)]++] = i;
        }
    }
};

#endif // FLAT_GRID_HPP
