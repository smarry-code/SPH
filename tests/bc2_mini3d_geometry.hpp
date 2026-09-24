#pragma once
// bc2_mini3d_geometry.hpp
// Уменьшенный кейс BC2 в полном 3D — для сквозного и регрессионного
// 3D-тестов (tests/test_e2e_3d.cpp, tests/test_regression.cpp при dim = 3).
//
// ЗАЧЕМ УМЕНЬШЕННЫЙ. Полный BC2 в 3D даже на сетке smoke (Δx = 0.8 мм) —
// это сотни тысяч частиц: сечение ударника радиусом 43.75 мм несёт ~9400
// узлов на слой. Поэтому толщины слоёв, материалы, плотности и скорости
// взяты из BC2 как есть (bc2::bc2_layers()), а радиусы всех слоёв уменьшены в
// одном отношении так, чтобы керамика имела R = 4 мм:
//     ударник 5.06 мм, керамика 4.00 мм, окно 2.94 мм.
// Это уже не опыт BC2, а его «тонкая» модель: боковая разгрузка доходит до
// оси почти сразу, и с приборным эталоном сверять нечего. Кейс проверяет
// ветку dim = 3 конвейера на тех же материалах, УРС и моделях прочности.
//
// РЕШЁТКА. Квадратная в плоскости (x, y): узлы x = (i + ½)·Δx, i ∈ [−n, n),
// то есть ось z проходит МЕЖДУ четырьмя центральными столбцами. Частица
// берётся, если x² + y² ≤ R². Такая решётка вместе с кругом переходит в себя
// при отражениях x → −x, y → −y и перестановке x ↔ y (группа квадрата), и
// тест пользуется этим: расчёт обязан сохранять ту же симметрию. Для этого
// геометрия сразу строит отображения частица -> её образ.
//
// По z — как в BC2 (bc2_geometry.hpp): слой режется на целое число ячеек,
// ближайшее к T/Δx, узлы по центрам ячеек, слои в контакте без зазоров.
// Масса — ρ·Δx·Δx·dz, длина сглаживания 2Δx.

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "param.hpp"
#include "particle_system.hpp"

#include "bc2_geometry.hpp"

namespace bc2m3
{

// Радиус керамики в уменьшенной модели и общий коэффициент для всех слоёв.
inline constexpr double kCeramicRadius = 4.0e-3;
inline constexpr double kBc2CeramicRadius = 0.0346;   // bc2_layers()[1].r_hi
inline constexpr double kRadiusScale = kCeramicRadius / kBc2CeramicRadius;

inline constexpr double kDxSmoke = 8.0e-4;            // как у smoke BC2

struct Geometry
{
    size_t ntotal = 0;
    double delta_x = 0.0;
    std::array<size_t, 4> layer_begin{};
    std::array<size_t, 4> layer_count{};
    std::array<size_t, 4> layer_columns{};    // столбцов (x, y) в сечении слоя
    std::array<size_t, 4> layer_nz{};
    std::array<double, 4> layer_radius{};
    std::array<double, 4> layer_dz{};

    // Образ частицы при симметриях решётки: сквозные индексы.
    std::vector<size_t> mirror_x;   // x → −x
    std::vector<size_t> mirror_y;   // y → −y
    std::vector<size_t> swap_xy;    // x ↔ y

    // Монитор: тыльная грань покоящейся керамики (граница с окном), один из
    // четырёх приосевых столбцов (x = y = +Δx/2). Там же, где у BC2 — на оси.
    size_t moni = 0;
};

inline Geometry build_geometry(ParticleSystem& sys, double delta_x)
{
    static_assert(dim == 3, "bc2_mini3d_geometry.hpp требует dim = 3");
    const auto layers = bc2::bc2_layers();

    Geometry g;
    g.delta_x = delta_x;

    // Столбцы сечения каждого слоя: (ix, iy) и их порядковый номер.
    std::array<std::vector<std::pair<int, int>>, 4> cols;
    std::array<std::map<std::pair<int, int>, size_t>, 4> col_index;
    size_t total = 0;
    for (size_t L = 0; L < layers.size(); ++L)
    {
        const double R = layers[L].r_hi * kRadiusScale;
        const int n = static_cast<int>(std::ceil(R / delta_x)) + 1;
        for (int ix = -n; ix < n; ++ix)
            for (int iy = -n; iy < n; ++iy)
            {
                const double x = (ix + 0.5) * delta_x;
                const double y = (iy + 0.5) * delta_x;
                if (x * x + y * y <= R * R)
                {
                    col_index[L][{ix, iy}] = cols[L].size();
                    cols[L].emplace_back(ix, iy);
                }
            }
        const double T = layers[L].z_hi - layers[L].z_lo;
        const size_t nz = bc2::cell_count(T, delta_x);
        g.layer_radius[L] = R;
        g.layer_columns[L] = cols[L].size();
        g.layer_nz[L] = nz;
        g.layer_dz[L] = T / static_cast<double>(nz);
        g.layer_begin[L] = total;
        g.layer_count[L] = nz * cols[L].size();
        total += g.layer_count[L];
    }
    g.ntotal = total;
    sys.resize(total);
    g.mirror_x.resize(total);
    g.mirror_y.resize(total);
    g.swap_xy.resize(total);

    // Сквозной индекс: начало слоя + iz·(столбцов) + номер столбца.
    auto index = [&](size_t L, size_t iz, int ix, int iy) {
        const auto it = col_index[L].find({ix, iy});
        if (it == col_index[L].end())
        {
            // Решётка и круг симметричны — образ существует всегда. Если нет,
            // сломано построение, и симметрийная проверка была бы ложной.
            throw std::logic_error("bc2_mini3d: у столбца нет симметричного образа");
        }
        return g.layer_begin[L] + iz * cols[L].size() + it->second;
    };

    for (size_t L = 0; L < layers.size(); ++L)
    {
        const bc2::LayerSpec& s = layers[L];
        const double dz = g.layer_dz[L];
        for (size_t iz = 0; iz < g.layer_nz[L]; ++iz)
        {
            const double z = s.z_lo + (0.5 + static_cast<double>(iz)) * dz;
            for (size_t c = 0; c < cols[L].size(); ++c)
            {
                const auto [ix, iy] = cols[L][c];
                const size_t k = g.layer_begin[L] + iz * cols[L].size() + c;
                sys.inzializ_index(k);
                sys.x[0][k]  = (ix + 0.5) * delta_x;
                sys.x[1][k]  = (iy + 0.5) * delta_x;
                sys.x[2][k]  = z;
                sys.vx[2][k] = s.vz;
                sys.itype[k] = s.itype;
                sys.rho[k]   = s.rho;
                sys.rho_0[k] = s.rho;
                sys.hsml[k]  = 2.0 * delta_x;
                sys.mass[k]  = s.rho * delta_x * delta_x * dz;
                g.mirror_x[k] = index(L, iz, -1 - ix, iy);
                g.mirror_y[k] = index(L, iz, ix, -1 - iy);
                g.swap_xy[k]  = index(L, iz, iy, ix);
            }
        }
    }

    // Тыльный слой покоящейся керамики (слой 2), столбец (0, 0).
    g.moni = index(2, g.layer_nz[2] - 1, 0, 0);
    return g;
}

} // namespace bc2m3
