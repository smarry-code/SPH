#pragma once
// bc2_geometry.hpp
// Построение геометрии кейса BC2 (удар по карбиду бора).
//
// Вынесено из tests/test_e2e_bc2.cpp (раздел «1. Геометрия кейса BC2») без
// изменений логики: постановку теперь использует и регрессионный тест
// tests/test_regression.cpp, и двух расходящихся копий быть не должно.
// Раскладка BC2 отличается от BC6 (tests/bc6_geometry.hpp): узлы по центрам
// ячеек, шаг подгоняется под толщину слоя отдельно по каждой оси.
//
#include <array>
#include <cmath>
#include <cstddef>

#include "param.hpp"
#include "particle_system.hpp"

namespace bc2
{

// Перенесена из закомментированного блока src/input.cpp:490-596 (внутри
// shock_tube). Здесь — ЕДИНСТВЕННЫЙ исполняемый экземпляр этой постановки;
// src/input.cpp намеренно не изменяется (ТЗ §1.4). Если постановка станет
// частью продового кода — переносить отсюда, а блок в input.cpp удалить.
//
// Исправления дефектов оригинала (ТЗ §1.3):
//   - добавлен двухпроходный подсчёт + sys.resize() ДО записи (в оригинале
//     запись в sys.x[0][ntotal] шла в контейнер нулевого размера — UB);
//   - накопление `r += delta_x` / `z += delta_x` в условии цикла заменено на
//     целочисленные счётчики: число частиц больше не зависит от накопленной
//     ошибки округления double на границе диапазона;
//   - `moni_particle = 12973` (глобал, которого больше нет) заменён на
//     геометрическое правило (см. pick_moni_particle);
//   - sys.axis заполняется после resize (std::vector<char>, bool-семантика).
//
// Приведение постановки к статье-источнику эталона (Grady & Moody,
// SAND96-0551, опыт BC2 = оригинальный тест CE6): узлы по центрам ячеек со
// своим шагом в каждом слое, контактная сборка без зазоров, монитор на оси,
// плотности и упругие константы B4C марки Eagle Picher. Подробности — в
// комментариях к bc2_layers, cell_count и pick_moni_particle.

struct LayerSpec
{
    double r_hi;      // верхняя граница r, включительно
    double z_lo;      // нижняя граница z, включительно
    double z_hi;      // верхняя граница z, включительно
    int    itype;
    double rho;
    double vz;        // vx[1], м/с
    const char* name;
};

struct Bc2Geometry
{
    size_t ntotal = 0;
    size_t moni_particle = 0;
    double delta_x = 0.0;
    std::array<size_t, 4> layer_begin{};
    std::array<size_t, 4> layer_count{};
    std::array<size_t, 4> layer_nr{};
    std::array<size_t, 4> layer_nz{};
    // Фактические шаги решётки внутри слоя. Отличаются от delta_x, потому что
    // слой разбивается на ЦЕЛОЕ число ячеек, а его номинальная толщина в общем
    // случае не кратна delta_x (слой 2: 3.917 мм при delta_x = 0.2 мм).
    std::array<double, 4> layer_dr{};
    std::array<double, 4> layer_dz{};
};

// Разбиение отрезка длиной L на целое число ячеек, ближайшее к L/dx.
// Узлы ставятся по ЦЕНТРАМ ячеек: lo + (i + 0.5)*d, i = 0..n-1, d = L/n.
// Зачем именно центры ячеек:
//   - слой воспроизводит свою НОМИНАЛЬНУЮ толщину точно (узловая схема
//     lo + i*dx отбрасывала хвост короче dx: слой 2 при dx = 0.2 мм давал
//     3.600 мм вместо 3.917 мм — на 8 % короче, то есть короче был и
//     нагружающий импульс);
//   - соседние слои получаются в КОНТАКТЕ без искусственных зазоров:
//     последний узел слоя и первый узел следующего отстоят ровно на шаг
//     решётки, то есть ровно так же, как два соседних узла внутри слоя.
inline size_t cell_count(double L, double dx)
{
    if (!(L > 0.0) || !(dx > 0.0)) return 1;
    const long n = std::lround(L / dx);
    return n < 1 ? 1u : static_cast<size_t>(n);
}

inline std::array<LayerSpec, 4> bc2_layers()
{
    // Опорные плоскости. Толщины и радиусы подтверждены статьёй-источником
    // эталона: Grady & Moody, SAND96-0551, «Shock Compression Profiles in
    // Ceramics» (1996), опыт BC2 (в статье — тест CE6).
    //
    // ВНИМАНИЕ: постановка ниже СОЗНАТЕЛЬНО расходится с закомментированным
    // блоком в src/input.cpp:490-596, из которого она когда-то была перенесена.
    // Блок в input.cpp не изменяется (ТЗ §1.4) и остаётся историческим
    // свидетельством, а не источником истины. Расхождения и их основание:
    //   - зазоры между слоями (+dx, +4dx, +5dx) УБРАНЫ: по статье сборка
    //     контактная, зазоров не задано, а прежние зазоры к тому же зависели
    //     от разрешения сетки (0.2 мм при dx = 0.2 мм против 0.8 мм при
    //     dx = 0.8 мм) — то есть постановка менялась вместе с сеткой;
    //   - плотности приведены к статье: подложка 626 -> 640, B4C 2510 -> 2517,
    //     LiF 2638 -> 2640 кг/м^3;
    //   - летящая пластина B4C даёт ровно 3.917 мм (см. cell_count), тогда как
    //     узловая схема оригинала обрезала её до 3.600 мм.
    const double Z1 = 0.006;
    const double Z2 = Z1 + 0.003917;   // 0.009917
    const double Z3 = Z2 + 0.009033;   // 0.018950
    const double Z4 = Z3 + 0.0254;     // 0.044350
    return std::array<LayerSpec, 4>{{
        // r_hi,     z_lo,  z_hi,  itype, rho,     vz,     name
        {0.04375,  0.0,   Z1,      8,    640.0, 2210.0, "ударник (полиуретан 626)"},
        {0.0346,   Z1,    Z2,     10,   2517.0, 2210.0, "керамика летящая (B4C)"},
        {0.0346,   Z2,    Z3,     10,   2517.0,    0.0, "керамика покоящаяся (B4C)"},
        {0.0254,   Z3,    Z4,      5,   2640.0,    0.0, "окно (LiF)"},
    }};
}

// Точка съёма эталонной кривой.
// Статья (SAND96-0551, BC2 = CE6) фиксирует место измерения: «at the center of
// the ceramic/window interface», то есть НА ОСИ. Поэтому берётся внутренний
// столбец решётки, ir = 0 (r = 0.5*dr), и последний узел по z в слое 3 —
// тыльная поверхность покоящейся керамики, граница с окном.
//
// Прежнее правило искало столбец, ближайший к r = 12.5 мм: так был устроен
// литерал moni_particle = 12973 из оригинального блока input.cpp. Это не
// «центр интерфейса», а точка на половине радиуса; статья такой точки не
// упоминает. Ось — самая шумная зона осесимметричного SPH, но выбор диктует
// постановка эксперимента, а не удобство численной кривой.
inline size_t pick_moni_particle(const Bc2Geometry& g)
{
    const size_t n_z = g.layer_nz[2];
    return g.layer_begin[2] + (n_z - 1);
}

inline Bc2Geometry build_bc2_geometry(ParticleSystem& sys, double delta_x)
{
    const auto layers = bc2_layers();

    Bc2Geometry g;
    g.delta_x = delta_x;

    // --- проход 1: только подсчёт (даёт размер ДО любой записи) -------------
    size_t total = 0;
    for (size_t L = 0; L < layers.size(); ++L)
    {
        const size_t n_r = cell_count(layers[L].r_hi, delta_x);
        const size_t n_z = cell_count(layers[L].z_hi - layers[L].z_lo, delta_x);
        g.layer_nr[L] = n_r;
        g.layer_nz[L] = n_z;
        g.layer_dr[L] = layers[L].r_hi / static_cast<double>(n_r);
        g.layer_dz[L] = (layers[L].z_hi - layers[L].z_lo) / static_cast<double>(n_z);
        g.layer_begin[L] = total;
        g.layer_count[L] = n_r * n_z;
        total += n_r * n_z;
    }
    g.ntotal = total;

    // --- resize ДО заполнения (исправление дефекта §1.3) --------------------
    sys.resize(total);

    // --- проход 2: заполнение ----------------------------------------------
    size_t k = 0;
    for (size_t L = 0; L < layers.size(); ++L)
    {
        const LayerSpec& s = layers[L];
        const double dr = g.layer_dr[L];
        const double dz = g.layer_dz[L];
        for (size_t ir = 0; ir < g.layer_nr[L]; ++ir)
        {
            const double r = (0.5 + static_cast<double>(ir)) * dr;
            for (size_t iz = 0; iz < g.layer_nz[L]; ++iz)
            {
                const double z = s.z_lo + (0.5 + static_cast<double>(iz)) * dz;
                sys.axis[k]  = (std::fabs(r) <= 0.5 * dr) ? 1 : 0;
                sys.hsml[k]  = 2.0 * delta_x;
                sys.u[k]     = 0.0;
                sys.p[k]     = 0.0;
                sys.x[0][k]  = r;
                sys.x[1][k]  = z;
                sys.itype[k] = s.itype;
                sys.rho[k]   = s.rho;
                sys.rho_0[k] = s.rho;
                sys.vx[0][k] = 0.0;
                sys.vx[1][k] = s.vz;
                // Плоская формула массы — из раскомментированной строки
                // оригинала. Два осесимметричных варианта там закомментированы
                // и намеренно не переносятся (менять постановку — не задача теста).
                // Шаги берутся фактические, по слою: при dr != dz != delta_x
                // произведение delta_x*delta_x дало бы массу, не отвечающую
                // занимаемой частицей ячейке.
                sys.mass[k]  = s.rho * dr * dz;
                ++k;
            }
        }
    }

    g.moni_particle = pick_moni_particle(g);
    return g;
}

} // namespace bc2
