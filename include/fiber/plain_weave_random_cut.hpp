// plain_weave_random_cut.hpp
// Ткань «случайная раскройка» (fiber/tz-random-cut-quilt.md) поверх генератора
// полотняного переплетения (include/fiber/plain_weave.hpp).
//
// Header-only. НЕ зависит от dim: только plain_weave.hpp и стандартная библиотека.
//
// ИДЕЯ. Прямоугольник W x H (центр (cx, cy)) замащивается диаграммой Вороного из N
// ячеек (одна итерация Ллойда). Каждая ячейка -- отдельный фрагмент полотняной ткани
// со своим поворотом сетки theta и сдвигом фазы (fx, fy). Рабочий многоугольник
// фрагмента -- ячейка, у которой биссекторы сдвинуты внутрь на gap/2 (границы
// прямоугольника не сдвигаются). Общих нитей у фрагментов нет: нити обрезаны по
// рабочему многоугольнику, фрагменты взаимодействуют контактом (его ищет ядро) и
// прошиты стежками BOND материала 1 через швы.
//
// МНОГОСЛОЙНЫЙ ПАКЕТ. L слоёв (layers), слой l на cz_l = cz - l*ply_pitch. Слои
// группы g = l / k (k = layers_per_cut) делят одну раскройку: одинаковые
// многоугольники, theta и фазу (копии со сдвигом по z). Раскройки групп независимы
// (один поток ГСЧ, расходуется по группам подряд, поэтому L = 1 даёт тот же
// результат, что однослойный режим). Все слои сшиты насквозь: точки стежков на
// швах ВСЕХ раскроек пакета, в каждой точке для каждой пары соседних слоёв (l, l+1)
// -- BOND между ближайшими в плане к точке бусинами слоя l и слоя l+1.
//
// ГЕОМЕТРИЯ ФРАГМЕНТА (центр c, u = (cos th, sin th), v = (-sin th, cos th)).
//   Основа j: c + (j + fx) p v + t u, вершины (пересечения с утком) на t = (m + fy) p.
//   Уток   i: c + (i + fy) p u + t v, вершины на t = (m + fx) p.
//   Вдоль нити локальная координата a = t - (фаза вдоль) * p, вершина w на a = w p.
//   z: основа cz + (-1)^j * zc(a), уток cz - (-1)^i * zc(a), zc -- профиль Curve
//   (Cos или Arcs) со знаком (-1)^w на полуволне w; в пересечении (i, j) основа на
//   (-1)^(i+j) A, уток на -(-1)^(i+j) A: разнесены на 2A.
//   Внутренние бусины -- вершины и по seg-1 бусин равных хорд между ними (полуволна
//   equal_chords считается один раз); концевые бусины -- точно на границе обрезки,
//   если до соседней бусины не меньше 0.2 h (h = p/seg).
//
// ПОРЯДОК НИТЕЙ (object_id сквозной): слой -> фрагмент -> основа (j по возрастанию)
// -> уток (i по возрастанию). Глобальный индекс бусины -- порядок записи блоков FIBER.
#pragma once

#include "fiber/plain_weave.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <ostream>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fiber
{

struct RandomCutParams
{
    PlainWeaveParams base;       // ткань: d, pitch, amp, seg, profile, cx/cy/cz, материал, clamp; nx/ny/overhang игнорируются
    int n_frag = 1;              // N: число фрагментов
    double size_x = 0.05;        // W, м
    double size_y = 0.05;        // H, м
    std::uint64_t seed = 1;      // зерно ГСЧ
    double gap = -1.0;           // зазор между фрагментами, м; < 0 -- по умолчанию d
    double rot_max_deg = 25.0;   // поворот сетки фрагмента равномерно в [-rot_max, +rot_max], градусы
    double stitch_step = 2.4e-3; // шаг стежков вдоль шва, м; 0 -- без стежков
    // Материал 1 стежка. ЗАГЛУШКИ: значения задаст заказчица (ТЗ §9 п.1).
    double stitch_E = 5.0e9;
    double stitch_nu = 0.3;
    double stitch_d = 0.2e-3;
    double stitch_rho = 1000.0;
    double stitch_eps_break = 0.2;
    double rc = 0.9e-3;          // радиус контакта только для диагностики (.cut.json)
    int layers = 1;              // L: число слоёв пакета
    int layers_per_cut = 0;      // k: слоёв на одну раскройку; 0 -- по умолчанию L (одна раскройка)
    double ply_pitch = -1.0;     // шаг слоёв по z, м; < 0 -- по умолчанию 2A + d + 1 мкм
};

// Вершины многоугольника и полуплоскости n.x <= c.
using RcPoint = std::array<double, 2>;
struct RcHalfPlane
{
    double nx, ny, c;
};

struct RandomCutFragment
{
    int id = 0;
    RcPoint center{{0.0, 0.0}};   // центр ячейки (после итерации Ллойда), начало решётки фрагмента
    double theta = 0.0;           // поворот сетки, рад
    double fx = 0.0, fy = 0.0;    // сдвиг фазы сетки (доли p)
    std::vector<RcPoint> cell;    // ячейка Вороного без зазора
    std::vector<RcPoint> polygon; // рабочий многоугольник (биссекторы сдвинуты на gap/2)
    std::vector<RcHalfPlane> hps; // полуплоскости рабочего многоугольника (прямоугольник + биссекторы)
    double area = 0.0;            // площадь рабочего многоугольника
};

struct RandomCutSeam
{
    int i = 0, j = 0;            // фрагменты, i < j
    RcPoint a{{0.0, 0.0}}, b{{0.0, 0.0}}; // концы отрезка шва (общая граница ячеек без зазора)
    double length = 0.0;
    std::vector<RcPoint> points; // точки стежков (середины ns равных частей)
};

// Одна раскройка (общая для группы слоёв).
struct RandomCutLayout
{
    std::vector<RandomCutFragment> frags;
    std::vector<RandomCutSeam> seams;
};

struct RandomCutStitch
{
    std::size_t a = 0, b = 0; // глобальные индексы бусин
    bool through = false;     // false -- стежок через шов в слое, true -- сквозной между слоями ply_a и ply_b = ply_a + 1
    int ply_a = 0, ply_b = 0;
    int frag_a = 0, frag_b = 0;
    int cut = 0;              // раскройка, которой принадлежит шов
    int seam = 0;             // номер шва в cuts[cut].seams
};

// Диапазон нитей и бусин фрагмента frag слоя ply (глобальные индексы).
struct RandomCutRange
{
    int ply = 0, frag = 0;
    std::size_t fiber_first = 0, fiber_count = 0;
    std::size_t bead_first = 0, bead_count = 0;
};

struct RandomCut
{
    RandomCutParams params;        // с подставленными значениями по умолчанию (gap, k, ply_pitch)
    PlainWeave W;                  // все нити всех слоёв; W.params.nx/ny = число нитей основы/утка
    std::vector<int> fiber_ply;    // [нить] -> слой
    std::vector<int> fiber_fragment; // [нить] -> фрагмент (внутри раскройки слоя)
    std::vector<RandomCutLayout> cuts;
    std::vector<int> ply_cut;      // [слой] -> раскройка
    std::vector<double> ply_cz;    // [слой] -> c_z
    std::vector<std::size_t> ply_bead_first; // [слой] -> первый глобальный индекс бусины
    std::vector<RandomCutRange> ranges;
    std::vector<RandomCutStitch> stitches;   // сначала в слоях, затем сквозные
    std::vector<std::pair<std::size_t, std::size_t>> stitch_pairs; // те же пары для BOND
    std::size_t n_stitches_in = 0, n_stitches_through = 0;
    std::vector<std::pair<std::size_t, std::size_t>> contacts; // диагностика: разные фрагменты одного слоя ближе rc
    double min_cross_dist = 0.0;   // мин. расстояние в плане между бусинами разных фрагментов одного слоя (inf, если одна)
};

namespace plain_weave_random_cut_detail
{

constexpr int max_fragments = 10000;
constexpr int max_layers = 10000;

// ГСЧ: mt19937_64, вещественное в [0, 1) вручную (не зависит от стандартной библиотеки).
struct Rng
{
    std::mt19937_64 g;
    explicit Rng(std::uint64_t s) : g(s) {}
    double uniform() { return static_cast<double>(g() >> 11) * (1.0 / 9007199254740992.0); }
};

inline double dist2(const RcPoint& a, const RcPoint& b)
{
    return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]);
}

// Прямоугольник [x0, x1] x [y0, y1] как полуплоскости.
inline std::vector<RcHalfPlane> rect_hps(double x0, double y0, double x1, double y1)
{
    return {{1.0, 0.0, x1}, {-1.0, 0.0, -x0}, {0.0, 1.0, y1}, {0.0, -1.0, -y0}};
}

// Полуплоскость «ближе к s_i, чем к s_j», граница сдвинута внутрь на shrink.
inline RcHalfPlane bisector(const RcPoint& si, const RcPoint& sj, double shrink)
{
    const double dx = sj[0] - si[0], dy = sj[1] - si[1];
    const double L = std::hypot(dx, dy);
    const double nx = dx / L, ny = dy / L;
    return {nx, ny, nx * 0.5 * (si[0] + sj[0]) + ny * 0.5 * (si[1] + sj[1]) - shrink};
}

// Sutherland-Hodgman: часть poly, где n.x <= c.
inline std::vector<RcPoint> clip(const std::vector<RcPoint>& poly, const RcHalfPlane& h)
{
    std::vector<RcPoint> out;
    const std::size_t m = poly.size();
    for (std::size_t q = 0; q < m; ++q)
    {
        const RcPoint& a = poly[q];
        const RcPoint& b = poly[(q + 1) % m];
        const double fa = h.nx * a[0] + h.ny * a[1] - h.c;
        const double fb = h.nx * b[0] + h.ny * b[1] - h.c;
        if (fa <= 0.0) out.push_back(a);
        if ((fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0))
        {
            const double t = fa / (fa - fb);
            out.push_back({a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])});
        }
    }
    return out;
}

inline std::vector<RcPoint> clip_all(double x0, double y0, double x1, double y1, const std::vector<RcHalfPlane>& hps)
{
    std::vector<RcPoint> poly = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    for (const auto& h : hps)
    {
        poly = clip(poly, h);
        if (poly.size() < 3) return {};
    }
    return poly;
}

inline double poly_area(const std::vector<RcPoint>& P)
{
    double s = 0.0;
    for (std::size_t q = 0; q < P.size(); ++q)
    {
        const RcPoint& a = P[q];
        const RcPoint& b = P[(q + 1) % P.size()];
        s += a[0] * b[1] - b[0] * a[1];
    }
    return 0.5 * std::fabs(s);
}

inline RcPoint poly_centroid(const std::vector<RcPoint>& P)
{
    double A = 0.0, cx = 0.0, cy = 0.0;
    for (std::size_t q = 0; q < P.size(); ++q)
    {
        const RcPoint& a = P[q];
        const RcPoint& b = P[(q + 1) % P.size()];
        const double cr = a[0] * b[1] - b[0] * a[1];
        A += cr;
        cx += (a[0] + b[0]) * cr;
        cy += (a[1] + b[1]) * cr;
    }
    A *= 0.5;
    return {cx / (6.0 * A), cy / (6.0 * A)};
}

// Интервал t, где p0 + t*dir внутри всех полуплоскостей; false -- пусто.
inline bool line_interval(const RcPoint& p0, const RcPoint& dir, const std::vector<RcHalfPlane>& hps, double& lo,
                          double& hi)
{
    lo = -std::numeric_limits<double>::infinity();
    hi = std::numeric_limits<double>::infinity();
    for (const auto& h : hps)
    {
        const double a = h.nx * dir[0] + h.ny * dir[1];
        const double b = h.c - (h.nx * p0[0] + h.ny * p0[1]);
        if (std::fabs(a) < 1e-14)
        {
            if (b < 0.0) return false;
        }
        else if (a > 0.0)
            hi = std::min(hi, b / a);
        else
            lo = std::max(lo, b / a);
    }
    return lo < hi;
}

// Ячейки Вороного без сдвига (для Ллойда и швов).
inline std::vector<RcHalfPlane> cell_hps(const std::vector<RcPoint>& seeds, int i, double shrink, int skip = -1)
{
    std::vector<RcHalfPlane> out;
    for (int j = 0; j < static_cast<int>(seeds.size()); ++j)
        if (j != i && j != skip) out.push_back(bisector(seeds[static_cast<std::size_t>(i)], seeds[static_cast<std::size_t>(j)], shrink));
    return out;
}

// Плоская сетка ячеек для поиска соседей в плане (x, y).
struct PlanGrid
{
    double cell = 1.0;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> map;
    long long ix_lo = 0, ix_hi = -1, iy_lo = 0, iy_hi = -1;

    static std::uint64_t key(long long ix, long long iy)
    {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(ix)) << 32) |
               static_cast<std::uint64_t>(static_cast<std::uint32_t>(iy));
    }
    long long ix(double x) const { return static_cast<long long>(std::floor(x / cell)); }

    PlanGrid(const std::vector<std::array<double, 3>>& pts, double cell_size) : cell(cell_size)
    {
        for (std::size_t q = 0; q < pts.size(); ++q)
        {
            const long long a = ix(pts[q][0]), b = ix(pts[q][1]);
            if (q == 0)
            {
                ix_lo = ix_hi = a;
                iy_lo = iy_hi = b;
            }
            ix_lo = std::min(ix_lo, a);
            ix_hi = std::max(ix_hi, a);
            iy_lo = std::min(iy_lo, b);
            iy_hi = std::max(iy_hi, b);
            map[key(a, b)].push_back(q);
        }
    }
    const std::vector<std::size_t>* at(long long a, long long b) const
    {
        const auto it = map.find(key(a, b));
        return it == map.end() ? nullptr : &it->second;
    }

    // Ближайшая в плане к pt точка, удовлетворяющая pred; при равенстве -- меньший индекс.
    // Возвращает pts.size()-подобное npos (SIZE_MAX), если таких нет.
    template <class Pred>
    std::size_t nearest(const std::vector<std::array<double, 3>>& pts, const RcPoint& pt, Pred pred) const
    {
        const long long a0 = ix(pt[0]), b0 = ix(pt[1]);
        double best = std::numeric_limits<double>::infinity();
        std::size_t bi = static_cast<std::size_t>(-1);
        const long long rmax = std::max({std::llabs(a0 - ix_lo), std::llabs(a0 - ix_hi), std::llabs(b0 - iy_lo),
                                         std::llabs(b0 - iy_hi)}) + 1;
        auto visit = [&](long long a, long long b) {
            const auto* v = at(a, b);
            if (!v) return;
            for (std::size_t q : *v)
            {
                if (!pred(q)) continue;
                const double d = (pts[q][0] - pt[0]) * (pts[q][0] - pt[0]) + (pts[q][1] - pt[1]) * (pts[q][1] - pt[1]);
                if (d < best || (d == best && q < bi))
                {
                    best = d;
                    bi = q;
                }
            }
        };
        for (long long r = 0; r <= rmax; ++r)
        {
            if (r == 0)
                visit(a0, b0);
            else
            {
                for (long long a = a0 - r; a <= a0 + r; ++a)
                {
                    visit(a, b0 - r);
                    visit(a, b0 + r);
                }
                for (long long b = b0 - r + 1; b <= b0 + r - 1; ++b)
                {
                    visit(a0 - r, b);
                    visit(a0 + r, b);
                }
            }
            // Точки кольца r + 1 не ближе r*cell.
            if (bi != static_cast<std::size_t>(-1) && static_cast<double>(r) * cell * static_cast<double>(r) * cell > best)
                break;
        }
        return bi;
    }
};

// Раскройка одной группы слоёв.
inline RandomCutLayout make_cut(Rng& rng, const RandomCutParams& P, double x0, double y0)
{
    const int N = P.n_frag;
    const double W = P.size_x, H = P.size_y;
    const double x1 = x0 + W, y1 = y0 + H;
    // 1. Центры с минимальным расстоянием.
    std::vector<RcPoint> seeds;
    double min_d = 0.55 * W / std::sqrt(static_cast<double>(N));
    int tries = 0;
    while (static_cast<int>(seeds.size()) < N)
    {
        const double ux = rng.uniform();
        const double uy = rng.uniform();
        const RcPoint s = {x0 + ux * W, y0 + uy * H};
        bool ok = true;
        for (const auto& t : seeds)
            if (dist2(s, t) < min_d * min_d)
            {
                ok = false;
                break;
            }
        if (ok)
        {
            seeds.push_back(s);
            tries = 0;
            continue;
        }
        if (++tries >= 10000)
        {
            min_d *= 0.97;
            tries = 0;
        }
    }
    // 2. Одна итерация Ллойда: центры -> центроиды ячеек.
    if (N > 1)
    {
        std::vector<RcPoint> moved(seeds.size());
        for (int i = 0; i < N; ++i)
        {
            const auto poly = clip_all(x0, y0, x1, y1, cell_hps(seeds, i, 0.0));
            moved[static_cast<std::size_t>(i)] = poly.empty() ? seeds[static_cast<std::size_t>(i)] : poly_centroid(poly);
        }
        seeds = moved;
    }
    else
        seeds[0] = {x0 + 0.5 * W, y0 + 0.5 * H};

    RandomCutLayout L;
    const double half_gap = 0.5 * P.gap;
    const double pi = std::acos(-1.0);
    for (int i = 0; i < N; ++i)
    {
        RandomCutFragment F;
        F.id = i;
        F.center = seeds[static_cast<std::size_t>(i)];
        F.cell = clip_all(x0, y0, x1, y1, cell_hps(seeds, i, 0.0));
        std::vector<RcHalfPlane> bis = cell_hps(seeds, i, half_gap);
        F.polygon = clip_all(x0, y0, x1, y1, bis);
        F.hps = rect_hps(x0, y0, x1, y1);
        F.hps.insert(F.hps.end(), bis.begin(), bis.end());
        F.area = F.polygon.empty() ? 0.0 : poly_area(F.polygon);
        L.frags.push_back(F);
    }
    // 3. Поворот и фаза сетки (после раскройки, по фрагментам подряд).
    for (auto& F : L.frags)
    {
        F.theta = (2.0 * rng.uniform() - 1.0) * P.rot_max_deg * pi / 180.0;
        F.fx = rng.uniform();
        F.fy = rng.uniform();
    }
    // 4. Швы: биссектор i-j внутри ячейки i (без зазора), полуплоскость самого j не учитывается.
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
        {
            const RcPoint& si = seeds[static_cast<std::size_t>(i)];
            const RcPoint& sj = seeds[static_cast<std::size_t>(j)];
            const double dx = sj[0] - si[0], dy = sj[1] - si[1];
            const double len = std::hypot(dx, dy);
            const RcPoint m = {0.5 * (si[0] + sj[0]), 0.5 * (si[1] + sj[1])};
            const RcPoint td = {-dy / len, dx / len};
            std::vector<RcHalfPlane> hp = rect_hps(x0, y0, x1, y1);
            const auto others = cell_hps(seeds, i, 0.0, j);
            hp.insert(hp.end(), others.begin(), others.end());
            double lo = 0.0, hi = 0.0;
            if (!line_interval(m, td, hp, lo, hi) || hi - lo <= 1e-9) continue;
            RandomCutSeam S;
            S.i = i;
            S.j = j;
            S.a = {m[0] + lo * td[0], m[1] + lo * td[1]};
            S.b = {m[0] + hi * td[0], m[1] + hi * td[1]};
            S.length = hi - lo;
            if (P.stitch_step > 0.0)
            {
                const long long ns = std::max<long long>(1, static_cast<long long>(std::floor(S.length / P.stitch_step)));
                for (long long s = 0; s < ns; ++s)
                {
                    const double t = lo + (static_cast<double>(s) + 0.5) * S.length / static_cast<double>(ns);
                    S.points.push_back({m[0] + t * td[0], m[1] + t * td[1]});
                }
            }
            L.seams.push_back(S);
        }
    return L;
}

// Нить фрагмента в плоскости группы (z относительно cz слоя).
struct PlaneFiber
{
    char direction = 'x';
    int index = 0;
    double offset = 0.0;
    int frag = 0;
    std::vector<std::array<double, 3>> pts;
    std::vector<std::size_t> vertex_bead;
};

} // namespace plain_weave_random_cut_detail

inline RandomCut build_plain_weave_random_cut(const RandomCutParams& Pin)
{
    using namespace plain_weave_detail;
    using namespace plain_weave_random_cut_detail;
    RandomCut R;
    RandomCutParams P = Pin;
    const PlainWeaveParams& q = P.base;
    const double p = q.pitch, d = q.d, A = q.amp;
    const int k = q.seg;

    // --- Проверки параметров.
    auto fin = [](double v, const char* name) {
        if (!std::isfinite(v))
            throw std::invalid_argument(std::string("plain_weave_random_cut: параметр ") + name + " должен быть конечным числом");
    };
    fin(q.d, "d");
    fin(q.pitch, "pitch");
    fin(q.amp, "amp");
    fin(q.cx, "cx");
    fin(q.cy, "cy");
    fin(q.cz, "cz");
    fin(q.E, "E");
    fin(q.nu, "nu");
    fin(q.rho, "rho");
    fin(q.eps_break, "eps_break");
    fin(q.mu_ff, "mu_ff");
    fin(q.mu_fs, "mu_fs");
    fin(P.size_x, "size_x");
    fin(P.size_y, "size_y");
    fin(P.gap, "gap");
    fin(P.rot_max_deg, "rot_max");
    fin(P.stitch_step, "stitch_step");
    fin(P.stitch_E, "stitch_E");
    fin(P.stitch_nu, "stitch_nu");
    fin(P.stitch_d, "stitch_d");
    fin(P.stitch_rho, "stitch_rho");
    fin(P.stitch_eps_break, "stitch_eps_break");
    fin(P.rc, "rc");
    fin(P.ply_pitch, "ply_pitch");
    if (q.clamp_edge_set) fin(q.clamp_edge, "clamp_edge");
    if (!(q.E > 0.0)) throw std::invalid_argument("plain_weave_random_cut: модуль E должен быть > 0");
    if (!(q.rho > 0.0)) throw std::invalid_argument("plain_weave_random_cut: плотность rho должна быть > 0");
    if (!(q.eps_break > 0.0)) throw std::invalid_argument("plain_weave_random_cut: eps_break должен быть > 0");
    if (q.mu_ff < 0.0 || q.mu_fs < 0.0) throw std::invalid_argument("plain_weave_random_cut: mu_ff, mu_fs должны быть >= 0");
    if (!(d > 0.0)) throw std::invalid_argument("plain_weave_random_cut: диаметр d должен быть > 0");
    if (!(p > 0.0)) throw std::invalid_argument("plain_weave_random_cut: шаг p должен быть > 0");
    if (!(A > 0.0)) throw std::invalid_argument("plain_weave_random_cut: амплитуда A должна быть > 0");
    if (k < 2 || k > plain_weave_max_seg)
        throw std::invalid_argument("plain_weave_random_cut: seg должно быть в [2, " + std::to_string(plain_weave_max_seg) + "]");
    if (P.n_frag < 1 || P.n_frag > max_fragments)
        throw std::invalid_argument("plain_weave_random_cut: число фрагментов N должно быть в [1, " +
                                    std::to_string(max_fragments) + "]");
    if (!(P.size_x > 0.0) || !(P.size_y > 0.0))
        throw std::invalid_argument("plain_weave_random_cut: размеры W, H должны быть > 0");
    if (P.gap < 0.0) P.gap = d;
    if (!(P.rot_max_deg >= 0.0 && P.rot_max_deg <= 180.0))
        throw std::invalid_argument("plain_weave_random_cut: rot_max должен быть в [0, 180] градусов");
    if (P.stitch_step < 0.0) throw std::invalid_argument("plain_weave_random_cut: stitch_step должен быть >= 0");
    if (P.stitch_step > 0.0 && (!(P.stitch_E > 0.0) || !(P.stitch_d > 0.0) || !(P.stitch_rho > 0.0) ||
                                !(P.stitch_eps_break > 0.0)))
        throw std::invalid_argument("plain_weave_random_cut: материал стежка: E, d, rho, eps_break должны быть > 0");
    if (!(P.rc > 0.0)) throw std::invalid_argument("plain_weave_random_cut: rc должен быть > 0");
    if (P.layers < 1 || P.layers > max_layers)
        throw std::invalid_argument("plain_weave_random_cut: число слоёв должно быть в [1, " + std::to_string(max_layers) + "]");
    if (P.layers_per_cut == 0) P.layers_per_cut = P.layers;
    if (P.layers_per_cut < 1 || P.layers_per_cut > P.layers)
        throw std::invalid_argument("plain_weave_random_cut: layers_per_cut должно быть в [1, layers]");
    if (P.ply_pitch < 0.0) P.ply_pitch = 2.0 * A + d + 1.0e-6;
    if (P.layers > 1 && !(P.ply_pitch > 0.0))
        throw std::invalid_argument("plain_weave_random_cut: ply_pitch должен быть > 0");
    {
        // Оценка числа бусин до построения: на слой 2*W*H/(p*h) внутренних + по 2 концевых на нить
        // (нитей не больше 2*(W+H)/p на фрагмент с запасом).
        const double h = p / k;
        const double per_ply = 2.0 * P.size_x * P.size_y / (p * h) +
                               4.0 * static_cast<double>(P.n_frag) * (P.size_x + P.size_y) / p * 2.0;
        const double est = per_ply * static_cast<double>(P.layers);
        if (!(est <= plain_weave_max_beads))
        {
            std::ostringstream s;
            s << "plain_weave_random_cut: оценка числа бусин " << est << " превышает предел " << plain_weave_max_beads;
            throw std::invalid_argument(s.str());
        }
    }

    R.params = P;
    PlainWeave& W = R.W;
    W.params = q;
    auto penalty = [&](bool bad, const std::string& msg) {
        if (!bad) return;
        if (!q.allow_penetration) throw std::invalid_argument(msg + " (разрешить: --allow-penetration)");
        const std::string w = "WARN " + msg;
        W.warnings.push_back(w);
        std::fprintf(stderr, "%s\n", w.c_str());
    };
    {
        std::ostringstream s1, s2, s3, s4;
        s1 << "plain_weave_random_cut: p < d (" << p << " < " << d << "): соседние параллельные нити перекрываются";
        s2 << "plain_weave_random_cut: 2A < d (" << 2.0 * A << " < " << d << "): нити проникают друг в друга в пересечении";
        s3 << "plain_weave_random_cut: gap < d (" << P.gap << " < " << d << "): нити соседних фрагментов проникают";
        s4 << "plain_weave_random_cut: ply_pitch < 2A + d (" << P.ply_pitch << " < " << 2.0 * A + d
           << "): слои проникают друг в друга";
        penalty(p < d, s1.str());
        penalty(2.0 * A < d, s2.str());
        penalty(P.gap < d, s3.str());
        penalty(P.layers > 1 && P.ply_pitch < 2.0 * A + d, s4.str());
    }

    const double x0 = q.cx - 0.5 * P.size_x, y0 = q.cy - 0.5 * P.size_y;
    const int ncut = (P.layers + P.layers_per_cut - 1) / P.layers_per_cut;

    // --- Раскройки групп (один поток ГСЧ).
    Rng rng(P.seed);
    for (int g = 0; g < ncut; ++g) R.cuts.push_back(make_cut(rng, P, x0, y0));

    // --- Полуволна и профиль.
    const Curve C{q.profile, p, A};
    std::vector<double> hs(static_cast<std::size_t>(k)), hz(static_cast<std::size_t>(k));
    {
        const std::vector<double> Ts = equal_chords(C, 0.0, 1.0, k);
        for (int m = 0; m < k; ++m)
        {
            double al, zz;
            C.at(Ts[static_cast<std::size_t>(m)], al, zz);
            hs[static_cast<std::size_t>(m)] = al;
            hz[static_cast<std::size_t>(m)] = zz;
        }
        hs[0] = 0.0;
        hz[0] = A;
    }
    const double h = p / k;
    // z профиля в точке a (координата вдоль нити, вершины в a = w p), без знака нити.
    auto zc_at = [&](double a) {
        const double w = std::floor(a / p);
        const double s = a - w * p; // [0, p)
        double t = s / p;
        if (q.profile == PlainWeaveProfile::Arcs)
        {
            double lo = 0.0, hi = 1.0;
            for (int it = 0; it < 200; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                if (mid == lo || mid == hi) break;
                double ss, zz;
                C.half(mid, ss, zz);
                if (ss < s)
                    lo = mid;
                else
                    hi = mid;
            }
            t = 0.5 * (lo + hi);
        }
        double ss, zz;
        C.half(t, ss, zz);
        const long long wi = static_cast<long long>(w);
        return (((wi % 2) + 2) % 2 == 1) ? -zz : zz;
    };
    auto parity_sign = [](long long n) { return (((n % 2) + 2) % 2 == 0) ? 1.0 : -1.0; };

    // --- Нити каждой раскройки в плоскости (z относительно cz слоя).
    std::vector<std::vector<PlaneFiber>> plane(static_cast<std::size_t>(ncut));
    for (int g = 0; g < ncut; ++g)
    {
        for (const auto& F : R.cuts[static_cast<std::size_t>(g)].frags)
        {
            if (F.polygon.empty()) continue;
            const RcPoint u = {std::cos(F.theta), std::sin(F.theta)};
            const RcPoint v = {-u[1], u[0]};
            for (int dir = 0; dir < 2; ++dir)
            {
                const bool warp = dir == 0;
                const RcPoint along = warp ? u : v;
                const RcPoint across = warp ? v : u;
                const double f_ac = warp ? F.fx : F.fy; // фаза поперёк (положение нити)
                const double f_al = warp ? F.fy : F.fx; // фаза вдоль (положение вершин)
                double smin = std::numeric_limits<double>::infinity(), smax = -smin;
                for (const auto& V : F.polygon)
                {
                    const double s = (V[0] - F.center[0]) * across[0] + (V[1] - F.center[1]) * across[1];
                    smin = std::min(smin, s);
                    smax = std::max(smax, s);
                }
                const long long j_lo = static_cast<long long>(std::floor(smin / p - f_ac)) - 1;
                const long long j_hi = static_cast<long long>(std::ceil(smax / p - f_ac)) + 1;
                for (long long j = j_lo; j <= j_hi; ++j)
                {
                    const double off = (static_cast<double>(j) + f_ac) * p;
                    const RcPoint p0 = {F.center[0] + off * across[0], F.center[1] + off * across[1]};
                    double t0 = 0.0, t1 = 0.0;
                    if (!line_interval(p0, along, F.hps, t0, t1) || t1 - t0 < 0.25 * h) continue;
                    const double a0 = t0 - f_al * p, a1 = t1 - f_al * p;
                    const double sg0 = warp ? parity_sign(j) : -parity_sign(j);
                    PlaneFiber PF;
                    PF.direction = warp ? 'x' : 'y';
                    PF.index = static_cast<int>(j);
                    PF.offset = off;
                    PF.frag = F.id;
                    auto bead = [&](double a, double zrel) -> std::array<double, 3> {
                        const double t = a + f_al * p;
                        return {p0[0] + t * along[0], p0[1] + t * along[1], sg0 * zrel};
                    };
                    std::vector<std::array<double, 3>> in;
                    std::vector<std::size_t> vb;
                    const long long w_lo = static_cast<long long>(std::floor(a0 / p));
                    const long long w_hi = static_cast<long long>(std::floor(a1 / p));
                    for (long long w = w_lo; w <= w_hi; ++w)
                    {
                        const double sw = parity_sign(w);
                        for (int m = 0; m < k; ++m)
                        {
                            const double a = static_cast<double>(w) * p + hs[static_cast<std::size_t>(m)];
                            if (a < a0 || a > a1) continue;
                            if (m == 0) vb.push_back(in.size());
                            in.push_back(bead(a, sw * hz[static_cast<std::size_t>(m)]));
                        }
                    }
                    auto chord3 = [](const std::array<double, 3>& x, const std::array<double, 3>& y) {
                        return std::sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
                                         (x[2] - y[2]) * (x[2] - y[2]));
                    };
                    const std::array<double, 3> e0 = bead(a0, zc_at(a0));
                    const std::array<double, 3> e1 = bead(a1, zc_at(a1));
                    if (in.empty() || chord3(e0, in.front()) >= 0.2 * h)
                    {
                        PF.pts.push_back(e0);
                        for (auto& x : vb) ++x;
                    }
                    PF.pts.insert(PF.pts.end(), in.begin(), in.end());
                    PF.vertex_bead = vb;
                    if (chord3(PF.pts.back(), e1) >= 0.2 * h) PF.pts.push_back(e1);
                    if (PF.pts.size() < 2) continue;
                    plane[static_cast<std::size_t>(g)].push_back(std::move(PF));
                }
            }
        }
        // Порядок: фрагмент -> основа -> уток (stable: внутри уже по возрастанию индекса).
        std::stable_sort(plane[static_cast<std::size_t>(g)].begin(), plane[static_cast<std::size_t>(g)].end(),
                         [](const PlaneFiber& a, const PlaneFiber& b) {
                             if (a.frag != b.frag) return a.frag < b.frag;
                             return a.direction < b.direction;
                         });
    }

    // --- Плоские массивы бусин раскроек (для поиска соседей) и метки фрагментов.
    std::vector<std::vector<std::array<double, 3>>> gpts(static_cast<std::size_t>(ncut));
    std::vector<std::vector<int>> gfrag(static_cast<std::size_t>(ncut));
    for (int g = 0; g < ncut; ++g)
        for (const auto& PF : plane[static_cast<std::size_t>(g)])
            for (const auto& X : PF.pts)
            {
                gpts[static_cast<std::size_t>(g)].push_back(X);
                gfrag[static_cast<std::size_t>(g)].push_back(PF.frag);
            }

    // --- Слои.
    int nwarp = 0, nweft = 0, obj = 0;
    std::size_t bead_run = 0;
    for (int l = 0; l < P.layers; ++l)
    {
        const int g = l / P.layers_per_cut;
        const double cz = q.cz - static_cast<double>(l) * P.ply_pitch;
        R.ply_cut.push_back(g);
        R.ply_cz.push_back(cz);
        R.ply_bead_first.push_back(bead_run);
        int cur_frag = -1;
        for (const auto& PF : plane[static_cast<std::size_t>(g)])
        {
            if (PF.frag != cur_frag)
            {
                RandomCutRange rg;
                rg.ply = l;
                rg.frag = PF.frag;
                rg.fiber_first = W.fibers.size();
                rg.bead_first = bead_run;
                R.ranges.push_back(rg);
                cur_frag = PF.frag;
            }
            PlainWeaveFiber F;
            F.object_id = obj++;
            F.direction = PF.direction;
            F.index = PF.index;
            F.offset = PF.offset;
            F.pts = PF.pts;
            for (auto& X : F.pts) X[2] += cz;
            F.vertex_bead = PF.vertex_bead;
            (PF.direction == 'x' ? nwarp : nweft)++;
            R.ranges.back().fiber_count++;
            R.ranges.back().bead_count += F.pts.size();
            bead_run += F.pts.size();
            W.fibers.push_back(std::move(F));
            R.fiber_ply.push_back(l);
            R.fiber_fragment.push_back(PF.frag);
        }
    }
    W.params.nx = nwarp;
    W.params.ny = nweft;
    W.params.overhang_x = W.params.overhang_y = 0.0;

    // --- Проверка зазора: мин. расстояние в плане между бусинами разных фрагментов (по раскройкам).
    R.min_cross_dist = std::numeric_limits<double>::infinity();
    for (int g = 0; g < ncut; ++g)
    {
        const auto& pts = gpts[static_cast<std::size_t>(g)];
        const auto& fr = gfrag[static_cast<std::size_t>(g)];
        const double cell = std::max(P.gap, 0.5 * p);
        const PlanGrid G(pts, cell);
        for (std::size_t a = 0; a < pts.size(); ++a)
        {
            const long long ia = G.ix(pts[a][0]), ib = G.ix(pts[a][1]);
            for (long long da = -1; da <= 1; ++da)
                for (long long db = -1; db <= 1; ++db)
                {
                    const auto* v = G.at(ia + da, ib + db);
                    if (!v) continue;
                    for (std::size_t b : *v)
                    {
                        if (b <= a || fr[a] == fr[b]) continue;
                        const double dd = std::hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]);
                        R.min_cross_dist = std::min(R.min_cross_dist, dd);
                    }
                }
        }
    }
    if (R.min_cross_dist < P.gap - 1e-9 * p)
    {
        std::ostringstream s;
        s << "plain_weave_random_cut: расстояние в плане между бусинами разных фрагментов " << R.min_cross_dist
          << " меньше зазора gap = " << P.gap;
        throw std::runtime_error(s.str());
    }

    // --- Стежки.
    if (P.stitch_step > 0.0)
    {
        std::vector<PlanGrid> grids;
        for (int g = 0; g < ncut; ++g) grids.emplace_back(gpts[static_cast<std::size_t>(g)], std::max(p, 1e-12));
        std::set<std::pair<std::size_t, std::size_t>> seen;
        auto add = [&](const RandomCutStitch& S) {
            const std::pair<std::size_t, std::size_t> key = {std::min(S.a, S.b), std::max(S.a, S.b)};
            if (S.a == S.b || !seen.insert(key).second) return false;
            R.stitches.push_back(S);
            return true;
        };
        // В слоях: через швы своей раскройки.
        for (int l = 0; l < P.layers; ++l)
        {
            const int g = R.ply_cut[static_cast<std::size_t>(l)];
            const auto& pts = gpts[static_cast<std::size_t>(g)];
            const auto& fr = gfrag[static_cast<std::size_t>(g)];
            const auto& cut = R.cuts[static_cast<std::size_t>(g)];
            for (std::size_t s = 0; s < cut.seams.size(); ++s)
            {
                const auto& S = cut.seams[s];
                for (const auto& pt : S.points)
                {
                    const std::size_t bi = grids[static_cast<std::size_t>(g)].nearest(pts, pt, [&](std::size_t x) { return fr[x] == S.i; });
                    const std::size_t bj = grids[static_cast<std::size_t>(g)].nearest(pts, pt, [&](std::size_t x) { return fr[x] == S.j; });
                    if (bi == static_cast<std::size_t>(-1) || bj == static_cast<std::size_t>(-1)) continue;
                    RandomCutStitch st;
                    st.a = R.ply_bead_first[static_cast<std::size_t>(l)] + bi;
                    st.b = R.ply_bead_first[static_cast<std::size_t>(l)] + bj;
                    st.ply_a = st.ply_b = l;
                    st.frag_a = S.i;
                    st.frag_b = S.j;
                    st.cut = g;
                    st.seam = static_cast<int>(s);
                    if (add(st)) ++R.n_stitches_in;
                }
            }
        }
        // Сквозные: точки швов всех раскроек, пары соседних слоёв.
        for (int l = 0; l + 1 < P.layers; ++l)
        {
            const int ga = R.ply_cut[static_cast<std::size_t>(l)], gb = R.ply_cut[static_cast<std::size_t>(l + 1)];
            for (int g = 0; g < ncut; ++g)
            {
                const auto& cut = R.cuts[static_cast<std::size_t>(g)];
                for (std::size_t s = 0; s < cut.seams.size(); ++s)
                    for (const auto& pt : cut.seams[s].points)
                    {
                        const std::size_t ba = grids[static_cast<std::size_t>(ga)].nearest(gpts[static_cast<std::size_t>(ga)], pt, [](std::size_t) { return true; });
                        const std::size_t bb = grids[static_cast<std::size_t>(gb)].nearest(gpts[static_cast<std::size_t>(gb)], pt, [](std::size_t) { return true; });
                        if (ba == static_cast<std::size_t>(-1) || bb == static_cast<std::size_t>(-1)) continue;
                        RandomCutStitch st;
                        st.through = true;
                        st.a = R.ply_bead_first[static_cast<std::size_t>(l)] + ba;
                        st.b = R.ply_bead_first[static_cast<std::size_t>(l + 1)] + bb;
                        st.ply_a = l;
                        st.ply_b = l + 1;
                        st.frag_a = gfrag[static_cast<std::size_t>(ga)][ba];
                        st.frag_b = gfrag[static_cast<std::size_t>(gb)][bb];
                        st.cut = g;
                        st.seam = static_cast<int>(s);
                        if (add(st)) ++R.n_stitches_through;
                    }
            }
        }
        for (const auto& S : R.stitches) R.stitch_pairs.emplace_back(S.a, S.b);
    }

    // --- Контакты диагностики: разные фрагменты одного слоя, 3D-расстояние < rc.
    for (int g = 0; g < ncut; ++g)
    {
        const auto& pts = gpts[static_cast<std::size_t>(g)];
        const auto& fr = gfrag[static_cast<std::size_t>(g)];
        const PlanGrid G(pts, P.rc);
        std::vector<std::pair<std::size_t, std::size_t>> loc;
        for (std::size_t a = 0; a < pts.size(); ++a)
        {
            const long long ia = G.ix(pts[a][0]), ib = G.ix(pts[a][1]);
            for (long long da = -1; da <= 1; ++da)
                for (long long db = -1; db <= 1; ++db)
                {
                    const auto* v = G.at(ia + da, ib + db);
                    if (!v) continue;
                    for (std::size_t b : *v)
                    {
                        if (b <= a || fr[a] == fr[b]) continue;
                        const double dx = pts[a][0] - pts[b][0], dy = pts[a][1] - pts[b][1], dz = pts[a][2] - pts[b][2];
                        if (dx * dx + dy * dy + dz * dz < P.rc * P.rc) loc.emplace_back(a, b);
                    }
                }
        }
        std::sort(loc.begin(), loc.end());
        for (int l = 0; l < P.layers; ++l)
            if (R.ply_cut[static_cast<std::size_t>(l)] == g)
                for (const auto& pr : loc)
                    R.contacts.emplace_back(R.ply_bead_first[static_cast<std::size_t>(l)] + pr.first,
                                            R.ply_bead_first[static_cast<std::size_t>(l)] + pr.second);
    }
    std::sort(R.contacts.begin(), R.contacts.end());

    // --- Зажим (только явный --clamp с --clamp-edge): боксы как в build_plain_weave, без требований к c_e.
    if (q.clamp != PlainWeaveClamp::None)
    {
        if (!q.clamp_edge_set)
            throw std::invalid_argument("plain_weave_random_cut: --clamp задан, но не задана граница зажима --clamp-edge (c_e)");
        const bool clamp_y = q.clamp == PlainWeaveClamp::Y || q.clamp == PlainWeaveClamp::XY;
        const bool clamp_x = q.clamp == PlainWeaveClamp::X || q.clamp == PlainWeaveClamp::XY;
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (const auto& F : W.fibers)
            for (const auto& X : F.pts)
                for (int a = 0; a < 3; ++a)
                {
                    lo[a] = std::min(lo[a], X[a]);
                    hi[a] = std::max(hi[a], X[a]);
                }
        for (int a = 0; a < 3; ++a)
        {
            lo[a] -= d;
            hi[a] += d;
        }
        const double ce = q.clamp_edge;
        auto add_box = [&](int axis, double blo, double bhi) {
            ClampBox b1, b2;
            for (int a = 0; a < 3; ++a)
                b1.lo[a] = b2.lo[a] = lo[a], b1.hi[a] = b2.hi[a] = hi[a];
            b1.lo[axis] = bhi;
            b2.hi[axis] = blo;
            W.clamp_boxes.push_back(b1); // верхний конец (чётный номер)
            W.clamp_boxes.push_back(b2); // нижний конец (нечётный)
        };
        if (clamp_x) add_box(0, q.cx - ce, q.cx + ce);
        if (clamp_y) add_box(1, q.cy - ce, q.cy + ce);
        std::size_t total = 0;
        for (auto& F : W.fibers)
        {
            for (const auto& X : F.pts)
                for (std::size_t bi = 0; bi < W.clamp_boxes.size(); ++bi)
                {
                    const auto& B = W.clamp_boxes[bi];
                    bool in = true;
                    for (int a = 0; a < 3; ++a)
                        if (X[a] < B.lo[a] || X[a] > B.hi[a]) in = false;
                    if (!in) continue;
                    if (bi % 2 == 0)
                        ++F.n_clamp_hi;
                    else
                        ++F.n_clamp_lo;
                    break;
                }
            F.clamped = F.n_clamp_lo + F.n_clamp_hi > 0;
            total += F.n_clamp_lo + F.n_clamp_hi;
        }
        if (total == 0)
            throw std::invalid_argument("plain_weave_random_cut: ни одна бусина не попала в боксы зажима (c_e слишком велик?)");
    }
    return R;
}

// Сумма площадей рабочих многоугольников раскройки g и длин её швов.
inline double random_cut_area(const RandomCut& R, int g)
{
    double s = 0.0;
    for (const auto& F : R.cuts[static_cast<std::size_t>(g)].frags) s += F.area;
    return s;
}
inline double random_cut_seam_length(const RandomCut& R, int g)
{
    double s = 0.0;
    for (const auto& S : R.cuts[static_cast<std::size_t>(g)].seams) s += S.length;
    return s;
}

// Дополнение .fib: стежки и строки заголовка о раскройке.
inline PlainWeaveFibExtra random_cut_fib_extra(const RandomCut& R)
{
    const auto& P = R.params;
    PlainWeaveFibExtra ex;
    if (!R.stitch_pairs.empty())
    {
        ex.pairs = &R.stitch_pairs;
        ex.stitch_E = P.stitch_E;
        ex.stitch_nu = P.stitch_nu;
        ex.stitch_d = P.stitch_d;
        ex.stitch_rho = P.stitch_rho;
        ex.stitch_eps_break = P.stitch_eps_break;
    }
    ex.window_area = P.size_x * P.size_y;
    std::ostringstream s1, s2, s3;
    s1 << std::setprecision(std::numeric_limits<double>::max_digits10);
    s2 << std::setprecision(std::numeric_limits<double>::max_digits10);
    s3 << std::setprecision(std::numeric_limits<double>::max_digits10);
    s1 << "случайная раскройка: N=" << P.n_frag << " size=" << P.size_x << "x" << P.size_y << " seed=" << P.seed
       << " gap=" << P.gap << " rot_max=" << P.rot_max_deg << " (nx/ny выше -- число нитей основы/утка, overhang не используется)";
    s2 << "слоёв " << P.layers << ", слоёв на раскройку " << P.layers_per_cut << ", раскроек " << R.cuts.size()
       << ", ply_pitch=" << P.ply_pitch;
    s3 << "стежки: шаг " << P.stitch_step << ", в слоях " << R.n_stitches_in << ", сквозных " << R.n_stitches_through
       << "; материал 1 (заглушка) E=" << P.stitch_E << " nu=" << P.stitch_nu << " d=" << P.stitch_d
       << " rho=" << P.stitch_rho << " eps_break=" << P.stitch_eps_break;
    ex.comments = {s1.str(), s2.str(), s3.str()};
    return ex;
}

// <out>.cut.json: раскройки, диапазоны фрагментов, стежки, контакты диагностики.
inline void write_random_cut_json(std::ostream& os, const RandomCut& R)
{
    const auto& P = R.params;
    const auto& q = P.base;
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    auto pt = [&](const RcPoint& a) { os << "[" << a[0] << ", " << a[1] << "]"; };
    os << "{\n";
    os << "  \"format\": \"random_cut_v1\",\n";
    os << "  \"size\": [" << P.size_x << ", " << P.size_y << "],\n";
    os << "  \"origin\": [" << q.cx - 0.5 * P.size_x << ", " << q.cy - 0.5 * P.size_y << "],\n";
    os << "  \"seed\": " << P.seed << ",\n";
    os << "  \"n_frag\": " << P.n_frag << ",\n";
    os << "  \"gap\": " << P.gap << ",\n";
    os << "  \"p\": " << q.pitch << ",\n";
    os << "  \"d\": " << q.d << ",\n";
    os << "  \"amp\": " << q.amp << ",\n";
    os << "  \"seg\": " << q.seg << ",\n";
    os << "  \"profile\": \"" << (q.profile == PlainWeaveProfile::Cos ? "cos" : "arcs") << "\",\n";
    os << "  \"rot_max_deg\": " << P.rot_max_deg << ",\n";
    os << "  \"stitch_step\": " << P.stitch_step << ",\n";
    os << "  \"rc\": " << P.rc << ",\n";
    os << "  \"layers\": " << P.layers << ",\n";
    os << "  \"layers_per_cut\": " << P.layers_per_cut << ",\n";
    os << "  \"ply_pitch\": " << P.ply_pitch << ",\n";
    os << "  \"n_beads\": " << (R.ply_bead_first.empty() ? 0 : [&] {
        std::size_t n = 0;
        for (const auto& F : R.W.fibers) n += F.pts.size();
        return n;
    }()) << ",\n";
    os << "  \"n_fibers\": " << R.W.fibers.size() << ",\n";
    os << "  \"plies\": [";
    for (std::size_t l = 0; l < R.ply_cut.size(); ++l)
        os << (l ? ", " : "") << "{\"cut\": " << R.ply_cut[l] << ", \"cz\": " << R.ply_cz[l]
           << ", \"bead_first\": " << R.ply_bead_first[l] << "}";
    os << "],\n";
    os << "  \"cuts\": [\n";
    for (std::size_t g = 0; g < R.cuts.size(); ++g)
    {
        const auto& cut = R.cuts[g];
        os << "    {\"fragments\": [\n";
        for (std::size_t f = 0; f < cut.frags.size(); ++f)
        {
            const auto& F = cut.frags[f];
            os << "      {\"id\": " << F.id << ", \"center\": ";
            pt(F.center);
            os << ", \"theta\": " << F.theta << ", \"phase\": [" << F.fx << ", " << F.fy << "], \"area\": " << F.area
               << ", \"polygon\": [";
            for (std::size_t v = 0; v < F.polygon.size(); ++v)
            {
                if (v) os << ", ";
                pt(F.polygon[v]);
            }
            os << "], \"cell\": [";
            for (std::size_t v = 0; v < F.cell.size(); ++v)
            {
                if (v) os << ", ";
                pt(F.cell[v]);
            }
            os << "]}" << (f + 1 < cut.frags.size() ? "," : "") << "\n";
        }
        os << "     ],\n     \"seams\": [\n";
        for (std::size_t s = 0; s < cut.seams.size(); ++s)
        {
            const auto& S = cut.seams[s];
            os << "      {\"i\": " << S.i << ", \"j\": " << S.j << ", \"a\": ";
            pt(S.a);
            os << ", \"b\": ";
            pt(S.b);
            os << ", \"length\": " << S.length << ", \"points\": [";
            for (std::size_t v = 0; v < S.points.size(); ++v)
            {
                if (v) os << ", ";
                pt(S.points[v]);
            }
            os << "]}" << (s + 1 < cut.seams.size() ? "," : "") << "\n";
        }
        os << "     ]}" << (g + 1 < R.cuts.size() ? "," : "") << "\n";
    }
    os << "  ],\n";
    os << "  \"ranges\": [\n";
    for (std::size_t r = 0; r < R.ranges.size(); ++r)
    {
        const auto& G = R.ranges[r];
        os << "    {\"ply\": " << G.ply << ", \"frag\": " << G.frag << ", \"fiber_first\": " << G.fiber_first
           << ", \"fiber_count\": " << G.fiber_count << ", \"bead_first\": " << G.bead_first
           << ", \"bead_count\": " << G.bead_count << "}" << (r + 1 < R.ranges.size() ? "," : "") << "\n";
    }
    os << "  ],\n";
    os << "  \"stitch_columns\": [\"a\", \"b\", \"through\", \"ply_a\", \"ply_b\", \"frag_a\", \"frag_b\", \"cut\", \"seam\"],\n";
    os << "  \"stitches\": [";
    for (std::size_t s = 0; s < R.stitches.size(); ++s)
    {
        const auto& S = R.stitches[s];
        os << (s ? ",\n    " : "\n    ") << "[" << S.a << ", " << S.b << ", " << (S.through ? 1 : 0) << ", " << S.ply_a
           << ", " << S.ply_b << ", " << S.frag_a << ", " << S.frag_b << ", " << S.cut << ", " << S.seam << "]";
    }
    os << (R.stitches.empty() ? "" : "\n  ") << "],\n";
    os << "  \"contacts\": [";
    for (std::size_t c = 0; c < R.contacts.size(); ++c)
        os << (c ? (c % 8 == 0 ? ",\n    " : ", ") : "\n    ") << "[" << R.contacts[c].first << ", " << R.contacts[c].second << "]";
    os << (R.contacts.empty() ? "" : "\n  ") << "]\n";
    os << "}\n";
}

inline void print_random_cut_summary(std::FILE* f, const RandomCut& R)
{
    const auto& P = R.params;
    const PlainWeaveSummary S = plain_weave_summary(R.W);
    std::size_t seams = 0;
    for (const auto& c : R.cuts) seams += c.seams.size();
    std::fprintf(f, "случайная раскройка: фрагментов %d, размер %.9g x %.9g м, seed %llu, зазор %.9g м\n", P.n_frag,
                 P.size_x, P.size_y, static_cast<unsigned long long>(P.seed), P.gap);
    std::fprintf(f, "слоёв %d, слоёв на раскройку %d, раскроек %zu, ply_pitch %.9g м\n", P.layers, P.layers_per_cut,
                 R.cuts.size(), P.ply_pitch);
    std::fprintf(f, "швов (по всем раскройкам): %zu\n", seams);
    std::fprintf(f, "нитей: %zu (основа %d, уток %d)\n", S.fibers, R.W.params.nx, R.W.params.ny);
    std::fprintf(f, "бусин: %zu, связей нитей: %zu, зажатых бусин: %zu\n", S.beads, S.bonds, S.clamped_beads);
    std::fprintf(f, "стежков: %zu (в слоях %zu, сквозных %zu)\n", R.stitches.size(), R.n_stitches_in,
                 R.n_stitches_through);
    std::fprintf(f, "контактов диагностики (rc = %.9g м): %zu\n", P.rc, R.contacts.size());
    std::fprintf(f, "l0: min %.9g м, max %.9g м\n", S.l0_min, S.l0_max);
    std::fprintf(f, "мин. расстояние в плане между фрагментами: %.9g м\n", R.min_cross_dist);
    for (std::size_t g = 0; g < R.cuts.size(); ++g)
    {
        std::fprintf(f, "раскройка %zu: площади фрагментов, м^2:", g);
        for (const auto& F : R.cuts[g].frags) std::fprintf(f, " %.6g", F.area);
        std::fprintf(f, "; сумма %.9g (прямоугольник %.9g)\n", random_cut_area(R, static_cast<int>(g)),
                     P.size_x * P.size_y);
    }
}

} // namespace fiber
