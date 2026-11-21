#pragma once
// plain_weave_chess.hpp -- «шахматная» ткань: пакет из k слоёв, каждый слой -- два подслоя кусочков
// с перекрытием (решение владельца 2026-11-21, experiments/2026-11-21-fabric-kT-w50-chess-ball1g-v550).
//
// ГЕОМЕТРИЯ.
//   * Полное полотно base (nx нитей основы вдоль x, ny нитей утка вдоль y, build_plain_weave) делится на
//     кусочки по линиям шва centered_patch_seams (plain_weave_quilt.hpp): ядро -- patch x patch нитей от центра,
//     у краёв остаток. Кусочек = ядро плюс по одной нити с каждой стороны, где есть сосед (7 x 7 при patch = 5):
//     крайняя нить кусочка ложится на крайнюю нить ядра соседа (перекрытие на одну нить).
//   * Шахматная доска: цвет кусочка c = (-1)^(mx + my) относительно центрального (c = +1). Слой k (с нуля) --
//     подслои 2k (к лицу) и 2k + 1 (к тылу) с шагом sub_pitch; c_z подслоя s = base.cz - s * sub_pitch.
//     Кусочек цвета c лежит на подслое 2k, если c * sign_k > 0, иначе на 2k + 1; sign_k = +1, а в режиме Flip у
//     нечётных k (слои 2, 4, ...) sign_k = -1 (чётность доски меняется от слоя к слою), в режиме Same -- всегда +1.
//     sub_pitch = D = 2A + d + eps: и между подслоями одного слоя (перекрытие -- одно пересечение над другим),
//     и между подслоем 2k + 1 и 2k + 2 (как шаг слоёв сплошного пакета).
//   * Углы: где сходятся четыре кусочка, квадрат 2 x 2 пересечений (нити по обе стороны обоих швов) входит
//     в расширения двух диагональных соседей одного подслоя. Пересечение такого квадрата принадлежит только
//     кусочку, в ядро которого входит; полосы перекрытия у внутренних кусочков -- patch - 2 нити.
//   * Нить кусочка -- участок нити полного полотна (те же бусины и координаты, z сдвинут на подслой):
//     внутри полотна обрезается на середине полуволны за крайним пересечением кусочка (нужно чётное seg),
//     у края полотна сохраняет выступ полотна (overhang_x, overhang_y; зажатые выступы -- в зажим).
//   * Стяжки: в каждом пересечении перекрытия (два кусочка одного слоя друг над другом) -- связь от нижней из
//     двух вершинных бусин верхнего кусочка к верхней вершинной бусине нижнего (длина sub_pitch - 2A = d + eps).
//   * Зажим: как у полотна (CLAMP_BOX по |y - cy| >= c_e и/или |x - cx| >= c_e), боксы -- по габариту слоя.
//
// ПОРЯДОК НИТЕЙ. Слой k: кусочки по mx (вдоль x), внутри -- по my; в кусочке сначала основа по j, затем уток
// по i. Бусины пакета -- подряд по слоям и нитям (как write_fabric_stack_fib в сцене Guo 2019).
#include "fiber/plain_weave.hpp"
#include "fiber/plain_weave_quilt.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fiber
{

enum class ChessMode
{
    Flip, // у слоя k + 1 чётность доски обратная
    Same  // у всех слоёв одна доска
};

inline const char* chess_mode_name(ChessMode m) { return m == ChessMode::Flip ? "flip" : "same"; }

struct ChessParams
{
    PlainWeaveParams base; // полное полотно; base.cz -- c_z подслоя 0 (лицевого)
    int patch = 5;         // ядро кусочка, нитей
    int layers = 1;        // число слоёв k
    double sub_pitch = 0.0; // шаг подслоёв по z, м (> 2A + d)
    ChessMode mode = ChessMode::Flip;
};

// Кусочки вдоль одной оси: ядра [lo, hi] в индексах нитей, нити у швов, номер центрального кусочка.
struct ChessAxis
{
    std::vector<std::pair<int, int>> core;
    std::vector<char> adj; // [нить] -- по одну из сторон линии шва
    int center = 0;

    int count() const { return static_cast<int>(core.size()); }
    std::pair<int, int> ext(int m) const
    {
        const auto& c = core[static_cast<std::size_t>(m)];
        return {m > 0 ? c.first - 1 : c.first, m + 1 < count() ? c.second + 1 : c.second};
    }
    int patch_of(int i) const
    {
        for (int m = 0; m < count(); ++m)
            if (core[static_cast<std::size_t>(m)].first <= i && i <= core[static_cast<std::size_t>(m)].second) return m;
        throw std::logic_error("plain_weave_chess: нить вне кусочков");
    }
};

struct ChessFabric
{
    ChessParams params;
    std::vector<PlainWeave> plies;   // слой k: нити всех его кусочков
    std::vector<double> sub_cz;      // c_z подслоёв (2 k)
    std::vector<std::pair<std::size_t, std::size_t>> stitch_pairs; // глобальные бусины (верхний кусочек, нижний)
    std::vector<std::array<int, 3>> stitch_cross;                  // (слой, i, j) стяжки
    ChessAxis ax, ay;                // ax -- вдоль x (индексы утка i, ny нитей), ay -- вдоль y (индексы основы j, nx)
    std::size_t pieces_per_layer = 0, stitches_per_layer = 0;
    double z_top = 0.0, z_bottom = 0.0; // огибающая пакета: оси крайних подслоёв +- (A + d/2)
};

namespace plain_weave_chess_detail
{

inline ChessAxis make_axis(double c, double pitch, int n, int patch)
{
    ChessAxis A;
    const std::vector<double> seams = centered_patch_seams(c, pitch, n, patch);
    A.adj.assign(static_cast<std::size_t>(n), 0);
    int lo = 0;
    for (double s : seams)
    {
        const double si = (s - c) / pitch + 0.5 * (n - 1); // полуцелое
        const int left = static_cast<int>(std::floor(si));
        if (std::fabs(si - (left + 0.5)) > 1e-6) throw std::logic_error("plain_weave_chess: шов не посередине между нитями");
        A.core.push_back({lo, left});
        A.adj[static_cast<std::size_t>(left)] = 1;
        A.adj[static_cast<std::size_t>(left + 1)] = 1;
        lo = left + 1;
    }
    A.core.push_back({lo, n - 1});
    A.center = A.patch_of((n - 1) / 2);
    return A;
}

// Входит ли пересечение (i вдоль x, j вдоль y) в кусочек (mx, my).
inline bool owns(const ChessAxis& ax, const ChessAxis& ay, int mx, int my, int i, int j)
{
    const auto cx = ax.core[static_cast<std::size_t>(mx)], cy = ay.core[static_cast<std::size_t>(my)];
    const auto ex = ax.ext(mx), ey = ay.ext(my);
    const bool in_cx = cx.first <= i && i <= cx.second, in_cy = cy.first <= j && j <= cy.second;
    const bool in_ex = ex.first <= i && i <= ex.second, in_ey = ey.first <= j && j <= ey.second;
    if (in_cx && in_cy) return true;
    if (!(in_ex && in_ey)) return false;
    if (ax.adj[static_cast<std::size_t>(i)] && ay.adj[static_cast<std::size_t>(j)]) return false; // угловой квадрат
    return (in_cx && in_ey) || (in_ex && in_cy);
}

inline int color(const ChessAxis& ax, const ChessAxis& ay, int mx, int my)
{
    return ((mx - ax.center) + (my - ay.center)) % 2 == 0 ? 1 : -1;
}

} // namespace plain_weave_chess_detail

// Подслой (0 -- лицевой) кусочка (mx, my) слоя k.
inline int chess_sublayer(const ChessFabric& C, int k, int mx, int my)
{
    const int c = plain_weave_chess_detail::color(C.ax, C.ay, mx, my);
    const int sign = (C.params.mode == ChessMode::Flip && k % 2 == 1) ? -1 : 1;
    return 2 * k + (c * sign > 0 ? 0 : 1);
}

inline ChessFabric build_plain_weave_chess(const ChessParams& P)
{
    using namespace plain_weave_chess_detail;
    const PlainWeaveParams& q = P.base;
    if (P.layers < 1) throw std::invalid_argument("plain_weave_chess: число слоёв должно быть >= 1");
    if (q.seg % 2 != 0) throw std::invalid_argument("plain_weave_chess: нужно чётное seg (обрез нити на середине полуволны)");
    if (!(P.sub_pitch >= 2.0 * q.amp + q.d))
        throw std::invalid_argument("plain_weave_chess: шаг подслоёв меньше 2A + d -- перекрытия входят друг в друга");
    ChessFabric C;
    C.params = P;
    C.ax = make_axis(q.cx, q.pitch, q.ny, P.patch); // вдоль x: нити утка i
    C.ay = make_axis(q.cy, q.pitch, q.nx, P.patch); // вдоль y: нити основы j
    const ChessAxis &ax = C.ax, &ay = C.ay;
    const int Mx = ax.count(), My = ay.count();
    C.pieces_per_layer = static_cast<std::size_t>(Mx) * static_cast<std::size_t>(My);

    const PlainWeave W0 = build_plain_weave(q); // полное полотно на подслое 0
    const int half = q.seg / 2;
    for (int s = 0; s < 2 * P.layers; ++s) C.sub_cz.push_back(q.cz - static_cast<double>(s) * P.sub_pitch);

    // Участок нити полного полотна по пересечениям [a, b] (n -- число пересечений нити), сдвиг по z.
    auto cut = [&](const PlainWeaveFiber& F, int a, int b, int n, double dz) {
        PlainWeaveFiber G;
        G.direction = F.direction;
        G.index = F.index;
        G.offset = F.offset;
        const std::size_t start = a == 0 ? 0 : F.vertex_bead[static_cast<std::size_t>(a)] - static_cast<std::size_t>(half);
        const std::size_t end = b == n - 1 ? F.pts.size() - 1
                                           : F.vertex_bead[static_cast<std::size_t>(b)] + static_cast<std::size_t>(half);
        for (std::size_t t = start; t <= end; ++t)
        {
            auto p = F.pts[t];
            p[2] += dz;
            G.pts.push_back(p);
        }
        for (int v = a; v <= b; ++v) G.vertex_bead.push_back(F.vertex_bead[static_cast<std::size_t>(v)] - start);
        G.n_clamp_lo = a == 0 ? F.n_clamp_lo : 0;
        G.n_clamp_hi = b == n - 1 ? F.n_clamp_hi : 0;
        G.clamped = F.clamped && (G.n_clamp_lo + G.n_clamp_hi) > 0;
        return G;
    };

    // (слой, mx, my, направление, индекс нити) -> (слой, нить в слое)
    std::map<std::tuple<int, int, int, char, int>, std::size_t> where;
    int oid = 0;
    for (int k = 0; k < P.layers; ++k)
    {
        PlainWeave L;
        L.params = q;
        L.params.cz = C.sub_cz[static_cast<std::size_t>(2 * k)];
        for (int mx = 0; mx < Mx; ++mx)
            for (int my = 0; my < My; ++my)
            {
                const int s = chess_sublayer(C, k, mx, my);
                const double dz = -static_cast<double>(s) * P.sub_pitch;
                const auto ex = ax.ext(mx), ey = ay.ext(my);
                for (int pass = 0; pass < 2; ++pass) // 0 -- основа (вдоль x), 1 -- уток (вдоль y)
                {
                    const auto across = pass == 0 ? ey : ex, along = pass == 0 ? ex : ey;
                    for (int t = across.first; t <= across.second; ++t)
                    {
                        int a = -1, b = -1;
                        for (int u = along.first; u <= along.second; ++u)
                        {
                            const bool own = pass == 0 ? owns(ax, ay, mx, my, u, t) : owns(ax, ay, mx, my, t, u);
                            if (!own) continue;
                            if (a < 0) a = u;
                            else if (u != b + 1)
                                throw std::logic_error("plain_weave_chess: пересечения нити кусочка не подряд");
                            b = u;
                        }
                        if (a < 0) continue;
                        const PlainWeaveFiber& F =
                            W0.fibers[static_cast<std::size_t>(pass == 0 ? t : q.nx + t)];
                        PlainWeaveFiber G = cut(F, a, b, pass == 0 ? q.ny : q.nx, dz);
                        G.object_id = oid++;
                        where[{k, mx, my, pass == 0 ? 'x' : 'y', t}] = L.fibers.size();
                        L.fibers.push_back(std::move(G));
                    }
                }
            }
        // Зажим: боксы полотна по габариту слоя; проверка «в боксе <=> зажатый выступ».
        if (q.clamp != PlainWeaveClamp::None)
        {
            double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
            for (const auto& F : L.fibers)
                for (const auto& p : F.pts)
                    for (int c = 0; c < 3; ++c)
                    {
                        lo[c] = std::min(lo[c], p[c]);
                        hi[c] = std::max(hi[c], p[c]);
                    }
            for (int c = 0; c < 3; ++c)
            {
                lo[c] -= q.d;
                hi[c] += q.d;
            }
            const bool cy = q.clamp == PlainWeaveClamp::Y || q.clamp == PlainWeaveClamp::XY;
            const bool cxx = q.clamp == PlainWeaveClamp::X || q.clamp == PlainWeaveClamp::XY;
            auto add_box = [&](int axis, double blo, double bhi) {
                ClampBox b1, b2;
                for (int c = 0; c < 3; ++c) b1.lo[c] = b2.lo[c] = lo[c], b1.hi[c] = b2.hi[c] = hi[c];
                b1.lo[axis] = bhi;
                b2.hi[axis] = blo;
                L.clamp_boxes.push_back(b1);
                L.clamp_boxes.push_back(b2);
            };
            if (cxx) add_box(0, q.cx - q.clamp_edge, q.cx + q.clamp_edge);
            if (cy) add_box(1, q.cy - q.clamp_edge, q.cy + q.clamp_edge);
            for (const auto& F : L.fibers)
            {
                std::size_t in = 0;
                for (const auto& p : F.pts)
                {
                    bool any = false;
                    for (const auto& B : L.clamp_boxes)
                    {
                        bool inside = true;
                        for (int c = 0; c < 3; ++c) inside = inside && p[c] >= B.lo[c] && p[c] <= B.hi[c];
                        any = any || inside;
                    }
                    if (any) ++in;
                }
                if (in != F.n_clamp_lo + F.n_clamp_hi)
                {
                    std::ostringstream m;
                    m << "plain_weave_chess: слой " << (k + 1) << ", нить " << F.direction << F.index << ": в боксах зажима "
                      << in << " бусин, ожидалось " << (F.n_clamp_lo + F.n_clamp_hi);
                    throw std::runtime_error(m.str());
                }
            }
        }
        C.plies.push_back(std::move(L));
    }

    // Глобальные индексы бусин: подряд по слоям и нитям.
    std::vector<std::vector<std::size_t>> first(C.plies.size());
    std::size_t g = 0;
    for (std::size_t k = 0; k < C.plies.size(); ++k)
        for (const auto& F : C.plies[k].fibers)
        {
            first[k].push_back(g);
            g += F.pts.size();
        }
    auto vertex = [&](int k, int mx, int my, char dir, int idx, int v, double& z) {
        const auto it = where.find({k, mx, my, dir, idx});
        if (it == where.end()) throw std::logic_error("plain_weave_chess: нет нити кусочка в пересечении");
        const PlainWeaveFiber& F = C.plies[static_cast<std::size_t>(k)].fibers[it->second];
        // вершины нити -- подряд с первого её пересечения
        const int v0 = dir == 'x' ? ax.ext(mx).first : ay.ext(my).first;
        int a = -1;
        for (int u = v0; a < 0; ++u)
            if (dir == 'x' ? owns(ax, ay, mx, my, u, idx) : owns(ax, ay, mx, my, idx, u)) a = u;
        const std::size_t loc = F.vertex_bead.at(static_cast<std::size_t>(v - a));
        z = F.pts[loc][2];
        return first[static_cast<std::size_t>(k)][it->second] + loc;
    };

    // Стяжки: пересечения перекрытия.
    for (int k = 0; k < P.layers; ++k)
    {
        std::size_t n_k = 0;
        for (int i = 0; i < q.ny; ++i)
            for (int j = 0; j < q.nx; ++j)
            {
                std::vector<std::pair<int, int>> ow;
                const int px = ax.patch_of(i), py = ay.patch_of(j);
                for (int mx = std::max(0, px - 1); mx <= std::min(Mx - 1, px + 1); ++mx)
                    for (int my = std::max(0, py - 1); my <= std::min(My - 1, py + 1); ++my)
                        if (owns(ax, ay, mx, my, i, j)) ow.push_back({mx, my});
                if (ow.size() < 2) continue;
                if (ow.size() > 2 ||
                    chess_sublayer(C, k, ow[0].first, ow[0].second) == chess_sublayer(C, k, ow[1].first, ow[1].second))
                {
                    std::ostringstream m;
                    m << "plain_weave_chess: пересечение (" << i << ", " << j << ") слоя " << (k + 1) << ": " << ow.size()
                      << " кусочка на одном подслое";
                    throw std::logic_error(m.str());
                }
                if (chess_sublayer(C, k, ow[0].first, ow[0].second) > chess_sublayer(C, k, ow[1].first, ow[1].second))
                    std::swap(ow[0], ow[1]);
                double zw, zf;
                const std::size_t uw = vertex(k, ow[0].first, ow[0].second, 'x', j, i, zw);
                const std::size_t uf = vertex(k, ow[0].first, ow[0].second, 'y', i, j, zf);
                const std::size_t up = zw < zf ? uw : uf; // нижняя нить верхнего кусочка
                const std::size_t lw = vertex(k, ow[1].first, ow[1].second, 'x', j, i, zw);
                const std::size_t lf = vertex(k, ow[1].first, ow[1].second, 'y', i, j, zf);
                const std::size_t dn = zw > zf ? lw : lf; // верхняя нить нижнего кусочка
                C.stitch_pairs.push_back({up, dn});
                C.stitch_cross.push_back({k, i, j});
                ++n_k;
            }
        if (k == 0) C.stitches_per_layer = n_k;
    }
    C.z_top = C.sub_cz.front() + q.amp + 0.5 * q.d;
    C.z_bottom = C.sub_cz.back() - q.amp - 0.5 * q.d;
    return C;
}

} // namespace fiber
