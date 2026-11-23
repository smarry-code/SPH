// tests/test_plain_weave_chess.cpp
// Тесты «шахматной» ткани (include/fiber/plain_weave_chess.hpp): малое полотно Guo 2019 (43 x 43 нити,
// p = 1,18 мм, d = 0,39 мм, A = 0,6 d, seg = 2, зажим |y| >= 25 мм, выступ 5,22 мм), кусочки 5 x 5.
//
// Разделы:
//   1 test_axes       -- кусочки по оси: 4 + 7*5 + 4 нити, центральный 19..23, нити у швов;
//   2 test_beads      -- нить кусочка == участок нити полного полотна, сдвинутый на свой подслой (побитово);
//   3 test_sublayers  -- подслои: flip / same, центральный кусочок на лицевом подслое;
//   4 test_owners     -- каждое пересечение -- у 1 или 2 кусочков; 2 -- на разных подслоях одного слоя;
//   5 test_stitches   -- число (864 на слой), вертикальны, длина D - 2A, верх -> низ;
//   6 test_clamp      -- зажатые выступы только у кусочков у края, в боксах ровно они;
//   7 test_distance   -- оси нитей разных нитей не ближе d (проникания нет), до 2 слоёв;
//   8 test_errors     -- нечётный seg, малый шаг подслоёв, |shift| >= patch -- исключения;
//   9 test_shift      -- сдвиг раскройки нечётных слоёв (shift = 2): кусочки 6 + 7*5 + 2, швы соседних слоёв
//                        не совпадают, владельцы и стяжки нечётного слоя, расстояние осей >= d;
//  10 test_through    -- стяжки между слоями: число = перебор пересечений, вертикальны, длина d + eps, соединяют
//                        подслой 2k + 1 с 2k + 2 (flip / same, сдвиг 0 / 2, 3 слоя).
#include "fiber/plain_weave.hpp"
#include "fiber/plain_weave_chess.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using fiber::ChessFabric;
using fiber::ChessMode;
using fiber::ChessParams;
using fiber::PlainWeaveParams;

static int failures = 0;

static void check(bool cond, const std::string& what)
{
    std::printf("  [%s] %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

// Малое полотно серии w50 (cases/guo2019/experiment/guo2019_scene.hpp, --fabric-n=43 --fabric-clamp-edge=25e-3
// --fabric-overhang-y=5.22e-3).
static PlainWeaveParams w50()
{
    PlainWeaveParams q;
    fiber::apply_plain_weave_preset("twaron_ct717", q);
    q.amp = 0.6 * q.d;
    q.seg = 2;
    q.nx = q.ny = 43;
    q.overhang_x = 0.5 * q.pitch;
    q.overhang_y = 5.22e-3;
    q.clamp = fiber::PlainWeaveClamp::Y;
    q.clamp_edge = 25e-3;
    q.clamp_edge_set = true;
    q.cx = q.cy = 0.0;
    q.cz = -(q.amp + 0.5 * q.d);
    return q;
}

static ChessParams params(int layers, ChessMode mode)
{
    ChessParams P;
    P.base = w50();
    P.patch = 5;
    P.layers = layers;
    P.sub_pitch = 2.0 * P.base.amp + P.base.d + 1e-6;
    P.mode = mode;
    return P;
}

static void test_axes()
{
    std::printf("1 test_axes\n");
    const ChessFabric C = fiber::build_plain_weave_chess(params(1, ChessMode::Flip));
    check(C.ax.count() == 9 && C.ay.count() == 9, "9 x 9 кусочков");
    check(C.ax.core.front() == std::make_pair(0, 3) && C.ax.core.back() == std::make_pair(39, 42), "края -- по 4 нити");
    check(C.ax.core[4] == std::make_pair(19, 23) && C.ax.center == 4, "центральный кусочек -- нити 19..23");
    int nadj = 0;
    for (char a : C.ax.adj) nadj += a ? 1 : 0;
    check(nadj == 16, "нитей у швов 16");
    check(C.ax.ext(0) == std::make_pair(0, 4) && C.ax.ext(4) == std::make_pair(18, 24), "расширение на нить, у края нет");
    check(C.pieces_per_layer == 81, "81 кусочек на слой");
}

static void test_beads()
{
    std::printf("2 test_beads\n");
    const ChessParams P = params(2, ChessMode::Flip);
    const ChessFabric C = fiber::build_plain_weave_chess(P);
    const fiber::PlainWeave W0 = fiber::build_plain_weave(P.base);
    bool same = true, vert = true;
    std::size_t nf = 0;
    for (std::size_t k = 0; k < C.plies.size(); ++k)
        for (const auto& F : C.plies[k].fibers)
        {
            ++nf;
            const auto& G = W0.fibers[static_cast<std::size_t>(F.direction == 'x' ? F.index : P.base.nx + F.index)];
            // ищем начало участка по x, y первой бусины
            std::size_t s0 = G.pts.size();
            for (std::size_t t = 0; t < G.pts.size(); ++t)
                if (G.pts[t][0] == F.pts[0][0] && G.pts[t][1] == F.pts[0][1]) { s0 = t; break; }
            if (s0 == G.pts.size() || s0 + F.pts.size() > G.pts.size()) { same = false; continue; }
            // подслой s по разности z; сдвиг генератора -- ровно -s * sub_pitch
            const double r = (F.pts[0][2] - G.pts[s0][2]) / P.sub_pitch;
            if (std::fabs(r - std::round(r)) > 1e-9) same = false;
            const double dz = -static_cast<double>(static_cast<long>(-std::round(r))) * P.sub_pitch;
            for (std::size_t t = 0; t < F.pts.size(); ++t)
            {
                const auto& a = F.pts[t];
                const auto& b = G.pts[s0 + t];
                if (a[0] != b[0] || a[1] != b[1] || a[2] != b[2] + dz) same = false;
            }
            for (std::size_t v = 0; v < F.vertex_bead.size(); ++v)
            {
                const auto& a = F.pts[F.vertex_bead[v]];
                if (std::fabs(std::fabs(a[2] - (P.base.cz + dz)) - P.base.amp) > 1e-12) vert = false;
            }
        }
    check(same, "бусины кусочков -- побитово участки нитей полотна со сдвигом на подслой");
    check(vert, "вершинные бусины -- на +-A от c_z подслоя");
    check(nf > 0, "нитей в пакете " + std::to_string(nf));
}

static void test_sublayers()
{
    std::printf("3 test_sublayers\n");
    const ChessFabric F = fiber::build_plain_weave_chess(params(3, ChessMode::Flip));
    const ChessFabric S = fiber::build_plain_weave_chess(params(3, ChessMode::Same));
    check(fiber::chess_sublayer(F, 0, 4, 4) == 0 && fiber::chess_sublayer(F, 0, 5, 4) == 1, "слой 1: центр на 0, сосед на 1");
    check(fiber::chess_sublayer(F, 1, 4, 4) == 3 && fiber::chess_sublayer(F, 1, 5, 4) == 2, "flip: слой 2 -- доска обратная");
    check(fiber::chess_sublayer(F, 2, 4, 4) == 4, "flip: слой 3 -- как слой 1");
    check(fiber::chess_sublayer(S, 1, 4, 4) == 2 && fiber::chess_sublayer(S, 1, 5, 4) == 3, "same: слой 2 -- та же доска");
    check(F.sub_cz.size() == 6 && std::fabs(F.sub_cz[1] - F.sub_cz[0] + F.params.sub_pitch) < 1e-15, "6 подслоёв с шагом D");
    const double A = F.params.base.amp, d = F.params.base.d;
    check(std::fabs(F.z_top) < 1e-15, "лицевая грань z = 0");
    check(std::fabs(F.z_bottom - (F.sub_cz.back() - A - 0.5 * d)) < 1e-15, "тыл -- нижний подслой");
}

static void test_owners()
{
    std::printf("4 test_owners\n");
    const ChessFabric C = fiber::build_plain_weave_chess(params(2, ChessMode::Flip));
    int n1 = 0, n2 = 0, bad = 0;
    for (int i = 0; i < 43; ++i)
        for (int j = 0; j < 43; ++j)
        {
            int n = 0;
            for (int mx = 0; mx < 9; ++mx)
                for (int my = 0; my < 9; ++my)
                    if (fiber::plain_weave_chess_detail::owns(C.ax, C.ay, mx, my, i, j)) ++n;
            if (n == 1) ++n1;
            else if (n == 2) ++n2;
            else ++bad;
        }
    check(bad == 0, "нет пересечений без кусочка и с тремя кусочками");
    check(n2 == 864 && n1 == 43 * 43 - 864, "в перекрытии 864 пересечения (16 нитей у швов x 27, по обеим осям)");
}

static void test_stitches()
{
    std::printf("5 test_stitches\n");
    const ChessParams P = params(2, ChessMode::Flip);
    const ChessFabric C = fiber::build_plain_weave_chess(P);
    check(C.stitches_per_layer == 864 && C.stitch_pairs.size() == 2 * 864, "864 стяжки на слой");
    std::vector<std::array<double, 3>> X;
    for (const auto& L : C.plies)
        for (const auto& F : L.fibers)
            for (const auto& p : F.pts) X.push_back(p);
    bool ok = true;
    const double len = P.sub_pitch - 2.0 * P.base.amp;
    for (const auto& pr : C.stitch_pairs)
    {
        const auto &a = X[pr.first], &b = X[pr.second];
        if (a[0] != b[0] || a[1] != b[1]) ok = false;
        if (std::fabs((a[2] - b[2]) - len) > 1e-12) ok = false;
    }
    check(ok, "стяжки вертикальны, сверху вниз, длина D - 2A = d + eps");
}

static void test_clamp()
{
    std::printf("6 test_clamp\n");
    const ChessParams P = params(1, ChessMode::Flip);
    const ChessFabric C = fiber::build_plain_weave_chess(P);
    const fiber::PlainWeave W0 = fiber::build_plain_weave(P.base);
    std::size_t ref = W0.fibers[static_cast<std::size_t>(P.base.nx)].n_clamp_lo; // выступ одного конца
    bool ok = true;
    std::size_t cl = 0, cfib = 0;
    for (const auto& F : C.plies[0].fibers)
    {
        if (F.direction == 'x' && (F.n_clamp_lo + F.n_clamp_hi) != 0) ok = false;
        if (F.direction == 'y' && F.clamped)
        {
            ++cfib;
            if (!(F.n_clamp_lo == ref || F.n_clamp_hi == ref)) ok = false;
        }
        cl += F.n_clamp_lo + F.n_clamp_hi;
    }
    // у края: каждая из 43 нитей утка -- в одном кусочке, у швов (16 нитей) -- в двух; два края
    check(ok, "зажаты только выступы нитей утка у края");
    check(cfib == 2 * (43 + 16), "зажатых нитей-кусочков 2 * (43 + 16) = 118");
    check(cl == cfib * ref, "зажатых бусин -- по выступу на нить");
    check(C.plies[0].clamp_boxes.size() == 2, "два бокса зажима");
}

// Наименьшее расстояние между осями разных нитей: точки хорд (по 4 на хорду), сетка ячеек.
static double min_axis_distance(const ChessFabric& C, double rc, double box)
{
    struct Pt { double x, y, z; int f; };
    std::vector<Pt> pts;
    int f = 0;
    for (const auto& L : C.plies)
        for (const auto& F : L.fibers)
        {
            for (std::size_t t = 0; t + 1 < F.pts.size(); ++t)
                for (int s = 0; s < 4; ++s)
                {
                    const double u = s / 4.0;
                    const auto &a = F.pts[t], &b = F.pts[t + 1];
                    const Pt p{a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1]), a[2] + u * (b[2] - a[2]), f};
                    if (std::fabs(p.x) <= box && std::fabs(p.y) <= box) pts.push_back(p);
                }
            ++f;
        }
    std::map<std::tuple<long, long, long>, std::vector<int>> grid;
    auto key = [&](const Pt& p) {
        return std::make_tuple(static_cast<long>(std::floor(p.x / rc)), static_cast<long>(std::floor(p.y / rc)),
                               static_cast<long>(std::floor(p.z / rc)));
    };
    for (int i = 0; i < static_cast<int>(pts.size()); ++i) grid[key(pts[static_cast<std::size_t>(i)])].push_back(i);
    double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i < static_cast<int>(pts.size()); ++i)
    {
        const Pt& p = pts[static_cast<std::size_t>(i)];
        const auto k = key(p);
        for (long dx = -1; dx <= 1; ++dx)
            for (long dy = -1; dy <= 1; ++dy)
                for (long dz = -1; dz <= 1; ++dz)
                {
                    const auto it = grid.find({std::get<0>(k) + dx, std::get<1>(k) + dy, std::get<2>(k) + dz});
                    if (it == grid.end()) continue;
                    for (int j : it->second)
                    {
                        const Pt& q = pts[static_cast<std::size_t>(j)];
                        if (q.f == p.f) continue;
                        best = std::min(best, std::sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) +
                                                        (p.z - q.z) * (p.z - q.z)));
                    }
                }
    }
    return best;
}

static void test_distance()
{
    std::printf("7 test_distance\n");
    for (ChessMode m : {ChessMode::Flip, ChessMode::Same})
    {
        const ChessParams P = params(2, m);
        const ChessFabric C = fiber::build_plain_weave_chess(P);
        const double dmin = min_axis_distance(C, P.base.d, 10e-3);
        char buf[160];
        std::snprintf(buf, sizeof buf, "%s: наименьшее расстояние осей разных нитей %.4f мм >= d = %.2f мм",
                      fiber::chess_mode_name(m), dmin * 1e3, P.base.d * 1e3);
        check(dmin >= P.base.d * (1.0 - 1e-9), buf);
    }
}

static void test_errors()
{
    std::printf("8 test_errors\n");
    bool thrown = false;
    try
    {
        ChessParams P = params(1, ChessMode::Flip);
        P.base.seg = 3;
        (void)fiber::build_plain_weave_chess(P);
    }
    catch (const std::invalid_argument&) { thrown = true; }
    check(thrown, "нечётный seg -- invalid_argument");
    thrown = false;
    try
    {
        ChessParams P = params(1, ChessMode::Flip);
        P.sub_pitch = 2.0 * P.base.amp + 0.5 * P.base.d;
        (void)fiber::build_plain_weave_chess(P);
    }
    catch (const std::invalid_argument&) { thrown = true; }
    check(thrown, "шаг подслоёв < 2A + d -- invalid_argument");
    thrown = false;
    try
    {
        ChessParams P = params(2, ChessMode::Flip);
        P.shift = 5;
        (void)fiber::build_plain_weave_chess(P);
    }
    catch (const std::invalid_argument&) { thrown = true; }
    check(thrown, "shift = patch -- invalid_argument");
}

static void test_shift()
{
    std::printf("9 test_shift\n");
    for (ChessMode m : {ChessMode::Flip, ChessMode::Same})
    {
        ChessParams P = params(2, m);
        P.shift = 2;
        const ChessFabric C = fiber::build_plain_weave_chess(P);
        const std::string tag = std::string(fiber::chess_mode_name(m)) + ": ";
        if (m == ChessMode::Flip)
        {
            const fiber::ChessAxis& a = C.ax_odd;
            check(a.count() == 9 && a.core.front() == std::make_pair(0, 5) && a.core.back() == std::make_pair(41, 42),
                  "нечётный слой: 9 кусочков, края 6 и 2 нити");
            check(a.core[a.center] == std::make_pair(21, 25), "нечётный слой: центральный кусочек -- нити 21..25");
            bool apart = true;
            for (std::size_t t = 0; t + 1 < a.core.size(); ++t)
                for (std::size_t u = 0; u + 1 < C.ax.core.size(); ++u)
                    if (a.core[t].second == C.ax.core[u].second) apart = false;
            check(apart, "линии шва слоёв 1 и 2 не совпадают (сдвиг 2 нити)");
            check(C.ay_odd.core == a.core, "по y сдвиг тот же");
            int n2 = 0, bad = 0;
            for (int i = 0; i < 43; ++i)
                for (int j = 0; j < 43; ++j)
                {
                    int n = 0;
                    for (int mx = 0; mx < a.count(); ++mx)
                        for (int my = 0; my < C.ay_odd.count(); ++my)
                            if (fiber::plain_weave_chess_detail::owns(a, C.ay_odd, mx, my, i, j)) ++n;
                    if (n == 2) ++n2;
                    else if (n != 1) ++bad;
                }
            check(bad == 0 && static_cast<std::size_t>(n2) == C.stitches_per_layer_odd,
                  "нечётный слой: пересечения у 1 или 2 кусочков, стяжек = пересечений перекрытия (" +
                      std::to_string(C.stitches_per_layer_odd) + ")");
        }
        check(C.stitch_pairs.size() == C.stitches_per_layer + C.stitches_per_layer_odd, tag + "стяжек -- сумма по слоям");
        std::vector<std::array<double, 3>> X;
        for (const auto& L : C.plies)
            for (const auto& F : L.fibers)
                for (const auto& p : F.pts) X.push_back(p);
        bool ok = true;
        const double len = P.sub_pitch - 2.0 * P.base.amp;
        for (const auto& pr : C.stitch_pairs)
        {
            const auto &a = X[pr.first], &b = X[pr.second];
            if (a[0] != b[0] || a[1] != b[1] || std::fabs((a[2] - b[2]) - len) > 1e-12) ok = false;
        }
        check(ok, tag + "стяжки вертикальны, длина d + eps");
        const double dmin = min_axis_distance(C, P.base.d, 10e-3);
        char buf[160];
        std::snprintf(buf, sizeof buf, "%sнаименьшее расстояние осей разных нитей %.4f мм >= d", tag.c_str(), dmin * 1e3);
        check(dmin >= P.base.d * (1.0 - 1e-9), buf);
    }
}

static void test_through()
{
    std::printf("10 test_through\n");
    for (ChessMode m : {ChessMode::Flip, ChessMode::Same})
        for (int sh : {0, 2})
        {
            ChessParams P = params(3, m);
            P.shift = sh;
            P.through = true;
            const ChessFabric C = fiber::build_plain_weave_chess(P);
            ChessParams P0 = P;
            P0.through = false;
            const ChessFabric C0 = fiber::build_plain_weave_chess(P0);
            char tag[64];
            std::snprintf(tag, sizeof tag, "%s, shift %d: ", fiber::chess_mode_name(m), sh);
            // перебор: пересечения, где у слоя k кусочек на 2k + 1, а у k + 1 -- на 2k + 2
            std::size_t want = 0;
            for (int k = 0; k + 1 < 3; ++k)
                for (int i = 0; i < 43; ++i)
                    for (int j = 0; j < 43; ++j)
                    {
                        bool up = false, dn = false;
                        for (int t = 0; t < 2; ++t)
                        {
                            const int kk = k + t;
                            const auto &ax = C.axis_x(kk), &ay = C.axis_y(kk);
                            for (int mx = 0; mx < ax.count(); ++mx)
                                for (int my = 0; my < ay.count(); ++my)
                                    if (fiber::plain_weave_chess_detail::owns(ax, ay, mx, my, i, j) &&
                                        fiber::chess_sublayer(C, kk, mx, my) == 2 * k + 1 + t)
                                        (t == 0 ? up : dn) = true;
                        }
                        if (up && dn) ++want;
                    }
            check(C.n_in == C0.stitch_pairs.size() && C.n_through == want && want > 0,
                  std::string(tag) + "стяжек между слоями " + std::to_string(C.n_through) + " (перебор " +
                      std::to_string(want) + "), в слоях как без них");
            std::vector<std::array<double, 3>> X;
            for (const auto& L : C.plies)
                for (const auto& F : L.fibers)
                    for (const auto& p : F.pts) X.push_back(p);
            bool ok = true;
            const double len = P.sub_pitch - 2.0 * P.base.amp;
            for (std::size_t t = C.n_in; t < C.stitch_pairs.size(); ++t)
            {
                const auto &a = X[C.stitch_pairs[t].first], &b = X[C.stitch_pairs[t].second];
                const int k = C.stitch_cross[t][0];
                const double za = C.sub_cz[static_cast<std::size_t>(2 * k + 1)] - P.base.amp;
                const double zb = C.sub_cz[static_cast<std::size_t>(2 * k + 2)] + P.base.amp;
                if (a[0] != b[0] || a[1] != b[1] || std::fabs((a[2] - b[2]) - len) > 1e-12 ||
                    std::fabs(a[2] - za) > 1e-9 || std::fabs(b[2] - zb) > 1e-9)
                    ok = false;
            }
            check(ok, std::string(tag) + "вертикальны, длина d + eps, от подслоя 2k + 1 к 2k + 2");
        }
}

int main()
{
    test_axes();
    test_beads();
    test_sublayers();
    test_owners();
    test_stitches();
    test_clamp();
    test_distance();
    test_shift();
    test_through();
    test_errors();
    std::printf("%s: %d ошибок\n", failures == 0 ? "OK" : "FAIL", failures);
    return failures == 0 ? 0 : 1;
}
