// tests/test_plain_weave_random_cut.cpp
// Тесты ткани «случайная раскройка» (include/fiber/plain_weave_random_cut.hpp,
// fiber/tz-random-cut-quilt.md, разделы 4 и 6), включая многослойный пакет.
//
// Разделы:
//   1 test_determinism   -- один seed -> побитово то же; другой seed -> другое;
//   2 test_areas         -- площади рабочих многоугольников + зазор x длины швов = W*H;
//   3 test_inside        -- бусины внутри рабочего многоугольника своего фрагмента;
//   4 test_gap           -- зазор между фрагментами; gap < d -> исключение / предупреждение;
//   5 test_crossings     -- |z_основы - z_утка| = 2A в пересечениях (Cos, Arcs, theta = 0);
//   6 test_l0            -- хорды внутри нитей равны по длине (1 %), концевые >= 0.2 h;
//   7 test_stitches      -- стежки в слое: число, пары из разных фрагментов, stitch_step = 0;
//   8 test_fib           -- .fib читается load_fiber_file: бусины, связи, материал 1;
//   9 test_single        -- N = 1;
//  10 test_stack         -- пакет слоёв: раскройки, копии, сквозные стежки, .fib;
//  11 test_clamp         -- зажим;
//  12 test_json          -- .cut.json и контакты диагностики.
#include "fiber/fiber_input.hpp"
#include "fiber/model.hpp"
#include "fiber/plain_weave.hpp"
#include "fiber/plain_weave_random_cut.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using fiber::PlainWeaveProfile;
using fiber::RandomCut;
using fiber::RandomCutParams;

static int failures = 0;

static void check(bool cond, const std::string& what)
{
    std::printf("  [%s] %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

// Базовая ткань: twaron_ct717 без зажима, N = 9, 50 x 50 мм, seed 3.
static RandomCutParams base_params(int N = 9, std::uint64_t seed = 3, PlainWeaveProfile prof = PlainWeaveProfile::Cos)
{
    RandomCutParams P;
    fiber::apply_plain_weave_preset("twaron_ct717", P.base);
    P.base.clamp = fiber::PlainWeaveClamp::None;
    P.base.clamp_edge_set = false;
    P.base.profile = prof;
    P.n_frag = N;
    P.size_x = 0.05;
    P.size_y = 0.05;
    P.seed = seed;
    return P;
}

// Исключение заданного типа при построении.
template <class Ex>
static bool throws(const RandomCutParams& P)
{
    try
    {
        fiber::build_plain_weave_random_cut(P);
    }
    catch (const Ex&)
    {
        return true;
    }
    catch (...)
    {
    }
    return false;
}

// Плоский массив бусин в порядке записи .fib.
struct Flat
{
    std::vector<std::array<double, 3>> pts;
    std::vector<int> ply, frag; // [бусина] -> слой, фрагмент
};

static Flat flatten(const RandomCut& R)
{
    Flat F;
    for (std::size_t f = 0; f < R.W.fibers.size(); ++f)
        for (const auto& X : R.W.fibers[f].pts)
        {
            F.pts.push_back(X);
            F.ply.push_back(R.fiber_ply[f]);
            F.frag.push_back(R.fiber_fragment[f]);
        }
    return F;
}

static std::size_t total_beads(const RandomCut& R)
{
    std::size_t n = 0;
    for (const auto& F : R.W.fibers) n += F.pts.size();
    return n;
}

static double chord(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

static double plan_dist(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return std::hypot(a[0] - b[0], a[1] - b[1]);
}

static double shoelace(const std::vector<fiber::RcPoint>& P)
{
    double s = 0.0;
    for (std::size_t q = 0; q < P.size(); ++q)
    {
        const auto& a = P[q];
        const auto& b = P[(q + 1) % P.size()];
        s += a[0] * b[1] - b[0] * a[1];
    }
    return 0.5 * std::fabs(s);
}

static bool same_fibers(const RandomCut& A, const RandomCut& B)
{
    if (A.W.fibers.size() != B.W.fibers.size()) return false;
    for (std::size_t f = 0; f < A.W.fibers.size(); ++f)
    {
        if (A.W.fibers[f].pts.size() != B.W.fibers[f].pts.size()) return false;
        for (std::size_t b = 0; b < A.W.fibers[f].pts.size(); ++b)
            if (A.W.fibers[f].pts[b] != B.W.fibers[f].pts[b]) return false;
    }
    return true;
}

static std::size_t seam_points(const RandomCut& R)
{
    std::size_t n = 0;
    for (const auto& c : R.cuts)
        for (const auto& S : c.seams) n += S.points.size();
    return n;
}

static const char* prof_name(PlainWeaveProfile p) { return p == PlainWeaveProfile::Cos ? "cos" : "arcs"; }

// --- 1. Детерминизм --------------------------------------------------------
static void test_determinism()
{
    std::printf("-- 1. детерминизм\n");
    const RandomCut A = fiber::build_plain_weave_random_cut(base_params());
    const RandomCut B = fiber::build_plain_weave_random_cut(base_params());
    check(same_fibers(A, B), "один seed: pts всех нитей совпадают побитово");
    check(A.stitch_pairs == B.stitch_pairs && !A.stitch_pairs.empty(), "один seed: stitch_pairs совпадают (и не пусты)");
    const RandomCut C = fiber::build_plain_weave_random_cut(base_params(9, 4));
    check(!same_fibers(A, C), "другой seed: нити другие");
    check(A.stitch_pairs != C.stitch_pairs, "другой seed: stitch_pairs другие");
}

// --- 2. Площади ------------------------------------------------------------
static void test_areas()
{
    std::printf("-- 2. площади\n");
    for (int N : {2, 9, 25})
    {
        const RandomCutParams P = base_params(N);
        const RandomCut R = fiber::build_plain_weave_random_cut(P);
        const double WH = P.size_x * P.size_y;
        double cells = 0.0;
        for (const auto& F : R.cuts[0].frags) cells += shoelace(F.cell);
        const double area = fiber::random_cut_area(R, 0);
        const double seams = fiber::random_cut_seam_length(R, 0);
        const double sum = area + R.params.gap * seams;
        std::printf("    N=%d: сумма ячеек %.9g, площадь+зазор*швы %.9g, W*H %.9g\n", N, cells, sum, WH);
        check(std::fabs(cells - WH) < 1e-9 * WH, "N=" + std::to_string(N) + ": сумма площадей ячеек = W*H (1e-9)");
        check(std::fabs(sum - WH) < 0.01 * WH, "N=" + std::to_string(N) + ": F.area + gap * длины швов = W*H (1 %)");
    }
}

// --- 3. Бусины внутри рабочего многоугольника ------------------------------
static void test_inside()
{
    std::printf("-- 3. бусины внутри рабочего многоугольника\n");
    for (auto prof : {PlainWeaveProfile::Cos, PlainWeaveProfile::Arcs})
    {
        const RandomCut R = fiber::build_plain_weave_random_cut(base_params(9, 3, prof));
        const Flat F = flatten(R);
        std::size_t bad = 0;
        double worst = -1e300;
        for (std::size_t b = 0; b < F.pts.size(); ++b)
        {
            const auto& frag = R.cuts[static_cast<std::size_t>(R.ply_cut[static_cast<std::size_t>(F.ply[b])])]
                                   .frags[static_cast<std::size_t>(F.frag[b])];
            for (const auto& h : frag.hps)
            {
                const double v = h.nx * F.pts[b][0] + h.ny * F.pts[b][1] - h.c;
                worst = std::max(worst, v);
                if (v > 1e-12) ++bad;
            }
        }
        std::printf("    макс. n.x - c = %.3e м\n", worst);
        check(bad == 0 && !F.pts.empty(), std::string(prof_name(prof)) + ": все бусины удовлетворяют hps (<= c + 1e-12)");
    }
}

// --- 4. Зазор --------------------------------------------------------------
static void test_gap()
{
    std::printf("-- 4. зазор между фрагментами\n");
    for (auto prof : {PlainWeaveProfile::Cos, PlainWeaveProfile::Arcs})
    {
        const RandomCutParams P = base_params(9, 3, prof);
        const RandomCut R = fiber::build_plain_weave_random_cut(P);
        const Flat F = flatten(R);
        // Перебор пар разных фрагментов одного слоя с отсечением по x (сортировка).
        std::vector<std::size_t> ord(F.pts.size());
        for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return F.pts[a][0] < F.pts[b][0]; });
        double dmin = 1e300;
        for (std::size_t ia = 0; ia < ord.size(); ++ia)
            for (std::size_t ib = ia + 1; ib < ord.size(); ++ib)
            {
                const std::size_t a = ord[ia], b = ord[ib];
                if (F.pts[b][0] - F.pts[a][0] > R.params.gap) break;
                if (F.ply[a] != F.ply[b] || F.frag[a] == F.frag[b]) continue;
                dmin = std::min(dmin, plan_dist(F.pts[a], F.pts[b]));
            }
        std::printf("    %s: мин. расстояние %.9g м (gap %.9g м, ядро: %.9g м)\n", prof_name(prof), dmin, R.params.gap,
                    R.min_cross_dist);
        check(dmin >= R.params.gap - 1e-9 * P.base.pitch,
              std::string(prof_name(prof)) + ": расстояние в плане между фрагментами >= gap - 1e-9 p");
        check(std::fabs(R.params.gap - P.base.d) < 1e-15, "gap по умолчанию = d");
    }
    RandomCutParams P = base_params();
    P.gap = 0.5 * P.base.d;
    check(throws<std::invalid_argument>(P), "gap < d без allow_penetration -> std::invalid_argument");
    P.base.allow_penetration = true;
    bool ok = true;
    RandomCut R;
    try
    {
        R = fiber::build_plain_weave_random_cut(P);
    }
    catch (const std::exception& e)
    {
        std::printf("    исключение: %s\n", e.what());
        ok = false;
    }
    check(ok, "с allow_penetration строится");
    check(ok && !R.W.warnings.empty(), "W.warnings непуст");
}

// --- 5. 2A в пересечениях --------------------------------------------------
static void crossing_case(const RandomCutParams& P, const std::string& name)
{
    const RandomCut R = fiber::build_plain_weave_random_cut(P);
    const double A = P.base.amp, p = P.base.pitch;
    double worst = 0.0;
    std::size_t matched = 0, unmatched_frag = 0;
    for (const auto& rg : R.ranges)
    {
        if (rg.ply != 0) continue;
        std::size_t m_frag = 0;
        for (std::size_t fa = rg.fiber_first; fa < rg.fiber_first + rg.fiber_count; ++fa)
        {
            const auto& Wa = R.W.fibers[fa];
            if (Wa.direction != 'x') continue;
            for (std::size_t fb = rg.fiber_first; fb < rg.fiber_first + rg.fiber_count; ++fb)
            {
                const auto& Wb = R.W.fibers[fb];
                if (Wb.direction != 'y') continue;
                for (std::size_t va : Wa.vertex_bead)
                    for (std::size_t vb : Wb.vertex_bead)
                    {
                        const auto &a = Wa.pts[va], &b = Wb.pts[vb];
                        if (plan_dist(a, b) > 1e-9 * p) continue;
                        ++matched;
                        ++m_frag;
                        worst = std::max(worst, std::fabs(std::fabs(a[2] - b[2]) - 2.0 * A));
                    }
            }
        }
        if (m_frag == 0) ++unmatched_frag;
    }
    std::printf("    %s: пересечений %zu, макс. ||dz| - 2A| = %.3e м\n", name.c_str(), matched, worst);
    check(matched > 0, name + ": пересечения найдены");
    check(unmatched_frag == 0, name + ": в каждом фрагменте есть пересечения");
    check(worst <= 1e-9 * A, name + ": в каждом пересечении |z_основы - z_утка| = 2A (1e-9 A)");
    check(2.0 * A >= P.base.d * (1.0 - 1e-12), name + ": 2A >= d (нет проникания в узле)");
}

static void test_crossings()
{
    std::printf("-- 5. 2A в пересечениях\n");
    crossing_case(base_params(9, 3, PlainWeaveProfile::Cos), "cos, N=9");
    crossing_case(base_params(9, 3, PlainWeaveProfile::Arcs), "arcs, N=9");
    crossing_case(base_params(5, 11, PlainWeaveProfile::Cos), "cos, N=5, seed 11");
    // theta = 0: сетка вдоль осей, положение нитей по фазе.
    RandomCutParams P = base_params(9, 3, PlainWeaveProfile::Cos);
    P.rot_max_deg = 0.0;
    crossing_case(P, "cos, rot_max=0");
    P.base.profile = PlainWeaveProfile::Arcs;
    crossing_case(P, "arcs, rot_max=0");
    const RandomCut R = fiber::build_plain_weave_random_cut(P);
    const double p = P.base.pitch;
    bool theta0 = true, pos_ok = true;
    for (const auto& F : R.cuts[0].frags) theta0 = theta0 && F.theta == 0.0;
    check(theta0, "rot_max=0: theta = 0 у всех фрагментов");
    for (std::size_t f = 0; f < R.W.fibers.size(); ++f)
    {
        const auto& Fb = R.W.fibers[f];
        const auto& frag = R.cuts[0].frags[static_cast<std::size_t>(R.fiber_fragment[f])];
        for (const auto& X : Fb.pts)
        {
            // Основа: y = c_y + (j + fx) p; уток: x = c_x + (i + fy) p.
            const double expect = Fb.direction == 'x' ? frag.center[1] + (Fb.index + frag.fx) * p
                                                      : frag.center[0] + (Fb.index + frag.fy) * p;
            const double got = Fb.direction == 'x' ? X[1] : X[0];
            if (std::fabs(got - expect) > 1e-12) pos_ok = false;
        }
    }
    check(pos_ok, "rot_max=0: основа на y = c_y + (j + fx) p, уток на x = c_x + (i + fy) p");
}

// --- 6. l0 -----------------------------------------------------------------
static void test_l0()
{
    std::printf("-- 6. длины хорд\n");
    for (auto prof : {PlainWeaveProfile::Cos, PlainWeaveProfile::Arcs})
    {
        const RandomCutParams P = base_params(9, 3, prof);
        const RandomCut R = fiber::build_plain_weave_random_cut(P);
        const double h = P.base.pitch / P.base.seg;
        double nominal = -1.0, worst_rel = 0.0, min_end = 1e300;
        std::size_t inner = 0, ends = 0;
        for (const auto& F : R.W.fibers)
        {
            const std::size_t n = F.pts.size();
            if (n < 2) continue;
            if (n >= 4)
            {
                if (nominal < 0.0) nominal = chord(F.pts[1], F.pts[2]);
                for (std::size_t b = 1; b + 2 < n; ++b)
                {
                    worst_rel = std::max(worst_rel, std::fabs(chord(F.pts[b], F.pts[b + 1]) - nominal) / nominal);
                    ++inner;
                }
            }
            min_end = std::min({min_end, chord(F.pts[0], F.pts[1]), chord(F.pts[n - 2], F.pts[n - 1])});
            ends += 2;
        }
        std::printf("    %s: номинал %.9g м, внутренних хорд %zu, макс. отклонение %.3e, мин. концевая %.4g h\n",
                    prof_name(prof), nominal, inner, worst_rel, min_end / h);
        const std::string nm = prof_name(prof);
        check(inner > 0 && worst_rel < 0.01, nm + ": внутренние хорды отличаются от номинала < 1 %");
        check(ends > 0 && min_end >= 0.2 * h - 1e-12, nm + ": концевые хорды >= 0.2 h");
    }
}

// --- 7. Стежки в слое ------------------------------------------------------
static void test_stitches()
{
    std::printf("-- 7. стежки в слое\n");
    const RandomCutParams P = base_params();
    const RandomCut R = fiber::build_plain_weave_random_cut(P);
    const Flat F = flatten(R);
    const std::size_t sp = seam_points(R);
    std::printf("    точек швов %zu, стежков %zu (в слоях %zu, сквозных %zu)\n", sp, R.stitches.size(), R.n_stitches_in,
                R.n_stitches_through);
    check(R.stitches.size() > 0 && R.stitches.size() <= sp, "0 < стежков <= число точек швов");
    check(R.n_stitches_in == R.stitches.size() && R.n_stitches_through == 0, "L = 1: все стежки в слое, сквозных нет");
    check(R.stitch_pairs.size() == R.stitches.size(), "stitch_pairs.size() = stitches.size()");
    std::set<std::pair<std::size_t, std::size_t>> uniq;
    bool diff_frag = true, close = true, in_range = true, same_ply = true, flags = true;
    const double lim = R.params.gap + 2.0 * P.base.pitch;
    double dmax = 0.0;
    for (const auto& S : R.stitches)
    {
        uniq.insert({std::min(S.a, S.b), std::max(S.a, S.b)});
        if (S.a >= F.pts.size() || S.b >= F.pts.size())
        {
            in_range = false;
            continue;
        }
        if (F.frag[S.a] == F.frag[S.b]) diff_frag = false;
        if (F.ply[S.a] != F.ply[S.b]) same_ply = false;
        if (F.frag[S.a] != S.frag_a || F.frag[S.b] != S.frag_b || S.through) flags = false;
        const double dd = plan_dist(F.pts[S.a], F.pts[S.b]);
        dmax = std::max(dmax, dd);
        if (dd > lim) close = false;
    }
    std::printf("    макс. расстояние в паре в плане %.9g м (предел %.9g м)\n", dmax, lim);
    check(in_range, "индексы стежков в пределах числа бусин");
    check(uniq.size() == R.stitches.size(), "пар-дубликатов нет");
    check(diff_frag, "бусины каждой пары из разных фрагментов (нет стежка внутри фрагмента)");
    check(same_ply && flags, "обе бусины в одном слое, frag_a/frag_b/through согласованы");
    check(close, "расстояние в паре в плане <= gap + 2 p");
    RandomCutParams P0 = P;
    P0.stitch_step = 0.0;
    const RandomCut R0 = fiber::build_plain_weave_random_cut(P0);
    check(R0.stitches.empty() && R0.stitch_pairs.empty(), "stitch_step = 0: стежков нет");
    check(same_fibers(R, R0), "stitch_step = 0: нити те же");
}

// --- 8. .fib ---------------------------------------------------------------
static void load_back(const RandomCut& R, const char* name, fiber::Model& model)
{
    namespace fs = std::filesystem;
    const fs::path fib = fs::temp_directory_path() / (std::string("test_pw_random_cut_") + name + ".fib");
    {
        std::ofstream os(fib);
        fiber::write_plain_weave_fib(os, R.W, "test_plain_weave_random_cut", fiber::random_cut_fib_extra(R));
    }
    fiber::load_fiber_file(fib.string(), model);
    std::error_code ec;
    fs::remove(fib, ec);
}

static void test_fib()
{
    std::printf("-- 8. .fib читается load_fiber_file\n");
    const RandomCut R = fiber::build_plain_weave_random_cut(base_params());
    fiber::Model model;
    load_back(R, "single", model);
    const std::size_t beads = total_beads(R);
    std::printf("    бусин %zu, нитей %zu, стежков %zu, связей в модели %zu\n", beads, R.W.fibers.size(),
                R.stitches.size(), model.topology.bonds.size());
    check(model.fiber.n == beads, "число бусин модели = сумма pts");
    check(model.topology.bonds.size() == beads - R.W.fibers.size() + R.stitches.size(),
          "число связей = (бусины - нити) + стежки");
    check(model.materials.size() >= 2, "материалов не меньше двух");
    if (model.materials.size() >= 2) check(model.materials[1].E == R.params.stitch_E, "материал 1 присутствует, E = stitch_E");
    std::size_t n_st = 0;
    for (const auto& b : model.topology.bonds)
        if (b.mat == 1) ++n_st;
    check(n_st == R.stitches.size(), "BOND материала 1 -- по числу стежков");
}

// --- 9. N = 1 --------------------------------------------------------------
static void test_single()
{
    std::printf("-- 9. N = 1\n");
    const RandomCut R = fiber::build_plain_weave_random_cut(base_params(1));
    check(R.stitches.empty() && R.stitch_pairs.empty(), "стежков нет");
    check(R.cuts.size() == 1 && R.cuts[0].frags.size() == 1 && R.cuts[0].seams.empty(),
          "одна раскройка, один фрагмент, швов нет");
    check(R.ranges.size() == 1 && !R.W.fibers.empty(), "один диапазон, нити построены");
    fiber::Model model;
    load_back(R, "n1", model);
    check(model.fiber.n == total_beads(R), "число бусин модели = сумма pts");
    check(model.topology.bonds.size() == total_beads(R) - R.W.fibers.size(), "связи = бусины - нити (стежков нет)");
    const double WH = 0.05 * 0.05;
    check(std::fabs(R.cuts[0].frags[0].area - WH) < 1e-9 * WH, "площадь единственного фрагмента = W*H");
}

// --- 10. Пакет -------------------------------------------------------------
static void stack_checks(const RandomCut& R, const std::string& name, bool single_cut)
{
    const Flat F = flatten(R);
    const int L = R.params.layers;
    const std::size_t nb = total_beads(R);
    std::set<std::pair<std::size_t, std::size_t>> uniq;
    bool adj = true, rng = true, vertical = true;
    std::size_t thr = 0;
    for (const auto& S : R.stitches)
    {
        uniq.insert({std::min(S.a, S.b), std::max(S.a, S.b)});
        if (S.a >= nb || S.b >= nb)
        {
            rng = false;
            continue;
        }
        if (S.through)
        {
            ++thr;
            if (std::abs(S.ply_a - S.ply_b) != 1) adj = false;
            if (F.ply[S.a] != S.ply_a || F.ply[S.b] != S.ply_b) rng = false;
            if (S.ply_a < 0 || S.ply_b >= L) rng = false;
            if (single_cut && (std::fabs(F.pts[S.a][0] - F.pts[S.b][0]) > 1e-15 || std::fabs(F.pts[S.a][1] - F.pts[S.b][1]) > 1e-15))
                vertical = false;
        }
        else if (F.ply[S.a] != F.ply[S.b] || S.ply_a != S.ply_b)
            rng = false;
    }
    check(thr == R.n_stitches_through && thr + R.n_stitches_in == R.stitches.size(), name + ": счётчики стежков согласованы");
    check(adj, name + ": сквозные стежки только между слоями (l, l+1)");
    check(rng, name + ": бусина a в диапазоне слоя ply_a, b -- слоя ply_b, индексы в пределах");
    check(uniq.size() == R.stitches.size(), name + ": пар-дубликатов нет");
    check(R.n_stitches_through <= static_cast<std::size_t>(L - 1) * seam_points(R),
          name + ": сквозных <= (L-1) * точек швов всех раскроек");
    if (single_cut)
        check(vertical && thr > 0, name + ": одна раскройка -- сквозные стежки строго вертикальны");
    else
        check(thr > 0, name + ": сквозные стежки есть");
    fiber::Model model;
    load_back(R, "stack", model);
    check(model.fiber.n == nb, name + ": .fib: число бусин модели = сумма pts");
    check(model.topology.bonds.size() == nb - R.W.fibers.size() + R.stitches.size(),
          name + ": .fib: связи = (бусины - нити) + stitches.size()");
}

static void test_stack()
{
    std::printf("-- 10. пакет слоёв\n");
    RandomCutParams P = base_params();
    P.layers = 4;
    P.layers_per_cut = 2;
    const RandomCut R = fiber::build_plain_weave_random_cut(P);
    const RandomCut S = fiber::build_plain_weave_random_cut(base_params());
    check((R.ply_cut == std::vector<int>{0, 0, 1, 1}) && R.cuts.size() == 2, "L=4, k=2: ply_cut = {0,0,1,1}, две раскройки");
    // Слои 0 и 1: копии со сдвигом z.
    check(R.ply_bead_first.size() == 4, "ply_bead_first по числу слоёв");
    const std::size_t n0 = R.ply_bead_first[1] - R.ply_bead_first[0];
    check(R.ply_bead_first[2] - R.ply_bead_first[1] == n0, "слои 0 и 1 содержат поровну бусин");
    const Flat F = flatten(R);
    bool xy_eq = true, dz_ok = true;
    double dz_err = 0.0;
    for (std::size_t b = 0; b < n0; ++b)
    {
        const auto &a = F.pts[R.ply_bead_first[0] + b], &c = F.pts[R.ply_bead_first[1] + b];
        if (a[0] != c[0] || a[1] != c[1]) xy_eq = false;
        const double e = std::fabs((a[2] - c[2]) - R.params.ply_pitch);
        dz_err = std::max(dz_err, e);
        if (e > 1e-15) dz_ok = false;
    }
    std::printf("    макс. отклонение dz от ply_pitch %.3e м\n", dz_err);
    check(xy_eq, "xy бусин слоя 1 побитово равны xy слоя 0");
    check(dz_ok, "z слоя 1 отличается от слоя 0 ровно на ply_pitch (1e-15)");
    bool centers_differ = false;
    for (std::size_t f = 0; f < R.cuts[0].frags.size(); ++f)
        if (R.cuts[0].frags[f].center != R.cuts[1].frags[f].center) centers_differ = true;
    check(centers_differ, "слои 0 и 2 -- разные раскройки (центры отличаются)");
    // Раскройка 0 пакета = однослойная с тем же seed.
    bool cut_same = R.cuts[0].frags.size() == S.cuts[0].frags.size();
    for (std::size_t f = 0; cut_same && f < S.cuts[0].frags.size(); ++f)
        cut_same = R.cuts[0].frags[f].center == S.cuts[0].frags[f].center && R.cuts[0].frags[f].theta == S.cuts[0].frags[f].theta &&
                   R.cuts[0].frags[f].fx == S.cuts[0].frags[f].fx && R.cuts[0].frags[f].fy == S.cuts[0].frags[f].fy;
    check(cut_same, "cuts[0] пакета = раскройка однослойного запуска (центры, theta, фаза побитово)");
    const Flat Fs = flatten(S);
    bool beads_same = n0 == Fs.pts.size();
    for (std::size_t b = 0; beads_same && b < n0; ++b) beads_same = F.pts[b] == Fs.pts[b];
    check(beads_same, "бусины слоя 0 пакета = бусины однослойного запуска");
    check(R.ply_cz.size() == 4 && R.ply_cz[0] == P.base.cz, "ply_cz[0] = cz");
    stack_checks(R, "L=4,k=2", false);

    // k = L (0 и 4): одна раскройка, одинаковые слои.
    for (int k : {0, 4})
    {
        RandomCutParams Q = base_params();
        Q.layers = 4;
        Q.layers_per_cut = k;
        const RandomCut T = fiber::build_plain_weave_random_cut(Q);
        const std::string nm = "L=4,k=" + std::to_string(k);
        check(T.cuts.size() == 1 && T.ply_cut == std::vector<int>({0, 0, 0, 0}), nm + ": одна раскройка, ply_cut = {0,0,0,0}");
        const Flat FT = flatten(T);
        const std::size_t m = T.ply_bead_first[1];
        bool same = FT.pts.size() == 4 * m;
        for (std::size_t l = 1; same && l < 4; ++l)
            for (std::size_t b = 0; b < m; ++b)
                if (FT.pts[b][0] != FT.pts[l * m + b][0] || FT.pts[b][1] != FT.pts[l * m + b][1]) same = false;
        check(same, nm + ": все 4 слоя одинаковы в плане");
        stack_checks(T, nm, true);
    }
    // L = 3, k = 2.
    {
        RandomCutParams Q = base_params();
        Q.layers = 3;
        Q.layers_per_cut = 2;
        const RandomCut T = fiber::build_plain_weave_random_cut(Q);
        check(T.cuts.size() == 2 && T.ply_cut == std::vector<int>({0, 0, 1}), "L=3,k=2: две раскройки, ply_cut = {0,0,1}");
        stack_checks(T, "L=3,k=2", false);
    }
    // Ошибка параметров.
    {
        RandomCutParams Q = base_params();
        Q.layers = 4;
        Q.layers_per_cut = 5;
        check(throws<std::invalid_argument>(Q), "layers_per_cut = 5 при L = 4 -> std::invalid_argument");
    }
}

// --- 11. Зажим -------------------------------------------------------------
static void test_clamp()
{
    std::printf("-- 11. зажим\n");
    RandomCutParams P = base_params();
    P.base.clamp = fiber::PlainWeaveClamp::Y;
    P.base.clamp_edge_set = false;
    check(throws<std::invalid_argument>(P), "clamp = Y без clamp_edge -> std::invalid_argument");
    P.base.clamp_edge = 0.02;
    P.base.clamp_edge_set = true;
    bool ok = true;
    RandomCut R;
    try
    {
        R = fiber::build_plain_weave_random_cut(P);
    }
    catch (const std::exception& e)
    {
        std::printf("    исключение: %s\n", e.what());
        ok = false;
    }
    check(ok, "clamp_edge = 0.02 строится");
    if (ok)
    {
        std::size_t clamped = 0;
        for (const auto& F : R.W.fibers) clamped += F.n_clamp_lo + F.n_clamp_hi;
        std::printf("    CLAMP_BOX: %zu, зажатых бусин: %zu\n", R.W.clamp_boxes.size(), clamped);
        check(R.W.clamp_boxes.size() == 2, "2 CLAMP_BOX");
        check(clamped > 0, "зажатых бусин > 0");
        std::ostringstream os;
        fiber::write_plain_weave_fib(os, R.W, "cmd", fiber::random_cut_fib_extra(R));
        const std::string txt = os.str();
        std::size_t boxes = 0, pos = 0;
        while ((pos = txt.find("\nCLAMP_BOX", pos)) != std::string::npos)
        {
            ++boxes;
            ++pos;
        }
        check(boxes == 2, ".fib содержит 2 строки CLAMP_BOX");
    }
    P.base.clamp_edge = 0.1;
    check(throws<std::invalid_argument>(P), "clamp_edge = 0.1 (ничего не зажато) -> std::invalid_argument");
}

// --- 12. JSON и контакты ---------------------------------------------------
static void test_json()
{
    std::printf("-- 12. .cut.json и контакты\n");
    for (int L : {1, 2})
    {
        RandomCutParams P = base_params();
        P.layers = L;
        P.layers_per_cut = 1;
        const RandomCut R = fiber::build_plain_weave_random_cut(P);
        std::ostringstream os;
        fiber::write_random_cut_json(os, R);
        const std::string s = os.str();
        const std::string nm = "L=" + std::to_string(L) + ": ";
        check(!s.empty() && s[0] == '{', nm + "JSON начинается с '{'");
        check(s.find("\"cuts\"") != std::string::npos && s.find("\"stitches\"") != std::string::npos &&
                  s.find("\"contacts\"") != std::string::npos,
              nm + "есть ключи cuts, stitches, contacts");
        const Flat F = flatten(R);
        bool pair_ok = true;
        for (const auto& c : R.contacts)
        {
            if (c.first >= F.pts.size() || c.second >= F.pts.size())
            {
                pair_ok = false;
                continue;
            }
            if (F.ply[c.first] != F.ply[c.second] || F.frag[c.first] == F.frag[c.second]) pair_ok = false;
            if (chord(F.pts[c.first], F.pts[c.second]) >= R.params.rc) pair_ok = false;
        }
        std::printf("    контактов %zu\n", R.contacts.size());
        check(!R.contacts.empty(), nm + "контактов > 0");
        check(pair_ok, nm + "контакты: разные фрагменты одного слоя, расстояние < rc");
    }
}

int main()
{
    test_determinism();
    test_areas();
    test_inside();
    test_gap();
    test_crossings();
    test_l0();
    test_stitches();
    test_fib();
    test_single();
    test_stack();
    test_clamp();
    test_json();
    if (failures == 0)
        std::printf("test_plain_weave_random_cut: все проверки пройдены\n");
    else
        std::printf("test_plain_weave_random_cut: ПРОВАЛ, %d проверок\n", failures);
    return failures == 0 ? 0 : 1;
}
