// tests/test_plain_weave_quilt.cpp
// Тесты ткани-«одеяла» из патчей и сетки стёжки (include/fiber/plain_weave_quilt.hpp,
// fiber/plan-quilt-fabric.md, разделы 3-5 и 10-12).
//
// Разделы:
//   1 test_equivalence   -- одеяло 34x34 по 3x3 == одна ткань 102x102 (сцена Guo 2019), побитово и .fib;
//   2 test_replication   -- бусины патча (a, b) == бусины патча (0, 0) + сдвиг на размер патча;
//   3 test_no_gaps       -- хорды через стыки не отличаются от обычных, l0_min/l0_max ткани;
//   4 test_labels        -- метки бусин по патчам, стыковая бусина -> меньший патч;
//   5 test_odd_patch     -- нечётный патч: исключение без флага, предупреждение с флагом;
//   6 test_window        -- окно 3x3 по 34 нити = 119.2 x 119.2 мм;
//   7 test_weave         -- |z_основы - z_утка| = 2A в пересечениях, нет проникания;
//   8 test_patches_file  -- формат .patches.txt;
//   9 test_stitches      -- сетка стёжки: число, вертикальность, every, нечётный seg, ошибки.
#include "fiber/plain_weave.hpp"
#include "fiber/plain_weave_quilt.hpp"
#include "guo2019_case.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using fiber::PlainWeave;
using fiber::PlainWeaveParams;
using fiber::PlainWeaveProfile;
using fiber::PlainWeaveQuilt;
using fiber::QuiltParams;

static int failures = 0;

static void check(bool cond, const std::string& what)
{
    std::printf("  [%s] %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

// Параметры ткани сцены Guo 2019 (tests/guo2019_scene.hpp, блок «ткань»), уровень base.
static PlainWeaveParams guo_params(double cz)
{
    const guo2019::LevelSpec L = guo2019::level_spec(guo2019::Level::Base);
    PlainWeaveParams q;
    fiber::apply_plain_weave_preset("twaron_ct717", q);
    q.amp = guo2019::kFabricAmpOverD * q.d;
    q.seg = L.seg;
    q.nx = L.nx;
    q.ny = L.ny;
    q.overhang_x = 0.5 * q.pitch;
    q.overhang_y = L.e_y;
    q.clamp = fiber::PlainWeaveClamp::Y;
    q.clamp_edge = L.c_e;
    q.clamp_edge_set = true;
    q.E = guo2019::kFabricE;
    q.nu = guo2019::kFabricNu;
    q.rho = guo2019::kFabricRho;
    q.eps_break = guo2019::kFabricEpsBreak;
    q.mu_ff = guo2019::kMuDefault;
    q.mu_fs = guo2019::kMuDefault;
    q.mu_fs_set = true;
    q.cx = 0.0;
    q.cy = 0.0;
    q.cz = cz;
    return q;
}

static QuiltParams guo_quilt(double cz = -(0.6 * 0.39e-3 + 0.5 * 0.39e-3))
{
    QuiltParams Q;
    Q.base = guo_params(cz);
    Q.patch_nx = guo2019::kFabricPatchBase;
    Q.patch_ny = guo2019::kFabricPatchBase;
    Q.tile_x = guo2019::kFabricTileBase;
    Q.tile_y = guo2019::kFabricTileBase;
    return Q;
}

// Малая ткань для проверок, где нужен перебор пар нитей.
static QuiltParams small_quilt(int seg, PlainWeaveProfile prof, int pnx, int pny, int tx, int ty)
{
    QuiltParams Q;
    fiber::apply_plain_weave_preset("test", Q.base);
    Q.base.amp = 0.6 * Q.base.d;
    Q.base.seg = seg;
    Q.base.profile = prof;
    Q.patch_nx = pnx;
    Q.patch_ny = pny;
    Q.tile_x = tx;
    Q.tile_y = ty;
    return Q;
}

static bool same_pts(const PlainWeave& A, const PlainWeave& B)
{
    if (A.fibers.size() != B.fibers.size()) return false;
    for (std::size_t f = 0; f < A.fibers.size(); ++f)
    {
        if (A.fibers[f].pts.size() != B.fibers[f].pts.size()) return false;
        if (A.fibers[f].vertex_bead != B.fibers[f].vertex_bead) return false;
        for (std::size_t b = 0; b < A.fibers[f].pts.size(); ++b)
            if (A.fibers[f].pts[b] != B.fibers[f].pts[b]) return false;
    }
    return true;
}

static std::size_t total_beads(const PlainWeave& W)
{
    std::size_t n = 0;
    for (const auto& F : W.fibers) n += F.pts.size();
    return n;
}

static double chord(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

// --- 1. Эквивалентность целой ткани ---------------------------------------
static void test_equivalence()
{
    std::printf("-- 1. эквивалентность целой ткани (Guo 2019 base, 34x34 по 3x3)\n");
    for (int variant = 0; variant < 2; ++variant)
    {
        const double cz = variant == 0 ? -(0.6 * 0.39e-3 + 0.5 * 0.39e-3) : -0.0123;
        const QuiltParams Q = guo_quilt(cz);
        const PlainWeaveQuilt R = fiber::build_plain_weave_quilt(Q);
        PlainWeaveParams q = guo_params(cz);
        const PlainWeave W = fiber::build_plain_weave(q);
        check(R.W.params.nx == 102 && R.W.params.ny == 102, "nx = ny = 102");
        check(same_pts(R.W, W), "fibers[*].pts и vertex_bead совпадают побитово с build_plain_weave");
        std::ostringstream a, b;
        fiber::write_plain_weave_fib(a, R.W, "cmd");
        fiber::write_plain_weave_fib(b, W, "cmd");
        check(a.str() == b.str(), "текст .fib совпадает");
        check(R.W.clamp_boxes.size() == W.clamp_boxes.size(), "число CLAMP_BOX совпадает");
    }
}

// --- 2. Тиражирование патча -----------------------------------------------
static void replication_case(const QuiltParams& Q, const char* name)
{
    const PlainWeaveQuilt R = fiber::build_plain_weave_quilt(Q);
    const PlainWeave& W = R.W;
    const int k = Q.base.seg;
    const double p = Q.base.pitch, tol = 1e-12 * p;
    std::size_t compared = 0, bad = 0;
    double maxerr = 0.0;
    for (std::size_t f = 0; f < W.fibers.size(); ++f)
    {
        const auto& F = W.fibers[f];
        const bool warp = F.direction == 'x';
        const int P = warp ? Q.patch_nx : Q.patch_ny;       // нитей поперёк в патче
        const int Palong = warp ? Q.patch_ny : Q.patch_nx;  // вершин вдоль нити в патче
        const int row = F.index / P;                        // ряд патчей поперёк нити
        const int loc = F.index - row * P;                  // номер нити в патче
        const std::size_t ref_f = warp ? static_cast<std::size_t>(loc) : static_cast<std::size_t>(W.params.nx) + static_cast<std::size_t>(loc);
        const auto& Fr = W.fibers[ref_f];
        const std::size_t b_lo = F.vertex_bead.front(), b_hi = F.vertex_bead.back(); // окно без выступов
        for (std::size_t b = b_lo; b <= b_hi; ++b)
        {
            const int a = static_cast<int>(b - b_lo) / (Palong * k); // патч вдоль нити
            const int aa = std::min(a, (warp ? Q.tile_x : Q.tile_y) - 1);
            const std::size_t br = b - static_cast<std::size_t>(aa) * Palong * k;
            if (br >= Fr.pts.size()) continue;
            // Только бусины, которые в (0,0) тоже внутри окна и не выступ.
            if (br < Fr.vertex_bead.front() || br > Fr.vertex_bead.back()) continue;
            double sx = 0, sy = 0;
            if (warp)
            {
                sx = aa * Palong * p; // вдоль x
                sy = row * P * p;     // по y
            }
            else
            {
                sy = aa * Palong * p;
                sx = row * P * p;
            }
            const double ex = std::fabs(F.pts[b][0] - (Fr.pts[br][0] + sx));
            const double ey = std::fabs(F.pts[b][1] - (Fr.pts[br][1] + sy));
            const double ez = std::fabs(F.pts[b][2] - Fr.pts[br][2]);
            maxerr = std::max({maxerr, ex, ey, ez});
            ++compared;
            if (ex > tol || ey > tol || ez > tol) ++bad;
        }
    }
    std::printf("    %s: сравнено бусин %zu, макс. отклонение %.3e м\n", name, compared, maxerr);
    check(compared > 0 && bad == 0, std::string(name) + ": бусины патча (a, b) = бусины патча (0, 0) + сдвиг");
}

static void test_replication()
{
    std::printf("-- 2. тиражирование патча\n");
    replication_case(guo_quilt(), "Guo 34x34, 3x3, cos");
    replication_case(small_quilt(4, PlainWeaveProfile::Arcs, 4, 6, 3, 2), "4x6, 3x2, arcs, k=4");
    replication_case(small_quilt(3, PlainWeaveProfile::Cos, 2, 2, 4, 3), "2x2, 4x3, cos, k=3");
}

// --- 3. Нет разрывов -------------------------------------------------------
static void test_no_gaps()
{
    std::printf("-- 3. нет разрывов на стыках\n");
    const QuiltParams Q = guo_quilt();
    const PlainWeaveQuilt R = fiber::build_plain_weave_quilt(Q);
    PlainWeaveParams q = guo_params(Q.base.cz);
    const PlainWeave W = fiber::build_plain_weave(q);
    const auto S1 = fiber::plain_weave_summary(R.W), S2 = fiber::plain_weave_summary(W);
    check(S1.l0_min == S2.l0_min && S1.l0_max == S2.l0_max, "l0_min/l0_max ткани из патчей = одной ткани");
    check(S1.bonds == S2.bonds && S1.beads == S2.beads, "число бусин и связей совпадает");
    // Хорды внутри окна (без выступов): все равны; пересекающие стык -- тоже.
    double lref = 0.0;
    {
        const auto& F = R.W.fibers[0];
        lref = chord(F.pts[F.vertex_bead[0]], F.pts[F.vertex_bead[0] + 1]);
    }
    std::size_t crossing = 0, bad = 0;
    for (std::size_t f = 0; f < R.W.fibers.size(); ++f)
    {
        const auto& F = R.W.fibers[f];
        for (std::size_t b = F.vertex_bead.front(); b < F.vertex_bead.back(); ++b)
        {
            const double l = chord(F.pts[b], F.pts[b + 1]);
            const bool cross = R.bead_patch[f][b] != R.bead_patch[f][b + 1];
            if (cross) ++crossing;
            if (std::fabs(l - lref) > 1e-7 * lref) ++bad;
        }
    }
    std::printf("    хорд через стыки %zu, хорд с длиной != обычной %zu (l0 = %.9g м)\n", crossing, bad, lref);
    check(crossing > 0, "стыки пересекаются хордами");
    check(bad == 0, "все хорды окна равны по длине, в том числе через стыки");
}

// --- 4. Метки --------------------------------------------------------------
static void test_labels()
{
    std::printf("-- 4. метки бусин по патчам\n");
    const QuiltParams Q = guo_quilt();
    const PlainWeaveQuilt R = fiber::build_plain_weave_quilt(Q);
    const int nP = Q.tile_x * Q.tile_y;
    std::vector<std::size_t> cnt(static_cast<std::size_t>(nP), 0);
    std::size_t sum = 0;
    bool range = true;
    for (const auto& v : R.bead_patch)
        for (int id : v)
        {
            if (id < 0 || id >= nP) range = false;
            else ++cnt[static_cast<std::size_t>(id)], ++sum;
        }
    check(range, "id патча в [0, tile_x*tile_y)");
    check(sum == total_beads(R.W), "сумма бусин по патчам = общее число бусин");
    // Внутренний патч (1,1): 34 нити основы и 34 утка по P*k бусин на нить.
    const std::size_t inner = static_cast<std::size_t>(Q.patch_nx) * Q.patch_ny * Q.base.seg * 2;
    check(cnt[static_cast<std::size_t>(1 * Q.tile_x + 1)] == inner, "внутренний патч (1,1): " + std::to_string(inner) + " бусин");
    check(R.seam_x.size() == 2 && R.seam_y.size() == 2, "стыков по 2 на ось");
    // Стыковая бусина (k = 2: на стыке есть бусина) -> патч с меньшим индексом.
    const double p = Q.base.pitch, tol = 1e-9 * p;
    bool seam_ok = true;
    std::size_t seam_beads = 0;
    for (std::size_t f = 0; f < R.W.fibers.size(); ++f)
    {
        const auto& F = R.W.fibers[f];
        const auto& seams = F.direction == 'x' ? R.seam_x : R.seam_y;
        const int ax = F.direction == 'x' ? 0 : 1;
        for (std::size_t b = 0; b < F.pts.size(); ++b)
            for (std::size_t t = 0; t < seams.size(); ++t)
                if (std::fabs(F.pts[b][ax] - seams[t]) < tol)
                {
                    ++seam_beads;
                    // вдоль оси нити ряд = t (меньший), у соседней (следующей) бусины -- t+1
                    const int id = R.bead_patch[f][b];
                    const int row_along = F.direction == 'x' ? id % Q.tile_x : id / Q.tile_x;
                    const int row_next = F.direction == 'x' ? R.bead_patch[f][b + 1] % Q.tile_x : R.bead_patch[f][b + 1] / Q.tile_x;
                    if (row_along != static_cast<int>(t) || row_next != static_cast<int>(t) + 1) seam_ok = false;
                }
    }
    check(seam_beads == 2 * 102 * 2 && seam_ok, "бусина на стыке отнесена к патчу с меньшим индексом (" + std::to_string(seam_beads) + " бусин)");
    // Выступы: крайние патчи.
    bool ovh = true;
    for (std::size_t f = 0; f < R.W.fibers.size(); ++f)
    {
        const auto& F = R.W.fibers[f];
        const int first = R.bead_patch[f].front(), last = R.bead_patch[f].back();
        if (F.direction == 'x')
        {
            if (first % Q.tile_x != 0 || last % Q.tile_x != Q.tile_x - 1) ovh = false;
        }
        else if (first / Q.tile_x != 0 || last / Q.tile_x != Q.tile_y - 1)
            ovh = false;
    }
    check(ovh, "выступы относятся к крайним патчам");
}

// --- 5. Нечётный патч ------------------------------------------------------
static void test_odd_patch()
{
    std::printf("-- 5. нечётный патч\n");
    QuiltParams Q = small_quilt(4, PlainWeaveProfile::Cos, 3, 3, 2, 2);
    bool thrown = false;
    try
    {
        fiber::build_plain_weave_quilt(Q);
    }
    catch (const std::invalid_argument&)
    {
        thrown = true;
    }
    check(thrown, "без allow_odd_patch -> std::invalid_argument");
    Q.allow_odd_patch = true;
    PlainWeaveQuilt R;
    bool ok = true;
    try
    {
        R = fiber::build_plain_weave_quilt(Q);
    }
    catch (const std::exception&)
    {
        ok = false;
    }
    check(ok, "с allow_odd_patch построение удаётся");
    bool warned = false;
    for (const auto& w : R.W.warnings)
        if (w.find("нечётно") != std::string::npos) warned = true;
    check(warned, "в warnings есть предупреждение");
    // Нечётный патч при tile = 1 разрешён без флага.
    Q = small_quilt(4, PlainWeaveProfile::Cos, 3, 5, 1, 1);
    ok = true;
    try
    {
        R = fiber::build_plain_weave_quilt(Q);
    }
    catch (const std::exception&)
    {
        ok = false;
    }
    check(ok && R.W.warnings.empty(), "нечётный патч при tile = 1: без ошибки и предупреждения");
    // Параметры < 1.
    Q = small_quilt(4, PlainWeaveProfile::Cos, 4, 4, 2, 2);
    Q.tile_y = 0;
    thrown = false;
    try
    {
        fiber::build_plain_weave_quilt(Q);
    }
    catch (const std::invalid_argument&)
    {
        thrown = true;
    }
    check(thrown, "tile_y = 0 -> std::invalid_argument");
}

// --- 6. Размер -------------------------------------------------------------
static void test_window()
{
    std::printf("-- 6. окно ткани\n");
    const PlainWeaveQuilt R = fiber::build_plain_weave_quilt(guo_quilt());
    double Lx = 0, Ly = 0;
    fiber::plain_weave_window(R.W.params, Lx, Ly);
    std::printf("    окно %.9g x %.9g мм\n", Lx * 1e3, Ly * 1e3);
    check(std::fabs(Lx - 101 * 1.18e-3) < 1e-12 && std::fabs(Ly - 101 * 1.18e-3) < 1e-12, "окно = 101 * 1.18 мм");
    check(std::fabs(Lx * 1e3 - 119.2) < 0.05 && std::fabs(Ly * 1e3 - 119.2) < 0.05, "окно ~ 119.2 x 119.2 мм");
}

// --- 7. Пересечения и проникание ------------------------------------------
static double seg_dist(const std::array<double, 3>& p1, const std::array<double, 3>& q1, const std::array<double, 3>& p2,
                       const std::array<double, 3>& q2)
{
    auto sub = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };
    auto dot = [](const std::array<double, 3>& a, const std::array<double, 3>& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
    auto clamp01 = [](double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); };
    const auto d1 = sub(q1, p1), d2 = sub(q2, p2), r = sub(p1, p2);
    const double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    double s, t;
    const double cc = dot(d1, r), b = dot(d1, d2);
    const double denom = a * e - b * b;
    s = denom > 1e-14 * a * e ? clamp01((b * f - cc * e) / denom) : 0.0;
    t = (b * s + f) / e;
    if (t < 0.0)
    {
        t = 0.0;
        s = clamp01(-cc / a);
    }
    else if (t > 1.0)
    {
        t = 1.0;
        s = clamp01((b - cc) / a);
    }
    std::array<double, 3> c1, c2;
    for (int i = 0; i < 3; ++i)
    {
        c1[i] = p1[i] + d1[i] * s;
        c2[i] = p2[i] + d2[i] * t;
    }
    return chord(c1, c2);
}

static void weave_case(const QuiltParams& Q, const char* name, bool pairwise)
{
    const PlainWeaveQuilt R = fiber::build_plain_weave_quilt(Q);
    const PlainWeave& W = R.W;
    const int nx = W.params.nx, ny = W.params.ny;
    const double A = Q.base.amp, d = Q.base.d;
    double worst = 0.0;
    for (int j = 0; j < nx; ++j)
        for (int i = 0; i < ny; ++i)
        {
            const auto& Fw = W.fibers[static_cast<std::size_t>(j)];
            const auto& Ff = W.fibers[static_cast<std::size_t>(nx + i)];
            const double zw = Fw.pts[Fw.vertex_bead[static_cast<std::size_t>(i)]][2];
            const double zf = Ff.pts[Ff.vertex_bead[static_cast<std::size_t>(j)]][2];
            worst = std::max(worst, std::fabs(std::fabs(zw - zf) - 2.0 * A));
        }
    std::printf("    %s: макс. ||zw - zf| - 2A| = %.3e м\n", name, worst);
    check(worst < 1e-12 * Q.base.pitch, std::string(name) + ": в каждом пересечении |z_основы - z_утка| = 2A");
    if (!pairwise) return;
    double dmin = 1e300;
    for (std::size_t a = 0; a < W.fibers.size(); ++a)
        for (std::size_t b = a + 1; b < W.fibers.size(); ++b)
        {
            const auto &Fa = W.fibers[a], &Fb = W.fibers[b];
            for (std::size_t i = 0; i + 1 < Fa.pts.size(); ++i)
                for (std::size_t j = 0; j + 1 < Fb.pts.size(); ++j)
                    dmin = std::min(dmin, seg_dist(Fa.pts[i], Fa.pts[i + 1], Fb.pts[j], Fb.pts[j + 1]));
        }
    std::printf("    %s: min расстояние сегмент-сегмент %.9g * d\n", name, dmin / d);
    check(dmin > 0.0, std::string(name) + ": оси не пересекаются");
    if (2.0 * A > d && Q.base.seg >= 4)
        check(dmin >= d * (1.0 - 1e-9), std::string(name) + ": при 2A > d нет проникания (min dist >= d)");
    else if (dmin < d * (1.0 - 1e-9))
        std::printf("    НАХОДКА (как в test_plain_weave): при k = 2 хорды сближаются у узла сильнее оси, min dist = %.4g*d; это свойство дискретизации, не квилта\n", dmin / d);
}

static void test_weave()
{
    std::printf("-- 7. пересечения и проникание\n");
    weave_case(guo_quilt(), "Guo 34x34, 3x3", false);
    weave_case(small_quilt(4, PlainWeaveProfile::Cos, 4, 4, 2, 2), "4x4, 2x2, cos, k=4", true);
    weave_case(small_quilt(4, PlainWeaveProfile::Arcs, 2, 4, 3, 2), "2x4, 3x2, arcs, k=4", true);
    weave_case(small_quilt(2, PlainWeaveProfile::Cos, 4, 4, 2, 2), "4x4, 2x2, cos, k=2", true);
}

// --- 8. .patches.txt -------------------------------------------------------
static void test_patches_file()
{
    std::printf("-- 8. формат .patches.txt\n");
    const QuiltParams Q = small_quilt(2, PlainWeaveProfile::Cos, 4, 4, 2, 2);
    const PlainWeaveQuilt R = fiber::build_plain_weave_quilt(Q);
    std::ostringstream os;
    fiber::write_plain_weave_quilt_patches(os, R);
    std::istringstream is(os.str());
    std::string line;
    std::size_t rows = 0, beads = 0;
    bool ok = true;
    while (std::getline(is, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        int oid, idx;
        char dir;
        std::size_t first, last, nr;
        ls >> oid >> dir >> idx >> first >> last >> nr;
        std::size_t covered = 0;
        for (std::size_t r = 0; r < nr; ++r)
        {
            std::string tok;
            ls >> tok;
            unsigned pid;
            std::size_t f1, l1;
            if (std::sscanf(tok.c_str(), "%u:%zu:%zu", &pid, &f1, &l1) != 3) ok = false;
            covered += l1 - f1 + 1;
        }
        if (covered != last - first + 1) ok = false;
        beads += covered;
        ++rows;
    }
    check(ok && rows == R.W.fibers.size() && beads == total_beads(R.W), "строка на нить, диапазоны покрывают все бусины");
}

// --- 9. Стёжка -------------------------------------------------------------
// Пакет как build_fabric_stack: копия q, cz сдвинут вниз на k*(2A + d + 1e-6).
static std::vector<PlainWeave> make_plies(const QuiltParams& Q, int n, PlainWeaveQuilt& first)
{
    std::vector<PlainWeave> plies;
    const double pitch = 2.0 * Q.base.amp + Q.base.d + 1e-6;
    for (int k = 0; k < n; ++k)
    {
        QuiltParams q = Q;
        q.base.cz = Q.base.cz - k * pitch;
        PlainWeaveQuilt R = fiber::build_plain_weave_quilt(q);
        if (k == 0) first = R;
        plies.push_back(std::move(R.W));
    }
    return plies;
}

static void stitch_case(const QuiltParams& Q, int nplies, int every, bool lx, bool ly, const char* name)
{
    PlainWeaveQuilt R0;
    const std::vector<PlainWeave> plies = make_plies(Q, nplies, R0);
    fiber::QuiltStitchParams sp;
    sp.every = every;
    sp.lines_x = lx;
    sp.lines_y = ly;
    const auto st = fiber::build_quilt_stitches(plies, R0.seam_x, R0.seam_y, sp);
    const int nx = R0.W.params.nx, ny = R0.W.params.ny;
    const std::size_t per_x = static_cast<std::size_t>((nx + every - 1) / every); // проколов на линии x (основа)
    const std::size_t per_y = static_cast<std::size_t>((ny + every - 1) / every);
    const std::size_t lines = (lx ? R0.seam_x.size() * per_x : 0) + (ly ? R0.seam_y.size() * per_y : 0);
    const std::size_t expect = lines * static_cast<std::size_t>(nplies - 1);
    std::printf("    %s: стежков %zu, ожидалось %zu\n", name, st.size(), expect);
    check(st.size() == expect, std::string(name) + ": число = (линии x проколы) x (N-1)");

    // Плоский массив бусин в порядке записи .fib слой за слоем.
    std::vector<std::array<double, 3>> all;
    std::vector<std::size_t> ply_of;
    for (std::size_t k = 0; k < plies.size(); ++k)
        for (const auto& F : plies[k].fibers)
            for (const auto& P : F.pts)
            {
                all.push_back(P);
                ply_of.push_back(k);
            }
    const double p = Q.base.pitch;
    bool vert = true, adj = true, onseam = true;
    std::set<std::pair<std::size_t, std::size_t>> uniq(st.begin(), st.end());
    for (const auto& pr : st)
    {
        if (pr.first >= all.size() || pr.second >= all.size())
        {
            vert = false;
            continue;
        }
        const auto &a = all[pr.first], &b = all[pr.second];
        if (std::fabs(a[0] - b[0]) >= 1e-12 || std::fabs(a[1] - b[1]) >= 1e-12) vert = false;
        if (ply_of[pr.second] != ply_of[pr.first] + 1) adj = false;
        const bool on_x = std::any_of(R0.seam_x.begin(), R0.seam_x.end(), [&](double s) { return std::fabs(a[0] - s) < 1e-9 * p; });
        const bool on_y = std::any_of(R0.seam_y.begin(), R0.seam_y.end(), [&](double s) { return std::fabs(a[1] - s) < 1e-9 * p; });
        if (!(on_x || on_y)) onseam = false;
        if (std::fabs(a[2] - b[2]) < 1e-12) vert = false; // слои разнесены по z
    }
    check(vert, std::string(name) + ": пары лежат на одной вертикали (|dx|, |dy| < 1e-12), слои разнесены по z");
    check(adj, std::string(name) + ": пары принадлежат соседним слоям (k, k+1)");
    check(onseam, std::string(name) + ": бусины на линиях шва");
    check(uniq.size() == st.size(), std::string(name) + ": пары не повторяются");
}

static void test_stitches()
{
    std::printf("-- 9. сетка стёжки\n");
    const QuiltParams Q = small_quilt(2, PlainWeaveProfile::Cos, 4, 6, 3, 2); // nx = 8, ny = 18; швов x: 2, y: 1
    stitch_case(Q, 2, 1, true, true, "N=2, every=1");
    stitch_case(Q, 4, 1, true, true, "N=4, every=1");
    stitch_case(Q, 4, 2, true, true, "N=4, every=2");
    stitch_case(Q, 2, 3, true, true, "N=2, every=3");
    stitch_case(Q, 3, 1, true, false, "N=3, только линии x");
    stitch_case(Q, 3, 1, false, true, "N=3, только линии y");
    const QuiltParams Q4 = small_quilt(4, PlainWeaveProfile::Arcs, 4, 4, 2, 2);
    stitch_case(Q4, 2, 1, true, true, "arcs k=4, N=2");

    // Точное число для every=1 и every=2 (прямая формула).
    {
        PlainWeaveQuilt R0;
        const auto plies = make_plies(Q, 4, R0);
        fiber::QuiltStitchParams a, b;
        b.every = 2;
        const auto s1 = fiber::build_quilt_stitches(plies, R0.seam_x, R0.seam_y, a);
        const auto s2 = fiber::build_quilt_stitches(plies, R0.seam_x, R0.seam_y, b);
        const std::size_t e1 = (2 * 8 + 1 * 18) * 3;
        const std::size_t e2 = (2 * 4 + 1 * 9) * 3;
        check(s1.size() == e1, "every=1: (2*8 + 1*18)*3 = " + std::to_string(e1));
        check(s2.size() == e2, "every=2: (2*4 + 1*9)*3 = " + std::to_string(e2));
    }
    // Один слой: стежков нет.
    {
        PlainWeaveQuilt R0;
        const auto plies = make_plies(Q, 1, R0);
        check(fiber::build_quilt_stitches(plies, R0.seam_x, R0.seam_y, {}).empty(), "один слой: стежков нет");
    }
    // Нечётный seg: на стыке нет бусины -> исключение; seams через quilt_seams.
    {
        const QuiltParams Qo = small_quilt(3, PlainWeaveProfile::Cos, 4, 4, 2, 2);
        PlainWeaveQuilt R0;
        const auto plies = make_plies(Qo, 2, R0);
        bool thrown = false;
        try
        {
            fiber::build_quilt_stitches(plies, R0.seam_x, R0.seam_y, {});
        }
        catch (const std::invalid_argument&)
        {
            thrown = true;
        }
        check(thrown, "нечётный seg: std::invalid_argument");
    }
    // quilt_seams согласован с build_plain_weave_quilt; нечётный патч без ограничений.
    {
        PlainWeaveQuilt R0;
        make_plies(Q, 1, R0);
        std::vector<double> sx, sy;
        fiber::quilt_seams(R0.W.params, Q.patch_nx, Q.patch_ny, Q.tile_x, Q.tile_y, sx, sy);
        check(sx == R0.seam_x && sy == R0.seam_y, "quilt_seams = seam_x/seam_y квилта");
        PlainWeaveParams full = Q.base;
        full.nx = 9;
        full.ny = 15; // патчи 3x5, плитки 3x3: нечётные
        fiber::quilt_seams(full, 3, 5, 3, 3, sx, sy);
        check(sx.size() == 2 && sy.size() == 2, "quilt_seams для нечётного патча: по 2 шва на ось");
        // шов лежит на середине полуволны между вершинами 4 и 5 (вдоль x)
        const double x4 = full.cx + (4 - 0.5 * (full.ny - 1)) * full.pitch;
        check(std::fabs(sx[0] - (x4 + 0.5 * full.pitch)) < 1e-12 * full.pitch, "шов между вершинами t*P-1 и t*P");
    }
    // Слои не на одной вертикали / разного состава -> исключение.
    {
        PlainWeaveQuilt R0;
        auto plies = make_plies(Q, 2, R0);
        QuiltParams Qs = Q;
        Qs.base.cx += 1e-5;
        plies[1] = fiber::build_plain_weave_quilt(Qs).W;
        bool thrown = false;
        try
        {
            fiber::build_quilt_stitches(plies, R0.seam_x, R0.seam_y, {});
        }
        catch (const std::invalid_argument&)
        {
            thrown = true;
        }
        check(thrown, "слои со сдвигом в плоскости: std::invalid_argument");
        Qs = small_quilt(2, PlainWeaveProfile::Cos, 4, 4, 2, 2);
        plies[1] = fiber::build_plain_weave_quilt(Qs).W;
        thrown = false;
        try
        {
            fiber::build_quilt_stitches(plies, R0.seam_x, R0.seam_y, {});
        }
        catch (const std::invalid_argument&)
        {
            thrown = true;
        }
        check(thrown, "слои разного состава: std::invalid_argument");
        fiber::QuiltStitchParams bad;
        bad.every = 0;
        thrown = false;
        try
        {
            fiber::build_quilt_stitches(plies, R0.seam_x, R0.seam_y, bad);
        }
        catch (const std::invalid_argument&)
        {
            thrown = true;
        }
        check(thrown, "every = 0: std::invalid_argument");
    }
}

int main()
{
    test_equivalence();
    test_replication();
    test_no_gaps();
    test_labels();
    test_odd_patch();
    test_window();
    test_weave();
    test_patches_file();
    test_stitches();
    if (failures == 0)
        std::printf("test_plain_weave_quilt: все проверки пройдены\n");
    else
        std::printf("test_plain_weave_quilt: ПРОВАЛ, %d проверок\n", failures);
    return failures == 0 ? 0 : 1;
}
