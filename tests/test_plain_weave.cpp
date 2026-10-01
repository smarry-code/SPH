// tests/test_plain_weave.cpp
// Тесты генератора полотняного переплетения (include/fiber/plain_weave.hpp,
// ТЗ weave/TZ_PLAIN_WEAVE_GENERATOR.md, §5). Этап 1: проверки 1-4, 8, 9
// (+ согласованность подсчёта зажатых бусин).
//
// Разделы:
//   1 test_count            -- число нитей/бусин/связей по формуле ТЗ;
//   2 test_weave            -- чередование над/под в пересечениях;
//   3 test_no_penetration   -- минимум расстояния сегмент-сегмент между нитями;
//   4 test_equal_chords     -- равные хорды (l0) внутри полуволны и выступа;
//   8 test_arcs_on_ellipses -- профиль arcs: точки на полуэллипсах;
//   9 test_param_errors     -- ошибки параметров;
//   + test_clamp_counts     -- число зажатых бусин в нечётном варианте.
//   5 test_roundtrip        -- запись .fib, load_fiber_file, сравнение с построенным;
//   6 test_density          -- поверхностная плотность (сводка = формула; Twaron);
//   7 test_contact          -- диагностика контакта на исходной геометрии;
//   7б test_contact_twaron  -- то же на геометрии Twaron, 12x12, k = 2/4, A = d/2, 0.6d.
#include "fiber/fiber_input.hpp"
#include "fiber/model.hpp"
#include "fiber/plain_weave.hpp"
#include "interactions/clamp.hpp"
#include "interactions/fiber_contact.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using fiber::PlainWeave;
using fiber::PlainWeaveFiber;
using fiber::PlainWeaveParams;
using fiber::PlainWeaveProfile;

static int failures = 0;

static void check(bool cond, const std::string& what)
{
    std::printf("  [%s] %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

struct Case
{
    std::string name;
    PlainWeaveParams q;
    PlainWeave W;
};

static PlainWeaveParams make_test_params()
{
    PlainWeaveParams q;
    fiber::apply_plain_weave_preset("test", q);
    return q;
}

static PlainWeaveParams make_odd_params(PlainWeaveProfile profile)
{
    PlainWeaveParams q = make_test_params();
    q.nx = 3;
    q.ny = 5;
    q.amp = 0.75 * q.d;
    q.seg = 3;
    q.overhang_x = 0.5 * q.pitch;
    q.overhang_y = 1.5 * q.pitch;
    q.profile = profile;
    q.clamp = fiber::PlainWeaveClamp::Y;
    // Середина выступа утка: y_{nx-1}-cy = (nx-1)/2*p, плюс 0.75 p.
    q.clamp_edge = 0.5 * (q.nx - 1) * q.pitch + 0.75 * q.pitch;
    q.clamp_edge_set = true;
    return q;
}

static std::string fiber_name(const PlainWeaveFiber& F)
{
    return "нить " + std::to_string(F.object_id) + " (" + F.direction + std::to_string(F.index) + ")";
}

static int overhang_segments(double e, double p, int k)
{
    return std::max(1, static_cast<int>(std::lround(e / (p / k))));
}

static double seg_len(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return std::sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2]));
}

// --- 1. Счёт ---------------------------------------------------------------
static void test_count(const Case& c)
{
    std::printf("-- 1. счёт (%s)\n", c.name.c_str());
    const auto& q = c.q;
    const auto& W = c.W;
    check(W.fibers.size() == static_cast<std::size_t>(q.nx + q.ny), "нитей nx+ny");
    bool ok = true;
    std::size_t beads = 0;
    for (const auto& F : W.fibers)
    {
        const bool warp = F.direction == 'x';
        const int nv = warp ? q.ny : q.nx; // пересечений вдоль нити
        const double e = warp ? fiber::plain_weave_overhang_x(q) : fiber::plain_weave_overhang_y(q);
        const std::size_t expect = static_cast<std::size_t>((nv - 1) * q.seg + 2 * overhang_segments(e, q.pitch, q.seg) + 1);
        beads += F.pts.size();
        if (F.pts.size() != expect)
        {
            ok = false;
            std::printf("    %s: бусин %zu, ожидалось %zu\n", fiber_name(F).c_str(), F.pts.size(), expect);
        }
        if (F.vertex_bead.size() != static_cast<std::size_t>(nv))
        {
            ok = false;
            std::printf("    %s: вершин %zu, ожидалось %d\n", fiber_name(F).c_str(), F.vertex_bead.size(), nv);
        }
    }
    check(ok, "бусин в каждой нити = N_half*k + 2*n_e + 1, вершин = число пересечений");
    const auto S = fiber::plain_weave_summary(W);
    check(S.beads == beads && S.bonds == beads - W.fibers.size(), "связей на одну меньше, чем бусин (на нить)");
    std::printf("    нитей %zu, бусин %zu (на нить: основа %zu, уток %zu), связей %zu\n", S.fibers, S.beads,
                W.fibers.front().pts.size(), W.fibers.back().pts.size(), S.bonds);
}

// --- 2. Переплетение -------------------------------------------------------
static void test_weave(const Case& c)
{
    std::printf("-- 2. переплетение (%s)\n", c.name.c_str());
    const auto& q = c.q;
    const auto& W = c.W;
    const double tol = 1e-12 * q.pitch;
    bool ok = true, ok_xy = true, ok_alt = true;
    double max_dz_err = 0.0;
    for (int j = 0; j < q.nx; ++j)
        for (int i = 0; i < q.ny; ++i)
        {
            const auto& Fw = W.fibers[j];
            const auto& Ff = W.fibers[q.nx + i];
            const auto& Pw = Fw.pts[Fw.vertex_bead[i]];
            const auto& Pf = Ff.pts[Ff.vertex_bead[j]];
            const double xi = q.cx + (i - 0.5 * (q.ny - 1)) * q.pitch;
            const double yj = q.cy + (j - 0.5 * (q.nx - 1)) * q.pitch;
            const double sgn = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            const double dz = Pw[2] - Pf[2];
            const double err = std::fabs(dz - 2.0 * q.amp * sgn);
            max_dz_err = std::max(max_dz_err, err);
            if (err > tol)
            {
                ok = false;
                std::printf("    (i=%d,j=%d): dz = %.17g, ожидалось %.17g\n", i, j, dz, 2.0 * q.amp * sgn);
            }
            if (Pw[0] != xi || Pw[1] != yj || Pf[0] != xi || Pf[1] != yj)
            {
                ok_xy = false;
                std::printf("    (i=%d,j=%d): xy основы (%.17g,%.17g), утка (%.17g,%.17g), ожидалось (%.17g,%.17g)\n", i, j,
                            Pw[0], Pw[1], Pf[0], Pf[1], xi, yj);
            }
        }
    check(ok, "dz = z_основы - z_утка = 2A*(-1)^(i+j) с допуском 1e-12*p");
    check(ok_xy, "xy бусин обеих нитей в пересечении совпадают с (x_i, y_j) точно");
    for (const auto& F : W.fibers)
        for (std::size_t v = 0; v + 1 < F.vertex_bead.size(); ++v)
        {
            const double za = F.pts[F.vertex_bead[v]][2] - q.cz;
            const double zb = F.pts[F.vertex_bead[v + 1]][2] - q.cz;
            if (!(za * zb < 0.0) || std::fabs(std::fabs(za) - q.amp) > 1e-12 * q.pitch)
            {
                ok_alt = false;
                std::printf("    %s: вершины %zu-%zu: z-cz = %.17g, %.17g\n", fiber_name(F).c_str(), v, v + 1, za, zb);
            }
        }
    check(ok_alt, "вдоль нити знак z-cz в вершинах чередуется, модуль = A");
    std::printf("    макс. |dz - 2A(-1)^(i+j)| = %.3e\n", max_dz_err);
}

// --- 3. Нет проникания -----------------------------------------------------
// Расстояние сегмент-сегмент (Ericson, Real-Time Collision Detection, 5.1.9),
// с отсечением параметров в [0, 1]; возвращает точки и параметры.
struct SegDist
{
    double dist, s, t;
    std::array<double, 3> c1, c2;
};

static SegDist segment_segment(const std::array<double, 3>& p1, const std::array<double, 3>& q1,
                               const std::array<double, 3>& p2, const std::array<double, 3>& q2)
{
    auto sub = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };
    auto dot = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    auto clamp01 = [](double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); };
    const auto d1 = sub(q1, p1), d2 = sub(q2, p2), r = sub(p1, p2);
    const double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    double s, t;
    if (a <= 1e-300 && e <= 1e-300)
        s = t = 0.0;
    else if (a <= 1e-300)
    {
        s = 0.0;
        t = clamp01(f / e);
    }
    else
    {
        const double cc = dot(d1, r);
        if (e <= 1e-300)
        {
            t = 0.0;
            s = clamp01(-cc / a);
        }
        else
        {
            const double b = dot(d1, d2);
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
        }
    }
    SegDist R;
    for (int k = 0; k < 3; ++k)
    {
        R.c1[k] = p1[k] + d1[k] * s;
        R.c2[k] = p2[k] + d2[k] * t;
    }
    const auto dv = sub(R.c1, R.c2);
    R.dist = std::sqrt(dot(dv, dv));
    R.s = s;
    R.t = t;
    return R;
}

static void test_no_penetration(const Case& c)
{
    std::printf("-- 3. нет проникания (%s)\n", c.name.c_str());
    const auto& q = c.q;
    const auto& W = c.W;
    double dmin = 1e300;
    SegDist best{};
    int fa = -1, fb = -1;
    std::size_t sa = 0, sb = 0;
    // Минимум вне окрестности узла: расстояние от точки на кривой до вершин-узлов
    // обеих нитей > p/(2k)... печатаем отдельно, где лежит минимум относительно узла.
    for (std::size_t a = 0; a < W.fibers.size(); ++a)
        for (std::size_t b = a + 1; b < W.fibers.size(); ++b)
        {
            const auto &Fa = W.fibers[a], &Fb = W.fibers[b];
            for (std::size_t i = 0; i + 1 < Fa.pts.size(); ++i)
                for (std::size_t j = 0; j + 1 < Fb.pts.size(); ++j)
                {
                    const SegDist r = segment_segment(Fa.pts[i], Fa.pts[i + 1], Fb.pts[j], Fb.pts[j + 1]);
                    if (r.dist < dmin)
                    {
                        dmin = r.dist;
                        best = r;
                        fa = static_cast<int>(a);
                        fb = static_cast<int>(b);
                        sa = i;
                        sb = j;
                    }
                }
        }
    std::printf("    минимум расстояния сегмент-сегмент: %.9g м = %.9g * d = %.9g * 2A\n", dmin, dmin / q.d,
                dmin / (2.0 * q.amp));
    std::printf("    достигается: нити %d и %d, сегменты %zu и %zu, параметры s=%.6g t=%.6g\n", fa, fb, sa, sb, best.s,
                best.t);
    std::printf("    точки: (%.9g, %.9g, %.9g) и (%.9g, %.9g, %.9g)\n", best.c1[0], best.c1[1], best.c1[2], best.c2[0],
                best.c2[1], best.c2[2]);
    // ТЗ §5 п.3: если минимум < 2A, допуск не ослабляется, это находка. При A = d/2
    // (2A = d) минимум < d возникает у самого узла: оси нитей в модели -- хорды
    // соседних бусин, а не гладкие кривые, и хорды сближаются меньше 2A. Это не
    // ошибка генератора, поэтому жёсткий гейт здесь -- отсутствие пересечения
    // (dmin > 0), строгое dmin >= d(1-1e-9) печатается как НАХОДКА (не FAIL).
    // При зазоре (2A > d) хорды обязаны держать d: здесь гейт строгий.
    check(dmin > 0.0, "min dist сегмент-сегмент > 0 (оси не пересекаются)");
    if (2.0 * q.amp > q.d) check(dmin >= q.d * (1.0 - 1e-9), "при 2A > d: min dist >= d*(1-1e-9)");
    if (dmin >= q.d * (1.0 - 1e-9))
        std::printf("    min dist >= d*(1-1e-9): проникания нет\n");
    else
        std::printf("    НАХОДКА: min dist = %.6g*d < d (хорды сближаются в узле сильнее оси, проникание %.3g%% от d); допуск не ослаблен\n",
                    dmin / q.d, 100.0 * (1.0 - dmin / q.d));
}

// --- 4. Равные хорды -------------------------------------------------------
static void test_equal_chords(const Case& c)
{
    std::printf("-- 4. равные хорды (%s)\n", c.name.c_str());
    const auto& q = c.q;
    const auto& W = c.W;
    bool ok = true;
    double worst = 0.0;
    auto group = [&](const PlainWeaveFiber& F, std::size_t b0, std::size_t b1, const char* what, std::size_t idx) {
        double lmin = 1e300, lmax = 0.0;
        for (std::size_t b = b0; b < b1; ++b)
        {
            const double l = seg_len(F.pts[b], F.pts[b + 1]);
            lmin = std::min(lmin, l);
            lmax = std::max(lmax, l);
        }
        const double rel = (lmax - lmin) / lmax;
        worst = std::max(worst, rel);
        if (rel > 1e-9)
        {
            ok = false;
            std::printf("    %s, %s %zu: l0 от %.17g до %.17g (отн. разброс %.3e), бусины %zu..%zu\n",
                        fiber_name(F).c_str(), what, idx, lmin, lmax, rel, b0, b1);
        }
    };
    for (const auto& F : W.fibers)
    {
        const auto& vb = F.vertex_bead;
        if (vb.front() > 0) group(F, 0, vb.front(), "выступ (начало)", 0);
        for (std::size_t w = 0; w + 1 < vb.size(); ++w) group(F, vb[w], vb[w + 1], "полуволна", w);
        if (vb.back() + 1 < F.pts.size()) group(F, vb.back(), F.pts.size() - 1, "выступ (конец)", 0);
    }
    check(ok, "l0 равны (отн. 1e-9) в каждой полуволне и в каждом выступе");
    const auto S = fiber::plain_weave_summary(W);
    std::printf("    min l0 = %.9g м, max l0 = %.9g м (p/k = %.9g), макс. отн. разброс в группе = %.3e\n", S.l0_min,
                S.l0_max, q.pitch / q.seg, worst);
}

// --- 8. Профиль arcs -------------------------------------------------------
// Нормировка: ((s - c)/(p/2))^2 + ((z - cz)/A)^2 = 1, c = x_0 + m*p -- центр
// ближайшей дуги (вдоль нити), допуск на |левая часть - 1| равен 1e-12.
static void test_arcs_on_ellipses(const Case& c)
{
    std::printf("-- 8. arcs: точки на полуэллипсах (%s)\n", c.name.c_str());
    const auto& q = c.q;
    const auto& W = c.W;
    bool ok = true, ok_vertex = true;
    double worst = 0.0;
    for (const auto& F : W.fibers)
    {
        const int axis = F.direction == 'x' ? 0 : 1;
        const double v0 = F.pts[F.vertex_bead[0]][axis];
        for (std::size_t b = 0; b < F.pts.size(); ++b)
        {
            const double s = F.pts[b][axis], z = F.pts[b][2];
            const double m = std::round((s - v0) / q.pitch);
            const double cc = v0 + m * q.pitch;
            const double v = std::pow((s - cc) / (0.5 * q.pitch), 2) + std::pow((z - q.cz) / q.amp, 2);
            worst = std::max(worst, std::fabs(v - 1.0));
            if (std::fabs(v - 1.0) > 1e-12)
            {
                ok = false;
                std::printf("    %s, бусина %zu: (%.17g, %.17g), невязка %.3e\n", fiber_name(F).c_str(), b, s, z, v - 1.0);
            }
        }
        for (std::size_t w = 0; w < F.vertex_bead.size(); ++w)
        {
            const auto& P = F.pts[F.vertex_bead[w]];
            if (std::fabs(std::fabs(P[2] - q.cz) - q.amp) > 1e-15 * q.pitch + 0.0)
            {
                ok_vertex = false;
                std::printf("    %s: вершина %zu не на вершине дуги: z-cz = %.17g\n", fiber_name(F).c_str(), w,
                            P[2] - q.cz);
            }
        }
    }
    check(ok, "каждая бусина на полуэллипсе (допуск 1e-12)");
    check(ok_vertex, "бусины стоят на вершинах дуг (|z-cz| = A) в пересечениях");
    std::printf("    макс. невязка уравнения эллипса = %.3e\n", worst);
}

// --- 9. Ошибки параметров --------------------------------------------------
template <class F> static bool throws_with_text(F&& f, std::string& msg)
{
    try
    {
        f();
    }
    catch (const std::exception& e)
    {
        msg = e.what();
        return !msg.empty();
    }
    return false;
}

static void test_param_errors()
{
    std::printf("-- 9. ошибки параметров\n");
    std::string msg;
    {
        PlainWeaveParams q = make_test_params();
        q.pitch = 0.9 * q.d;
        check(throws_with_text([&] { fiber::build_plain_weave(q); }, msg), "p < d -> исключение");
        std::printf("    текст: %s\n", msg.c_str());
        q.allow_penetration = true;
        bool no_throw = true;
        PlainWeave W;
        try
        {
            W = fiber::build_plain_weave(q);
        }
        catch (...)
        {
            no_throw = false;
        }
        check(no_throw && !W.warnings.empty() && W.warnings[0].rfind("WARN", 0) == 0,
              "p < d с allow_penetration -> успех и WARN в warnings");
    }
    {
        PlainWeaveParams q = make_test_params();
        q.amp = 0.4 * q.d; // 2A < d
        check(throws_with_text([&] { fiber::build_plain_weave(q); }, msg), "2A < d -> исключение");
        std::printf("    текст: %s\n", msg.c_str());
        q.allow_penetration = true;
        bool no_throw = true;
        PlainWeave W;
        try
        {
            W = fiber::build_plain_weave(q);
        }
        catch (...)
        {
            no_throw = false;
        }
        check(no_throw && !W.warnings.empty() && W.warnings[0].rfind("WARN", 0) == 0,
              "2A < d с allow_penetration -> успех и WARN в warnings");
    }
    {
        PlainWeaveParams q = make_test_params();
        q.clamp = fiber::PlainWeaveClamp::Y; // c_e не задан
        check(throws_with_text([&] { fiber::build_plain_weave(q); }, msg), "clamp=Y без clamp_edge -> исключение");
        std::printf("    текст: %s\n", msg.c_str());
        q.allow_penetration = true;
        check(throws_with_text([&] { fiber::build_plain_weave(q); }, msg), "то же с allow_penetration -> исключение");
    }
    {
        // Тест: nx = ny = 4, p = 1 мм, d = 0.5 мм, e = p/2: y_{nx-1}-cy = 1.5p,
        // диапазон c_e: [1.5p + d/2, 1.5p + p/2] = [1.75p, 2p] (в p: d/2 = 0.25p).
        PlainWeaveParams q = make_test_params();
        q.clamp = fiber::PlainWeaveClamp::Y;
        q.clamp_edge_set = true;
        const double lo = 1.5 * q.pitch + 0.5 * q.d, hi = 1.5 * q.pitch + plain_weave_overhang_y(q);
        const double vals[3] = {lo - 0.01 * q.pitch, hi + 0.01 * q.pitch, 0.0};
        const char* names[3] = {"c_e ниже диапазона (в зажим попадают бусины основы)", "c_e выше диапазона (нет зажатой бусины)",
                                "c_e = 0"};
        for (int a = 0; a < 3; ++a)
        {
            q.clamp_edge = vals[a];
            q.allow_penetration = false;
            check(throws_with_text([&] { fiber::build_plain_weave(q); }, msg), std::string(names[a]) + " -> исключение");
            q.allow_penetration = true;
            check(throws_with_text([&] { fiber::build_plain_weave(q); }, msg),
                  std::string(names[a]) + " с allow_penetration -> исключение");
        }
        std::printf("    текст: %s\n", msg.c_str());
        q.allow_penetration = false;
        q.clamp_edge = lo;
        bool ok_lo = true;
        try
        {
            fiber::build_plain_weave(q);
        }
        catch (...)
        {
            ok_lo = false;
        }
        check(ok_lo, "c_e на нижней границе диапазона допустим");
        q.clamp_edge = 0.5 * (lo + hi);
        bool ok_mid = true;
        try
        {
            fiber::build_plain_weave(q);
        }
        catch (...)
        {
            ok_mid = false;
        }
        check(ok_mid, "c_e посередине диапазона допустим");
    }
    {
        // Зажим по x: диапазон по x_{ny-1}-cx и e_x.
        PlainWeaveParams q = make_test_params();
        q.clamp = fiber::PlainWeaveClamp::X;
        q.clamp_edge_set = true;
        q.clamp_edge = 1.5 * q.pitch + plain_weave_overhang_x(q) + 0.01 * q.pitch;
        check(throws_with_text([&] { fiber::build_plain_weave(q); }, msg), "clamp=X, c_e выше диапазона -> исключение");
    }
}

// --- Зажатые бусины (нечётный вариант) --------------------------------------
static void test_clamp_counts(const Case& c)
{
    std::printf("-- зажим (%s)\n", c.name.c_str());
    const auto& q = c.q;
    const auto& W = c.W;
    check(W.clamp_boxes.size() == 2, "два бокса для clamp=Y");
    std::size_t total = 0, counted = 0;
    bool warp_free = true, weft_ok = true;
    const double b_hi = q.cy + q.clamp_edge, b_lo = q.cy - q.clamp_edge;
    for (const auto& F : W.fibers)
    {
        total += static_cast<std::size_t>(F.n_clamp_lo + F.n_clamp_hi);
        if (F.direction == 'x')
        {
            if (F.clamped || F.n_clamp_lo != 0 || F.n_clamp_hi != 0) warp_free = false;
            continue;
        }
        if (!F.clamped || F.n_clamp_lo < 1 || F.n_clamp_hi < 1) weft_ok = false;
        for (const auto& P : F.pts)
            if (P[1] <= b_lo || P[1] >= b_hi) ++counted;
    }
    check(counted == total, "число зажатых бусин = sum(n_clamp_lo + n_clamp_hi)");
    check(warp_free, "ни одна бусина основы не зажата");
    check(weft_ok, "у каждой нити утка зажато хотя бы по одной бусине на каждом конце");
    std::printf("    зажатых бусин всего %zu; у нити утка: lo=%d hi=%d\n", total, W.fibers.back().n_clamp_lo,
                W.fibers.back().n_clamp_hi);
    const auto S = fiber::plain_weave_summary(W);
    check(S.clamped_beads == total, "сводка: clamped_beads совпадает");
    check(W.warnings.empty(), "нет WARN о бусинах на границе c_e");
}

static void warn(const std::string& what)
{
    std::printf("  WARN : %s\n", what.c_str());
}

// Временный каталог (образец: e2e::ScopedWorkDir, но без смены cwd).
struct TempDir
{
    std::filesystem::path dir;
    TempDir()
    {
        static int counter = 0;
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        dir = std::filesystem::temp_directory_path() /
              ("plain_weave_test_" + std::to_string(static_cast<long long>(stamp)) + "_" + std::to_string(counter++));
        std::filesystem::create_directories(dir);
    }
    ~TempDir()
    {
        std::error_code ec;
        std::filesystem::remove_all(dir, ec);
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
};

struct LayoutRow
{
    int object_id = 0;
    char direction = 'x';
    int index = 0;
    int clamped = 0;
    std::size_t first = 0, last = 0;
    int n_lo = 0, n_hi = 0;
};

// Записать .fib и раскладку, загрузить fib в Model; раскладка читается обратно из текста.
static std::vector<LayoutRow> write_and_load(const PlainWeave& W, fiber::Model& model)
{
    TempDir tmp;
    const auto fib = tmp.dir / "weave.fib";
    {
        std::ofstream os(fib);
        fiber::write_plain_weave_fib(os, W, "test_plain_weave");
    }
    std::ostringstream lay;
    fiber::write_plain_weave_layout(lay, W);
    fiber::load_fiber_file(fib.string(), model);
    std::vector<LayoutRow> rows;
    std::istringstream is(lay.str());
    std::string line;
    while (std::getline(is, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        LayoutRow r;
        double offset;
        ls >> r.object_id >> r.direction >> r.index >> offset >> r.clamped >> r.first >> r.last >> r.n_lo >> r.n_hi;
        rows.push_back(r);
    }
    return rows;
}

static bool within_ulp(double a, double b)
{
    return a == b || a == std::nextafter(b, a) || b == std::nextafter(a, b);
}

// --- 5. Круговой путь через парсер -------------------------------------------
static void test_roundtrip(const std::string& name, const PlainWeaveParams& q)
{
    std::printf("-- 5. круговой путь через парсер (%s)\n", name.c_str());
    const PlainWeave W = fiber::build_plain_weave(q);
    fiber::Model model;
    const auto rows = write_and_load(W, model);
    const auto& fb = model.fiber;
    const std::size_t nf = W.fibers.size();
    std::size_t nbeads = 0;
    for (const auto& F : W.fibers) nbeads += F.pts.size();
    check(rows.size() == nf, "раскладка: строка на каждую нить");
    check(fb.size() == nbeads, "число бусин в модели совпадает с построенным");

    // Координаты, object_id, масса.
    bool ok_x = true, ok_id = true, ok_layout = true;
    std::size_t first = 0;
    double sum_l0 = 0.0;
    for (std::size_t f = 0; f < nf && f < rows.size(); ++f)
    {
        const auto& F = W.fibers[f];
        const auto& R = rows[f];
        if (R.first != first || R.last != first + F.pts.size() - 1) ok_layout = false;
        for (std::size_t b = 0; b < F.pts.size(); ++b)
        {
            const std::size_t g = R.first + b;
            if (g >= fb.size()) { ok_x = false; continue; }
            for (std::size_t dd = 0; dd < 3; ++dd)
                if (!within_ulp(fb.x[dd][g], F.pts[b][dd]))
                {
                    if (ok_x)
                        std::printf("    %s бусина %zu: координата %zu: файл %.17g, построено %.17g\n",
                                    fiber_name(F).c_str(), b, dd, fb.x[dd][g], F.pts[b][dd]);
                    ok_x = false;
                }
            if (fb.object_id[g] != R.object_id) ok_id = false;
            if (b + 1 < F.pts.size()) sum_l0 += seg_len(F.pts[b], F.pts[b + 1]);
        }
        first += F.pts.size();
    }
    check(ok_layout, "раскладка: first_bead..last_bead подряд по нитям");
    check(ok_x, "координаты бусин совпадают (до 1 ulp)");
    check(ok_id, "object_id бусин = object_id нити из раскладки");

    std::size_t expect_bonds = 0;
    for (const auto& F : W.fibers) expect_bonds += F.pts.size() - 1;
    check(model.topology.bonds.size() == expect_bonds, "число связей = сумма (бусин - 1)");

    double msum = 0.0;
    for (double m : fb.mass) msum += m;
    const double pi = std::acos(-1.0);
    const double mexp = q.rho * pi * q.d * q.d / 4.0 * sum_l0;
    check(std::fabs(msum - mexp) <= 1e-12 * mexp, "сумма масс = rho*pi*d^2/4*sum(l0) (отн. 1e-12)");
    std::printf("    sum m = %.12g кг, формула %.12g кг, отн. разница %.3g\n", msum, mexp,
                std::fabs(msum - mexp) / mexp);

    check(model.materials.size() >= 1, "материал 0 загружен");
    if (!model.materials.empty())
    {
        check(model.materials[0].mu_ff == q.mu_ff && model.materials[0].mu_fs == q.mu_fs,
              "mu_ff, mu_fs материала равны заданным (" + std::to_string(q.mu_ff) + ", " + std::to_string(q.mu_fs) +
                  ")");
    }

    // Зажим.
    const auto& mask = model.clamps.fiber_mask;
    if (q.clamp == fiber::PlainWeaveClamp::None)
    {
        check(mask.empty(), "clamp=None: fiber_mask пуст (empty)");
        return;
    }
    check(mask.size() == fb.size(), "clamp: размер fiber_mask == nfiber");
    if (mask.size() != fb.size()) return;
    const bool cy_on = q.clamp == fiber::PlainWeaveClamp::Y || q.clamp == fiber::PlainWeaveClamp::XY;
    const bool cx_on = q.clamp == fiber::PlainWeaveClamp::X || q.clamp == fiber::PlainWeaveClamp::XY;
    std::vector<std::uint8_t> expect(fb.size(), 0);
    std::size_t layout_total = 0, nonzero = 0;
    bool ok_ends = true, ok_none_x = true, ok_value = true, ok_mask = true;
    first = 0;
    for (std::size_t f = 0; f < nf; ++f)
    {
        const auto& F = W.fibers[f];
        const bool yfib = F.direction == 'y';
        const bool on = (yfib && cy_on) || (!yfib && cx_on);
        const int axis = yfib ? 1 : 0;
        const double c = yfib ? q.cy : q.cx;
        int lo = 0, hi = 0;
        for (std::size_t b = 0; b < F.pts.size(); ++b)
        {
            const double v = F.pts[b][axis];
            if (on && (v <= c - q.clamp_edge || v >= c + q.clamp_edge))
            {
                expect[first + b] = interactions::kClampAll;
                (v <= c - q.clamp_edge ? lo : hi)++;
            }
        }
        if (on && (lo < 1 || hi < 1)) ok_ends = false;
        if (!on && (rows[f].n_lo != 0 || rows[f].n_hi != 0)) ok_none_x = false;
        layout_total += static_cast<std::size_t>(rows[f].n_lo + rows[f].n_hi);
        first += F.pts.size();
    }
    for (std::size_t g = 0; g < mask.size(); ++g)
    {
        if (mask[g] != 0) ++nonzero;
        if (mask[g] != 0 && mask[g] != interactions::kClampAll) ok_value = false;
        if (mask[g] != expect[g]) ok_mask = false;
    }
    // Нити, не зажимаемые вовсе, не содержат зажатых бусин.
    first = 0;
    for (const auto& F : W.fibers)
    {
        const bool yfib = F.direction == 'y';
        if (!((yfib && cy_on) || (!yfib && cx_on)))
            for (std::size_t b = 0; b < F.pts.size(); ++b)
                if (mask[first + b] != 0) ok_none_x = false;
        first += F.pts.size();
    }
    check(ok_mask, "маска == ожидаемая (строится из PlainWeave независимо от парсера)");
    check(ok_value, "ненулевые значения маски == kClampAll");
    check(nonzero == layout_total, "число зажатых бусин = sum(n_clamp_lo + n_clamp_hi) раскладки");
    check(ok_ends, "у каждой зажимаемой нити >= 1 зажатая бусина на каждом конце");
    check(ok_none_x, "у незажимаемых нитей зажатых бусин нет");
    std::printf("    зажатых бусин %zu\n", nonzero);
}

// --- 6. Плотность -------------------------------------------------------------
static void test_density(const std::string& name, const PlainWeaveParams& q)
{
    std::printf("-- 6. плотность (%s)\n", name.c_str());
    const PlainWeave W = fiber::build_plain_weave(q);
    fiber::Model model;
    write_and_load(W, model);
    const auto S = fiber::plain_weave_summary(W);
    double msum = 0.0;
    for (double m : model.fiber.mass) msum += m;
    const double Lx = (q.ny - 1) * q.pitch, Ly = (q.nx - 1) * q.pitch;
    const double rho_model = msum / (Lx * Ly);
    check(std::fabs(S.areal_density - rho_model) <= 1e-12 * rho_model, "areal_density сводки = sum m (модель) / (Lx*Ly)");
    std::printf("    areal_density = %.9g кг/м^2, из модели %.9g\n", S.areal_density, rho_model);
}

static void test_density_twaron()
{
    std::printf("-- 6. плотность (twaron_ct717, только генерация)\n");
    PlainWeaveParams q;
    fiber::apply_plain_weave_preset("twaron_ct717", q);
    const PlainWeave W = fiber::build_plain_weave(q);
    const auto S = fiber::plain_weave_summary(W);
    check(std::isfinite(S.areal_density) && S.areal_density > 0.0, "плотность конечна и положительна");
    const double pi = std::acos(-1.0);
    const double Aw = pi * q.d * q.d / 4.0;
    // Номинал: две системы нитей, по одной нити на шаг p в каждой.
    const double nominal = 2.0 * q.rho * Aw / q.pitch;
    // Длина нити на полуволну (хорды) и коэффициент удлинения l_полуволны / p.
    const auto& F = W.fibers[0];
    double half_len = 0.0;
    for (int s = 0; s < q.seg; ++s)
        half_len += seg_len(F.pts[F.vertex_bead[1] - s], F.pts[F.vertex_bead[1] - s - 1]);
    const double crimp = half_len / q.pitch;
    std::printf("    бусин %zu, нитей %zu, плотность %.6g кг/м^2 (ожидание ТЗ 0.29 +- 0.03)\n", S.beads, S.fibers,
                S.areal_density);
    std::printf("    разложение: номинал 2*rho*(pi d^2/4)/p = %.6g кг/м^2\n", nominal);
    std::printf("    коэффициент удлинения на полуволну (хорды) = %.6g (l0 ~ %.4g м, k = %d)\n", crimp, half_len / q.seg,
                q.seg);
    std::printf("    номинал * удлинение = %.6g кг/м^2\n", nominal * crimp);
    // Вклад выступов: длина нитей между крайними пересечениями (окно) и за ними.
    const double Lx = (q.ny - 1) * q.pitch, Ly = (q.nx - 1) * q.pitch;
    double l_win = 0.0, l_over = 0.0;
    for (const auto& Fi : W.fibers)
        for (std::size_t b = 0; b + 1 < Fi.pts.size(); ++b)
        {
            const double l = seg_len(Fi.pts[b], Fi.pts[b + 1]);
            (b >= Fi.vertex_bead.front() && b < Fi.vertex_bead.back() ? l_win : l_over) += l;
        }
    std::printf("    e_x = %.4g м (основа, %d нитей), e_y = %.4g м (уток, %d нитей); окно %.4g x %.4g м\n",
                fiber::plain_weave_overhang_x(q), q.nx, fiber::plain_weave_overhang_y(q), q.ny, Lx, Ly);
    std::printf("    длина нитей: в окне %.5g м, в выступах %.5g м; плотность окна %.6g, вклад выступов %.6g кг/м^2\n",
                l_win, l_over, q.rho * Aw * l_win / (Lx * Ly), q.rho * Aw * l_over / (Lx * Ly));
    const double diff = S.areal_density - 0.29;
    if (std::fabs(diff) > 0.03)
        warn("плотность Twaron " + std::to_string(S.areal_density) + " вне 0.29 +- 0.03 (находка, не FAIL)");
}

// --- 7. Контакт на исходной геометрии -----------------------------------------
struct ContactMeasure
{
    std::size_t pair_count = 0, segment_pairs = 0, contacts = 0;
    double max_pen = 0.0;
    double per_crossing = 0.0;
    bool finite = true;
    std::size_t crossing_contacts = 0; // оценка: контакты сегментов нитей разных направлений
    std::size_t parallel_contacts = 0; // оценка: контакты между параллельными нитями
};

static ContactMeasure measure_contact(const PlainWeaveParams& q)
{
    const PlainWeave W = fiber::build_plain_weave(q);
    fiber::Model model;
    write_and_load(W, model);
    rhs::CpuExecutor exec(2);
    interactions::FiberContactSearchResult search;
    interactions::build_fiber_contact_search(exec, model.fiber, model.topology, model.materials, search);
    const auto diag = interactions::compute_fiber_contact_diagnostics(exec, model.fiber, model.materials, search);
    ContactMeasure m;
    m.pair_count = diag.pair_count;
    m.segment_pairs = diag.segment_pair_count;
    m.contacts = diag.contact_count;
    m.max_pen = diag.max_penetration;
    m.per_crossing = static_cast<double>(diag.contact_count) / (static_cast<double>(q.nx) * q.ny);
    // Характеристика: контакты между нитями одного направления (соседние параллельные)
    // и разных (пересечения). object_id < nx -- основа.
    for (const auto& prim : search.primitives)
    {
        const auto g = interactions::compute_fiber_contact_geometry(model.fiber, model.materials, prim);
        if (!g.contact || g.rejected_by_ownership || g.degenerate) continue;
        const bool wa = model.fiber.object_id[prim.a0] < q.nx, wb = model.fiber.object_id[prim.b0] < q.nx;
        (wa == wb ? m.parallel_contacts : m.crossing_contacts)++;
    }
    m.finite = std::isfinite(m.max_pen) && std::isfinite(m.per_crossing) && std::isfinite(diag.min_dt > 0 ? 1.0 : 0.0);
    return m;
}

static void report_contact(const char* tag, const PlainWeaveParams& q, const std::string& label)
{
    const ContactMeasure m = measure_contact(q);
    std::printf("    %s %s: pair_count=%zu segment_pair_count=%zu contact_count=%zu max_pen=%.6g м (%.4g d) "
                "контактов на пересечение=%.4g\n",
                tag, label.c_str(), m.pair_count, m.segment_pairs, m.contacts, m.max_pen, m.max_pen / q.d,
                m.per_crossing);
    std::printf("      из них: разные направления (пересечения) %zu, параллельные нити %zu\n", m.crossing_contacts,
                m.parallel_contacts);
    check(m.finite, std::string(tag) + " " + label + ": диагностика конечна (нет NaN/Inf)");
    if (m.per_crossing > 1.0)
        warn(std::string(tag) + " " + label + ": " + std::to_string(m.per_crossing) +
             " контактов на пересечение > 1 (седло, ТЗ §2)");
}

static void test_contact()
{
    std::printf("-- 7. контакт на исходной геометрии (пресет test)\n");
    for (double fa : {0.5, 0.6})
    {
        PlainWeaveParams q = make_test_params();
        q.amp = fa * q.d;
        report_contact("7:", q, "A=" + std::to_string(fa) + "d");
    }
}

static void test_contact_twaron()
{
    std::printf("-- 7б. контакт на геометрии twaron_ct717 (12x12, clamp=None)\n");
    for (int k : {2, 4})
        for (double fa : {0.5, 0.6})
        {
            PlainWeaveParams q;
            fiber::apply_plain_weave_preset("twaron_ct717", q);
            q.nx = q.ny = 12;
            q.clamp = fiber::PlainWeaveClamp::None;
            q.clamp_edge_set = false;
            q.seg = k;
            q.amp = fa * q.d;
            const ContactMeasure m = measure_contact(q);
            std::printf("7б: k=%d A=%.1fd contacts=%zu per_crossing=%.4g max_pen/d=%.4g (pair_count=%zu "
                        "segment_pairs=%zu; пересечения %zu, параллельные %zu)\n",
                        k, fa, m.contacts, m.per_crossing, m.max_pen / q.d, m.pair_count, m.segment_pairs,
                        m.crossing_contacts, m.parallel_contacts);
            check(m.finite, "7б: k=" + std::to_string(k) + " A=" + std::to_string(fa) + "d: конечно (нет NaN/Inf)");
            if (m.per_crossing > 1.0)
                warn("7б: k=" + std::to_string(k) + " A=" + std::to_string(fa) + "d: " +
                     std::to_string(m.per_crossing) + " контактов на пересечение > 1 (седло)");
        }
}

static void run_case(const Case& c)
{
    std::printf("=== вариант %s ===\n", c.name.c_str());
    test_count(c);
    test_weave(c);
    test_no_penetration(c);
    test_equal_chords(c);
    if (c.q.profile == PlainWeaveProfile::Arcs) test_arcs_on_ellipses(c);
    if (c.q.clamp != fiber::PlainWeaveClamp::None) test_clamp_counts(c);
}

int main()
{
    std::vector<Case> cases;
    auto add = [&](const std::string& name, const PlainWeaveParams& q) {
        Case c;
        c.name = name;
        c.q = q;
        c.W = fiber::build_plain_weave(q);
        cases.push_back(std::move(c));
    };
    add("test (cos)", make_test_params());
    add("нечётный (cos, 3x5, A=0.75d, k=3, clamp=Y)", make_odd_params(PlainWeaveProfile::Cos));
    {
        PlainWeaveParams q = make_test_params();
        q.profile = PlainWeaveProfile::Arcs;
        add("test (arcs)", q);
    }
    add("нечётный (arcs, 3x5, A=0.75d, k=3, clamp=Y)", make_odd_params(PlainWeaveProfile::Arcs));

    for (const auto& c : cases) run_case(c);
    test_param_errors();

    // 5, 6: прогоны на test и нечётном варианте, варианты трения, clamp=XY.
    {
        PlainWeaveParams q = make_test_params();
        q.mu_ff = 0.3;
        q.mu_fs = 0.2;
        q.mu_fs_set = true;
        test_roundtrip("test, mu_ff=0.3 mu_fs=0.2", q);
        q.mu_ff = 0.0;
        test_roundtrip("test, только mu_fs=0.2 (mu_ff=0 явно)", q);
        test_roundtrip("test, без трения", make_test_params());
        PlainWeaveParams o = make_odd_params(PlainWeaveProfile::Cos);
        o.mu_ff = 0.3;
        o.mu_fs = 0.2;
        test_roundtrip("нечётный cos, clamp=Y", o);
        test_roundtrip("нечётный arcs, clamp=Y", make_odd_params(PlainWeaveProfile::Arcs));
        PlainWeaveParams xy = make_test_params();
        xy.clamp = fiber::PlainWeaveClamp::XY;
        // Диапазон c_e для 4x4, e = p/2: [1.5p + d/2, 1.5p + p/2] = [1.75, 2.0] мм; середина.
        xy.clamp_edge = 1.5 * xy.pitch + 0.375 * xy.pitch;
        xy.clamp_edge_set = true;
        test_roundtrip("test, clamp=XY", xy);
        test_density("test", make_test_params());
        test_density("нечётный cos, clamp=Y", make_odd_params(PlainWeaveProfile::Cos));
    }
    test_density_twaron();
    test_contact();
    test_contact_twaron();

    std::printf("\n%s (failures = %d)\n", failures == 0 ? "ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ" : "ЕСТЬ ОШИБКИ", failures);
    return failures == 0 ? 0 : 1;
}
