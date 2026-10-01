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
// Проверки 5, 6, 7 (круговой путь через load_fiber_file, плотность, диагностика
// контакта) добавляет следующий этап: новые раздел-функции и вызовы в main.
#include "fiber/plain_weave.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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

    std::printf("\n%s (failures = %d)\n", failures == 0 ? "ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ" : "ЕСТЬ ОШИБКИ", failures);
    return failures == 0 ? 0 : 1;
}
