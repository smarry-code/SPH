// plain_weave.hpp
// Генератор геометрии тканого полотна ПОЛОТНЯНОГО переплетения (plain weave)
// в формате файла нитей проекта (.fib, fiber::load_fiber_file).
// ТЗ: weave/TZ_PLAIN_WEAVE_GENERATOR.md (этап 1: геометрия, запись, сводка).
//
// Header-only. НЕ зависит от dim: не включает particle_system.hpp, clamp.hpp,
// model.hpp -- только стандартную библиотеку. Бокс закрепления -- шесть чисел.
//
// ГЕОМЕТРИЯ (ТЗ §3). Плоскость ткани x-y, толщина по z, центр (cx, cy, cz).
//   Основа j = 0..nx-1 идёт вдоль x на y_j = cy + (j - (nx-1)/2)*p,
//   уток   i = 0..ny-1 идёт вдоль y на x_i = cx + (i - (ny-1)/2)*p.
//   Основа j: z = cz + A*(-1)^j * cos(pi*(x - x_0)/p),
//   уток   i: z = cz - A*(-1)^i * cos(pi*(y - y_0)/p).
//   В пересечении (i, j) основа на cz + A*(-1)^(i+j), уток на cz - A*(-1)^(i+j).
//   Профиль "arcs": ось собрана из полуэллипсов с центрами в вершинах
//   (x_0 + m*p), полуось вдоль нити p/2, по z -- A; на стыке z = cz.
//
// ДИСКРЕТИЗАЦИЯ. Бусины стоят ТОЧНО в каждом пересечении (вершине): координаты
// вершины берутся из формул позиций, а не из накопленной суммы. Между соседними
// вершинами ровно k сегментов РАВНОЙ ДЛИНЫ ХОРДЫ (трактовка ТЗ §3.4: проверка 4
// требует равных l0; "равные дуги" хорды не уравнивают). Алгоритм: на полуволне
// (параметр T в [0, 1]) методом стрельбы подбирается длина хорды L: от вершины
// k-1 раз берётся точка кривой на евклидовом расстоянии L (бисекция по
// параметру), L подгоняется бисекцией так, чтобы расстояние от последней точки
// до конечной вершины тоже было L. Профиль полуволны считается один раз и
// переиспользуется со сдвигом и знаком. Выступ e делится на
// ne = round(e/(p/k)) (при e > 0 не меньше 1; при e == 0 выступа нет, ne = 0:
// ни точек выступа, ни дубля крайней бусины) равных хорд от крайней вершины до точки кривой с
// координатой вдоль нити ровно x_0 - e (x_end + e); выступ симметричен
// относительно вершины, поэтому считается один раз на направление.
//
// ЗАЖИМ (ТЗ §3.5): по границе c_e пишутся два CLAMP_BOX на зажимаемое
// направление; b_hi = c + c_e и b_lo = c - c_e вычисляются один раз и
// используются и в боксе, и в подсчёте зажатых бусин (замкнутая проверка).
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fiber
{

enum class PlainWeaveProfile
{
    Cos,
    Arcs
};

enum class PlainWeaveClamp
{
    None,
    X,
    Y,
    XY
};

struct PlainWeaveParams
{
    double d = 0.5e-3;     // диаметр нити, м
    double pitch = 1.0e-3; // p: шаг между осями параллельных нитей, м
    double amp = 0.25e-3;  // A: амплитуда волнистости по z, м
    int nx = 4;            // число нитей основы (вдоль x)
    int ny = 4;            // число нитей утка (вдоль y)
    int seg = 4;           // k: сегментов на полуволну, >= 2
    // Выступ за крайнюю поперечную нить; отрицательное значение = "не задано",
    // тогда p/2 (конец нити на середине полуволны).
    double overhang_x = -1.0;
    double overhang_y = -1.0;
    PlainWeaveProfile profile = PlainWeaveProfile::Cos;
    double cx = 0.0, cy = 0.0, cz = 0.0;
    PlainWeaveClamp clamp = PlainWeaveClamp::None;
    double clamp_edge = 0.0;       // c_e, м
    bool clamp_edge_set = false;   // c_e задан (пресетом или явно)
    // Материал (по умолчанию -- как пресет test, демо "4 слоя").
    double E = 5.0e7;
    double nu = 0.3;
    double rho = 1200.0;
    double eps_break = 0.15;
    double mu_ff = 0.0;
    double mu_fs = 0.0;
    bool mu_fs_set = false; // для CLI; запись mu_fs решает значение (см. write_plain_weave_fib)
    bool allow_penetration = false;
};

// Выступ по направлению: не задан -> p/2.
inline double plain_weave_overhang_x(const PlainWeaveParams& q)
{
    return q.overhang_x < 0.0 ? 0.5 * q.pitch : q.overhang_x;
}
inline double plain_weave_overhang_y(const PlainWeaveParams& q)
{
    return q.overhang_y < 0.0 ? 0.5 * q.pitch : q.overhang_y;
}

// Три пресета ТЗ §4.3. CLI применяет пресет, затем явные параметры.
inline void apply_plain_weave_preset(const std::string& name, PlainWeaveParams& q)
{
    if (name == "test")
    {
        q.d = 0.5e-3;
        q.pitch = 1.0e-3;
        q.amp = 0.5 * q.d;
        q.nx = 4;
        q.ny = 4;
        q.seg = 4;
        q.overhang_x = q.overhang_y = -1.0;
        q.profile = PlainWeaveProfile::Cos;
        q.clamp = PlainWeaveClamp::None;
        q.clamp_edge_set = false;
        q.E = 5.0e7;
        q.nu = 0.3;
        q.rho = 1200.0;
        q.eps_break = 0.15;
        q.mu_ff = 0.0;
        q.mu_fs = 0.0;
    }
    else if (name == "circle_draw")
    {
        // Скрипт Circle_draw.txt: полуэллипсы 0.1 и 0.2 мм, шаг центров 0.4 мм.
        // d = 2A = 0.2 мм (касание, ТЗ §10 п.2) -- отступление от подписи
        // s_d = 0.4 мм в скрипте: при A = 0.1 мм такие нити пересекались бы.
        q.d = 0.2e-3;
        q.pitch = 0.4e-3;
        q.amp = 0.1e-3;
        q.nx = 132;
        q.ny = 132;
        q.seg = 4;
        q.overhang_x = q.overhang_y = -1.0; // p/2
        q.profile = PlainWeaveProfile::Arcs;
        q.clamp = PlainWeaveClamp::None;
        q.clamp_edge_set = false;
        q.E = 5.0e7;
        q.nu = 0.3;
        q.rho = 1200.0;
        q.eps_break = 0.15;
        q.mu_ff = 0.0;
        q.mu_fs = 0.0;
    }
    else if (name == "twaron_ct717")
    {
        q.d = 0.39e-3;     // оценка: круг площадью 0.12 мм^2
        q.pitch = 1.18e-3; // оценка: 8.5 нит/см
        q.amp = 0.5 * q.d;
        q.nx = 102;
        q.ny = 102; // окно 120 мм
        q.seg = 4;
        q.profile = PlainWeaveProfile::Cos;
        q.E = 8.5e10;       // оценка: середина динамического диапазона 78-97 ГПа
        q.nu = 0.3;
        q.rho = 1440.0;     // оценка: справочная плотность арамида, не из статьи
        q.eps_break = 0.03; // оценка: eps_b ~ 2 U_T / sigma_T
        q.mu_ff = 0.0;
        q.mu_fs = 0.0;
        q.overhang_x = -1.0; // p/2 (с учётом --pitch, если задан)
        q.overhang_y = 5.0e-3; // оценка: ширина зажима в статье не задана
        q.clamp = PlainWeaveClamp::Y;
        q.clamp_edge = 60.0e-3; // оценка: кромка окна статьи
        q.clamp_edge_set = true;
    }
    else
    {
        throw std::invalid_argument("apply_plain_weave_preset: неизвестный пресет '" + name +
                                    "' (допустимы: test, circle_draw, twaron_ct717)");
    }
}

struct PlainWeaveFiber
{
    int object_id = 0;
    char direction = 'x'; // 'x' основа, 'y' уток
    int index = 0;
    double offset = 0.0; // y_j для основы, x_i для утка
    bool clamped = false;
    std::size_t n_clamp_lo = 0; // зажатых бусин у начала нити
    std::size_t n_clamp_hi = 0; // зажатых бусин у конца нити
    std::vector<std::array<double, 3>> pts;
    std::vector<std::size_t> vertex_bead; // индексы бусин в пересечениях по порядку вдоль нити
};

struct ClampBox
{
    double lo[3] = {0, 0, 0};
    double hi[3] = {0, 0, 0};
};

struct PlainWeave
{
    PlainWeaveParams params;
    std::vector<PlainWeaveFiber> fibers; // сначала основа (nx), затем уток (ny)
    std::vector<ClampBox> clamp_boxes;
    std::vector<std::string> warnings;
};

namespace plain_weave_detail
{

// Предохранители: общее число бусин (оценка до построения), сегментов на
// полуволну и хорд одного выступа.
constexpr double plain_weave_max_beads = 5.0e7;
constexpr int plain_weave_max_seg = 10000;
constexpr double plain_weave_max_overhang_chords = 1.0e7;

// Кривая полуволны с локальным параметром T (вершины в целых T).
// Локально: вершина T=0 на z=+A, T=1 на z=-A; на следующей полуволне знак меняется.
struct Curve
{
    PlainWeaveProfile type;
    double p, A;

    // Точка внутри одной полуволны, t in [0, 1]: (s вдоль нити от вершины, z).
    void half(double t, double& s, double& z) const
    {
        if (type == PlainWeaveProfile::Cos)
        {
            const double pi = std::acos(-1.0);
            s = p * t;
            z = A * std::cos(pi * t);
        }
        else
        {
            // Полуэллипсы с центрами в вершинах, угол эллипса по каждой половинке.
            const double pi = std::acos(-1.0);
            if (t <= 0.5)
            {
                const double th = pi * t; // 0..pi/2
                s = 0.5 * p * std::sin(th);
                z = A * std::cos(th);
            }
            else
            {
                const double ph = pi * (1.0 - t); // pi/2..0
                s = p - 0.5 * p * std::sin(ph);
                z = -A * std::cos(ph);
            }
        }
    }

    // Общая точка кривой: along = w*p + s, z со знаком (-1)^w.
    void at(double T, double& along, double& z) const
    {
        const double w = std::floor(T);
        const double t = T - w;
        double s, zz;
        half(t, s, zz);
        along = w * p + s;
        const long long wi = static_cast<long long>(w);
        const bool odd = ((wi % 2) + 2) % 2 == 1;
        z = odd ? -zz : zz;
    }

    double dist(double Ta, double Tb) const
    {
        double a1, z1, a2, z2;
        at(Ta, a1, z1);
        at(Tb, a2, z2);
        return std::hypot(a2 - a1, z2 - z1);
    }
};

// Следующая точка на расстоянии L от точки T в направлении dir (+1/-1).
// Возвращает false, если за пределом Tb.
inline bool next_at_distance(const Curve& C, double T, double L, double dir, double Tb, double& out)
{
    const double step = 1.0 / 4096.0;
    double lo = T, hi = T;
    bool found = false;
    for (int guard = 0; guard < 4096 * 64; ++guard)
    {
        hi = lo + dir * step;
        if ((hi - Tb) * dir > 0.0)
        {
            // Проверим сам предел Tb: если и там расстояние < L, нет решения.
            if (C.dist(T, Tb) < L) return false;
            hi = Tb;
            found = true;
            break;
        }
        if (C.dist(T, hi) >= L)
        {
            found = true;
            break;
        }
        lo = hi;
    }
    if (!found) return false;
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (mid == lo || mid == hi) break;
        if (C.dist(T, mid) >= L)
            hi = mid;
        else
            lo = mid;
    }
    out = hi;
    return true;
}

// Параметры T точек n равных хорд от Ta до Tb (включая концы), стрельба по L.
inline std::vector<double> equal_chords(const Curve& C, double Ta, double Tb, int n)
{
    const double dir = Tb > Ta ? 1.0 : -1.0;
    // Верхняя оценка: длина ломаной по мелкой сетке.
    double len = 0.0;
    {
        const int M = 4096;
        double prev = Ta;
        for (int i = 1; i <= M; ++i)
        {
            const double cur = Ta + (Tb - Ta) * i / M;
            len += C.dist(prev, cur);
            prev = cur;
        }
    }
    auto walk = [&](double L, std::vector<double>* out) -> double {
        double T = Ta;
        if (out)
        {
            out->clear();
            out->push_back(Ta);
        }
        for (int s = 1; s < n; ++s)
        {
            double Tn;
            if (!next_at_distance(C, T, L, dir, Tb, Tn)) return -L; // L слишком велико
            T = Tn;
            if (out) out->push_back(T);
        }
        return C.dist(T, Tb) - L;
    };
    double lo = 0.0, hi = len + 1e-300;
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (mid == lo || mid == hi) break;
        if (walk(mid, nullptr) > 0.0)
            lo = mid;
        else
            hi = mid;
    }
    std::vector<double> Ts;
    walk(0.5 * (lo + hi), &Ts);
    Ts.push_back(Tb);
    // Стрельба обязана дать ровно n хорд (n + 1 точек); иначе -- сбой алгоритма.
    if (Ts.size() != static_cast<std::size_t>(n) + 1)
        throw std::runtime_error("plain_weave: не удалось построить равные хорды (стрельба дала " +
                                 std::to_string(Ts.size()) + " точек вместо " + std::to_string(n + 1) + ")");
    return Ts;
}

// Профиль выступа e от вершины со знаком +A: (u_m, z_m), m = 0..n_e; u_0 = 0,
// z_0 = A, u_{n_e} = e точно.
struct Overhang
{
    std::vector<double> u, z;
};

// Число хорд выступа e: 0 при e == 0 (выступа нет), иначе max(1, round(e/(p/k))).
// Вызывающий гарантирует e/(p/k) <= plain_weave_max_overhang_chords.
inline long long overhang_chords(double e, double p, int k)
{
    if (e == 0.0) return 0;
    return std::max<long long>(1, std::llround(e / (p / k)));
}

inline Overhang make_overhang(const Curve& C, double e, int k)
{
    if (e == 0.0)
    {
        // Нулевой выступ: единственная «точка» -- сама вершина, новых бусин нет.
        Overhang o0;
        o0.u.push_back(0.0);
        o0.z.push_back(C.A);
        return o0;
    }
    const int ne = static_cast<int>(overhang_chords(e, C.p, k));
    // T_end: along(T_end) = e (along монотонно растёт по T).
    double lo = 0.0, hi = e / C.p + 2.0;
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        double al, zz;
        C.at(mid, al, zz);
        if (al < e)
            lo = mid;
        else
            hi = mid;
    }
    const double Tend = 0.5 * (lo + hi);
    const std::vector<double> Ts = equal_chords(C, 0.0, Tend, ne);
    Overhang o;
    for (std::size_t m = 0; m < Ts.size(); ++m)
    {
        double al, zz;
        C.at(Ts[m], al, zz);
        o.u.push_back(al);
        o.z.push_back(zz);
    }
    o.u.front() = 0.0;
    o.z.front() = C.A;
    o.u.back() = e;
    return o;
}

} // namespace plain_weave_detail

inline PlainWeave build_plain_weave(const PlainWeaveParams& q)
{
    using namespace plain_weave_detail;
    PlainWeave W;
    W.params = q;
    const double p = q.pitch, d = q.d, A = q.amp;
    const int k = q.seg;

    // Конечность и допустимость всех вещественных параметров (генератор не
    // должен писать файл, который парсер или физика отвергнут).
    {
        auto fin = [](double v, const char* name) {
            if (!std::isfinite(v))
                throw std::invalid_argument(std::string("plain_weave: параметр ") + name + " должен быть конечным числом");
        };
        fin(q.d, "d");
        fin(q.pitch, "pitch");
        fin(q.amp, "amp");
        fin(q.overhang_x, "overhang_x");
        fin(q.overhang_y, "overhang_y");
        fin(q.cx, "cx");
        fin(q.cy, "cy");
        fin(q.cz, "cz");
        if (q.clamp_edge_set) fin(q.clamp_edge, "clamp_edge");
        fin(q.E, "E");
        fin(q.nu, "nu");
        fin(q.rho, "rho");
        fin(q.eps_break, "eps_break");
        fin(q.mu_ff, "mu_ff");
        fin(q.mu_fs, "mu_fs");
        if (!(q.E > 0.0)) throw std::invalid_argument("plain_weave: модуль E должен быть > 0");
        if (!(q.rho > 0.0)) throw std::invalid_argument("plain_weave: плотность rho должна быть > 0");
        if (!(q.eps_break > 0.0)) throw std::invalid_argument("plain_weave: eps_break должен быть > 0");
        if (q.mu_ff < 0.0) throw std::invalid_argument("plain_weave: mu_ff должен быть >= 0");
        if (q.mu_fs < 0.0) throw std::invalid_argument("plain_weave: mu_fs должен быть >= 0");
    }

    if (!(d > 0.0)) throw std::invalid_argument("plain_weave: диаметр d должен быть > 0");
    if (!(p > 0.0)) throw std::invalid_argument("plain_weave: шаг p должен быть > 0");
    if (!(A > 0.0)) throw std::invalid_argument("plain_weave: амплитуда A должна быть > 0");
    if (k < 2) throw std::invalid_argument("plain_weave: число сегментов на полуволну k должно быть >= 2");
    if (k > plain_weave_max_seg)
        throw std::invalid_argument("plain_weave: seg слишком велико (максимум " + std::to_string(plain_weave_max_seg) + ")");
    if (q.nx < 1 || q.ny < 1) throw std::invalid_argument("plain_weave: nx, ny должны быть >= 1");
    const double ex = plain_weave_overhang_x(q), ey = plain_weave_overhang_y(q);
    if (ex < 0.0 || ey < 0.0) throw std::invalid_argument("plain_weave: выступы должны быть >= 0");
    // Размер: диапазон перед приведением к целым и оценка числа бусин ДО построения.
    {
        const double chx = ex / (p / k), chy = ey / (p / k);
        if (!(chx <= plain_weave_max_overhang_chords) || !(chy <= plain_weave_max_overhang_chords))
            throw std::invalid_argument("plain_weave: выступ слишком велик относительно шага (больше 1e7 хорд)");
        const double nex = static_cast<double>(overhang_chords(ex, p, k));
        const double ney = static_cast<double>(overhang_chords(ey, p, k));
        const double bx = static_cast<double>(q.nx) * (static_cast<double>(q.ny - 1) * k + 1.0 + 2.0 * nex);
        const double by = static_cast<double>(q.ny) * (static_cast<double>(q.nx - 1) * k + 1.0 + 2.0 * ney);
        if (!(bx + by <= plain_weave_max_beads))
        {
            std::ostringstream s;
            s << "plain_weave: оценка числа бусин " << bx + by << " превышает предел " << plain_weave_max_beads;
            throw std::invalid_argument(s.str());
        }
    }

    auto penalty = [&](bool bad, const std::string& msg) {
        if (!bad) return;
        if (!q.allow_penetration) throw std::invalid_argument(msg + " (разрешить: --allow-penetration)");
        const std::string w = "WARN " + msg;
        W.warnings.push_back(w);
        std::fprintf(stderr, "%s\n", w.c_str());
    };
    {
        std::ostringstream s1, s2;
        s1 << "plain_weave: p < d (" << p << " < " << d << "): соседние параллельные нити перекрываются";
        s2 << "plain_weave: 2A < d (" << 2.0 * A << " < " << d << "): нити проникают друг в друга в пересечении";
        penalty(p < d, s1.str());
        penalty(2.0 * A < d, s2.str());
    }

    // Положения нитей.
    auto coord = [&](double c, int idx, int n) { return c + (idx - 0.5 * (n - 1)) * p; };

    // Границы зажима (ТЗ §3.5).
    const bool clamp_y = q.clamp == PlainWeaveClamp::Y || q.clamp == PlainWeaveClamp::XY;
    const bool clamp_x = q.clamp == PlainWeaveClamp::X || q.clamp == PlainWeaveClamp::XY;
    if (q.clamp != PlainWeaveClamp::None)
    {
        if (!q.clamp_edge_set)
            throw std::invalid_argument("plain_weave: --clamp задан, но не задана граница зажима --clamp-edge (c_e)");
        const double ce = q.clamp_edge;
        if (clamp_y)
        {
            const double yl = coord(q.cy, q.nx - 1, q.nx) - q.cy;
            if (ce < yl + 0.5 * d || ce > yl + ey)
            {
                std::ostringstream s;
                s << "plain_weave: c_e = " << ce << " вне допустимого диапазона для зажима по y: ["
                  << yl + 0.5 * d << ", " << yl + ey << "]";
                throw std::invalid_argument(s.str());
            }
        }
        if (clamp_x)
        {
            const double xl = coord(q.cx, q.ny - 1, q.ny) - q.cx;
            if (ce < xl + 0.5 * d || ce > xl + ex)
            {
                std::ostringstream s;
                s << "plain_weave: c_e = " << ce << " вне допустимого диапазона для зажима по x: ["
                  << xl + 0.5 * d << ", " << xl + ex << "]";
                throw std::invalid_argument(s.str());
            }
        }
    }

    const Curve C{q.profile, p, A};
    // Полуволна: локальные (s_m, z_m), m = 0..k-1 (m = k -- следующая вершина).
    std::vector<double> hs(k), hz(k);
    {
        const std::vector<double> Ts = equal_chords(C, 0.0, 1.0, k);
        for (int m = 0; m < k; ++m)
        {
            double al, zz;
            C.at(Ts[m], al, zz);
            hs[m] = al;
            hz[m] = zz;
        }
        hs[0] = 0.0;
        hz[0] = A;
    }
    const Overhang ovx = make_overhang(C, ex, k); // выступ нитей основы (вдоль x)
    const Overhang ovy = make_overhang(C, ey, k); // выступ нитей утка (вдоль y)

    // Нить вдоль направления dir: nv вершин с координатами along(idx), sigma0 -- знак вершины 0.
    auto build_fiber = [&](char dir, int idx, int object_id, double offset, int nv, double c_al,
                           int n_al, double sigma0, const Overhang& ov) {
        PlainWeaveFiber F;
        F.object_id = object_id;
        F.direction = dir;
        F.index = idx;
        F.offset = offset;
        auto push = [&](double along, double z) {
            if (dir == 'x')
                F.pts.push_back({along, offset, z});
            else
                F.pts.push_back({offset, along, z});
        };
        const double v0 = coord(c_al, 0, n_al);
        const int ne = static_cast<int>(ov.u.size()) - 1; // 0 при нулевом выступе
        for (int m = ne; m >= 1; --m) push(v0 - ov.u[m], q.cz + sigma0 * ov.z[m]);
        const int Nhalf = nv - 1;
        for (int w = 0; w < Nhalf; ++w)
        {
            const double sg = (w % 2 == 0) ? sigma0 : -sigma0;
            const double vw = coord(c_al, w, n_al);
            F.vertex_bead.push_back(F.pts.size());
            push(vw, q.cz + sg * A);
            for (int m = 1; m < k; ++m) push(vw + hs[m], q.cz + sg * hz[m]);
        }
        const double sgN = (Nhalf % 2 == 0) ? sigma0 : -sigma0;
        const double vN = coord(c_al, Nhalf, n_al);
        F.vertex_bead.push_back(F.pts.size());
        push(vN, q.cz + sgN * A);
        for (int m = 1; m <= ne; ++m) push(vN + ov.u[m], q.cz + sgN * ov.z[m]);
        return F;
    };

    for (int j = 0; j < q.nx; ++j)
        W.fibers.push_back(build_fiber('x', j, j, coord(q.cy, j, q.nx), q.ny, q.cx, q.ny,
                                       (j % 2 == 0) ? 1.0 : -1.0, ovx));
    for (int i = 0; i < q.ny; ++i)
        W.fibers.push_back(build_fiber('y', i, q.nx + i, coord(q.cx, i, q.ny), q.nx, q.cy, q.nx,
                                       (i % 2 == 0) ? -1.0 : 1.0, ovy));

    // Зажим.
    if (q.clamp != PlainWeaveClamp::None)
    {
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (const auto& F : W.fibers)
            for (const auto& P : F.pts)
                for (int a = 0; a < 3; ++a)
                {
                    lo[a] = std::min(lo[a], P[a]);
                    hi[a] = std::max(hi[a], P[a]);
                }
        for (int a = 0; a < 3; ++a)
        {
            lo[a] -= d;
            hi[a] += d;
        }
        const double ce = q.clamp_edge;
        const double bx_hi = q.cx + ce, bx_lo = q.cx - ce;
        const double by_hi = q.cy + ce, by_lo = q.cy - ce;
        auto add_box = [&](int axis, double blo, double bhi) {
            ClampBox b1, b2;
            for (int a = 0; a < 3; ++a)
                b1.lo[a] = b2.lo[a] = lo[a], b1.hi[a] = b2.hi[a] = hi[a];
            b1.lo[axis] = bhi; // верхний конец: [b_hi, hi]
            b2.hi[axis] = blo; // нижний конец: [lo, b_lo]
            W.clamp_boxes.push_back(b1);
            W.clamp_boxes.push_back(b2);
        };
        if (clamp_x) add_box(0, bx_lo, bx_hi);
        if (clamp_y) add_box(1, by_lo, by_hi);

        auto in_box = [](const ClampBox& b, const std::array<double, 3>& P) {
            for (int a = 0; a < 3; ++a)
                if (P[a] < b.lo[a] || P[a] > b.hi[a]) return false;
            return true;
        };
        for (auto& F : W.fibers)
        {
            const bool clamped = (F.direction == 'y' && clamp_y) || (F.direction == 'x' && clamp_x);
            F.clamped = clamped;
            const int axis = F.direction == 'y' ? 1 : 0;
            const double blo = axis == 1 ? by_lo : bx_lo, bhi = axis == 1 ? by_hi : bx_hi;
            for (std::size_t b = 0; b < F.pts.size(); ++b)
            {
                const auto& P = F.pts[b];
                bool any = false;
                for (const auto& box : W.clamp_boxes) any = any || in_box(box, P);
                const bool lo_hit = clamped && P[axis] <= blo;
                const bool hi_hit = clamped && P[axis] >= bhi;
                if (lo_hit) ++F.n_clamp_lo;
                if (hi_hit) ++F.n_clamp_hi;
                if (any != (lo_hit || hi_hit))
                {
                    std::ostringstream s;
                    s << "plain_weave: нить " << F.object_id << " (" << F.direction << F.index << "), бусина " << b
                      << " (" << P[0] << ", " << P[1] << ", " << P[2] << "): "
                      << (any ? "попала в бокс закрепления, хотя не должна" : "не попала в бокс, хотя должна");
                    throw std::runtime_error(s.str());
                }
                if (clamped && std::fabs(std::fabs(P[axis] - (axis == 1 ? q.cy : q.cx)) - ce) < 1e-9 * p)
                {
                    std::ostringstream s;
                    s << "WARN нить " << F.object_id << " бусина " << b << " на границе c_e (" << P[0] << ", "
                      << P[1] << ", " << P[2] << "): результат зависит от округления, сдвиньте c_e";
                    W.warnings.push_back(s.str());
                    std::fprintf(stderr, "%s\n", s.str().c_str());
                }
            }
            if (clamped && (F.n_clamp_lo < 1 || F.n_clamp_hi < 1))
            {
                std::ostringstream s;
                s << "plain_weave: у зажимаемой нити " << F.object_id << " нет зажатой бусины на конце (lo="
                  << F.n_clamp_lo << ", hi=" << F.n_clamp_hi << ")";
                throw std::runtime_error(s.str());
            }
        }
    }
    return W;
}

// Окно полотна для плотности: Lx = x_{ny-1} - x_0 = (ny-1)p, Ly = (nx-1)p --
// между крайними нитями. При n = 1 по направлению окно между нитями вырождено
// (0), поэтому берётся ширина одной ячейки p.
inline void plain_weave_window(const PlainWeaveParams& q, double& Lx, double& Ly)
{
    Lx = q.ny > 1 ? (static_cast<double>(q.ny) - 1.0) * q.pitch : q.pitch;
    Ly = q.nx > 1 ? (static_cast<double>(q.nx) - 1.0) * q.pitch : q.pitch;
}

struct PlainWeaveSummary
{
    std::size_t fibers = 0, beads = 0, bonds = 0, clamped_beads = 0;
    double l0_min = 0.0, l0_max = 0.0;
    double thickness = 0.0;    // 2A + d
    double areal_density = 0.0; // кг/м^2: sum m / (L_x * L_y)
    double cost_ratio = 0.0;    // оценка: бусин / 340 (демо "4 слоя")
};

// Плотность: Σm / (L_x*L_y), Σm = rho*pi*d^2/4*Σl0 (как load_fiber_file),
// L_x = x_{ny-1}-x_0 = (ny-1)p, L_y = (nx-1)p -- окно между крайними нитями.
// При n = 1 размер окна 0, тогда по этому направлению берётся p.
inline PlainWeaveSummary plain_weave_summary(const PlainWeave& W)
{
    const auto& q = W.params;
    PlainWeaveSummary S;
    S.fibers = W.fibers.size();
    double sum_l0 = 0.0;
    S.l0_min = std::numeric_limits<double>::infinity();
    S.l0_max = 0.0;
    for (const auto& F : W.fibers)
    {
        S.beads += F.pts.size();
        S.clamped_beads += F.n_clamp_lo + F.n_clamp_hi;
        for (std::size_t b = 0; b + 1 < F.pts.size(); ++b)
        {
            const auto &a = F.pts[b], &c = F.pts[b + 1];
            const double l = std::sqrt((c[0] - a[0]) * (c[0] - a[0]) + (c[1] - a[1]) * (c[1] - a[1]) +
                                       (c[2] - a[2]) * (c[2] - a[2]));
            sum_l0 += l;
            S.l0_min = std::min(S.l0_min, l);
            S.l0_max = std::max(S.l0_max, l);
        }
    }
    S.bonds = S.beads - S.fibers;
    if (S.bonds == 0) S.l0_min = 0.0;
    S.thickness = 2.0 * q.amp + q.d;
    const double pi = std::acos(-1.0);
    const double mass = q.rho * pi * q.d * q.d / 4.0 * sum_l0;
    double Lx = 0.0, Ly = 0.0;
    plain_weave_window(q, Lx, Ly);
    S.areal_density = mass / (Lx * Ly);
    S.cost_ratio = static_cast<double>(S.beads) / 340.0;
    return S;
}

inline void print_plain_weave_summary(std::FILE* f, const PlainWeave& W)
{
    const PlainWeaveSummary S = plain_weave_summary(W);
    std::fprintf(f, "нитей: %zu (основа %d, уток %d)\n", S.fibers, W.params.nx, W.params.ny);
    std::fprintf(f, "бусин: %zu, связей: %zu, зажатых бусин: %zu\n", S.beads, S.bonds, S.clamped_beads);
    std::fprintf(f, "l0: min %.9g м, max %.9g м\n", S.l0_min, S.l0_max);
    std::fprintf(f, "толщина полотна 2A+d: %.9g м\n", S.thickness);
    std::fprintf(f, "поверхностная плотность: %.9g кг/м^2\n", S.areal_density);
    std::fprintf(f, "оценка стоимости (бусин/340, демо \"4 слоя\"): %.4g\n", S.cost_ratio);
}

// Необязательное дополнение .fib (случайная раскройка, plain_weave_random_cut.hpp):
// стежки BOND материала 1 и уточнения заголовка. По умолчанию пусто -- вывод прежний.
struct PlainWeaveFibExtra
{
    // Материал 1 стежков (пишется, если pairs != nullptr): формат как MATERIAL 0.
    double stitch_E = 0.0, stitch_nu = 0.0, stitch_d = 0.0, stitch_rho = 0.0, stitch_eps_break = 0.0;
    const std::vector<std::pair<std::size_t, std::size_t>>* pairs = nullptr; // BOND i j 1 (глобальные индексы)
    std::vector<std::string> comments; // доп. строки комментария заголовка (без '#')
    double window_area = 0.0;          // > 0: площадь окна для поверхностной плотности вместо plain_weave_window
};

// Запись .fib (формат: src/fiber_input.cpp). Правило записи трения (ТЗ §3.6):
// если mu_fs != 0 пишутся оба поля (mu_ff = 0 явно, если не задан); иначе при
// mu_ff != 0 -- только mu_ff; если оба 0, необязательные поля не пишутся.
inline void write_plain_weave_fib(std::ostream& os, const PlainWeave& W, const std::string& command_line,
                                  const PlainWeaveFibExtra& extra = PlainWeaveFibExtra{})
{
    const auto& q = W.params;
    PlainWeaveSummary S = plain_weave_summary(W);
    if (extra.window_area > 0.0)
    {
        double Lx = 0.0, Ly = 0.0;
        plain_weave_window(q, Lx, Ly);
        S.areal_density = S.areal_density * (Lx * Ly) / extra.window_area;
    }
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "# Полотняное переплетение (plain weave), сгенерировано plain_weave.hpp\n";
    os << "# команда: " << command_line << "\n";
    os << "# d=" << q.d << " pitch=" << q.pitch << " amp=" << q.amp << " nx=" << q.nx << " ny=" << q.ny
       << " seg=" << q.seg << "\n";
    os << "# overhang_x=" << plain_weave_overhang_x(q) << " overhang_y=" << plain_weave_overhang_y(q)
       << " profile=" << (q.profile == PlainWeaveProfile::Cos ? "cos" : "arcs") << "\n";
    os << "# center=" << q.cx << "," << q.cy << "," << q.cz << " allow_penetration=" << (q.allow_penetration ? 1 : 0)
       << "\n";
    os << "# E=" << q.E << " nu=" << q.nu << " rho=" << q.rho << " eps_break=" << q.eps_break << " mu_ff=" << q.mu_ff
       << " mu_fs=" << q.mu_fs << "\n";
    const char* cl = q.clamp == PlainWeaveClamp::None ? "none"
                     : q.clamp == PlainWeaveClamp::X  ? "x"
                     : q.clamp == PlainWeaveClamp::Y  ? "y"
                                                      : "xy";
    os << "# clamp=" << cl;
    if (q.clamp != PlainWeaveClamp::None) os << " clamp_edge=" << q.clamp_edge;
    os << "\n";
    os << "# итого: нитей " << S.fibers << ", бусин " << S.beads << ", связей " << S.bonds << ", зажатых бусин "
       << S.clamped_beads << ", поверхностная плотность " << S.areal_density << " кг/м^2\n";
    for (const auto& c : extra.comments) os << "# " << c << "\n";
    for (const auto& w : W.warnings) os << "# " << w << "\n";

    os << "MATERIAL 0 " << q.E << " " << q.nu << " " << q.d << " " << q.rho << " " << q.eps_break;
    if (q.mu_fs != 0.0)
        os << " " << q.mu_ff << " " << q.mu_fs;
    else if (q.mu_ff != 0.0)
        os << " " << q.mu_ff;
    os << "\n";
    if (extra.pairs)
        os << "MATERIAL 1 " << extra.stitch_E << " " << extra.stitch_nu << " " << extra.stitch_d << " "
           << extra.stitch_rho << " " << extra.stitch_eps_break << "\n";

    for (const auto& F : W.fibers)
    {
        os << "FIBER 0 " << F.pts.size() << "\n";
        for (const auto& P : F.pts) os << P[0] << " " << P[1] << " " << P[2] << "\n";
    }
    if (extra.pairs)
    {
        os << "# стёжка: " << extra.pairs->size() << " стежков (материал 1)\n";
        for (const auto& pr : *extra.pairs) os << "BOND " << pr.first << " " << pr.second << " 1\n";
    }
    if (!W.clamp_boxes.empty())
    {
        os << "# закрепление: --clamp=" << cl << ", c_e = " << q.clamp_edge << ", ожидаемое число зажатых бусин "
           << S.clamped_beads << " (маска по умолчанию xyz)\n";
        for (const auto& b : W.clamp_boxes)
            os << "CLAMP_BOX " << b.lo[0] << " " << b.lo[1] << " " << b.lo[2] << " " << b.hi[0] << " " << b.hi[1]
               << " " << b.hi[2] << "\n";
    }
}

// Раскладка: строка на нить; first/last -- глобальные индексы бусин подряд по блокам.
inline void write_plain_weave_layout(std::ostream& os, const PlainWeave& W)
{
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "# object_id direction index offset clamped first_bead last_bead n_clamp_lo n_clamp_hi\n";
    std::size_t first = 0;
    for (const auto& F : W.fibers)
    {
        const std::size_t last = first + F.pts.size() - 1;
        os << F.object_id << " " << F.direction << " " << F.index << " " << F.offset << " " << (F.clamped ? 1 : 0)
           << " " << first << " " << last << " " << F.n_clamp_lo << " " << F.n_clamp_hi << "\n";
        first = last + 1;
    }
}

} // namespace fiber
