#pragma once
// bc6_geometry.hpp
// Построение геометрии кейса BC6 по правилу §2.3 постановки.
//
// Файл НЕ зависит от кейса BC2 и не переиспользует его построитель: правила
// раскладки у них разные. BC2 ставит узлы по центрам ячеек и подгоняет шаг
// под толщину слоя отдельно по каждой оси, чтобы толщина воспроизводилась
// точно. BC6 держит шаг решётки РАВНЫМ Δx по обеим осям и мирится с тем, что
// размах тела короче паспортной толщины не более чем на Δx: так требует
// постановка, потому что фиксированный зазор удара 1.000 мм должен быть
// одинаков при любом разрешении, а это возможно только при едином шаге.
//
// ПРАВИЛО РАСКЛАДКИ (§2.3):
//     r   = (i_r + 1/2)·Δx          смещение на полшага удерживает частицы
//     z   = z0 + i_z·Δx             строго внутри r > 0
//     n_z = floor(T/Δx) + 1
//     n_r = число i_r, для которых (i_r + 1/2)·Δx <= R
// Локальный индекс внутри тела: i_r·n_z + i_z, i_r — внешний цикл; сквозная
// нумерация с нуля, тела создаются в порядке 1 → 4.
//
// РАССТАНОВКА ТЕЛ ВДОЛЬ ОСИ. Тела в контакте (1↔2 и 3↔4) продолжают единую
// решётку: первый слой следующего тела на Δx дальше последнего слоя
// предыдущего. Между телами 2 и 3 — фиксированный зазор 1.000 мм.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <sstream>
#include <string>

#include "param.hpp"
#include "particle_system.hpp"

#include "bc6_case.hpp"
#include "e2e_harness.hpp"

namespace bc6
{

// Число слоёв частиц вдоль оси на протяжении T. Допуск 1e-9 снимает случай,
// когда T/Δx математически целое, но в double чуть меньше целого.
inline size_t layers_along_z(double T, double dx)
{
    if (!(T > 0.0) || !(dx > 0.0)) return 1;
    const long n = static_cast<long>(std::floor(T / dx + 1e-9)) + 1;
    return n < 1 ? 1u : static_cast<size_t>(n);
}

// Число столбцов частиц по радиусу: i_r от 0 до наибольшего, при котором
// (i_r + 1/2)·Δx <= R. Замкнутая форма — floor(R/Δx + 1/2).
inline size_t columns_along_r(double R, double dx)
{
    if (!(R > 0.0) || !(dx > 0.0)) return 1;
    const long n = static_cast<long>(std::floor(R / dx + 0.5 + 1e-9));
    return n < 1 ? 1u : static_cast<size_t>(n);
}

struct Geometry
{
    double delta_x = 0.0;
    size_t ntotal  = 0;
    std::array<size_t, 4> begin{}, count{}, n_r{}, n_z{};
    std::array<double, 4> z_first{}, z_last{};   // координаты крайних СЛОЁВ, м
    double impact_gap = 0.0;                     // z_first[2] − z_last[1];
                                                 // сборка встык => равен Δx

    size_t moni_axis   = 0;    // основной монитор: тыльный слой образца, i_r = 0
    size_t moni_ctl    = 0;    // контрольный: тот же слой, r ≈ 12.5 мм
    size_t moni_ctl_ir = 0;
    double moni_ctl_r  = 0.0;
    double moni_axis_r = 0.0;

    // Фактический размах тела вдоль оси (короче паспортной толщины не более
    // чем на Δx — см. §2.2 постановки).
    double extent(size_t b) const { return static_cast<double>(n_z[b] - 1) * delta_x; }
};

// Планирование без выделения частиц: нужно режиму --plot-from, где расчёта
// нет вовсе, а подпись графика обязана назвать разрешение и число частиц.
inline Geometry plan_geometry(double delta_x)
{
    Geometry g;
    g.delta_x = delta_x;

    double z0 = 0.0;
    size_t total = 0;
    for (size_t b = 0; b < kBodies.size(); ++b)
    {
        const BodySpec& s = kBodies[b];
        g.n_r[b] = columns_along_r(s.R, delta_x);
        g.n_z[b] = layers_along_z(s.T, delta_x);
        g.begin[b] = total;
        g.count[b] = g.n_r[b] * g.n_z[b];
        total += g.count[b];

        g.z_first[b] = z0;
        g.z_last[b]  = z0 + static_cast<double>(g.n_z[b] - 1) * delta_x;
        // Сборка встык: плоскость удара раскладывается так же, как
        // склеенные границы — решётка продолжается с шагом Δx. Поле
        // impact_after больше не меняет раскладку и осталось признаком
        // того, какая из границ является плоскостью удара.
        (void)s.impact_after;
        z0 = g.z_last[b] + delta_x;
    }
    g.ntotal = total;
    g.impact_gap = g.z_first[2] - g.z_last[1];

    // Мониторы — геометрически, не литералом индекса.
    const size_t nz = g.n_z[kSampleBody];
    const size_t nr = g.n_r[kSampleBody];
    g.moni_axis   = g.begin[kSampleBody] + (nz - 1);
    g.moni_axis_r = 0.5 * delta_x;

    // Контрольный столбец: наименьший i_r, у которого (i_r + 1/2)·Δx >= 12.5 мм.
    // Ограничение по n_r обязательно: в smoke-режиме у образца всего 48
    // столбцов, и литерал i_r = 62 из §2.4 там попросту не существует.
    long ir = static_cast<long>(std::ceil(kMoniCtlRadius / delta_x - 0.5 - 1e-9));
    if (ir < 0) ir = 0;
    if (ir > static_cast<long>(nr) - 1) ir = static_cast<long>(nr) - 1;
    g.moni_ctl_ir = static_cast<size_t>(ir);
    g.moni_ctl_r  = (static_cast<double>(g.moni_ctl_ir) + 0.5) * delta_x;
    g.moni_ctl    = g.begin[kSampleBody] + g.moni_ctl_ir * nz + (nz - 1);

    return g;
}

// Построение с заполнением частиц. Возвращает ту же схему, что plan_geometry.
inline Geometry build_geometry(ParticleSystem& sys, double delta_x)
{
    const Geometry g = plan_geometry(delta_x);

    // resize ДО первой записи: иначе запись идёт мимо выделенной памяти.
    sys.resize(g.ntotal);

    size_t k = 0;
    for (size_t b = 0; b < kBodies.size(); ++b)
    {
        const BodySpec& s = kBodies[b];
        // Масса — по ПЛОСКОЙ формуле ρ₀·Δx², а не по кольцевой 2πr·ρ₀·Δx²
        // (§3 постановки). Так же сделано в BC2; менять конвенцию в рамках
        // одного кейса нельзя, иначе кейсы станут несопоставимы.
        const double mass = s.rho0 * delta_x * delta_x;
        for (size_t ir = 0; ir < g.n_r[b]; ++ir)
        {
            const double r = (static_cast<double>(ir) + 0.5) * delta_x;
            for (size_t iz = 0; iz < g.n_z[b]; ++iz)
            {
                const double z = g.z_first[b] + static_cast<double>(iz) * delta_x;
                sys.axis[k]  = (ir == 0) ? 1 : 0;
                // h = 2Δx при построении; на первом же шаге h_upgrade при
                // sle = 1 заменит её на 1.3·(m/ρ)^(1/dim) ≈ 1.3·Δx.
                sys.hsml[k]  = 2.0 * delta_x;
                sys.x[0][k]  = r;
                sys.x[1][k]  = z;
                sys.vx[0][k] = 0.0;
                sys.vx[1][k] = s.vz;
                sys.itype[k] = s.itype;
                // body[] = itype[] (fiber/PLAN.md §3.1) -- см. то же
                // обоснование в tests/bc2_geometry.hpp: без этой строки
                // Solid-Solid-Contact (body[i]!=body[j]) никогда не срабатывает.
                sys.body[k]  = s.itype;
                sys.rho[k]   = s.rho0;
                sys.rho_0[k] = s.rho0;
                sys.mass[k]  = mass;
                // Сборка ненапряжена и покоится (§3). resize уже обнулил эти
                // поля, но задать их явно дешевле, чем полагаться на это.
                sys.p[k]    = 0.0;
                sys.u[k]    = 0.0;
                sys.D[k]    = 0.0;
                sys.e_pl[k] = 0.0;
                sys.sxx[k]  = 0.0; sys.syy[k] = 0.0; sys.szz[k] = 0.0;
                sys.sxy[k]  = 0.0; sys.syz[k] = 0.0; sys.sxz[k] = 0.0;
                ++k;
            }
        }
    }
    return g;
}

// ---------------------------------------------------------------------------
// Ожидаемое время прихода фронта на монитор (§5.3).
// Пробег берётся по ФАКТИЧЕСКОМУ размаху частиц образца, а не по паспортной
// толщине: при Δx = 0.2 мм это 10.4 мм вместо 10.526, и разница в 9 нс
// сопоставима с допуском ±0.05 мкс.
// kTImpact равен нулю (сборка встык). Слагаемое оставлено явно: оно
// показывает, из чего состоит ожидание, и вернёт смысл, если геометрию
// когда-нибудь снова разведут.
inline double t_front_expected(const Geometry& g)
{
    return kTImpact + g.extent(kSampleBody) / kClB4C;
}

// ---------------------------------------------------------------------------
inline void report_geometry(const Geometry& g)
{
    std::printf("\n[§2] Геометрия\n");
    for (size_t b = 0; b < kBodies.size(); ++b)
    {
        std::printf("  info : тело %zu, %-22s n_r=%4zu n_z=%4zu частиц=%6zu "
                    "z=%.4f…%.4f мм размах=%.4f мм\n",
                    b + 1, kBodies[b].name, g.n_r[b], g.n_z[b], g.count[b],
                    g.z_first[b] * 1e3, g.z_last[b] * 1e3, g.extent(b) * 1e3);
    }
    std::printf("  info : всего частиц %zu; Δx = %.4f мм; шаг на плоскости удара %.6f мм\n",
                g.ntotal, g.delta_x * 1e3, g.impact_gap * 1e3);
    std::printf("  info : монитор основной  #%zu (i_r = 0,  r = %.3f мм)\n",
                g.moni_axis, g.moni_axis_r * 1e3);
    std::printf("  info : монитор контрольный #%zu (i_r = %zu, r = %.3f мм)\n",
                g.moni_ctl, g.moni_ctl_ir, g.moni_ctl_r * 1e3);
    std::printf("  info : удар при t = %.4f мкс (сборка встык); фронт ожидается при %.4f мкс\n",
                kTImpact * 1e6, t_front_expected(g) * 1e6);
}

// Проверки §6, зависящие только от геометрии. Работают в обоих режимах.
//
// Числа частиц ВЫЧИСЛЯЮТСЯ по формулам §2.3 (см. layers_along_z /
// columns_along_r), а литералы ниже только сверяются с вычисленным — как и
// требует §2.4 постановки. Это замки от регрессии, а не источник истины.
inline void assert_control_numbers(const Geometry& g, bool full)
{
    using e2e::check;
    using e2e::note;

    std::printf("\n[§6] Контрольные числа геометрии\n");

    // --- не зависит от режима ---
    check(std::fabs(g.impact_gap - g.delta_x) < 1e-12,
          "сборка встык: слои на плоскости удара расходятся ровно на Δx");
    check(std::fabs((g.z_first[1] - g.z_last[0]) - g.delta_x) < 1e-12,
          "тела 1 и 2 в контакте: слои расходятся ровно на Δx");
    check(std::fabs((g.z_first[3] - g.z_last[2]) - g.delta_x) < 1e-12,
          "тела 3 и 4 в контакте: слои расходятся ровно на Δx");

    // Ряд радиусов невозрастающий, но НЕ строго убывающий (§2.5): тела 1-2
    // имеют одинаковый радиус, и тела 3-4 тоже. Построитель обязан это
    // выдерживать.
    bool non_increasing = true, strictly_decreasing = true;
    for (size_t b = 1; b < kBodies.size(); ++b)
    {
        if (kBodies[b].R > kBodies[b - 1].R + 1e-15) non_increasing = false;
        if (!(kBodies[b].R < kBodies[b - 1].R - 1e-15)) strictly_decreasing = false;
    }
    check(non_increasing, "ряд радиусов невозрастающий");
    check(!strictly_decreasing,
          "ряд радиусов НЕ строго убывающий (43.75 → 43.75 → 38.05 → 38.05) — "
          "построитель не предполагает строгого убывания");

    // Мониторы.
    const size_t lo = g.begin[kSampleBody];
    const size_t hi = lo + g.count[kSampleBody];
    const size_t nz = g.n_z[kSampleBody];
    check(g.moni_axis >= lo && g.moni_axis < hi,
          "основной монитор лежит внутри образца (тело 3)");
    check(g.moni_ctl >= lo && g.moni_ctl < hi,
          "контрольный монитор лежит внутри образца (тело 3)");
    check((g.moni_axis - lo) % nz == nz - 1,
          "основной монитор на ТЫЛЬНОМ слое образца");
    check((g.moni_ctl - lo) % nz == nz - 1,
          "контрольный монитор на ТЫЛЬНОМ слое образца");
    check(g.moni_axis == lo + (nz - 1),
          "основной монитор в приосевом столбце (i_r = 0)");

    // --- по режимам: замки от регрессии ---
    const std::array<size_t, 4> expect = full
        ? std::array<size_t, 4>{8979, 2409, 10070, 23180}
        : std::array<size_t, 4>{605, 165, 672, 1488};
    const size_t expect_total = full ? 44638u : 2930u;

    for (size_t b = 0; b < 4; ++b)
    {
        std::ostringstream os;
        os << "тело " << (b + 1) << ": частиц " << g.count[b]
           << " (ожидалось " << expect[b] << ")";
        check(g.count[b] == expect[b], os.str());
    }
    {
        std::ostringstream os;
        os << "всего частиц " << g.ntotal << " (ожидалось " << expect_total << ")";
        check(g.ntotal == expect_total, os.str());
    }

    if (full)
    {
        check(g.moni_axis == 11440, "сквозной индекс основного монитора = 11440 (§2.4)");
        check(g.moni_ctl == 14726, "сквозной индекс контрольного монитора = 14726 (§2.4)");
        check(g.moni_ctl_ir == 62, "контрольный столбец i_r = 62 (§2.4)");
        check(std::fabs(g.moni_ctl_r - 12.5e-3) < 1e-9,
              "контрольный монитор ровно на r = 12.5 мм");
        // Фактические размахи (§2.2): 8.0 / 2.0 / 10.4 / 24.2 мм.
        const std::array<double, 4> ext{8.0e-3, 2.0e-3, 10.4e-3, 24.2e-3};
        for (size_t b = 0; b < 4; ++b)
        {
            std::ostringstream os;
            os << "тело " << (b + 1) << ": размах " << (g.extent(b) * 1e3)
               << " мм (ожидалось " << (ext[b] * 1e3) << ")";
            check(std::fabs(g.extent(b) - ext[b]) < 1e-12, os.str());
        }
    }
    else
    {
        // В smoke контрольный столбец до 12.5 мм не дотягивает: у образца
        // всего 48 столбцов, ближайший центр — 12.4 мм. Это не дефект, а
        // следствие огрубления сетки, и проверять надо именно ограничение.
        std::ostringstream os;
        os << "контрольный монитор смещён на r = " << (g.moni_ctl_r * 1e3)
           << " мм (ближайший столбец при Δx = " << (g.delta_x * 1e3) << " мм)";
        note(os.str());
        check(g.moni_ctl_ir < g.n_r[kSampleBody],
              "контрольный столбец ограничен числом столбцов образца");
    }
}

} // namespace bc6
