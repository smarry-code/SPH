#pragma once
// jc_taylor.hpp
// Тест Тейлора (E4) и сходимость по разрешению (E5) для модели Джонсона-Кука:
// цилиндр STEEL 4340 (itype 12) Ø 7.62 мм, L0 = 25.4 мм, 200 м/с в жёсткую
// стенку. План: fiber/jc-implementation-plan.md §5.3 E4/E5.
//
// ВСЕ числа постановки и построение геометрии — здесь; прогон, проверки и
// режимы командной строки — в tests/test_e2e_jc_taylor.cpp.
//
// ---------------------------------------------------------------------------
// КАК ЗАДАНА СТЕНКА
// ---------------------------------------------------------------------------
// Жёсткой стенки в решателе нет. Взята точная симметричная замена: ДВА
// одинаковых цилиндра летят навстречу друг другу со скоростями ±V0 вдоль z и
// касаются в плоскости z = 0. По симметрии z -> -z в этой плоскости v_z = 0 и
// нет касательного напряжения: ровно то, что даёт гладкая жёсткая стенка.
// Анализируется половина z > 0 (первая в нумерации частиц).
//
// Это ОДИН сплошной материал, а не два тела с контактом:
//   * контактные силы (contact_force_type 1..3) работают на парах разных
//     itype (ext_force.cpp, is_contact_material_pair), два одинаковых тела из
//     одного материала контакта не видят. SPH-пары при этом не фильтруются
//     по телам (sph_same_body_only = false), и градиент давления через
//     плоскость z = 0 передаёт нагрузку так же, как внутри тела;
//   * при сжатии (всё время нагружения) взаимопроникновения нет: давление в
//     слое у плоскости положительно и растёт, поэтому частицы половин не
//     переходят через z = 0 (проверяется в тесте: min z половины > 0);
//   * единственное отличие от настоящей стенки — растяжение. Стенка не держит
//     растяжения, а сплошная среда держит. Это сказывается лишь в конце
//     (отскок торца), когда нормальное напряжение на торце уже около нуля;
//     тест выводит минимальное давление в приторцевом слое, чтобы это было
//     видно.
//
// Рассмотренные и отвергнутые альтернативы:
//   * Clamp (include/interactions/clamp.hpp, слой полностью закреплённых
//     частиц как стенка). Закреплённые частицы остаются обычными частицами
//     SPH: у них копится девиатор от радиального скольжения торца, то есть
//     получается стенка с трением (до предела текучести стенки) и к тому же
//     сжимаемая. Тейлор-эксперимент идёт о гладкую стенку; трение меняет
//     грибовидный торец заметно сильнее, чем неточность стенки из
//     симметрии. Выигрыш Clamp — вдвое меньше частиц, но не оправдывает
//     систематической ошибки.
//   * Зеркальные призраки у z = 0 (по аналогии с призраками у оси) — требует
//     правки src/, запрещено.
//
// ---------------------------------------------------------------------------
// ОСЕСИММЕТРИЯ (SPH_DIM=2, axis_sim): x[0] = r, x[1] = z; sxx = S_rr,
// syy = S_zz, szz = S_theta; масса — плоская rho0*dx^2 (как BC6, §3).
// Частицы: r = (ir + 1/2)*dx, z = ±(iz + 1/2)*dx. Решётка симметрична
// относительно z = 0 (первый слой каждой половины на dx/2 от плоскости).
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>

#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"

namespace jct
{

// ===========================================================================
// Исходные данные опыта
// ===========================================================================
inline constexpr int    kItype = 12;          // STEEL 4340 (make_jc_library)
inline constexpr double kRho0  = 7830.0;      // кг/м³
inline constexpr double kR0    = 3.81e-3;     // радиус, м (Ø 7.62 мм)
inline constexpr double kL0    = 25.4e-3;     // длина одной половины, м
inline constexpr double kV0    = 200.0;       // скорость удара, м/с

// Число ячеек по радиусу. n_z = L0/dx получается целым при n_r, кратном 3
// (L0/R0 = 20/3): 6 -> 40 слоёв, 12 -> 80, 24 -> 160, 48 -> 320.
// Замеры (2 потока OpenMP, g++ -O2, ноутбук): 6 — 1 с, 12 — 3 с, 24 — 30 с
// (7680 частиц, ~4300 шагов, 0.9 мкс на частицу-шаг); 48 — оценка 5-8 мин.
inline constexpr size_t kNrE4     = 24;       // E4: dx = 0.159 мм, 160 слоёв на половину
inline constexpr size_t kNrConv   = 12;       // E5: грубый уровень h (тонкий h/2 = 24 = уровень E4)
inline constexpr size_t kNrSmoke  = 6;        // smoke: dx = 0.635 мм, 40 слоёв

// Время и критерий остановки.
inline constexpr double kTMax        = 150.0e-6;  // потолок времени расчёта, с
inline constexpr double kStopKeFrac  = 0.01;      // стоп, когда KE половины < 1 % от начальной
inline constexpr size_t kProbeEvery  = 25;        // период опроса наблюдателем, шагов
inline constexpr size_t kSmokeSteps  = 300;       // длина smoke-прогона, шагов

// Допуски приёмки (план E4/E5).
inline constexpr double kLfLo       = 0.60;       // L_final/L0 в (0.60, 0.95)
inline constexpr double kLfHi       = 0.95;
inline constexpr double kRadiusTol  = 0.01;       // допуск убывания радиуса торца, доля R0
inline constexpr double kEplFaceL   = 0.05;       // максимум e_pl не дальше kEplFaceL*L0 (1.27 мм) от z = 0:
                                                  // «несколько h» при h ~ 1.3 dx; физическое расстояние, не от dx
                                                  // (замер: z_max = 0.34 / 0.69 мм при n_r = 12 / 24)
inline constexpr double kEnergyTol  = 0.02;       // баланс энергии, доля E_mech0
inline constexpr double kAutodynTol = 0.05;       // сверка с Autodyn: L_final, R_max
inline constexpr double kConvTol    = 0.03;       // E5: |L(h) - L(h/2)|/L(h/2)

// ===========================================================================
// Геометрия
// ===========================================================================
struct Geometry
{
    double dx = 0.0;
    size_t n_r = 0, n_z = 0;      // на ПОЛОВИНУ
    size_t nhalf = 0;             // частиц в половине z > 0: индексы [0, nhalf)
    size_t ntotal = 0;            // обе половины

    // Сквозной индекс частицы половины z > 0.
    size_t index(size_t ir, size_t iz) const { return ir * n_z + iz; }
    double L0() const { return static_cast<double>(n_z) * dx; }
    double R0() const { return static_cast<double>(n_r) * dx; }
};

// Планирование без выделения частиц. dx = R0/n_r (радиус точный); n_z — ближайшее
// целое к L0/dx. ok = false, если длина половины при этом отличается от L0
// больше чем на 0,1 % (для n_r кратного 3 отличие нулевое).
inline Geometry plan_geometry(size_t n_r, bool* ok = nullptr)
{
    Geometry g;
    g.n_r = n_r < 1 ? 1 : n_r;
    g.dx = kR0 / static_cast<double>(g.n_r);
    g.n_z = static_cast<size_t>(std::floor(kL0 / g.dx + 0.5));
    if (g.n_z < 1) g.n_z = 1;
    g.nhalf = g.n_r * g.n_z;
    g.ntotal = 2 * g.nhalf;
    if (ok) *ok = std::fabs(g.L0() - kL0) <= 1.0e-3 * kL0;
    return g;
}

// Размер решётки по заданному dx (ближайшее целое число ячеек по радиусу).
inline size_t nr_from_dx(double dx)
{
    const long n = static_cast<long>(std::floor(kR0 / dx + 0.5));
    return n < 1 ? 1u : static_cast<size_t>(n);
}

inline Geometry build_geometry(ParticleSystem& sys, size_t n_r)
{
    const Geometry g = plan_geometry(n_r);
    sys.resize(g.ntotal);
    const double mass = kRho0 * g.dx * g.dx;
    for (size_t half = 0; half < 2; ++half)
    {
        const double sgn = (half == 0) ? 1.0 : -1.0;   // половина z > 0 летит к стенке (-z)
        for (size_t ir = 0; ir < g.n_r; ++ir)
        {
            const double r = (static_cast<double>(ir) + 0.5) * g.dx;
            for (size_t iz = 0; iz < g.n_z; ++iz)
            {
                const size_t k = half * g.nhalf + g.index(ir, iz);
                const double z = sgn * (static_cast<double>(iz) + 0.5) * g.dx;
                sys.axis[k] = (ir == 0) ? 1 : 0;
                sys.hsml[k] = 2.0 * g.dx;       // h_upgrade (sle = 1) на первом шаге пересчитает
                sys.x[0][k] = r;
                sys.x[1][k] = z;
                sys.vx[0][k] = 0.0;
                sys.vx[1][k] = -sgn * kV0;
                sys.itype[k] = kItype;
                sys.body[k]  = kItype;
                sys.rho[k]   = kRho0;
                sys.rho_0[k] = kRho0;
                sys.mass[k]  = mass;
                sys.p[k] = 0.0; sys.u[k] = 0.0; sys.D[k] = 0.0; sys.e_pl[k] = 0.0;
                sys.sxx[k] = sys.syy[k] = sys.szz[k] = 0.0;
                sys.sxy[k] = sys.syz[k] = sys.sxz[k] = 0.0;
            }
        }
    }
    return g;
}

// ===========================================================================
// Конфигурация решателя (параметры BC6, tests/data/bc6_initial_conditions.txt:
// pa_sph = 2, nnps = 2, sle = 1, skf = 1, CFL = 0.1, XSPH 0.25, alpha = beta = 1)
// ===========================================================================
inline void apply_config(SimulationConfig& cfg, double t_total)
{
    cfg.apply_defaults();
    cfg.pa_sph = 2;
    cfg.nnps = 2;
    cfg.sle = 1;
    cfg.skf = 1;
    cfg.CFL = 0.1;
    cfg.epsilon = 0.25;
    cfg.alpha = 1.0;
    cfg.beta = 1.0;
    cfg.etq = 0.01;
    cfg.contact_force_type = 0;
    cfg.t_total = t_total;
    // Кадры конвейера не нужны (свои кадры пишет наблюдатель): один кадр t = 0.
    cfg.save_t = 1.0;
    cfg.print_step = 100000000;
}

// ===========================================================================
// Измерения по половине z > 0
// ===========================================================================
struct Measure
{
    double ke = 0.0;            // кинетическая энергия половины, Дж (по плоской массе)
    double L = 0.0;             // длина по осевому столбцу (ir <= 1)
    double z_face = 0.0;        // положение торца на оси
    double R_face = 0.0;        // радиус торца: max r в приторцевом слое + dx/2
    double e_pl_axis = 0.0;     // e_pl осевой частицы у торца (ir = 0, iz = 0)
    double e_pl_max = 0.0;      // максимум e_pl по половине
    double e_pl_max_r = 0.0, e_pl_max_z = 0.0;
    double T_max = 0.0;
    double p_min_face = 0.0;    // минимум давления в слое у торца
    double z_min = 0.0;         // наименьший z частицы половины (проникание через плоскость < 0)
    double vz_com = 0.0;        // масс-взвешенная осевая скорость половины
    // Полная энергия половины sum m (|v|^2/2 + u) с разными весами (диагностика баланса):
    // e_flat -- без веса (как energy_log), e_r0 -- вес r0/R0 (масса кольца ~ 2 pi r0 m),
    // e_r -- вес r/R0 по текущему радиусу.
    double e_flat = 0.0, e_r0 = 0.0, e_r = 0.0;
};

inline Measure measure_half(const ParticleSystem& sys, const Geometry& g)
{
    Measure m;
    double z_lo = 1e30;
    for (size_t k = 0; k < g.nhalf; ++k) z_lo = std::min(z_lo, sys.x[1][k]);
    m.z_min = z_lo;

    double zb = -1e30, zf = 1e30, mass = 0.0, mvz = 0.0;
    m.p_min_face = 1e30;
    m.R_face = 0.0;
    m.e_pl_max = -1.0;
    for (size_t k = 0; k < g.nhalf; ++k)
    {
        const double vr = sys.vx[0][k], vz = sys.vx[1][k];
        m.ke += 0.5 * sys.mass[k] * (vr * vr + vz * vz);
        mass += sys.mass[k];
        mvz += sys.mass[k] * vz;
        const size_t ir = k / g.n_z;
        const double e_k = sys.mass[k] * (0.5 * (vr * vr + vz * vz) + sys.u[k]);
        m.e_flat += e_k;
        m.e_r0 += e_k * (static_cast<double>(ir) + 0.5) * g.dx / kR0;
        m.e_r += e_k * sys.x[0][k] / kR0;
        if (ir <= 1)
        {
            zb = std::max(zb, sys.x[1][k]);
            zf = std::min(zf, sys.x[1][k]);
        }
        if (sys.x[1][k] < z_lo + 1.01 * g.dx)
        {
            m.R_face = std::max(m.R_face, sys.x[0][k]);
            m.p_min_face = std::min(m.p_min_face, sys.p[k]);
        }
        if (sys.e_pl[k] > m.e_pl_max)
        {
            m.e_pl_max = sys.e_pl[k];
            m.e_pl_max_r = sys.x[0][k];
            m.e_pl_max_z = sys.x[1][k];
        }
        m.T_max = std::max(m.T_max, sys.T[k]);
    }
    m.R_face += 0.5 * g.dx;
    m.z_face = zf;
    m.L = zb - zf + g.dx;
    m.e_pl_axis = sys.e_pl[g.index(0, 0)];
    m.vz_com = mass > 0.0 ? mvz / mass : 0.0;
    return m;
}

} // namespace jct
