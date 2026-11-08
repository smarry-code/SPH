// test_axisym.cpp
// Unit-тест осесимметричной ветки силового расчёта (axis_sim && dim == 2).
//
// ЗАЧЕМ. В этой ветке нет однозначной записи о том, какая переменная какую
// компоненту тензора хранит, и это уже приводило к ошибке. Кинематика
// заполняется по индексам координат (0 -> exx, 1 -> eyy), поэтому фактическое
// соглашение таково:
//
//     sxx = S_rr      syy = S_zz      szz = S_theta      sxy = S_rz
//     exx = eps_rr    eyy = eps_zz    ezz = eps_theta = v_r / r
//
// То есть ОСЕВАЯ компонента лежит в syy, а szz — ОКРУЖНАЯ. Часть формул
// осесимметрии была записана в противоположном соглашении (szz как осевая,
// окружная не хранится и равна -(sxx + szz)), и расхождение проявилось сразу
// в четырёх местах: радиальный источник уравнения движения, осевой источник,
// невычисляемая ezz и девиаторная работа в энергии.
//
// Тест закрепляет соглашение проверками, которые следуют из симметрии задачи
// и не зависят от точности SPH-аппроксимации:
//
//   A1  одноосное сжатие вдоль оси: dsxx == dszz (S_rr = S_theta по симметрии)
//   A2  то же: dsxx + dsyy + dszz == 0 (девиатор остаётся бесследовым)
//   A3  то же: dsyy == -2*dsxx (осевая компонента — именно syy)
//   B1  радиальное расширение: dsxx == dszz
//   B2  то же: сумма == 0                     <- ловит ezz = 0
//   C   однородное объёмное сжатие: все три производные нулевые  <- ловит ezz = 0
//   D   покой с однородным девиатором одноосного сжатия: dv_r/dt == 0
//                                               <- ловит (2*sxx + szz) вместо (sxx - szz)
//   E   покой с однородным сдвигом: dv_z/dt == S_rz/(rho r)
//                                               <- ловит sxz (нулевую при dim = 2) вместо sxy
//
// Проверки A1/B1/C/D — точные тождества схемы, поэтому допуск задан
// относительно масштаба самих производных, а не абсолютным числом.
// fiber/PLAN.md, этап Ф4-B2, шаг 1: probe_rates() ниже больше НЕ вызывает
// старый путь (direct_find/con_density/int_force) напрямую — конвейер
// заменён на тот же набор вызовов, что src/single_step.cpp использует для
// продового расчёта (Solid derived -> GhostSet -> SolidExtended -> поиск ->
// Solid-SPH kernel -> reduce -> finalize), БЕЗ стадии Solid-Solid-Contact:
// тест проверяет только соглашение об индексах девиатора в осесимметричной
// ветке finalize_solid, а не контакт. Смысл и допуски проверок A1-E2 не
// изменились — сравниваются те же величины (dsxx/dsyy/dszz/dvr/dvz), теперь
// читаемые из interactions::SolidFinal вместо выходных массивов int_force.
#include "interactions/solid_sph.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <span>
#include <string>
#include <utility>
#include <vector>

static_assert(dim == 2, "тест осесимметрии требует сборки с SPH_DIM=2");
static_assert(axis_sim, "тест осесимметрии требует axis_sim == true");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

// Проверка «a == b» с допуском, отнесённым к масштабу сравниваемых величин.
// Абсолютный порог здесь не годится: производные девиатора имеют порядок
// 2*G*eps ~ 1e14 Па/с, и любое фиксированное число либо всё пропустит, либо
// всё завалит.
bool near_rel(double a, double b, double scale, double tol)
{
    return std::fabs(a - b) <= tol * std::max(std::fabs(scale), 1.0);
}

// --- Расчётная площадка -----------------------------------------------------
// Регулярная решётка в меридиональной полуплоскости. Пробная частица берётся
// в глубине: до любой границы не меньше четырёх слоёв, то есть её носитель
// (2h = 2.6*dx) заполнен целиком и краевые эффекты в неё не достают.
constexpr size_t kNr = 14;
constexpr size_t kNz = 14;
constexpr double kDx   = 1.0e-4;     // 0.1 мм
constexpr double kRho0 = 2506.0;     // B4C
constexpr int    kItype = 10;        // B4C, JH-2 (модуль сдвига 199.4 ГПа)
constexpr size_t kProbeIr = 7;
constexpr size_t kProbeIz = 7;

size_t flat_index(size_t ir, size_t iz) { return ir * kNz + iz; }

struct Rates
{
    double dsxx = 0.0, dsyy = 0.0, dszz = 0.0;
    double dvr  = 0.0, dvz  = 0.0;
    double scale = 0.0;              // масштаб задачи для относительных допусков
};

using VelocityField = std::function<void(double r, double z, double& vr, double& vz)>;

// Однократный проход конвейера: зеркальные частицы -> соседи -> неразрывность
// -> внутренние силы. Возвращает производные на пробной частице.
//
// Начальный девиатор задаётся снаружи (deviator), потому что проверка D
// требует ненулевого S при нулевых скоростях: источник (S_rr - S_theta)/r
// иначе тождественно нулевой и ошибку в нём не увидеть.
Rates probe_rates(const VelocityField& field,
                  double s_rr = 0.0, double s_zz = 0.0, double s_theta = 0.0,
                  double s_rz = 0.0)
{
    const size_t n = kNr * kNz;

    ParticleSystem sys;
    // Новый путь не дописывает зеркальные частицы в sys (они живут только в
    // interactions::SolidExtended, построенном из GhostSet) — в отличие от
    // старого virt_part.cpp, запас под них здесь не нужен.
    sys.resize(n);

    for (size_t ir = 0; ir < kNr; ++ir)
    {
        for (size_t iz = 0; iz < kNz; ++iz)
        {
            const size_t k = flat_index(ir, iz);
            const double r = (static_cast<double>(ir) + 0.5) * kDx;
            const double z = static_cast<double>(iz) * kDx;

            double vr = 0.0, vz = 0.0;
            field(r, z, vr, vz);

            sys.x[0][k]  = r;
            sys.x[1][k]  = z;
            sys.vx[0][k] = vr;
            sys.vx[1][k] = vz;
            sys.itype[k] = kItype;
            sys.rho[k]   = kRho0;
            sys.rho_0[k] = kRho0;     // rho == rho_0 => mu = 0 => p = 0
            sys.mass[k]  = kRho0 * kDx * kDx;
            sys.hsml[k]  = 1.3 * kDx;
            sys.axis[k]  = (ir == 0) ? 1 : 0;
            sys.p[k]     = 0.0;
            sys.u[k]     = 0.0;
            sys.D[k]     = 0.0;
            sys.e_pl[k]  = 0.0;
            sys.sxx[k]   = s_rr;      // соглашение: sxx = S_rr
            sys.syy[k]   = s_zz;      //             syy = S_zz
            sys.szz[k]   = s_theta;   //             szz = S_theta
            sys.sxy[k]   = s_rz;      //             sxy = S_rz
            sys.syz[k]   = 0.0;
            sys.sxz[k]   = 0.0;
        }
    }

    SimulationConfig cfg;
    cfg.pa_sph          = 2;
    cfg.skf             = 1;          // кубический сплайн
    cfg.density_method  = 1;
    cfg.maxn            = n;
    cfg.max_interaction = n * 200;     // не используется новым путём (нет max_interaction-буферов), оставлено для минимальности diff

    // Тот же набор вызовов, что src/single_step.cpp (новый путь, fiber/PLAN.md
    // Ф4-B1), без стадии Solid-Solid-Contact (тест не про контакт).
    rhs::CpuExecutor exec;

    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec, cfg.mats, sys, n, derived);

    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, n, ghosts);

    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, n, derived, ghosts, ext);

    interactions::SolidSphSearchResult search;
    interactions::build_solid_sph_search(exec, cfg, ext, search);

    rhs::IndexedContributions<interactions::SolidDelta> sph_buffers;
    interactions::compute_solid_sph_contributions(exec, cfg, ext, search.primitives, sph_buffers);

    const std::vector<interactions::SolidDelta> reduced = rhs::reduce_contributions(exec, sph_buffers, n);

    std::vector<interactions::SolidFinal> finalres;
    interactions::finalize_solid(exec, cfg.mats, sys, n, reduced, derived.p, finalres);

    const size_t probe = flat_index(kProbeIr, kProbeIz);
    const interactions::SolidFinal& f = finalres[probe];
    Rates out;
    out.dsxx = f.dsxx;
    out.dsyy = f.dsyy;
    out.dszz = f.dszz;
    out.dvr  = f.dv[0];
    out.dvz  = f.dv[1];
    out.scale = std::max({std::fabs(out.dsxx), std::fabs(out.dsyy), std::fabs(out.dszz)});
    return out;
}

// --- A. Одноосное сжатие вдоль оси -----------------------------------------
// v_r = 0, v_z = -k*z. Тогда eps_rr = eps_theta = 0, eps_zz = -k.
// Девиатор обязан расти как S_rr = S_theta = +2Gk/3, S_zz = -4Gk/3.
void test_uniaxial_compression()
{
    std::printf("\n[A] Одноосное сжатие вдоль оси (v_r = 0, v_z = -k z)\n");
    constexpr double k = 1.0e5;       // 1/с
    const Rates r = probe_rates([](double, double z, double& vr, double& vz) {
        vr = 0.0;
        vz = -k * z;
    });

    std::printf("        dsxx = %.6e  dsyy = %.6e  dszz = %.6e  (Па/с)\n",
                r.dsxx, r.dsyy, r.dszz);

    check(near_rel(r.dsxx, r.dszz, r.scale, 1e-12),
          "A1 dsxx == dszz: S_rr и S_theta растут одинаково (осевая симметрия)");
    check(near_rel(r.dsxx + r.dsyy + r.dszz, 0.0, r.scale, 1e-12),
          "A2 dsxx + dsyy + dszz == 0: девиатор остаётся бесследовым");
    check(near_rel(r.dsyy, -2.0 * r.dsxx, r.scale, 1e-12),
          "A3 dsyy == -2*dsxx: осевая компонента девиатора лежит в syy");
    check(r.dsyy < 0.0, "A4 dsyy < 0: сжатие вдоль оси даёт отрицательную S_zz");
}

// --- B. Радиальное расширение ----------------------------------------------
// v_r = a*r, v_z = 0. Тогда eps_rr = eps_theta = a, eps_zz = 0 — и девиатор
// снова обязан быть симметричным по (rr, theta) и бесследовым.
//
// Это ключевой случай для ezz: окружная скорость деформации здесь РАВНА a,
// а не нулю, и если её не вычислять, dszz расходится с dsxx ровно на 2*G*a.
void test_radial_expansion()
{
    std::printf("\n[B] Радиальное расширение (v_r = a r, v_z = 0)\n");
    constexpr double a = 1.0e5;       // 1/с
    const Rates r = probe_rates([](double rr, double, double& vr, double& vz) {
        vr = a * rr;
        vz = 0.0;
    });

    std::printf("        dsxx = %.6e  dsyy = %.6e  dszz = %.6e  (Па/с)\n",
                r.dsxx, r.dsyy, r.dszz);

    // Точного равенства здесь быть не может, и это свойство схемы, а не
    // дефект: eps_rr собирается SPH-суммой по парам, а eps_theta = v_r/r
    // берётся аналитически. Измеренное расхождение на этой решётке 3.1 %.
    // Допуск 10 % оставляет тридцатикратный запас до ошибки, которую проверка
    // ловит: при ezz = 0 величина dszz равна -4Ga/3, то есть расходится
    // с dsxx почти на 300 % и меняет знак.
    check(near_rel(r.dsxx, r.dszz, r.scale, 0.10),
          "B1 dsxx ~= dszz: eps_theta = v_r/r равна eps_rr, значит и девиаторы равны");
    check(near_rel(r.dsxx + r.dsyy + r.dszz, 0.0, r.scale, 1e-6),
          "B2 dsxx + dsyy + dszz == 0: след не дрейфует при v_r != 0");
}

// --- C. Однородное объёмное сжатие -----------------------------------------
// v = -c*(r, z): все три главные скорости деформации равны -c, форма не
// меняется, девиатор обязан остаться нулевым.
void test_isotropic_compression()
{
    std::printf("\n[C] Однородное объёмное сжатие (v = -c (r, z))\n");
    constexpr double c = 1.0e5;       // 1/с
    const Rates r = probe_rates([](double rr, double z, double& vr, double& vz) {
        vr = -c * rr;
        vz = -c * z;
    });

    // Масштаб для допуска — не сами производные (они должны быть нулевыми), а
    // характерная величина 2*G*c, с которой они сравниваются.
    constexpr double kShearB4C = 199.4e9;
    const double scale = 2.0 * kShearB4C * c;
    std::printf("        dsxx = %.6e  dsyy = %.6e  dszz = %.6e  (масштаб 2Gc = %.3e)\n",
                r.dsxx, r.dsyy, r.dszz, scale);

    // Остаточные 0.4-0.7 % от 2Gc — та же дискретизационная погрешность, что
    // и в B: численные eps_rr, eps_zz против аналитической eps_theta. При
    // ezz = 0 величина dszz равна 2Gc, то есть все 100 % масштаба.
    check(near_rel(r.dsxx, 0.0, scale, 0.05), "C1 dsxx ~= 0: форма не меняется");
    check(near_rel(r.dsyy, 0.0, scale, 0.05), "C2 dsyy ~= 0");
    check(near_rel(r.dszz, 0.0, scale, 0.05), "C3 dszz ~= 0: окружная компонента тоже");
}

// --- D. Осевой источник уравнения движения ---------------------------------
// Покой, однородный бесследовый девиатор состояния одноосного сжатия:
// S_rr = S_theta = S, S_zz = -2S. Источник (S_rr - S_theta)/r при таком
// состоянии тождественно нулевой, значит радиального ускорения быть не должно.
// Запись (2*S_rr + S_theta)/r дала бы 3S/r — при S = 1 ГПа и r = 0.75 мм это
// 4e12 Па/м, то есть ускорение порядка 1.6e9 м/с².
void test_axial_source_at_rest()
{
    std::printf("\n[D] Покой + однородный девиатор одноосного сжатия\n");
    constexpr double S = 1.0e9;       // Па
    const Rates r = probe_rates([](double, double, double& vr, double& vz) {
        vr = 0.0;
        vz = 0.0;
    }, /*s_rr=*/S, /*s_zz=*/-2.0 * S, /*s_theta=*/S);

    // Масштаб — ускорение, которое дал бы ошибочный источник 3S/(rho*r).
    const double r_probe = (static_cast<double>(kProbeIr) + 0.5) * kDx;
    const double scale = 3.0 * S / (kRho0 * r_probe);
    std::printf("        dv_r/dt = %.6e  dv_z/dt = %.6e  (масштаб ошибки 3S/(rho r) = %.3e м/с²)\n",
                r.dvr, r.dvz, scale);

    check(near_rel(r.dvr, 0.0, scale, 1e-3),
          "D1 dv_r/dt == 0: источник (S_rr - S_theta)/r на этом состоянии нулевой");
    check(near_rel(r.dvz, 0.0, scale, 1e-3),
          "D2 dv_z/dt == 0: однородное поле напряжений не создаёт осевой силы");
}

// --- E. Осевой источник от сдвига ------------------------------------------
// Покой, однородный сдвиг S_rz = T, остальные компоненты нулевые. Осевое
// уравнение движения в цилиндрических координатах содержит источник
// sigma_rz / r, поэтому радиальная сила отсутствует, а осевое ускорение
// обязано быть T/(rho*r) — ненулевым.
//
// Проверка нужна потому, что этот член брал компоненту sxz, которой при
// dim == 2 просто нет: она тождественно нулевая, и весь источник молча
// пропадал. Сдвиг S_rz в плоскости меридиана хранится в sxy.
void test_shear_source_at_rest()
{
    std::printf("\n[E] Покой + однородный сдвиг S_rz\n");
    constexpr double T = 1.0e9;       // Па
    const Rates r = probe_rates([](double, double, double& vr, double& vz) {
        vr = 0.0;
        vz = 0.0;
    }, /*s_rr=*/0.0, /*s_zz=*/0.0, /*s_theta=*/0.0, /*s_rz=*/T);

    const double r_probe = (static_cast<double>(kProbeIr) + 0.5) * kDx;
    const double expected = T / (kRho0 * r_probe);
    std::printf("        dv_z/dt = %.6e, ожидалось T/(rho r) = %.6e м/с² (%.1f %%)\n",
                r.dvz, expected, 100.0 * r.dvz / expected);

    // Допуск широкий намеренно. Источник собирается как сумма ПО ПАРАМ с
    // весом w[k], а самопара (j == i) в списке соседей отсутствует, поэтому
    // сумма недобирает самовклад ядра m/rho * W(0) = dx^2 * 15/(7 pi h^2) *
    // 2/3, что при h = 1.3 dx даёт 0.269. Ожидаемая доля 0.731, измеренная
    // 0.741 — разница остаётся на дискретизацию. Это свойство реализации
    // (нормировка Шепарда не применяется), общее для обоих осевых источников,
    // и закреплять его точным числом здесь не следует: проверка существует
    // ради того, чтобы член не исчез совсем, а ошибка давала ровно ноль.
    check(near_rel(r.dvz, expected, expected, 0.5),
          "E1 dv_z/dt ~ T/(rho r): осевой источник от сдвига присутствует");
    check(near_rel(r.dvr, 0.0, expected, 1e-3),
          "E2 dv_r/dt == 0: однородный сдвиг не создаёт радиальной силы");
}

// --- F. Поиск пар и сборка примитивов при разном числе потоков ---------------
// Осесимметричная решётка с зеркальными (ghost) частицами у оси: массивы
// поиска (primitives, real_pairs, contact_pairs) обязаны совпасть ПОБАЙТНО при
// 1 и 7 потоках (сборка по индексам: счёт -> префиксная сумма -> запись).
// Решётка больше основной площадки, чтобы блоков сборки было несколько.
void test_search_threads()
{
    std::printf("F. поиск и сборка примитивов: 1 и 7 потоков\n");
    constexpr size_t nr = 60, nz = 60;
    const size_t n = nr * nz;

    ParticleSystem sys;
    sys.resize(n);
    for (size_t ir = 0; ir < nr; ++ir)
        for (size_t iz = 0; iz < nz; ++iz)
        {
            const size_t k = ir * nz + iz;
            sys.x[0][k]  = (static_cast<double>(ir) + 0.5) * kDx;
            sys.x[1][k]  = static_cast<double>(iz) * kDx;
            sys.itype[k] = kItype;
            sys.rho[k]   = kRho0;
            sys.rho_0[k] = kRho0;
            sys.mass[k]  = kRho0 * kDx * kDx;
            sys.hsml[k]  = 1.3 * kDx;
            sys.axis[k]  = (ir == 0) ? 1 : 0;
        }

    SimulationConfig cfg;
    cfg.pa_sph          = 2;
    cfg.skf             = 1;
    cfg.density_method  = 1;
    cfg.maxn            = n;
    cfg.max_interaction = n * 200;

    rhs::CpuExecutor exec1(1);
    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec1, cfg.mats, sys, n, derived);
    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, n, ghosts);
    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, n, derived, ghosts, ext);

    interactions::SolidSphSearchResult s1, s7;
    interactions::build_solid_sph_search(exec1, cfg, ext, s1);
    interactions::build_solid_sph_search(rhs::CpuExecutor(7), cfg, ext, s7);

    size_t ghost_prims = 0;
    for (const auto& pr : s1.primitives) if (pr.j >= ext.nreal) ++ghost_prims;
    std::printf("  пар: %zu, primitives: %zu (с ghost: %zu), real_pairs: %zu, блоков: %zu\n",
                s1.pairs.size(), s1.primitives.size(), ghost_prims, s1.real_pairs.size(),
                s1.blocks.size() - 1);
    check(ghost_prims > 0, "F1 в сцене есть примитивы с ghost-частицами");
    check(s1.blocks.size() > 2, "F2 блоков сборки несколько");

    // Порядок пар у find_pairs_parallel зависит от числа потоков, поэтому
    // каждый результат сверяется с последовательной сборкой (прежний цикл с
    // push_back) ПО СВОИМ ЖЕ парам, а между потоками -- как множество.
    auto reference = [&](const interactions::SolidSphSearchResult& s,
                         std::vector<interactions::SolidSphPrimitive>& prims,
                         std::vector<std::pair<rhs::index_t, rhs::index_t>>& real) {
        for (const auto& pr : s.pairs)
        {
            const rhs::index_t i = pr.first, j = pr.second;
            double dx[dim];
            double r2 = 0.0;
            for (size_t d = 0; d < dim; ++d) { dx[d] = ext.x[d][i] - ext.x[d][j]; r2 += dx[d] * dx[d]; }
            const double r = std::sqrt(r2);
            const double mh = 0.5 * (ext.hsml[i] + ext.hsml[j]);
            interactions::SolidSphPrimitive prim;
            prim.i = i; prim.j = j;
            double w = 0.0;
            double dwdx_buf[dim];
            kernel(cfg.skf, r, std::span<const double>(dx, dim), mh, w, dwdx_buf);
            prim.w = w;
            for (size_t d = 0; d < dim; ++d) prim.dwdx[d] = dwdx_buf[d];
            prims.push_back(prim);
            if (j < ext.nreal) real.emplace_back(i, j);
        }
    };
    auto same_as_reference = [&](const interactions::SolidSphSearchResult& s) {
        std::vector<interactions::SolidSphPrimitive> prims;
        std::vector<std::pair<rhs::index_t, rhs::index_t>> real;
        reference(s, prims, real);
        return s.primitives.size() == prims.size() &&
               std::memcmp(s.primitives.data(), prims.data(), prims.size() * sizeof(prims[0])) == 0 &&
               s.real_pairs == real;
    };
    check(same_as_reference(s1), "F3 1 поток: primitives (побайтно) и real_pairs == последовательная сборка");
    check(same_as_reference(s7), "F4 7 потоков: primitives (побайтно) и real_pairs == последовательная сборка");
    check(s1.contact_pairs.empty() && s7.contact_pairs.empty(), "F5 одно тело: contact_pairs пуст");
    auto sorted_pairs = [](std::vector<std::pair<rhs::index_t, rhs::index_t>> v) {
        std::sort(v.begin(), v.end());
        return v;
    };
    check(s1.primitives.size() == s7.primitives.size() &&
              sorted_pairs(s1.real_pairs) == sorted_pairs(s7.real_pairs),
          "F6 множества real_pairs и число примитивов совпадают при 1 и 7 потоках");
}

// --- G. Работа давления на окружной деформации ------------------------------
// du/dt = -(p/rho) div v, а div v в осесимметрии = плоская дивергенция + v_r/r.
// Плоская часть приходит из парной суммы (SolidDelta::du), окружную
// finalize_solid добавляет сам: -(p/rho) v_r/r. Одна частица, парные суммы
// нулевые, девиатор нулевой -- du должно быть ровно этим членом.
void test_hoop_pressure_work()
{
    std::printf("\n[G] Работа давления на окружной деформации -(p/rho) v_r/r\n");
    ParticleSystem sys;
    sys.resize(1);
    constexpr double r = 3.0e-3, vr = 50.0, rho = 2600.0, P = 2.0e9;
    sys.x[0][0] = r;     sys.x[1][0] = 1.0e-3;
    sys.vx[0][0] = vr;   sys.vx[1][0] = -10.0;
    sys.rho[0] = rho;    sys.rho_0[0] = kRho0;
    sys.mass[0] = kRho0 * kDx * kDx;
    sys.itype[0] = kItype;

    rhs::CpuExecutor exec;
    SimulationConfig cfg;
    const std::vector<interactions::SolidDelta> reduced(1);
    std::vector<interactions::SolidFinal> f0, fp;
    interactions::finalize_solid(exec, cfg.mats, sys, 1, reduced, std::vector<double>{0.0}, f0);
    interactions::finalize_solid(exec, cfg.mats, sys, 1, reduced, std::vector<double>{P}, fp);

    const double expect = -P / rho * vr / r;
    std::printf("        du(p=0) = %.6e  du(p) = %.6e  ожидается %.6e (Дж/(кг*с))\n", f0[0].du, fp[0].du, expect);
    check(f0[0].du == 0.0, "G1 p = 0: du = 0 (нет парных сумм и девиатора)");
    check(near_rel(fp[0].du, expect, expect, 1e-14), "G2 du = -(p/rho) v_r/r");
    check(fp[0].drho == f0[0].drho, "G3 drho от p не зависит");
}

// --- G4. Масса кольца для учёта энергии ------------------------------------
// energy_mass = 2*pi*r0*m, где r0 -- начальный радиус (init_reference_radius),
// а не текущий: масса кольца постоянна, даже когда частица смещается по r.
void test_energy_mass_ring()
{
    std::printf("\n[G4] energy_mass = 2*pi*r0*m (r0 фиксирован)\n");
    ParticleSystem sys;
    sys.resize(2);
    sys.x[0][0] = 1.0e-3;  sys.x[0][1] = 4.0e-3;
    sys.mass[0] = 2.0;     sys.mass[1] = 3.0;
    sys.init_reference_radius(2);
    sys.x[0][0] = 2.5e-3;  sys.x[0][1] = 7.0e-3; // сдвиг: r0 меняться не должен

    constexpr double pi = 3.14159265358979323846;
    const double e0 = 2.0 * pi * 1.0e-3 * 2.0, e1 = 2.0 * pi * 4.0e-3 * 3.0;
    check(sys.r0[0] == 1.0e-3 && sys.r0[1] == 4.0e-3, "G4a r0 не изменился после сдвига x");
    check(near_rel(energy_mass(sys, 0), e0, e0, 1e-14), "G4b energy_mass(0) = 2*pi*r0*m");
    check(near_rel(energy_mass(sys, 1), e1, e1, 1e-14), "G4c energy_mass(1) = 2*pi*r0*m");
}

// --- G5. Парный баланс энергии при гомологическом расширении (диагностика) --
// Сырой reduced (без finalize_solid): dv[d], du. Состояние: однородное p > 0
// (rho = 1.01*rho_0), v_r = k*r, v_z = 0, S = 0, u = 0, mass "плоская". Решётка
// по r растянута: r = lam*r0 (плотность НЕ пересчитывается -- проверяется
// правая часть, не физическое состояние).
// Мощность на частицу: w_i = m_i*(v_i.dv_i + 0.5*du_i) (0.5 -- как в finalize_solid).
// Суммы с весами: 1 (P_flat), 2*pi*r0 (P_r0, = energy_mass), 2*pi*r (P_r).
// Масштаб S = sum m*|v_r|*|dv_r|. "Ось" = ir < 3 (носители с ghost-частицами).
struct PairEnergy
{
    double flat_all = 0.0, flat_in = 0.0, flat_axis = 0.0;
    double r0_all = 0.0, r0_in = 0.0;
    double r_all = 0.0, r_in = 0.0;
    double scale_all = 0.0, scale_in = 0.0;
    double p_min = 0.0, p_max = 0.0;
};

PairEnergy pair_energy_balance(double lam, int pair_mass = 0)
{
    constexpr size_t nr = 24, nz = 24;
    constexpr double k = 1.0e4;
    constexpr double pi = 3.14159265358979323846;
    const size_t n = nr * nz;

    ParticleSystem sys;
    sys.resize(n);
    for (size_t ir = 0; ir < nr; ++ir)
        for (size_t iz = 0; iz < nz; ++iz)
        {
            const size_t i = ir * nz + iz;
            sys.x[0][i]  = (static_cast<double>(ir) + 0.5) * kDx;
            sys.x[1][i]  = static_cast<double>(iz) * kDx;
            sys.itype[i] = kItype;
            sys.rho[i]   = 1.01 * kRho0;
            sys.rho_0[i] = kRho0;
            sys.mass[i]  = kRho0 * kDx * kDx;
            sys.hsml[i]  = 1.3 * kDx;
            sys.axis[i]  = (ir == 0) ? 1 : 0;
        }
    sys.init_reference_radius(n);           // r0 = (ir+0.5)*dx
    for (size_t i = 0; i < n; ++i)
    {
        sys.x[0][i]  = lam * sys.r0[i];
        sys.vx[0][i] = k * sys.x[0][i];
        sys.vx[1][i] = 0.0;
    }

    SimulationConfig cfg;
    cfg.pa_sph = 2; cfg.skf = 1; cfg.density_method = 1;
    cfg.axisym_pair_mass = pair_mass;
    cfg.maxn = n; cfg.max_interaction = n * 200;

    rhs::CpuExecutor exec;
    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec, cfg.mats, sys, n, derived);
    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, n, ghosts);
    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, n, derived, ghosts, ext);
    interactions::SolidSphSearchResult search;
    interactions::build_solid_sph_search(exec, cfg, ext, search);
    rhs::IndexedContributions<interactions::SolidDelta> buf;
    interactions::compute_solid_sph_contributions(exec, cfg, ext, search.primitives, buf);
    const std::vector<interactions::SolidDelta> red = rhs::reduce_contributions(exec, buf, n);

    PairEnergy o;
    o.p_min = o.p_max = derived.p[0];
    for (size_t i = 0; i < n; ++i)
    {
        o.p_min = std::min(o.p_min, derived.p[i]);
        o.p_max = std::max(o.p_max, derived.p[i]);
    }
    for (size_t ir = 0; ir < nr; ++ir)
        for (size_t iz = 0; iz < nz; ++iz)
        {
            const size_t i = ir * nz + iz;
            const double m = sys.mass[i];
            const double vdv = sys.vx[0][i] * red[i].dv[0] + sys.vx[1][i] * red[i].dv[1];
            const double w = m * (vdv + 0.5 * red[i].du);
            const double sc = m * std::fabs(sys.vx[0][i]) * std::fabs(red[i].dv[0]);
            o.flat_all += w;
            o.r0_all   += 2.0 * pi * sys.r0[i] * w;
            o.r_all    += 2.0 * pi * sys.x[0][i] * w;
            o.scale_all += sc;
            if (ir < 3) o.flat_axis += w;
            const bool inner = ir >= 3 && ir + 3 < nr && iz >= 3 && iz + 3 < nz;
            if (inner)
            {
                o.flat_in += w;
                o.r0_in   += 2.0 * pi * sys.r0[i] * w;
                o.r_in    += 2.0 * pi * sys.x[0][i] * w;
                o.scale_in += sc;
            }
        }
    return o;
}

void test_pair_energy_balance()
{
    std::printf("\n[G5] парный баланс энергии, гомологическое расширение v_r = k r (диагностика)\n");
    constexpr double pi = 3.14159265358979323846;
    for (const int pm : {0, 1})
    for (const double lam : {1.0, 1.2})
    {
        const PairEnergy o = pair_energy_balance(lam, pm);
        // P_r0, P_r несут лишний множитель 2*pi*r (размерность) -- делим на
        // 2*pi*12*dx (середина решётки), чтобы сравнивать с P_flat/S.
        const double wref = 2.0 * pi * 12.0 * kDx;
        std::printf("  [%s] lam = %.1f: p in [%.4e, %.4e] Pa\n", pm ? "ring" : "flat", lam, o.p_min, o.p_max);
        std::printf("    все:        P_flat/S = %+.3e, P_r0/(wref*S) = %+.3e, P_r/(wref*S) = %+.3e  (S = %.3e)\n",
                    o.flat_all / o.scale_all, o.r0_all / (wref * o.scale_all),
                    o.r_all / (wref * o.scale_all), o.scale_all);
        std::printf("    внутренние: P/S_all  = %+.3e, P_r0/(wref*S) = %+.3e, P_r/(wref*S) = %+.3e  (S = %.3e)\n",
                    o.flat_in / o.scale_all, o.r0_in / (wref * o.scale_all),
                    o.r_in / (wref * o.scale_all), o.scale_in);
        std::printf("    вклад зоны оси (ir<3) в P_flat(все): %+.3e от S\n", o.flat_axis / o.scale_all);
        check(o.p_min > 0.0 && std::fabs(o.p_max - o.p_min) <= 1e-9 * o.p_max,
              "G5 исходное p однородно и > 0");
        // G5a: ворот нет -- см. отчёт (ghost/hoop у оси не дают парного баланса).
    }
}

// --- G6. Чистая парная консервативность: одна пара вне оси, режим flat/ring --
// Две частицы (r_i < r_j, соседи в ядре, ghost нет: r > 2h), разные массы, плотности,
// скорости (сближение -- искусственная вязкость активна). Вклад пары, взвешенный
// физическим весом w = 2*pi*r0*m:  B = sum_{i,j} w_a*(v_a . dv_a + 0.5*du_a).
// Ring: B ~ 0 (машинный ноль относительно масштаба sum |w v.dv|); flat: B != 0
// (дисбаланс m_i m_j (r0_i - r0_j)/2 * h*(v_i+v_j) -- печатается).
double pair_balance_one(int mode, double& scale)
{
    constexpr double pi = 3.14159265358979323846;
    ParticleSystem sys;
    sys.resize(2);
    const double r_i = 10.0 * kDx, r_j = r_i + 1.1 * kDx;
    const double rr[2] = {r_i, r_j};
    const double zz[2] = {0.0, 0.4 * kDx};
    const double vr[2] = {+5.0, -3.0}, vz[2] = {1.5, -0.5};
    const double mm[2] = {kRho0 * kDx * kDx, 1.7 * kRho0 * kDx * kDx};
    const double rh[2] = {1.01 * kRho0, 1.03 * kRho0};
    for (size_t a = 0; a < 2; ++a)
    {
        sys.x[0][a] = rr[a]; sys.x[1][a] = zz[a];
        sys.vx[0][a] = vr[a]; sys.vx[1][a] = vz[a];
        sys.itype[a] = kItype; sys.rho[a] = rh[a]; sys.rho_0[a] = kRho0;
        sys.mass[a] = mm[a]; sys.hsml[a] = 1.3 * kDx; sys.axis[a] = 0;
    }
    sys.init_reference_radius(2);
    // r0 != r: рассинхронизируем (решётка деформирована), r0_i/r0_j независимы от x
    sys.r0[0] = 0.93 * rr[0]; sys.r0[1] = 1.08 * rr[1];

    SimulationConfig cfg;
    cfg.pa_sph = 2; cfg.skf = 1; cfg.density_method = 1;
    cfg.axisym_pair_mass = mode;
    cfg.maxn = 2; cfg.max_interaction = 400;

    rhs::CpuExecutor exec;
    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec, cfg.mats, sys, 2, derived);
    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, 2, ghosts);
    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, 2, derived, ghosts, ext);
    interactions::SolidSphSearchResult search;
    interactions::build_solid_sph_search(exec, cfg, ext, search);
    rhs::IndexedContributions<interactions::SolidDelta> buf;
    interactions::compute_solid_sph_contributions(exec, cfg, ext, search.primitives, buf);
    const std::vector<interactions::SolidDelta> red = rhs::reduce_contributions(exec, buf, 2);

    check(ghosts.source_index.empty() && search.primitives.size() == 1, "G6 одна пара без ghost");
    double B = 0.0;
    scale = 0.0;
    for (size_t a = 0; a < 2; ++a)
    {
        const double w = 2.0 * pi * sys.r0[a] * sys.mass[a];
        const double vdv = sys.vx[0][a] * red[a].dv[0] + sys.vx[1][a] * red[a].dv[1];
        B += w * (vdv + 0.5 * red[a].du);
        scale += w * std::fabs(vdv);
    }
    return B;
}

void test_pair_conservation()
{
    std::printf("\n[G6] парная консервативность (одна пара, вес 2*pi*r0*m)\n");
    double s_flat = 0.0, s_ring = 0.0;
    const double b_flat = pair_balance_one(0, s_flat);
    const double b_ring = pair_balance_one(1, s_ring);
    std::printf("  flat: B/scale = %+.3e  (scale = %.3e)\n", b_flat / s_flat, s_flat);
    std::printf("  ring: B/scale = %+.3e  (scale = %.3e)\n", b_ring / s_ring, s_ring);
    check(std::fabs(b_ring) / s_ring < 1e-12, "G6a ring: баланс пары < 1e-12");
    check(std::fabs(b_flat) / s_flat > 1e-6, "G6b flat: баланс пары != 0 при r0_i != r0_j");
}

} // namespace

int main()
{
    std::printf("=== Unit-тест осесимметричной ветки (axis_sim, dim = 2) ===\n");

    test_uniaxial_compression();
    test_radial_expansion();
    test_isotropic_compression();
    test_axial_source_at_rest();
    test_shear_source_at_rest();
    test_search_threads();
    test_hoop_pressure_work();
    test_energy_mass_ring();
    test_pair_energy_balance();
    test_pair_conservation();

    std::printf("\n=== %s (провалов: %d) ===\n",
                failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", failures);
    return failures == 0 ? 0 : 1;
}
