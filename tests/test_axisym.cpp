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
#include "density.hpp"
#include "direct_find.hpp"
#include "int_force.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "virt_part.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
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
    sys.resize(n * 2);               // вторая половина — запас под зеркальные

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
    cfg.nnps            = 1;          // direct_find: без сетки, без её настроек
    cfg.pa_sph          = 2;
    cfg.skf             = 1;          // кубический сплайн
    cfg.density_method  = 1;
    cfg.maxn            = n * 2;
    cfg.max_interaction = n * 200;

    size_t nvirt = 0;
    virt_part(sys, n, nvirt, cfg.maxn);
    const size_t nall = n + nvirt;

    size_t niac = 0;
    std::vector<size_t> pair_i(cfg.max_interaction), pair_j(cfg.max_interaction);
    std::vector<double> w(cfg.max_interaction, 0.0);
    std::vector<std::vector<double>> dwdx(dim, std::vector<double>(cfg.max_interaction, 0.0));
    std::vector<int> countiac(nall, 0);
    direct_find(cfg, sys, nall, niac, pair_i, pair_j, w, dwdx, countiac);

    PairGraph graph;
    graph.num_colors = 0;             // atomic-ветвь потребителей пар

    std::vector<double> drho(nall, 0.0);
    con_density(cfg, sys, nall, niac, pair_i, pair_j, dwdx, drho, graph);

    std::vector<double> eta(nall, 0.0), t(nall, 0.0);
    std::vector<double> tdsdt(nall, 0.0), dedt(nall, 0.0);
    std::vector<double> dsxx(nall, 0.0), dsyy(nall, 0.0), dszz(nall, 0.0);
    std::vector<double> dsxy(nall, 0.0), dsyz(nall, 0.0), dsxz(nall, 0.0);
    std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(nall, 0.0));
    WorkBuffers wb;

    int_force(cfg, sys, nall, niac, eta, pair_i, pair_j, dwdx, w, t, dvxdt,
              tdsdt, dedt, dsxx, dsyy, dszz, dsxy, dsyz, dsxz, drho, wb, graph);

    const size_t probe = flat_index(kProbeIr, kProbeIz);
    Rates out;
    out.dsxx = dsxx[probe];
    out.dsyy = dsyy[probe];
    out.dszz = dszz[probe];
    out.dvr  = dvxdt[0][probe];
    out.dvz  = dvxdt[1][probe];
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

} // namespace

int main()
{
    std::printf("=== Unit-тест осесимметричной ветки (axis_sim, dim = 2) ===\n");

    test_uniaxial_compression();
    test_radial_expansion();
    test_isotropic_compression();
    test_axial_source_at_rest();
    test_shear_source_at_rest();

    std::printf("\n=== %s (провалов: %d) ===\n",
                failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", failures);
    return failures == 0 ? 0 : 1;
}
