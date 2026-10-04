// test_kinematics_3d.cpp
// Unit-тест трёхмерной ветки силового расчёта (dim == 3).
//
// ЗАЧЕМ. Сквозные тесты и регрессия считают только dim = 2: осесимметричная
// ветка проверена test_axisym, а ветки `dim == 3` в int_force не проверял
// никто. Между тем нити (fiber/PLAN.md) живут только в 3D, и Solid, с которым
// они контактируют, будет считаться именно этими ветками. В них шесть
// компонент девиатора, три компоненты поворота и девять слагаемых импульса —
// перепутать xz с yz или знак при повороте легко, а заметить в сквозном
// прогоне трудно: ошибка остаётся конечной и ничего не роняет.
//
// КАК. Кубическая решётка 11×11×11 частиц B4C; поле скоростей и напряжений
// задаётся аналитически, конвейер проходит один раз (соседи -> неразрывность
// -> внутренние силы), производные снимаются на центральной частице. Её
// носитель (2h = 2.6·dx) заполнен целиком: до границы пять слоёв.
//
// ПРОВЕРКИ опираются на свойства, которые не зависят от точности SPH-
// аппроксимации, — симметрию решётки и алгебру тензоров:
//   A  простой сдвиг по трём парам осей: ненулевыми выходят ровно та
//      скорость деформации и тот поворот, которые положено, а три случая
//      переходят друг в друга перестановкой осей;
//   B  объёмное сжатие не рождает девиатора; одноосное сжатие даёт
//      бесследовый девиатор, симметричный по двум поперечным осям;
//   C  жёсткое вращение вокруг каждой оси при произвольном бесследовом
//      девиаторе: скорости деформации нулевые, производная девиатора — ровно
//      поворот R·S − S·R (производная Яуманна), след и S:S не меняются;
//   D  покой при однородном напряжении: ускорения нет;
//   E  линейное поле одной компоненты напряжения σ_ab = k·x_b: ускорение
//      появляется только по оси a — ловит компоненту девиатора, попавшую не в
//      то уравнение движения.
// fiber/PLAN.md, этап Ф4-B2, шаг 1: probe() ниже больше НЕ вызывает старый
// путь (direct_find/con_density/int_force) напрямую — конвейер заменён на
// тот же набор вызовов, что src/single_step.cpp использует для продового
// расчёта (Solid derived -> GhostSet -> SolidExtended -> поиск -> Solid-SPH
// kernel -> reduce -> finalize), БЕЗ стадии Solid-Solid-Contact (тест не про
// контакт). Смысл проверок A-E не изменился; e/r/ds/a читаются из
// interactions::SolidFinal (шим-поля exx..eyz/rxy..ryz и dsxx..dsxz), а не
// из sys, куда их писал старый int_force.cpp напрямую.
#include "interactions/solid_sph.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

static_assert(dim == 3, "test_kinematics_3d требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

// |a − b| <= tol·scale. Масштаб задаёт каждая проверка сама: производные
// девиатора порядка 2Gγ ~ 1e16 Па/с, и абсолютный порог здесь не годится.
bool near(double a, double b, double scale, double tol)
{
    return std::fabs(a - b) <= tol * std::fabs(scale);
}

// --- Расчётная площадка -----------------------------------------------------
constexpr size_t kN     = 11;
constexpr size_t kProbe = 5;          // центр решётки
constexpr double kDx    = 1.0e-4;     // 0.1 мм
constexpr double kRho0  = 2506.0;     // B4C
constexpr int    kItype = 10;         // B4C, JH-2

size_t flat_index(size_t ix, size_t iy, size_t iz) { return (ix * kN + iy) * kN + iz; }

// Координата узла относительно пробной частицы: поля задаются в этих
// координатах, поэтому на пробной частице линейные поля равны нулю.
double coord(size_t i) { return (static_cast<double>(i) - static_cast<double>(kProbe)) * kDx; }

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;

using VelocityField = std::function<Vec3(const Vec3&)>;
// Полное напряжение (девиатор, p = 0 при rho = rho_0) в точке.
using StressField = std::function<Mat3(const Vec3&)>;

struct Probe
{
    Mat3 e{};      // скорость деформации (симметричная)
    Mat3 r{};      // поворот в соглашении кода: r_ab = ½(∂v_a/∂x_b − ∂v_b/∂x_a)
    Mat3 ds{};     // производная девиатора
    Vec3 a{};      // ускорение
    double G = 0.0;
};

Mat3 zero_stress(const Vec3&) { return Mat3{}; }

Probe probe(const VelocityField& vfield, const StressField& sfield = zero_stress)
{
    const size_t n = kN * kN * kN;
    ParticleSystem sys;
    sys.resize(n);

    for (size_t ix = 0; ix < kN; ++ix)
        for (size_t iy = 0; iy < kN; ++iy)
            for (size_t iz = 0; iz < kN; ++iz)
            {
                const size_t k = flat_index(ix, iy, iz);
                const Vec3 x{coord(ix), coord(iy), coord(iz)};
                const Vec3 v = vfield(x);
                const Mat3 s = sfield(x);
                sys.inzializ_index(k);
                for (size_t d = 0; d < 3; ++d)
                {
                    sys.x[d][k]  = x[d];
                    sys.vx[d][k] = v[d];
                }
                sys.itype[k] = kItype;
                sys.rho[k]   = kRho0;
                sys.rho_0[k] = kRho0;    // rho == rho_0 => p = 0
                sys.mass[k]  = kRho0 * kDx * kDx * kDx;
                sys.hsml[k]  = 1.3 * kDx;
                sys.sxx[k] = s[0][0]; sys.syy[k] = s[1][1]; sys.szz[k] = s[2][2];
                sys.sxy[k] = s[0][1]; sys.syz[k] = s[1][2]; sys.sxz[k] = s[0][2];
            }

    SimulationConfig cfg;
    cfg.pa_sph          = 2;
    cfg.skf             = 1;          // кубический сплайн
    cfg.density_method  = 1;
    cfg.maxn            = n;
    cfg.max_interaction = n * 200;     // не используется новым путём, оставлено для минимальности diff

    // Тот же набор вызовов, что src/single_step.cpp (новый путь, fiber/PLAN.md
    // Ф4-B1), без стадии Solid-Solid-Contact (тест не про контакт). В 3D
    // GhostSet всегда пуст (axis_sim && dim==2 ложно при dim==3).
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
    interactions::finalize_solid(exec, cfg.mats, sys, n, reduced, finalres);

    const size_t p = flat_index(kProbe, kProbe, kProbe);
    const interactions::SolidFinal& f = finalres[p];
    Probe out;
    out.e = {{{f.exx, f.exy, f.exz},
              {f.exy, f.eyy, f.eyz},
              {f.exz, f.eyz, f.ezz}}};
    out.r = {{{0.0, f.rxy, f.rxz},
              {-f.rxy, 0.0, f.ryz},
              {-f.rxz, -f.ryz, 0.0}}};
    out.ds = {{{f.dsxx, f.dsxy, f.dsxz},
               {f.dsxy, f.dsyy, f.dsyz},
               {f.dsxz, f.dsyz, f.dszz}}};
    out.a = {f.dv[0], f.dv[1], f.dv[2]};
    out.G = cfg.mats.const_JH_2[static_cast<size_t>(kItype - 10)][1];
    return out;
}

double max_abs(const Mat3& m)
{
    double s = 0.0;
    for (const auto& row : m) for (double v : row) s = std::max(s, std::fabs(v));
    return s;
}

Mat3 mul(const Mat3& a, const Mat3& b)
{
    Mat3 c{};
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 3; ++k) c[i][j] += a[i][k] * b[k][j];
    return c;
}

const char* kAxis = "xyz";

std::string comp(size_t a, size_t b) { return std::string("s") + kAxis[a] + kAxis[b]; }

// --- A. Простой сдвиг -------------------------------------------------------
// v_a = γ·x_b (a ≠ b). Ожидается: e_ab = γ/2, r_ab = γ/2 (соглашение кода),
// остальные компоненты e и r нулевые; ds_ab = 2G·e_ab. Три сдвига — (x по y),
// (y по z), (z по x) — переводятся друг в друга циклической перестановкой осей,
// и решётка к ней симметрична: значения обязаны совпасть до округления.
void test_simple_shear()
{
    std::printf("\n[A] Простой сдвиг v_a = γ x_b по трём парам осей\n");
    constexpr double gamma = 1.0e5;   // 1/с
    const std::array<std::array<size_t, 2>, 3> pairs = {{{0, 1}, {1, 2}, {2, 0}}};
    std::array<Probe, 3> res;
    for (size_t c = 0; c < 3; ++c)
    {
        const size_t a = pairs[c][0], b = pairs[c][1];
        res[c] = probe([a, b](const Vec3& x) {
            Vec3 v{};
            v[a] = gamma * x[b];
            return v;
        });
        const Probe& r = res[c];
        const double scale = 0.5 * gamma;
        std::printf("        v_%c = γ %c: e_%c%c = %.6e, r_%c%c = %.6e (γ/2 = %.3e)\n",
                    kAxis[a], kAxis[b], kAxis[a], kAxis[b], r.e[a][b], kAxis[a], kAxis[b],
                    r.r[a][b], scale);

        // Дискретизация (нормировка ядра на решётке) даёт единицы процентов.
        check(near(r.e[a][b], scale, scale, 0.05) && near(r.r[a][b], scale, scale, 0.05),
              std::string("A1 v_") + kAxis[a] + " = γ " + kAxis[b] + ": e_" + kAxis[a] +
                  kAxis[b] + " ≈ γ/2 и r_" + kAxis[a] + kAxis[b] + " ≈ γ/2");
        bool others_zero = true;
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
            {
                const bool target = (i == a && j == b) || (i == b && j == a);
                if (!target)
                {
                    others_zero = others_zero && near(r.e[i][j], 0.0, scale, 1e-10) &&
                                  near(r.r[i][j], 0.0, scale, 1e-10);
                }
            }
        check(others_zero, "A2   остальные компоненты e и r нулевые");
        check(near(r.ds[a][b], 2.0 * r.G * r.e[a][b], 2.0 * r.G * scale, 1e-10),
              "A3   d" + comp(a, b) + "/dt == 2G e_" + kAxis[a] + kAxis[b] +
                  " (без начального девиатора поворот не вносит вклада)");
    }
    const double scale = 0.5 * gamma;
    check(near(res[0].e[0][1], res[1].e[1][2], scale, 1e-12) &&
              near(res[0].e[0][1], res[2].e[2][0], scale, 1e-12) &&
              near(res[0].r[0][1], res[1].r[1][2], scale, 1e-12) &&
              near(res[0].r[0][1], res[2].r[2][0], scale, 1e-12),
          "A4 три сдвига совпадают при циклической перестановке осей x→y→z→x");
}

// --- B. Сжатие ----------------------------------------------------------------
// Объёмное сжатие v = −c·x не меняет формы: девиатор не рождается.
// Одноосное v_a = −k·x_a: девиатор бесследовый, по двум поперечным осям
// одинаковый, по оси сжатия — минус удвоенный поперечный; три оси сжатия
// переходят друг в друга перестановкой.
void test_compression()
{
    std::printf("\n[B] Сжатие: объёмное и одноосное\n");
    constexpr double c = 1.0e5;
    const Probe iso = probe([](const Vec3& x) { return Vec3{-c * x[0], -c * x[1], -c * x[2]}; });
    const double scale_iso = 2.0 * iso.G * c;
    std::printf("        объёмное: ds = (%.3e %.3e %.3e), масштаб 2Gc = %.3e\n",
                iso.ds[0][0], iso.ds[1][1], iso.ds[2][2], scale_iso);
    check(max_abs(iso.ds) <= 1e-10 * scale_iso, "B1 объёмное сжатие не рождает девиатора");

    std::array<double, 3> along{}, across{};
    for (size_t a = 0; a < 3; ++a)
    {
        const Probe r = probe([a](const Vec3& x) {
            Vec3 v{};
            v[a] = -c * x[a];
            return v;
        });
        const double scale = 2.0 * r.G * c;
        const size_t b = (a + 1) % 3, d = (a + 2) % 3;
        const double tr = r.ds[0][0] + r.ds[1][1] + r.ds[2][2];
        std::printf("        сжатие по %c: ds_%c%c = %.6e, поперечные %.6e %.6e, след %.2e\n",
                    kAxis[a], kAxis[a], kAxis[a], r.ds[a][a], r.ds[b][b], r.ds[d][d], tr);
        check(near(tr, 0.0, scale, 1e-10),
              std::string("B2 сжатие по ") + kAxis[a] + ": девиатор бесследовый");
        check(near(r.ds[b][b], r.ds[d][d], scale, 1e-12) &&
                  near(r.ds[a][a], -2.0 * r.ds[b][b], scale, 1e-10) && r.ds[a][a] < 0.0,
              std::string("B3 сжатие по ") + kAxis[a] + ": поперечные компоненты равны, " +
                  "осевая = −2·поперечная < 0");
        along[a] = r.ds[a][a];
        across[a] = r.ds[b][b];
    }
    const double scale = 2.0 * iso.G * c;
    check(near(along[0], along[1], scale, 1e-12) && near(along[0], along[2], scale, 1e-12) &&
              near(across[0], across[1], scale, 1e-12) && near(across[0], across[2], scale, 1e-12),
          "B4 три оси сжатия дают одинаковый отклик");
}

// --- C. Жёсткое вращение ------------------------------------------------------
// v = ω × x — вращение без деформации. Скорости деформации нулевые, а девиатор
// обязан лишь поворачиваться: dS/dt = R·S − S·R, где R — поворот в соглашении
// кода. Это и есть определение производной Яуманна в int_force:
//   dS_ab/dt = Σ_c (S_ac R_bc + S_cb R_ac) = (R·S − S·R)_ab.
// Ожидание собирается из R, который вычислил сам код: проверяется алгебра
// поворотных членов, а не точность SPH-градиента. Следствия поворота — след
// девиатора и S:S не меняются — проверяются отдельно: при ошибке знака в
// одной компоненте они и называют её.
void test_rigid_rotation()
{
    std::printf("\n[C] Жёсткое вращение при произвольном девиаторе\n");
    constexpr double omega = 1.0e5;   // рад/с
    // Бесследовый девиатор общего вида, все компоненты разные и ненулевые.
    const Mat3 S = {{{1.0e9, 0.3e9, -0.5e9},
                     {0.3e9, -0.4e9, 0.7e9},
                     {-0.5e9, 0.7e9, -0.6e9}}};
    const auto uniform = [S](const Vec3&) { return S; };

    for (size_t axis = 0; axis < 3; ++axis)
    {
        const Probe r = probe([axis](const Vec3& x) {
            Vec3 w{};
            w[axis] = omega;
            return Vec3{w[1] * x[2] - w[2] * x[1], w[2] * x[0] - w[0] * x[2],
                        w[0] * x[1] - w[1] * x[0]};
        }, uniform);

        const Mat3 rs = mul(r.r, S), sr = mul(S, r.r);
        Mat3 expected{};
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j) expected[i][j] = rs[i][j] - sr[i][j];

        const double scale = max_abs(r.r) * max_abs(S);
        const std::string tag = std::string(" (вращение вокруг ") + kAxis[axis] + ")";
        std::printf("        вокруг %c: |r| = %.3e, масштаб |r||S| = %.3e\n", kAxis[axis],
                    max_abs(r.r), scale);

        check(max_abs(r.e) <= 1e-10 * omega, "C1 скорости деформации нулевые" + tag);
        bool all_ok = true;
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = i; j < 3; ++j)
            {
                const bool ok = near(r.ds[i][j], expected[i][j], scale, 1e-9);
                if (!ok)
                {
                    std::printf("        d%s/dt = %.6e, ожидалось (R·S − S·R) = %.6e\n",
                                comp(i, j).c_str(), r.ds[i][j], expected[i][j]);
                }
                all_ok = all_ok && ok;
            }
        check(all_ok, "C2 dS/dt == R·S − S·R: девиатор только поворачивается" + tag);
        const double tr = r.ds[0][0] + r.ds[1][1] + r.ds[2][2];
        check(near(tr, 0.0, scale, 1e-9),
              "C3 след dS/dt == 0: поворот не делает девиатор шаровым" + tag);
        double sds = 0.0;
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j) sds += S[i][j] * r.ds[i][j];
        check(near(sds, 0.0, scale * max_abs(S), 1e-9),
              "C4 S:dS/dt == 0: интенсивность напряжений при повороте постоянна" + tag);
    }
}

// --- D. Покой при однородном напряжении ----------------------------------------
// Дивергенция однородного поля нулевая — ускорения быть не должно. Каждая из
// шести компонент задаётся отдельно.
void test_uniform_stress_at_rest()
{
    std::printf("\n[D] Покой + однородная компонента напряжения\n");
    constexpr double s0 = 1.0e9;
    const double scale = s0 / (kRho0 * kDx);   // ускорение от градиента s0 на шаге решётки
    for (size_t a = 0; a < 3; ++a)
        for (size_t b = a; b < 3; ++b)
        {
            const Probe r = probe([](const Vec3&) { return Vec3{}; }, [a, b](const Vec3&) {
                Mat3 s{};
                s[a][b] = s[b][a] = s0;
                return s;
            });
            const double amax = std::max({std::fabs(r.a[0]), std::fabs(r.a[1]), std::fabs(r.a[2])});
            check(amax <= 1e-10 * scale,
                  "D  однородная " + comp(a, b) + " не создаёт ускорения");
        }
}

// --- E. Линейное поле одной компоненты напряжения ------------------------------
// σ_ab = k·x_b (и симметричная σ_ba): дивергенция даёт ускорение
// k/ρ по оси a; для a ≠ b слагаемое σ_ba = k·x_b по оси b не дифференцируется
// (производная по x_a равна нулю). Значит ускорение — только по оси a.
// Величина проверяется грубо (дискретизация), направление — строго.
void test_linear_stress()
{
    std::printf("\n[E] Линейное поле одной компоненты напряжения σ_ab = k x_b\n");
    constexpr double k = 1.0e13;      // Па/м
    const double expected = k / kRho0;
    for (size_t a = 0; a < 3; ++a)
        for (size_t b = 0; b < 3; ++b)
        {
            const Probe r = probe([](const Vec3&) { return Vec3{}; }, [a, b](const Vec3& x) {
                Mat3 s{};
                s[a][b] = s[b][a] = k * x[b];
                return s;
            });
            const size_t o1 = (a + 1) % 3, o2 = (a + 2) % 3;
            const std::string name = comp(a, b) + " = k " + kAxis[b];
            std::printf("        %-12s a = (%.4e %.4e %.4e), ожидалось %.4e по %c\n",
                        name.c_str(), r.a[0], r.a[1], r.a[2], expected, kAxis[a]);
            check(near(r.a[a], expected, expected, 0.05) &&
                      near(r.a[o1], 0.0, expected, 1e-9) && near(r.a[o2], 0.0, expected, 1e-9),
                  "E  " + name + ": ускорение k/ρ только по оси " + kAxis[a]);
        }
}

} // namespace

int main()
{
    std::printf("=== Unit-тест трёхмерной ветки (dim = 3) ===\n");

    test_simple_shear();
    test_compression();
    test_rigid_rotation();
    test_uniform_stress_at_rest();
    test_linear_stress();

    std::printf("\n=== %s (провалов: %d) ===\n", failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", failures);
    return failures == 0 ? 0 : 1;
}
