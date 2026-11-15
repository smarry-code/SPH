// test_art_stress.cpp
// Unit-тест искусственного напряжения Монагана (include/interactions/artificial_stress.hpp,
// член в compute_solid_sph_contributions). Собирается для SPH_DIM = 1, 2, 3.
//
// Чистые функции (в любой сборке, все размерности):
//   A1  1D: растяжение sigma > 0 -> R = -eps*sigma, сжатие -> 0;
//   A2  2D: тензор Q diag(l1, l2) Q^T -> R = Q diag(r(l1), r(l2)) Q^T (1e-12),
//       в том числе изотропный и чисто сдвиговый;
//   A3  3D (Якоби): то же для 1000 случайных поворотов и главных значений обоих знаков;
//   A4  всестороннее сжатие -> R = 0 точно;
//   A5  pair_factor: f(r = dp) = 1, f^n, W(dp) <= 0 -> 0.
// Конвейер (облако частиц STEEL 4340, растяжение p < 0, случайный девиатор):
//   P1  art_stress_eps = 0: результат побитово равен расчёту без build_art_stress;
//   P2  eps = 0.3: член сохраняет импульс  sum m_i (dv_on - dv_off) = 0 (в осесимметрии -- по z);
//   P3  ... и энергию  sum m_i v_i.(dv_on - dv_off) + sum m_i (du_on - du_off) = 0;
//   P4  сжатие (p > 0, s = 0): R = 0 у всех, результат побитово равен eps = 0;
//   P5  член ненулевой при растяжении (иначе P2-P3 проверяли бы пустоту).
#include "interactions/artificial_stress.hpp"
#include "interactions/solid_sph.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace
{
int failures = 0;

void check(bool cond, const char* what)
{
    if (!cond)
    {
        std::printf("FAIL: %s\n", what);
        ++failures;
    }
}

bool close(double a, double b, double tol)
{
    return std::fabs(a - b) <= tol * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

constexpr double kEps = 0.3;

void test_pure()
{
    // A1
    check(close(art_stress::r_1d(kEps, 2.0e8, -1.0e8), -kEps * 3.0e8, 1e-15), "A1 растяжение 1D");
    check(art_stress::r_1d(kEps, -2.0e8, 1.0e8) == 0.0, "A1 сжатие 1D -> 0");

    // A2: 2D
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    double worst2 = 0.0;
    for (int t = 0; t < 1000; ++t)
    {
        const double th = 3.14159265358979 * U(rng);
        const double l1 = 1e9 * U(rng), l2 = 1e9 * U(rng);
        const double c = std::cos(th), s = std::sin(th);
        const double a = l1 * c * c + l2 * s * s, b = l1 * s * s + l2 * c * c, xy = (l1 - l2) * c * s;
        const double r1 = art_stress::principal_r(kEps, l1), r2 = art_stress::principal_r(kEps, l2);
        const double exx = r1 * c * c + r2 * s * s, eyy = r1 * s * s + r2 * c * c, exy = (r1 - r2) * c * s;
        double rxx, ryy, rxy;
        art_stress::r_2d(kEps, a, b, xy, rxx, ryy, rxy);
        worst2 = std::max({worst2, std::fabs(rxx - exx), std::fabs(ryy - eyy), std::fabs(rxy - exy)});
    }
    check(worst2 <= 1e-12 * 1e9 * kEps * 4, "A2 2D: R = Q diag(r) Q^T");
    {
        double rxx, ryy, rxy;
        art_stress::r_2d(kEps, 5e8, 5e8, 0.0, rxx, ryy, rxy);
        check(close(rxx, -kEps * 5e8, 1e-15) && close(ryy, -kEps * 5e8, 1e-15) && rxy == 0.0, "A2 изотропное растяжение");
        art_stress::r_2d(kEps, 0.0, 0.0, 4e8, rxx, ryy, rxy);
        // главные +4e8 (ось 45°) и -4e8: R = -eps*4e8 * e1 e1^T, e1 = (1,1)/sqrt2
        check(close(rxx, -0.5 * kEps * 4e8, 1e-14) && close(ryy, -0.5 * kEps * 4e8, 1e-14)
                  && close(rxy, -0.5 * kEps * 4e8, 1e-14),
              "A2 чистый сдвиг");
    }

    // A3: 3D
    double worst3 = 0.0;
    for (int t = 0; t < 1000; ++t)
    {
        // случайный поворот: ортонормировка Грама-Шмидта
        double q[3][3];
        for (auto& row : q) for (double& v : row) v = U(rng);
        for (int i = 0; i < 3; ++i)
        {
            for (int k = 0; k < i; ++k)
            {
                double d = 0.0;
                for (int m = 0; m < 3; ++m) d += q[i][m] * q[k][m];
                for (int m = 0; m < 3; ++m) q[i][m] -= d * q[k][m];
            }
            double nn = 0.0;
            for (int m = 0; m < 3; ++m) nn += q[i][m] * q[i][m];
            nn = std::sqrt(nn);
            for (int m = 0; m < 3; ++m) q[i][m] /= nn;
        }
        const double l[3] = {1e9 * U(rng), 1e9 * U(rng), 1e9 * U(rng)};
        const double r[3] = {art_stress::principal_r(kEps, l[0]), art_stress::principal_r(kEps, l[1]),
                             art_stress::principal_r(kEps, l[2])};
        auto build = [&](const double* lam, int i, int j) {
            return lam[0] * q[0][i] * q[0][j] + lam[1] * q[1][i] * q[1][j] + lam[2] * q[2][i] * q[2][j];
        };
        art_stress::Sym3 sig;
        sig.xx = build(l, 0, 0); sig.yy = build(l, 1, 1); sig.zz = build(l, 2, 2);
        sig.xy = build(l, 0, 1); sig.yz = build(l, 1, 2); sig.xz = build(l, 0, 2);
        const art_stress::Sym3 R = art_stress::r_3d(kEps, sig);
        worst3 = std::max({worst3, std::fabs(R.xx - build(r, 0, 0)), std::fabs(R.yy - build(r, 1, 1)),
                           std::fabs(R.zz - build(r, 2, 2)), std::fabs(R.xy - build(r, 0, 1)),
                           std::fabs(R.yz - build(r, 1, 2)), std::fabs(R.xz - build(r, 0, 2))});
    }
    std::printf("A2 max err = %.3e, A3 max err = %.3e (Па)\n", worst2, worst3);
    check(worst3 <= 1e-12 * 1e9 * kEps * 6, "A3 3D Якоби: R = Q diag(r) Q^T");

    // A4
    {
        art_stress::Sym3 sig;
        sig.xx = -1e9; sig.yy = -2e9; sig.zz = -3e8; sig.xy = 1e8; sig.yz = -5e7; sig.xz = 2e7;
        const art_stress::Sym3 R = art_stress::r_3d(kEps, sig);
        check(R.xx == 0.0 && R.yy == 0.0 && R.zz == 0.0 && R.xy == 0.0 && R.yz == 0.0 && R.xz == 0.0,
              "A4 всестороннее сжатие -> R = 0");
    }

    // A5
    check(art_stress::pair_factor(0.7, 0.7, 4.0) == 1.0, "A5 f(r = dp) = 1");
    check(close(art_stress::pair_factor(0.35, 0.7, 4.0), 0.0625, 1e-15), "A5 f^n");
    check(art_stress::pair_factor(0.35, 0.0, 4.0) == 0.0, "A5 W(dp) = 0 -> 0");
}

struct Result
{
    std::vector<interactions::SolidFinal> fin;
};

// Облако частиц 4340 на сетке с шумом; tension -- растяжение (rho < rho0, случайный s),
// иначе сжатие (rho > rho0, s = 0).
ParticleSystem make_cloud(bool tension, std::size_t& n)
{
    constexpr std::size_t kSide = (dim == 1) ? 40 : (dim == 2 ? 14 : 7);
    constexpr double dx = 1.0e-4, rho0 = 7830.0;
    n = 1;
    for (std::size_t d = 0; d < dim; ++d) n *= kSide;
    ParticleSystem sys;
    sys.resize(n);
    std::mt19937_64 rng(777);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    for (std::size_t k = 0; k < n; ++k)
    {
        sys.inzializ_index(k);
        std::size_t rem = k;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const std::size_t id = rem % kSide;
            rem /= kSide;
            // В осесимметрии облако вдали от оси: без ghost-пар (они несимметричны).
            const double off = (d == 0) ? 100.0 * dx : 0.0;
            sys.x[d][k] = off + (static_cast<double>(id) + 0.15 * U(rng)) * dx;
            sys.vx[d][k] = 30.0 * U(rng);
        }
        sys.itype[k] = 12; // STEEL 4340, Linear EOS
        sys.rho_0[k] = rho0;
        sys.rho[k] = tension ? rho0 * (0.99 + 0.002 * U(rng)) : rho0 * (1.01 + 0.002 * U(rng));
        double vol = 1.0;
        for (std::size_t d = 0; d < dim; ++d) vol *= dx;
        sys.mass[k] = rho0 * vol;
        sys.hsml[k] = 1.3 * dx;
        if (tension)
        {
            sys.sxx[k] = 3e8 * U(rng); sys.syy[k] = 3e8 * U(rng); sys.szz[k] = -(sys.sxx[k] + sys.syy[k]);
            if constexpr (dim >= 2) sys.sxy[k] = 2e8 * U(rng);
            if constexpr (dim == 3) { sys.syz[k] = 2e8 * U(rng); sys.sxz[k] = 2e8 * U(rng); }
        }
    }
    return sys;
}

Result run(const ParticleSystem& sys, std::size_t n, double eps, bool call_art, bool* any_r = nullptr)
{
    SimulationConfig cfg;
    cfg.apply_defaults();
    cfg.pa_sph = 2;
    cfg.skf = 1;
    cfg.density_method = 1;
    cfg.art_stress_eps = eps;
    rhs::CpuExecutor exec;
    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec, cfg.mats, sys, n, derived);
    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, n, ghosts);
    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, n, derived, ghosts, ext);
    if (call_art) interactions::build_art_stress(exec, cfg, sys, ghosts, ext);
    if (any_r)
    {
        *any_r = false;
        for (std::size_t k = 0; k < ext.art_rxx.size(); ++k)
            if (ext.art_rxx[k] != 0.0 || ext.art_ryy[k] != 0.0 || ext.art_rzz[k] != 0.0 || ext.art_rxy[k] != 0.0
                || ext.art_ryz[k] != 0.0 || ext.art_rxz[k] != 0.0)
                *any_r = true;
    }
    interactions::SolidSphSearchResult search;
    interactions::build_solid_sph_search(exec, cfg, ext, search);
    rhs::IndexedContributions<interactions::SolidDelta> buf;
    interactions::compute_solid_sph_contributions(exec, cfg, ext, search.primitives, buf);
    const std::vector<interactions::SolidDelta> reduced = rhs::reduce_contributions(exec, buf, n);
    Result r;
    interactions::finalize_solid(exec, cfg.mats, sys, n, reduced, derived.p, r.fin);
    return r;
}

bool bitwise_equal(const Result& a, const Result& b)
{
    if (a.fin.size() != b.fin.size()) return false;
    for (std::size_t i = 0; i < a.fin.size(); ++i)
    {
        if (a.fin[i].du != b.fin[i].du || a.fin[i].drho != b.fin[i].drho) return false;
        for (std::size_t d = 0; d < dim; ++d)
            if (a.fin[i].dv[d] != b.fin[i].dv[d]) return false;
    }
    return true;
}

void test_pipeline()
{
    std::size_t n = 0;
    const ParticleSystem sys = make_cloud(true, n);
    const Result none = run(sys, n, 0.0, false);
    const Result off = run(sys, n, 0.0, true);
    check(bitwise_equal(none, off), "P1 eps = 0 побитово равен расчёту без члена");

    bool any_r = false;
    const Result on = run(sys, n, kEps, true, &any_r);
    check(any_r, "P5 при растяжении R != 0");

    double mom[dim] = {}, mom_scale = 0.0, en = 0.0, en_scale = 0.0, ddv_max = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double m = sys.mass[i];
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double ddv = on.fin[i].dv[d] - off.fin[i].dv[d];
            mom[d] += m * ddv;
            mom_scale += std::fabs(m * ddv);
            en += m * sys.vx[d][i] * ddv;
            en_scale += std::fabs(m * sys.vx[d][i] * ddv);
            ddv_max = std::max(ddv_max, std::fabs(ddv));
        }
        const double ddu = on.fin[i].du - off.fin[i].du;
        en += m * ddu;
        en_scale += std::fabs(m * ddu);
    }
    // В осесимметрии осевой член (R_rr - R_tt)/r радиальный импульс не сохраняет
    // (как и у s): проверяется только осевая компонента z.
    const std::size_t d_mom0 = (axis_sim && dim == 2) ? 1 : 0;
    double mom_err = 0.0;
    for (std::size_t d = d_mom0; d < dim; ++d) mom_err = std::max(mom_err, std::fabs(mom[d]));
    std::printf("P: max|ddv| = %.3e м/с^2, импульс %.3e / %.3e, энергия %.3e / %.3e\n", ddv_max, mom_err, mom_scale,
                std::fabs(en), en_scale);
    check(ddv_max > 0.0, "P5 член меняет ускорения");
    // Разность on - off содержит округление основного члена (~1e-16 от |dv|): допуск 1e-9.
    check(mom_err <= 1e-9 * mom_scale, "P2 импульс сохраняется");
    check(std::fabs(en) <= 1e-9 * en_scale, "P3 энергия сохраняется");

    // P4: сжатие
    std::size_t n2 = 0;
    const ParticleSystem sys2 = make_cloud(false, n2);
    bool any_r2 = true;
    const Result c_off = run(sys2, n2, 0.0, true);
    const Result c_on = run(sys2, n2, kEps, true, &any_r2);
    check(!any_r2, "P4 сжатие: R = 0 у всех частиц");
    check(bitwise_equal(c_off, c_on), "P4 сжатие: результат побитово равен eps = 0");
}
} // namespace

int main()
{
    std::printf("test_art_stress, dim = %zu, axis_sim = %d\n", static_cast<std::size_t>(dim), axis_sim ? 1 : 0);
    test_pure();
    test_pipeline();
    if (failures == 0) std::printf("OK\n");
    return failures == 0 ? 0 : 1;
}
