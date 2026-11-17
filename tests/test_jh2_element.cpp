// test_jh2_element.cpp
// Одноэлементный приёмочный расчёт модели керамики JH-2 (B4C, itype 10):
// одномерная деформация сжатия со скоростью 1e5 1/с до mu = 0.1 и разгрузка
// (jh2/jh2-fix-plan.md, этап E, п. 2; качественно -- рис. 8 J&H 1999).
// Одна частица, без SPH-соседей. На каждом шаге те же формулы, что в коде:
//   rho *= exp(-exx*dt); p, c -- urs_cem_JH2;
//   девиатор (dim == 1): dsxx = 2G*(exx + drho/(3 rho)),
//   dsyy = dszz = 2G*drho/(3 rho); затем prochnost_cem_JH2.
// Соглашение знаков кода: p > 0 при сжатии, exx < 0 при сжатии,
// sxx < 0 при сжатии, продольное напряжение sigma_x = -p + sxx (сжатие < 0).
// Собирается с SPH_DIM=1. CSV: --csv=<файл>
// (t, mu, p, sxx, sigma_x, sigma_eq, D, e_pl, dp).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include "param.hpp"
#include "materials.hpp"
#include "particle_system.hpp"
#include "eos.hpp"
#include "ppochnost.hpp"

namespace
{
int failures = 0;

void check(bool ok, const char* name, double got, double ref)
{
    std::printf("%s %s: получено=%.9e ожидалось=%.9e\n", ok ? "OK  " : "FAIL", name, got, ref);
    if (!ok)
        ++failures;
}

// Независимый эталон предела текучести sigma_eq (формулы (3)-(6) J&H 1999).
double ref_strength(const std::vector<double>& c, double P_n, double D, double eps_i)
{
    const double s_hel = c[3], T_n = c[4], A = c[5], N = c[6], B = c[7], M = c[8], sfmax = c[9], C = c[10];
    double s0 = 0.0, sf = 0.0;
    if (P_n <= 0.0)
    {
        if (P_n + T_n > 0.0)
            s0 = s_hel * A * std::pow(P_n + T_n, N);
    }
    else
    {
        s0 = s_hel * A * std::pow(P_n + T_n, N);
        sf = s_hel * B * std::pow(P_n, M);
    }
    if (sf > sfmax * s_hel)
        sf = sfmax * s_hel;
    double s = s0 - D * (s0 - sf);
    if (eps_i > 0.0)
        s *= (1.0 + C * std::log(eps_i));
    return s;
}

double hydro(const std::vector<double>& c, double mu)
{
    return c[13] * mu + c[14] * mu * mu + c[15] * mu * mu * mu;
}
} // namespace

int main(int argc, char** argv)
{
    std::string csv_path;
    for (int k = 1; k < argc; ++k)
        if (std::strncmp(argv[k], "--csv=", 6) == 0)
            csv_path = argv[k] + 6;

    Materials mats;
    const auto& c = mats.const_JH_2[0];
    const double rho0 = c[0], G = c[1], P_hel = c[2], K1 = c[13];
    const double edot = 1.0e5;               // 1/с
    const double eps_i = (2.0 / 3.0) * edot; // эквивалентная скорость в коде (exx = edot)
    const double mu_max = 0.1;
    const int nsteps_load = 4000;
    const double dt = std::log(1.0 + mu_max) / edot / nsteps_load;

    // --- ожидаемая точка HEL из условия пересечения (бисекция по mu) ---
    // упругая ветка: p = K(mu), sigma_eq = 2G*ln(1+mu) (девиатор 1D); текучесть
    // при sigma_eq = sigma_y(P*, D = 0).
    auto f = [&](double mu) { return 2.0 * G * std::log(1.0 + mu) - ref_strength(c, hydro(c, mu) / P_hel, 0.0, eps_i); };
    double lo = 1e-6, hi = 0.2;
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        (f(mid) > 0.0 ? hi : lo) = mid;
    }
    const double mu_hel = 0.5 * (lo + hi);
    const double sx_hel = hydro(c, mu_hel) + (4.0 / 3.0) * G * std::log(1.0 + mu_hel);
    std::printf("HEL ожидаемая: mu_HEL=%.6f |sigma_x|_HEL=%.4f ГПа\n", mu_hel, sx_hel * 1e-9);

    ParticleSystem sys;
    sys.resize(1);
    sys.inzializ_index(0);
    sys.itype[0] = 10;
    sys.rho_0[0] = rho0;
    sys.rho[0] = rho0;
    sys.p[0] = 0.0;
    sys.D[0] = 0.0;
    sys.dp[0] = 0.0;
    sys.e_pl[0] = 0.0;
    sys.sxx[0] = sys.syy[0] = sys.szz[0] = sys.sxy[0] = sys.syz[0] = sys.sxz[0] = 0.0;
    sys.eyy[0] = sys.ezz[0] = sys.exy[0] = sys.exz[0] = sys.eyz[0] = 0.0;

    FILE* fcsv = nullptr;
    if (!csv_path.empty())
    {
        fcsv = std::fopen(csv_path.c_str(), "w");
        if (fcsv)
            std::fprintf(fcsv, "# B4C, 1D, edot=%.3g 1/s, dt=%.6e; p>0 при сжатии, sigma_x=-p+sxx (сжатие < 0), sigma_eq=1.5|sxx|\n"
                               "t,mu,p,sxx,sigma_x,sigma_eq,D,e_pl,dp\n",
                         edot, dt);
    }

    double t = 0.0, mu = 0.0;
    bool yielded = false, unloading = false;
    double mu_yield = 0, sx_yield = 0, mu_D05 = -1, mu_D1 = -1, max_dp = 0, p_unload_005 = -1, hyd_unload_005 = 0;
    double prev_D = 0, prev_ep = 0, prev_dp = 0, max_el_dev = 0;
    bool mono_ok = true, range_ok = true, dp_pos_ok = true, soft_ok = true;
    double soft_p = 0, soft_hyd = 0, soft_dp = 0;
    int unload_checks = 0;
    long n = 0;
    const long nmax = 3L * nsteps_load + 10;

    while (n < nmax)
    {
        const double exx = unloading ? edot : -edot;
        sys.exx[0] = exx;
        const double rho_old = sys.rho[0];
        sys.rho[0] = rho_old * std::exp(-exx * dt);
        const double drho = -rho_old * exx;
        mu = sys.rho[0] / rho0 - 1.0;
        urs_cem_JH2(mats, 10, sys.rho[0], rho0, sys.dp[0], sys.c[0], sys.p[0], sys.D[0]);
        // девиатор по формулам int_force/solid_sph для dim == 1
        sys.sxx[0] += 2.0 * G * (exx + drho / (3.0 * rho_old)) * dt;
        sys.syy[0] += 2.0 * G * (drho / (3.0 * rho_old)) * dt;
        sys.szz[0] += 2.0 * G * (drho / (3.0 * rho_old)) * dt;
        prochnost_cem_JH2(mats, sys, 0, dt, drho);
        t += dt;
        ++n;

        const double sxx = sys.sxx[0], p = sys.p[0], D = sys.D[0], ep = sys.e_pl[0], dp = sys.dp[0];
        const double sigma_x = -p + sxx, sigma_eq = 1.5 * std::fabs(sxx);
        if (fcsv)
            std::fprintf(fcsv, "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", t, mu, p, sxx, sigma_x, sigma_eq, D, ep, dp);

        if (!unloading)
        {
            if (!yielded)
            {
                if (ep > 0.0)
                {
                    yielded = true;
                    mu_yield = mu;
                    sx_yield = std::fabs(sigma_x);
                }
                else
                {
                    // упругий участок: p = K(mu), sxx = -(4/3) G ln(1+mu)
                    const double dev = std::fabs(sxx + (4.0 / 3.0) * G * std::log(1.0 + mu)) / std::fabs(sxx);
                    const double dph = std::fabs(p - hydro(c, mu)) / std::fabs(p);
                    if (dev > max_el_dev)
                        max_el_dev = dev;
                    if (dph > max_el_dev)
                        max_el_dev = dph;
                    if (n == 1)
                    {
                        const double modulus = std::fabs(sigma_x) / mu;
                        const double ref = K1 + 4.0 * G / 3.0;
                        check(std::fabs(modulus - ref) <= 1e-3 * ref, "упругий модуль |sigma_x|/mu на первом шаге = K+4G/3", modulus, ref);
                    }
                }
            }
            if (mu_D05 < 0 && D >= 0.5)
                mu_D05 = mu;
            if (mu_D1 < 0 && D >= 1.0)
                mu_D1 = mu;
            if (mu >= mu_max)
                unloading = true;
        }
        else
        {
            // разгрузка: при D = 1 давление выше гидростаты на dp
            if (D >= 1.0 && mu > 0.0)
            {
                const double hyd = hydro(c, mu);
                ++unload_checks;
                if (!(p > hyd && std::fabs(p - hyd - dp) <= 1e-9 * std::fabs(p) && dp > 0.0))
                    soft_ok = false;
                soft_p = p;
                soft_hyd = hyd;
                soft_dp = dp;
            }
            if (p_unload_005 < 0 && mu <= 0.05)
            {
                p_unload_005 = p;
                hyd_unload_005 = hydro(c, mu);
            }
            if (p <= 0.0 || mu <= 0.0)
                break;
        }

        if (D < 0.0 || D > 1.0)
            range_ok = false;
        if (ep < 0.0 || ep < prev_ep || dp < prev_dp || D < prev_D)
            mono_ok = false;
        if (D > 1e-6 && !(dp > 0.0))
            dp_pos_ok = false;
        if (dp > max_dp)
            max_dp = dp;
        prev_D = D;
        prev_ep = ep;
        prev_dp = dp;
    }
    if (fcsv)
        std::fclose(fcsv);

    std::printf("шагов=%ld, dt=%.3e с\n", n, dt);
    std::printf("HEL (расчёт): mu=%.6f |sigma_x|=%.4f ГПа; D=0.5 при mu=%.6f; D=1 при mu=%.6f\n",
                mu_yield, sx_yield * 1e-9, mu_D05, mu_D1);
    std::printf("max dp=%.6e Па; p(разгрузка, mu=0.05)=%.6e Па, гидростата=%.6e Па, разность=%.6e Па\n",
                max_dp, p_unload_005, hyd_unload_005, p_unload_005 - hyd_unload_005);
    std::printf("последняя точка разгрузки с D=1: p=%.6e гидростата=%.6e dp=%.6e\n", soft_p, soft_hyd, soft_dp);

    check(max_el_dev < 1e-6, "упругий участок: совпадение с гидростатой и (4/3)G ln(1+mu)", max_el_dev, 0.0);
    check(yielded && std::fabs(sx_yield - sx_hel) <= 0.05 * sx_hel, "излом на HEL: |sigma_x| в пределах 5 %", sx_yield, sx_hel);
    check(yielded && std::fabs(mu_yield - mu_hel) <= 0.05 * mu_hel + 2.0 * mu_max / nsteps_load,
          "mu_HEL совпадает с условием пересечения", mu_yield, mu_hel);
    check(mu_D05 > 0 && mu_D1 > 0 && mu_D1 > mu_hel && mu_D1 - mu_hel < 0.03,
          "D растёт 0 -> 1 на малом интервале mu после HEL", mu_D1 - mu_yield, 0.03);
    check(range_ok, "D in [0,1]", 0, 0);
    check(mono_ok, "e_pl >= 0 и не убывает; D, dp не убывают", 0, 0);
    check(dp_pos_ok, "dp > 0 при D > 0", max_dp, 0);
    check(unload_checks > 0 && soft_ok, "разгрузка при D = 1: p = гидростата + dp > гидростаты", unload_checks, 1);

    std::printf("%s (failures=%d)\n", failures == 0 ? "ALL OK" : "FAILED", failures);
    return failures != 0;
}
