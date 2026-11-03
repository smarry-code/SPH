// tests/test_jc_eos.cpp
// Unit-тесты eos_linear / eos_shock (fiber/jc-implementation-plan.md, §3.3, U15-U17).
#include "eos.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool same(double a, double b) { return !(a < b) && !(a > b); }

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1e-300, std::fabs(b)); }

// Независимая реализация адиабаты.
double pH_ref(const ShockEos& e, double rho0, double mu)
{
    if (mu <= 0.0) return rho0 * e.c0 * e.c0 * mu;
    const double d = 1.0 - (e.s1 - 1.0) * mu;
    return rho0 * e.c0 * e.c0 * mu * (1.0 + mu) / (d * d);
}

double eH_ref(const ShockEos& e, double rho0, double mu)
{
    if (mu <= 0.0) return 0.0;
    return 0.5 * pH_ref(e, rho0, mu) / rho0 * mu / (1.0 + mu);
}

void test_u15()
{
    std::printf("U15 Linear (4340)\n");
    const LinearEos e{159.0e9};
    const double rho0 = 7830.0;
    double p = 1.0, c = 0.0;
    eos_linear(e, rho0, rho0, p, c);
    check(same(p, 0.0), "p(rho0) == 0 exactly");
    check(near_rel(c, std::sqrt(e.K / rho0), 1e-14), "c = sqrt(K/rho0)");
    eos_linear(e, rho0 * 1.01, rho0, p, c);
    check(near_rel(p, 1.59e9, 1e-12), "p(mu=0.01) = 1.59e9");
    eos_linear(e, rho0 * 0.99, rho0, p, c);
    check(near_rel(p, -1.59e9, 1e-12), "p(mu=-0.01) = -1.59e9");
}

struct ShockCase
{
    const char* name;
    double rho0;
    ShockEos e;
};

const ShockCase kCases[] = {
    {"1006", 7896.0, {2.17, 4569.0, 1.49}},
    {"7039", 2770.0, {2.0, 5328.0, 1.338}},
};

void test_u16()
{
    std::printf("U16 Shock\n");
    for (const auto& k : kCases)
    {
        std::printf(" material %s\n", k.name);
        double p = 1.0, c = 0.0;
        eos_shock(k.e, k.rho0, k.rho0, 0.0, p, c);
        check(same(p, 0.0), "p(rho0, u=0) == 0 exactly");
        check(c > 0.0 && near_rel(c, k.e.c0, 1e-12), "c(rho0,u=0) = c0");

        bool ok_u0 = true, ok_h = true, ok_t = true;
        for (double mu : {0.01, 0.05, 0.1, 0.2, 0.3})
        {
            const double rho = k.rho0 * (1.0 + mu);
            const double pH = pH_ref(k.e, k.rho0, mu);
            const double eH = eH_ref(k.e, k.rho0, mu);
            eos_shock(k.e, rho, k.rho0, 0.0, p, c);
            ok_u0 = ok_u0 && near_rel(p, pH - k.e.gamma0 * k.rho0 * eH, 1e-12);
            // Точное тождество: на u = e_H(rho) давление лежит на адиабате.
            eos_shock(k.e, rho, k.rho0, eH, p, c);
            ok_h = ok_h && near_rel(p, pH, 1e-12);
        }
        for (double mu : {-0.01, -0.05, -0.1})
        {
            eos_shock(k.e, k.rho0 * (1.0 + mu), k.rho0, 0.0, p, c);
            ok_t = ok_t && near_rel(p, k.rho0 * k.e.c0 * k.e.c0 * mu, 1e-12);
        }
        check(ok_u0, "u=0, mu>0: p = p_H - G0*rho0*e_H");
        check(ok_h, "u=e_H, mu>0: p = p_H (on the Hugoniot)");
        check(ok_t, "u=0, mu<0: p = rho0*c0^2*mu (linear)");
    }
}

void test_u17()
{
    std::printf("U17 Shock sound speed vs finite difference\n");
    for (const auto& k : kCases)
    {
        std::printf(" material %s\n", k.name);
        bool ok = true;
        for (double mu : {-0.1, -0.03, 0.02, 0.05, 0.1, 0.2, 0.3})
        {
            for (double u : {0.0, 1.0e5, 1.0e6, 5.0e6})
            {
                const double rho = k.rho0 * (1.0 + mu);
                double p = 0.0, c = 0.0, pp = 0.0, pm = 0.0, cc = 0.0;
                eos_shock(k.e, rho, k.rho0, u, p, c);
                const double h = 1e-5 * rho;
                eos_shock(k.e, rho + h, k.rho0, u, pp, cc);
                eos_shock(k.e, rho - h, k.rho0, u, pm, cc);
                const double dpdrho = (pp - pm) / (2.0 * h);
                const double c2_fd = dpdrho + p / (rho * rho) * (k.e.gamma0 * k.rho0);
                const bool good = near_rel(c * c, c2_fd, 1e-6);
                if (!good) std::printf("    mu=%g u=%g c2=%.9g fd=%.9g\n", mu, u, c * c, c2_fd);
                ok = ok && good;
            }
        }
        check(ok, "c^2 matches central difference (rel 1e-6)");
    }
    // Защита: c <= 0 -> c0 (сильное растяжение + отрицательная энергия).
    double p, c;
    eos_shock(kCases[0].e, 0.5 * kCases[0].rho0, kCases[0].rho0, -1.0e9, p, c);
    check(std::isfinite(c) && c > 0.0, "fallback keeps c finite and positive");
    // Защита сингулярности D -> 0: конечные p, c.
    const double mu_sing = 1.0 / (kCases[0].e.s1 - 1.0);
    eos_shock(kCases[0].e, kCases[0].rho0 * (1.0 + mu_sing), kCases[0].rho0, 0.0, p, c);
    check(std::isfinite(p) && std::isfinite(c), "singularity D->0 gives finite p, c");
}

} // namespace

int main()
{
    test_u15();
    test_u16();
    test_u17();
    std::printf("%s (%d failures)\n", failures == 0 ? "ALL PASSED" : "FAILED", failures);
    return failures;
}
