// test_jh2.cpp
// Тесты модели керамики JH-2 (src/ppochnost.cpp): приращение пластической
// деформации, sigma_i в плоском случае, разрыхление (dp) и повреждение (D).
// Обоснование -- jh2/jh2-code-review.md, план -- jh2/jh2-fix-plan.md.
// Собирается с SPH_DIM=1 (одномерная деформация: sxx -- единственная
// ненулевая компонента девиатора, sigma_eq = 1.5*|sxx|, J2 = 0.75*sxx^2).
// Все константы берутся из Materials (const_JH_2[0] -- B4C, itype 10).
#include <cmath>
#include <cstdio>
#include <vector>
#include "param.hpp"
#include "materials.hpp"
#include "particle_system.hpp"
#include "ppochnost.hpp"

namespace
{
int failures = 0;

void check(bool ok, const char* name, double got, double ref)
{
    const double rel = (ref != 0.0) ? std::fabs(got - ref) / std::fabs(ref) : std::fabs(got);
    std::printf("%s %s: код=%.12e эталон=%.12e отн.откл.=%.3e\n", ok ? "OK  " : "FAIL", name, got, ref, rel);
    if (!ok)
        ++failures;
}

bool close_rel(double a, double b, double tol)
{
    return std::fabs(a - b) <= tol * std::fabs(b);
}

// Эталонный предел текучести sigma_eq JH-2 (формулы (3)-(6) J&H 1999) для
// заданных P*, D и eps_i (эквивалентная скорость деформации в единицах кода).
// Написан независимо от ppochnost.cpp.
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

// Одна частица B4C, 1D, P* = P_n, rho = (1+mu)*rho_0, скорость деформации exx = eps_rate.
void setup(ParticleSystem& sys, const Materials& mats, double P_n, double mu, double eps_rate)
{
    const auto& c = mats.const_JH_2[0];
    sys.resize(1);
    sys.inzializ_index(0);
    sys.itype[0] = 10;
    sys.rho_0[0] = c[0];
    sys.rho[0] = c[0] * (1.0 + mu);
    sys.p[0] = P_n * c[2];
    sys.exx[0] = eps_rate;
    sys.D[0] = 0.0;
    sys.dp[0] = 0.0;
    sys.e_pl[0] = 0.0;
}
} // namespace

int main()
{
    Materials mats;
    const auto& c = mats.const_JH_2[0];
    const double G = c[1], T_n = c[4], K1 = c[13];
    const double D1 = c[11], D2 = c[12];
    const double rate = 1.0e4;               // eps_dot, 1/с
    const double eps_i = (2.0 / 3.0) * rate; // эквивалентная скорость в коде при exx = eps_dot
    const double dt = 1.0e-9;

    // --- T1: приращение пластической деформации ---
    {
        ParticleSystem sys;
        setup(sys, mats, 1.0, 0.0, rate);
        const double sy = ref_strength(c, 1.0, 0.0, eps_i);
        const double s_prev = sy / 1.5;
        const double trial = s_prev + 2.0 * G * (2.0 / 3.0) * rate * dt;
        sys.sxx[0] = trial;
        const double drho = -sys.rho[0] * rate;
        prochnost_cem_JH2(mats, sys, 0, dt, drho);
        const double ref = (2.0 / 3.0) * rate * dt;
        check(close_rel(sys.e_pl[0], ref, 1e-9), "T1 приращение e_pl", sys.e_pl[0], ref);
        // после возврата напряжение на поверхности текучести
        check(close_rel(1.5 * sys.sxx[0], sy, 1e-12), "T1 sigma_eq после возврата", 1.5 * sys.sxx[0], sy);
    }

    // --- T2: sigma_i в плоском случае ---
    {
        const double v = sigma_i_general(2, false, 1.0, 2.0, 0.5, 0.0, 0.0);
        check(std::fabs(v - 7.25) < 1e-12, "T2 sigma_i плоский", v, 7.25);
        const double va = sigma_i_general(2, true, 1.0, 2.0, 0.5, 0.0, 0.0);
        check(std::fabs(va - 7.25) < 1e-12, "T2 sigma_i осесимметричный", va, 7.25);
        const double v1 = sigma_i_general(1, true, 2.0, 0.0, 0.0, 0.0, 0.0);
        check(std::fabs(v1 - 3.0) < 1e-12, "T2 sigma_i 1D", v1, 3.0);
        const double v3 = sigma_i_general(3, false, 1.0, 2.0, 0.5, 0.25, 0.125);
        const double r3 = 1.0 + 2.0 + 4.0 + 0.25 + 0.015625 + 0.0625;
        check(std::fabs(v3 - r3) < 1e-12, "T2 sigma_i 3D", v3, r3);
    }

    // --- T3: разрыхление при росте D (формула (14) J&H 1999) ---
    {
        ParticleSystem sys;
        const double mu = 0.01;
        setup(sys, mats, 1.0, mu, rate);
        const double ef = D1 * std::pow(1.0 + T_n, D2);
        const double dD = 0.1;
        const double d_ep = dD * ef;
        const double s_old = ref_strength(c, 1.0, 0.0, eps_i);
        sys.sxx[0] = (s_old + 3.0 * G * d_ep) / 1.5;
        prochnost_cem_JH2(mats, sys, 0, dt, 0.0);
        const double s_new = ref_strength(c, 1.0, dD, eps_i);
        const double dU = (s_old * s_old - s_new * s_new) / (6.0 * G);
        const double K1mu = K1 * mu;
        const double dp_ref = -K1mu + std::sqrt(K1mu * K1mu + 2.0 * K1 * dU);
        check(close_rel(sys.D[0], dD, 1e-9), "T3 D", sys.D[0], dD);
        check(dp_ref > 0.0 && close_rel(sys.dp[0], dp_ref, 1e-6), "T3 dp", sys.dp[0], dp_ref);
        check(close_rel(sys.e_pl[0], d_ep, 1e-9), "T3 e_pl", sys.e_pl[0], d_ep);
    }

    // --- T4: нет разрыхления без роста D ---
    {
        ParticleSystem sys;
        setup(sys, mats, 1.0, 0.01, rate);
        sys.D[0] = 0.3;
        sys.dp[0] = 5.0e6;
        const double sy = ref_strength(c, 1.0, 0.3, eps_i);
        sys.sxx[0] = 0.5 * sy / 1.5; // упругий шаг
        const double s0 = sys.sxx[0];
        prochnost_cem_JH2(mats, sys, 0, dt, 0.0);
        check(sys.D[0] == 0.3 && sys.dp[0] == 5.0e6 && sys.e_pl[0] == 0.0 && sys.sxx[0] == s0,
              "T4a упругий шаг: D, dp, e_pl, sxx не меняются", sys.dp[0], 5.0e6);

        ParticleSystem s2;
        setup(s2, mats, 1.0, 0.01, rate);
        s2.D[0] = 1.0;
        s2.dp[0] = 5.0e6;
        const double sy1 = ref_strength(c, 1.0, 1.0, eps_i);
        s2.sxx[0] = 2.0 * sy1 / 1.5; // пластика при D = 1
        prochnost_cem_JH2(mats, s2, 0, dt, 0.0);
        check(s2.D[0] == 1.0 && s2.dp[0] == 5.0e6, "T4b D = 1: D и dp не меняются", s2.dp[0], 5.0e6);
        check(close_rel(1.5 * s2.sxx[0], sy1, 1e-12), "T4b sigma_eq на пределе остаточной прочности",
              1.5 * s2.sxx[0], sy1);
    }

    // --- T5: цельная керамика, пластика пренебрежимо мала: dp растёт ровно
    //     на ΔU(D_old -> D_new), а не на пробную работу пластики ---
    {
        ParticleSystem sys;
        const double mu = 0.01;
        setup(sys, mats, 1.0, mu, rate);
        const double ef = D1 * std::pow(1.0 + T_n, D2);
        const double dD = 1.0e-7;
        const double s_old = ref_strength(c, 1.0, 0.0, eps_i);
        sys.sxx[0] = (s_old + 3.0 * G * dD * ef) / 1.5;
        prochnost_cem_JH2(mats, sys, 0, dt, 0.0);
        const double s_new = ref_strength(c, 1.0, dD, eps_i);
        const double dU = (s_old - s_new) * (s_old + s_new) / (6.0 * G);
        const double K1mu = K1 * mu;
        const double dp_ref = -K1mu + std::sqrt(K1mu * K1mu + 2.0 * K1 * dU);
        check(dp_ref > 0.0 && close_rel(sys.dp[0], dp_ref, 1e-5), "T5 dp при малом D", sys.dp[0], dp_ref);
    }

    // --- T6 (D1): ограничение eps* = max(eps*, 1) за флагом jh2_clamp_low_rate ---
    {
        const double eps_low = 0.5; // eps* < 1: множитель 1 + C ln(eps*) < 1
        const double s_off = jh2_strength(c, 1.0, 0.0, eps_low, false);
        const double s_on = jh2_strength(c, 1.0, 0.0, eps_low, true);
        const double s_1 = ref_strength(c, 1.0, 0.0, 0.0); // без множителя скорости
        check(!mats.jh2_clamp_low_rate, "T6 флаг по умолчанию выключен", 0.0, 0.0);
        check(close_rel(s_off, ref_strength(c, 1.0, 0.0, eps_low), 1e-12) && s_off < s_1,
              "T6 флаг off, eps*<1: множитель < 1", s_off, ref_strength(c, 1.0, 0.0, eps_low));
        check(close_rel(s_on, s_1, 1e-12), "T6 флаг on, eps*<1: множитель = 1", s_on, s_1);
        const double s_hi = jh2_strength(c, 1.0, 0.0, 3.0, true);
        check(close_rel(s_hi, ref_strength(c, 1.0, 0.0, 3.0), 1e-12), "T6 флаг on, eps*>1: без изменений", s_hi,
              ref_strength(c, 1.0, 0.0, 3.0));
    }

    // --- T7 (D2): нижняя граница eps_f (jh2_efmin) ---
    {
        Materials m2;
        const double ef = D1 * std::pow(1.0 + T_n, D2);
        m2.jh2_efmin = 10.0 * ef;
        ParticleSystem sys;
        setup(sys, m2, 1.0, 0.01, rate);
        const double d_ep = 0.1 * ef;
        const double s_old = ref_strength(c, 1.0, 0.0, eps_i);
        sys.sxx[0] = (s_old + 3.0 * G * d_ep) / 1.5;
        prochnost_cem_JH2(m2, sys, 0, dt, 0.0);
        check(close_rel(sys.D[0], d_ep / m2.jh2_efmin, 1e-9), "T7 D при efmin = 10 ef", sys.D[0], d_ep / m2.jh2_efmin);
        // efmin = 0 (по умолчанию): D = d_ep / ef, как в оригинале
        ParticleSystem s0;
        setup(s0, mats, 1.0, 0.01, rate);
        s0.sxx[0] = (s_old + 3.0 * G * d_ep) / 1.5;
        prochnost_cem_JH2(mats, s0, 0, dt, 0.0);
        check(close_rel(s0.D[0], 0.1, 1e-9), "T7 D при efmin = 0", s0.D[0], 0.1);
    }

    std::printf("%s (failures=%d)\n", failures == 0 ? "ALL OK" : "FAILED", failures);
    return failures != 0;
}
