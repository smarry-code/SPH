// eos.cpp
// Определения уравнений состояния (УРС) и функций давления.
// Объявления — в eos.hpp.
#include "eos.hpp"

#include <algorithm>

void p_gas(double rho, double u, double& p, double& c)
{
    // Адиабатическое уравнение состояния газа
    const double gamma = 1.4;
    p = (gamma - 1.0) * rho * u;
    if (p >= 0)
    {
        c = std::sqrt(gamma * p / rho);
    }
    else
    {
        c = 100;
    }
}

void p_art_water(double rho, double& p, double& c)
{
    // Упрощённое уравнение состояния для воды
    const double B = 1000.0; // Модуль объёмной упругости
    p = B * (std::pow(rho / 1000.0, 7.0) - 1.0);
    c = std::sqrt(7.0 * B * std::pow(rho / 1000.0, 6.0) / 1000.0);
}

void eos_mi_grunaizen(const Materials& mats, size_t nmat, double rho, double rho_0, double u, double& p, double& c)
{
    // Материальные константы (ранее — глобалы param.hpp): aa, zz, gamm, n_Me.
    const auto& aa = mats.aa;
    const auto& zz = mats.zz;
    const auto& gamm = mats.gamm;
    const size_t n_Me = mats.n_Me;
    double w = rho / rho_0;
    double A1 = aa[nmat - n_Me] * pow(w, zz[nmat - n_Me]) * (zz[nmat - n_Me] - gamm[nmat - n_Me]) / (zz[nmat - n_Me] - 1);
    double A2 = aa[nmat - n_Me] * gamm[nmat - n_Me];
    double A3 = w * aa[nmat - n_Me] * zz[nmat - n_Me] * (gamm[nmat - n_Me] - 1) / (zz[nmat - n_Me] - 1);
    p = A1 - A2 + A3 + u * rho * (gamm[nmat - n_Me] - 1);
    double px = aa[nmat - n_Me] * (pow(w, zz[nmat - n_Me]) - 1);
    double C = (p - px) / double(rho);
    double zz1 = zz[nmat - n_Me] - 1.0;
    double aa1 = aa[nmat - n_Me] * zz[nmat - n_Me] / rho_0;
    double as2 = aa1 * pow(w, zz1) + C * gamm[nmat - n_Me];
    if (as2 > 0)
    {
        c = sqrt(as2);
    }
    else
    {
        c = 100.0;
    }
}
void eos_ua_teta(const Materials& mats, size_t nmat, double rho, double rho_0, double &p, double &c)
{
    // Материальные константы (ранее — глобалы param.hpp): aa, zz, n_Me.
    const auto& aa = mats.aa;
    const auto& zz = mats.zz;
    const size_t n_Me = mats.n_Me;
    double a1 = aa[nmat - n_Me];
    double z1 = zz[nmat - n_Me];
    double az1 = aa[nmat - n_Me] * zz[nmat - n_Me] / rho_0;
    double z11 = zz[nmat - n_Me] - 1.0;
    p = a1 * (pow((rho / rho_0), z1) - 1.0);
    double as2 = az1 * pow((rho / rho_0), z11);
    if (as2 > 0)
    {
        c = sqrt(as2);
    }
    else
    {
        c = 100.;
    }
}
void urs_cem_JH2(const Materials& mats, size_t nmat, double qo, double qo0, double dp, double &as, double &pv, double &D)
{
    // Материальные константы (ранее — глобал param.hpp): const_JH_2.
    const auto& const_JH_2 = mats.const_JH_2;
    size_t mat = nmat - 10;
    double K1 = const_JH_2[mat][13];
    double K2 = const_JH_2[mat][14];
    double K3 = const_JH_2[mat][15];
    double T = const_JH_2[mat][16];
    double mu = (qo / qo0) - 1;
    if (D == 0)
    {
        if (mu >= 0)
        {
            pv = K1 * mu + K2 * pow(mu, 2) + K3 * pow(mu, 3);
            as = sqrt((K1 / qo0) + 2 * (K2 / qo0) * mu + 3 * (K3 / qo0) * pow(mu, 2));
        }
        else
        {
            pv = K1 * mu;
            if (pv < -T)
            {
                //pv = 0.0;
                pv = -T;
            }
            as = sqrt(K1 / qo0);
        }
    }
    else
    {
        if (mu >= 0)
        {
            pv = K1 * mu + K2 * pow(mu, 2) + K3 * pow(mu, 3) + dp;
            as = sqrt((K1 / qo0) + 2 * (K2 / qo0) * mu + 3 * (K3 / qo0) * pow(mu, 2));
        }
        else
        {
            pv = K1 * mu;
            if (pv < (-T * (1 - D)))
            {
                //pv = 0.0;
                pv = -T * (1 - D);
            }
            as = sqrt(K1 / qo0);
        }
    }
}

void eos_linear(const LinearEos& eos, double rho, double rho0, double& p, double& c)
{
    // Autodyn Linear: энергия не входит в давление.
    const double mu = rho / rho0 - 1.0;
    p = eos.K * mu;
    c = std::sqrt(eos.K / rho0);
}

void eos_shock(const ShockEos& eos, double rho, double rho0, double u, double& p, double& c)
{
    // Ми-Грюнайзен по ударной адиабате Us = c0 + S1·up (как в Autodyn).
    // Обозначения: μ = ρ/ρ0 − 1, A = ρ0·c0², D = 1 − (S1−1)μ, Γ = Γ0·ρ0.
    //   μ > 0: p_H = A·μ(1+μ)/D²,  e_H = ½·(p_H/ρ0)·μ/(1+μ) = ½·c0²·μ²/D²;
    //   μ ≤ 0: p_H = A·μ,          e_H = 0;
    //   p = p_H + Γ·(u − e_H).
    // Производные по μ (dμ/dρ = 1/ρ0):
    //   dp_H/dμ = A·[(1+2μ)·D + 2(S1−1)μ(1+μ)]/D³   (μ ≤ 0: A),
    //   de_H/dμ = c0²·μ/D³                           (μ ≤ 0: 0),
    //   ∂p/∂ρ|_u = (dp_H/dμ − Γ·de_H/dμ)/ρ0,   ∂p/∂u = Γ.
    //   c² = ∂p/∂ρ|_u + (p/ρ²)·∂p/∂u.
    // Сингулярность: D → 0 при сжатии μ → 1/(S1−1) (за пределами физического
    // диапазона). D ограничивается снизу значением kDmin, чтобы p и c оставались
    // конечными (давление в этой области — насыщенное, не бесконечное).
    constexpr double kDmin = 1.0e-3;
    const double mu = rho / rho0 - 1.0;
    const double A = rho0 * eos.c0 * eos.c0;
    const double gam = eos.gamma0 * rho0;

    double pH, eH, dpH, deH;
    if (mu > 0.0)
    {
        const double D = std::max(1.0 - (eos.s1 - 1.0) * mu, kDmin);
        const double D2 = D * D;
        const double D3 = D2 * D;
        pH = A * mu * (1.0 + mu) / D2;
        eH = 0.5 * eos.c0 * eos.c0 * mu * mu / D2;
        dpH = A * ((1.0 + 2.0 * mu) * D + 2.0 * (eos.s1 - 1.0) * mu * (1.0 + mu)) / D3;
        deH = eos.c0 * eos.c0 * mu / D3;
    }
    else
    {
        pH = A * mu;
        eH = 0.0;
        dpH = A;
        deH = 0.0;
    }
    p = pH + gam * (u - eH);

    const double dpdrho = (dpH - gam * deH) / rho0;
    const double c2 = dpdrho + (p / (rho * rho)) * gam;
    c = (std::isfinite(c2) && c2 > 0.0) ? std::sqrt(c2) : eos.c0;
}

void eos_jc(const JcMaterial& m, double rho, double rho0, double u, double D, double& p, double& c)
{
    if (m.eos == EosKind::Linear)
        eos_linear(m.lin, rho, rho0, p, c);
    else
        eos_shock(m.shock, rho, rho0, u, p, c);
    // Разрушенная частица: растягивающее давление обнуляется (аналог Hydro Pmin = 0).
    if (D >= 1.0 && p < 0.0) p = 0.0;
}
