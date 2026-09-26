// eos.cpp
// Определения уравнений состояния (УРС) и функций давления.
// Объявления — в eos.hpp.
#include "eos.hpp"

void p_gas(real_t rho, real_t u, real_t& p, real_t& c)
{
    // Адиабатическое уравнение состояния газа
    const real_t gamma = real_t(1.4);
    p = (gamma - real_t(1.0)) * rho * u;
    if (p >= 0)
    {
        c = std::sqrt(gamma * p / rho);
    }
    else
    {
        c = 100;
    }
}

void p_art_water(real_t rho, real_t& p, real_t& c)
{
    // Упрощённое уравнение состояния для воды
    const real_t B = real_t(1000.0); // Модуль объёмной упругости
    p = B * (std::pow(rho / real_t(1000.0), real_t(7.0)) - real_t(1.0));
    c = std::sqrt(real_t(7.0) * B * std::pow(rho / real_t(1000.0), real_t(6.0)) / real_t(1000.0));
}

void eos_mi_grunaizen(const Materials& mats, size_t nmat, real_t rho, real_t rho_0, real_t u, real_t& p, real_t& c)
{
    // Материальные константы (ранее — глобалы param.hpp): aa, zz, gamm, n_Me.
    const auto& aa = mats.aa;
    const auto& zz = mats.zz;
    const auto& gamm = mats.gamm;
    const size_t n_Me = mats.n_Me;
    real_t w = rho / rho_0;
    real_t A1 = aa[nmat - n_Me] * pow(w, zz[nmat - n_Me]) * (zz[nmat - n_Me] - gamm[nmat - n_Me]) / (zz[nmat - n_Me] - 1);
    real_t A2 = aa[nmat - n_Me] * gamm[nmat - n_Me];
    real_t A3 = w * aa[nmat - n_Me] * zz[nmat - n_Me] * (gamm[nmat - n_Me] - 1) / (zz[nmat - n_Me] - 1);
    p = A1 - A2 + A3 + u * rho * (gamm[nmat - n_Me] - 1);
    real_t px = aa[nmat - n_Me] * (pow(w, zz[nmat - n_Me]) - 1);
    real_t C = (p - px) / real_t(rho);
    real_t zz1 = zz[nmat - n_Me] - real_t(1.0);
    real_t aa1 = aa[nmat - n_Me] * zz[nmat - n_Me] / rho_0;
    real_t as2 = aa1 * pow(w, zz1) + C * gamm[nmat - n_Me];
    if (as2 > 0)
    {
        c = sqrt(as2);
    }
    else
    {
        c = real_t(100.0);
    }
}
void eos_ua_teta(const Materials& mats, size_t nmat, real_t rho, real_t rho_0, real_t &p, real_t &c)
{
    // Материальные константы (ранее — глобалы param.hpp): aa, zz, n_Me.
    const auto& aa = mats.aa;
    const auto& zz = mats.zz;
    const size_t n_Me = mats.n_Me;
    real_t a1 = aa[nmat - n_Me];
    real_t z1 = zz[nmat - n_Me];
    real_t az1 = aa[nmat - n_Me] * zz[nmat - n_Me] / rho_0;
    real_t z11 = zz[nmat - n_Me] - real_t(1.0);
    p = a1 * (pow((rho / rho_0), z1) - real_t(1.0));
    real_t as2 = az1 * pow((rho / rho_0), z11);
    if (as2 > 0)
    {
        c = sqrt(as2);
    }
    else
    {
        c = real_t(100.);
    }
}
void urs_cem_JH2(const Materials& mats, size_t nmat, real_t qo, real_t qo0, real_t dp, real_t &as, real_t &pv, real_t &D)
{
    // Материальные константы (ранее — глобал param.hpp): const_JH_2.
    const auto& const_JH_2 = mats.const_JH_2;
    size_t mat = nmat - 10;
    real_t K1 = const_JH_2[mat][13];
    real_t K2 = const_JH_2[mat][14];
    real_t K3 = const_JH_2[mat][15];
    real_t T = const_JH_2[mat][16];
    real_t mu = (qo / qo0) - 1;
    if (D == 0)
    {
        if (mu >= 0)
        {
            pv = K1 * mu + K2 * std::pow(mu, real_t(2)) + K3 * std::pow(mu, real_t(3));
            as = sqrt((K1 / qo0) + 2 * (K2 / qo0) * mu + 3 * (K3 / qo0) * std::pow(mu, real_t(2)));
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
            pv = K1 * mu + K2 * std::pow(mu, real_t(2)) + K3 * std::pow(mu, real_t(3)) + dp;
            as = sqrt((K1 / qo0) + 2 * (K2 / qo0) * mu + 3 * (K3 / qo0) * std::pow(mu, real_t(2)));
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
