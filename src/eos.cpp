// eos.cpp
// Определения уравнений состояния (УРС) и функций давления.
// Объявления — в eos.hpp.
#include "eos.hpp"

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
