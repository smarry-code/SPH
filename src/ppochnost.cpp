// ppochnost.cpp
// Определения моделей прочности (Мизес, JH-2) и интенсивности напряжений.
// Объявления — в ppochnost.hpp.
#include "ppochnost.hpp"
#include <stdexcept>

double sigma_i(double &sxx, double &syy, double &sxy, double &syz, double &sxz)
{
    return sigma_i_general(static_cast<int>(dim), axis_sim, sxx, syy, sxy, syz, sxz);
}

double sigma_i_general(int dim_, bool axis_, double sxx, double syy, double sxy, double syz, double sxz)
{
    if (dim_ == 1) // одномерная деформация
    {
        return 0.75 * sxx * sxx;
        // return (sxx * sxx + syy * syy + sxx * syy);
        // return (sxx - syy) * (sxx - syy) / 3.0;
    }
    else if (dim_ == 2)
    {

        if (axis_ == true)
        {
            return (sxx * sxx + sxx * syy + syy * syy + sxy * sxy); // осесимметричный случай
        }
        else
        {
            // Плоский случай: szz = -(sxx + syy) (след девиатора равен нулю),
            // J2 = (sxx^2 + syy^2 + szz^2)/2 + sxy^2 -- та же формула, что и
            // для осесимметричного случая. Прежнее sxx^2 + sxy^2 верно лишь
            // при syy = -sxx (чистый сдвиг).
            return (sxx * sxx + sxx * syy + syy * syy + sxy * sxy);
        }
    }
    else if (dim_ == 3) // трехмерный случай
    {
        return (sxx * sxx + sxx * syy + syy * syy + sxy * sxy + sxz * sxz + syz * syz);
    }
    throw std::logic_error("sigma_i_general: неподдерживаемая размерность dim_");
}
void prochnost_mises(const Materials& mats, size_t nmat, double &sxx, double &syy, double &szz, double &sxy, double &syz, double &sxz)
{
    // Материальные константы (ранее — глобалы param.hpp): G, sigy, n_Me.
    const auto& G = mats.G;
    const auto& sigy = mats.sigy;
    const size_t n_Me = mats.n_Me;
    double gem = G[nmat - n_Me];
    double ge2m = 2.0 * G[nmat - n_Me];
    double smiz2m = pow(sigy[nmat - n_Me], 2) / 3.0;
    double smizm = sqrt(pow(sigy[nmat - n_Me], 2) / 3.0);
    double ripl = sigma_i(sxx, syy, sxy, syz, sxz);
    bool flag_e = true;
    if (ripl > smiz2m)
    {
        double corr = smizm / sqrt(ripl);
        sxx *= corr;
        syy *= corr;
        szz *= corr;
        sxy *= corr;
        syz *= corr;
        sxz *= corr;
        flag_e = false;
    }
}
void prochnost_cem_JH2(const Materials& mats, ParticleSystem& sys, size_t i, double dt, double drho)
{
    // Материальные константы (ранее — глобал param.hpp): const_JH_2.
    const auto& const_JH_2 = mats.const_JH_2;
    size_t nmat = static_cast<size_t>(sys.itype[i] - 10);
    double gem = const_JH_2[nmat][1];
    double P_hel = const_JH_2[nmat][2];
    double T_n = const_JH_2[nmat][4];
    double s_hel = const_JH_2[nmat][3];
    double A = const_JH_2[nmat][5];
    double N = const_JH_2[nmat][6];
    double B = const_JH_2[nmat][7];
    double M = const_JH_2[nmat][8];
    double s_f_max = const_JH_2[nmat][9];
    double C = const_JH_2[nmat][10];
    double D1 = const_JH_2[nmat][11];
    double D2 = const_JH_2[nmat][12];
    double K1 = const_JH_2[nmat][13];
    double beta_JH = 1.0;
    double P_n = sys.p[i] / P_hel;
    double eps_i = (sqrt(2) * sqrt((sys.exx[i] - sys.eyy[i]) * (sys.exx[i] - sys.eyy[i]) + (sys.eyy[i] - sys.ezz[i]) * (sys.eyy[i] - sys.ezz[i]) +
                                   (sys.ezz[i] - sys.exx[i]) * (sys.ezz[i] - sys.exx[i]) + 6 * (sys.exy[i] * sys.exy[i] + sys.exz[i] * sys.exz[i] + sys.eyz[i] * sys.eyz[i]))) / 3;
    double sigt0 = 0.0;
    double sigtf = 0.0;
    if (P_n <= 0)
    {
        if ((P_n + T_n) > 0)
        {
            sigt0 = s_hel * A * pow(P_n + T_n, N);
            // sigtf = s_hel*B*pow(fabs(P_n), M);
        }
    }
    else
    {
        sigt0 = s_hel * A * pow(P_n + T_n, N);
        sigtf = s_hel * B * pow(P_n, M);
    }
    if (sigtf > s_f_max * s_hel)
    {
        sigtf = s_f_max * s_hel;
    }
    double sigt = (sigt0 - sys.D[i] * (sigt0 - sigtf));
    if (eps_i > 0)
    {
        sigt *= (1 + C * log(eps_i));
    }
    double smiz2m = pow(sigt, 2) / 3.0;
    double smizm = sqrt(pow(sigt, 2) / 3.0);
    double sxx = sys.sxx[i];
    double syy = sys.syy[i];
    double sxy = sys.sxy[i];
    double sxz = sys.sxz[i];
    double syz = sys.syz[i];
    // Пробное (до радиального возврата) напряжение: sys.sxx..sxz -- девиатор.
    // ripl = sqrt(J2_пробное), seq_tr = sqrt(3*J2_пробное) -- эквивалентное
    // напряжение (в 1D seq = 1.5*|sxx|); sigt -- предел текучести в тех же
    // единицах (sigma_eq), возврат при seq_tr > sigt (sqrt(J2) > sigt/sqrt(3)).
    double ripl = sqrt(sigma_i(sxx, syy, sxy, syz, sxz));
    const double seq_tr = sqrt(3.0) * ripl;
    const double D_old = sys.D[i];
    // drho (скорость изменения плотности) JH-2 больше не использует:
    // приращение пластической деформации находится из радиального возврата.
    (void)drho;
    (void)dt;
    bool flag_e = true;
    if (ripl > smizm)
    {
        double corr = smizm / ripl;
        sys.sxx[i] *= corr;
        sys.syy[i] *= corr;
        sys.szz[i] *= corr;
        sys.sxy[i] *= corr;
        sys.syz[i] *= corr;
        sys.sxz[i] *= corr;
        flag_e = false;
    }
    double mu = (sys.rho[i] / sys.rho_0[i]) - 1;
    if (mu >= 0)
    {
        double u0 = (2*ripl) / (6 * gem);
        double u1 = (2*sigma_i(sys.sxx[i], sys.syy[i], sys.sxy[i], sys.syz[i], sys.sxz[i])) / (6 * gem);
        double du = u0 - u1;
        if (du > 0)
        {
            double ds_1 = pow((K1 * mu + sys.dp[i]), 2);
            double ds_2 = 2 * beta_JH * K1 * du;
            if (fabs(ds_1) >= fabs(ds_2))
            {
                sys.dp[i] = -K1 * mu + sqrt(ds_1 + ds_2);
            }
        }
    }
    if (flag_e == false && D_old < 1)
    {
        // Радиальный возврат (Wilkins): приращение эквивалентной пластической
        // деформации за шаг d_ep = (sigma_eq_пробное - sigma_y)/(3G), без
        // тензорных компонент и вращательных членов (jh2/jh2-code-review.md, п. 1).
        double d_ep = (seq_tr - sigt) / (3.0 * gem);
        if ((P_n + T_n) > 0)
        {
            double znam = D1 * pow(P_n + T_n, D2);
            sys.e_pl[i] += d_ep;
            sys.D[i] += d_ep / znam;
            if (sys.D[i] > 1)
            {
                sys.D[i] = 1.0;
            }
        }
        else
        {
            sys.D[i] = 1.0;
        }
    }
}
