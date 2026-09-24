// ppochnost.cpp
// Определения моделей прочности (Мизес, JH-2) и интенсивности напряжений.
// Объявления — в ppochnost.hpp.
#include "ppochnost.hpp"

double sigma_i(double &sxx, double &syy, double &sxy, double &syz, double &sxz)
{
    if (dim == 1) // одномерная деформация
    {
        return 0.75 * sxx * sxx;
        // return (sxx * sxx + syy * syy + sxx * syy);
        // return (sxx - syy) * (sxx - syy) / 3.0;
    }
    else if (dim == 2)
    {
        
        if (axis_sim == true)
        {
            return (sxx * sxx + sxx * syy + syy * syy + sxy * sxy); // осесимметричный случай
        }
        else
        {
            return (sxx * sxx + sxy * sxy); // плоский случай
        }
    }
    else if (dim == 3) // трехмерный случай
    {
        return (sxx * sxx + sxx * syy + syy * syy + sxy * sxy + sxz * sxz + syz * syz);
    }
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
    double szz = sys.szz[i];
    double sxy = sys.sxy[i];
    double sxz = sys.sxz[i];
    double syz = sys.syz[i];
    double ripl = sqrt(sigma_i(sxx, syy, sxy, syz, sxz));
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
    if ((flag_e == false) && (sys.D[i] < 1))
    {
        double dsx = (sys.sxx[i] - sxx) / dt;
        double dsy = (sys.syy[i] - syy) / dt;
        double dsz = (sys.szz[i] - szz) / dt;
        double dsxy = (sys.sxy[i] - sxy) / dt;
        double dsyz = (sys.syz[i] - syz) / dt;
        double dsxz = (sys.sxz[i] - sxz) / dt;
        double eps_e_x = 0.0, eps_e_y = 0.0, eps_e_z = 0.0, eps_e_xy = 0.0, eps_e_yz = 0.0, eps_e_xz = 0.0;
        //double eps_vol = (sys.exx[i] + sys.eyy[i] + sys.ezz[i]) / 3.0;
        double eps_vol = drho / (3.0 * sys.rho[i]);
        if (dim > 1)
        {
            eps_e_x = ((dsx - 2 * (sys.sxy[i] * sys.rxy[i] + sys.sxz[i] * sys.rxz[i])) / gem) + eps_vol;
            eps_e_y = ((dsy - 2 * (sys.syz[i] * sys.ryz[i] - sys.sxy[i] * sys.rxy[i])) / gem) + eps_vol;
            // Поворотная часть — та же, что у dszz в int_force (−2(sxz·rxz +
            // syz·ryz)); знак при syz·ryz исправлен вместе с ней.
            eps_e_z = ((dsz - (-2 * (sys.sxz[i] * sys.rxz[i] + sys.syz[i] * sys.ryz[i]))) / gem) + eps_vol;
            eps_e_xy = (dsxy - (sys.rxy[i] * (sys.syy[i] - sys.sxx[i]) + sys.sxz[i] * sys.ryz[i] + sys.syz[i] * sys.rxz[i])) / gem;
            eps_e_yz = (dsyz - (-sys.sxy[i] * sys.rxz[i] - sys.sxz[i] * sys.rxy[i] + sys.ryz[i] * (sys.szz[i] - sys.syy[i]))) / gem;
            eps_e_xz = (dsxz - (sys.rxz[i] * (sys.szz[i] - sys.sxx[i]) - sys.sxy[i] * sys.ryz[i] + sys.syz[i] * sys.rxy[i])) / gem;
        }
        else if (dim == 1)
        {
            eps_e_x = (dsx / gem) + eps_vol;
        }
        double eps_p_x = sys.exx[i] - eps_e_x;
        double eps_p_y = sys.eyy[i] - eps_e_y;
        double eps_p_z = sys.ezz[i] - eps_e_z;
        double eps_p_xy = sys.exy[i] - eps_e_xy;
        double eps_p_yz = sys.eyz[i] - eps_e_yz;
        double eps_p_xz = sys.exz[i] - eps_e_xz;
        double eps_i_f = (sqrt(2) * sqrt((eps_p_x - eps_p_y) * (eps_p_x - eps_p_y) + (eps_p_y - eps_p_z) * (eps_p_y - eps_p_z) +
                                         (eps_p_z - eps_p_x) * (eps_p_z - eps_p_x) + 6 * (eps_p_xy * eps_p_xy + eps_p_xz * eps_p_xz + eps_p_yz * eps_p_yz))) / 3;
        if ((P_n + T_n) > 0)
        {
            double znam = D1 * pow(P_n + T_n, D2);
            sys.e_pl[i] += eps_i_f * dt;
            sys.D[i] += ((eps_i_f * dt) / znam);
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
