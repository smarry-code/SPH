// ppochnost.cpp
// Определения моделей прочности (Мизес, JH-2) и интенсивности напряжений.
// Объявления — в ppochnost.hpp.
#include "ppochnost.hpp"

real_t sigma_i(real_t &sxx, real_t &syy, real_t &sxy, real_t &syz, real_t &sxz)
{
    if (dim == 1) // одномерная деформация
    {
        return real_t(0.75) * sxx * sxx;
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
void prochnost_mises(const Materials& mats, size_t nmat, real_t &sxx, real_t &syy, real_t &szz, real_t &sxy, real_t &syz, real_t &sxz)
{
    // Материальные константы (ранее — глобалы param.hpp): G, sigy, n_Me.
    const auto& G = mats.G;
    const auto& sigy = mats.sigy;
    const size_t n_Me = mats.n_Me;
    real_t gem = G[nmat - n_Me];
    real_t ge2m = real_t(2.0) * G[nmat - n_Me];
    real_t smiz2m = std::pow(sigy[nmat - n_Me], real_t(2)) / real_t(3.0);
    real_t smizm = sqrt(std::pow(sigy[nmat - n_Me], real_t(2)) / real_t(3.0));
    real_t ripl = sigma_i(sxx, syy, sxy, syz, sxz);
    bool flag_e = true;
    if (ripl > smiz2m)
    {
        real_t corr = smizm / sqrt(ripl);
        sxx *= corr;
        syy *= corr;
        szz *= corr;
        sxy *= corr;
        syz *= corr;
        sxz *= corr;
        flag_e = false;
    }
}
void prochnost_cem_JH2(const Materials& mats, ParticleSystem& sys, size_t i, real_t dt, real_t drho)
{
    // Материальные константы (ранее — глобал param.hpp): const_JH_2.
    const auto& const_JH_2 = mats.const_JH_2;
    size_t nmat = static_cast<size_t>(sys.itype[i] - 10);
    real_t gem = const_JH_2[nmat][1];
    real_t P_hel = const_JH_2[nmat][2];
    real_t T_n = const_JH_2[nmat][4];
    real_t s_hel = const_JH_2[nmat][3];
    real_t A = const_JH_2[nmat][5];
    real_t N = const_JH_2[nmat][6];
    real_t B = const_JH_2[nmat][7];
    real_t M = const_JH_2[nmat][8];
    real_t s_f_max = const_JH_2[nmat][9];
    real_t C = const_JH_2[nmat][10];
    real_t D1 = const_JH_2[nmat][11];
    real_t D2 = const_JH_2[nmat][12];
    real_t K1 = const_JH_2[nmat][13];
    real_t beta_JH = real_t(1.0);
    real_t P_n = sys.p[i] / P_hel;
    real_t eps_i = (std::sqrt(real_t(2)) * sqrt((sys.exx[i] - sys.eyy[i]) * (sys.exx[i] - sys.eyy[i]) + (sys.eyy[i] - sys.ezz[i]) * (sys.eyy[i] - sys.ezz[i]) +
                                   (sys.ezz[i] - sys.exx[i]) * (sys.ezz[i] - sys.exx[i]) + 6 * (sys.exy[i] * sys.exy[i] + sys.exz[i] * sys.exz[i] + sys.eyz[i] * sys.eyz[i]))) / 3;
    real_t sigt0 = 0;
    real_t sigtf = 0;
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
    real_t sigt = (sigt0 - sys.D[i] * (sigt0 - sigtf));
    if (eps_i > 0)
    {
        sigt *= (1 + C * log(eps_i));
    }
    real_t smiz2m = std::pow(sigt, real_t(2)) / real_t(3.0);
    real_t smizm = sqrt(std::pow(sigt, real_t(2)) / real_t(3.0));
    real_t sxx = sys.sxx[i];
    real_t syy = sys.syy[i];
    real_t szz = sys.szz[i];
    real_t sxy = sys.sxy[i];
    real_t sxz = sys.sxz[i];
    real_t syz = sys.syz[i];
    real_t ripl = sqrt(sigma_i(sxx, syy, sxy, syz, sxz));
    bool flag_e = true;
    if (ripl > smizm)
    {
        real_t corr = smizm / ripl;
        sys.sxx[i] *= corr;
        sys.syy[i] *= corr;
        sys.szz[i] *= corr;
        sys.sxy[i] *= corr;
        sys.syz[i] *= corr;
        sys.sxz[i] *= corr;
        flag_e = false;
    }
    real_t mu = (sys.rho[i] / sys.rho_0[i]) - 1;
    if (mu >= 0)
    {
        real_t u0 = (2*ripl) / (6 * gem);
        real_t u1 = (2*sigma_i(sys.sxx[i], sys.syy[i], sys.sxy[i], sys.syz[i], sys.sxz[i])) / (6 * gem);
        real_t du = u0 - u1;
        if (du > 0)
        {
            real_t ds_1 = std::pow((K1 * mu + sys.dp[i]), real_t(2));
            real_t ds_2 = 2 * beta_JH * K1 * du;
            if (fabs(ds_1) >= fabs(ds_2))
            {
                sys.dp[i] = -K1 * mu + sqrt(ds_1 + ds_2);
            }
        }
    }
    if ((flag_e == false) && (sys.D[i] < 1))
    {
        real_t dsx = (sys.sxx[i] - sxx) / dt;
        real_t dsy = (sys.syy[i] - syy) / dt;
        real_t dsz = (sys.szz[i] - szz) / dt;
        real_t dsxy = (sys.sxy[i] - sxy) / dt;
        real_t dsyz = (sys.syz[i] - syz) / dt;
        real_t dsxz = (sys.sxz[i] - sxz) / dt;
        real_t eps_e_x = 0, eps_e_y = 0, eps_e_z = 0, eps_e_xy = 0, eps_e_yz = 0, eps_e_xz = 0;
        //real_t eps_vol = (sys.exx[i] + sys.eyy[i] + sys.ezz[i]) / 3.0;
        real_t eps_vol = drho / (real_t(3.0) * sys.rho[i]);
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
        real_t eps_p_x = sys.exx[i] - eps_e_x;
        real_t eps_p_y = sys.eyy[i] - eps_e_y;
        real_t eps_p_z = sys.ezz[i] - eps_e_z;
        real_t eps_p_xy = sys.exy[i] - eps_e_xy;
        real_t eps_p_yz = sys.eyz[i] - eps_e_yz;
        real_t eps_p_xz = sys.exz[i] - eps_e_xz;
        real_t eps_i_f = (std::sqrt(real_t(2)) * sqrt((eps_p_x - eps_p_y) * (eps_p_x - eps_p_y) + (eps_p_y - eps_p_z) * (eps_p_y - eps_p_z) +
                                         (eps_p_z - eps_p_x) * (eps_p_z - eps_p_x) + 6 * (eps_p_xy * eps_p_xy + eps_p_xz * eps_p_xz + eps_p_yz * eps_p_yz))) / 3;
        if ((P_n + T_n) > 0)
        {
            real_t znam = D1 * pow(P_n + T_n, D2);
            sys.e_pl[i] += eps_i_f * dt;
            sys.D[i] += ((eps_i_f * dt) / znam);
            if (sys.D[i] > 1)
            {
                sys.D[i] = real_t(1.0);
            }
        }
        else
        {
            sys.D[i] = real_t(1.0);
        }
    }
}
