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
void prochnost_mises(size_t nmat, double &sxx, double &syy, double &szz, double &sxy, double &syz, double &sxz)
{
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
void prochnost_cem_JH2(Particle &Elem, double dt, double drho)
{
    size_t nmat = Elem.itype - 10;
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
    double P_n = Elem.p / P_hel;
    double eps_i = (sqrt(2) * sqrt((Elem.exx - Elem.eyy) * (Elem.exx - Elem.eyy) + (Elem.eyy - Elem.ezz) * (Elem.eyy - Elem.ezz) +
                                   (Elem.ezz - Elem.exx) * (Elem.ezz - Elem.exx) + 6 * (Elem.exy * Elem.exy + Elem.exz * Elem.exz + Elem.eyz * Elem.eyz))) / 3;
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
    double sigt = (sigt0 - Elem.D * (sigt0 - sigtf));
    if (eps_i > 0)
    {
        sigt *= (1 + C * log(eps_i));
    }
    double smiz2m = pow(sigt, 2) / 3.0;
    double smizm = sqrt(pow(sigt, 2) / 3.0);
    double sxx = Elem.sxx;
    double syy = Elem.syy;
    double szz = Elem.szz;
    double sxy = Elem.sxy;
    double sxz = Elem.sxz;
    double syz = Elem.syz;
    double ripl = sqrt(sigma_i(sxx, syy, sxy, syz, sxz));
    bool flag_e = true;
    if (ripl > smizm)
    {
        double corr = smizm / ripl;
        Elem.sxx *= corr;
        Elem.syy *= corr;
        Elem.szz *= corr;
        Elem.sxy *= corr;
        Elem.syz *= corr;
        Elem.sxz *= corr;
        flag_e = false;
    }
    double mu = (Elem.rho / Elem.rho_0) - 1;
    if (mu >= 0)
    {
        double u0 = (2*ripl) / (6 * gem);
        double u1 = (2*sigma_i(Elem.sxx, Elem.syy, Elem.sxy, Elem.syz, Elem.sxz)) / (6 * gem);
        double du = u0 - u1;
        if (du > 0)
        {
            double ds_1 = pow((K1 * mu + Elem.dp), 2);
            double ds_2 = 2 * beta_JH * K1 * du;
            if (fabs(ds_1) >= fabs(ds_2))
            {
                Elem.dp = -K1 * mu + sqrt(ds_1 + ds_2);
            }
        }
    }
    if ((flag_e == false) && (Elem.D < 1))
    {
        double dsx = (Elem.sxx - sxx) / dt;
        double dsy = (Elem.syy - syy) / dt;
        double dsz = (Elem.szz - szz) / dt;
        double dsxy = (Elem.sxy - sxy) / dt;
        double dsyz = (Elem.syz - syz) / dt;
        double dsxz = (Elem.sxz - sxz) / dt;
        double eps_e_x = 0.0, eps_e_y = 0.0, eps_e_z = 0.0, eps_e_xy = 0.0, eps_e_yz = 0.0, eps_e_xz = 0.0;
        //double eps_vol = (Elem.exx + Elem.eyy + Elem.ezz) / 3.0;
        double eps_vol = drho / (3.0 * Elem.rho);
        if (dim > 1)
        {
            eps_e_x = ((dsx - 2 * (Elem.sxy * Elem.rxy + Elem.sxz * Elem.rxz)) / gem) + eps_vol;
            eps_e_y = ((dsy - 2 * (Elem.syz * Elem.ryz - Elem.sxy * Elem.rxy)) / gem) + eps_vol;
            eps_e_z = ((dsz - (-2 * (Elem.sxz * Elem.rxz - Elem.syz * Elem.ryz))) / gem) + eps_vol;
            eps_e_xy = (dsxy - (Elem.rxy * (Elem.syy - Elem.sxx) + Elem.sxz * Elem.ryz + Elem.syz * Elem.rxz)) / gem;
            eps_e_yz = (dsyz - (-Elem.sxy * Elem.rxz - Elem.sxz * Elem.rxy + Elem.ryz * (Elem.szz - Elem.syy))) / gem;
            eps_e_xz = (dsxz - (Elem.rxz * (Elem.szz - Elem.sxx) - Elem.sxy * Elem.ryz + Elem.syz * Elem.rxy)) / gem;
        }
        else if (dim == 1)
        {
            eps_e_x = (dsx / gem) + eps_vol;
        }
        double eps_p_x = Elem.exx - eps_e_x;
        double eps_p_y = Elem.eyy - eps_e_y;
        double eps_p_z = Elem.ezz - eps_e_z;
        double eps_p_xy = Elem.exy - eps_e_xy;
        double eps_p_yz = Elem.eyz - eps_e_yz;
        double eps_p_xz = Elem.exz - eps_e_xz;
        double eps_i_f = (sqrt(2) * sqrt((eps_p_x - eps_p_y) * (eps_p_x - eps_p_y) + (eps_p_y - eps_p_z) * (eps_p_y - eps_p_z) +
                                         (eps_p_z - eps_p_x) * (eps_p_z - eps_p_x) + 6 * (eps_p_xy * eps_p_xy + eps_p_xz * eps_p_xz + eps_p_yz * eps_p_yz))) / 3;
        if ((P_n + T_n) > 0)
        {
            double znam = D1 * pow(P_n + T_n, D2);
            Elem.e_pl += eps_i_f * dt;
            Elem.D += ((eps_i_f * dt) / znam);
            if (Elem.D > 1)
            {
                Elem.D = 1.0;
            }
        }
        else
        {
            Elem.D = 1.0;
        }
    }
}