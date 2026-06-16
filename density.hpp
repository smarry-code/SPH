// density.cpp

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

void eos_mi_grunaizen(size_t nmat, double rho, double rho_0, double u, double& p, double& c)
{
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
void eos_ua_teta(size_t nmat, double rho, double rho_0, double &p, double &c)
{
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
void urs_cem_JH2(size_t nmat, double qo, double qo0, double dp, double &as, double &pv, double &D)
{
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
// Функция вычисления плотности через SPH-суммирование
void sum_density(std::vector<Particle>& Elem, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j, const std::vector<double>& w)
{
    // Инициализация временных массивов
    std::vector<double> wi(maxn, 0.0);
    std::vector<double> hv(dim, 0.0);
    double r = 0.0;

    // Вычисление самоплотности для каждой частицы
    for (size_t i = 0; i < ntotal; ++i)
    {
        double selfdens = 0.0;
        double tdwdx[dim] = {};
        kernel(r, hv, Elem[i].hsml, selfdens, tdwdx);
        wi[i] = selfdens * Elem[i].mass / Elem[i].rho;
    }



    // Расчет нормализационного знаменателя
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k] - 1; // Переводим индексы в C++ (начинаются с 0)
        size_t j = pair_j[k] - 1;

        wi[i] += Elem[j].mass / Elem[j].rho * w[k];
        wi[j] += Elem[i].mass / Elem[i].rho * w[k];
    }

    // Инициализация плотности
    for (size_t i = 0; i < ntotal; ++i)
    {
        Elem[i].rho = 0.0;
    }
    // Добавление самовклада
    for (size_t i = 0; i < ntotal; ++i)
    {
        double selfdens = 0.0;
        double tdwdx[dim] = {};
        kernel(r, hv, Elem[i].hsml, selfdens, tdwdx);
        Elem[i].rho += selfdens * Elem[i].mass;
    }

    // Обновление плотности
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k] - 1;
        size_t j = pair_j[k] - 1;

        Elem[i].rho += Elem[j].mass * w[k];
        Elem[j].rho += Elem[i].mass * w[k];
    }

    // Нормализация плотности
    if (nor_density)
    {
        for (size_t i = 0; i < ntotal; ++i)
        {
            if (wi[i] > 1e-10)   // Защита от деления на ноль
            {
                Elem[i].rho /= wi[i];
            }
        }
    }
}

// Функция вычисления изменения плотности через уравнение непрерывности
void con_density(std::vector<Particle>& Elem, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
                 const std::vector<std::vector<double>>& dwdx, std::vector<double>& drhodt)
{
    // Инициализация производной плотности
    std::fill(drhodt.begin(), drhodt.begin() + ntotal, 0.0);
    /*for (size_t i = 0; i < ntotal; ++i)
    {
        // Вычисление давления через уравнение состояния
        if (std::abs(itype[i]) == 1)
        {
            p_gas(rho[i], u[i], p[i], c[i]);
        }
        else if (std::abs(itype[i]) == 2)
        {
            p_art_water(rho[i], p[i], c[i]);
        }
        else if (std::abs(itype[i]) >= 3)
        {
            eos_mi_grunaizen(itype[i], rho[i], rho_0[i], u[i], p[i], c[i]);
        }
    } */
    // Обработка всех пар взаимодействий
    #pragma omp parallel for schedule(auto)
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k] - 1;
        size_t j = pair_j[k] - 1;
        std::vector<double> dvx(dim);
        for (size_t d = 0; d < dim; ++d)
        {
            dvx[d] = Elem[i].vx[d] - Elem[j].vx[d];
            if (Elem[i].itype != Elem[j].itype)
            {
                dvx[d] = 2 * (Elem[i].vx[d] - (Elem[i].vx[d] + Elem[j].vx[d]) * 0.5);
            }
        }
        double vcc = 0.0;
        for (size_t d = 0; d < dim; ++d)
        {
            vcc += dvx[d] * dwdx[d][k];
        }

        // drhodt[i] += mass[j] * vcc;
        // drhodt[j] += mass[i] * vcc;
        double q_i = (Elem[i].rho / Elem[j].rho) * Elem[j].mass * vcc;
        double q_j = (Elem[j].rho / Elem[i].rho) * Elem[i].mass * vcc;
#pragma omp atomic update
        drhodt[i] += q_i;
#pragma omp atomic update
        drhodt[j] += q_j;
        if (axis_sim == true && dim == 2)
        {
            // drhodt[i] += (Elem[i].rho / Elem[j].rho) * Elem[j].mass * (Elem[j].vx[0] /  Elem[j].x[0]) * w[k];
            // drhodt[j] += (Elem[j].rho / Elem[i].rho) * Elem[i].mass * (Elem[i].vx[0] /  Elem[i].x[0]) * w[k];
            double ax_i = (Elem[i].vx[0] / Elem[i].x[0]);
            double ax_j = (Elem[j].vx[0] / Elem[j].x[0]);
#pragma omp atomic update
            drhodt[i] -= ax_i;
#pragma omp atomic update
            drhodt[j] -= ax_j;
        }
    }
}

