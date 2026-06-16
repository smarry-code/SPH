// ext_force.cpp
void virtual_symm_2d (std::vector<Particle>& Elem, size_t i, size_t j, std::vector<std::vector<double>>& dvxdt)
{
    constexpr double rr0 = 1.25e-5;  // Пороговое расстояние
    constexpr double dd = 1.0e-2;    // Сила отталкивания
    constexpr int p1 = 12;           // Степень для отталкивающей силы
    constexpr int p2 = 4;            // Степень для притягивающей силы
    if (axis_sim == true && dim == 2)
    {
        if (Elem[i].hsml >= Elem[i].x[0])
        {
            Particle buf_elem = Elem[j];
            buf_elem.x[0] = -Elem[j].x[0];
            buf_elem.vx[0] = -Elem[j].vx[0];
            buf_elem.x[0] = -Elem[j].x[0];
            buf_elem.sxy = -Elem[j].sxy;
            std::vector<double> dx(dim);
            double rr = 0.0;
            for (size_t d = 0; d < dim; ++d)
            {
                dx[d] = Elem[i].x[d] - buf_elem.x[d];
                rr += dx[d] * dx[d];
            }
            rr = std::sqrt(rr);
            if (rr < rr0)
            {
                double inv_rr = 1.0 / rr;
                double force = ((std::pow(rr0 * inv_rr, p1)) - (std::pow(rr0 * inv_rr, p2))) * inv_rr * inv_rr;

                // Применение силы
                for (size_t d = 0; d < dim; ++d)
                {
                    double f = dd * dx[d] * force;
                    dvxdt[d][i] += f;
                    dvxdt[d][j] -= f;
                }
            }
        }
    }
}

// Основная функция ext_force
void ext_force(std::vector<Particle>& Elem, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j, std::vector<std::vector<double>>& dvxdt)
{

    // Параметры для силы отталкивания
    constexpr double rr0 = 1.25e-5;  // Пороговое расстояние
    constexpr double dd = 1.0e-2;    // Сила отталкивания
    constexpr int p1 = 12;           // Степень для отталкивающей силы
    constexpr int p2 = 4;            // Степень для притягивающей силы

    // Обнуление массива ускорений
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            dvxdt[d][i] = 0.0;
        }
    }

    // Учёт самогравитации
    if (self_gravity)
    {
        // Гравитация действует в последней размерности (например, z в 3D)
        size_t gravity_dim = dim - 1;
        for (size_t i = 0; i < ntotal; ++i)
        {
            dvxdt[gravity_dim][i] = -9.8;  // Ускорение свободного падения
        }
    }

    // Обработка всех пар взаимодействий
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k] - 1;  // Переводим в C++ индексацию
        size_t j = pair_j[k] - 1;

        // Проверка: одна частица реальная (itype > 0), другая виртуальная (itype < 0)
        bool is_real_virtual_pair = (Elem[i].itype > 0 && Elem[j].itype < 0);
        //bool is_contact_pair_ = (itype[i] != itype[j]);
        if (is_real_virtual_pair)
        {
            // Вычисление разности координат
            std::vector<double> dx(dim);
            double rr = 0.0;

            for (size_t d = 0; d < dim; ++d)
            {
                dx[d] = Elem[i].x[d] - Elem[j].x[d];
                rr += dx[d] * dx[d];
            }

            rr = std::sqrt(rr);

            // Проверка на близость частиц
            if (rr < rr0)
            {
                double inv_rr = 1.0 / rr;
                double force = ((std::pow(rr0 * inv_rr, p1)) - (std::pow(rr0 * inv_rr, p2))) * inv_rr * inv_rr;

                // Применение силы
                for (size_t d = 0; d < dim; ++d)
                {
                    double f = dd * dx[d] * force;
                    dvxdt[d][i] += f;
                    dvxdt[d][j] -= f;
                }
            }
        }
        double A = 1e5;
        double n1 = 6;
        double n2 = 4;
        bool is_contact_pair_ = false;
        if (is_contact_pair_)
        {
            // Вычисление разности координат
            std::vector<double> dx(dim);
            double rr = 0.0;

            for (size_t d = 0; d < dim; ++d)
            {
                dx[d] = Elem[i].x[d] - Elem[j].x[d];
                rr += dx[d] * dx[d];
            }
            rr = std::sqrt(rr);
            double pe = 0.5 * (Elem[i].hsml + Elem[j].hsml) / rr;
            if (pe >= 1)
            {
                double inv_rr = 1.0 / rr;
                double force = A * ((std::pow(pe, n1)) - (std::pow(pe, n2))) * inv_rr * inv_rr;
                for (size_t d = 0; d < dim; ++d)
                {
                    double f = dx[d] * force;
                    dvxdt[d][i] += f;
                    dvxdt[d][j] -= f;
                }
            }
        }
    }

    // Проверка на NaN и бесконечность
    bool has_nan_or_inf = false;

    // Проверка массива dvxdt
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            if (std::isnan(dvxdt[d][i]) || std::isinf(dvxdt[d][i]))
            {
                std::cerr << "ERROR: dvxdt[" << d << "][" << i << "] = "
                          << dvxdt[d][i] << " (NaN or Inf detected)" << std::endl;
                has_nan_or_inf = true;
            }
        }
    }

    // Если обнаружены ошибки, выбрасываем исключение
    if (has_nan_or_inf)
    {
        throw std::runtime_error("NaN or Inf detected in ext_force");
    }
}
