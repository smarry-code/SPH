// kernel.cpp

// Функция вычисления ядра и его производной
void kernel(double r, const std::vector<double>& dx, double mhsml, double& w, double* dwdx)
{
    // Инициализация выходных значений
    w = 0.0;

    for (size_t d = 0; d < dim; d++)
    {
        dwdx[d] = 0.0;
    }
    // Вычисление относительного расстояния
    double q = r / mhsml;
    double factor = 0.0;

    // Выбор ядра в зависимости от skf
    if (skf == 1)    // Кубический сплайн (Monaghan 1985)
    {
        // Определение нормировочного фактора в зависимости от размерности
        if (dim == 1)
        {
            factor = 1.0 / mhsml;
        }
        else if (dim == 2)
        {
            factor = 15.0 / (7.0 * pi * mhsml * mhsml);
        }
        else if (dim == 3)
        {
            factor = 3.0 / (2.0 * pi * mhsml * mhsml * mhsml);
        }
        else
        {
            throw std::runtime_error("Неверная размерность: dim = " + std::to_string(dim));
        }

        // Вычисление ядра и его производной
        if (q >= 0.0 && q <= 1.0)
        {
            w = factor * (2.0/3.0 - q*q + q*q*q / 2.0);
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = factor * (-2.0 + 1.5 * q) * dx[d] / (mhsml * mhsml);
            }
        }
        else if (q > 1.0 && q <= 2.0)
        {
            double t = 2.0 - q;
            w = factor * (1.0/6.0) * std::pow(2.0 - q, 3);
            //w = factor * (1.0/6.0) * t * t * t;
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = -factor * 0.5 * std::pow(2.0 - q, 2) * dx[d] / (r * mhsml);
                //dwdx[d] = -factor * 0.5 * t * t * dx[d] / (r * mhsml);
            }
        }
        // Иначе w и dwdx остаются нулевыми

    }
    else if (skf == 2)      // Гауссово ядро
    {
        factor = 1.0 / (std::pow(mhsml, dim) * std::pow(pi, dim/2.0));

        if (q >= 0.0 && q <= 3.0)
        {
            w = factor * std::exp(-q*q);
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = w * (-2.0 * dx[d]) / (mhsml * mhsml);
            }
        }
        // Иначе w и dwdx остаются нулевыми

    }
    else if (skf == 3)      // Квинтическое ядро (Morris 1997)
    {
        if (dim == 1)
        {
            factor = 1.0 / (120.0 * mhsml);
        }
        else if (dim == 2)
        {
            factor = 7.0 / (478.0 * pi * mhsml * mhsml);
        }
        else if (dim == 3)
        {
            factor = 1.0 / (120.0 * pi * mhsml * mhsml * mhsml);
        }
        else
        {
            throw std::runtime_error("Неверная размерность: dim = " + std::to_string(dim));
        }

        if (q >= 0.0 && q <= 1.0)
        {
            w = factor * (std::pow(3.0 - q, 5) - 6.0 * std::pow(2.0 - q, 5) + 15.0 * std::pow(1.0 - q, 5));
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = factor * (-120.0 + 120.0*q - 50.0*q*q) * dx[d] / (mhsml * mhsml);
            }
        }
        else if (q > 1.0 && q <= 2.0)
        {
            w = factor * (std::pow(3.0 - q, 5) - 6.0 * std::pow(2.0 - q, 5));
            for (size_t d = 0; d < dim; ++d)
            {
                if (r > 1e-8)
                {
                    dwdx[d] = factor * (-5.0 * std::pow(3.0 - q, 4) + 30.0 * std::pow(2.0 - q, 4)) * dx[d] / (r * mhsml);
                }
            }
        }
        else if (q > 2.0 && q <= 3.0)
        {
            w = factor * std::pow(3.0 - q, 5);
            for (size_t d = 0; d < dim; ++d)
            {
                if (r > 1e-8)
                {
                    dwdx[d] = -factor * 5.0 * std::pow(3.0 - q, 4) * dx[d] / (r * mhsml);
                }
            }
        }
        // Иначе w и dwdx остаются нулевыми
    }
}
