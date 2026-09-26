// kernel.cpp
// Определение функции вычисления ядра сглаживания и его производной.
// Объявление — в kernel.hpp.
#include "kernel.hpp"

void kernel(int skf, real_t r, std::span<const real_t> dx, real_t mhsml, real_t& w, real_t* dwdx)
{
    // Инициализация выходных значений
    w = 0;

    for (size_t d = 0; d < dim; d++)
    {
        dwdx[d] = 0;
    }
    // Вычисление относительного расстояния
    real_t q = r / mhsml;
    real_t factor = 0;

    // Выбор ядра в зависимости от skf
    if (skf == 1)    // Кубический сплайн (Monaghan 1985)
    {
        // Определение нормировочного фактора в зависимости от размерности
        if (dim == 1)
        {
            factor = real_t(1.0) / mhsml;
        }
        else if (dim == 2)
        {
            factor = real_t(15.0) / (real_t(7.0) * pi * mhsml * mhsml);
        }
        else if (dim == 3)
        {
            factor = real_t(3.0) / (real_t(2.0) * pi * mhsml * mhsml * mhsml);
        }
        else
        {
            throw std::runtime_error("Неверная размерность: dim = " + std::to_string(dim));
        }

        // Вычисление ядра и его производной
        if (q >= 0 && q <= real_t(1.0))
        {
            w = factor * (real_t(2.0)/real_t(3.0) - q*q + q*q*q / real_t(2.0));
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = factor * (-real_t(2.0) + real_t(1.5) * q) * dx[d] / (mhsml * mhsml);
            }
        }
        else if (q > real_t(1.0) && q <= real_t(2.0))
        {
            real_t t = real_t(2.0) - q;
            w = factor * (real_t(1.0)/real_t(6.0)) * ((real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q));
            //w = factor * (1.0/6.0) * t * t * t;
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = -factor * real_t(0.5) * ((real_t(2.0) - q)*(real_t(2.0) - q)) * dx[d] / (r * mhsml);
                //dwdx[d] = -factor * 0.5 * t * t * dx[d] / (r * mhsml);
            }
        }
        // Иначе w и dwdx остаются нулевыми

    }
    else if (skf == 2)      // Гауссово ядро
    {
        factor = real_t(1.0) / (std::pow(mhsml, static_cast<real_t>(dim)) * std::pow(pi, dim/real_t(2.0)));

        if (q >= 0 && q <= real_t(3.0))
        {
            w = factor * std::exp(-q*q);
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = w * (-real_t(2.0) * dx[d]) / (mhsml * mhsml);
            }
        }
        // Иначе w и dwdx остаются нулевыми

    }
    else if (skf == 3)      // Квинтическое ядро (Morris 1997)
    {
        if (dim == 1)
        {
            factor = real_t(1.0) / (real_t(120.0) * mhsml);
        }
        else if (dim == 2)
        {
            factor = real_t(7.0) / (real_t(478.0) * pi * mhsml * mhsml);
        }
        else if (dim == 3)
        {
            factor = real_t(1.0) / (real_t(120.0) * pi * mhsml * mhsml * mhsml);
        }
        else
        {
            throw std::runtime_error("Неверная размерность: dim = " + std::to_string(dim));
        }

        if (q >= 0 && q <= real_t(1.0))
        {
            w = factor * (((real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)) - real_t(6.0) * ((real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)) + real_t(15.0) * ((real_t(1.0) - q)*(real_t(1.0) - q)*(real_t(1.0) - q)*(real_t(1.0) - q)*(real_t(1.0) - q)));
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = factor * (-real_t(120.0) + real_t(120.0)*q - real_t(50.0)*q*q) * dx[d] / (mhsml * mhsml);
            }
        }
        else if (q > real_t(1.0) && q <= real_t(2.0))
        {
            w = factor * (((real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)) - real_t(6.0) * ((real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)));
            for (size_t d = 0; d < dim; ++d)
            {
                if (r > real_t(1e-8))
                {
                    dwdx[d] = factor * (-real_t(5.0) * ((real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)) + real_t(30.0) * ((real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q))) * dx[d] / (r * mhsml);
                }
            }
        }
        else if (q > real_t(2.0) && q <= real_t(3.0))
        {
            w = factor * ((real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q));
            for (size_t d = 0; d < dim; ++d)
            {
                if (r > real_t(1e-8))
                {
                    dwdx[d] = -factor * real_t(5.0) * ((real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)) * dx[d] / (r * mhsml);
                }
            }
        }
        // Иначе w и dwdx остаются нулевыми
    }
}
