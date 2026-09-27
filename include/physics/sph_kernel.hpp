#pragma once
// physics/sph_kernel.hpp
// Ядро сглаживания (W) и его градиент (dW/dx) для skf = 1 (кубический
// сплайн, Monaghan 1985), 2 (Гаусс) и 3 (квинтическое, Morris 1997) --
// fiber/CUDA_PLAN.md §3.2, этап Г2: "функции ядра сглаживания (W, ∇W)".
//
// Построчная копия арифметики kernel.cpp (Ф0-Ф8): порядок операций и типы
// промежуточных величин не менялись -- требование побитной регрессии (план
// §3.2). Единственное содержательное отличие от kernel.cpp -- сигнатура
// (const real_t* dx вместо std::span<const real_t>, чтобы функция была
// SPH_HD: std::span как параметр физической функции не нужен, а сырой
// указатель -- тот же приём, что и у SolidView/DeviceBuffer в плане §3.2)
// и отсутствие исключения на недопустимой размерности: dim -- constexpr
// (param.hpp, SPH_DIM), задаётся при сборке и всегда равен 1, 2 или 3 для
// всех целей проекта, поэтому ветка "иначе" в kernel.cpp была недостижима
// уже в исходном коде -- её отсутствие здесь не меняет поведение ни для
// одной существующей сборки, а бросать исключение из SPH_HD-функции нельзя
// (план §3.2, "Правила кода").
//
// kernel.cpp/kernel.hpp сохраняют старую сигнатуру (std::span, есть другие
// потребители в будущих частях Г2) и делегируют сюда без изменения
// арифметики.
#include <cmath>
#include <cstddef>
#include "param.hpp"       // dim, pi
#include "physics/hd.hpp"  // SPH_HD

namespace physics
{

// Вычисление ядра W и его производной dW/dx.
// skf -- тип ядра; r -- расстояние; dx -- компоненты вектора x_ij (длина
// dim); mhsml -- средняя длина сглаживания; w -- значение ядра; dwdx[d] --
// градиент ядра (длина dim).
SPH_HD inline void sph_kernel_w_dwdx(int skf, real_t r, const real_t* dx, real_t mhsml, real_t& w, real_t* dwdx)
{
    // Инициализация выходных значений
    w = 0;

    for (std::size_t d = 0; d < dim; d++)
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

        // Вычисление ядра и его производной
        if (q >= 0 && q <= real_t(1.0))
        {
            w = factor * (real_t(2.0)/real_t(3.0) - q*q + q*q*q / real_t(2.0));
            for (std::size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = factor * (-real_t(2.0) + real_t(1.5) * q) * dx[d] / (mhsml * mhsml);
            }
        }
        else if (q > real_t(1.0) && q <= real_t(2.0))
        {
            w = factor * (real_t(1.0)/real_t(6.0)) * ((real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q));
            for (std::size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = -factor * real_t(0.5) * ((real_t(2.0) - q)*(real_t(2.0) - q)) * dx[d] / (r * mhsml);
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
            for (std::size_t d = 0; d < dim; ++d)
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

        if (q >= 0 && q <= real_t(1.0))
        {
            w = factor * (((real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)) - real_t(6.0) * ((real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)) + real_t(15.0) * ((real_t(1.0) - q)*(real_t(1.0) - q)*(real_t(1.0) - q)*(real_t(1.0) - q)*(real_t(1.0) - q)));
            for (std::size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = factor * (-real_t(120.0) + real_t(120.0)*q - real_t(50.0)*q*q) * dx[d] / (mhsml * mhsml);
            }
        }
        else if (q > real_t(1.0) && q <= real_t(2.0))
        {
            w = factor * (((real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)*(real_t(3.0) - q)) - real_t(6.0) * ((real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)*(real_t(2.0) - q)));
            for (std::size_t d = 0; d < dim; ++d)
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
            for (std::size_t d = 0; d < dim; ++d)
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

} // namespace physics
