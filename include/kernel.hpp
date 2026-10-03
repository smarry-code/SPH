#pragma once
// kernel.hpp
// Вычисление ядра сглаживания и его производной.
// Определение в заголовке (а не в kernel.cpp): цикл сборки примитивов
// build_solid_sph_search выбирает skf один раз вне цикла и вызывает
// kernel_cubic/gauss/quintic напрямую. Сами ядра -- [[gnu::noinline]]: при
// встраивании в цикл сборки проход записи на 1 потоке замедлялся на ~20 %
// (A: 40 мс против 33 мс с вызовом, speedup/prims_parallel.md).
//
// Параметр dx принимает std::span<const double> — это допускает передачу
// как std::vector<double>, так и std::array<double, dim> (и плоских буферов SoA в будущем),
// устраняя принудительное использование std::vector у вызывающих сторон.
// skf (тип ядра: 1 - кубический сплайн, 2 - Гаусс, 3 - квинтическое) передаётся
// явно (Interface Segregation): ядро — низкоуровневый примитив, ему нужен только
// тип ядра, а не весь SimulationConfig. Ранее читался из глобала skf (param.hpp).
#include <span>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include "param.hpp"

// Кубический сплайн (Monaghan 1985).
[[gnu::noinline]] inline void kernel_cubic(double r, std::span<const double> dx, double mhsml, double& w, double* dwdx)
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
        w = factor * (1.0/6.0) * ((2.0 - q)*(2.0 - q)*(2.0 - q));
        //w = factor * (1.0/6.0) * t * t * t;
        for (size_t d = 0; d < dim; ++d)
        {
            dwdx[d] = -factor * 0.5 * ((2.0 - q)*(2.0 - q)) * dx[d] / (r * mhsml);
            //dwdx[d] = -factor * 0.5 * t * t * dx[d] / (r * mhsml);
        }
    }
    // Иначе w и dwdx остаются нулевыми
}

// Гауссово ядро.
[[gnu::noinline]] inline void kernel_gauss(double r, std::span<const double> dx, double mhsml, double& w, double* dwdx)
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

// Квинтическое ядро (Morris 1997).
[[gnu::noinline]] inline void kernel_quintic(double r, std::span<const double> dx, double mhsml, double& w, double* dwdx)
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
        w = factor * (((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)) - 6.0 * ((2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q)) + 15.0 * ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)));
        for (size_t d = 0; d < dim; ++d)
        {
            dwdx[d] = factor * (-120.0 + 120.0*q - 50.0*q*q) * dx[d] / (mhsml * mhsml);
        }
    }
    else if (q > 1.0 && q <= 2.0)
    {
        w = factor * (((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)) - 6.0 * ((2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q)));
        for (size_t d = 0; d < dim; ++d)
        {
            if (r > 1e-8)
            {
                dwdx[d] = factor * (-5.0 * ((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)) + 30.0 * ((2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q))) * dx[d] / (r * mhsml);
            }
        }
    }
    else if (q > 2.0 && q <= 3.0)
    {
        w = factor * ((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q));
        for (size_t d = 0; d < dim; ++d)
        {
            if (r > 1e-8)
            {
                dwdx[d] = -factor * 5.0 * ((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)) * dx[d] / (r * mhsml);
            }
        }
    }
    // Иначе w и dwdx остаются нулевыми
}

// Вычисление ядра W и его производной dW/dx.
// skf — тип ядра; r — расстояние; dx — компоненты вектора x_ij;
// mhsml — средняя длина сглаживания; w — значение ядра; dwdx[d] — градиент ядра.
// При неизвестном skf w и dwdx нулевые.
inline void kernel(int skf, double r, std::span<const double> dx, double mhsml, double& w, double* dwdx)
{
    if (skf == 1) kernel_cubic(r, dx, mhsml, w, dwdx);
    else if (skf == 2) kernel_gauss(r, dx, mhsml, w, dwdx);
    else if (skf == 3) kernel_quintic(r, dx, mhsml, w, dwdx);
    else
    {
        w = 0.0;
        for (size_t d = 0; d < dim; d++) dwdx[d] = 0.0;
    }
}
