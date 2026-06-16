// h_upgrade.cpp

// Основная функция h_upgrade
void h_upgrade(std::vector<Particle>& Elem, double dt, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
               const std::vector<std::vector<double>>& dwdx)
{
    // Метод 0: Длины сглаживания не обновляются
    if (sle == 0)
    {
        return;
    }

    // Метод 1: h = fac * (mass/rho)^(1/dim)
    else if (sle == 1)
    {
        constexpr double fac = 1.3;  // Коэффициент масштабирования
        for (size_t i = 0; i < ntotal; ++i)
        {
            Elem[i].hsml = fac * std::pow(Elem[i].mass / Elem[i].rho, 1.0 / dim);
        }
    }

    // Метод 2: dh/dt = (h/rho) * (dvx/dx + dvy/dy + dvz/dz)
    else if (sle == 2)
    {
        std::vector<double> vcc(maxn, 0.0);
        std::vector<double> dhsml(maxn, 0.0);
        // Вычисление дивергенции скорости
        for (size_t k = 0; k < niac; ++k)
        {
            size_t i = pair_i[k] - 1;  // Переводим в C++ индексацию
            size_t j = pair_j[k] - 1;

            std::vector<double> dvx(dim);
            for (size_t d = 0; d < dim; ++d)
            {
                dvx[d] = Elem[i].vx[d] - Elem[j].vx[d];
            }

            double hvcc = 0.0;
            for (size_t d = 0; d < dim; ++d)
            {
                hvcc += dvx[d] * dwdx[d][k];
            }

            vcc[i] += Elem[j].mass * hvcc / Elem[j].rho;
            vcc[j] += Elem[i].mass * hvcc / Elem[i].rho;
        }

        // Обновление длины сглаживания
        for (size_t i = 0; i < ntotal; ++i)
        {
            dhsml[i] = (Elem[i].hsml / dim) * vcc[i];
            Elem[i].hsml += dt * dhsml[i];
            if (Elem[i].hsml <= 0.0)
            {
                Elem[i].hsml -= dt * dhsml[i];
            }
        }
    }
    else
    {
        throw std::runtime_error("Unknown smoothing length evolution method (sle)");
    }

    // Проверка на NaN и бесконечность
    bool has_nan_or_inf = false;

    for (size_t i = 0; i < ntotal; ++i)
    {
        if (std::isnan(Elem[i].hsml) || std::isinf(Elem[i].hsml))
        {
            std::cerr << "ERROR: hsml[" << i << "] = "
                      << Elem[i].hsml << " (NaN or Inf detected)" << std::endl;
            has_nan_or_inf = true;
        }
    }

    // Если обнаружены ошибки, выбрасываем исключение
    if (has_nan_or_inf)
    {
        throw std::runtime_error("NaN or Inf detected in h_upgrade");
    }
}
