// art_visc.cpp
void art_visc(std::vector<Particle>& Elem, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
              const std::vector<std::vector<double>>& dwdx, std::vector<std::vector<double>>& dvxdt, std::vector<double>& dedt)
{
    // Инициализация массивов
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            dvxdt[d][i] = 0.0;
        }
        dedt[i] = 0.0;
    }

    // Основной цикл по парам взаимодействий
#pragma omp parallel for schedule(auto)    
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k] - 1;  // Переводим индексы в C++ (начинаются с 0)
        size_t j = pair_j[k] - 1;
        double mhsml = (Elem[i].hsml + Elem[j].hsml) / 2.0;  // Средняя длина сглаживания
        // Вычисление относительной скорости и расстояния
        double vr = 0.0;
        double rr = 0.0;
        std::vector<double> dvx(dim);
        std::vector<double> dx(dim);

        for (size_t d = 0; d < dim; ++d)
        {
            dvx[d] = Elem[i].vx[d] - Elem[j].vx[d];
            dx[d] = Elem[i].x[d] - Elem[j].x[d];
            vr += dvx[d] * dx[d];
            rr += dx[d] * dx[d];
        }
        // Применяем искусственную вязкость только для сближающихся частиц
        if (vr < 0.0)
        {
            double muv = mhsml * vr / (rr + mhsml * mhsml * etq * etq); // Вязкость

            // Средние значения звука и плотности
            double mc = (Elem[i].c + Elem[j].c) / 2.0;
            double mrho = (Elem[i].rho + Elem[j].rho) / 2.0;
            double piv = (beta * muv - alpha * mc) * muv / mrho; // Искусственное давление

            // Обновление ускорения и внутренней энергии
            for (size_t d = 0; d < dim; ++d)
            {
                double h = -piv * dwdx[d][k]; // Сила на единицу массы
                double dvxdt_i = Elem[j].mass * h;
                double dvxdt_j = Elem[i].mass * h;
                double dedt_i = Elem[j].mass * dvx[d] * h;
                double dedt_j = Elem[i].mass * dvx[d] * h;
                // Применение сил симметрично
#pragma omp atomic update
                dvxdt[d][i] += dvxdt_i;
#pragma omp atomic update
                dvxdt[d][j] -= dvxdt_j;
                // Изменение энергии из-за работы вязкости
#pragma omp atomic update
                dedt[i] -= dedt_i;
#pragma omp atomic update
                dedt[j] -= dedt_j;
            }
        }
    }
    // Финальная корректировка энергии
    for (size_t i = 0; i < ntotal; ++i)
    {
        dedt[i] = 0.5 * dedt[i];
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

    // Проверка массива dedt
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (std::isnan(dedt[i]) || std::isinf(dedt[i]))
        {
            std::cerr << "ERROR: dedt[" << i << "] = "
                      << dedt[i] << " (NaN or Inf detected)" << std::endl;
            has_nan_or_inf = true;
        }
    }

    // Если обнаружены ошибки, выбрасываем исключение
    if (has_nan_or_inf)
    {
        throw std::runtime_error("NaN or Inf detected in art_visc");
    }
}
