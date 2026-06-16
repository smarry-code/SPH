// av_vel.cpp

// Основная функция av_vel
void av_vel(std::vector<Particle>& Elem, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i,
            const std::vector<size_t>& pair_j, const std::vector<double>& w, std::vector<std::vector<double>>& av)
{


    // Обнуление массива средней скорости
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            av[d][i] = 0.0;
        }
    }
    // Основной цикл по парам взаимодействий
#pragma omp parallel for schedule(auto)      
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k] - 1;  // Переводим индексы в C++ (начинаются с 0)
        size_t j = pair_j[k] - 1;
        std::vector<double> dvx(dim);
        for (size_t d = 0; d < dim; ++d)
        {
            dvx[d] = Elem[i].vx[d] - Elem[j].vx[d];
            double rho_sum = Elem[i].rho + Elem[j].rho;
            double av_i = 2.0 * Elem[j].mass * dvx[d] / rho_sum * w[k];
            double av_j = 2.0 * Elem[i].mass * dvx[d] / rho_sum * w[k];
#pragma omp atomic update
            av[d][i] -= av_i;
#pragma omp atomic update
            av[d][j] += av_j;
        }
    }

    // Масштабирование средней скорости
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            av[d][i] *= epsilon;
        }
    }

    // Проверка на NaN и бесконечность
    bool has_nan_or_inf = false;
    // Проверка массива av
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            if (std::isnan(av[d][i]) || std::isinf(av[d][i]))
            {
                std::cerr << "ERROR: av[" << d << "][" << i << "] = "
                          << av[d][i] << " (NaN or Inf detected)" << std::endl;
                has_nan_or_inf = true;
            }
        }
    }
    // Если обнаружены ошибки, выбрасываем исключение
    if (has_nan_or_inf)
    {
        throw std::runtime_error("NaN or Inf detected in av_vel");
    }
}
