// av_vel.cpp
// Определение функции усреднения скорости для предотвращения проникания.
// Переведено на SoA: горячий цикл читает поля частиц через плоские
// __restrict__-указатели (gather перед циклом). scatter не нужен — выход
// в плоский массив av.
// Объявление — в av_vel.hpp.
#include "av_vel.hpp"

// Основная функция av_vel
void av_vel(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i,
            const std::vector<size_t>& pair_j, const std::vector<double>& w, std::vector<std::vector<double>>& av, const PairGraph& graph)
{


    // Обнуление массива средней скорости
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            av[d][i] = 0.0;
        }
    }

    // Мост AoS -> SoA: плоские __restrict__-указатели на читаемые поля.
    const double* __restrict__ prho  = sys.rho.data();
    const double* __restrict__ pmass = sys.mass.data();
    const double* __restrict__ pvx[dim];
    for (size_t d = 0; d < dim; ++d)
    {
        pvx[d] = sys.vx[d].data();
    }

    // Основной цикл по парам взаимодействий
    auto body = [&]<class Acc>(size_t i, size_t j, size_t k)
    {
                std::array<double, dim> dvx{};
                for (size_t d = 0; d < dim; ++d)
                {
                    dvx[d] = pvx[d][i] - pvx[d][j];
                    double rho_sum = prho[i] + prho[j];
                    double av_i = 2.0 * pmass[j] * dvx[d] / rho_sum * w[k];
                    double av_j = 2.0 * pmass[i] * dvx[d] / rho_sum * w[k];
                    Acc::add(av[d][i], (-(av_i)));
                    Acc::add(av[d][j], av_j);
                }
    
    };
    // Единая параллельная область: одна команда потоков на весь цикл,
    // между цветами — только неявный барьер #pragma omp for (без fork-join).
    // Ранее for_each_color создавал отдельный parallel-регион на каждый цвет
    // (num_colors fork-join за вызов). Пары одного цвета не разделяют частиц
    // => конкурентной записи нет => atomic не нужен.
#pragma omp parallel
    {
        if (graph.num_colors > 0)
        {
            for (int _col = 0; _col < graph.num_colors; ++_col)
            {
#pragma omp for schedule(static)
                for (size_t _p = graph.color_offset[static_cast<size_t>(_col)]; _p < graph.color_offset[static_cast<size_t>(_col) + 1]; ++_p)
                {
                    size_t k = graph.pairs_by_color[_p];
                    body.template operator()<AccPlain>(pair_i[k], pair_j[k], k);
                }
            }
        }
        else
        {
#pragma omp for schedule(auto)
            for (size_t k = 0; k < niac; ++k)
            {
                body.template operator()<AccAtomic>(pair_i[k], pair_j[k], k);
            }
        }
    }

    // Масштабирование средней скорости
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            av[d][i] *= cfg.epsilon;
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
