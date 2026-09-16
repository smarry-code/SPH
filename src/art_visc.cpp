// art_visc.cpp
// Искусственная вязкость. Переведено на SoA: горячий цикл читает
// поля частиц через плоские __restrict__-указатели в ParticleSystem (gather
// перед циклом). scatter не нужен — выход идёт в плоские массивы dvxdt/dedt.
// Объявление — в art_visc.hpp.
#include "art_visc.hpp"

void art_visc(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
              const std::vector<std::vector<double>>& dwdx, std::vector<std::vector<double>>& dvxdt, std::vector<double>& dedt, const PairGraph& graph)
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

    // Мост AoS -> SoA: плоские __restrict__-указатели на читаемые поля.
    const double* __restrict__ phsml = sys.hsml.data();
    const double* __restrict__ pc    = sys.c.data();
    const double* __restrict__ prho  = sys.rho.data();
    const double* __restrict__ pmass = sys.mass.data();
    const double* __restrict__ pvx[dim];
    const double* __restrict__ px[dim];
    for (size_t d = 0; d < dim; ++d)
    {
        pvx[d] = sys.vx[d].data();
        px[d] = sys.x[d].data();
    }

    // Основной цикл по парам взаимодействий
    auto body = [&]<class Acc>(size_t i, size_t j, size_t k)
    {
                double mhsml = (phsml[i] + phsml[j]) / 2.0;  // Средняя длина сглаживания
                // Вычисление относительной скорости и расстояния
                double vr = 0.0;
                double rr = 0.0;
                std::array<double, dim> dvx{};
                std::array<double, dim> dx{};

                for (size_t d = 0; d < dim; ++d)
                {
                    dvx[d] = pvx[d][i] - pvx[d][j];
                    dx[d] = px[d][i] - px[d][j];
                    vr += dvx[d] * dx[d];
                    rr += dx[d] * dx[d];
                }
                // Применяем искусственную вязкость только для сближающихся частиц
                if (vr < 0.0)
                {
                    double muv = mhsml * vr / (rr + mhsml * mhsml * cfg.etq * cfg.etq); // Вязкость

                    // Средние значения звука и плотности
                    double mc = (pc[i] + pc[j]) / 2.0;
                    double mrho = (prho[i] + prho[j]) / 2.0;
                    double piv = (cfg.beta * muv - cfg.alpha * mc) * muv / mrho; // Искусственное давление

                    // Обновление ускорения и внутренней энергии
                    for (size_t d = 0; d < dim; ++d)
                    {
                        double h = -piv * dwdx[d][k]; // Сила на единицу массы
                        double dvxdt_i = pmass[j] * h;
                        double dvxdt_j = pmass[i] * h;
                        double dedt_i = pmass[j] * dvx[d] * h;
                        double dedt_j = pmass[i] * dvx[d] * h;
                        // Применение сил симметрично
                        Acc::add(dvxdt[d][i], dvxdt_i);
                        Acc::add(dvxdt[d][j], (-(dvxdt_j)));
                        // Изменение энергии из-за работы вязкости
                        Acc::add(dedt[i], (-(dedt_i)));
                        Acc::add(dedt[j], (-(dedt_j)));
                    }
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
