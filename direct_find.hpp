// direct_find.cpp
#include <algorithm>
#include "kernel.hpp"

// Основная функция direct_find
void direct_find(std::vector<Particle>& Elem, size_t ntotal, size_t& niac, std::vector<size_t>& pair_i,
                 std::vector<size_t>& pair_j, std::vector<double>& w, std::vector<std::vector<double>>& dwdx, std::vector<int>& countiac)
{
    // Определение scale_k на основе типа ядра
    int scale_k = 0;
    if (skf == 1) scale_k = 2;  // Кубический сплайн
    else if (skf == 2) scale_k = 3;  // Гауссово ядро
    else if (skf == 3) scale_k = 3;  // Квинтическое ядро

    // Инициализация массивов
    countiac.assign(ntotal, 0);  // Обнуляем счетчики
    niac = 0;

    // Основной цикл по всем парам частиц
    for (size_t i = 0; i < ntotal - 1; ++i)
    {
        for (size_t j = i + 1; j < ntotal; ++j)
        {
            // Вычисление разности координат
            std::vector<double> dxiac(dim);
            double driac = 0.0;

            for (size_t d = 0; d < dim; ++d)
            {
                dxiac[d] = Elem[i].x[d] - Elem[j].x[d];
                driac += dxiac[d] * dxiac[d];
            }
            double mhsml = (Elem[i].hsml + Elem[j].hsml) / 2.0;
            double r = std::sqrt(driac);
            // Проверка на взаимодействие
            if (r < scale_k * mhsml)
            {
                if (niac < max_interaction)
                {
                    // Сохраняем пару
                    pair_i[niac] = i + 1;  // Fortran-стиль индексации
                    pair_j[niac] = j + 1;

                    // Обновляем счетчики
                    countiac[i]++;
                    countiac[j]++;
                    // Вычисление ядра и его производной
                    double tdwdx[dim];
                    kernel(r, dxiac, mhsml, w[niac], tdwdx);
                    for (size_t d = 0; d < dim; ++d)
                    {
                        dwdx[d][niac] = tdwdx[d];
                    }
                    niac++;
                }
                else
                {
                    throw std::runtime_error("Превышено максимальное количество взаимодействий");
                }
            }
        }
    }

    // Статистика
    /*if (itimestep % print_step == 0 && int_stat)
    {
        int sumiac = 0;
        int maxiac = 0;
        int miniac = 1000;
        int noiac = 0;
        size_t maxp = 0;
        size_t minp = 0;

        for (size_t i = 0; i < ntotal; ++i)
        {
            sumiac += countiac[i];
            if (countiac[i] > maxiac)
            {
                maxiac = countiac[i];
                maxp = i + 1;
            }
            if (countiac[i] < miniac)
            {
                miniac = countiac[i];
                minp = i + 1;
            }
            if (countiac[i] == 0) noiac++;
        }
        std::ofstream file("log.txt", std::ios::app);
        file << " >> Статистика: взаимодействия на частицу:\n";
        file << "*** Частица " << maxp << " максимальные взаимодействия: " << maxiac << "\n";
        file << "*** Частица " << minp << " минимальные взаимодействия: " << miniac << "\n";
        file << "*** Среднее: " << double(sumiac)/double(ntotal) << "\n";
        file << "*** Общее число пар: " << niac << "\n";
        file << "*** Частицы без взаимодействий: " << noiac << "\n";
        file.close();
    }*/
}

