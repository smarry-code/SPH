// direct_find.cpp
// Определение функции прямого поиска соседних частиц.
// Объявление — в direct_find.hpp.
#include "direct_find.hpp"

// Основная функция direct_find
void direct_find(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t& niac, std::vector<size_t>& pair_i,
                 std::vector<size_t>& pair_j, std::vector<double>& w, std::vector<std::vector<double>>& dwdx, std::vector<int>& countiac)
{
    // Определение scale_k на основе типа ядра
    int scale_k = 0;
    if (cfg.skf == 1) scale_k = 2;  // Кубический сплайн
    else if (cfg.skf == 2) scale_k = 3;  // Гауссово ядро
    else if (cfg.skf == 3) scale_k = 3;  // Квинтическое ядро

    // Инициализация массивов
    countiac.assign(ntotal, 0);  // Обнуляем счетчики
    niac = 0;

    // Основной цикл по всем парам частиц
    for (size_t i = 0; i < ntotal - 1; ++i)
    {
        for (size_t j = i + 1; j < ntotal; ++j)
        {
            // Вычисление разности координат
            std::array<double, dim> dxiac{};
            double driac = 0.0;

            for (size_t d = 0; d < dim; ++d)
            {
                dxiac[d] = sys.x[d][i] - sys.x[d][j];
                driac += dxiac[d] * dxiac[d];
            }
            double mhsml = (sys.hsml[i] + sys.hsml[j]) / 2.0;
            double r = std::sqrt(driac);
            // Проверка на взаимодействие
            if (r < scale_k * mhsml)
            {
                if (niac < cfg.max_interaction)
                {
                    // Сохраняем пару (0-базовая индексация)
                    pair_i[niac] = i;
                    pair_j[niac] = j;

                    // Обновляем счетчики
                    countiac[i]++;
                    countiac[j]++;
                    // Вычисление ядра и его производной
                    double tdwdx[dim];
                    kernel(cfg.skf, r, dxiac, mhsml, w[niac], tdwdx);
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
