#include <map>
#include <unordered_map>
#include <string>
#include <cmath>
#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>
#include <omp.h>
namespace
{
// Статическая переменная для хранения времени начала выполнения программы
const auto program_start = std::chrono::high_resolution_clock::now();
}
void time_elapsed(double& elapsed)
{
    // Получаем текущее время
    auto now = std::chrono::high_resolution_clock::now();
    // Вычисляем разницу между текущим временем и временем старта программы
    std::chrono::duration<double> duration = now - program_start;
    // Сохраняем результат в переданную переменную
    elapsed = duration.count();
}

// Вспомогательная функция: комбинирование хешей (из boost/hash_combine)
template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    seed ^= std::hash<T> {}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Хеш для std::tuple
struct TupleHash
{
    template <typename... Args>
    std::size_t operator()(const std::tuple<Args...>& t) const
    {
        return hash_impl(t, std::index_sequence_for<Args...> {});
    }

private:
    template <typename Tuple, std::size_t... Is>
    std::size_t hash_impl(const Tuple& t, std::index_sequence<Is...>) const
    {
        std::size_t seed = 0;
        (hash_combine(seed, std::get<Is>(t)), ...); // C++17 fold expression
        return seed;
    }
};

using CellKey = std::tuple<int, int, int>;
using CellMap = std::unordered_map<CellKey, std::vector<size_t>, TupleHash>;

void link_list(std::vector<Particle>& Elem, int itimestep, size_t ntotal, size_t& niac, std::vector<size_t>& pair_i,
               std::vector<size_t>& pair_j, std::vector<double>& w, std::vector<std::vector<double>>& dwdx)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    double s1, s2 = 0.0;
    time_elapsed(s1);
    for (size_t i = 0; i < pair_i.size(); i++)
    {
        pair_i[i] = 0;
        pair_j[i] = 0;
    }
    int scale_k = 0;
    if (skf == 1) scale_k = 2;  // Кубический сплайн
    else if (skf == 2) scale_k = 3;  // Гауссово ядро
    else if (skf == 3) scale_k = 3;  // Квинтическое ядро
    double max_hmsl = 0.0;
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (Elem[i].hsml > max_hmsl) max_hmsl = Elem[i].hsml;
    }
    double dx = scale_k * max_hmsl;
    CellMap Cell;
    Cell.reserve(ntotal);
    std::vector<CellKey> particle_cell(ntotal);
    for (size_t i = 0; i < ntotal; ++i)
    {
        int ix = static_cast<int>(std::floor(Elem[i].x[0] / dx));
        int iy = (dim >= 2) ? static_cast<int>(std::floor(Elem[i].x[1] / dx)) : 0;
        int iz = (dim >= 3) ? static_cast<int>(std::floor(Elem[i].x[2] / dx)) : 0;

        CellKey key{ix, iy, iz};
        Cell[key].push_back(i);
        particle_cell[i] = key;
    }
    niac = 0;
    time_elapsed(s2);
    //std::wcout << "Затраченное время на последовательную секцию = " << s2 - s1 << " секунд" << std::endl;
    time_elapsed(s1);
    std::vector<double> thread_time(omp_get_max_threads(), 0.0);
    #pragma omp parallel
    {
        auto tid = omp_get_thread_num();
        double start = omp_get_wtime();  // Начало замера
        #pragma omp for schedule(auto)
        for (size_t i = 0; i < ntotal; ++i)
        {
            auto [ix, iy, iz] = particle_cell[i];
            CellKey kl {ix, iy, iz};
            int fac = static_cast<int>(std::round(scale_k * Elem[i].hsml / dx));
            int ky_n = 0;
            int ky_k = 0;
            int kz_n = 0;
            int kz_k = 0;
            if (dim >= 2)
            {
                ky_n = -fac;
                ky_k = fac;
            }
            if (dim >= 3)
            {
                kz_n = -fac;
                kz_k = fac;
            }
            //if (fac < 1) fac = 1;
            //if (fac == 1)
            if (fac < 1)
            {
                for (auto j : Cell[kl])
                {
                    if (j > i)
                    {
                        // Вычисление разности координат
                        std::vector<double> dxiac(dim);
                        double driac = 0.0;

                        for (size_t d = 0; d < dim; ++d)
                        {
                            dxiac[d] = Elem[i].x[d] - Elem[j].x[d];
                            driac += dxiac[d] * dxiac[d];
                        }
                        double mhsml = 0.5 * (Elem[i].hsml + Elem[j].hsml);
                        // Проверка на взаимодействие
                        double r = std::sqrt(driac);
                        if (driac < ((scale_k * mhsml) * (scale_k * mhsml)))
                        {
                            size_t idx;
                            #pragma omp atomic capture
                            {
                                idx = niac;
                                niac++;
                            }
                            if (idx < max_interaction)
                            {
                                double tdwdx[dim];
                                double w_buf = 0.0;
                                kernel(r, dxiac, mhsml,w_buf, tdwdx);
                                //countiac[i]++;
                                pair_i[idx] = i + 1;
                                pair_j[idx] = j + 1;
                                w[idx] = w_buf;
                                for (size_t d = 0; d < dim; ++d)
                                {
                                    dwdx[d][idx] = tdwdx[d];
                                    tdwdx[d] = 0.0;
                                }                               
                            }
                            else
                            {
                                std::cout << "Превышен  max_interaction, niac = " << niac << std::endl;
                            }
                        }
                    }
                }
            }
            else
            {
                for (int kx = -fac; kx <= fac; ++kx)
                {
                    for (int ky = ky_n; ky <= ky_k; ++ky)
                    {
                        for (int kz = kz_n; kz <= kz_k; ++kz)
                        {
                            CellKey neighbor_key{ix + kx, iy + ky, iz + kz};
                            auto it = Cell.find(neighbor_key);
                            if (it != Cell.end())
                            {
                                for (size_t j : it->second)
                                {
                                    if (j > i)
                                    {
                                        // Вычисление разности координат
                                        std::vector<double> dxiac(dim);
                                        double driac = 0.0;

                                        for (size_t d = 0; d < dim; ++d)
                                        {
                                            dxiac[d] = Elem[i].x[d] - Elem[j].x[d];
                                            driac += dxiac[d] * dxiac[d];
                                        }
                                        double mhsml = 0.5 * (Elem[i].hsml + Elem[j].hsml);
                                        double r = std::sqrt(driac);
                                        // Проверка на взаимодействие
                                        if (driac < ((scale_k * mhsml) * (scale_k * mhsml)))
                                        {
                                            size_t idx;
                                            #pragma omp atomic capture
                                            {
                                                idx = niac;
                                                niac++;
                                            }
                                            if (idx < max_interaction)
                                            {
                                                double tdwdx[dim];
                                                double w_buf = 0.0;
                                                kernel(r, dxiac, mhsml, w_buf, tdwdx);
                                                //countiac[i]++;
                                                pair_i[idx] = i + 1;
                                                pair_j[idx] = j + 1;
                                                w[idx] = w_buf;
                                                for (size_t d = 0; d < dim; ++d)
                                                {
                                                    dwdx[d][idx] = tdwdx[d];
                                                    tdwdx[d] = 0.0;
                                                }                                               
                                            }
                                            /*else
                                            {
                                                std::cout << "Превышен  max_interaction, niac = " << niac << std::endl;
                                            }*/
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            double end = omp_get_wtime();    // Конец замера
            // Сохраняем время в общий массив
            thread_time[tid] = end - start;
        }


    }
    /*for (int t = 0; t < omp_get_max_threads(); ++t)
    {
        std::cout << "Поток " << t << ": " << thread_time[t] << " сек\n";
    }
    time_elapsed(s2);
    std::wcout << "Затраченное время на параллельную секцию = " << s2 - s1 << " секунд" << std::endl;
    std::wcout << "Общее число взаимодействий частиц - " << niac << std::endl;*/
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
        //std::ios::sync_with_stdio(false);
        std::ofstream file("log.txt", std::ios::app);
        file << " >> Статистика: взаимодействия на частицу:\n";
        file << "*** Частица " << maxp << " максимальные взаимодействия: " << maxiac << "\n";
        file << "*** Частица " << minp << " минимальные взаимодействия: " << miniac << "\n";
        file << "*** Среднее: " << double(sumiac)/double(ntotal) << "\n";
        file << "*** Общее число пар: " << niac << "\n";
        file << "*** Частицы без взаимодействий: " << noiac << "\n";
        file.close();
    } */
}
