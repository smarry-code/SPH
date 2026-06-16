#include <map>
#include <unordered_map>
#include <string>
#include <cmath>
#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>

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
                 std::vector<size_t>& pair_j, std::vector<double>& w, std::vector<std::vector<double>>& dwdx, std::vector<int>& countiac, double dx)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    int scale_k = 0;
    if (skf == 1) scale_k = 2;  // Кубический сплайн
    else if (skf == 2) scale_k = 3;  // Гауссово ядро
    else if (skf == 3) scale_k = 3;  // Квинтическое ядро
    countiac.assign(ntotal, 0);  // Обнуляем счетчики
    niac = 0;
    double max_hmsl = 0.0;
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (Elem[i].hsml > max_hmsl) max_hmsl = Elem[i].hsml;
    }
    dx = 2 * max_hmsl;
    CellMap Cell;
    for (size_t i = 0; i < ntotal; ++i)
    {
        std::string kl = "";
        if (dim == 1)
        {
            auto i_x = static_cast<int>(std::trunc(Elem[i].x[0] / dx));
            //auto i_x = int(Elem[i].x[0] / dx);
            kl += std::to_string(i_x);
        }
        else if (dim == 2)
        {
            auto i_x = static_cast<int>(std::trunc(Elem[i].x[0] / dx));
            auto i_y = static_cast<int>(std::trunc(Elem[i].x[1] / dx));
            //auto i_x = int(Elem[i].x[0] / dx);
            //auto i_y = int(Elem[i].x[1] / dx);
            kl += std::to_string(i_x) + "_" + std::to_string(i_y);
        }
        else if (dim == 3)
        {
            auto i_x = static_cast<int>(std::trunc(Elem[i].x[0] / dx));
            auto i_y = static_cast<int>(std::trunc(Elem[i].x[1] / dx));
            auto i_z = static_cast<int>(std::trunc(Elem[i].x[2] / dx));
            //auto i_x = int(Elem[i].x[0] / dx);
            //auto i_y = int(Elem[i].x[1] / dx);
            //auto i_z = int(Elem[i].x[2] / dx);
            kl += std::to_string(i_x) + "_" + std::to_string(i_y) + "_" + std::to_string(i_z);
        }
        Cell[kl].push_back(i);
    }
    for (size_t i = 0; i < ntotal; ++i)
    {
        std::string kl = "";
        int i_x, i_y, i_z = 0;
        if (dim == 1)
        {
            i_x = static_cast<int>(std::trunc(Elem[i].x[0] / dx));
            //i_x = int(Elem[i].x[0] / dx);
            kl += std::to_string(i_x);
        }
        else if (dim == 2)
        {
            i_x = static_cast<int>(std::trunc(Elem[i].x[0] / dx));
            i_y = static_cast<int>(std::trunc(Elem[i].x[1] / dx));
            //i_x = int(Elem[i].x[0] / dx);
            //i_y = int(Elem[i].x[1] / dx);
            kl += std::to_string(i_x) + "_" + std::to_string(i_y);
        }
        else if (dim == 3)
        {
            i_x = static_cast<int>(std::trunc(Elem[i].x[0] / dx));
            i_y = static_cast<int>(std::trunc(Elem[i].x[1] / dx));
            i_z = static_cast<int>(std::trunc(Elem[i].x[2] / dx));
            //i_x = int(Elem[i].x[0] / dx);
            //i_y = int(Elem[i].x[1] / dx);
            //i_z = int(Elem[i].x[2] / dx);
            kl += std::to_string(i_x) + "_" + std::to_string(i_y) + "_" + std::to_string(i_z);
        }
        int fac = static_cast<int>(std::round(scale_k * Elem[i].hsml / dx));
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
                    //double mhsml = Elem[i].hsml;
                    //if (Elem[j].hsml > Elem[i].hsml) mhsml = Elem[j].hsml;
                    double mhsml = 0.5 * (Elem[i].hsml + Elem[j].hsml);
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
        }
        else
        {
            std::vector <std::string> kl_1;
            //for (int k = -(fac - 1); k <= fac - 1; k ++)
            for (int k = -fac; k <= fac; k ++)
            {
                std::string kl_buf = "";
                int i_x_1, i_y_1, i_z_1 = 0;
                if (dim == 1)
                {
                    i_x_1 = i_x + k;
                    kl_buf += std::to_string(i_x_1);
                    kl_1.push_back(kl_buf);
                    kl_buf = "";
                }
                else if (dim == 2)
                {
                    i_x_1 = i_x + k;
                    //for (int k_y = -(fac - 1); k_y <= fac - 1; k_y++)
                    for (int k_y = -fac; k_y <= fac; k_y++)
                    {
                        i_y_1 = i_y + k_y;
                        kl_buf += std::to_string(i_x_1) + "_" + std::to_string(i_y_1);
                        kl_1.push_back(kl_buf);
                        kl_buf = "";
                    }
                }
                else if (dim == 3)
                {
                    i_x_1 = i_x + k;
                    //for (int k_y = -(fac - 1); k_y <= fac - 1; k_y++)
                    for (int k_y = -fac; k_y <= fac; k_y++)
                    {
                        i_y_1 = i_y + k_y;
                        //for (int k_z = -(fac - 1); k_z <= fac - 1; k_z++)
                        for (int k_z = -fac; k_z <= fac; k_z++)
                        {
                             i_z_1 = i_z + k_z;
                             kl_buf += std::to_string(i_x_1) + "_" + std::to_string(i_y_1) + "_" + std::to_string(i_z_1);
                             kl_1.push_back(kl_buf);
                             kl_buf = "";
                        }
                    }
                }
            }
            for (auto k : kl_1)
            {
                for (auto j : Cell[k])
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
            }
        }
    }
    // Статистика
    if (itimestep % print_step == 0 && int_stat)
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
    }
}
