// ext_force.cpp
// Определения функций внешних сил (отталкивание от границы, самогравитация).
// Объявления — в ext_force.hpp.
#include "ext_force.hpp"

void virtual_symm_2d (ParticleSystem& sys, size_t i, size_t j, std::vector<std::vector<double>>& dvxdt)
{
    constexpr double rr0 = 1.25e-5;  // Пороговое расстояние
    constexpr double dd = 1.0e-2;    // Сила отталкивания
    constexpr int p1 = 12;           // Степень для отталкивающей силы
    constexpr int p2 = 4;            // Степень для притягивающей силы
    if (axis_sim == true && dim == 2)
    {
        if (sys.hsml[i] >= sys.x[0][i])
        {
            // Зеркальный образ частицы j: координата x[0] отражается (-> -x[0][j]).
            // Работаем напрямую с полями sys, без создания локальной Particle.
            std::array<double, dim> dx{};
            double rr = 0.0;
            dx[0] = sys.x[0][i] - (-sys.x[0][j]);
            rr += dx[0] * dx[0];
            for (size_t d = 1; d < dim; ++d)
            {
                dx[d] = sys.x[d][i] - sys.x[d][j];
                rr += dx[d] * dx[d];
            }
            rr = std::sqrt(rr);
            if (rr < rr0)
            {
                double inv_rr = 1.0 / rr;
                double force = ((std::pow(rr0 * inv_rr, p1)) - (std::pow(rr0 * inv_rr, p2))) * inv_rr * inv_rr;

                // Применение силы
                for (size_t d = 0; d < dim; ++d)
                {
                    double f = dd * dx[d] * force;
                    dvxdt[d][i] += f;
                    dvxdt[d][j] -= f;
                }
            }
        }
    }
}

// Основная функция ext_force
void ext_force(ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j, std::vector<std::vector<double>>& dvxdt)
{

    // Параметры для силы отталкивания
    constexpr double rr0 = 1.25e-5;  // Пороговое расстояние
    constexpr double dd = 1.0e-2;    // Сила отталкивания
    constexpr int p1 = 12;           // Степень для отталкивающей силы
    constexpr int p2 = 4;            // Степень для притягивающей силы

    // Обнуление массива ускорений
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            dvxdt[d][i] = 0.0;
        }
    }

    // Учёт самогравитации
    if (self_gravity)
    {
        // Гравитация действует в последней размерности (например, z в 3D)
        size_t gravity_dim = dim - 1;
        for (size_t i = 0; i < ntotal; ++i)
        {
            dvxdt[gravity_dim][i] = -9.8;  // Ускорение свободного падения
        }
    }

    // Обработка всех пар взаимодействий
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k];  // 0-базовая индексация
        size_t j = pair_j[k];

        // Проверка: одна частица реальная (itype > 0), другая виртуальная (itype < 0)
        bool is_real_virtual_pair = (sys.itype[i] > 0 && sys.itype[j] < 0);
        //bool is_contact_pair_ = (itype[i] != itype[j]);
        if (is_real_virtual_pair)
        {
            // Вычисление разности координат
            std::array<double, dim> dx{};
            double rr = 0.0;

            for (size_t d = 0; d < dim; ++d)
            {
                dx[d] = sys.x[d][i] - sys.x[d][j];
                rr += dx[d] * dx[d];
            }

            rr = std::sqrt(rr);

            // Проверка на близость частиц
            if (rr < rr0)
            {
                double inv_rr = 1.0 / rr;
                double force = ((std::pow(rr0 * inv_rr, p1)) - (std::pow(rr0 * inv_rr, p2))) * inv_rr * inv_rr;

                // Применение силы
                for (size_t d = 0; d < dim; ++d)
                {
                    double f = dd * dx[d] * force;
                    dvxdt[d][i] += f;
                    dvxdt[d][j] -= f;
                }
            }
        }
        double A = 1e5;
        double n1 = 6;
        double n2 = 4;
        bool is_contact_pair_ = false;
        if (is_contact_pair_)
        {
            // Вычисление разности координат
            std::array<double, dim> dx{};
            double rr = 0.0;

            for (size_t d = 0; d < dim; ++d)
            {
                dx[d] = sys.x[d][i] - sys.x[d][j];
                rr += dx[d] * dx[d];
            }
            rr = std::sqrt(rr);
            double pe = 0.5 * (sys.hsml[i] + sys.hsml[j]) / rr;
            if (pe >= 1)
            {
                double inv_rr = 1.0 / rr;
                double force = A * ((std::pow(pe, n1)) - (std::pow(pe, n2))) * inv_rr * inv_rr;
                for (size_t d = 0; d < dim; ++d)
                {
                    double f = dx[d] * force;
                    dvxdt[d][i] += f;
                    dvxdt[d][j] -= f;
                }
            }
        }
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

    // Если обнаружены ошибки, выбрасываем исключение
    if (has_nan_or_inf)
    {
        throw std::runtime_error("NaN or Inf detected in ext_force");
    }
}
