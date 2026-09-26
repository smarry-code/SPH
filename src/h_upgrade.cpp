// h_upgrade.cpp
// Определение функции обновления длины сглаживания.
// Объявление — в h_upgrade.hpp.
#include "h_upgrade.hpp"

// Основная функция h_upgrade
void h_upgrade(const SimulationConfig& cfg, ParticleSystem& sys, real_t dt, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
               const std::vector<std::vector<real_t>>& dwdx)
{
    // Метод 0: Длины сглаживания не обновляются
    if (cfg.sle == 0)
    {
        return;
    }

    // Метод 1: h = fac * (mass/rho)^(1/dim)
    else if (cfg.sle == 1)
    {
        constexpr real_t fac = real_t(1.3);         // Коэффициент масштабирования
        constexpr real_t rho_floor = real_t(1e-12); // защита плотности от 0/отрицательных/NaN
        for (size_t i = 0; i < ntotal; ++i)
        {
            sys.hsml[i] = fac * std::pow(sys.mass[i] / sys.rho[i], real_t(1.0) / dim);
        }
    }

    // Метод 2: dh/dt = (h/rho) * (dvx/dx + dvy/dy + dvz/dz)
    else if (cfg.sle == 2)
    {
        std::vector<real_t> vcc(ntotal, real_t(0.0));
        std::vector<real_t> dhsml(ntotal, real_t(0.0));
        // Вычисление дивергенции скорости
        for (size_t k = 0; k < niac; ++k)
        {
            size_t i = pair_i[k];  // 0-базовая индексация
            size_t j = pair_j[k];

            std::array<real_t, dim> dvx{};
            for (size_t d = 0; d < dim; ++d)
            {
                dvx[d] = sys.vx[d][i] - sys.vx[d][j];
            }

            real_t hvcc = 0;
            for (size_t d = 0; d < dim; ++d)
            {
                hvcc += dvx[d] * dwdx[d][k];
            }

            vcc[i] += sys.mass[j] * hvcc / sys.rho[j];
            vcc[j] += sys.mass[i] * hvcc / sys.rho[i];
        }

        // Обновление длины сглаживания
        for (size_t i = 0; i < ntotal; ++i)
        {
            dhsml[i] = (sys.hsml[i] / dim) * vcc[i];
            sys.hsml[i] += dt * dhsml[i];
            if (sys.hsml[i] <= 0)
            {
                sys.hsml[i] -= dt * dhsml[i];
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
        if (std::isnan(sys.hsml[i]) || std::isinf(sys.hsml[i]))
        {
            std::cerr << "ERROR: hsml[" << i << "] = "
                      << sys.hsml[i] << " (NaN or Inf detected)" << std::endl;
            has_nan_or_inf = true;
        }
    }

    // Если обнаружены ошибки, выбрасываем исключение
    if (has_nan_or_inf)
    {
        throw std::runtime_error("NaN or Inf detected in h_upgrade");
    }
}
