// density.cpp
// Определения функций расчёта плотности (суммирование, уравнение непрерывности).
// Объявления — в density.hpp.
#include "density.hpp"

#include <stdexcept>
#include <string>

// Функция вычисления плотности через SPH-суммирование
void sum_density(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j, const std::vector<double>& w)
{
    // Инициализация временных массивов
    std::vector<double> wi(ntotal, 0.0);
    std::array<double, dim> hv{};
    double r = 0.0;

    // Вычисление самоплотности для каждой частицы
    for (size_t i = 0; i < ntotal; ++i)
    {
        double selfdens = 0.0;
        double tdwdx[dim] = {};
        kernel(cfg.skf, r, hv, sys.hsml[i], selfdens, tdwdx);
        wi[i] = selfdens * sys.mass[i] / sys.rho[i];
    }



    // Расчет нормализационного знаменателя
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k]; // 0-базовая индексация
        size_t j = pair_j[k];

        wi[i] += sys.mass[j] / sys.rho[j] * w[k];
        wi[j] += sys.mass[i] / sys.rho[i] * w[k];
    }

    // Инициализация плотности
    for (size_t i = 0; i < ntotal; ++i)
    {
        sys.rho[i] = 0.0;
    }
    // Добавление самовклада
    for (size_t i = 0; i < ntotal; ++i)
    {
        double selfdens = 0.0;
        double tdwdx[dim] = {};
        kernel(cfg.skf, r, hv, sys.hsml[i], selfdens, tdwdx);
        sys.rho[i] += selfdens * sys.mass[i];
    }

    // Обновление плотности
    for (size_t k = 0; k < niac; ++k)
    {
        size_t i = pair_i[k]; // 0-базовая индексация
        size_t j = pair_j[k];

        sys.rho[i] += sys.mass[j] * w[k];
        sys.rho[j] += sys.mass[i] * w[k];
    }

    // Нормализация плотности
    if (nor_density)
    {
        for (size_t i = 0; i < ntotal; ++i)
        {
            if (wi[i] > 1e-10)   // Защита от деления на ноль
            {
                sys.rho[i] /= wi[i];
            }
        }
    }
}

// Функция вычисления изменения плотности через уравнение непрерывности
void con_density(const SimulationConfig& cfg, const ParticleSystem& sys, size_t ntotal, size_t niac, const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
                 const std::vector<std::vector<double>>& dwdx, std::vector<double>& drhodt, const PairGraph& graph)
{
    // Локальные псевдонимы параметров плотности (тело функции и лямбды body
    // не изменяются): density_method, xi, epsilon_h ранее — глобалы param.hpp.
    const int density_method = cfg.density_method;
    const double xi = cfg.xi;
    const double epsilon_h = cfg.epsilon_h;
    // Инициализация производной плотности
    std::fill_n(drhodt.begin(), ntotal, 0.0);
    /*for (size_t i = 0; i < ntotal; ++i)
    {
        // Вычисление давления через уравнение состояния
        if (std::abs(itype[i]) == 1)
        {
            p_gas(rho[i], u[i], p[i], c[i]);
        }
        else if (std::abs(itype[i]) == 2)
        {
            p_art_water(rho[i], p[i], c[i]);
        }
        else if (std::abs(itype[i]) >= 3)
        {
            eos_mi_grunaizen(itype[i], rho[i], rho_0[i], u[i], p[i], c[i]);
        }
    } */
    // Мост AoS -> SoA: плоские __restrict__-указатели на читаемые поля.
    const double* __restrict__ pmass = sys.mass.data();
    const double* __restrict__ prho  = sys.rho.data();
    const double* __restrict__ phsml = sys.hsml.data();
    const double* __restrict__ pc    = sys.c.data();
    const double* __restrict__ pvx[dim];
    const double* __restrict__ px[dim];
    for (size_t d = 0; d < dim; ++d)
    {
        pvx[d] = sys.vx[d].data();
        px[d] = sys.x[d].data();
    }

    // Обработка всех пар взаимодействий
    auto body = [&]<class Acc>(size_t i, size_t j, size_t k)
    {
                std::array<double, dim> dvx{};
                for (size_t d = 0; d < dim; ++d)
                {
                    dvx[d] = pvx[d][i] - pvx[d][j];
                }
                double vcc = 0.0;
                for (size_t d = 0; d < dim; ++d)
                {
                    vcc += dvx[d] * dwdx[d][k];
                }

                if (density_method == 1)
                {
                    Acc::add(drhodt[i], pmass[j] * vcc);
                    Acc::add(drhodt[j], pmass[i] * vcc);
                }
                else if (density_method == 2)
                {
                    double q_i = (prho[i] / prho[j]) * pmass[j] * vcc;
                    double q_j = (prho[j] / prho[i]) * pmass[i] * vcc;
                    Acc::add(drhodt[i], q_i);
                    Acc::add(drhodt[j], q_j);
                }
                else if (density_method == 3)
                {
                    // Схема с искусственной диффузией плотности (Bonet-Lok / Molteni-Colagrossi)
                    // База: как density_method == 1 (простая форма)
                    double q_i = pmass[j] * vcc;
                    double q_j = pmass[i] * vcc;

                    // Диффузионный член D_i (формулы 15-16)
                    // h_ij = (h_i + h_j) / 2
                    double h_ij = 0.5 * (phsml[i] + phsml[j]);
            
                    // Вычисление x_i - x_j и |x_i - x_j|^2
                    double rr = 0.0;
                    std::array<double, dim> dx{};
                    for (size_t d = 0; d < dim; ++d)
                    {
                        dx[d] = px[d][i] - px[d][j];
                        rr += dx[d] * dx[d];
                    }
            
                    // Удельные объёмы v = 1/rho
                    double vi = 1.0 / prho[i];
                    double vj = 1.0 / prho[j];
            
                    // Знаменатель с регуляризатором: |x_ij|^2 + epsilon_h * h_ij^2
                    double denom = rr + epsilon_h * h_ij * h_ij;
            
                    // psi_ij = 2 * (vi/vj - 1) * x_ij / denom
                    // Коэффициенты для psi_ij и psi_ji (антисимметрия: x_ji = -x_ij)
                    double coef_i = 2.0 * (vi / vj - 1.0) / denom;
                    double coef_j = 2.0 * (vj / vi - 1.0) / denom;
            
                    // Скалярное произведение psi_ij · grad_W_ij
                    double psi_dot_dwdx_i = 0.0;
                    double psi_dot_dwdx_j = 0.0;
                    for (size_t d = 0; d < dim; ++d)
                    {
                        psi_dot_dwdx_i += (coef_i * dx[d]) * dwdx[d][k];
                        psi_dot_dwdx_j += (coef_j * (-dx[d])) * dwdx[d][k];
                    }
            
                    // Объёмы частиц: dV = m * v
                    double dVi = pmass[i] * vi;
                    double dVj = pmass[j] * vj;
                    double c_0 = 0.5 * (pc[i] + pc[j]);
                    // c0 — локальные скорости звука (прагматичный выбор)
                    // D_i имеет размерность 1/T (т.к. уравнение в 1.md записано для log(v))
                    double D_i = xi * h_ij * c_0 * psi_dot_dwdx_i * dVj;
                    double D_j = xi * h_ij * c_0 * psi_dot_dwdx_j * dVi;

                    // Умножение на -ρ_i для согласования размерностей: [drhodt] = [ρ]/T
                    Acc::add(drhodt[i], q_i - prho[i] * D_i);
                    Acc::add(drhodt[j], q_j - prho[j] * D_j);
                }
                else
                {
                    // Защита от деления на ноль для нормализованного варианта
                    double rho_i = (std::max)(prho[i], 1e-10);
                    double rho_j = (std::max)(prho[j], 1e-10);
                    // Fallback: нормализованная форма по умолчанию
                    std::cerr << "Предупреждение: неизвестный density_method = " << density_method
                              << ", используется вариант 2 (нормализованный)" << std::endl;
                    double q_i = (rho_i / rho_j) * pmass[j] * vcc;
                    double q_j = (rho_j / rho_i) * pmass[i] * vcc;
                    Acc::add(drhodt[i], q_i);
                    Acc::add(drhodt[j], q_j);
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

    // Осесимметричная поправка к уравнению непрерывности.
    // В цилиндрических координатах при осевой симметрии
    //     div(v) = dv_r/dr + v_r/r + dv_z/dz,
    // тогда как парная сумма выше даёт только плоскую часть (dv_r/dr + dv_z/dz).
    // Недостающий член: drho/dt += -rho_i * v_r_i / r_i — он относится к ЧАСТИЦЕ,
    // а не к паре, поэтому добавляется здесь, один раз на частицу.
    // Та же форма и то же место (пер-частичный цикл) — у осевого члена энергии
    // в int_force.cpp:620-624.
    if (axis_sim == true && dim == 2)
    {
        for (size_t i = 0; i < ntotal; ++i)
        {
            // Здесь ntotal = ntotal + nvirt: диапазон включает зеркальные
            // частицы, у которых r < 0 ПО ПОСТРОЕНИЮ, и для них член корректен
            // (v_r и r меняют знак вместе, отношение сохраняется). Поэтому
            // проверяется только истинная особенность r == 0: деление дало бы
            // Inf, который разошёлся бы по полям и проявился много позже и
            // далеко от причины. Клампы, удерживавшие осевые частицы
            // (time_integration.cpp, ветка sys.axis), закомментированы,
            // поэтому это состояние достижимо. Падаем в точке возникновения.
            if (px[0][i] == 0.0)
            {
                throw std::runtime_error(
                    "con_density: частица " + std::to_string(i) + " с r = "
                    + std::to_string(px[0][i]) + " == 0 (itype "
                    + std::to_string(sys.itype[i])
                    + "): осевой член неразрывности не определён");
            }
            drhodt[i] -= prho[i] * pvx[0][i] / px[0][i];
        }
    }
}
