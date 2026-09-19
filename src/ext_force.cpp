// ext_force.cpp
// Определения функций внешних сил (отталкивание от границы, самогравитация,
// контакт «частица-частица» — part_0.md, ADR-007). Объявления — в ext_force.hpp.
#include "ext_force.hpp"
#include "materials.hpp"

namespace
{

// Признак контактной пары (решение D2=A): разные материалы, обе частицы реальные.
// Ограничение: два тела из одного материала контакта не увидят — в проекте нет
// идентификатора тела (part id из k-файла затирается номером материала,
// input.cpp:230 и 786-807). Смена критерия — правка только этой функции.
//
// nreal — число реальных частиц; зеркальные лежат за ними, в [nreal, ntotal),
// и наследуют положительный itype оригинала (virt_part.cpp -> inzializ_2d_symm_index),
// поэтому по одному материалу от реальных неотличимы. Контакт им запрещён:
//  - призрак представляет то же тело на угле theta+pi, и пара «реальная + призрак»
//    дублирует отталкивание, уже описанное прямой парой «реальная + реальная»;
//  - реакция, записанная в слот призрака, теряется: интегрирование идёт по
//    [0, nreal), так что импульс пары не сохраняется;
//  - счётчики wb (contact_pair_count, contact_max_delta, min_dt_contact)
//    засорялись бы фантомными парами, а min_dt_contact ещё и урезал бы шаг
//    всего расчёта по несуществующему контакту.
// Штрафной контакт определён между дискретными телами, а призрак — не тело,
// а приём замыкания носителя ядра в континуальных SPH-суммах.
bool is_contact_material_pair(const ParticleSystem& sys, size_t i, size_t j, size_t nreal)
{
    return i < nreal && j < nreal
        && sys.itype[i] > 0 && sys.itype[j] > 0 && sys.itype[i] != sys.itype[j];
}

// Локальный CFL-шаг пары: тот же критерий, что compute_cfl_timestep
// (time_integration.cpp:13-27), суженный на пару. Гарантирует
// dt_applied <= dt_ref, поэтому жёсткость от dt_ref не образует петли
// «жёсткость -> шаг -> жёсткость».
double pair_dt_ref(const SimulationConfig& cfg, const ParticleSystem& sys, size_t k)
{
    double v2 = 0.0;
    for (size_t d = 0; d < dim; ++d) v2 += sys.vx[d][k] * sys.vx[d][k];
    const double denom = sys.c[k] + std::sqrt(v2);
    return denom > 0.0 ? cfg.CFL * sys.hsml[k] / denom : 1e30;
}

} // namespace

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
void ext_force(const SimulationConfig& cfg, ParticleSystem& sys, double dt,
               size_t ntotal, size_t nreal, size_t niac,
               const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
               std::vector<std::vector<double>>& dvxdt, WorkBuffers& wb)
{
    const int contact_type = cfg.contact_force_type;
    // Защита от рассогласования вызывающей стороны: nreal больше полного
    // диапазона означало бы, что «реальными» объявлены несуществующие индексы.
    if (nreal > ntotal)
    {
        throw std::runtime_error("ext_force: nreal > ntotal");
    }

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
    wb.reset_contact();

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
        if (contact_type != contact::Off && is_contact_material_pair(sys, i, j, nreal))
        {
            std::array<double, dim> dx{};      // x_i - x_j
            double r2 = 0.0;
            for (size_t d = 0; d < dim; ++d) { dx[d] = sys.x[d][i] - sys.x[d][j]; r2 += dx[d]*dx[d]; }
            const double r = std::sqrt(r2);
            const double p = 0.5 * (sys.hsml[i] + sys.hsml[j]) - r;    // (2.1)

            if (r > 0.0 && p > 0.0)
            {
                double dot = 0.0, v2 = 0.0;
                for (size_t d = 0; d < dim; ++d)
                {
                    const double dv = sys.vx[d][i] - sys.vx[d][j];
                    dot += dx[d] * dv;
                    v2  += dv * dv;
                }
                contact::PairGeom g;
                g.r = r;  g.p = p;
                g.pdot  = -dot / r;                 // > 0 при сближении
                g.v_rel = std::sqrt(v2);

                const double Ri = 0.5 * sys.hsml[i], Rj = 0.5 * sys.hsml[j];
                const double Gi = shear_modulus_by_itype(cfg.mats, sys.itype[i]);
                const double Gj = shear_modulus_by_itype(cfg.mats, sys.itype[j]);

                contact::PairProps pr;
                pr.m_red  = sys.mass[i] * sys.mass[j] / (sys.mass[i] + sys.mass[j]);
                pr.c_pair = std::min(sys.c[i], sys.c[j]);
                pr.g_eff  = (Gi > 0.0 && Gj > 0.0) ? Gi * Gj / (Gi + Gj) : 0.0;
                pr.r_eff  = Ri * Rj / (Ri + Rj);
                pr.dt_ref = std::min(pair_dt_ref(cfg, sys, i), pair_dt_ref(cfg, sys, j));

                const contact::Result res =
                    contact::penalty_force(contact_type, g, pr, cfg.contact_penalty_scale, dt);

                ++wb.contact_pair_count;
                if (p > wb.contact_max_delta) wb.contact_max_delta = p;
                if (res.dt_limit < wb.min_dt_contact) wb.min_dt_contact = res.dt_limit;

                if (res.fn > 0.0)
                {
                    const double inv_r = 1.0 / r;
                    for (size_t d = 0; d < dim; ++d)
                    {
                        const double f = res.fn * dx[d] * inv_r;    // вдоль n = (x_i - x_j)/r
                        dvxdt[d][i] += f / sys.mass[i];             // dvxdt — ускорения
                        dvxdt[d][j] -= f / sys.mass[j];             // точный баланс импульса
                    }
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
