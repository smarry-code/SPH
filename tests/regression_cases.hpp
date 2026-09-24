#pragma once
// regression_cases.hpp
// Подготовка кейсов для регрессионного теста (tests/test_regression.cpp):
// BC2 и BC6 при dim = 2, уменьшенный BC2 (tests/bc2_mini3d_geometry.hpp) при
// dim = 3. Кейс другой размерности в сборке просто недоступен.
//
// Постановка ровно та же, что у smoke-режима сквозных тестов: та же геометрия
// (tests/bc2_geometry.hpp, tests/bc6_geometry.hpp), тот же шаг решётки
// Δx = 8e-4, те же конфиги tests/data/bc{2,6}_initial_conditions.txt, тот же
// порядок подготовки, что в prepare_run() из tests/test_e2e_bc6.cpp:
//   умолчания -> разбор конфига -> геометрия -> запас под зеркальные частицы
//   -> УРС (заполняет c) -> начальный шаг по Куранту.
// Отличие одно: t_total здесь не задаёт длину прогона. Останов делает предел
// шагов сторожа (RunGuard::set_step_limit), поэтому t_total остаётся из
// конфига — он заведомо больше, чем проходят регрессионные прогоны.

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>

#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "input.hpp"
#include "int_force.hpp"          // urs()
#include "time_integration.hpp"   // compute_cfl_timestep()

// Геометрии пишут координаты по индексам своей размерности, поэтому в сборку
// попадают только нужные: SPH_DIM — тот же макрос, из которого берётся dim.
#if SPH_DIM == 2
#include "bc2_geometry.hpp"
#include "bc6_case.hpp"
#include "bc6_geometry.hpp"
#elif SPH_DIM == 3
#include "bc2_mini3d_geometry.hpp"
#endif

namespace reg
{

enum class Case { BC2, BC6, BC2_MINI3D };

inline const char* case_name(Case c)
{
    switch (c)
    {
        case Case::BC2: return "bc2";
        case Case::BC6: return "bc6";
        default:        return "bc2mini3d";
    }
}

// Размерность, в которой кейс определён.
inline size_t case_dim(Case c) { return c == Case::BC2_MINI3D ? 3u : 2u; }

// Шаг решётки регрессии — smoke-разрешение сквозных тестов: 0.8 мм у всех
// кейсов. На нём 2152 (BC2), 2930 (BC6) и 3648 (уменьшенный BC2 в 3D)
// частиц: прогон занимает секунды, эталон — порядка мегабайта на точку.
inline double case_delta_x(Case) { return 8e-4; }

// Уменьшенный BC2 берёт конфиг BC2: те же материалы и параметры решателя.
inline std::filesystem::path default_config(Case c, const std::filesystem::path& data_dir)
{
    return data_dir / (c == Case::BC6 ? "bc6_initial_conditions.txt"
                                      : "bc2_initial_conditions.txt");
}

struct Prepared
{
    size_t ntotal = 0;
    size_t moni = 0;
    double delta_x = 0.0;
    double dt = 0.0;          // начальный шаг
};

inline Prepared prepare_case(Case c, int contact_force_type,
                             const std::filesystem::path& config,
                             SimulationConfig& cfg, ParticleSystem& sys)
{
    cfg.apply_defaults();
    // parse_sph_config при неоткрывшемся файле молча оставляет умолчания —
    // без этой проверки регрессия сверяла бы «не тот» кейс.
    if (!std::filesystem::exists(config))
    {
        throw std::runtime_error("не найден конфиг: " + config.string());
    }
    parse_sph_config(cfg, config.string());
    cfg.contact_force_type = contact_force_type;

    if (case_dim(c) != dim)
    {
        throw std::runtime_error(std::string("кейс ") + case_name(c) + " определён при dim = " +
                                 std::to_string(case_dim(c)) + ", а тест собран с dim = " +
                                 std::to_string(dim));
    }

    Prepared p;
    p.delta_x = case_delta_x(c);
#if SPH_DIM == 2
    if (c == Case::BC2)
    {
        const bc2::Bc2Geometry g = bc2::build_bc2_geometry(sys, p.delta_x);
        p.ntotal = g.ntotal;
        p.moni = g.moni_particle;
    }
    else
    {
        const bc6::Geometry g = bc6::build_geometry(sys, p.delta_x);
        p.ntotal = g.ntotal;
        p.moni = g.moni_axis;
    }
#elif SPH_DIM == 3
    {
        const bc2m3::Geometry g = bc2m3::build_geometry(sys, p.delta_x);
        p.ntotal = g.ntotal;
        p.moni = g.moni;
    }
#endif
    cfg.moni_particle = p.moni;

    // Запас под зеркальные частицы — как в main.cpp и в сквозных тестах.
    cfg.maxn = static_cast<size_t>(static_cast<double>(p.ntotal) * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    if (sys.size() < cfg.maxn)
    {
        const size_t prev = sys.size();
        sys.resize(cfg.maxn);
        for (size_t i = prev; i < sys.size(); ++i) sys.inzializ_index(i);
    }

    // c заполняет УРС, и только потом считается шаг.
    for (size_t i = 0; i < p.ntotal; ++i) urs(cfg.mats, sys, i);
    p.dt = compute_cfl_timestep(cfg, sys, p.ntotal);
    return p;
}

} // namespace reg
