#pragma once
// interactions/solid_contact.hpp
// Solid-Solid-Contact: primitive/kernel для новой схемы RHS (fiber/SPH.md,
// раздел "Solid-Solid-Contact"; fiber/PLAN.md, этап Ф4, часть A).
//
// Физика силы (типы 1-3 -- pinball/Swegle/импульсный) ПОЛНОСТЬЮ переиспользует
// include/contact_force.hpp (contact::penalty_force -> pinball_force/
// swegle_force/impulse_force) -- ТУ ЖЕ функцию, что вызывает старый путь
// (src/ext_force.cpp). Формулы не продублированы и не изменены ни на бит.
//
// Геометрические помощники particle_radius()/pair_dt_ref() в src/ext_force.cpp
// объявлены в анонимном namespace той единицы трансляции и не экспортируются,
// поэтому здесь они воспроизведены заново под другими именами
// (solid_contact_radius / solid_contact_pair_dt_ref) -- сознательное,
// временное дублирование мелкой геометрической обвязки (fiber/PLAN.md, этап
// Ф4, пункт 5: "либо переиспользуемая функция... либо продублируй формулы
// с комментарием, что дублирование временное"). Для самой физики силы выбран
// первый вариант (реюз), для этой обвязки -- второй. Уйдёт в часть Б вместе
// с удалением src/ext_force.cpp.
//
// ОТЛИЧИЕ ОТ старого пути: признак контактной пары -- body[i] != body[j]
// (fiber/SPH.md, "Solid-Solid-Contact"; fiber/PLAN.md §3.1), а НЕ
// itype[i] != itype[j], как в ext_force.cpp::is_contact_material_pair.
// Для встроенных генераторов body[] заполняется равным itype
// (particle_system.hpp), поэтому на BC2/BC6 критерии дают одинаковый набор
// пар; расхождение возможно только для нескольких тел одного материала,
// которых текущие тесты не строят.
//
// fiber/CUDA_PLAN.md §3.2, этап Г2 (часть "Solid-Contact"): геометрия пары,
// PairGeom/PairProps и раскладка penalty_force по осям вынесены в
// physics::compute_solid_contact_pair (physics/solid_contact_pair.hpp) --
// SPH_HD-функция, общая для CPU и (в будущих частях Г2/Г4) GPU. Здесь --
// только сборка входа (в т.ч. модуль сдвига G из POD physics::MaterialTable,
// physics/material_table.hpp, вместо Materials/std::vector) и раскладка
// выхода в SolidDelta/SolidContactDiagnostics.
//
// solid_contact_radius/solid_contact_pair_dt_ref СОХРАНЕНЫ как тонкие
// обёртки над physics::solid_contact_particle_radius/
// solid_contact_particle_dt_ref -- под ПРЕЖНИМИ именем и сигнатурой: их
// напрямую вызывает include/interactions/fiber_solid_contact.hpp (Fiber-
// Solid-Contact), которого эта часть Г2 не касается.
#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include "contact_force.hpp"
#include "interactions/solid_sph.hpp" // SolidDelta
#include "materials.hpp"
#include "particle_system.hpp"
#include "physics/material_table.hpp"
#include "physics/solid_contact_pair.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"
#include "simulation_config.hpp"

namespace interactions
{

constexpr std::uint32_t kSolidContactInteractionType = 11;

// Тонкая обёртка над physics::solid_contact_particle_radius -- прежние имя и
// сигнатура (см. комментарий вверху файла: используется fiber_solid_contact.hpp).
inline real_t solid_contact_radius(const ParticleSystem& sys, std::size_t k)
{
    return physics::solid_contact_particle_radius(sys.mass[k], sys.rho[k], sys.rho_0[k]);
}

// Тонкая обёртка над physics::solid_contact_particle_dt_ref -- прежние имя и
// сигнатура. c -- значения скорости звука из SolidDerived (EOS этого RHS),
// а не sys.c (которое в новом пути не пишется).
inline real_t solid_contact_pair_dt_ref(const SimulationConfig& cfg, const ParticleSystem& sys,
                                         const std::vector<real_t>& c, std::size_t k)
{
    real_t vx[dim];
    for (std::size_t d = 0; d < dim; ++d) vx[d] = sys.vx[d][k];
    return physics::solid_contact_particle_dt_ref(static_cast<real_t>(cfg.CFL), sys.hsml[k], vx, c[k]);
}

// ----------------------------------------------------------------------
// compute_solid_contact_contributions -- kernel, эмитирует contributions
// (только dv[dim], остальные поля SolidDelta нулевые) в тот же reduce, что
// и Solid-SPH (kSolidContactInteractionType отличает interaction_type для
// детерминированного режима).
//
// real_pairs -- ПОДМНОЖЕСТВО (i,j), i<j, ОБА реальные, взятое из ТОГО ЖЕ
// поиска, что построил Solid-SPH (build_solid_sph_search) -- "Одна и та же
// построенная spatial structure может использоваться несколькими interaction
// types" (SPH.md, "Выполнение RHS"). Радиус этого поиска (scale_k*hsml)
// ШИРЕ фактического контактного диапазона (R_i+R_j) -- ТОЧНО как в старом
// пути, где ext_force фильтрует те же niac-пары от density/int_force, а не
// строит собственный más узкий список.
//
// dt_applied -- ФАКТИЧЕСКИЙ применяемый шаг интегрирования (не dt_ref пары):
// у типа 3 (импульсный) сила зависит от него напрямую, как в
// src/ext_force.cpp (параметр dt).
// ----------------------------------------------------------------------
inline void compute_solid_contact_contributions(
    const rhs::CpuExecutor& exec, const SimulationConfig& cfg, const Materials& mats,
    const ParticleSystem& sys, const std::vector<real_t>& c, real_t dt_applied,
    const std::vector<std::pair<rhs::index_t, rhs::index_t>>& real_pairs,
    std::vector<rhs::ContributionBuffer<SolidDelta>>& thread_buffers)
{
    const std::size_t nworkers = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const int contact_type = cfg.contact_force_type;
    if (contact_type == contact::Off || real_pairs.empty()) return;

    // Таблица модулей сдвига -- один раз на вызов (Materials за время RHS не
    // меняется), а не на пару: physics::compute_solid_contact_pair сама
    // Materials/std::vector не видит (правила SPH_HD-кода, план §3.2).
    const physics::MaterialTable mat_table = physics::build_material_table(mats);
    const real_t penalty_scale = static_cast<real_t>(cfg.contact_penalty_scale);
    const real_t impulse_gamma = static_cast<real_t>(cfg.contact_impulse_gamma);

    exec.parallel_for(real_pairs.size(), [&](std::size_t k) {
        const rhs::index_t i = real_pairs[k].first;
        const rhs::index_t j = real_pairs[k].second;
        if (sys.body[i] == sys.body[j]) return; // fiber/PLAN.md §3.1: критерий body[]

        // Сборка входа -- геометрия и материалы пары; сама математика (r, p,
        // PairGeom/PairProps, penalty_force, разложение силы по осям) --
        // physics::compute_solid_contact_pair (physics/solid_contact_pair.hpp).
        physics::SolidContactPairInput pair_in{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            pair_in.x_i[d] = sys.x[d][i];   pair_in.x_j[d] = sys.x[d][j];
            pair_in.vx_i[d] = sys.vx[d][i]; pair_in.vx_j[d] = sys.vx[d][j];
        }
        pair_in.mass_i = sys.mass[i]; pair_in.mass_j = sys.mass[j];
        pair_in.Ri = solid_contact_radius(sys, i);
        pair_in.Rj = solid_contact_radius(sys, j);
        pair_in.c_i = c[i]; pair_in.c_j = c[j];
        pair_in.G_i = physics::material_shear_modulus(mat_table, sys.itype[i]);
        pair_in.G_j = physics::material_shear_modulus(mat_table, sys.itype[j]);
        pair_in.dt_ref_i = solid_contact_pair_dt_ref(cfg, sys, c, i);
        pair_in.dt_ref_j = solid_contact_pair_dt_ref(cfg, sys, c, j);
        pair_in.contact_type = contact_type;
        pair_in.penalty_scale = penalty_scale;
        pair_in.dt_applied = dt_applied;
        pair_in.impulse_gamma = impulse_gamma;

        physics::SolidContactPairOutput pair_out;
        physics::compute_solid_contact_pair(pair_in, pair_out);
        if (!pair_out.contact) return;   // r<=0 || p<=0 -- как в исходнике
        if (pair_out.fn <= 0) return;    // res.fn<=0 -- как в исходнике

        SolidDelta di{}, dj{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            di.dv[d] = pair_out.dv_i[d];
            dj.dv[d] = pair_out.dv_j[d];
        }
        int w_id = rhs::CpuExecutor::worker_id();
        if (w_id < 0 || static_cast<std::size_t>(w_id) >= thread_buffers.size()) w_id = 0;
        rhs::ContributionWriter<SolidDelta> writer(thread_buffers[static_cast<std::size_t>(w_id)], kSolidContactInteractionType);
        writer.emit(i, di, static_cast<std::uint64_t>(k), 0);
        writer.emit(j, dj, static_cast<std::uint64_t>(k), 1);
    });
}

// ----------------------------------------------------------------------
// Диагностика (fiber/PLAN.md: "диагностика contact_pair_count,
// contact_max_delta -- через parallel_reduce"; dt_limit -- как сейчас в
// contact_force.hpp/ADR-007). Отдельный проход через CpuExecutor::
// parallel_reduce, буквально по формулировке плана; независим от
// compute_solid_contact_contributions (та же геометрия пересчитывается).
// Вызывается КАЖДЫЙ шаг из src/single_step.cpp (fiber/PLAN.md, этап
// Ф4-B1) -- ADR-007 требует min_dt_contact ограничивающим продовый шаг
// сейчас же, а не по отдельному запросу; цена второго прохода по real_pairs
// приемлема (см. src/single_step.cpp, комментарий у вызова).
// ----------------------------------------------------------------------
struct SolidContactDiagnostics
{
    std::size_t pair_count = 0;
    real_t max_delta = 0;
    real_t min_dt = contact::kNoLimit;
};

inline SolidContactDiagnostics compute_solid_contact_diagnostics(
    const rhs::CpuExecutor& exec, const SimulationConfig& cfg, const Materials& mats,
    const ParticleSystem& sys, const std::vector<real_t>& c, real_t dt_applied,
    const std::vector<std::pair<rhs::index_t, rhs::index_t>>& real_pairs)
{
    const SolidContactDiagnostics zero{};
    const int contact_type = cfg.contact_force_type;
    if (contact_type == contact::Off || real_pairs.empty()) return zero;

    const physics::MaterialTable mat_table = physics::build_material_table(mats);
    const real_t penalty_scale = static_cast<real_t>(cfg.contact_penalty_scale);
    const real_t impulse_gamma = static_cast<real_t>(cfg.contact_impulse_gamma);

    auto combine = [](SolidContactDiagnostics a, SolidContactDiagnostics b) {
        SolidContactDiagnostics r;
        r.pair_count = a.pair_count + b.pair_count;
        r.max_delta = std::max(a.max_delta, b.max_delta);
        r.min_dt = std::min(a.min_dt, b.min_dt);
        return r;
    };

    return exec.parallel_reduce(
        real_pairs.size(), zero,
        [&](std::size_t k) -> SolidContactDiagnostics {
            const rhs::index_t i = real_pairs[k].first;
            const rhs::index_t j = real_pairs[k].second;
            if (sys.body[i] == sys.body[j]) return zero;

            // ТА ЖЕ physics::compute_solid_contact_pair, что и в
            // compute_solid_contact_contributions -- геометрия пары считается
            // заново (как и в исходнике: "та же геометрия пересчитывается",
            // см. комментарий выше), поле out.contact/out.p/out.dt_limit
            // заполняются ДО проверки "fn<=0" внутри самой функции, поэтому
            // диагностика (в отличие от contributions) НЕ фильтрует по fn --
            // ровно как исходный код, который не проверял res.fn здесь.
            physics::SolidContactPairInput pair_in{};
            for (std::size_t d = 0; d < dim; ++d)
            {
                pair_in.x_i[d] = sys.x[d][i];   pair_in.x_j[d] = sys.x[d][j];
                pair_in.vx_i[d] = sys.vx[d][i]; pair_in.vx_j[d] = sys.vx[d][j];
            }
            pair_in.mass_i = sys.mass[i]; pair_in.mass_j = sys.mass[j];
            pair_in.Ri = solid_contact_radius(sys, i);
            pair_in.Rj = solid_contact_radius(sys, j);
            pair_in.c_i = c[i]; pair_in.c_j = c[j];
            pair_in.G_i = physics::material_shear_modulus(mat_table, sys.itype[i]);
            pair_in.G_j = physics::material_shear_modulus(mat_table, sys.itype[j]);
            pair_in.dt_ref_i = solid_contact_pair_dt_ref(cfg, sys, c, i);
            pair_in.dt_ref_j = solid_contact_pair_dt_ref(cfg, sys, c, j);
            pair_in.contact_type = contact_type;
            pair_in.penalty_scale = penalty_scale;
            pair_in.dt_applied = dt_applied;
            pair_in.impulse_gamma = impulse_gamma;

            physics::SolidContactPairOutput pair_out;
            physics::compute_solid_contact_pair(pair_in, pair_out);
            if (!pair_out.contact) return zero; // r<=0 || p<=0 -- как в исходнике

            SolidContactDiagnostics d1;
            d1.pair_count = 1;
            d1.max_delta = pair_out.p;
            d1.min_dt = pair_out.dt_limit;
            return d1;
        },
        combine);
}

} // namespace interactions
