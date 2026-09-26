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
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "contact_force.hpp"
#include "interactions/solid_sph.hpp" // SolidDelta
#include "materials.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"
#include "simulation_config.hpp"

namespace interactions
{

constexpr std::uint32_t kSolidContactInteractionType = 11;

// Копия particle_radius (src/ext_force.cpp, анонимный namespace) -- дословно
// та же формула: R = 0.5*(m/max(rho,rho_0))^(1/dim).
inline real_t solid_contact_radius(const ParticleSystem& sys, std::size_t k)
{
    const real_t rho = sys.rho[k];
    if (!(rho > 0)) return real_t(0.0);
    const real_t rho_eff = std::max(rho, sys.rho_0[k]);
    const real_t volume = sys.mass[k] / rho_eff;
    return real_t(0.5) * std::pow(volume, real_t(1.0) / static_cast<real_t>(dim));
}

// Копия pair_dt_ref (src/ext_force.cpp, анонимный namespace). c -- значения
// скорости звука из SolidDerived (EOS этого RHS), а не sys.c (которое в новом
// пути не пишется).
inline real_t solid_contact_pair_dt_ref(const SimulationConfig& cfg, const ParticleSystem& sys,
                                         const std::vector<real_t>& c, std::size_t k)
{
    real_t v2 = 0;
    for (std::size_t d = 0; d < dim; ++d) v2 += sys.vx[d][k] * sys.vx[d][k];
    const real_t denom = c[k] + std::sqrt(v2);
    return denom > 0 ? static_cast<real_t>(cfg.CFL) * sys.hsml[k] / denom : real_t(1e30);
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

    exec.parallel_for(real_pairs.size(), [&](std::size_t k) {
        const rhs::index_t i = real_pairs[k].first;
        const rhs::index_t j = real_pairs[k].second;
        if (sys.body[i] == sys.body[j]) return; // fiber/PLAN.md §3.1: критерий body[]

        real_t dx[dim];
        real_t r2 = 0;
        for (std::size_t d = 0; d < dim; ++d) { dx[d] = sys.x[d][i] - sys.x[d][j]; r2 += dx[d] * dx[d]; }
        const real_t r = std::sqrt(r2);
        const real_t Ri = solid_contact_radius(sys, i);
        const real_t Rj = solid_contact_radius(sys, j);
        const real_t p = (Ri + Rj) - r;
        if (!(r > 0) || !(p > 0)) return;

        real_t dot = 0, v2 = 0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const real_t dv = sys.vx[d][i] - sys.vx[d][j];
            dot += dx[d] * dv;
            v2 += dv * dv;
        }
        contact::PairGeom g;
        g.r = r; g.p = p; g.pdot = -dot / r; g.v_rel = std::sqrt(v2);

        const real_t Gi = shear_modulus_by_itype(mats, sys.itype[i]);
        const real_t Gj = shear_modulus_by_itype(mats, sys.itype[j]);

        contact::PairProps pr;
        pr.m_red = sys.mass[i] * sys.mass[j] / (sys.mass[i] + sys.mass[j]);
        pr.c_pair = std::min(c[i], c[j]);
        pr.g_eff = (Gi > 0 && Gj > 0) ? Gi * Gj / (Gi + Gj) : real_t(0.0);
        pr.r_eff = Ri * Rj / (Ri + Rj);
        pr.dt_ref = std::min(solid_contact_pair_dt_ref(cfg, sys, c, i), solid_contact_pair_dt_ref(cfg, sys, c, j));

        const contact::Result res = contact::penalty_force(contact_type, g, pr, static_cast<real_t>(cfg.contact_penalty_scale),
                                                             dt_applied, static_cast<real_t>(cfg.contact_impulse_gamma));
        if (res.fn <= 0) return;

        const real_t inv_r = real_t(1.0) / r;
        SolidDelta di{}, dj{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            const real_t f = res.fn * dx[d] * inv_r;
            di.dv[d] = f / sys.mass[i];
            dj.dv[d] = -f / sys.mass[j];
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

            real_t dx[dim];
            real_t r2 = 0;
            for (std::size_t d = 0; d < dim; ++d) { dx[d] = sys.x[d][i] - sys.x[d][j]; r2 += dx[d] * dx[d]; }
            const real_t r = std::sqrt(r2);
            const real_t Ri = solid_contact_radius(sys, i);
            const real_t Rj = solid_contact_radius(sys, j);
            const real_t p = (Ri + Rj) - r;
            if (!(r > 0) || !(p > 0)) return zero;

            real_t dot = 0;
            for (std::size_t d = 0; d < dim; ++d) dot += dx[d] * (sys.vx[d][i] - sys.vx[d][j]);
            contact::PairGeom g;
            g.r = r; g.p = p; g.pdot = -dot / r;

            const real_t Gi = shear_modulus_by_itype(mats, sys.itype[i]);
            const real_t Gj = shear_modulus_by_itype(mats, sys.itype[j]);
            contact::PairProps pr;
            pr.m_red = sys.mass[i] * sys.mass[j] / (sys.mass[i] + sys.mass[j]);
            pr.c_pair = std::min(c[i], c[j]);
            pr.g_eff = (Gi > 0 && Gj > 0) ? Gi * Gj / (Gi + Gj) : real_t(0.0);
            pr.r_eff = Ri * Rj / (Ri + Rj);
            pr.dt_ref = std::min(solid_contact_pair_dt_ref(cfg, sys, c, i), solid_contact_pair_dt_ref(cfg, sys, c, j));

            const contact::Result res = contact::penalty_force(contact_type, g, pr, static_cast<real_t>(cfg.contact_penalty_scale),
                                                                 dt_applied, static_cast<real_t>(cfg.contact_impulse_gamma));
            SolidContactDiagnostics d1;
            d1.pair_count = 1;
            d1.max_delta = p;
            d1.min_dt = res.dt_limit;
            return d1;
        },
        combine);
}

} // namespace interactions
