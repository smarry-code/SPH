#pragma once
// interactions/solid_contact.hpp
// Solid-Solid-Contact: primitive/kernel для новой схемы RHS (notes/fiber-model/fiber-problem-statement.md,
// раздел "Solid-Solid-Contact"; notes/fiber-model/fiber-plan.md, этап Ф4, часть A).
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
// временное дублирование мелкой геометрической обвязки (notes/fiber-model/fiber-plan.md, этап
// Ф4, пункт 5: "либо переиспользуемая функция... либо продублируй формулы
// с комментарием, что дублирование временное"). Для самой физики силы выбран
// первый вариант (реюз), для этой обвязки -- второй. Уйдёт в часть Б вместе
// с удалением src/ext_force.cpp.
//
// ОТЛИЧИЕ ОТ старого пути: признак контактной пары -- body[i] != body[j]
// (notes/fiber-model/fiber-problem-statement.md, "Solid-Solid-Contact"; notes/fiber-model/fiber-plan.md §3.1), а НЕ
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
inline double solid_contact_radius(const ParticleSystem& sys, std::size_t k)
{
    const double rho = sys.rho[k];
    if (!(rho > 0.0)) return 0.0;
    const double rho_eff = std::max(rho, sys.rho_0[k]);
    const double volume = sys.mass[k] / rho_eff;
    return 0.5 * std::pow(volume, 1.0 / static_cast<double>(dim));
}

// Копия pair_dt_ref (src/ext_force.cpp, анонимный namespace). c -- значения
// скорости звука из SolidDerived (EOS этого RHS), а не sys.c (которое в новом
// пути не пишется).
inline double solid_contact_pair_dt_ref(const SimulationConfig& cfg, const ParticleSystem& sys,
                                         const std::vector<double>& c, std::size_t k)
{
    double v2 = 0.0;
    for (std::size_t d = 0; d < dim; ++d) v2 += sys.vx[d][k] * sys.vx[d][k];
    const double denom = c[k] + std::sqrt(v2);
    return denom > 0.0 ? cfg.CFL * sys.hsml[k] / denom : 1e30;
}

// ----------------------------------------------------------------------
// compute_solid_contact_contributions -- kernel, эмитирует contributions
// (только dv[dim], остальные поля SolidDelta нулевые) в тот же reduce, что
// и Solid-SPH (kSolidContactInteractionType отличает interaction_type для
// детерминированного режима).
//
// contact_pairs -- межтельные (body[i] != body[j]) пары (i,j), i<j, ОБА
// реальные: подмножество real_pairs, взятое из ТОГО ЖЕ
// поиска, что построил Solid-SPH (build_solid_sph_search) -- "Одна и та же
// построенная spatial structure может использоваться несколькими interaction
// types" (fiber-problem-statement.md, "Выполнение RHS"). Радиус этого поиска (scale_k*hsml)
// ШИРЕ фактического контактного диапазона (R_i+R_j) -- ТОЧНО как в старом
// пути, где ext_force фильтрует те же niac-пары от density/int_force, а не
// строит собственный más узкий список. Исходный индекс пары в real_pairs (k)
// идёт в meta вкладов (primitive_id), как раньше.
//
// dt_applied -- ФАКТИЧЕСКИЙ применяемый шаг интегрирования (не dt_ref пары):
// у типа 3 (импульсный) сила зависит от него напрямую, как в
// src/ext_force.cpp (параметр dt).
// ----------------------------------------------------------------------
inline void compute_solid_contact_contributions(
    const rhs::CpuExecutor& exec, const SimulationConfig& cfg, const Materials& mats,
    const ParticleSystem& sys, const std::vector<double>& c, double dt_applied,
    const std::vector<SolidContactPair>& contact_pairs,
    std::vector<rhs::ContributionBuffer<SolidDelta>>& thread_buffers)
{
    const std::size_t nworkers = exec.det_segments(); // детерминированно: не зависит от числа потоков
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const int contact_type = cfg.contact_force_type;
    if (contact_type == contact::Off || contact_pairs.empty()) return;

    exec.parallel_for_segments(contact_pairs.size(), [&](std::size_t seg, std::size_t kk) {
        const rhs::index_t i = contact_pairs[kk].i;
        const rhs::index_t j = contact_pairs[kk].j;
        const std::uint64_t k = contact_pairs[kk].k; // исходный индекс в real_pairs
        // notes/fiber-model/fiber-plan.md §3.1: критерий body[]; список уже межтельный, проверка -- защита
        if (sys.body[i] == sys.body[j]) return;

        double dx[dim];
        double r2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d) { dx[d] = sys.x[d][i] - sys.x[d][j]; r2 += dx[d] * dx[d]; }
        const double r = std::sqrt(r2);
        const double Ri = solid_contact_radius(sys, i);
        const double Rj = solid_contact_radius(sys, j);
        const double p = (Ri + Rj) - r;
        if (!(r > 0.0) || !(p > 0.0)) return;

        double dot = 0.0, v2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double dv = sys.vx[d][i] - sys.vx[d][j];
            dot += dx[d] * dv;
            v2 += dv * dv;
        }
        contact::PairGeom g;
        g.r = r; g.p = p; g.pdot = -dot / r; g.v_rel = std::sqrt(v2);

        const double Gi = shear_modulus_by_itype(mats, sys.itype[i]);
        const double Gj = shear_modulus_by_itype(mats, sys.itype[j]);

        contact::PairProps pr;
        pr.m_red = sys.mass[i] * sys.mass[j] / (sys.mass[i] + sys.mass[j]);
        pr.c_pair = std::min(c[i], c[j]);
        pr.g_eff = (Gi > 0.0 && Gj > 0.0) ? Gi * Gj / (Gi + Gj) : 0.0;
        pr.r_eff = Ri * Rj / (Ri + Rj);
        pr.dt_ref = std::min(solid_contact_pair_dt_ref(cfg, sys, c, i), solid_contact_pair_dt_ref(cfg, sys, c, j));

        const contact::Result res = contact::penalty_force(contact_type, g, pr, cfg.contact_penalty_scale,
                                                             dt_applied, cfg.contact_impulse_gamma);
        if (res.fn <= 0.0) return;

        const double inv_r = 1.0 / r;
        SolidDelta di{}, dj{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double f = res.fn * dx[d] * inv_r;
            di.dv[d] = f / sys.mass[i];
            dj.dv[d] = -f / sys.mass[j];
        }
        rhs::ContributionWriter<SolidDelta> writer(thread_buffers[seg], kSolidContactInteractionType);
        writer.emit(i, di, k, 0);
        writer.emit(j, dj, k, 1);
    });
}

// ----------------------------------------------------------------------
// Диагностика (notes/fiber-model/fiber-plan.md: "диагностика contact_pair_count,
// contact_max_delta -- через parallel_reduce"; dt_limit -- как сейчас в
// contact_force.hpp/ADR-007). Отдельный проход через CpuExecutor::
// parallel_reduce, буквально по формулировке плана; независим от
// compute_solid_contact_contributions (та же геометрия пересчитывается).
// Идёт по contact_pairs -- межтельным парам того же поиска, а не по всем
// real_pairs (см. build_solid_sph_search).
// Вызывается КАЖДЫЙ шаг из src/single_step.cpp (notes/fiber-model/fiber-plan.md, этап
// Ф4-B1) -- ADR-007 требует min_dt_contact ограничивающим продовый шаг
// сейчас же, а не по отдельному запросу; цена второго прохода по
// contact_pairs приемлема (см. src/single_step.cpp, комментарий у вызова).
// ----------------------------------------------------------------------
struct SolidContactDiagnostics
{
    std::size_t pair_count = 0;
    double max_delta = 0.0;
    double min_dt = contact::kNoLimit;
};

inline SolidContactDiagnostics compute_solid_contact_diagnostics(
    const rhs::CpuExecutor& exec, const SimulationConfig& cfg, const Materials& mats,
    const ParticleSystem& sys, const std::vector<double>& c, double dt_applied,
    const std::vector<SolidContactPair>& contact_pairs)
{
    const SolidContactDiagnostics zero{};
    const int contact_type = cfg.contact_force_type;
    if (contact_type == contact::Off || contact_pairs.empty()) return zero;

    auto combine = [](SolidContactDiagnostics a, SolidContactDiagnostics b) {
        SolidContactDiagnostics r;
        r.pair_count = a.pair_count + b.pair_count;
        r.max_delta = std::max(a.max_delta, b.max_delta);
        r.min_dt = std::min(a.min_dt, b.min_dt);
        return r;
    };

    return exec.parallel_reduce(
        contact_pairs.size(), zero,
        [&](std::size_t kk) -> SolidContactDiagnostics {
            const rhs::index_t i = contact_pairs[kk].i;
            const rhs::index_t j = contact_pairs[kk].j;
            if (sys.body[i] == sys.body[j]) return zero; // защита: список уже межтельный

            double dx[dim];
            double r2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d) { dx[d] = sys.x[d][i] - sys.x[d][j]; r2 += dx[d] * dx[d]; }
            const double r = std::sqrt(r2);
            const double Ri = solid_contact_radius(sys, i);
            const double Rj = solid_contact_radius(sys, j);
            const double p = (Ri + Rj) - r;
            if (!(r > 0.0) || !(p > 0.0)) return zero;

            double dot = 0.0;
            for (std::size_t d = 0; d < dim; ++d) dot += dx[d] * (sys.vx[d][i] - sys.vx[d][j]);
            contact::PairGeom g;
            g.r = r; g.p = p; g.pdot = -dot / r;

            const double Gi = shear_modulus_by_itype(mats, sys.itype[i]);
            const double Gj = shear_modulus_by_itype(mats, sys.itype[j]);
            contact::PairProps pr;
            pr.m_red = sys.mass[i] * sys.mass[j] / (sys.mass[i] + sys.mass[j]);
            pr.c_pair = std::min(c[i], c[j]);
            pr.g_eff = (Gi > 0.0 && Gj > 0.0) ? Gi * Gj / (Gi + Gj) : 0.0;
            pr.r_eff = Ri * Rj / (Ri + Rj);
            pr.dt_ref = std::min(solid_contact_pair_dt_ref(cfg, sys, c, i), solid_contact_pair_dt_ref(cfg, sys, c, j));

            const contact::Result res = contact::penalty_force(contact_type, g, pr, cfg.contact_penalty_scale,
                                                                 dt_applied, cfg.contact_impulse_gamma);
            SolidContactDiagnostics d1;
            d1.pair_count = 1;
            d1.max_delta = p;
            d1.min_dt = res.dt_limit;
            return d1;
        },
        combine);
}

} // namespace interactions
