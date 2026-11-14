// time_integration.cpp
// Определения функций временной интеграции.
// Объявления — в time_integration.hpp.
#include "time_integration.hpp"
#include "checkpoint.hpp"
#include "output.hpp"
#include "fiber/fiber_output.hpp"
#include "interactions/clamp.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_friction.hpp"
#include "interactions/fiber_normal_projection.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "interactions/postprocessing.hpp"
#include "ode/ode_solver.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include <fstream>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <limits>
#ifdef SPH_CONTACT_STATS
#include "interactions/fiber_contact_stats.hpp"
#include <memory>
#endif
#ifdef SPH_ENERGY_LOG
#include "energy_log.hpp"
#include <memory>
#endif

#ifdef _WIN32
#define NOMINMAX // std::min/std::max (rhs/contributions.hpp, interactions/fiber_bond.hpp) конфликтуют с макросами min/max windows.h
#include <windows.h>
#endif

double compute_cfl_timestep(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal)
{
    double min_dt = 1e6;                      // шаг по куранту
    for (size_t i = 0; i < ntotal; ++i)
    {
        double v = 0.0;
        for (size_t d = 0; d < dim; ++d)
        {
             v +=  sys.vx[d][i] * sys.vx[d][i];
        }
        v = std::sqrt(v);  // Скорость частицы
        double dt = cfg.CFL * sys.hsml[i] / (sys.c[i] + v);
        if (dt < min_dt) min_dt = dt;
    }
    return min_dt;
}

namespace detail
{
// Связь цикла с чекпоинтами (include/checkpoint.hpp). nullptr вместо указателя на
// структуру -- цикл побитово прежний и без накладных расходов.
//   restart       -- восстановленное состояние цикла (после checkpoint::read в обёртке;
//                    Model уже перезаписана): цикл продолжается с него, а не с нуля.
//   on_checkpoint -- запись чекпоинта (обёртка замыкает cfg и Model и зовёт checkpoint::write);
//                    пуст -- чекпоинты не пишутся.
struct CheckpointHooks
{
    const checkpoint::Payload* restart = nullptr;
    std::function<void(checkpoint::Payload&&)> on_checkpoint;
};

// Копия вектора истории в буфер решателя: пустой источник -- пропуск (EulerSolver),
// размер больше буфера -- ошибка (чекпоинт от другой сцены).
inline void restore_history_vec(std::vector<double>& dst, const std::vector<double>& src)
{
    if (src.empty()) return;
    if (src.size() > dst.size()) throw std::runtime_error("рестарт: размер истории решателя в чекпоинте больше буфера");
    for (size_t i = 0; i < src.size(); ++i) dst[i] = src[i];
}

// ----------------------------------------------------------------------
// run_time_integration<Solver> -- ЕДИНСТВЕННАЯ реализация главного цикла по
// времени (fiber/PLAN.md, этап Ф8, раздел E «Унификация time_integration»).
// Обе публичные перегрузки time_integration() ниже -- тонкие обёртки над этим
// шаблоном: разница между "Solid-only" (старая сигнатура, ParticleSystem&) и
// "Solid+Fiber" (fiber::Model&) -- только в том, какие fb/topology/materials
// передаются (пустые по умолчанию -> nfiber == 0 -> вся Fiber-физика ниже
// no-op по существующим `if (nfiber > 0)`-веткам, как и раньше в перегрузке
// для fiber::Model). Никакой копии ParticleSystem не возникает — оба
// вызывающих места передают его по ссылке.
//
// Solver -- ode::PredictorCorrectorSolver (по умолчанию, fiber/PLAN.md §Ф8)
// или ode::EulerSolver (include/ode/ode_solver.hpp); выбор — по
// cfg.ode_solver, ОДИН раз перед циклом while, на стороне вызывающих
// перегрузок. Интерфейс обоих решателей единообразен (predict()/advance()),
// поэтому тело цикла не содержит if/switch по типу решателя (fiber/PLAN.md
// §3.9: «Дизайн — выбор через if/switch по enum на этапе вызова» — здесь этот
// выбор целиком сведён к выбору ШАБЛОННОГО ПАРАМЕТРА при инстанцировании).
//
// Порядок одного шага (fiber/SPH.md, «Описание алгоритмов и структур
// данных» / PLAN.md §3.9): solver.predict() -> RHS (single_step + Fiber-Bond/
// Fiber-Fiber-Contact/Fiber-Solid-Contact, ОДНО вычисление на шаг) -> новый
// [1] проекция закрепления -> dt -> единая стадия трения нить-нить и Fiber-Solid
// (fiber/fiber_friction_projection.md, fiber/fiber_solid_friction_implementation.md;
// внутри неё [2] c_k = c_s = 0 у полностью зажатых и [3] повторная проекция
// после прибавления трения) -> solver.advance() -> Postprocessing (apply_material_strength, apply_fiber_bond_constitutive, разрыв
// связей Fiber-Bond) -> вывод. update_smoothing_length -- ИСКЛЮЧЕНИЕ, см.
// большой комментарий в include/interactions/postprocessing.hpp (шапка
// файла) -- почему он остаётся МЕЖДУ RHS и вычислением нового dt, а не после
// advance().
//
// Точки [1]-[3] -- стадия constraint кинематического закрепления
// (interactions/clamp.hpp; fiber/SPH.md, «Clamp (кинематическое
// закрепление)»): выполняются ТОЛЬКО при непустом clamps, иначе цикл
// побитово прежний.
// ----------------------------------------------------------------------
template <class Solver>
void run_time_integration(const SimulationConfig& cfg, ParticleSystem& sys, fiber::FiberState& fb,
                           fiber::FiberTopology& topology, std::vector<fiber::FiberMaterial>& materials,
                           size_t ntotal, double& dt, double& time, double /*delta_x*/, RunGuard* guard,
                           const interactions::ClampSet& clamps, interactions::ClampReport* clamp_report,
                           const CheckpointHooks* ckpt = nullptr,
                           interactions::NormalProjectionStats* proj_stats = nullptr)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    const size_t nfiber = fb.size();

    // Проекция скорости нормального контакта (законы 1-3, fiber/TZ_LAW3_PROJECTION.md):
    // стадия между выбором dt и стадией трения. Оба ключа 0 -- proj_on == false, ни одна
    // строка ниже не выполняется, результат побитово прежний.
    validate_fiber_contact_projection(cfg);
    const bool proj_ff = cfg.fiber_fiber_contact_predictor == 1;
    const bool proj_fs = cfg.fiber_solid_contact_predictor == 1;
    const bool proj_on = nfiber > 0 && (proj_ff || proj_fs);

    const double t_total = cfg.t_total;
    const size_t print_step = cfg.print_step;
    const double save_t = cfg.save_t;
    const size_t moni_particle = cfg.moni_particle;

    // Solid RHS-буферы охватывают реальные И виртуальные частицы: single_step
    // вызывает узлы конвейера с ntotal + nvirt, и int_force/art_visc/av_vel
    // индексируют эти массивы всем диапазоном. Число зеркальных частиц
    // определяется каждый шаг заново, поэтому берётся запас maxn (main.cpp:
    // maxn = ntotal * 1.4 «под возможные виртуальные частицы»). Интегрируются
    // по-прежнему только реальные частицы [0, ntotal).
    const size_t nbuf = (cfg.maxn > ntotal) ? cfg.maxn : ntotal;

    Solver solver;
    solver.reset(nbuf, ntotal, nfiber);

    // переиспользуемые рабочие буферы (живут всю симуляцию, передаются в single_step)
    WorkBuffers wb;

    // --- Fiber RHS-инфраструктура (fiber/PLAN.md, этапы Ф3/Ф5/Ф6/Ф7):
    // per-thread буферы contributions/событий, переиспользуются между шагами
    // ("память только растёт", как и остальные буферы RhsWorkspace). Не
    // относятся к состоянию решателя (не f^{n-1}) -- это инфраструктура
    // ВЫЧИСЛЕНИЯ ОДНОГО RHS, а не история интегрирования по времени.
    rhs::CpuExecutor fiber_exec;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_thread_buffers;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_contact_thread_buffers;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_friction_thread_buffers;
    // Mass splitting стадии трения: c_k (число активных слотов на бусину) и
    // буферы подсчёта -- переиспользуются между шагами (память только растёт).
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> fiber_friction_count_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidCountDelta>> solid_friction_count_buffers;
    std::vector<double> fiber_friction_bead_count; // c_k (оба вида контактов)
    std::vector<double> solid_friction_count;      // c_s
    // Вклады стадии трения в Solid (interaction_type 13) и слоты Fiber-Solid.
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> solid_friction_thread_buffers;
    std::vector<interactions::FiberSolidFrictionContact> friction_solid_contacts;
    // Контакты нить-нить текущего шага и слоты стадии трения (fiber/
    // fiber_friction_projection.md, §3, §5): живут на уровне функции, потому
    // что стадия трения выполняется ПОСЛЕ вычисления dt, вне блока RHS.
    // Переиспользуются между шагами (память только растёт).
    interactions::FiberContactSearchResult fiber_contact_search;
    std::vector<interactions::FiberFrictionContact> friction_contacts;
    // Списки активных слотов трения (режим списка: ядра сбрасывают только слоты прошлого шага, стадия
    // трения обходит только активные; см. FrictionActiveRef) и рабочая память их сбора по отрезкам.
    // Списки живут между шагами вместе с массивами слотов; при рестарте из чекпоинта слоты не
    // восстанавливаются (не сериализуются) -- пустые списки согласованы со свежими слотами.
    std::vector<interactions::FrictionActiveRef> fric_ff_active;
    std::vector<interactions::FrictionActiveRef> fric_fs_active;
    std::vector<rhs::CpuExecutor::SegPart<interactions::FrictionActiveRef>> fric_ff_parts;
    std::vector<rhs::CpuExecutor::SegPart<interactions::FrictionActiveRef>> fric_fs_parts;
    // Нормальные слоты и рабочая память стадии проекции (заполняются ядрами только при
    // predictor = 1 соответствующего вида; иначе остаются пустыми) и таблица тёплого старта.
    std::vector<interactions::FiberNormalContact> normal_ff;
    std::vector<interactions::FiberSolidNormalContact> normal_fs;
    std::vector<rhs::CpuExecutor::SegPart<interactions::FiberNormalContact>> normal_ff_parts;      // сбор по отрезкам
    std::vector<rhs::CpuExecutor::SegPart<interactions::FiberSolidNormalContact>> normal_fs_parts;
    interactions::NormalProjectionWorkspace proj_ws;
    interactions::NormalWarmStart proj_warm;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fiber_solid_fiber_buffers;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fiber_solid_solid_buffers;
    std::vector<std::vector<interactions::FiberBreakEvent>> fiber_break_thread_buffers;
    // Результаты reduce и его рабочая память (память только растёт): одна
    // ReduceScratch на все вызовы шага, они выполняются строго последовательно.
    rhs::ReduceScratch fiber_reduce_scratch;
    std::vector<fiber::FiberDelta> fiber_bond_delta;
    std::vector<fiber::FiberDelta> fiber_contact_delta;
    std::vector<fiber::FiberDelta> fiber_solid_fiber_delta;
    std::vector<interactions::SolidDelta> fiber_solid_solid_delta;
    std::vector<fiber::FiberDelta> friction_delta;
    std::vector<interactions::SolidDelta> friction_solid_delta;
    std::vector<interactions::FiberCountDelta> fiber_friction_count_delta;
    std::vector<interactions::SolidCountDelta> solid_friction_count_delta;
    // Поиск Fiber-Solid-Contact текущего шага: primitives читаются dt_limit,
    // рядом живёт рабочая память поиска (индекс, пары, копии координат).
    interactions::FiberSolidContactSearchResult fiber_solid_search;
    // Законы нормальной силы контактов нитей (interactions/fiber_contact_law.hpp):
    // 0 = Герц (прежний путь); dt_bead нужен только при типе != 0.
    const interactions::FiberContactLaw fiber_solid_law{cfg.fiber_solid_contact_type,
                                                        cfg.fiber_solid_contact_penalty_scale,
                                                        cfg.fiber_solid_contact_impulse_gamma,
                                                        cfg.fiber_solid_contact_impulse_push_tau,
                                                        cfg.fiber_contact_impulse_push_depth,
                                                        cfg.fiber_contact_impulse_push_vmax};
    const interactions::FiberContactLaw fiber_fiber_law{cfg.fiber_fiber_contact_type,
                                                        cfg.fiber_fiber_contact_penalty_scale,
                                                        cfg.fiber_fiber_contact_impulse_gamma,
                                                        cfg.fiber_fiber_contact_impulse_push_tau,
                                                        cfg.fiber_contact_impulse_push_depth,
                                                        cfg.fiber_contact_impulse_push_vmax};
    const bool fiber_need_bead_dt_ref = interactions::fiber_contact_law_is_penalty(fiber_solid_law) ||
                                        interactions::fiber_contact_law_is_penalty(fiber_fiber_law);
    // Кэш опорных пределов шага нитей (bond_dt_limit, bead_dt_ref): пересчёт только при смене топологии.
    // Изначально invalid (в т.ч. после рестарта: топология восстановлена до цикла, кэш локален).
    interactions::FiberDtRefCache dt_cache;

    // Файл мониторинга открывается ОДИН раз (дозапись; до цикла log.txt уже
    // содержит строки main.cpp/input.cpp); std::endl сбрасывает буфер каждый шаг.
    std::ofstream log_file("log.txt", std::ios::app);

    // Кинематическое закрепление (interactions/clamp.hpp): проверки, списки
    // зажатых и рабочая память -- ОДИН раз до цикла (в шаге -- без аллокаций).
    const bool has_clamps = !clamps.empty();
    interactions::ClampWorkspace clamp_ws;
    if (has_clamps)
    {
        interactions::validate_clamps(clamps, sys, ntotal, fb);
        interactions::prepare_clamp_workspace(clamps, ntotal, nfiber, clamp_ws);
        if (clamp_report != nullptr) interactions::clamp_report_counts(clamp_ws, sys, *clamp_report);
    }

    // Начальное eps_prev неупругих связей = текущая деформация (иначе связь,
    // созданная при l != l0, даёт на первом шаге фиктивный скачок Δε и w_v).
    // Для упругих материалов no-op.
    // При рестарте eps_prev берётся из чекпоинта (BondState восстановлен), не пересчитывается.
    const checkpoint::Payload* const restart = (ckpt != nullptr) ? ckpt->restart : nullptr;
    if (nfiber > 0 && restart == nullptr) interactions::init_fiber_bond_state(fb, topology, materials);

    size_t current_ts = 0;
    double next_output_time = 0.0;

    // Рестарт: состояние цикла, история решателя, счётчики сторожа. dt и time --
    // ссылки вызывающего, значения из чекпоинта перекрывают пересчитанные снаружи.
    if (restart != nullptr)
    {
        const auto& L = restart->loop;
        time = L.time;
        dt = L.dt;
        current_ts = L.step;
        next_output_time = L.next_output_time;
        const auto& H = restart->hist;
        restore_history_vec(solver.solid_rhs.du, H.du);
        restore_history_vec(solver.solid_rhs.drho, H.drho);
        restore_history_vec(solver.solid_rhs.dsxx, H.dsxx);
        restore_history_vec(solver.solid_rhs.dsyy, H.dsyy);
        restore_history_vec(solver.solid_rhs.dszz, H.dszz);
        restore_history_vec(solver.solid_rhs.dsxy, H.dsxy);
        restore_history_vec(solver.solid_rhs.dsyz, H.dsyz);
        restore_history_vec(solver.solid_rhs.dsxz, H.dsxz);
        for (size_t d = 0; d < dim; ++d)
        {
            restore_history_vec(solver.solid_rhs.dvx[d], H.dvx[d]);
            restore_history_vec(solver.fiber_rhs.dv[d], H.fiber_dv[d]);
        }
        solver.set_first_step(H.first_step);
        if (guard != nullptr) guard->restore_progress(L.guard_steps, L.guard_dt_first, L.guard_dt_min);
        proj_warm.set(restart->proj_warm, restart->proj_h_prev);
    }

    // Периодические и финальные чекпоинты. Выключено -- ckpt_on == false, в цикле одна ветка.
    const bool ckpt_on = ckpt != nullptr && static_cast<bool>(ckpt->on_checkpoint) && cfg.checkpoint_every_steps > 0;
    size_t last_ckpt_ts = current_ts; // шаг последнего записанного (или восстановленного) чекпоинта

#ifdef SPH_ENERGY_LOG
    // Лог механической энергии (include/energy_log.hpp): запасы -- каждые
    // energy_log_every шагов, каналы потерь -- накапливаются на КАЖДОМ шаге.
    // При energy_log_every == 0 elog пуст, и все врезки ниже -- одна ветка.
    // Врезки только читают состояние: траектория побитово та же, что без лога.
    std::unique_ptr<energy::EnergyLogger> elog;
    energy::EnergyChannels echan;
    interactions::FiberEnergyTap ftap; // приращения энергии нитей (postprocessing.hpp)
    std::vector<double> sigma_clip;    // σ_m < 0 на отсечке по связям шага (RHS-ядро -> clip_loss)
    // Проход «только демпфер» закона 4 для канала D_contact_fiber (рабочая память).
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> elog_damp_ff_buffers, elog_damp_fs_fbuf;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> elog_damp_fs_sbuf;
    std::vector<interactions::FiberFrictionContact> elog_damp_ff_slots;
    std::vector<interactions::FiberSolidFrictionContact> elog_damp_fs_slots;
    std::vector<fiber::FiberDelta> elog_damp_ff_delta, elog_damp_fs_fdelta;
    std::vector<interactions::SolidDelta> elog_damp_fs_sdelta;
    // Стадия проекции: ускорения до стадии (для разности -> канал D_contact_fiber).
    std::vector<double> proj_dv_fiber_before[dim], proj_dv_solid_before[dim];
    if (cfg.energy_log_every > 0)
    {
        elog = std::make_unique<energy::EnergyLogger>(cfg.energy_log_file, cfg.energy_log_every, restart != nullptr);
        solver.u_clamp_acc = &echan.u_clamp;
        ftap.w_p = &echan.w_p_fiber;
        ftap.w_v = &echan.w_v_fiber;
        ftap.brk = &echan.brk;
        ftap.brk_contact = &echan.brk_contact;
        ftap.clip = &echan.clip;
        ftap.sigma_clip = &sigma_clip;
        ftap.contact_potential = [&]() {
            double uff = 0.0, ufs = 0.0;
            energy::compute_hertz_potential(fiber_exec, cfg, sys, ntotal, fb, topology, materials, elog->workspace(),
                                            uff, ufs);
            return uff + ufs;
        };
        wb.energy_log = true;
        if (restart != nullptr)
        {
            // Рестарт: CSV дописывается без заголовка и строки шага 0; каналы и база E_mech0 -- из чекпоинта.
            if (restart->has_energy)
            {
                echan.from_array(restart->energy.chan);
                elog->restore_e0(restart->energy.e_mech0, restart->energy.have_e0);
            }
        }
        else
            elog->write_row(0, time, dt,
                            energy::compute_energy_stocks(fiber_exec, cfg, sys, ntotal, fb, topology, materials,
                                                          elog->workspace()),
                            echan);
    }
#endif

    // Снимок живого состояния решателя/цикла для checkpoint::write (вызывается редко).
    auto make_payload = [&]() {
        checkpoint::Payload pl;
        pl.loop.time = time;
        pl.loop.dt = dt;
        pl.loop.step = current_ts;
        pl.loop.next_output_time = next_output_time;
        pl.proj_warm = proj_warm.entries();
        pl.proj_h_prev = proj_warm.h_prev();
        if (guard != nullptr)
        {
            pl.loop.guard_steps = guard->steps();
            pl.loop.guard_dt_first = guard->dt_first();
            pl.loop.guard_dt_min = guard->dt_min();
        }
        auto& H = pl.hist;
        H.first_step = solver.first_step();
        auto cp = [ntotal](std::vector<double>& dst, const std::vector<double>& src) {
            dst.assign(src.begin(), src.begin() + static_cast<std::ptrdiff_t>(ntotal));
        };
        cp(H.du, solver.solid_rhs.du);
        cp(H.drho, solver.solid_rhs.drho);
        cp(H.dsxx, solver.solid_rhs.dsxx);
        cp(H.dsyy, solver.solid_rhs.dsyy);
        cp(H.dszz, solver.solid_rhs.dszz);
        cp(H.dsxy, solver.solid_rhs.dsxy);
        cp(H.dsyz, solver.solid_rhs.dsyz);
        cp(H.dsxz, solver.solid_rhs.dsxz);
        for (size_t d = 0; d < dim; ++d)
        {
            cp(H.dvx[d], solver.solid_rhs.dvx[d]);
            H.fiber_dv[d].assign(solver.fiber_rhs.dv[d].begin(),
                                 solver.fiber_rhs.dv[d].begin() + static_cast<std::ptrdiff_t>(nfiber));
        }
#ifdef SPH_ENERGY_LOG
        if (elog)
        {
            pl.has_energy = true;
            echan.to_array(pl.energy.chan);
            pl.energy.e_mech0 = elog->e_mech0();
            pl.energy.have_e0 = elog->have_e0();
        }
#endif
        return pl;
    };
    auto write_checkpoint = [&]() {
        ckpt->on_checkpoint(make_payload());
        last_ckpt_ts = current_ts;
    };

#ifdef SPH_CONTACT_STATS
    // Диагностика контактов нитей (interactions/fiber_contact_stats.hpp): без
    // contact_hold_n только читает состояние; сводка пишется в деструкторе.
    std::unique_ptr<interactions::contact_stats::Collector> cstats;
    if (cfg.contact_stats_every > 0 && nfiber > 0)
    {
        if (proj_on)
            throw std::runtime_error(
                "contact_stats не поддерживается в режиме проекции нормального контакта (fiber_*_contact_predictor = 1)");
        cstats = std::make_unique<interactions::contact_stats::Collector>(
            cfg, cfg.contact_stats_every, cfg.contact_stats_file, cfg.contact_stats_summary, cfg.contact_hold_n);
    }
#endif

    // Режим списка слотов трения включён на всех шагах расчёта; контактная статистика читает плотные
    // массивы слотов (и меняет их на месте) -- с ней список отключается.
    bool fric_list_mode = true;
#ifdef SPH_CONTACT_STATS
    if (cstats) fric_list_mode = false;
#endif

    // Основной цикл по времени
    while (time < t_total)
    {
        current_ts += 1;
        // Вывод прогресса
        if (current_ts % print_step == 0)
        {
            std::wcout << "______________________________________" << std::endl;
            std::wcout << "  Текущий щаг по времени = " << dt
                      << "    текущее время расчета = " << time + dt
                      << "    текущий цикл = " << current_ts << std::endl;
            if (cfg.contact_force_type != 0)
            {
                const bool limited_by_contact = (wb.min_dt_contact < 1e30) &&
                    (wb.min_dt_contact <= compute_cfl_timestep(cfg, sys, ntotal));
                std::wcout << "  Шаг ограничен " << (limited_by_contact ? "КОНТАКТОМ" : "CFL")
                          << "    contact_pair_count = " << wb.contact_pair_count
                          << "    contact_max_delta = " << wb.contact_max_delta << std::endl;
            }
            std::wcout << "______________________________________" << std::endl;
        }

        // OdeSolver: предиктор (State::save + полушаг по f^{n-1}; нет-оп у
        // Эйлера и на первом шаге PredictorCorrector).
        solver.predict(fiber_exec, sys, ntotal, fb, nfiber, dt);

        // RHS: Solid (SPH+контакт, как раньше) -- пишет f^n в solver.solid_rhs.
        single_step(cfg, sys, current_ts, dt, ntotal, solver.solid_rhs.t, solver.solid_rhs.tdsdt, solver.solid_rhs.dx,
                    solver.solid_rhs.dvx, solver.solid_rhs.du, solver.solid_rhs.ds, solver.solid_rhs.drho,
                    solver.solid_rhs.av, solver.solid_rhs.dsxx, solver.solid_rhs.dsyy, solver.solid_rhs.dszz,
                    solver.solid_rhs.dsxy, solver.solid_rhs.dsyz, solver.solid_rhs.dsxz, wb);

        // Postprocessing, этап Ф7 (fiber/PLAN.md §3.8): sle=1 -- пересчёт h.
        // МЕСТО ВЫЗОВА -- сразу после single_step(), ДО compute_cfl_timestep():
        // см. большой комментарий у interactions::update_smoothing_length
        // (include/interactions/postprocessing.hpp) -- почему это обязательно
        // для бит-в-бит совпадения с прежним поведением (h_upgrade раньше
        // вызывался отсюда же, изнутри single_step.cpp).
        interactions::update_smoothing_length(cfg, sys, dt, ntotal);

        // RHS: Fiber-Bond (fiber/PLAN.md, этап Ф3) — readonly search по
        // topology.bonds, contributions -> reduce (fiber/SPH.md, «Основные
        // требования к interaction kernel»).
        // RHS: Fiber-Fiber-Contact (fiber/PLAN.md, этап Ф5) — readonly
        // spatial search по бусинам РАЗНЫХ нитей (interactions/fiber_contact.hpp),
        // тот же FiberDelta reduce, СКЛАДЫВАЕТСЯ с вкладом Fiber-Bond (обе
        // физики независимо производят ускорение бусин — fiber/PLAN.md §3.3:
        // «FiberDelta: dv[dim]», сумма нескольких interaction types в общий dv).
        double fiber_dt_limit = std::numeric_limits<double>::max();
        bool has_active_friction_contacts = false;      // нить-нить
        bool has_active_solid_friction_contacts = false; // Fiber-Solid
        if (nfiber > 0)
        {
#ifdef SPH_ENERGY_LOG
            interactions::compute_fiber_bond(fiber_exec, fb, topology, materials, fiber_thread_buffers,
                                             fiber_bond_delta, fiber_reduce_scratch, elog ? &sigma_clip : nullptr);
#else
            interactions::compute_fiber_bond(fiber_exec, fb, topology, materials, fiber_thread_buffers,
                                             fiber_bond_delta, fiber_reduce_scratch);
#endif
            dt_cache.ensure(fiber_exec, fb, topology, materials, fiber_need_bead_dt_ref);

            // L_max (текущая макс. длина связи) -- один раз за шаг для обоих поисков (fb.x, topology не меняются).
            const double fiber_L_max = interactions::fiber_max_current_bond_length(fiber_exec, fb, topology);
            // Режим отсечки: в primitives остаются только пары, не прошедшие fiber_contact_surely_apart; ядро
            // и предел шага работают по виду (индекс k полного списка -- в prim_k), результат побитово тот же.
            interactions::FiberContactBuildOptions fiber_build_opt;
            fiber_build_opt.cull_apart = true;
            fiber_build_opt.l_max = &fiber_L_max;
            interactions::build_fiber_contact_search(fiber_exec, fb, topology, materials, fiber_contact_search,
                                                     fiber_build_opt);
            const interactions::FiberContactPrimView fiber_pview =
                interactions::fiber_contact_view(fiber_contact_search);
            // dt здесь ещё шаг ПРОШЛОГО цикла (новый считается после RHS) --
            // тип 3 строит силу по нему, как Solid-Solid.
            interactions::FiberContactStepParams fiber_fiber_params;
            fiber_fiber_params.law = fiber_fiber_law;
            fiber_fiber_params.dt_prev = dt;
            fiber_fiber_params.cfg = &cfg;
            fiber_fiber_params.bead_dt_ref = &dt_cache.bead_dt_ref;
            fiber_fiber_params.pair_dt_limit = &fiber_contact_search.pair_dt_limit;
            fiber_fiber_params.normal_ff = proj_ff ? &normal_ff : nullptr; // проекция: слоты вместо силы
            fiber_fiber_params.normal_ff_parts = &normal_ff_parts;
            // Закон 0: предел шага отдаёт само ядро (тот же результат, что compute_fiber_contact_dt_limit).
            double fiber_fiber_law0_dt = std::numeric_limits<double>::max();
            if (fiber_fiber_law.type == 0) fiber_fiber_params.law0_dt_min_out = &fiber_fiber_law0_dt;
            if (fric_list_mode)
            {
                fiber_fiber_params.fric_active = &fric_ff_active;
                fiber_fiber_params.fric_active_parts = &fric_ff_parts;
            }
            // Только нормальная сила; заполняет слоты трения friction_contacts.
            has_active_friction_contacts = interactions::compute_fiber_contact_contributions(
                fiber_exec, fb, materials, fiber_pview, fiber_contact_thread_buffers,
                friction_contacts, fiber_fiber_params);
            rhs::reduce_contributions_into(fiber_exec, fiber_contact_thread_buffers, fb.size(), fiber_contact_delta,
                                           fiber_reduce_scratch);

            // Поэлементные циклы сборки RHS -- параллельно (каждый элемент пишется один раз, то же выражение).
            for (size_t d = 0; d < dim; ++d)
                fiber_exec.parallel_for(nfiber, [&](std::size_t i) {
                    solver.fiber_rhs.dv[d][i] = fiber_bond_delta[i].dv[d] + fiber_contact_delta[i].dv[d];
                });

            fiber_dt_limit = dt_cache.bond_dt_limit;
            const double fiber_contact_dt_limit =
                (fiber_fiber_law.type == 0)
                    ? fiber_fiber_law0_dt
                    : interactions::compute_fiber_contact_dt_limit(fiber_exec, fb, materials, fiber_pview,
                                                                   fiber_fiber_params);
            if (fiber_contact_dt_limit < fiber_dt_limit) fiber_dt_limit = fiber_contact_dt_limit;

            // RHS: Fiber-Solid-Contact (fiber/PLAN.md, этап Ф6) -- механический
            // контакт бусины нити с реальной Solid-частицей. Вклад в Solid идёт
            // ЧЕРЕЗ ТОТ ЖЕ МЕХАНИЗМ, что и Solid-SPH/Solid-Solid-Contact внутри
            // single_step (SolidDelta + rhs::reduce_contributions) — здесь он
            // лишь ДОБАВЛЯЕТСЯ к уже посчитанному single_step() dvx
            // (interactions/fiber_solid_contact.hpp, шапка файла, решения 1-6).
            // Вклад в Fiber складывается с Fiber-Bond/Fiber-Fiber-Contact в тот
            // же dv, как и они между собой.
            // dt здесь ещё шаг ПРОШЛОГО цикла (новый считается после RHS) --
            // тип 3 строит силу по нему, как Solid-Solid.
            interactions::FiberContactStepParams fiber_solid_params;
            fiber_solid_params.law = fiber_solid_law;
            fiber_solid_params.dt_prev = dt;
            fiber_solid_params.cfg = &cfg;
            fiber_solid_params.bead_dt_ref = &dt_cache.bead_dt_ref;
            fiber_solid_params.pair_dt_limit = &fiber_solid_search.pair_dt_limit;
            fiber_solid_params.normal_fs = proj_fs ? &normal_fs : nullptr; // проекция: слоты вместо силы
            fiber_solid_params.normal_fs_parts = &normal_fs_parts;
            double fiber_solid_law0_dt = std::numeric_limits<double>::max();
            if (fiber_solid_law.type == 0) fiber_solid_params.law0_dt_min_out = &fiber_solid_law0_dt;
            if (fric_list_mode)
            {
                fiber_solid_params.fric_active = &fric_fs_active;
                fiber_solid_params.fric_active_parts = &fric_fs_parts;
            }
            has_active_solid_friction_contacts = interactions::compute_fiber_solid_contact(
                fiber_exec, cfg.mats, sys, ntotal, fb, topology, materials, fiber_solid_fiber_buffers,
                fiber_solid_solid_buffers, fiber_solid_search, fiber_solid_fiber_delta, fiber_solid_solid_delta,
                fiber_reduce_scratch, &friction_solid_contacts, fiber_solid_params, &fiber_L_max);

            for (size_t d = 0; d < dim; ++d)
                fiber_exec.parallel_for(nfiber,
                                        [&](std::size_t i) { solver.fiber_rhs.dv[d][i] += fiber_solid_fiber_delta[i].dv[d]; });
            for (size_t d = 0; d < dim; ++d)
                fiber_exec.parallel_for(
                    ntotal, [&](std::size_t i) { solver.solid_rhs.dvx[d][i] += fiber_solid_solid_delta[i].dv[d]; });

            const double fiber_solid_dt_limit =
                (fiber_solid_law.type == 0)
                    ? fiber_solid_law0_dt
                    : interactions::compute_fiber_solid_contact_dt_limit(fiber_exec, cfg.mats, sys, fb, materials,
                                                                         fiber_solid_search.primitives,
                                                                         fiber_solid_params);
            if (fiber_solid_dt_limit < fiber_dt_limit) fiber_dt_limit = fiber_solid_dt_limit;
#ifdef SPH_ENERGY_LOG
            // Ускорение непотенциальной части контактов нитей -- для канала
            // D_contact_fiber: законы 1-3 -- вся сила; закон 4 -- только демпфер
            // (отдельный проход ядра с damper_only; упругая часть -- в U_hertz_*);
            // у Герца (0) канал не нужен.
            if (elog && (fiber_fiber_law.type != 0 || fiber_solid_law.type != 0))
            {
                const int tff = fiber_fiber_law.type, tfs = fiber_solid_law.type;
                const bool d4ff = tff == interactions::kFiberLawHertzDamped;
                const bool d4fs = tfs == interactions::kFiberLawHertzDamped;
                if (d4ff)
                {
                    interactions::FiberContactStepParams p = fiber_fiber_params;
                    p.damper_only = true;
                    p.pair_dt_limit = nullptr;
                    p.law0_dt_min_out = nullptr;
                    p.fric_active = nullptr;
                    p.fric_active_parts = nullptr;
                    interactions::compute_fiber_contact_contributions(fiber_exec, fb, materials, fiber_pview,
                                                                      elog_damp_ff_buffers, elog_damp_ff_slots, p);
                    rhs::reduce_contributions_into(fiber_exec, elog_damp_ff_buffers, fb.size(), elog_damp_ff_delta,
                                                   fiber_reduce_scratch);
                }
                if (d4fs)
                {
                    interactions::FiberContactStepParams p = fiber_solid_params;
                    p.damper_only = true;
                    p.pair_dt_limit = nullptr;
                    p.law0_dt_min_out = nullptr;
                    p.fric_active = nullptr;
                    p.fric_active_parts = nullptr;
                    interactions::compute_fiber_solid_contact_contributions(
                        fiber_exec, cfg.mats, sys, fb, materials, fiber_solid_search.primitives, elog_damp_fs_fbuf,
                        elog_damp_fs_sbuf, elog_damp_fs_slots, p);
                    rhs::reduce_contributions_into(fiber_exec, elog_damp_fs_fbuf, fb.size(), elog_damp_fs_fdelta,
                                                   fiber_reduce_scratch);
                    rhs::reduce_contributions_into(fiber_exec, elog_damp_fs_sbuf, ntotal, elog_damp_fs_sdelta,
                                                   fiber_reduce_scratch);
                }
                const bool ff = tff != 0, fs = tfs != 0;
                for (size_t d = 0; d < dim; ++d)
                {
                    elog->contact_fiber_acc[d].assign(nfiber, 0.0);
                    for (size_t i = 0; i < nfiber; ++i)
                        elog->contact_fiber_acc[d][i] =
                            (ff ? (d4ff ? elog_damp_ff_delta[i].dv[d] : fiber_contact_delta[i].dv[d]) : 0.0) +
                            (fs ? (d4fs ? elog_damp_fs_fdelta[i].dv[d] : fiber_solid_fiber_delta[i].dv[d]) : 0.0);
                }
                elog->contact_solid_acc.resize(dim);
                for (size_t d = 0; d < dim; ++d)
                {
                    elog->contact_solid_acc[d].assign(ntotal, 0.0);
                    if (fs)
                        for (size_t i = 0; i < ntotal; ++i)
                            elog->contact_solid_acc[d][i] =
                                d4fs ? elog_damp_fs_sdelta[i].dv[d] : fiber_solid_solid_delta[i].dv[d];
                }
                elog->has_contact_fiber = true;
            }
#endif
#ifdef SPH_CONTACT_STATS
            if (cstats)
                cstats->after_rhs(static_cast<long long>(current_ts), cfg.mats, sys, fb, topology, materials,
                                  fiber_contact_search, friction_contacts, fiber_fiber_law, fiber_solid_search,
                                  friction_solid_contacts, fiber_solid_law, dt, has_active_friction_contacts,
                                  has_active_solid_friction_contacts);
#endif
        }

        // [1] Закрепление: RHS шага полностью сведён (все Reduce и сложения
        // выше завершены) -- проекция a_{i,d} <- 0 (у Solid ещё и XSPH av)
        // для зажатых компонент; реакция -- до обнуления. Стоит ДО dt и ДО
        // стадии трения: пробная скорость зажатой частицы в трении = 0.
        double clamp_solid_force[dim] = {};
        double clamp_fiber_force[dim] = {};
        if (has_clamps)
        {
            interactions::apply_clamp_projection(fiber_exec, clamp_ws.solid, sys.mass, solver.solid_rhs.dvx,
                                                 &solver.solid_rhs.av, clamp_ws, clamp_solid_force,
                                                 &clamp_ws.acc_solid, /*acc_add=*/false);
            interactions::apply_clamp_projection(fiber_exec, clamp_ws.fiber, fb.mass, solver.fiber_rhs.dv,
                                                 static_cast<std::vector<std::vector<double>>*>(nullptr), clamp_ws,
                                                 clamp_fiber_force, &clamp_ws.acc_fiber, /*acc_add=*/false);
#ifdef SPH_ENERGY_LOG
            // Ускорение контакта в зажатых компонентах тоже обнуляется -- иначе
            // в канале контакта появилась бы работа неподвижной опоры.
            if (wb.has_contact_dv)
            {
                double unused_force[dim] = {};
                interactions::apply_clamp_projection(fiber_exec, clamp_ws.solid, sys.mass, wb.contact_dv,
                                                     static_cast<std::vector<std::vector<double>>*>(nullptr),
                                                     clamp_ws, unused_force);
            }
            // То же для ускорения контактов нитей с законом 1-3 (канал D_contact_fiber).
            if (elog && elog->has_contact_fiber)
            {
                double unused_force[dim] = {};
                interactions::apply_clamp_projection(fiber_exec, clamp_ws.solid, sys.mass, elog->contact_solid_acc,
                                                     static_cast<std::vector<std::vector<double>>*>(nullptr),
                                                     clamp_ws, unused_force);
                interactions::apply_clamp_projection(fiber_exec, clamp_ws.fiber, fb.mass, elog->contact_fiber_acc,
                                                     static_cast<std::vector<std::vector<double>>*>(nullptr),
                                                     clamp_ws, unused_force);
            }
#endif
        }

        // Общий шаг (fiber/SPH.md, «Шаг по времени»): min(CFL, контакт Solid, dt_limit связей).
        dt = compute_cfl_timestep(cfg, sys, ntotal);
        if (wb.min_dt_contact < dt) dt = wb.min_dt_contact;
        if (fiber_dt_limit < dt) dt = fiber_dt_limit;
        if (cfg.dt_scale < 1.0) dt *= cfg.dt_scale; // проверка сходимости по шагу

        // Стадия нормального контакта законов 1-3 как ПРОЕКЦИЯ СКОРОСТИ
        // (include/interactions/fiber_normal_projection.hpp): ядра RHS для видов с
        // predictor = 1 не дали силы, а заполнили нормальные слоты. Стадия стоит после
        // выбора dt (нужен h) и ДО трения: добавляет Δv/h в solver.fiber_rhs.dv и
        // solver.solid_rhs.dvx (после зажимов [1]); затем слоты трения проецируемых видов
        // заполняются из нормальных (N = P/h), слоты трения остальных видов не трогаются.
        // Работа стадии идёт в канал D_contact_fiber лога энергии.
        if (proj_on)
        {
            const auto fupd = solver.fiber_velocity_update(fb, dt);
            const auto supd = solver.solid_velocity_update(sys, dt);
            const double h_proj = fupd.h;
#ifdef SPH_ENERGY_LOG
            if (elog)
            {
                for (size_t d = 0; d < dim; ++d)
                {
                    proj_dv_fiber_before[d].assign(solver.fiber_rhs.dv[d].begin(),
                                                   solver.fiber_rhs.dv[d].begin() + static_cast<std::ptrdiff_t>(nfiber));
                    proj_dv_solid_before[d].assign(solver.solid_rhs.dvx[d].begin(),
                                                   solver.solid_rhs.dvx[d].begin() + static_cast<std::ptrdiff_t>(ntotal));
                }
            }
#endif
            proj_warm.apply(fiber_exec, normal_ff, normal_fs, h_proj, cfg.fiber_contact_projection_warm == 1);
            interactions::compute_normal_counts(fiber_exec, nfiber, ntotal, normal_ff, normal_fs, proj_ws);
            if (has_clamps)
                interactions::apply_clamp_friction_counts(fiber_exec, clamp_ws, proj_ws.bead_count,
                                                          proj_ws.solid_count);
            interactions::NormalProjectionStats stage_stats;
            interactions::run_normal_projection(fiber_exec, fb, fupd, solver.fiber_rhs.dv, normal_ff, sys, supd,
                                                solver.solid_rhs.dvx, normal_fs, cfg.fiber_contact_projection_iters,
                                                proj_ws, stage_stats);
            if (proj_stats != nullptr) *proj_stats = stage_stats;
            // [3] Вклад стадии в зажатые компоненты -- часть реакции: повторная проекция.
            if (has_clamps)
            {
                interactions::apply_clamp_projection(fiber_exec, clamp_ws.solid, sys.mass, solver.solid_rhs.dvx,
                                                     static_cast<std::vector<std::vector<double>>*>(nullptr),
                                                     clamp_ws, clamp_solid_force, &clamp_ws.acc_solid,
                                                     /*acc_add=*/true);
                interactions::apply_clamp_projection(fiber_exec, clamp_ws.fiber, fb.mass, solver.fiber_rhs.dv,
                                                     static_cast<std::vector<std::vector<double>>*>(nullptr),
                                                     clamp_ws, clamp_fiber_force, &clamp_ws.acc_fiber,
                                                     /*acc_add=*/true);
            }
#ifdef SPH_ENERGY_LOG
            // Фактически применённое ускорение стадии -- в канал D_contact_fiber (после проекции
            // зажимов). Для проецируемых видов ядра RHS силу не давали, поэтому прибавка
            // к уже снятому ускорению остальных контактов нитей.
            if (elog && elog->has_contact_fiber)
            {
                for (size_t d = 0; d < dim; ++d)
                {
                    for (size_t i = 0; i < nfiber; ++i)
                        elog->contact_fiber_acc[d][i] += solver.fiber_rhs.dv[d][i] - proj_dv_fiber_before[d][i];
                    for (size_t i = 0; i < ntotal; ++i)
                        elog->contact_solid_acc[d][i] += solver.solid_rhs.dvx[d][i] - proj_dv_solid_before[d][i];
                }
            }
#endif
            proj_warm.rebuild(normal_ff, normal_fs, h_proj);
            // Слоты трения проецируемых видов из нормальных; остальные вид -- как дало ядро.
            if (proj_ff)
                has_active_friction_contacts = interactions::fill_friction_slots_from_normal_ff(
                    fiber_exec, h_proj, normal_ff, friction_contacts, fric_list_mode ? &fric_ff_active : nullptr);
            if (proj_fs)
                has_active_solid_friction_contacts = interactions::fill_friction_slots_from_normal_fs(
                    fiber_exec, h_proj, normal_fs, friction_solid_contacts, fric_list_mode ? &fric_fs_active : nullptr);
        }

        // Единая стадия трения (fiber/fiber_friction_projection.md, §5;
        // fiber/fiber_solid_friction_implementation.md, §5): проекция на
        // уровне скоростей по формуле обновления скорости решателя для ОБОИХ
        // видов контактов (нить-нить и Fiber-Solid) ОДНИМ проходом Якоби, с
        // mass splitting (c_k по обоим видам, c_s для частиц Solid).
        // Добавляет ускорение трения в solver.fiber_rhs.dv и solver.solid_rhs.dvx
        // (после этого там ПОЛНОЕ ускорение шага -- оно же станет f^{n-1} для
        // предиктора PC). При mu_ff=mu_fs=0 у всех материалов активных слотов
        // нет -- стадия не выполняется вообще (побитово прежний результат).
        if (nfiber > 0 && (has_active_friction_contacts || has_active_solid_friction_contacts))
        {
            const auto fupd = solver.fiber_velocity_update(fb, dt);
            const auto supd = solver.solid_velocity_update(sys, dt);
#ifdef SPH_ENERGY_LOG
            // Трение: снимок ускорений до стадии; после неё (с повторной
            // проекцией зажимов [3]) разность -- фактически применённое
            // ускорение трения, его работа считается перед advance().
            if (elog)
            {
                for (size_t d = 0; d < dim; ++d)
                    elog->friction_fiber[d].assign(solver.fiber_rhs.dv[d].begin(),
                                                   solver.fiber_rhs.dv[d].begin() + static_cast<std::ptrdiff_t>(nfiber));
                elog->friction_solid.resize(dim);
                for (size_t d = 0; d < dim; ++d)
                    elog->friction_solid[d].assign(solver.solid_rhs.dvx[d].begin(),
                                                   solver.solid_rhs.dvx[d].begin() + static_cast<std::ptrdiff_t>(ntotal));
            }
#endif
            interactions::compute_friction_counts(fiber_exec, nfiber, ntotal, friction_contacts,
                                                   friction_solid_contacts, fiber_friction_count_buffers,
                                                   solid_friction_count_buffers, fiber_friction_bead_count,
                                                   solid_friction_count, fiber_friction_count_delta,
                                                   solid_friction_count_delta, fiber_reduce_scratch,
                                                   fric_list_mode ? &fric_ff_active : nullptr,
                                                   fric_list_mode ? &fric_fs_active : nullptr,
                                                   fric_list_mode ? nullptr : fiber_contact_search.seg_bound);
            // [2] Бесконечная масса полностью зажатых в mass splitting: c = 0.
            if (has_clamps)
                interactions::apply_clamp_friction_counts(fiber_exec, clamp_ws, fiber_friction_bead_count,
                                                          solid_friction_count);
            interactions::compute_friction_contributions(
                fiber_exec, fb, fupd, solver.fiber_rhs.dv, friction_contacts, sys, supd, solver.solid_rhs.dvx,
                friction_solid_contacts, fiber_friction_bead_count, solid_friction_count, /*solid_kinematic=*/false,
                fiber_friction_thread_buffers, solid_friction_thread_buffers,
                fric_list_mode ? &fric_ff_active : nullptr, fric_list_mode ? &fric_fs_active : nullptr,
                fric_list_mode ? nullptr : fiber_contact_search.seg_bound);
            rhs::reduce_contributions_into(fiber_exec, fiber_friction_thread_buffers, nfiber, friction_delta,
                                           fiber_reduce_scratch);
            for (size_t d = 0; d < dim; ++d)
                fiber_exec.parallel_for(nfiber, [&](std::size_t i) { solver.fiber_rhs.dv[d][i] += friction_delta[i].dv[d]; });
            if (has_active_solid_friction_contacts)
            {
                rhs::reduce_contributions_into(fiber_exec, solid_friction_thread_buffers, ntotal, friction_solid_delta,
                                               fiber_reduce_scratch);
                for (size_t d = 0; d < dim; ++d)
                    fiber_exec.parallel_for(
                        ntotal, [&](std::size_t i) { solver.solid_rhs.dvx[d][i] += friction_solid_delta[i].dv[d]; });
            }
            // [3] Вклад трения в зажатую частицу -- часть реакции: повторная
            // проекция (без XSPH). Вклад трения в свободные частицы от
            // зажатых не трогается.
            if (has_clamps)
            {
                interactions::apply_clamp_projection(fiber_exec, clamp_ws.solid, sys.mass, solver.solid_rhs.dvx,
                                                     static_cast<std::vector<std::vector<double>>*>(nullptr),
                                                     clamp_ws, clamp_solid_force, &clamp_ws.acc_solid,
                                                     /*acc_add=*/true);
                interactions::apply_clamp_projection(fiber_exec, clamp_ws.fiber, fb.mass, solver.fiber_rhs.dv,
                                                     static_cast<std::vector<std::vector<double>>*>(nullptr),
                                                     clamp_ws, clamp_fiber_force, &clamp_ws.acc_fiber,
                                                     /*acc_add=*/true);
            }
#ifdef SPH_ENERGY_LOG
            if (elog)
            {
                for (size_t d = 0; d < dim; ++d)
                    for (size_t i = 0; i < nfiber; ++i)
                        elog->friction_fiber[d][i] = solver.fiber_rhs.dv[d][i] - elog->friction_fiber[d][i];
                for (size_t d = 0; d < dim; ++d)
                    for (size_t i = 0; i < ntotal; ++i)
                        elog->friction_solid[d][i] = solver.solid_rhs.dvx[d][i] - elog->friction_solid[d][i];
                elog->has_friction = true;
            }
#endif
        }

        // Остаток условия Синьорини по итоговым скоростям шага (после трения и зажимов [3]) -- только
        // диагностика для монитора сцены: трение может дать нормальную составляющую в соседних контактах.
        if (proj_on && proj_stats != nullptr && proj_stats->enabled)
        {
            const auto fupd = solver.fiber_velocity_update(fb, dt);
            const auto supd = solver.solid_velocity_update(sys, dt);
            interactions::normal_projection_post_residual(fiber_exec, fb, fupd, solver.fiber_rhs.dv, normal_ff, sys,
                                                          supd, solver.solid_rhs.dvx, normal_fs, proj_ws, *proj_stats);
        }

        // Реакция опоры шага ([1] + [3], в этом порядке) -- в отчёт сцене.
        if (has_clamps && clamp_report != nullptr)
        {
            for (size_t d = 0; d < dim; ++d)
            {
                clamp_report->solid_force[d] = clamp_solid_force[d];
                clamp_report->fiber_force[d] = clamp_fiber_force[d];
            }
            // Сумма модулей и максимум полной реакции частиц ([1] + [3]); диагностика, на счёт не влияет.
            interactions::clamp_abs_summary(fiber_exec, clamp_ws.acc_solid, clamp_ws.solid.size(), clamp_ws,
                                            clamp_report->solid_force_abs, clamp_report->solid_force_max);
            interactions::clamp_abs_summary(fiber_exec, clamp_ws.acc_fiber, clamp_ws.fiber.size(), clamp_ws,
                                            clamp_report->fiber_force_abs, clamp_report->fiber_force_max);
        }

#ifdef SPH_ENERGY_LOG
        // Каналы шага: работа вклада на средней скорости шага (energy_log.hpp,
        // work_mid_*), потеря = -работа:
        //   трение -- ускорение трения шага (Fiber и Solid);
        //   контакт Solid-Solid -- его ускорение (после проекции зажимов);
        //   XSPH -- поправка av (у Эйлера av нет);
        //   тепло ИВ -- sum m*du_av*h (справочно, уже в U_int).
        if (elog)
        {
            const auto supd = solver.solid_velocity_update(sys, dt);
            if (elog->has_friction)
            {
                const auto fupd = solver.fiber_velocity_update(fb, dt);
                echan.friction -= energy::work_mid_fiber(fiber_exec, fb, fupd, solver.fiber_rhs.dv,
                                                         elog->friction_fiber, nfiber) +
                                  energy::work_mid_solid(fiber_exec, sys, supd, solver.solid_rhs.dvx, ntotal,
                                                         elog->friction_solid, false);
                elog->has_friction = false;
            }
            if (elog->has_contact_fiber)
            {
                const auto fupd = solver.fiber_velocity_update(fb, dt);
                echan.contact_fiber -= energy::work_mid_fiber(fiber_exec, fb, fupd, solver.fiber_rhs.dv,
                                                              elog->contact_fiber_acc, nfiber) +
                                       energy::work_mid_solid(fiber_exec, sys, supd, solver.solid_rhs.dvx, ntotal,
                                                              elog->contact_solid_acc, false);
                elog->has_contact_fiber = false;
            }
            if (wb.has_contact_dv && ntotal > 0)
                echan.contact_ss -=
                    energy::work_mid_solid(fiber_exec, sys, supd, solver.solid_rhs.dvx, ntotal, wb.contact_dv, false);
            if (supd.av != nullptr && ntotal > 0)
                echan.xsph -= energy::work_mid_solid(fiber_exec, sys, supd, solver.solid_rhs.dvx, ntotal, *supd.av, true);
            if (wb.du_av.size() >= ntotal)
                echan.q_av += energy::mass_weighted_sum(fiber_exec, sys.mass, wb.du_av, ntotal) * supd.h;
        }
#endif
        // OdeSolver: State::advance -- обновление на полный шаг (Solid: u/rho/s/
        // v/x с правилами, специфичными для Solid (u>=0, XSPH), — внутри
        // решателя, fiber/PLAN.md §3.9; Fiber: только x, v).
        solver.advance(fiber_exec, sys, ntotal, fb, nfiber, dt);
#ifdef SPH_CONTACT_STATS
        if (cstats) cstats->after_advance(time + dt, dt, fb, sys, friction_contacts, friction_solid_contacts);
#endif

        // Postprocessing, этап Ф7 (fiber/PLAN.md §3.8) -- ПОСЛЕ шага
        // ODESolver (predictor-corrector/Euler Solid И Fiber уже завершены
        // выше):
        //   1) Мизес/JH-2 -- interactions::apply_material_strength, буквально
        //      прежний код, читает финальные sys.sxx.. (место вызова не
        //      изменилось);
        //   2) разрыв связей Fiber-Bond -- interactions::
        //      run_fiber_topology_postprocessing: поиск событий (readonly
        //      parallel_for по topology.bonds, ТЕКУЩИЕ fb.x после advance()
        //      выше) -> применение строго последовательно (prev/next,
        //      новый object_id j-стороне, topology.bonds/version). Fiber-Bond
        //      kernel и Fiber-Fiber-Contact search строят primitives заново
        //      из topology/object_id на КАЖДОМ вызове RHS -- отдельной
        //      инвалидации кэша не требуется (fiber/SPH.md, «Postprocessing и
        //      изменение topology»).
#ifdef SPH_ENERGY_LOG
        // Лог энергии: пластическая работа Мизеса (Solid) и каналы нитей --
        // пластическая/вязкая работа связей и запас порвавшихся связей (FiberEnergyTap).
        interactions::apply_material_strength(cfg, sys, ntotal, dt, solver.solid_rhs.drho,
                                              elog ? &echan.w_plastic : nullptr);
        if (nfiber > 0)
        {
            interactions::run_fiber_topology_postprocessing(fiber_exec, fb, topology, materials, dt,
                                                            fiber_break_thread_buffers, elog ? &ftap : nullptr);
        }
#else
        interactions::apply_material_strength(cfg, sys, ntotal, dt, solver.solid_rhs.drho);
        if (nfiber > 0)
        {
            // Внутренние переменные связей (пластика/Максвелл/диссипация) --
            // с dt этого шага, по x^{n+1}, ПЕРЕД поиском разрыва: внутри обёртки
            // (apply_fiber_bond_constitutive). Для упругих материалов no-op.
            interactions::run_fiber_topology_postprocessing(fiber_exec, fb, topology, materials, dt,
                                                            fiber_break_thread_buffers);
        }
#endif
        if (time >= next_output_time)
        {
            output(sys, ntotal, time, current_ts, dt);
            if (nfiber > 0) fiber::fiber_output(fb, topology, materials, current_ts, time);
            next_output_time += save_t;
        }
        time += dt;

        // Вывод информации о мониторинге
        // Пишется ОСЕВАЯ компонента. В осесимметричной постановке ось удара —
        // x[1] (z), и именно осевую скорость поверхности измеряет интерферометр,
        // с которым сверяется эталон. Осевая компонента — последняя координата:
        // z в 3D, z (ось симметрии) в осесимметричном 2D, x в 1D.
        const size_t moni_axis = dim - 1;
        log_file << time << " " << sys.vx[moni_axis][moni_particle] << " " << sys.itype[moni_particle] << std::endl;

#ifdef SPH_ENERGY_LOG
        if (elog && elog->due(current_ts))
            elog->write_row(current_ts, time, dt,
                            energy::compute_energy_stocks(fiber_exec, cfg, sys, ntotal, fb, topology, materials,
                                                          elog->workspace()),
                            echan);
#endif

        // Сторож времени/застоя (run_guard.hpp). Без него цикл при коллапсе
        // шага крутится бесконечно: время задачи практически не растёт, а
        // признаков зависания для ОС нет. Причина обрыва остаётся в guard.
        // Чекпоинт -- ПОСЛЕ should_stop: наблюдатель шага (монитор) уже отработал, а
        // счётчики сторожа включают этот шаг. Периодический (current_ts % N == 0) и
        // финальный -- перед выходом по любой причине останова сторожа.
        const bool stop = guard != nullptr && guard->should_stop(dt, time, t_total);
        if (ckpt_on && (stop || current_ts % cfg.checkpoint_every_steps == 0)) write_checkpoint();
        if (stop)
        {
            // Штатный останов по числу шагов (StepLimit) аварией не является.
            if (guard->failed()) std::cerr << guard->message() << std::endl;
            break;
        }
    }
    // Нормальное завершение (time >= t_total): финальный чекпоинт, если последний шаг не записан.
    if (ckpt_on && current_ts != last_ckpt_ts) write_checkpoint();
}

} // namespace detail

// Основная функция временной интеграции (fiber/PLAN.md, этап Ф8: тонкая
// обёртка над detail::run_time_integration с пустыми Fiber-структурами по
// умолчанию -- FiberState::n == 0, FiberTopology::bonds.empty(), поэтому вся
// Fiber-физика внутри общего цикла no-op, как и раньше в перегрузке для
// fiber::Model). Выбор решателя -- по cfg.ode_solver (0 = PredictorCorrector,
// по умолчанию; 1 = Euler), инстанцирует run_time_integration нужным типом.
void time_integration(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, double& dt, double& time, double delta_x,
                      RunGuard* guard)
{
    fiber::FiberState empty_fb;
    fiber::FiberTopology empty_topology;
    std::vector<fiber::FiberMaterial> empty_materials;
    static const interactions::ClampSet kNoClamps; // Solid-only: зажимов нет
    if (cfg.ode_solver == 1)
        detail::run_time_integration<ode::EulerSolver>(cfg, sys, empty_fb, empty_topology, empty_materials, ntotal, dt,
                                                        time, delta_x, guard, kNoClamps, nullptr);
    else
        detail::run_time_integration<ode::PredictorCorrectorSolver>(cfg, sys, empty_fb, empty_topology,
                                                                     empty_materials, ntotal, dt, time, delta_x, guard,
                                                                     kNoClamps, nullptr);
}

// Перегрузка для fiber::Model (этап Ф2: fiber/PLAN.md пункт 7 — сигнатура;
// этап Ф3+: физика Fiber-Bond/Fiber-Fiber-Contact/Fiber-Solid-Contact/
// Postprocessing). Тонкая обёртка над detail::run_time_integration — та же
// функция, что и у старой сигнатуры выше, различаются только передаваемые
// model.fiber/model.topology/model.materials вместо пустых структур.
//
// Чекпоинты/рестарт (include/checkpoint.hpp, настройки cfg.checkpoint_* / restart_*):
//   - cfg.restart_file непуст: до цикла checkpoint::read перезаписывает Model, dt, time
//     и состояние цикла. cfg вызывающего АВТОРИТЕТЕН для запуска (t_total, save_t и
//     т.п. -- новые значения); физические отличия от сохранённого без cfg.restart_force
//     -- исключение std::runtime_error со списком полей. ntotal должен совпадать с сохранённым.
//   - cfg.checkpoint_file непуст и checkpoint_every_steps > 0: запись каждые N шагов
//     и финальная (останов сторожа / конец расчёта).
void time_integration(const SimulationConfig& cfg, fiber::Model& model, size_t ntotal, double& dt, double& time,
                      double delta_x, RunGuard* guard, const checkpoint::UserState* user_state)
{
    detail::CheckpointHooks hooks;
    checkpoint::Payload restored;
    const bool want_write = !cfg.checkpoint_file.empty() && cfg.checkpoint_every_steps > 0;
    if (!cfg.restart_file.empty())
    {
        SimulationConfig eff = cfg; // read() на выходе отдаёт действующие настройки (= новые)
        size_t nt = ntotal;
        checkpoint::RestoreOptions opt;
        opt.force = cfg.restart_force;
        // read() ничего не меняет в model при ошибке (проверки до фиксации), поэтому при
        // битом основном файле (CRC, обрезка) безопасно пробуем ротационную копию .prev.
        try
        {
            checkpoint::read(cfg.restart_file, eff, model, nt, restored, opt);
        }
        catch (const std::exception& e)
        {
            const std::string prev = cfg.restart_file + ".prev";
            if (!std::ifstream(prev).good()) throw;
            std::cerr << "ПРЕДУПРЕЖДЕНИЕ: не удалось прочитать " << cfg.restart_file << " (" << e.what()
                      << "); пробую " << prev << std::endl;
            eff = cfg;
            nt = ntotal;
            checkpoint::read(prev, eff, model, nt, restored, opt);
        }
        hooks.restart = &restored;
        if (user_state != nullptr && user_state->restore) user_state->restore(restored.user_state, restored.loop.guard_steps);
    }
    if (want_write)
    {
        hooks.on_checkpoint = [&cfg, &model, ntotal, user_state](checkpoint::Payload&& pl) {
            if (user_state != nullptr && user_state->save) pl.user_state = user_state->save();
            checkpoint::write(cfg.checkpoint_file, cfg, model, ntotal, pl);
        };
    }
    const detail::CheckpointHooks* h = (hooks.restart != nullptr || want_write) ? &hooks : nullptr;
    if (cfg.ode_solver == 1)
        detail::run_time_integration<ode::EulerSolver>(cfg, model.solid, model.fiber, model.topology, model.materials,
                                                        ntotal, dt, time, delta_x, guard, model.clamps,
                                                        &model.clamp_report, h, &model.proj_stats);
    else
        detail::run_time_integration<ode::PredictorCorrectorSolver>(cfg, model.solid, model.fiber, model.topology,
                                                                     model.materials, ntotal, dt, time, delta_x, guard,
                                                                     model.clamps, &model.clamp_report, h,
                                                                     &model.proj_stats);
}
