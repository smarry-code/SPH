#pragma once
// single_step.hpp
// Объявления функции печати текущего времени (time_print) и оркестратора
// одного шага симуляции (single_step). Определения — в single_step.cpp.
#include <cstddef>
#include <vector>
#include "pair_graph.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "diagnostics.hpp"
#include "buffers.hpp"

// Печать текущего системного времени (используется в main.cpp).
void time_print();

// Оркестратор одного шага (fiber/PLAN.md, этап Ф4-B1: точка переключения
// продового расчёта на новый путь RHS -- include/interactions/solid_sph.hpp
// + include/interactions/solid_contact.hpp). Обе перегрузки time_integration
// (ParticleSystem&, fiber::Model&) вызывают ЭТУ функцию, поэтому изменение
// её тела переключает оба цикла одновременно (src/time_integration.cpp).
// Временный шим совместимости (до этапа Ф7, fiber/PLAN.md §3.6): p, c,
// exx..eyz, rxy..ryz по-прежнему пишутся в sys, чтобы prochnost_mises/
// prochnost_cem_JH2 (src/ppochnost.cpp, вызываются в time_integration.cpp
// ПОСЛЕ single_step) и compute_cfl_timestep (читает sys.c) продолжали
// работать без изменений.
void single_step(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, double dt, size_t ntotal, std::vector<double>& t, std::vector<double>& tdsdt,
                 std::vector<std::vector<double>>& dx, std::vector<std::vector<double>>& dvx, std::vector<double>& du, std::vector<double>& ds,
                 std::vector<double>& drho, std::vector<std::vector<double>>& av, std::vector<double>& dsxx, std::vector<double>& dsyy,
                 std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, WorkBuffers& wb,
                 PairGraph& graph);

// Старый путь (direct_find/link_list -> con_density -> int_force -> art_visc
// -> ext_force -> h_upgrade -> av_vel), СОХРАНЁН НЕТРОНУТЫМ и компилируемым
// для части Б2 (сравнение/откат), но НЕ вызывается ни из одной перегрузки
// time_integration (fiber/PLAN.md, этап Ф4-B1). Сигнатура и тело --
// побитная копия прежнего single_step().
// TODO(Ф4-B2): удалить вместе со старым путём (direct_find, link_list,
// int_force, art_visc, av_vel, ext_force, virt_part) после независимого
// ревью части B1.
void single_step_legacy(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, double dt, size_t ntotal, std::vector<double>& t, std::vector<double>& tdsdt,
                 std::vector<std::vector<double>>& dx, std::vector<std::vector<double>>& dvx, std::vector<double>& du, std::vector<double>& ds,
                 std::vector<double>& drho, std::vector<std::vector<double>>& av, std::vector<double>& dsxx, std::vector<double>& dsyy,
                 std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, WorkBuffers& wb,
                 PairGraph& graph);
