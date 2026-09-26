#pragma once
// single_step.hpp
// Объявления функции печати текущего времени (time_print) и оркестратора
// одного шага симуляции (single_step). Определения — в single_step.cpp.
//
// fiber/PLAN.md, этап Ф4-B2: старый путь (direct_find/link_list ->
// con_density -> int_force -> art_visc -> ext_force -> h_upgrade -> av_vel,
// прежде — single_step_legacy) удалён вместе с graph-параметром: новый путь
// (interactions/solid_sph.hpp + interactions/solid_contact.hpp) не строит
// раскраску пар (rhs::reduce_contributions суммирует по шардам
// owner-computes, а не по цветам графа), поэтому PairGraph больше не нужен
// ни одному потребителю single_step.
#include <cstddef>
#include <vector>
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "diagnostics.hpp"
#include "buffers.hpp"

// Печать текущего системного времени (используется в main.cpp).
void time_print();

// Оркестратор одного шага (fiber/PLAN.md, этап Ф4: новый путь RHS для Solid,
// include/interactions/solid_sph.hpp + include/interactions/solid_contact.hpp).
// Обе перегрузки time_integration (ParticleSystem&, fiber::Model&) вызывают
// ЭТУ функцию, поэтому изменение её тела переключает оба цикла одновременно
// (src/time_integration.cpp). Временный шим совместимости (до этапа Ф7,
// fiber/PLAN.md §3.6): p, c, exx..eyz, rxy..ryz по-прежнему пишутся в sys,
// чтобы prochnost_mises/prochnost_cem_JH2 (src/ppochnost.cpp, вызываются в
// time_integration.cpp ПОСЛЕ single_step) и compute_cfl_timestep (читает
// sys.c) продолжали работать без изменений.
void single_step(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, real_t dt, size_t ntotal, std::vector<real_t>& t, std::vector<real_t>& tdsdt,
                 std::vector<std::vector<real_t>>& dx, std::vector<std::vector<real_t>>& dvx, std::vector<real_t>& du, std::vector<real_t>& ds,
                 std::vector<real_t>& drho, std::vector<std::vector<real_t>>& av, std::vector<real_t>& dsxx, std::vector<real_t>& dsyy,
                 std::vector<real_t>& dszz, std::vector<real_t>& dsxy, std::vector<real_t>& dsyz, std::vector<real_t>& dsxz, WorkBuffers& wb);
