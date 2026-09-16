#pragma once
// int_force.hpp
// Объявления функций внутренних сил SPH (urs — диспетчер УРС, int_force —
// основной силовой расчёт). Определения — в int_force.cpp.
//
// Горячие циклы (body1/body2) переведены с AoS-доступа на плоские
// __restrict__-указатели в ParticleSystem (SoA): потоковые чтения/записи в
// непрерывные массивы -> кэш-дружелюбность, доказательство неалиасинга для
// автовекторизации.
#include <cstddef>
#include <vector>
#include <array>
#include <cassert>
#include "buffers.hpp"
#include "pair_graph.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "eos.hpp"

// Диспетчеризация уравнения состояния по типу материала (itype) через switch
// с default — управляемый отказ.
void urs(const Materials& mats, ParticleSystem& sys, size_t i);

// Основная функция расчёта внутренних сил SPH.
void int_force(const SimulationConfig& cfg, ParticleSystem& sys, size_t ntotal, size_t niac,
               const std::vector<double>& eta, const std::vector<size_t>& pair_i,
               const std::vector<size_t>& pair_j, const std::vector<std::vector<double>>& dwdx,
               const std::vector<double>& w, std::vector<double>& t,
               std::vector<std::vector<double>>& dvxdt, std::vector<double>& tdsdt, std::vector<double>& dedt, std::vector<double>& dsxx, std::vector<double>& dsyy,
               std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, const std::vector<double>& drho,
               WorkBuffers& wb, const PairGraph& graph);
