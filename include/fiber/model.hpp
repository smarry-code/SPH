#pragma once
// fiber/model.hpp
// Model — объединение всех compute-структур одного расчёта (fiber/PLAN.md
// §3.1: «Model объединяет SolidState, FiberState, FiberTopology и таблицу
// материалов»).
//
// SolidState — существующий ParticleSystem (particle_system.hpp), без
// переименования: Ф2 не меняет архитектуру Solid, только добавляет данные
// нитей рядом с ним. Model владеет ParticleSystem по значению — так и тесты
// модели данных нитей, и будущий вызывающий код (time_integration) работают
// с одним объектом, не задумываясь о времени жизни отдельных частей.
#include "fiber/fiber_system.hpp"
#include "particle_system.hpp"

namespace fiber
{

struct Model
{
    ParticleSystem solid;
    FiberState fiber;
    FiberTopology topology;
    std::vector<FiberMaterial> materials;
};

} // namespace fiber
