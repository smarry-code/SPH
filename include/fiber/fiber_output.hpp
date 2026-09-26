#pragma once
// fiber/fiber_output.hpp
// Вывод состояния нитей (fiber/PLAN.md, этап Ф2, пункт 4), аналог output()
// (output.hpp) для compute-структуры Fiber. Определение — в
// src/fiber_output.cpp. Формат файла — в комментарии над fiber_output().
//
// Вызывается из time_integration рядом с output() (та же частота save_t),
// только когда нити есть.
#include <cstddef>
#include "fiber/model.hpp"

namespace fiber
{

// Пишет fiber_output_<itimestep>.txt (5 цифр, дополнено нулями слева, как
// sph_output_<itimestep>.txt у output()): по каждой частице object_id,
// координаты, скорости; по каждой связи — l0, текущую длину l, l_break и
// относительное удлинение (l-l0)/l0.
void fiber_output(const FiberState& fiber, const FiberTopology& topology, std::size_t itimestep, double time);
void fiber_output(const Model& model, std::size_t itimestep, double time);

} // namespace fiber
