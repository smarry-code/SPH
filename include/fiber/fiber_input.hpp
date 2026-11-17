#pragma once
// fiber/fiber_input.hpp
// Разбор файла ломаных нитей (fiber/PLAN.md, этап Ф2, пункт 3). Формат файла
// подробно описан в комментарии над load_fiber_file (src/fiber_input.cpp).
// Объявление здесь, определение — в src/fiber_input.cpp (общий конвейер,
// как остальные модули проекта: заголовки в include/, исходники в src/).
#include <string>
#include "fiber/model.hpp"

namespace fiber
{

// Разбирает файл нитей filename в model.fiber/model.topology/model.materials
// (существующее содержимое model.solid не трогается). Бросает
// std::runtime_error:
//   - если сборка сделана не с SPH_DIM=3 (check_fiber_requires_dim3);
//   - если файл не открывается или содержит ошибку формата;
//   - если FIBER ссылается на неопределённый MATERIAL;
//   - если нить содержит меньше двух точек или нулевой/отрицательный сегмент;
//   - если BOND (связь-стежок) ссылается на бусину вне диапазона, на неопределённый
//     материал, соединяет бусину с собой или совпадающие бусины.
void load_fiber_file(const std::string& filename, Model& model);

} // namespace fiber
