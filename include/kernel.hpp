#pragma once
// kernel.hpp
// Объявление функции вычисления ядра сглаживания и его производной.
// Определение — в kernel.cpp.
//
// Параметр dx принимает std::span<const real_t> — это допускает передачу
// как std::vector<real_t>, так и std::array<real_t, dim> (и плоских буферов SoA в будущем),
// устраняя принудительное использование std::vector у вызывающих сторон.
// skf (тип ядра: 1 - кубический сплайн, 2 - Гаусс, 3 - квинтическое) передаётся
// явно (Interface Segregation): ядро — низкоуровневый примитив, ему нужен только
// тип ядра, а не весь SimulationConfig. Ранее читался из глобала skf (param.hpp).
#include <span>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include "param.hpp"

// Вычисление ядра W и его производной dW/dx.
// skf — тип ядра; r — расстояние; dx — компоненты вектора x_ij;
// mhsml — средняя длина сглаживания; w — значение ядра; dwdx[d] — градиент ядра.
void kernel(int skf, real_t r, std::span<const real_t> dx, real_t mhsml, real_t& w, real_t* dwdx);
