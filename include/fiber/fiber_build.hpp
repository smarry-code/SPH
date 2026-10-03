#pragma once
// fiber/fiber_build.hpp
// Общая сборка модели нитей из готовых ломаных и разбор директивы MATERIAL.
// Используется разбором файла нитей (load_fiber_file, src/fiber_input.cpp) и
// импортом LS-DYNA .k (сборщик нитей из *ELEMENT_BEAM). Определения — в
// src/fiber_input.cpp.
#include <array>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "fiber/model.hpp"
#include "interactions/clamp.hpp"
#include "param.hpp"

namespace fiber
{

// Ломаная нити в недеформированной конфигурации: индекс материала и
// последовательность точек (x, y, z).
struct PolylineSpec
{
    int material_id = -1;
    std::vector<std::array<double, 3>> points;
};

// Бокс закрепления: бусины нитей внутри [lo, hi] получают маску mask
// (биты компонент, interactions::kClampAll == закреплена целиком).
// Не называть ClampBox: fiber::ClampBox занято в fiber/plain_weave.hpp.
struct ClampBoxSpec
{
    double lo[dim];
    double hi[dim];
    std::uint8_t mask = interactions::kClampAll;
};

// Собирает model.materials/fiber/topology (и clamps.fiber_mask) из ломаных.
// model.solid не трогается. Для каждой ломаной:
//   - бусина на каждую точку; скорость нулевая; object_id = номер ломаной;
//     material_id = polyline.material_id;
//   - масса бусины rho * A * (l_prev + l_next) / 2, A = pi*d^2/4, для
//     концевой бусины берётся только один прилежащий сегмент (второе
//     слагаемое отсутствует, а не делится пополам дважды);
//   - связь на каждый сегмент: l0 = длина сегмента,
//     l_break = (1 + eps_break) * l0; prev/next по цепочке (-1 на концах).
// Маска CLAMP: если boxes непуст, fiber_mask заполняется нулями и
// размечается interactions::mark_clamp_box по каждому боксу; иначе
// fiber_mask очищается (пустой, не заполненный нулями).
// Увеличивает model.topology.version.
// Бросает std::runtime_error: материал ломаной вне materials (индекс вне
// диапазона), меньше двух точек, нулевой/отрицательный сегмент (в тексте
// последней ошибки префикс who: для .fib - "fiber::load_fiber_file", для .k - свой).
void build_model_from_polylines(const std::vector<PolylineSpec>& polylines,
                                const std::vector<FiberMaterial>& materials,
                                const std::vector<ClampBoxSpec>& boxes, Model& model,
                                const std::string& who = "fiber::load_fiber_file");

// Разбор директивы "MATERIAL <id> <E> <nu> <d> <rho> <eps_break> [mu_ff [mu_fs]]".
// iss позиционирован сразу после слова MATERIAL; line_no — номер строки для
// сообщений; who — префикс сообщений об ошибках (например
// "fiber::load_fiber_file"). Записывает материал в materials[id]
// (с расширением массива) и отмечает material_defined[id] = 1; повторное
// определение id перезаписывает прежнее. Бросает std::runtime_error при
// неверном формате и отрицательных mu_ff/mu_fs.
void parse_material_directive(std::istringstream& iss, std::size_t line_no, const std::string& who,
                              std::vector<FiberMaterial>& materials, std::vector<char>& material_defined);

} // namespace fiber
