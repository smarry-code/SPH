#ifndef RHS_TYPES_HPP
#define RHS_TYPES_HPP

#include <cstddef>

// ============================================================================
// rhs/types.hpp — общие типы CPU backend вычисления RHS (notes/fiber-model/fiber-problem-statement.md, раздел
// "Организация хранения частиц" и notes/fiber-model/fiber-plan.md §3). Этап Ф1: только
// инфраструктура, никакой физики (Fiber/Solid) здесь и не подразумевается.
// ============================================================================

namespace rhs
{

// Стабильный в рамках одного вычисления RHS индекс частицы внутри своей
// compute-структуры (fiber-problem-statement.md: "index_t particle_id").
using index_t = std::size_t;

} // namespace rhs

#endif // RHS_TYPES_HPP
