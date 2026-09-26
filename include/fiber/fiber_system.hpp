#pragma once
// fiber/fiber_system.hpp
// Модель данных нити (fiber/PLAN.md, этап Ф2, раздел «3.1. Данные»):
// FiberState (переменные состояния бусин), FiberTopology (связи и соседство
// по нити), FiberMaterial (свойства материала нити).
//
// Все поля — real_t, хранение — Structure of Arrays, частицы всех нитей
// одного типа лежат в общих массивах (fiber/SPH.md, «Организация хранения
// частиц»), по образцу ParticleSystem (particle_system.hpp).
//
// Компилируется при ЛЮБОМ SPH_DIM (заголовок — часть общего конвейера
// src/*.cpp, см. CMakeLists.txt: SPH_PIPELINE_SOURCES собирается один раз на
// каждую цель со своим SPH_DIM). Физически нить требует dim == 3
// (fiber/SPH.md: «Нити поддерживаются только в трёхмерной постановке»,
// fiber/PLAN.md Р9) — эту границу проверяет check_fiber_requires_dim3()
// НИЖЕ, а не сам факт компиляции структур: FiberState::x[dim]/vx[dim]
// зеркалирует ParticleSystem::x[dim]/vx[dim] и безопасен при любом dim,
// пока код обращается только к индексам [0, dim).
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#include "param.hpp"

namespace fiber
{

// Проверка границы Р9. Бросает исключение, если код, реально использующий
// модель нитей (разбор файла нитей, построение FiberState/FiberTopology,
// вывод), вызван в сборке с SPH_DIM, отличным от 3. Вызывается ДО любого
// обращения к данным нити — как из parse_sph_config (как только в конфиге
// обнаружен путь к файлу нитей), так и из load_fiber_file/fiber_output
// повторно (defence in depth), чтобы проверка срабатывала «на этапе разбора
// конфигурации/подготовки расчёта, а не потом в середине шага»
// (fiber/PLAN.md, этап Ф2, пункт 5).
//
// dim — compile-time constexpr (param.hpp, ADR-002), поэтому ветка
// dim != 3 недостижима в бинарнике, собранном с SPH_DIM=3: её проверяет
// отдельная лёгкая цель test_fiber_dim_guard, собранная БЕЗ SPH_DIM
// (dim = 1 по умолчанию) — см. CMakeLists.txt и
// tests/test_fiber_dim_guard.cpp.
inline void check_fiber_requires_dim3()
{
    if constexpr (dim != 3)
    {
        throw std::runtime_error(
            "Нити (fiber_file) заданы в конфигурации, но сборка выполнена с "
            "SPH_DIM=" + std::to_string(dim) +
            " -- нити поддерживаются только при SPH_DIM=3 (fiber/SPH.md, "
            "раздел \"Fiber\"; fiber/PLAN.md, решение Р9).");
    }
}

// --- FiberState -------------------------------------------------------------
// Переменные состояния бусин нити: координаты, скорости, постоянная масса
// (из ненагруженной конфигурации, fiber/SPH.md), номер нити (меняется при
// разрыве, этап Ф7) и индекс материала (не меняется при разрыве).
struct FiberState
{
    std::vector<real_t> x[dim];
    std::vector<real_t> vx[dim];
    std::vector<real_t> mass;
    std::vector<int> object_id;
    std::vector<int> material_id;

    std::size_t n = 0;

    void resize(std::size_t new_n)
    {
        n = new_n;
        for (std::size_t d = 0; d < dim; ++d)
        {
            x[d].resize(n);
            vx[d].resize(n);
        }
        mass.resize(n);
        object_id.resize(n);
        material_id.resize(n);
    }

    std::size_t size() const { return n; }
};

// --- FiberTopology -----------------------------------------------------------
// Связь топологических соседей (fiber/SPH.md, Fiber-Bond, «Постановка
// проекта: только растяжение»): пара бусин i<j по цепочке нити, длина в
// ненагруженной конфигурации l0 и порог разрыва l_break = (1+eps_break)*l0.
struct FiberBond
{
    std::size_t i = 0;
    std::size_t j = 0;
    real_t l0 = 0.0;
    real_t l_break = 0.0;
};

// prev[]/next[] — соседи по нити, -1 на конце (частица — конец нити, если
// prev == -1 или next == -1). bonds — каждая связь один раз. version —
// topology_version (fiber/PLAN.md §3.1): растёт при изменении топологии
// (разрыв связи, этап Ф7); используется, чтобы инвалидировать кэшированные
// от старой топологии структуры (Postprocessing, fiber/SPH.md).
struct FiberTopology
{
    std::vector<long long> prev;
    std::vector<long long> next;
    std::vector<FiberBond> bonds;
    unsigned long long version = 0;

    void resize(std::size_t new_n)
    {
        prev.assign(new_n, -1);
        next.assign(new_n, -1);
    }

    std::size_t size() const { return prev.size(); }
};

// --- FiberDelta ---------------------------------------------------------------
// Изменение состояния Fiber-частицы, суммируемое на стадии Reduce
// (fiber/PLAN.md §3.3: «FiberDelta: dv[dim]»; fiber/SPH.md, «Основные
// требования к interaction kernel», п.6: тип Delta определяется
// compute-структурой target частицы). Единственное поле — ускорение;
// interaction kernels (Fiber-Bond, позже Fiber-Fiber/Fiber-Solid) пишут в
// него через rhs::ContributionWriter<FiberDelta>, суммирование —
// rhs::reduce_contributions (operator+= ниже).
struct FiberDelta
{
    real_t dv[dim] = {};

    FiberDelta& operator+=(const FiberDelta& other)
    {
        for (std::size_t d = 0; d < dim; ++d) dv[d] += other.dv[d];
        return *this;
    }
};

// --- FiberMaterial -----------------------------------------------------------
// Свойства материала нити (fiber/PLAN.md §3.1; fiber/SPH.md, Fiber-Bond и
// Fiber-Solid-Contact): модуль Юнга, коэффициент Пуассона, диаметр,
// плотность, критическое удлинение разрыва.
struct FiberMaterial
{
    real_t E = 0.0;
    real_t nu = 0.0;
    real_t d = 0.0;
    real_t rho = 0.0;
    real_t eps_break = 0.0;
};

} // namespace fiber
