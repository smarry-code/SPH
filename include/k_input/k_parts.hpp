// k_parts.hpp
// Таблица деталей (sidecar-файл) для импорта LS-DYNA .k (этап K3,
// readk/k-import-plan.md, разд. 3 Р2): привязка PID из .k к физике проекта.
//
// Формат (текст, поля через пробелы, '#' - комментарий до конца строки, одна
// директива на строку, директивы регистрозависимы как в .fib):
//   SOLID_PART <pid> <itype> <rho> [hsml=<h>] [v=<vx> [<vy> [<vz>]]]
//   MATERIAL <id> <E> <nu> <d> <rho> <eps_break> [mu_ff [mu_fs]]   (как в .fib)
//   FIBER_PART <pid> <material_id>
// Именованные параметры SOLID_PART - в любом порядке, каждый не более раза.
// Любая ошибка - std::runtime_error "<файл>:<строка>: сообщение".
#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "k_input/k_reader.hpp"

namespace k_input
{

struct SolidPart
{
    long long pid = 0;
    int itype = 0;                    // 1..15 (12..15 - Джонсон-Кук)
    double rho = 0.0;                 // > 0
    bool has_hsml = false;
    double hsml = 0.0;                // > 0, если has_hsml
    bool has_v = false;
    int nv = 0;                       // сколько компонент задано в v= (1..3)
    std::array<double, 3> v{0.0, 0.0, 0.0};
    std::size_t line = 0;             // строка в файле таблицы (1-based)
};

struct FiberPart
{
    long long pid = 0;
    int material_id = 0;
    std::size_t line = 0;
};

struct KPartsTable
{
    std::string file;                         // имя файла (или virtual_name) для сообщений
    std::vector<SolidPart> solid;
    std::vector<FiberPart> fiber;
    std::vector<fiber::FiberMaterial> materials;   // по индексу id
    std::vector<char> material_defined;            // 1, если id объявлен

    // PID -> индекс в solid / fiber.
    std::unordered_map<long long, std::size_t> solid_index;
    std::unordered_map<long long, std::size_t> fiber_index;

    const SolidPart* find_solid(long long pid) const;
    const FiberPart* find_fiber(long long pid) const;
};

KPartsTable read_k_parts_file(const std::string& path);
KPartsTable read_k_parts_from_string(const std::string& text, const std::string& virtual_name);

// Сверка таблицы с колодой. Ошибки (runtime_error с именем файла таблицы):
// PID из deck.sph без SOLID_PART; при fiber_mode - PID из deck.beams без
// FIBER_PART. Предупреждения (в warnings): PID таблицы, которого нет ни в
// sph, ни в beams, ни в solids, ни в *PART колоды. PID из *ELEMENT_SOLID
// проверяются при solid_elements = true (K9): PID без SOLID_PART - ошибка.
void check_parts_against_deck(const KPartsTable& table, const KDeck& deck, bool fiber_mode,
                              std::vector<std::string>& warnings, bool solid_elements = false);

} // namespace k_input
