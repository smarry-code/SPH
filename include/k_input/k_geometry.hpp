// k_geometry.hpp
// Сборка начальной геометрии из LS-DYNA .k (этап K4, readk/k-import-plan.md,
// Р3-Р6, Р9): KDeck + таблица деталей -> ParticleSystem (режим solid).
// Любая ошибка - std::runtime_error с файлом и строкой.
#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "fiber/model.hpp"
#include "k_input/k_parts.hpp"
#include "k_input/k_reader.hpp"
#include "materials.hpp"
#include "particle_system.hpp"

namespace k_input
{

enum class KMode { Solid, SolidFiber };

struct KImportOptions
{
    std::string k_file;
    KMode mode = KMode::Solid;
    // Пустой parts_file - "режим без таблицы" (наследие parse_k_file): заполняются
    // только x, mass, itype = body = PID; rho, hsml и прочее - нули, urs не
    // вызывается. В режиме SolidFiber таблица обязательна.
    std::string parts_file;
    double length_scale = 1.0;    // множитель координат (Р5)
    double mass_scale = 1.0;      // множитель масс (Р5)
    long long monitor_nid = -1;   // NID мониторной частицы (Р4); < 0 - не задан
};

struct KImportSummary
{
    std::size_t nodes = 0;               // узлов в колоде
    std::size_t sph_particles = 0;
    std::array<double, 3> bbox_min{0.0, 0.0, 0.0};   // габарит SPH-частиц, м
    std::array<double, 3> bbox_max{0.0, 0.0, 0.0};
    double dx_median = 0.0;              // медиана (m/rho)^(1/dim); 0 без таблицы
    std::vector<std::string> warnings;
    std::map<std::string, std::size_t> skipped;   // пропущенные ключевые слова
    // Поля для нитей добавит K5.
};

// Solid-часть: по записи deck.sph на частицу, в порядке файла. table == nullptr -
// режим без таблицы. moni_particle меняется только при opts.monitor_nid >= 0.
void build_solid_from_deck(const KDeck& deck, const KPartsTable* table, const KImportOptions& opts,
                           const Materials& mats, ParticleSystem& sys, std::size_t& moni_particle,
                           KImportSummary& summary);

// Точка входа: читает колоду и таблицу, сверяет их, пишет в model.solid.
KImportSummary load_k_geometry(const KImportOptions& opts, const Materials& mats, fiber::Model& model,
                               std::size_t& moni_particle);

} // namespace k_input
