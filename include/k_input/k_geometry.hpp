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

#include "fiber/fiber_build.hpp"
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
    // Нити (режим solid_fiber).
    std::size_t fibers = 0, beads = 0, bonds = 0;
    double seg_min = 0.0, seg_max = 0.0;   // длина сегмента l0, м
    double ld_min = 0.0, ld_max = 0.0;     // l0 / d
    // Закрепления из *BOUNDARY_SPC_* (K8).
    std::size_t clamped_particles = 0, clamped_beads = 0;
};

// Solid-часть: по записи deck.sph на частицу, в порядке файла. table == nullptr -
// режим без таблицы. moni_particle меняется только при opts.monitor_nid >= 0.
void build_solid_from_deck(const KDeck& deck, const KPartsTable* table, const KImportOptions& opts,
                           const Materials& mats, ParticleSystem& sys, std::size_t& moni_particle,
                           KImportSummary& summary);

// Нити из *ELEMENT_BEAM (Р1, Р9): нить - простая незамкнутая цепочка балок с
// одним PID; порядок нитей - по первой балке в файле, обход от конца с меньшим
// NID. Требует SPH_DIM=3 и FIBER_PART для каждого PID балок. Пишет в
// model.materials/fiber/topology (model.solid не трогает).
void build_fibers_from_deck(const KDeck& deck, const KPartsTable& table, const KImportOptions& opts,
                            fiber::Model& model, KImportSummary& summary,
                            std::vector<long long>* bead_nids = nullptr);

// Закрепления *BOUNDARY_SPC_NODE/_SET (K8): DOFX/DOFY/DOFZ != 0 -> бит компоненты
// (только c < dim) в model.clamps.solid_mask (частица i <-> deck.sph[i]) и
// fiber_mask (бусина k <-> bead_nids[k]). CID != 0 - ошибка; повторы - по OR;
// без единого бита маска остаётся пустой; непустая fiber_mask объединяется по OR.
// Компонента >= dim и узлы без роли - предупреждения в summary.
void apply_spc_clamps(const KDeck& deck, const std::vector<long long>& bead_nids, fiber::Model& model,
                      KImportSummary& summary);

// Точка входа: читает колоду и таблицу, сверяет их, пишет в model.solid.
KImportSummary load_k_geometry(const KImportOptions& opts, const Materials& mats, fiber::Model& model,
                               std::size_t& moni_particle);

} // namespace k_input
