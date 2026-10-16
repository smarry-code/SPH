#pragma once
// interactions/solid_workspace.hpp
// Постоянная рабочая память шага Solid (single_step, src/single_step.cpp).
//
// Раньше single_step создавал все эти структуры локально на каждом вызове, и
// каждый шаг выделял/освобождал сотни мегабайт (glibc отдаёт крупные блоки
// через mmap/munmap -> порядка 4e5 page fault'ов на шаг; speedup/base.md
// §4.2). Теперь они живут между шагами: все стадии (compute_solid_derived,
// build_ghost_set, build_solid_extended, build_solid_sph_search,
// compute_solid_*_contributions, reduce_contributions_into, finalize_solid)
// заполняют их через assign/resize/clear, поэтому ёмкость только растёт и
// на шагах после первого новых аллокаций почти нет. Размеры каждый раз
// задаются заново (ntotal и число ghost могут меняться от шага к шагу) --
// содержимое с предыдущего шага не используется.
//
// Владелец -- WorkBuffers::solid_ws (buffers.hpp), создаётся лениво в
// single_step при первом вызове.
#include <vector>

#include "interactions/solid_contact.hpp"
#include "interactions/solid_sph.hpp"
#include "rhs/contributions.hpp"

namespace interactions
{

struct SolidStepWorkspace
{
    SolidDerived derived;
    GhostSet ghosts;
    SolidExtended ext;
    SolidSphSearchResult search; // вместе с SpatialIndex и буфером пар поиска

    // Вклады по interaction type; в reduce идут в этом порядке (сначала
    // Solid-SPH, затем Solid-Solid-Contact). Solid-SPH пишет по заранее
    // известным позициям в постоянный массив; контакт пишет вклад только при
    // fn > 0 (число заранее неизвестно) -- в буферы потоков.
    rhs::IndexedContributions<SolidDelta> sph_contribs;
    std::vector<rhs::ContributionBuffer<SolidDelta>> contact_buffers;
    std::vector<SolidDelta> reduced;
    rhs::ReduceScratch reduce_scratch;
    std::vector<SolidFinal> finalres;
#ifdef SPH_ENERGY_LOG
    std::vector<SolidDelta> contact_reduced; // только вклады контакта (лог энергии)
#endif
};

} // namespace interactions
