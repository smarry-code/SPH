// test_k_clamp.cpp
// Unit-тесты закреплений из *BOUNDARY_SPC_NODE/_SET (readk/k-import-plan.md, этап K8):
//   - SPC_NODE на SPH-частицах: solid_mask по индексу частицы (порядок *ELEMENT_SPH);
//   - SPC_SET через *SET_NODE_LIST и *SET_NODE_LIST_GENERATE;
//   - SPC на узлах нити (solid_fiber): fiber_mask по порядку бусин;
//   - OR при повторных SPC; пустая маска без SPC и при нулевых DOF;
//   - узел без роли - предупреждение; CID != 0 - ошибка; fixed-формат (ширина 10).
// SPH_DIM=3.
#include "fiber/model.hpp"
#include "k_input/k_geometry.hpp"
#include "simulation_config.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

using namespace k_input;

template <class F>
std::string error_of(F&& f)
{
    try { f(); }
    catch (const std::exception& e) { return e.what(); }
    return std::string();
}

bool has(const std::string& s, const std::string& sub) { return s.find(sub) != std::string::npos; }

const char* kTable =
    "SOLID_PART 1 5 1000\n"
    "MATERIAL 0 2.0e9 0.35 1.0e-4 1140 0.02\n"
    "FIBER_PART 2 0\n";

// SPH-частицы в порядке файла: NID 3, 1, 4, 2 -> индексы 0..3. Узел 9 - без роли.
// Нить PID 2: цепочка NID 10-5-7-8, обход от 8 -> бусины 0..3 = NID 8, 7, 5, 10.
const char* kGeom =
    "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n4,3,0,0\n9,50,0,0\n"
    "10,0,5,0\n5,1,5,0\n7,2.5,5,0\n8,4,5,0\n"
    "*PART\nA\n1\n*PART\nF\n2\n"
    "*ELEMENT_SPH\n3,1,1.0\n1,1,1.0\n4,1,1.0\n2,1,1.0\n"
    "*ELEMENT_BEAM\n1,2,10,5\n2,2,5,7\n3,2,7,8\n";

struct Built
{
    fiber::Model model;
    KImportSummary sum;
};

Built build(const std::string& tail, bool fibers = false, const std::string& geom = kGeom)
{
    Built b;
    const KDeck deck = read_k_deck_from_string(geom + tail, "mem.k");
    const KPartsTable t = read_k_parts_from_string(kTable, "t.parts");
    KImportOptions o;
    o.k_file = "mem.k";
    o.mode = fibers ? KMode::SolidFiber : KMode::Solid;
    std::size_t moni = 0;
    build_solid_from_deck(deck, nullptr, o, Materials{}, b.model.solid, moni, b.sum);
    std::vector<long long> beads;
    if (fibers) build_fibers_from_deck(deck, t, o, b.model, b.sum, &beads);
    apply_spc_clamps(deck, beads, b.model, b.sum);
    return b;
}

std::string geom_wo_beams()
{
    std::string g = kGeom;
    g.erase(g.find("*ELEMENT_BEAM"));
    return g;
}

bool has_warning(const Built& b, const std::string& sub)
{
    for (const auto& w : b.sum.warnings)
        if (has(w, sub)) return true;
    return false;
}

void test_node()
{
    std::printf("[spc_node]\n");
    Built b = build("*BOUNDARY_SPC_NODE\n1,0,1,1,1,0,0,0\n3,0,1,0,0\n4,0,0,1,0\n", false, geom_wo_beams());
    const auto& m = b.model.clamps.solid_mask;
    // NID 3 -> частица 0 (x), NID 1 -> 1 (xyz), NID 4 -> 2 (y), NID 2 -> 3 (свободна).
    check(m.size() == 4, "solid_mask полного размера");
    check(m.size() == 4 && m[0] == 1 && m[1] == 7 && m[2] == 2 && m[3] == 0, "биты по индексам частиц");
    check(b.model.clamps.fiber_mask.empty(), "fiber_mask пуста");
    check(b.sum.clamped_particles == 3 && b.sum.clamped_beads == 0, "сводка: 3 частицы");
}

void test_set()
{
    std::printf("[spc_set]\n");
    Built b = build("*SET_NODE_LIST\n5\n1,2\n*SET_NODE_LIST_GENERATE\n6\n3,4\n"
                    "*BOUNDARY_SPC_SET\n5,0,0,1,0\n6,0,1,0,0\n",
                    false, geom_wo_beams());
    const auto& m = b.model.clamps.solid_mask;
    // NID 1,2 -> частицы 1,3 (y); NID 3,4 -> частицы 0,2 (x).
    check(m.size() == 4 && m[0] == 1 && m[1] == 2 && m[2] == 1 && m[3] == 2, "LIST и GENERATE");
}

void test_fiber()
{
    std::printf("[fiber]\n");
    Built b = build("*BOUNDARY_SPC_NODE\n8,0,1,1,1\n10,0,0,0,1\n1,0,1,0,0\n", true);
    const auto& f = b.model.clamps.fiber_mask;
    check(f.size() == 4 && f[0] == 7 && f[1] == 0 && f[2] == 0 && f[3] == 4, "fiber_mask: NID 8 -> бусина 0, NID 10 -> 3");
    check(b.model.clamps.solid_mask.size() == 4 && b.model.clamps.solid_mask[1] == 1, "solid_mask: NID 1 -> частица 1");
    check(b.sum.clamped_particles == 1 && b.sum.clamped_beads == 2, "сводка: 1 частица, 2 бусины");
}

void test_or_empty()
{
    std::printf("[or/empty]\n");
    Built b = build("*SET_NODE_LIST\n5\n1\n*BOUNDARY_SPC_NODE\n1,0,1\n1,0,0,1\n"
                    "*BOUNDARY_SPC_SET\n5,0,0,0,1\n",
                    false, geom_wo_beams());
    check(b.model.clamps.solid_mask.size() == 4 && b.model.clamps.solid_mask[1] == 7, "повторы объединены по OR");
    Built e = build("", true);
    check(e.model.clamps.empty(), "без SPC маски пусты");
    Built z = build("*BOUNDARY_SPC_NODE\n1,0,0,0,0,1,1,1\n", false, geom_wo_beams());
    check(z.model.clamps.empty(), "нулевые DOFX..Z (вращательные игнорируются) -> пусто");
    check(z.sum.clamped_particles == 0, "сводка: 0");
}

void test_errors()
{
    std::printf("[warn/errors]\n");
    Built b = build("*BOUNDARY_SPC_NODE\n9,0,1,1,1\n1,0,1\n", false, geom_wo_beams());
    check(has_warning(b, "1 узлов"), "узел без роли: предупреждение");
    check(b.model.clamps.solid_mask.size() == 4 && b.model.clamps.solid_mask[1] == 1, "остальные SPC применены");
    // Узел нити в режиме solid - тоже без роли.
    Built s = build("*BOUNDARY_SPC_NODE\n8,0,1\n", false);
    check(has_warning(s, "1 узлов") && s.model.clamps.empty(), "узел нити в режиме solid: без роли");

    std::string e = error_of([] { build("*BOUNDARY_SPC_NODE\n1,5,1,1,1\n", false, geom_wo_beams()); });
    check(has(e, "CID=5") && has(e, "не поддерживается") && has(e, "mem.k"), "CID != 0 (NODE): " + e);
    e = error_of([] { build("*SET_NODE_LIST\n5\n1\n*BOUNDARY_SPC_SET\n5,2,1,1,1\n", false, geom_wo_beams()); });
    check(has(e, "CID=2"), "CID != 0 (SET): " + e);
}

void test_fixed()
{
    std::printf("[fixed]\n");
    // Ширина 10: NID=3, CID пусто, DOFX=0, DOFY=1, DOFZ=1.
    Built b = build("*BOUNDARY_SPC_NODE\n         3                   0         1         1\n", false, geom_wo_beams());
    check(b.model.clamps.solid_mask.size() == 4 && b.model.clamps.solid_mask[0] == 6, "fixed-формат: NID 3 -> y|z");
}

} // namespace

int main()
{
    test_node();
    test_set();
    test_fiber();
    test_or_empty();
    test_errors();
    test_fixed();
    std::printf(failures ? "FAILED: %d\n" : "all ok\n", failures);
    return failures ? 1 : 0;
}
