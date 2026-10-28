// test_k_fiber.cpp
// Unit-тесты сборщика нитей из *ELEMENT_BEAM (readk/k-import-plan.md, этап K5,
// Р1, Р8, Р9) и эквивалентность .fib <-> .k на небольшом тканом полотне:
//   - цепочки с перемешанным порядком балок и обратными балками дают ту же
//     ломаную, что упорядоченная запись; обход от конца с меньшим NID;
//   - ветвление, кольцо, N1 == N2, дубликат ребра, узел и в SPH, и в балке,
//     PID без FIBER_PART, разные PID в цепочке, нет нитей -> ошибки с файлом/строкой;
//   - режим solid пропускает балки с предупреждением;
//   - object_id, material_id и порядок нитей (Р9);
//   - .fib (plain weave) и .k дают побитово одинаковую модель нитей.
// SPH_DIM=3 (нити - только 3D).
#include "fiber/fiber_input.hpp"
#include "fiber/plain_weave.hpp"
#include "k_input/k_geometry.hpp"
#include "simulation_config.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
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

bool bit_eq(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// Таблица: SPH-деталь 1; нити PID 2 -> материал 0, PID 3 -> материал 1.
const char* kTable =
    "SOLID_PART 1 5 1000\n"
    "MATERIAL 0 2.0e9 0.35 1.0e-4 1140 0.02\n"
    "MATERIAL 1 1.0e9 0.30 2.0e-4 900 0.05\n"
    "FIBER_PART 2 0\nFIBER_PART 3 1\n";

// Один SPH-узел (NID 1) вдали + заданные узлы/балки; *PART 1..4.
std::string deck_text(const std::string& nodes, const std::string& beams)
{
    return "*NODE\n1,100,0,0\n" + nodes + "*PART\nA\n1\n*PART\nB\n2\n*PART\nC\n3\n*PART\nD\n4\n"
           "*ELEMENT_SPH\n1,1,1.0\n*ELEMENT_BEAM\n" + beams;
}

struct Built
{
    fiber::Model model;
    KImportSummary sum;
};

Built build(const std::string& text, const std::string& table = kTable)
{
    Built b;
    const KDeck deck = read_k_deck_from_string(text, "mem.k");
    const KPartsTable t = read_k_parts_from_string(table, "t.parts");
    KImportOptions o;
    o.k_file = "mem.k";
    o.mode = KMode::SolidFiber;
    build_fibers_from_deck(deck, t, o, b.model, b.sum);
    return b;
}

bool same_fibers(const fiber::Model& a, const fiber::Model& b)
{
    if (a.fiber.size() != b.fiber.size() || a.topology.bonds.size() != b.topology.bonds.size()
        || a.materials.size() != b.materials.size())
        return false;
    for (std::size_t i = 0; i < a.fiber.size(); ++i)
    {
        for (std::size_t d = 0; d < dim; ++d)
            if (!bit_eq(a.fiber.x[d][i], b.fiber.x[d][i]) || !bit_eq(a.fiber.vx[d][i], b.fiber.vx[d][i])) return false;
        if (!bit_eq(a.fiber.mass[i], b.fiber.mass[i]) || a.fiber.object_id[i] != b.fiber.object_id[i]
            || a.fiber.material_id[i] != b.fiber.material_id[i] || a.topology.prev[i] != b.topology.prev[i]
            || a.topology.next[i] != b.topology.next[i])
            return false;
    }
    for (std::size_t k = 0; k < a.topology.bonds.size(); ++k)
    {
        const auto& x = a.topology.bonds[k];
        const auto& y = b.topology.bonds[k];
        if (x.i != y.i || x.j != y.j || !bit_eq(x.l0, y.l0) || !bit_eq(x.l_break, y.l_break)) return false;
    }
    for (std::size_t m = 0; m < a.materials.size(); ++m)
    {
        const auto& x = a.materials[m];
        const auto& y = b.materials[m];
        if (!bit_eq(x.E, y.E) || !bit_eq(x.nu, y.nu) || !bit_eq(x.d, y.d) || !bit_eq(x.rho, y.rho)
            || !bit_eq(x.eps_break, y.eps_break) || !bit_eq(x.mu_ff, y.mu_ff) || !bit_eq(x.mu_fs, y.mu_fs))
            return false;
    }
    return true;
}

const char* kChainNodes = "10,0,0,0\n5,1,0,0\n7,2.5,0,0\n8,4,0,0\n";

void test_order()
{
    std::printf("[order]\n");
    // Цепочка 10-5-7-8; обход от меньшего NID-конца (8): точки x = 4, 2.5, 1, 0.
    Built a = build(deck_text(kChainNodes, "1,2,10,5\n2,2,5,7\n3,2,7,8\n"));
    check(a.model.fiber.size() == 4 && a.sum.fibers == 1 && a.sum.beads == 4 && a.sum.bonds == 3,
          "одна нить из 4 бусин, 3 связи");
    check(a.model.fiber.x[0][0] == 4.0 && a.model.fiber.x[0][3] == 0.0, "старт с конца с меньшим NID (8)");
    Built b = build(deck_text(kChainNodes, "3,2,8,7\n1,2,5,10\n2,2,5,7\n"));
    check(same_fibers(a.model, b.model), "перемешанный порядок и обратные балки: та же ломаная");
    check(a.sum.seg_min == 1.0 && a.sum.seg_max == 1.5 && a.sum.ld_min > 0.0, "сводка: длины сегментов");

    // Р9: нити по первой балке в файле, object_id/material_id.
    const std::string nodes2 = "11,0,0,0\n12,1,0,0\n13,2,0,0\n21,0,5,0\n22,1,5,0\n";
    Built c = build(deck_text(nodes2, "1,3,22,21\n2,2,11,12\n3,2,12,13\n"));
    check(c.sum.fibers == 2 && c.model.fiber.size() == 5, "две нити, 5 бусин");
    check(c.model.fiber.object_id[0] == 0 && c.model.fiber.object_id[1] == 0 && c.model.fiber.object_id[2] == 1,
          "object_id = номер нити (первая нить - PID 3)");
    check(c.model.fiber.material_id[0] == 1 && c.model.fiber.material_id[1] == 1 && c.model.fiber.material_id[2] == 0,
          "material_id из FIBER_PART");
    check(c.model.fiber.x[1][0] == 5.0 && c.model.fiber.x[0][0] == 0.0 && c.model.fiber.x[0][1] == 1.0,
          "первая нить - от NID 21 (меньший), затем нить PID 2 от NID 11");
    check(c.model.fiber.x[0][2] == 0.0 && c.model.fiber.x[0][4] == 2.0, "вторая нить 11-12-13");
}

void expect_err(const std::string& text, const std::string& sub, const std::string& what,
                const std::string& table = kTable)
{
    const std::string e = error_of([&] { build(text, table); });
    check(!e.empty() && has(e, sub) && has(e, "mem.k:"), what + ": " + e);
}

void test_errors()
{
    std::printf("[errors]\n");
    const std::string nodes = "2,1,0,0\n3,2,0,0\n4,3,0,0\n5,4,0,0\n";
    expect_err(deck_text(nodes, "1,2,2,3\n2,2,3,4\n3,2,3,5\n"), "ветвление", "ветвление");
    expect_err(deck_text(nodes, "1,2,2,3\n2,2,3,4\n3,2,4,2\n"), "кольцо", "кольцо");
    expect_err(deck_text(nodes, "1,2,2,2\n"), "N1 == N2", "N1 == N2");
    expect_err(deck_text(nodes, "1,2,2,3\n2,2,3,2\n"), "дубликат ребра", "дубликат ребра");
    expect_err(deck_text(nodes, "1,2,2,3\n2,3,3,4\n"), "разными PID", "разные PID в цепочке");
    expect_err(deck_text(nodes, "1,4,2,3\n"), "FIBER_PART", "PID балки без FIBER_PART");
    // Узел и в SPH, и в балке (NID 1 - SPH).
    expect_err(deck_text(nodes, "1,2,1,2\n"), "одновременно", "узел SPH и балки");
    // Нет нитей.
    std::string e = error_of([&] { build(deck_text(nodes, "")); });
    check(has(e, "нет ни одной нити") && has(e, "mem.k"), "нет балок: " + e);
    // Нулевой сегмент: префикс не от .fib.
    e = error_of([&] { build(deck_text("2,1,0,0\n3,1,0,0\n", "1,2,2,3\n")); });
    check(has(e, "нулевая") && has(e, "mem.k") && !has(e, "load_fiber_file"), "нулевой сегмент: нейтральный префикс: " + e);
}

void test_solid_mode()
{
    std::printf("[solid mode]\n");
    const KDeck deck = read_k_deck_from_string(deck_text(kChainNodes, "1,2,10,5\n2,2,5,7\n"), "mem.k");
    const KPartsTable t = read_k_parts_from_string("SOLID_PART 1 5 1000\n", "t.parts");
    SimulationConfig cfg;
    cfg.apply_defaults();
    ParticleSystem sys;
    std::size_t moni = 0;
    KImportSummary s;
    KImportOptions o;
    o.mode = KMode::Solid;
    build_solid_from_deck(deck, &t, o, cfg.mats, sys, moni, s);
    bool warn = false;
    for (const auto& w : s.warnings)
        if (has(w, "ELEMENT_BEAM") && has(w, "PID: 2")) warn = true;
    check(sys.size() == 1 && warn, "solid: балки пропущены с предупреждением");
}

// ---- Эквивалентность .fib <-> .k ---------------------------------------
void test_equivalence()
{
    std::printf("[equivalence .fib <-> .k]\n");
    fiber::PlainWeaveParams q;
    q.nx = 3;
    q.ny = 2;
    q.seg = 3;
    q.clamp = fiber::PlainWeaveClamp::None;
    const fiber::PlainWeave W = fiber::build_plain_weave(q);

    std::ostringstream fib;
    fiber::write_plain_weave_fib(fib, W, "test_k_fiber");
    {
        std::ofstream("__k_fiber_w.fib") << fib.str();
    }
    fiber::Model ma;
    fiber::load_fiber_file("__k_fiber_w.fib", ma);

    // MATERIAL - та же строка, что в .fib.
    std::string mat_line;
    {
        std::istringstream in(fib.str());
        std::string ln;
        while (std::getline(in, ln))
            if (ln.rfind("MATERIAL", 0) == 0) mat_line = ln;
    }

    std::ostringstream k;
    char buf[160];
    k << "*NODE\n1,100,0,0\n";
    long long nid = 1000;
    std::ostringstream beams;
    long long eid = 1;
    for (const auto& F : W.fibers)
    {
        const long long base = nid;
        for (const auto& P : F.pts)
        {
            std::snprintf(buf, sizeof buf, "%lld,%.17g,%.17g,%.17g\n", nid++, P[0], P[1], P[2]);
            k << buf;
        }
        for (std::size_t j = 0; j + 1 < F.pts.size(); ++j)
        {
            // Половина балок записана в обратном направлении.
            const long long a = base + static_cast<long long>(j), b = a + 1;
            beams << eid++ << ",2," << ((j % 2) ? b : a) << "," << ((j % 2) ? a : b) << "\n";
        }
    }
    k << "*PART\nS\n1\n*PART\nF\n2\n*ELEMENT_SPH\n1,1,1.0\n*ELEMENT_BEAM\n" << beams.str();
    { std::ofstream("__k_fiber_w.k") << k.str(); }
    { std::ofstream("__k_fiber_w.parts") << "SOLID_PART 1 5 1000\n" << mat_line << "\nFIBER_PART 2 0\n"; }

    SimulationConfig cfg;
    cfg.apply_defaults();
    fiber::Model mb;
    std::size_t moni = 0;
    KImportOptions o;
    o.k_file = "__k_fiber_w.k";
    o.parts_file = "__k_fiber_w.parts";
    o.mode = KMode::SolidFiber;
    const KImportSummary s = load_k_geometry(o, cfg.mats, mb, moni);
    check(ma.fiber.size() > 0 && mb.solid.size() == 1, "модель собрана, Solid из 1 частицы");
    check(s.fibers == W.fibers.size() && s.beads == ma.fiber.size() && s.bonds == ma.topology.bonds.size(),
          "сводка нитей совпадает с .fib");
    check(same_fibers(ma, mb), "fiber.x, mass, object_id, material_id, prev/next, bonds побитово");

    std::remove("__k_fiber_w.fib");
    std::remove("__k_fiber_w.k");
    std::remove("__k_fiber_w.parts");
}

} // namespace

int main()
{
    try
    {
        test_order();
        test_errors();
        test_solid_mode();
        test_equivalence();
    }
    catch (const std::exception& e)
    {
        std::printf("исключение: %s\n", e.what());
        return 1;
    }
    if (failures)
    {
        std::printf("\nFAILED: %d\n", failures);
        return 1;
    }
    std::printf("\nвсе проверки пройдены\n");
    return 0;
}
