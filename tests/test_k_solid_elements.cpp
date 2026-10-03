// test_k_solid_elements.cpp
// Unit-тесты *ELEMENT_SOLID -> частицы (readk/k-import-plan.md, этап K9):
//   - куб 2x2x2 из 8 гексаэдров: V = a^3, центры в центрах ячеек;
//   - тетраэдр: V = |det|/6, центр - среднее вершин; клин (N4=N3, N8=N7) - объём призмы;
//   - 10-узловой тетраэдр (середины игнорируются);
//   - режим Skip: предупреждение и 0 частиц-солидов;
//   - нулевой объём, PID без SOLID_PART, нет таблицы, узел и в SPH -> ошибки;
//   - SPH + солиды: порядок частиц и размер solid_mask при SPC на SPH-узле.
// SPH_DIM=3.
#include "fiber/model.hpp"
#include "k_input/k_geometry.hpp"
#include "param.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <cmath>
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

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(std::fabs(b), 1e-300); }

const char* kTable = "SOLID_PART 1 5 1000\nSOLID_PART 2 10 2510 hsml=0.5\n";

struct Result
{
    fiber::Model model;
    std::size_t moni = 0;
    KImportSummary sum;
    const ParticleSystem& sys() const { return model.solid; }
};

Result build(const std::string& deck_text, bool particles = true, const char* table = kTable)
{
    Result r;
    const KDeck deck = read_k_deck_from_string(deck_text, "mem.k");
    SimulationConfig cfg;
    cfg.apply_defaults();
    KImportOptions o;
    o.k_file = "mem.k";
    if (particles) o.solid_elements = KSolidElements::Particles;
    if (table)
    {
        const KPartsTable t = read_k_parts_from_string(table, "t.parts");
        check_parts_against_deck(t, deck, false, r.sum.warnings, particles);
        build_solid_from_deck(deck, &t, o, cfg.mats, r.model.solid, r.moni, r.sum);
    }
    else
        build_solid_from_deck(deck, nullptr, o, cfg.mats, r.model.solid, r.moni, r.sum);
    apply_spc_clamps(deck, {}, r.model, r.sum);
    return r;
}

const char* kParts = "*PART\nA\n1\n*PART\nB\n2\n";

std::string cube_deck(double a)
{
    std::string s = "*NODE\n";
    auto id = [](int i, int j, int k) { return 1 + i + 3 * j + 9 * k; };
    for (int k = 0; k < 3; ++k)
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                s += std::to_string(id(i, j, k)) + "," + std::to_string(i * a) + "," + std::to_string(j * a) + ","
                     + std::to_string(k * a) + "\n";
    s += kParts;
    s += "*ELEMENT_SOLID\n";
    int eid = 1;
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i)
                s += std::to_string(eid++) + ",1," + std::to_string(id(i, j, k)) + "," + std::to_string(id(i + 1, j, k))
                     + "," + std::to_string(id(i + 1, j + 1, k)) + "," + std::to_string(id(i, j + 1, k)) + ","
                     + std::to_string(id(i, j, k + 1)) + "," + std::to_string(id(i + 1, j, k + 1)) + ","
                     + std::to_string(id(i + 1, j + 1, k + 1)) + "," + std::to_string(id(i, j + 1, k + 1)) + "\n";
    return s;
}

void test_cube()
{
    std::printf("[cube]\n");
    const double a = 0.375;   // степень двойки: координаты в тексте точные
    Result r = build(cube_deck(a));
    const ParticleSystem& s = r.sys();
    check(s.size() == 8 && r.sum.solid_element_particles == 8 && r.sum.sph_particles == 0, "8 частиц");
    bool ok_m = true, ok_x = true;
    for (std::size_t e = 0; e < s.size() && e < 8; ++e)
    {
        ok_m = ok_m && near_rel(s.mass[e], 1000.0 * a * a * a, 1e-13) && s.itype[e] == 5 && s.body[e] == 1
               && s.rho[e] == 1000.0;
        const int i = static_cast<int>(e) & 1, j = (static_cast<int>(e) >> 1) & 1, k = (static_cast<int>(e) >> 2) & 1;
        const double c[3] = {(i + 0.5) * a, (j + 0.5) * a, (k + 0.5) * a};
        for (int d = 0; d < 3; ++d) ok_x = ok_x && std::fabs(s.x[d][e] - c[d]) < 1e-14;
    }
    check(ok_m, "mass = rho a^3 (1e-13), itype, body, rho");
    check(ok_x, "центры в центрах ячеек, порядок файла");
    check(near_rel(s.hsml[0], 2.0 * a, 1e-12), "hsml = 2 V^(1/3)");
    check(near_rel(r.sum.dx_median, a, 1e-12), "dx_median = a");
    check(near_rel(r.sum.bbox_max[0], 1.5 * a, 1e-12) && near_rel(r.sum.bbox_min[0], 0.5 * a, 1e-12),
          "габарит по частицам");
}

void test_shapes()
{
    std::printf("[shapes]\n");
    // Тетраэдр.
    Result t = build(std::string("*NODE\n1,0,0,0\n2,2,0,0\n3,0,3,0\n4,1,1,4\n") + kParts
                     + "*ELEMENT_SOLID\n1,2,1,2,3,4\n");
    const double V = 2.0 * 3.0 * 4.0 / 6.0;
    check(t.sys().size() == 1 && near_rel(t.sys().mass[0], 2510.0 * V, 1e-14), "тетраэдр: V = det/6");
    check(std::fabs(t.sys().x[0][0] - 0.75) < 1e-14 && std::fabs(t.sys().x[1][0] - 1.0) < 1e-14
              && std::fabs(t.sys().x[2][0] - 1.0) < 1e-14,
          "центр - среднее вершин");
    check(t.sys().hsml[0] == 0.5, "hsml из таблицы");
    // Клин: треугольник 1,2,3 внизу, 5,6,7 вверху; N4=N3, N8=N7. Объём = 0.5*2*3*h, h=2.
    Result w = build(std::string("*NODE\n1,0,0,0\n2,2,0,0\n3,0,3,0\n5,0,0,2\n6,2,0,2\n7,0,3,2\n") + kParts
                     + "*ELEMENT_SOLID\n1,1,1,2,3,3,5,6,7,7\n");
    check(w.sys().size() == 1 && near_rel(w.sys().mass[0], 1000.0 * 6.0, 1e-13), "клин: объём призмы");
    check(std::fabs(w.sys().x[2][0] - 1.0) < 1e-13 && std::fabs(w.sys().x[0][0] - 2.0 / 3.0) < 1e-13,
          "клин: центр масс");
    // 10 узлов: середины (произвольные) игнорируются.
    Result q = build(std::string("*NODE\n1,0,0,0\n2,2,0,0\n3,0,3,0\n4,1,1,4\n5,9,9,9\n") + kParts
                     + "*ELEMENT_SOLID (ten nodes format)\n1,2\n1,2,3,4,5,5,5,5,5,5\n");
    check(q.sys().size() == 1 && near_rel(q.sys().mass[0], 2510.0 * V, 1e-14), "10 узлов: углы N1-N4");
}

void test_skip_and_errors()
{
    std::printf("[skip/errors]\n");
    const std::string tet = std::string("*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n") + kParts;
    Result s = build(tet + "*ELEMENT_SOLID\n1,1,1,2,3,4\n", false);
    bool w = false;
    for (const auto& x : s.sum.warnings) w = w || has(x, "ELEMENT_SOLID");
    check(s.sys().size() == 0 && s.sum.solid_element_particles == 0 && w, "Skip: предупреждение, 0 частиц");

    std::string e = error_of([&] {
        build(std::string("*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,1,1,0\n") + kParts + "*ELEMENT_SOLID\n1,1,1,2,3,4\n");
    });
    check(has(e, "объём") && has(e, "EID 1") && has(e, "mem.k"), "нулевой объём: " + e);
    e = error_of([&] { build(tet + "*ELEMENT_SOLID\n1,2,1,2,3,4\n", true, "SOLID_PART 1 5 1000\n"); });
    check(has(e, "PID 2") && has(e, "SOLID_PART"), "PID без SOLID_PART: " + e);
    e = error_of([&] { build(tet + "*ELEMENT_SOLID\n1,1,1,2,3,4\n", true, nullptr); });
    check(has(e, "таблиц"), "нет таблицы: " + e);
    e = error_of([&] { build(tet + "*ELEMENT_SPH\n1,1,1.0\n*ELEMENT_SOLID\n1,1,1,2,3,4\n"); });
    check(has(e, "ELEMENT_SPH"), "узел и в SPH, и в солиде: " + e);
}

void test_mixed()
{
    std::printf("[mixed]\n");
    const std::string d = std::string("*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n10,5,5,5\n11,6,5,5\n") + kParts
                          + "*ELEMENT_SPH\n11,1,2.0\n10,1,3.0\n*ELEMENT_SOLID\n1,1,1,2,3,4\n"
                            "*BOUNDARY_SPC_NODE\n10,0,0,1\n1,0,1,1,1\n";
    Result r = build(d);
    const ParticleSystem& s = r.sys();
    check(s.size() == 3 && r.sum.sph_particles == 2 && r.sum.solid_element_particles == 1, "2 SPH + 1 солид");
    check(s.mass[0] == 2.0 && s.mass[1] == 3.0 && near_rel(s.mass[2], 1000.0 / 6.0, 1e-14),
          "SPH первыми, затем солид");
    const auto& m = r.model.clamps.solid_mask;
    check(m.size() == 3 && m[0] == 0 && m[1] == 2 && m[2] == 0, "solid_mask размера 3, SPC на SPH-узле 10");
    bool w = false;
    for (const auto& x : r.sum.warnings) w = w || has(x, "ни SPH-частицей");
    check(w, "SPC на узле солида: предупреждение без роли");
    check(near_rel(r.sum.bbox_max[0], 6.0, 1e-14) && r.sum.bbox_min[0] < 0.3, "габарит по всем частицам");
}

} // namespace

int main()
{
    test_cube();
    test_shapes();
    test_skip_and_errors();
    test_mixed();
    std::printf(failures ? "FAILED: %d\n" : "all ok\n", failures);
    return failures ? 1 : 0;
}
