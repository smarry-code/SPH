// test_fiber_model.cpp
// Unit-тесты модели данных нитей (fiber/PLAN.md, этап Ф2):
//   - разбор файла ломаных (load_fiber_file): число частиц, топология
//     prev/next, l0/l_break на связях, масса бусин по формуле fiber/SPH.md;
//   - ошибки формата файла нитей;
//   - fiber::fiber_output — файл создаётся и содержит ожидаемые значения;
//   - fiber::check_fiber_requires_dim3() не бросает при dim == 3 (ветка
//     throw при dim != 3 недостижима здесь компилятивно — её проверяет
//     отдельная цель test_fiber_dim_guard, собранная без SPH_DIM=3);
//   - body[] в ParticleSystem: из PID k-файла (parse_k_file) и из itype[]
//     во встроенном генераторе (shock_tube).
//
// Собирается с SPH_DIM=3 (CMakeLists.txt): нити (fiber/SPH.md) поддержаны
// только в этой размерности.
#include "fiber/fiber_input.hpp"
#include "fiber/fiber_output.hpp"
#include "fiber/model.hpp"
#include "input.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_model требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near(double a, double b, double tol = 1e-9)
{
    return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

// Пишет содержимое content в файл filename (перезаписывает).
void write_file(const std::string& filename, const std::string& content)
{
    std::ofstream out(filename);
    out << content;
}

// --- Фикстура: две ломаные, разное число точек и материалы -----------------
// Нить A (материал 0): 4 точки на оси x с равным шагом L = 0.01 -> 3 связи.
// Нить B (материал 1): 2 точки, несимметричный отрезок -> 1 связь.
const char* kFixturePath = "__fiber_model_fixture.fib";
const double kL = 0.01;              // шаг нити A
const double kMatA_d = 1.0e-3, kMatA_rho = 1200.0, kMatA_eps = 0.05;
const double kMatB_d = 2.0e-3, kMatB_rho = 900.0,  kMatB_eps = 0.10;

std::string fixture_content()
{
    std::ostringstream s;
    s << "# fixture: две нити, разные материалы и число точек\n";
    s << "MATERIAL 0 2.0e9 0.35 " << kMatA_d << " " << kMatA_rho << " " << kMatA_eps << "\n";
    s << "MATERIAL 1 1.5e9 0.30 " << kMatB_d << " " << kMatB_rho << " " << kMatB_eps << "\n";
    s << "FIBER 0 4\n";
    s << "0.0 0.0 0.0\n";
    s << kL << " 0.0 0.0\n";
    s << (2 * kL) << " 0.0 0.0\n";
    s << (3 * kL) << " 0.0 0.0\n";
    s << "FIBER 1 2\n";
    s << "0.0 0.01 0.0\n";
    s << "0.0 0.02 0.005\n";
    return s.str();
}

void test_topology_and_counts()
{
    std::printf("test_topology_and_counts:\n");
    write_file(kFixturePath, fixture_content());

    fiber::Model model;
    fiber::load_fiber_file(kFixturePath, model);

    check(model.fiber.size() == 6, "6 бусин (4 + 2)");
    check(model.topology.bonds.size() == 4, "4 связи (3 + 1)");
    check(model.materials.size() == 2, "2 материала в таблице");

    // Нить A: индексы 0..3. Концы (0, 3) — ровно один сосед -1; внутренние
    // (1, 2) — оба соседа не -1.
    check(model.topology.prev[0] == -1 && model.topology.next[0] != -1, "бусина 0 — конец (prev=-1)");
    check(model.topology.prev[3] != -1 && model.topology.next[3] == -1, "бусина 3 — конец (next=-1)");
    check(model.topology.prev[1] != -1 && model.topology.next[1] != -1, "бусина 1 — внутренняя");
    check(model.topology.prev[2] != -1 && model.topology.next[2] != -1, "бусина 2 — внутренняя");
    check(model.topology.next[0] == 1 && model.topology.prev[1] == 0, "связность 0-1");
    check(model.topology.next[1] == 2 && model.topology.prev[2] == 1, "связность 1-2");
    check(model.topology.next[2] == 3 && model.topology.prev[3] == 2, "связность 2-3");

    // Нить B: индексы 4..5, оба — концы (n=2, внутренних бусин нет).
    check(model.topology.prev[4] == -1 && model.topology.next[4] == 5, "бусина 4 — начало нити B");
    check(model.topology.prev[5] == 4 && model.topology.next[5] == -1, "бусина 5 — конец нити B");

    // object_id: нить A -> 0, нить B -> 1 (порядок блоков FIBER в файле).
    check(model.fiber.object_id[0] == 0 && model.fiber.object_id[3] == 0, "object_id нити A = 0");
    check(model.fiber.object_id[4] == 1 && model.fiber.object_id[5] == 1, "object_id нити B = 1");
    check(model.fiber.material_id[0] == 0, "material_id нити A = 0");
    check(model.fiber.material_id[4] == 1, "material_id нити B = 1");
}

void test_masses()
{
    std::printf("test_masses:\n");
    fiber::Model model;
    fiber::load_fiber_file(kFixturePath, model);

    const double A = pi * kMatA_d * kMatA_d / 4.0;
    const double mass_end = kMatA_rho * A * kL / 2.0;
    const double mass_mid = kMatA_rho * A * kL;

    check(near(model.fiber.mass[0], mass_end), "масса конца нити A = rho*A*L/2");
    check(near(model.fiber.mass[1], mass_mid), "масса внутренней бусины A = rho*A*L");
    check(near(model.fiber.mass[2], mass_mid), "масса внутренней бусины A (2) = rho*A*L");
    check(near(model.fiber.mass[3], mass_end), "масса конца нити A (3) = rho*A*L/2");

    // Нить B: единственный сегмент длиной sqrt(0.01^2+0.005^2), обе бусины —
    // концы этого сегмента, каждая получает половину его массы.
    const double lB = std::sqrt(0.01 * 0.01 + 0.005 * 0.005);
    const double AB = pi * kMatB_d * kMatB_d / 4.0;
    const double massB_end = kMatB_rho * AB * lB / 2.0;
    check(near(model.fiber.mass[4], massB_end), "масса бусины 4 нити B");
    check(near(model.fiber.mass[5], massB_end), "масса бусины 5 нити B");
}

void test_l0_and_l_break()
{
    std::printf("test_l0_and_l_break:\n");
    fiber::Model model;
    fiber::load_fiber_file(kFixturePath, model);

    // Первые 3 связи принадлежат нити A (материал 0, шаг L, eps_break = 0.05).
    for (size_t k = 0; k < 3; ++k)
    {
        const fiber::FiberBond& b = model.topology.bonds[k];
        check(near(b.l0, kL), "l0 связи нити A = L (k=" + std::to_string(k) + ")");
        check(near(b.l_break, (1.0 + kMatA_eps) * kL), "l_break = (1+eps_break)*l0 (нить A, k=" + std::to_string(k) + ")");
        check(b.i == k && b.j == k + 1, "индексы связи (нить A, k=" + std::to_string(k) + ")");
    }

    // 4-я связь — нить B (материал 1, eps_break = 0.10), несимметричный отрезок.
    const double lB = std::sqrt(0.01 * 0.01 + 0.005 * 0.005);
    const fiber::FiberBond& bB = model.topology.bonds[3];
    check(near(bB.l0, lB), "l0 связи нити B = |p1-p0|");
    check(near(bB.l_break, (1.0 + kMatB_eps) * lB), "l_break нити B = (1+eps_break)*l0");
    check(bB.i == 4 && bB.j == 5, "индексы связи нити B");
}

void test_fiber_output_writes_expected_values()
{
    std::printf("test_fiber_output_writes_expected_values:\n");
    fiber::Model model;
    fiber::load_fiber_file(kFixturePath, model);

    const size_t itimestep = 7;
    fiber::fiber_output(model, itimestep, 1.25e-6);

    std::ifstream in("fiber_output_00007.txt");
    check(in.is_open(), "файл fiber_output_00007.txt создан");
    if (!in.is_open()) return;

    std::string header_line;
    std::getline(in, header_line);
    std::istringstream hs(header_line);
    std::string kw; double t; size_t step, n, nb;
    hs >> kw >> t >> kw >> step >> kw >> n >> kw >> nb;
    check(near(t, 1.25e-6, 1e-6), "TIME записано верно");
    check(step == itimestep, "STEP записано верно");
    check(n == 6, "N_PARTICLES == 6 в заголовке");
    check(nb == 4, "N_BONDS == 4 в заголовке");

    std::string line;
    std::getline(in, line);              // "# id object_id ..."
    std::getline(in, line);
    check(line == "PARTICLES", "маркер блока PARTICLES");

    size_t particle_lines = 0;
    while (true)
    {
        if (!std::getline(in, line)) break;
        if (line == "# i j l0 l l_break elongation") break;
        ++particle_lines;
    }
    check(particle_lines == 6, "6 строк в блоке PARTICLES");

    std::getline(in, line);
    check(line == "BONDS", "маркер блока BONDS");

    // Первая связь: свежезагруженная модель, l должно совпадать с l0 -> удлинение ~0.
    std::getline(in, line);
    std::istringstream bs(line);
    size_t bi, bj; double l0, l, lbreak, elong;
    bs >> bi >> bj >> l0 >> l >> lbreak >> elong;
    check(bi == 0 && bj == 1, "первая связь в выводе — (0,1)");
    check(near(l, l0, 1e-4), "l == l0 для только что загруженной нити (без деформации)");
    check(std::fabs(elong) < 1e-4, "удлинение ~0 для только что загруженной нити");
}

void test_dim3_guard_does_not_throw()
{
    std::printf("test_dim3_guard_does_not_throw:\n");
    bool threw = false;
    try { fiber::check_fiber_requires_dim3(); }
    catch (const std::runtime_error&) { threw = true; }
    check(!threw, "check_fiber_requires_dim3() не бросает при dim == 3");
}

void test_malformed_files_throw()
{
    std::printf("test_malformed_files_throw:\n");

    // Ссылка на неопределённый материал.
    {
        write_file("__fiber_bad_material.fib", "FIBER 0 2\n0 0 0\n1 0 0\n");
        fiber::Model model;
        bool threw = false;
        try { fiber::load_fiber_file("__fiber_bad_material.fib", model); }
        catch (const std::runtime_error&) { threw = true; }
        check(threw, "FIBER на неопределённый MATERIAL бросает исключение");
    }

    // Меньше двух точек.
    {
        write_file("__fiber_bad_points.fib", "MATERIAL 0 1e9 0.3 1e-3 1000 0.05\nFIBER 0 1\n0 0 0\n");
        fiber::Model model;
        bool threw = false;
        try { fiber::load_fiber_file("__fiber_bad_points.fib", model); }
        catch (const std::runtime_error&) { threw = true; }
        check(threw, "нить с одной точкой бросает исключение");
    }

    // Несуществующий файл.
    {
        fiber::Model model;
        bool threw = false;
        try { fiber::load_fiber_file("__fiber_does_not_exist.fib", model); }
        catch (const std::runtime_error&) { threw = true; }
        check(threw, "отсутствующий файл бросает исключение");
    }
}

// --- body[] --------------------------------------------------------------
void test_body_from_k_file()
{
    std::printf("test_body_from_k_file:\n");
    const char* kfile = "__fiber_model_body.k";
    write_file(kfile,
        "*NODE\n"
        "1,0.0,0.0,0.0\n"
        "2,1.0,0.0,0.0\n"
        "3,2.0,0.0,0.0\n"
        "*ELEMENT_SPH\n"
        "1,5,1.0\n"
        "2,7,2.0\n"
        "3,5,3.0\n");

    SimulationConfig cfg;
    cfg.apply_defaults();
    ParticleSystem sys;
    parse_k_file(cfg, kfile, sys);

    check(sys.size() == 3, "3 частицы разобраны из k-файла");
    const int expected_pid[3] = {5, 7, 5};
    for (size_t i = 0; i < 3 && i < sys.size(); ++i)
    {
        check(sys.body[i] == expected_pid[i], "body[" + std::to_string(i) + "] == PID k-файла");
    }
}

void test_body_from_shock_tube()
{
    std::printf("test_body_from_shock_tube:\n");
    SimulationConfig cfg;
    cfg.apply_defaults();
    ParticleSystem sys;
    size_t ntotal = 0;
    double dx = 0.0;
    shock_tube(cfg, sys, ntotal, dx);

    check(ntotal > 0, "shock_tube сгенерировал частицы");
    bool all_match = true;
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (sys.body[i] != sys.itype[i]) { all_match = false; break; }
    }
    check(all_match, "body[i] == itype[i] для всех частиц shock_tube");
}

} // namespace

int main()
{
    test_topology_and_counts();
    test_masses();
    test_l0_and_l_break();
    test_fiber_output_writes_expected_values();
    test_dim3_guard_does_not_throw();
    test_malformed_files_throw();
    test_body_from_k_file();
    test_body_from_shock_tube();

    if (failures == 0)
    {
        std::printf("test_fiber_model: OK\n");
    }
    else
    {
        std::printf("test_fiber_model: %d failure(s)\n", failures);
    }
    return failures;
}
