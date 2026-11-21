// tests/test_fiber_stitch.cpp
// Связи-стежки (notes/fabrics/plan-quilt-fabric.md, разделы 10-12, этап 2a): директива
// BOND в .fib, собственный материал связи (FiberBond::mat), только растяжение,
// разрыв связи вне цепочки (apply_fiber_break_events не трогает prev/next и
// object_id). Стиль -- по образцу tests/test_fiber_bond.cpp и
// tests/test_fiber_model.cpp (фикстура .fib на диске). Собирается с SPH_DIM=3.
// Round-trip чекпоинта со стежком -- в tests/test_checkpoint.cpp.
#include "fiber/fiber_input.hpp"
#include "fiber/fiber_system.hpp"
#include "fiber/model.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/postprocessing.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_stitch требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b)); }

void write_file(const std::string& filename, const std::string& content)
{
    std::ofstream out(filename);
    out << content;
}

// Материалы: 0 -- нить, 1 -- стежок (другие E, d, eps_break).
const double kE0 = 2.0e9, kD0 = 1.0e-3, kRho0 = 1200.0, kEps0 = 0.05;
const double kE1 = 5.0e8, kD1 = 3.0e-4, kRho1 = 1000.0, kEps1 = 0.2;
const double kL = 0.01;    // длина нити внутри слоя
const double kGap = 0.002; // расстояние между слоями (по z)
const char* kPath = "__fiber_stitch_fixture.fib";

// Два «слоя» по одной нити из 2 бусин: бусины 0,1 -- слой 1 (z = 0), 2,3 -- слой 2
// (z = kGap). Стежок в тестах: 1 -- 3 (вертикальная связь).
std::string fixture(const std::string& bonds)
{
    std::ostringstream s;
    s.precision(17);
    s << "MATERIAL 0 " << kE0 << " 0.35 " << kD0 << " " << kRho0 << " " << kEps0 << "\n";
    s << "MATERIAL 1 " << kE1 << " 0.30 " << kD1 << " " << kRho1 << " " << kEps1 << "\n";
    s << "FIBER 0 2\n0 0 0\n" << kL << " 0 0\n";
    s << "FIBER 0 2\n0 0 " << kGap << "\n" << kL << " 0 " << kGap << "\n";
    s << bonds;
    return s.str();
}

// Возвращает текст ошибки загрузки ("" -- успех).
std::string load_error(const std::string& content)
{
    write_file(kPath, content);
    fiber::Model m;
    try
    {
        fiber::load_fiber_file(kPath, m);
    }
    catch (const std::runtime_error& e)
    {
        return e.what();
    }
    return "";
}

// Сдвиг второго слоя (бусины 2, 3) по z.
void shift_layer2(fiber::Model& m, double dz)
{
    m.fiber.x[2][2] += dz;
    m.fiber.x[2][3] += dz;
}

void test_parsing()
{
    std::printf("test_parsing:\n");
    write_file(kPath, fixture("BOND 1 3 1\n"));
    fiber::Model m;
    fiber::load_fiber_file(kPath, m);

    check(m.topology.bonds.size() == 3, "3 связи: 2 цепочечные + 1 стежок");
    const fiber::FiberBond& c0 = m.topology.bonds[0];
    const fiber::FiberBond& c1 = m.topology.bonds[1];
    const fiber::FiberBond& s = m.topology.bonds[2];
    check(c0.mat == -1 && c1.mat == -1, "цепочечные связи: mat == -1");
    check(fiber::bond_material_id(c0, m.fiber) == 0 && fiber::bond_material_id(c1, m.fiber) == 0,
          "цепочечные связи берут материал бусины (0)");
    check(s.i == 1 && s.j == 3 && s.mat == 1, "стежок: i=1, j=3, mat=1");
    check(near_rel(s.l0, kGap, 1e-14), "l0 стежка = расстояние между бусинами");
    check(near_rel(s.l_break, (1.0 + kEps1) * kGap, 1e-14), "l_break = (1 + eps_break материала 1) * l0");
    check(fiber::bond_material_id(s, m.fiber) == 1, "материал стежка = 1");
    check(fiber::is_chain_bond(c0, m.topology) && fiber::is_chain_bond(c1, m.topology),
          "цепочечные связи -- в цепочке");
    check(!fiber::is_chain_bond(s, m.topology), "стежок вне цепочки");
    check(m.topology.next[0] == 1 && m.topology.prev[1] == 0 && m.topology.next[1] == -1 &&
              m.topology.prev[2] == -1 && m.topology.next[2] == 3 && m.topology.prev[3] == 2,
          "стежок не меняет prev/next");
    check(m.fiber.object_id[0] == 0 && m.fiber.object_id[3] == 1, "object_id нитей без изменений");

    // Без BOND -- только цепочечные связи, массы те же.
    write_file(kPath, fixture(""));
    fiber::Model m0;
    fiber::load_fiber_file(kPath, m0);
    check(m0.topology.bonds.size() == 2, "без BOND: 2 связи");
    check(m0.fiber.mass == m.fiber.mass, "масса бусин не зависит от BOND");
}

void test_parser_errors()
{
    std::printf("test_parser_errors:\n");
    std::string e = load_error(fixture("BOND 1 4 1\n"));
    check(!e.empty() && e.find("диапазона") != std::string::npos, "индекс вне диапазона -> ошибка");
    check(e.find("строке") != std::string::npos, "в сообщении есть номер строки");
    e = load_error(fixture("BOND -1 3 1\n"));
    check(!e.empty(), "отрицательный индекс -> ошибка");
    e = load_error(fixture("BOND 2 2 1\n"));
    check(!e.empty() && e.find("саму с собой") != std::string::npos, "i == j -> ошибка");
    e = load_error(fixture("BOND 1 3 7\n"));
    check(!e.empty() && e.find("неопределённый материал") != std::string::npos, "неопределённый материал -> ошибка");
    e = load_error(fixture("BOND 1 3\n"));
    check(!e.empty(), "не хватает полей -> ошибка");
    e = load_error(fixture("BOND 1 3 1 5\n"));
    check(!e.empty() && e.find("лишние") != std::string::npos, "лишний токен -> ошибка");
    {
        // Совпадающие бусины: l0 == 0.
        const std::string s =
            "MATERIAL 0 2e9 0.35 1e-3 1200 0.05\nFIBER 0 2\n0 0 0\n0.01 0 0\nFIBER 0 2\n0 0 0\n0.01 0 0\nBOND 0 2 0\n";
        e = load_error(s);
        check(!e.empty() && e.find("l0") != std::string::npos, "совпадающие бусины (l0 = 0) -> ошибка");
    }
    e = load_error(fixture("BOND 1 3 1\n"));
    check(e.empty(), "корректный BOND загружается");
}

void test_force_and_compression()
{
    std::printf("test_force_and_compression:\n");
    write_file(kPath, fixture("BOND 1 3 1\n"));
    fiber::Model m;
    fiber::load_fiber_file(kPath, m);
    const double l0 = m.topology.bonds[2].l0;
    const double A1 = pi * kD1 * kD1 / 4.0;

    rhs::CpuExecutor exec(1);
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;

    // Растяжение: слой 2 уходит вверх на 10% l0 (цепочечные связи по x не меняются).
    const double dz = 0.1 * l0;
    shift_layer2(m, dz);
    std::vector<fiber::FiberDelta> out = interactions::compute_fiber_bond(exec, m.fiber, m.topology, m.materials, buffers);
    const double F = kE1 * A1 * dz / l0; // E1*A1*(l-l0)/l0 по материалу стежка
    check(near_rel(out[1].dv[2] * m.fiber.mass[1], F, 1e-9), "растяжение: сила на бусину 1 = E1*A1*(l-l0)/l0 (по z)");
    check(near_rel(out[3].dv[2] * m.fiber.mass[3], -F, 1e-9), "растяжение: сила на бусину 3 равна и противоположна");
    check(out[0].dv[2] == 0.0 && out[2].dv[2] == 0.0, "бусины без стежка не ускоряются по z");
    check(out[1].dv[0] == 0.0 && out[1].dv[1] == 0.0, "сила стежка направлена вдоль связи");

    // Контроль: по материалу бусины (0) сила была бы другой.
    const double A0 = pi * kD0 * kD0 / 4.0;
    check(!near_rel(F, kE0 * A0 * dz / l0, 1e-3), "сила по материалу 1 отличается от силы по материалу 0");

    // Сжатие: слой 2 приближен -> стежок силы не даёт.
    fiber::Model mc;
    fiber::load_fiber_file(kPath, mc);
    shift_layer2(mc, -0.3 * l0);
    out = interactions::compute_fiber_bond(exec, mc.fiber, mc.topology, mc.materials, buffers);
    bool zero = true;
    for (std::size_t k = 0; k < 4; ++k)
        for (std::size_t d = 0; d < dim; ++d) zero = zero && out[k].dv[d] == 0.0;
    check(zero, "сжатие: стежок силы не даёт");

    const double dt = interactions::compute_fiber_bond_dt_limit(exec, m.fiber, m.topology, m.materials);
    check(dt > 0.0 && std::isfinite(dt), "dt_limit конечен и положителен");
}

void test_break()
{
    std::printf("test_break:\n");
    write_file(kPath, fixture("BOND 1 3 1\n"));
    fiber::Model m;
    fiber::load_fiber_file(kPath, m);
    const double l0 = m.topology.bonds[2].l0;

    rhs::CpuExecutor exec(1);
    std::vector<std::vector<interactions::FiberBreakEvent>> tb;

    // Ниже порога (1 + eps_break1)*l0: событий нет.
    shift_layer2(m, 0.5 * kEps1 * l0);
    auto ev = interactions::find_fiber_break_events(exec, m.fiber, m.topology, tb);
    check(ev.empty(), "ниже порога разрыва событий нет");

    // Выше порога стежка; цепочечные связи по x не меняются, их порог (0.05) не достигнут.
    shift_layer2(m, 1.5 * kEps1 * l0);
    const std::vector<long long> prev0 = m.topology.prev, next0 = m.topology.next;
    const std::vector<int> obj0 = m.fiber.object_id;
    const unsigned long long ver0 = m.topology.version;
    ev = interactions::find_fiber_break_events(exec, m.fiber, m.topology, tb);
    check(ev.size() == 1 && ev[0].bond_index == 2, "рвётся ровно стежок (связь 2)");
    interactions::apply_fiber_break_events(m.fiber, m.topology, ev);
    check(m.topology.bonds.size() == 2, "стежок удалён из bonds");
    check(m.topology.version == ver0 + 1, "topology.version вырос");
    check(m.topology.prev == prev0 && m.topology.next == next0, "prev/next не изменились");
    check(m.fiber.object_id == obj0, "object_id не изменились");
    check(m.topology.bonds[0].mat == -1 && m.topology.bonds[1].mat == -1, "цепочечные связи остались с mat == -1");
    check(fiber::is_chain_bond(m.topology.bonds[0], m.topology) && fiber::is_chain_bond(m.topology.bonds[1], m.topology),
          "цепочки целы");
}

} // namespace

int main()
{
    test_parsing();
    test_parser_errors();
    test_force_and_compression();
    test_break();
    if (failures != 0)
    {
        std::printf("FAILED: %d\n", failures);
        return 1;
    }
    std::printf("OK\n");
    return 0;
}
