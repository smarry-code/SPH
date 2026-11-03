// test_k_parts.cpp
// Unit-тесты таблицы деталей (readk/k-import-plan.md, этап K3, Р2):
//   - пример из Р2 разбирается, значения совпадают;
//   - hsml=/v= в любом порядке, 1-3 компоненты v;
//   - каждая ошибка с файлом и номером строки;
//   - MATERIAL после FIBER_PART допустим;
//   - сверка с колодой: ошибки и предупреждение.
#include "k_input/k_parts.hpp"
#include "k_input/k_reader.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

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

KPartsTable parse(const std::string& t) { return read_k_parts_from_string(t, "t.parts"); }

const char* kExample =
    "# SOLID_PART <PID> <itype> <rho> [hsml=<h>] [v=<vx> <vy> <vz>]\n"
    "SOLID_PART 1 7  320  v=1546 0 0   # ударник\n"
    "SOLID_PART 2 10 2510              # керамика\n"
    "SOLID_PART 3 10 2510\n"
    "SOLID_PART 4 5  2638\n"
    "MATERIAL 0 2.0e9 0.35 1.0e-4 1140 0.02\n"
    "FIBER_PART 5 0\n";

void test_example()
{
    std::printf("[example]\n");
    KPartsTable t = parse(kExample);
    check(t.solid.size() == 4 && t.fiber.size() == 1, "4 SOLID_PART, 1 FIBER_PART");
    const SolidPart* p1 = t.find_solid(1);
    check(p1 && p1->itype == 7 && p1->rho == 320.0 && !p1->has_hsml, "PID 1: itype=7 rho=320");
    check(p1 && p1->has_v && p1->nv == 3 && p1->v[0] == 1546.0 && p1->v[1] == 0.0 && p1->v[2] == 0.0,
          "PID 1: v=1546 0 0");
    check(p1 && p1->line == 2, "PID 1: строка 2");
    const SolidPart* p4 = t.find_solid(4);
    check(p4 && p4->itype == 5 && p4->rho == 2638.0 && !p4->has_v && p4->nv == 0, "PID 4 без v");
    // itype 12-15 (Джонсон-Кук) принимаются.
    check(parse("SOLID_PART 7 12 7830\n").find_solid(7) != nullptr, "itype = 12 принят");
    check(parse("SOLID_PART 7 15 2700\n").find_solid(7) != nullptr, "itype = 15 принят");
    check(t.find_solid(5) == nullptr && t.find_fiber(5) && t.find_fiber(5)->material_id == 0,
          "PID 5 - нить с материалом 0");
    check(t.find_fiber(1) == nullptr && t.find_solid(99) == nullptr, "поиск отсутствующего PID");
    check(t.materials.size() == 1 && t.material_defined[0] && t.materials[0].E == 2.0e9
              && t.materials[0].rho == 1140.0 && t.materials[0].eps_break == 0.02,
          "MATERIAL 0 разобран");

    t = parse("SOLID_PART 1 3 100 v=5 hsml=0.01\nSOLID_PART 2 3 100 hsml=2e-3 v=1 2\n"
              "SOLID_PART 3 3 100 v= 4 5 6\n");
    check(t.solid[0].has_hsml && t.solid[0].hsml == 0.01 && t.solid[0].nv == 1 && t.solid[0].v[0] == 5.0,
          "v=5 hsml=0.01 (v первым)");
    check(t.solid[1].hsml == 2e-3 && t.solid[1].nv == 2 && t.solid[1].v[1] == 2.0, "hsml= затем v=1 2");
    check(t.solid[2].nv == 3 && t.solid[2].v[2] == 6.0, "v= с пробелом после '='");

    // MATERIAL после FIBER_PART допустим.
    t = parse("FIBER_PART 5 3\nMATERIAL 3 2.0e9 0.35 1.0e-4 1140 0.02 0.1 0.2\n");
    check(t.fiber.size() == 1 && t.materials[3].mu_ff == 0.1 && t.materials[3].mu_fs == 0.2,
          "MATERIAL после FIBER_PART допустим");
    t = parse("\n# только комментарии\n\n");
    check(t.solid.empty() && t.fiber.empty(), "пустая таблица");
}

void expect_error(const std::string& text, const std::string& sub, const std::string& line,
                  const std::string& what)
{
    const std::string e = error_of([&] { parse(text); });
    check(!e.empty() && has(e, sub) && has(e, "t.parts:" + line), what + ": " + e);
}

void test_errors()
{
    std::printf("[errors]\n");
    expect_error("# c\nBOGUS 1 2\n", "BOGUS", "2", "неизвестная директива");
    expect_error("solid_part 1 3 1\n", "solid_part", "1", "регистрозависимость");
    expect_error("\nSOLID_PART 1 0 100\n", "itype", "2", "itype = 0");
    expect_error("SOLID_PART 1 16 100\n", "itype", "1", "itype = 16");
    expect_error("SOLID_PART 1 3 0\n", "rho", "1", "rho = 0");
    expect_error("SOLID_PART 1 3 -5\n", "rho", "1", "rho < 0");
    expect_error("SOLID_PART 1 3 100 hsml=0\n", "hsml", "1", "hsml = 0");
    expect_error("SOLID_PART 1 3 100 hsml=-1\n", "hsml", "1", "hsml < 0");
    expect_error("SOLID_PART 1 3 100\nSOLID_PART 1 4 100\n", "повтор PID 1", "2", "повтор SOLID_PART");
    expect_error("MATERIAL 0 1 0.3 1 1 0.1\nSOLID_PART 5 3 1\nFIBER_PART 5 0\n", "повтор PID 5", "3",
                 "PID в SOLID_PART и FIBER_PART");
    expect_error("FIBER_PART 5 0\nSOLID_PART 5 3 1\n", "повтор PID 5", "2", "FIBER_PART, затем SOLID_PART");
    expect_error("\nFIBER_PART 5 7\nMATERIAL 0 1 0.3 1 1 0.1\n", "необъявленный MATERIAL 7", "2",
                 "необъявленный MATERIAL");
    expect_error("SOLID_PART 1 3 nan\n", "rho", "1", "rho = nan");
    expect_error("SOLID_PART 1 3 inf\n", "rho", "1", "rho = inf");
    expect_error("SOLID_PART 1 3 1e999\n", "rho", "1", "rho переполнение");
    expect_error("SOLID_PART 1 3 abc\n", "rho", "1", "rho нечисловое");
    expect_error("SOLID_PART 1 3 100 v=1 x\n", "'x'", "1", "v= нечисловое");
    expect_error("SOLID_PART 1 3 100 v=1 2 3 4\n", "'4'", "1", "v= четыре числа");
    expect_error("SOLID_PART 1 3 100 v=\n", "v=", "1", "v= без значений");
    expect_error("SOLID_PART 1 3 100 v=1 v=2\n", "повторяется", "1", "v= дважды");
    expect_error("SOLID_PART 1 3 100 hsml=1 hsml=2\n", "повторяется", "1", "hsml= дважды");
    expect_error("SOLID_PART 1 3 100 foo\n", "foo", "1", "лишний токен");
    expect_error("SOLID_PART 1 3 100 hsml=abc\n", "hsml", "1", "hsml нечисловое");
    expect_error("SOLID_PART 1 3\n", "SOLID_PART", "1", "мало полей");
    expect_error("SOLID_PART 1.5 3 100\n", "PID", "1", "PID дробный");
    expect_error("FIBER_PART 5\n", "FIBER_PART", "1", "FIBER_PART мало полей");
    expect_error("FIBER_PART 5 -1\n", "material_id", "1", "material_id < 0");
    expect_error("MATERIAL 0 1 0.3\n", "MATERIAL", "1", "MATERIAL неполный");
    expect_error("MATERIAL 0 1 0.3 1 1 0.1 -1\n", "MATERIAL", "1", "mu_ff < 0");

    const std::string e = error_of([] { read_k_parts_file("no_such_file.parts"); });
    check(has(e, "no_such_file.parts"), "нет файла: " + e);
}

KDeck make_deck()
{
    return read_k_deck_from_string(
        "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,2,0\n"
        "*PART\nA\n1\n*PART\nB\n2\n*PART\nC\n3\n"
        "*ELEMENT_SPH\n1,1,1.0\n2,1,1.0\n"
        "*ELEMENT_BEAM\n1,2,3,4\n",
        "deck.k");
}

void test_deck()
{
    std::printf("[deck]\n");
    const KDeck deck = make_deck();
    std::vector<std::string> w;

    KPartsTable t = parse("SOLID_PART 1 3 100\nMATERIAL 0 1 0.3 1 1 0.1\nFIBER_PART 2 0\n");
    check(error_of([&] { check_parts_against_deck(t, deck, true, w); }).empty() && w.empty(),
          "согласованная таблица: без ошибок и предупреждений");

    // Режим solid: балки без FIBER_PART не мешают.
    t = parse("SOLID_PART 1 3 100\n");
    w.clear();
    check(error_of([&] { check_parts_against_deck(t, deck, false, w); }).empty(), "solid: FIBER_PART не нужен");
    std::string e = error_of([&] { check_parts_against_deck(t, deck, true, w); });
    check(has(e, "PID 2") && has(e, "FIBER_PART") && has(e, "t.parts"), "solid_fiber: нет FIBER_PART: " + e);

    t = parse("MATERIAL 0 1 0.3 1 1 0.1\nFIBER_PART 2 0\n");
    e = error_of([&] { check_parts_against_deck(t, deck, false, w); });
    check(has(e, "PID 1") && has(e, "SOLID_PART") && has(e, "deck.k:"), "SPH PID без SOLID_PART: " + e);

    // Лишние PID таблицы - предупреждения; PID из *PART без элементов - известен.
    t = parse("SOLID_PART 1 3 100\nSOLID_PART 3 3 100\nSOLID_PART 99 3 100\n"
              "MATERIAL 0 1 0.3 1 1 0.1\nFIBER_PART 2 0\nFIBER_PART 98 0\n");
    w.clear();
    check(error_of([&] { check_parts_against_deck(t, deck, true, w); }).empty(), "лишние PID - не ошибка");
    check(w.size() == 2 && has(w[0], "PID 99") && has(w[0], "t.parts:3") && has(w[1], "PID 98")
              && has(w[1], "t.parts:6"),
          "предупреждения о PID 99 и 98");
}

} // namespace

int main()
{
    try
    {
        test_example();
        test_errors();
        test_deck();
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
