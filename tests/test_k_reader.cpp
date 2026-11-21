// test_k_reader.cpp
// Unit-тесты слоя чтения LS-DYNA .k (readk/k-import-plan.md, этап K2):
//   - комментарии, регистр, опции и текст в скобках в имени ключевого слова;
//   - одна и та же модель в free, fixed, long (LONG=Y и суффикс '+') даёт
//     одинаковый KDeck;
//   - пустые поля = умолчания, D-экспонента, мусор в числе = ошибка;
//   - *INCLUDE (относительный путь, *INCLUDE_PATH, цикл, нет файла), *END;
//   - &параметр, *INCLUDE_TRANSFORM, *ELEMENT_BEAM_THICKNESS -> ошибка;
//   - номер строки в ошибках; висячие ссылки и дубликаты; узел в SPH и балке;
//   - сводка пропущенных ключевых слов; *PART в стиле LS-PrePost.
//
// Тест НЕ зависит от корневого 1.k (производительность - этап K4).
#include "k_input/k_cards.hpp"
#include "k_input/k_reader.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef SPH_TEST_DATA_DIR
#error "SPH_TEST_DATA_DIR должен быть задан CMake"
#endif

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

const std::string kDir = std::string(SPH_TEST_DATA_DIR) + "/k/";

using namespace k_input;

// Текст исключения или пустая строка, если не брошено.
template <class F>
std::string error_of(F&& f)
{
    try { f(); }
    catch (const std::exception& e) { return e.what(); }
    return std::string();
}

bool has(const std::string& s, const std::string& sub)
{
    return s.find(sub) != std::string::npos;
}

// Подпись KDeck без источников: сравнение моделей из разных форматов.
std::string signature(const KDeck& d)
{
    std::ostringstream o;
    char b[512];
    o << "title=" << d.title << "\n";
    for (const KNode& n : d.nodes)
    {
        std::snprintf(b, sizeof b, "N %lld %.17g %.17g %.17g\n", n.nid, n.x, n.y, n.z);
        o << b;
    }
    for (const KSph& s : d.sph)
    {
        std::snprintf(b, sizeof b, "S %lld %lld %.17g %d\n", s.nid, s.pid, s.mass, s.is_volume ? 1 : 0);
        o << b;
    }
    for (const KBeam& e : d.beams)
        o << "B " << e.eid << " " << e.pid << " " << e.n1 << " " << e.n2 << " " << e.n3 << "\n";
    for (const KSolid& e : d.solids)
    {
        o << "E " << e.eid << " " << e.pid << " " << e.nn;
        for (int k = 0; k < e.nn; ++k) o << " " << e.n[static_cast<std::size_t>(k)];
        o << "\n";
    }
    for (const KPart& p : d.parts) o << "P " << p.pid << " [" << p.title << "]\n";
    for (const KSpcNode& c : d.spc_nodes)
        o << "C " << c.nid << " " << c.dofx << c.dofy << c.dofz << "\n";
    for (const KSpcSet& c : d.spc_sets)
        o << "CS " << c.nsid << " " << c.dofx << c.dofy << c.dofz << "\n";
    for (const auto& kv : d.node_sets)
    {
        o << "SET " << kv.first << ":";
        for (long long n : kv.second.nids) o << " " << n;
        o << "\n";
    }
    for (const auto& kv : d.skipped) o << "SKIP " << kv.first << " " << kv.second << "\n";
    return o.str();
}

void write_file(const std::string& name, const std::string& content)
{
    std::ofstream out(name, std::ios::binary);
    out << content;
}

// Минимальная корректная заготовка для проверок ссылок.
const char* kBase =
    "*PART\nP\n1\n*PART\nQ\n2\n"
    "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n";

KDeck parse(const std::string& text, bool validate = true)
{
    return read_k_deck_from_string(text, "mem.k", validate);
}

// --- 1. Эталонная модель в разных форматах ---------------------------------
void test_formats()
{
    std::printf("[formats]\n");
    KDeck df = read_k_deck(kDir + "model_free.k");
    check(df.title == "Small test model", "TITLE прочитан");
    check(df.nodes.size() == 16, "free: 16 узлов");
    check(df.sph.size() == 4, "free: NEND 1..3 развёрнут + 1 частица = 4");
    check(df.sph[0].nid == 1 && df.sph[2].nid == 3 && df.sph[2].mass == 0.5 && df.sph[2].pid == 1,
          "free: NEND включительно, тот же PID и MASS");
    check(df.sph[3].mass == -0.125, "free: отрицательная MASS сохранена как есть");
    check(df.beams.size() == 2 && df.beams[0].n3 == 13 && df.beams[1].n3 == 0,
          "free: балки, пустой N3 = 0");
    check(df.solids.size() == 1 && df.solids[0].nn == 8, "free: солид 8 узлов");
    check(df.parts.size() == 2 && df.parts[0].title == "Impactor" && df.parts[1].pid == 2,
          "free: *PART заголовок + PID (стиль LS-PrePost)");
    check(df.spc_nodes.size() == 1 && df.spc_nodes[0].dofx == 1 && df.spc_nodes[0].dofz == 1,
          "free: BOUNDARY_SPC_NODE");
    check(df.spc_sets.size() == 1 && df.spc_sets[0].nsid == 5 && df.spc_sets[0].dofx == 1 &&
          df.spc_sets[0].dofy == 0, "free: BOUNDARY_SPC_SET, пустые DOF = 0");
    check(df.node_sets.at(5).nids.size() == 10, "free: SET_NODE_LIST 10 узлов");
    check(df.node_sets.at(6).nids.size() == 7 && df.node_sets.at(6).nids.front() == 20 &&
          df.node_sets.at(6).nids.back() == 27, "free: SET_NODE_LIST_GENERATE включительно");
    check(df.skipped.count("MAT_ELASTIC") && df.skipped.at("MAT_ELASTIC") == 1, "сводка пропущенных: MAT_ELASTIC=1");
    check(df.find_node(99) == nullptr, "строки после *END не читаются");
    const KNode* n12 = df.find_node(12);
    check(n12 && n12->x == 1.0 && n12->z == -2.0, "D-экспонента 1.0D+0 и -2.0");
    const KNode* n11 = df.find_node(11);
    check(n11 && n11->y == 0.0, "пустое поле = умолчание 0.0");

    const std::string sf = signature(df);
    check(signature(read_k_deck(kDir + "model_fixed.k")) == sf, "fixed даёт тот же KDeck, что free");
    check(signature(read_k_deck(kDir + "model_long.k")) == sf, "long (LONG=Y) даёт тот же KDeck");
    check(signature(read_k_deck(kDir + "model_plus.k")) == sf, "long через суффикс '+' даёт тот же KDeck");
    check(signature(read_k_deck_raw(kDir + "model_fixed.k")) == sf, "read_k_deck_raw то же без проверки");

    // фикс. формат: ширины берутся из таблицы
    check(k_input::field_offset(k_input::kNodeCard, 1, false) == 8 &&
          k_input::field_offset(k_input::kNodeCard, 4, false) == 56 &&
          k_input::field_offset(k_input::kNodeCard, 1, true) == 20 &&
          k_input::field_offset(k_input::kNodeCard, 4, true) == 80,
          "k_cards: смещения standard (8,56) и long (20,80)");
}

// --- 2. Синтаксис ключевых слов и полей -------------------------------------
void test_syntax()
{
    std::printf("[syntax]\n");
    KDeck d = parse(std::string(kBase) +
        "$ комментарий\n*element_sph\r\n1,1,0.5\r\n$ ещё\r\n2,1,0.5\r\n"
        "*ELEMENT_SPH_VOLUME\n3,2,-2.0\n");
    check(d.sph.size() == 3, "CRLF, комментарии, регистр: 3 частицы");
    check(!d.sph[0].is_volume && d.sph[2].is_volume && d.sph[2].mass == -2.0, "_VOLUME -> is_volume");

    d = parse(std::string(kBase) + "*ELEMENT_SOLID (ten nodes format)\n"
        "7,1\n1,2,3,1,2,3,1,2,3,1\n");
    check(d.solids.size() == 1 && d.solids[0].nn == 10 && d.solids[0].eid == 7 && d.solids[0].n[9] == 1,
          "(ten nodes format): две карточки, 10 узлов");
    d = parse(std::string(kBase) + "*Element_Solid  (ten nodes format)\n"
        "7       1\n       1       2       3       1       2       3       1       2       3       1\n");
    check(d.solids.size() == 1 && d.solids[0].nn == 10, "(ten nodes format) в fixed");
    d = parse(std::string(kBase) + "*ELEMENT_SOLID\n7,1,1,2,3,1\n");
    check(d.solids[0].nn == 4, "солид с 4 узлами");

    d = parse(std::string(kBase) + "*KEYWORD MEMORY=100m NCPU=2\n*NODE+\n"
        "                   4                 1.0                 2.0                 3.0\n", false);
    check(d.find_node(4) && d.find_node(4)->z == 3.0, "*KEYWORD с прочими параметрами и *NODE+");
    d = parse("*KEYWORD LONG=Y\n*NODE-\n4,1,2,3\n       5             0.5\n", false);
    check(d.nodes.size() == 2 && d.nodes[1].nid == 5 && d.nodes[1].x == 0.5, "LONG=Y + '-' -> standard-ширины блока");
    d = parse("*KEYWORD LONG=S\n*NODE\n4,1,2,3\n", false);
    check(d.nodes.size() == 1, "LONG=S");

    d = parse("*NODE\n1,1e3,1.5E+2,1.0d-3\n2, ,5 ,\n", false);
    check(d.nodes[0].x == 1000.0 && d.nodes[0].y == 150.0 && d.nodes[0].z == 1.0e-3,
          "вещественные: 1e3, 1.5E+2, 1.0d-3");
    check(d.nodes[1].x == 0.0 && d.nodes[1].y == 5.0 && d.nodes[1].z == 0.0,
          "free: пустые поля и недостающий хвост = умолчания");

    d = parse("*TITLE\nМодель 1\nвторая строка\n*NODE\n1,0,0,0\n", false);
    check(d.title == "Модель 1", "*TITLE: следующая строка");

    d = parse(std::string(kBase) + "*PART_CONTACT\nx\n*PARTICLE_BLAST\n1\n*MAT_ELASTIC\n*MAT_ELASTIC\n", false);
    check(d.parts.size() == 2 && d.skipped.at("PART_CONTACT") == 1 && d.skipped.at("PARTICLE_BLAST") == 1 &&
          d.skipped.at("MAT_ELASTIC") == 2, "*PART_xxx/*PARTICLE_* пропускаются, не путаются с *PART");
    d = parse("*PART\nname\n5,1,2,3\nextra card line\n  1.0 2.0\n*NODE\n1,0,0,0\n", false);
    check(d.parts.size() == 1 && d.parts[0].pid == 5 && d.nodes.size() == 1,
          "*PART: лишние карточки до следующего * пропущены");
    d = parse("*PART\n\n6\n", false);
    check(d.parts.size() == 1 && d.parts[0].pid == 6 && d.parts[0].title.empty(), "*PART с пустым заголовком");
}

// --- 3. Ошибки ---------------------------------------------------------------
void test_errors()
{
    std::printf("[errors]\n");
    std::string e = error_of([] { parse("*NODE\n1,0,0,0\n2,0,&px,0\n", false); });
    check(has(e, "mem.k:3") && has(e, "PARAMETER"), "&параметр -> ошибка с номером строки 3: " + e);
    e = error_of([] { parse("*NODE\n1,0,0,0\n2,abc,0,0\n", false); });
    check(has(e, "mem.k:3"), "мусор в числовом поле -> ошибка: " + e);
    e = error_of([] { parse("*NODE\n1,0,1.5x,0\n", false); });
    check(!e.empty(), "хвост '1.5x' -> ошибка");
    e = error_of([] { parse("*NODE\n1.5,0,0,0\n", false); });
    check(!e.empty(), "нецелый NID -> ошибка");
    e = error_of([] { parse("*NODE\n,0,0,0\n", false); });
    check(has(e, "mem.k:2") && has(e, "NID"), "обязательное поле пусто -> ошибка: " + e);
    e = error_of([] { parse("*NODE\n       1\n", false); });
    check(e.empty(), "fixed: пустые X/Y/Z допустимы");
    e = error_of([] { parse("*INCLUDE_TRANSFORM\nfile.k\n", false); });
    check(has(e, "mem.k:1") && has(e, "INCLUDE_TRANSFORM"), "*INCLUDE_TRANSFORM -> ошибка: " + e);
    e = error_of([] { parse("*NODE\n1,0,0,0\n*ELEMENT_BEAM_THICKNESS\n1,1,1,2\n", false); });
    check(has(e, "mem.k:3") && has(e, "ELEMENT_BEAM_THICKNESS"), "*ELEMENT_BEAM_THICKNESS -> ошибка: " + e);
    e = error_of([] { parse("*ELEMENT_SOLID_ORTHO\n", false); });
    check(!e.empty(), "*ELEMENT_SOLID_ORTHO -> ошибка");
    e = error_of([] { parse("*ELEMENT_SPH_GENERAL\n", false); });
    check(!e.empty(), "*ELEMENT_SPH_xxx кроме _VOLUME -> ошибка");
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_SPH\n1,1,0.5,0\n2,1,0.5,1\n", false); });
    check(has(e, "NEND") && has(e, "mem.k:"), "NEND < NID -> ошибка: " + e);
    e = error_of([] { parse("*PART\nname\n", false); });
    check(!e.empty(), "*PART без карточки PID -> ошибка");
    e = error_of([] { parse(std::string(kBase) + "*SET_NODE_LIST_GENERATE\n1\n5,2\n", false); });
    check(!e.empty(), "диапазон набора с концом < начала -> ошибка");
    e = error_of([] { read_k_deck(kDir + "no_such.k"); });
    check(!e.empty(), "несуществующий главный файл -> ошибка");

    // ссылки и дубликаты (Р8)
    e = error_of([] { parse(std::string(kBase) + "*NODE\n2,5,5,5\n"); });
    check(has(e, "mem.k:") && has(e, "дубликат NID 2"), "дубликат NID в *NODE: " + e);
    e = error_of([] { parse(std::string(kBase) + "*PART\nR\n2\n"); });
    check(has(e, "дубликат PID 2"), "дубликат PID: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_SPH\n1,1,1.0\n1,1,1.0\n"); });
    check(has(e, "дубликат NID 1"), "дубликат NID в *ELEMENT_SPH: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_BEAM\n1,2,1,2\n1,2,2,3\n"); });
    check(has(e, "дубликат EID 1"), "дубликат EID балок: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_SOLID\n1,1,1,2,3,1\n1,1,1,2,3,1\n"); });
    check(has(e, "дубликат EID 1"), "дубликат EID солидов: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_SPH\n9,1,1.0\n"); });
    check(has(e, "узел 9"), "SPH ссылается на несуществующий узел: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_BEAM\n1,2,1,9\n"); });
    check(has(e, "узел 9"), "балка: несуществующий узел: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_SOLID\n1,1,1,2,3,9\n"); });
    check(has(e, "узел 9"), "солид: несуществующий узел: " + e);
    e = error_of([] { parse(std::string(kBase) + "*BOUNDARY_SPC_NODE\n9,,1,1,1\n"); });
    check(has(e, "узел 9"), "SPC_NODE: несуществующий узел: " + e);
    e = error_of([] { parse(std::string(kBase) + "*BOUNDARY_SPC_SET\n7,,1,1,1\n"); });
    check(has(e, "набор узлов 7"), "SPC_SET: неизвестный SID: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_SPH\n1,3,1.0\n"); });
    check(has(e, "PID 3"), "SPH: PID без *PART: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_BEAM\n1,3,1,2\n"); });
    check(has(e, "PID 3"), "балка: PID без *PART: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_SPH\n2,1,1.0\n*ELEMENT_BEAM\n1,2,1,2\n"); });
    check(has(e, "одновременно") && has(e, "mem.k:"), "узел и в SPH, и в балке: " + e);
    e = error_of([] { parse(std::string(kBase) + "*ELEMENT_SPH\n1,1,1.0\n*ELEMENT_BEAM\n1,2,2,3\n*BOUNDARY_SPC_NODE\n3,,1\n"); });
    check(e.empty(), "корректная модель проходит проверку: " + e);
    // порядок блоков произволен: элементы до *NODE
    e = error_of([] { parse("*ELEMENT_SPH\n1,1,1.0\n*NODE\n1,0,0,0\n*PART\nP\n1\n"); });
    check(e.empty(), "ссылки разрешаются после чтения всего файла: " + e);
}

// --- 4. Include --------------------------------------------------------------
void test_include()
{
    std::printf("[include]\n");
    KDeck d = read_k_deck(kDir + "inc_main.k");
    check(d.nodes.size() == 2 && d.sph.size() == 2, "include: относительный путь + *INCLUDE_PATH, *END только в include");
    check(d.find_node(99) == nullptr, "*END в include-файле прекращает чтение только его");
    check(d.files.size() == 3, "include: источники файлов учтены");
    check(d.where(d.nodes[1].src).find("inc_in_path.k:2") != std::string::npos, "источник узла из include: файл и строка");

    std::string e = error_of([] { read_k_deck(kDir + "cycle_a.k"); });
    check(has(e, "цикл") && has(e, "cycle_"), "цикл include -> ошибка: " + e);
    e = error_of([] { read_k_deck(kDir + "missing_inc.k"); });
    check(has(e, "missing_inc.k:2") && has(e, "no_such_file.k"), "нет include-файла -> ошибка с файлом и строкой: " + e);

    // продолжение ' +' и *INCLUDE_PATH_RELATIVE: временные файлы в cwd теста
    write_file("__k_inc_a.k", "*NODE\n1,0,0,0\n");
    write_file("__k_inc_main.k",
               "*INCLUDE\n__k_inc_ +\na.k\n");
    d = read_k_deck_raw("__k_inc_main.k");
    check(d.nodes.size() == 1, "*INCLUDE: имя файла с продолжением ' +'");
    write_file("__k_inc_main2.k", "*INCLUDE_PATH_RELATIVE\n.\n*INCLUDE\n__k_inc_a.k\n");
    d = read_k_deck_raw("__k_inc_main2.k");
    check(d.nodes.size() == 1, "*INCLUDE_PATH_RELATIVE");
    write_file("__k_bad.k", "*NODE\n1,0,0,0\n2,0,zzz,0\n");
    write_file("__k_inc_main3.k", "*INCLUDE\n__k_bad.k\n");
    e = error_of([] { read_k_deck_raw("__k_inc_main3.k"); });
    check(has(e, "__k_bad.k:3"), "ошибка в include-файле называет include-файл и его строку: " + e);
    write_file("__k_inc_main4.k", "*INCLUDE_TRANSFORM\n__k_inc_a.k\n");
    e = error_of([] { read_k_deck_raw("__k_inc_main4.k"); });
    check(has(e, "__k_inc_main4.k:1"), "*INCLUDE_TRANSFORM из файла: " + e);
    for (const char* f : {"__k_inc_a.k", "__k_inc_main.k", "__k_inc_main2.k", "__k_bad.k",
                          "__k_inc_main3.k", "__k_inc_main4.k"})
        std::remove(f);
}

} // namespace

int main()
{
    test_formats();
    test_syntax();
    test_errors();
    test_include();
    if (failures)
    {
        std::printf("\nFAILED: %d\n", failures);
        return 1;
    }
    std::printf("\nвсе проверки пройдены\n");
    return 0;
}
