// test_fiber_input_inelastic.cpp
// Т15 (notes/fiber-model/fiber-plasticity-design.md, раздел 8): директивы неупругости связи
// в .fib -- PLASTIC, PLASTIC_TABLE, RATE, VISCO, FAILURE.
//   - все директивы корректно попадают в поля FiberMaterial;
//   - порядок директив/ключей свободен, директивы до и после MATERIAL/FIBER;
//   - файл без директив даёт упругий материал (inelastic() == false);
//   - каждая ошибка -- std::runtime_error с понятным текстом (проверяется
//     подстрока и номер строки);
//   - вывод fiber_output: колонки eps_p sigma w_p w_v только при наличии
//     неупругих материалов, иначе формат прежний (6 колонок).
// Собирается с SPH_DIM=3 (нити, notes/fiber-model/fiber-problem-statement.md).
#include "fiber/fiber_input.hpp"
#include "fiber/fiber_output.hpp"
#include "fiber/model.hpp"
#include "param.hpp"

#include "e2e_harness.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_input_inelastic требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

const char* kPath = "__fiber_inelastic_fixture.fib";

const char* kBase =
    "MATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.30\n"
    "FIBER 0 3\n0 0 0\n0.01 0 0\n0.02 0 0\n";

void write_file(const std::string& content)
{
    std::ofstream out(kPath);
    out << content;
}

fiber::Model load(const std::string& content)
{
    write_file(content);
    fiber::Model m;
    fiber::load_fiber_file(kPath, m);
    return m;
}

// Ожидаем runtime_error, в тексте которого есть все подстроки (через '|').
void expect_error(const std::string& name, const std::string& content, const std::string& needles)
{
    try
    {
        load(content);
        check(false, name + ": ожидалось исключение, его нет");
    }
    catch (const std::runtime_error& e)
    {
        const std::string msg = e.what();
        bool all = true;
        std::size_t pos = 0;
        while (pos <= needles.size())
        {
            const std::size_t bar = needles.find('|', pos);
            const std::string n = needles.substr(pos, bar == std::string::npos ? std::string::npos : bar - pos);
            if (msg.find(n) == std::string::npos) all = false;
            if (bar == std::string::npos) break;
            pos = bar + 1;
        }
        check(all, name + ": '" + msg + "'");
    }
}

bool near(double a, double b) { return std::fabs(a - b) <= 1e-12 * std::max(1.0, std::fabs(b)); }

void test_no_directives_is_elastic()
{
    std::printf("test_no_directives_is_elastic:\n");
    const fiber::Model m = load(kBase);
    const fiber::FiberMaterial& mat = m.materials[0];
    check(!mat.inelastic(), "без директив материал упругий");
    check(!mat.plastic.enabled() && mat.rate.model == fiber::RateModel::None, "plastic/rate выключены");
    check(mat.visco.eta == 0 && mat.visco.xi == 0 && mat.visco.E1 == 0 && mat.visco.tau == 0, "visco нулевое");
    check(mat.failure.eps_p == 0 && mat.failure.sigma_u == 0, "failure выключено");
}

void test_all_directives()
{
    std::printf("test_all_directives:\n");
    // Директивы ДО MATERIAL/FIBER, ключи в произвольном порядке/регистре.
    const fiber::Model m = load(std::string(
                                    "PLASTIC 0 H=0.4e9 SY=60e6   # линейное\n"
                                    "RATE 0 CS P=5 C=40\n"
                                    "VISCO 0 E1=0.5e9 TAU=1e-5 ETA=3.0\n"
                                    "FAILURE 0 SIGMA_U=90e6 EPS_P=0.25\n") +
                                kBase);
    const fiber::FiberMaterial& mat = m.materials[0];
    check(mat.plastic.type == fiber::HardeningType::Linear && near(mat.plastic.sy, 60e6) && near(mat.plastic.H, 0.4e9),
          "PLASTIC линейное: SY, H");
    check(mat.rate.model == fiber::RateModel::CowperSymonds && near(mat.rate.a, 40) && near(mat.rate.b, 5),
          "RATE CS: C=40, P=5");
    check(near(mat.visco.eta, 3.0) && near(mat.visco.E1, 0.5e9) && near(mat.visco.tau, 1e-5) && mat.visco.xi == 0,
          "VISCO: ETA, E1, TAU");
    check(near(mat.failure.eps_p, 0.25) && near(mat.failure.sigma_u, 90e6), "FAILURE: EPS_P, SIGMA_U");
    check(mat.inelastic() && mat.plastic.enabled(), "материал неупругий");
    check(near(mat.E, 2.0e9) && near(mat.eps_break, 0.30), "поля MATERIAL не затронуты");

    // Идеальная пластичность (H не задан) и степенное упрочнение, директивы ПОСЛЕ FIBER.
    const fiber::Model m2 = load(std::string(kBase) + "PLASTIC 0 SY=50e6\n");
    check(m2.materials[0].plastic.type == fiber::HardeningType::Linear && m2.materials[0].plastic.H == 0.0,
          "PLASTIC без H: идеальная пластичность");
    const fiber::Model m3 = load(std::string(kBase) + "PLASTIC 0 SY=60e6 K=300e6 N=0.4\nRATE 0 JC C=0.05 EDOT0=1e-3\n");
    const fiber::FiberMaterial& p3 = m3.materials[0];
    check(p3.plastic.type == fiber::HardeningType::Power && near(p3.plastic.K, 300e6) && near(p3.plastic.n, 0.4),
          "PLASTIC степенное: K, N");
    check(p3.rate.model == fiber::RateModel::JohnsonCook && near(p3.rate.a, 0.05) && near(p3.rate.b, 1e-3),
          "RATE JC: C, EDOT0");
    const fiber::Model m4 = load(std::string(kBase) + "PLASTIC 0 SY=1e6\nRATE 0 PERZYNA ETA=1e6 M=1.0\nVISCO 0 XI=0.05\n");
    const fiber::FiberMaterial& p4 = m4.materials[0];
    check(p4.rate.model == fiber::RateModel::Perzyna && near(p4.rate.a, 1e6) && near(p4.rate.b, 1.0),
          "RATE PERZYNA: ETA, M");
    check(near(p4.visco.xi, 0.05) && p4.visco.eta == 0.0, "VISCO XI");
}

void test_table_and_second_material()
{
    std::printf("test_table_and_second_material:\n");
    const fiber::Model m = load(
        "MATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.30\n"
        "MATERIAL 1 1.0e9 0.35 1.0e-4 1140.0 0.30\n"
        "PLASTIC_TABLE 1 4   # таблица\n"
        "0.000  60e6\n"
        "0.020  72e6   # комментарий\n"
        "\n"
        "0.100  95e6\n"
        "0.300  95e6\n"
        "RATE 1 CS C=40 P=5\n"
        "FIBER 0 2\n0 0 0\n0.01 0 0\n"
        "FIBER 1 2\n0 0.01 0\n0.01 0.01 0\n");
    check(!m.materials[0].inelastic(), "материал 0 остался упругим");
    const fiber::PlasticParams& p = m.materials[1].plastic;
    check(p.type == fiber::HardeningType::Table && p.table_eps.size() == 4 && p.table_sigma.size() == 4 && p.enabled(),
          "таблица: 4 точки");
    check(near(p.table_eps[1], 0.02) && near(p.table_sigma[1], 72e6) && near(p.table_eps[3], 0.3) &&
              near(p.table_sigma[3], 95e6),
          "значения таблицы");
    check(m.materials[1].rate.model == fiber::RateModel::CowperSymonds, "RATE к таблице");
}

void test_errors()
{
    std::printf("test_errors:\n");
    const std::string B = kBase;
    // PLASTIC
    expect_error("неизвестный ключ", B + "PLASTIC 0 SY=60e6 FOO=1\n", "PLASTIC|неизвестный ключ 'FOO'|строке 6");
    expect_error("токен без '='", B + "PLASTIC 0 SY 60e6\n", "ожидалось КЛЮЧ=значение");
    expect_error("значение не число", B + "PLASTIC 0 SY=abc\n", "не конечное число");
    expect_error("нет SY", B + "PLASTIC 0 H=1e9\n", "не задан обязательный ключ SY");
    expect_error("SY<=0", B + "PLASTIC 0 SY=0\n", "SY должен быть > 0");
    expect_error("H<0", B + "PLASTIC 0 SY=1e6 H=-1\n", "H должен быть >= 0");
    expect_error("K<0", B + "PLASTIC 0 SY=1e6 K=-1 N=0.5\n", "K должен быть >= 0");
    expect_error("N>1", B + "PLASTIC 0 SY=1e6 K=1e6 N=1.5\n", "N должен быть в (0, 1]");
    expect_error("N=0", B + "PLASTIC 0 SY=1e6 K=1e6 N=0\n", "N должен быть в (0, 1]");
    expect_error("K без N", B + "PLASTIC 0 SY=1e6 K=1e6\n", "оба ключа K и N");
    expect_error("H и K", B + "PLASTIC 0 SY=1e6 H=1e6 K=1e6 N=0.5\n", "несовместимы");
    expect_error("повтор ключа", B + "PLASTIC 0 SY=1e6 SY=2e6\n", "повторяется");
    expect_error("повтор PLASTIC", B + "PLASTIC 0 SY=1e6\nPLASTIC 0 SY=2e6\n", "повторная директива|первая в строке 6");
    expect_error("нет id", B + "PLASTIC SY=1e6\n", "ожидался id материала");
    expect_error("неопределённый материал", B + "PLASTIC 3 SY=1e6\n", "неопределённый материал 3");
    // PLASTIC_TABLE
    expect_error("PLASTIC после таблицы", B + "PLASTIC_TABLE 0 2\n0 1e6\n0.1 2e6\nPLASTIC 0 SY=1e6\n", "PLASTIC и PLASTIC_TABLE");
    expect_error("таблица после PLASTIC", B + "PLASTIC 0 SY=1e6\nPLASTIC_TABLE 0 2\n0 1e6\n0.1 2e6\n", "PLASTIC и PLASTIC_TABLE");
    expect_error("повтор таблицы", B + "PLASTIC_TABLE 0 2\n0 1e6\n0.1 2e6\nPLASTIC_TABLE 0 2\n0 1e6\n0.1 2e6\n", "повторная директива");
    expect_error("N<2", B + "PLASTIC_TABLE 0 1\n0 1e6\n", "число строк таблицы должно быть >= 2");
    expect_error("N не число", B + "PLASTIC_TABLE 0 abc\n", "ожидалось число строк");
    expect_error("первая eps_p != 0", B + "PLASTIC_TABLE 0 2\n0.01 1e6\n0.1 2e6\n", "первая точка таблицы|eps_p = 0");
    expect_error("eps_p не возрастает", B + "PLASTIC_TABLE 0 3\n0 1e6\n0.1 2e6\n0.1 3e6\n", "строго возрастать");
    expect_error("sigma убывает", B + "PLASTIC_TABLE 0 3\n0 3e6\n0.1 2e6\n0.2 4e6\n", "не должен убывать");
    expect_error("sigma<=0", B + "PLASTIC_TABLE 0 2\n0 0\n0.1 2e6\n", "sigma_s должен быть > 0");
    expect_error("мало строк", B + "PLASTIC_TABLE 0 3\n0 1e6\n0.1 2e6\n", "конец файла");
    expect_error("мало строк (директива вместо строки)", B + "PLASTIC_TABLE 0 3\n0 1e6\n0.1 2e6\nPLASTIC 0 SY=1\n",
                 "строка таблицы 3 из 3");
    expect_error("лишний токен в строке таблицы", B + "PLASTIC_TABLE 0 2\n0 1e6 7\n0.1 2e6\n", "лишний токен");
    // RATE
    expect_error("RATE без PLASTIC", B + "RATE 0 CS C=40 P=5\n", "RATE|требует PLASTIC");
    expect_error("RATE неизвестная модель", B + "PLASTIC 0 SY=1e6\nRATE 0 XYZ C=1\n", "неизвестная модель 'XYZ'");
    expect_error("RATE без модели", B + "PLASTIC 0 SY=1e6\nRATE 0\n", "ожидалась модель");
    expect_error("RATE CS ключ JC", B + "PLASTIC 0 SY=1e6\nRATE 0 CS C=1 EDOT0=2\n", "неизвестный ключ 'EDOT0'");
    expect_error("RATE CS без P", B + "PLASTIC 0 SY=1e6\nRATE 0 CS C=1\n", "не задан обязательный ключ P");
    expect_error("RATE CS C<=0", B + "PLASTIC 0 SY=1e6\nRATE 0 CS C=0 P=5\n", "C должен быть > 0");
    expect_error("RATE CS P<1", B + "PLASTIC 0 SY=1e6\nRATE 0 CS C=1 P=0.5\n", "P должен быть >= 1");
    expect_error("RATE JC C<=0", B + "PLASTIC 0 SY=1e6\nRATE 0 JC C=0 EDOT0=1\n", "C должен быть > 0");
    expect_error("RATE JC EDOT0<=0", B + "PLASTIC 0 SY=1e6\nRATE 0 JC C=1 EDOT0=0\n", "EDOT0 должен быть > 0");
    expect_error("RATE PERZYNA ETA<0", B + "PLASTIC 0 SY=1e6\nRATE 0 PERZYNA ETA=-1 M=1\n", "ETA должен быть > 0");
    expect_error("RATE PERZYNA ETA=0", B + "PLASTIC 0 SY=1e6\nRATE 0 PERZYNA ETA=0 M=1\n", "ETA должен быть > 0");
    expect_error("RATE PERZYNA M<=0", B + "PLASTIC 0 SY=1e6\nRATE 0 PERZYNA ETA=1 M=0\n", "M должен быть > 0");
    expect_error("повтор RATE", B + "PLASTIC 0 SY=1e6\nRATE 0 CS C=1 P=2\nRATE 0 JC C=1 EDOT0=1\n", "повторная директива");
    // VISCO
    expect_error("VISCO пусто", B + "VISCO 0\n", "не задан ни один ключ");
    expect_error("VISCO неизвестный ключ", B + "VISCO 0 FOO=1\n", "неизвестный ключ 'FOO'");
    expect_error("VISCO ETA<0", B + "VISCO 0 ETA=-1\n", "ETA должен быть >= 0");
    expect_error("VISCO XI<0", B + "VISCO 0 XI=-1\n", "XI должен быть >= 0");
    expect_error("VISCO E1<0", B + "VISCO 0 E1=-1 TAU=1\n", "E1 должен быть >= 0");
    expect_error("VISCO TAU<=0", B + "VISCO 0 E1=1e6 TAU=0\n", "TAU должен быть > 0");
    expect_error("VISCO E1 без TAU", B + "VISCO 0 E1=1e6\n", "E1 > 0 требует TAU");
    expect_error("VISCO TAU без E1", B + "VISCO 0 TAU=1e-5\n", "TAU задан без E1");
    expect_error("VISCO ETA и XI", B + "VISCO 0 ETA=1 XI=0.1\n", "ETA и XI одновременно");
    expect_error("повтор VISCO", B + "VISCO 0 ETA=1\nVISCO 0 XI=0.1\n", "повторная директива");
    // FAILURE
    expect_error("FAILURE пусто", B + "FAILURE 0\n", "не задан ни один ключ");
    expect_error("FAILURE EPS_P<=0", B + "FAILURE 0 EPS_P=0\n", "EPS_P должен быть > 0");
    expect_error("FAILURE SIGMA_U<=0", B + "FAILURE 0 SIGMA_U=-5\n", "SIGMA_U должен быть > 0");
    expect_error("FAILURE ключ", B + "FAILURE 0 SY=5\n", "неизвестный ключ 'SY'");
    expect_error("повтор FAILURE", B + "FAILURE 0 EPS_P=0.1\nFAILURE 0 EPS_P=0.2\n", "повторная директива");
    // Неизвестная директива ведёт себя как раньше
    expect_error("неизвестная директива", B + "PLASTICITY 0 SY=1\n", "неизвестное ключевое слово 'PLASTICITY'");
}

std::vector<std::string> read_lines(const std::string& path)
{
    std::ifstream in(path);
    std::vector<std::string> v;
    std::string l;
    while (std::getline(in, l)) v.push_back(l);
    return v;
}

std::vector<std::string> tokens(const std::string& l)
{
    std::istringstream s(l);
    std::vector<std::string> t;
    std::string x;
    while (s >> x) t.push_back(x);
    return t;
}

void test_output_columns()
{
    std::printf("test_output_columns:\n");
    // Упругий: прежний формат (6 колонок).
    fiber::Model me = load(kBase);
    fiber::fiber_output(me, 1, 0.0);
    {
        const auto lines = read_lines("fiber_output_00001.txt");
        std::size_t b = 0;
        while (b < lines.size() && lines[b] != "BONDS") ++b;
        check(b > 0 && lines[b - 1] == "# i j l0 l l_break elongation", "упругий: заголовок BONDS прежний");
        check(b + 1 < lines.size() && tokens(lines[b + 1]).size() == 6, "упругий: 6 колонок в строке связи");
    }
    // Неупругий: +eps_p sigma w_p w_v; растянем вторую связь и проверим sigma/eps_p.
    fiber::Model mi = load(std::string(kBase) + "PLASTIC 0 SY=60e6 H=0.4e9\n");
    const double E = mi.materials[0].E;
    const double A = pi * mi.materials[0].d * mi.materials[0].d / 4.0;
    mi.fiber.x[0][2] = 0.0202; // связь 1 растянута на 1% от l0 (x1 = 0.01, x2 = 0.0202)
    const double l1 = mi.fiber.x[0][2] - mi.fiber.x[0][1];
    const double eps1 = l1 / mi.topology.bonds[1].l0 - 1.0;
    mi.topology.bonds[1].st.eps_p = 0.002;
    mi.topology.bonds[1].st.w_p = 5.0;
    mi.topology.bonds[1].st.w_v = 2.0;
    fiber::fiber_output(mi, 2, 0.0);
    const auto lines = read_lines("fiber_output_00002.txt");
    std::size_t b = 0;
    while (b < lines.size() && lines[b] != "BONDS") ++b;
    check(b > 0 && lines[b - 1] == "# i j l0 l l_break elongation eps_p sigma w_p w_v", "неупругий: заголовок BONDS расширен");
    const auto t0 = tokens(lines[b + 1]);
    const auto t1 = tokens(lines[b + 2]);
    check(t0.size() == 10 && t1.size() == 10, "неупругий: 10 колонок");
    if (t1.size() == 10)
    {
        check(near(std::stod(t1[6]), 0.002), "eps_p = состояние связи");
        check(std::fabs(std::stod(t1[7]) - E * (eps1 - 0.002)) < 1e-6 * E * eps1, "sigma = E (ε - ε_p) (без вязкости, v = 0)");
        check(near(std::stod(t1[8]), 5.0 * A * mi.topology.bonds[1].l0), "w_p в Дж = w_p·A·l0");
        check(near(std::stod(t1[9]), 2.0 * A * mi.topology.bonds[1].l0), "w_v в Дж = w_v·A·l0");
    }
    if (t0.size() == 10)
        check(std::stod(t0[7]) == 0.0 && std::stod(t0[6]) == 0.0, "ненагруженная связь: sigma = eps_p = 0");
}

// --- Доработки парсера (ревью этапа 3) ---------------------------------------
std::string to_crlf(const std::string& s)
{
    std::string r;
    for (char ch : s)
    {
        if (ch == '\n') r += '\r';
        r += ch;
    }
    return r;
}

// Читает CRLF-файл; true, если загрузился без исключения.
bool load_binary(const std::string& content, fiber::Model& m)
{
    {
        std::ofstream out(kPath, std::ios::binary);
        out << content;
    }
    try
    {
        fiber::load_fiber_file(kPath, m);
    }
    catch (const std::exception& e)
    {
        std::printf("    %s\n", e.what());
        return false;
    }
    return true;
}

void test_parser_edge_cases()
{
    std::printf("test_parser_edge_cases:\n");
    const std::string B = kBase;
    // Повторный MATERIAL после PLASTIC: параметры неупругости сохраняются и
    // применяются к ПОСЛЕДНЕМУ определению материала.
    {
        const fiber::Model m = load(
            "MATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.30\n"
            "PLASTIC 0 SY=60e6 H=0.4e9\n"
            "MATERIAL 0 3.0e9 0.30 2.0e-4 1000.0 0.20\n"
            "FIBER 0 2\n0 0 0\n0.01 0 0\n");
        const fiber::FiberMaterial& mat = m.materials[0];
        check(near(mat.E, 3.0e9) && near(mat.d, 2.0e-4) && near(mat.eps_break, 0.20),
              "повторный MATERIAL: действует последнее определение");
        check(mat.plastic.enabled() && near(mat.plastic.sy, 60e6) && near(mat.plastic.H, 0.4e9),
              "повторный MATERIAL: PLASTIC сохранён");
    }
    // FIBER требует уже определённый MATERIAL, директива неупругости этого не меняет.
    expect_error("FIBER до MATERIAL (директива раньше)",
                 "PLASTIC 0 SY=1e6\nFIBER 0 2\n0 0 0\n0.01 0 0\nMATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.30\n",
                 "FIBER ссылается на неопределённый материал 0");
    {
        // Директива до MATERIAL, FIBER после MATERIAL -- нормально.
        const fiber::Model m = load(
            "PLASTIC 0 SY=1e6\nMATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.30\nFIBER 0 2\n0 0 0\n0.01 0 0\n");
        check(m.materials[0].plastic.enabled(), "директива до MATERIAL применяется");
    }
    // Директива для пропущенного id: MATERIAL 2 создаёт слоты 0 и 1 неопределёнными.
    expect_error("директива для пропущенного id",
                 "MATERIAL 2 2.0e9 0.35 1.0e-4 1140.0 0.30\nFIBER 2 2\n0 0 0\n0.01 0 0\nPLASTIC 1 SY=1e6\n",
                 "неопределённый материал 1");
    // Частичный разбор чисел и нечисловые значения.
    expect_error("SY=1e9abc", B + "PLASTIC 0 SY=1e9abc\n", "не конечное число");
    expect_error("SY=nan", B + "PLASTIC 0 SY=nan\n", "не конечное число");
    expect_error("SY=inf", B + "PLASTIC 0 SY=inf\n", "не конечное число");
    expect_error("SY=1e999", B + "PLASTIC 0 SY=1e999\n", "не конечное число");
    expect_error("SY= (пусто)", B + "PLASTIC 0 SY=\n", "ожидалось КЛЮЧ=значение");
    expect_error("=1 (нет ключа)", B + "PLASTIC 0 =1\n", "ожидалось КЛЮЧ=значение");
    // Регистр ключей и названия модели нечувствителен; ключевое слово директивы -- нет.
    {
        const fiber::Model m = load(B + "PLASTIC 0 sy=60e6 h=0.4e9\nRATE 0 cs c=40 p=5\n");
        check(near(m.materials[0].plastic.sy, 60e6) && near(m.materials[0].plastic.H, 0.4e9),
              "ключи sy=, h= в нижнем регистре");
        check(m.materials[0].rate.model == fiber::RateModel::CowperSymonds && near(m.materials[0].rate.a, 40) &&
                  near(m.materials[0].rate.b, 5),
              "модель cs и ключи c=, p= в нижнем регистре");
    }
    expect_error("ключевое слово директивы чувствительно к регистру", B + "plastic 0 SY=1e6\n",
                 "неизвестное ключевое слово 'plastic'");
    // PLASTIC_TABLE: лишний токен в заголовке и нецелое число строк.
    expect_error("PLASTIC_TABLE лишний токен в заголовке", B + "PLASTIC_TABLE 0 2 extra\n0 1e6\n0.1 2e6\n",
                 "лишний токен 'extra'");
    expect_error("PLASTIC_TABLE нецелое N", B + "PLASTIC_TABLE 0 2.5\n0 1e6\n0.1 2e6\n", "ожидалось число строк");
    // VISCO с нулевыми ETA/XI допустимо, но неупругости нет (inelastic() == false;
    // задокументировано в fib-file.md).
    {
        const fiber::Model m = load(B + "VISCO 0 ETA=0\n");
        check(!m.materials[0].inelastic(), "VISCO ETA=0 -> inelastic() == false");
        const fiber::Model m2 = load(B + "VISCO 0 XI=0\n");
        check(!m2.materials[0].inelastic(), "VISCO XI=0 -> inelastic() == false");
        const fiber::Model m3 = load(B + "VISCO 0 ETA=0 XI=0\n");
        check(!m3.materials[0].inelastic(), "VISCO ETA=0 XI=0 -> inelastic() == false");
    }
    // FAILURE сам по себе не делает материал неупругим, но включает критерий разрыва.
    {
        const fiber::Model m = load(B + "FAILURE 0 SIGMA_U=50e6\n");
        check(!m.materials[0].inelastic() && m.materials[0].has_failure(),
              "FAILURE без PLASTIC/VISCO: inelastic() == false, has_failure() == true");
    }
    // CRLF-файл.
    {
        fiber::Model m;
        const bool ok = load_binary(to_crlf(B + "PLASTIC 0 SY=60e6 H=0.4e9   # c\nRATE 0 CS C=40 P=5\n"
                                                "VISCO 0 E1=0.5e9 TAU=1e-5\nFAILURE 0 EPS_P=0.2\n"),
                                    m);
        check(ok, "CRLF-файл читается");
        if (ok)
        {
            const fiber::FiberMaterial& mat = m.materials[0];
            check(near(mat.plastic.sy, 60e6) && near(mat.rate.a, 40) && near(mat.visco.E1, 0.5e9) &&
                      near(mat.failure.eps_p, 0.2),
                  "CRLF: все директивы разобраны");
        }
        fiber::Model m2;
        const bool ok2 = load_binary(to_crlf(B + "PLASTIC_TABLE 0 3\n0 1e6\n0.1 2e6\n0.2 3e6\n"), m2);
        check(ok2 && m2.materials[0].plastic.table_eps.size() == 3, "CRLF: PLASTIC_TABLE разобрана");
    }
}

// --- Доработки вывода (ревью этапа 3) ----------------------------------------
std::size_t bonds_line(const std::vector<std::string>& lines)
{
    std::size_t b = 0;
    while (b < lines.size() && lines[b] != "BONDS") ++b;
    return b;
}

void test_output_extended()
{
    std::printf("test_output_extended:\n");
    // (1) Неупругость только через VISCO ETA, ненулевые скорости (ветка eps_dot).
    {
        fiber::Model m = load(std::string(kBase) + "VISCO 0 ETA=50\n");
        check(m.materials[0].inelastic(), "VISCO ETA > 0 -> неупругий");
        const double E = m.materials[0].E;
        const double l0 = m.topology.bonds[1].l0;
        m.fiber.x[0][2] = 0.0202; // связь 1 растянута на 2%
        m.fiber.vx[0][1] = 1.0;
        m.fiber.vx[0][2] = 11.0; // относительная скорость вдоль связи 10 м/с
        const double l = m.fiber.x[0][2] - m.fiber.x[0][1];
        const double eps = l / l0 - 1.0;
        const double eps_dot = 10.0 / l0;
        const double sigma_exp = E * eps + 50.0 * eps_dot;
        fiber::fiber_output(m, 3, 0.0);
        const auto lines = read_lines("fiber_output_00003.txt");
        const std::size_t b = bonds_line(lines);
        const auto t1 = tokens(lines[b + 2]);
        check(t1.size() == 10, "VISCO-only: 10 колонок");
        if (t1.size() == 10)
        {
            const double sigma = std::stod(t1[7]);
            check(std::fabs(sigma - sigma_exp) <= 1e-12 * sigma_exp, "sigma = E eps + eta eps_dot (ветка eps_dot): " + t1[7]);
            check(std::stod(t1[6]) == 0.0 && std::stod(t1[8]) == 0.0, "eps_p = w_p = 0 без пластики");
        }
        // Быстрая разгрузка: вязкий член отрицателен и велик -> sigma = max(sigma, 0) = 0.
        m.fiber.vx[0][2] = -1.0e6;
        fiber::fiber_output(m, 4, 0.0);
        const auto l2 = read_lines("fiber_output_00004.txt");
        const auto t2 = tokens(l2[bonds_line(l2) + 2]);
        check(t2.size() == 10 && std::stod(t2[7]) == 0.0, "быстрая разгрузка: sigma = 0");
    }
    // (2) Смешанная модель: материал 0 упругий, материал 1 неупругий.
    {
        fiber::Model m = load(
            "MATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.30\n"
            "MATERIAL 1 1.0e9 0.35 1.0e-4 1140.0 0.30\n"
            "PLASTIC 1 SY=60e6 H=0.4e9\n"
            "FIBER 0 3\n0 0 0\n0.01 0 0\n0.02 0 0\n"
            "FIBER 1 3\n0 0.05 0\n0.01 0.05 0\n0.02 0.05 0\n");
        check(!m.materials[0].inelastic() && m.materials[1].inelastic(), "смешанная модель: 0 упругий, 1 неупругий");
        m.fiber.x[0][1] = 0.0099; // упругая связь 0 сжата
        m.fiber.x[0][2] = 0.0202; // упругая связь 1 растянута на 3%
        m.fiber.x[0][5] = 0.0204; // неупругая связь (4-5) растянута на 4%
        m.topology.bonds[3].st.eps_p = 0.004;
        m.topology.bonds[3].st.w_p = 3.0;
        const double A = pi * 1.0e-4 * 1.0e-4 / 4.0;
        fiber::fiber_output(m, 5, 0.0);
        const auto lines = read_lines("fiber_output_00005.txt");
        const std::size_t b = bonds_line(lines);
        check(lines[b - 1] == "# i j l0 l l_break elongation eps_p sigma w_p w_v",
              "смешанная модель: расширенный заголовок");
        const auto t0 = tokens(lines[b + 1]); // упругая, сжата
        const auto t1 = tokens(lines[b + 2]); // упругая, растянута
        const auto t3 = tokens(lines[b + 4]); // неупругая
        check(t0.size() == 10 && t1.size() == 10 && t3.size() == 10, "смешанная модель: у всех связей 10 колонок");
        if (t0.size() == 10 && t1.size() == 10 && t3.size() == 10)
        {
            check(std::stod(t0[7]) == 0.0 && std::stod(t0[6]) == 0.0, "упругая сжатая связь: sigma = 0");
            const double eps1 = (0.0202 - 0.0099) / 0.01 - 1.0;
            check(std::fabs(std::stod(t1[7]) - 2.0e9 * eps1) <= 1e-12 * 2.0e9 * eps1 && std::stod(t1[6]) == 0.0 &&
                      std::stod(t1[8]) == 0.0 && std::stod(t1[9]) == 0.0,
                  "упругая растянутая связь: sigma = E eps, eps_p = w = 0");
            const double eps3 = (0.0204 - 0.01) / 0.01 - 1.0;
            check(near(std::stod(t3[6]), 0.004) &&
                      std::fabs(std::stod(t3[7]) - 1.0e9 * (eps3 - 0.004)) <= 1e-9 * 1.0e9 &&
                      std::fabs(std::stod(t3[8]) - 3.0 * A * 0.01) <= 1e-12 * 3.0 * A * 0.01,
                  "неупругая связь: eps_p, sigma = E(eps - eps_p), w_p");
        }
    }
    // (3) eps_p > eps (связь провисла): sigma = 0; точность новых колонок.
    {
        fiber::Model m = load(std::string(kBase) + "PLASTIC 0 SY=60e6 H=0.4e9\n");
        m.fiber.x[0][2] = 0.0203; // связь 1: eps = 0.03
        m.topology.bonds[1].st.eps_p = 0.05;
        m.topology.bonds[0].st.eps_p = 1.23456789012345e-5;
        fiber::fiber_output(m, 6, 0.0);
        const auto lines = read_lines("fiber_output_00006.txt");
        const std::size_t b = bonds_line(lines);
        const auto t0 = tokens(lines[b + 1]);
        const auto t1 = tokens(lines[b + 2]);
        check(t1.size() == 10 && std::stod(t1[7]) == 0.0 && near(std::stod(t1[6]), 0.05),
              "eps < eps_p: sigma = 0, eps_p выведено");
        check(t0.size() == 10 && std::fabs(std::stod(t0[6]) - 1.23456789012345e-5) <= 1e-15 * 1.23456789012345e-5,
              "eps_p выведен с 17 значащими цифрами");
        // Прежние колонки -- с точностью потока по умолчанию (6 значащих цифр).
        check(t0.size() == 10 && t0[3].size() <= 8, "прежние колонки прежней точности: l = " + t0[3]);
    }
}

} // namespace

int main()
{
    // Все файлы (фикстура .fib и fiber_output_*.txt) -- во временном каталоге.
    e2e::ScopedWorkDir wd("sph_fiber_inelastic", "parser_output");
    test_no_directives_is_elastic();
    test_all_directives();
    test_table_and_second_material();
    test_errors();
    test_output_columns();
    test_parser_edge_cases();
    test_output_extended();
    if (failures == 0) std::printf("test_fiber_input_inelastic: OK\n");
    else std::printf("test_fiber_input_inelastic: %d failure(s)\n", failures);
    return failures;
}
