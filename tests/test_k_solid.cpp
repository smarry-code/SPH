// test_k_solid.cpp
// Unit-тесты Solid-сборщика импорта LS-DYNA .k (readk/k-import-plan.md, этап K4,
// Р3-Р6): NEND, MASS<0, _VOLUME, единицы, Р6 для своей dim (SPH_DIM = 1/2/3),
// hsml по умолчанию и из таблицы, urs выполнен, body = PID, k_monitor_nid,
// v= с лишними компонентами, PID без SOLID_PART, пропуск балок, режим без
// таблицы, сквозной load_k_geometry через временные файлы.
#include "k_input/k_geometry.hpp"
#include "param.hpp"
#include "simulation_config.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
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

constexpr int kDim = static_cast<int>(dim);

// Узлы вдоль оси x; элементы: 1..2 (NEND), 3 (MASS<0), 4 (_VOLUME, PID 2).
const char* kDeckText =
    "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n4,3,0,0\n"
    "*PART\nA\n1\n*PART\nB\n2\n"
    "*ELEMENT_SPH\n1,1,0.5,2\n3,1,-2.0\n"
    "*ELEMENT_SPH_VOLUME\n4,2,4.0\n";

std::string v_for_dim()
{
    return kDim == 1 ? "v=10" : kDim == 2 ? "v=10 20" : "v=10 20 30";
}

std::string table_text()
{
    return "SOLID_PART 1 5 1000 " + v_for_dim() + "\nSOLID_PART 2 10 2510 hsml=0.5\n";
}

struct Result
{
    ParticleSystem sys;
    std::size_t moni = 777;
    KImportSummary sum;
};

Result build(const std::string& deck_text, const std::string* table_src, KImportOptions opts = {})
{
    Result r;
    const KDeck deck = read_k_deck_from_string(deck_text, "mem.k");
    SimulationConfig cfg;
    cfg.apply_defaults();
    if (table_src)
    {
        const KPartsTable t = read_k_parts_from_string(*table_src, "t.parts");
        build_solid_from_deck(deck, &t, opts, cfg.mats, r.sys, r.moni, r.sum);
    }
    else
        build_solid_from_deck(deck, nullptr, opts, cfg.mats, r.sys, r.moni, r.sum);
    return r;
}

bool near(double a, double b) { return std::fabs(a - b) <= 1e-12 * std::fmax(1.0, std::fabs(b)); }

void test_basic()
{
    std::printf("[basic dim=%d]\n", kDim);
    const std::string tab = table_text();
    Result r = build(kDeckText, &tab);
    const ParticleSystem& s = r.sys;
    check(s.size() == 4, "NEND: 4 частицы");
    check(r.sum.sph_particles == 4 && r.sum.nodes == 4, "сводка: 4 узла, 4 частицы");
    check(s.x[0][0] == 0.0 && s.x[0][1] == 1.0 && s.x[0][3] == 3.0, "x по порядку файла");
    check(r.sum.bbox_min[0] == 0.0 && r.sum.bbox_max[0] == 3.0, "габарит по x");
    check(s.mass[0] == 0.5 && s.mass[1] == 0.5, "MASS>0 -> масса");
    const double vol_dim = 1.0;  // length_scale = 1
    check(near(s.mass[2], 2.0 * vol_dim * 1000.0), "MASS<0 -> |MASS|*rho");
    check(near(s.mass[3], 4.0 * vol_dim * 2510.0), "_VOLUME -> |MASS|*rho (PID 2)");
    check(s.itype[0] == 5 && s.itype[3] == 10, "itype из таблицы");
    check(s.body[0] == 1 && s.body[2] == 1 && s.body[3] == 2, "body = PID");
    check(s.rho[0] == 1000.0 && s.rho_0[0] == 1000.0 && s.rho[3] == 2510.0, "rho = rho_0 из таблицы");
    const double hs = 2.0 * std::pow(0.5 / 1000.0, 1.0 / kDim);
    check(near(s.hsml[0], hs), "hsml по умолчанию 2*(m/rho)^(1/dim)");
    check(s.hsml[3] == 0.5, "hsml из таблицы");
    check(s.vx[0][0] == 10.0 && s.vx[0][3] == 0.0, "vx из v=");
    if (kDim >= 2) check(s.vx[1][0] == 20.0, "vy из v=");
    if (kDim >= 3) check(s.vx[2][0] == 30.0, "vz из v=");
    bool urs_ok = true;
    for (std::size_t i = 0; i < s.size(); ++i)
        if (!(s.c[i] > 0.0) || !std::isfinite(s.p[i])) urs_ok = false;
    check(urs_ok, "urs выполнен: c > 0 (itype 5 и 10)");
    check(r.sum.dx_median > 0.0, "dx_median > 0");
    check(r.moni == 777, "moni_particle не тронут без monitor_nid");
}

void test_units()
{
    std::printf("[units dim=%d]\n", kDim);
    const std::string tab = table_text();
    KImportOptions o;
    o.length_scale = 1e-3;
    o.mass_scale = 2.0;
    Result r = build(kDeckText, &tab, o);
    check(near(r.sys.x[0][1], 1e-3) && near(r.sys.x[0][3], 3e-3), "length_scale на координаты");
    check(near(r.sum.bbox_max[0], 3e-3), "габарит в метрах");
    check(near(r.sys.mass[0], 1.0), "mass_scale на MASS>0");
    check(near(r.sys.mass[2], 2.0 * std::pow(1e-3, kDim) * 1000.0), "length_scale^dim на объём");
    check(near(r.sys.hsml[0], 2.0 * std::pow(1.0 / 1000.0, 1.0 / kDim)), "hsml по масштабированной массе");
}

void test_dim()
{
    std::printf("[dim rule dim=%d]\n", kDim);
    const std::string tab = table_text();
    // Лишняя координата: dim=1 -> y, dim=2 -> z.
    const char* bad = kDim == 1 ? "*NODE\n1,0,1,0\n2,1,0,0\n" : "*NODE\n1,0,0,1\n2,1,0,0\n";
    const std::string dtext = std::string(bad) + "*PART\nA\n1\n*ELEMENT_SPH\n1,1,1.0\n2,1,1.0\n";
    const std::string t1 = "SOLID_PART 1 5 1000\n";
    if (kDim < 3)
    {
        const std::string e = error_of([&] { build(dtext, &t1); });
        check(has(e, "mem.k:") && has(e, "NID 1") && has(e, "dim=" + std::to_string(kDim)), "ненулевая лишняя координата: " + e);
        // Малая добавка в пределах допуска - успех.
        const char* ok = kDim == 1 ? "*NODE\n1,0,1e-12,0\n2,5,0,0\n" : "*NODE\n1,0,0,1e-12\n2,5,0,0\n";
        const std::string otext = std::string(ok) + "*PART\nA\n1\n*ELEMENT_SPH\n1,1,1.0\n2,1,1.0\n";
        check(error_of([&] { build(otext, &t1); }).empty(), "лишняя координата в допуске - успех");
    }
    else
    {
        const std::string dtext3 = "*NODE\n1,0,0,1\n2,1,1,0\n*PART\nA\n1\n*ELEMENT_SPH\n1,1,1.0\n2,1,1.0\n";
        Result r = build(dtext3, &t1);
        check(r.sys.x[2][0] == 1.0 && r.sys.x[1][1] == 1.0, "dim=3: все три координаты");
    }
    // v= с лишними компонентами.
    const std::string tz = "SOLID_PART 1 5 1000 v=1 0 0\nSOLID_PART 2 10 2510\n";
    check(error_of([&] { build(kDeckText, &tz); }).empty(), "v= с нулевыми лишними компонентами - успех");
    const std::string tnz = "SOLID_PART 1 5 1000 v=1 2 3\nSOLID_PART 2 10 2510\n";
    const std::string e = error_of([&] { build(kDeckText, &tnz); });
    if (kDim < 3)
        check(has(e, "t.parts:1") && has(e, "v="), "v= с ненулевой лишней компонентой -> ошибка: " + e);
    else
        check(e.empty(), "dim=3: v=1 2 3 допустимо");
}

void test_errors_monitor_skip()
{
    std::printf("[errors dim=%d]\n", kDim);
    const std::string tab = table_text();
    KImportOptions o;
    o.monitor_nid = 3;
    Result r = build(kDeckText, &tab, o);
    check(r.moni == 2, "k_monitor_nid=3 -> индекс 2");
    o.monitor_nid = 99;
    std::string e = error_of([&] { build(kDeckText, &tab, o); });
    check(has(e, "99") && has(e, "monitor"), "нет такого NID -> ошибка: " + e);

    const std::string no2 = "SOLID_PART 1 5 1000\n";
    e = error_of([&] { build(kDeckText, &no2); });
    check(has(e, "PID 2") && has(e, "SOLID_PART"), "PID без SOLID_PART -> ошибка: " + e);

    e = error_of([&] { build("*NODE\n1,0,0,0\n*PART\nA\n1\n*ELEMENT_SPH\n1,1,0\n", &no2); });
    check(has(e, "MASS = 0") && has(e, "mem.k:"), "MASS = 0 -> ошибка: " + e);

    const std::string deck_b = std::string(kDeckText)
        + "*NODE\n6,0,0,0\n7,0.5,0,0\n*PART\nF\n3\n*ELEMENT_BEAM\n1,3,6,7\n";
    r = build(deck_b, &tab);
    bool beam_warn = false, solid_warn = false;
    for (const auto& w : r.sum.warnings)
        if (has(w, "ELEMENT_BEAM") && has(w, "PID: 3") && has(w, "1 элементов")) beam_warn = true;
    check(beam_warn && r.sys.size() == 4, "балки пропущены с предупреждением");
    r = build(std::string(kDeckText) + "*NODE\n8,0,0,0\n*ELEMENT_SOLID\n1,1,1,1,1,1\n", &tab);
    for (const auto& w : r.sum.warnings)
        if (has(w, "ELEMENT_SOLID")) solid_warn = true;
    check(solid_warn, "*ELEMENT_SOLID пропущены с предупреждением");

    // Режим без таблицы.
    r = build("*NODE\n1,0,0,0\n2,1,0,0\n*PART\nA\n1\n*PART\nB\n2\n*ELEMENT_SPH\n1,1,0.5\n2,2,0.25\n", nullptr);
    check(r.sys.size() == 2 && r.sys.itype[1] == 2 && r.sys.body[1] == 2 && r.sys.mass[0] == 0.5, "без таблицы: itype=body=PID");
    check(r.sys.rho[0] == 0.0 && r.sys.hsml[0] == 0.0 && r.sys.c[0] == 0.0 && r.sum.dx_median == 0.0,
          "без таблицы: rho/hsml/c нули, urs не вызван");
    e = error_of([&] { build(kDeckText, nullptr); });
    check(has(e, "MASS-объём"), "без таблицы: MASS<0 -> ошибка: " + e);
}

void test_load()
{
    std::printf("[load dim=%d]\n", kDim);
    {
        std::ofstream("__k_solid.k") << kDeckText;
        std::ofstream("__k_solid.parts") << table_text();
    }
    SimulationConfig cfg;
    cfg.apply_defaults();
    fiber::Model model;
    std::size_t moni = 0;
    KImportOptions o;
    o.k_file = "__k_solid.k";
    o.parts_file = "__k_solid.parts";
    o.monitor_nid = 4;
    KImportSummary s = load_k_geometry(o, cfg.mats, model, moni);
    check(model.solid.size() == 4 && moni == 3 && s.sph_particles == 4, "load_k_geometry: 4 частицы, монитор 3");

    { std::ofstream("__k_solid_bad.parts") << "SOLID_PART 1 5 1000\n"; }
    o.parts_file = "__k_solid_bad.parts";
    const std::string e = error_of([&] { load_k_geometry(o, cfg.mats, model, moni); });
    check(has(e, "PID 2") && has(e, "SOLID_PART"), "load: PID без SOLID_PART: " + e);

    o.parts_file.clear();
    o.mode = KMode::SolidFiber;
    check(!error_of([&] { load_k_geometry(o, cfg.mats, model, moni); }).empty(), "solid_fiber без таблицы -> ошибка");
    std::remove("__k_solid.k");
    std::remove("__k_solid.parts");
    std::remove("__k_solid_bad.parts");
}

} // namespace

int main()
{
    try
    {
        test_basic();
        test_units();
        test_dim();
        test_errors_monitor_skip();
        test_load();
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
