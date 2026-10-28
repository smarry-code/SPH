// test_k_e2e.cpp
// Сквозной тест подключения импорта .k к input()/time_integration (этап K6,
// readk/k-import-plan.md, Р4, Р5, Р10, раздел 4.4): Initial_Conditions.txt с
// ключами k_* -> input() -> несколько шагов time_integration(Model).
// input() читает Initial_Conditions.txt из текущего каталога, а цикл пишет
// кадры и log.txt туда же, поэтому тест работает во временном каталоге и
// восстанавливает cwd. SPH_DIM=3.
#include "fiber/fiber_input.hpp"
#include "input.hpp"
#include "run_guard.hpp"
#include "time_integration.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool bit_eq(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

std::string read_all(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("не открыть " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void write_file(const fs::path& p, const std::string& text)
{
    std::ofstream out(p, std::ios::binary);
    out << text;
}

// Временный рабочий каталог с восстановлением cwd.
class WorkDir
{
public:
    WorkDir()
    {
        old_ = fs::current_path();
        dir_ = fs::temp_directory_path() / "sph_k_e2e_work";
        std::error_code ec;
        fs::remove_all(dir_, ec);
        fs::create_directories(dir_);
        fs::current_path(dir_);
    }
    ~WorkDir()
    {
        std::error_code ec;
        fs::current_path(old_, ec);
        fs::remove_all(dir_, ec);
    }

private:
    fs::path old_, dir_;
};

// Сетка 4x4x2 SPH-узлов (NID 1..32, PID 1) с шагом 1 мм и две нити по 4 узла
// (PID 2, NID 101..104 и 201..204), в стороне от сетки.
std::string make_deck(bool with_beams)
{
    std::ostringstream k;
    k << "*KEYWORD\n*PART\nplate\n1\n*PART\nfiber\n2\n*NODE\n";
    int nid = 1;
    for (int iz = 0; iz < 2; ++iz)
        for (int iy = 0; iy < 4; ++iy)
            for (int ix = 0; ix < 4; ++ix)
                k << nid++ << "," << ix * 1.0e-3 << "," << iy * 1.0e-3 << "," << iz * 1.0e-3 << "\n";
    for (int f = 0; f < 2; ++f)
        for (int j = 0; j < 4; ++j)
            k << (100 * (f + 1) + 1 + j) << "," << (0.02 + j * 1.0e-3) << "," << f * 0.01 << ",0.0\n";
    k << "*ELEMENT_SPH\n";
    for (int n = 1; n <= 32; ++n) k << n << ",1,2.51e-6\n";
    if (with_beams)
    {
        k << "*ELEMENT_BEAM\n";
        int e = 1;
        for (int f = 0; f < 2; ++f)
            for (int j = 0; j < 3; ++j)
                k << e++ << ",2," << (100 * (f + 1) + 1 + j) << "," << (100 * (f + 1) + 2 + j) << "\n";
    }
    k << "*END\n";
    return k.str();
}

const char* kTable =
    "SOLID_PART 1 10 2510 hsml=0.0013\n"
    "MATERIAL 0 2.0e9 0.35 1.0e-4 1140 0.02\n"
    "FIBER_PART 2 0\n";

const char* kFib =
    "MATERIAL 0 2.0e9 0.35 1.0e-4 1140 0.02\n"
    "FIBER 0 3\n0.05 0.0 0.0\n0.051 0.0 0.0\n0.052 0.0 0.0\n"
    "FIBER 0 2\n0.05 0.01 0.0\n0.05 0.011 0.0\n";

// Initial_Conditions.txt: общая часть из tests/data + ключи k_*.
void write_ic(const std::string& extra)
{
    std::string base = read_all(std::string(SPH_TEST_DATA_DIR) + "/guo2019_initial_conditions.txt");
    write_file("Initial_Conditions.txt", base + "\n" + extra + "\n");
}

template <class F>
std::string error_of(F&& f)
{
    try { f(); }
    catch (const std::exception& e) { return e.what(); }
    return std::string();
}

bool solid_finite(const ParticleSystem& s, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < dim; ++d)
            if (!std::isfinite(s.x[d][i]) || !std::isfinite(s.vx[d][i])) return false;
    return true;
}

bool fiber_finite(const fiber::FiberState& f)
{
    for (std::size_t k = 0; k < f.size(); ++k)
        for (std::size_t d = 0; d < dim; ++d)
            if (!std::isfinite(f.x[d][k]) || !std::isfinite(f.vx[d][k])) return false;
    return true;
}

// Несколько шагов основного цикла; как main: расширение model.solid до maxn.
bool run_steps(SimulationConfig& cfg, fiber::Model& model, std::size_t ntotal, std::size_t steps)
{
    cfg.maxn = static_cast<std::size_t>(ntotal * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    const std::size_t prev = model.solid.size();
    if (prev < cfg.maxn)
    {
        model.solid.resize(cfg.maxn);
        for (std::size_t i = prev; i < model.solid.size(); ++i) model.solid.inzializ_index(i);
    }
    cfg.t_total = 1.0;
    double dt = compute_cfl_timestep(cfg, model.solid, ntotal);
    double time = 0.0;
    RunGuard guard(0.0, 0.0, 0);
    guard.set_step_limit(steps);
    time_integration(cfg, model, ntotal, dt, time, 0.0, &guard);
    return guard.reason() == RunGuard::Reason::StepLimit && guard.steps() == steps;
}

void test_solid()
{
    std::printf("[k_mode=solid]\n");
    write_file("deck.k", make_deck(false));
    write_file("deck.parts", kTable);
    write_ic("k_file = deck.k\nk_mode = solid\nk_parts_file = deck.parts\nk_monitor_nid = 5\n");
    SimulationConfig cfg;
    cfg.apply_defaults();
    fiber::Model model;
    std::size_t ntotal = 0;
    double dx = 0.0;
    input(cfg, model, ntotal, dx);
    check(ntotal == 32 && model.solid.size() == 32, "загружено 32 частицы");
    check(model.fiber.size() == 0, "нитей нет");
    check(model.solid.itype[0] == 10 && model.solid.rho[0] == 2510.0, "itype 10, rho 2510");
    check(model.solid.c[0] > 0.0, "urs выполнен (c > 0)");
    check(dx > 0.0 && std::fabs(dx - 1.0e-3) < 1.0e-4, "dx_median ~ 1 мм");
    check(cfg.moni_particle == 4, "k_monitor_nid -> cfg.moni_particle (NID 5 -> индекс 4)");
    check(run_steps(cfg, model, ntotal, 5), "5 шагов выполнены");
    check(solid_finite(model.solid, ntotal), "после шагов всё конечно");
}

void test_solid_fiber()
{
    std::printf("[k_mode=solid_fiber]\n");
    write_file("deck.k", make_deck(true));
    write_file("deck.parts", kTable);
    write_ic("k_file = deck.k\nk_mode = solid_fiber\nk_parts_file = deck.parts\n");
    SimulationConfig cfg;
    cfg.apply_defaults();
    fiber::Model model;
    std::size_t ntotal = 0;
    double dx = 0.0;
    input(cfg, model, ntotal, dx);
    check(ntotal == 32, "32 частицы Solid");
    check(model.fiber.size() == 8 && model.topology.bonds.size() == 6, "2 нити: 8 бусин, 6 связей");
    check(run_steps(cfg, model, ntotal, 5), "5 шагов выполнены");
    check(solid_finite(model.solid, ntotal) && fiber_finite(model.fiber), "после шагов всё конечно");
}

void test_k_plus_fib()
{
    std::printf("[.k (solid) + fiber_file]\n");
    write_file("deck.k", make_deck(false));
    write_file("deck.parts", kTable);
    write_file("a.fib", kFib);

    write_ic("k_file = deck.k\nk_parts_file = deck.parts\n");
    SimulationConfig cfg0;
    cfg0.apply_defaults();
    fiber::Model base;
    std::size_t n0 = 0;
    double dx0 = 0.0;
    input(cfg0, base, n0, dx0);

    write_ic("k_file = deck.k\nk_parts_file = deck.parts\nfiber_file = a.fib\n");
    SimulationConfig cfg;
    cfg.apply_defaults();
    fiber::Model model;
    std::size_t n = 0;
    double dx = 0.0;
    input(cfg, model, n, dx);

    fiber::Model direct;
    fiber::load_fiber_file("a.fib", direct);
    bool same_fiber = model.fiber.size() == direct.fiber.size() &&
                      model.topology.bonds.size() == direct.topology.bonds.size();
    if (same_fiber)
    {
        for (std::size_t k = 0; k < direct.fiber.size(); ++k)
        {
            same_fiber = same_fiber && bit_eq(model.fiber.mass[k], direct.fiber.mass[k]) &&
                         model.fiber.object_id[k] == direct.fiber.object_id[k];
            for (std::size_t d = 0; d < dim; ++d)
                same_fiber = same_fiber && bit_eq(model.fiber.x[d][k], direct.fiber.x[d][k]);
        }
        for (std::size_t b = 0; b < direct.topology.bonds.size(); ++b)
            same_fiber = same_fiber && bit_eq(model.topology.bonds[b].l0, direct.topology.bonds[b].l0);
    }
    check(same_fiber && model.fiber.size() == 5, "model.fiber совпадает с load_fiber_file (5 бусин)");
    bool same_solid = n == n0 && model.solid.size() == base.solid.size();
    for (std::size_t i = 0; same_solid && i < n0; ++i)
    {
        same_solid = bit_eq(model.solid.mass[i], base.solid.mass[i]) &&
                     bit_eq(model.solid.hsml[i], base.solid.hsml[i]) && model.solid.itype[i] == base.solid.itype[i];
        for (std::size_t d = 0; d < dim; ++d)
            same_solid = same_solid && bit_eq(model.solid.x[d][i], base.solid.x[d][i]);
    }
    check(same_solid, "Solid совпадает с импортом без fiber_file");
    check(run_steps(cfg, model, n, 3), "3 шага выполнены");
}

void test_errors()
{
    std::printf("[ошибки конфигурации]\n");
    write_file("deck.k", make_deck(true));
    write_file("deck.parts", kTable);
    write_file("a.fib", kFib);
    auto run_ic = [&](const std::string& extra)
    {
        write_ic(extra);
        SimulationConfig cfg;
        cfg.apply_defaults();
        fiber::Model model;
        std::size_t n = 0;
        double dx = 0.0;
        input(cfg, model, n, dx);
    };
    const std::string base = "k_file = deck.k\nk_parts_file = deck.parts\n";
    std::string e = error_of([&] { run_ic(base + "k_mode = solid_fiber\nfiber_file = a.fib\n"); });
    check(e.find("два источника нитей") != std::string::npos, "solid_fiber + fiber_file -> ошибка: " + e);
    e = error_of([&] { run_ic(base + "k_mode = bogus\n"); });
    check(e.find("k_mode") != std::string::npos, "неверный k_mode -> ошибка: " + e);
    e = error_of([&] { run_ic(base + "k_length_scale = 0\n"); });
    check(e.find("k_length_scale") != std::string::npos, "k_length_scale=0 -> ошибка: " + e);
    e = error_of([&] { run_ic(base + "k_mass_scale = abc\n"); });
    check(e.find("k_mass_scale") != std::string::npos, "k_mass_scale=abc -> ошибка: " + e);
    e = error_of([&] { run_ic(base + "k_solid_elements = x\n"); });
    check(e.find("k_solid_elements") != std::string::npos, "неверный k_solid_elements -> ошибка: " + e);
    e = error_of([&] { run_ic(base + "k_monitor_nid = 1.5\n"); });
    check(e.find("k_monitor_nid") != std::string::npos, "нецелый k_monitor_nid -> ошибка: " + e);
}

} // namespace

int main()
{
    try
    {
        WorkDir wd;
        test_solid();
        test_solid_fiber();
        test_k_plus_fib();
        test_errors();
    }
    catch (const std::exception& e)
    {
        std::printf("  FAIL: исключение: %s\n", e.what());
        ++failures;
    }
    std::printf("%s (%d)\n", failures == 0 ? "PASS" : "FAIL", failures);
    return failures == 0 ? 0 : 1;
}
