// test_sph_body_isolation.cpp
// Опция ядра sph_same_body_only (изоляция SPH-пар Solid по телам;
// ТЗ эксперименты Guo 2019, §4.2).
//
// Два куба Al (itype 4) 4x4x4 частиц, dx = 1 мм, hsml = 2 dx, сложены гранью
// к грани по непрерывной решётке; тело 1 и тело 2.
//
// ПРОВЕРКИ:
//   1  поиск пар (build_solid_sph_search): с флагом пар между кубами РОВНО
//      ноль, без флага -- больше нуля; множества пар внутри каждого куба
//      совпадают (отсортированные списки (i, j)) при флаге и без него;
//   2  flag = true при contact_force_type = 0 -> std::runtime_error;
//   3  20 шагов time_integration (PC, contact_force_type = 1, кубы летят друг
//      на друга со скоростью 10 м/с):
//        A  флаг выключен, оба куба body = 1;
//        B  флаг включён,  оба куба body = 1 (ни одна пара не отвергается)
//           -- x, v, rho, p, u совпадают с A ПОБИТНО: флаг сам по себе не
//           вносит численных побочных эффектов;
//        C  флаг включён, тела 1 и 2 -- значения конечны, результат
//           отличается от A (печатается max |dv|).
//
// ПОБИТОВАЯ идентичность режима «флаг выключен» с кодом ДО изменения здесь
// не проверяется: её обеспечивают эталоны `ctest -L regression`.
//
// Для побитового сравнения (A против B) число потоков OpenMP фиксируется
// равным 1 (параллельная редукция с ULP-расхождениями иначе даёт шум).
#include "input.hpp"
#include "int_force.hpp"
#include "interactions/solid_sph.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/cpu_executor.hpp"
#include "run_guard.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"

#include "e2e_harness.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

static_assert(dim == 3, "test_sph_body_isolation требует сборки с SPH_DIM=3");

#ifndef SPH_TEST_DATA_DIR
#define SPH_TEST_DATA_DIR "tests/data"
#endif

namespace fs = std::filesystem;
using e2e::check;
using e2e::note;

namespace
{

constexpr std::size_t kSide = 4;                 // частиц на ребре куба
constexpr double kDx = 1.0e-3;                   // 1 мм
constexpr double kRhoAl = 2700.0;
constexpr int kItypeAl = 4;
constexpr double kApproachSpeed = 10.0;          // м/с, каждый куб
constexpr std::size_t kCubeN = kSide * kSide * kSide;
constexpr std::size_t kN = 2 * kCubeN;
constexpr std::size_t kSteps = 20;

// Конфигурация: как в других сквозных тестах (разбор файла BC2 даёт рабочие
// skf/density_method/pa_sph и материалы), затем точечные переопределения.
SimulationConfig make_config(bool same_body_only, int contact_type)
{
    SimulationConfig cfg;
    cfg.apply_defaults();
    const fs::path config_path = fs::path(SPH_TEST_DATA_DIR) / "bc2_initial_conditions.txt";
    if (fs::exists(config_path))
    {
        e2e::ScopedSilence silence(true);
        parse_sph_config(cfg, config_path.string());
    }
    cfg.ode_solver = 0;                           // PredictorCorrector
    cfg.skf = 1;                                  // кубический сплайн, носитель 2h
    cfg.contact_force_type = contact_type;
    cfg.sph_same_body_only = same_body_only;
    cfg.num_threads = 1;
    cfg.maxn = static_cast<std::size_t>(1.4 * static_cast<double>(kN)) + 1;
    cfg.max_interaction = 100 * cfg.maxn;
    cfg.moni_particle = 0;
    cfg.save_t = 1.0e30;                          // кадр только на первом шаге (во временный каталог)
    cfg.print_step = 1000000000;
    return cfg;
}

// Два куба рядом по оси x: частицы 0..63 -- первый, 64..127 -- второй.
// Шаг решётки dx и между кубами (непрерывная решётка 8 x 4 x 4).
ParticleSystem make_system(const SimulationConfig& cfg, int body_a, int body_b, double speed)
{
    ParticleSystem sys;
    sys.resize(cfg.maxn);
    for (std::size_t i = 0; i < cfg.maxn; ++i) sys.inzializ_index(i);

    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t ix = 0; ix < kSide; ++ix)
            for (std::size_t iy = 0; iy < kSide; ++iy)
                for (std::size_t iz = 0; iz < kSide; ++iz)
                {
                    const std::size_t k = c * kCubeN + (ix * kSide + iy) * kSide + iz;
                    sys.x[0][k] = (static_cast<double>(c * kSide + ix) + 0.5) * kDx;
                    sys.x[1][k] = (static_cast<double>(iy) + 0.5) * kDx;
                    sys.x[2][k] = (static_cast<double>(iz) + 0.5) * kDx;
                    sys.vx[0][k] = (c == 0) ? speed : -speed;
                    sys.itype[k] = kItypeAl;
                    sys.body[k] = (c == 0) ? body_a : body_b;
                    sys.rho[k] = kRhoAl;
                    sys.rho_0[k] = kRhoAl;
                    sys.mass[k] = kRhoAl * kDx * kDx * kDx;
                    sys.hsml[k] = 2.0 * kDx;
                }
    for (std::size_t i = 0; i < kN; ++i) urs(cfg.mats, sys, i);
    return sys;
}

using PairList = std::vector<std::pair<std::size_t, std::size_t>>;

struct SearchOut
{
    PairList cross;   // пары между телами
    PairList inside;  // пары внутри тела (i<j после нормализации, отсортированы)
    PairList contact; // real_pairs как есть (в порядке поиска): по ним работает Solid-Solid-Contact
};

// Поиск пар через ту же цепочку, что single_step (без ghost: в 3D их нет).
// Исключение build_solid_sph_search пробрасывается вызывающему.
SearchOut run_search(const SimulationConfig& cfg, const ParticleSystem& sys)
{
    rhs::CpuExecutor exec;
    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec, cfg.mats, sys, kN, derived);
    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, kN, ghosts);
    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, kN, derived, ghosts, ext);
    interactions::SolidSphSearchResult search;
    interactions::build_solid_sph_search(exec, cfg, ext, search);

    SearchOut out;
    for (const auto& rp : search.real_pairs)
        out.contact.emplace_back(static_cast<std::size_t>(rp.first), static_cast<std::size_t>(rp.second));
    for (const auto& pr : search.primitives)
    {
        const std::size_t i = pr.i, j = pr.j;
        if (j >= kN) continue; // ghost-пар в 3D нет, защита
        const auto p = std::make_pair(std::min(i, j), std::max(i, j));
        if (sys.body[i] != sys.body[j]) out.cross.push_back(p);
        else out.inside.push_back(p);
    }
    std::sort(out.cross.begin(), out.cross.end());
    std::sort(out.inside.begin(), out.inside.end());
    return out;
}

// Прямой подсчёт пар между кубами по номерам (не по body): cube = индекс / kCubeN.
std::size_t count_cross_by_index(const SimulationConfig& cfg, const ParticleSystem& sys)
{
    const SearchOut so = run_search(cfg, sys);
    std::size_t n = 0;
    for (const auto& p : so.cross) if ((p.first / kCubeN) != (p.second / kCubeN)) ++n;
    for (const auto& p : so.inside) if ((p.first / kCubeN) != (p.second / kCubeN)) ++n;
    return n;
}

void test_search()
{
    std::printf("test_search:\n");
    SimulationConfig cfg_off = make_config(false, 1);
    SimulationConfig cfg_on = make_config(true, 1);
    const ParticleSystem sys = make_system(cfg_off, 1, 2, 0.0);

    const SearchOut off = run_search(cfg_off, sys);
    const SearchOut on = run_search(cfg_on, sys);
    const std::size_t cross_off = count_cross_by_index(cfg_off, sys);
    const std::size_t cross_on = count_cross_by_index(cfg_on, sys);

    note("пар между кубами: без флага " + std::to_string(cross_off) + ", с флагом " + std::to_string(cross_on));
    note("пар внутри кубов: без флага " + std::to_string(off.inside.size()) + ", с флагом " +
         std::to_string(on.inside.size()));

    check(cross_off > 0, "без флага пары между кубами есть");
    check(cross_on == 0, "с флагом пар между кубами ровно ноль");
    check(off.cross.size() == cross_off && on.cross.empty(),
          "разбиение по body согласуется с разбиением по номеру куба");
    check(!on.inside.empty(), "пары внутри кубов сохранены");
    check(on.inside == off.inside, "множества пар внутри каждого куба совпадают при флаге и без него");
    // Контакт Solid-Solid берёт пары из real_pairs того же поиска: изоляция
    // убирает только SPH-primitive, список контактных пар обязан остаться прежним
    // (иначе разные тела не взаимодействуют вовсе -- найдено smoke-тестом Guo 2019).
    std::size_t contact_cross_on = 0;
    for (const auto& p : on.contact) if (sys.body[p.first] != sys.body[p.second]) ++contact_cross_on;
    note("контактных пар (real_pairs) между кубами с флагом: " + std::to_string(contact_cross_on));
    check(contact_cross_on == cross_off, "с флагом пары между кубами остаются в real_pairs для контакта");
    // Порядок пар у find_pairs_parallel зависит от потоков -- сравниваются множества.
    PairList cs_on = on.contact, cs_off = off.contact;
    std::sort(cs_on.begin(), cs_on.end());
    std::sort(cs_off.begin(), cs_off.end());
    check(cs_on == cs_off, "множества real_pairs с флагом и без него совпадают");
}

void test_exception()
{
    std::printf("test_exception:\n");
    SimulationConfig cfg = make_config(true, 0);
    const ParticleSystem sys = make_system(cfg, 1, 2, 0.0);
    bool threw = false;
    try
    {
        (void)run_search(cfg, sys);
    }
    catch (const std::runtime_error& ex)
    {
        threw = true;
        note(std::string("сообщение: ") + ex.what());
    }
    check(threw, "sph_same_body_only = 1 при contact_force_type = 0 бросает std::runtime_error");

    // Контрольный случай: тот же набор, флаг выключен -- исключения нет.
    SimulationConfig cfg_ok = make_config(false, 0);
    bool threw_ok = false;
    try
    {
        (void)run_search(cfg_ok, sys);
    }
    catch (const std::exception&)
    {
        threw_ok = true;
    }
    check(!threw_ok, "флаг выключен и contact_force_type = 0 -- исключения нет");
}

struct RunResult
{
    ParticleSystem sys;
    bool threw = false;
    std::string what;
    std::size_t steps = 0;
};

RunResult run_case(bool same_body_only, int body_a, int body_b, const std::string& tag)
{
    SimulationConfig cfg = make_config(same_body_only, 1);
    RunResult res;
    res.sys = make_system(cfg, body_a, body_b, kApproachSpeed);
    // Скорости: куб A едет в +x, куб B -- в -x (внутри make_system).
    double dt = compute_cfl_timestep(cfg, res.sys, kN);
    double time = 0.0;
    RunGuard guard(0.0, 0.0, 0);
    guard.set_step_limit(kSteps);
    {
        e2e::ScopedWorkDir wd("sph_body_isolation", tag);
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(cfg, res.sys, kN, dt, time, kDx, &guard);
        }
        catch (const std::exception& ex)
        {
            res.threw = true;
            res.what = ex.what();
        }
    }
    res.steps = guard.steps();
    return res;
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b, std::size_t n)
{
    return std::memcmp(a.data(), b.data(), n * sizeof(double)) == 0;
}

bool all_finite(const ParticleSystem& s, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d)
            if (!std::isfinite(s.x[d][i]) || !std::isfinite(s.vx[d][i])) return false;
        if (!std::isfinite(s.rho[i]) || !std::isfinite(s.p[i]) || !std::isfinite(s.u[i])) return false;
    }
    return true;
}

double max_dv(const ParticleSystem& a, const ParticleSystem& b, std::size_t n)
{
    double m = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t i = 0; i < n; ++i) m = std::max(m, std::fabs(a.vx[d][i] - b.vx[d][i]));
    return m;
}

void test_integration()
{
    std::printf("test_integration:\n");
#ifdef _OPENMP
    omp_set_num_threads(1); // для побитового сравнения A и B
#endif

    const RunResult A = run_case(false, 1, 1, "A");
    const RunResult B = run_case(true, 1, 1, "B");
    const RunResult C = run_case(true, 1, 2, "C");

    check(!A.threw, "прогон A (флаг выключен) без исключений" + (A.threw ? ": " + A.what : std::string()));
    check(!B.threw, "прогон B (флаг включён, одно тело) без исключений" + (B.threw ? ": " + B.what : std::string()));
    check(!C.threw, "прогон C (флаг включён, два тела) без исключений" + (C.threw ? ": " + C.what : std::string()));
    if (A.threw || B.threw || C.threw) return;

    check(A.steps == kSteps && B.steps == kSteps && C.steps == kSteps,
          "все прогоны сделали ровно " + std::to_string(kSteps) + " шагов (A=" + std::to_string(A.steps) +
              ", B=" + std::to_string(B.steps) + ", C=" + std::to_string(C.steps) + ")");

    bool bits = true;
    for (std::size_t d = 0; d < dim; ++d)
        bits = bits && same_bits(A.sys.x[d], B.sys.x[d], kN) && same_bits(A.sys.vx[d], B.sys.vx[d], kN);
    bits = bits && same_bits(A.sys.rho, B.sys.rho, kN) && same_bits(A.sys.p, B.sys.p, kN) &&
           same_bits(A.sys.u, B.sys.u, kN);
    check(bits, "B (флаг, одно тело) побитно совпадает с A (флаг выключен): x, v, rho, p, u");

    check(all_finite(C.sys, kN), "C: все значения конечны");
    const double dvAC = max_dv(A.sys, C.sys, kN);
    note("max |dv| между A и C = " + std::to_string(dvAC) + " м/с");
    check(dvAC > 0.0, "C отличается от A (пары между телами исключены из SPH-сумм)");
}

} // namespace

int main()
{
    std::printf("test_sph_body_isolation (sph_same_body_only)\n");
    test_search();
    test_exception();
    test_integration();

    if (e2e::failures == 0) std::printf("\nВСЕ ПРОВЕРКИ ПРОЙДЕНЫ\n");
    else std::printf("\nПРОВАЛЕНО проверок: %d\n", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
