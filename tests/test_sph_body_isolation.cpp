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
//           отличается от A (печатается max |dv|);
//   4  список межтельных пар search.contact_pairs: это подпоследовательность
//      real_pairs с body[i] != body[j] (тот же порядок, k -- индекс в
//      real_pairs); compute_solid_contact_diagnostics по нему при
//      CpuExecutor(1), (2), (7) даёт ровно те же pair_count, max_delta, min_dt,
//      что последовательный эталон по ВСЕМ real_pairs (кубы с проникновением).
//   5  сборка примитивов build_solid_sph_search параллельна (счёт -> префиксная
//      сумма -> запись по индексам): primitives, real_pairs, contact_pairs
//      ПОБАЙТНО равны последовательному эталону (копия прежнего цикла с
//      push_back) при CpuExecutor(1), (2), (7), с изоляцией и без, skf = 1..3;
//      повторные вызовы в один SolidSphSearchResult не перевыделяют буферы
//      (data()/capacity() постоянны, пока размер не превышает достигнутый
//      максимум), blocks выровнен по кэш-линии.
//
// ПОБИТОВАЯ идентичность режима «флаг выключен» с кодом ДО изменения здесь
// не проверяется: её обеспечивают эталоны `ctest -L regression`.
//
// Для побитового сравнения (A против B) число потоков OpenMP фиксируется
// равным 1 (параллельная редукция с ULP-расхождениями иначе даёт шум).
#include "input.hpp"
#include "int_force.hpp"
#include "interactions/solid_contact.hpp"
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
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

static_assert(dim == 3, "test_sph_body_isolation требует сборки с SPH_DIM=3");

#ifndef SPH_CASES_DIR
// Задаётся из CMake (функция sph_case). Фолбэк -- запуск из корня репозитория.
#define SPH_CASES_DIR "cases"
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
    const fs::path config_path = fs::path(SPH_CASES_DIR) / "bc2" / "regression" / "initial_conditions.txt";
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

// Список контактных пар и диагностика контакта по нему (parallel_reduce с
// порцией n/(8*nthreads)): сцена с проникновением -- второй куб сдвинут
// навстречу первому на 0.3 dx (радиус контакта R = dx/2, зазор на решётке 0).
void test_contact_pairs()
{
    std::printf("test_contact_pairs:\n");
    SimulationConfig cfg = make_config(true, 1);
    ParticleSystem sys = make_system(cfg, 1, 2, kApproachSpeed);
    for (std::size_t k = kCubeN; k < kN; ++k)
    {
        sys.x[0][k] -= 0.3 * kDx;
        // неоднородность скоростей: разные pdot, а значит и разные dt_limit пар
        sys.vx[1][k] += 0.5 * std::sin(static_cast<double>(k));
    }

    rhs::CpuExecutor exec1;
    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec1, cfg.mats, sys, kN, derived);
    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, kN, ghosts);
    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, kN, derived, ghosts, ext);
    interactions::SolidSphSearchResult search;
    interactions::build_solid_sph_search(exec1, cfg, ext, search);

    // contact_pairs == подпоследовательность real_pairs с body[i] != body[j]
    std::vector<interactions::SolidContactPair> expected;
    for (std::size_t k = 0; k < search.real_pairs.size(); ++k)
    {
        const rhs::index_t i = search.real_pairs[k].first, j = search.real_pairs[k].second;
        if (sys.body[i] != sys.body[j])
            expected.push_back({i, j, static_cast<std::uint64_t>(k)});
    }
    note("real_pairs: " + std::to_string(search.real_pairs.size()) + ", contact_pairs: " +
         std::to_string(search.contact_pairs.size()));
    check(!search.contact_pairs.empty(), "межтельные пары найдены");
    bool same_list = search.contact_pairs.size() == expected.size();
    for (std::size_t q = 0; same_list && q < expected.size(); ++q)
    {
        const auto& a = search.contact_pairs[q];
        const auto& b = expected[q];
        same_list = a.i == b.i && a.j == b.j && a.k == b.k;
    }
    check(same_list, "contact_pairs == подпоследовательность real_pairs с body[i] != body[j], k -- индекс в real_pairs");
    bool k_ok = true;
    for (const auto& cp : search.contact_pairs)
        k_ok = k_ok && cp.k < search.real_pairs.size() && search.real_pairs[cp.k].first == cp.i &&
               search.real_pairs[cp.k].second == cp.j;
    check(k_ok, "k каждой контактной пары указывает на ту же пару в real_pairs");

    // Эталон: последовательный проход по ВСЕМ real_pairs (диагностика сама
    // отфильтрует пары одного тела по body), CpuExecutor(1) -- без OpenMP-порций.
    std::vector<interactions::SolidContactPair> all;
    all.reserve(search.real_pairs.size());
    for (std::size_t k = 0; k < search.real_pairs.size(); ++k)
        all.push_back({search.real_pairs[k].first, search.real_pairs[k].second, static_cast<std::uint64_t>(k)});
    const double dt = 1.0e-7;
    const interactions::SolidContactDiagnostics ref = interactions::compute_solid_contact_diagnostics(
        rhs::CpuExecutor(1), cfg, cfg.mats, sys, derived.c, dt, all);
    note("эталон: pair_count = " + std::to_string(ref.pair_count) + ", max_delta = " +
         std::to_string(ref.max_delta) + ", min_dt = " + std::to_string(ref.min_dt));
    check(ref.pair_count > 0, "в сцене есть контактные пары с проникновением (pair_count > 0)");
    check(ref.max_delta > 0.0, "max_delta > 0");

    for (int nt : {1, 2, 7})
    {
        const interactions::SolidContactDiagnostics got = interactions::compute_solid_contact_diagnostics(
            rhs::CpuExecutor(nt), cfg, cfg.mats, sys, derived.c, dt, search.contact_pairs);
        check(got.pair_count == ref.pair_count && got.max_delta == ref.max_delta && got.min_dt == ref.min_dt,
              "диагностика по contact_pairs при CpuExecutor(" + std::to_string(nt) +
                  ") == эталон по всем real_pairs (pair_count, max_delta, min_dt)");
    }
}

// --- Параллельная сборка примитивов: эталон и память -------------------------
// Куб side^3 на тело, два тела рядом по x (непрерывная решётка 2*side x side x
// side); больше kSide, чтобы блоков сборки (n / 2048 пар) было несколько.
ParticleSystem make_lattice(const SimulationConfig& cfg, std::size_t side)
{
    const std::size_t cube = side * side * side;
    ParticleSystem sys;
    sys.resize(2 * cube);
    for (std::size_t i = 0; i < 2 * cube; ++i) sys.inzializ_index(i);
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t ix = 0; ix < side; ++ix)
            for (std::size_t iy = 0; iy < side; ++iy)
                for (std::size_t iz = 0; iz < side; ++iz)
                {
                    const std::size_t k = c * cube + (ix * side + iy) * side + iz;
                    sys.x[0][k] = (static_cast<double>(c * side + ix) + 0.5) * kDx;
                    sys.x[1][k] = (static_cast<double>(iy) + 0.5) * kDx;
                    sys.x[2][k] = (static_cast<double>(iz) + 0.5) * kDx;
                    sys.itype[k] = kItypeAl;
                    sys.body[k] = static_cast<int>(c) + 1;
                    sys.rho[k] = kRhoAl;
                    sys.rho_0[k] = kRhoAl;
                    sys.mass[k] = kRhoAl * kDx * kDx * kDx;
                    sys.hsml[k] = 2.0 * kDx;
                }
    for (std::size_t i = 0; i < 2 * cube; ++i) urs(cfg.mats, sys, i);
    return sys;
}

struct SolidScene
{
    interactions::SolidDerived derived;
    interactions::GhostSet ghosts;
    interactions::SolidExtended ext;
};

void make_scene(const SimulationConfig& cfg, const ParticleSystem& sys, std::size_t n, SolidScene& sc)
{
    rhs::CpuExecutor exec(1);
    interactions::compute_solid_derived(exec, cfg.mats, sys, n, sc.derived);
    interactions::build_ghost_set(sys, n, sc.ghosts);
    interactions::build_solid_extended(sys, n, sc.derived, sc.ghosts, sc.ext);
}

// Дословная копия прежнего ядра (src/kernel.cpp @ 06e0b99, до переноса в
// kernel.hpp): эталон сборки не зависит от нового inline-ядра, так что
// test_parallel_build проверяет и дословность kernel_cubic/gauss/quintic.
static void kernel_06e0b99(int skf, double r, std::span<const double> dx, double mhsml, double& w, double* dwdx)
{
    // Инициализация выходных значений
    w = 0.0;

    for (size_t d = 0; d < dim; d++)
    {
        dwdx[d] = 0.0;
    }
    // Вычисление относительного расстояния
    double q = r / mhsml;
    double factor = 0.0;

    // Выбор ядра в зависимости от skf
    if (skf == 1)    // Кубический сплайн (Monaghan 1985)
    {
        // Определение нормировочного фактора в зависимости от размерности
        if (dim == 1)
        {
            factor = 1.0 / mhsml;
        }
        else if (dim == 2)
        {
            factor = 15.0 / (7.0 * pi * mhsml * mhsml);
        }
        else if (dim == 3)
        {
            factor = 3.0 / (2.0 * pi * mhsml * mhsml * mhsml);
        }
        else
        {
            throw std::runtime_error("Неверная размерность: dim = " + std::to_string(dim));
        }

        // Вычисление ядра и его производной
        if (q >= 0.0 && q <= 1.0)
        {
            w = factor * (2.0/3.0 - q*q + q*q*q / 2.0);
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = factor * (-2.0 + 1.5 * q) * dx[d] / (mhsml * mhsml);
            }
        }
        else if (q > 1.0 && q <= 2.0)
        {
            double t = 2.0 - q;
            w = factor * (1.0/6.0) * ((2.0 - q)*(2.0 - q)*(2.0 - q));
            //w = factor * (1.0/6.0) * t * t * t;
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = -factor * 0.5 * ((2.0 - q)*(2.0 - q)) * dx[d] / (r * mhsml);
                //dwdx[d] = -factor * 0.5 * t * t * dx[d] / (r * mhsml);
            }
        }
        // Иначе w и dwdx остаются нулевыми

    }
    else if (skf == 2)      // Гауссово ядро
    {
        factor = 1.0 / (std::pow(mhsml, dim) * std::pow(pi, dim/2.0));

        if (q >= 0.0 && q <= 3.0)
        {
            w = factor * std::exp(-q*q);
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = w * (-2.0 * dx[d]) / (mhsml * mhsml);
            }
        }
        // Иначе w и dwdx остаются нулевыми

    }
    else if (skf == 3)      // Квинтическое ядро (Morris 1997)
    {
        if (dim == 1)
        {
            factor = 1.0 / (120.0 * mhsml);
        }
        else if (dim == 2)
        {
            factor = 7.0 / (478.0 * pi * mhsml * mhsml);
        }
        else if (dim == 3)
        {
            factor = 1.0 / (120.0 * pi * mhsml * mhsml * mhsml);
        }
        else
        {
            throw std::runtime_error("Неверная размерность: dim = " + std::to_string(dim));
        }

        if (q >= 0.0 && q <= 1.0)
        {
            w = factor * (((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)) - 6.0 * ((2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q)) + 15.0 * ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)));
            for (size_t d = 0; d < dim; ++d)
            {
                dwdx[d] = factor * (-120.0 + 120.0*q - 50.0*q*q) * dx[d] / (mhsml * mhsml);
            }
        }
        else if (q > 1.0 && q <= 2.0)
        {
            w = factor * (((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)) - 6.0 * ((2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q)));
            for (size_t d = 0; d < dim; ++d)
            {
                if (r > 1e-8)
                {
                    dwdx[d] = factor * (-5.0 * ((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)) + 30.0 * ((2.0 - q)*(2.0 - q)*(2.0 - q)*(2.0 - q))) * dx[d] / (r * mhsml);
                }
            }
        }
        else if (q > 2.0 && q <= 3.0)
        {
            w = factor * ((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q));
            for (size_t d = 0; d < dim; ++d)
            {
                if (r > 1e-8)
                {
                    dwdx[d] = -factor * 5.0 * ((3.0 - q)*(3.0 - q)*(3.0 - q)*(3.0 - q)) * dx[d] / (r * mhsml);
                }
            }
        }
        // Иначе w и dwdx остаются нулевыми
    }
}

// Эталон: прежний последовательный цикл сборки (build_solid_sph_search до
// параллельной версии), от уже найденных пар, с прежним ядром.
struct RefBuild
{
    std::vector<interactions::SolidSphPrimitive> primitives;
    std::vector<std::pair<rhs::index_t, rhs::index_t>> real_pairs;
    std::vector<interactions::SolidContactPair> contact_pairs;
};

RefBuild reference_build(const SimulationConfig& cfg, const interactions::SolidExtended& ext,
                         const std::vector<std::pair<rhs::index_t, rhs::index_t>>& pairs)
{
    RefBuild out;
    const std::size_t nreal = ext.nreal;
    const bool same_body_only = cfg.sph_same_body_only;
    auto add_real_pair = [&](rhs::index_t i, rhs::index_t j) {
        if (ext.body[i] != ext.body[j])
            out.contact_pairs.push_back({i, j, static_cast<std::uint64_t>(out.real_pairs.size())});
        out.real_pairs.emplace_back(i, j);
    };
    for (const auto& pr : pairs)
    {
        const rhs::index_t i = pr.first, j = pr.second;
        if (same_body_only && j < nreal && ext.body[i] != ext.body[j])
        {
            add_real_pair(i, j);
            continue;
        }
        double dx[dim];
        double r2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d) { dx[d] = ext.x[d][i] - ext.x[d][j]; r2 += dx[d] * dx[d]; }
        const double r = std::sqrt(r2);
        const double mh = 0.5 * (ext.hsml[i] + ext.hsml[j]);

        interactions::SolidSphPrimitive prim;
        prim.i = i; prim.j = j;
        double w = 0.0;
        double dwdx_buf[dim];
        kernel_06e0b99(cfg.skf, r, std::span<const double>(dx, dim), mh, w, dwdx_buf);
        prim.w = w;
        for (std::size_t d = 0; d < dim; ++d) prim.dwdx[d] = dwdx_buf[d];
        out.primitives.push_back(prim);
        if (j < nreal) add_real_pair(i, j);
    }
    return out;
}

// У SolidSphPrimitive нет паддинга (size_t, size_t, double, double[dim]):
// побайтное сравнение по всему элементу корректно.
static_assert(sizeof(interactions::SolidSphPrimitive) ==
                  2 * sizeof(rhs::index_t) + (1 + dim) * sizeof(double),
              "SolidSphPrimitive с паддингом: memcmp по элементу недопустим");

bool same_prims(const std::vector<interactions::SolidSphPrimitive>& a,
                const std::vector<interactions::SolidSphPrimitive>& b)
{
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(a[0])) == 0);
}

bool same_contact(const std::vector<interactions::SolidContactPair>& a,
                  const std::vector<interactions::SolidContactPair>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t q = 0; q < a.size(); ++q)
        if (a[q].i != b[q].i || a[q].j != b[q].j || a[q].k != b[q].k) return false;
    return true;
}

bool equals_reference(const interactions::SolidSphSearchResult& got, const RefBuild& ref)
{
    return same_prims(got.primitives, ref.primitives) && got.real_pairs == ref.real_pairs &&
           same_contact(got.contact_pairs, ref.contact_pairs);
}

void test_parallel_build()
{
    std::printf("test_parallel_build:\n");
    constexpr std::size_t kSideBig = 8;
    for (int skf : {0, 1, 2, 3}) // 0 -- неизвестный skf: ядро даёт нули (ветка default)
        for (bool same_body : {false, true})
        {
            SimulationConfig cfg = make_config(same_body, 1);
            cfg.skf = skf;
            const ParticleSystem sys = make_lattice(cfg, kSideBig);
            const std::size_t n = 2 * kSideBig * kSideBig * kSideBig;
            SolidScene sc;
            make_scene(cfg, sys, n, sc);
            for (int nt : {1, 2, 7})
            {
                rhs::CpuExecutor exec(nt);
                interactions::SolidSphSearchResult search;
                interactions::build_solid_sph_search(exec, cfg, sc.ext, search);
                const RefBuild ref = reference_build(cfg, sc.ext, search.pairs);
                const std::string tag = "skf=" + std::to_string(skf) + ", same_body_only=" +
                                        std::to_string(same_body) + ", потоков " + std::to_string(nt);
                if (skf == 1 && nt == 1)
                    note("пар: " + std::to_string(search.pairs.size()) + ", primitives: " +
                         std::to_string(ref.primitives.size()) + ", real_pairs: " +
                         std::to_string(ref.real_pairs.size()) + ", contact_pairs: " +
                         std::to_string(ref.contact_pairs.size()));
                check(!ref.primitives.empty() && !ref.contact_pairs.empty(), "эталон непуст: " + tag);
                check(equals_reference(search, ref), "primitives/real_pairs/contact_pairs побайтно == эталон: " + tag);
            }
        }
}

void test_build_memory()
{
    std::printf("test_build_memory:\n");
    SimulationConfig cfg = make_config(true, 1);
    rhs::CpuExecutor exec(7);
    interactions::SolidSphSearchResult search;

    auto run = [&](std::size_t side, const std::string& tag) {
        const ParticleSystem sys = make_lattice(cfg, side);
        SolidScene sc;
        make_scene(cfg, sys, 2 * side * side * side, sc);
        interactions::build_solid_sph_search(exec, cfg, sc.ext, search);
        const RefBuild ref = reference_build(cfg, sc.ext, search.pairs);
        check(equals_reference(search, ref), "результат == эталон: " + tag);
    };
    struct Snap
    {
        const void* p[4];
        std::size_t c[4];
        bool operator==(const Snap&) const = default;
    };
    auto snap = [&]() {
        return Snap{{search.primitives.data(), search.real_pairs.data(), search.contact_pairs.data(),
                     search.blocks.data()},
                    {search.primitives.capacity(), search.real_pairs.capacity(),
                     search.contact_pairs.capacity(), search.blocks.capacity()}};
    };

    run(8, "сцена 8, вызов 1");
    Snap s0 = snap();
    run(8, "сцена 8, вызов 2");
    check(snap() == s0, "повторный вызов: data()/capacity() primitives, real_pairs, contact_pairs, blocks не изменились");
    run(8, "сцена 8, вызов 3");
    check(snap() == s0, "третий вызов: data()/capacity() не изменились");

    run(10, "сцена 10 (больше)");
    s0 = snap();
    run(6, "сцена 6 (меньше)");
    check(snap() == s0, "меньшая сцена: data()/capacity() не изменились");
    run(8, "сцена 8 (снова)");
    check(snap() == s0, "сцена 8 после максимума: data()/capacity() не изменились");
    run(10, "сцена 10 (снова максимум)");
    check(snap() == s0, "снова максимум: data()/capacity() не изменились");

    check(!search.blocks.empty() && reinterpret_cast<std::uintptr_t>(search.blocks.data()) % 64 == 0,
          "blocks.data() выровнен по 64 байта");
    check(search.blocks.size() > 2, "в сцене несколько блоков сборки (blocks.size() = " +
                                       std::to_string(search.blocks.size()) + ")");
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
    test_contact_pairs();
    test_parallel_build();
    test_build_memory();
    test_exception();
    test_integration();

    if (e2e::failures == 0) std::printf("\nВСЕ ПРОВЕРКИ ПРОЙДЕНЫ\n");
    else std::printf("\nПРОВАЛЕНО проверок: %d\n", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
