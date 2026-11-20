// test_regression.cpp
// Регрессионный тест конвейера на кейсах BC2 и BC6.
//
// ЗАЧЕМ. Сквозные тесты test_e2e_bc2/bc6 проверяют живучесть конвейера и
// сверяют форму кривой с ПРИБОРНЫМ эталоном в грубых воротах (амплитуда в
// пределах ×2). Ошибка, меняющая, скажем, осевой член импульса на 5 %, ими не
// ловится. План нитей (notes/fiber-model/fiber-plan.md) требует на этапах Ф2–Ф8, чтобы «BC2/BC6
// не изменились», прежде всего на Ф4 — переносе Solid на новую схему RHS.
// Этот тест делает такое утверждение проверяемым.
//
// ЧТО ДЕЛАЕТ. Кейс (BC2 или BC6, smoke-разрешение, заданный
// contact_force_type) считается с нуля ровно N шагов — останов по пределу
// шагов сторожа, а не по t_total, чтобы при чуть разошедшемся dt сравнивались
// состояния после одинакового числа шагов. Затем все поля всех частиц
// снимаются прямо из памяти и сверяются с эталоном из
// cases/<кейс>/regression/snapshots/ по относительным нормам с допуском
// (tests/regression_snapshot.hpp).
//
// ДВЕ КОНТРОЛЬНЫЕ ТОЧКИ, два независимых прогона:
//   N = 1  — одно вычисление RHS и один шаг предиктор-корректора. Это «стенд
//            эквивалентности RHS на одном состоянии» из рисков Ф4
//            (notes/fiber-model/fiber-plan.md §5). Допуск тесный.
//   N = K  — накопление по времени: контакт, JH-2, Мизес, h, XSPH. Допуск
//            шире: расхождения на уровне округления растут от шага к шагу.
//
// ЭТАЛОН генерируется ЭТИМ ЖЕ тестом (--emit-golden) при OMP_NUM_THREADS=1.
// Готовая цель: cmake --build build --target regression_golden. Ни один
// зарегистрированный ctest-тест в дерево исходников не пишет.
//
// ВАЖНО ПРО РАЗМЕРНОСТЬ: программа собирается дважды и компилирует свои копии
// src/*.cpp (решение E1, ADR-002):
//   test_regression    -DSPH_DIM=2  кейсы bc2, bc6 (осесимметрия), K = 100;
//   test_regression_3d -DSPH_DIM=3  кейс bc2mini3d (уменьшенный BC2 в 3D,
//                                   tests/bc2_mini3d_geometry.hpp), K = 50.
// Состав полей эталона зависит от размерности (tests/regression_snapshot.hpp),
// размерность пишется в шапку эталона, и эталон чужой размерности не читается.

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"
#include "run_guard.hpp"

#include "e2e_harness.hpp"
#include "regression_cases.hpp"
#include "regression_snapshot.hpp"

// Одна и та же программа собирается дважды: test_regression (SPH_DIM=2,
// кейсы BC2 и BC6 в осесимметричной постановке) и test_regression_3d
// (SPH_DIM=3, уменьшенный BC2 в полном 3D, cases/bc2/regression/bc2_mini3d_geometry.hpp).
static_assert(dim == 2 || dim == 3, "test_regression собирается с -DSPH_DIM=2 или 3");
static_assert(dim != 2 || axis_sim, "кейсы BC2/BC6 требуют axis_sim = true");

#ifndef SPH_CASES_DIR
// Задаётся из CMake (функция sph_case). Фолбэк -- запуск из корня репозитория.
#define SPH_CASES_DIR "cases"
#endif

namespace fs = std::filesystem;

using e2e::check;
using e2e::note;

namespace
{

// ===========================================================================
// 1. Контрольные точки и допуски
// ===========================================================================

constexpr size_t kStepsShort = 1;
// 2D: K = 100 — наибольшее K <= 100, при котором разброс между числами
// потоков остаётся <= 1e-6 по всем полям (критерий калибровки). Фактически
// он на этом горизонте не растёт вовсе: от N = 1 до N = 100 — те же ~1e-13.
// 3D: K = 50 — бюджет времени. Шаг в 3D стоит ~0.18 с (3648 частиц), и
// четыре типа контакта по N = 1 и N = 50 укладываются примерно в минуту.
constexpr size_t kStepsLong  = (dim == 3) ? 50 : 100;

// Допуски. Нормы — в шапке tests/regression_snapshot.hpp.
//
// КАЛИБРОВКА (2026-09-24, Windows 11, clang++ из LLVM, Release -O3, 12
// потоков). Однопоточный эталон против прогонов при OMP_NUM_THREADS = 1, 2,
// 4, 12 (по 3 повтора), все 8 кейсов, N = 1, 25, 50, 100. Допуск
// tol = max(100·s, 1e-11), округлённый вверх до ряда 1-2-5, где s — худший
// наблюдаемый разброс поля по всем кейсам. Пол 1e-11 перекрывает ошибку
// представления %.12e в файле эталона (~5e-13).
//
// Что показала калибровка:
//   * s почти целиком — ошибка представления в файле: однопоточный прогон
//     против однопоточного эталона даёт те же 1e-13..5e-13, что и
//     многопоточный. Собственный шум потоков виден только в e_pl, D, dt;
//   * худшие s: N = 1 — hsml 4.6e-13; N = 100 — D 1.4e-12 (Linf, частица у
//     порога повреждения JH-2), e_pl 9.3e-13, v1 4.0e-13;
//   * роста с числом шагов практически нет — горизонт K = 100 не хаотичен.
//
// Допуск не расширять вслепую. Правка, которая лишь переставляет слагаемые,
// даёт шум того же порядка, что потоки, и проходит. Если такая правка всё же
// не прошла, это повод искать ошибку, а не поднимать число. Переносить
// допуск на другую машину или компилятор можно только после перекалибровки.
#if SPH_DIM == 2
//                             x0     x1     v0     v1     rho    p      u      c      hsml   e_pl   D      sxx    syy    szz    sxy
constexpr std::array<double, reg::kFieldCount> kInfShort = {2e-11, 2e-11, 5e-11, 5e-11, 2e-11, 1e-11, 1e-11, 1e-11, 5e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 2e-11};
constexpr std::array<double, reg::kFieldCount> kL2Short  = {1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 5e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 2e-11};
constexpr std::array<double, reg::kFieldCount> kInfLong  = {2e-11, 2e-11, 5e-11, 5e-11, 2e-11, 5e-11, 5e-11, 1e-11, 5e-11, 1e-10, 2e-10, 5e-11, 5e-11, 5e-11, 5e-11};
constexpr std::array<double, reg::kFieldCount> kL2Long   = {2e-11, 2e-11, 2e-11, 2e-11, 2e-11, 2e-11, 5e-11, 1e-11, 5e-11, 2e-11, 1e-11, 2e-11, 2e-11, 2e-11, 2e-11};
#elif SPH_DIM == 3
// 3D, уменьшенный BC2 (КАЛИБРОВКА 2026-09-25, та же машина и сборка, та же
// процедура: OMP_NUM_THREADS = 1, 2, 4, 12, по 3 повтора, 4 типа контакта,
// N = 1, 25, 50; код уже с исправленным знаком syz·ryz в dszz).
// Картина та же, что в 2D: s почти целиком — ошибка представления %.12e
// (однопоточный прогон даёт те же числа); шум потоков заметен только в
// e_pl, v2 и dt. Худшие s: N = 1 — hsml 1.7e-13 при v2 2.2e-13; N = 50 —
// hsml 4.8e-13, e_pl 3.6e-13, p 3.3e-13. Роста с числом шагов нет.
//                             x0     x1     x2     v0     v1     v2     rho    p      u      c      hsml   e_pl   D      sxx    syy    szz    sxy    syz    sxz
constexpr std::array<double, reg::kFieldCount> kInfShort = {2e-11, 2e-11, 2e-11, 1e-11, 1e-11, 5e-11, 2e-11, 1e-11, 2e-11, 1e-11, 2e-11, 1e-11, 1e-11, 2e-11, 2e-11, 1e-11, 1e-11, 2e-11, 2e-11};
constexpr std::array<double, reg::kFieldCount> kL2Short  = {1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 2e-11, 1e-11, 2e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 2e-11, 2e-11};
constexpr std::array<double, reg::kFieldCount> kInfLong  = {2e-11, 2e-11, 2e-11, 1e-11, 1e-11, 5e-11, 2e-11, 5e-11, 1e-11, 1e-11, 5e-11, 5e-11, 1e-11, 5e-11, 5e-11, 2e-11, 1e-11, 5e-11, 5e-11};
constexpr std::array<double, reg::kFieldCount> kL2Long   = {2e-11, 2e-11, 1e-11, 1e-11, 1e-11, 2e-11, 1e-11, 2e-11, 1e-11, 1e-11, 5e-11, 1e-11, 1e-11, 5e-11, 5e-11, 2e-11, 2e-11, 2e-11, 2e-11};
#endif

// Для произвольного --steps берётся таблица N = K: калибровка показала, что
// разброс от числа шагов почти не зависит.
reg::Tolerances tolerances_for(size_t steps)
{
    reg::Tolerances t;
    const bool shortrun = (steps <= kStepsShort);
    for (size_t k = 0; k < reg::kFieldCount; ++k)
    {
        t.field[k] = shortrun ? reg::Tolerance{kInfShort[k], kL2Short[k]}
                              : reg::Tolerance{kInfLong[k], kL2Long[k]};
    }
    // t и dt: наблюдаемый разброс не больше 6.3e-15 (dt, N = 100) — пол.
    t.time_rel = 1e-11;
    t.dt_rel   = 1e-11;
    return t;
}


std::string golden_name(reg::Case c, int contact, size_t steps)
{
    return std::string(reg::case_name(c)) + "_contact" + std::to_string(contact) + "_step" +
           std::to_string(steps) + ".txt";
}

int omp_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// ===========================================================================
// 2. Прогон
// ===========================================================================

struct RunSettings
{
    reg::Case kase = reg::Case::BC2;
    int    contact = 0;
    size_t steps = kStepsShort;
    fs::path config;
    bool   silence = true;
    bool   keep_workdir = false;
    double wall_limit_s = 300.0;
};

struct RunResult
{
    bool ok = false;
    std::string error;
    reg::Snapshot snap;
    double elapsed_s = 0.0;
    std::string workdir;
};

RunResult run_steps(const RunSettings& rs)
{
    RunResult res;
    e2e::ScopedWorkDir wd("sph_regression", std::string(reg::case_name(rs.kase)) + "_c" +
                                                std::to_string(rs.contact) + "_s" +
                                                std::to_string(rs.steps));
    if (rs.keep_workdir) wd.keep();
    res.workdir = wd.dir().string();

    SimulationConfig cfg;
    ParticleSystem sys;
    reg::Prepared prep;
    try
    {
        prep = reg::prepare_case(rs.kase, rs.contact, rs.config, cfg, sys);
    }
    catch (const std::exception& e)
    {
        res.error = std::string("подготовка прогона не удалась: ") + e.what();
        return res;
    }
    const std::vector<double> mass0(sys.mass.begin(),
                                    sys.mass.begin() + static_cast<std::ptrdiff_t>(prep.ntotal));

    // Детектор застоя выключен (0): его оценка «сколько шагов осталось до
    // t_total» здесь бессмысленна — прогон заведомо обрывается раньше.
    RunGuard guard(rs.wall_limit_s, 0.0, 1);
    guard.set_step_limit(rs.steps);

    double dt = prep.dt;
    double time = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    try
    {
        e2e::ScopedSilence silence(rs.silence);
        time_integration(cfg, sys, prep.ntotal, dt, time, prep.delta_x, &guard);
    }
    catch (const std::exception& e)
    {
        res.error = std::string("исключение в конвейере: ") + e.what() + " (на шаге " +
                    std::to_string(guard.steps()) + ")";
        return res;
    }
    res.elapsed_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    if (guard.reason() != RunGuard::Reason::StepLimit)
    {
        std::ostringstream os;
        os << "прогон не дошёл до предела шагов: сделано " << guard.steps() << " из "
           << rs.steps;
        if (guard.failed()) os << "\n" << guard.message();
        else os << " (достигнут t_total = " << cfg.t_total << " — увеличьте его в конфиге)";
        res.error = os.str();
        return res;
    }

    res.snap = reg::capture(sys, prep.ntotal);
    res.snap.case_name = reg::case_name(rs.kase);
    res.snap.contact_force_type = rs.contact;
    res.snap.delta_x = prep.delta_x;
    res.snap.steps = guard.steps();
    res.snap.time = time;
    res.snap.dt = guard.dt_last();
    res.snap.threads = omp_threads();

    // Инварианты, которых нет в эталоне: масса постоянна, а при dim = 2 syz и
    // sxz тождественно нулевые (при dim = 3 они входят в эталон как поля).
    // Нарушение — дефект, а не расхождение.
    size_t mass_changed = 0, shear_nonzero = 0;
    for (size_t i = 0; i < prep.ntotal; ++i)
    {
        if (sys.mass[i] != mass0[i]) ++mass_changed;
        if (dim == 2 && (sys.syz[i] != 0.0 || sys.sxz[i] != 0.0)) ++shear_nonzero;
    }
    check(mass_changed == 0, "масса частиц не изменилась за прогон (" +
                                 std::to_string(mass_changed) + " изменённых)");
    if (dim == 2)
    {
        check(shear_nonzero == 0, "syz и sxz нулевые при dim = 2 (" +
                                      std::to_string(shear_nonzero) + " ненулевых)");
    }

    res.ok = true;
    return res;
}

// ===========================================================================
// 3. Самопроверка компаратора
// ===========================================================================
// Без расчёта: синтетический снимок с масштабами реального кейса. Доказывает,
// что сверка (а) пропускает эталон против самого себя и после записи/чтения
// файла, (б) ловит возмущение выше допуска и называет именно ту частицу и то
// поле, (в) не ловит возмущение ниже допуска, (г) ловит NaN, расхождение
// itype, ntotal и числа шагов. Тест, который не падает на заведомой поломке,
// регрессию не охраняет.

reg::Snapshot synthetic_snapshot(size_t n)
{
    reg::Snapshot s;
    s.case_name = "bc2";
    s.contact_force_type = 0;
    s.delta_x = 8e-4;
    s.steps = kStepsLong;
    s.time = 4.4e-7;
    s.dt = 8.7621e-9;
    s.threads = 1;
    // Характерные масштабы полей кейса BC2, СИ.
    std::array<double, reg::kFieldCount> amp{};
    for (size_t d = 0; d < dim; ++d)
    {
        amp[reg::kX0 + d] = 0.03 + 0.01 * static_cast<double>(d);
        amp[reg::kV0 + d] = (d + 1 == dim) ? 2210.0 : 50.0;
    }
    amp[reg::kRho] = 2517.0; amp[reg::kP] = 3e9; amp[reg::kU] = 1e5; amp[reg::kC] = 1.3e4;
    amp[reg::kHsml] = 1.6e-3; amp[reg::kEpl] = 0.05; amp[reg::kD] = 0.3;
    for (size_t k = reg::kSxx; k < reg::kFieldCount; ++k)
        amp[k] = 1e9 * (1.0 + 0.5 * static_cast<double>(k - reg::kSxx));
    s.itype.resize(n);
    for (size_t k = 0; k < reg::kFieldCount; ++k) s.f[k].resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        s.itype[i] = (i < n / 2) ? 8 : 10;
        for (size_t k = 0; k < reg::kFieldCount; ++k)
        {
            const double phase = 0.37 * static_cast<double>(i) + 1.3 * static_cast<double>(k);
            s.f[k][i] = amp[k] * (0.6 + 0.4 * std::sin(phase));
        }
    }
    return s;
}

void comparator_selfcheck()
{
    std::printf("\n[самопроверка] Компаратор снимков\n");
    const size_t n = 64;
    const size_t victim = 17;
    const reg::Snapshot ref = synthetic_snapshot(n);

    for (const size_t steps : {kStepsShort, kStepsLong})
    {
        const reg::Tolerances tol = tolerances_for(steps);
        const std::string tag = " [допуски N = " + std::to_string(steps) + "]";

        check(reg::compare(ref, ref, tol).ok(), "эталон против самого себя проходит" + tag);

        for (const size_t k : {size_t(reg::kP), size_t(reg::kV1), size_t(reg::kD)})
        {
            double max_ref = 0.0;
            for (double v : ref.f[k]) max_ref = std::max(max_ref, std::fabs(v));
            const double scale = std::max(max_ref, reg::kFieldFloor[k]);
            const std::string fname = reg::kFieldName[k];

            reg::Snapshot big = ref;
            big.f[k][victim] += 10.0 * tol.field[k].inf * scale;
            const reg::Report rb = reg::compare(ref, big, tol);
            bool others_ok = true;
            for (size_t j = 0; j < reg::kFieldCount; ++j)
                if (j != k) others_ok = others_ok && rb.field[j].ok;
            check(!rb.ok() && !rb.field[k].ok && rb.field[k].worst == victim && others_ok,
                  "возмущение 10·допуск в поле " + fname + " частицы " +
                      std::to_string(victim) + " ловится и указывает на неё" + tag);

            reg::Snapshot small = ref;
            small.f[k][victim] += 0.1 * tol.field[k].inf * scale;
            check(reg::compare(ref, small, tol).ok(),
                  "возмущение 0,1·допуск в поле " + fname + " не ловится" + tag);
        }
    }

    const reg::Tolerances tol = tolerances_for(kStepsShort);
    {
        reg::Snapshot bad = ref;
        bad.f[reg::kV1][victim] = std::nan("");
        const reg::Report r = reg::compare(ref, bad, tol);
        check(!r.ok() && !r.field[reg::kV1].ok, "NaN в расчёте ловится");
    }
    {
        reg::Snapshot bad = ref;
        bad.itype[victim] = 5;
        const reg::Report r = reg::compare(ref, bad, tol);
        check(!r.comparable && r.first_itype_mismatch == victim,
              "расхождение itype делает снимки несравнимыми и называет частицу");
    }
    {
        reg::Snapshot bad = ref;
        bad.itype.pop_back();
        for (auto& v : bad.f) v.pop_back();
        check(!reg::compare(ref, bad, tol).comparable, "расхождение ntotal ловится");
    }
    {
        reg::Snapshot bad = ref;
        bad.steps += 1;
        check(!reg::compare(ref, bad, tol).ok(), "расхождение числа шагов ловится");
    }

    // Запись и чтение: формат теряет не больше ошибки представления %.12e,
    // то есть сверка с самым тесным допуском после круга обязана пройти.
    e2e::ScopedWorkDir wd("sph_regression", "selfcheck");
    const std::string path = (wd.dir() / "roundtrip.txt").string();
    std::string err;
    reg::Snapshot back;
    const bool io_ok = reg::write_snapshot(ref, path, err) && reg::read_snapshot(path, back, err);
    check(io_ok, "эталон записывается и читается" + (err.empty() ? "" : " (" + err + ")"));
    if (io_ok)
    {
        check(back.case_name == ref.case_name && back.steps == ref.steps &&
                  back.contact_force_type == ref.contact_force_type &&
                  back.time == ref.time && back.dt == ref.dt && back.delta_x == ref.delta_x,
              "метаданные переживают запись и чтение точно");
        check(reg::compare(ref, back, tol).ok(),
              "поля после записи и чтения проходят сверку с тесным допуском");
    }
}

// ===========================================================================
// 4. Режимы
// ===========================================================================

struct Options
{
    bool selfcheck_only = false;
    bool have_case = false;
    reg::Case kase = reg::Case::BC2;
    int contact = 0;
    std::vector<size_t> steps{kStepsShort, kStepsLong};
    fs::path emit_dir;
    fs::path golden_dir;   // пусто -- cases/<кейс>/regression/snapshots (reg::default_golden_dir)
    fs::path config;
    std::vector<std::pair<std::string, std::string>> meta;
    bool verbose = false;
    bool keep_workdir = false;
    double wall_limit_s = 300.0;
};

void run_checkpoint(const Options& o, size_t steps)
{
    const std::string name = golden_name(o.kase, o.contact, steps);
    std::printf("\n[N = %zu] %s\n", steps, name.c_str());

    RunSettings rs;
    rs.kase = o.kase;
    rs.contact = o.contact;
    rs.steps = steps;
    rs.config = o.config;
    rs.silence = !o.verbose;
    rs.keep_workdir = o.keep_workdir;
    rs.wall_limit_s = o.wall_limit_s;

    const RunResult r = run_steps(rs);
    if (o.keep_workdir) note("каталог прогона сохранён: " + r.workdir);
    if (!r.ok)
    {
        check(false, "прогон не состоялся: " + r.error);
        return;
    }
    {
        std::ostringstream os;
        os << "ntotal = " << r.snap.size() << ", шагов = " << r.snap.steps << ", t = "
           << r.snap.time << ", dt = " << r.snap.dt << ", потоков OpenMP = " << r.snap.threads
           << ", " << e2e::ru_num(r.elapsed_s, 2) << " с";
        note(os.str());
    }

    if (!o.emit_dir.empty())
    {
        reg::Snapshot s = r.snap;
        for (const auto& [k, v] : o.meta) s.extra[k] = v;
        std::error_code ec;
        fs::create_directories(o.emit_dir, ec);
        const fs::path out = o.emit_dir / name;
        std::string err;
        check(reg::write_snapshot(s, out.string(), err),
              "эталон записан: " + out.string() + (err.empty() ? "" : " (" + err + ")"));
        return;
    }

    const fs::path gpath = o.golden_dir / name;
    reg::Snapshot ref;
    std::string err;
    if (!reg::read_snapshot(gpath.string(), ref, err))
    {
        check(false, err + "\n         сгенерировать: cmake --build build --target regression_golden");
        return;
    }
    const reg::Tolerances tol = tolerances_for(steps);
    const reg::Report rep = reg::compare(ref, r.snap, tol);
    reg::print_report(rep, ref, r.snap, tol);
    check(rep.ok(), "состояние после " + std::to_string(steps) + " " +
                        e2e::ru_plural(steps, "шага", "шагов", "шагов") +
                        " совпадает с эталоном в пределах допусков");
}

void print_usage()
{
    std::printf(
        "Использование: test_regression --case=КЕЙС [--contact-type=N] [ключи]\n"
        "  --case=КЕЙС        bc2 | bc6 (сборка dim = 2), bc2mini3d (сборка dim = 3,\n"
        "                     test_regression_3d); сетка smoke сквозных тестов\n"
        "  --contact-type=N   contact_force_type, N = 0..3 (по умолчанию 0)\n"
        "  --steps=S          одна контрольная точка на S шагов вместо пары\n"
        "                     N = %zu и N = %zu\n"
        "  --emit-golden=DIR  записать эталоны в DIR вместо сверки (эталон\n"
        "                     снимается при OMP_NUM_THREADS=1; готовая цель:\n"
        "                     cmake --build build --target regression_golden)\n"
        "  --meta=КЛЮЧ=ЗНАЧ   строка шапки эталона (commit, date, ...), повторяемый\n"
        "  --golden-dir=DIR   откуда читать эталоны (по умолчанию\n"
        "                     cases/<кейс>/regression/snapshots)\n"
        "  --config=ФАЙЛ      другой конфиг решателя (проверка чувствительности:\n"
        "                     копия конфига с изменённым параметром обязана\n"
        "                     провалить сверку)\n"
        "  --selfcheck        только самопроверка компаратора, без расчёта\n"
        "  --verbose          не глушить вывод конвейера\n"
        "  --keep-workdir     не удалять временный каталог прогона\n"
        "  --time-limit=СЕК   бюджет времени на один прогон (по умолчанию 300)\n",
        kStepsShort, kStepsLong);
}

bool parse_args(int argc, char** argv, Options& o, int& exit_code)
{
    auto value = [](const std::string& a, const char* key) { return a.substr(std::strlen(key)); };
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        try
        {
            if (a == "--selfcheck") o.selfcheck_only = true;
            else if (a == "--verbose") o.verbose = true;
            else if (a == "--keep-workdir") o.keep_workdir = true;
            else if (a.rfind("--case=", 0) == 0)
            {
                const std::string v = value(a, "--case=");
                if (v == "bc2") o.kase = reg::Case::BC2;
                else if (v == "bc6") o.kase = reg::Case::BC6;
                else if (v == "bc2mini3d") o.kase = reg::Case::BC2_MINI3D;
                else throw std::invalid_argument("кейс должен быть bc2, bc6 или bc2mini3d");
                if (reg::case_dim(o.kase) != dim)
                {
                    throw std::invalid_argument("кейс определён при dim = " +
                                                std::to_string(reg::case_dim(o.kase)) +
                                                ", а программа собрана с dim = " +
                                                std::to_string(dim));
                }
                o.have_case = true;
            }
            else if (a.rfind("--contact-type=", 0) == 0)
            {
                o.contact = std::stoi(value(a, "--contact-type="));
                if (o.contact < 0 || o.contact > 3) throw std::invalid_argument("N = 0..3");
            }
            else if (a.rfind("--steps=", 0) == 0)
            {
                const long long s = std::stoll(value(a, "--steps="));
                if (s < 1) throw std::invalid_argument("S >= 1");
                o.steps = {static_cast<size_t>(s)};
            }
            else if (a.rfind("--emit-golden=", 0) == 0) o.emit_dir = value(a, "--emit-golden=");
            else if (a.rfind("--golden-dir=", 0) == 0) o.golden_dir = value(a, "--golden-dir=");
            else if (a.rfind("--config=", 0) == 0) o.config = value(a, "--config=");
            else if (a.rfind("--meta=", 0) == 0)
            {
                const std::string kv = value(a, "--meta=");
                const size_t eq = kv.find('=');
                if (eq == std::string::npos || eq == 0) throw std::invalid_argument("КЛЮЧ=ЗНАЧ");
                o.meta.emplace_back(kv.substr(0, eq), kv.substr(eq + 1));
            }
            else if (a.rfind("--time-limit=", 0) == 0)
            {
                o.wall_limit_s = std::max(0.0, std::stod(value(a, "--time-limit=")));
            }
            else if (a == "-h" || a == "--help")
            {
                print_usage();
                exit_code = 0;
                return false;
            }
            else throw std::invalid_argument("неизвестный аргумент");
        }
        catch (const std::exception& e)
        {
            std::fprintf(stderr, "Неверный аргумент %s: %s\n", a.c_str(), e.what());
            print_usage();
            exit_code = 2;
            return false;
        }
    }
    if (!o.selfcheck_only && !o.have_case)
    {
        std::fprintf(stderr, "Нужен --case=КЕЙС или --selfcheck\n");
        print_usage();
        exit_code = 2;
        return false;
    }
    if (o.config.empty()) o.config = reg::default_config(o.kase, SPH_CASES_DIR);
    if (o.golden_dir.empty()) o.golden_dir = reg::default_golden_dir(o.kase, SPH_CASES_DIR);
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    Options o;
    int exit_code = 0;
    if (!parse_args(argc, argv, o, exit_code)) return exit_code;

    const fs::path cwd = fs::current_path();
    const e2e::DirSnapshot before = e2e::snapshot_dir(cwd);

    if (o.selfcheck_only)
    {
        std::printf("=== регрессия: самопроверка компаратора ===\n");
    }
    else
    {
        std::printf("=== регрессия %s, contact_force_type = %d, %s ===\n",
                    reg::case_name(o.kase), o.contact,
                    o.emit_dir.empty() ? "сверка с эталоном" : "ГЕНЕРАЦИЯ эталона");
        note("конфиг: " + o.config.string());
        if (!o.emit_dir.empty() && omp_threads() != 1)
        {
            e2e::warn("эталон снимается при " + std::to_string(omp_threads()) +
                      " потоках OpenMP; воспроизводимый эталон — при OMP_NUM_THREADS=1");
        }
    }

    // Самопроверка идёт и в обычном режиме, до прогона: сверка, не способная
    // упасть, не должна выдавать зелёный результат.
    comparator_selfcheck();

    if (!o.selfcheck_only)
    {
        for (const size_t s : o.steps) run_checkpoint(o, s);
    }

    e2e::check_no_side_effects(cwd, before);

    std::printf("\n=== регрессия: %s (провалов: %d) ===\n",
                e2e::failures == 0 ? "ПРОЙДЕНА" : "ПРОВАЛЕНА", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
