#pragma once
// tests/clamp_integration_fiber.hpp
// Интеграционные тесты кинематического закрепления бусин Fiber через РЕАЛЬНЫЙ
// цикл time_integration(cfg, fiber::Model&, ...) (fiber/TZ_CLAMP_SOLID_FIBER.md,
// §3, §8, §10: F1-F5; Clamp/FIBER_CLAMP_SPEC.md, §8).
//   F1 -- закреплённые концы нити неподвижны побитово (PC и Euler), свободные
//         бусины движутся, кадры fiber_output без NaN/Inf;
//   F2 -- то же при трении нить-нить (mu_ff > 0) со скользящей по зажиму нитью;
//         пара полностью закреплённых нитей (c_k = 0) -- без NaN, реакция конечна;
//   F3 -- предварительно натянутая струна: период основной моды против
//         T = 2L*sqrt(mu/N), три уровня сгущения, энергия;
//   F4 -- баланс импульса: P^{n+1} - P^n = h*R_fiber (Fiber без Solid);
//   F5 -- разрыв связи на краю зажима: бусина остаётся закреплённой.
// Файл включается в tests/test_clamp.cpp ПОСЛЕ определения check()/failures/
// same_bits() и вызывается из main через run_clamp_fiber_integration_tests().
//
// Сцена без Solid: model.solid -- одна фиктивная частица (ntotal = 0): цикл
// пишет log.txt с sys.vx[dim-1][moni_particle], поэтому вектор не должен быть
// пустым. Все прогоны идут во временном каталоге (e2e::ScopedWorkDir): цикл
// пишет log.txt и fiber_output_*.txt в текущий каталог.
#include "fiber/fiber_output.hpp"
#include "fiber/model.hpp"
#include "interactions/clamp.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_friction.hpp"
#include "ode/ode_solver.hpp"
#include "rhs/cpu_executor.hpp"
#include "run_guard.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"

#include "e2e_harness.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace clamp_fiber_it
{

// ---------------------------------------------------------------------------
// Материал и сцена
// ---------------------------------------------------------------------------
// E = 1e7 Па, d = 1 мм, rho = 1000 кг/м^3: осевая скорость звука 100 м/с,
// EA = 7.85 Н. Шаг dt_limit связей ~ 0.8*l0/(sqrt(2)*100) ~ 6e-5 с при l0 = 1 см.
inline fiber::FiberMaterial make_mat(double mu_ff, double eps_break)
{
    fiber::FiberMaterial m;
    m.E = 1.0e7;
    m.nu = 0.3;
    m.d = 1.0e-3;
    m.rho = 1000.0;
    m.eps_break = eps_break;
    m.mu_ff = mu_ff;
    return m;
}

inline double mat_area(const fiber::FiberMaterial& m) { return interactions::fiber_cross_section_area(m); }

// Добавляет прямую нить из n бусин: x_k = origin + k*spacing*dir, связи (k,k+1)
// с l0, l_break = (1+eps_break)*l0; масса по ненагруженной конфигурации
// m = rho*A*l0 (у концов половина), как в src/fiber_input.cpp.
// Возвращает индекс первой бусины.
inline std::size_t add_chain(fiber::Model& model, std::size_t n, const double (&origin)[3], const double (&dir)[3],
                             double spacing, double l0, int material_id, int object_id)
{
    const std::size_t first = model.fiber.size();
    model.fiber.resize(first + n);
    model.topology.prev.resize(first + n, -1);
    model.topology.next.resize(first + n, -1);
    const fiber::FiberMaterial& mat = model.materials[static_cast<std::size_t>(material_id)];
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t g = first + k;
        for (std::size_t d = 0; d < 3; ++d)
        {
            model.fiber.x[d][g] = origin[d] + static_cast<double>(k) * spacing * dir[d];
            model.fiber.vx[d][g] = 0.0;
        }
        const double seg = (k == 0 || k + 1 == n) ? 0.5 : 1.0;
        model.fiber.mass[g] = mat.rho * mat_area(mat) * l0 * seg;
        model.fiber.object_id[g] = object_id;
        model.fiber.material_id[g] = material_id;
        if (k > 0) model.topology.prev[g] = static_cast<long long>(g - 1);
        if (k + 1 < n) model.topology.next[g] = static_cast<long long>(g + 1);
    }
    for (std::size_t k = 0; k + 1 < n; ++k)
    {
        fiber::FiberBond b;
        b.i = first + k;
        b.j = first + k + 1;
        b.l0 = l0;
        b.l_break = (1.0 + mat.eps_break) * l0;
        model.topology.bonds.push_back(b);
    }
    return first;
}

inline void init_model(fiber::Model& model)
{
    model.solid.resize(1); // фиктивная частица, ntotal = 0
    model.solid.inzializ_index(0);
    model.fiber.resize(0);
    model.topology.resize(0);
}

inline void set_mask(fiber::Model& model, std::size_t first, std::size_t count, std::uint8_t mk = interactions::kClampAll)
{
    if (model.clamps.fiber_mask.size() != model.fiber.size()) model.clamps.fiber_mask.resize(model.fiber.size(), 0);
    for (std::size_t k = 0; k < count; ++k) model.clamps.fiber_mask[first + k] = mk;
}

inline SimulationConfig make_cfg(int ode_solver)
{
    SimulationConfig cfg;
    cfg.apply_defaults();
    cfg.maxn = 1;
    cfg.max_interaction = 100;
    cfg.moni_particle = 0;
    cfg.ode_solver = ode_solver; // 0 = PC, 1 = Euler
    cfg.t_total = 1.0e30;        // останов -- по числу шагов (RunGuard) или явно
    return cfg;
}

inline double bond_dt_limit(const fiber::Model& m)
{
    double best = 1e30;
    for (const auto& b : m.topology.bonds)
    {
        const fiber::FiberMaterial& mat = m.materials[static_cast<std::size_t>(m.fiber.material_id[b.i])];
        const double k = mat.E * mat_area(mat) / b.l0;
        const double mi = m.fiber.mass[b.i], mj = m.fiber.mass[b.j];
        best = std::min(best, 0.8 / std::sqrt(k / (mi * mj / (mi + mj))));
    }
    return best;
}

// ---------------------------------------------------------------------------
// Кадры fiber_output
// ---------------------------------------------------------------------------
struct Frame
{
    double time = 0.0;
    std::vector<double> x[3], v[3];
    bool finite = true; // все числа кадра конечны
};

inline bool parse_number(const std::string& tok, double& out)
{
    char* end = nullptr;
    out = std::strtod(tok.c_str(), &end); // понимает nan/inf
    return end != tok.c_str() && *end == '\0';
}

// Читает все fiber_output_*.txt текущего каталога по возрастанию номера шага.
// Любой нечисловой/неконечный токен в данных -> finite = false.
inline std::vector<Frame> read_frames()
{
    namespace fs = std::filesystem;
    std::vector<std::string> names;
    for (const auto& e : fs::directory_iterator(fs::current_path()))
    {
        const std::string n = e.path().filename().string();
        if (n.rfind("fiber_output_", 0) == 0) names.push_back(n);
    }
    std::sort(names.begin(), names.end());
    std::vector<Frame> frames;
    for (const auto& n : names)
    {
        std::ifstream f(n);
        Frame fr;
        std::string line;
        bool in_particles = false, in_bonds = false;
        while (std::getline(f, line))
        {
            if (line.rfind("TIME", 0) == 0)
            {
                std::istringstream is(line);
                std::string a, t;
                is >> a >> t;
                double tv = 0.0;
                if (!parse_number(t, tv) || !std::isfinite(tv)) fr.finite = false;
                fr.time = tv;
                continue;
            }
            if (line.empty() || line[0] == '#') continue;
            if (line == "PARTICLES") { in_particles = true; in_bonds = false; continue; }
            if (line == "BONDS") { in_particles = false; in_bonds = true; continue; }
            std::istringstream is(line);
            std::vector<std::string> tok;
            std::string t;
            while (is >> t) tok.push_back(t);
            for (const auto& s : tok)
            {
                double val;
                if (!parse_number(s, val) || !std::isfinite(val)) fr.finite = false;
            }
            if (in_particles && tok.size() == 9)
            {
                for (int d = 0; d < 3; ++d)
                {
                    double xv = 0.0, vv = 0.0;
                    parse_number(tok[3 + d], xv);
                    parse_number(tok[6 + d], vv);
                    fr.x[d].push_back(xv);
                    fr.v[d].push_back(vv);
                }
            }
            (void)in_bonds;
        }
        frames.push_back(std::move(fr));
    }
    return frames;
}

// ---------------------------------------------------------------------------
// Прогон: временный каталог, тишина, time_integration; dt/time сохраняются
// между вызовами (время и dt -- сквозные).
// ---------------------------------------------------------------------------
struct Runner
{
    e2e::ScopedWorkDir wd;
    double dt = 1.0e-9; // стартовое значение; цикл перезаписывает его вычисленным шагом
    double time = 0.0;
    std::size_t steps_done = 0;
    bool threw = false;
    std::string what;

    explicit Runner(const std::string& tag) : wd("sph_clamp_fiber", tag) {}

    // Делает ровно steps шагов. save_t -- период кадров fiber_output.
    void run(SimulationConfig cfg, fiber::Model& model, std::size_t steps, double save_t)
    {
        cfg.save_t = save_t;
        cfg.print_step = 1000000000;
        RunGuard guard(0.0, 0.0, 0);
        guard.set_step_limit(steps);
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(cfg, model, 0, dt, time, 0.0, &guard);
        }
        catch (const std::exception& ex)
        {
            threw = true;
            what = ex.what();
        }
        steps_done += guard.steps();
    }

    // Прогон до времени t_end (шаги не ограничены, кроме страховки).
    void run_until(SimulationConfig cfg, fiber::Model& model, double t_end, double save_t)
    {
        cfg.t_total = t_end;
        cfg.save_t = save_t;
        cfg.print_step = 1000000000;
        RunGuard guard(0.0, 0.0, 0);
        guard.set_step_limit(5000000);
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(cfg, model, 0, dt, time, 0.0, &guard);
        }
        catch (const std::exception& ex)
        {
            threw = true;
            what = ex.what();
        }
        steps_done += guard.steps();
    }
};

inline bool bits_eq(double a, double b) { return ::same_bits(a, b); }

inline bool all_finite(const fiber::Model& m)
{
    for (std::size_t k = 0; k < m.fiber.size(); ++k)
        for (std::size_t d = 0; d < 3; ++d)
            if (!std::isfinite(m.fiber.x[d][k]) || !std::isfinite(m.fiber.vx[d][k])) return false;
    return true;
}

inline bool report_finite(const fiber::Model& m)
{
    for (std::size_t d = 0; d < 3; ++d)
        if (!std::isfinite(m.clamp_report.fiber_force[d])) return false;
    return true;
}

// Закреплённые бусины (полностью): x побитово == x0, v == 0.0 (проверка по всем компонентам).
inline bool clamped_still(const fiber::Model& m, const std::vector<double> (&x0)[3])
{
    bool ok = true;
    for (std::size_t k = 0; k < m.fiber.size(); ++k)
    {
        if (m.clamps.fiber_mask[k] != interactions::kClampAll) continue;
        for (std::size_t d = 0; d < 3; ++d)
            if (!bits_eq(m.fiber.x[d][k], x0[d][k]) || m.fiber.vx[d][k] != 0.0) ok = false;
    }
    return ok;
}

inline void save_x(const fiber::Model& m, std::vector<double> (&x0)[3])
{
    for (std::size_t d = 0; d < 3; ++d) x0[d] = m.fiber.x[d];
}

inline bool frames_clean(const std::vector<Frame>& fr)
{
    for (const auto& f : fr)
        if (!f.finite) return false;
    return !fr.empty();
}

// Число шагов сцен F1/F2/F4. ode::EulerSolver (x += dt*v^n, v += dt*a) --
// явный Эйлер, на шаге dt_limit связей (ω*dt до ~1.1 для осевых мод) его
// амплитуда растёт ~1.5 раза за шаг: длинный прогон уходит в разнос (а
// пространственный поиск контактов на разлетевшейся нити -- в гигантскую
// сетку). Поэтому для Euler -- короткий прогон, для PC -- длинный.
inline std::size_t steps_for(int solver) { return solver == 0 ? 3000 : 30; }

inline const char* solver_name(int s) { return s == 0 ? "PC" : "Euler"; }

// ---------------------------------------------------------------------------
// Сцена F1: одна прямая нить вдоль x, 21 бусина, l0 = 1 см (без натяжения).
// Закреплены полностью бусины 0..2 и 18..20. Всем свободным бусинам 3..17 задана
// поперечная (z) скорость 0.5 м/с: связи у краёв зажима вытягиваются сразу
// (реакция опоры ненулевая с первых шагов), дальше нить колеблется как струна.
// eps_break = 10 -- разрыва нет.
// ---------------------------------------------------------------------------
constexpr std::size_t kF1N = 21;
constexpr double kF1L0 = 0.01;

inline void build_f1(fiber::Model& m)
{
    init_model(m);
    m.materials = {make_mat(0.0, 10.0)};
    const double o[3] = {0, 0, 0}, dr[3] = {1, 0, 0};
    add_chain(m, kF1N, o, dr, kF1L0, kF1L0, 0, 0);
    set_mask(m, 0, 3);
    set_mask(m, kF1N - 3, 3);
    for (std::size_t k = 3; k < kF1N - 3; ++k) m.fiber.vx[2][k] = 0.5;
}

inline void test_f1()
{
    std::printf("F1: закреплённые концы нити неподвижны (реальный цикл):\n");
    for (int solver = 0; solver < 2; ++solver)
    {
        const auto t0 = std::chrono::steady_clock::now();
        const std::string tag = std::string(" [") + solver_name(solver) + "]";
        fiber::Model m;
        build_f1(m);
        std::vector<double> x0[3];
        save_x(m, x0);
        const std::size_t steps = steps_for(solver);
        Runner r("f1");
        r.run(make_cfg(solver), m, steps, bond_dt_limit(m) * static_cast<double>(steps) / 30.0);
        check(!r.threw, "F1: time_integration не бросил исключение" + tag + (r.threw ? ": " + r.what : ""));
        check(r.steps_done == steps, "F1: сделано ровно " + std::to_string(steps) + " шагов" + tag);
        check(clamped_still(m, x0), "F1: закреплённые x побитово == x0, v == 0.0" + tag);
        double maxdisp = 0.0, maxv = 0.0;
        for (std::size_t k = 3; k < kF1N - 3; ++k)
        {
            for (std::size_t d = 0; d < 3; ++d)
            {
                maxdisp = std::max(maxdisp, std::fabs(m.fiber.x[d][k] - x0[d][k]));
                maxv = std::max(maxv, std::fabs(m.fiber.vx[d][k]));
            }
        }
        check(maxdisp > 1e-4 && maxv > 1e-3, "F1: свободные бусины движутся (max смещение " +
                                                  std::to_string(maxdisp) + " м, max |v| " + std::to_string(maxv) + ")" + tag);
        check(all_finite(m) && report_finite(m), "F1: конечное итоговое состояние и реакция" + tag);
        const auto frames = read_frames();
        check(frames_clean(frames), "F1: кадры fiber_output без NaN/Inf (" + std::to_string(frames.size()) + " кадров)" + tag);
        // Закреплённые бусины и в кадрах: v == 0, x == x0 (вывод 6 знаков -- сравнение по значению вывода).
        bool frames_still = true;
        for (const auto& f : frames)
            for (std::size_t k = 0; k < kF1N; ++k)
            {
                if (m.clamps.fiber_mask[k] != interactions::kClampAll) continue;
                for (std::size_t d = 0; d < 3; ++d)
                    if (f.v[d][k] != 0.0 || std::fabs(f.x[d][k] - x0[d][k]) > 1e-6 * (1.0 + std::fabs(x0[d][k])))
                        frames_still = false;
            }
        check(frames_still, "F1: закреплённые бусины неподвижны во всех кадрах" + tag);
        const double R = std::sqrt(m.clamp_report.fiber_force[0] * m.clamp_report.fiber_force[0] +
                                   m.clamp_report.fiber_force[1] * m.clamp_report.fiber_force[1] +
                                   m.clamp_report.fiber_force[2] * m.clamp_report.fiber_force[2]);
        std::printf("    info: %s, шагов %zu, max смещение %.4g м, |R_last| = %.4g Н, %.2f с\n", solver_name(solver),
                    r.steps_done, maxdisp, R,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
}

// ---------------------------------------------------------------------------
// F2. Сцена A: нить F1 (вдоль x, закреплены 0..2 и 18..20) + свободная нить B
// вдоль y (11 бусин, l0 = 1 см), переброшенная «канатом» через зажатый сегмент
// (1,2) в точке x = 1.5 см, z = 0.9 мм (зазор меньше d = 1 мм -- контакт с
// самого начала). Скорость B: vy = 0.3 м/с (скольжение вдоль зажима), крайние
// бусины B (k = 0..2 и 8..10) тянутся вниз vz = -0.2 м/с: натяжение каната
// прижимает его к зажиму на всё время скольжения. mu_ff = 0.5 у обеих нитей.
// Сравнение с прогоном mu_ff = 0: импульс B вдоль y (среднее vy) при трении
// должен быть заметно меньше (Euler -- короткий прогон, только инварианты).
// ---------------------------------------------------------------------------
inline void build_f2(fiber::Model& m, double mu)
{
    init_model(m);
    m.materials = {make_mat(mu, 10.0)};
    const double o[3] = {0, 0, 0}, dr[3] = {1, 0, 0};
    add_chain(m, kF1N, o, dr, kF1L0, kF1L0, 0, 0);
    set_mask(m, 0, 3);
    set_mask(m, kF1N - 3, 3);
    for (std::size_t k = 3; k < kF1N - 3; ++k) m.fiber.vx[2][k] = 0.5;
    const double ob[3] = {0.015, -0.05, 0.9e-3}, db[3] = {0, 1, 0};
    const std::size_t b0 = add_chain(m, 11, ob, db, kF1L0, kF1L0, 0, 1);
    for (std::size_t k = 0; k < 11; ++k)
    {
        m.fiber.vx[1][b0 + k] = 0.3;
        m.fiber.vx[2][b0 + k] = (k <= 2 || k >= 8) ? -0.2 : 0.0;
    }
    // маска должна быть размера nfiber (set_mask выше сделал размер по первой нити -- подгоняем)
    m.clamps.fiber_mask.resize(m.fiber.size(), 0);
}

inline double mean_vy_B(const fiber::Model& m)
{
    double s = 0.0;
    for (std::size_t k = kF1N; k < m.fiber.size(); ++k) s += m.fiber.vx[1][k];
    return s / static_cast<double>(m.fiber.size() - kF1N);
}

// ---------------------------------------------------------------------------
// Прямая проверка стадии трения для подсцены F2 «две нити касаются с проникновением»
// (нить A вдоль x, нить B вдоль y, зазор 0.9 мм < d, mu_ff = 0.5). Повторяет
// последовательность src/time_integration.cpp для нить-нить: поиск контактов ->
// слоты трения -> dv = Bond + Contact -> [1] проекция закрепления -> пробная
// скорость по FiberVelocityUpdate -> compute_friction_counts -> [2]
// apply_clamp_friction_counts -> compute_friction_contributions -> reduce.
// Чтобы проверка не была вырожденной (нулевые скорости дают нулевое трение
// тривиально), к dv ДО проекции добавляется касательная нагрузка (A: +x, B: +y),
// которая на закреплённых бусинах должна быть снята проекцией [1].
// ---------------------------------------------------------------------------
struct F2FrictionProbe
{
    bool has_active = false;       // compute_fiber_contact_contributions вернул true
    std::size_t active_slots = 0;  // активные слоты трения (mu*N > 0)
    std::size_t force_slots = 0;   // слоты с has_force (mu*N > 0 и |u_t| > 0)
    double max_raw_dv_clamped = 0.0; // max |dv| закреплённых бусин ДО проекции (нагрузка + нормаль)
    double max_ut = 0.0;           // max |u_t| пробной относительной скорости по слотам
    double max_delta = 0.0;        // max |friction_delta| по всем бусинам после reduce
    double max_delta_A = 0.0, max_delta_B = 0.0;
};

inline F2FrictionProbe probe_f2_friction(fiber::Model& m, std::size_t nA, bool project)
{
    F2FrictionProbe res;
    const std::size_t nf = m.fiber.size();
    rhs::CpuExecutor exec;
    interactions::FiberContactSearchResult search;
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> bond_buf, contact_buf, fr_fiber_buf;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> fr_solid_buf;
    std::vector<rhs::ContributionBuffer<interactions::FiberCountDelta>> cnt_fiber_buf;
    std::vector<rhs::ContributionBuffer<interactions::SolidCountDelta>> cnt_solid_buf;
    std::vector<interactions::FiberFrictionContact> slots;
    std::vector<interactions::FiberSolidFrictionContact> no_solid_slots;
    std::vector<double> bead_count, solid_count;

    const auto bond = interactions::compute_fiber_bond(exec, m.fiber, m.topology, m.materials, bond_buf);
    interactions::build_fiber_contact_search(exec, m.fiber, m.topology, m.materials, search);
    res.has_active = interactions::compute_fiber_contact_contributions(exec, m.fiber, m.materials, search.primitives,
                                                                        contact_buf, slots);
    const auto contact = rhs::reduce_contributions(exec, contact_buf, nf);
    std::vector<double> dv[3];
    for (std::size_t d = 0; d < 3; ++d)
    {
        dv[d].assign(nf, 0.0);
        for (std::size_t k = 0; k < nf; ++k)
            dv[d][k] = bond[k].dv[d] + contact[k].dv[d] + (k < nA ? (d == 0 ? 5.0e3 : 0.0) : (d == 1 ? 5.0e3 : 0.0));
    }
    for (std::size_t k = 0; k < nf; ++k)
        if (m.clamps.fiber_mask[k] == interactions::kClampAll)
            for (std::size_t d = 0; d < 3; ++d) res.max_raw_dv_clamped = std::max(res.max_raw_dv_clamped, std::fabs(dv[d][k]));

    interactions::ClampWorkspace ws;
    interactions::validate_clamps(m.clamps, m.solid, 0, m.fiber);
    interactions::prepare_clamp_workspace(m.clamps, 0, nf, ws);
    if (project)
    {
        double R[3] = {};
        interactions::apply_clamp_projection(exec, ws.fiber, m.fiber.mass, dv,
                                             static_cast<std::vector<std::vector<double>>*>(nullptr), ws, R);
    }

    const double h = 1.0e-5;
    ode::FiberVelocityUpdate upd{};
    for (std::size_t d = 0; d < 3; ++d) upd.base[d] = &m.fiber.vx[d];
    upd.h = h;
    ode::SolidVelocityUpdate supd{};
    for (std::size_t d = 0; d < 3; ++d) supd.base[d] = &m.solid.vx[d];
    supd.h = 0.0;
    supd.av = nullptr;
    std::vector<std::vector<double>> zero_dvx(3);

    interactions::compute_friction_counts(exec, nf, 0, slots, no_solid_slots, cnt_fiber_buf, cnt_solid_buf, bead_count,
                                           solid_count);
    interactions::apply_clamp_friction_counts(exec, ws, bead_count, solid_count); // [2] (c_k = 0 у полностью зажатых)

    for (const auto& c : slots)
    {
        if (!c.active) continue;
        ++res.active_slots;
        interactions::FiberFrictionImpulse imp;
        interactions::fiber_friction_impulse(m.fiber, upd, dv, c, imp, bead_count.data());
        if (imp.has_force) ++res.force_slots;
        res.max_ut = std::max(res.max_ut, imp.ut_norm);
    }
    interactions::compute_friction_contributions(exec, m.fiber, upd, dv, slots, m.solid, supd, zero_dvx, no_solid_slots,
                                                  bead_count, solid_count, false, fr_fiber_buf, fr_solid_buf);
    const auto fr = rhs::reduce_contributions(exec, fr_fiber_buf, nf);
    for (std::size_t k = 0; k < nf; ++k)
        for (std::size_t d = 0; d < 3; ++d)
        {
            const double a = std::fabs(fr[k].dv[d]);
            res.max_delta = std::max(res.max_delta, a);
            (k < nA ? res.max_delta_A : res.max_delta_B) = std::max(k < nA ? res.max_delta_A : res.max_delta_B, a);
        }
    return res;
}

inline void build_f2_pair(fiber::Model& m, std::size_t maskA, std::size_t maskB_count)
{
    init_model(m);
    m.materials = {make_mat(0.5, 10.0)};
    const double oa[3] = {0, 0, 0}, da[3] = {1, 0, 0};
    const double ob[3] = {0.05, -0.05, 0.9e-3}, db[3] = {0, 1, 0};
    add_chain(m, 11, oa, da, kF1L0, kF1L0, 0, 0);
    add_chain(m, 11, ob, db, kF1L0, kF1L0, 0, 1);
    m.clamps.fiber_mask.assign(m.fiber.size(), 0);
    set_mask(m, 0, maskA);
    set_mask(m, 11, maskB_count);
}

inline void check_f2_friction_stage()
{
    // 1) Боевая последовательность, обе нити полностью закреплены: активный слот есть,
    //    пробная u_t == 0 (проекция [1] сняла нагрузку), ни одного слота с силой, вклады ровно 0.0.
    fiber::Model m;
    build_f2_pair(m, 11, 11);
    const F2FrictionProbe p = probe_f2_friction(m, 11, true);
    std::printf("    info: трение зажатая-зажатая: активных слотов %zu, слотов с силой %zu, max|dv_raw| закреплённых до проекции = %.3e, "
                "max|u_t| = %.3e, max|friction_delta| = %.3e\n",
                p.active_slots, p.force_slots, p.max_raw_dv_clamped, p.max_ut, p.max_delta);
    check(p.has_active && p.active_slots > 0, "F2: подсцена -- есть активный слот трения (mu*N > 0)");
    check(p.max_raw_dv_clamped > 1.0, "F2: подсцена -- до проекции на закреплённых бусинах ненулевое ускорение (проверка не вырождена)");
    check(p.force_slots == 0 && p.max_ut == 0.0, "F2: стадия трения -- у контакта двух полностью закреплённых нитей нет силы (has_force == false, u_t == 0)");
    check(p.max_delta == 0.0, "F2: стадия трения -- вклад трения на КАЖДОЙ бусине ровно 0.0 (до проекции [3], не маскируется суммой)");

    // 2) Положительный контроль А: проекция [1] пропущена -- тот же слот видит пробную u_t != 0 и даёт силу.
    fiber::Model mc;
    build_f2_pair(mc, 11, 11);
    const F2FrictionProbe pc = probe_f2_friction(mc, 11, false);
    std::printf("    info: контроль (без проекции [1]): слотов с силой %zu, max|u_t| = %.3e, max|friction_delta| = %.3e\n",
                pc.force_slots, pc.max_ut, pc.max_delta);
    check(pc.force_slots > 0 && pc.max_delta > 0.0,
          "F2: контроль -- без проекции [1] тот же слот даёт ненулевую силу трения (проверка чувствительна)");

    // 3) Положительный контроль Б: A закреплена, B свободна и движется -- трение есть и действует на B.
    fiber::Model mf;
    build_f2_pair(mf, 11, 0);
    for (std::size_t k = 11; k < 22; ++k) mf.fiber.vx[1][k] = 1.0;
    const F2FrictionProbe pf = probe_f2_friction(mf, 11, true);
    std::printf("    info: контроль (A закреплена, B свободна и движется): слотов с силой %zu, max|u_t| = %.3e, "
                "max|delta| B = %.3e\n", pf.force_slots, pf.max_ut, pf.max_delta_B);
    check(pf.force_slots > 0 && pf.max_delta_B > 0.0,
          "F2: контроль -- у закреплённой A и свободной движущейся B слот трения даёт ненулевую силу на B");
}

inline void test_f2()
{
    std::printf("F2: закрепление + трение нить-нить (mu_ff = 0.5):\n");
    check_f2_friction_stage();
    double vyB[2][2] = {}; // [solver][mu 0 / 0.5]
    for (int solver = 0; solver < 2; ++solver)
    {
        for (int use_mu = 0; use_mu < 2; ++use_mu)
        {
            const auto t0 = std::chrono::steady_clock::now();
            const std::string tag = std::string(" [") + solver_name(solver) + ", mu_ff=" + (use_mu ? "0.5" : "0") + "]";
            fiber::Model m;
            build_f2(m, use_mu ? 0.5 : 0.0);
            std::vector<double> x0[3];
            save_x(m, x0);
            const std::size_t steps = steps_for(solver);
            Runner r("f2");
            r.run(make_cfg(solver), m, steps, bond_dt_limit(m) * static_cast<double>(steps) / 30.0);
            check(!r.threw, "F2: time_integration не бросил исключение" + tag + (r.threw ? ": " + r.what : ""));
            check(clamped_still(m, x0), "F2: закреплённые x побитово == x0, v == 0.0" + tag);
            check(all_finite(m) && report_finite(m), "F2: нет NaN/Inf в состоянии и реакции" + tag);
            check(frames_clean(read_frames()), "F2: кадры без NaN/Inf" + tag);
            vyB[solver][use_mu] = mean_vy_B(m);
            std::printf("    info:%s mean vy(B) = %.5g, %.2f с\n", tag.c_str(), vyB[solver][use_mu],
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        }
    }
    // Эффект трения (PC, 100 шагов -- до того, как канат выйдет из режима проскальзывания и
    // динамика усложнится): среднее vy нити B с mu_ff = 0.5 заметно меньше, чем при mu_ff = 0.
    {
        double v[2];
        for (int use_mu = 0; use_mu < 2; ++use_mu)
        {
            fiber::Model m;
            build_f2(m, use_mu ? 0.5 : 0.0);
            Runner r("f2fr");
            r.run(make_cfg(0), m, 100, 1e30);
            v[use_mu] = mean_vy_B(m);
        }
        check(v[1] < v[0] - 0.01, "F2: трение о зажатый сегмент замедляет скользящую нить (PC, 100 шагов): mean vy = " +
                                      std::to_string(v[1]) + " < " + std::to_string(v[0]) + " - 0.01");
    }

    // Справка + проверка: баланс импульса по y при трении о ПОЛНОСТЬЮ закреплённый сегмент:
    // единственная сила по y на нить B -- трение, и оно должно попасть в R_fiber (точка [3]).
    {
        fiber::Model m;
        build_f2(m, 0.5);
        Runner r("f2bal");
        double sum_dPy = 0.0, sum_hRy = 0.0;
        for (int s = 0; s < 100; ++s)
        {
            double P0 = 0.0, P1 = 0.0;
            for (std::size_t k = 0; k < m.fiber.size(); ++k) P0 += m.fiber.mass[k] * m.fiber.vx[1][k];
            r.run(make_cfg(0), m, 1, 1e30);
            for (std::size_t k = 0; k < m.fiber.size(); ++k) P1 += m.fiber.mass[k] * m.fiber.vx[1][k];
            sum_dPy += P1 - P0;
            sum_hRy += 0.5 * r.dt * m.clamp_report.fiber_force[1];
        }
        std::printf("    info: баланс P_y с трением о зажим (PC, 100 шагов): sum dP_y = %.6g, sum h*R_y = %.6g\n",
                    sum_dPy, sum_hRy);
        // По y на B действуют только трение и (малая, симметричная) нормальная составляющая каната;
        // нить A по y не движется. Допуск 1e-3 от суммарного импульса, потерянного B на трение.
        check(std::fabs(sum_dPy - sum_hRy) <= 1e-3 * std::fabs(sum_dPy) + 1e-15,
              "F2: импульс трения о зажим входит в R_fiber: sum dP_y == sum h*R_y (1e-3)");
    }

    // Подсцена: две ПОЛНОСТЬЮ закреплённые нити касаются с проникновением, mu_ff = 0.5.
    // Контактная пара имеет активный слот трения, но c_k = 0 у обеих сторон
    // (бесконечная масса) -- нет деления на ноль; нити неподвижны, реакция
    // конечна, а нормальные силы контакта взаимно уравновешены (сумма реакций
    // по всем зажатым бусинам = 0, Fiber-Bond при ненагруженной нити не даёт силы).
    for (int solver = 0; solver < 2; ++solver)
    {
        const std::string tag = std::string(" [") + solver_name(solver) + "]";
        fiber::Model m;
        init_model(m);
        m.materials = {make_mat(0.5, 10.0)};
        const double oa[3] = {0, 0, 0}, da[3] = {1, 0, 0};
        const double ob[3] = {0.05, -0.05, 0.9e-3}, db[3] = {0, 1, 0}; // пересечение над сегментом, зазор 0.9 мм < d
        add_chain(m, 11, oa, da, kF1L0, kF1L0, 0, 0);
        add_chain(m, 11, ob, db, kF1L0, kF1L0, 0, 1);
        set_mask(m, 0, m.fiber.size());
        std::vector<double> x0[3];
        save_x(m, x0);
        // dt_limit не запрашивается у контакта: берём шаг связей, цикл сам ограничит.
        // Контроль: та же сцена, закреплена только нить A -- контакт реально создаёт силу.
        double Rc = 0.0;
        {
            fiber::Model mc;
            init_model(mc);
            mc.materials = {make_mat(0.5, 10.0)};
            add_chain(mc, 11, oa, da, kF1L0, kF1L0, 0, 0);
            add_chain(mc, 11, ob, db, kF1L0, kF1L0, 0, 1);
            set_mask(mc, 0, 11);
            mc.clamps.fiber_mask.resize(mc.fiber.size(), 0);
            Runner rc("f2cc_ctl");
            rc.run(make_cfg(solver), mc, 1, 1e30);
            const double* Q = mc.clamp_report.fiber_force;
            Rc = std::sqrt(Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2]);
        }
        check(Rc > 1e-4, "F2: контроль -- контакт создаёт силу на закреплённой A, |R_A| = " + std::to_string(Rc) + tag);
        Runner r("f2cc");
        r.run(make_cfg(solver), m, solver == 0 ? 300 : 30, 1e30);
        check(!r.threw, "F2: полностью закреплённые нити -- нет исключения" + tag + (r.threw ? ": " + r.what : ""));
        check(clamped_still(m, x0), "F2: обе нити неподвижны побитово, v == 0.0" + tag);
        check(all_finite(m) && report_finite(m), "F2: состояние и clamp_report конечны (нет деления на c_k = 0)" + tag);
        const double* R = m.clamp_report.fiber_force;
        const double Rn = std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
        // Нормальная сила контакта -- по оси z, равные и противоположные на нитях: сумма ~ 0.
        // Проверяем, что контакт вообще был (реакция по отдельной нити ненулевая) -- через
        // дополнительный расчёт: суммарная реакция мала по сравнению с силой на одной нити
        // (оценка сверху через EA*1e-3 = 7.85e-3 Н -- масштаб силы контакта Герца здесь заведомо выше).
        check(Rn <= 1e-12 * Rc, "F2: суммарная реакция полностью закреплённых нитей ~ 0 (третий закон), |R| = " +
                                    std::to_string(Rn) + " <= 1e-12*|R_A| (" + std::to_string(Rc) + ")" + tag);
    }
}

// ---------------------------------------------------------------------------
// F3: предварительно натянутая струна. Концы закреплены, бусины стоят с
// шагом l = L/nseg, связи l0 = l/(1+eps) => натяжение N = E*A*eps (закон
// связи: F = E*A*(l/l0 - 1)). Погонная масса в натянутом состоянии
// mu = m/l = rho*A*l0/l = rho*A/(1+eps). Аналитический период основной
// моды T = 2L*sqrt(mu/N). Начальное смещение z_k = a*sin(pi*k/nseg), v = 0.
// Параметры: L = 0.4 м, eps = 0.02, a = 4e-4 м (нелинейная добавка к N
// ~ (pi*a/L)^2/4/eps ~ 1e-4 относительно), сгущение nseg = 10, 20, 40.
// Период -- по моментам нулей z средней бусины в кадрах fiber_output (64 кадра
// за период, линейная интерполяция нуля; 3 периода, ~6 нулей; кадры печатаются
// с 6 значащими цифрами). Энергия -- из итогового состояния модели.
// Допуск: погрешность периода на самом мелком уровне <= 2 %, а по уровням не
// возрастает (с запасом 1e-5 на шум вывода); относительный дрейф полной
// энергии за 3 периода <= 1e-3 (PC).
// ---------------------------------------------------------------------------
struct F3Result
{
    double T_meas = 0.0, T_ref = 0.0, err = 0.0, drift = 0.0;
    int ncross = 0;
};

inline double total_energy(const fiber::Model& m)
{
    double ke = 0.0, pe = 0.0;
    for (std::size_t k = 0; k < m.fiber.size(); ++k)
        for (std::size_t d = 0; d < 3; ++d) ke += 0.5 * m.fiber.mass[k] * m.fiber.vx[d][k] * m.fiber.vx[d][k];
    for (const auto& b : m.topology.bonds)
    {
        const fiber::FiberMaterial& mat = m.materials[static_cast<std::size_t>(m.fiber.material_id[b.i])];
        double l2 = 0.0;
        for (std::size_t d = 0; d < 3; ++d)
        {
            const double dx = m.fiber.x[d][b.j] - m.fiber.x[d][b.i];
            l2 += dx * dx;
        }
        const double l = std::sqrt(l2);
        if (l > b.l0) pe += mat.E * mat_area(mat) * (l - b.l0) * (l - b.l0) / (2.0 * b.l0); // интеграл F(l) dl
    }
    return ke + pe;
}

inline F3Result run_f3(std::size_t nseg, int solver)
{
    const double L = 0.4, eps = 0.02, amp = 4e-4;
    const double l = L / static_cast<double>(nseg), l0 = l / (1.0 + eps);
    fiber::Model m;
    init_model(m);
    m.materials = {make_mat(0.0, 10.0)};
    const double o[3] = {0, 0, 0}, dr[3] = {1, 0, 0};
    add_chain(m, nseg + 1, o, dr, l, l0, 0, 0);
    for (std::size_t k = 0; k <= nseg; ++k)
        m.fiber.x[2][k] = amp * std::sin(pi * static_cast<double>(k) / static_cast<double>(nseg));
    set_mask(m, 0, 1);
    set_mask(m, nseg, 1);

    const fiber::FiberMaterial& mat = m.materials[0];
    const double A = mat_area(mat);
    const double N = mat.E * A * eps;
    const double mu = mat.rho * A * l0 / l;
    F3Result res;
    res.T_ref = 2.0 * L * std::sqrt(mu / N);
    const double E0 = total_energy(m);

    Runner r("f3");
    r.run_until(make_cfg(solver), m, 3.1 * res.T_ref, res.T_ref / 64.0);
    check(!r.threw, "F3: нет исключения (nseg = " + std::to_string(nseg) + ")");
    const auto frames = read_frames();
    check(frames_clean(frames), "F3: кадры без NaN/Inf (nseg = " + std::to_string(nseg) + ")");
    const std::size_t mid = nseg / 2;
    std::vector<double> tc;
    for (std::size_t i = 1; i < frames.size(); ++i)
    {
        const double z0 = frames[i - 1].x[2][mid], z1 = frames[i].x[2][mid];
        if ((z0 < 0.0) != (z1 < 0.0))
            tc.push_back(frames[i - 1].time + (frames[i].time - frames[i - 1].time) * z0 / (z0 - z1));
    }
    res.ncross = static_cast<int>(tc.size());
    if (tc.size() >= 3)
        res.T_meas = 2.0 * (tc.back() - tc.front()) / static_cast<double>(tc.size() - 1);
    res.err = res.T_meas > 0.0 ? std::fabs(res.T_meas / res.T_ref - 1.0) : 1.0;
    res.drift = std::fabs(total_energy(m) / E0 - 1.0);
    return res;
}

inline void test_f3()
{
    std::printf("F3: струна с преднатяжением, период основной моды:\n");
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t nsegs[3] = {10, 20, 40};
    F3Result res[3];
    for (int i = 0; i < 3; ++i)
    {
        res[i] = run_f3(nsegs[i], 0);
        std::printf("    info: nseg = %zu (бусин %zu): T_ref = %.6g с, T_meas = %.6g с, ошибка = %.3e, "
                    "нулей %d, дрейф энергии %.3e\n",
                    nsegs[i], nsegs[i] + 1, res[i].T_ref, res[i].T_meas, res[i].err, res[i].ncross, res[i].drift);
        check(res[i].ncross >= 5, "F3: найдено >= 5 нулей средней бусины (nseg = " + std::to_string(nsegs[i]) + ")");
        check(res[i].drift <= 1e-3, "F3: дрейф полной энергии <= 1e-3 (nseg = " + std::to_string(nsegs[i]) + ")");
    }
    check(res[2].err <= 0.02, "F3: ошибка периода на самом мелком уровне <= 2 %");
    check(res[1].err <= res[0].err + 1e-5 && res[2].err <= res[1].err + 1e-5,
          "F3: ошибка периода не растёт при сгущении (10 -> 20 -> 40)");
    std::printf("    info: F3 всего %.2f с\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}

// ---------------------------------------------------------------------------
// F4: баланс импульса на сцене F1 без Solid. Внутренние силы (связи) попарно
// уравновешены, значит для свободных бусин sum m*a = -sum_{зажатые} m*a_raw = R,
// а для шага: P^{n+1} - P^n = h * R_fiber. Подход: time_integration с
// set_step_limit(1) вызывается подряд на одной модели; каждый вызов создаёт
// новый решатель => ПЕРВЫЙ шаг: PC имеет h = dt/2 (v += dt/2 * a), Euler h = dt,
// где dt -- значение, которое time_integration записывает в свой аргумент.
// Критерий ТЗ -- чистая относительная ошибка 1e-12 на каждом шаге (без абсолютных
// добавок), режим «от покоя», см. run_f4. Число шагов: PC 200, Euler 30 (см. steps_for).
// ---------------------------------------------------------------------------
// Сцена F4: геометрия F1, скорости нулевые, свободные бусины смещены по z
// профилем 1.5e-3*sin(pi*k/20) (связи у зажима натянуты, реакция ненулевая с
// первого шага). Старт из покоя нужен для точности: P растёт от нуля вместе с
// dP, и округление P (~eps*sum(m|v|)) не затмевает dP -- иначе относительная
// ошибка 1e-12 недостижима физически.
inline void build_f4(fiber::Model& m)
{
    build_f1(m);
    for (std::size_t k = 0; k < kF1N; ++k)
    {
        m.fiber.vx[2][k] = 0.0;
        m.fiber.x[2][k] = 1.5e-3 * std::sin(pi * static_cast<double>(k) / 20.0);
    }
}

// Два режима (оба на каждом решателе):
//   reset=true  -- КРИТЕРИЙ ТЗ. Перед каждым шагом v := 0 у всех бусин (P^n = 0 ТОЧНО),
//       свободные бусины случайно смещены (LCG, до 3 % l0), чтобы реакция была
//       ненулевой и разной на каждом шаге. Тогда dP = P^{n+1} - 0 вычисляется без
//       вычитания накопленного импульса, и на КАЖДОМ шаге требуется чистая
//       относительная ошибка (без каких-либо абсолютных добавок):
//           |dP - h R| <= 1e-12 * max(|dP|, h|R|),
//       при этом на каждом шаге утверждается h|R| > 0 (проверка не вырождена).
//   reset=false -- сквозная динамика из состояния build_f4, ТОЛЬКО справка (info) и
//       проверка ненулевой реакции. Критерий 1e-12 здесь НЕ применяется: P^{n+1} - P^n
//       -- разность двух накопленных сумм, и её ошибка -- округление v^{n+1}
//       (~eps*sum m|v|, а не ошибка зажима), поэтому относительная ошибка dP
//       определяется отношением sum m|v| / |dP|, а не корректностью R. Проверку
//       точности зажима несёт режим reset.
inline void run_f4(int solver, bool reset)
{
    const std::string tag = std::string(" [") + solver_name(solver) + (reset ? ", от покоя" : ", сквозной") + "]";
    fiber::Model m;
    build_f4(m);
    Runner r("f4");
    const SimulationConfig cfg = make_cfg(solver);
    double max_rel = 0.0, max_hR = 0.0, min_hR = 1e300, max_info = 0.0;
    int n_bad = 0, n_zero_R = 0;
    bool ok_run = true;
    const int nsteps = static_cast<int>(steps_for(solver) / (solver == 0 ? 15 : 1)); // PC 200, Euler 30
    std::uint64_t lcg = 12345;
    auto rnd = [&lcg]()
    {
        lcg = lcg * 6364136223846793005ull + 1442695040888963407ull;
        return static_cast<double>(lcg >> 11) / 9007199254740992.0 - 0.5; // [-0.5, 0.5)
    };
    std::vector<double> xbase[3];
    save_x(m, xbase);
    for (int s = 0; s < nsteps; ++s)
    {
        if (reset)
        {
            for (std::size_t k = 0; k < m.fiber.size(); ++k)
                for (std::size_t d = 0; d < 3; ++d)
                {
                    m.fiber.vx[d][k] = 0.0;
                    // Положение свободных бусин ПЕРЕЗАДАЁТСЯ на каждом шаге (а не накапливается):
                    // база (build_f4) + линейный сдвиг по x от +1e-3 (бусина 3) до -1e-3
                    // (последняя свободная) -- оба крайних звена (у зажимов) гарантированно
                    // растянуты (связь работает только на растяжение), поэтому h|R| > 0 на
                    // каждом шаге -- + случайное смещение до 3 % l0.
                    if (m.clamps.fiber_mask[k] == 0)
                    {
                        const double frac = (static_cast<double>(k) - 3.0) / static_cast<double>(kF1N - 7); // 0..1
                        m.fiber.x[d][k] = xbase[d][k] + (d == 0 ? 1e-3 * (1.0 - 2.0 * frac) : 0.0) + 6e-4 * rnd();
                    }
                }
        }
        double P0[3] = {}, P1[3] = {};
        for (std::size_t k = 0; k < m.fiber.size(); ++k)
            for (std::size_t d = 0; d < 3; ++d) P0[d] += m.fiber.mass[k] * m.fiber.vx[d][k];
        r.run(cfg, m, 1, 1e30);
        if (r.threw) { ok_run = false; break; }
        const double h = solver == 0 ? 0.5 * r.dt : r.dt;
        double Pscale = 0.0;
        for (std::size_t k = 0; k < m.fiber.size(); ++k)
            for (std::size_t d = 0; d < 3; ++d)
            {
                P1[d] += m.fiber.mass[k] * m.fiber.vx[d][k];
                Pscale += m.fiber.mass[k] * std::fabs(m.fiber.vx[d][k]);
            }
        double dP2 = 0.0, hR2 = 0.0, diff2 = 0.0;
        for (std::size_t d = 0; d < 3; ++d)
        {
            const double dP = P1[d] - P0[d];
            const double hR = h * m.clamp_report.fiber_force[d];
            dP2 += dP * dP;
            hR2 += hR * hR;
            diff2 += (dP - hR) * (dP - hR);
        }
        const double ndP = std::sqrt(dP2), nhR = std::sqrt(hR2);
        const double den = std::max(ndP, nhR);
        const double err = std::sqrt(diff2);
        max_hR = std::max(max_hR, nhR);
        min_hR = std::min(min_hR, nhR);
        if (!(nhR > 0.0)) ++n_zero_R;
        if (reset)
        {
            // Чистая относительная проверка на КАЖДОМ шаге; нулевая реакция -- тоже нарушение.
            if (!(nhR > 0.0) || !(err <= 1e-12 * den)) ++n_bad;
            if (den > 0.0) max_rel = std::max(max_rel, err / den);
        }
        else if (Pscale > 0.0)
            max_info = std::max(max_info, err / Pscale); // справка: ошибка относительно sum m|v|
    }
    check(ok_run, "F4: нет исключений" + tag);
    check(max_hR > 0.0, "F4: реакция ненулевая (max h|R| = " + std::to_string(max_hR) + ")" + tag);
    if (reset)
    {
        check(n_zero_R == 0 && min_hR > 0.0,
              "F4: h|R| > 0 на каждом шаге (min h|R| = " + std::to_string(min_hR) + ")" + tag);
        check(n_bad == 0 && max_rel <= 1e-12,
              "F4: на КАЖДОМ шаге |dP - hR| <= 1e-12*max(|dP|, h|R|); max относительная ошибка = " +
                  std::to_string(max_rel) + ", нарушений " + std::to_string(n_bad) + " из " + std::to_string(nsteps) + tag);
    }
    std::printf("    info:%s %d шагов, h|R| in [%.3e, %.3e], max rel = %.3e, нарушений %d%s", tag.c_str(), nsteps, min_hR,
                max_hR, max_rel, n_bad, reset ? "\n" : "");
    if (!reset) std::printf(", справка: max |dP-hR|/sum(m|v|) = %.3e (округление P, не критерий)\n", max_info);
}

inline void test_f4()
{
    std::printf("F4: баланс импульса P^{n+1} - P^n = h*R:\n");
    for (int solver = 0; solver < 2; ++solver)
        for (int reset = 1; reset >= 0; --reset) run_f4(solver, reset != 0);
}

// ---------------------------------------------------------------------------
// F5: разрыв связи на краю зажима. Нить из 11 бусин вдоль x, l0 = 1 см,
// eps_break = 1e-3; закреплены полностью бусины 0..2; бусинам 3..10 задана
// скорость +x 0.5 м/с и сдвиг +5 % l0 (вытягивание). Связь (2,3) -- на краю зажима -- рвётся
// по нормальному критерию l > l_break в Postprocessing. Проверяется: связь
// исчезла из topology.bonds, object_id бусины 3 изменился, маска бусины 2 не
// изменилась, её x побитово == x0, v == 0, а после разрыва реакция опоры = 0
// (бусина 2 остаётся закреплённой, но связей со свободной частью нет).
// ---------------------------------------------------------------------------
inline void test_f5()
{
    std::printf("F5: разрыв связи на краю зажима:\n");
    for (int solver = 0; solver < 2; ++solver)
    {
        const std::string tag = std::string(" [") + solver_name(solver) + "]";
        fiber::Model m;
        init_model(m);
        m.materials = {make_mat(0.0, 1.0e-3)};
        const double o[3] = {0, 0, 0}, dr[3] = {1, 0, 0};
        add_chain(m, 11, o, dr, kF1L0, kF1L0, 0, 0);
        set_mask(m, 0, 3);
        // Свободная часть сдвинута на 5 % l0 (связь (2,3) натянута сверх l_break = 1.001*l0
        // уже на первом шаге) и движется с vx = 0.5 м/с дальше.
        for (std::size_t k = 3; k < 11; ++k)
        {
            m.fiber.x[0][k] += 0.05 * kF1L0;
            m.fiber.vx[0][k] = 0.5;
        }
        const std::vector<std::uint8_t> mask0 = m.clamps.fiber_mask;
        std::vector<double> x0[3];
        save_x(m, x0);
        const std::size_t nb0 = m.topology.bonds.size();
        const int id3_0 = m.fiber.object_id[3];

        Runner r("f5");
        r.run(make_cfg(solver), m, 1, 1e30); // один шаг: связь (2,3) растянута сильно -- рвётся
        check(!r.threw, "F5: нет исключения" + tag);
        const double Rx_first = m.clamp_report.fiber_force[0];
        check(Rx_first != 0.0 && std::isfinite(Rx_first),
              "F5: реакция на первом шаге (до разрыва) ненулевая, R_x = " + std::to_string(Rx_first) + tag);
        bool still_has_bond23 = false;
        for (const auto& b : m.topology.bonds)
            if (b.i == 2 && b.j == 3) still_has_bond23 = true;
        check(!still_has_bond23 && m.topology.bonds.size() < nb0, "F5: связь (2,3) разорвана по l > l_break" + tag);
        check(m.fiber.object_id[3] != id3_0 && m.fiber.object_id[3] != m.fiber.object_id[2],
              "F5: у бусины 3 новый object_id, отличный от object_id бусины 2" + tag);
        check(m.clamps.fiber_mask == mask0, "F5: маска не изменилась после разрыва" + tag);
        check(clamped_still(m, x0), "F5: закреплённые x побитово == x0, v == 0.0 после разрыва" + tag);

        r.run(make_cfg(solver), m, 100, 1e30);
        check(!r.threw, "F5: нет исключения при продолжении расчёта" + tag);
        check(m.clamps.fiber_mask == mask0 && clamped_still(m, x0),
              "F5: после 100 шагов за разрывом бусины по-прежнему неподвижны, маска прежняя" + tag);
        check(m.fiber.object_id[2] == m.fiber.object_id[0] && m.fiber.object_id[3] != m.fiber.object_id[2],
              "F5: зажатая часть остаётся одним object_id, свободная отделена" + tag);
        const double* R = m.clamp_report.fiber_force;
        check(R[0] == 0.0 && R[1] == 0.0 && R[2] == 0.0,
              "F5: реакция после разрыва ровно 0 (нет связей со свободной частью)" + tag);
        check(all_finite(m), "F5: состояние конечно" + tag);
    }
}

} // namespace clamp_fiber_it

inline void run_clamp_fiber_integration_tests()
{
    using namespace clamp_fiber_it;
    const auto t0 = std::chrono::steady_clock::now();
    auto lap = [&](const char* name, void (*f)())
    {
        const auto a = std::chrono::steady_clock::now();
        f();
        std::printf("    время %s: %.2f с\n", name, std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count());
    };
    lap("F1", test_f1);
    lap("F2", test_f2);
    lap("F3", test_f3);
    lap("F4", test_f4);
    lap("F5", test_f5);
    std::printf("    время F1-F5: %.2f с\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}
