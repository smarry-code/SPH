// test_e2e_guo2019.cpp
// Smoke-тест сцен Guo 2019 (ТЗ Guo 2019 §9): шар + ткань + Al-пластина, уровень
// smoke, V_0 = 100 м/с, t_total = 20 мкс, кадры каждые 2 мкс. Четыре регистрации
// в ctest: --config=A | 1T | 1T_A | A_1T | 2T | 4T.
//
// Прогон идёт одним непрерывным time_integration(cfg, fiber::Model&, ...) с
// монитором GuoShotMonitor в RunGuard::set_step_observer (tests/guo2019_monitor.hpp),
// во временном каталоге (ScopedWorkDir), вывод конвейера заглушён.
//
// Для 2T и 4T (пакет из 2 и 4 слоёв) дополнительно S1m, S2m, S14, S15 (ТЗ 2T/4T §6.1).
//
// Проверки S1-S13 -- жёсткие (e2e::check), справочные величины -- note/warn
// (ТЗ Guo 2019 §9). Проверки не ослабляются, чтобы тест позеленел.

#include "guo2019_monitor.hpp"
#include "guo2019_scene.hpp"

#include "e2e_curve.hpp"
#include "e2e_fiber_frames.hpp"
#include "e2e_frames.hpp"
#include "e2e_harness.hpp"

#include "run_guard.hpp"
#include "time_integration.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using e2e::check;
using e2e::note;
using e2e::warn;
using namespace guo2019;

namespace
{

// Дорогая диагностика монитора (контакты, туннелирование) -- раз в kDiagEvery шагов (ТЗ Guo 2019 §6.1: стоимость монитора <= 5 %).
constexpr std::size_t kDiagEvery = 10;

std::string fmt(double v, int prec = 4)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*e", prec, v);
    return buf;
}

std::string fmtg(double v, int prec = 6)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", prec, v);
    return buf;
}

const char* reason_name(RunGuard::Reason r)
{
    switch (r)
    {
    case RunGuard::Reason::None: return "None";
    case RunGuard::Reason::WallClock: return "WallClock";
    case RunGuard::Reason::Stalled: return "Stalled";
    case RunGuard::Reason::StepLimit: return "StepLimit";
    }
    return "?";
}

// Параметры времени и вывода smoke-прогона (ТЗ Guo 2019 §9).
void setup_smoke_cfg(GuoScene& sc)
{
    sc.cfg.t_total = kSmokeTTotal;
    sc.cfg.save_t = kSmokeSaveT;
    sc.cfg.print_step = 1000000000;
}

GuoScene make_scene(Config config, const fs::path& workdir, bool body_isolation = true)
{
    e2e::ScopedSilence silence(true);
    return build_guo_scene(config, Level::Smoke, kSmokeV0, workdir, kMuDefault, body_isolation);
}

// ---------------------------------------------------------------------------
// Чтение кадров Solid: зажатые частицы на x_0 (с точностью печати) и v == 0
// ---------------------------------------------------------------------------
struct FrameClampResult
{
    std::size_t frames = 0, rows_checked = 0, violations = 0, bad_rowcount = 0;
    std::string first;
};

FrameClampResult check_solid_frames_clamp(const fs::path& dir, const std::vector<std::uint8_t>& mask,
                                          const std::vector<double> (&x0)[3], std::size_t ntotal)
{
    FrameClampResult res;
    std::vector<fs::path> files;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(dir, ec))
    {
        std::size_t ts = 0;
        if (e.path().extension() == ".txt" && e2e::frame_timestep(e.path().filename().string(), ts))
            files.push_back(e.path());
    }
    std::sort(files.begin(), files.end());
    for (const auto& p : files)
    {
        ++res.frames;
        std::ifstream in(p);
        std::string line;
        std::size_t row = 0;
        while (std::getline(in, line))
        {
            if (line.empty() || !(line[0] >= '0' && line[0] <= '9')) continue;
            if (row < ntotal && row < mask.size() && mask[row] != 0)
            {
                std::istringstream is(line);
                double v[9];
                bool ok = true;
                for (int k = 0; k < 9 && ok; ++k) ok = static_cast<bool>(is >> v[k]);
                if (ok)
                {
                    ++res.rows_checked;
                    for (int d = 0; d < 3; ++d)
                    {
                        if (!(mask[row] & (1u << d))) continue;
                        const double x = v[3 + d], vel = v[6 + d];
                        // кадр печатает 6 значащих цифр: допуск 1e-5 отн.
                        const double tol = 1.0e-5 * std::fabs(x0[d][row]) + 1.0e-18;
                        if (!(std::fabs(x - x0[d][row]) <= tol) || !(vel == 0.0))
                        {
                            ++res.violations;
                            if (res.first.empty())
                                res.first = p.filename().string() + ", частица " + std::to_string(row) + ", компонента " +
                                            std::to_string(d);
                        }
                    }
                }
            }
            ++row;
        }
        if (row != ntotal) ++res.bad_rowcount;
    }
    return res;
}

FrameClampResult check_fiber_frames_clamp(const fs::path& dir, const std::vector<std::uint8_t>& mask,
                                          const std::vector<double> (&x0)[3])
{
    FrameClampResult res;
    for (const auto& p : e2e::list_fiber_frames(dir))
    {
        ++res.frames;
        e2e::FiberFrame fr;
        std::string err;
        if (!e2e::read_fiber_frame(p, fr, &err))
        {
            ++res.violations;
            if (res.first.empty()) res.first = err;
            continue;
        }
        for (std::size_t k = 0; k < fr.id.size(); ++k)
        {
            const std::size_t id = fr.id[k];
            if (id >= mask.size() || mask[id] == 0) continue;
            ++res.rows_checked;
            for (int d = 0; d < 3; ++d)
            {
                if (!(mask[id] & (1u << d))) continue;
                const double tol = 1.0e-5 * std::fabs(x0[d][id]) + 1.0e-18;
                if (!(std::fabs(fr.x[k][d] - x0[d][id]) <= tol) || !(fr.v[k][d] == 0.0))
                {
                    ++res.violations;
                    if (res.first.empty())
                        res.first = p.filename().string() + ", бусина " + std::to_string(id) + ", компонента " +
                                    std::to_string(d);
                }
            }
        }
    }
    return res;
}

std::size_t count_lines(const fs::path& p)
{
    std::ifstream in(p);
    std::string line;
    std::size_t n = 0;
    while (std::getline(in, line)) ++n;
    return n;
}

// ---------------------------------------------------------------------------
// S3: останов монитора на заданном шаге (только конфигурация 1T)
// ---------------------------------------------------------------------------
void test_stop_request()
{
    std::printf("\n[S3] Останов request_stop на шаге 7 (конфигурация 1T)\n");
    e2e::ScopedWorkDir wd("sph_e2e_guo2019", "1T_stop");
    GuoScene sc = make_scene(Config::T1, wd.dir());
    setup_smoke_cfg(sc);
    double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal), time = 0.0;
    RunGuard guard(100.0, 2.0e6, 200);
    std::size_t calls = 0;
    guard.set_step_observer([&](double, double) {
        if (++calls == 7) request_stop(guard);
    });
    bool threw = false;
    std::string what;
    {
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(sc.cfg, sc.model, sc.ntotal, dt, time, sc.spec.dx_s, &guard);
        }
        catch (const std::exception& ex)
        {
            threw = true;
            what = ex.what();
        }
    }
    check(!threw, "прогон с request_stop не бросил исключение" + (threw ? ": " + what : std::string()));
    check(guard.steps() == 7, "ровно 7 шагов (факт " + std::to_string(guard.steps()) + ")");
    check(guard.reason() == RunGuard::Reason::StepLimit,
          std::string("причина останова StepLimit (факт ") + reason_name(guard.reason()) + ")");
    check(!guard.failed(), "останов не считается аварией");
    const std::size_t lines = count_lines(wd.dir() / "log.txt");
    check(lines == 7, "в log.txt ровно 7 строк (факт " + std::to_string(lines) + ")");
}

// ---------------------------------------------------------------------------
// S11: изоляция тел без контакта -- исключение
// ---------------------------------------------------------------------------
void test_isolation_without_contact(Config config)
{
    std::printf("\n[S11] sph_same_body_only = 1 при contact_force_type = 0\n");
    e2e::ScopedWorkDir wd("sph_e2e_guo2019", std::string(config_name(config)) + "_s11");
    GuoScene sc = make_scene(config, wd.dir());
    setup_smoke_cfg(sc);
    sc.cfg.sph_same_body_only = true;
    sc.cfg.contact_force_type = 0;
    double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal), time = 0.0;
    RunGuard guard(100.0, 2.0e6, 200);
    guard.set_step_limit(2);
    bool threw_rt = false, threw_other = false;
    std::string what;
    {
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(sc.cfg, sc.model, sc.ntotal, dt, time, sc.spec.dx_s, &guard);
        }
        catch (const std::runtime_error& ex)
        {
            threw_rt = true;
            what = ex.what();
        }
        catch (...)
        {
            threw_other = true;
        }
    }
    check(threw_rt && !threw_other, "исключение std::runtime_error при старте" + (threw_rt ? ": " + what : std::string()));
}

// ---------------------------------------------------------------------------
// S9: нить не прошла сквозь целую пластину
// ---------------------------------------------------------------------------
void check_s9(const GuoScene& sc, const MonitorStats& st)
{
    if (!(sc.config == Config::T1_A || sc.config == Config::A_T1)) return;
    std::printf("\n[S9] Порядок слоёв\n");
    if (!(st.m_ej == 0.0))
    {
        note("m_ej = " + fmtg(st.m_ej) + " кг на t_end > 0: пластина пробита, S9 не применяется");
        return;
    }
    const auto& s = sc.model.solid;
    const auto& fb = sc.model.fiber;
    const std::size_t nx = sc.plate_nx, ny = sc.plate_ny;
    const double dxx = sc.dx_p_axis[0], dyy = sc.dx_p_axis[1];
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> zmin(nx * ny, inf), zmax(nx * ny, -inf);
    auto cell = [&](double x, double y, long& ix, long& iy) {
        ix = static_cast<long>(std::floor((x + sc.plate_half_x) / dxx));
        iy = static_cast<long>(std::floor((y + sc.plate_half_y) / dyy));
        return ix >= 0 && iy >= 0 && ix < static_cast<long>(nx) && iy < static_cast<long>(ny);
    };
    for (std::size_t i = sc.plate_begin; i < sc.plate_end; ++i)
    {
        long ix, iy;
        if (!cell(s.x[0][i], s.x[1][i], ix, iy)) continue;
        const std::size_t c = static_cast<std::size_t>(ix) * ny + static_cast<std::size_t>(iy);
        zmin[c] = std::min(zmin[c], s.x[2][i]);
        zmax[c] = std::max(zmax[c], s.x[2][i]);
    }
    const bool fabric_above = sc.config == Config::T1_A; // 1T/A: ткань спереди (выше пластины); A/1T -- сзади
    std::size_t in_window = 0, skipped = 0, bad = 0;
    for (std::size_t k = 0; k < fb.size(); ++k)
    {
        const double bx = fb.x[0][k], by = fb.x[1][k], bz = fb.x[2][k];
        if (std::fabs(bx) > sc.plate_half_x || std::fabs(by) > sc.plate_clamp_edge) continue;
        long ix, iy;
        if (!cell(bx, by, ix, iy)) continue;
        const std::size_t c = static_cast<std::size_t>(ix) * ny + static_cast<std::size_t>(iy);
        if (!(zmin[c] <= zmax[c]))
        {
            ++skipped;
            continue;
        }
        ++in_window;
        if (fabric_above ? (bz < zmin[c]) : (bz > zmax[c])) ++bad;
    }
    note("бусин в окне: " + std::to_string(in_window) + ", ячеек без частиц пластины пропущено: " +
         std::to_string(skipped));
    check(bad == 0, "ни одна бусина не оказалась по другую сторону пластины (нарушений " + std::to_string(bad) + ")");
}

// ---------------------------------------------------------------------------
// Основной прогон и проверки S1-S10, S12, S13
// ---------------------------------------------------------------------------
void run_main_case(Config config, bool keep_workdir, double time_limit, bool body_isolation)
{
    const std::string cname = config_name(config);
    e2e::ScopedWorkDir wd("sph_e2e_guo2019", cname);
    if (keep_workdir)
    {
        wd.keep();
        note("рабочий каталог сохранён: " + wd.dir().string());
    }
    GuoScene sc = make_scene(config, wd.dir(), body_isolation);
    setup_smoke_cfg(sc);
    const SimulationConfig& cfg = sc.cfg;
    fiber::Model& model = sc.model;

    std::printf("\n[Сцена] конфигурация %s\n", cname.c_str());
    print_scene_summary(stdout, sc);

    // ------------------------------------------------------------- S1 (до прогона)
    std::printf("\n[S1] Числа частиц, нитей, масса шара\n");
    check(sc.n_ball == sc.expected.n_ball, "частиц шара " + std::to_string(sc.n_ball) + " == формула " +
                                               std::to_string(sc.expected.n_ball));
    check(sc.n_plate == sc.expected.n_plate, "частиц пластины " + std::to_string(sc.n_plate) + " == формула " +
                                                 std::to_string(sc.expected.n_plate));
    check(sc.ntotal == sc.expected.ntotal, "Solid всего " + std::to_string(sc.ntotal) + " == формула " +
                                               std::to_string(sc.expected.ntotal));
    check(sc.actual_beads == sc.expected.beads, "бусин " + std::to_string(sc.actual_beads) + " == формула " +
                                                    std::to_string(sc.expected.beads));
    check(sc.actual_bonds == sc.expected.bonds, "связей " + std::to_string(sc.actual_bonds) + " == формула " +
                                                    std::to_string(sc.expected.bonds));
    {
        const double rel = std::fabs(sc.ball_mass - sc.ball_mass_target) / sc.ball_mass_target;
        check(rel <= kTolMassRel, "масса шара " + fmtg(sc.ball_mass, 12) + " кг против заданной " +
                                      fmtg(sc.ball_mass_target, 12) + " (отн. " + fmt(rel, 2) + " <= 1e-12)");
    }
    if (sc.has_fabric)
    {
        note("поверхностная плотность ткани в окне (связи целиком внутри): " + fmtg(sc.areal_density_window, 5) +
             " кг/м^2, окно мало (smoke), с статьёй 0.28 не сравнивается; summary (все бусины): " +
             fmtg(sc.areal_density_summary, 5) + " кг/м^2");
        // Зажатые бусины -- только у нитей вдоль y, поровну на концах (по факту fiber_mask)
        const auto& fm = model.clamps.fiber_mask;
        std::size_t bad_dir = 0, bad_balance = 0, bad_mask = 0, total = 0;
        for (const auto& y : sc.yarns)
        {
            const std::size_t mid = y.bead_begin + (y.bead_end - y.bead_begin) / 2;
            std::size_t lo = 0, hi = 0;
            for (std::size_t k = y.bead_begin; k < y.bead_end; ++k)
            {
                if (k >= fm.size() || fm[k] == 0) continue;
                if (fm[k] != interactions::kClampAll) ++bad_mask;
                (k < mid ? lo : hi) += 1;
            }
            total += lo + hi;
            if (y.direction == 'x' && lo + hi != 0) ++bad_dir;
            if (y.direction == 'y' && (lo != hi || lo != y.n_clamp_lo || hi != y.n_clamp_hi)) ++bad_balance;
        }
        check(bad_dir == 0, "зажатых бусин у нитей вдоль x нет (нарушений " + std::to_string(bad_dir) + ")");
        check(bad_balance == 0, "у каждой нити вдоль y зажато поровну на концах и столько, сколько дал генератор "
                                "(нарушений " + std::to_string(bad_balance) + ")");
        check(bad_mask == 0, "маска зажатых бусин kClampAll у всех (нарушений " + std::to_string(bad_mask) + ")");
        check(total == sc.expected.clamped_beads && total == sc.actual_clamped_beads,
              "зажатых бусин по fiber_mask " + std::to_string(total) + " == формула " +
                  std::to_string(sc.expected.clamped_beads));
    }
    if (sc.n_plies >= 2)
    {
        const int N = sc.n_plies;
        const int nx = sc.W.params.nx, ny = sc.W.params.ny;
        std::printf("\n[S1m] Пакет слоёв: %d слоёв\n", N);
        const auto& fm = model.clamps.fiber_mask;
        std::size_t sum_beads = 0, sum_bonds = 0, sum_clamped = 0, yarn_total = 0, clamped_first = 0;
        for (int k = 0; k < N; ++k)
        {
            const auto& pl = sc.ply_layout[static_cast<std::size_t>(k)];
            const auto& pe = sc.ply_expected[static_cast<std::size_t>(k)];
            const std::string tag = "слой " + std::to_string(k + 1) + ": ";
            const std::size_t beads = pl.bead_end - pl.bead_begin;
            std::size_t bonds = 0;
            for (const auto& b : model.topology.bonds)
                if (b.i >= pl.bead_begin && b.i < pl.bead_end) ++bonds;
            std::size_t clamped = 0, bad_dir = 0, bad_bal = 0, bad_ply = 0, bad_oid = 0;
            for (std::size_t y = pl.yarn_begin; y < pl.yarn_end; ++y)
            {
                const auto& Y = sc.yarns[y];
                if (Y.ply != k + 1) ++bad_ply;
                // object_id нитей слоя k (с единицы) -- ровно [(k-1)(nx+ny); k(nx+ny)), по порядку
                if (Y.object_id != k * (nx + ny) + static_cast<int>(y - pl.yarn_begin)) ++bad_oid;
                const std::size_t mid = Y.bead_begin + (Y.bead_end - Y.bead_begin) / 2;
                std::size_t lo = 0, hi = 0;
                for (std::size_t b = Y.bead_begin; b < Y.bead_end; ++b)
                    if (b < fm.size() && fm[b] != 0) (b < mid ? lo : hi) += 1;
                clamped += lo + hi;
                if (Y.direction == 'x' && lo + hi != 0) ++bad_dir;
                if (Y.direction == 'y' && (lo != hi || lo != Y.n_clamp_lo || hi != Y.n_clamp_hi)) ++bad_bal;
            }
            if (k == 0) clamped_first = clamped;
            check(beads == pe.beads, tag + "бусин " + std::to_string(beads) + " == формула " + std::to_string(pe.beads));
            check(bonds == pe.bonds, tag + "связей " + std::to_string(bonds) + " == формула " + std::to_string(pe.bonds));
            check(clamped == pe.clamped_beads && clamped == clamped_first,
                  tag + "зажатых бусин " + std::to_string(clamped) + " == формула " + std::to_string(pe.clamped_beads) +
                      " и == слою 1 (" + std::to_string(clamped_first) + ")");
            check(bad_dir == 0 && bad_bal == 0, tag + "зажаты только нити вдоль y, поровну на концах (нарушений " +
                                                    std::to_string(bad_dir + bad_bal) + ")");
            check(pl.yarn_end - pl.yarn_begin == static_cast<std::size_t>(nx + ny) && bad_ply == 0 && bad_oid == 0,
                  tag + "нитей " + std::to_string(pl.yarn_end - pl.yarn_begin) + " == nx+ny, object_id в диапазоне [" +
                      std::to_string(k * (nx + ny)) + "; " + std::to_string((k + 1) * (nx + ny)) + ") (нарушений " +
                      std::to_string(bad_ply + bad_oid) + ")");
            note(tag + "поверхностная плотность: в окне " + fmtg(pe.areal_density_window, 5) + " кг/м^2, summary " +
                 fmtg(pe.areal_density_summary, 5) + " кг/м^2");
            sum_beads += beads;
            sum_bonds += bonds;
            sum_clamped += clamped;
            yarn_total += pl.yarn_end - pl.yarn_begin;
        }
        check(sum_beads == sc.actual_beads && sum_bonds == sc.actual_bonds,
              "суммы по слоям: бусин " + std::to_string(sum_beads) + ", связей " + std::to_string(sum_bonds) +
                  " == всего в модели " + std::to_string(sc.actual_beads) + ", " + std::to_string(sc.actual_bonds));
        check(sum_clamped == sc.expected.clamped_beads && yarn_total == sc.yarns.size(),
              "сумма зажатых по слоям " + std::to_string(sum_clamped) + " == формула " +
                  std::to_string(sc.expected.clamped_beads) + " (сверка с model.clamp_report -- после прогона), нитей " +
                  std::to_string(yarn_total) + " == " + std::to_string(sc.yarns.size()));
        note("плотность, сумма по слоям: summary " + fmtg(sc.areal_density_summary, 5) + " кг/м^2, в окне " +
             fmtg(sc.areal_density_window, 5) + " кг/м^2");
    }
    if (const auto es = elastic_by_itype(cfg.mats, kBallItype))
        note("сталь (itype 3): E = " + fmtg(es->E * 1e-9, 4) + " ГПа, nu = " + fmtg(es->nu, 4));

    // ------------------------------------------------------------- S2 (до прогона)
    std::printf("\n[S2] Стартовая укладка\n");
    {
        const BallGap bg = ball_target_gap(sc);
        check(std::fabs(bg.analytic - kBallGap0) <= kTolGap, "зазор шар--мишень (аналитический) " +
                                                                  fmt(bg.analytic, 12) + " м == g_0 = " + fmt(kBallGap0) +
                                                                  " +- 1e-12 (откл. " + fmt(bg.analytic - kBallGap0, 2) + ")");
        note("зазор шар--мишень по граням частиц " + fmt(bg.particle_min, 6) + " м, по контактным радиусам " +
             fmt(bg.contact_min, 6) + " м");
    }
    if (sc.has_fabric)
    {
        const GapReport gr = fabric_solid_min_gap(sc);
        auto report = [&](const char* who, const GapSide& g) {
            if (!g.found)
            {
                note(std::string("зазор ткань--") + who + ": пар в окне " + fmt(gr.margin, 2) +
                     " м нет (зазор больше окна)");
                return;
            }
            char loc[512];
            std::snprintf(loc, sizeof loc,
                          "зазор ткань--%s: min %.6e м (расстояние %.6e, r_f %.4e, r_s %.4e); связь %zu, частица %zu; "
                          "точка оси (%.5e, %.5e, %.5e), центр частицы (%.5e, %.5e, %.5e)",
                          who, g.gap, g.dist, g.r_f, g.r_s, g.bond, g.particle, g.seg_point[0], g.seg_point[1],
                          g.seg_point[2], g.part_pos[0], g.part_pos[1], g.part_pos[2]);
            note(loc);
            check(g.gap >= 0.0, std::string("нет перекрытия ткань--") + who + ": зазор " + fmt(g.gap, 4) + " м >= 0");
        };
        report("шар", gr.ball);
        if (sc.has_plate) report("пластина", gr.plate);
    }

    if (sc.n_plies >= 2)
    {
        std::printf("\n[S2m] Зазор между слоями (до прогона)\n");
        const PlyGapReport pg = ply_vertex_gaps(sc);
        check(pg.valid && pg.n_pairs > 0, "проверено пар (пара слоёв, пересечение): " + std::to_string(pg.n_pairs));
        check(std::fabs(pg.min_gap - kEpsPly) <= kTolGap && std::fabs(pg.max_gap - kEpsPly) <= kTolGap,
              "зазор слоёв по вершинным бусинам: min " + fmt(pg.min_gap, 12) + " м, max " + fmt(pg.max_gap, 12) +
                  " м == eps_ply = " + fmt(kEpsPly, 3) + " +- 1e-12 (откл. min " + fmt(pg.min_gap - kEpsPly, 2) +
                  ", max " + fmt(pg.max_gap - kEpsPly, 2) + ")");
    }

    // ------------------------------------------------------------------ прогон
    // Снимок x_0 для проверки кадров (S6)
    std::vector<double> sx0[3], fx0[3];
    for (int d = 0; d < 3; ++d)
    {
        sx0[d].assign(model.solid.x[d].begin(), model.solid.x[d].begin() + static_cast<std::ptrdiff_t>(sc.ntotal));
        fx0[d].assign(model.fiber.x[d].begin(), model.fiber.x[d].end());
    }
    const std::vector<std::uint8_t> solid_mask = model.clamps.solid_mask;
    const std::vector<std::uint8_t> fiber_mask = model.clamps.fiber_mask;

    double dt = compute_cfl_timestep(sc.cfg, model.solid, sc.ntotal), time = 0.0;
    const double dt_init = dt;
    RunGuard guard(time_limit, 2.0e6, 200);
    MonitorOptions mo;
    mo.csv_path = (wd.dir() / "monitor.csv").string();
    mo.guard = &guard;
    mo.n_diag = kDiagEvery;
    GuoShotMonitor monitor(sc, mo);
    GuoShotMonitor::FFSplitCheck split_step1, split_last;
    bool have_split_step1 = false;
    guard.set_step_observer([&](double d, double t) {
        monitor.on_step(d, t);
        if (sc.n_plies >= 2 && !have_split_step1)
        {
            // S15 на шаге 1: оба способа счёта по одному поиску на состоянии после шага 1
            split_step1 = monitor.ff_split_check();
            have_split_step1 = true;
        }
    });

    bool threw = false;
    std::string what;
    const auto w0 = std::chrono::steady_clock::now();
    {
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(sc.cfg, model, sc.ntotal, dt, time, sc.spec.dx_s, &guard);
        }
        catch (const std::exception& ex)
        {
            threw = true;
            what = ex.what();
        }
    }
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
    monitor.finalize();
    const MonitorStats& st = monitor.stats();
    if (sc.n_plies >= 2) split_last = monitor.ff_split_check(); // состояние на последнем шаге

    std::printf("\n[Прогон] %zu шагов, t = %s с, %.2f с стенных, начальный dt (CFL) = %s с, сторож: %s\n", st.steps,
                fmt(time).c_str(), wall, fmt(dt_init).c_str(), reason_name(guard.reason()));
    check(!threw, "нет исключений в time_integration" + (threw ? ": " + what : std::string()));
    note("монитор: " + fmtg(st.monitor_wall_s, 4) + " с = " + fmtg(100.0 * st.monitor_wall_s / std::max(wall, 1e-12), 3) +
         " % времени прогона");

    // ---------------------------------------------------------------------- S1 (после прогона)
    std::printf("\n[S1] Зажатые после прогона (model.clamp_report)\n");
    check(model.clamp_report.solid_clamped == sc.expected.clamped_plate,
          "зажатых частиц Solid " + std::to_string(model.clamp_report.solid_clamped) + " == формула " +
              std::to_string(sc.expected.clamped_plate));
    check(model.clamp_report.fiber_clamped == sc.expected.clamped_beads,
          "зажатых бусин " + std::to_string(model.clamp_report.fiber_clamped) + " == формула " +
              std::to_string(sc.expected.clamped_beads));
    if (sc.n_plies >= 2)
    {
        std::size_t sum = 0;
        for (const auto& pe : sc.ply_expected) sum += pe.clamped_beads;
        check(model.clamp_report.fiber_clamped == sum,
              "S1m: model.clamp_report.fiber_clamped " + std::to_string(model.clamp_report.fiber_clamped) +
                  " == сумма по слоям " + std::to_string(sum));
    }

    // ---------------------------------------------------------------------- S2 (первый шаг)
    std::printf("\n[S2] Первый шаг\n");
    if (sc.has_plate)
        check(st.plate_vmax_step1 == 0.0, "скорости частиц пластины после шага 1 побитово 0 (max |v| = " +
                                              fmt(st.plate_vmax_step1, 3) + ")");
    if (sc.has_fabric)
    {
        note("максимальная |v| бусин после шага 1: " + fmt(st.bead_vmax_step1, 3) + " м/с");
        if (st.bead_vmax_step1 > 0.0)
            warn("бусины движутся на первом шаге (max |v| = " + fmt(st.bead_vmax_step1, 3) +
                 " м/с): ложные контакты Fiber-Fiber или контакт с Solid");
    }

    if (sc.n_plies >= 2)
    {
        std::printf("\n[S2m] Контакты нить--нить на исходной укладке и после шага 1\n");
        check(st.ff_intra_initial == 0 && st.ff_inter_initial == 0,
              "исходная укладка: контактов внутри слоя " + std::to_string(st.ff_intra_initial) + ", между слоями " +
                  std::to_string(st.ff_inter_initial) + " (ожидание 0 и 0)");
        check(st.ff_intra_step1 == 0 && st.ff_inter_step1 == 0,
              "после шага 1: контактов внутри слоя " + std::to_string(st.ff_intra_step1) + ", между слоями " +
                  std::to_string(st.ff_inter_step1) + " (ожидание 0 и 0)");

        std::printf("\n[S14] Порядок слоёв\n");
        check(st.ply_order_violations_step1 == 0,
              "шаг 1: нарушений порядка слоёв " + std::to_string(st.ply_order_violations_step1));
        check(st.ply_order_violations == 0,
              "максимум нарушений за все шаги диагностики " + std::to_string(st.ply_order_violations) +
                  (st.t_ply_order_first >= 0.0 ? " (первое при t = " + fmt(st.t_ply_order_first, 3) + " с)"
                                               : std::string()));

        std::printf("\n[S15] Послойный счёт контактов против compute_fiber_contact_diagnostics\n");
        auto s15 = [&](const char* when, const GuoShotMonitor::FFSplitCheck& r) {
            check(r.match(), std::string(when) + ": диагностика ядра " + std::to_string(r.old_count) +
                                 " контактов, max delta/d " + fmtg(r.old_delta_d, 9) + "; послойно внутри " +
                                 std::to_string(r.intra) + " + между " + std::to_string(r.inter) + " = " +
                                 std::to_string(r.intra + r.inter) + ", max delta/d " +
                                 fmtg(std::max(r.delta_intra_d, r.delta_inter_d), 9) + " (внутри " +
                                 fmtg(r.delta_intra_d, 6) + ", между " + fmtg(r.delta_inter_d, 6) + ")");
        };
        check(have_split_step1, "проверка S15 на шаге 1 выполнена");
        s15("шаг 1", split_step1);
        s15("последний шаг", split_last);
    }

    // ---------------------------------------------------------------------- S3
    std::printf("\n[S3] Сторож и дошли ли до t_total\n");
    check(!guard.failed(), std::string("сторож не оборвал расчёт аварийно (причина ") + reason_name(guard.reason()) + ")");
    check(time >= sc.cfg.t_total || st.stopped_by_monitor,
          "расчёт дошёл до t_total = " + fmt(sc.cfg.t_total) + " (t = " + fmt(time) + ") или штатно остановлен монитором");

    // ---------------------------------------------------------------------- S4
    std::printf("\n[S4] Кадры\n");
    const e2e::FrameScan fs_scan = e2e::scan_frames(wd.dir(), e2e::FrameProbe{});
    const e2e::FrameCheck& fc = fs_scan.health;
    check(fc.frames >= 2, "кадров Solid: " + std::to_string(fc.frames));
    check(fc.rows == fc.frames * sc.ntotal && fc.short_rows == 0,
          "строк в кадрах Solid " + std::to_string(fc.rows) + " == кадров x ntotal = " +
              std::to_string(fc.frames * sc.ntotal));
    check(fc.bad_finite == 0 && fc.bad_rho == 0 && fc.bad_mass == 0,
          "Solid: нет NaN/Inf, rho > 0, m > 0 (не-конечных " + std::to_string(fc.bad_finite) + ", rho<=0 " +
              std::to_string(fc.bad_rho) + ", m<=0 " + std::to_string(fc.bad_mass) + ")" +
              (fc.first_problem.empty() ? std::string() : " " + fc.first_problem));
    const e2e::FiberFrameScan ff_scan = e2e::scan_fiber_frames(wd.dir());
    if (sc.has_fabric)
    {
        check(ff_scan.frames == fc.frames, "кадров нитей " + std::to_string(ff_scan.frames) + " == кадров Solid " +
                                                std::to_string(fc.frames));
        check(ff_scan.health, "нити: все кадры разобраны, нет NaN/Inf в x, v (не-конечных " +
                                  std::to_string(ff_scan.nonfinite) + ")" +
                                  (ff_scan.first_problem.empty() ? std::string() : " " + ff_scan.first_problem));
        bool counts_ok = true;
        for (std::size_t n : ff_scan.n_particles) counts_ok = counts_ok && n == sc.actual_beads;
        check(counts_ok, "число бусин в каждом кадре нитей равно " + std::to_string(sc.actual_beads));
    }
    else
        check(ff_scan.frames == 0, "ткани нет: кадров нитей не пишется (факт " + std::to_string(ff_scan.frames) + ")");

    // ---------------------------------------------------------------------- S5
    std::printf("\n[S5] Шаг по времени\n");
    {
        const double dtf = guard.dt_first(), dtm = guard.dt_min();
        note("dt_first = " + fmt(dtf) + " с, dt_min = " + fmt(dtm) + " с, шагов " + std::to_string(guard.steps()));
        check(dtf > 0.0 && dtm > 0.0 && dtf / dtm <= kDtRatioMax,
              "dt_first/dt_min = " + fmtg(dtm > 0.0 ? dtf / dtm : 0.0, 4) + " <= " + fmtg(kDtRatioMax, 3));
        const double expected = sc.cfg.t_total / dtf;
        check(static_cast<double>(guard.steps()) < kStepsFactor * expected,
              "шагов " + std::to_string(guard.steps()) + " < " + fmtg(kStepsFactor * expected, 6) + " (10x от оценки " +
                  fmtg(expected, 6) + " по первому dt)");
    }

    // ---------------------------------------------------------------------- S6
    std::printf("\n[S6] Неподвижность зажимов\n");
    check(st.clamp_violations == 0,
          "на каждом шаге зажатые x == x_0 побитово и v == 0.0 (нарушений " + std::to_string(st.clamp_violations) + ")" +
              (st.clamp_first_violation.empty() ? std::string() : ": " + st.clamp_first_violation));
    {
        const FrameClampResult fr = check_solid_frames_clamp(wd.dir(), solid_mask, sx0, sc.ntotal);
        check(fr.violations == 0 && fr.bad_rowcount == 0,
              "кадры Solid: зажатые частицы на x_0 (допуск печати 1e-5) и v = 0 (проверено строк " +
                  std::to_string(fr.rows_checked) + ", нарушений " + std::to_string(fr.violations) + ")" +
                  (fr.first.empty() ? std::string() : ": " + fr.first));
        if (sc.has_plate) check(fr.rows_checked > 0, "в кадрах Solid есть зажатые частицы пластины");
        if (sc.has_fabric)
        {
            const FrameClampResult ft = check_fiber_frames_clamp(wd.dir(), fiber_mask, fx0);
            check(ft.violations == 0 && ft.rows_checked > 0,
                  "кадры нитей: зажатые бусины на x_0 (допуск печати 1e-5) и v = 0 (проверено строк " +
                      std::to_string(ft.rows_checked) + ", нарушений " + std::to_string(ft.violations) + ")" +
                      (ft.first.empty() ? std::string() : ": " + ft.first));
        }
    }

    // ---------------------------------------------------------------------- S7
    std::printf("\n[S7] Массы\n");
    check(st.mass_violations == 0, "суммы масс Solid и нитей равны стартовым побитово на каждом шаге (нарушений " +
                                       std::to_string(st.mass_violations) + ")");

    // ---------------------------------------------------------------------- S8
    std::printf("\n[S8] Контакт с мишенью\n");
    check(std::fabs(st.v_z_end) < sc.V0 - kMinDecel, "|v_z(t_end)| = " + fmtg(std::fabs(st.v_z_end), 8) + " < V_0 - 0.5 = " +
                                                         fmtg(sc.V0 - kMinDecel, 6) + " м/с");

    // ---------------------------------------------------------------------- S9
    check_s9(sc, st);

    // ---------------------------------------------------------------------- S12, S13
    std::printf("\n[S12] Энергобаланс\n");
    {
        const double rmax = st.E_tot_max / st.E_tot_0, rend = st.E_tot_end / st.E_tot_0;
        check(rmax <= kEnergyRatioFail, "E_tot_max/E_tot_0 = " + fmtg(rmax, 8) + " <= " + fmtg(kEnergyRatioFail, 4) +
                                            " (E_tot_end/E_tot_0 = " + fmtg(rend, 8) + ")");
        if (rmax > kEnergyRatioWarn && rmax <= kEnergyRatioFail)
            warn("E_tot_max/E_tot_0 = " + fmtg(rmax, 8) + " в (1.005; 1.02]");
        monitor.print_energy_decomposition(stdout);
    }
    std::printf("\n[S13] Туннелирование\n");
    check(st.n_tunnel == 0, "n_tunnel = " + std::to_string(st.n_tunnel) + " (нить не прошла сквозь шар)");

    // ---------------------------------------------------------------------- справочные
    std::printf("\n[Справка]\n");
    if (sc.has_fabric)
    {
        note("max delta/d нить--нить = " + fmtg(st.max_delta_ff_d, 5) + ", нить--Solid = " + fmtg(st.max_delta_fs_d, 5));
        if (st.max_delta_ff_d > kDeltaOverDWarn) warn("max delta/d нить--нить = " + fmtg(st.max_delta_ff_d, 5) + " > 0.3");
        if (st.max_delta_fs_d > kDeltaOverDWarn) warn("max delta/d нить--Solid = " + fmtg(st.max_delta_fs_d, 5) + " > 0.3");
        note("контакты Fiber-Fiber: на исходной укладке " + std::to_string(st.ff_contacts_initial) + " (delta/d " +
             fmtg(st.ff_delta_d_initial, 4) + "), после шага 1 " + std::to_string(st.ff_contacts_step1) + " (delta/d " +
             fmtg(st.ff_delta_d_step1, 4) + "), максимум за прогон " + std::to_string(st.ff_contacts_max));
        note("contacts_per_crossing: максимум " + fmtg(st.contacts_per_crossing, 5) + ", после шага 1 " +
             fmtg(st.contacts_per_crossing_step1, 5));
        if (sc.n_plies >= 2)
        {
            note("contacts_per_crossing_inter: максимум " + fmtg(st.contacts_per_crossing_inter, 5) + ", после шага 1 " +
                 fmtg(st.contacts_per_crossing_inter_step1, 5) + "; max delta/d между слоями " +
                 fmtg(st.max_delta_ff_inter_d, 5));
            if (st.contacts_per_crossing_inter > 1.0)
                warn("contacts_per_crossing_inter = " + fmtg(st.contacts_per_crossing_inter, 5) + " > 1");
            if (st.max_delta_ff_inter_d > 0.3)
                warn("max delta/d нить--нить между слоями = " + fmtg(st.max_delta_ff_inter_d, 5) + " > 0.3");
            for (int k = 0; k < sc.n_plies; ++k)
            {
                const std::size_t q = static_cast<std::size_t>(k);
                note("слой " + std::to_string(k + 1) + ": порвано вдоль y " + std::to_string(st.n_broken_v_p[q]) +
                     ", вдоль x " + std::to_string(st.n_broken_h_p[q]) + ", t первого разрыва " +
                     (st.t_first_break_p[q] >= 0.0 ? fmt(st.t_first_break_p[q], 3) + " с" : std::string("нет")) +
                     ", прогиб слоя на t_end " + fmt(st.w_f_end_p[q], 3) + " м, втягивание кромок на t_end " +
                     fmt(st.pullin_end_p[q], 3) + " м");
            }
        }
        if (st.contacts_per_crossing > 1.0)
            warn("contacts_per_crossing = " + fmtg(st.contacts_per_crossing, 5) + " > 1 контакта на пересечение");
        if (st.ff_contacts_initial > 0 || st.ff_contacts_step1 > 0)
            note("стартовые контакты Fiber-Fiber не равны нулю (ожидание при A = 0.6 d: 0)");
        note("разрывов связей: " + std::to_string(st.n_bonds_0 - st.n_bonds_end) + ", нитей порвано: вдоль y " +
             std::to_string(st.n_broken_v) + ", вдоль x " + std::to_string(st.n_broken_h) + "; w_f_max = " +
             fmt(st.w_f_max, 3) + " м");
        note("реакции: F_clamp_s_z_max = " + fmtg(st.F_s_z_max, 5) + " Н (t = " + fmt(st.F_s_z_max_t, 3) +
             "), F_clamp_f_z_max = " + fmtg(st.F_f_z_max, 5) + " Н (t = " + fmt(st.F_f_z_max_t, 3) + ")");
    }
    note("w_p_max = " + fmt(st.w_p_max, 3) + " м, m_ej = " + fmtg(st.m_ej, 5) + " кг, z_c = " + fmt(st.z_c_end, 6) +
         " м, v_z = " + fmtg(st.v_z_end, 8) + " м/с");
    note(std::string("исход: ") + outcome_name(st.outcome) + ", V_R = " +
         ((st.outcome == Outcome::Perforated || st.outcome == Outcome::NotPerforated) ? fmtg(st.V_R, 6) : "не определена") +
         " (на smoke 20 мкс не проверяется)");
    monitor.write_shot_txt((wd.dir() / "shot.txt").string(), wall, reason_name(guard.reason()));
    note("wall прогона " + fmtg(wall, 4) + " с");
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    bool have_config = false;
    Config config = Config::A;
    bool keep_workdir = false;
    bool body_isolation = true; // --no-body-isolation: ДИАГНОСТИКА (ТЗ Guo 2019 §8.4 п. 3), не режим приёмки
    double time_limit = 100.0;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "-h" || a == "--help")
        {
            std::printf("Использование: test_e2e_guo2019 --config=A|1T|1T_A|A_1T|2T|4T [--keep-workdir] [--time-limit=<с>]\n");
            return 0;
        }
        else if (a == "--keep-workdir") keep_workdir = true;
        else if (a == "--no-body-isolation") body_isolation = false;
        else if (a.rfind("--config=", 0) == 0)
        {
            if (!parse_config(a.substr(std::strlen("--config=")), config))
            {
                std::fprintf(stderr, "Неверная конфигурация: %s (допустимо A|1T|1T_A|A_1T|2T|4T)\n", a.c_str());
                return 2;
            }
            have_config = true;
        }
        else if (a.rfind("--time-limit=", 0) == 0)
        {
            try { time_limit = std::stod(a.substr(std::strlen("--time-limit="))); }
            catch (const std::exception&)
            {
                std::fprintf(stderr, "Неверное значение --time-limit: %s\n", a.c_str());
                return 2;
            }
            if (!(time_limit > 0.0))
            {
                std::fprintf(stderr, "--time-limit должен быть положительным\n");
                return 2;
            }
        }
        else
        {
            std::fprintf(stderr, "Неизвестный аргумент: %s\n", a.c_str());
            return 2;
        }
    }
    if (!have_config)
    {
        std::fprintf(stderr, "Нужен --config=A|1T|1T_A|A_1T|2T|4T\n");
        return 2;
    }

    const fs::path cwd = fs::current_path();
    const e2e::DirSnapshot snap = e2e::snapshot_dir(cwd);

    std::printf("=== e2e Guo2019 %s: smoke, V0 = %.0f м/с, t_total = %.0f мкс ===\n", config_name(config), kSmokeV0,
                kSmokeTTotal * 1e6);
    try
    {
        run_main_case(config, keep_workdir, time_limit, body_isolation);
        if (config == Config::T1) test_stop_request();
        test_isolation_without_contact(config);
    }
    catch (const std::exception& ex)
    {
        check(false, std::string("необработанное исключение: ") + ex.what());
    }
    // S10: каталог процесса вернулся, артефакты конвейера в нём не появились.
    std::printf("\n[S10] Рабочий каталог процесса\n");
    check(fs::current_path() == cwd, "текущий каталог после теста тот же: " + cwd.string());
    e2e::check_no_side_effects(cwd, snap);

    std::printf("\n=== e2e Guo2019 %s: %s (провалов: %d) ===\n", config_name(config),
                e2e::failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
