// fabric_checks.cpp
// Диагностические проверочные сцены модели ткани Guo 2019 (шаг 2 диагностики). Не тест: add_test() не добавлен
// (кроме короткого smoke-прогона, см. CMakeLists.txt). Ядро (src/, include/) не меняется: сцены собираются из
// штатных частей -- генератор include/fiber/plain_weave.hpp, load_fiber_file, Clamp, time_integration().
//
// Подкоманды:
//   fabric_checks strip   -- A2: растяжение полоски ткани (или такой же полоски ПРЯМЫХ нитей) с заданной скоростью
//                            деформации; сила по реакции неподвижного зажима, деформация по перемещению подвижного.
//   fabric_checks capstan -- B1: тест Эйлера (капстан): нить огибает неподвижный цилиндр на угол theta,
//                            цилиндр -- частицы Solid (Clamp) либо зажатые тонкие нити (для mu_ff).
//   fabric_checks drum    -- B2: свободные колебания натянутого квадратного лоскута ткани после поперечного толчка
//                            (затухание от трения-проекции; для D_friction нужна сборка с -DSPH_ENERGY_LOG=ON).
//
// ПОДВИЖНЫЙ ЗАЖИМ без правки ядра. Clamp держит бусину неподвижной (проекция ускорения в 0), а validate_clamps
// требует v(0) = 0, поэтому «зажим с постоянной скоростью» собран иначе: группа бусин ТЯЖЁЛАЯ (масса бусины
// увеличена до M, без маски Clamp) и стартует со скоростью v. Изменение скорости за прогон dv = ∫F dt / M;
// M_total = 1e5 кг (по умолчанию) даёт dv/v < 0.5 % даже при 10 1/с и F ~ 2.5 кН (контролируется:
// «dv_heavy_over_v» в summary.txt). Тяжёлая группа жёсткая: все её бусины имеют одну скорость, связи внутри группы не
// растягиваются. Неподвижный зажим -- штатный Clamp, и его реакция (model.clamp_report.fiber_force) -- это сила
// на неподвижном конце (как на силоизмерительном стержне).
//
// Результаты пишутся в --out=<каталог>/<tag>/ (CSV + summary.txt); расчёт идёт во временном каталоге внутри.

#include "e2e_harness.hpp"
#include "guo2019_case.hpp"

#include "fiber/fiber_input.hpp"
#include "fiber/model.hpp"
#include "fiber/plain_weave.hpp"
#include "input.hpp"
#include "int_force.hpp"
#include "interactions/clamp.hpp"
#include "interactions/fiber_bond.hpp"
#include "materials.hpp"
#include "output.hpp"
#include "particle_system.hpp"
#include "run_guard.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace
{

constexpr double kPi = 3.14159265358979323846;

// Материал нити Guo (как tests/guo2019_case.hpp): twaron_ct717, E = 85 ГПа, rho = 1440, nu = 0.3.
constexpr double kD = 0.39e-3;
constexpr double kPitch = 1.18e-3;

// ---------------------------------------------------------------------------
// Командная строка: ключи --name=value
// ---------------------------------------------------------------------------
struct Args
{
    std::string cmd;
    std::map<std::string, std::string> kv;
    bool has(const std::string& k) const { return kv.count(k) != 0; }
    std::string str(const std::string& k, const std::string& def) const
    {
        const auto it = kv.find(k);
        return it == kv.end() ? def : it->second;
    }
    double num(const std::string& k, double def) const
    {
        const auto it = kv.find(k);
        if (it == kv.end()) return def;
        std::size_t pos = 0;
        const double v = std::stod(it->second, &pos);
        if (pos != it->second.size() || !std::isfinite(v)) throw std::runtime_error("неверное значение --" + k + "=" + it->second);
        return v;
    }
};

Args parse(int argc, char** argv)
{
    Args a;
    if (argc < 2) throw std::runtime_error("нужна подкоманда: strip | capstan | drum");
    a.cmd = argv[1];
    for (int i = 2; i < argc; ++i)
    {
        const std::string s = argv[i];
        if (s.rfind("--", 0) != 0) throw std::runtime_error("неизвестный аргумент: " + s);
        const std::size_t eq = s.find('=');
        if (eq == std::string::npos) a.kv[s.substr(2)] = "1";
        else a.kv[s.substr(2, eq - 2)] = s.substr(eq + 1);
    }
    return a;
}

void print_usage()
{
    std::printf(
        "Использование: fabric_checks <подкоманда> [--ключ=значение ...]\n"
        "  strip   --kind=weave|straight --amp-d=0.6 [--seg=2] --material=elastic|plastic --rate=300 [--mu=0.1]\n"
        "          [--eps-max=0.14] [--nx=31 --ny=10] [--fib-inelastic=<файл.fib>] --out=<каталог> --tag=<имя>\n"
        "  capstan --drum=solid|fiber --theta-deg=180 --mu=0.3 [--T1=<Н>] [--v=0.05] [--t-end=<с>] --out=<каталог> --tag=<имя>\n"
        "  drum    --mu=0.1 --stretch=0.07 --v0=3 [--n=18] [--t-end=1e-3] [--energy-log=N] --out=<каталог> --tag=<имя>\n"
        "  общие:  [--threads=N] [--time-limit=<с>]\n");
}

// ---------------------------------------------------------------------------
// Общая обвязка прогона
// ---------------------------------------------------------------------------
struct Scene
{
    SimulationConfig cfg;
    fiber::Model model;
    std::size_t ntotal = 0;
    double delta_x = 0.0;
};

SimulationConfig make_cfg(std::size_t ntotal)
{
    const fs::path p = fs::path(SPH_CASES_DIR) / "guo2019" / "experiment" / "initial_conditions.txt";
    if (!fs::exists(p)) throw std::runtime_error("нет файла конфигурации " + p.string());
    SimulationConfig cfg;
    cfg.apply_defaults();
    {
        e2e::ScopedSilence silence(true);
        parse_sph_config(cfg, p.string());
    }
    guo2019::apply_guo_materials(cfg.mats);
    cfg.sph_same_body_only = true;
    cfg.maxn = std::max<std::size_t>(static_cast<std::size_t>(1.4 * static_cast<double>(ntotal)) + 1, 1);
    cfg.max_interaction = 100 * cfg.maxn;
    cfg.moni_particle = 0;
    cfg.save_t = 1.0e30;
    cfg.print_step = 1000000000;
    return cfg;
}

// Сцена без Solid: фиктивная Solid-частица, ntotal = 0 (как tests/clamp_e2e_fiber.hpp).
void init_fiber_only(Scene& sc)
{
    sc.cfg = make_cfg(0);
    sc.ntotal = 0;
    sc.model.solid.resize(sc.cfg.maxn);
    for (std::size_t i = 0; i < sc.cfg.maxn; ++i) sc.model.solid.inzializ_index(i);
    sc.model.fiber.resize(0);
    sc.model.topology.resize(0);
}

void load_fib_text(Scene& sc, const std::string& text, const fs::path& dir)
{
    const fs::path p = dir / "scene.fib";
    {
        std::ofstream f(p);
        f << text;
    }
    fiber::load_fiber_file(p.string(), sc.model);
}

// Директивы неупругости (PLASTIC/RATE/VISCO/FAILURE) вставляются сразу после первой строки "MATERIAL 0".
std::string insert_inelastic(const std::string& fib, const std::string& directives_file)
{
    std::ifstream in(directives_file);
    if (!in) throw std::runtime_error("не открыть " + directives_file);
    std::string body, ln;
    while (std::getline(in, ln))
    {
        if (!ln.empty() && ln.back() == '\r') ln.pop_back();
        body += ln + "\n";
    }
    const std::string text = "\n" + fib; // MATERIAL 0 может стоять в первой строке
    const std::size_t m = text.find("\nMATERIAL 0 ");
    if (m == std::string::npos) throw std::runtime_error("в .fib нет строки MATERIAL 0");
    const std::size_t e = text.find('\n', m + 1);
    return text.substr(1, e) + body + text.substr(e + 1);
}

struct RunInfo
{
    std::size_t steps = 0;
    double wall_s = 0.0, dt_first = 0.0, dt_min = 0.0, t_reached = 0.0;
    bool failed = false;
};

// Один непрерывный прогон; рабочий каталог -- workdir (log.txt, кадр шага 0 пишет конвейер в CWD).
RunInfo run_scene(Scene& sc, double t_end, double time_limit, const fs::path& workdir, std::function<void(double, double, RunGuard&)> obs,
                  std::size_t energy_every = 0)
{
    RunInfo ri;
    sc.cfg.t_total = t_end;
    sc.cfg.save_t = 1.0e30;
    sc.cfg.print_step = 1000000000;
#ifdef SPH_ENERGY_LOG
    sc.cfg.energy_log_every = energy_every;
    sc.cfg.energy_log_file = (workdir / "energy_log.csv").string();
#else
    (void)energy_every;
#endif
    fs::create_directories(workdir);
    const fs::path old = fs::current_path();
    fs::current_path(workdir);
    RunGuard guard(time_limit, 0.0, 0);
    guard.set_step_observer([&](double dt, double t) { obs(dt, t, guard); });
    double dt = sc.ntotal > 0 ? compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal) : 1.0e-9;
    double time = 0.0;
    const auto w0 = std::chrono::steady_clock::now();
    try
    {
        e2e::ScopedSilence silence(true);
        time_integration(sc.cfg, sc.model, sc.ntotal, dt, time, sc.delta_x, &guard);
    }
    catch (...)
    {
        fs::current_path(old);
        throw;
    }
    ri.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
    fs::current_path(old);
    ri.steps = guard.steps();
    ri.dt_first = guard.dt_first();
    ri.dt_min = guard.dt_min();
    ri.t_reached = time;
    ri.failed = guard.failed();
    // служебные кадры шага 0 не нужны
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(workdir, ec))
    {
        const std::string n = e.path().filename().string();
        if (n.rfind("sph_output_", 0) == 0 || n.rfind("fiber_output_", 0) == 0 || n == "sph_gp.csv") fs::remove(e.path(), ec);
    }
    return ri;
}

void request_stop(RunGuard& g) { g.set_step_limit(g.steps() + 1); }

double area_of(const fiber::FiberMaterial& m) { return interactions::fiber_cross_section_area(m); }

// Растяжение связи (l - l0)/l0.
double bond_strain(const fiber::Model& m, const fiber::FiberBond& b)
{
    double l2 = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        const double r = m.fiber.x[d][b.j] - m.fiber.x[d][b.i];
        l2 += r * r;
    }
    return (std::sqrt(l2) - b.l0) / b.l0;
}

// Материал нити Guo по умолчанию -> параметры генератора.
void set_guo_material(fiber::PlainWeaveParams& q, bool plastic, double mu)
{
    q.E = guo2019::kFabricE;
    q.nu = guo2019::kFabricNu;
    q.rho = guo2019::kFabricRho;
    q.eps_break = plastic ? 0.06 : guo2019::kFabricEpsBreak;
    q.mu_ff = mu;
    q.mu_fs = mu;
    q.mu_fs_set = true;
}

std::string material_line(double E, double d, double rho, double eps_break, double mu, int id)
{
    char buf[256];
    std::snprintf(buf, sizeof buf, "MATERIAL %d %.17g 0.3 %.17g %.17g %.17g %.17g %.17g\n", id, E, d, rho, eps_break, mu, mu);
    return buf;
}

void write_kv(const fs::path& p, const std::vector<std::pair<std::string, std::string>>& kv)
{
    std::ofstream f(p);
    for (const auto& e : kv) f << e.first << "=" << e.second << "\n";
}

std::string g9(double v)
{
    char b[64];
    std::snprintf(b, sizeof b, "%.9g", v);
    return b;
}

void set_threads(const Args& a)
{
#ifdef _OPENMP
    if (a.has("threads")) omp_set_num_threads(static_cast<int>(a.num("threads", 1)));
#else
    (void)a;
#endif
}

// ===========================================================================
// A2: полоска
// ===========================================================================
int cmd_strip(const Args& a)
{
    set_threads(a);
    const bool weave = a.str("kind", "weave") == "weave";
    const double amp_d = a.num("amp-d", 0.6);
    const int seg = static_cast<int>(a.num("seg", 2)); // сегментов на полуволну (base = 2)
    const bool plastic = a.str("material", "elastic") == "plastic";
    const double rate = a.num("rate", 300.0);
    const double mu = a.num("mu", 0.1);
    const double eps_max = a.num("eps-max", 0.14);
    const int nx = static_cast<int>(a.num("nx", 31)); // нитей x (поперечных), задают длину зажатых нитей y
    const int ny = static_cast<int>(a.num("ny", 10)); // зажатых нитей y (нагружаемых)
    const double M_total = a.num("M-total", 1.0e5); // кг на всю тяжёлую группу
    const fs::path out = fs::absolute(a.str("out", "out")) / a.str("tag", "strip");
    const std::string inel = a.str("fib-inelastic", "");
    if (plastic && inel.empty()) throw std::runtime_error("--material=plastic требует --fib-inelastic=<файл>");
    fs::create_directories(out);

    // ---- геометрия и материал
    const double eps_break = plastic ? 0.06 : guo2019::kFabricEpsBreak;
    std::string fib;
    const double e_over = 3.0e-3; // выступ зажатых нитей за последнее пересечение
    const double y_last = 0.5 * static_cast<double>(nx - 1) * kPitch;
    const double c_e = y_last + 1.0e-3;
    fiber::PlainWeave W;
    if (weave)
    {
        fiber::PlainWeaveParams q;
        fiber::apply_plain_weave_preset("twaron_ct717", q);
        q.amp = amp_d * q.d;
        q.seg = seg;
        q.nx = nx;
        q.ny = ny;
        q.overhang_y = e_over;
        q.clamp = fiber::PlainWeaveClamp::Y;
        q.clamp_edge = c_e;
        q.clamp_edge_set = true;
        set_guo_material(q, plastic, mu);
        q.cx = q.cy = q.cz = 0.0;
        W = fiber::build_plain_weave(q);
        std::ostringstream os;
        fiber::write_plain_weave_fib(os, W, "fabric_checks strip");
        fib = os.str();
    }
    else
    {
        // прямые нити: ny нитей вдоль y на тех же x, что у основы утка; шаг бусин ~ как у ткани (p/2)
        std::ostringstream os;
        os.precision(17);
        os << material_line(guo2019::kFabricE, kD, guo2019::kFabricRho, eps_break, mu, 0);
        const double y_end = y_last + e_over;
        const int nb = static_cast<int>(std::lround(2.0 * y_end / (0.5 * kPitch))) + 1;
        for (int i = 0; i < ny; ++i)
        {
            const double x = (static_cast<double>(i) - 0.5 * static_cast<double>(ny - 1)) * kPitch;
            os << "FIBER 0 " << nb << "\n";
            for (int k = 0; k < nb; ++k)
                os << x << " " << (-y_end + 2.0 * y_end * static_cast<double>(k) / static_cast<double>(nb - 1)) << " 0\n";
        }
        os << "CLAMP_BOX -1 " << -1.0 << " -1 1 " << -c_e << " 1\n";
        os << "CLAMP_BOX -1 " << c_e << " -1 1 1 1\n";
        fib = os.str();
    }
    if (plastic) fib = insert_inelastic(fib, inel);

    Scene sc;
    init_fiber_only(sc);
    const fs::path work = out / "run";
    fs::create_directories(work);
    load_fib_text(sc, fib, work);
    fiber::Model& m = sc.model;
    const std::size_t nf = m.fiber.size();
    if (m.clamps.fiber_mask.size() != nf) throw std::runtime_error("нет зажимов в сцене полоски");

    // ---- подвижный конец (y > 0): тяжёлые бусины со скоростью v, без маски Clamp
    std::vector<std::size_t> hi, lo;
    for (std::size_t k = 0; k < nf; ++k)
        if (m.clamps.fiber_mask[k] != 0)
        {
            if (m.fiber.x[1][k] > 0.0) hi.push_back(k);
            else lo.push_back(k);
        }
    if (hi.empty() || lo.empty()) throw std::runtime_error("не нашлись зажатые бусины с обоих концов");
    double y_hi_in = 1e9, y_lo_in = -1e9;
    for (std::size_t k : hi) y_hi_in = std::min(y_hi_in, m.fiber.x[1][k]);
    for (std::size_t k : lo) y_lo_in = std::max(y_lo_in, m.fiber.x[1][k]);
    const double L0 = y_hi_in - y_lo_in;
    const double v = rate * L0;
    const double Mb = M_total / static_cast<double>(hi.size());
    std::vector<double> y0_hi(hi.size());
    for (std::size_t q2 = 0; q2 < hi.size(); ++q2)
    {
        const std::size_t k = hi[q2];
        m.clamps.fiber_mask[k] = 0;
        m.fiber.mass[k] = Mb;
        m.fiber.vx[1][k] = v;
        y0_hi[q2] = m.fiber.x[1][k];
    }
    // длина нитей на самом деле стыкуется между зажимами: число нагружаемых нитей = ny
    const fiber::FiberMaterial& mat0 = m.materials.at(0);
    const double EA = mat0.E * area_of(mat0);
    const double nEA = static_cast<double>(ny) * EA; // жёсткость прямых нитей по деформации, Н (F = nEA * eps)

    const double t_end = eps_max / rate;
    const int nbin = 3000;
    const double dtb = t_end / nbin;
    struct Row { double t, eps, F; std::size_t nb; };
    std::vector<Row> rows;
    double acc_F = 0.0, acc_eps = 0.0, acc_t = 0.0;
    std::size_t acc_n = 0;
    double next_bin = dtb;
    double Fmax = 0.0, eps_Fmax = 0.0, eps_first_break = -1.0;
    const std::size_t bonds0 = m.topology.bonds.size();
    std::size_t nb_prev = bonds0;
    double vh_end = v;
    auto obs = [&](double /*dt*/, double t, RunGuard& g) {
        double u = 0.0;
        for (std::size_t q2 = 0; q2 < hi.size(); ++q2) u += m.fiber.x[1][hi[q2]] - y0_hi[q2];
        u /= static_cast<double>(hi.size());
        const double eps = u / L0;
        const double F = -m.clamp_report.fiber_force[1];
        acc_F += F;
        acc_eps += eps;
        acc_t += t;
        ++acc_n;
        const std::size_t nbn = m.topology.bonds.size();
        if (nbn < nb_prev && eps_first_break < 0.0) eps_first_break = eps;
        nb_prev = nbn;
        vh_end = m.fiber.vx[1][hi[0]];
        if (t >= next_bin)
        {
            Row r{acc_t / static_cast<double>(acc_n), acc_eps / static_cast<double>(acc_n), acc_F / static_cast<double>(acc_n), nbn};
            rows.push_back(r);
            if (r.F > Fmax)
            {
                Fmax = r.F;
                eps_Fmax = r.eps;
            }
            acc_F = acc_eps = acc_t = 0.0;
            acc_n = 0;
            next_bin += dtb;
            // всё порвалось: сила упала ниже 2 % максимума после заметного растяжения
            if (rows.size() > 50 && Fmax > 2.0e-4 * nEA && r.F < 0.02 * Fmax && eps > 0.03) request_stop(g);
        }
    };
    std::printf("[strip %s] kind=%s amp/d=%.2f %s rate=%g 1/c L0=%.4f мм v=%.3f м/с, нитей %d, бусин %zu, связей %zu, nEA=%.4g Н\n",
                a.str("tag", "strip").c_str(), weave ? "weave" : "straight", amp_d, plastic ? "plastic" : "elastic", rate, L0 * 1e3, v, ny,
                nf, bonds0, nEA);
    std::fflush(stdout);
    const RunInfo ri = run_scene(sc, t_end, a.num("time-limit", 0.0), work, obs);
    {
        std::ofstream f(out / "strip.csv");
        f.precision(10);
        f << "t,eps,F,F_over_nEA,bonds\n";
        for (const Row& r : rows) f << r.t << "," << r.eps << "," << r.F << "," << r.F / nEA << "," << r.nb << "\n";
    }
    write_kv(out / "summary.txt",
             {{"kind", weave ? "weave" : "straight"},
              {"amp_d", g9(amp_d)},
              {"seg", std::to_string(seg)},
              {"material", plastic ? "plastic" : "elastic"},
              {"rate", g9(rate)},
              {"mu", g9(mu)},
              {"nx", std::to_string(nx)},
              {"n_loaded", std::to_string(ny)},
              {"L0", g9(L0)},
              {"v", g9(v)},
              {"EA", g9(EA)},
              {"nEA", g9(nEA)},
              {"beads", std::to_string(nf)},
              {"bonds0", std::to_string(bonds0)},
              {"steps", std::to_string(ri.steps)},
              {"dt_first", g9(ri.dt_first)},
              {"dt_min", g9(ri.dt_min)},
              {"wall_s", g9(ri.wall_s)},
              {"t_reached", g9(ri.t_reached)},
              {"t_end", g9(t_end)},
              {"F_max", g9(Fmax)},
              {"eps_at_Fmax", g9(eps_Fmax)},
              {"eps_first_break", g9(eps_first_break)},
              {"dv_heavy_over_v", g9(std::fabs(vh_end - v) / v)},
              {"failed_guard", ri.failed ? "1" : "0"}});
    std::printf("[strip %s] шагов %zu, wall %.1f с, dt=%.3g с, F_max=%.4g Н при eps=%.4f, первый разрыв eps=%.4f, dv/v=%.2e\n",
                a.str("tag", "strip").c_str(), ri.steps, ri.wall_s, ri.dt_first, Fmax, eps_Fmax, eps_first_break,
                std::fabs(vh_end - v) / v);
    std::error_code ec;
    fs::remove(work / "scene.fib", ec);
    return 0;
}

// ===========================================================================
// B1: капстан
// ===========================================================================
int cmd_capstan(const Args& a)
{
    set_threads(a);
    const bool solid = a.str("drum", "solid") == "solid";
    const double theta = a.num("theta-deg", 180.0) * kPi / 180.0;
    const double mu = a.num("mu", 0.3);
    const double R_c = a.num("R-c", 2.0e-3);     // радиус цилиндра (поверхность решётки/нитей), м
    const double ds = a.num("ds", 0.2e-3);       // целевой шаг бусин нити
    const double dxs = a.num("dx-s", 0.1e-3);    // шаг решётки Solid / толщина тонкой нити барабана
    const double Ltail = a.num("tail", 6.0e-3);
    const double v = a.num("v", 0.1);
    const double M_heavy = a.num("M-heavy", 1.0);
    const fs::path out = fs::absolute(a.str("out", "out")) / a.str("tag", "capstan");
    fs::create_directories(out);
    const double e_theory = std::exp(mu * theta);
    double T1 = a.has("T1") ? a.num("T1", 1.0) : std::min(8.0, std::max(1.0, 100.0 / e_theory));

    const double E = guo2019::kFabricE, rho = guo2019::kFabricRho;
    const double A = kPi * kD * kD / 4.0;
    const double eps1 = T1 / (E * A);

    // ---- радиус траектории центра нити
    // R_c -- радиус контактной поверхности цилиндра: центры частиц Solid / оси тонких нитей лежат на R_c - dx/2
    const double R_t = R_c + 0.5 * kD;

    // ---- точки нити: вход по касательной, дуга, выход по касательной; плоскость z = z(alpha) (подъём при theta > pi)
    const double a0 = -0.5 * kPi;
    const double rise = theta > kPi + 1e-9 ? 0.8e-3 : 0.0;
    const std::size_t nA = std::max<std::size_t>(2, static_cast<std::size_t>(std::lround(theta * R_t / ds)));
    const double dA = theta / static_cast<double>(nA);
    const double chord = 2.0 * R_t * std::sin(0.5 * dA);
    std::vector<std::array<double, 3>> pts;
    auto zof = [&](double al) { return rise * ((al - a0) / theta - 0.5); };
    const std::size_t nT = std::max<std::size_t>(2, static_cast<std::size_t>(std::lround(Ltail / chord)));
    {
        const double tx = -std::sin(a0), ty = std::cos(a0);
        const double sx = R_t * std::cos(a0), sy = R_t * std::sin(a0);
        for (std::size_t k = nT; k >= 1; --k) pts.push_back({sx - tx * chord * static_cast<double>(k), sy - ty * chord * static_cast<double>(k), zof(a0)});
        for (std::size_t j = 0; j <= nA; ++j)
        {
            const double al = a0 + dA * static_cast<double>(j);
            pts.push_back({R_t * std::cos(al), R_t * std::sin(al), zof(al)});
        }
        const double a1 = a0 + theta;
        const double ex = R_t * std::cos(a1), ey = R_t * std::sin(a1);
        const double ux = -std::sin(a1), uy = std::cos(a1);
        for (std::size_t k = 1; k <= nT; ++k) pts.push_back({ex + ux * chord * static_cast<double>(k), ey + uy * chord * static_cast<double>(k), zof(a1)});
    }
    const std::size_t nbeads_thread = pts.size();
    const double a_end = a0 + theta;
    const double dir_end[3] = {-std::sin(a_end), std::cos(a_end), 0.0};

    // ---- .fib: материал нити 0, (материал барабана 1), нить, (нити барабана)
    std::ostringstream os;
    os.precision(17);
    os << material_line(E, kD, rho, 0.03, mu, 0);
    const double zspan = 0.5 * rise + 0.3e-3; // половина высоты барабана по z (контакт при расстоянии < 0.25 мм)
    std::size_t n_drum_fibers = 0;
    if (!solid)
    {
        os << material_line(E, dxs, rho, 0.03, mu, 1);
    }
    os << "FIBER 0 " << pts.size() << "\n";
    for (const auto& p : pts) os << p[0] << " " << p[1] << " " << p[2] << "\n";
    const double r_drum_y = R_c - 0.5 * dxs; // оси тонких нитей барабана
    if (!solid)
    {
        const std::size_t Ny = static_cast<std::size_t>(std::lround(2.0 * kPi * r_drum_y / dxs));
        const std::size_t nz = static_cast<std::size_t>(std::lround(2.0 * zspan / 0.3e-3)) + 1;
        for (std::size_t i = 0; i < Ny; ++i)
        {
            const double al = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(Ny);
            os << "FIBER 1 " << nz << "\n";
            for (std::size_t k = 0; k < nz; ++k)
                os << r_drum_y * std::cos(al) << " " << r_drum_y * std::sin(al) << " "
                   << (-zspan + 2.0 * zspan * static_cast<double>(k) / static_cast<double>(nz - 1)) << "\n";
            ++n_drum_fibers;
        }
    }

    // ---- Solid-цилиндр: слои радиусов R_c - 0.5 dx, ... (3 слоя), кольца по z
    Scene sc;
    std::size_t n_solid = 0;
    struct SP { double x, y, z; };
    std::vector<SP> sp;
    if (solid)
    {
        const std::size_t nzl = static_cast<std::size_t>(std::lround(2.0 * zspan / dxs)) + 1;
        for (int layer = 0; layer < 1; ++layer)
        {
            const double r = R_c - 0.5 * dxs - static_cast<double>(layer) * dxs;
            if (r <= 0.0) break;
            const std::size_t Nr = static_cast<std::size_t>(std::lround(2.0 * kPi * r / dxs));
            for (std::size_t i = 0; i < Nr; ++i)
            {
                const double al = 2.0 * kPi * (static_cast<double>(i) + 0.5 * static_cast<double>(layer % 2)) / static_cast<double>(Nr);
                for (std::size_t k = 0; k < nzl; ++k)
                    sp.push_back({r * std::cos(al), r * std::sin(al), -0.5 * (static_cast<double>(nzl - 1)) * dxs + static_cast<double>(k) * dxs});
            }
        }
        n_solid = sp.size();
    }
    sc.ntotal = n_solid;
    sc.cfg = make_cfg(n_solid);
    sc.delta_x = solid ? dxs : 0.0;
    // Цилиндр неподвижен (Clamp), от его материала нужна только жёсткость контакта. Сталь Guo (aa = 1.45e11) даёт
    // c ~ 4.3 км/с и CFL-шаг по hsml = 2 dx в ~10 раз меньше шага нитей. Поэтому aa и G уменьшены в 18 и 10 раз:
    // E_цилиндра ~ 2e10 Па (контактное проникание при N ~ 1 Н порядка 2 мкм), c ~ 0.4 км/с, шаг задают нити.
    sc.cfg.mats.aa[guo2019::kSlotSteel] = 1.44e9;
    sc.cfg.mats.G[guo2019::kSlotSteel] = 8.0e9;
    {
        ParticleSystem& s = sc.model.solid;
        s.resize(sc.cfg.maxn);
        for (std::size_t i = 0; i < sc.cfg.maxn; ++i) s.inzializ_index(i);
        const double rho_s = guo2019::kRhoSteelRef;
        for (std::size_t i = 0; i < n_solid; ++i)
        {
            s.x[0][i] = sp[i].x;
            s.x[1][i] = sp[i].y;
            s.x[2][i] = sp[i].z;
            s.itype[i] = guo2019::kBallItype;
            s.body[i] = static_cast<int>(i) + 1; // у каждой частицы своё тело: SPH-суммы не считаются (sph_same_body_only), цилиндр неподвижен
            s.rho[i] = rho_s;
            s.rho_0[i] = rho_s;
            s.mass[i] = rho_s * dxs * dxs * dxs;
            s.hsml[i] = 2.0 * dxs;
        }
        for (std::size_t i = 0; i < n_solid; ++i) urs(sc.cfg.mats, s, i);
    }
    sc.model.fiber.resize(0);
    sc.model.topology.resize(0);
    const fs::path work = out / "run";
    fs::create_directories(work);
    load_fib_text(sc, os.str(), work);
    fiber::Model& m = sc.model;
    if (m.fiber.size() < nbeads_thread) throw std::runtime_error("капстан: неверное число бусин");

    // ---- начальное натяжение T1: l0 = l/(1+eps1) у всех связей нити (барабан-нити остаются ненапряжёнными)
    for (auto& b : m.topology.bonds)
        if (b.i < nbeads_thread && b.j < nbeads_thread)
        {
            double l2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double r = m.fiber.x[d][b.j] - m.fiber.x[d][b.i];
                l2 += r * r;
            }
            b.l0 = std::sqrt(l2) / (1.0 + eps1);
            b.l_break = (1.0 + m.materials[0].eps_break) * b.l0;
        }
    // ---- зажимы: первая бусина нити (T1), весь барабан-нити и Solid-цилиндр
    m.clamps.fiber_mask.assign(m.fiber.size(), 0);
    m.clamps.fiber_mask[0] = interactions::kClampAll;
    for (std::size_t k = nbeads_thread; k < m.fiber.size(); ++k) m.clamps.fiber_mask[k] = interactions::kClampAll;
    m.clamps.solid_mask.assign(sc.ntotal, 0);
    for (std::size_t i = 0; i < n_solid; ++i) m.clamps.solid_mask[i] = interactions::kClampAll;
    if (sc.ntotal == 0) m.clamps.solid_mask.clear();
    // ---- тяжёлая последняя бусина со скоростью v вдоль касательной выхода
    const std::size_t kh = nbeads_thread - 1;
    m.fiber.mass[kh] = M_heavy;

    // ---- наблюдатель: T1 (связи у зажима), T2 (связи у тяжёлой бусины), бины по времени
    const double t_end = a.num("t-end", 6.0e-4);
    const int nbin = 1500;
    const double dtb = t_end / nbin;
    const std::size_t nav = 5;
    struct Row { double t, T1, T2, Tarc_mid, pen; };
    std::vector<Row> rows;
    std::vector<std::size_t> bidx_lo, bidx_hi; // индексы связей
    {
        for (std::size_t q2 = 0; q2 < m.topology.bonds.size(); ++q2)
        {
            const auto& b = m.topology.bonds[q2];
            for (std::size_t s = 1; s <= nav; ++s)
            {
                if (b.i == s && b.j == s + 1) bidx_lo.push_back(q2);
                if (b.i == kh - 1 - s && b.j == kh - s) bidx_hi.push_back(q2);
            }
        }
    }
    auto tension = [&](const std::vector<std::size_t>& ids) {
        double s = 0.0;
        for (std::size_t q2 : ids) s += E * A * std::max(bond_strain(m, m.topology.bonds[q2]), 0.0);
        return s / static_cast<double>(ids.size());
    };
    // ---- фаза успокоения (кинетическое демпфирование): предварительное натяжение T1 неравновесно (контакт с цилиндром
    //      податлив: проникание в несколько мкм снимает заметную часть натяжения), и без успокоения нить звенела бы.
    //      Наблюдатель обнуляет скорости бусин нити, когда кинетическая энергия падает (пик пройден), -- это только
    //      подготовка начального состояния; фаза тяги идёт без вмешательства в состояние.
    {
        const double t_settle = a.num("t-settle", 1.0e-4);
        double ke_prev = 0.0;
        auto settle_obs = [&](double /*dt*/, double /*t*/, RunGuard& /*g*/) {
            double ke = 0.0;
            for (std::size_t k = 0; k < nbeads_thread; ++k)
                for (std::size_t d = 0; d < 3; ++d) ke += 0.5 * m.fiber.mass[k] * m.fiber.vx[d][k] * m.fiber.vx[d][k];
            if (ke < ke_prev)
            {
                for (std::size_t k = 0; k < nbeads_thread; ++k)
                    for (std::size_t d = 0; d < 3; ++d) m.fiber.vx[d][k] = 0.0;
                ke = 0.0;
            }
            ke_prev = ke;
        };
        if (t_settle > 0.0)
        {
            // трение на время успокоения выключено: иначе прилипание фиксирует неравновесное распределение натяжения
            std::vector<std::pair<double, double>> mu_save;
            for (auto& mt : m.materials) { mu_save.emplace_back(mt.mu_ff, mt.mu_fs); mt.mu_ff = 0.0; mt.mu_fs = 0.0; }
            const RunInfo rs = run_scene(sc, t_settle, a.num("time-limit", 0.0), work, settle_obs);
            for (std::size_t q3 = 0; q3 < m.materials.size(); ++q3) { m.materials[q3].mu_ff = mu_save[q3].first; m.materials[q3].mu_fs = mu_save[q3].second; }
            for (std::size_t k = 0; k < nbeads_thread; ++k)
                for (std::size_t d = 0; d < 3; ++d) m.fiber.vx[d][k] = 0.0;
            double pen = 0.0;
            std::size_t npen = 0;
            for (std::size_t k = nT; k + nT < nbeads_thread; ++k)
            {
                pen += R_t - std::hypot(m.fiber.x[0][k], m.fiber.x[1][k]);
                ++npen;
            }
            {
                double sa = 0.0, st = 0.0;
                std::size_t na = 0, nt = 0;
                for (const auto& b : m.topology.bonds)
                    if (b.i < nbeads_thread && b.j < nbeads_thread)
                    {
                        const double tt = E * A * bond_strain(m, b);
                        if (b.i >= nT && b.j <= nT + nA) { sa += tt; ++na; }
                        else { st += tt; ++nt; }
                    }
                std::printf("[capstan %s] после успокоения: среднее натяжение (со знаком) дуги %.3f Н, хвостов %.3f Н\n", a.str("tag", "capstan").c_str(),
                            sa / static_cast<double>(std::max<std::size_t>(na, 1)), st / static_cast<double>(std::max<std::size_t>(nt, 1)));
            }
            std::printf("[capstan %s] успокоение: %zu шагов, wall %.1f с, среднее проникание бусин дуги %.3g мкм, T1=%.3f T2=%.3f Н\n",
                        a.str("tag", "capstan").c_str(), rs.steps, rs.wall_s, 1e6 * pen / static_cast<double>(std::max<std::size_t>(npen, 1)),
                        tension(bidx_lo), tension(bidx_hi));
        }
    }
    for (std::size_t d = 0; d < 3; ++d) m.fiber.vx[d][kh] = v * dir_end[d];
    double a1 = 0.0, a2 = 0.0, at = 0.0;
    std::size_t an = 0;
    double next_bin = dtb;
    auto obs = [&](double /*dt*/, double t, RunGuard& /*g*/) {
        a1 += tension(bidx_lo);
        a2 += tension(bidx_hi);
        at += t;
        ++an;
        if (t >= next_bin)
        {
            rows.push_back({at / static_cast<double>(an), a1 / static_cast<double>(an), a2 / static_cast<double>(an), 0.0, 0.0});
            a1 = a2 = at = 0.0;
            an = 0;
            next_bin += dtb;
        }
    };
    std::printf("[capstan %s] drum=%s theta=%.0f deg mu=%.2f exp(mu*theta)=%.4f T1=%.3f Н v=%.3f м/с R_t=%.3f мм бусин нити %zu (ds=%.3f мм), Solid %zu, нитей барабана %zu\n",
                a.str("tag", "capstan").c_str(), solid ? "solid" : "fiber", theta * 180.0 / kPi, mu, e_theory, T1, v, R_t * 1e3, nbeads_thread,
                chord * 1e3, n_solid, n_drum_fibers);
    std::fflush(stdout);
    const RunInfo ri = run_scene(sc, t_end, a.num("time-limit", 0.0), work, obs);
    {
        std::ofstream f(out / "capstan.csv");
        f.precision(10);
        f << "t,T1,T2,ratio\n";
        for (const Row& r : rows) f << r.t << "," << r.T1 << "," << r.T2 << "," << (r.T1 > 0 ? r.T2 / r.T1 : 0.0) << "\n";
    }
    // плато: среднее отношение за последние 25 % прогона
    double sum_r = 0.0;
    std::size_t nr = 0;
    double T1e = 0, T2e = 0;
    for (std::size_t i = rows.size() * 3 / 4; i < rows.size(); ++i)
        if (rows[i].T1 > 0.0)
        {
            sum_r += rows[i].T2 / rows[i].T1;
            T1e += rows[i].T1;
            T2e += rows[i].T2;
            ++nr;
        }
    const double ratio = nr ? sum_r / static_cast<double>(nr) : 0.0;
    if (nr) { T1e /= static_cast<double>(nr); T2e /= static_cast<double>(nr); }
    // дискретная теория для цепочки с шагом dA: (1 + 2 mu sin(dA/2))^nA / (1 - ...)? -- справочно: (1+mu*2 sin(dA/2))^n
    const double e_disc = std::pow(1.0 + 2.0 * mu * std::sin(0.5 * dA), static_cast<double>(nA));
    write_kv(out / "summary.txt",
             {{"drum", solid ? "solid" : "fiber"},
              {"theta_deg", g9(theta * 180.0 / kPi)},
              {"mu", g9(mu)},
              {"exp_mu_theta", g9(e_theory)},
              {"exp_discrete", g9(e_disc)},
              {"ratio_plateau", g9(ratio)},
              {"dev_pct", g9(100.0 * (ratio / e_theory - 1.0))},
              {"T1_plateau", g9(T1e)},
              {"T2_plateau", g9(T2e)},
              {"T1_init", g9(T1)},
              {"v", g9(v)},
              {"R_t", g9(R_t)},
              {"chord", g9(chord)},
              {"dalpha", g9(dA)},
              {"beads_thread", std::to_string(nbeads_thread)},
              {"solid", std::to_string(n_solid)},
              {"drum_fibers", std::to_string(n_drum_fibers)},
              {"steps", std::to_string(ri.steps)},
              {"dt_first", g9(ri.dt_first)},
              {"dt_min", g9(ri.dt_min)},
              {"wall_s", g9(ri.wall_s)},
              {"t_end", g9(t_end)}});
    std::printf("[capstan %s] T2/T1=%.4f (плато) против exp(mu*theta)=%.4f (откл. %+.1f %%), дискр. %.4f; T1=%.3f T2=%.3f Н; шагов %zu, dt=%.3g, wall %.1f с\n",
                a.str("tag", "capstan").c_str(), ratio, e_theory, 100.0 * (ratio / e_theory - 1.0), e_disc, T1e, T2e, ri.steps, ri.dt_first, ri.wall_s);
    std::error_code ec;
    fs::remove(work / "scene.fib", ec);
    return 0;
}

// ===========================================================================
// B2: колебания лоскута
// ===========================================================================
int cmd_drum(const Args& a)
{
    set_threads(a);
    const double mu = a.num("mu", 0.1);
    const double stretch = a.num("stretch", 0.02);
    const double v0 = a.num("v0", 3.0);
    const double r0 = a.num("r0", 3.0e-3);
    const double E_mod = guo2019::kFabricE;
    const int n = static_cast<int>(a.num("n", 18));
    const double amp_d = a.num("amp-d", 0.6);
    const double t_end = a.num("t-end", 1.0e-3);
    const std::size_t elog = static_cast<std::size_t>(a.num("energy-log", 0.0));
    const fs::path out = fs::absolute(a.str("out", "out")) / a.str("tag", "drum");
    fs::create_directories(out);

    fiber::PlainWeaveParams q;
    fiber::apply_plain_weave_preset("twaron_ct717", q);
    q.amp = amp_d * q.d;
    q.seg = 2;
    q.nx = n;
    q.ny = n;
    const double last = 0.5 * static_cast<double>(n - 1) * q.pitch;
    q.overhang_x = q.overhang_y = 0.5 * q.pitch;
    q.clamp = fiber::PlainWeaveClamp::XY;
    q.clamp_edge = last + 0.35e-3;
    q.clamp_edge_set = true;
    set_guo_material(q, false, mu);
    q.cx = q.cy = q.cz = 0.0;
    const fiber::PlainWeave W = fiber::build_plain_weave(q);
    std::ostringstream os;
    fiber::write_plain_weave_fib(os, W, "fabric_checks drum");

    Scene sc;
    init_fiber_only(sc);
    const fs::path work = out / "run";
    fs::create_directories(work);
    load_fib_text(sc, os.str(), work);
    fiber::Model& m = sc.model;
    const std::size_t nf = m.fiber.size();

    // раскладка перекрёстков: основа (вдоль x) j, уток (вдоль y) i; бусина перекрёстка по vertex_bead
    std::vector<std::size_t> first(W.fibers.size());
    {
        std::size_t f0 = 0;
        for (std::size_t k = 0; k < W.fibers.size(); ++k)
        {
            first[k] = f0;
            f0 += W.fibers[k].pts.size();
        }
    }
    struct Cross { std::size_t w, f; };
    std::vector<Cross> cross;
    for (std::size_t j = 0; j < static_cast<std::size_t>(q.nx); ++j)
        for (std::size_t i = 0; i < static_cast<std::size_t>(q.ny); ++i)
            cross.push_back({first[j] + W.fibers[j].vertex_bead.at(i), first[static_cast<std::size_t>(q.nx) + i] + W.fibers[static_cast<std::size_t>(q.nx) + i].vertex_bead.at(j)});

    // ---- предварительное растяжение: x, y всех бусин масштабируются на (1+stretch), z не меняется;
    //      l0 связей остаются исходными, поэтому натяжение = stretch (в проекции)
    for (std::size_t k = 0; k < nf; ++k)
    {
        m.fiber.x[0][k] *= (1.0 + stretch);
        m.fiber.x[1][k] *= (1.0 + stretch);
    }
    // ---- успокоение предварительно растянутого состояния (вязкое демпфирование скоростей, трение выключено): масштабирование
    //      даёт неравновесное состояние, его релаксация без успокоения шла бы с диссипацией на трении и затмила бы колебания
    {
        const double t_settle = a.num("t-settle", 4.0e-4);
        const double gamma = a.num("damp", 0.02); // доля скорости, гасимая на каждом шаге (динамическая релаксация)
        auto settle_obs = [&](double, double, RunGuard&) {
            for (std::size_t k = 0; k < nf; ++k)
                for (std::size_t d = 0; d < 3; ++d) m.fiber.vx[d][k] *= 1.0 - gamma;
        };
        std::vector<std::pair<double, double>> mu_save;
        for (auto& mt : m.materials) { mu_save.emplace_back(mt.mu_ff, mt.mu_fs); mt.mu_ff = 0.0; mt.mu_fs = 0.0; }
        const RunInfo rs = run_scene(sc, t_settle, a.num("time-limit", 0.0), work, settle_obs);
        for (std::size_t q3 = 0; q3 < m.materials.size(); ++q3) { m.materials[q3].mu_ff = mu_save[q3].first; m.materials[q3].mu_fs = mu_save[q3].second; }
        for (std::size_t k = 0; k < nf; ++k)
            for (std::size_t d = 0; d < 3; ++d) m.fiber.vx[d][k] = 0.0;
        double Ub = 0.0;
        for (const auto& bnd : m.topology.bonds)
        {
            const double e = bond_strain(m, bnd);
            if (e > 0.0) Ub += 0.5 * E_mod * area_of(m.materials[0]) * e * e * bnd.l0;
        }
        std::printf("[drum %s] успокоение: %zu шагов, wall %.1f с, U_bond=%.4g Дж\n", a.str("tag", "drum").c_str(), rs.steps, rs.wall_s, Ub);
    }
    // ---- толчок: v_z = v0 exp(-(r/r0)^2) у свободных бусин
    for (std::size_t k = 0; k < nf; ++k)
    {
        if (m.clamps.fiber_mask.size() == nf && m.clamps.fiber_mask[k] != 0) continue;
        const double r2 = m.fiber.x[0][k] * m.fiber.x[0][k] + m.fiber.x[1][k] * m.fiber.x[1][k];
        m.fiber.vx[2][k] = v0 * std::exp(-r2 / (r0 * r0));
    }
    std::size_t kc = 0;
    {
        double best = 1e9;
        for (std::size_t k = 0; k < nf; ++k)
        {
            const double r2 = m.fiber.x[0][k] * m.fiber.x[0][k] + m.fiber.x[1][k] * m.fiber.x[1][k];
            if (r2 < best) { best = r2; kc = k; }
        }
    }
    const double zc0 = m.fiber.x[2][kc];
    std::vector<std::array<double, 3>> rel0(cross.size());
    for (std::size_t c = 0; c < cross.size(); ++c)
        for (std::size_t d = 0; d < 3; ++d) rel0[c][d] = m.fiber.x[d][cross[c].w] - m.fiber.x[d][cross[c].f];
    const fiber::FiberMaterial& mat0 = m.materials.at(0);
    const double EA = mat0.E * area_of(mat0);
    (void)EA;
    double KE0 = 0.0;
    for (std::size_t k = 0; k < nf; ++k)
        for (std::size_t d = 0; d < 3; ++d) KE0 += 0.5 * m.fiber.mass[k] * m.fiber.vx[d][k] * m.fiber.vx[d][k];

    struct Row { double t, KE, zc, zmax, slip_rms, slip_max, frac_slip; };
    std::vector<Row> rows;
    const int nbin = 4000;
    const double dtb = t_end / nbin;
    double next_bin = 0.0;
    auto obs = [&](double /*dt*/, double t, RunGuard& /*g*/) {
        if (t < next_bin) return;
        next_bin += dtb;
        double KE = 0.0;
        double zmax = 0.0;
        for (std::size_t k = 0; k < nf; ++k)
        {
            for (std::size_t d = 0; d < 3; ++d) KE += 0.5 * m.fiber.mass[k] * m.fiber.vx[d][k] * m.fiber.vx[d][k];
        }
        // относительное тангенциальное (в плоскости x-y) смещение пары нитей на перекрёстке
        double s2 = 0.0, smax = 0.0;
        std::size_t nslip = 0;
        for (std::size_t c = 0; c < cross.size(); ++c)
        {
            double dx2 = 0.0;
            for (std::size_t d = 0; d < 2; ++d)
            {
                const double r = (m.fiber.x[d][cross[c].w] - m.fiber.x[d][cross[c].f]) - rel0[c][d];
                dx2 += r * r;
            }
            const double sl = std::sqrt(dx2);
            s2 += sl * sl;
            smax = std::max(smax, sl);
            if (sl > 0.05 * kD) ++nslip;
        }
        zmax = std::fabs(m.fiber.x[2][kc] - zc0);
        rows.push_back({t, KE, m.fiber.x[2][kc] - zc0, zmax, std::sqrt(s2 / static_cast<double>(cross.size())), smax,
                        static_cast<double>(nslip) / static_cast<double>(cross.size())});
    };
    std::printf("[drum %s] mu=%.2f stretch=%.3f v0=%.2f м/с лоскут %dx%d, бусин %zu, перекрёстков %zu, KE0=%.4g Дж\n", a.str("tag", "drum").c_str(),
                mu, stretch, v0, n, n, nf, cross.size(), KE0);
    std::fflush(stdout);
    const RunInfo ri = run_scene(sc, t_end, a.num("time-limit", 0.0), work, obs, elog);
    {
        std::ofstream f(out / "drum.csv");
        f.precision(10);
        f << "t,KE,zc,zabs,slip_rms,slip_max,frac_slip\n";
        for (const Row& r : rows)
            f << r.t << "," << r.KE << "," << r.zc << "," << r.zmax << "," << r.slip_rms << "," << r.slip_max << "," << r.frac_slip << "\n";
    }
    write_kv(out / "summary.txt",
             {{"mu", g9(mu)}, {"stretch", g9(stretch)}, {"v0", g9(v0)}, {"n", std::to_string(n)}, {"beads", std::to_string(nf)},
              {"crossings", std::to_string(cross.size())}, {"KE0", g9(KE0)}, {"steps", std::to_string(ri.steps)},
              {"dt_first", g9(ri.dt_first)}, {"dt_min", g9(ri.dt_min)}, {"wall_s", g9(ri.wall_s)}, {"t_end", g9(t_end)},
              {"t_reached", g9(ri.t_reached)}});
    std::printf("[drum %s] шагов %zu, dt=%.3g с, wall %.1f с, t=%.3g\n", a.str("tag", "drum").c_str(), ri.steps, ri.dt_first, ri.wall_s,
                ri.t_reached);
    std::error_code ec;
    fs::remove(work / "scene.fib", ec);
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))
        {
            print_usage();
            return 0;
        }
        const Args a = parse(argc, argv);
        if (a.cmd == "strip") return cmd_strip(a);
        if (a.cmd == "capstan") return cmd_capstan(a);
        if (a.cmd == "drum") return cmd_drum(a);
        std::fprintf(stderr, "Неизвестная подкоманда: %s\n", a.cmd.c_str());
        print_usage();
        return 2;
    }
    catch (const std::exception& ex)
    {
        std::fprintf(stderr, "ОШИБКА: %s\n", ex.what());
        return 1;
    }
}
