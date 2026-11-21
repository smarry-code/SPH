#pragma once
// ceramic_shot.hpp
// Общий драйвер экспериментальных утилит для 2D-опытов с керамикой (bc2_shot, bc6_shot).
// По образцу guo2019_shot: постоянный каталог прогона, серия по Δx и типу контакта,
// число потоков, сторож времени, чекпоинты и рестарт, итог в shot.txt, сводка runs.csv.
//
// Кейс отдаёт CaseHooks: построение геометрии (номер частицы-монитора, ключи для shot.txt)
// и разбор кривой монитора против эталона (ключи для shot.txt, подпись графика).
// Порядок подготовки прогона -- тот же, что у сквозных тестов (prepare_run() в
// tests/test_e2e_bc6.cpp): умолчания -> разбор конфига -> геометрия -> запас под
// зеркальные частицы -> УРС -> начальный шаг по Куранту.
//
// Счёт идёт через time_integration(cfg, fiber::Model&, ...) без нитей: это та же реализация,
// что у перегрузки для ParticleSystem (src/time_integration.cpp), но с чекпоинтами.
//
// Каталог прогона: <out>/<кейс>_dx<мкм>_c<тип>[_<tag>]/ -- в нём же CWD прогона, поэтому
// файлы конвейера (log.txt, кадры sph_output_*.txt) пишутся туда.
//   log.txt                      монитор конвейера, строка на шаг (дозапись)
//   sph_output_*.txt             кадры раз в save_t
//   curve.txt                    кривая монитора без повторов после рестарта
//   compare.svg                  расчёт против приборного эталона
//   shot.txt                     итог: параметры, исход, метрики сверки (ключ=значение)
//   initial_conditions.used.txt  копия применённого файла параметров решателя
//   scene_id.txt, checkpoint.ckpt[.prev]   при --checkpoint-every
// Сводка серии -- <out>/runs.csv (строка на прогон, дозапись).
//
// Коды выхода: 0 -- все прогоны дошли до t_total или до --max-steps; 1 -- хотя бы один оборван
// (сторож времени/застоя, исключение); 2 -- ошибка ключей. Сверка с эталоном воротами НЕ является: это отчёт.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <iostream>
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
#include "input.hpp"
#include "int_force.hpp"          // urs()
#include "time_integration.hpp"
#include "run_guard.hpp"
#include "fiber/model.hpp"

#include "e2e_harness.hpp"
#include "e2e_curve.hpp"
#include "e2e_plot.hpp"

namespace ceramic
{

namespace fs = std::filesystem;
using KV = std::vector<std::pair<std::string, std::string>>;

inline std::string fmtg(double v, int prec = 9)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*g", prec, v);
    return buf;
}

// Итог построения геометрии.
struct Built
{
    size_t ntotal = 0;
    size_t moni = 0;   // частица-монитор (cfg.moni_particle)
    KV kv;             // ключи для shot.txt (слои, число частиц, место монитора)
};

// Разбор кривой. Заполняет ключи для shot.txt и подпись графика; lab.ref_shift --
// сдвиг эталона при отрисовке (0 -- сырая шкала).
using AnalyzeFn = std::function<void(const e2e::Curve& calc, const e2e::Curve& ref, double delta_x,
                                     size_t ntotal, KV& kv, e2e::PlotLabels& lab)>;

struct CaseHooks
{
    std::string name;          // "bc2" -- префикс каталога прогона
    std::string title;         // для справки и shot.txt
    fs::path setup_dir;        // cases/<кейс>/experiment
    double dx_default = 2e-4;  // полный кейс
    std::function<Built(ParticleSystem&, double)> build;
    AnalyzeFn analyze;
};

struct Options
{
    fs::path out;
    fs::path ic;                       // пусто -- setup_dir/initial_conditions.txt
    std::vector<double> dx_list;       // пусто -- dx_default
    std::vector<int> contact_list;     // пусто -- из конфига
    double t_total = 0.0;              // 0 -- из конфига
    double save_t = 0.0;               // 0 -- из конфига
    int threads = 0;                   // 0 -- OMP_NUM_THREADS / все
    double time_limit = 0.0;           // 0 -- без лимита
    long long max_steps = 0;           // 0 -- без лимита
    size_t checkpoint_every = 0;       // 0 -- без чекпоинтов
    bool restart = false;
    bool force = false;                // --restart --force: физические отличия cfg допускаются
    bool overwrite = false;            // стереть прежний каталог прогона
    bool quiet = false;                // глушить вывод конвейера
    std::string tag;
};

inline const char* reason_name(RunGuard::Reason r)
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

inline void print_usage(const CaseHooks& h)
{
    std::printf(
        "%s_shot -- %s (экспериментальная постановка %s)\n"
        "Использование: %s_shot [ключи]\n"
        "  --out=<каталог>          корень результатов (по умолчанию runs/%s)\n"
        "  --contact-type=N         contact_force_type 0..3 (по умолчанию -- из конфига)\n"
        "  --contact-list=N,M,...   серия по типам контакта\n"
        "  --delta-x=М              шаг решётки, м (по умолчанию %g -- полный кейс)\n"
        "  --dx-list=М,М,...        серия по шагу решётки (сходимость)\n"
        "  --t-total=С              время расчёта (по умолчанию -- из конфига)\n"
        "  --save-t=С               интервал кадров (по умолчанию -- из конфига)\n"
        "  --ic=<файл>              другой файл параметров решателя\n"
        "  --threads=N              число потоков OpenMP\n"
        "  --time-limit=С           бюджет стенного времени на прогон, 0 -- без лимита (по умолчанию)\n"
        "  --max-steps=N            останов после N шагов (StepLimit)\n"
        "  --checkpoint-every=N     чекпоинт каждые N шагов (checkpoint.ckpt в каталоге прогона)\n"
        "  --restart [--force]      продолжить прогон из checkpoint.ckpt (те же ключи сцены)\n"
        "  --overwrite              стереть прежний каталог прогона и считать заново\n"
        "  --tag=<строка>           суффикс имени каталога прогона\n"
        "  --quiet                  не печатать вывод конвейера\n"
        "Каталог прогона: <out>/%s_dx<мкм>_c<тип>[_<tag>]/ -- log.txt, кадры, curve.txt, compare.svg, shot.txt.\n"
        "Коды выхода: 0 -- все прогоны дошли до t_total (или --max-steps), 1 -- оборван хотя бы один, 2 -- ошибка ключей.\n",
        h.name.c_str(), h.title.c_str(), h.setup_dir.string().c_str(), h.name.c_str(), h.name.c_str(),
        h.dx_default, h.name.c_str());
}

inline bool parse_double(const std::string& s, double& out)
{
    try
    {
        size_t pos = 0;
        out = std::stod(s, &pos);
        return pos == s.size() && std::isfinite(out);
    }
    catch (const std::exception&) { return false; }
}

inline bool parse_ll(const std::string& s, long long& out)
{
    try
    {
        size_t pos = 0;
        out = std::stoll(s, &pos);
        return pos == s.size();
    }
    catch (const std::exception&) { return false; }
}

template <class T, class F>
bool parse_list(const std::string& s, std::vector<T>& out, F parse_one)
{
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        T v{};
        if (!parse_one(tok, v)) return false;
        out.push_back(v);
    }
    return !out.empty();
}

// 0 -- разобрано, 1 -- справка, 2 -- ошибка.
inline int parse_args(int argc, char** argv, const CaseHooks& h, Options& o)
{
    o.out = fs::path("runs") / h.name;
    auto bad = [](const std::string& a) {
        std::fprintf(stderr, "ошибка: неверный ключ или значение: %s\n", a.c_str());
        return 2;
    };
    auto pdx = [](const std::string& t, double& v) { return parse_double(t, v) && v > 0.0; };
    auto pct = [](const std::string& t, int& v) {
        long long n = 0;
        if (!parse_ll(t, n) || n < 0 || n > 3) return false;
        v = static_cast<int>(n);
        return true;
    };
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto val = [&](const char* key) { return a.substr(std::char_traits<char>::length(key)); };
        auto starts = [&](const char* key) { return a.rfind(key, 0) == 0; };
        double d = 0.0;
        long long n = 0;
        if (a == "-h" || a == "--help") return 1;
        else if (starts("--out=")) o.out = val("--out=");
        else if (starts("--ic="))
        {
            o.ic = val("--ic=");
            if (!fs::exists(o.ic)) return bad(a + " (файл не найден)");
        }
        else if (starts("--contact-type="))
        {
            int c = 0;
            if (!pct(val("--contact-type="), c)) return bad(a);
            o.contact_list = {c};
        }
        else if (starts("--contact-list="))
        {
            o.contact_list.clear();
            if (!parse_list<int>(val("--contact-list="), o.contact_list, pct)) return bad(a);
        }
        else if (starts("--delta-x="))
        {
            if (!pdx(val("--delta-x="), d)) return bad(a);
            o.dx_list = {d};
        }
        else if (starts("--dx-list="))
        {
            o.dx_list.clear();
            if (!parse_list<double>(val("--dx-list="), o.dx_list, pdx)) return bad(a);
        }
        else if (starts("--t-total="))
        {
            if (!pdx(val("--t-total="), o.t_total)) return bad(a);
        }
        else if (starts("--save-t="))
        {
            if (!pdx(val("--save-t="), o.save_t)) return bad(a);
        }
        else if (starts("--threads="))
        {
            if (!parse_ll(val("--threads="), n) || n < 1 || n > 4096) return bad(a);
            o.threads = static_cast<int>(n);
        }
        else if (starts("--time-limit="))
        {
            if (!parse_double(val("--time-limit="), d) || d < 0.0) return bad(a);
            o.time_limit = d;
        }
        else if (starts("--max-steps="))
        {
            if (!parse_ll(val("--max-steps="), n) || n < 1) return bad(a);
            o.max_steps = n;
        }
        else if (starts("--checkpoint-every="))
        {
            if (!parse_ll(val("--checkpoint-every="), n) || n < 1) return bad(a);
            o.checkpoint_every = static_cast<size_t>(n);
        }
        else if (a == "--restart") o.restart = true;
        else if (a == "--force") o.force = true;
        else if (a == "--overwrite") o.overwrite = true;
        else if (a == "--quiet") o.quiet = true;
        else if (starts("--tag="))
        {
            o.tag = val("--tag=");
            if (o.tag.find_first_of("/\\ ") != std::string::npos) return bad(a);
        }
        else return bad(a);
    }
    if (o.restart && o.overwrite) return bad("--restart вместе с --overwrite");
    if (o.force && !o.restart) return bad("--force без --restart");
    if (o.dx_list.empty()) o.dx_list = {h.dx_default};
    if (o.ic.empty()) o.ic = h.setup_dir / "initial_conditions.txt";
    return 0;
}

struct CwdGuard
{
    fs::path old;
    explicit CwdGuard(const fs::path& to) : old(fs::current_path()) { fs::current_path(to); }
    ~CwdGuard()
    {
        std::error_code ec;
        fs::current_path(old, ec);
    }
    CwdGuard(const CwdGuard&) = delete;
    CwdGuard& operator=(const CwdGuard&) = delete;
};

// Кривая из log.txt без повторов: после рестарта конвейер дописывает в log.txt шаги с момента
// чекпоинта, и шаги между чекпоинтом и обрывом оказываются в файле дважды. Остаётся последняя
// запись: точки, время которых не меньше времени очередной, отбрасываются.
inline e2e::Curve monotone_curve(const e2e::Curve& c)
{
    e2e::Curve out;
    for (size_t i = 0; i < c.t.size(); ++i)
    {
        while (!out.t.empty() && out.t.back() >= c.t[i])
        {
            out.t.pop_back();
            out.v.pop_back();
            out.itype.pop_back();
        }
        out.t.push_back(c.t[i]);
        out.v.push_back(c.v[i]);
        out.itype.push_back(c.itype[i]);
    }
    return out;
}

inline std::string run_dir_name(const CaseHooks& h, double dx, int contact, const std::string& tag)
{
    const long long um = std::llround(dx * 1e6);
    std::string s = h.name + "_dx" + std::to_string(um) + "_c" + std::to_string(contact);
    if (!tag.empty()) s += "_" + tag;
    return s;
}

struct RunSummary
{
    std::string dir;
    double dx = 0.0;
    int contact = 0;
    size_t ntotal = 0;
    std::string outcome;   // completed | stopped (--max-steps) | guard | error
    double t_reached = 0.0;
    size_t steps = 0;
    double wall = 0.0;
};

inline void append_runs_csv(const fs::path& out, const CaseHooks& h, const RunSummary& r)
{
    const fs::path p = out / "runs.csv";
    const bool fresh = !fs::exists(p);
    std::ofstream f(p, std::ios::app);
    if (fresh) f << "case,dir,delta_x,contact,ntotal,outcome,t_reached,steps,wall_s\n";
    f << h.name << ',' << r.dir << ',' << fmtg(r.dx) << ',' << r.contact << ',' << r.ntotal << ','
      << r.outcome << ',' << fmtg(r.t_reached) << ',' << r.steps << ',' << fmtg(r.wall, 6) << '\n';
}

inline void write_kv(const fs::path& p, const KV& kv)
{
    std::ofstream f(p, std::ios::trunc);
    for (const auto& [k, v] : kv) f << k << '=' << v << '\n';
}

// Один прогон. Возвращает true, если дошёл до t_total.
inline bool run_one(const CaseHooks& h, const Options& o, double dx, int contact_opt, const e2e::Curve& ref,
                    RunSummary& sum)
{
    // Конфиг читается до каталога: тип контакта по умолчанию берётся из него.
    SimulationConfig cfg;
    cfg.apply_defaults();
    if (!fs::exists(o.ic)) throw std::runtime_error("не найден файл параметров решателя: " + o.ic.string());
    // parse_sph_config молча оставляет умолчания, если файл не открылся, -- путь
    // берётся абсолютным, потому что ниже меняется CWD.
    const fs::path ic_abs = fs::absolute(o.ic);
    parse_sph_config(cfg, ic_abs.string());
    const int contact = contact_opt >= 0 ? contact_opt : cfg.contact_force_type;
    cfg.contact_force_type = contact;
    if (o.t_total > 0.0) cfg.t_total = o.t_total;
    if (o.save_t > 0.0) cfg.save_t = o.save_t;
    if (o.threads > 0)
    {
        cfg.num_threads = o.threads;
#ifdef _OPENMP
        omp_set_num_threads(o.threads);
#endif
    }
    else
    {
#ifdef _OPENMP
        cfg.num_threads = omp_get_max_threads();
#endif
    }

    const std::string name = run_dir_name(h, dx, contact, o.tag);
    const fs::path sub = fs::absolute(o.out / name);
    sum.dir = name;
    sum.dx = dx;
    sum.contact = contact;

    if (!o.restart && fs::exists(sub / "log.txt"))
    {
        if (!o.overwrite)
            throw std::runtime_error("каталог прогона уже содержит расчёт: " + sub.string() +
                                     " (--restart -- продолжить, --overwrite -- считать заново)");
        fs::remove_all(sub);
    }
    if (o.restart && !fs::exists(sub))
        throw std::runtime_error("--restart: нет каталога прогона " + sub.string());
    fs::create_directories(sub);
    CwdGuard cwd(sub);

    std::printf("\n=== %s: delta_x = %g м, contact_force_type = %d, t_total = %g с -> %s\n", h.name.c_str(), dx,
                contact, cfg.t_total, sub.string().c_str());
    std::fflush(stdout);

    fs::copy_file(ic_abs, sub / "initial_conditions.used.txt", fs::copy_options::overwrite_existing);

    fiber::Model model;
    ParticleSystem& sys = model.solid;
    const Built b = h.build(sys, dx);
    const size_t ntotal = b.ntotal;
    sum.ntotal = ntotal;
    cfg.moni_particle = b.moni;

    // Запас под зеркальные частицы -- как в main.cpp и в сквозных тестах.
    cfg.maxn = static_cast<size_t>(static_cast<double>(ntotal) * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    if (sys.size() < cfg.maxn)
    {
        const size_t prev = sys.size();
        sys.resize(cfg.maxn);
        for (size_t i = prev; i < sys.size(); ++i) sys.inzializ_index(i);
    }
    for (size_t i = 0; i < ntotal; ++i) urs(cfg.mats, sys, i);
    double dt = compute_cfl_timestep(cfg, sys, ntotal);
    double time = 0.0;

    // Сцена для сверки при рестарте: геометрия в cfg чекпоинта не входит.
    std::ostringstream sid;
    sid << "case=" << h.name << " delta_x=" << fmtg(dx, 17) << " ntotal=" << ntotal << " moni=" << b.moni
        << " contact=" << contact << "\n";
    const fs::path ckpt = sub / "checkpoint.ckpt";
    if (o.restart)
    {
        const fs::path id_path = sub / "scene_id.txt";
        if (fs::exists(id_path))
        {
            std::ifstream idf(id_path);
            const std::string saved((std::istreambuf_iterator<char>(idf)), std::istreambuf_iterator<char>());
            if (saved != sid.str())
                throw std::runtime_error("--restart: сцена не совпадает с сохранённой:\n  сохранено: " + saved +
                                         "  сейчас:    " + sid.str());
        }
        const fs::path prev = ckpt.string() + ".prev";
        if (!fs::exists(ckpt) && !fs::exists(prev))
            throw std::runtime_error("--restart: нет чекпоинта " + ckpt.string() + " (и .prev)");
        cfg.restart_file = (fs::exists(ckpt) ? ckpt : prev).string();
        cfg.restart_force = o.force;
    }
    else
    {
        // log.txt конвейер дописывает; новый прогон начинается с пустого файла (с BOM, как в main.cpp).
        std::ofstream lf("log.txt", std::ios::trunc);
        lf << "\xEF\xBB\xBF";
    }
    if (o.checkpoint_every > 0)
    {
        cfg.checkpoint_file = ckpt.string();
        cfg.checkpoint_every_steps = o.checkpoint_every;
        if (!o.restart)
        {
            std::ofstream idf(sub / "scene_id.txt", std::ios::trunc);
            idf << sid.str();
        }
    }

    RunGuard guard(o.time_limit, 2.0e6, 200);
    if (o.max_steps > 0) guard.set_step_limit(static_cast<size_t>(o.max_steps));

    std::string error;
    const auto t0 = std::chrono::steady_clock::now();
    try
    {
        e2e::ScopedSilence silence(o.quiet);
        time_integration(cfg, model, ntotal, dt, time, dx, &guard);
    }
    catch (const std::exception& e)
    {
        error = e.what();
    }
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::string err;
    const e2e::Curve curve = monotone_curve(e2e::read_curve(sub / "log.txt", err));
    if (error.empty() && !err.empty()) error = err;

    std::string outcome = "completed";
    if (!error.empty()) outcome = "error";
    else if (guard.failed()) outcome = "guard";
    else if (guard.triggered()) outcome = "stopped";   // --max-steps: штатный останов
    const double t_reached = curve.empty() ? time : curve.t.back();

    e2e::PlotSeries s;
    s.curve = curve;
    s.type = contact;
    s.aborted = outcome != "completed";
    s.steps = curve.size();
    std::string werr;
    e2e::write_curve_dump(s, sub / "curve.txt", h.name.c_str(), werr);

    KV kv;
    kv.emplace_back("case", h.name);
    kv.emplace_back("title", h.title);
    kv.emplace_back("setup", h.setup_dir.string());
    kv.emplace_back("ic", ic_abs.string());
    kv.emplace_back("delta_x", fmtg(dx));
    kv.emplace_back("contact_force_type", std::to_string(contact));
    kv.emplace_back("t_total", fmtg(cfg.t_total));
    kv.emplace_back("save_t", fmtg(cfg.save_t));
    kv.emplace_back("threads", std::to_string(cfg.num_threads));
    kv.emplace_back("ntotal", std::to_string(ntotal));
    kv.emplace_back("moni_particle", std::to_string(b.moni));
    for (const auto& p : b.kv) kv.push_back(p);
    kv.emplace_back("outcome", outcome);
    kv.emplace_back("guard_reason", reason_name(guard.reason()));
    if (guard.triggered()) kv.emplace_back("guard_message", guard.message());
    if (!error.empty()) kv.emplace_back("error", error);
    kv.emplace_back("restart", o.restart ? "1" : "0");
    kv.emplace_back("t_reached", fmtg(t_reached));
    kv.emplace_back("steps", std::to_string(curve.size()));
    kv.emplace_back("dt_last", fmtg(dt));
    kv.emplace_back("wall_s", fmtg(wall, 6));

    e2e::PlotLabels lab;
    if (h.analyze && !curve.empty()) h.analyze(curve, ref, dx, ntotal, kv, lab);
    if (!curve.empty())
    {
        if (lab.title.empty()) lab.title = h.title;
        lab.ref_legend = "эталон (прибор)";
        e2e::emit_plot({s}, ref, lab, sub / "compare.svg");
    }
    write_kv(sub / "shot.txt", kv);

    sum.outcome = outcome;
    sum.t_reached = t_reached;
    sum.steps = curve.size();
    sum.wall = wall;

    std::printf("=== %s: %s, t = %g с, шагов %zu, %.1f с%s%s\n", name.c_str(), outcome.c_str(), t_reached,
                curve.size(), wall, error.empty() ? "" : ": ", error.c_str());
    if (guard.triggered()) std::printf("    сторож: %s\n", guard.message().c_str());
    std::fflush(stdout);
    return outcome == "completed" || outcome == "stopped";
}

inline int run_main(int argc, char** argv, const CaseHooks& h)
{
    Options o;
    const int pa = parse_args(argc, argv, h, o);
    if (pa == 1)
    {
        print_usage(h);
        return 0;
    }
    if (pa != 0)
    {
        print_usage(h);
        return 2;
    }

    std::string err;
    const fs::path ref_path = h.setup_dir / "reference.txt";
    const e2e::Curve ref = e2e::read_curve(ref_path, err);
    if (!err.empty()) std::fprintf(stderr, "предупреждение: эталон не прочитан: %s\n", err.c_str());

    std::error_code ec;
    fs::create_directories(o.out, ec);

    std::vector<int> contacts = o.contact_list;
    if (contacts.empty()) contacts = {-1};   // из конфига
    bool all_ok = true;
    for (double dx : o.dx_list)
    {
        for (int c : contacts)
        {
            RunSummary sum;
            try
            {
                all_ok = run_one(h, o, dx, c, ref, sum) && all_ok;
            }
            catch (const std::exception& e)
            {
                std::fprintf(stderr, "ошибка прогона (delta_x = %g, contact = %d): %s\n", dx, c, e.what());
                sum.outcome = "error";
                all_ok = false;
            }
            if (!sum.dir.empty()) append_runs_csv(o.out, h, sum);
        }
    }
    return all_ok ? 0 : 1;
}

} // namespace ceramic
