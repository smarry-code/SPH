// guo2019_shot.cpp
// Серийные выстрелы Guo 2019 (ТЗ Guo 2019 §8): один или несколько выстрелов стального
// шара по мишени A / 1T / 1T_A / A_1T, сводка shots.csv, ballistic limit V_BL (vbl.txt)
// и замер стоимости (ТЗ Guo 2019 §7). Не тест: add_test() не добавлен.
//
// Использование:
//   guo2019_shot --config=A|1T|1T_A|A_1T|<k>T --level=base|fine|smoke|scaled
//                (--v0=<м/с> | --paper-shots | --v0-list=a,b,c)
//                --out=<каталог> [--frames=off|camera] [--t-max=<с>] [--time-limit=<с>]
//                [--max-steps=N] [--no-body-isolation] [--mu=<μ_ff=μ_fs>] [--mu-ff=<x>] [--mu-fs=<x>] [--repeat=N]
//                [--bracket|--no-bracket] [--n-mon=N] [--n-diag=N] [--threads=N]
//                [--probe-frame] [--plate-nz=N] [--save-t=<с>] [--ball-d=<м>] [--ball-mass=<кг>] [--ball-dx=<м>]
//                [--fiber-inelastic=<файл>] [--eps-break=<x>] [--energy-log=N]
//                [--checkpoint] [--checkpoint-every=<N шагов>] [--restart [--force]]
//   <k>T -- пакет из k слоёв ткани без пластины, k = 1..64 (1T, 2T, 4T -- как раньше; прочие k без данных статьи,
//   t_max по умолчанию 600 мкс, оценка). --ball-d / --ball-mass / --ball-dx заменяют диаметр, массу и шаг решётки шара
//   уровня (плотность выводится из массы); к имени каталога выстрела добавляется _bD<мм>m<г>dx<мм>, в shot.txt --
//   ball_D=, ball_m=, ball_dx=, ball_rho=. При --checkpoint рядом пишется scene_id.txt; --restart с другим числом слоёв
//   или другим шаром отклоняется (старые чекпоинты без scene_id.txt не проверяются).
//   --mu задаёт ОБА коэффициента (μ_ff нить-нить и μ_fs нить-Solid); --mu-ff / --mu-fs переопределяют только свой
//   коэффициент и имеют приоритет над --mu независимо от порядка ключей. Не заданный из пары берётся из --mu
//   (по умолчанию kMuDefault). Без новых ключей поведение прежнее. Значения пишутся в shot.txt (mu_ff=, mu_fs=).
//   --energy-log=N -- energy_log.csv каждые N шагов в подкаталоге выстрела (только сборка с -DSPH_ENERGY_LOG=ON;
//   include/energy_log.hpp).
//   --checkpoint / --checkpoint-every=N -- чекпоинт <каталог выстрела>/checkpoint.ckpt (+ .prev) каждые N шагов
//   интегрирования (N >= 1, по умолчанию 10000) и финальный при останове; --restart -- продолжить выстрел из него
//   (нет файла -- ошибка); --force -- разрешить менять физические настройки при --restart (только с --restart).
//   Состояние монитора хранится в чекпоинте (блок USER, формат v3): после --restart вердикт, максимумы (n_tunnel,
//   max_delta_*) и monitor.csv -- как при непрерывном счёте; из чекпоинта v1/v2 -- только с шага рестарта.
//   --fiber-inelastic -- файл с директивами PLASTIC/PLASTIC_TABLE/RATE/VISCO/FAILURE для материала 0 (без MATERIAL/FIBER);
//   разбирается штатным парсером .fib и переносится в материалы сцены; --eps-break переопределяет eps_break
//   (и l_break всех связей = (1+eps_break) l0).
//
// Каждый выстрел -- свежая сцена (build_guo_scene) в своём подкаталоге
// <out>/<config>_<level>[_nz<N>][_bD<мм>m<г>dx<мм>]_v<V0>_r<i>/; файлы конвейера (log.txt, кадры) пишутся в CWD,
// поэтому на время выстрела делается chdir туда. Код сцены и монитора -- tests/guo2019_*.hpp.

#include "guo2019_monitor.hpp"
#include "guo2019_scene.hpp"

#include "checkpoint.hpp"
#include "e2e_harness.hpp"

#include "fiber/fiber_input.hpp"
#include "fiber/fiber_output.hpp"
#include "fiber/fiber_plasticity.hpp"
#include "output.hpp"
#include "run_guard.hpp"
#include "time_integration.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;
using namespace guo2019;

namespace
{

using KV = std::vector<std::pair<std::string, std::string>>;

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

std::string fmtg(double v, int prec = 9)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", prec, v);
    return buf;
}

std::string fmt_v0(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.1f", v);
    return buf;
}

// ---------------------------------------------------------------------------
// Параметры командной строки
// ---------------------------------------------------------------------------
struct Options
{
    Config config = Config::A;
    int plies = 0; // число слоёв для Config::TK (для 1T, 2T, 4T тоже заполняется)
    Level level = Level::Base;
    bool have_config = false, have_level = false;
    bool have_v0 = false, paper_shots = false, have_list = false;
    std::vector<double> v0_list; // --v0 или --v0-list
    fs::path out;
    bool frames_camera = false;
    double save_t = -1.0; // > 0 -- интервал кадров (--save-t, для анимации); иначе по --frames
    double t_max = -1.0; // < 0 -- default_t_max(config)
    double time_limit = 0.0;
    std::size_t max_steps = 0;
    bool body_isolation = true;
    double mu = kMuDefault;
    double mu_ff = -1.0; // < 0 -- не задан, берётся mu
    double mu_fs = -1.0; // < 0 -- не задан, берётся mu
    int repeat = 1;
    bool repeat_set = false;
    int bracket = -1; // -1 -- по умолчанию (включён с --paper-shots), 0 -- выкл., 1 -- вкл.
    std::size_t n_mon = 1, n_diag = 10;
    int threads = 0; // 0 -- как в cfg по умолчанию
    int probe_frame = -1; // -1 -- по умолчанию (включён при --max-steps)
    int plate_nz = 0;     // > 0 -- переопределить n_z пластины уровня (сравнение стоимости)
    double ball_d = 0.0, ball_m = 0.0, ball_dx = 0.0; // > 0 -- переопределить D_s, m_s, dx_s шара уровня
    fs::path fiber_inelastic; // непусто -- файл директив неупругости нити (материал 0)
    double eps_break = -1.0;  // > 0 -- переопределить eps_break нити
    double fiber_E = -1.0;    // > 0 -- переопределить модуль Юнга нити, Па (все материалы ткани)
    // Тип 3: выталкивание при глубоком проникании (оба контакта нитей) и gamma; множитель шага.
    double impulse_push_tau = 0.0, impulse_push_depth = 0.1, impulse_gamma = -1.0, impulse_push_vmax = 0.0;
    double dt_scale = 1.0;
    std::size_t energy_log = 0; // > 0 -- energy_log.csv каждые N шагов (только сборка с SPH_ENERGY_LOG)
    // Законы нормальной силы контактов нитей (fiber_contact_law.hpp): -1 -- как в
    // guo2019_initial_conditions.txt (по умолчанию 0 = Герц), иначе 0..3.
    int fs_contact = -1;
    int ff_contact = -1;
    // Нормальная сила законов 1-3 как проекция скорости: -1 -- не задано (берётся из конфига).
    int ff_predictor = -1, fs_predictor = -1;
    int proj_iters = -1, proj_warm = -1;
    std::size_t contact_stats = 0; // > 0 -- contact_stats.csv каждые N шагов (только сборка с SPH_CONTACT_STATS)
    bool hold_n = false;           // контрфакт «удержание N» (только сборка с SPH_CONTACT_STATS)
    bool checkpoint = false;       // писать чекпоинты (--checkpoint или --checkpoint-every)
    std::size_t checkpoint_every = 10000; // период чекпоинта, шаги интегрирования
    bool restart = false;          // продолжить выстрел из checkpoint.ckpt
    bool force = false;            // --force: разрешить менять физические настройки при рестарте
};

bool parse_double(const std::string& s, double& out)
{
    try
    {
        std::size_t pos = 0;
        out = std::stod(s, &pos);
        return pos == s.size() && std::isfinite(out);
    }
    catch (const std::exception&)
    {
        return false;
    }
}

bool parse_ll(const std::string& s, long long& out)
{
    try
    {
        std::size_t pos = 0;
        out = std::stoll(s, &pos);
        return pos == s.size();
    }
    catch (const std::exception&)
    {
        return false;
    }
}

void print_usage()
{
    std::printf("Использование: guo2019_shot --config=A|1T|1T_A|A_1T|<k>T --level=base|fine|smoke|scaled\n"
                "                   (--v0=<м/с> | --paper-shots | --v0-list=a,b,c)\n"
                "                   --out=<каталог> [--frames=off|camera] [--t-max=<с>] [--time-limit=<с>]\n"
                "                   [--max-steps=N] [--no-body-isolation] [--mu=<μ>] [--mu-ff=<μ>] [--mu-fs=<μ>]\n"
                "                   [--repeat=N]\n"
                "                   [--bracket|--no-bracket] [--n-mon=N] [--n-diag=N] [--threads=N]\n"
                "                   [--probe-frame|--no-probe-frame] [--plate-nz=N] [--save-t=<с>]\n"
                "                   [--ball-d=<м>] [--ball-mass=<кг>] [--ball-dx=<м>]  (шар: диаметр, масса, шаг решётки; каждый > 0)\n"
                "                   [--fs-contact=0..4] [--ff-contact=0..4]  (закон контакта нитей, 0 = Герц, 4 = Герц + демпфер)\n"
                "                   [--ff-contact-predictor=0|1] [--fs-contact-predictor=0|1] [--contact-proj-iters=1..100] [--contact-proj-warm=0|1]\n"
                "                   (нормальная сила законов 1-3 как проекция скорости; по умолчанию 0)\n"
                "                   [--contact-stats=N] [--hold-n]  (диагностика контактов, сборка с -DSPH_CONTACT_STATS=ON)\n"
                "                   [--fiber-inelastic=<файл.fib>] [--eps-break=<x>] [--fiber-E=<Па>]\n"
                "                   [--energy-log=N]  (лог энергии, сборка с -DSPH_ENERGY_LOG=ON)\n"
                "                   [--impulse-push-tau=<с>] [--impulse-push-depth=<доля r>] [--impulse-push-vmax=<м/с>] [--impulse-gamma=<g>]  (тип 3)\n"
                "                   [--dt-scale=<s из (0;1]>]  (множитель шага, сходимость)\n"
                "                   [--checkpoint] [--checkpoint-every=<N>] [--restart [--force]]\n"
                "  <k>T: k слоёв ткани без пластины, k = 1..64 (1T, 2T, 4T -- как раньше; прочие k без данных статьи, t_max 600 мкс).\n"
                "  --ball-*: к имени каталога выстрела добавляется _bD<мм>m<г>dx<мм>; в shot.txt -- ball_D, ball_m, ball_dx, ball_rho;\n"
                "  --restart с другим числом слоёв или шаром отклоняется (по scene_id.txt, пишется при --checkpoint).\n"
                "  --checkpoint: чекпоинт <каталог выстрела>/checkpoint.ckpt (+ .prev); --checkpoint-every=N: период в шагах\n"
                "  интегрирования (N >= 1, по умолчанию 10000; ключ сам включает чекпоинты); финальный чекпоинт пишется при останове.\n"
                "  --restart: продолжить выстрел из checkpoint.ckpt (если файла нет -- ошибка; при битом основном берётся .prev).\n"
                "  Состояние монитора хранится в чекпоинте (v3): вердикт, максимумы и monitor.csv после --restart полные\n"
                "  (из чекпоинта v1/v2 -- только с шага рестарта).\n"
                "  --force: при --restart разрешить менять физические настройки (без него любое отличие -- отказ со списком полей).\n"
                "  --force без --restart -- ошибка. Риски --force:\n%s"
                "  --mu задаёт оба коэффициента трения (нить-нить и нить-Solid); --mu-ff / --mu-fs переопределяют только\n"
                "  свой и имеют приоритет над --mu независимо от порядка ключей; не заданный из пары берётся из --mu\n"
                "  (по умолчанию %g).\n", checkpoint::force_risks_text().c_str(), kMuDefault);
}

// 0 -- разобрано, 1 -- справка, 2 -- ошибка
int parse_args(int argc, char** argv, Options& o)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto val = [&](const char* key) -> std::string { return a.substr(std::strlen(key)); };
        auto bad = [&](const char* what) {
            std::fprintf(stderr, "Неверное значение %s: %s\n", what, a.c_str());
            return 2;
        };
        double d = 0.0;
        long long n = 0;
        if (a == "-h" || a == "--help") return 1;
        else if (a.rfind("--config=", 0) == 0)
        {
            std::string cerr_msg;
            if (!parse_config(val("--config="), o.config, o.plies, &cerr_msg))
            {
                if (!cerr_msg.empty()) std::fprintf(stderr, "%s\n", cerr_msg.c_str());
                return bad("--config");
            }
            o.have_config = true;
        }
        else if (a.rfind("--level=", 0) == 0)
        {
            if (!parse_level(val("--level="), o.level)) return bad("--level");
            o.have_level = true;
        }
        else if (a.rfind("--v0-list=", 0) == 0)
        {
            std::stringstream ss(val("--v0-list="));
            std::string tok;
            while (std::getline(ss, tok, ','))
            {
                if (!parse_double(tok, d) || !(d > 0.0)) return bad("--v0-list");
                o.v0_list.push_back(d);
            }
            if (o.v0_list.empty()) return bad("--v0-list");
            o.have_list = true;
        }
        else if (a.rfind("--v0=", 0) == 0)
        {
            if (!parse_double(val("--v0="), d) || !(d > 0.0)) return bad("--v0");
            o.v0_list.push_back(d);
            o.have_v0 = true;
        }
        else if (a == "--paper-shots") o.paper_shots = true;
        else if (a.rfind("--out=", 0) == 0) o.out = val("--out=");
        else if (a.rfind("--frames=", 0) == 0)
        {
            const std::string f = val("--frames=");
            if (f == "off") o.frames_camera = false;
            else if (f == "camera") o.frames_camera = true;
            else return bad("--frames");
        }
        else if (a.rfind("--save-t=", 0) == 0)
        {
            if (!parse_double(val("--save-t="), d) || !(d > 0.0)) return bad("--save-t");
            o.save_t = d;
        }
        else if (a.rfind("--t-max=", 0) == 0)
        {
            if (!parse_double(val("--t-max="), d) || !(d > 0.0)) return bad("--t-max");
            o.t_max = d;
        }
        else if (a.rfind("--time-limit=", 0) == 0)
        {
            if (!parse_double(val("--time-limit="), d) || d < 0.0) return bad("--time-limit");
            o.time_limit = d;
        }
        else if (a.rfind("--max-steps=", 0) == 0)
        {
            if (!parse_ll(val("--max-steps="), n) || n < 1) return bad("--max-steps");
            o.max_steps = static_cast<std::size_t>(n);
        }
        else if (a == "--checkpoint") o.checkpoint = true;
        else if (a.rfind("--checkpoint-every=", 0) == 0)
        {
            if (!parse_ll(val("--checkpoint-every="), n) || n < 1) return bad("--checkpoint-every");
            o.checkpoint_every = static_cast<std::size_t>(n);
            o.checkpoint = true;
        }
        else if (a == "--restart") o.restart = true;
        else if (a == "--force") o.force = true;
        else if (a == "--no-body-isolation") o.body_isolation = false;
        else if (a.rfind("--mu=", 0) == 0)
        {
            if (!parse_double(val("--mu="), d) || d < 0.0) return bad("--mu");
            o.mu = d;
        }
        else if (a.rfind("--mu-ff=", 0) == 0)
        {
            if (!parse_double(val("--mu-ff="), d) || d < 0.0) return bad("--mu-ff");
            o.mu_ff = d;
        }
        else if (a.rfind("--mu-fs=", 0) == 0)
        {
            if (!parse_double(val("--mu-fs="), d) || d < 0.0) return bad("--mu-fs");
            o.mu_fs = d;
        }
        else if (a.rfind("--repeat=", 0) == 0)
        {
            if (!parse_ll(val("--repeat="), n) || n < 1 || n > 100) return bad("--repeat");
            o.repeat = static_cast<int>(n);
            o.repeat_set = true;
        }
        else if (a == "--bracket") o.bracket = 1;
        else if (a == "--no-bracket") o.bracket = 0;
        else if (a.rfind("--n-mon=", 0) == 0)
        {
            if (!parse_ll(val("--n-mon="), n) || n < 1) return bad("--n-mon");
            o.n_mon = static_cast<std::size_t>(n);
        }
        else if (a.rfind("--n-diag=", 0) == 0)
        {
            if (!parse_ll(val("--n-diag="), n) || n < 1) return bad("--n-diag");
            o.n_diag = static_cast<std::size_t>(n);
        }
        else if (a.rfind("--threads=", 0) == 0)
        {
            if (!parse_ll(val("--threads="), n) || n < 1 || n > 4096) return bad("--threads");
            o.threads = static_cast<int>(n);
        }
        else if (a.rfind("--plate-nz=", 0) == 0)
        {
            if (!parse_ll(val("--plate-nz="), n) || n < 1 || n > 64) return bad("--plate-nz");
            o.plate_nz = static_cast<int>(n);
        }
        else if (a.rfind("--ball-d=", 0) == 0)
        {
            if (!parse_double(val("--ball-d="), d) || !(d > 0.0)) return bad("--ball-d");
            o.ball_d = d;
        }
        else if (a.rfind("--ball-mass=", 0) == 0)
        {
            if (!parse_double(val("--ball-mass="), d) || !(d > 0.0)) return bad("--ball-mass");
            o.ball_m = d;
        }
        else if (a.rfind("--ball-dx=", 0) == 0)
        {
            if (!parse_double(val("--ball-dx="), d) || !(d > 0.0)) return bad("--ball-dx");
            o.ball_dx = d;
        }
        else if (a.rfind("--fs-contact=", 0) == 0)
        {
            if (!parse_ll(val("--fs-contact="), n) || n < 0 || n > 4) return bad("--fs-contact");
            o.fs_contact = static_cast<int>(n);
        }
        else if (a.rfind("--ff-contact=", 0) == 0)
        {
            if (!parse_ll(val("--ff-contact="), n) || n < 0 || n > 4) return bad("--ff-contact");
            o.ff_contact = static_cast<int>(n);
        }
        else if (a.rfind("--ff-contact-predictor=", 0) == 0)
        {
            if (!parse_ll(val("--ff-contact-predictor="), n) || n < 0 || n > 1) return bad("--ff-contact-predictor");
            o.ff_predictor = static_cast<int>(n);
        }
        else if (a.rfind("--fs-contact-predictor=", 0) == 0)
        {
            if (!parse_ll(val("--fs-contact-predictor="), n) || n < 0 || n > 1) return bad("--fs-contact-predictor");
            o.fs_predictor = static_cast<int>(n);
        }
        else if (a.rfind("--contact-proj-iters=", 0) == 0)
        {
            if (!parse_ll(val("--contact-proj-iters="), n) || n < 1 || n > 100) return bad("--contact-proj-iters");
            o.proj_iters = static_cast<int>(n);
        }
        else if (a.rfind("--contact-proj-warm=", 0) == 0)
        {
            if (!parse_ll(val("--contact-proj-warm="), n) || n < 0 || n > 1) return bad("--contact-proj-warm");
            o.proj_warm = static_cast<int>(n);
        }
        else if (a.rfind("--contact-stats=", 0) == 0 || a == "--hold-n")
        {
#ifndef SPH_CONTACT_STATS
            std::fprintf(stderr, "%s: сборка без SPH_CONTACT_STATS (cmake -DSPH_CONTACT_STATS=ON)\n", a.c_str());
            return 2;
#else
            if (a == "--hold-n")
                o.hold_n = true;
            else
            {
                if (!parse_ll(val("--contact-stats="), n) || n < 1) return bad("--contact-stats");
                o.contact_stats = static_cast<std::size_t>(n);
            }
#endif
        }
        else if (a.rfind("--fiber-inelastic=", 0) == 0)
        {
            o.fiber_inelastic = val("--fiber-inelastic=");
            if (o.fiber_inelastic.empty() || !fs::exists(o.fiber_inelastic)) return bad("--fiber-inelastic (файл не найден)");
            o.fiber_inelastic = fs::absolute(o.fiber_inelastic);
        }
        else if (a.rfind("--fiber-E=", 0) == 0)
        {
            if (!parse_double(val("--fiber-E="), d) || !(d > 0.0)) return bad("--fiber-E");
            o.fiber_E = d;
        }
        else if (a.rfind("--impulse-push-tau=", 0) == 0)
        {
            if (!parse_double(val("--impulse-push-tau="), d) || d < 0.0) return bad("--impulse-push-tau");
            o.impulse_push_tau = d;
        }
        else if (a.rfind("--impulse-push-depth=", 0) == 0)
        {
            if (!parse_double(val("--impulse-push-depth="), d) || d < 0.0 || d >= 1.0) return bad("--impulse-push-depth");
            o.impulse_push_depth = d;
        }
        else if (a.rfind("--impulse-push-vmax=", 0) == 0)
        {
            if (!parse_double(val("--impulse-push-vmax="), d) || d < 0.0) return bad("--impulse-push-vmax");
            o.impulse_push_vmax = d;
        }
        else if (a.rfind("--impulse-gamma=", 0) == 0)
        {
            if (!parse_double(val("--impulse-gamma="), d) || !(d > 0.0) || d > 1.6) return bad("--impulse-gamma");
            o.impulse_gamma = d;
        }
        else if (a.rfind("--dt-scale=", 0) == 0)
        {
            if (!parse_double(val("--dt-scale="), d) || !(d > 0.0) || d > 1.0) return bad("--dt-scale");
            o.dt_scale = d;
        }
        else if (a.rfind("--eps-break=", 0) == 0)
        {
            if (!parse_double(val("--eps-break="), d) || !(d > 0.0)) return bad("--eps-break");
            o.eps_break = d;
        }
        else if (a.rfind("--energy-log=", 0) == 0)
        {
            if (!parse_ll(val("--energy-log="), n) || n < 1) return bad("--energy-log");
#ifndef SPH_ENERGY_LOG
            std::fprintf(stderr, "--energy-log: сборка без SPH_ENERGY_LOG (cmake -DSPH_ENERGY_LOG=ON)\n");
            return 2;
#else
            o.energy_log = static_cast<std::size_t>(n);
#endif
        }
        else if (a == "--probe-frame") o.probe_frame = 1;
        else if (a == "--no-probe-frame") o.probe_frame = 0;
        else
        {
            std::fprintf(stderr, "Неизвестный аргумент: %s\n", a.c_str());
            return 2;
        }
    }
    if (!o.have_config) { std::fprintf(stderr, "Нужен --config=A|1T|1T_A|A_1T|<k>T\n"); return 2; }
    if (!o.have_level) { std::fprintf(stderr, "Нужен --level=base|fine|smoke|scaled\n"); return 2; }
    if (o.out.empty()) { std::fprintf(stderr, "Нужен --out=<каталог>\n"); return 2; }
    const int modes = (o.have_v0 ? 1 : 0) + (o.paper_shots ? 1 : 0) + (o.have_list ? 1 : 0);
    if (modes != 1)
    {
        std::fprintf(stderr, "Нужно ровно одно из --v0=, --paper-shots, --v0-list=\n");
        return 2;
    }
    if (o.force && !o.restart) { std::fprintf(stderr, "--force допустим только вместе с --restart\n"); return 2; }
    if (o.bracket < 0) o.bracket = o.paper_shots ? 1 : 0;
    if (!o.repeat_set) o.repeat = o.paper_shots ? 2 : 1; // ТЗ Guo 2019 §8.3: --repeat=2 для выстрелов статьи
    if (o.probe_frame < 0) o.probe_frame = o.max_steps > 0 ? 1 : 0;
    return 0;
}

// ---------------------------------------------------------------------------
// Запись одного выстрела
// ---------------------------------------------------------------------------
struct ShotRecord
{
    double V0 = 0.0;
    int rep = 0;
    bool ok = false; // выстрел отработал и записал shot.txt
    std::string error;
    Outcome outcome = Outcome::Undetermined;
    double V_R = 0.0, e = 0.0, m_s = 0.0;
    fs::path dir;
    KV kv; // содержимое shot.txt (в порядке записи)
};

struct CwdGuard
{
    fs::path old;
    explicit CwdGuard(const fs::path& to) : old(fs::current_path()) { fs::current_path(to); }
    ~CwdGuard()
    {
        std::error_code ec;
        fs::current_path(old, ec);
    }
};

struct FrameSizes
{
    std::size_t n_sph = 0, n_fiber = 0;
    std::uintmax_t bytes_sph = 0, bytes_fiber = 0;
};

FrameSizes scan_frame_sizes(const fs::path& dir)
{
    FrameSizes fz;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(dir, ec))
    {
        const std::string n = e.path().filename().string();
        std::error_code ec2;
        if (n.rfind("sph_output_", 0) == 0)
        {
            ++fz.n_sph;
            fz.bytes_sph += fs::file_size(e.path(), ec2);
        }
        else if (n.rfind("fiber_output_", 0) == 0)
        {
            ++fz.n_fiber;
            fz.bytes_fiber += fs::file_size(e.path(), ec2);
        }
    }
    return fz;
}

KV parse_kv_line(const std::string& line)
{
    KV out;
    std::stringstream ss(line);
    std::string tok;
    while (ss >> tok)
    {
        const std::size_t eq = tok.find('=');
        if (eq == std::string::npos) continue;
        out.emplace_back(tok.substr(0, eq), tok.substr(eq + 1));
    }
    return out;
}

std::string kv_get(const KV& kv, const std::string& key)
{
    for (const auto& p : kv)
        if (p.first == key) return p.second;
    return std::string();
}

bool g_first_scene_printed = false;
bool g_first_scene_printed_push_warn = false; // предупреждение о --impulse-push-tau -- один раз

// Неупругость нити и eps_break (после build_guo_scene). Возвращает пары ключ=значение для shot.txt.
// Файл директив разбирается штатным load_fiber_file: временный .fib = MATERIAL 0 сцены + директивы файла + заглушка FIBER
// (парсер требует хотя бы одну нить); plastic/rate/visco/failure копируются во все материалы сцены.
KV apply_fiber_overrides(const Options& opt, GuoScene& sc, const fs::path& workdir)
{
    KV info;
    // Модуль Юнга -- первым: временный материал директив неупругости ниже строится уже с ним.
    // Сцена (геометрия, массы, l0) от E не зависит; E читается только в счёте.
    if (opt.fiber_E > 0.0)
        for (fiber::FiberMaterial& m : sc.model.materials) m.E = opt.fiber_E;
    if (!opt.fiber_inelastic.empty())
    {
        std::ifstream in(opt.fiber_inelastic);
        if (!in) throw std::runtime_error("не открыть " + opt.fiber_inelastic.string());
        std::string body, ln;
        while (std::getline(in, ln))
        {
            if (!ln.empty() && ln.back() == '\r') ln.pop_back();
            const std::size_t k = ln.find_first_not_of(" \t");
            if (k != std::string::npos && (ln.compare(k, 8, "MATERIAL") == 0 || ln.compare(k, 5, "FIBER") == 0))
                throw std::runtime_error("--fiber-inelastic: файл должен содержать только директивы PLASTIC/RATE/VISCO/FAILURE "
                                         "(без MATERIAL/FIBER): " + ln);
            body += ln + "\n";
        }
        const fiber::FiberMaterial& m0 = sc.model.materials.at(0);
        char buf[512];
        const fs::path p = workdir / "fiber_inelastic_material.fib";
        {
            std::ofstream f(p);
            std::snprintf(buf, sizeof buf, "MATERIAL 0 %.17g %.17g %.17g %.17g %.17g\n", m0.E, m0.nu, m0.d, m0.rho,
                          m0.eps_break);
            f << buf << body << "FIBER 0 2\n0 0 0\n1e-3 0 0\n";
        }
        fiber::Model tmp;
        fiber::load_fiber_file(p.string(), tmp);
        if (tmp.materials.empty() || !(tmp.materials.at(0).inelastic() || tmp.materials.at(0).has_failure()))
            throw std::runtime_error("--fiber-inelastic: директивы не включили неупругость/критерий разрыва материала 0");
        for (fiber::FiberMaterial& m : sc.model.materials)
        {
            m.plastic = tmp.materials[0].plastic;
            m.rate = tmp.materials[0].rate;
            m.visco = tmp.materials[0].visco;
            m.failure = tmp.materials[0].failure;
            m.damage = tmp.materials[0].damage;
        }
        const fiber::FiberMaterial& m = sc.model.materials.at(0);
        info.emplace_back("fiber_inelastic", opt.fiber_inelastic.filename().string());
        info.emplace_back("pl_sy", fmtg(m.plastic.sy));
        info.emplace_back("pl_H", fmtg(m.plastic.H));
        info.emplace_back("rate_model", std::to_string(static_cast<long long>(m.rate.model)));
        info.emplace_back("visco_eta", fmtg(m.visco.eta));
        info.emplace_back("fail_eps_p", fmtg(m.failure.eps_p));
        info.emplace_back("fail_sigma_u", fmtg(m.failure.sigma_u));
        if (m.damage.enabled())
        {
            info.emplace_back("damage_deps_f", fmtg(m.damage.deps_f));
            info.emplace_back("damage_gf", fmtg(m.damage.gf));
        }
    }
    if (opt.eps_break > 0.0)
    {
        for (fiber::FiberMaterial& m : sc.model.materials) m.eps_break = opt.eps_break;
        // l_break лежит в связи: (1 + eps_break) l0, как при загрузке .fib
        for (fiber::FiberBond& b : sc.model.topology.bonds) b.l_break = (1.0 + opt.eps_break) * b.l0;
    }
    info.emplace_back("fiber_eps_break", fmtg(sc.model.materials.at(0).eps_break));
    info.emplace_back("fiber_E", fmtg(sc.model.materials.at(0).E));
    info.emplace_back("impulse_push_tau", fmtg(opt.impulse_push_tau));
    info.emplace_back("impulse_push_depth", fmtg(opt.impulse_push_depth));
    info.emplace_back("impulse_push_vmax", fmtg(opt.impulse_push_vmax));
    info.emplace_back("impulse_gamma", fmtg(opt.impulse_gamma > 0.0 ? opt.impulse_gamma : 1.0));
    info.emplace_back("dt_scale", fmtg(opt.dt_scale));
    // Ключи проекции нормального контакта пишутся в shot.txt отдельно (по sc.cfg после применения CLI,
    // см. run_shot): здесь cfg выстрела ещё не заполнен, а сама функция зовётся только при переопределениях нити.
    return info;
}

ShotRecord run_shot(const Options& opt, double V0, int rep, const std::string& note)
{
    ShotRecord rec;
    rec.V0 = V0;
    rec.rep = rep;
    const std::string cname = config_name_str(opt.config, opt.plies);
    const std::string lname = level_name(opt.level);
    const std::string nz_tag = opt.plate_nz > 0 ? "_nz" + std::to_string(opt.plate_nz) : std::string();
    // Шар не по умолчанию -- свой каталог, чтобы не пересечься с обычными выстрелами
    std::string ball_tag;
    if (opt.ball_d > 0.0 || opt.ball_m > 0.0 || opt.ball_dx > 0.0)
    {
        const LevelSpec L0 = level_spec(opt.level);
        const double bd = opt.ball_d > 0.0 ? opt.ball_d : L0.D_s, bm = opt.ball_m > 0.0 ? opt.ball_m : L0.m_s,
                     bx = opt.ball_dx > 0.0 ? opt.ball_dx : L0.dx_s;
        ball_tag = "_bD" + fmtg(bd * 1e3, 6) + "m" + fmtg(bm * 1e3, 6) + "dx" + fmtg(bx * 1e3, 6);
    }
    const fs::path sub =
        fs::absolute(opt.out) / (cname + "_" + lname + nz_tag + ball_tag + "_v" + fmt_v0(V0) + "_r" + std::to_string(rep));
    rec.dir = sub;
    try
    {
        fs::create_directories(sub);
        CwdGuard cwd(sub);

        // ---- сцена
        const auto tb0 = std::chrono::steady_clock::now();
        GuoScene sc;
        {
            e2e::ScopedSilence silence(true);
            sc = build_guo_scene(opt.config, opt.level, V0, sub, opt.mu, opt.body_isolation, opt.plate_nz,
                                  opt.plies, opt.ball_d, opt.ball_m, opt.ball_dx);
        }
        const double build_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - tb0).count();
        rec.m_s = sc.ball_mass;
        // Раздельное трение: сцена построена с --mu для обоих; --mu-ff/--mu-fs переопределяют материалы ткани.
        const double mu_ff = opt.mu_ff >= 0.0 ? opt.mu_ff : opt.mu;
        const double mu_fs = opt.mu_fs >= 0.0 ? opt.mu_fs : opt.mu;
        if (sc.has_fabric && (opt.mu_ff >= 0.0 || opt.mu_fs >= 0.0))
            for (fiber::FiberMaterial& m : sc.model.materials)
            {
                m.mu_ff = mu_ff;
                m.mu_fs = mu_fs;
            }
        KV fiber_info;
        if (sc.has_fabric && (!opt.fiber_inelastic.empty() || opt.eps_break > 0.0 || opt.fiber_E > 0.0))
        {
            fiber_info = apply_fiber_overrides(opt, sc, sub);
            std::printf("[нить] неупругость/eps_break:");
            for (const auto& kvp : fiber_info) std::printf(" %s=%s", kvp.first.c_str(), kvp.second.c_str());
            std::printf("\n");
        }
        if (!g_first_scene_printed)
        {
            g_first_scene_printed = true;
            print_scene_summary(stdout, sc);
        }

        const double t_max = opt.t_max > 0.0 ? opt.t_max : default_t_max(opt.config, opt.level);
        sc.cfg.t_total = t_max;
        if (opt.fs_contact >= 0) sc.cfg.fiber_solid_contact_type = opt.fs_contact;
        if (opt.ff_contact >= 0) sc.cfg.fiber_fiber_contact_type = opt.ff_contact;
        sc.cfg.fiber_solid_contact_impulse_push_tau = sc.cfg.fiber_fiber_contact_impulse_push_tau = opt.impulse_push_tau;
        if (opt.impulse_push_tau > 0.0 && !g_first_scene_printed_push_warn)
        {
            g_first_scene_printed_push_warn = true;
            std::fprintf(stderr, "ПРЕДУПРЕЖДЕНИЕ: --impulse-push-tau включает выталкивание закона 3 -- источник энергии "
                                 "(выбросы до -37 %% E0, ложные разрывы; fiber_contact_law.hpp, FiberContactLaw::push_tau). "
                                 "Только для воспроизведения серий 2026-10-20/21.\n");
        }
        sc.cfg.fiber_contact_impulse_push_depth = opt.impulse_push_depth;
        sc.cfg.fiber_contact_impulse_push_vmax = opt.impulse_push_vmax;
        if (opt.impulse_gamma > 0.0)
            sc.cfg.fiber_solid_contact_impulse_gamma = sc.cfg.fiber_fiber_contact_impulse_gamma = opt.impulse_gamma;
        sc.cfg.dt_scale = opt.dt_scale;
        if (opt.ff_predictor >= 0) sc.cfg.fiber_fiber_contact_predictor = opt.ff_predictor;
        if (opt.fs_predictor >= 0) sc.cfg.fiber_solid_contact_predictor = opt.fs_predictor;
        if (opt.proj_iters >= 0) sc.cfg.fiber_contact_projection_iters = opt.proj_iters;
        if (opt.proj_warm >= 0) sc.cfg.fiber_contact_projection_warm = opt.proj_warm;
        validate_fiber_contact_projection(sc.cfg); // ошибка настроек -- исключение выстрела (ОШИБКА)
#ifdef SPH_CONTACT_STATS
        // Статистика контактов -- в подкаталог выстрела (CWD).
        sc.cfg.contact_stats_every = opt.contact_stats;
        sc.cfg.contact_hold_n = opt.hold_n;
#endif
        // Кадр на шаге 0 пишется всегда (next_output_time стартует с 0), поэтому при frames=off
        // save_t огромен: после него кадров больше не будет.
        sc.cfg.save_t = opt.save_t > 0.0 ? opt.save_t : (opt.frames_camera ? kBaseSaveT : 1.0e30);
        sc.cfg.print_step = 1000000000;
#ifdef SPH_ENERGY_LOG
        // Лог энергии -- в подкаталог выстрела (CWD), include/energy_log.hpp.
        sc.cfg.energy_log_every = opt.energy_log;
#endif
        // cfg.num_threads конвейер не читает (его применяет только main.cpp), поэтому
        // число потоков задаётся напрямую: CpuExecutor берёт omp_get_max_threads().
        if (opt.threads > 0)
        {
            sc.cfg.num_threads = opt.threads;
#ifdef _OPENMP
            omp_set_num_threads(opt.threads);
#endif
        }

        // Чекпоинты/рестарт: файл -- в каталоге выстрела (абсолютный путь, CWD = sub). При --restart
        // dt/time ниже перекрываются состоянием цикла из чекпоинта (time_integration).
        const fs::path ckpt_path = sub / "checkpoint.ckpt";
        // Идентификатор сцены (число слоёв и шар): параметры шара и слоёв не входят в cfg чекпоинта, поэтому
        // сверяются отдельным текстовым файлом; формат чекпоинта не меняется, старые чекпоинты без файла допускаются.
        const std::string scene_id = "config=" + config_name_str(opt.config, sc.n_plies) + " plies=" + std::to_string(sc.n_plies) +
                                     " ball_D=" + fmtg(sc.spec.D_s, 17) + " ball_m=" + fmtg(sc.spec.m_s, 17) +
                                     " ball_dx=" + fmtg(sc.spec.dx_s, 17) + "\n";
        const fs::path id_path = sub / "scene_id.txt";
        if (opt.checkpoint && !opt.restart)
        {
            std::ofstream idf(id_path);
            idf << scene_id;
        }
        if (opt.restart && fs::exists(id_path))
        {
            std::ifstream idf(id_path);
            const std::string saved((std::istreambuf_iterator<char>(idf)), std::istreambuf_iterator<char>());
            if (saved != scene_id)
                throw std::runtime_error("--restart: сцена не совпадает с сохранённой (число слоёв/шар):\n  сохранено: " + saved +
                                         "  сейчас:    " + scene_id);
        }
        if (opt.checkpoint)
        {
            sc.cfg.checkpoint_file = ckpt_path.string();
            sc.cfg.checkpoint_every_steps = opt.checkpoint_every;
        }
        if (opt.restart)
        {
            const fs::path prev = ckpt_path.string() + ".prev";
            if (!fs::exists(ckpt_path) && !fs::exists(prev))
                throw std::runtime_error("--restart: нет чекпоинта " + ckpt_path.string() + " (и .prev) в каталоге выстрела");
            sc.cfg.restart_file = (fs::exists(ckpt_path) ? ckpt_path : prev).string();
            sc.cfg.restart_force = opt.force;
        }

        double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal) * sc.cfg.dt_scale, time = 0.0;
        RunGuard guard(opt.time_limit, 2.0e6, 200); // time_limit <= 0 -- без лимита по стенному времени
        if (opt.max_steps > 0) guard.set_step_limit(opt.max_steps);
        MonitorOptions mo;
        mo.n_mon = opt.n_mon;
        mo.n_diag = opt.n_diag;
        mo.t_max = t_max;
        mo.csv_path = (sub / "monitor.csv").string();
        mo.guard = &guard;
        mo.restart = opt.restart;
        GuoShotMonitor monitor(sc, mo);
        guard.set_step_observer([&](double d, double t) { monitor.on_step(d, t); });
        // Состояние монитора -- в блоке USER чекпоинта (только при --checkpoint / --restart).
        checkpoint::UserState mon_state;
        mon_state.save = [&monitor] { return monitor.save_state(); };
        mon_state.restore = [&monitor](const std::string& b, std::size_t n) { monitor.restore_state(b, n); };
        const checkpoint::UserState* mon_state_ptr = (opt.checkpoint || opt.restart) ? &mon_state : nullptr;

        // ---- прогон
        const auto w0 = std::chrono::steady_clock::now();
        {
            e2e::ScopedSilence silence(true);
            time_integration(sc.cfg, sc.model, sc.ntotal, dt, time, sc.spec.dx_s, &guard, mon_state_ptr);
        }
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
        monitor.finalize();
        const MonitorStats& st = monitor.stats();

        // ---- кадры: кадр шага 0 есть всегда; --probe-frame добавляет кадр конечного состояния
        FrameSizes fz0 = scan_frame_sizes(sub);
        FrameSizes probe;
        if (opt.probe_frame == 1)
        {
            e2e::ScopedSilence silence(true);
            output(sc.model.solid, sc.ntotal, time, 99999, dt);
            if (sc.has_fabric) fiber::fiber_output(sc.model, 99999, time);
            const FrameSizes fz1 = scan_frame_sizes(sub);
            probe.bytes_sph = fz1.bytes_sph - fz0.bytes_sph;
            probe.bytes_fiber = fz1.bytes_fiber - fz0.bytes_fiber;
        }
        const FrameSizes fz = scan_frame_sizes(sub);
        // размер одного кадра: из пробного, иначе из средних по записанным кадрам
        const double frame_sph = opt.probe_frame == 1 ? static_cast<double>(probe.bytes_sph)
                                 : (fz0.n_sph ? static_cast<double>(fz0.bytes_sph) / fz0.n_sph : 0.0);
        const double frame_fiber = opt.probe_frame == 1 ? static_cast<double>(probe.bytes_fiber)
                                   : (fz0.n_fiber ? static_cast<double>(fz0.bytes_fiber) / fz0.n_fiber : 0.0);
        if (!opt.frames_camera && !(opt.save_t > 0.0) && opt.probe_frame != 1)
        {
            // кадры выключены: служебный кадр шага 0 после замера размера удаляем
            std::error_code ec;
            for (const auto& e : fs::directory_iterator(sub, ec))
            {
                const std::string n = e.path().filename().string();
                if (n.rfind("sph_output_", 0) == 0 || n.rfind("fiber_output_", 0) == 0) fs::remove(e.path());
            }
        }

        // StepLimit от request_stop монитора -- штатный останов по исходу или t_max, не от --max-steps
        const std::string greason =
            (st.stopped_by_monitor && guard.reason() == RunGuard::Reason::StepLimit) ? "Monitor" : reason_name(guard.reason());
        monitor.write_shot_txt((sub / "shot.txt").string(), wall, greason);

        // ---- стоимость: дописываем ключи в shot.txt
        const double steps = static_cast<double>(st.steps);
        const double per_step = steps > 0.0 ? wall / steps : 0.0;
        const double dt_mean = steps > 0.0 ? st.t_end / steps : 0.0;
        const double steps_fc = dt_mean > 0.0 ? t_max / dt_mean : 0.0;
        const double wall_fc_h = steps_fc * per_step / 3600.0;
        const double mon_share = wall > 0.0 ? 100.0 * st.monitor_wall_s / wall : 0.0;
        std::string line;
        {
            std::ifstream in(sub / "shot.txt");
            std::getline(in, line);
        }
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r' || line.back() == ' ')) line.pop_back();
        auto add = [&](const std::string& k, const std::string& v) { line += " " + k + "=" + v; };
        add("rep", std::to_string(rep));
        add("t_max", fmtg(t_max));
        add("ntotal", std::to_string(sc.ntotal));
        add("n_ball", std::to_string(sc.n_ball));
        add("n_plate", std::to_string(sc.n_plate));
        add("plate_nz", std::to_string(sc.plate_nz));
        add("n_beads", std::to_string(sc.actual_beads));
        add("n_bonds", std::to_string(sc.actual_bonds));
        add("n_mon", std::to_string(opt.n_mon));
        add("n_diag", std::to_string(opt.n_diag));
        add("scene_build_s", fmtg(build_s, 6));
        add("wall_per_step_s", fmtg(per_step, 6));
        add("dt_mean", fmtg(dt_mean));
        add("steps_forecast", fmtg(steps_fc, 8));
        add("wall_forecast_h", fmtg(wall_fc_h, 6));
        add("monitor_wall_s", fmtg(st.monitor_wall_s, 6));
        add("monitor_share_pct", fmtg(mon_share, 4));
        // Только при --restart: restored -- вердикт и максимумы за весь выстрел; partial -- только с шага рестарта.
        if (opt.restart) add("monitor_state", monitor.state_restored() ? "restored" : "partial");
        add("frames_sph_n", std::to_string(fz.n_sph));
        add("frames_fiber_n", std::to_string(fz.n_fiber));
        add("frames_bytes_total", std::to_string(fz.bytes_sph + fz.bytes_fiber));
        add("frame_sph_bytes", fmtg(frame_sph, 10));
        add("frame_fiber_bytes", fmtg(frame_fiber, 10));
        for (const auto& kvp : fiber_info) add(kvp.first, kvp.second);
        // Режим проекции (любой вид контакта): ключи из ИТОГОВОЙ конфигурации sc.cfg (после CLI); при
        // выключенном режиме строка shot.txt не меняется. Место в строке -- прежнее (перед mu_ff).
        if (sc.cfg.fiber_fiber_contact_predictor == 1 || sc.cfg.fiber_solid_contact_predictor == 1)
        {
            add("ff_contact_predictor", std::to_string(sc.cfg.fiber_fiber_contact_predictor == 1 ? 1 : 0));
            add("fs_contact_predictor", std::to_string(sc.cfg.fiber_solid_contact_predictor == 1 ? 1 : 0));
            add("contact_proj_iters", std::to_string(sc.cfg.fiber_contact_projection_iters));
            add("contact_proj_warm", std::to_string(sc.cfg.fiber_contact_projection_warm));
        }
        add("mu_ff", fmtg(mu_ff));
        add("mu_fs", fmtg(mu_fs));
        if (opt.ball_d > 0.0 || opt.ball_m > 0.0 || opt.ball_dx > 0.0) // шар не по умолчанию
        {
            add("ball_D", fmtg(sc.spec.D_s));
            add("ball_m", fmtg(sc.spec.m_s));
            add("ball_dx", fmtg(sc.spec.dx_s));
            add("ball_rho", fmtg(sc.rho_s));
        }
        {
            std::ofstream outf(sub / "shot.txt");
            outf << line << "\n";
        }

        rec.kv = parse_kv_line(line);
        rec.outcome = st.outcome;
        rec.V_R = st.V_R;
        rec.e = monitor.shot_energy();
        rec.ok = true;

        // ---- вывод в консоль
        const bool det = st.outcome == Outcome::Perforated || st.outcome == Outcome::NotPerforated;
        std::printf("[выстрел] %s %s%s V0=%s r=%d: %s, V_R=%s, e=%s Дж, шагов=%zu, wall=%.2f с (%.4f с/шаг), "
                    "dt_first=%.3e с, монитор %.2f %%, сторож=%s%s\n",
                    cname.c_str(), lname.c_str(), nz_tag.c_str(), fmt_v0(V0).c_str(), rep, outcome_name(st.outcome),
                    det ? fmtg(st.V_R, 6).c_str() : "-", det ? fmtg(rec.e, 6).c_str() : "-", st.steps, wall, per_step,
                    st.dt_first, mon_share, greason.c_str(), note.c_str());
        std::printf("          стоимость: ntotal=%zu n_ball=%zu n_plate=%zu бусин=%zu связей=%zu; сцена %.2f с; dt_first=%.3e "
                    "dt_min=%.3e dt_mean=%.3e; прогноз: %.3g шагов до t_max=%.0f мкс, %.3g ч; кадр: Solid %.3g МБ, нити "
                    "%.3g МБ; кадров %zu+%zu\n",
                    sc.ntotal, sc.n_ball, sc.n_plate, sc.actual_beads, sc.actual_bonds, build_s, st.dt_first, st.dt_min,
                    dt_mean, steps_fc, t_max * 1e6, wall_fc_h, frame_sph / 1048576.0, frame_fiber / 1048576.0, fz.n_sph,
                    fz.n_fiber);
        if (det || st.steps > 0) monitor.print_energy_decomposition(stdout);
    }
    catch (const std::exception& ex)
    {
        rec.ok = false;
        rec.error = ex.what();
        std::printf("[выстрел] %s %s V0=%s r=%d: ОШИБКА: %s\n", cname.c_str(), lname.c_str(), fmt_v0(V0).c_str(), rep,
                    ex.what());
    }
    return rec;
}

// ---------------------------------------------------------------------------
// V_BL (ТЗ Guo 2019 §8.3)
// ---------------------------------------------------------------------------
struct VblResult
{
    bool has_n = false, has_p = false;
    double max_n = 0.0, min_p = 0.0;
    std::size_t n_det = 0, n_undet = 0, n_used = 0;
    bool defined = false;
    double e_bar = 0.0, V_BL = 0.0;
    double m_s = 0.0;
};

VblResult compute_vbl(const std::vector<ShotRecord>& recs)
{
    VblResult r;
    for (const auto& s : recs)
    {
        if (!s.ok) continue;
        r.m_s = s.m_s;
        if (s.outcome == Outcome::NotPerforated)
        {
            ++r.n_det;
            if (!r.has_n || s.V0 > r.max_n) r.max_n = s.V0;
            r.has_n = true;
        }
        else if (s.outcome == Outcome::Perforated)
        {
            ++r.n_det;
            if (!r.has_p || s.V0 < r.min_p) r.min_p = s.V0;
            r.has_p = true;
        }
        else
            ++r.n_undet;
    }
    if (r.has_n && r.has_p && r.m_s > 0.0)
    {
        double sum = 0.0;
        for (const auto& s : recs)
            if (s.ok && s.outcome == Outcome::Perforated && s.V0 > r.max_n && s.V0 <= r.max_n + 15.0)
            {
                sum += s.e;
                ++r.n_used;
            }
        if (r.n_used > 0)
        {
            r.e_bar = sum / static_cast<double>(r.n_used);
            r.V_BL = std::sqrt(2.0 * r.e_bar / r.m_s);
            r.defined = true;
        }
    }
    return r;
}

void write_vbl(const Options& opt, const std::vector<ShotRecord>& recs, std::FILE* extra)
{
    const VblResult v = compute_vbl(recs);
    const PaperVBL pv = paper_vbl(opt.config);
    auto emit = [&](std::FILE* f) {
        std::fprintf(f, "config=%s level=%s mu=%g body_isolation=%d\n", config_name_str(opt.config, opt.plies).c_str(), level_name(opt.level),
                     opt.mu, opt.body_isolation ? 1 : 0);
        std::fprintf(f, "выстрелов с определённым исходом: %zu, не определён (исключены): %zu\n", v.n_det, v.n_undet);
        if (v.has_n) std::fprintf(f, "наибольшая «н» (не пробит): %.1f м/с\n", v.max_n);
        else std::fprintf(f, "наибольшая «н» (не пробит): нет\n");
        if (v.has_p) std::fprintf(f, "наименьшая «п» (пробит): %.1f м/с\n", v.min_p);
        else std::fprintf(f, "наименьшая «п» (пробит): нет\n");
        if (v.has_n && v.has_p)
            std::fprintf(f, "интервал [н; п] = [%.1f; %.1f] м/с%s\n", v.max_n, v.min_p,
                         v.min_p < v.max_n ? " (немонотонный: «п» ниже «н»)" : "");
        if (v.defined)
        {
            std::fprintf(f, "ē = %.6g Дж по %zu пробивным выстрелам в пределах 15 м/с над %.1f, m_s = %.6g кг\n", v.e_bar,
                         v.n_used, v.max_n, v.m_s);
            std::fprintf(f, "V_BL = %.2f м/с (статья: %.0f +%.0f/-%.0f)\n", v.V_BL, pv.value, pv.plus, pv.minus);
        }
        else
            std::fprintf(f, "V_BL не определён (%s)\n",
                         !v.has_p ? "нет пробитых выстрелов"
                         : !v.has_n ? "нет непробитых выстрелов"
                                    : "нет пробитых выстрелов в пределах 15 м/с над наибольшей «н»");
        // 2T/4T: V_BL и e_BL статьи (ТЗ 2T/4T §4.5); для остальных конфигураций формат vbl.txt прежний
        if (opt.config == Config::T2 || opt.config == Config::T4)
            std::fprintf(f, "V_BL статьи = %.0f м/с, e_BL статьи = %.1f Дж\n", pv.value, paper_e_bl(opt.config));
    };
    std::FILE* f = std::fopen((fs::absolute(opt.out) / "vbl.txt").string().c_str(), "w");
    if (f)
    {
        emit(f);
        std::fclose(f);
    }
    if (extra) emit(extra);
}

void write_shots_csv(const Options& opt, const std::vector<ShotRecord>& recs)
{
    const ShotRecord* first = nullptr;
    for (const auto& s : recs)
        if (s.ok) { first = &s; break; }
    if (!first) return;
    std::FILE* f = std::fopen((fs::absolute(opt.out) / "shots.csv").string().c_str(), "w");
    if (!f) return;
    for (std::size_t k = 0; k < first->kv.size(); ++k) std::fprintf(f, "%s%s", k ? "," : "", first->kv[k].first.c_str());
    std::fprintf(f, "\n");
    for (const auto& s : recs)
    {
        if (!s.ok) continue;
        for (std::size_t k = 0; k < first->kv.size(); ++k)
            std::fprintf(f, "%s%s", k ? "," : "", kv_get(s.kv, first->kv[k].first).c_str());
        std::fprintf(f, "\n");
    }
    std::fclose(f);
}

bool has_tried(const std::vector<ShotRecord>& recs, double v)
{
    for (const auto& s : recs)
        if (std::fabs(s.V0 - v) < 1.0e-6) return true;
    return false;
}

std::string paper_note(const std::vector<PaperShot>& ps, const ShotRecord& r)
{
    std::string out;
    bool any = false;
    for (const auto& p : ps)
        if (std::fabs(p.V0 - r.V0) < 1.0e-6)
        {
            any = true;
            out += std::string(out.empty() ? "" : "; ") + (p.perforated ? "п" : "н") +
                   (p.perforated ? " V_R=" + fmtg(p.V_R, 4) : std::string());
        }
    if (!any) return std::string();
    const bool det = r.outcome == Outcome::Perforated || r.outcome == Outcome::NotPerforated;
    const char* model = !det ? "?" : (r.outcome == Outcome::Perforated ? "п" : "н");
    // совпадение исхода с хотя бы одним опытом при этой скорости
    bool match = false;
    for (const auto& p : ps)
        if (std::fabs(p.V0 - r.V0) < 1.0e-6 && det && p.perforated == (r.outcome == Outcome::Perforated)) match = true;
    return std::string(" | опыт: ") + out + " | модель: " + model + (r.outcome == Outcome::Perforated ? " V_R=" + fmtg(r.V_R, 4) : "") +
           (!det ? " | ИСХОД НЕ ОПРЕДЕЛЁН" : match ? " | ИСХОД СОВПАЛ" : " | ИСХОД НЕ СОВПАЛ");
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    Options opt;
    const int pr = parse_args(argc, argv, opt);
    if (pr == 1)
    {
        print_usage();
        return 0;
    }
    if (pr != 0)
    {
        print_usage();
        return 2;
    }

    std::error_code ec;
    fs::create_directories(opt.out, ec);
    if (ec)
    {
        std::fprintf(stderr, "Не удалось создать каталог %s: %s\n", opt.out.string().c_str(), ec.message().c_str());
        return 2;
    }

    const std::vector<PaperShot> paper = paper_shots(opt.config);
    std::vector<double> v_first;
    if (opt.paper_shots)
    {
        for (const auto& p : paper)
            if (std::find_if(v_first.begin(), v_first.end(), [&](double v) { return std::fabs(v - p.V0) < 1e-6; }) ==
                v_first.end())
                v_first.push_back(p.V0); // дубликаты скоростей статьи (78 и 78) считаются один раз
    }
    else
        v_first = opt.v0_list;

    std::printf("=== guo2019_shot: config=%s level=%s, выстрелов в серии %zu, mu=%g, изоляция тел %s, кадры %s, "
                "t_max=%.0f мкс, time-limit=%s, n_mon=%zu, n_diag=%zu, repeat=%d, bracket=%s, probe-frame=%s ===\n",
                config_name_str(opt.config, opt.plies).c_str(), level_name(opt.level), v_first.size(), opt.mu,
                opt.body_isolation ? "вкл." : "ВЫКЛ.", opt.frames_camera ? "camera" : "off",
                (opt.t_max > 0.0 ? opt.t_max : default_t_max(opt.config, opt.level)) * 1e6,
                opt.time_limit > 0.0 ? fmtg(opt.time_limit, 6).c_str() : "нет", opt.n_mon, opt.n_diag, opt.repeat,
                opt.bracket ? "да" : "нет", opt.probe_frame ? "да" : "нет");

    std::vector<ShotRecord> recs;
    int n_err = 0;
    auto shoot = [&](double v, int rep, const std::string& tag) {
        ShotRecord r = run_shot(opt, v, rep, tag);
        if (!r.ok) ++n_err;
        const bool is_paper = opt.paper_shots && rep == 0;
        if (r.ok && is_paper)
        {
            const std::string pn = paper_note(paper, r);
            if (!pn.empty()) std::printf("          сравнение с опытом%s\n", pn.c_str());
        }
        recs.push_back(std::move(r));
        write_shots_csv(opt, recs);
        write_vbl(opt, recs, nullptr);
    };

    // ---- 1. основная серия
    for (double v : v_first) shoot(v, 0, "");

    // ---- 2. расширение бисекцией (ТЗ Guo 2019 §8.2): до 6 выстрелов, пока интервал [н; п] > 5 м/с
    if (opt.bracket)
    {
        std::printf("\n--- bracket: уточнение границы «пробит / не пробит» (до 6 дополнительных выстрелов) ---\n");
        for (int extra = 0; extra < 6; ++extra)
        {
            const VblResult v = compute_vbl(recs);
            if (v.n_det == 0)
            {
                std::printf("bracket: нет выстрелов с определённым исходом, расширять не от чего\n");
                break;
            }
            double next = 0.0;
            const char* why = "";
            if (!v.has_n)
            {
                next = v.min_p - 20.0; // все «п»: идём вниз
                why = "все пробиты: вниз на 20 м/с";
                if (!(next > 0.0)) break;
            }
            else
            {
                // lo -- наибольшая «н»; hi -- наименьшая «п» выше неё
                const double lo = v.max_n;
                double hi = std::numeric_limits<double>::infinity();
                for (const auto& s : recs)
                    if (s.ok && s.outcome == Outcome::Perforated && s.V0 > lo) hi = std::min(hi, s.V0);
                if (!std::isfinite(hi))
                {
                    next = lo + 20.0; // все «н»: идём вверх
                    why = "все не пробиты: вверх на 20 м/с";
                }
                else if (hi - lo <= 5.0)
                {
                    std::printf("bracket: интервал [%.1f; %.1f] <= 5 м/с, достаточно\n", lo, hi);
                    break;
                }
                else
                {
                    next = std::round(0.5 * (lo + hi) * 10.0) / 10.0;
                    why = "бисекция";
                }
            }
            if (has_tried(recs, next))
            {
                std::printf("bracket: V0=%.1f уже считался (исход не определён), останов\n", next);
                break;
            }
            std::printf("bracket %d/6: V0=%.1f (%s)\n", extra + 1, next, why);
            shoot(next, 0, " [bracket]");
        }
    }

    // ---- 3. повторы около V_BL модели (ТЗ Guo 2019 §8.3)
    if (opt.repeat > 1)
    {
        const VblResult v = compute_vbl(recs);
        if (v.defined)
        {
            std::vector<double> near;
            for (const auto& s : recs)
                if (s.ok && s.rep == 0 && std::fabs(s.V0 - v.V_BL) <= 5.0 &&
                    std::find_if(near.begin(), near.end(), [&](double x) { return std::fabs(x - s.V0) < 1e-6; }) == near.end())
                    near.push_back(s.V0);
            std::printf("\n--- повторы: %zu скоростей в пределах 5 м/с от V_BL модели %.2f м/с, по %d доп. выстрелу ---\n",
                        near.size(), v.V_BL, opt.repeat - 1);
            for (double x : near)
                for (int r = 1; r < opt.repeat; ++r) shoot(x, r, " [повтор]");
        }
        else
            std::printf("\n--- повторы пропущены: V_BL модели не определён ---\n");
    }

    std::printf("\n=== Итог: %s %s ===\n", config_name_str(opt.config, opt.plies).c_str(), level_name(opt.level));
    write_shots_csv(opt, recs);
    write_vbl(opt, recs, stdout);
    // Послойные разрывы последнего пробивного выстрела (ТЗ 2T/4T §4.5), только для пакета из >= 2 слоёв
    if (fabric_plies(opt.config, opt.plies) >= 2)
    {
        const ShotRecord* last = nullptr;
        for (const auto& s : recs)
            if (s.ok && s.outcome == Outcome::Perforated) last = &s;
        if (!last)
            std::printf("\nПослойные разрывы: пробивных выстрелов в серии нет\n");
        else
        {
            const int np = std::atoi(kv_get(last->kv, "n_plies").c_str());
            std::printf("\nПослойные разрывы последнего пробивного выстрела (V0=%.1f м/с, r=%d; слой 1 -- лицевой):\n", last->V0,
                        last->rep);
            for (int k = 1; k <= np; ++k)
            {
                const std::string q = "_p" + std::to_string(k);
                std::printf("  слой %d: порвано вдоль y (зажатые) %s, вдоль x (свободные) %s, t_первого разрыва %s с, прогиб "
                            "слоя при нём %s м\n",
                            k, kv_get(last->kv, "n_broken_v" + q).c_str(), kv_get(last->kv, "n_broken_h" + q).c_str(),
                            kv_get(last->kv, "t_first_break" + q).c_str(), kv_get(last->kv, "w_f_at_first_break" + q).c_str());
            }
        }
    }
    std::printf("shots.csv и vbl.txt: %s\n", fs::absolute(opt.out).string().c_str());
    if (n_err > 0) std::printf("ВНИМАНИЕ: выстрелов с исключением: %d\n", n_err);
    return 0;
}
