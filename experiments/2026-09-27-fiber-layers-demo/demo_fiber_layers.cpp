// experiments/2026-09-27-fiber-layers-demo/demo_fiber_layers.cpp
//
// Демонстрационный расчёт: менее массивный снаряд из SPH-частиц (алюминий)
// летит в сеть из 4 слоёв нитей (направления x,y,x,y сверху вниз по ходу
// снаряда, каждый слой -- несколько ПАРАЛЛЕЛЬНЫХ нитей, как ткань). Два
// сценария с ОДИНАКОВОЙ сетью и ОДИНАКОВЫМ снарядом, различающиеся только
// начальной скоростью (и, при необходимости, t_total): "pierce" -- снаряд
// прорывает все 4 слоя, "stop" -- снаряд останавливается/отскакивает.
//
// По образцу experiments/2026-09-26-fiber-net-demo/demo_fiber_net.cpp: та же
// архитектура (write_net_file -> fiber::load_fiber_file -> псевдо-анкеры ->
// SimulationConfig -> chdir -> time_integration(cfg, fiber::Model&, ...)).
// Это НЕ тест: add_test() не добавляется (см. CMakeLists.txt).
//
// Запуск: demo_fiber_layers [output_base_dir] [pierce|stop]
//   output_base_dir -- каталог, ВНУТРИ которого создаётся подкаталог сценария
//   (output_base_dir/pierce или output_base_dir/stop); туда пишутся net.fib,
//   net_layout.txt (метаданные слоёв для анализа в Python) и кадры
//   sph_output_*.txt/fiber_output_*.txt. По умолчанию --
//   experiments/2026-09-27-fiber-layers-demo/data, сценарий -- pierce.
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "fiber/fiber_input.hpp"
#include "fiber/model.hpp"
#include "int_force.hpp" // urs()
#include "materials.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "run_guard.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"

#ifdef _WIN32
#include <windows.h>
#endif

static_assert(dim == 3, "demo_fiber_layers требует сборки с SPH_DIM=3 (нити -- только в 3D, notes/fiber-model/fiber-problem-statement.md)");

namespace fs = std::filesystem;

namespace
{

// ---------------------------------------------------------------------
// Сеть -- 4 слоя по z, направления (сверху вниз, по ходу снаряда вдоль -z):
// x, y, x, y. Каждый слой -- fibers_per_layer ПАРАЛЛЕЛЬНЫХ прямых нитей,
// разнесённых по перпендикулярной в плоскости слоя оси на fiber_spacing
// (> d, иначе соседи в слое касались бы уже на старте). Соседние слои по z
// разнесены РОВНО на d (касание без проникания на старте -- как в старом
// демо). Все нити центрированы по (x=0,y=0) в своей паре координат --
// траектория снаряда проходит через центр сети.
//
// FIBER-блоки пишутся по слоям подряд (все 5 нитей слоя 0, затем слоя 1, и
// т.д.), поэтому object_id бусины в кадре 0 (fiber::load_fiber_file
// присваивает object_id = порядковый номер FIBER-блока, см. src/
// fiber_input.cpp) напрямую даёт номер слоя: layer = object_id_at_t0 /
// fibers_per_layer. net_layout.txt дублирует это явно для Python-анализа,
// который читает object_id ИЗ КАДРА 0 (после разрывов object_id меняется --
// notes/fiber-model/fiber-problem-statement.md, разрыв связи -- но id бусины (индекс в FiberState) не
// меняется НИКОГДА, поэтому связка "id бусины из кадра 0 -> слой" остаётся
// верной все кадры).
// ---------------------------------------------------------------------
void write_net_file(const std::string& fib_path, const std::string& layout_path, double E, double nu, double d,
                     double rho, double eps_break, double l0, int n_beads, int fibers_per_layer,
                     double fiber_spacing, double layer_gap, int n_layers, double mu_ff, double mu_fs)
{
    std::ofstream f(fib_path);
    if (!f.is_open())
    {
        throw std::runtime_error("demo_fiber_layers: не удалось создать файл сети '" + fib_path + "'");
    }
    std::ofstream lay(layout_path);
    if (!lay.is_open())
    {
        throw std::runtime_error("demo_fiber_layers: не удалось создать файл разметки '" + layout_path + "'");
    }

    f << "# Демо-сеть: 4 слоя нитей x,y,x,y (experiments/2026-09-27-fiber-layers-demo)\n";
    f << "# Сгенерировано demo_fiber_layers.cpp -- см. README.md рядом за обоснованием параметров.\n";
    f << "MATERIAL 0 " << E << " " << nu << " " << d << " " << rho << " " << eps_break;
    if (mu_ff > 0.0 || mu_fs > 0.0)
    {
        f << " " << mu_ff;
        if (mu_fs > 0.0)
        {
            f << " " << mu_fs;
        }
    }
    f << "\n";

    lay << "# layer_index direction fiber_index_in_layer object_id_t0 offset z_layer\n";

    const int half = n_beads / 2; // n_beads нечётное -> бусины -half..0..half, есть центральная бусина
    const int half_fibers = (fibers_per_layer - 1) / 2;

    int fid = 0;
    for (int layer = 0; layer < n_layers; ++layer)
    {
        const bool along_x = (layer % 2 == 0); // 0,2 -- вдоль x; 1,3 -- вдоль y
        // z слоя: слои центрированы вокруг z=0, шаг между соседними -- ровно
        // layer_gap (касание без проникания на старте).
        const double z_layer = (static_cast<double>(n_layers - 1) / 2.0 - static_cast<double>(layer)) * layer_gap;

        for (int fl = 0; fl < fibers_per_layer; ++fl)
        {
            const double offset = static_cast<double>(fl - half_fibers) * fiber_spacing;

            f << "FIBER 0 " << n_beads << "\n";
            for (int k = -half; k <= half; ++k)
            {
                const double s = static_cast<double>(k) * l0;
                if (along_x)
                    f << s << " " << offset << " " << z_layer << "\n";
                else
                    f << offset << " " << s << " " << z_layer << "\n";
            }

            lay << layer << " " << (along_x ? "x" : "y") << " " << fl << " " << fid << " " << offset << " "
                << z_layer << "\n";
            ++fid;
        }
    }

    f.close();
    lay.close();
}

// ---------------------------------------------------------------------
// Снаряд -- сплошной куб SPH-частиц (алюминий, itype=body=4), решётка по
// образцу demo_fiber_net.cpp::build_projectile_cube() (узлы
// x=(i+0.5)*dx), без отбора по радиусу -- куб, а не шар.
// ---------------------------------------------------------------------
std::size_t build_projectile_cube(ParticleSystem& sys, int n_side, double dx, double rho0, double cx, double cy,
                                   double cz, double vz0)
{
    const int half = n_side / 2;
    const std::size_t ntot = static_cast<std::size_t>(n_side) * static_cast<std::size_t>(n_side) *
                              static_cast<std::size_t>(n_side);
    sys.resize(ntot);

    std::size_t k = 0;
    for (int ix = -half; ix < half; ++ix)
        for (int iy = -half; iy < half; ++iy)
            for (int iz = -half; iz < half; ++iz)
            {
                sys.inzializ_index(k);
                sys.x[0][k] = cx + (static_cast<double>(ix) + 0.5) * dx;
                sys.x[1][k] = cy + (static_cast<double>(iy) + 0.5) * dx;
                sys.x[2][k] = cz + (static_cast<double>(iz) + 0.5) * dx;
                sys.vx[2][k] = vz0;
                sys.itype[k] = 4; // алюминий (materials.hpp)
                sys.body[k]  = 4; // body[]=itype[] -- notes/fiber-model/fiber-plan.md §3.1
                sys.rho[k]   = rho0;
                sys.rho_0[k] = rho0;
                sys.hsml[k]  = 2.0 * dx;
                sys.mass[k]  = rho0 * dx * dx * dx;
                ++k;
            }
    return ntot;
}

} // namespace

int main(int argc, char** argv)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // UTF-8 в консоли Windows (как main.cpp/time_integration.cpp)
#endif

    std::string output_base_arg = "experiments/2026-09-27-fiber-layers-demo/data";
    std::string scenario = "pierce";
    if (argc > 1) output_base_arg = argv[1];
    if (argc > 2) scenario = argv[2];
    if (scenario != "pierce" && scenario != "stop")
    {
        std::cerr << "demo_fiber_layers: сценарий должен быть 'pierce' или 'stop', получено '" << scenario << "'\n";
        return 1;
    }

    const fs::path out_dir = fs::absolute(fs::path(output_base_arg) / scenario);
    fs::create_directories(out_dir);
    const fs::path net_path = out_dir / "net.fib";
    const fs::path layout_path = out_dir / "net_layout.txt";

    // Переменные окружения:
    // - SPH_DEMO_MU_FF: коэффициент трения нить-нить (если не задана, 0)
    // - SPH_DEMO_MU_FS: коэффициент трения нить-Solid (если не задана, 0)
    // - SPH_DEMO_BUDGET_S: временной бюджет в секундах для RunGuard (если не задана, 900)
    double mu_ff = 0.0;
    const char* mu_ff_env = std::getenv("SPH_DEMO_MU_FF");
    if (mu_ff_env != nullptr)
    {
        try
        {
            mu_ff = std::stod(mu_ff_env);
        }
        catch (const std::exception&)
        {
            std::fprintf(stderr, "demo_fiber_layers: ошибка разбора SPH_DEMO_MU_FF='%s'\n", mu_ff_env);
            return 1;
        }
        if (mu_ff < 0.0)
        {
            std::fprintf(stderr, "demo_fiber_layers: SPH_DEMO_MU_FF должен быть >= 0, получено %.6e\n", mu_ff);
            return 1;
        }
    }

    double mu_fs = 0.0;
    const char* mu_fs_env = std::getenv("SPH_DEMO_MU_FS");
    if (mu_fs_env != nullptr)
    {
        try
        {
            mu_fs = std::stod(mu_fs_env);
        }
        catch (const std::exception&)
        {
            std::fprintf(stderr, "demo_fiber_layers: ошибка разбора SPH_DEMO_MU_FS='%s'\n", mu_fs_env);
            return 1;
        }
        if (mu_fs < 0.0)
        {
            std::fprintf(stderr, "demo_fiber_layers: SPH_DEMO_MU_FS должен быть >= 0, получено %.6e\n", mu_fs);
            return 1;
        }
    }

    double budget_s = 900.0;
    const char* budget_env = std::getenv("SPH_DEMO_BUDGET_S");
    if (budget_env != nullptr)
    {
        try
        {
            budget_s = std::stod(budget_env);
        }
        catch (const std::exception&)
        {
            std::fprintf(stderr, "demo_fiber_layers: ошибка разбора SPH_DEMO_BUDGET_S='%s'\n", budget_env);
            return 1;
        }
        if (budget_s <= 0.0)
        {
            std::fprintf(stderr, "demo_fiber_layers: SPH_DEMO_BUDGET_S должен быть > 0, получено %.6e\n", budget_s);
            return 1;
        }
    }

    // =====================================================================
    // Параметры сети -- ОДИНАКОВЫЕ для обоих сценариев (README.md -- откуда
    // взялись числа: оценка энергоёмкости сети против кинетической энергии
    // снаряда).
    // =====================================================================
    const double E_f = 5.0e7, nu_f = 0.3, d_f = 1.5e-3, rho_f = 1200.0;
    const double eps_break = 0.15; // подобрано так, чтобы порог pierce/stop лежал в единицах-десятках м/с
    const double l0 = 2.0e-3;      // шаг бусин вдоль нити
    const int n_beads = 17;        // 16 сегментов * 2 мм = 32 мм на нить (~4x снаряда, 7.8 мм)

    const int n_layers = 4;
    const int fibers_per_layer = 5;         // ориентир 4-6
    const double fiber_spacing = 3.0e-3;     // > d_f (1.5 мм) -- соседи в слое не касаются на старте
    const double layer_gap = d_f;            // соседние слои касаются без проникания (как в demo_fiber_net)

    const double anchor_mass_factor = 50.0; // псевдо-анкеры -- как в demo_fiber_net

    write_net_file(net_path.string(), layout_path.string(), E_f, nu_f, d_f, rho_f, eps_break, l0, n_beads,
                   fibers_per_layer, fiber_spacing, layer_gap, n_layers, mu_ff, mu_fs);

    // --- Снаряд: алюминий, itype=body=4, МЕНЕЕ МАССИВНЫЙ, чем в
    // demo_fiber_net (куб 8^3, dx=1.2мм, ~2.4 г): куб 6^3=216 частиц (>=100),
    // dx=1.3мм (крупнее, чем в demo_fiber_net -- смягчает ужесточение CFL от
    // меньшего размера), масса ~1.28 г (~54% прежней).
    const int n_side = 6;
    const double dx_p = 1.3e-3;
    const double rho_s = 2700.0; // плотность алюминия
    const double half_extent = (n_side / 2) * dx_p;
    const double gap = 4.0e-3; // зазор снаряд -> верхний слой на старте

    // Верхний слой (layer=0, направление x) на z = +1.5*d_f (см. write_net_file).
    const double z_top_layer = (static_cast<double>(n_layers - 1) / 2.0) * layer_gap;
    const double cz = z_top_layer + half_extent + gap;

    // --- Скорость и длительность сценария -- ЕДИНСТВЕННОЕ различие между
    // pierce и stop (net и снаряд идентичны). Числа -- README.md.
    double v0 = 0.0, t_total = 0.0, save_t = 0.0;
    if (scenario == "pierce")
    {
        v0 = 35.0;
        t_total = 6.0e-4;
        save_t = 8.0e-6; // ~75 кадров
    }
    else // stop
    {
        // v0=15 м/с ещё пробивал сеть (остаточная vz~-6.1 м/с, чуть ниже
        // нижнего слоя); v0=10 м/с при t_total=1.0e-3 выходил на плато
        // vz~-0.7..-0.8 м/с (почти остановка, ползучее продавливание только
        // верхнего слоя, z снаряда не опускался ниже -2.29 мм -- нижний
        // слой на z=-2.25 мм даже не задет) -- калибровочные прогоны, см.
        // README.md. Взято v0=10 м/с с увеличенным t_total, чтобы плато vz
        // успело разрешиться в остановку/отскок, а не обрывалось на полпути.
        v0 = 10.0;
        t_total = 1.8e-3;
        save_t = 2.2e-5; // ~82 кадра
    }
    // Необязательные argv[3]/argv[4]/argv[5] -- переопределение
    // t_total/save_t/v0 для калибровочных прогонов, не влияет на геометрию
    // сети/материалы.
    if (argc > 3) t_total = std::stod(argv[3]);
    if (argc > 4) save_t = std::stod(argv[4]);
    if (argc > 5) v0 = std::stod(argv[5]);

    fiber::Model model;
    const std::size_t ntotal = build_projectile_cube(model.solid, n_side, dx_p, rho_s, 0.0, 0.0, cz, -v0);

    fiber::load_fiber_file(net_path.string(), model);

    // Псевдо-анкеры -- концевые бусины: prev==-1 || next==-1
    // (fiber/fiber_system.hpp, FiberTopology).
    std::size_t n_anchors = 0;
    for (std::size_t i = 0; i < model.fiber.size(); ++i)
    {
        if (model.topology.prev[i] == -1 || model.topology.next[i] == -1)
        {
            model.fiber.mass[i] *= anchor_mass_factor;
            ++n_anchors;
        }
    }

    // =====================================================================
    // SimulationConfig -- значения по умолчанию там, где нет причин их менять.
    // =====================================================================
    SimulationConfig cfg;
    cfg.apply_defaults();
    cfg.fiber_file = net_path.string(); // только для протокола -- парсинг уже сделан выше
    cfg.t_total = t_total;
    cfg.save_t = save_t;
    cfg.contact_force_type = 0; // единственное Solid-тело -- Solid-Solid-Contact не нужен
    // SPH_DEMO_CFL -- число Куранта вместо значения по умолчанию (проверка
    // сходимости баланса энергии по шагу по времени).
    if (const char* cfl_env = std::getenv("SPH_DEMO_CFL"))
    {
        cfg.CFL = std::stod(cfl_env);
        std::printf("demo_fiber_layers: CFL=%.4f\n", cfg.CFL);
    }
    cfg.moni_particle = 0;
    cfg.print_step = 2000; // много шагов ожидается -- реже пишем прогресс

    cfg.maxn = static_cast<std::size_t>(static_cast<double>(ntotal) * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    {
        const std::size_t prev = model.solid.size();
        model.solid.resize(cfg.maxn);
        for (std::size_t i = prev; i < model.solid.size(); ++i) model.solid.inzializ_index(i);
    }
    for (std::size_t i = 0; i < ntotal; ++i) urs(cfg.mats, model.solid, i);

    const double mass_total = rho_s * dx_p * dx_p * dx_p * static_cast<double>(ntotal);
    std::printf(
        "demo_fiber_layers: scenario=%s v0=%.3f m/s t_total=%.4e save_t=%.4e\n", scenario.c_str(), v0, t_total,
        save_t);
    std::printf(
        "demo_fiber_layers: snar'ad ntotal=%zu (%d^3, dx=%.4e), massa=%.6e kg, niti busin=%zu (%zu ankerov), "
        "sv'azey=%zu, sloev=%d x %d niti\n",
        ntotal, n_side, dx_p, mass_total, model.fiber.size(), n_anchors, model.topology.bonds.size(), n_layers,
        fibers_per_layer);
    std::printf("demo_fiber_layers: mu_ff=%.6e, mu_fs=%.6e, budget=%.3f s\n", mu_ff, mu_fs, budget_s);
    std::printf("demo_fiber_layers: output_dir=%s\n", out_dir.string().c_str());

    // Смена рабочего каталога -- output()/fiber_output()/log.txt/sph_gp.csv
    // жёстко пишут в CWD.
    fs::current_path(out_dir);

    double dt = compute_cfl_timestep(cfg, model.solid, ntotal);
    double time = 0.0;
    RunGuard guard(budget_s, 0.0, 1);

    const std::size_t bonds_before = model.topology.bonds.size();
    time_integration(cfg, model, ntotal, dt, time, dx_p, &guard);

    std::printf("demo_fiber_layers: t=%.6e iz %.6e, shagov=%zu, sv'azey ostalos'=%zu iz %zu\n", time, cfg.t_total,
                guard.steps(), model.topology.bonds.size(), bonds_before);
    if (guard.triggered())
    {
        std::printf("%s\n", guard.message().c_str());
        if (guard.failed()) return 1;
    }
    return 0;
}
