// experiments/2026-09-26-fiber-net-demo/demo_fiber_net.cpp
//
// Демонстрационный расчёт: снаряд из SPH-частиц (алюминий) пробивает
// крестообразную сеть из двух нитей. Постановка -- fiber/DEMO_VISUALIZATION_TASK.md.
//
// В отличие от e2e-тестов (tests/test_e2e_fiber_contact.cpp,
// tests/test_e2e_fiber_solid.cpp), которые интегрируют локальным
// PredictorCorrector-драйвером мимо Solid-подсистемы, это НЕ тест: сцена
// прогоняется через ПРОДОВЫЙ цикл time_integration(cfg, fiber::Model&, ...)
// (include/time_integration.hpp) -- тот же путь, что main.cpp, только с
// непустой fiber::Model. Это же служит дополнительной проверкой продового
// пути на сцене нового размера (~500 Solid-частиц + сеть), которая раньше
// не прогонялась (постановка, раздел 4).
//
// Сеть задаётся файлом формата fiber_input.cpp (src/fiber_input.cpp,
// комментарий в шапке -- полное описание формата) и загружается
// fiber::load_fiber_file() -- продовый путь загрузки (SimulationConfig::
// fiber_file, см. постановку, раздел 3, способ 2 -- предпочтительный).
//
// Вывод кадров: time_integration сама вызывает output() (src/output.cpp,
// Solid) и fiber::fiber_output() (src/fiber_output.cpp, нити) с одинаковой
// периодичностью cfg.save_t (src/time_integration.cpp -- вызов добавлен
// туда же, рядом с output(), см. постановку, раздел 4). Драйвер здесь этот
// вызов не дублирует и не оборачивает -- только готовит сцену и передаёт её
// в общий цикл.
//
// Запуск: demo_fiber_net [output_dir]
//   output_dir -- каталог для net.fib и кадров (sph_output_*.txt,
//   fiber_output_*.txt, log.txt, sph_gp.csv); создаётся, если не существует.
//   По умолчанию -- experiments/2026-09-26-fiber-net-demo/data/run ОТНОСИТЕЛЬНО
//   текущего рабочего каталога запуска. time_integration/output()/
//   fiber_output() жёстко пишут в текущий рабочий каталог процесса, поэтому
//   драйвер переходит в output_dir (fs::current_path) ПОСЛЕ того, как все
//   относительные пути (в т.ч. output_dir сам по себе) резолвятся в
//   абсолютные.
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

static_assert(dim == 3, "demo_fiber_net требует сборки с SPH_DIM=3 (нити -- только в 3D, fiber/SPH.md)");

namespace fs = std::filesystem;

namespace
{

// ---------------------------------------------------------------------
// Сеть -- две прямые нити крестом: P вдоль оси x на высоте z=zP, Q вдоль
// оси y на высоте z=zQ, обе центрированы в (0,0) -- траектория снаряда
// проходит через точку их проекционного пересечения (постановка, раздел 2).
// Пишется в формате fiber_input.cpp (MATERIAL/FIBER + точки ломаной),
// читается фактическим продовым парсером fiber::load_fiber_file().
// ---------------------------------------------------------------------
void write_net_file(const std::string& path, double E, double nu, double d, double rho, double eps_break, double l0,
                     int n_beads, double zP, double zQ, double mu_ff, double mu_fs)
{
    std::ofstream f(path);
    if (!f.is_open())
    {
        throw std::runtime_error("demo_fiber_net: не удалось создать файл сети '" + path + "'");
    }

    f << "# Демо-сеть: крест из двух нитей (fiber/DEMO_VISUALIZATION_TASK.md)\n";
    f << "# Сгенерировано demo_fiber_net.cpp -- см. README.md рядом за обоснованием параметров.\n";
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

    const int half = n_beads / 2; // n_beads нечётное -> бусины -half..0..half, есть центральная

    f << "FIBER 0 " << n_beads << "\n";
    for (int k = -half; k <= half; ++k)
    {
        f << (static_cast<double>(k) * l0) << " 0.0 " << zP << "\n";
    }

    f << "FIBER 0 " << n_beads << "\n";
    for (int k = -half; k <= half; ++k)
    {
        f << "0.0 " << (static_cast<double>(k) * l0) << " " << zQ << "\n";
    }

    f.close();
}

// ---------------------------------------------------------------------
// Снаряд -- сплошной куб SPH-частиц (алюминий, itype=body=4), решётка по
// образцу tests/bc2_mini3d_geometry.hpp::build_geometry() (узлы
// x=(i+0.5)*dx), без отбора по радиусу -- куб, а не шар (постановка,
// раздел 3). Обязательно заполняет sys.body[k] (=itype[k]) -- без этого
// Solid-Solid-Contact/фильтры по телу не работают (fiber/PLAN.md §3.1).
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
                sys.body[k]  = 4; // body[]=itype[] -- fiber/PLAN.md §3.1
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

    // Каталог вывода -- аргумент командной строки или значение по умолчанию.
    // Резолвится в АБСОЛЮТНЫЙ путь ДО смены рабочего каталога (иначе
    // относительный output_dir потерял бы смысл после chdir ниже).
    std::string output_dir_arg = "experiments/2026-09-26-fiber-net-demo/data/run";
    if (argc > 1) output_dir_arg = argv[1];

    const fs::path out_dir = fs::absolute(output_dir_arg);
    fs::create_directories(out_dir);
    const fs::path net_path = out_dir / "net.fib";

    // Переменные окружения:
    // - SPH_DEMO_MU_FF: коэффициент трения нить-нить (если не задана, 0)
    // - SPH_DEMO_MU_FS: коэффициент трения нить-Solid (если не задана, 0)
    // - SPH_DEMO_BUDGET_S: временной бюджет в секундах для RunGuard (если не задана, 600)
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
            std::fprintf(stderr, "demo_fiber_net: ошибка разбора SPH_DEMO_MU_FF='%s'\n", mu_ff_env);
            return 1;
        }
        if (mu_ff < 0.0)
        {
            std::fprintf(stderr, "demo_fiber_net: SPH_DEMO_MU_FF должен быть >= 0, получено %.6e\n", mu_ff);
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
            std::fprintf(stderr, "demo_fiber_net: ошибка разбора SPH_DEMO_MU_FS='%s'\n", mu_fs_env);
            return 1;
        }
        if (mu_fs < 0.0)
        {
            std::fprintf(stderr, "demo_fiber_net: SPH_DEMO_MU_FS должен быть >= 0, получено %.6e\n", mu_fs);
            return 1;
        }
    }

    double budget_s = 600.0;
    const char* budget_env = std::getenv("SPH_DEMO_BUDGET_S");
    if (budget_env != nullptr)
    {
        try
        {
            budget_s = std::stod(budget_env);
        }
        catch (const std::exception&)
        {
            std::fprintf(stderr, "demo_fiber_net: ошибка разбора SPH_DEMO_BUDGET_S='%s'\n", budget_env);
            return 1;
        }
        if (budget_s <= 0.0)
        {
            std::fprintf(stderr, "demo_fiber_net: SPH_DEMO_BUDGET_S должен быть > 0, получено %.6e\n", budget_s);
            return 1;
        }
    }

    // =====================================================================
    // Параметры сцены (fiber/DEMO_VISUALIZATION_TASK.md, разделы 1-2).
    // Итоговые числа и обоснование выбора -- см. README.md рядом с этим файлом.
    // =====================================================================

    // --- Нить: материал того же порядка, что test_e2e_fiber_contact.cpp ---
    const double E_f = 2.0e7, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    const double eps_break = 0.05; // подобран экспериментально -- см. README
    const double l0 = 2.0e-3;      // шаг бусин вдоль нити
    const int n_beads = 21;        // 20 сегментов * 2 мм = 40 мм на нить (>>снаряда)

    // D0 = zP - zQ = d_f: нити касаются без проникания, контактная сила на
    // старте нулевая; Fiber-Fiber-Contact включается только от удара снаряда.
    const double zP = 0.5e-3;  // нить P (вдоль x) чуть выше плоскости z=0
    const double zQ = -0.5e-3; // нить Q (вдоль y) чуть ниже

    // Псевдо-анкеры: утяжеление концевых бусин, чтобы сеть не улетала целиком
    // от удара (постановка, раздел 2, "открытый вопрос реализации").
    const double anchor_mass_factor = 50.0;

    write_net_file(net_path.string(), E_f, nu_f, d_f, rho_f, eps_break, l0, n_beads, zP, zQ, mu_ff, mu_fs);

    // --- Снаряд: алюминий, itype=body=4, как build_plate() в
    // tests/test_e2e_fiber_solid.cpp -----------------------------------------
    const int n_side = 8;        // 8^3 = 512 частиц (в диапазоне 300-800)
    const double dx_p = 1.2e-3;  // шаг решётки снаряда
    const double rho_s = 2700.0; // плотность алюминия (как в тестах Ф6)
    const double v0 = 30.0;      // м/с, вдоль -z (десятки м/с)
    const double half_extent = (n_side / 2) * dx_p;
    const double gap = 4.0e-3; // зазор снаряд -> верхняя нить на старте
    const double cz = zP + half_extent + gap;

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
    // SimulationConfig -- значения по умолчанию там, где нет причин их менять
    // (постановка, раздел 3).
    // =====================================================================
    SimulationConfig cfg;
    cfg.apply_defaults();
    cfg.fiber_file = net_path.string(); // только для протокола -- парсинг уже сделан выше
    cfg.t_total = 5.0e-4;               // ~3.7x времени подлёта снаряда до сети
    cfg.save_t  = 8.0e-6;               // ~62 кадра за весь прогон (30-100, постановка §3)
    cfg.contact_force_type = 0;         // единственное Solid-тело -- Solid-Solid-Contact не нужен
    cfg.moni_particle = 0;

    cfg.maxn = static_cast<std::size_t>(static_cast<double>(ntotal) * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    {
        // Хвост [ntotal, maxn) -- под возможные виртуальные частицы
        // (main.cpp: "запас нужен под возможные виртуальные частицы").
        const std::size_t prev = model.solid.size();
        model.solid.resize(cfg.maxn);
        for (std::size_t i = prev; i < model.solid.size(); ++i) model.solid.inzializ_index(i);
    }
    for (std::size_t i = 0; i < ntotal; ++i) urs(cfg.mats, model.solid, i);

    std::printf("demo_fiber_net: snar'ad ntotal=%zu (%d^3), niti busin=%zu (%zu ankerov), sv'azey=%zu\n", ntotal,
                n_side, model.fiber.size(), n_anchors, model.topology.bonds.size());
    std::printf("demo_fiber_net: mu_ff=%.6e, mu_fs=%.6e, budget=%.3f s\n", mu_ff, mu_fs, budget_s);
    std::printf("demo_fiber_net: output_dir=%s\n", out_dir.string().c_str());

    // Смена рабочего каталога -- output()/fiber_output()/log.txt/sph_gp.csv
    // жёстко пишут в CWD (постановка, раздел 5).
    fs::current_path(out_dir);

    double dt = compute_cfl_timestep(cfg, model.solid, ntotal);
    double time = 0.0;
    // Сторож -- только бюджет по времени, без слежения за застоем:
    // число шагов оценено заранее (см. README), ложного срабатывания не
    // ожидается; если оценка неверна -- лучше явный, читаемый обрыв по
    // wall-clock, чем зависший процесс.
    RunGuard guard(budget_s, 0.0, 1);

    const std::size_t bonds_before = model.topology.bonds.size();
    time_integration(cfg, model, ntotal, dt, time, dx_p, &guard);

    std::printf("demo_fiber_net: t=%.6e iz %.6e, shagov=%zu, sv'azey ostalos'=%zu iz %zu\n", time, cfg.t_total,
                guard.steps(), model.topology.bonds.size(), bonds_before);
    if (guard.triggered())
    {
        std::printf("%s\n", guard.message().c_str());
        if (guard.failed()) return 1;
    }
    return 0;
}
