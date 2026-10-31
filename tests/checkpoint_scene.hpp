#pragma once
// tests/checkpoint_scene.hpp
// Небольшая переиспользуемая сцена для тестов чекпоинта/рестарта (ветка
// fiber_checkpoint_restart). Код скопирован из tests/test_k_e2e.cpp (make_deck,
// write_ic, загрузка через input()), сам test_k_e2e.cpp не менялся.
//
// СОСТАВ СЦЕНЫ.
//   Solid: пластина 4x4x2 SPH-узлов, 32 частицы, шаг 1 мм (itype 10, rho 2510),
//          импорт из .k-колоды (k_file = deck.k, k_parts_file = deck.parts).
//   Fiber: два пучка по 4 бусины (шаг 1 мм, L = 3 мм) = 8 бусин и 6 связей,
//          из .fib (fiber_file = scene.fib). Материал упругопластичный Twaron:
//          PLASTIC + FAILURE + DAMAGE (образец experiments/2026-10-26-guo2019-
//          damage-summary/twaron_s25_dmg0.00415.fib), eps_break -- аргумент.
//   Нити зажаты CLAMP_BOX у левого конца (x = 0.050 м), правый конец разгоняется
//          линейным профилем скорости vx = kPullSpeed * (x - 0.050) / L, так что
//          связи растягиваются с первого шага. При малом eps_break связи рвутся
//          (проверка разрыва и повреждённости), при большом -- нет.
//   Пластина и нити разнесены (пластина при x, y in [0, 3] мм, нити при
//          x in [50, 53] мм), поэтому контакт Solid-Fiber не участвует.
//
// ЧТО ЧИТАЕТ И ПИШЕТ. input() берёт Initial_Conditions.txt из текущего каталога,
// а цикл по времени пишет кадры и log.txt туда же. Поэтому build_checkpoint_scene
// создаёт workdir, переходит в него (cwd не восстанавливает -- это делает вызывающая
// сторона, например через e2e::ScopedWorkDir) и сразу готовит состояние к
// time_integration: maxn/max_interaction, расширенный model.solid, dt по Куранту.
// Шаги выполняет вызывающая сторона: time_integration(sc.cfg, sc.model,
// sc.ntotal, sc.dt, time, sc.delta_x, &guard).

#include "fiber/fiber_input.hpp"
#include "fiber/model.hpp"
#include "input.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

// Готовая к запуску сцена. dt и delta_x передаются в time_integration по ссылке
// и по значению соответственно; dt меняется во время прогона (это ожидаемо).
struct CkptScene
{
    SimulationConfig cfg;
    fiber::Model model;
    std::size_t ntotal = 0;
    double dt = 0.0;
    double delta_x = 0.0;
};

namespace ckpt_scene_detail
{

// Скорость вытягивания правого конца нитей, м/с. dt ~ 1.4e-8 с, поэтому за 40
// шагов (~0.5 мкс) при 60 м/с связи растягиваются примерно на 1 %.
constexpr double kPullSpeed = 60.0;
// Положение зажима (левый конец) и длина нити вдоль x, м.
constexpr double kRootX = 0.050;
constexpr double kFiberLen = 0.003;

// Таблица частей для импорта .k (как в test_k_e2e.cpp, kTable).
inline constexpr const char* kParts =
    "SOLID_PART 1 10 2510 hsml=0.0013\n"
    "MATERIAL 0 2.0e9 0.35 1.0e-4 1140 0.02\n";

inline std::string read_all(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("не открыть " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

inline void write_file(const std::filesystem::path& p, const std::string& text)
{
    std::ofstream out(p, std::ios::binary);
    if (!out) throw std::runtime_error("не создать " + p.string());
    out << text;
}

// Сетка 4x4x2 SPH-узлов (NID 1..32, PID 1), шаг 1 мм. Нити в .k не входят:
// они берутся из .fib (fiber_file).
inline std::string make_solid_deck()
{
    std::ostringstream k;
    k << "*KEYWORD\n*PART\nplate\n1\n*NODE\n";
    int nid = 1;
    for (int iz = 0; iz < 2; ++iz)
        for (int iy = 0; iy < 4; ++iy)
            for (int ix = 0; ix < 4; ++ix)
                k << nid++ << "," << ix * 1.0e-3 << "," << iy * 1.0e-3 << "," << iz * 1.0e-3 << "\n";
    k << "*ELEMENT_SPH\n";
    for (int n = 1; n <= 32; ++n) k << n << ",1,2.51e-6\n";
    k << "*END\n";
    return k.str();
}

// Файл нитей: материал 0 (Twaron, E = 80 ГПа как в образце с повреждённостью),
// два пучка по 4 бусины вдоль x (y = 0 и y = 10 мм), зажим у левого конца.
inline std::string make_fiber_file(double eps_break)
{
    std::ostringstream os;
    os << std::setprecision(12);
    os << "# Сцена чекпоинт-теста (tests/checkpoint_scene.hpp): Twaron, пластичность + повреждённость\n";
    os << "MATERIAL 0 80.0e9 0.35 1.0e-4 1140 " << eps_break << "\n";
    os << "PLASTIC 0 SY=2.04e9 H=100e9\n";
    os << "FAILURE 0 SIGMA_U=2.5e9\n";
    os << "DAMAGE 0 DEPS_F=0.00415\n";
    for (int f = 0; f < 2; ++f)
    {
        const double y = f * 0.01;
        os << "FIBER 0 4\n";
        for (int j = 0; j < 4; ++j) os << (kRootX + j * 1.0e-3) << " " << y << " 0.0\n";
        os << "CLAMP_BOX " << kRootX - 0.5e-3 << " " << y - 0.5e-3 << " -0.5e-3 "
           << kRootX + 0.5e-3 << " " << y + 0.5e-3 << " 0.5e-3\n";
    }
    return os.str();
}

// Initial_Conditions.txt: общая часть из tests/data + ключи импорта.
inline void write_ic(const std::string& extra)
{
    std::string base = read_all(std::string(SPH_TEST_DATA_DIR) + "/guo2019_initial_conditions.txt");
    write_file("Initial_Conditions.txt", base + "\n" + extra + "\n");
}

} // namespace ckpt_scene_detail

// Собирает сцену в каталоге workdir (создаётся, cwd переводится туда).
//   eps_break -- критическое удлинение связи нити (MATERIAL), малое значение
//                форсирует разрыв; ~0.5 даёт разрывов не ожидать.
//   ode_solver -- 0 = предиктор-корректор, 1 = Эйлер (cfg.ode_solver).
// Бросает std::runtime_error при сбое загрузки и std::invalid_argument при
// неверном ode_solver.
inline CkptScene build_checkpoint_scene(const std::string& workdir, double eps_break, int ode_solver)
{
    namespace fs = std::filesystem;
    using namespace ckpt_scene_detail;

    if (ode_solver != 0 && ode_solver != 1)
        throw std::invalid_argument("build_checkpoint_scene: ode_solver вне 0..1");

    fs::create_directories(workdir);
    fs::current_path(workdir);
    write_file("deck.k", make_solid_deck());
    write_file("deck.parts", kParts);
    write_file("scene.fib", make_fiber_file(eps_break));
    write_ic("k_file = deck.k\nk_parts_file = deck.parts\nfiber_file = scene.fib\n");

    CkptScene sc;
    sc.cfg.apply_defaults();
    input(sc.cfg, sc.model, sc.ntotal, sc.delta_x);
    sc.cfg.ode_solver = ode_solver;

    if (sc.ntotal != 32 || sc.model.fiber.size() != 8 || sc.model.topology.bonds.size() != 6)
        throw std::runtime_error("build_checkpoint_scene: состав сцены не совпал с ожидаемым");
    if (sc.model.clamps.empty())
        throw std::runtime_error("build_checkpoint_scene: зажимы нитей не загружены");

    // Скорости правых бусин: линейный профиль от зажима. Бусины внутри CLAMP_BOX
    // зажаты проекцией и своих скоростей не получают.
    fiber::FiberState& fb = sc.model.fiber;
    for (std::size_t k = 0; k < fb.size(); ++k)
    {
        const double x = fb.x[0][k];
        if (x > kRootX + 0.5e-3)
            fb.vx[0][k] = kPullSpeed * (x - kRootX) / kFiberLen;
    }

    // Как run_steps в test_k_e2e.cpp: расширение model.solid до maxn и dt по Куранту.
    sc.cfg.maxn = static_cast<std::size_t>(sc.ntotal * 1.4);
    sc.cfg.max_interaction = 100 * sc.cfg.maxn;
    const std::size_t prev = sc.model.solid.size();
    if (prev < sc.cfg.maxn)
    {
        sc.model.solid.resize(sc.cfg.maxn);
        for (std::size_t i = prev; i < sc.model.solid.size(); ++i) sc.model.solid.inzializ_index(i);
    }
    sc.cfg.t_total = 1.0;
    sc.dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal);
    return sc;
}
