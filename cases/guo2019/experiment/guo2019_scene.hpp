#pragma once
// guo2019_scene.hpp
// ЭКСПЕРИМЕНТАЛЬНАЯ КОПИЯ (cases/README.md, «Роли постановок»): ответвлена от
// cases/guo2019/regression/ в a9c3964; меняется свободно под расчёты
// (guo2019_shot, fiber_t0_dump, fabric_checks). Тесты её не используют, кроме
// проверок самих этих утилит. Отличия: diff -r cases/guo2019/regression cases/guo2019/experiment.
// Сцена «шар + ткань из нитей + Al-пластина» для опытов Guo et al. 2019 (ТЗ Guo
// 2019 §3, §4.1, §5). Одна функция build_guo_scene() собирает fiber::Model и
// SimulationConfig по имени конфигурации (A, 1T, 1T_A, A_1T, 2T, 4T) и уровню
// разрешения (smoke, base, fine). От ctest не зависит: бросает исключения, а не
// ведёт счётчик failures. Все числа -- cases/guo2019/experiment/guo2019_case.hpp.
//
// ПОРЯДОК СБОРКИ Model (ТЗ Guo 2019 §4.1):
//   1) частицы Solid (шар, затем пластина) -> model.solid;
//   2) fiber::load_fiber_file (НЕ трогает model.solid; заполняет model.fiber,
//      model.topology, model.materials и model.clamps.fiber_mask из директив
//      CLAMP_BOX, которые записал генератор; при отсутствии CLAMP_BOX очищает
//      fiber_mask);
//   3) model.clamps.solid_mask размера ровно ntotal (а не model.solid.size():
//      при maxn = 1.4 ntotal массивы длиннее, и validate_clamps бросил бы).
//
// ИНДЕКСЫ. Solid: шар [0, n_ball), пластина [n_ball, ntotal); внутри пластины
// индекс n_ball + (ix*ny + iy)*nz + iz, iz = 0 -- ТЫЛЬНЫЙ (нижний по z) слой.
// Бусины ткани: подряд по нитям в порядке генератора (сначала основа вдоль x,
// затем уток вдоль y); у каждой нити диапазон [bead_begin, bead_end).
// ПАКЕТ СЛОЁВ (2T, 4T; ТЗ 2T/4T §3.2, §4.2): блоки FIBER всех слоёв идут в одном
// .fib по слоям (слой 1 -- лицевой, по нему бьёт шар), внутри слоя в порядке
// генератора. object_id нити -- глобальный: (k - 1)(nx + ny) + локальный id.
// Слои одинаковые, без сдвига в плоскости, шаг по z Delta_ply = 2A + d + kEpsPly.
//
// НАПРАВЛЕНИЯ НИТЕЙ. В генераторе «основа» (direction 'x', nx нитей) идёт вдоль
// x, «уток» (direction 'y', ny нитей) -- вдоль y. Зажим clamp = Y пишет два
// CLAMP_BOX по оси y, и зажатой оказывается нить с direction == 'y': это
// вертикальные (зажатые) нити статьи, идущие вдоль y (ТЗ Guo 2019 §3.1).
// Горизонтальные ('x') свободны, их концы у x = +-(60 мм) не закреплены.
//
// hsml ЧАСТИЦ SOLID. Шар: 2 dx_s. Пластина: 2 dx_p, dx_p = t/n_z (номинальный
// шаг, общий для обеих осей; фактические шаги по x и y отличаются от него
// малой долей, см. dx_p_axis). Это начальное значение, как у bc2/bc6
// (hsml = 2 dx): при sle = 1 цикл сам пересчитывает h = 1.3 (m/rho)^(1/3).
// Масса частицы пластины rho dx_x dx_y dx_z, контактный радиус
// solid_contact_radius = 0.5 (dx_x dx_y dx_z)^(1/3), поэтому при неравных шагах
// он чуть отличается от dx_z/2 (это влияет на зазор ткань--пластина в 1T/A).

#include "e2e_harness.hpp"
#include "guo2019_case.hpp"
#include "lattice.hpp"

#include "fiber/fiber_input.hpp"
#include "fiber/model.hpp"
#include "fiber/plain_weave.hpp"
#include "fiber/plain_weave_quilt.hpp"
#include "fiber/plain_weave_random_cut.hpp"
#include "input.hpp"
#include "int_force.hpp"
#include "interactions/clamp.hpp"
#include "interactions/solid_contact.hpp"
#include "materials.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef SPH_CASES_DIR
// Задаётся из CMake (функция sph_case). Фолбэк -- запуск из корня репозитория.
#define SPH_CASES_DIR "cases"
#endif

static_assert(dim == 3, "guo2019_scene.hpp: сцена Guo 2019 требует SPH_DIM=3");

namespace guo2019
{

// Раскладка одной нити в общих массивах бусин (диапазон полуоткрытый).
struct YarnLayout
{
    int ply = 1;                    // номер слоя пакета, 1 -- лицевой
    char direction = 'x';           // 'x' горизонтальная (свободная), 'y' вертикальная (зажатая)
    int index = 0;                  // номер нити внутри своего направления
    int object_id = 0;              // глобальный номер нити (порядковый блок FIBER в .fib); при одном слое равен id генератора
    std::size_t bead_begin = 0;     // первая бусина нити
    std::size_t bead_end = 0;       // за последней бусиной
    std::size_t n_clamp_lo = 0;     // зажатых бусин у начала нити
    std::size_t n_clamp_hi = 0;     // зажатых бусин у конца нити
    bool clamped = false;
};

// Раскладка слоя пакета в общих массивах (диапазоны полуоткрытые).
struct PlyLayout
{
    std::size_t bead_begin = 0, bead_end = 0; // бусины слоя в model.fiber
    std::size_t yarn_begin = 0, yarn_end = 0; // нити слоя в GuoScene::yarns
    int object_id_begin = 0;                  // object_id первой нити слоя
};

// Ожидаемые числа слоя (из PlainWeave, ТЗ 2T/4T §3.2) и его поверхностная плотность.
struct PlyExpected
{
    std::size_t beads = 0, bonds = 0, clamped_beads = 0;
    double areal_density_window = 0.0;  // кг/м^2, связи слоя целиком в окне
    double areal_density_summary = 0.0; // plain_weave_summary слоя
    std::size_t bonds_in_window = 0;
};

struct GuoScene
{
    // --- основное ---------------------------------------------------------
    Config config = Config::A;
    Level level = Level::Smoke;
    LevelSpec spec{};
    double V0 = 0.0;                // м/с, шар летит вдоль -z
    double mu = 0.0;                // mu_ff = mu_fs
    bool body_isolation = true;
    bool has_fabric = false, has_plate = false;
    SimulationConfig cfg;
    fiber::Model model;
    fiber::PlainWeave W;            // ткань; при пакете -- слой 1 (лицевой); пусто, если has_fabric == false
    int n_plies = 0;                // число слоёв ткани (fabric_plies(config, plies))
    std::vector<fiber::PlainWeave> plies; // слои, слой 1 первым
    std::vector<double> ply_cz;     // c_z слоёв
    double ply_pitch = 0.0;         // Delta_ply = 2A + d + kEpsPly
    std::vector<PlyLayout> ply_layout;
    std::vector<PlyExpected> ply_expected;
    std::filesystem::path fib_path, layout_path;
    std::size_t n_stitch = 0;       // число связей-стежков (0 -- без стёжки)
    std::vector<std::pair<std::size_t, std::size_t>> stitch_pairs; // пары глобальных бусин стежков
    std::vector<std::array<double, 3>> stitch_r0; // начальные векторы пар x_j - x_i (по .fib, t = 0; не из чекпоинта)
    // Мозаика (MosaicSpec): ткань -- случайная раскройка, стежки через швы в слоях и сквозные между слоями.
    bool mosaic = false;
    std::vector<char> stitch_through; // [стежок] 1 -- сквозной между слоями; пусто без мозаики (стежки сетки -- как раньше)
    std::size_t n_stitch_in = 0, n_stitch_through = 0;
    int mosaic_cuts = 0;                // число раскроек пакета
    double mosaic_min_cross_dist = 0.0; // мин. расстояние в плане между бусинами разных фрагментов слоя, м
    double mosaic_size_x = 0.0, mosaic_size_y = 0.0, mosaic_gap = 0.0, mosaic_stitch_step = 0.0;
    bool mosaic_own_field = false;      // поле или c_e заданы явно (MosaicSpec::size_x/y, clamp_edge)
    std::filesystem::path cut_json_path;

    // --- диапазоны индексов Solid -----------------------------------------
    std::size_t ntotal = 0;
    std::size_t n_ball = 0, n_plate = 0;
    std::size_t ball_begin = 0, ball_end = 0;     // [0, n_ball)
    std::size_t plate_begin = 0, plate_end = 0;   // [n_ball, ntotal)
    std::vector<std::size_t> plate_back_layer;    // индексы тылового (нижнего по z) слоя пластины
    std::size_t moni_particle = 0;                // частица шара, ближайшая к его центру

    // --- шар --------------------------------------------------------------
    double rho_s = 0.0;             // плотность частиц шара, кг/м^3
    double ball_mass_target = 0.0;  // заданная масса шара, кг
    double ball_mass = 0.0;         // фактическая сумма масс частиц шара
    double R_eff = 0.0;             // max|r_k| + dx_s/2
    double z_c0 = 0.0;              // начальный z центра шара = z_f + g_0 + R_eff

    // --- геометрия мишени (ТЗ Guo 2019 §3.5) --------------------------------
    double z_f = kZFront;           // лицевая поверхность мишени
    double z_b = 0.0;               // тыльная поверхность мишени
    double eps_g = kEpsG;
    double fabric_cz = 0.0;         // c_z центра полотна (если есть)
    double fabric_amp = 0.0;        // A
    double fabric_d = 0.0;          // d
    double plate_z_top = 0.0, plate_z_bot = 0.0; // грани пластины (если есть)
    double plate_half_x = 0.0, plate_half_y = 0.0, plate_clamp_edge = 0.0;
    std::size_t plate_nx = 0, plate_ny = 0, plate_nz = 0;
    double plate_t = 0.0;                      // толщина пластины уровня, м
    double dx_p = 0.0;                         // номинальный шаг t/n_z
    double dx_p_axis[3] = {0, 0, 0};           // фактические шаги x, y, z
    double dx_p_mismatch = 0.0;                // (max - min)/min по осям
    double plate_mass_particle = 0.0;          // rho dx dy dz

    // --- ожидаемые числа по формулам ----------------------------------------
    struct Expected
    {
        std::size_t n_ball = 0;        // узлы решётки внутри радиуса
        std::size_t n_plate = 0;       // nx ny nz (cell_count)
        std::size_t ntotal = 0;
        std::size_t beads = 0;         // сумма pts.size()
        std::size_t bonds = 0;         // сумма (pts.size() - 1)
        std::size_t clamped_beads = 0; // сумма n_clamp_lo + n_clamp_hi
        std::size_t clamped_plate = 0; // nx nz (число рядов с |y| >= edge)
        std::size_t clamped_rows = 0;  // число рядов решётки y с |y| >= edge
    } expected;

    // Фактические числа из построенной Model (для сверки с expected).
    std::size_t actual_beads = 0, actual_bonds = 0;
    std::size_t actual_clamped_beads = 0, actual_clamped_plate = 0;

    // --- раскладка нитей ---------------------------------------------------
    std::vector<YarnLayout> yarns;

    // --- проверка узлов пластины на границе зажима -------------------------------
    std::size_t plate_nodes_on_edge = 0;      // узлов с |y| ровно равным plate_clamp_edge
    double plate_min_edge_dist = std::numeric_limits<double>::infinity(); // min ||y| - edge|

    // --- поверхностная плотность ткани -------------------------------------------
    double areal_density_window = 0.0;  // кг/м^2, только связи с обоими концами в окне
    double areal_density_summary = 0.0; // plain_weave_summary (все бусины, включая выступы)
    double window_Lx = 0.0, window_Ly = 0.0;
    std::size_t bonds_in_window = 0;
};

// ===========================================================================
// Шар: решётка с отбором по радиусу
// ===========================================================================
namespace detail
{

// Обход узлов (i+1/2) dx от центра; f(a, b, c, r) для узлов с r <= D/2 - dx/2 + tol.
template <class F>
void for_each_ball_node(double D, double dx, F&& f)
{
    const long M = static_cast<long>(std::ceil(0.5 * D / dx)) + 1;
    const double lim = 0.5 * D - 0.5 * dx + kBallRadiusTolDx * dx;
    for (long ix = -M; ix < M; ++ix)
        for (long iy = -M; iy < M; ++iy)
            for (long iz = -M; iz < M; ++iz)
            {
                const double a = (static_cast<double>(ix) + 0.5) * dx;
                const double b = (static_cast<double>(iy) + 0.5) * dx;
                const double c = (static_cast<double>(iz) + 0.5) * dx;
                const double r = std::sqrt(a * a + b * b + c * c);
                if (r <= lim) f(a, b, c, r);
            }
}

} // namespace detail

// Число узлов решётки шара внутри радиуса -- ожидаемое n_ball по формуле.
inline std::size_t ball_lattice_count(double D, double dx)
{
    std::size_t n = 0;
    detail::for_each_ball_node(D, dx, [&](double, double, double, double) { ++n; });
    return n;
}

// Координаты узлов пластины (общие для построения и для формул ожидаемых чисел).
inline double plate_node(double lo, double d, std::size_t i) { return lo + (static_cast<double>(i) + 0.5) * d; }

// ===========================================================================
// Пакет слоёв ткани (ТЗ 2T/4T §4.2)
// ===========================================================================
// Слой k (с нуля) -- копия q1 с cz = q1.cz - k * ply_pitch, построенная build_plain_weave.
// c_z считается умножением, без накопления суммы. n_plies == 1 даёт ровно build_plain_weave(q1).
inline std::vector<fiber::PlainWeave> build_fabric_stack(const fiber::PlainWeaveParams& q1, int n_plies,
                                                         double ply_pitch)
{
    if (n_plies < 1) throw std::invalid_argument("guo2019: build_fabric_stack: n_plies < 1");
    std::vector<fiber::PlainWeave> out;
    out.reserve(static_cast<std::size_t>(n_plies));
    for (int k = 0; k < n_plies; ++k)
    {
        fiber::PlainWeaveParams q = q1;
        q.cz = q1.cz - static_cast<double>(k) * ply_pitch;
        out.push_back(fiber::build_plain_weave(q));
    }
    return out;
}

namespace detail
{

// Материалы слоёв обязаны совпадать: в .fib пакета одна строка MATERIAL 0.
inline bool same_fabric_material(const fiber::PlainWeaveParams& a, const fiber::PlainWeaveParams& b)
{
    return a.E == b.E && a.nu == b.nu && a.d == b.d && a.rho == b.rho && a.eps_break == b.eps_break &&
           a.mu_ff == b.mu_ff && a.mu_fs == b.mu_fs;
}

} // namespace detail

// Стёжка слоёв пакета (этап 2в): сетка стежков на линиях шва «одеяла».
struct StitchSpec
{
    bool enabled = false;
    bool seam = false;              // true -- шов-строчка в плоскости слоя (>= 1 слоя), false -- сетка между слоями (>= 2)
    int tile = kFabricTileBase;     // сетка патчей tile x tile; nx и ny делятся на tile
    int patch = 0;                  // > 0: кусочки по patch нитей от центра, у краёв остаток (вместо tile)
    bool lines_x = true, lines_y = true;
    int every = 1;
    double E = 8.5e10, nu = 0.3, d = 0.2e-3, rho = 1440.0, eps_break = 0.03; // материал 1 стежка
};

// Мозаичная ткань: случайная раскройка (include/fiber/plain_weave_random_cut.hpp) вместо сплошного
// полотна, решение владельца 2026-11-13. Поле -- габарит сплошной ткани уровня: по x (ny - 1) p + 2 e_x
// (свободные концы по p/2), по y (nx - 1) p + 2 e_y (выступы в зажим); зажим -- боксы раскройки |y| >= c_e
// для всех бусин независимо от направления нити. Слои группы layers_per_cut делят одну раскройку.
// Стежки: через швы в каждом слое и сквозные между соседними слоями, материал 1 -- из StitchSpec
// (E, nu, d, rho, eps_break); --stitch=grid|seam вместе с мозаикой не допускается.
struct MosaicSpec
{
    bool enabled = false;
    int n_frag = 1;            // N фрагментов
    std::uint64_t seed = 1;
    int layers_per_cut = 0;    // k; 0 -- одна раскройка на весь пакет
    double gap = -1.0;         // зазор между фрагментами, м; < 0 -- d
    double rot_max_deg = 25.0; // поворот сетки фрагмента в [-rot_max, rot_max], градусы
    double stitch_step = -1.0; // шаг стежков вдоль шва, м; < 0 -- шаг нитей p (прокол на каждой нити)
    // Своё поле (решение владельца 2026-11-15): size_x x size_y, м, по центру (0, 0), и граница зажима c_e, м
    // (зажаты бусины |y| >= c_e). < 0 -- габарит сплошной ткани уровня и c_e уровня. Окно плотности при своём
    // поле -- свободная часть без полушага с краёв: (size_x - p) x (2 c_e - p).
    double size_x = -1.0, size_y = -1.0;
    double clamp_edge = -1.0;
};

// Дополнение .fib: MATERIAL 1 и строки BOND i j 1 (глобальные индексы бусин).
struct StitchFibExtra
{
    const StitchSpec* spec = nullptr;
    const std::vector<std::pair<std::size_t, std::size_t>>* pairs = nullptr;
};

// .fib пакета слоёв: одна MATERIAL 0, блоки FIBER всех слоёв по порядку, затем все CLAMP_BOX
// (правило полей mu_ff/mu_fs -- как у write_plain_weave_fib, include/fiber/plain_weave.hpp).
inline void write_fabric_stack_fib(std::ostream& os, const std::vector<fiber::PlainWeave>& plies,
                                   const std::string& command_line, const StitchFibExtra& stitch = StitchFibExtra{})
{
    if (plies.empty()) throw std::invalid_argument("guo2019: write_fabric_stack_fib: нет слоёв");
    const auto& q = plies.front().params;
    for (const auto& W : plies)
        if (!detail::same_fabric_material(q, W.params))
            throw std::runtime_error("guo2019: write_fabric_stack_fib: материалы слоёв различаются");

    std::size_t beads = 0, bonds = 0, clamped = 0, nfib = 0;
    for (const auto& W : plies)
    {
        const fiber::PlainWeaveSummary S = fiber::plain_weave_summary(W);
        beads += S.beads;
        bonds += S.bonds;
        clamped += S.clamped_beads;
        nfib += S.fibers;
    }
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "# Пакет слоёв полотняного переплетения (plain weave), сгенерировано plain_weave.hpp\n";
    os << "# команда: " << command_line << "\n";
    os << "# слоёв: " << plies.size() << ", Delta_ply = "
       << (plies.size() > 1 ? plies.front().params.cz - plies[1].params.cz : 0.0) << "\n";
    for (std::size_t k = 0; k < plies.size(); ++k)
        os << "# слой " << (k + 1) << ": c_z = " << plies[k].params.cz << "\n";
    os << "# d=" << q.d << " pitch=" << q.pitch << " amp=" << q.amp << " nx=" << q.nx << " ny=" << q.ny
       << " seg=" << q.seg << "\n";
    os << "# overhang_x=" << fiber::plain_weave_overhang_x(q) << " overhang_y=" << fiber::plain_weave_overhang_y(q)
       << " profile=" << (q.profile == fiber::PlainWeaveProfile::Cos ? "cos" : "arcs") << "\n";
    os << "# center=" << q.cx << "," << q.cy << " allow_penetration=" << (q.allow_penetration ? 1 : 0) << "\n";
    os << "# E=" << q.E << " nu=" << q.nu << " rho=" << q.rho << " eps_break=" << q.eps_break << " mu_ff=" << q.mu_ff
       << " mu_fs=" << q.mu_fs << "\n";
    const char* cl = q.clamp == fiber::PlainWeaveClamp::None ? "none"
                     : q.clamp == fiber::PlainWeaveClamp::X   ? "x"
                     : q.clamp == fiber::PlainWeaveClamp::Y   ? "y"
                                                              : "xy";
    os << "# clamp=" << cl;
    if (q.clamp != fiber::PlainWeaveClamp::None) os << " clamp_edge=" << q.clamp_edge;
    os << "\n";
    os << "# итого: слоёв " << plies.size() << ", нитей " << nfib << ", бусин " << beads << ", связей " << bonds
       << ", зажатых бусин " << clamped << "\n";
    for (std::size_t k = 0; k < plies.size(); ++k)
        for (const auto& w : plies[k].warnings) os << "# слой " << (k + 1) << ": " << w << "\n";

    os << "MATERIAL 0 " << q.E << " " << q.nu << " " << q.d << " " << q.rho << " " << q.eps_break;
    if (q.mu_fs != 0.0)
        os << " " << q.mu_ff << " " << q.mu_fs;
    else if (q.mu_ff != 0.0)
        os << " " << q.mu_ff;
    os << "\n";
    if (stitch.spec && stitch.pairs)
        os << "MATERIAL 1 " << stitch.spec->E << " " << stitch.spec->nu << " " << stitch.spec->d << " "
           << stitch.spec->rho << " " << stitch.spec->eps_break << "\n";

    for (const auto& W : plies)
        for (const auto& F : W.fibers)
        {
            os << "FIBER 0 " << F.pts.size() << "\n";
            for (const auto& P : F.pts) os << P[0] << " " << P[1] << " " << P[2] << "\n";
        }
    if (stitch.spec && stitch.pairs)
    {
        os << "# стёжка: " << stitch.pairs->size() << " стежков (материал 1)\n";
        for (const auto& pr : *stitch.pairs) os << "BOND " << pr.first << " " << pr.second << " 1\n";
    }
    for (std::size_t k = 0; k < plies.size(); ++k)
    {
        if (plies[k].clamp_boxes.empty()) continue;
        os << "# закрепление слоя " << (k + 1) << ": --clamp=" << cl << ", c_e = " << q.clamp_edge << "\n";
        for (const auto& b : plies[k].clamp_boxes)
            os << "CLAMP_BOX " << b.lo[0] << " " << b.lo[1] << " " << b.lo[2] << " " << b.hi[0] << " " << b.hi[1]
               << " " << b.hi[2] << "\n";
    }
}

// Раскладка пакета: строка на нить, первый столбец ply, object_id/first_bead/last_bead глобальные.
inline void write_fabric_stack_layout(std::ostream& os, const std::vector<fiber::PlainWeave>& plies)
{
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "# ply object_id direction index offset clamped first_bead last_bead n_clamp_lo n_clamp_hi\n";
    std::size_t first = 0;
    int oid = 0;
    for (std::size_t k = 0; k < plies.size(); ++k)
        for (const auto& F : plies[k].fibers)
        {
            const std::size_t last = first + F.pts.size() - 1;
            os << (k + 1) << " " << oid << " " << F.direction << " " << F.index << " " << F.offset << " "
               << (F.clamped ? 1 : 0) << " " << first << " " << last << " " << F.n_clamp_lo << " " << F.n_clamp_hi
               << "\n";
            first = last + 1;
            ++oid;
        }
}

// ===========================================================================
// Сборка сцены
// ===========================================================================
inline GuoScene build_guo_scene(Config config, Level level, double V0, const std::filesystem::path& workdir,
                                double mu = kMuDefault, bool body_isolation = true, int plate_nz_override = 0,
                                int plies = 0, double ball_d_override = 0.0, double ball_m_override = 0.0,
                                double ball_dx_override = 0.0, int fabric_n_override = 0,
                                double fabric_clamp_edge_override = 0.0, double fabric_overhang_y_override = 0.0,
                                const StitchSpec& stitch = StitchSpec{}, const MosaicSpec& mosaic = MosaicSpec{})
{
    // plies -- число слоёв для Config::TK; ball_*_override > 0 заменяют D_s, m_s, dx_s шара уровня
    // (решение владельца 2026-11-09); плотность по-прежнему выводится из массы: m_s / (n_ball dx_s^3).
    // fabric_*_override > 0 заменяют габариты ткани уровня (решение владельца 2026-11-14, малое полотно):
    // число нитей nx = ny, границу зажима c_e и выступ зажатых нитей e_y. Только для мишеней без пластины:
    // пластина своих габаритов не меняет, и окно ткани с ней разошлось бы.
    if (config == Config::TK && (plies < 1 || plies > kMaxPlies))
        throw std::runtime_error("guo2019: для конфигурации kT число слоёв должно быть от 1 до " + std::to_string(kMaxPlies));
    GuoScene sc;
    sc.config = config;
    sc.level = level;
    sc.spec = level_spec(level);
    if (plate_nz_override > 0) sc.spec.n_z = plate_nz_override; // только для сравнения стоимости
    if (ball_d_override > 0.0) sc.spec.D_s = ball_d_override;
    if (ball_m_override > 0.0) sc.spec.m_s = ball_m_override;
    if (ball_dx_override > 0.0) sc.spec.dx_s = ball_dx_override;
    if (fabric_n_override > 0 || fabric_clamp_edge_override > 0.0 || fabric_overhang_y_override > 0.0)
    {
        if (!config_has_fabric(config) || config_has_plate(config))
            throw std::runtime_error("guo2019: габариты ткани (--fabric-*) меняются только для мишеней без пластины (kT)");
        if (fabric_n_override > 0) sc.spec.nx = sc.spec.ny = fabric_n_override;
        if (fabric_clamp_edge_override > 0.0) sc.spec.c_e = fabric_clamp_edge_override;
        if (fabric_overhang_y_override > 0.0) sc.spec.e_y = fabric_overhang_y_override;
    }
    sc.V0 = V0;
    sc.mu = mu;
    sc.body_isolation = body_isolation;
    sc.has_fabric = config_has_fabric(config);
    sc.has_plate = config_has_plate(config);
    sc.n_plies = fabric_plies(config, plies);
    const LevelSpec& L = sc.spec;
    const double plate_t = L.t;
    sc.plate_t = plate_t;

    // ---------------------------------------------------------------- шар
    struct Node { double a, b, c; };
    std::vector<Node> ball;
    double max_r = 0.0;
    detail::for_each_ball_node(L.D_s, L.dx_s, [&](double a, double b, double c, double r) {
        ball.push_back({a, b, c});
        max_r = std::max(max_r, r);
    });
    if (ball.empty()) throw std::runtime_error("guo2019: решётка шара пуста (dx_s слишком велик)");
    sc.n_ball = ball.size();
    sc.ball_begin = 0;
    sc.ball_end = sc.n_ball;
    sc.ball_mass_target = L.m_s;
    const double dxs3 = L.dx_s * L.dx_s * L.dx_s;
    sc.rho_s = L.m_s / (static_cast<double>(sc.n_ball) * dxs3); // масса ровно m_s
    sc.R_eff = max_r + 0.5 * L.dx_s;
    sc.z_c0 = sc.z_f + L.g_0 + sc.R_eff;

    // -------------------------------------------------------------- ткань
    std::optional<fiber::RandomCut> rcut; // мозаика: раскройка (пишется в .fib ниже)
    if (mosaic.enabled && !sc.has_fabric) throw std::invalid_argument("guo2019: мозаика требует конфигурации с тканью");
    if (mosaic.enabled && stitch.enabled)
        throw std::invalid_argument("guo2019: мозаика несовместима с --stitch=grid|seam (стежки мозаики строятся по швам раскройки)");
    if (sc.has_fabric)
    {
        fiber::PlainWeaveParams q;
        fiber::apply_plain_weave_preset("twaron_ct717", q);
        q.amp = kFabricAmpOverD * q.d;                 // A = 0.6 d, решение §2 п.1
        q.seg = L.seg;
        q.nx = L.nx;
        q.ny = L.ny;
        q.overhang_x = 0.5 * q.pitch;                  // свободные нити (вдоль x): конец на середине полуволны
        q.overhang_y = L.e_y;                          // зажатые нити (вдоль y): выступ в зажим
        q.clamp = fiber::PlainWeaveClamp::Y;
        q.clamp_edge = L.c_e;
        q.clamp_edge_set = true;
        q.E = kFabricE;
        q.nu = kFabricNu;
        q.rho = kFabricRho;
        q.eps_break = kFabricEpsBreak;
        q.mu_ff = mu;
        q.mu_fs = mu;
        q.mu_fs_set = true;
        q.cx = 0.0;
        q.cy = 0.0;
        const double Aamp = q.amp, d = q.d;
        q.cz = (config == Config::A_T1) ? (-plate_t - kEpsG - Aamp - 0.5 * d) // A/1T
                                        : -(Aamp + 0.5 * d);                  // 1T, 1T/A
        sc.fabric_cz = q.cz;
        sc.fabric_amp = Aamp;
        sc.fabric_d = d;
        sc.ply_pitch = 2.0 * Aamp + d + kEpsPly; // Delta_ply (ТЗ 2T/4T §3.2)
        if (mosaic.enabled)
        {
            if ((mosaic.size_x > 0.0) != (mosaic.size_y > 0.0))
                throw std::invalid_argument("guo2019: мозаика: размер поля задаётся парой size_x, size_y");
            const double fx = mosaic.size_x > 0.0 ? mosaic.size_x
                                                  : static_cast<double>(q.ny - 1) * q.pitch + 2.0 * q.overhang_x;
            const double fy = mosaic.size_y > 0.0 ? mosaic.size_y
                                                  : static_cast<double>(q.nx - 1) * q.pitch + 2.0 * q.overhang_y;
            if (mosaic.clamp_edge > 0.0) q.clamp_edge = mosaic.clamp_edge;
            if (!(q.clamp_edge < 0.5 * fy))
                throw std::invalid_argument("guo2019: мозаика: граница зажима c_e должна быть меньше size_y / 2");
            sc.mosaic_own_field = mosaic.size_x > 0.0 || mosaic.clamp_edge > 0.0;
            fiber::RandomCutParams P;
            P.base = q;
            P.n_frag = mosaic.n_frag;
            P.size_x = fx;
            P.size_y = fy;
            P.seed = mosaic.seed;
            P.gap = mosaic.gap;
            P.rot_max_deg = mosaic.rot_max_deg;
            P.stitch_step = mosaic.stitch_step < 0.0 ? q.pitch : mosaic.stitch_step;
            P.stitch_E = stitch.E;
            P.stitch_nu = stitch.nu;
            P.stitch_d = stitch.d;
            P.stitch_rho = stitch.rho;
            P.stitch_eps_break = stitch.eps_break;
            P.layers = sc.n_plies;
            P.layers_per_cut = mosaic.layers_per_cut;
            P.ply_pitch = sc.ply_pitch;
            rcut = fiber::build_plain_weave_random_cut(P);
            const fiber::RandomCut& R = *rcut;
            // Слои -- куски полотна раскройки (нити в порядке записи); параметры -- сплошной ткани уровня
            // (окно, плотность summary), c_z -- свой у слоя.
            sc.plies.assign(static_cast<std::size_t>(sc.n_plies), fiber::PlainWeave{});
            for (int k = 0; k < sc.n_plies; ++k)
            {
                sc.plies[static_cast<std::size_t>(k)].params = q;
                sc.plies[static_cast<std::size_t>(k)].params.cz = R.ply_cz[static_cast<std::size_t>(k)];
            }
            for (std::size_t f = 0; f < R.W.fibers.size(); ++f)
                sc.plies[static_cast<std::size_t>(R.fiber_ply[f])].fibers.push_back(R.W.fibers[f]);
            sc.mosaic = true;
            sc.stitch_pairs = R.stitch_pairs;
            sc.stitch_through.clear();
            for (const auto& st : R.stitches) sc.stitch_through.push_back(st.through ? 1 : 0);
            if (sc.stitch_through.size() != sc.stitch_pairs.size())
                throw std::logic_error("guo2019: мозаика: число стежков и пар BOND различается");
            sc.n_stitch = sc.stitch_pairs.size();
            sc.n_stitch_in = R.n_stitches_in;
            sc.n_stitch_through = R.n_stitches_through;
            sc.mosaic_cuts = static_cast<int>(R.cuts.size());
            sc.mosaic_min_cross_dist = R.min_cross_dist;
            sc.mosaic_size_x = R.params.size_x;
            sc.mosaic_size_y = R.params.size_y;
            sc.mosaic_gap = R.params.gap;
            sc.mosaic_stitch_step = R.params.stitch_step;
        }
        else // Генератор сам проверяет диапазон c_e и бросает исключение.
            sc.plies = build_fabric_stack(q, sc.n_plies, sc.ply_pitch);
        sc.W = sc.plies.front(); // слой 1 (лицевой): на него опирается остальной код
        for (const auto& P : sc.plies) sc.ply_cz.push_back(P.params.cz);
        if (stitch.enabled)
        {
            if (stitch.seam ? sc.n_plies < 1 : sc.n_plies < 2)
                throw std::invalid_argument(stitch.seam ? "guo2019: шов-строчка требует >= 1 слоя ткани"
                                                        : "guo2019: стёжка требует пакета из >= 2 слоёв");
            std::vector<double> sx, sy;
            if (stitch.patch > 0)
            {
                // как в quilt_seams: x-линии -- по ny нитям вдоль y, y-линии -- по nx (у kT nx = ny)
                sx = fiber::centered_patch_seams(q.cx, q.pitch, q.ny, stitch.patch);
                sy = fiber::centered_patch_seams(q.cy, q.pitch, q.nx, stitch.patch);
            }
            else
            {
                if (stitch.tile < 1 || q.nx % stitch.tile != 0 || q.ny % stitch.tile != 0)
                    throw std::invalid_argument("guo2019: стёжка: nx и ny должны делиться на stitch-tile");
                fiber::quilt_seams(q, q.nx / stitch.tile, q.ny / stitch.tile, stitch.tile, stitch.tile, sx, sy);
            }
            fiber::QuiltStitchParams sp;
            sp.every = stitch.every;
            sp.lines_x = stitch.lines_x;
            sp.lines_y = stitch.lines_y;
            sc.stitch_pairs = stitch.seam ? fiber::build_quilt_seam_stitches(sc.plies, sx, sy, sp)
                                          : fiber::build_quilt_stitches(sc.plies, sx, sy, sp);
            sc.n_stitch = sc.stitch_pairs.size();
        }
    }

    // ----------------------------------------------------- грани и тыл мишени
    {
        const double A = sc.fabric_amp, d = sc.fabric_d;
        switch (config)
        {
        case Config::T1:
            sc.z_b = -(2.0 * A + d);
            break;
        case Config::A:
            sc.plate_z_top = 0.0;
            sc.plate_z_bot = -plate_t;
            sc.z_b = -plate_t;
            break;
        case Config::T1_A:
        {
            const double shift = 2.0 * A + d + kEpsG; // пластина сдвинута под ткань
            sc.plate_z_top = -shift;
            sc.plate_z_bot = -shift - plate_t;
            sc.z_b = sc.plate_z_bot;
            break;
        }
        case Config::A_T1:
            sc.plate_z_top = 0.0;
            sc.plate_z_bot = -plate_t;
            sc.z_b = -plate_t - kEpsG - 2.0 * A - d;
            break;
        case Config::T2:
        case Config::T4:
        case Config::TK:
            // Пакет без пластины: z_b = c_z,N - A - d/2 = -(N (2A + d) + (N - 1) eps_ply)
            sc.z_b = -(static_cast<double>(sc.n_plies) * (2.0 * A + d) + static_cast<double>(sc.n_plies - 1) * kEpsPly);
            break;
        }
    }

    // ----------------------------------------------------------- пластина
    if (sc.has_plate)
    {
        sc.plate_half_x = L.plate_half_x;
        sc.plate_half_y = L.plate_half_y;
        sc.plate_clamp_edge = L.plate_clamp_edge;
        sc.plate_nz = static_cast<std::size_t>(L.n_z);
        sc.dx_p = plate_t / static_cast<double>(L.n_z);
        // Правило lattice.hpp (cell_count, как у BC2): целое число ячеек, ближайшее к L/dx, шаг подгоняется.
        sc.plate_nx = lattice::cell_count(2.0 * sc.plate_half_x, sc.dx_p);
        sc.plate_ny = lattice::cell_count(2.0 * sc.plate_half_y, sc.dx_p);
        sc.dx_p_axis[0] = 2.0 * sc.plate_half_x / static_cast<double>(sc.plate_nx);
        sc.dx_p_axis[1] = 2.0 * sc.plate_half_y / static_cast<double>(sc.plate_ny);
        sc.dx_p_axis[2] = plate_t / static_cast<double>(sc.plate_nz);
        const double dmax = *std::max_element(sc.dx_p_axis, sc.dx_p_axis + 3);
        const double dmin = *std::min_element(sc.dx_p_axis, sc.dx_p_axis + 3);
        sc.dx_p_mismatch = (dmax - dmin) / dmin;
        sc.plate_mass_particle = kPlateRho * sc.dx_p_axis[0] * sc.dx_p_axis[1] * sc.dx_p_axis[2];
        sc.n_plate = sc.plate_nx * sc.plate_ny * sc.plate_nz;
    }
    sc.ntotal = sc.n_ball + sc.n_plate;
    sc.plate_begin = sc.n_ball;
    sc.plate_end = sc.ntotal;

    // ----------------------------------------------------- конфиг решателя
    {
        const std::filesystem::path p = std::filesystem::path(SPH_CASES_DIR) / "guo2019" / "experiment" / "initial_conditions.txt";
        if (!std::filesystem::exists(p)) throw std::runtime_error("guo2019: нет файла конфигурации " + p.string());
        sc.cfg.apply_defaults();
        {
            e2e::ScopedSilence silence(true);
            parse_sph_config(sc.cfg, p.string());
        }
        apply_guo_materials(sc.cfg.mats); // сталь в слоте 3 -- ДО urs()
        sc.cfg.sph_same_body_only = body_isolation;
        sc.cfg.maxn = static_cast<std::size_t>(1.4 * static_cast<double>(sc.ntotal)) + 1;
        sc.cfg.max_interaction = 100 * sc.cfg.maxn;
    }

    // ------------------------------------------------------ частицы Solid
    ParticleSystem& s = sc.model.solid;
    s.resize(sc.cfg.maxn);
    for (std::size_t i = 0; i < sc.cfg.maxn; ++i) s.inzializ_index(i);
    {
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < sc.n_ball; ++k)
        {
            s.x[0][k] = ball[k].a;
            s.x[1][k] = ball[k].b;
            s.x[2][k] = sc.z_c0 + ball[k].c;
            s.vx[2][k] = -V0;
            s.itype[k] = kBallItype;
            s.body[k] = kBallBody;
            s.rho[k] = sc.rho_s;
            s.rho_0[k] = sc.rho_s;
            s.mass[k] = sc.rho_s * dxs3;
            s.hsml[k] = 2.0 * L.dx_s;
            sc.ball_mass += s.mass[k];
            const double r = std::sqrt(ball[k].a * ball[k].a + ball[k].b * ball[k].b + ball[k].c * ball[k].c);
            if (r < best)
            {
                best = r;
                sc.moni_particle = k;
            }
        }
    }
    if (sc.has_plate)
    {
        const std::size_t nx = sc.plate_nx, ny = sc.plate_ny, nz = sc.plate_nz;
        const double dxx = sc.dx_p_axis[0], dyy = sc.dx_p_axis[1], dzz = sc.dx_p_axis[2];
        for (std::size_t ix = 0; ix < nx; ++ix)
            for (std::size_t iy = 0; iy < ny; ++iy)
                for (std::size_t iz = 0; iz < nz; ++iz)
                {
                    const std::size_t k = sc.n_ball + (ix * ny + iy) * nz + iz;
                    s.x[0][k] = plate_node(-sc.plate_half_x, dxx, ix);
                    s.x[1][k] = plate_node(-sc.plate_half_y, dyy, iy);
                    s.x[2][k] = plate_node(sc.plate_z_bot, dzz, iz); // iz = 0 -- тыльный слой
                    s.itype[k] = kPlateItype;
                    s.body[k] = kPlateBody;
                    s.rho[k] = kPlateRho;
                    s.rho_0[k] = kPlateRho;
                    s.mass[k] = sc.plate_mass_particle;
                    s.hsml[k] = 2.0 * sc.dx_p;
                    if (iz == 0) sc.plate_back_layer.push_back(k);
                    // Узлы на границе зажима: |y| ровно равно edge, и минимальное расстояние до неё.
                    const double ay = std::fabs(s.x[1][k]);
                    if (ay == sc.plate_clamp_edge) ++sc.plate_nodes_on_edge;
                    sc.plate_min_edge_dist = std::min(sc.plate_min_edge_dist, std::fabs(ay - sc.plate_clamp_edge));
                }
        // Ожидаемое число зажатых частиц пластины -- по формуле решётки (ряды y).
        for (std::size_t iy = 0; iy < ny; ++iy)
        {
            const double y = plate_node(-sc.plate_half_y, dyy, iy);
            if (y <= -sc.plate_clamp_edge || y >= sc.plate_clamp_edge) ++sc.expected.clamped_rows;
        }
        sc.expected.clamped_plate = nx * nz * sc.expected.clamped_rows;
    }
    for (std::size_t i = 0; i < sc.ntotal; ++i) urs(sc.cfg.mats, s, i); // cfg.mats со сталью уже применён

    sc.cfg.moni_particle = sc.moni_particle;

    // ----------------------------------------------------- нити (после Solid)
    if (sc.has_fabric)
    {
        std::filesystem::create_directories(workdir);
        sc.fib_path = workdir / "guo2019_fabric.fib";
        sc.layout_path = workdir / "guo2019_fabric_layout.txt";
        {
            std::ofstream f(sc.fib_path);
            if (!f) throw std::runtime_error("guo2019: не удалось открыть " + sc.fib_path.string());
            if (rcut)
            {
                fiber::PlainWeaveFibExtra ex = fiber::random_cut_fib_extra(*rcut);
                std::ostringstream c3; // материал стежка задан сценой, а не заглушкой генератора
                c3 << std::setprecision(17) << "стежки: шаг " << rcut->params.stitch_step << ", в слоях " << rcut->n_stitches_in
                   << ", сквозных " << rcut->n_stitches_through << "; материал 1 E=" << stitch.E << " nu=" << stitch.nu
                   << " d=" << stitch.d << " rho=" << stitch.rho << " eps_break=" << stitch.eps_break;
                if (ex.comments.size() >= 3) ex.comments[2] = c3.str();
                fiber::write_plain_weave_fib(f, rcut->W, "guo2019_scene build_guo_scene (мозаика)", ex);
            }
            else if (sc.n_plies == 1 && !stitch.enabled) // шов-строчка для 1T идёт через писатель пакета (MATERIAL 1, BOND)
                fiber::write_plain_weave_fib(f, sc.W, "guo2019_scene build_guo_scene");
            else
            {
                StitchFibExtra ex;
                if (stitch.enabled)
                {
                    ex.spec = &stitch;
                    ex.pairs = &sc.stitch_pairs;
                }
                write_fabric_stack_fib(f, sc.plies, "guo2019_scene build_guo_scene", ex);
            }
        }
        {
            std::ofstream f(sc.layout_path);
            if (!f) throw std::runtime_error("guo2019: не удалось открыть " + sc.layout_path.string());
            if (sc.n_plies == 1 && !rcut)
                fiber::write_plain_weave_layout(f, sc.W);
            else
                write_fabric_stack_layout(f, sc.plies);
        }
        if (rcut)
        {
            sc.cut_json_path = workdir / "guo2019_fabric.cut.json";
            std::ofstream f(sc.cut_json_path);
            if (!f) throw std::runtime_error("guo2019: не удалось открыть " + sc.cut_json_path.string());
            fiber::write_random_cut_json(f, *rcut);
        }
        fiber::load_fiber_file(sc.fib_path.string(), sc.model); // заполняет model.clamps.fiber_mask из CLAMP_BOX
        sc.stitch_r0.clear();
        sc.stitch_r0.reserve(sc.stitch_pairs.size());
        for (const auto& pr : sc.stitch_pairs)
        {
            const auto& x = sc.model.fiber.x;
            sc.stitch_r0.push_back({x[0][pr.second] - x[0][pr.first], x[1][pr.second] - x[1][pr.first],
                                    x[2][pr.second] - x[2][pr.first]});
        }

        std::size_t first = 0;
        int oid = 0; // глобальный object_id: (k - 1)(nx + ny) + локальный id
        for (int k = 0; k < sc.n_plies; ++k)
        {
            PlyLayout pl;
            pl.bead_begin = first;
            pl.yarn_begin = sc.yarns.size();
            pl.object_id_begin = oid;
            PlyExpected pe;
            for (const auto& F : sc.plies[static_cast<std::size_t>(k)].fibers)
            {
                YarnLayout y;
                y.ply = k + 1;
                y.direction = F.direction;
                y.index = F.index;
                y.object_id = oid++;
                y.bead_begin = first;
                y.bead_end = first + F.pts.size();
                y.n_clamp_lo = F.n_clamp_lo;
                y.n_clamp_hi = F.n_clamp_hi;
                y.clamped = F.clamped;
                sc.yarns.push_back(y);
                first = y.bead_end;
                pe.beads += F.pts.size();
                pe.bonds += F.pts.size() - 1;
                pe.clamped_beads += F.n_clamp_lo + F.n_clamp_hi;
            }
            pl.bead_end = first;
            pl.yarn_end = sc.yarns.size();
            sc.ply_layout.push_back(pl);
            sc.ply_expected.push_back(pe);
            sc.expected.beads += pe.beads;
            sc.expected.bonds += pe.bonds;
            sc.expected.clamped_beads += pe.clamped_beads;
        }
    }
    else
    {
        sc.model.fiber.resize(0);
        sc.model.topology.resize(0);
    }

    // -------------------------------------------------- закрепление пластины
    sc.model.clamps.solid_mask.assign(sc.ntotal, 0);
    if (sc.has_plate)
    {
        rhs::CpuExecutor exec(1);
        const double big = 1.0; // 1 м -- с запасом за габариты пластины по x и z
        const double edge = sc.plate_clamp_edge;
        const double lo_a[dim] = {-big, -big, -big}, hi_a[dim] = {big, -edge, big}; // y <= -edge
        const double lo_b[dim] = {-big, edge, -big}, hi_b[dim] = {big, big, big};   // y >= edge
        interactions::mark_clamp_box(exec, sc.model.solid.x, sc.ntotal, lo_a, hi_a, interactions::kClampAll,
                                     sc.model.clamps.solid_mask);
        interactions::mark_clamp_box(exec, sc.model.solid.x, sc.ntotal, lo_b, hi_b, interactions::kClampAll,
                                     sc.model.clamps.solid_mask);
    }
    for (std::size_t i = sc.ball_begin; i < sc.ball_end; ++i) sc.model.clamps.solid_mask[i] = 0; // шар не зажимается

    // ------------------------------------------------ ожидаемые числа, факт
    sc.expected.n_ball = ball_lattice_count(L.D_s, L.dx_s);
    sc.expected.n_plate = sc.has_plate ? sc.plate_nx * sc.plate_ny * sc.plate_nz : 0;
    sc.expected.ntotal = sc.expected.n_ball + sc.expected.n_plate;
    sc.actual_beads = sc.model.fiber.size();
    sc.actual_bonds = sc.model.topology.bonds.size() - sc.n_stitch; // стежки в «связи ткани» не входят
    for (auto m : sc.model.clamps.fiber_mask) if (m != 0) ++sc.actual_clamped_beads;
    for (std::size_t i = sc.plate_begin; i < sc.plate_end; ++i)
        if (sc.model.clamps.solid_mask[i] != 0) ++sc.actual_clamped_plate;

    // -------------------------------- поверхностная плотность ткани (в окне)
    if (sc.has_fabric)
    {
        const auto& q = sc.W.params;
        fiber::plain_weave_window(q, sc.window_Lx, sc.window_Ly);
        if (sc.mosaic && sc.mosaic_own_field) // своё поле мозаики: свободная часть без полушага с краёв
        {
            sc.window_Lx = sc.mosaic_size_x - q.pitch;
            sc.window_Ly = 2.0 * q.clamp_edge - q.pitch;
        }
        const double tol = 1.0e-9 * q.pitch;
        const double hx = 0.5 * sc.window_Lx + tol, hy = 0.5 * sc.window_Ly + tol;
        const auto& fb = sc.model.fiber;
        auto inside = [&](std::size_t i) {
            return std::fabs(fb.x[0][i] - q.cx) <= hx && std::fabs(fb.x[1][i] - q.cy) <= hy;
        };
        const double pi = std::acos(-1.0);
        double mass = 0.0;
        std::vector<double> ply_mass(static_cast<std::size_t>(sc.n_plies), 0.0);
        for (const auto& b : sc.model.topology.bonds)
        {
            if (b.mat >= 0) continue; // стежок
            if (!inside(b.i) || !inside(b.j)) continue;
            std::size_t kply = 0; // слой связи по бусине b.i (число слоёв произвольное)
            while (kply + 1 < sc.ply_layout.size() && b.i >= sc.ply_layout[kply].bead_end) ++kply;
            const auto& mat = sc.model.materials[static_cast<std::size_t>(fb.material_id[b.i])];
            const double dm = mat.rho * pi * mat.d * mat.d / 4.0 * b.l0;
            mass += dm;
            ply_mass[kply] += dm;
            ++sc.bonds_in_window;
            ++sc.ply_expected[kply].bonds_in_window;
        }
        sc.areal_density_window = mass / (sc.window_Lx * sc.window_Ly);
        sc.areal_density_summary = 0.0;
        for (std::size_t k = 0; k < sc.ply_expected.size(); ++k)
        {
            sc.ply_expected[k].areal_density_window = ply_mass[k] / (sc.window_Lx * sc.window_Ly);
            sc.ply_expected[k].areal_density_summary = fiber::plain_weave_summary(sc.plies[k]).areal_density;
            sc.areal_density_summary += sc.ply_expected[k].areal_density_summary;
        }
    }

    return sc;
}

// ===========================================================================
// S2m: зазоры между слоями пакета (ТЗ 2T/4T §4.2)
// ===========================================================================
struct PlyGapReport
{
    bool valid = false; // есть хотя бы одна пара соседних слоёв
    double min_gap = std::numeric_limits<double>::infinity();
    double max_gap = -std::numeric_limits<double>::infinity();
    int min_ply = 0, min_i = 0, min_j = 0; // где достигнут минимум: пара слоёв (min_ply, min_ply + 1), пересечение (i, j)
    int max_ply = 0, max_i = 0, max_j = 0;
    std::size_t n_pairs = 0;               // число проверенных (пара слоёв, пересечение)
};

// Для пары соседних слоёв (k, k+1) и каждого пересечения (i, j): нижняя по z из двух вершинных бусин
// слоя k и верхняя из двух вершинных бусин слоя k+1; зазор = разность z минус d. Вершины: у основы j
// (вдоль x) вершина i, у утка i (вдоль y) вершина j. По координатам model.fiber (начальное состояние).
inline PlyGapReport ply_vertex_gaps(const GuoScene& sc)
{
    PlyGapReport r;
    if (sc.n_plies < 2 || sc.mosaic) return r; // мозаика: сетки пересечений nx x ny нет
    r.valid = true;
    const auto& fb = sc.model.fiber;
    const int nx = sc.W.params.nx, ny = sc.W.params.ny;
    const double d = sc.fabric_d;
    auto bead = [&](int ply, char dir, int idx, int vert) {
        // нити слоя: сначала nx основы (j = 0..nx-1), затем ny утка (i = 0..ny-1)
        const std::size_t loc = (dir == 'x') ? static_cast<std::size_t>(idx) : static_cast<std::size_t>(nx + idx);
        const std::size_t yi = sc.ply_layout[static_cast<std::size_t>(ply)].yarn_begin + loc;
        const auto& F = sc.plies[static_cast<std::size_t>(ply)].fibers[loc];
        return sc.yarns[yi].bead_begin + F.vertex_bead[static_cast<std::size_t>(vert)];
    };
    for (int k = 0; k + 1 < sc.n_plies; ++k)
        for (int j = 0; j < nx; ++j)
            for (int i = 0; i < ny; ++i)
            {
                const double za = fb.x[2][bead(k, 'x', j, i)], zb = fb.x[2][bead(k, 'y', i, j)];
                const double zc = fb.x[2][bead(k + 1, 'x', j, i)], zd = fb.x[2][bead(k + 1, 'y', i, j)];
                const double gap = std::min(za, zb) - std::max(zc, zd) - d;
                ++r.n_pairs;
                if (gap < r.min_gap)
                {
                    r.min_gap = gap;
                    r.min_ply = k + 1;
                    r.min_i = i;
                    r.min_j = j;
                }
                if (gap > r.max_gap)
                {
                    r.max_gap = gap;
                    r.max_ply = k + 1;
                    r.max_i = i;
                    r.max_j = j;
                }
            }
    return r;
}

// ===========================================================================
// S2: геометрия укладки без прогона
// ===========================================================================
struct GapSide
{
    bool found = false;       // найдена хотя бы одна пара в окне поиска
    double gap = std::numeric_limits<double>::infinity(); // min (dist - r_f - r_s), м
    double dist = 0.0, r_f = 0.0, r_s = 0.0;
    std::size_t bond = 0;     // индекс в model.topology.bonds
    std::size_t particle = 0; // глобальный индекс Solid-частицы
    double seg_point[3] = {0, 0, 0}; // ближайшая точка оси сегмента
    double part_pos[3] = {0, 0, 0};  // центр частицы
};

struct GapReport
{
    GapSide ball, plate;
    double margin = 0.0; // окно поиска: пары с зазором больше margin не рассматриваются
};

namespace detail
{

inline double point_segment_dist(const double p[3], const double a[3], const double b[3], double q[3])
{
    double ab[3], ap[3];
    double ab2 = 0.0, t = 0.0;
    for (int d = 0; d < 3; ++d)
    {
        ab[d] = b[d] - a[d];
        ap[d] = p[d] - a[d];
        ab2 += ab[d] * ab[d];
        t += ap[d] * ab[d];
    }
    t = ab2 > 0.0 ? std::min(1.0, std::max(0.0, t / ab2)) : 0.0;
    double d2 = 0.0;
    for (int d = 0; d < 3; ++d)
    {
        q[d] = a[d] + t * ab[d];
        const double e = p[d] - q[d];
        d2 += e * e;
    }
    return std::sqrt(d2);
}

inline std::uint64_t cell_key(long long ix, long long iy, long long iz)
{
    constexpr long long off = 1LL << 20;
    constexpr std::uint64_t m = (1ULL << 21) - 1ULL;
    return (static_cast<std::uint64_t>(ix + off) & m) | ((static_cast<std::uint64_t>(iy + off) & m) << 21) |
           ((static_cast<std::uint64_t>(iz + off) & m) << 42);
}

// Минимум зазора по связям нити и Solid-частицам из [i0, i1); равномерная сетка, O(N).
inline GapSide min_gap_range(const GuoScene& sc, std::size_t i0, std::size_t i1, double margin)
{
    GapSide res;
    const auto& s = sc.model.solid;
    const auto& fb = sc.model.fiber;
    const auto& bonds = sc.model.topology.bonds;
    if (i1 <= i0 || bonds.empty()) return res;

    double rs_max = 0.0, org[3] = {1e300, 1e300, 1e300};
    for (std::size_t k = i0; k < i1; ++k)
    {
        rs_max = std::max(rs_max, interactions::solid_contact_radius(s, k));
        for (int d = 0; d < 3; ++d) org[d] = std::min(org[d], s.x[d][k]);
    }
    double rf_max = 0.0;
    for (const auto& m : sc.model.materials) rf_max = std::max(rf_max, 0.5 * m.d);
    const double W = rf_max + rs_max + margin; // размер ячейки = радиус окна
    if (!(W > 0.0)) return res;

    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> grid;
    grid.reserve((i1 - i0) / 4 + 16);
    auto cell = [&](double x, int d) { return static_cast<long long>(std::floor((x - org[d]) / W)); };
    for (std::size_t k = i0; k < i1; ++k)
        grid[cell_key(cell(s.x[0][k], 0), cell(s.x[1][k], 1), cell(s.x[2][k], 2))].push_back(
            static_cast<std::uint32_t>(k));

    for (std::size_t bi = 0; bi < bonds.size(); ++bi)
    {
        if (bonds[bi].mat >= 0) continue; // стежок -- не сегмент нити
        double a[3], b[3], lo[3], hi[3];
        for (int d = 0; d < 3; ++d)
        {
            a[d] = fb.x[d][bonds[bi].i];
            b[d] = fb.x[d][bonds[bi].j];
            lo[d] = std::min(a[d], b[d]) - W;
            hi[d] = std::max(a[d], b[d]) + W;
        }
        const double rf = 0.5 * sc.model.materials[static_cast<std::size_t>(fb.material_id[bonds[bi].i])].d;
        const long long cx0 = cell(lo[0], 0), cx1 = cell(hi[0], 0);
        const long long cy0 = cell(lo[1], 1), cy1 = cell(hi[1], 1);
        const long long cz0 = cell(lo[2], 2), cz1 = cell(hi[2], 2);
        for (long long cx = cx0; cx <= cx1; ++cx)
            for (long long cy = cy0; cy <= cy1; ++cy)
                for (long long cz = cz0; cz <= cz1; ++cz)
                {
                    const auto it = grid.find(cell_key(cx, cy, cz));
                    if (it == grid.end()) continue;
                    for (std::uint32_t k : it->second)
                    {
                        const double p[3] = {s.x[0][k], s.x[1][k], s.x[2][k]};
                        double q[3];
                        const double dist = point_segment_dist(p, a, b, q);
                        const double rs = interactions::solid_contact_radius(s, k);
                        const double gap = dist - rf - rs;
                        if (gap < res.gap)
                        {
                            res.found = true;
                            res.gap = gap;
                            res.dist = dist;
                            res.r_f = rf;
                            res.r_s = rs;
                            res.bond = bi;
                            res.particle = k;
                            for (int d = 0; d < 3; ++d)
                            {
                                res.seg_point[d] = q[d];
                                res.part_pos[d] = p[d];
                            }
                        }
                    }
                }
    }
    return res;
}

} // namespace detail

// S2 (ТЗ Guo 2019 §3.5, §9): min по связям нити и Solid-частицам из
// (расстояние от центра частицы до оси сегмента) - (r_f + r_s), отдельно для
// шара и пластины. Верхней границы у зазора нет; ищутся пары с зазором меньше
// margin (по умолчанию 2 мм): если в окне пар нет, found == false (зазор
// больше margin). Без ткани отчёт пуст.
inline GapReport fabric_solid_min_gap(const GuoScene& sc, double margin = 2.0e-3)
{
    GapReport r;
    r.margin = margin;
    if (!sc.has_fabric) return r;
    r.ball = detail::min_gap_range(sc, sc.ball_begin, sc.ball_end, margin);
    r.plate = detail::min_gap_range(sc, sc.plate_begin, sc.plate_end, margin);
    return r;
}

// Зазор шар--мишень (ТЗ Guo 2019 §3.5, §9 S2).
struct BallGap
{
    double analytic = 0.0;     // z_c0 - R_eff - z_f; обязано равняться g_0 с точностью kTolGap
    double particle_min = 0.0; // min_k (z_k - dx_s/2) - z_f: по граням ячеек частиц шара (>= analytic)
    double contact_min = 0.0;  // min_k (z_k - r_s,k) - z_f: по контактным радиусам частиц
};

inline BallGap ball_target_gap(const GuoScene& sc)
{
    BallGap g;
    g.analytic = sc.z_c0 - sc.R_eff - sc.z_f;
    g.particle_min = std::numeric_limits<double>::infinity();
    g.contact_min = std::numeric_limits<double>::infinity();
    const auto& s = sc.model.solid;
    for (std::size_t k = sc.ball_begin; k < sc.ball_end; ++k)
    {
        g.particle_min = std::min(g.particle_min, s.x[2][k] - 0.5 * sc.spec.dx_s - sc.z_f);
        g.contact_min = std::min(g.contact_min, s.x[2][k] - interactions::solid_contact_radius(s, k) - sc.z_f);
    }
    return g;
}

// ===========================================================================
// Печать сводки сцены
// ===========================================================================
inline void print_scene_summary(std::FILE* f, const GuoScene& sc)
{
    auto cmp = [&](const char* name, std::size_t actual, std::size_t expected) {
        std::fprintf(f, "  %-34s %zu (формула %zu) %s\n", name, actual, expected, actual == expected ? "OK" : "РАСХОЖДЕНИЕ");
    };
    const LevelSpec& L = sc.spec;
    std::fprintf(f, "Сцена Guo 2019: конфигурация %s, уровень %s, V0 = %.3f м/с, mu = %.3f, изоляция тел = %d\n",
                 config_name_str(sc.config, sc.n_plies).c_str(), level_name(sc.level), sc.V0, sc.mu, sc.body_isolation ? 1 : 0);
    std::fprintf(f, "Частицы и нити (факт против формулы):\n");
    cmp("частиц шара", sc.n_ball, sc.expected.n_ball);
    cmp("частиц пластины", sc.n_plate, sc.expected.n_plate);
    cmp("Solid всего", sc.ntotal, sc.expected.ntotal);
    if (sc.has_fabric)
    {
        cmp("бусин", sc.actual_beads, sc.expected.beads);
        cmp("связей", sc.actual_bonds, sc.expected.bonds);
        cmp("зажатых бусин", sc.actual_clamped_beads, sc.expected.clamped_beads);
    }
    if (sc.has_plate) cmp("зажатых частиц пластины", sc.actual_clamped_plate, sc.expected.clamped_plate);

    std::fprintf(f, "Шар: D_s = %.6g м, dx_s = %.6g м, m_s = %.9g кг (факт %.9g, отн. откл. %.3e)\n", L.D_s, L.dx_s,
                 sc.ball_mass_target, sc.ball_mass, (sc.ball_mass - sc.ball_mass_target) / sc.ball_mass_target);
    std::fprintf(f, "  rho_s = %.6g кг/м^3 (сплошной шар %.0f x %.3g = %.0f, откл. %.2f %%)\n", sc.rho_s, kRhoSteelRef, L.ball_rho_factor,
                 kRhoSteelRef * L.ball_rho_factor,
                 100.0 * (sc.rho_s - kRhoSteelRef * L.ball_rho_factor) / (kRhoSteelRef * L.ball_rho_factor));
    std::fprintf(f, "  R_eff = %.6g м против D_s/2 = %.6g м (разность %.3e м); z_c0 = %.9g м, g_0 = %.3g м\n", sc.R_eff,
                 0.5 * L.D_s, sc.R_eff - 0.5 * L.D_s, sc.z_c0, L.g_0);
    if (const auto es = elastic_by_itype(sc.cfg.mats, kBallItype))
        std::fprintf(f, "  сталь (слот itype %d): E = %.4g Па, nu = %.4f (ожидание ~%.0f ГПа, ~%.2f)\n", kBallItype, es->E,
                     es->nu, kSteelERef * 1e-9, kSteelNuRef);
    else
        std::fprintf(f, "  сталь: elastic_by_itype не определён\n");
    std::fprintf(f, "  eps_g = %.3e м\n", sc.eps_g);
    std::fprintf(f, "Мишень: z_f = %.9g м, z_b = %.9g м\n", sc.z_f, sc.z_b);
    if (sc.has_plate)
    {
        std::fprintf(f, "Пластина: грани z in [%.9g; %.9g] м, %zu x %zu x %zu узлов, зажим |y| >= %.4g м\n",
                     sc.plate_z_bot, sc.plate_z_top, sc.plate_nx, sc.plate_ny, sc.plate_nz, sc.plate_clamp_edge);
        std::fprintf(f, "  шаги dx = %.6g, dy = %.6g, dz = %.6g м (номинал t/n_z = %.6g); расхождение по осям %.3f %%%s\n",
                     sc.dx_p_axis[0], sc.dx_p_axis[1], sc.dx_p_axis[2], sc.dx_p, 100.0 * sc.dx_p_mismatch,
                     sc.dx_p_mismatch > kPlateDxMismatchWarn ? " -- БОЛЬШЕ 2 %, находка" : "");
        std::fprintf(f, "  узлов на |y| == edge ровно: %zu, min ||y| - edge| = %.6g м; рядов в зажиме по y: %zu\n",
                     sc.plate_nodes_on_edge, sc.plate_min_edge_dist, sc.expected.clamped_rows);
    }
    if (sc.has_fabric)
    {
        const auto& q = sc.W.params;
        std::fprintf(f, "Ткань: d = %.4g м, p = %.4g м, A = %.4g м, c_z = %.9g м, k = %d, %d x %d нитей, зажим |y| >= %.4g м\n",
                     q.d, q.pitch, q.amp, q.cz, q.seg, q.nx, q.ny, q.clamp_edge);
        std::size_t ny_cl = 0;
        for (const auto& y : sc.yarns) if (y.clamped) ++ny_cl;
        std::fprintf(f, "  нитей: %zu, зажатых (вдоль y): %zu\n", sc.yarns.size(), ny_cl);
        std::fprintf(f, "  поверхностная плотность в окне (только связи целиком внутри, %zu связей, окно %.4g x %.4g м): %.5f кг/м^2\n",
                     sc.bonds_in_window, sc.window_Lx, sc.window_Ly, sc.areal_density_window);
        std::fprintf(f, "  plain_weave_summary (все бусины, включая выступы): %.5f кг/м^2; статья %.2f, ожидание модели %.3f (k = 2) / %.3f (k = 4)\n",
                     sc.areal_density_summary, kArealDensityPaper, kArealDensityExpectK2, kArealDensityExpectK4);
        std::fprintf(f, "  толщина полотна 2A + d = %.4g м\n", 2.0 * q.amp + q.d);
        if (sc.mosaic)
            std::fprintf(f, "  мозаика: поле %.6g x %.6g м, раскроек %d, gap = %.4g м, мин. расстояние между фрагментами %.4g м; "
                         "стежков %zu (в слоях %zu, сквозных %zu), шаг вдоль шва %.4g м\n",
                         sc.mosaic_size_x, sc.mosaic_size_y, sc.mosaic_cuts, sc.mosaic_gap, sc.mosaic_min_cross_dist, sc.n_stitch,
                         sc.n_stitch_in, sc.n_stitch_through, sc.mosaic_stitch_step);
    }
    if (sc.has_fabric && sc.n_plies >= 2)
    {
        std::fprintf(f, "Пакет слоёв: %d, Delta_ply = 2A + d + eps_ply = %.9g м (eps_ply = %.3e м)\n", sc.n_plies, sc.ply_pitch,
                     kEpsPly);
        std::size_t sb = 0, sn = 0, scl = 0;
        double sum_w = 0.0, sum_s = 0.0;
        for (int k = 0; k < sc.n_plies; ++k)
        {
            const auto& pl = sc.ply_layout[static_cast<std::size_t>(k)];
            const auto& pe = sc.ply_expected[static_cast<std::size_t>(k)];
            std::fprintf(f, "  слой %d: c_z = %.9g м, бусин / связей / зажатых = %zu / %zu / %zu, нити %zu..%zu (object_id с %d), "
                         "плотность в окне %.5f кг/м^2 (summary %.5f)\n",
                         k + 1, sc.ply_cz[static_cast<std::size_t>(k)], pe.beads, pe.bonds, pe.clamped_beads, pl.yarn_begin,
                         pl.yarn_end - 1, pl.object_id_begin, pe.areal_density_window, pe.areal_density_summary);
            sb += pe.beads;
            sn += pe.bonds;
            scl += pe.clamped_beads;
            sum_w += pe.areal_density_window;
            sum_s += pe.areal_density_summary;
        }
        std::fprintf(f, "  сумма: бусин / связей / зажатых = %zu / %zu / %zu; плотность в окне %.5f кг/м^2, summary %.5f "
                     "(статья %.2f / %.2f для 2T / 4T)\n",
                     sb, sn, scl, sum_w, sum_s, 2.0 * kArealDensityPaper, 4.0 * kArealDensityPaper);
        const PlyGapReport pg = ply_vertex_gaps(sc);
        if (!pg.valid)
            std::fprintf(f, "  зазор слоёв по вершинным бусинам: не считается (мозаика)\n");
        else
        std::fprintf(f, "  зазор слоёв по вершинным бусинам (мин / макс; ожидание eps_ply = %.3e м): %.12e / %.12e м по %zu "
                     "пересечениям; мин в паре слоёв %d, пересечение (%d, %d); макс в паре %d, (%d, %d)\n",
                     kEpsPly, pg.min_gap, pg.max_gap, pg.n_pairs, pg.min_ply, pg.min_i, pg.min_j, pg.max_ply, pg.max_i,
                     pg.max_j);
    }
    const BallGap bg = ball_target_gap(sc);
    std::fprintf(f, "Зазор шар--мишень: аналитический %.12e м (g_0 = %.3e, откл. %.2e); по граням частиц %.6e м; по контактным радиусам %.6e м\n",
                 bg.analytic, L.g_0, bg.analytic - L.g_0, bg.particle_min, bg.contact_min);
}

} // namespace guo2019
