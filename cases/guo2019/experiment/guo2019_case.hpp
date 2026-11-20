#pragma once
// guo2019_case.hpp
// ЭКСПЕРИМЕНТАЛЬНАЯ КОПИЯ (cases/README.md, «Роли постановок»): ответвлена от
// cases/guo2019/regression/ в a9c3964; меняется свободно под расчёты
// (guo2019_shot, fiber_t0_dump, fabric_checks). Тесты её не используют, кроме
// проверок самих этих утилит. Отличия: diff -r cases/guo2019/regression cases/guo2019/experiment.
// ВСЕ числа постановки опытов Guo et al., IJIE 127 (2019) 1-16: удар стального
// шара по мишеням A (Al-пластина), 1T (один слой ткани Twaron CT 717), 1T/A
// (ткань спереди, Al сзади), A/1T (Al спереди, ткань сзади), 2T и 4T (два и
// четыре несклеенных слоя ткани; ТЗ 2T/4T -- notes/guo2019/guo2019-2t-4t-tz.md).
// Постановка: ТЗ Guo 2019 §3 (геометрия и материалы), §7 (уровни разрешения),
// §8 (серии выстрелов), §9 (допуски smoke).
//
// Здесь и только здесь живут числа самого опыта и допуски приёмки. Построение
// сцены -- cases/guo2019/experiment/guo2019_scene.hpp. Параметры РЕШАТЕЛЯ лежат отдельно, в
// cases/guo2019/experiment/initial_conditions.txt, и разбираются боевым
// parse_sph_config().
//
// ЕДИНИЦЫ: метры, секунды, м/с, кг, кг/м^3, Па. Постановка записана в мм/мкс/г,
// поэтому у каждой константы в комментарии стоит исходное значение и пометка
// источника: «статья» (Guo 2019, текст), «рис.» (снято с рисунка, +-1-2 м/с),
// «оценка» (нет в статье, принято исполнителем ТЗ; пересчитывается и
// фиксируется в отчёте), «решение §N» (принято в ТЗ Guo 2019).

#include "materials.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace guo2019
{

// ===========================================================================
// Конфигурации мишени и уровни разрешения
// ===========================================================================
enum class Config
{
    A,    // одна Al-пластина (опорная конфигурация для нормировки e/(e_A + e_1T))
    T1,   // один слой ткани, имя "1T"
    T1_A, // ткань спереди, Al сзади, имя "1T_A" (в статье 1T/A)
    A_T1, // Al спереди, ткань сзади, имя "A_1T" (в статье A/1T)
    T2,   // два слоя ткани без пластины, имя "2T"
    T4,   // четыре слоя ткани без пластины, имя "4T"
    TK    // k слоёв ткани без пластины, k = 3, 5..kMaxPlies (1, 2, 4 -- это T1, T2, T4); имя "<k>T"; k хранится отдельно
};

// Предел числа слоёв ткани для Config::TK (решение владельца 2026-11-09).
inline constexpr int kMaxPlies = 64;

inline const char* config_name(Config c)
{
    switch (c)
    {
    case Config::A: return "A";
    case Config::T1: return "1T";
    case Config::T1_A: return "1T_A";
    case Config::A_T1: return "A_1T";
    case Config::T2: return "2T";
    case Config::T4: return "4T";
    case Config::TK: return "kT"; // число слоёв неизвестно: используйте config_name_str(c, plies)
    }
    return "?";
}

// Имя с учётом числа слоёв (для TK -- "<plies>T"; для остальных конфигураций plies игнорируется).
inline std::string config_name_str(Config c, int plies)
{
    return c == Config::TK ? std::to_string(plies) + "T" : std::string(config_name(c));
}

// true -- имя распознано (A, 1T, 1T_A, A_1T; также допускаются 1T/A и A/1T; "<k>T", k = 1..kMaxPlies).
// plies -- число слоёв ткани для "<k>T" (0 для остальных имён). "1T", "2T", "4T" дают T1, T2, T4 (как раньше),
// прочие k -- Config::TK. Причина отказа по диапазону k -- в err.
inline bool parse_config(const std::string& s, Config& out, int& plies, std::string* err = nullptr)
{
    plies = 0;
    if (s.size() >= 2 && s.back() == 'T' && s[0] >= '1' && s[0] <= '9' &&
        s.find_first_not_of("0123456789") == s.size() - 1)
    {
        const std::string digits = s.substr(0, s.size() - 1);
        const int k = digits.size() > 4 ? kMaxPlies + 1 : std::stoi(digits);
        if (k > kMaxPlies)
        {
            if (err) *err = "число слоёв ткани должно быть от 1 до " + std::to_string(kMaxPlies) + ": " + s;
            return false;
        }
        plies = k;
        out = k == 1 ? Config::T1 : k == 2 ? Config::T2 : k == 4 ? Config::T4 : Config::TK;
        return true;
    }
    if (s == "A") out = Config::A;
    else if (s == "1T") out = Config::T1;
    else if (s == "1T_A" || s == "1T/A") out = Config::T1_A;
    else if (s == "A_1T" || s == "A/1T") out = Config::A_T1;
    else return false;
    return true;
}

// Прежний интерфейс: только именованные конфигурации (TK без числа слоёв смысла не имеет, поэтому не принимается).
inline bool parse_config(const std::string& s, Config& out)
{
    int plies = 0;
    return parse_config(s, out, plies) && out != Config::TK;
}

// Число слоёв ткани в мишени: 0 у A, 1 у 1T, 1T_A, A_1T, 2 у 2T, 4 у 4T.
// Для Config::TK число слоёв берётся из plies (иначе plies игнорируется).
inline int fabric_plies(Config c, int plies = 0)
{
    switch (c)
    {
    case Config::A: return 0;
    case Config::T1:
    case Config::T1_A:
    case Config::A_T1: return 1;
    case Config::T2: return 2;
    case Config::T4: return 4;
    case Config::TK: return plies;
    }
    return 0;
}

inline bool config_has_fabric(Config c) { return c == Config::TK || fabric_plies(c) > 0; }
inline bool config_has_plate(Config c)
{
    return c != Config::T1 && c != Config::T2 && c != Config::T4 && c != Config::TK;
}

enum class Level
{
    Smoke, // ctest, секунды (ТЗ Guo 2019 §7)
    Base,  // геометрия статьи, guo2019_shot
    Fine,  // сходимость по разрешению
    Scaled // уменьшенная постановка (решение владельца 2026-10-01, отступление от ТЗ)
};

inline const char* level_name(Level l)
{
    switch (l)
    {
    case Level::Smoke: return "smoke";
    case Level::Base: return "base";
    case Level::Fine: return "fine";
    case Level::Scaled: return "scaled";
    }
    return "?";
}

inline bool parse_level(const std::string& s, Level& out)
{
    if (s == "smoke") out = Level::Smoke;
    else if (s == "base") out = Level::Base;
    else if (s == "fine") out = Level::Fine;
    else if (s == "scaled") out = Level::Scaled;
    else return false;
    return true;
}

// ===========================================================================
// Ударник: стальной шар (ТЗ Guo 2019 §3.2)
// ===========================================================================
inline constexpr double kBallD = 12.0e-3;     // 12 мм, статья
inline constexpr double kBallMass = 7.0e-3;   // 7.000 г, статья
inline constexpr double kBallGap0 = 0.5e-3;   // g_0 = 0.5 мм между R_eff и лицевой поверхностью мишени, оценка
inline constexpr double kRhoSteelRef = 7736.0; // кг/м^3 сплошного шара (7 г / (pi/6 12^3 мм^3)); только для печати
inline constexpr int kBallItype = 3;          // слот тантала, константы стали задаёт apply_guo_materials, решение §13 п.4
inline constexpr int kBallBody = 1;           // тело 1
// Smoke: шар уменьшен до D_s = 4 мм, масса по объёму: m = m_s (4/12)^3, ТЗ Guo 2019 §7.
inline constexpr double kBallDSmoke = 4.0e-3; // 4 мм, оценка (smoke)

// Шаг решётки шара по уровням, ТЗ Guo 2019 §7 (все -- оценки разрешения).
inline constexpr double kBallDxSmoke = 0.8e-3; // 0.8 мм
inline constexpr double kBallDxBase = 0.6e-3;  // 0.6 мм (около 4.2 тыс. частиц)
inline constexpr double kBallDxFine = 0.4e-3;  // 0.4 мм (около 14 тыс. частиц)

// Уровень scaled: уменьшенная постановка, решение владельца 2026-10-01, ОТСТУПЛЕНИЕ ОТ ТЗ.
inline constexpr double kBallDScaled = 6.0e-3;      // D_s = 6 мм (было 12 мм, масштаб 1/2), решение владельца 2026-10-01, отступление от ТЗ
inline constexpr double kBallDxScaled = 0.3e-3;     // 0.3 мм = D_s/20 (было 0.6 мм у base, масштаб 1/2; та же относительная
                                                    // точность), решение владельца 2026-10-01, отступление от ТЗ
inline constexpr double kBallRhoFactorScaled = 2.0; // плотность шара x2 (было x1), решение владельца 2026-10-01, отступление от ТЗ
// m_s = 2 * 7 г * (1/2)^3 = 1.75 г (было 7 г), решение владельца 2026-10-01, отступление от ТЗ
inline constexpr double kBallMassScaled = kBallRhoFactorScaled * kBallMass * (kBallDScaled / kBallD) * (kBallDScaled / kBallD) * (kBallDScaled / kBallD);

// Допуск отбора узлов решётки шара по радиусу (в долях dx): |r| <= D/2 - dx/2 + tol*dx,
// чтобы узлы, лежащие на границе в точной арифметике, отбирались одинаково на любой платформе.
inline constexpr double kBallRadiusTolDx = 1.0e-12;

// ===========================================================================
// Материалы: переопределения cfg.mats (ТЗ Guo 2019 §3.2, §3.4)
// ===========================================================================
// Слот массивов Materials = itype - Materials::n_Me (n_Me = 3).
inline constexpr std::size_t kSlotSteel = static_cast<std::size_t>(kBallItype) - Materials::n_Me; // 0
inline constexpr int kPlateItype = 4; // алюминий, реестр
inline constexpr std::size_t kSlotAl = static_cast<std::size_t>(kPlateItype) - Materials::n_Me;   // 1

inline constexpr double kSteelG = 80.0e9;     // G[0] = 80 ГПа, оценка (справочная для стали)
inline constexpr double kSteelSigy = 2.0e9;   // sigy[0] = 2.0 ГПа, оценка: закалённый шарикоподшипниковый шар
inline constexpr double kSteelZz = 4.0;       // zz[0] = 4.0, оценка
inline constexpr double kSteelAa = 4.0e10;    // aa[0] = 4.0e10 Па (K = aa zz = 160 ГПа), оценка
inline constexpr double kSteelGamm = 2.0;     // gamm[0] = 2.0, оценка
inline constexpr double kAlSigy = 300.0e6;    // sigy[1] = 300 МПа вместо 240 МПа реестра, статья (sigma_y ~ 300 МПа)
// G[1] = 26 ГПа из реестра не трогается (статья: E = 67 ГПа, nu = 0.32 -> G = 25.4 ГПа).

// Ожидаемые упругие константы шара из elastic_by_itype(cfg.mats, 3): E ~ 206 ГПа, nu ~ 0.29 (только для печати).
inline constexpr double kSteelERef = 206.0e9;
inline constexpr double kSteelNuRef = 0.29;

// Слот 3 в сценах Guo 2019 означает СТАЛЬ (в реестре materials.hpp -- тантал).
// Слот 4 (Al): меняется только предел текучести. materials.hpp не меняется.
inline void apply_guo_materials(Materials& m)
{
    m.G[kSlotSteel] = kSteelG;
    m.sigy[kSlotSteel] = kSteelSigy;
    m.zz[kSlotSteel] = kSteelZz;
    m.aa[kSlotSteel] = kSteelAa;
    m.gamm[kSlotSteel] = kSteelGamm;
    m.sigy[kSlotAl] = kAlSigy;
}

// ===========================================================================
// Ткань: 1 слой Twaron CT 717 (ТЗ Guo 2019 §3.3)
// ===========================================================================
// d = 0.39 мм (оценка: круг площадью 0.12 мм^2) и p = 1.18 мм (оценка: 8.5 нит/см) --
// значения пресета twaron_ct717 генератора; сцена их НЕ дублирует, а берёт из пресета.
// Здесь они -- справочные (для печати и сравнения с пресетом).
inline constexpr double kFabricDPreset = 0.39e-3;     // 0.39 мм, оценка
inline constexpr double kFabricPitchPreset = 1.18e-3; // 1.18 мм, оценка
inline constexpr double kFabricAmpOverD = 0.6;        // A = 0.6 d = 0.234 мм на всех уровнях, решение §2 п.1 (при d/2 -- 4 стартовых контакта на пересечение)
inline constexpr double kFabricE = 85.0e9;            // E = 85 ГПа, оценка: середина динамического диапазона 78-97 ГПа
inline constexpr double kFabricNu = 0.3;              // оценка
inline constexpr double kFabricRho = 1440.0;          // 1440 кг/м^3, оценка: справочная плотность арамида
inline constexpr double kFabricEpsBreak = 0.03;       // eps_b = 3 %, оценка: eps_b ~ 2 U_T / sigma_T
inline constexpr double kMuDefault = 0.3;             // mu_ff = mu_fs = 0.3, оценка (в статье трения нет; проверено на демо проекта, решение §13 п.8)
inline constexpr double kFabricOverhangYBase = 5.0e-3; // e_y = c = 5 мм (ширина зажима), оценка; выступ вдоль x (свободные нити) = p/2
inline constexpr double kFabricClampEdgeBase = 60.0e-3; // c_e = 60 мм: бусины зажатых нитей (вдоль y) с |y| >= 60 мм неподвижны, статья (верх и низ зажаты)
inline constexpr int kFabricNBase = 102;              // nx = ny = 102 нити (окно 120 мм), пресет
inline constexpr int kFabricPatchBase = 34;           // патч 34 x 34 нити для «одеяла» (plain_weave_quilt.hpp); сцену не меняет
inline constexpr int kFabricTileBase = 3;             // сетка патчей 3 x 3
static_assert(kFabricPatchBase * kFabricTileBase == kFabricNBase, "патч x сетка должны давать kFabricNBase нитей");
inline constexpr int kFabricSegBase = 2;             // k = 2 сегмента на полуволну, base, ТЗ Guo 2019 §7
inline constexpr int kFabricSegFine = 4;              // k = 4, fine
// smoke: 12 x 12 нитей с d, p пресета, k = 2, A = 0.6 d.
inline constexpr int kFabricNSmoke = 12;
inline constexpr int kFabricSegSmoke = 2;
inline constexpr double kFabricOverhangYSmoke = 2.0e-3;  // e_y = 2 мм, ТЗ Guo 2019 §7
inline constexpr double kFabricClampEdgeSmoke = 7.0e-3;  // c_e = 7 мм: зажим ткани совпадает с зажимом пластины |y| >= 7 мм;
                                                         // допустимый диапазон генератора 6.69-8.49 мм
// Поверхностная плотность ткани: статья 0.28 кг/м^2; ожидание модели по коду генератора --
// 0.317 (k = 2) и 0.320 (k = 4) кг/м^2 (только для печати, не проверка).
inline constexpr double kArealDensityPaper = 0.28;
inline constexpr double kArealDensityExpectK2 = 0.317;
inline constexpr double kArealDensityExpectK4 = 0.320;

// ===========================================================================
// Пластина: Al 6061-T6 (ТЗ Guo 2019 §3.4)
// ===========================================================================
inline constexpr double kPlateT = 1.6e-3;   // t = 1.6 мм, статья
inline constexpr double kPlateRho = 2700.0; // 2700 кг/м^3, справочная для 6061 (согласуется с «+6.5 %» статьи)
inline constexpr int kPlateBody = 2;        // тело 2
// base/fine: x in [-60; 60] мм (свободные кромки), y in [-65; 65] мм (по 5 мм в зажим), зажим |y| >= 60 мм.
inline constexpr double kPlateHalfXBase = 60.0e-3;       // 60 мм, статья (окно 120 мм)
inline constexpr double kPlateHalfYBase = 65.0e-3;       // 65 мм, оценка (ширина зажима 5 мм не задана)
inline constexpr double kPlateClampEdgeBase = 60.0e-3;   // 60 мм, статья
// smoke: окно 14 x 14 мм, зажим 2 мм: x in [-7; 7], y in [-9; 9], зажим |y| >= 7 мм (ТЗ Guo 2019 §7).
inline constexpr double kPlateHalfXSmoke = 7.0e-3;       // 7 мм, ТЗ §7
inline constexpr double kPlateHalfYSmoke = 9.0e-3;       // 9 мм, ТЗ §7
inline constexpr double kPlateClampEdgeSmoke = 7.0e-3;   // 7 мм, ТЗ §7
// Слоёв по толщине n_z (dx_p = t/n_z): smoke 2 (0.8 мм), base 3 (0.533 мм), fine 4 (0.4 мм); оценки разрешения.
inline constexpr int kPlateNzSmoke = 2;
inline constexpr int kPlateNzBase = 3;
inline constexpr int kPlateNzFine = 4;
// scaled: решение владельца 2026-10-01, отступление от ТЗ.
inline constexpr double kPlateTScaled = 0.8e-3;           // t = 0.8 мм (было 1.6 мм, масштаб 1/2; t/D = 0.133 как в статье)
inline constexpr int kPlateNzScaled = 3;                  // n_z = 3 (dx_p = 0.2667 мм; те же слои, что у base)
inline constexpr double kPlateHalfXScaled = 15.0e-3;      // 15 мм (было 60 мм, масштаб 1/4)
inline constexpr double kPlateClampWidthScaled = 1.25e-3; // ширина зажима 1.25 мм (было 5 мм, масштаб 1/4)
inline constexpr double kPlateClampEdgeScaled = 15.0e-3;  // зажим |y| >= 15 мм (было 60 мм, масштаб 1/4)
inline constexpr double kPlateHalfYScaled = kPlateClampEdgeScaled + kPlateClampWidthScaled; // 16.25 мм (было 65 мм, масштаб 1/4)
inline constexpr int kFabricNScaled = 26;                 // nx = ny = 26 (было 102): окно между крайними нитями 25 x 1.18 = 29.5 мм (масштаб ~1/4);
                                                          // ткань twaron НЕ масштабируется (d, p, A, материал как в пресете)
inline constexpr int kFabricSegScaled = 2;                // k = 2, как у base
inline constexpr double kFabricOverhangYScaled = 1.25e-3; // e_y = 1.25 мм (было 5 мм, масштаб 1/4), решение владельца 2026-10-01
inline constexpr double kFabricClampEdgeScaled = 15.0e-3; // c_e = 15 мм (было 60 мм, масштаб 1/4); диапазон генератора [14.945; 16.0] мм
inline constexpr double kTMaxScaleScaled = 0.5;           // t_max x0.5 (A 300 мкс, 1T 500, 1T_A и A_1T 400), оценка: времена процесса
                                                          // ~ размеру шара при той же скорости; решение владельца 2026-10-01
inline constexpr double kWpRadiusDefault = 30.0e-3;       // круг для w_p: r <= 30 мм, ТЗ Guo 2019 §6.1
inline constexpr double kWpRadiusScaled = 7.5e-3;         // 7.5 мм (было 30 мм, масштаб 1/4), решение владельца 2026-10-01

// Расхождение шагов по осям больше этого -- находка, её нужно описать (ТЗ Guo 2019 §3.4).
inline constexpr double kPlateDxMismatchWarn = 0.02; // 2 %

// ===========================================================================
// Укладка (ТЗ Guo 2019 §3.5)
// ===========================================================================
inline constexpr double kEpsG = 1.0e-6; // eps_g = 1 мкм: добавочный зазор ткань--пластина (страховка от округления), оценка
inline constexpr double kZFront = 0.0;  // z_f: лицевая поверхность мишени
// Зазор между соседними слоями пакета ткани (2T, 4T): отдельное имя, чтобы менять независимо от kEpsG;
// значение то же, 1 мкм, оценка (решение 2T/4T §12 п.2).
inline constexpr double kEpsPly = kEpsG;

// ===========================================================================
// Описатель уровня
// ===========================================================================
struct LevelSpec
{
    double dx_s;                // шаг решётки шара, м
    double D_s;                 // диаметр шара, м
    double m_s;                 // масса шара, кг
    int n_z;                    // слоёв пластины по толщине
    double plate_half_x;        // полуразмер пластины по x (свободные кромки), м
    double plate_half_y;        // полуразмер по y (включая зажим), м
    double plate_clamp_edge;    // зажаты частицы пластины с |y| >= plate_clamp_edge, м
    int nx, ny;                 // число нитей ткани: основа (вдоль x), уток (вдоль y)
    int seg;                    // k: сегментов нити на полуволну
    double e_y;                 // выступ нитей вдоль y за крайнюю нить, м
    double c_e;                 // граница зажима ткани |y| >= c_e, м
    double t;                   // толщина пластины, м
    double t_max_scale;         // множитель к default_t_max(config)
    double wp_radius;           // радиус круга для w_p, м
    double ball_rho_factor;     // плотность шара относительно сплошной стали (только печать: m_s её уже учитывает)
    double g_0;                 // зазор шар--мишень, м
};

inline LevelSpec level_spec(Level l)
{
    LevelSpec s{};
    s.t = kPlateT;
    s.t_max_scale = 1.0;
    s.wp_radius = kWpRadiusDefault;
    s.ball_rho_factor = 1.0;
    s.g_0 = kBallGap0;
    switch (l)
    {
    case Level::Smoke:
        s.dx_s = kBallDxSmoke;
        s.D_s = kBallDSmoke;
        // масса по объёму: m_s (4/12)^3
        s.m_s = kBallMass * (kBallDSmoke / kBallD) * (kBallDSmoke / kBallD) * (kBallDSmoke / kBallD);
        s.n_z = kPlateNzSmoke;
        s.plate_half_x = kPlateHalfXSmoke;
        s.plate_half_y = kPlateHalfYSmoke;
        s.plate_clamp_edge = kPlateClampEdgeSmoke;
        s.nx = s.ny = kFabricNSmoke;
        s.seg = kFabricSegSmoke;
        s.e_y = kFabricOverhangYSmoke;
        s.c_e = kFabricClampEdgeSmoke;
        break;
    case Level::Base:
    case Level::Fine:
        s.dx_s = (l == Level::Base) ? kBallDxBase : kBallDxFine;
        s.D_s = kBallD;
        s.m_s = kBallMass;
        s.n_z = (l == Level::Base) ? kPlateNzBase : kPlateNzFine;
        s.plate_half_x = kPlateHalfXBase;
        s.plate_half_y = kPlateHalfYBase;
        s.plate_clamp_edge = kPlateClampEdgeBase;
        s.nx = s.ny = kFabricNBase;
        s.seg = (l == Level::Base) ? kFabricSegBase : kFabricSegFine;
        s.e_y = kFabricOverhangYBase;
        s.c_e = kFabricClampEdgeBase;
        break;
    case Level::Scaled:
        s.dx_s = kBallDxScaled;
        s.D_s = kBallDScaled;
        s.m_s = kBallMassScaled;
        s.ball_rho_factor = kBallRhoFactorScaled;
        s.t = kPlateTScaled;
        s.n_z = kPlateNzScaled;
        s.plate_half_x = kPlateHalfXScaled;
        s.plate_half_y = kPlateHalfYScaled;
        s.plate_clamp_edge = kPlateClampEdgeScaled;
        s.nx = s.ny = kFabricNScaled;
        s.seg = kFabricSegScaled;
        s.e_y = kFabricOverhangYScaled;
        s.c_e = kFabricClampEdgeScaled;
        s.t_max_scale = kTMaxScaleScaled;
        s.wp_radius = kWpRadiusScaled;
        break;
    }
    return s;
}

// ===========================================================================
// Время расчёта и кадры
// ===========================================================================
// t_max по умолчанию для серий (ТЗ Guo 2019 §8.1, все оценки): потолок времени до исхода.
inline double default_t_max(Config c)
{
    switch (c)
    {
    case Config::A: return 600.0e-6;    // 600 мкс
    case Config::T1: return 1000.0e-6;  // 1000 мкс
    case Config::T1_A:
    case Config::A_T1: return 800.0e-6; // 800 мкс
    case Config::T2: return 800.0e-6;   // 800 мкс, оценка (2T 4T §3: V_BL в 1.9 раза выше, чем у 1T)
    case Config::T4: return 600.0e-6;   // 600 мкс, оценка (в Fig. 9 4T при 269 м/с прорван примерно к 200 мкс)
    case Config::TK: return 600.0e-6;   // 600 мкс, оценка (решение владельца 2026-11-09; по аналогии с 4T)
    }
    return 800.0e-6;
}

// t_max с учётом уровня: scaled -- в t_max_scale раз короче.
inline double default_t_max(Config c, Level l)
{
    return default_t_max(c) * level_spec(l).t_max_scale;
}

// Smoke (ТЗ Guo 2019 §9): V_0 = 100 м/с, t_total = 20 мкс, кадр раз в 2 мкс.
inline constexpr double kSmokeV0 = 100.0;        // м/с
inline constexpr double kSmokeTTotal = 20.0e-6;  // 20 мкс
inline constexpr double kSmokeSaveT = 2.0e-6;    // 2 мкс
// Серии на base/fine: шаг камеры статьи 40 000 кадр/с = 25 мкс (ТЗ Guo 2019 §7).
inline constexpr double kBaseSaveT = 25.0e-6;

// ===========================================================================
// Критерий исхода выстрела (ТЗ Guo 2019 §6.2)
// ===========================================================================
inline constexpr double kOutcomeWindow = 20.0e-6; // Delta t_ff = 20 мкс: окно свободного полёта / непрерывного v_z >= 0, оценка
inline constexpr double kOutcomeDv = 0.1;         // порог |Delta v_z| на окне, м/с, оценка

// ===========================================================================
// Выстрелы статьи (ТЗ Guo 2019 §8.2). V_0 сняты с рисунков (рис.), +-1-2 м/с.
// ===========================================================================
struct PaperShot
{
    double V0;       // м/с
    bool perforated; // п (пробит) / н (не пробит)
    double V_R;      // остаточная скорость опыта, м/с (0 для «н»)
};

inline std::vector<PaperShot> paper_shots(Config c)
{
    switch (c)
    {
    case Config::A:
        return {{121, false, 0}, {126, false, 0}, {130, false, 0}, {134, true, 40}, {135, true, 44}};
    case Config::T1:
        return {{59, false, 0}, {64, false, 0}, {70, false, 0}, {78, true, 15}, {78, true, 20}};
    case Config::T1_A:
        return {{149, false, 0}, {161, false, 0}, {175, false, 0}, {182, true, 45},
                {183, true, 67}, {184, false, 0}, {192, true, 60}};
    case Config::A_T1:
        return {{120, false, 0}, {145, false, 0}, {151, true, 29}, {152, true, 25}, {171, true, 72}};
    // 2T и 4T: рис. 10, векторная разметка PDF, +-0.3 м/с (ТЗ 2T/4T §3.1)
    case Config::T2:
        return {{137.6, false, 0}, {143.7, false, 0}, {147.7, true, 29.0}, {148.6, true, 10.2}, {153.5, true, 41.6}};
    case Config::T4:
        return {{247.4, false, 0}, {252.4, false, 0}, {257.7, false, 0}, {264.5, false, 0},
                {265.8, true, 83.0}, {267.2, true, 27.7}, {268.7, true, 86.2}};
    case Config::TK: return {}; // в статье нет
    }
    return {};
}

// Баллистический предел статьи V_BL, м/с, с отклонениями «+plus / -minus» (статья).
struct PaperVBL
{
    double value, plus, minus;
};

inline PaperVBL paper_vbl(Config c)
{
    switch (c)
    {
    case Config::A: return {128.0, 4.0, 1.0};
    case Config::T1: return {76.0, 3.0, 3.0};
    case Config::T1_A: return {178.0, 6.0, 8.0};
    case Config::A_T1: return {148.0, 3.0, 3.0};
    case Config::T2: return {147.0, 3.0, 3.0};
    case Config::T4: return {258.0, 8.0, 5.0};
    case Config::TK: return {0, 0, 0}; // в статье нет
    }
    return {0, 0, 0};
}

// Поглощённая энергия при V_BL статьи, Дж (статья, §4.2.2): 1T 20.2, 2T 75.6, 4T 233.0; для остальных 0 (не задана).
inline double paper_e_bl(Config c)
{
    switch (c)
    {
    case Config::T1: return 20.2;
    case Config::T2: return 75.6;
    case Config::T4: return 233.0;
    default: return 0.0;
    }
}

// Нормировки на 1T и механизм для 2T, 4T (статья §4.2.2-§4.2.3, Fig. 12-17); только для печати, не проверки.
inline constexpr double kPaperVblRatio2T = 1.9;    // V_BL(2T)/V_BL(1T), статья
inline constexpr double kPaperVblRatio4T = 3.4;    // V_BL(4T)/V_BL(1T), статья
inline constexpr double kPaperERatio2T = 3.7;      // e(2T)/e(1T) при V_BL, статья
inline constexpr double kPaperERatio4T = 11.5;     // e(4T)/e(1T) при V_BL, статья
inline constexpr double kPaperSpecific2T = 1.87;   // (e/n)(2T)/e(1T), статья, Fig. 13
inline constexpr double kPaperSpecific4T = 2.88;   // (e/n)(4T)/e(1T), статья, Fig. 13
inline constexpr double kPaperVblExp = 0.887;      // V_BL ~ n^0.887, подгонка статьи, Fig. 12a
inline constexpr double kPaperEExp = 1.773;        // e ~ n^1.773, подгонка статьи, Fig. 12b
// Порванные вертикальные (зажатые) нити по слоям, слой 1 -- лицевой, статья Fig. 16.
inline constexpr int kPaperBrokenV4T[4] = {5, 6, 7, 18};
inline constexpr int kPaperBrokenV1T[1] = {3};
inline constexpr int kPaperUnravel1T = 3;          // распущенные вертикальные нити у свободных кромок, 1T, статья Fig. 17
inline constexpr int kPaperUnravel4TPly1 = 9;      // то же, слой 1 у 4T, статья Fig. 17
inline constexpr double kPaperDeflMinM = 25.0e-3;  // максимальный прогиб перед разрушением 25-27 мм (1T, 2T, 4T), статья Fig. 14
inline constexpr double kPaperDeflMaxM = 27.0e-3;

// Нормировка e/(e_A + e_1T), статья: 1.43 для 1T/A, 0.99 для A/1T; для A и 1T -- 0 (не определена).
inline constexpr double kPaperNormT1_A = 1.43;
inline constexpr double kPaperNormA_T1 = 0.99;
inline double paper_norm(Config c)
{
    return c == Config::T1_A ? kPaperNormT1_A : (c == Config::A_T1 ? kPaperNormA_T1 : 0.0);
}

// ===========================================================================
// Допуски приёмки (ТЗ Guo 2019 §9)
// ===========================================================================
inline constexpr double kTolGap = 1.0e-12;        // S2: зазор шар--мишень = g_0 +- 1e-12 м
inline constexpr double kTolMassRel = 1.0e-12;    // S1: масса шара равна заданной, отн.
inline constexpr double kEnergyRatioFail = 1.02;  // S12: E_tot_max/E_tot_0 <= 1.02 (падение при превышении)
inline constexpr double kEnergyRatioWarn = 1.005; // порог WARN, как energy_balance.py проекта
inline constexpr double kDeltaOverDWarn = 0.3;    // WARN при max delta/d > 0.3 (friction-report.md §10 п.9)
inline constexpr double kDtRatioMax = 100.0;      // S5: dt_first/dt_min <= 100
inline constexpr double kStepsFactor = 10.0;      // S5: шагов меньше, чем в 10 раз от оценки по первому dt
inline constexpr double kMinDecel = 0.5;          // S8: |v_z(t_end)| < V_0 - 0.5 м/с

} // namespace guo2019
