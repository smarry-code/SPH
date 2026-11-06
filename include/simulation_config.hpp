#pragma once
#include <iostream>
#include <stdexcept>
#include <string>
// simulation_config.hpp
// Инкапсуляция runtime-параметров симуляции (подэтап 5.7 этапа 6.3
// «Архитектурная чистота»).
//
// Ранее параметры расчёта (pa_sph, nnps, skf, CFL, alpha/beta, ...) были
// неконстантными глобалами в param.hpp и пронизывали весь конвейер как
// неявное глобальное состояние. Теперь они собраны в один конфиг-объект:
//   - создаётся в main();
//   - наполняется в parse_sph_config()/input() (значения из Initial_Conditions.txt);
//   - передаётся по const SimulationConfig& во все узлы расчёта.
// Мутабельные поля (пишутся после старта):
//   - maxn / max_interaction — устанавливаются в main после input (ntotal*1.4);
//   - moni_particle — пишется в input (shock_tube/parse_k_file), читается в time_integration.
//
// ЧТО ОСТАЁТСЯ В param.hpp (намеренно, ADR-2 этапа 6.3):
//   - dim (constexpr) — определяет РАЗМЕР C-массивов ParticleSystem (x[dim], vx[dim]),
//     это compile-time характеристика структуры данных, а не runtime-параметр;
//   - логические флаги (visc, axis_sim, summation_density, ...) — compile-time
//     переключатели веток компиляции; не являются mutable состоянием.
// Эти constexpr-сущности корректны с точки зрения ODR и не блокируют 5.2-В.
#include <cstddef>
#include <string>
#include "materials.hpp"

struct SimulationConfig
{
    // --- SPH-метод и алгоритмы (загружаются из Initial_Conditions.txt) ---
    int pa_sph = 2;                     // Метод аппроксимации давления
    int nnps   = 1;                     // Метод поиска соседних частиц
    int sle    = 1;                     // Алгоритм обновления длины сглаживания
    int skf    = 1;                     // Тип ядра сглаживания
    int density_method          = 1;    // Метод расчёта плотности (1/2/3)
    int coloring_niac_threshold = 0;    // Порог graph coloring vs atomic

    // --- ODESolver (fiber/PLAN.md, этап Ф8; §3.9) ---
    // 0 - PredictorCorrector (по умолчанию, схема проекта, на неё
    //     откалиброваны эталоны BC2/BC6), 1 - явный Euler (эталон для тестов
    //     порядка сходимости, неустойчив на незатухающих упругих связях --
    //     не для рабочих расчётов).
    int ode_solver = 0;

    // --- Контакт «частица-частица» (part_0.md, ADR-007) ---
    // 0 - выключен (по умолчанию), 1 - пинбол Белицко, 2 - Свигл, 3 - импульсный.
    int    contact_force_type   = 0;
    double contact_penalty_scale = 1.0;   // масштабный коэффициент Kp (типы 1-2)
    // Доля нормальной скорости сближения, снимаемая за шаг (только тип 3).
    // 1.0 — гасить сближение пары за шаг полностью (буквальная (2.10));
    // меньше — оставлять паре (1-gamma)*pdot. Рабочий диапазон (0; 1.6].
    double contact_impulse_gamma = 1.0;

    // --- Закон нормальной силы контактов нитей (fiber/SPH.md, «Законы
    //     нормальной силы»; include/interactions/fiber_contact_law.hpp) ---
    // 0 = ГЕРЦ (прежнее поведение), НЕ «выключено» -- в отличие от
    // contact_force_type. 1 - пинбол, 2 - Свигл, 3 - импульсный, 4 - Герц + демпфер
    // (Kp -- масштаб демпфера).
    int    fiber_solid_contact_type          = 0;   // Fiber-Solid; 0 = Герц, не выключено
    int    fiber_fiber_contact_type          = 0;   // Fiber-Fiber; 0 = Герц, не выключено
    double fiber_solid_contact_penalty_scale = 1.0; // Kp (типы 1, 2)
    double fiber_solid_contact_impulse_gamma = 1.0; // gamma (тип 3), рабочий диапазон (0; 1.6]
    double fiber_fiber_contact_penalty_scale = 1.0; // Kp (типы 1, 2)
    double fiber_fiber_contact_impulse_gamma = 1.0; // gamma (тип 3), рабочий диапазон (0; 1.6]
    // Выталкивание типа 3 при глубоком проникании (fiber/plan-fiber-contact-laws12-fix.md;
    // FiberContactLaw::push_tau): к pdot добавляется max(p - depth*r, 0)/tau, tau в с.
    // tau = 0 -- выключено (побитово прежний тип 3); depth -- доля суммы радиусов пары.
    double fiber_solid_contact_impulse_push_tau = 0.0;
    double fiber_fiber_contact_impulse_push_tau = 0.0;
    double fiber_contact_impulse_push_depth = 0.1;
    double fiber_contact_impulse_push_vmax = 0.0; // предел скорости выталкивания, м/с; 0 -- без предела
    // Нормальная сила законов 1-3 как проекция скорости (fiber/TZ_LAW3_PROJECTION.md):
    // 0 -- прежняя явная сила (побитово прежнее поведение), 1 -- проекция.
    int    fiber_fiber_contact_predictor = 0; // Fiber-Fiber, 0/1
    int    fiber_solid_contact_predictor = 0; // Fiber-Solid, 0/1
    int    fiber_contact_projection_iters = 8; // число итераций проекции, 1..100
    int    fiber_contact_projection_warm  = 1; // тёплый старт импульсов, 0/1
    // Множитель общего шага по времени (проверка сходимости по шагу); 1 -- без изменений,
    // допустимо (0; 1].
    double dt_scale = 1.0;

    // --- Изоляция SPH-сумм Solid по телам ---
    // При true SPH-пары реальных частиц разных тел (sys.body) отвергаются;
    // тела взаимодействуют только через Solid-Solid-Contact
    // (contact_force_type != 0). По умолчанию false (прежнее поведение).
    bool sph_same_body_only = false;

    // --- Шаг по времени / вывод ---
    double CFL       = 0.1;             // Число Куранта
    double t_total   = 3e-6;            // Полное время расчёта
    double save_t    = 2.5e-8;          // Интервал сохранения кадров
    size_t print_step = 50;             // Шаг вывода прогресса на экран

    // --- Вязкость / диффузия / усреднение ---
    double alpha = 1.0, beta = 1.0, etq = 0.01;          // Псевдовязкость
    double epsilon = 0.3;                                 // Усреднение скоростей (Monaghan 1992)
    double xi        = 0.001;                             // Интенсивность диффузии плотности
    double epsilon_h = 0.01;                              // Регуляризатор знаменателя psi_ij

    // --- Параллелизм (N3: вместо хардкода omp_set_num_threads(7) в main) ---
    int num_threads = 7;

    // --- Runtime-размеры (НЕ constexpr; устанавливаются в main после input) ---
    size_t maxn = 0;
    size_t max_interaction = 0;

    // --- Мутабельное runtime-состояние ---
    size_t moni_particle = 0;           // Пишется в input, читается в time_integration

    // --- Материалы (константны после инициализации; ADR-3 этапа 6.3) ---
    Materials mats;

    // --- Нити (fiber/PLAN.md, этап Ф2) ---
    // Путь к файлу ломаных нитей (см. формат в fiber_input.cpp). Пусто (по
    // умолчанию) — нитей в расчёте нет. Если задан, а сборка сделана не с
    // SPH_DIM=3, parse_sph_config бросает std::runtime_error немедленно
    // (fiber/SPH.md: "Нити поддерживаются только в трёхмерной постановке").
    std::string fiber_file;

    // --- Импорт LS-DYNA .k (readk/k-import-plan.md, этап K6) ---
    // k_file пуст (по умолчанию) -- импорта нет. Непустой k_file имеет
    // приоритет над встроенными генераторами (shocktube/shearcavity).
    std::string k_file;
    std::string k_mode = "solid";           // "solid" | "solid_fiber"
    std::string k_parts_file;               // таблица деталей (Р2); пусто -- режим без таблицы
    double k_length_scale = 1.0;            // множитель координат (Р5), > 0
    double k_mass_scale   = 1.0;            // множитель масс (Р5), > 0
    long long k_monitor_nid = -1;           // NID мониторной частицы; < 0 -- не задан
    std::string k_solid_elements = "skip";  // "skip" | "particles" (*ELEMENT_SOLID)

    // --- Чекпоинты / рестарт (include/checkpoint.hpp) ---
    // Настройки запуска, не физика: в чекпоинт не сериализуются и не сравниваются.
    // Пустые/0 -- выключено, цикл побитово прежний.
    std::string checkpoint_file;            // путь чекпоинта; пусто -- не писать
    size_t checkpoint_every_steps = 0;      // период записи, в шагах интегрирования; 0 -- выкл.
    std::string restart_file;               // чекпоинт для продолжения; пусто -- старт с нуля
    bool restart_force = false;             // --force: разрешить менять физические настройки при рестарте

#ifdef SPH_ENERGY_LOG
    // --- Лог механической энергии (include/energy_log.hpp) ---
    // Только в сборке с -DSPH_ENERGY_LOG=ON. energy_log_every = k > 0 -- строка
    // лога каждые k шагов (после шага интегратора), 0 -- лог выключен.
    size_t energy_log_every = 0;
    std::string energy_log_file = "energy_log.csv";
#endif

#ifdef SPH_CONTACT_STATS
    // --- Диагностика контактов нитей (interactions/fiber_contact_stats.hpp) ---
    // Только в сборке с -DSPH_CONTACT_STATS=ON. contact_stats_every = k > 0 --
    // строка CSV каждые k шагов и сводка в конце; contact_hold_n -- контрфакт
    // «удержание N» (меняет физику, только для диагностики).
    size_t contact_stats_every = 0;
    std::string contact_stats_file = "contact_stats.csv";
    std::string contact_stats_summary = "contact_stats_summary.txt";
    bool contact_hold_n = false;
#endif

    // Установка значений по умолчанию (перенос значений из бывших param.hpp).
    // Вызывается из main ДО parse_sph_config: явно фиксирует контракт
    // инициализации конфига (независимо от default member initializers).
    // parse_sph_config переопределяет поля, явно заданные в Initial_Conditions.txt.
    void apply_defaults();
};

// inline-определение: безопасно при включении в несколько TU (5.2-В).
inline void SimulationConfig::apply_defaults()
{
    pa_sph = 2;  nnps = 1;  sle = 1;  skf = 1;
    density_method = 1;  coloring_niac_threshold = 0;
    ode_solver = 0;
    sph_same_body_only = false;
    contact_force_type = 0;  contact_penalty_scale = 1.0;  contact_impulse_gamma = 1.0;
    // Контакты нитей: 0 = Герц (не «выключено»).
    fiber_solid_contact_type = 0;  fiber_solid_contact_penalty_scale = 1.0;  fiber_solid_contact_impulse_gamma = 1.0;
    fiber_fiber_contact_type = 0;  fiber_fiber_contact_penalty_scale = 1.0;  fiber_fiber_contact_impulse_gamma = 1.0;
    fiber_solid_contact_impulse_push_tau = 0.0;  fiber_fiber_contact_impulse_push_tau = 0.0;
    fiber_contact_impulse_push_depth = 0.1;  fiber_contact_impulse_push_vmax = 0.0;  dt_scale = 1.0;
    fiber_fiber_contact_predictor = 0;  fiber_solid_contact_predictor = 0;
    fiber_contact_projection_iters = 8;  fiber_contact_projection_warm = 1;
    CFL = 0.1;  t_total = 3e-6;  save_t = 2.5e-8;  print_step = 50;
    alpha = 1.0;  beta = 1.0;  etq = 0.01;  epsilon = 0.3;  xi = 0.001;  epsilon_h = 0.01;
    num_threads = 7;
    maxn = 0;  max_interaction = 0;  moni_particle = 0;
    fiber_file.clear();
    k_file.clear();  k_mode = "solid";  k_parts_file.clear();
    k_length_scale = 1.0;  k_mass_scale = 1.0;  k_monitor_nid = -1;
    k_solid_elements = "skip";
#ifdef SPH_ENERGY_LOG
    energy_log_every = 0;  energy_log_file = "energy_log.csv";
#endif
}

// Проверка настроек проекции скорости нормального контакта (fiber/TZ_LAW3_PROJECTION.md).
// Бросает std::runtime_error (по-русски). Обе проекции выключены -- ничего не делает.
// Вызывается после разбора конфига/аргументов (input.cpp, guo2019_shot, time_integration).
inline void validate_fiber_contact_projection(const SimulationConfig& cfg)
{
    const int pf = cfg.fiber_fiber_contact_predictor, ps = cfg.fiber_solid_contact_predictor;
    if (pf == 0 && ps == 0) return;
    if (pf != 0 && pf != 1)
        throw std::runtime_error("fiber_fiber_contact_predictor: допустимы 0 и 1, получено " + std::to_string(pf));
    if (ps != 0 && ps != 1)
        throw std::runtime_error("fiber_solid_contact_predictor: допустимы 0 и 1, получено " + std::to_string(ps));
    if (cfg.fiber_contact_projection_iters < 1 || cfg.fiber_contact_projection_iters > 100)
        throw std::runtime_error("fiber_contact_projection_iters: допустимо 1..100, получено " +
                                 std::to_string(cfg.fiber_contact_projection_iters));
    if (cfg.fiber_contact_projection_warm != 0 && cfg.fiber_contact_projection_warm != 1)
        throw std::runtime_error("fiber_contact_projection_warm: допустимы 0 и 1, получено " +
                                 std::to_string(cfg.fiber_contact_projection_warm));
    bool gamma_warned = false;
    auto check_kind = [&](const char* name, int predictor, int type, double push_tau, double gamma)
    {
        if (predictor != 1) return;
        if (type < 1 || type > 3)
            throw std::runtime_error(std::string(name) + "_contact_predictor = 1 требует закона контакта 1, 2 или 3, а " +
                                     name + "_contact_type = " + std::to_string(type));
        if (push_tau > 0.0)
            throw std::runtime_error(std::string(name) + "_contact_impulse_push_tau > 0 запрещён в режиме проекции (" +
                                     name + "_contact_predictor = 1)");
        if (gamma != 1.0 && !gamma_warned)
        {
            gamma_warned = true;
            std::cerr << "ПРЕДУПРЕЖДЕНИЕ: " << name << "_contact_impulse_gamma = " << gamma
                      << " != 1: γ игнорируется в режиме проекции\n";
        }
    };
    check_kind("fiber_fiber", pf, cfg.fiber_fiber_contact_type, cfg.fiber_fiber_contact_impulse_push_tau,
               cfg.fiber_fiber_contact_impulse_gamma);
    check_kind("fiber_solid", ps, cfg.fiber_solid_contact_type, cfg.fiber_solid_contact_impulse_push_tau,
               cfg.fiber_solid_contact_impulse_gamma);
}
