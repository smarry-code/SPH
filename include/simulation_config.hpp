#pragma once
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

    // --- Контакт «частица-частица» (part_0.md, ADR-007) ---
    // 0 - выключен (по умолчанию), 1 - пинбол Белицко, 2 - Свигл, 3 - импульсный.
    int    contact_force_type   = 0;
    double contact_penalty_scale = 1.0;   // масштабный коэффициент Kp (типы 1-2)
    // Доля нормальной скорости сближения, снимаемая за шаг (только тип 3).
    // 1.0 — гасить сближение пары за шаг полностью (буквальная (2.10));
    // меньше — оставлять паре (1-gamma)*pdot. Рабочий диапазон (0; 1.6].
    double contact_impulse_gamma = 1.0;

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
    contact_force_type = 0;  contact_penalty_scale = 1.0;  contact_impulse_gamma = 1.0;
    CFL = 0.1;  t_total = 3e-6;  save_t = 2.5e-8;  print_step = 50;
    alpha = 1.0;  beta = 1.0;  etq = 0.01;  epsilon = 0.3;  xi = 0.001;  epsilon_h = 0.01;
    num_threads = 7;
    maxn = 0;  max_interaction = 0;  moni_particle = 0;
}
