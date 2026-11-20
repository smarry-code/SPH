# ADR-003: Инкапсуляция параметров в `SimulationConfig`

**Статус:** Принято
**Дата:** Подэтап 5.7 (этап 6.3 «Архитектурная чистота»)

## Контекст

Ранее параметры расчёта (`pa_sph`, `nnps`, `skf`, `CFL`, `alpha`/`beta`,
`t_total`, ...) были неконстантными глобалами в
[`param.hpp`](../../../include/param.hpp) и пронизывали весь конвейер как
неявное глобальное состояние. Это приводило к:

- скрытым зависимостям между модулями;
- невозможности изолированного тестирования;
- риску гонок при параллельной записи в глобалы;
- нарушению ODR при включении заголовков в несколько единиц трансляции.

## Решение

Параметры собраны в один конфиг-объект
[`SimulationConfig`](../../../include/simulation_config.hpp:31):

- создаётся в [`run_sph()`](../../../src/main.cpp:65) (вызывается из `main()`);
- наполняется значениями по умолчанию (`apply_defaults()`) и из
  [`Initial_Conditions.txt`](../../03-reference/config/initial-conditions.md)
  (`parse_sph_config()`);
- передаётся по `const SimulationConfig&` во все узлы расчёта.

Мутабельные поля (пишутся после старта) явно выделены: `maxn`/`max_interaction`
устанавливаются в `main` после `input()`; `moni_particle` пишется в `input` и
читается в `time_integration`.

## Последствия

**Положительные:**

- Устранено неконстантное глобальное состояние;
- Явный контракт инициализации конфига;
- Изолируемое тестирование модулей с собственным `SimulationConfig`.

**Отрицательные:**

- Увеличение числа параметров в сигнатурах функций;
- Часть полей остаётся мутабельной (`maxn`, `moni_particle`).

## Связанные файлы

- [`simulation_config.hpp`](../../../include/simulation_config.hpp:31): определение структуры.
- [`input.cpp`](../../../src/input.cpp:110): `parse_sph_config()`.

## Актуальное состояние (2026-11)

- Конфиг создаётся в `run_sph()` ([`main.cpp:65`](../../../src/main.cpp:65)); `main()`
  только вызывает её и перехватывает исключения импорта. `apply_defaults()` вызывается
  до `parse_sph_config()` ([`input.cpp:110`](../../../src/input.cpp:110)).
- Структура сильно выросла: помимо исходных полей в ней типы контактов и их параметры,
  `ode_solver`, ключи импорта `.k`, `fiber_file`, `jc_heating`/`jc_chi`, `xsph_mode`,
  `energy_log_*` и др. Полный перечень — в
  [справке по `SimulationConfig`](../../03-reference/config/simulation-config.md) и
  [`Initial_Conditions.txt`](../../03-reference/config/initial-conditions.md).
- Часть констант метода (`visc`, `axis_sim`, `summation_density` и др.) остаётся
  `constexpr` в [`param.hpp`](../../../include/param.hpp) (см. ADR-002).
