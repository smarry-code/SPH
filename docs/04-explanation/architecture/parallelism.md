# Параллелизм: OpenMP, `CpuExecutor` и reduce по шардам

> Схема: `CpuExecutor` + `ContributionWriter`/буферы вкладов + reduce по шардам
> (owner-computes), [ADR-008](../decisions/adr-008-rhs-architecture-fiber.md). Ключ
> `coloring_niac_threshold` читается из конфигурации и пишется в чекпоинт, но
> **ни на что не влияет**.

## Проблема конкурирующей записи

В парных циклах SPH несколько пар могут писать в одну частицу: пара
$k_1=(3,7)$ и пара $k_2=(3,9)$ обе дают вклад в ускорение частицы 3. Одновременный
`+=` из разных потоков — гонка данных. В прежнем коде это решалось
`#pragma omp atomic update`, затем раскраской графа пар. Для нескольких
compute-структур с разнородными взаимодействиями (Solid-SPH, Solid–Solid-контакт,
Fiber-Bond, Fiber–Fiber, Fiber–Solid; примитив может объединять до четырёх
частиц) раскраска неприменима, поэтому принят другой подход.

## Решение: вклады и reduce по шардам

Взаимодействие разделено на стадии (см. [Обзор конвейера](pipeline-overview.md)):

```mermaid
flowchart LR
    A["Поиск<br/>(readonly)"] --> B["Примитивы<br/>(readonly)"]
    B --> C["Ядро: вклады<br/>ContributionWriter"]
    C --> D["Reduce по шардам<br/>(owner-computes)"]
    D --> E["Finalize<br/>parallel_for по частицам"]
```

- **Ядро взаимодействия** не пишет в общий RHS. Оно формирует типизированные
  вклады `Contribution<Delta>` (цель, `Delta`, метаданные `interaction_type`,
  `primitive_id`, `slot`) в буфер своего блока/потока:
  [`ContributionBuffer`](../../../include/rhs/contributions.hpp:56) (число вкладов
  заранее неизвестно) или
  [`IndexedContributions`](../../../include/rhs/contributions.hpp:93)
  (число вкладов известно заранее, позиции — префиксная сумма; так работает
  Solid-SPH). Вклады пишутся через
  [`ContributionWriter`](../../../include/rhs/contributions.hpp:149)
  (`emit(target, delta, primitive_id[, slot])`). `Delta` — `SolidDelta`
  ([`solid_sph.hpp`](../../../include/interactions/solid_sph.hpp:99)) или
  `fiber::FiberDelta`.
- **Reduce** — [`rhs::reduce_contributions_into`](../../../include/rhs/contributions.hpp:353):
  диапазон частиц $[0, N)$ режется на непрерывные шарды (по умолчанию
  $4\times$ число потоков, не больше $N$); ссылки на вклады раскладываются по
  шардам схемой count → prefix sum → scatter (сами вклады не копируются), затем
  каждый шард суммирует свой диапазон независимо от остальных. Шард —
  единственный владелец своих частиц, поэтому `atomic`/`mutex` не нужны.
  Можно свёртывать логически склеенные наборы (например, Solid-SPH +
  контакт Solid–Solid в `single_step`).
- **Finalize** — отдельный `parallel_for` по частицам вне ядра.

Физический код (`interactions/*.hpp`) не содержит `#pragma omp`, `mutex` и
`atomic`: параллелизм только через `CpuExecutor`.

## Детерминизм

Результат не зависит от числа потоков:

- параллельные запросы `SpatialIndex` и ядра режут диапазон на фиксированное число
  отрезков `CpuExecutor::kDetSegments = 64`
  ([`parallel_for_segments`](../../../include/rhs/cpu_executor.hpp)); отрезок
  целиком обрабатывает один поток по возрастанию индекса, вклады кладутся в буфер
  с номером **отрезка**, а не потока;
- порядок суммирования на каждую частицу в reduce: шард → номер буфера →
  позиция в буфере; поэтому суммирование совпадает с последовательным при
  любом числе потоков;
- `parallel_reduce` объединяет частичные результаты тех же отрезков по порядку,
  поэтому сумма `double` воспроизводима;
- `ReduceOptions::deterministic` дополнительно сортирует вклады шарда по
  `(target, interaction_type, primitive_id, slot)` (по умолчанию `false`).

## `CpuExecutor`

[`rhs::CpuExecutor`](../../../include/rhs/cpu_executor.hpp:41) — единая абстракция
параллельности (один объект на вызов `single_step` и один `fiber_exec` на цикл),
реализация на OpenMP:

| Примитив | Назначение |
|---|---|
| `parallel_for(n, grain, f)` | независимые итерации, `schedule(dynamic, grain)` |
| `parallel_for_segments(n, f)` | детерминированное разбиение на 64 отрезка, `f(seg, i)` |
| `parallel_reduce(n, init, map, combine)` | свёртка с фиксированным порядком объединения |
| `parallel_scan(in, out, n)` | исключающая префиксная сумма |
| `parallel_sort(first, last, less)` | блочная сортировка + `inplace_merge` |
| `make_task_group()` | `run(f)`/`wait()` на `std::async` (для независимых ветвей) |
| `num_workers()`, `worker_id()` | число потоков, номер текущего потока |

Без `-fopenmp` (нет макроса `_OPENMP`) executor работает последовательно.
`CpuExecutor(0)` берёт `omp_get_max_threads()`.

## Рабочая память

Буферы вкладов, массивы пар, индексы и `ReduceScratch` живут между шагами
(`SolidStepWorkspace` в `WorkBuffers::solid_ws`, буферы нитей в
`run_time_integration`) и **только растут**: на шагах после первого новых
аллокаций почти нет. Общий именованный реестр буферов —
[`RhsWorkspace`](../../../include/rhs/rhs_workspace.hpp) (`get<T>(key)`).

## Настройка числа потоков

Число потоков OpenMP задаётся в
[`SimulationConfig::num_threads`](../../../include/simulation_config.hpp:93)
(по умолчанию 7; ключа в `Initial_Conditions.txt` нет, значение меняется только в коде; в чекпоинте сверяется) и устанавливается в [`run_sph()`](../../../src/main.cpp:67) до чтения конфигурации через
`omp_set_num_threads()`; `CpuExecutor` подхватывает это значение через
`omp_get_max_threads()`.

## Профилирование

CMake-опция `SPH_PROFILE` (по умолчанию `OFF`) сохранена, но блоков
`#ifdef SPH_PROFILE` в коде сейчас нет; `force_timing.csv` не пишется.
Время всего расчёта печатает `main()`. Лог энергии (`SPH_ENERGY_LOG`) и
диагностика контактов нитей (`SPH_CONTACT_STATS`) к параллелизму отношения не имеют.

## Тестирование

- `executor_test` — [`tests/test_executor.cpp`](../../../tests/test_executor.cpp)
  (`for`/`reduce`/`scan`/`sort`);
- `contributions_test` — [`tests/test_contributions.cpp`](../../../tests/test_contributions.cpp)
  (reduce, детерминизм при 1/2/7 потоках);
- `spatial_index_test` — [`tests/test_spatial_index.cpp`](../../../tests/test_spatial_index.cpp).

## Связанные материалы

- [ADR-008: Архитектура RHS и нити](../decisions/adr-008-rhs-architecture-fiber.md)
- [Поиск соседних частиц](neighbor-search.md)
- Файлы: [`cpu_executor.hpp`](../../../include/rhs/cpu_executor.hpp),
  [`contributions.hpp`](../../../include/rhs/contributions.hpp)
