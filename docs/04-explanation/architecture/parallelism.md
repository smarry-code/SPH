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
  [`ContributionBuffer`](../../../include/rhs/contributions.hpp:60) (число вкладов
  заранее неизвестно) или
  [`IndexedContributions`](../../../include/rhs/contributions.hpp:97)
  (число вкладов известно заранее, позиции — префиксная сумма; так работает
  Solid-SPH). Вклады пишутся через
  [`ContributionWriter`](../../../include/rhs/contributions.hpp:153)
  (`emit(target, delta, primitive_id[, slot])`). `Delta` — `SolidDelta`
  ([`solid_sph.hpp`](../../../include/interactions/solid_sph.hpp:100)) или
  `fiber::FiberDelta`.
- **Reduce** — [`rhs::reduce_contributions_into`](../../../include/rhs/contributions.hpp:377):
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

**Контракт детерминизма** (шапка [`contributions.hpp`](../../../include/rhs/contributions.hpp)).
При `deterministic = false` (так зовут все ядра расчёта) сортировки нет: порядок суммирования на
частицу — номер буфера, затем позиция в буфере. Независимость от числа потоков держится на
производителе: буфер $g$ заполняет ровно один отрезок $g$ фиксированного разбиения, а внутри отрезка
вклады идут по возрастанию номера primitive. Ядро, которое пишет вклады иначе (`worker_id()`,
динамическая раздача), обязано включать `deterministic`. Следствие: обход **подмножества** primitive
(только активных слотов, только выживших после отсечки) в порядке возрастания номера с **прежним
номером отрезка** сохраняет порядок сложения на частицу и потому результат побитово. На этом
построены новые примитивы ниже и стадии контакта нитей (трение по активным слотам, проекция
нормального контакта).

## `CpuExecutor`

[`rhs::CpuExecutor`](../../../include/rhs/cpu_executor.hpp:41) — единая абстракция
параллельности (один объект на вызов `single_step` и один `fiber_exec` на цикл),
реализация на OpenMP:

| Примитив | Назначение |
|---|---|
| `parallel_for(n, grain, f)` | независимые итерации, `schedule(dynamic, grain)` |
| `parallel_for_segments(n, f)` | детерминированное разбиение на 64 отрезка, `f(seg, i)` |
| `parallel_reduce(n, init, map, combine)` | свёртка с фиксированным порядком объединения |
| `seg_begin(n, s)` | начало отрезка $s$: $\lfloor n/64\rfloor s + \lfloor (n \bmod 64)\,s/64\rfloor$, без переполнения $n\cdot s$ |
| `seg_of_index(k, n)` | номер отрезка, в котором индекс $k$ стоит при `parallel_for_segments(n)` |
| `parallel_for_seg_groups(m, seg_of, f)`, `parallel_reduce_seg_groups(...)` | обход/свёртка **подмножества** индексов (`j = 0..m-1` по возрастанию), у каждого известен номер исходного отрезка `seg_of(j)`; вклады ложатся в те же буферы и в том же порядке, что при обходе полного диапазона с пропуском невыбранных |
| `seg_group_bounds`, `parallel_for_seg_bounds(bound, f)`, `parallel_reduce_seg_bounds(...)` | то же по готовым границам `bound[0..64]` (границы считаются один раз и используются в нескольких проходах) |
| `SegPart<T>`, `prepare_segment_parts`, `concat_segments(parts, out)` | сбор результатов по отрезкам: отрезок $s$ пишет в `parts[s].v`, затем `concat_segments` склеивает по порядку отрезков (смещения — префиксная сумма, копирование параллельно); `SegPart` с `alignas(64)` |
| `parallel_scan(in, out, n)` | исключающая префиксная сумма |
| `parallel_sort(first, last, less)` | блочная сортировка + `inplace_merge` |
| `make_task_group()` | `run(f)`/`wait()` на `std::async` (для независимых ветвей) |
| `num_workers()`, `worker_id()` | число потоков, номер текущего потока |

Без `-fopenmp` (нет макроса `_OPENMP`) executor работает последовательно.
`CpuExecutor(0)` берёт `omp_get_max_threads()`.

Примитивы «по подмножеству» нужны, чтобы не обходить плотные массивы целиком: контакт нитей
(`cull_apart`, список активных слотов трения, сжатый массив нормальных слотов проекции) работает
только с нужными элементами, но суммирует в том же порядке, что полный обход (если невыбранные
элементы дают нейтральный элемент свёртки: $0.0$ в сумме, целый счётчик, `max` с нулём).

### Выравнивание `alignas(64)`

`ContributionBuffer` ([`contributions.hpp`](../../../include/rhs/contributions.hpp:60)) и
`CpuExecutor::SegPart` выровнены по 64 байтам: буферы лежат в `std::vector` по одному на отрезок,
соседние отрезки заполняют разные потоки, и без выравнивания заголовки `vector` (24 Б) делили бы
кэш-линию (ложное разделение при каждом `push_back`). На порядок вкладов и результат выравнивание
не влияет.

## Что выполняется параллельно

Помимо ядер и reduce, через `CpuExecutor` идут (результат побитово совпадает с последовательным):

- **Циклы решателя** ([`ode_solver.hpp`](../../../include/ode/ode_solver.hpp)): `predict`/`correct`/`advance`
  предиктор-корректора, шаг Эйлера, `pc_finalize_velocity_position*` (включая режим `xsph_mode = position`),
  `euler_step_position_velocity`, `euler_step_scalar` — `parallel_for` по частицам (каждый элемент пишется
  один раз тем же выражением).
- **`SpatialIndex::build(exec, ...)`**: копирование координат, bbox (`parallel_reduce`), ячейка каждой
  частицы, подсчёт по сегментам и раскладка — внутреннее состояние равно последовательному `build`
  побитово (см. [Поиск соседей](neighbor-search.md)). Вызывается из `build_solid_sph_search`,
  `build_fiber_contact_search`, `build_fiber_solid_contact_search`.
- **`build_solid_extended(exec, ...)`** ([`solid_sph.hpp`](../../../include/interactions/solid_sph.hpp:340)):
  копирование реальных и зеркальных частиц; максимум $h$ для радиуса поиска — `parallel_reduce` с `max`
  ([`solid_sph.hpp`](../../../include/interactions/solid_sph.hpp:532)).
- **`Diagnostics::validate(exec, ...)`** ([`diagnostics.hpp`](../../../include/diagnostics.hpp)):
  параллельная проверка «есть ли нечисло» (булево ИЛИ — точно при любом порядке); сообщения печатает прежний
  последовательный проход, только если нарушение найдено.
- **Контакт нитей**: `compute_fiber_bead_dt_ref`, `FiberDtRefCache`, фильтр `cull_apart`, трение по
  активным слотам, стадия проекции нормального контакта — см. [Контакт нитей](../models/05-fiber-contact.md).

**Остаётся последовательным:** накопление `u_clamp_acc` (энергия, внесённая отсечкой $u\ge0$, только при
`SPH_ENERGY_LOG`) — сумма `double`, порядок слагаемых должен быть фиксированным, поэтому в сборке с
`SPH_ENERGY_LOG` поэлементный цикл `correct_solid`/`step_solid` (скалярные поля Solid) идёт серийно,
а без неё — через `exec.parallel_for` (комментарий в `ode_solver.hpp`). Серийна и
проверка конечности плотности/массы в `compute_solid_derived`: исключение нельзя бросать из тела
`parallel for` OpenMP.

## Рабочая память

Буферы вкладов, массивы пар, индексы и `ReduceScratch` живут между шагами
(`SolidStepWorkspace` в `WorkBuffers::solid_ws`, буферы нитей в
`run_time_integration`) и **только растут**: на шагах после первого новых
аллокаций почти нет. Общий именованный реестр буферов —
[`RhsWorkspace`](../../../include/rhs/rhs_workspace.hpp) (`get<T>(key)`).

## Настройка числа потоков

Число потоков OpenMP задаётся в
[`SimulationConfig::num_threads`](../../../include/simulation_config.hpp:104)
(по умолчанию 7; ключа в `Initial_Conditions.txt` нет, значение меняется только в коде; в чекпоинте
свободно меняется при рестарте) и устанавливается в [`run_sph()`](../../../src/main.cpp:67) до чтения конфигурации через
`omp_set_num_threads()`; `CpuExecutor` подхватывает это значение через
`omp_get_max_threads()`. Конвейер сам `cfg.num_threads` не читает: утилиты (`guo2019_shot --threads=N`)
вызывают `omp_set_num_threads` напрямую. `parallel_for(n, f)` без `grain` берёт порцию $\max(1, n/64)$.

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
- `spatial_index_test` — [`tests/test_spatial_index.cpp`](../../../tests/test_spatial_index.cpp);
- `fiber_contact_projection_test` (метка `3d`) — проекция нормального контакта, в т. ч. побитовое
  совпадение при 1 и 8 потоках (Т7), сжатые слоты против полного массива при 1/2/7/8 (Т8); `checkpoint_test` — независимость результата рестарта от числа потоков.

## Связанные материалы

- [ADR-008: Архитектура RHS и нити](../decisions/adr-008-rhs-architecture-fiber.md)
- [Поиск соседних частиц](neighbor-search.md)
- Файлы: [`cpu_executor.hpp`](../../../include/rhs/cpu_executor.hpp),
  [`contributions.hpp`](../../../include/rhs/contributions.hpp)
