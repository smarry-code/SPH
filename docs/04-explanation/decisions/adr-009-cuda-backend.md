# ADR-009: CUDA-бэкенд

**Статус:** Принято
**Дата:** 2026-09-27

## Контекст

Проект переносит конвейер Solid + Fiber ([SPH.md](../../../fiber/SPH.md),
[IMPLEMENTATION_REPORT.md](../../../fiber/IMPLEMENTATION_REPORT.md)) на
GPU — план [`fiber/CUDA_PLAN.md`](../../../fiber/CUDA_PLAN.md). Архитектура
Ф0–Ф8 (SoA, primitives-без-указателей, contributions без атомарных
операций, count→scan→fill) изначально проектировалась с расчётом на второй
backend (ADR-008), поэтому перенос — не переписывание, а добавление
CUDA-реализации тех же формул рядом с CPU-путём.

Первый вопрос — каким компилятором собирать `.cu`-файлы в проекте, который
уже собирается Clang 22 (MSVC ABI, Ninja, C++20, `-Weverything`). Три
кандидата: Clang как компилятор CUDA (`enable_language(CUDA)` с
`CMAKE_CUDA_COMPILER=clang++`), `nvcc` с MSVC `cl` как хост-компилятором на
Windows, `nvcc` с Clang как хостом на Linux. Спайк вне репозитория
(2026-09-27) собрал пробную программу — ядро с device-лямбдой,
`cub::DeviceRadixSort::SortPairs`, `cub::DeviceScan::ExclusiveSum`,
`cub::DeviceReduce::Min`, `pow/exp/sqrt/log` на GPU против CPU, включение
`precision.hpp` и `interactions/solid_sph.hpp` целиком — в fp64 и fp32,
всеми тремя связками.

Оборудование машины разработки: NVIDIA GeForce RTX 4050 Laptop (Ada
Lovelace, `sm_89`, 6 ГБ, драйвер 617.14, режим WDDM), CUDA Toolkit 13.4
(рядом установлен 12.6), MSVC Build Tools VS 2026 (toolset 14.51, `cl`
19.51), CMake 4.3.2.

## Решение

**Компилятор `.cu` — `nvcc`; Clang как компилятор CUDA не используется
(решение владельца, Р-Г1).** Обе связки технически работают (таблица
ниже), но выбрана `nvcc` + Clang: `.cu`-файлы компилирует `nvcc` 13.4, весь
остальной проект (`.cpp`, тесты) — Clang 22, в одной сборке
(`CXX=clang++`, `CUDA=nvcc`). Следствия:

* На Windows у `nvcc` единственный поддерживаемый хост-компилятор для
  хост-части `.cu` — MSVC `cl` (`-ccbin clang` там не поддерживается);
  сборка идёт из окружения `vcvars64.bat` (скрипт
  [`tools/cuda_env.bat`](../../../tools/cuda_env.bat) и пресет `cuda`,
  [`CMakePresets.json`](../../../CMakePresets.json)). Объектники `cl` и
  Clang линкуются вместе без проблем — оба используют MSVC ABI и одну STL.
* На Linux хост-компилятором `nvcc` по умолчанию тоже становится Clang
  (`CMAKE_CUDA_HOST_COMPILER=clang++`); если версия `nvcc` на конкретной
  машине не поддерживает установленный Clang как хост — откат на GCC.
* Заголовки, попадающие в `.cu` (`precision.hpp`, будущий
  `include/physics/*`, представления состояния), обязаны собираться и
  Clang, и `cl`: без clang-специфичных расширений, конструкции под
  `#ifdef __clang__`.

### Проверка связок (спайк 2026-09-27, повторена при написании `test_cuda_toolchain`)

| Связка | Результат | Обязательные флаги/условия |
|---|---|---|
| Clang 22 как компилятор CUDA + SDK 13.4 | работает; математика GPU/CPU: 4e-16 (fp64), 2e-7 (fp32) | `-DNOMINMAX` (иначе макросы `min`/`max` из `windows.h` ломают CUB); окружение MSVC не нужно |
| `nvcc` 13.4 + MSVC 14.51 (выбрано) | работает, те же числа | `-Xcompiler=/Zc:preprocessor` (иначе CCCL останавливается с `#error`), `/Zc:__cplusplus`, `--extended-lambda`; из окружения `vcvars64.bat` |
| `nvcc` 12.6 + MSVC 14.51 | **не работает** | `nvcc` 12.6 отклоняет MSVC 14.51 как неподдерживаемый (`host_config.h`, C1189); с `-allow-unsupported-compiler` версийная проверка снимается, но `cudafe++` падает с `ACCESS_VIOLATION` (0xC0000005) на препроцессированном файле — не версийный гейт, а реальная несовместимость. Минимальная версия CUDA в проекте (12.4, Р-Г1а) остаётся требованием ко второй, Linux-машине с другим тулчейном, а не к этой |
| Clang 22 + SDK 12.6 | не работает | ошибка в libcu++ (`aligned_storage`) при компиляции хост-части |

Числа `test_cuda_toolchain`/`test_cuda_toolchain_fp32` при интеграции в
CMake (этап Г0, а не только внешний спайк): математика GPU/CPU —
**4.26e-16 отн.** (fp64, допуск 1e-15), **2.27e-07 отн.** (fp32, допуск
1e-6); `cub::DeviceRadixSort::SortPairs` (со стабильностью, n=65536),
`cub::DeviceScan::ExclusiveSum` (n=16384), `cub::DeviceReduce::Min`
(n=16384) — точное совпадение с CPU в обеих точностях; повторный запуск
сортировки в том же процессе даёт тот же результат (детерминизм GPU-пути,
Р-Г5).

**Архитектуры устройства — `SPH_CUDA_ARCHITECTURES`, по умолчанию
`89-real;75-virtual` (Р-Г1а).** Машинный код для карты разработки
(`sm_89`) плюс PTX архитектуры Turing (75) — младшей, которую поддерживает
CUDA 13: PTX X запускается на любой карте с compute capability ≥ X, поэтому
эта же сборка запустится и на второй машине владельца (Linux, NVIDIA,
модель неизвестна). Минимальная версия CUDA в проекте — **12.4**, а не 13:
CUDA 13 не поддерживает карты старше Turing, а на Linux-машине может
оказаться, например, Pascal или Volta. Linux-пресет `cuda-linux`
использует `SPH_CUDA_ARCHITECTURES=native` (CMake сам определяет карту
сборочной машины).

**Остальные решения плана (Р-Г2…Р-Г10)** приняты вместе с планом и здесь
только перечислены — подробности и обоснование в
[`fiber/CUDA_PLAN.md`](../../../fiber/CUDA_PLAN.md), §2:

* Р-Г2 — одна математика на оба backend: физика одного primitive/частицы —
  `SPH_HD`-функции (`__host__ __device__`), CPU-kernels переводятся на них
  до написания GPU-кода, под побитной регрессией (этап Г2).
* Р-Г3 — CPU-путь остаётся эталоном и основной сборкой: `SPH_ENABLE_CUDA`
  по умолчанию `OFF`, при `OFF` сборка и 62 CPU-теста не меняются.
* Р-Г4 — состояние живёт на GPU весь расчёт, копирование на хост — только
  для вывода и диагностики.
* Р-Г5 — reduce без атомарных операций, тот же канонический порядок
  `(interaction_type, primitive_id, slot)`, что и в SPH.md; GPU-путь
  детерминирован от запуска к запуску.
* Р-Г6 — точность через `real_t`/`SPH_PRECISION` (этап ГF); **рабочий
  режим GPU — fp32**, fp64 остаётся эталоном корректности (решение
  владельца, 2026-09-27: на BC6 с контактом типов 2–3 fp32 расходится с
  fp64 на пороговой физике, см. [FP32_REPORT.md](../../../fiber/FP32_REPORT.md)).
* Р-Г7 — порядок этапов: 3D Solid → нити → 2D осесимметричный Solid.
* Р-Г8 — ограничения CPU-версии (радиус поиска нитей по `l0`, `sle=2` не
  поддерживается, нет трения и самоконтакта нити) переносятся как есть.
* Р-Г9 — индексы на GPU — `uint32_t` (на CPU остаётся `size_t`).
* Р-Г10 — процесс: код → независимое ревью → сборка → полный прогон
  тестов, как в Ф0–Ф8.

## Каркас сборки (этап Г0)

* [`cmake/SphCuda.cmake`](../../../cmake/SphCuda.cmake): опция
  `SPH_ENABLE_CUDA` (по умолчанию `OFF`, при `OFF` файл сразу
  `return()`), `find_package(CUDAToolkit 12.4 REQUIRED)`,
  `enable_language(CUDA)` с проверками (CUDA < 12.4, `CMAKE_CUDA_COMPILER`
  указывает на Clang, на Windows `cl` не найден в PATH — все три дают
  понятный `FATAL_ERROR` со ссылкой на `tools/cuda_env.bat`), кэш-переменная
  `SPH_CUDA_ARCHITECTURES`, функция `sph_cuda_configure_target()` с общими
  флагами (`--extended-lambda`, `-DNOMINMAX`, на Windows —
  `-Xcompiler=/Zc:preprocessor,/Zc:__cplusplus,/utf-8,/EHsc,/W4`) и
  заготовка `SPH_CUDA_STRICT_FP` (`-fmad=false`, план §5, используется с
  Г4).
* [`tools/cuda_env.bat`](../../../tools/cuda_env.bat): находит
  `vcvars64.bat` (через `vswhere`, с фиксированным путём BuildTools как
  запасным вариантом), включает окружение MSVC, выполняет переданную
  команду и пробрасывает её код возврата.
* [`CMakePresets.json`](../../../CMakePresets.json): пресет `cuda`
  (наследует `clang-windows`, добавляет `SPH_ENABLE_CUDA=ON`,
  `CMAKE_CUDA_COMPILER=nvcc`, каталог `build-cuda`) и `cuda-linux`
  (Ninja, `nvcc`, хост Clang, `SPH_CUDA_ARCHITECTURES=native`) —
  непроверен, второй машины нет.
* [`tests/cuda/test_cuda_toolchain.cu`](../../../tests/cuda/test_cuda_toolchain.cu):
  первый CUDA-тест, цели `test_cuda_toolchain` (метка `cuda`) и
  `test_cuda_toolchain_fp32` (метки `cuda;fp32`), `SKIP_RETURN_CODE 77`
  при отсутствии GPU.
* [`bench/cuda_microbench.cu`](../../../bench/cuda_microbench.cu): не тест
  CTest, отдельная цель.

### Числа микробенчмарка (RTX 4050 Laptop, sm_89, драйвер 617.14)

| Метрика | Значение |
|---|---|
| Пропускная способность памяти (copy, чтение+запись) | 169.8 ГБ/с |
| Пиковая арифметика FP32 (цепочки FMA) | 9356 GFLOP/s |
| Пиковая арифметика FP64 (цепочки FMA) | 170.6 GFLOP/s |
| Отношение FP32/FP64 | ≈ 54.8× |

Согласуется с ожиданием плана (§1.1): FP64 на этой карте на два порядка
медленнее FP32, а не «даром» — рабочий режим GPU оправданно выбран fp32
(Р-Г6), fp64 остаётся эталоном корректности, а не целью по скорости.
Пропускная способность (169.8 ГБ/с) — около 88 % паспортных ~192 ГБ/с,
разумно для однопроходного `cudaMemcpy`-подобного копирования без
дополнительной настройки под конкретную карту.

## Проверка критериев готовности Г0

```text
cmake --build build                              # SPH_ENABLE_CUDA=OFF (по умолчанию)
ctest --test-dir build -E full -j 6               # 62/62 CPU-тестов, без изменений

tools\cuda_env.bat cmake --preset cuda
tools\cuda_env.bat cmake --build build-cuda -j 8
tools\cuda_env.bat ctest --test-dir build-cuda -L cuda --output-on-failure
                                                   # test_cuda_toolchain, test_cuda_toolchain_fp32 — зелёные
ctest --test-dir build-cuda -E full -LE cuda -j 6 # те же 62 CPU-теста, в build-cuda
```

Оба прогона дали 62/62 и 2/2 соответственно (числа сверки — выше).

## Последствия

**Положительные:**

* Оба технически рабочих варианта компилятора зафиксированы числами, а не
  предположением — решение владельца (`nvcc`) опирается на измеренное
  поведение, а не только на предпочтение;
* `SPH_ENABLE_CUDA=OFF` по умолчанию гарантирует, что CPU-сборка и её 62
  теста не видят этого этапа вообще — риск регресса CPU-пути от
  интеграции CUDA в CMake минимален;
* `test_cuda_toolchain` — недорогая (< 1 c на обе точности) защита от
  скрытой поломки тулчейна при будущих обновлениях CUDA/драйвера/MSVC/
  Clang, задолго до того, как будет что-то реальное в `include/physics/*`
  (Г2+);
* Понятные `FATAL_ERROR` в `cmake/SphCuda.cmake` (версия CUDA, Clang как
  компилятор CUDA, отсутствие `cl` в PATH) экономят время на диагностику
  типичных ошибок конфигурации.

**Отрицательные:**

* Связка `nvcc` + `cl` привязывает сборку CUDA-целей на Windows к
  окружению MSVC (`tools/cuda_env.bat`) — обычный `cmake --build` без него
  не соберёт `build-cuda`; альтернатива (Clang как компилятор CUDA) не
  имела бы этого ограничения, но не выбрана владельцем;
* CUDA 12.6, установленная рядом на машине разработки, не годится для
  проверки с текущим MSVC 14.51 (`cudafe++` падает) — сборка на CUDA 12.x
  остаётся проверенной только по документации NVIDIA и ожидается рабочей
  на менее новом MSVC (например, toolset VS 2022), но не подтверждена
  измерением на этой машине;
* Пресет `cuda-linux` и вторая машина (Р-Г1а) не проверены — архитектура и
  версия CUDA там неизвестны до первого доступа к машине.

## Связанные файлы

- [`fiber/CUDA_PLAN.md`](../../../fiber/CUDA_PLAN.md): план переноса на
  GPU, §1.1 (тулчейн), §2 (решения Р-Г1…Р-Г10), §4 (этап Г0).
- [`cmake/SphCuda.cmake`](../../../cmake/SphCuda.cmake): опция
  `SPH_ENABLE_CUDA`, подключение `nvcc`, архитектуры, флаги.
- [`tools/cuda_env.bat`](../../../tools/cuda_env.bat): окружение MSVC для
  сборки CUDA-целей на Windows.
- [`CMakePresets.json`](../../../CMakePresets.json): пресеты `cuda`,
  `cuda-linux`.
- [`tests/cuda/test_cuda_toolchain.cu`](../../../tests/cuda/test_cuda_toolchain.cu),
  [`bench/cuda_microbench.cu`](../../../bench/cuda_microbench.cu).
- [`docs/02-how-to/02-build-cuda.md`](../../02-how-to/02-build-cuda.md):
  инструкция по сборке и запуску.
- [ADR-008](adr-008-rhs-architecture-fiber.md): архитектура RHS, на
  которой строится GPU-перенос (SoA, primitives, contributions без
  атомарных операций).
- [FP32_REPORT.md](../../../fiber/FP32_REPORT.md): решение о fp32 как
  рабочем режиме (Р-Г6) и его физические ограничения.
