# Опции CMake

Опции сборки определены в [`CMakeLists.txt`](../../../CMakeLists.txt:36) и
управляют компиляцией проекта SPH.

| Опция | По умолчанию | Описание |
|---|---|---|
| `SPH_USE_OPENMP` | `ON` | Включить OpenMP (`-fopenmp`) |
| `SPH_NATIVE_ARCH` | `OFF` | Использовать `-march=native` (выкл. для воспроизводимого baseline) |
| `SPH_OLEVEL` | `2` | Уровень оптимизации Clang/GCC (2 или 3) |
| `SPH_PROFILE` | `OFF` | Определяет макрос `SPH_PROFILE` для цели `SPH`. Блоков `#ifdef SPH_PROFILE` в текущем коде нет, `force_timing.csv` не создаётся |
| `SPH_DIM` | `1` | Размерность исполняемого `SPH` (1, 2 или 3). Тестовые цели задают свою размерность сами. Нити (`fiber_file`) и `k_mode=solid_fiber` требуют `3` |
| `SPH_CONTACT_STATS` | `OFF` | Диагностика контактов нитей (`#ifdef SPH_CONTACT_STATS`, [`fiber_contact_stats.hpp`](../../../include/interactions/fiber_contact_stats.hpp)); включает ключи `guo2019_shot --contact-stats` и `--hold-n` (определяется только для `guo2019_shot`, не для `SPH`) |
| `SPH_ENERGY_LOG` | `OFF` | Лог механической энергии (`#ifdef SPH_ENERGY_LOG`, [`energy_log.hpp`](../../../include/energy_log.hpp)); включает ключи `energy_log_every`/`energy_log_file` в `Initial_Conditions.txt` и `--energy-log` у `guo2019_shot`/`fabric_checks`. Определяется для `SPH`, `guo2019_shot`, `fabric_checks` |
| `SPH_E2E_SLOW_TESTS` | `OFF` | Регистрировать в ctest медленные полные прогоны BC2 и BC6 (метка `slow`) |
| `SPH_BC2_RUN_TIME_LIMIT` | `900` | Бюджет сторожа на один полный прогон BC2, с (`0` — без ограничения); общий для тестов и цели `bc2_plot` |
| `SPH_BC6_RUN_TIME_LIMIT` | `1800` | То же для BC6 |
| `SPH_BC2_PLOT`, `SPH_BC6_PLOT` | `contact/e2e/bc{2,6}_monitor_velocity.svg` | Куда цели `bc2_plot` / `bc6_plot` пишут график (путь к файлу) |

## Примеры

```bash
# Базовая Release-сборка
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release

# Нативная архитектура + -O3
cmake -S . -B build -G Ninja -DSPH_NATIVE_ARCH=ON -DSPH_OLEVEL=3

# Сборка с профилированием
cmake -S . -B build -G Ninja -DSPH_PROFILE=ON

# Трёхмерный исполняемый SPH (например, для импорта .k с нитями)
cmake -S . -B build -G Ninja -DSPH_DIM=3
```

## Целевые стандарты

- Язык: **C++20** (`CMAKE_CXX_STANDARD=20`, без расширений);
- Тип сборки по умолчанию — `Release` (если не задан `CMAKE_BUILD_TYPE`); флаги Release переопределены на `-O${SPH_OLEVEL}` без `NDEBUG`, так что `assert` остаются активными;
- Предупреждения: Clang — `-Weverything` (с подавлениями
  `-Wno-missing-prototypes`, `-Wno-c++98-compat*`), GCC — `-Wall -Wextra -Wpedantic`, MSVC — `/W4`; `-fexceptions` для Clang/GCC;
- Windows: определяется макрос `_WIN32`, UTF-8 консоль включается в коде.

## Пресеты

[`CMakePresets.json`](../../../CMakePresets.json) (генератор Ninja, каталог `build/`, компилятор LLVM в `C:/Program Files/LLVM/bin`) содержит три конфигурационных пресета; `SPH_DIM` и диагностические опции в них не заданы (по умолчанию `SPH_DIM=1`, `OFF`).

| Пресет | Что задаёт |
|---|---|
| `clang-windows` | Release, `SPH_USE_OPENMP=ON`, `SPH_NATIVE_ARCH=OFF`, `SPH_OLEVEL=2` — baseline |
| `clang-windows-debug` | наследует `clang-windows`, `CMAKE_BUILD_TYPE=Debug` |
| `clang-windows-native` | наследует `clang-windows`, `SPH_NATIVE_ARCH=ON`, `SPH_OLEVEL=3` |

Для каждого есть одноимённый build-пресет: `cmake --preset clang-windows`, затем `cmake --build --preset clang-windows`.

## Дополнительные цели

Помимо `SPH` в `CMakeLists.txt` собираются утилиты `guo2019_shot`, `fabric_checks`, `plain_weave_gen` (см. [параметры запуска](../cli.md)), демонстрации `demo_fiber_net`, `demo_fiber_layers` и `fiber_t0_dump` (только последняя с `EXCLUDE_FROM_ALL`), тесты (`enable_testing()`) и служебные цели `bc2_plot`, `bc6_plot`, `regression_golden` (скрипты в `cmake/`).
