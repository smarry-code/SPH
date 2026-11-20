# Быстрый старт

Это руководство проведёт вас от пустого репозитория до запущенного расчёта.
Следуя ему шаг за шагом, вы соберёте проект и получите первый набор выходных
данных.

## Предварительные требования

| Инструмент | Минимальная версия | Назначение |
|---|---|---|
| **Clang** (LLVM) | 14+ | Компилятор C++20 |
| **CMake** | 3.21+ | Система сборки |
| **Ninja** (рекомендуется) | — | Генератор сборки |
| **OpenMP** (libomp) | — | Параллелизм расчёта |
| **Python 3** | 3.8+ | Сборка сайта документации (необязательно) |

!!! note "Windows"
    На Windows используется LLVM/Clang с `libomp`. Прямой флаг `-fopenmp`
    применяется в [`CMakeLists.txt`](../../CMakeLists.txt:1616) вместо
    `find_package(OpenMP)` — это надёжнее для дистрибутива LLVM под Windows.
    UTF-8 консоль включается вызовом `SetConsoleOutputCP(65001)` в коде.

## Шаг 1. Клонирование

```bash
git clone <url-репозитория> SPH
cd SPH
```

## Шаг 2. Сборка

Базовая конфигурация (`Release`, `-O2`, OpenMP включён) воспроизводит baseline:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

После сборки исполняемый файл `SPH` (на Windows — `SPH.exe`) окажется в
`build/`.

Те же конфигурации заданы пресетами [`CMakePresets.json`](../../CMakePresets.json)
(Ninja, каталог `build/`): `clang-windows` (Release, `-O2`), `clang-windows-debug`
(Debug) и `clang-windows-native` (`-march=native`, `-O3`). В пресетах прописаны пути
`C:/Program Files/LLVM/bin/clang.exe` и `clang++.exe`; на другой машине используйте
командные строки выше.

```bash
cmake --preset clang-windows
cmake --build --preset clang-windows
```

### Опции сборки

Опции CMake (определены в [`CMakeLists.txt`](../../CMakeLists.txt:36)):

| Опция | По умолчанию | Описание |
|---|---|---|
| `SPH_USE_OPENMP` | `ON` | Включить OpenMP |
| `SPH_NATIVE_ARCH` | `OFF` | Использовать `-march=native` (выкл. для воспроизводимого baseline) |
| `SPH_OLEVEL` | `2` | Уровень оптимизации (2 или 3) |
| `SPH_PROFILE` | `OFF` | Определяет макрос `SPH_PROFILE` для цели `SPH`; в текущих `src/` и `include/` блоков `#ifdef SPH_PROFILE` нет, так что опция ничего не меняет |
| `SPH_DIM` | `1` | Размерность исполняемого `SPH` (1, 2 или 3); нити и импорт `.k` в режиме `solid_fiber` требуют 3 |
| `SPH_ENERGY_LOG` | `OFF` | Лог механической энергии ([`include/energy_log.hpp`](../../include/energy_log.hpp)); применяется к `SPH`, `guo2019_shot`, `fabric_checks` |
| `SPH_CONTACT_STATS` | `OFF` | Диагностика контактов нитей (только `guo2019_shot`) |
| `SPH_E2E_SLOW_TESTS` | `OFF` | Регистрировать медленные тесты (метка `slow`), см. [запуск тестов](../02-how-to/01-run-tests.md) |

Пример включения нативной архитектуры и `-O3`:

```bash
cmake -S . -B build -G Ninja -DSPH_NATIVE_ARCH=ON -DSPH_OLEVEL=3
```

## Шаг 3. Запуск первого расчёта

По умолчанию сборка настроена на тест **«ударная труба»** (флаг
[`shocktube = true`](../../include/param.hpp:283)) и одномерную размерность
([`dim = SPH_DIM`](../../include/param.hpp:16), по умолчанию `SPH_DIM=1`).

```bash
cd build
./SPH
```

Программа:

1. Читает параметры из [`Initial_Conditions.txt`](../03-reference/config/initial-conditions.md)
   (функция [`parse_sph_config()`](../../src/input.cpp:110));
2. Генерирует начальную конфигурацию 5000 частиц (ударная труба);
3. Выполняет временной цикл интегрирования с выводом прогресса в консоль;
4. Сохраняет кадры состояния в файлы `sph_output_*.txt`.

!!! tip "Каталог запуска"
    Входной файл `Initial_Conditions.txt` и выходные файлы создаются в
    **текущем рабочем каталоге**. Запускайте `SPH` из каталога, где разрешена
    запись (например, `build/` или отдельный рабочий каталог сценария).
    Исполняемый файл `SPH` не принимает аргументов командной строки. Если
    `Initial_Conditions.txt` в каталоге запуска нет, ошибки не будет: параметры
    останутся значениями по умолчанию. Скопируйте файл из корня репозитория в
    рабочий каталог.

## Шаг 4. Где искать результаты

После завершения расчёта в рабочем каталоге появятся:

| Файл | Описание |
|---|---|
| `sph_output_NNNNN.txt` | Кадры состояния частиц (см. [формат вывода](../03-reference/file-formats/dat-output.md)) |
| `sph_gp.csv` | Накопительный дамп координат/скоростей всех частиц по времени |
| `log.txt` | Журнал расчёта (включая мониторируемую частицу) |

## Шаг 5. Проверка работоспособности (тесты)

В проекте настроены CTest-тесты (регистрируются после [`enable_testing()`](../../CMakeLists.txt:72)):

```bash
ctest --test-dir build --output-on-failure
```

Набор большой, от юнит-тестов до сквозных прогонов. Для примера:

- `spatial_index_test` — сеточный поиск соседей (`SpatialIndex`);
- `executor_test`, `contributions_test` — исполнитель `CpuExecutor` и редукция вкладов;
- `particle_system_test` — мост AoS↔SoA;
- `contact_test` — контактные силы «частица-частица».

Полный список, метки (`3d`, `e2e`, `slow`, `regression`, `jc`) и запуск по меткам —
в [руководстве по тестам](../02-how-to/01-run-tests.md). Сквозные тесты и регрессия
считают минуты и используют несколько потоков OpenMP.

## Шаг 6. Просмотр документации (необязательно)

Локальный сайт документации:

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

Откройте <http://127.0.0.1:8000>.

---

## Что дальше?

- [Спецификации входных/выходных файлов](../03-reference/file-formats/dat-output.md)
- [Обзор архитектуры расчётного конвейера](../04-explanation/architecture/pipeline-overview.md)
- [Краткое руководство по форматам](../03-reference/config/initial-conditions.md)
