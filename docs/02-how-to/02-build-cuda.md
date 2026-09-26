# Сборка и запуск CUDA-бэкенда

CUDA-бэкенд (см. [`fiber/CUDA_PLAN.md`](../../../fiber/CUDA_PLAN.md)) — опция
сборки, выключенная по умолчанию: обычная CPU-сборка (пресет
`clang-windows`, каталог `build`) не видит ни языка CUDA, ни новых целей и
не меняется этим этапом. Детали устройства — [ADR-009](../04-explanation/decisions/adr-009-cuda-backend.md).

## Требования

* CUDA Toolkit **12.4 или новее** (проверено на 13.4; см. ADR-009 про 12.6).
* На Windows — компилятор `.cu` это **`nvcc`**, а его хост-компилятор — **MSVC
  `cl`** (решение Р-Г1: связка `nvcc` + Clang, `-ccbin clang` на Windows не
  поддерживается). `cl` не попадает в PATH сам по себе — нужно окружение
  `vcvars64.bat`.
* Остальной проект (`.cpp`, тесты) по-прежнему собирает Clang — как в
  CPU-сборке.

## Windows: пресет `cuda`

`cl.exe` есть только внутри окружения MSVC, поэтому конфигурация и сборка
идут через обёртку [`tools/cuda_env.bat`](../../../tools/cuda_env.bat): она
находит `vcvars64.bat` (через `vswhere`, с фиксированным путём BuildTools
как запасным вариантом), включает окружение MSVC и запускает переданную
команду в нём.

```bash
tools\cuda_env.bat cmake --preset cuda
tools\cuda_env.bat cmake --build build-cuda -j 8
tools\cuda_env.bat ctest --test-dir build-cuda -L cuda --output-on-failure
```

Пресет `cuda` (`CMakePresets.json`) наследует `clang-windows` (Clang, C++20,
OpenMP, `-O2`) и добавляет `SPH_ENABLE_CUDA=ON`, `CMAKE_CUDA_COMPILER=nvcc`,
каталог сборки `build-cuda` — отдельный от CPU-сборки `build`.

!!! note "Почему CPU-тесты в build-cuda тоже запускаются через cuda_env.bat"
    Сама сборка и конфигурация требуют окружения MSVC (иначе `nvcc` не
    находит `cl`), но уже собранные `.exe` от окружения не зависят — их
    можно запускать и без `tools\cuda_env.bat`, например:
    `ctest --test-dir build-cuda -E full -LE cuda -j 6`.

### Что проверяют CUDA-тесты этапа Г0

`test_cuda_toolchain` (метка `cuda`, `SPH_PRECISION=64`) и
`test_cuda_toolchain_fp32` (метки `cuda;fp32`, `SPH_PRECISION=32`) — первый
CUDA-тест проекта ([`tests/cuda/test_cuda_toolchain.cu`](../../../tests/cuda/test_cuda_toolchain.cu)).
Он ловит поломку тулчейна (обновление CUDA/драйвера/MSVC/Clang), а не физику:
device-лямбда, `cub::DeviceRadixSort::SortPairs` (стабильность),
`cub::DeviceScan::ExclusiveSum`, `cub::DeviceReduce::Min`, сверка
`pow/exp/sqrt/log` в `real_t` между GPU и CPU, включение заголовков проекта
(`precision.hpp`, `param.hpp`) через `cl`. Печатает имя устройства, `sm_XY`
и объём памяти. Если `cudaGetDeviceCount` не находит рабочий GPU, тест
завершается кодом 77, и CTest показывает его как пропущенный (`Not Run`), а
не проваленный (`SKIP_RETURN_CODE 77` в `CMakeLists.txt`).

### Полный набор в сборке с CUDA

Критерий готовности Г0 — оба режима зелёные в одной и той же сборке
`build-cuda`:

```bash
tools\cuda_env.bat ctest --test-dir build-cuda -L cuda --output-on-failure
ctest --test-dir build-cuda -E full -LE cuda -j 6 --output-on-failure
```

Второй набор — те же 62 CPU-теста, что и в обычной сборке `build`; CUDA не
меняет их поведение (Р-Г3: CPU-путь остаётся эталоном и основной сборкой).

## Другая видеокарта

Архитектуры задаются кэш-переменной `SPH_CUDA_ARCHITECTURES`
(`CMAKE_CUDA_ARCHITECTURES`, `cmake/SphCuda.cmake`). По умолчанию —
`89-real;75-virtual`: машинный код под RTX 4050 Laptop (карта разработки,
`sm_89`) и PTX архитектуры Turing (75) — младшей, которую поддерживает
CUDA 13; PTX X запускается на любой карте с compute capability ≥ X, поэтому
одна и та же сборка `build-cuda` запускается и на другой карте (Р-Г1а).

Чтобы собрать под конкретную карту явно (например, если целевая карта старше
Turing и нужна сборка только под неё):

```bash
tools\cuda_env.bat cmake --preset cuda -DSPH_CUDA_ARCHITECTURES=75-real
```

Значение `native` определяет архитектуру карты, установленной на самой
сборочной машине (CMake ≥ 3.24) — так настроен пресет `cuda-linux` (ниже).

## Linux: пресет `cuda-linux`

```bash
cmake --preset cuda-linux
cmake --build build-cuda -j 8
ctest --test-dir build-cuda -L cuda --output-on-failure
```

На Linux хост-компилятором `nvcc` по умолчанию тоже назначается Clang
(`CMAKE_CUDA_HOST_COMPILER=clang++`, решение Р-Г1) — окружения вроде
`vcvars64.bat` там не нужно, `cmake/SphCuda.cmake` не требует отдельного
скрипта-обёртки для Linux. `SPH_CUDA_ARCHITECTURES=native` подбирает
архитектуру под карту сборочной машины. Если конкретная связка `nvcc` +
Clang на машине не собирается, откат на GCC — передать явно:

```bash
cmake --preset cuda-linux -DCMAKE_CUDA_HOST_COMPILER=g++
```

!!! warning "Пресет не проверен"
    Вторая машина владельца (Linux, карта NVIDIA неизвестной модели,
    решение Р-Г1а) на момент этапа Г0 недоступна. Пресет `cuda-linux`
    создан по тем же принципам, что и `cuda`, но живой прогон
    `test_cuda_toolchain` там не выполнялся — см. ADR-009.

## Микробенчмарк тулчейна

`cuda_microbench` — не тест CTest, а отдельная цель: пропускная способность
памяти и пиковая арифметика FP32/FP64 (цепочки FMA) на установленной карте.
Числа нужны заранее, до переноса физики (Г2+), чтобы знать потолок
ускорения каждого режима точности (Р-Г6: рабочий режим GPU — fp32).

```bash
tools\cuda_env.bat cmake --build build-cuda --target cuda_microbench
build-cuda\cuda_microbench.exe
```

Числа на карте разработки (RTX 4050 Laptop) — в
[ADR-009](../04-explanation/decisions/adr-009-cuda-backend.md).

## Опции сборки

| Переменная | По умолчанию | Назначение |
|---|---|---|
| `SPH_ENABLE_CUDA` | `OFF` | Включает язык CUDA, цели `test_cuda_toolchain*` и `cuda_microbench` |
| `SPH_CUDA_ARCHITECTURES` | `89-real;75-virtual` | `CMAKE_CUDA_ARCHITECTURES`; `native` — определить по карте сборочной машины |
| `SPH_CUDA_STRICT_FP` | `OFF` | Заготовка строгого режима FP (`-fmad=false`) для разбора расхождений CPU/GPU (план §5); используется начиная с Г4 |

## Понятные ошибки конфигурации

`cmake/SphCuda.cmake` останавливает конфигурацию с объяснением, если:

* CUDA Toolkit найден версией меньше 12.4 (Р-Г1а: CUDA 13 не поддерживает
  карты старше Turing);
* `CMAKE_CUDA_COMPILER` указывает на Clang (решение Р-Г1: компилятор `.cu` в
  проекте — только `nvcc`; Clang как компилятор CUDA технически работает,
  но не используется);
* сборка на Windows запущена не из окружения MSVC (`cl.exe` не найден в
  PATH) — вместо невнятной ошибки CMake сразу называется причина и
  подсказывается `tools\cuda_env.bat`.
