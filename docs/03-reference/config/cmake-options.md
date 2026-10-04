# Опции CMake

Опции сборки определены в [`CMakeLists.txt`](../../../CMakeLists.txt:19) и
управляют компиляцией проекта SPH.

| Опция | По умолчанию | Описание |
|---|---|---|
| `SPH_USE_OPENMP` | `ON` | Включить OpenMP (`-fopenmp`) |
| `SPH_NATIVE_ARCH` | `OFF` | Использовать `-march=native` (выкл. для воспроизводимого baseline) |
| `SPH_OLEVEL` | `2` | Уровень оптимизации Clang/GCC (2 или 3) |
| `SPH_PROFILE` | `OFF` | Включить замеры фаз (`#ifdef SPH_PROFILE`, запись в `force_timing.csv`) |
| `SPH_DIM` | `1` | Размерность исполняемого `SPH` (1, 2 или 3). Тестовые цели задают свою размерность сами. Нити (`fiber_file`) и `k_mode=solid_fiber` требуют `3` |

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
- Базовые флаги Clang: `-O2 -fexceptions -Weverything` (с подавлениями
  `-Wno-missing-prototypes`, `-Wno-c++98-compat*`);
- Windows: определяется макрос `_WIN32`, UTF-8 консоль включается в коде.
