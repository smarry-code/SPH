# ============================================================================
# CMake-сборка проекта SPH
# Этап 3.0 рефакторинга: переход с прямой компиляции clang++ (tasks.json)
# на систему сборки CMake. Точно воспроизводит baseline-конфигурацию:
#   Clang / C++20 / -O2 / -fopenmp / -fexceptions / -Weverything (+ подавления)
# Без -march=native по умолчанию (чтобы baseline совпадал с прежней сборкой).
# ============================================================================

cmake_minimum_required(VERSION 3.16)

project(SPH LANGUAGES CXX)

# --- AppleClang + OpenMP: системный clang не понимает голый -fopenmp,
#     нужен -Xpreprocessor -fopenmp и явные include/lib из Homebrew libomp
#     (fix(build), см. contact/01-plan.md, шаг 0, решение D5).
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE SPH_LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT SPH_LIBOMP_PREFIX)
        set(SPH_LIBOMP_PREFIX "/opt/homebrew/opt/libomp")
    endif()
    set(SPH_APPLECLANG_OMP_COMPILE_OPTIONS -Xpreprocessor -fopenmp -I${SPH_LIBOMP_PREFIX}/include)
    set(SPH_APPLECLANG_OMP_LINK_OPTIONS -L${SPH_LIBOMP_PREFIX}/lib -lomp)
endif()

# --- Стандарт: C++20 (как в прежней сборке -std=c++20) ----------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Опции ------------------------------------------------------------------
option(SPH_USE_OPENMP  "Включить OpenMP"            ON)
option(SPH_NATIVE_ARCH "Использовать -march=native"  OFF)   # OFF: baseline без -march
set(SPH_OLEVEL "2" CACHE STRING "Уровень оптимизации (2|3)") # baseline = -O2
set_property(CACHE SPH_OLEVEL PROPERTY STRINGS 2 3)
option(SPH_PROFILE     "Включить профилировочные замеры (#ifdef SPH_PROFILE)" OFF)  # OFF: вне горячего пути

# --- Конфигурация по умолчанию ----------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Исходники модулей расчётного конвейера (подэтап 5.2-В):
#     объявления — в include/*.hpp, определения — в src/*.cpp (main.cpp в src/).
#     Автосборка всех src/*.cpp (точки входа main.cpp и модулей).
file(GLOB SPH_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
# --- Исполняемый модуль (имя SPH.exe на Windows — как прежде) ----------------
add_executable(SPH ${SPH_SOURCES})
set_target_properties(SPH PROPERTIES OUTPUT_NAME "SPH")

# --- Тесты (лежат в корне проекта; подключают заголовки из include/) --------
enable_testing()
add_executable(test_flat_grid test_flat_grid.cpp)
target_include_directories(test_flat_grid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME flat_grid_test COMMAND test_flat_grid)

add_executable(test_pair_graph test_pair_graph.cpp)
target_include_directories(test_pair_graph PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME pair_graph_test COMMAND test_pair_graph)
# pair_graph.hpp использует OpenMP (спекулятивная блочно-параллельная окраска, фаза B1).
# Без -fopenmp omp_get_max_threads()==1 и параллельный путь не тестируется.
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_pair_graph PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_pair_graph PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_pair_graph PRIVATE -fopenmp)
    target_link_options(test_pair_graph PRIVATE -fopenmp)
endif()

# --- Тест моста AoS <-> SoA (Этап 4.3.1): тождество gather/scatter. ---
add_executable(test_particle_system test_particle_system.cpp)
target_include_directories(test_particle_system PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME particle_system_test COMMAND test_particle_system)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(test_particle_system PRIVATE -Wno-unused-parameter)
endif()

# --- Тест контактных сил (part_0.md, ADR-007): аналитика типов 1-3,
#     баланс импульса, ограничитель шага. Линкуется с src/ext_force.cpp,
#     чтобы покрыть и оркестратор.
add_executable(test_contact test_contact.cpp src/ext_force.cpp)
target_include_directories(test_contact PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME contact_test COMMAND test_contact)

# --- Бенчмарк фазы B (Этап 1): сравнение серийного (B2) и параллельного (B1) путей.
add_executable(bench_color bench_color.cpp)
target_include_directories(bench_color PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(bench_color PRIVATE -O2 -Wno-unused-parameter)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(bench_color PRIVATE -O2)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(bench_color PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(bench_color PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(bench_color PRIVATE -fopenmp)
    target_link_options(bench_color PRIVATE -fopenmp)
endif()

# --- Предупреждения: воспроизводим -Weverything (Clang) с прежними подавлениями
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(SPH PRIVATE
        -Weverything
        -Wno-missing-prototypes
        -Wno-c++98-compat
        -Wno-c++98-compat-pedantic
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(SPH PRIVATE -Wall -Wextra -Wpedantic)
elseif(MSVC)
    target_compile_options(SPH PRIVATE /W4)
endif()

# --- Оптимизация: переопределяем дефолтные Release-флаги CMake (-O3 -DNDEBUG),
#     чтобы ТОЧНО воспроизвести прежнюю baseline-сборку tasks.json
#     (только -O2, без NDEBUG — assert остаются активными, как раньше).
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "-O${SPH_OLEVEL}" CACHE STRING "Release C++ flags" FORCE)
endif()

# --- Исключения (-fexceptions) — как в прежней сборке -----------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(SPH PRIVATE -fexceptions)
endif()

# --- OpenMP -----------------------------------------------------------------
# Прямой -fopenmp для Clang на Windows (надёжнее find_package(OpenMP),
# который не всегда находит libomp в дистрибутиве LLVM/clang).
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(SPH PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(SPH PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(SPH PRIVATE -fopenmp)
    target_link_options(SPH PRIVATE -fopenmp)
endif()

# --- -march=native — ВЫКЛЮЧЕН по умолчанию (baseline); включается опцией
#     как предпосылка SIMD-активации (Этап 4.4) -------------------------------
if(SPH_NATIVE_ARCH AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(SPH PRIVATE -march=native)
endif()

# --- Windows: UTF-8 консоль (см. main.cpp SetConsoleOutputCP) ----------------
if(WIN32)
    target_compile_definitions(SPH PRIVATE _WIN32)
endif()

# --- Профилирование: включает блоки #ifdef SPH_PROFILE в pair_graph/int_force
#     (тайминг фаз build и параллельных циклов, запись в force_timing.csv) ----
if(SPH_PROFILE)
    target_compile_definitions(SPH PRIVATE SPH_PROFILE)
endif()

# --- Каталоги include: заголовки — в include/, исходники — в src/ ------------
target_include_directories(SPH PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
