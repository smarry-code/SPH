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
option(SPH_CONTACT_STATS "Диагностика контактов нитей (#ifdef SPH_CONTACT_STATS, include/interactions/fiber_contact_stats.hpp)" OFF)
option(SPH_ENERGY_LOG  "Лог механической энергии (#ifdef SPH_ENERGY_LOG, include/energy_log.hpp)" OFF)  # OFF: кода лога в сборке нет
# Размерность основного исполняемого SPH (1|2|3). Тестовые цели задают свою
# SPH_DIM сами и от этой опции не зависят. Нити и импорт .k в режиме
# solid_fiber требуют 3.
set(SPH_DIM "1" CACHE STRING "Размерность исполняемого SPH (1|2|3)")
set_property(CACHE SPH_DIM PROPERTY STRINGS 1 2 3)
if(NOT SPH_DIM MATCHES "^[123]$")
    message(FATAL_ERROR "SPH_DIM должен быть 1, 2 или 3, получено: ${SPH_DIM}")
endif()

# --- Конфигурация по умолчанию ----------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Исходники модулей расчётного конвейера (подэтап 5.2-В):
#     объявления — в include/*.hpp, определения — в src/*.cpp (main.cpp в src/).
#     Автосборка всех src/*.cpp (точки входа main.cpp и модулей).
file(GLOB SPH_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
# Модули конвейера БЕЗ точки входа: переиспользуются целями тестов, которые
# приносят собственный main() (ТЗ contact/04-e2e-test-brief.md §7 — общий
# список исходников в переменной, без второго file(GLOB)).
set(SPH_PIPELINE_SOURCES ${SPH_SOURCES})
list(REMOVE_ITEM SPH_PIPELINE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
# --- Исполняемый модуль (имя SPH.exe на Windows — как прежде) ----------------
add_executable(SPH ${SPH_SOURCES})
set_target_properties(SPH PROPERTIES OUTPUT_NAME "SPH")
target_compile_definitions(SPH PRIVATE SPH_DIM=${SPH_DIM})

# --- Тесты (исходники в tests/, данные в tests/data/; заголовки из include/) ---
enable_testing()

# --- Инфраструктура нитей, этап Ф1 (fiber/PLAN.md): CpuExecutor,
#     ContributionWriter/reduce, SpatialIndex, RhsWorkspace (include/rhs/).
#     Чистая инфраструктура без физики — не зависит от SPH_DIM/param.hpp,
#     размерность задаётся шаблонным параметром SpatialIndex<Dim>, поэтому
#     цели собираются один раз с параметрами по умолчанию.
add_executable(test_executor tests/test_executor.cpp)
target_include_directories(test_executor PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME executor_test COMMAND test_executor)
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_executor PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_executor PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_executor PRIVATE -fopenmp)
    target_link_options(test_executor PRIVATE -fopenmp)
endif()

add_executable(test_spatial_index tests/test_spatial_index.cpp)
target_include_directories(test_spatial_index PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME spatial_index_test COMMAND test_spatial_index)
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_spatial_index PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_spatial_index PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_spatial_index PRIVATE -fopenmp)
    target_link_options(test_spatial_index PRIVATE -fopenmp)
endif()

add_executable(test_contributions tests/test_contributions.cpp)
target_include_directories(test_contributions PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME contributions_test COMMAND test_contributions)
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_contributions PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_contributions PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_contributions PRIVATE -fopenmp)
    target_link_options(test_contributions PRIVATE -fopenmp)
endif()

# --- Тест моста AoS <-> SoA (Этап 4.3.1): тождество gather/scatter. ---
add_executable(test_particle_system tests/test_particle_system.cpp)
target_include_directories(test_particle_system PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME particle_system_test COMMAND test_particle_system)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(test_particle_system PRIVATE -Wno-unused-parameter)
endif()

# --- Тест контактных сил (part_0.md, ADR-007): аналитика типов 1-3,
#     баланс импульса, ограничитель шага. Линкуется с src/ext_force.cpp,
#     чтобы покрыть и оркестратор.
add_executable(test_contact tests/test_contact.cpp src/ext_force.cpp)
target_include_directories(test_contact PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME contact_test COMMAND test_contact)

# --- Тест осесимметричной ветки силового расчёта ----------------------------
# Закрепляет соглашение об индексах девиатора при axis_sim && dim == 2
# (sxx = S_rr, syy = S_zz, szz = S_theta) и вытекающие из симметрии тождества:
# бесследовость, равенство S_rr и S_theta, нулевой радиальный источник на
# состоянии одноосного сжатия. Линкует конвейер своими объектниками с
# SPH_DIM=2 — по той же причине, что и e2e-цели ниже (ADR-002).
add_executable(test_axisym tests/test_axisym.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_axisym PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_axisym PRIVATE SPH_DIM=2)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_axisym PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_axisym PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_axisym PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_axisym PRIVATE -fopenmp)
    target_link_options(test_axisym PRIVATE -fopenmp)
endif()
add_test(NAME axisym_test COMMAND test_axisym)

# --- Тест трёхмерной ветки силового расчёта (dim == 3) ----------------------
# По образцу test_axisym, но для ветки dim == 3, которую не проверяет ни
# один сквозной тест: простой сдвиг, сжатие, жёсткое вращение (производная
# Яуманна), покой при однородном и линейном поле напряжений. Нити
# (fiber/PLAN.md) живут только в 3D, и Solid рядом с ними считается этой
# веткой. Свои объектники конвейера с SPH_DIM=3 — по той же причине (ADR-002).
add_executable(test_kinematics_3d tests/test_kinematics_3d.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_kinematics_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_kinematics_3d PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_kinematics_3d PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_kinematics_3d PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_kinematics_3d PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_kinematics_3d PRIVATE -fopenmp)
    target_link_options(test_kinematics_3d PRIVATE -fopenmp)
endif()
add_test(NAME kinematics_3d_test COMMAND test_kinematics_3d)
set_tests_properties(kinematics_3d_test PROPERTIES LABELS "3d")

# --- Модель данных нитей (fiber/PLAN.md, этап Ф2): FiberState, FiberTopology,
#     FiberMaterial, Model, разбор файла ломаных (fiber_input.cpp), масса/l0/
#     l_break, body[] в Solid. Нити только при SPH_DIM=3 (Р9) — собственные
#     объектники конвейера с SPH_DIM=3, по образцу test_kinematics_3d/
#     test_e2e_3d (решение E1, adr-002-dim-constexpr.md).
add_executable(test_fiber_model tests/test_fiber_model.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_model PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_model PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_model PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_model PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_model PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_model PRIVATE -fopenmp)
    target_link_options(test_fiber_model PRIVATE -fopenmp)
endif()
add_test(NAME fiber_model_test COMMAND test_fiber_model)
set_tests_properties(fiber_model_test PROPERTIES LABELS "3d")

# --- Чтение LS-DYNA .k (readk/k-import-plan.md, этап K2): KReader + KDeck
#     (include/k_input/, src/k_reader.cpp - попадает в SPH_SOURCES через GLOB).
#     Слой формата не зависит от размерности; SPH_DIM=1 - как в плане (разд. 5).
add_executable(test_k_reader tests/test_k_reader.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_reader PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_reader PRIVATE
    SPH_DIM=1
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_reader PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_reader PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_reader PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_reader PRIVATE -fopenmp)
    target_link_options(test_k_reader PRIVATE -fopenmp)
endif()
add_test(NAME k_reader_test COMMAND test_k_reader)

# --- Таблица деталей для импорта .k (readk/k-import-plan.md, этап K3, Р2):
#     include/k_input/k_parts.hpp, src/k_parts.cpp. MATERIAL нитей - код .fib,
#     поэтому SPH_DIM=3 (нити - только 3D, Р9).
add_executable(test_k_parts tests/test_k_parts.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_parts PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_parts PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_parts PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_parts PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_parts PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_parts PRIVATE -fopenmp)
    target_link_options(test_k_parts PRIVATE -fopenmp)
endif()
add_test(NAME k_parts_test COMMAND test_k_parts)
set_tests_properties(k_parts_test PROPERTIES LABELS "3d")

# --- Solid-сборщик импорта .k (readk/k-import-plan.md, этап K4): src/k_geometry.cpp.
#     test_k_solid собирается для SPH_DIM 1/2/3 (правило Р6 для каждой dim).

add_executable(test_k_solid_1d tests/test_k_solid.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_solid_1d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_solid_1d PRIVATE SPH_DIM=1)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_solid_1d PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_solid_1d PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_solid_1d PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_solid_1d PRIVATE -fopenmp)
    target_link_options(test_k_solid_1d PRIVATE -fopenmp)
endif()
add_test(NAME k_solid_1d_test COMMAND test_k_solid_1d)

add_executable(test_k_solid_2d tests/test_k_solid.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_solid_2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_solid_2d PRIVATE SPH_DIM=2)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_solid_2d PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_solid_2d PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_solid_2d PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_solid_2d PRIVATE -fopenmp)
    target_link_options(test_k_solid_2d PRIVATE -fopenmp)
endif()
add_test(NAME k_solid_2d_test COMMAND test_k_solid_2d)

add_executable(test_k_solid_3d tests/test_k_solid.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_solid_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_solid_3d PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_solid_3d PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_solid_3d PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_solid_3d PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_solid_3d PRIVATE -fopenmp)
    target_link_options(test_k_solid_3d PRIVATE -fopenmp)
endif()
add_test(NAME k_solid_3d_test COMMAND test_k_solid_3d)
set_tests_properties(k_solid_3d_test PROPERTIES LABELS "3d")

# Эквивалентность старого parse_k_file (копия-эталон в тесте) и нового импорта на 1.k.

add_executable(test_k_legacy_equiv tests/test_k_legacy_equiv.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_legacy_equiv PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_legacy_equiv PRIVATE SPH_DIM=3
    SPH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_legacy_equiv PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_legacy_equiv PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_legacy_equiv PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_legacy_equiv PRIVATE -fopenmp)
    target_link_options(test_k_legacy_equiv PRIVATE -fopenmp)
endif()
add_test(NAME k_legacy_equiv_test COMMAND test_k_legacy_equiv)
set_tests_properties(k_legacy_equiv_test PROPERTIES LABELS "3d")

# --- Нити из *ELEMENT_BEAM (readk/k-import-plan.md, этап K5) и эквивалентность .fib <-> .k.
add_executable(test_k_fiber tests/test_k_fiber.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_fiber PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_fiber PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_fiber PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_fiber PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_fiber PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_fiber PRIVATE -fopenmp)
    target_link_options(test_k_fiber PRIVATE -fopenmp)
endif()
add_test(NAME k_fiber_test COMMAND test_k_fiber)
set_tests_properties(k_fiber_test PROPERTIES LABELS "3d")

# --- Сквозной тест подключения импорта .k к input()/time_integration (этап K6).
add_executable(test_k_e2e tests/test_k_e2e.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_e2e PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_e2e PRIVATE SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_e2e PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_e2e PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_e2e PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_e2e PRIVATE -fopenmp)
    target_link_options(test_k_e2e PRIVATE -fopenmp)
endif()
add_test(NAME k_e2e_test COMMAND test_k_e2e)
set_tests_properties(k_e2e_test PROPERTIES LABELS "3d")

# --- Чекпоинт/рестарт: сцена tests/checkpoint_scene.hpp (Solid + нити из .fib с
#     пластичностью и повреждённостью). Заготовка: тесты будут добавляться в
#     tests/test_checkpoint.cpp. Вторая цель -- та же программа со SPH_ENERGY_LOG.
add_executable(test_checkpoint tests/test_checkpoint.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_checkpoint PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_checkpoint PRIVATE SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_checkpoint PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_checkpoint PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_checkpoint PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_checkpoint PRIVATE -fopenmp)
    target_link_options(test_checkpoint PRIVATE -fopenmp)
endif()
add_test(NAME checkpoint_test COMMAND test_checkpoint)
set_tests_properties(checkpoint_test PROPERTIES LABELS "3d" TIMEOUT 600
    ENVIRONMENT OMP_NUM_THREADS=1)

add_executable(test_checkpoint_energy tests/test_checkpoint.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_checkpoint_energy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_checkpoint_energy PRIVATE SPH_DIM=3
    SPH_ENERGY_LOG
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_checkpoint_energy PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_checkpoint_energy PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_checkpoint_energy PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_checkpoint_energy PRIVATE -fopenmp)
    target_link_options(test_checkpoint_energy PRIVATE -fopenmp)
endif()
add_test(NAME checkpoint_energy_test COMMAND test_checkpoint_energy)
set_tests_properties(checkpoint_energy_test PROPERTIES LABELS "3d" TIMEOUT 600
    ENVIRONMENT OMP_NUM_THREADS=1)

# --- Закрепления из *BOUNDARY_SPC_* (readk/k-import-plan.md, этап K8).
add_executable(test_k_clamp tests/test_k_clamp.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_clamp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_clamp PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_clamp PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_clamp PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_clamp PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_clamp PRIVATE -fopenmp)
    target_link_options(test_k_clamp PRIVATE -fopenmp)
endif()
add_test(NAME k_clamp_test COMMAND test_k_clamp)
set_tests_properties(k_clamp_test PROPERTIES LABELS "3d")

# --- *ELEMENT_SOLID -> частицы (readk/k-import-plan.md, этап K9).
add_executable(test_k_solid_elements tests/test_k_solid_elements.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_k_solid_elements PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_k_solid_elements PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_solid_elements PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_k_solid_elements PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_k_solid_elements PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_k_solid_elements PRIVATE -fopenmp)
    target_link_options(test_k_solid_elements PRIVATE -fopenmp)
endif()
add_test(NAME k_solid_elements_test COMMAND test_k_solid_elements)
set_tests_properties(k_solid_elements_test PROPERTIES LABELS "3d")

# --- Fiber-Bond (fiber/PLAN.md, этап Ф3): primitive/kernel силы растяжения
#     топологических соседей (include/interactions/fiber_bond.hpp), dt_limit
#     связей, предиктор-корректор для fiber::FiberState. Нити — только при
#     SPH_DIM=3 (Р9), по образцу test_fiber_model.
add_executable(test_fiber_bond tests/test_fiber_bond.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_bond PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_bond PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_bond PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_bond PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_bond PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_bond PRIVATE -fopenmp)
    target_link_options(test_fiber_bond PRIVATE -fopenmp)
endif()
add_test(NAME fiber_bond_test COMMAND test_fiber_bond)
set_tests_properties(fiber_bond_test PROPERTIES LABELS "3d")

# --- Пластичность нити, этап 1 (fiber/fiber-plasticity-design.md): чистая
#     функция bond_constitutive_update, тесты Т1-Т7. Без решателя и SPH_DIM.
add_executable(test_fiber_plasticity tests/test_fiber_plasticity.cpp)
target_include_directories(test_fiber_plasticity PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(test_fiber_plasticity PRIVATE -Weverything -Wno-missing-prototypes -Wno-padded
        -Wno-unsafe-buffer-usage -Wno-unsafe-buffer-usage-in-libc-call -Wno-c++98-compat -Wno-c++98-compat-pedantic)
endif()
add_test(NAME fiber_plasticity_test COMMAND test_fiber_plasticity)

# --- Пластичность нити, этап 2 (fiber/fiber-plasticity-design.md): поля
#     состояния FiberBond, RHS-ядро, проход Postprocessing; Т0, Т8-Т11 +
#     полный цикл пластики. По образцу test_fiber_bond.
add_executable(test_fiber_plasticity_integration tests/test_fiber_plasticity_integration.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_plasticity_integration PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_plasticity_integration PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_plasticity_integration PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_plasticity_integration PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_plasticity_integration PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_plasticity_integration PRIVATE -fopenmp)
    target_link_options(test_fiber_plasticity_integration PRIVATE -fopenmp)
endif()
add_test(NAME fiber_plasticity_integration_test COMMAND test_fiber_plasticity_integration)
set_tests_properties(fiber_plasticity_integration_test PROPERTIES LABELS "3d")

# --- Пластичность нити, этап 5 (валидация): Т12 двухволновая картина в
#     стержне из бусин (Карман-Рахматулин), Т13 независимость от дискретизации
#     (10 и 40 бусин, с RATE и без).
add_executable(test_fiber_plasticity_validation tests/test_fiber_plasticity_validation.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_plasticity_validation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_plasticity_validation PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_plasticity_validation PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_plasticity_validation PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_plasticity_validation PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_plasticity_validation PRIVATE -fopenmp)
    target_link_options(test_fiber_plasticity_validation PRIVATE -fopenmp)
endif()
add_test(NAME fiber_plasticity_validation_test COMMAND test_fiber_plasticity_validation)
set_tests_properties(fiber_plasticity_validation_test PROPERTIES LABELS "3d" TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=2")

# --- Т16: сцена Guo 2019, нить упругая против упруго-пластичной (smoke для
#     ctest; без аргументов). С аргументами (--level/--config/--v0/--t-max/--out)
#     -- длинный прогон с таблицей summary.txt (experiments/...-fiber-plasticity-guo).
add_executable(test_fiber_plasticity_guo tests/test_fiber_plasticity_guo.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_plasticity_guo PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/tests)
target_compile_definitions(test_fiber_plasticity_guo PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_plasticity_guo PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_plasticity_guo PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_plasticity_guo PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_plasticity_guo PRIVATE -fopenmp)
    target_link_options(test_fiber_plasticity_guo PRIVATE -fopenmp)
endif()
add_test(NAME fiber_plasticity_guo_smoke COMMAND test_fiber_plasticity_guo)
set_tests_properties(fiber_plasticity_guo_smoke PROPERTIES LABELS "3d;e2e" TIMEOUT 240
    ENVIRONMENT "OMP_NUM_THREADS=2")

# --- Т15: директивы неупругости связи в .fib (PLASTIC, PLASTIC_TABLE, RATE,
#     VISCO, FAILURE), ошибки разбора и колонки eps_p/sigma в выводе BONDS.
add_executable(test_fiber_input_inelastic tests/test_fiber_input_inelastic.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_input_inelastic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_input_inelastic PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_input_inelastic PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_input_inelastic PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_input_inelastic PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_input_inelastic PRIVATE -fopenmp)
    target_link_options(test_fiber_input_inelastic PRIVATE -fopenmp)
endif()
add_test(NAME fiber_input_inelastic_test COMMAND test_fiber_input_inelastic)
set_tests_properties(fiber_input_inelastic_test PROPERTIES LABELS "3d")

# --- ODESolver (fiber/PLAN.md, этап Ф8): интерфейс OdeSolver, PredictorCorrector
#     и Euler (include/ode/ode_solver.hpp) -- первый порядок сходимости Эйлера
#     на задаче с аналитическим решением, ограниченность энергии нити без
#     контактов на PredictorCorrector. Нити -- только при SPH_DIM=3 (Р9), по
#     образцу test_fiber_bond.
add_executable(test_ode_solver tests/test_ode_solver.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_ode_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_ode_solver PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_ode_solver PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_ode_solver PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_ode_solver PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_ode_solver PRIVATE -fopenmp)
    target_link_options(test_ode_solver PRIVATE -fopenmp)
endif()
add_test(NAME ode_solver_test COMMAND test_ode_solver)
set_tests_properties(ode_solver_test PROPERTIES LABELS "3d")

# --- Postprocessing: разрыв связей Fiber-Bond (fiber/PLAN.md, этап Ф7):
#     каркас "поиск событий -> применение" (include/interactions/
#     postprocessing.hpp), критерий l > l_break, новый object_id
#     отделившейся части, взаимодействие через Fiber-Fiber-Contact после
#     разрыва. Нити -- только при SPH_DIM=3 (Р9), по образцу test_fiber_bond.
add_executable(test_fiber_break tests/test_fiber_break.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_break PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_break PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_break PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_break PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_break PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_break PRIVATE -fopenmp)
    target_link_options(test_fiber_break PRIVATE -fopenmp)
endif()
add_test(NAME fiber_break_test COMMAND test_fiber_break)
set_tests_properties(fiber_break_test PROPERTIES LABELS "3d")

# --- Fiber-Fiber-Contact (fiber/PLAN.md, этап Ф5): primitive/kernel
#     пространственного контакта бусин разных нитей
#     (include/interactions/fiber_contact.hpp), dt_limit контакта,
#     диагностика "пар на пересечение". Нити — только при SPH_DIM=3 (Р9), по
#     образцу test_fiber_bond.
add_executable(test_fiber_contact tests/test_fiber_contact.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_contact PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_contact PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_contact PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_contact PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_contact PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_contact PRIVATE -fopenmp)
    target_link_options(test_fiber_contact PRIVATE -fopenmp)
endif()
add_test(NAME fiber_contact_test COMMAND test_fiber_contact)
set_tests_properties(fiber_contact_test PROPERTIES LABELS "3d")

# --- Трение нить-нить: проекция на уровне скоростей (стадия между dt и
#     advance(); include/interactions/fiber_friction.hpp,
#     fiber/fiber_friction_projection.md, §8). Нити -- только SPH_DIM=3.
add_executable(test_fiber_friction tests/test_fiber_friction.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_friction PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_friction PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_friction PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_friction PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_friction PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_friction PRIVATE -fopenmp)
    target_link_options(test_fiber_friction PRIVATE -fopenmp)
endif()
add_test(NAME fiber_friction_test COMMAND test_fiber_friction)
set_tests_properties(fiber_friction_test PROPERTIES LABELS "3d")

# --- Кинематическое закрепление Solid/Fiber (include/interactions/clamp.hpp,
#     fiber/TZ_CLAMP_SOLID_FIBER.md, §10). Регистрируется при 1, 2 и 4 потоках.
add_executable(test_clamp tests/test_clamp.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_clamp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_clamp PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_clamp PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_clamp PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_clamp PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_clamp PRIVATE -fopenmp)
    target_link_options(test_clamp PRIVATE -fopenmp)
endif()
add_test(NAME clamp_test_t1 COMMAND test_clamp)
set_tests_properties(clamp_test_t1 PROPERTIES LABELS "3d" ENVIRONMENT "OMP_NUM_THREADS=1")
add_test(NAME clamp_test_t2 COMMAND test_clamp)
set_tests_properties(clamp_test_t2 PROPERTIES LABELS "3d" ENVIRONMENT "OMP_NUM_THREADS=2")
add_test(NAME clamp_test_t4 COMMAND test_clamp)
set_tests_properties(clamp_test_t4 PROPERTIES LABELS "3d" ENVIRONMENT "OMP_NUM_THREADS=4")

# --- E2E-критерий готовности Ф5 (fiber/PLAN.md: «две нити, летящие навстречу
#     крест-накрест, отскакивают без взаимного проникания»).
add_executable(test_e2e_fiber_contact tests/test_e2e_fiber_contact.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_fiber_contact PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_e2e_fiber_contact PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_fiber_contact PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_fiber_contact PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_fiber_contact PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_fiber_contact PRIVATE -fopenmp)
    target_link_options(test_e2e_fiber_contact PRIVATE -fopenmp)
endif()
add_test(NAME e2e_fiber_contact_test COMMAND test_e2e_fiber_contact)
set_tests_properties(e2e_fiber_contact_test PROPERTIES LABELS "3d")

# --- Fiber-Solid-Contact (fiber/PLAN.md, этап Ф6): primitive/kernel
#     механического контакта бусины нити с реальной Solid-частицей
#     (include/interactions/fiber_solid_contact.hpp), elastic_by_itype
#     (include/materials.hpp), dt_limit контакта. Нити — только при
#     SPH_DIM=3 (Р9), по образцу test_fiber_contact.
add_executable(test_fiber_solid tests/test_fiber_solid.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_solid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_solid PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_solid PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_solid PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_solid PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_solid PRIVATE -fopenmp)
    target_link_options(test_fiber_solid PRIVATE -fopenmp)
endif()
add_test(NAME fiber_solid_test COMMAND test_fiber_solid)
set_tests_properties(fiber_solid_test PROPERTIES LABELS "3d")

# --- Законы pinball/Swegle/импульсный для Fiber-Fiber (fiber/plan-fiber-contact-types.md, этап 2):
#     ядро fiber_contact.hpp, предел шага и диагностика.
add_executable(test_fiber_contact_types tests/test_fiber_contact_types.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_contact_types PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_contact_types PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_contact_types PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_contact_types PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_contact_types PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_contact_types PRIVATE -fopenmp)
    target_link_options(test_fiber_contact_types PRIVATE -fopenmp)
endif()
add_test(NAME fiber_contact_types_test COMMAND test_fiber_contact_types)
set_tests_properties(fiber_contact_types_test PROPERTIES LABELS "3d")

# --- Законы pinball/Swegle/импульсный для Fiber-Solid (fiber/plan-fiber-contact-types.md, этап 1):
#     include/interactions/fiber_contact_law.hpp, ядро fiber_solid_contact.hpp.
add_executable(test_fiber_solid_types tests/test_fiber_solid_types.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_solid_types PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_solid_types PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_solid_types PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_solid_types PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_solid_types PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_solid_types PRIVATE -fopenmp)
    target_link_options(test_fiber_solid_types PRIVATE -fopenmp)
endif()
add_test(NAME fiber_solid_types_test COMMAND test_fiber_solid_types)
set_tests_properties(fiber_solid_types_test PROPERTIES LABELS "3d")

# --- E2E-критерий готовности Ф6 (fiber/PLAN.md: «smoke одного Solid в 3D,
#     затем нить, падающая на пластину, без проникания»).
add_executable(test_e2e_fiber_solid tests/test_e2e_fiber_solid.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_fiber_solid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_e2e_fiber_solid PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_fiber_solid PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_fiber_solid PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_fiber_solid PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_fiber_solid PRIVATE -fopenmp)
    target_link_options(test_e2e_fiber_solid PRIVATE -fopenmp)
endif()
add_test(NAME e2e_fiber_solid_test COMMAND test_e2e_fiber_solid)
set_tests_properties(e2e_fiber_solid_test PROPERTIES LABELS "3d" TIMEOUT 300)

# --- Опция ядра sph_same_body_only: изоляция SPH-пар Solid по телам
#     (ТЗ эксперименты Guo 2019, §4.2).
add_executable(test_sph_body_isolation tests/test_sph_body_isolation.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_sph_body_isolation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_sph_body_isolation PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_sph_body_isolation PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_sph_body_isolation PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_sph_body_isolation PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_sph_body_isolation PRIVATE -fopenmp)
    target_link_options(test_sph_body_isolation PRIVATE -fopenmp)
endif()
add_test(NAME sph_body_isolation_test COMMAND test_sph_body_isolation)
set_tests_properties(sph_body_isolation_test PROPERTIES LABELS "3d" TIMEOUT 120)

# --- Самопроверка читателя кадров нитей (tests/e2e_fiber_frames.hpp).
add_executable(test_e2e_fiber_frames tests/test_e2e_fiber_frames.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_fiber_frames PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_e2e_fiber_frames PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_fiber_frames PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_fiber_frames PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_fiber_frames PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_fiber_frames PRIVATE -fopenmp)
    target_link_options(test_e2e_fiber_frames PRIVATE -fopenmp)
endif()
add_test(NAME e2e_fiber_frames_selfcheck COMMAND test_e2e_fiber_frames)
set_tests_properties(e2e_fiber_frames_selfcheck PROPERTIES LABELS "3d;e2e" TIMEOUT 60)

# --- Smoke-тест сцен Guo 2019 (ТЗ Guo 2019 §9): шар + ткань + Al-пластина,
#     по одной регистрации на конфигурацию мишени.
add_executable(test_e2e_guo2019 tests/test_e2e_guo2019.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_guo2019 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_e2e_guo2019 PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_guo2019 PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_guo2019 PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_guo2019 PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_guo2019 PRIVATE -fopenmp)
    target_link_options(test_e2e_guo2019 PRIVATE -fopenmp)
endif()
foreach(_guo_cfg A 1T 1T_A A_1T 2T 4T)
    add_test(NAME e2e_guo2019_smoke_${_guo_cfg}
        COMMAND test_e2e_guo2019 --config=${_guo_cfg} --time-limit=100)
    # 2 потока: при ctest -jN 8 потоков на тест переподписывают ядра, и
    # smoke упирается в --time-limit (замер: 1T_A 8 с при 8 потоках один,
    # >100 с вчетвером при -j4; 15 с при 2 потоках).
    set_tests_properties(e2e_guo2019_smoke_${_guo_cfg} PROPERTIES LABELS "3d;e2e" TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=2")
endforeach()

# --- Т0 пластичности нити: дамп состояния нитей с %.17g для побитового
#     сравнения до/после (tests/fiber_t0_dump.cpp), без add_test. Инструмент:
#     EXCLUDE_FROM_ALL, собирается явно (--target fiber_t0_dump).
add_executable(fiber_t0_dump EXCLUDE_FROM_ALL tests/fiber_t0_dump.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(fiber_t0_dump PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(fiber_t0_dump PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(fiber_t0_dump PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(fiber_t0_dump PRIVATE -fopenmp)
    target_link_options(fiber_t0_dump PRIVATE -fopenmp)
endif()

# --- Лог механической энергии (include/energy_log.hpp) на smoke Guo 2019:
#     формат CSV, сверка запасов, замыкание баланса, лог не возмущает расчёт.
#     Собирается с SPH_ENERGY_LOG всегда, независимо от опции.
add_executable(test_energy_log tests/test_energy_log.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_energy_log PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_energy_log PRIVATE
    SPH_DIM=3
    SPH_ENERGY_LOG
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_energy_log PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_energy_log PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_energy_log PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_energy_log PRIVATE -fopenmp)
    target_link_options(test_energy_log PRIVATE -fopenmp)
endif()
add_test(NAME energy_log_test COMMAND test_energy_log --time-limit=200)
set_tests_properties(energy_log_test PROPERTIES LABELS "3d;e2e" TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=2")

# --- Программа серийных выстрелов Guo 2019 (ТЗ Guo 2019 §8), без add_test.
add_executable(guo2019_shot tools/guo2019/guo2019_shot.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(guo2019_shot PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/tests)
target_compile_definitions(guo2019_shot PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(guo2019_shot PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(guo2019_shot PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(guo2019_shot PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(guo2019_shot PRIVATE -fopenmp)
    target_link_options(guo2019_shot PRIVATE -fopenmp)
endif()

# --- CLI чекпоинтов guo2019_shot: коды ошибок и сравнение раздельного прогона с непрерывным
#     (tests/checkpoint_cli_test.sh; требует собранный guo2019_shot).
add_test(NAME checkpoint_cli_test COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/checkpoint_cli_test.sh $<TARGET_FILE:guo2019_shot>)
set_tests_properties(checkpoint_cli_test PROPERTIES LABELS "3d" TIMEOUT 600 ENVIRONMENT OMP_NUM_THREADS=1)

# --- Проверочные сцены модели ткани Guo 2019 (шаг 2 диагностики): растяжение полоски, капстан, колебания
#     лоскута (tools/fabric_checks/fabric_checks.cpp). Не тест; add_test -- только короткий smoke (ниже).
add_executable(fabric_checks tools/fabric_checks/fabric_checks.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(fabric_checks PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/tests)
target_compile_definitions(fabric_checks PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(fabric_checks PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(fabric_checks PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(fabric_checks PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(fabric_checks PRIVATE -fopenmp)
    target_link_options(fabric_checks PRIVATE -fopenmp)
endif()
# Короткие smoke-прогоны (секунды): полоска прямых нитей и капстан на 10 мкс -- сцены собираются и считаются без исключений.
add_test(NAME fabric_checks_smoke_strip
    COMMAND fabric_checks strip --kind=straight --material=elastic --rate=300 --eps-max=0.004
            --out=${CMAKE_CURRENT_BINARY_DIR}/fabric_checks_smoke --tag=strip --threads=2)
add_test(NAME fabric_checks_smoke_capstan
    COMMAND fabric_checks capstan --drum=solid --theta-deg=90 --mu=0.3 --t-settle=1e-5 --t-end=1e-5
            --out=${CMAKE_CURRENT_BINARY_DIR}/fabric_checks_smoke --tag=capstan --threads=2)
set_tests_properties(fabric_checks_smoke_strip fabric_checks_smoke_capstan PROPERTIES LABELS "3d;e2e" TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=2")

# --- Генератор полотняного переплетения (weave/TZ_PLAIN_WEAVE_GENERATOR.md):
#     include/fiber/plain_weave.hpp (header-only), проверки геометрии.
add_executable(test_plain_weave tests/test_plain_weave.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_plain_weave PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_plain_weave PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_plain_weave PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_plain_weave PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_plain_weave PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_plain_weave PRIVATE -fopenmp)
    target_link_options(test_plain_weave PRIVATE -fopenmp)
endif()
add_test(NAME plain_weave_test COMMAND test_plain_weave)
set_tests_properties(plain_weave_test PROPERTIES LABELS "3d" TIMEOUT 120)

# --- CLI генератор полотняного переплетения (weave/TZ_PLAIN_WEAVE_GENERATOR.md, §4.2):
#     инструмент для создания геометрии с командной строки.
add_executable(plain_weave_gen tools/plain_weave/plain_weave_gen.cpp)
target_include_directories(plain_weave_gen PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(plain_weave_gen PRIVATE -fexceptions)
endif()

# CLI-тест: cmake -P скрипт создаёт уникальный каталог в системном temp и удаляет его.
add_test(NAME plain_weave_gen_cli_test
         COMMAND ${CMAKE_COMMAND} -DGEN=$<TARGET_FILE:plain_weave_gen>
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/plain_weave_cli_test.cmake)
set_tests_properties(plain_weave_gen_cli_test PROPERTIES LABELS "3d" TIMEOUT 60)

# --- Ткань-«одеяло» из патчей и сетка стёжки слоёв (include/fiber/plain_weave_quilt.hpp).
add_executable(test_plain_weave_quilt tests/test_plain_weave_quilt.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_plain_weave_quilt PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_plain_weave_quilt PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_plain_weave_quilt PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_plain_weave_quilt PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_plain_weave_quilt PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_plain_weave_quilt PRIVATE -fopenmp)
    target_link_options(test_plain_weave_quilt PRIVATE -fopenmp)
endif()
add_test(NAME plain_weave_quilt_test COMMAND test_plain_weave_quilt)
set_tests_properties(plain_weave_quilt_test PROPERTIES LABELS "3d" TIMEOUT 120)

# --- Эквивалентность Solid-цикла в новой перегрузке time_integration
#     (fiber::Model&) со старой (ParticleSystem&) при пустых нитях
#     (fiber/PLAN.md, этап Ф3). Новая перегрузка не переиспользует старый
#     цикл, а содержит его построчную копию (src/time_integration.cpp,
#     TODO(Ф8)); до этого теста копия не была упражнена НИ ОДНИМ тестом с
#     непустым Solid (test_fiber_bond.cpp проверяет только Fiber-физику при
#     ntotal == 0). Кейс — уменьшенная 3D-геометрия BC2
#     (tests/bc2_mini3d_geometry.hpp), та же, что использует test_e2e_3d.
add_executable(test_fiber_time_integration_equivalence
    tests/test_fiber_time_integration_equivalence.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_time_integration_equivalence PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_time_integration_equivalence PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_time_integration_equivalence PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_time_integration_equivalence PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_time_integration_equivalence PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_time_integration_equivalence PRIVATE -fopenmp)
    target_link_options(test_fiber_time_integration_equivalence PRIVATE -fopenmp)
endif()
add_test(NAME fiber_time_integration_equivalence_test COMMAND test_fiber_time_integration_equivalence)
set_tests_properties(fiber_time_integration_equivalence_test PROPERTIES LABELS "3d" TIMEOUT 300)

# --- Границы SPH_DIM != 3 для нитей (fiber/PLAN.md Ф2, пункт 5):
#     fiber::check_fiber_requires_dim3() обязана бросать исключение. dim —
#     compile-time constexpr (ADR-002), поэтому ветку dim != 3 нельзя
#     исполнить в цели, собранной с SPH_DIM=3 (там dim == 3 всегда,
#     см. test_fiber_model выше). Отдельная лёгкая цель БЕЗ SPH_DIM
#     (dim = 1 по умолчанию, param.hpp) проверяет именно эту ветку;
#     check_fiber_requires_dim3() — header-only, конвейер не нужен.
add_executable(test_fiber_dim_guard tests/test_fiber_dim_guard.cpp)
target_include_directories(test_fiber_dim_guard PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME fiber_dim_guard_test COMMAND test_fiber_dim_guard)

# --- End-to-end тест на кейсе BC2 (contact/04-e2e-test-brief.md) ------------
# Цель линкует ВЕСЬ конвейер (src/*.cpp кроме main.cpp) и компилирует его
# СВОИМИ объектниками с SPH_DIM=2: объектники цели SPH собраны с dim=1 и
# несовместимы по ODR (layout ParticleSystem зависит от dim).
# См. решение E1 и docs/04-explanation/decisions/adr-002-dim-constexpr.md.
option(SPH_E2E_SLOW_TESTS "Регистрировать в ctest медленные полные прогоны BC2 и BC6 (метка slow)" OFF)

add_executable(test_e2e_bc2 tests/test_e2e_bc2.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_bc2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# SPH_TEST_DATA_DIR — АБСОЛЮТНЫЙ путь от каталога исходников, а не от рабочего:
# тест запускается из каталога сборки (ctest) и из временного каталога прогона.
target_compile_definitions(test_e2e_bc2 PRIVATE
    SPH_DIM=2
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_bc2 PRIVATE -fexceptions)
endif()
# Конвейер использует OpenMP (interactions/solid_sph.hpp, rhs/cpu_executor.hpp).
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_bc2 PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_bc2 PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_bc2 PRIVATE -fopenmp)
    target_link_options(test_e2e_bc2 PRIVATE -fopenmp)
endif()

# Smoke-режим (ТЗ §5): в прогоне ctest по умолчанию, цель по времени < 30 с.
# Два рубежа защиты от зависания, именно в таком порядке:
#   --time-limit — внутренний сторож (run_guard.hpp). Срабатывает ПЕРВЫМ и даёт
#                  разбор: не хватило времени или расчёт перестал продвигаться.
#   TIMEOUT      — жёсткий рубеж CTest на случай, если откажет сам сторож
#                  (например, зависание внутри одного шага, куда он не заглядывает).
#
# Регистрируется ПО ТЕСТУ НА ТИП КОНТАКТА (0..3). Прежде один тест делал
# базовый прогон и следом сценарии §6: обрыв любого из них прекращал разбор
# остальных, а в отчёте ctest всё сводилось к одной строке «провален», по
# которой нельзя сказать, какой именно тип сломался. Отдельные тесты падают
# независимо, идут параллельно под `ctest -j` и называют виновника именем.
foreach(_type RANGE 0 3)
    add_test(NAME e2e_bc2_smoke_contact${_type}
             COMMAND test_e2e_bc2 --contact-type=${_type} --time-limit=120)
    set_tests_properties(e2e_bc2_smoke_contact${_type} PROPERTIES LABELS "e2e" TIMEOUT 300)
endforeach()

# Бюджет сторожа на ОДИН полный прогон BC2, секунды. Общий для тестов полного
# режима и для цели bc2_plot, чтобы они не разъезжались.
set(SPH_BC2_RUN_TIME_LIMIT 900 CACHE STRING
    "Бюджет сторожа на один полный прогон BC2, с (0 — без ограничения)")

# Полный режим (ТЗ §5): минуты, поэтому по умолчанию НЕ регистрируется.
# Включение:  cmake -S . -B build -DSPH_E2E_SLOW_TESTS=ON && ctest --test-dir build -L slow
# (голая метка LABELS "slow" из ТЗ не исключает тест из прогона по умолчанию —
#  CTest не умеет «зарегистрирован, но по умолчанию не запускается», поэтому
#  регистрация гейтится опцией сборки.)
if(SPH_E2E_SLOW_TESTS)
    # По тесту на тип контакта, как и в smoke. Теперь каждый тест — РОВНО ОДИН
    # прогон, поэтому жёсткий рубеж CTest покрывает один бюджет сторожа с
    # запасом (300 с сторожа < 900 с TIMEOUT) и обещанный порядок сохраняется:
    # сторож срабатывает первым и говорит, чего не хватило.
    # Бюджет сторожа. 300 с, стоявшие здесь раньше, полному прогону НЕ хватает:
    # на машине разработки он идёт около 4,3 шага/с, то есть 1531 шагу нужно
    # ~355 с, а сценарию типа 2 до обрыва — ~390 с. Прежде это было незаметно,
    # потому что провал списывался на сам сценарий; сторож же честно писал
    # «шаг здоров — не хватило именно времени».
    #
    # Поднятие бюджета не ослабляет защиту: коллапс шага ловит ВТОРОЙ механизм
    # того же сторожа — детектор застоя (потолок оценки оставшихся шагов и
    # терпение в 200 шагов, run_guard.hpp), и он от стенных часов не зависит.
    # Бюджет времени — это backstop на случай зависания, и он обязан покрывать
    # самую медленную поддерживаемую машину, а не калиброваться под быструю.
    foreach(_type RANGE 0 3)
        add_test(NAME e2e_bc2_full_contact${_type}
                 COMMAND test_e2e_bc2 --full --contact-type=${_type}
                         --time-limit=${SPH_BC2_RUN_TIME_LIMIT})
        set_tests_properties(e2e_bc2_full_contact${_type}
                             PROPERTIES LABELS "slow;e2e" TIMEOUT 1800)
    endforeach()
endif()

# Перестройка графика для отчёта contact/06-e2e-case.md:
#     cmake --build build --target bc2_plot
# ЧЕТЫРЕ полных прогона (по одному на contact_force_type), каждый выгружает
# свою кривую, затем пятый шаг — быстрый, без расчёта — собирает из выгрузок
# SVG с пятью кривыми: эталон и типы 0..3. В файл выгружается ровно тот объект
# Curve, который проверил компаратор, поэтому график по-прежнему не может
# разойтись с приведёнными в отчёте метриками.
#
# Почему отдельная цель, а не часть прогона тестов: запись в дерево исходников —
# побочный эффект. Тест обязан его не иметь (он это сам и проверяет, §7
# check_no_side_effects), обязан работать на checkout только для чтения и
# допускать параллельный ctest. Поэтому --emit-plot выключен по умолчанию,
# ни одна зарегистрированная ctest-цель его не передаёт, а обновление картинки
# остаётся осознанным действием владельца.
set(SPH_BC2_PLOT "${CMAKE_CURRENT_SOURCE_DIR}/contact/e2e/bc2_monitor_velocity.svg"
    CACHE FILEPATH "Куда цель bc2_plot записывает график BC2")
# Прогон идёт через cmake -P: цель должна выжить при неудаче теста по существу
# расчёта, если график при этом записан (подробности — в cmake/bc2_plot.cmake).
add_custom_target(bc2_plot
    COMMAND ${CMAKE_COMMAND}
            -DSPH_PLOT_EXE=$<TARGET_FILE:test_e2e_bc2>
            "-DSPH_PLOT_OUT=${SPH_BC2_PLOT}"
            -DSPH_PLOT_TIME_LIMIT=${SPH_BC2_RUN_TIME_LIMIT}
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/bc2_plot.cmake"
    DEPENDS test_e2e_bc2
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Полный прогон BC2 и перестройка ${SPH_BC2_PLOT}"
    USES_TERMINAL VERBATIM)

# --- End-to-end тест на кейсе BC6 (contact/06-e2e-bc6.md) -------------------
# Откольный опыт на B4C: ударник и окно PMMA, подложка — пена 160. Схема та
# же, что у BC2 (свои объектники конвейера с SPH_DIM=2 по решению E1, по тесту
# на contact_force_type, smoke по умолчанию, полный режим под опцией), но
# кейс независим: свои данные (tests/bc6_case.hpp), своя геометрия
# (tests/bc6_geometry.hpp), свой конфиг и свой эталон в tests/data/.
# Общая обвязка вынесена в tests/e2e_*.hpp; tests/test_e2e_bc2.cpp НЕ изменён
# и продолжает жить со своими копиями.
add_executable(test_e2e_bc6 tests/test_e2e_bc6.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_bc6 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_e2e_bc6 PRIVATE
    SPH_DIM=2
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_bc6 PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_bc6 PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_bc6 PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_bc6 PRIVATE -fopenmp)
    target_link_options(test_e2e_bc6 PRIVATE -fopenmp)
endif()

# Smoke-режим: delta_x = 8e-4, t_total = 1e-6. Физические проверки там
# пропускаются с явным предупреждением — при таком шаге носитель ядра
# (2h = 3.2 мм) шире зазора удара 1.0 мм, тела взаимодействуют с t = 0.
# Smoke проверяет живучесть конвейера, геометрию и самопроверки инструментов.
foreach(_type RANGE 0 3)
    add_test(NAME e2e_bc6_smoke_contact${_type}
             COMMAND test_e2e_bc6 --contact-type=${_type} --time-limit=120)
    set_tests_properties(e2e_bc6_smoke_contact${_type} PROPERTIES LABELS "e2e" TIMEOUT 300)
endforeach()

# Бюджет сторожа на ОДИН полный прогон BC6, секунды. Общий для тестов полного
# режима и для цели bc6_plot, чтобы они не разъезжались.
#
# Значение поставлено ПО ИЗМЕРЕНИЮ, а не по оценке, и измерение это пришлось
# делать дважды. Первая попытка исходила из скорости BC2 и дала 1800 с; прогон
# не уложился, пройдя 1216 шагов из ~1590. Но виноват был не кейс: замер шёл на
# сборке Debug (-O0). На Release (-O3) полный прогон занимает 471 с и проходит
# все 1614 шагов, то есть 3.4 шага/с.
#
# 1800 — примерно четырёхкратный запас над измеренными 471 с. Это НЕ оценка
# длительности, а backstop на случай зависания, и он обязан покрывать самую
# медленную поддерживаемую машину, включая сборку без оптимизации. Коллапс шага
# ловит второй механизм сторожа (детектор застоя), от стенных часов не
# зависящий, поэтому запас по времени защиту не ослабляет.
#
# Отдельно стоит знать, куда уходит время: single_step на каждом шаге выделяет
# и обнуляет pair_i/pair_j/w/dwdx размером max_interaction = 100*maxn
# (src/single_step.cpp:114-118). Для BC6 это 250 МБ на шаг при действительном
# числе пар примерно в 13 раз меньше ёмкости. Переделка — отдельная задача.
set(SPH_BC6_RUN_TIME_LIMIT 1800 CACHE STRING
    "Бюджет сторожа на один полный прогон BC6, с (0 — без ограничения)")

if(SPH_E2E_SLOW_TESTS)
    foreach(_type RANGE 0 3)
        add_test(NAME e2e_bc6_full_contact${_type}
                 COMMAND test_e2e_bc6 --full --contact-type=${_type}
                         --time-limit=${SPH_BC6_RUN_TIME_LIMIT})
        set_tests_properties(e2e_bc6_full_contact${_type}
                             PROPERTIES LABELS "slow;e2e" TIMEOUT 3600)
    endforeach()
endif()

# Перестройка графика для отчёта contact/06-e2e-bc6.md:
#     cmake --build build --target bc6_plot
# Устройство и обоснование — те же, что у bc2_plot выше.
set(SPH_BC6_PLOT "${CMAKE_CURRENT_SOURCE_DIR}/contact/e2e/bc6_monitor_velocity.svg"
    CACHE FILEPATH "Куда цель bc6_plot записывает график BC6")
add_custom_target(bc6_plot
    COMMAND ${CMAKE_COMMAND}
            -DSPH_PLOT_EXE=$<TARGET_FILE:test_e2e_bc6>
            "-DSPH_PLOT_OUT=${SPH_BC6_PLOT}"
            -DSPH_PLOT_TIME_LIMIT=${SPH_BC6_RUN_TIME_LIMIT}
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/bc6_plot.cmake"
    DEPENDS test_e2e_bc6
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Полный прогон BC6 и перестройка ${SPH_BC6_PLOT}"
    USES_TERMINAL VERBATIM)

# --- Сквозной 3D-тест на уменьшенном BC2 (tests/test_e2e_3d.cpp) -----------
# Материалы, толщины и скорости BC2, радиусы уменьшены до 4 мм у керамики
# (tests/bc2_mini3d_geometry.hpp). Проверяет ветку dim = 3 целиком: здоровье,
# сохранение массы, импульса и энергии, симметрию решётки, приход волны.
# Свои объектники конвейера с SPH_DIM=3 (решение E1).
add_executable(test_e2e_3d tests/test_e2e_3d.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_e2e_3d PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_3d PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_3d PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_3d PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_3d PRIVATE -fopenmp)
    target_link_options(test_e2e_3d PRIVATE -fopenmp)
endif()
# Smoke: ровно 100 шагов, по тесту на тип контакта.
foreach(_type RANGE 0 3)
    add_test(NAME e2e_3d_smoke_contact${_type}
             COMMAND test_e2e_3d --contact-type=${_type} --time-limit=240)
    set_tests_properties(e2e_3d_smoke_contact${_type} PROPERTIES LABELS "e2e;3d" TIMEOUT 300)
endforeach()
# Баланс энергии без XSPH: поправка av прибавляется к скорости и уносит
# энергию, поэтому сохранение полной энергии проверяется при epsilon = 0.
add_test(NAME e2e_3d_energy_noxsph COMMAND test_e2e_3d --no-xsph --time-limit=240)
set_tests_properties(e2e_3d_energy_noxsph PROPERTIES LABELS "e2e;3d" TIMEOUT 300)
# Полный прогон до t_total (3.5 мкс) — под той же опцией, что полные BC2/BC6.
if(SPH_E2E_SLOW_TESTS)
    foreach(_type RANGE 0 3)
        add_test(NAME e2e_3d_full_contact${_type}
                 COMMAND test_e2e_3d --full --contact-type=${_type} --time-limit=1800)
        set_tests_properties(e2e_3d_full_contact${_type}
                             PROPERTIES LABELS "slow;e2e;3d" TIMEOUT 3600)
    endforeach()
endif()

# --- E2E-тесты закрепления Solid и Fiber (tests/test_e2e_clamp.cpp) ---------
# Постановка -- Clamp/CLAMP_TESTS.md: сцены с аналитическим ответом (реакция
# опоры, отражение волны, частота), только через боевой time_integration.
# Один запуск -- один вариант (--case, --solver). Прогоны при OMP_NUM_THREADS=1
# (остальной цикл не обязан быть побитово независим от числа потоков); ES1 и
# EF1 -- дополнительно при 4 потоках с теми же допусками (§1.1).
# ES4 (жёсткая стенка с контактом Solid-Solid) не регистрируется: SPH-суммы
# Solid не изолируют тела по body[], постановка откладывает ES4 до изоляции.
add_executable(test_e2e_clamp tests/test_e2e_clamp.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_clamp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_e2e_clamp PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_clamp PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_clamp PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_clamp PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_clamp PRIVATE -fopenmp)
    target_link_options(test_e2e_clamp PRIVATE -fopenmp)
endif()

# Оба решателя: Solid-сцены, EF1, EF2 (§1.1). Только PC: EF3-EF6, EK1.
foreach(_case ES1 ES2 ES3 EF1 EF2)
    foreach(_solver pc euler)
        add_test(NAME e2e_clamp_${_case}_${_solver}_smoke
                 COMMAND test_e2e_clamp --case=${_case} --solver=${_solver} --time-limit=240)
        set_tests_properties(e2e_clamp_${_case}_${_solver}_smoke PROPERTIES
                             LABELS "e2e;3d" TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1")
    endforeach()
endforeach()
foreach(_case EF3 EF4 EF5 EF6 EK1)
    add_test(NAME e2e_clamp_${_case}_pc_smoke
             COMMAND test_e2e_clamp --case=${_case} --solver=pc --time-limit=240)
    set_tests_properties(e2e_clamp_${_case}_pc_smoke PROPERTIES
                         LABELS "e2e;3d" TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1")
endforeach()
# EF1-b: разрыв у кромки зажима (вариант EF1 с v0 = 1.5 eps_break c), только PC.
add_test(NAME e2e_clamp_EF1b_pc_smoke
         COMMAND test_e2e_clamp --case=EF1 --variant=b --solver=pc --time-limit=240)
set_tests_properties(e2e_clamp_EF1b_pc_smoke PROPERTIES
                     LABELS "e2e;3d" TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1")
foreach(_case ES1 EF1)
    add_test(NAME e2e_clamp_${_case}_pc_smoke_t4
             COMMAND test_e2e_clamp --case=${_case} --solver=pc --time-limit=240)
    set_tests_properties(e2e_clamp_${_case}_pc_smoke_t4 PROPERTIES
                         LABELS "e2e;3d" TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4")
endforeach()
if(SPH_E2E_SLOW_TESTS)
    foreach(_case ES1 ES2 ES3 EF1 EF2 EF3)
        add_test(NAME e2e_clamp_${_case}_pc_full
                 COMMAND test_e2e_clamp --case=${_case} --solver=pc --full --time-limit=3000)
        set_tests_properties(e2e_clamp_${_case}_pc_full PROPERTIES
                             LABELS "slow;e2e;3d" TIMEOUT 3600 ENVIRONMENT "OMP_NUM_THREADS=1")
    endforeach()
endif()

# --- Регрессионный тест на кейсах BC2 и BC6 (tests/test_regression.cpp) ----
# Сквозные тесты выше сверяют расчёт с ПРИБОРНЫМ эталоном в грубых воротах и
# ловят только кратную поломку. Регрессия сверяет все поля каждой частицы после
# N = 1 и N = K шагов с эталоном, снятым этим же кодом (tests/data/regression/),
# и ловит любое изменение расчёта сверх допуска. Сборка — как у e2e: свои
# объектники конвейера с SPH_DIM=2 (решение E1).
add_executable(test_regression tests/test_regression.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_regression PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_regression PRIVATE
    SPH_DIM=2
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_regression PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_regression PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_regression PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_regression PRIVATE -fopenmp)
    target_link_options(test_regression PRIVATE -fopenmp)
endif()

add_test(NAME regression_selfcheck COMMAND test_regression --selfcheck)
set_tests_properties(regression_selfcheck PROPERTIES LABELS "regression" TIMEOUT 60)
# По тесту на кейс и тип контакта — по той же причине, что у e2e: провал
# одного не скрывает остальные, и виновник называется по имени.
foreach(_case bc2 bc6)
    foreach(_type RANGE 0 3)
        add_test(NAME regression_${_case}_contact${_type}
                 COMMAND test_regression --case=${_case} --contact-type=${_type}
                         --time-limit=240)
        set_tests_properties(regression_${_case}_contact${_type}
                             PROPERTIES LABELS "regression" TIMEOUT 300)
    endforeach()
endforeach()

# Та же программа в 3D: кейс bc2mini3d (уменьшенный BC2, 3648 частиц), K = 50.
# Нужна, чтобы перенос Solid на новую схему RHS (fiber/PLAN.md, Ф4) не сломал
# ветку dim = 3, в которой будут жить нити.
add_executable(test_regression_3d tests/test_regression.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_regression_3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_regression_3d PRIVATE
    SPH_DIM=3
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_regression_3d PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_regression_3d PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_regression_3d PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_regression_3d PRIVATE -fopenmp)
    target_link_options(test_regression_3d PRIVATE -fopenmp)
endif()
add_test(NAME regression_selfcheck_3d COMMAND test_regression_3d --selfcheck)
set_tests_properties(regression_selfcheck_3d PROPERTIES LABELS "regression;3d" TIMEOUT 60)
foreach(_type RANGE 0 3)
    add_test(NAME regression_bc2mini3d_contact${_type}
             COMMAND test_regression_3d --case=bc2mini3d --contact-type=${_type}
                     --time-limit=240)
    set_tests_properties(regression_bc2mini3d_contact${_type}
                         PROPERTIES LABELS "regression;3d" TIMEOUT 300)
endforeach()

# Перегенерация эталонов — осознанное действие владельца, а не часть прогона
# тестов (запись в дерево исходников):
#     cmake --build build --target regression_golden
# Скрипт ставит OMP_NUM_THREADS=1, обходит обе программы (2D и 3D) и вписывает
# в шапку commit и дату.
add_custom_target(regression_golden
    COMMAND ${CMAKE_COMMAND}
            -DSPH_REG_EXE=$<TARGET_FILE:test_regression>
            -DSPH_REG_EXE_3D=$<TARGET_FILE:test_regression_3d>
            "-DSPH_REG_OUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/regression"
            "-DSPH_REG_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/regression_golden.cmake"
    DEPENDS test_regression test_regression_3d
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Перегенерация регрессионных эталонов в tests/data/regression"
    USES_TERMINAL VERBATIM)

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

# --- Профилирование: включает блоки #ifdef SPH_PROFILE (тайминг фаз, запись
#     в force_timing.csv) ----
if(SPH_PROFILE)
    target_compile_definitions(SPH PRIVATE SPH_PROFILE)
endif()

# --- Лог механической энергии: включает блоки #ifdef SPH_ENERGY_LOG и ключи
#     конфига energy_log_every/energy_log_file. Без опции кода лога в сборке нет.
if(SPH_ENERGY_LOG)
    target_compile_definitions(SPH PRIVATE SPH_ENERGY_LOG)
    target_compile_definitions(guo2019_shot PRIVATE SPH_ENERGY_LOG)
    target_compile_definitions(fabric_checks PRIVATE SPH_ENERGY_LOG)
endif()

# --- Диагностика контактов нитей: блоки #ifdef SPH_CONTACT_STATS и ключи
#     guo2019_shot --contact-stats/--hold-n. Без опции кода статистики нет.
if(SPH_CONTACT_STATS)
    target_compile_definitions(guo2019_shot PRIVATE SPH_CONTACT_STATS)
endif()

# --- Каталоги include: заголовки — в include/, исходники — в src/ ------------
target_include_directories(SPH PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# --- Демонстрационный расчёт: снаряд из SPH-частиц пробивает крестообразную
#     сеть из двух нитей (fiber/DEMO_VISUALIZATION_TASK.md). НЕ тест: строит
#     сцену и прогоняет её через продовый time_integration(cfg, fiber::Model&,
#     ...), чтобы получить кадры sph_output_*.txt/fiber_output_*.txt для
#     внешней Python-визуализации (experiments/2026-09-26-fiber-net-demo/).
#     add_test() сознательно не добавляется — это демо-артефакт, не
#     регрессионный тест (постановка, раздел 8 «Вне скоупа»).
add_executable(demo_fiber_net experiments/2026-09-26-fiber-net-demo/demo_fiber_net.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(demo_fiber_net PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(demo_fiber_net PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(demo_fiber_net PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(demo_fiber_net PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(demo_fiber_net PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(demo_fiber_net PRIVATE -fopenmp)
    target_link_options(demo_fiber_net PRIVATE -fopenmp)
endif()

# --- Демонстрационный расчёт: менее массивный снаряд пробивает/не пробивает
#     сеть из 4 слоёв нитей (x,y,x,y), experiments/2026-09-27-fiber-layers-demo/.
#     НЕ тест: та же архитектура, что demo_fiber_net (продовый
#     time_integration(cfg, fiber::Model&, ...)), два сценария (pierce/stop)
#     через argv, add_test() сознательно не добавляется.
add_executable(demo_fiber_layers experiments/2026-09-27-fiber-layers-demo/demo_fiber_layers.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(demo_fiber_layers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(demo_fiber_layers PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(demo_fiber_layers PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(demo_fiber_layers PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(demo_fiber_layers PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(demo_fiber_layers PRIVATE -fopenmp)
    target_link_options(demo_fiber_layers PRIVATE -fopenmp)
endif()

# --- Повреждённость связи с разупрочнением (DAMAGE; fiber/plan-fiber-damage.md):
#     закон на одной связи, GF, разрыв по D >= 1, разбор .fib, детерминизм.
add_executable(test_fiber_damage tests/test_fiber_damage.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_fiber_damage PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_fiber_damage PRIVATE SPH_DIM=3)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_damage PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_fiber_damage PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_fiber_damage PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_fiber_damage PRIVATE -fopenmp)
    target_link_options(test_fiber_damage PRIVATE -fopenmp)
endif()
add_test(NAME fiber_damage_test COMMAND test_fiber_damage)
set_tests_properties(fiber_damage_test PROPERTIES LABELS "3d")
