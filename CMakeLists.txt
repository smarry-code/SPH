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
# Модули конвейера БЕЗ точки входа: переиспользуются целями тестов, которые
# приносят собственный main() (ТЗ contact/04-e2e-test-brief.md §7 — общий
# список исходников в переменной, без второго file(GLOB)).
set(SPH_PIPELINE_SOURCES ${SPH_SOURCES})
list(REMOVE_ITEM SPH_PIPELINE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
# --- Исполняемый модуль (имя SPH.exe на Windows — как прежде) ----------------
add_executable(SPH ${SPH_SOURCES})
set_target_properties(SPH PROPERTIES OUTPUT_NAME "SPH")

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

# --- Регрессия fp32 (fiber/CUDA_PLAN.md, этап Г0-FP32) -----------------------
# Та же программа и те же кейсы, но конвейер собран с SPH_PRECISION=32
# (include/precision.hpp: real_t = float). Свои эталоны —
# tests/data/regression_fp32/ (снимаются этой же программой, цель
# regression_golden_fp32), свои допуски (tests/test_regression.cpp, блок
# SPH_PRECISION == 32). fp64-эталоны и fp64-цели выше не затрагиваются.
#
# -Wdouble-promotion / -Wimplicit-float-conversion: любое неявное повышение
# float -> double в конвейере означает, что часть «fp32»-расчёта на деле идёт
# в double (литерал 2.0 вместо real_t(2.0) и т.п.). Предупреждения, а не
# ошибки: тестовая обвязка (геометрия, снимок) законно работает в double.
option(SPH_FP32_WARNINGS "Предупреждения о повышении float -> double в fp32-целях" ON)
# src/input.cpp — разбор конфигурации и встроенные генераторы начальных
# условий: значения считаются в double и записываются в поля real_t. Сужение
# там и есть граница ввода (include/precision.hpp), а не утечка double в расчёт,
# поэтому для этого файла оба предупреждения выключены.
if(SPH_FP32_WARNINGS AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set_source_files_properties(src/input.cpp PROPERTIES
        COMPILE_OPTIONS "-Wno-double-promotion;-Wno-implicit-float-conversion")
endif()
foreach(_dim 2 3)
    if(_dim EQUAL 2)
        set(_t test_regression_fp32)
    else()
        set(_t test_regression_3d_fp32)
    endif()
    add_executable(${_t} tests/test_regression.cpp ${SPH_PIPELINE_SOURCES})
    target_include_directories(${_t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_compile_definitions(${_t} PRIVATE
        SPH_DIM=${_dim}
        SPH_PRECISION=32
        SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(${_t} PRIVATE -fexceptions)
        if(SPH_FP32_WARNINGS)
            target_compile_options(${_t} PRIVATE -Wdouble-promotion)
            if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                target_compile_options(${_t} PRIVATE -Wimplicit-float-conversion)
            endif()
        endif()
    endif()
    if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        target_compile_options(${_t} PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
        target_link_options(${_t} PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
    elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(${_t} PRIVATE -fopenmp)
        target_link_options(${_t} PRIVATE -fopenmp)
    endif()
endforeach()

add_test(NAME regression_fp32_selfcheck COMMAND test_regression_fp32 --selfcheck)
set_tests_properties(regression_fp32_selfcheck PROPERTIES LABELS "regression;fp32" TIMEOUT 60)
foreach(_case bc2 bc6)
    foreach(_type RANGE 0 3)
        add_test(NAME regression_fp32_${_case}_contact${_type}
                 COMMAND test_regression_fp32 --case=${_case} --contact-type=${_type}
                         --time-limit=240)
        set_tests_properties(regression_fp32_${_case}_contact${_type}
                             PROPERTIES LABELS "regression;fp32" TIMEOUT 300)
    endforeach()
endforeach()
add_test(NAME regression_fp32_selfcheck_3d COMMAND test_regression_3d_fp32 --selfcheck)
set_tests_properties(regression_fp32_selfcheck_3d PROPERTIES LABELS "regression;fp32;3d" TIMEOUT 60)
foreach(_type RANGE 0 3)
    add_test(NAME regression_fp32_bc2mini3d_contact${_type}
             COMMAND test_regression_3d_fp32 --case=bc2mini3d --contact-type=${_type}
                     --time-limit=240)
    set_tests_properties(regression_fp32_bc2mini3d_contact${_type}
                         PROPERTIES LABELS "regression;fp32;3d" TIMEOUT 300)
endforeach()

# Перегенерация fp32-эталонов — тот же скрипт, что у regression_golden
# (OMP_NUM_THREADS=1, commit и дата в шапке), другие программы и каталог.
#     cmake --build build --target regression_golden_fp32
add_custom_target(regression_golden_fp32
    COMMAND ${CMAKE_COMMAND}
            -DSPH_REG_EXE=$<TARGET_FILE:test_regression_fp32>
            -DSPH_REG_EXE_3D=$<TARGET_FILE:test_regression_3d_fp32>
            "-DSPH_REG_OUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/regression_fp32"
            "-DSPH_REG_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/regression_golden.cmake"
    DEPENDS test_regression_fp32 test_regression_3d_fp32
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Перегенерация fp32-эталонов в tests/data/regression_fp32"
    USES_TERMINAL VERBATIM)

# Сквозной BC2 в fp32: тот же tests/test_e2e_bc2.cpp (приборный эталон, грубые
# ворота формы кривой), конвейер с SPH_PRECISION=32. Отвечает на вопрос
# «остаётся ли fp32-расчёт физикой», которого регрессия (сверка fp32 с самим
# собой) не задаёт. Smoke — в обычном прогоне, полный — под SPH_E2E_SLOW_TESTS.
add_executable(test_e2e_bc2_fp32 tests/test_e2e_bc2.cpp ${SPH_PIPELINE_SOURCES})
target_include_directories(test_e2e_bc2_fp32 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_e2e_bc2_fp32 PRIVATE
    SPH_DIM=2
    SPH_PRECISION=32
    SPH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_bc2_fp32 PRIVATE -fexceptions)
endif()
if(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_compile_options(test_e2e_bc2_fp32 PRIVATE ${SPH_APPLECLANG_OMP_COMPILE_OPTIONS})
    target_link_options(test_e2e_bc2_fp32 PRIVATE ${SPH_APPLECLANG_OMP_LINK_OPTIONS})
elseif(SPH_USE_OPENMP AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_e2e_bc2_fp32 PRIVATE -fopenmp)
    target_link_options(test_e2e_bc2_fp32 PRIVATE -fopenmp)
endif()
foreach(_type RANGE 0 3)
    add_test(NAME e2e_fp32_bc2_smoke_contact${_type}
             COMMAND test_e2e_bc2_fp32 --contact-type=${_type} --time-limit=120)
    set_tests_properties(e2e_fp32_bc2_smoke_contact${_type} PROPERTIES LABELS "e2e;fp32" TIMEOUT 300)
endforeach()
if(SPH_E2E_SLOW_TESTS)
    foreach(_type RANGE 0 3)
        add_test(NAME e2e_fp32_bc2_full_contact${_type}
                 COMMAND test_e2e_bc2_fp32 --full --contact-type=${_type}
                         --time-limit=${SPH_BC2_RUN_TIME_LIMIT})
        set_tests_properties(e2e_fp32_bc2_full_contact${_type}
                             PROPERTIES LABELS "slow;e2e;fp32" TIMEOUT 1800)
    endforeach()
endif()

# --- Точность расчёта основной программы (include/precision.hpp): 64 — double
#     (по умолчанию, эталонный расчёт), 32 — float.
set(SPH_PRECISION "64" CACHE STRING "Точность расчёта цели SPH: 64 (double) | 32 (float)")
set_property(CACHE SPH_PRECISION PROPERTY STRINGS 64 32)
target_compile_definitions(SPH PRIVATE SPH_PRECISION=${SPH_PRECISION})

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

# --- Каталоги include: заголовки — в include/, исходники — в src/ ------------
target_include_directories(SPH PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
