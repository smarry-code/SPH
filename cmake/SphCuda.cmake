# ============================================================================
# CUDA-бэкенд (fiber/CUDA_PLAN.md, этап Г0): опция SPH_ENABLE_CUDA,
# подключение nvcc, архитектуры устройства и флаги компиляции .cu.
# Подключается из CMakeLists.txt после project(... LANGUAGES CXX) — язык
# CUDA включается здесь, только когда опция ON, и не трогает CPU-цели.
#
# Решение Р-Г1 (fiber/CUDA_PLAN.md §2, спайк 2026-09-27): .cu компилирует
# ТОЛЬКО nvcc, весь остальной проект — Clang, в одной сборке. На Windows у
# nvcc единственный поддерживаемый хост-компилятор — MSVC cl (флаг
# -ccbin clang там не поддерживается); Clang как КОМПИЛЯТОР CUDA технически
# работает (спайк это подтвердил), но в проекте не используется — решение
# владельца зафиксировано в ADR-009.
# ============================================================================

option(SPH_ENABLE_CUDA "Собирать CUDA-бэкенд (nvcc) в дополнение к CPU-конвейеру" OFF)

if(NOT SPH_ENABLE_CUDA)
    return()
endif()

# --- Архитектуры устройства (решение Р-Г1а). По умолчанию — машинный код
#     для карты разработки (RTX 4050 Laptop, sm_89) и PTX архитектуры Turing
#     (75) — младшей, которую поддерживает CUDA 13: PTX X запускается на
#     любой карте с compute capability >= X, поэтому один и тот же
#     build-cuda запустится и на другой карте (вторая машина владельца —
#     Linux с NVIDIA неизвестной модели). Пресет cuda-linux переопределяет
#     это на "native": CMake >= 3.24 сам определяет карту сборочной машины.
set(SPH_CUDA_ARCHITECTURES "89-real;75-virtual" CACHE STRING
    "Архитектуры CUDA (CMAKE_CUDA_ARCHITECTURES); native — определить по карте сборочной машины")

# --- Заготовка строгого режима FP (план §5 «Как сравнивать CPU и GPU»):
#     отключает слияние умножения со сложением (FMA) у nvcc, чтобы при
#     разборе расхождений CPU/GPU свести их к ULP библиотечных функций.
#     Используется начиная с Г4 (SPH_HD-функции ещё не существуют в Г0).
option(SPH_CUDA_STRICT_FP "Собирать CUDA-код с -fmad=false (разбор расхождений CPU/GPU, §5 плана)" OFF)

# --- Windows: nvcc компилирует хост-часть .cu ТОЛЬКО через MSVC cl и ищет
#     его исключительно через PATH (без vcvars nvcc падает с "Cannot find
#     compiler 'cl.exe' in PATH", а enable_language(CUDA) — невнятной
#     ошибкой CMake). Проверяем сами и даём понятную инструкцию (Г0: скрипт
#     tools/cuda_env.bat и пресет cuda).
if(WIN32)
    find_program(SPH_CUDA_HOST_CL_EXE NAMES cl)
    if(NOT SPH_CUDA_HOST_CL_EXE)
        message(FATAL_ERROR
            "SPH_ENABLE_CUDA=ON: компилятор MSVC cl.exe не найден в PATH. "
            "nvcc на Windows компилирует хост-часть .cu только через cl и "
            "ищет его через PATH (fiber/CUDA_PLAN.md, решение Р-Г1). "
            "Запустите конфигурацию из окружения MSVC:\n"
            "    tools\\cuda_env.bat cmake --preset cuda")
    endif()
endif()

# --- Компилятор CUDA не должен оказаться Clang (решение Р-Г1: единственный
#     используемый компилятор .cu в проекте — nvcc). Проверяем ДО
#     enable_language(CUDA): если пользователь явно передал
#     -DCMAKE_CUDA_COMPILER=clang++, ошибка должна быть понятной, а не
#     свалиться в глубине CUDA-специфичной конфигурации CMake.
if(DEFINED CMAKE_CUDA_COMPILER AND CMAKE_CUDA_COMPILER MATCHES "clang")
    message(FATAL_ERROR
        "SPH_ENABLE_CUDA=ON: CMAKE_CUDA_COMPILER указывает на Clang "
        "(${CMAKE_CUDA_COMPILER}). Решение Р-Г1 (fiber/CUDA_PLAN.md, ADR-009): "
        "компилятор .cu в проекте — nvcc; Clang как компилятор CUDA не "
        "используется. Уберите -DCMAKE_CUDA_COMPILER или укажите nvcc.")
endif()

enable_language(CUDA)

# Повторная проверка ПОСЛЕ enable_language: если CMake сам выбрал Clang
# (например, nvcc не нашёлся, а clang++ в PATH есть) — тоже понятная ошибка,
# хотя откатить enable_language(CUDA) уже нельзя, конфигурация всё равно
# останавливается здесь же.
if(CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR
        "SPH_ENABLE_CUDA=ON: CMake выбрал Clang компилятором CUDA "
        "(${CMAKE_CUDA_COMPILER}). Решение Р-Г1: нужен nvcc. Укажите его "
        "явно: -DCMAKE_CUDA_COMPILER=nvcc (пресет cuda уже это делает).")
endif()

# --- CUDAToolkit: пути CCCL (CUB/Thrust, в CUDA 13 — include/cccl)
#     добавляются им сам. Минимальная версия — 12.4 (Р-Г1а): CUDA 13 не
#     поддерживает карты старше Turing (sm_75), а вторая машина владельца
#     (Linux) может оказаться со старой картой. Сборка на 12.x проверяется
#     на машине разработки, где CUDA 12.6 установлена рядом (Г0, отчёт —
#     в ADR-009).
find_package(CUDAToolkit 12.4 REQUIRED)
if(CUDAToolkit_VERSION VERSION_LESS 12.4)
    message(FATAL_ERROR
        "SPH_ENABLE_CUDA=ON: найден CUDA Toolkit ${CUDAToolkit_VERSION}, "
        "нужен >= 12.4 (fiber/CUDA_PLAN.md, решение Р-Г1а).")
endif()

set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES "${SPH_CUDA_ARCHITECTURES}")

# --- Хост-компилятор nvcc на Linux: по умолчанию тот же Clang, что и для
#     остального проекта (Р-Г1); если конкретная связка nvcc+Clang на
#     машине не собирается, откат на GCC — передаётся явно
#     (-DCMAKE_CUDA_HOST_COMPILER=g++). На Windows CMAKE_CUDA_HOST_COMPILER
#     НЕ задаётся: там у nvcc единственный хост — MSVC cl, найденный через
#     PATH (проверка выше); явное указание сбивает автоопределение
#     тулсета MSVC самим nvcc.
if(NOT WIN32 AND NOT CMAKE_CUDA_HOST_COMPILER)
    set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE FILEPATH
        "Хост-компилятор nvcc для .cu (по умолчанию — тот же Clang, что и CXX)")
endif()

# --- sph_cuda_configure_target(target): общие флаги CUDA-цели. Вызывается
#     для каждой CUDA-цели (первая — test_cuda_toolchain и cuda_microbench,
#     Г0; далее — src/cuda/*.cu, Г3+).
#   --extended-lambda            device-лямбды (их используют CUB/Thrust и
#                                 будущий include/physics/*, Р-Г2).
#   -DNOMINMAX                   windows.h определяет макросы min/max,
#                                 которые ломают CUB (план §1.1/§3.2);
#                                 SPH_HD-код использует sph_min/sph_max.
#   -Xcompiler=/Zc:preprocessor  без него CCCL останавливается с
#                                 "#error: traditional preprocessor" под cl.
#   -Xcompiler=/Zc:__cplusplus   иначе __cplusplus == 199711L под cl.
#   -Xcompiler=/utf-8,/EHsc,/W4  как у остального проекта под MSVC-хостом.
function(sph_cuda_configure_target target)
    set_target_properties(${target} PROPERTIES
        CUDA_STANDARD 20
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION OFF)
    target_compile_options(${target} PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
        $<$<COMPILE_LANGUAGE:CUDA>:-DNOMINMAX>)
    if(WIN32)
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor,/Zc:__cplusplus,/utf-8,/EHsc,/W4>)
    endif()
    if(SPH_CUDA_STRICT_FP)
        target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>)
    endif()
    target_link_libraries(${target} PRIVATE CUDA::cudart_static)
endfunction()

message(STATUS
    "SPH: CUDA backend ON — компилятор ${CMAKE_CUDA_COMPILER} "
    "(CUDAToolkit ${CUDAToolkit_VERSION}), архитектуры ${CMAKE_CUDA_ARCHITECTURES}")
