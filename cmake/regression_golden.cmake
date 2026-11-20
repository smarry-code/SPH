# regression_golden.cmake
# Обёртка цели regression_golden: перегенерация эталонов регрессионного теста
# (tests/test_regression.cpp) по всем кейсам — BC2 и BC6 (dim = 2) и
# уменьшенный BC2 в 3D (dim = 3), contact_force_type 0..3, контрольные точки
# N = 1 и N = K.
#
# Почему отдельная цель, а не ключ теста. Запись в дерево исходников — побочный
# эффект, которого у теста быть не должно. Перегенерация эталона — осознанное
# действие владельца: только при намеренном изменении физики, с причиной в
# сообщении коммита (docs/02-how-to/01-run-tests.md, раздел «Регрессия»).
#
# Почему OMP_NUM_THREADS=1. При нескольких потоках результат плавает на уровне
# округления (порядок пар в поиске соседей раздаётся через atomic). Эталон,
# снятый в один поток, воспроизводим: повторная генерация на той же сборке
# даёт те же файлы, и их дифф показывает только настоящие изменения.
#
# Ожидаемые переменные (-D): SPH_REG_EXE (2D), SPH_REG_EXE_3D, SPH_REG_CASES_DIR (каталог cases/),
# SPH_REG_SOURCE_DIR.

if(NOT SPH_REG_EXE OR NOT SPH_REG_EXE_3D OR NOT SPH_REG_CASES_DIR OR NOT SPH_REG_SOURCE_DIR)
    message(FATAL_ERROR
        "regression_golden.cmake: нужны -DSPH_REG_EXE, -DSPH_REG_EXE_3D, -DSPH_REG_CASES_DIR, "
        "-DSPH_REG_SOURCE_DIR")
endif()

set(ENV{OMP_NUM_THREADS} 1)

execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY "${SPH_REG_SOURCE_DIR}"
    OUTPUT_VARIABLE _commit OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _git_code ERROR_QUIET)
if(NOT _git_code EQUAL 0 OR NOT _commit)
    set(_commit "unknown")
endif()
execute_process(
    COMMAND git status --porcelain --untracked-files=no
    WORKING_DIRECTORY "${SPH_REG_SOURCE_DIR}"
    OUTPUT_VARIABLE _dirty OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(_dirty)
    # Эталон снят с незакоммиченного кода: commit в шапке его не описывает.
    set(_commit "${_commit}+dirty")
endif()
string(TIMESTAMP _date "%Y-%m-%d" UTC)

set(_bad "")
# Кейс -> программа его размерности.
set(_exe_bc2 "${SPH_REG_EXE}")
set(_exe_bc6 "${SPH_REG_EXE}")
set(_exe_bc2mini3d "${SPH_REG_EXE_3D}")
# Кейс -> каталог постановки (уменьшенный BC2 живёт в каталоге BC2).
set(_dir_bc2 "bc2")
set(_dir_bc6 "bc6")
set(_dir_bc2mini3d "bc2")
foreach(_case bc2 bc6 bc2mini3d)
    set(_out "${SPH_REG_CASES_DIR}/${_dir_${_case}}/regression/snapshots")
    file(MAKE_DIRECTORY "${_out}")
    foreach(_type RANGE 0 3)
        message(STATUS "regression_golden: ${_case}, contact_force_type = ${_type}")
        execute_process(
            COMMAND "${_exe_${_case}}" "--case=${_case}" "--contact-type=${_type}"
                    "--emit-golden=${_out}"
                    "--meta=commit=${_commit}" "--meta=date=${_date}"
                    "--meta=generated=cmake --build build --target regression_golden"
            RESULT_VARIABLE _code)
        if(NOT _code EQUAL 0)
            list(APPEND _bad "${_case}/contact${_type}(код ${_code})")
        endif()
    endforeach()
endforeach()

if(_bad)
    string(REPLACE ";" ", " _bad_list "${_bad}")
    message(FATAL_ERROR "regression_golden: генерация не удалась для ${_bad_list}")
endif()
message(STATUS "regression_golden: эталоны записаны в ${SPH_REG_CASES_DIR}/<кейс>/regression/snapshots (commit ${_commit}, ${_date})")
