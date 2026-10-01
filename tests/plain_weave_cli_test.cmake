# CLI-тест plain_weave_gen (cmake -P): -DGEN=<путь к plain_weave_gen>.
# Выходные файлы -- в уникальном каталоге системного temp, удаляется в конце.
if(NOT GEN)
    message(FATAL_ERROR "не задан -DGEN=<plain_weave_gen>")
endif()
if(DEFINED ENV{TMPDIR} AND NOT "$ENV{TMPDIR}" STREQUAL "")
    set(_tmp "$ENV{TMPDIR}")
else()
    set(_tmp "/tmp")
endif()
string(RANDOM LENGTH 12 ALPHABET "abcdefghijklmnopqrstuvwxyz0123456789" _rnd)
string(TIMESTAMP _ts "%Y%m%d%H%M%S")
set(_dir "${_tmp}/plain_weave_cli_test_${_ts}_${_rnd}")
file(MAKE_DIRECTORY "${_dir}")

set(_fail "")
function(_run expected_rc)
    execute_process(COMMAND ${GEN} ${ARGN} RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
    if(NOT rc EQUAL expected_rc)
        set(_fail "${_fail}\nвызов '${ARGN}': код ${rc}, ожидался ${expected_rc}" PARENT_SCOPE)
    endif()
endfunction()

_run(0 --preset=test --out=${_dir}/test.fib)
foreach(f test.fib test.fib.layout.txt)
    if(NOT EXISTS "${_dir}/${f}")
        set(_fail "${_fail}\nне создан файл ${f}")
    endif()
endforeach()
# Нулевой выступ допустим.
_run(0 --preset=test --overhang-x=0 --out=${_dir}/oh0.fib)
# p < d без --allow-penetration -- ошибка построения (1).
_run(1 --preset=test --pitch=1e-4 --out=${_dir}/bad.fib)
# Неизвестный ключ и нечисловые значения -- ошибка использования (2).
_run(2 --preset=test --no-such-key=1 --out=${_dir}/bad2.fib)
_run(2 --preset=test --amp=inf --out=${_dir}/bad3.fib)
_run(2 --preset=test --rho=1e999 --out=${_dir}/bad4.fib)
# Невалидное значение материала -- ошибка построения (1).
_run(1 --preset=test --mu-ff=-1 --out=${_dir}/bad5.fib)

file(REMOVE_RECURSE "${_dir}")
if(NOT "${_fail}" STREQUAL "")
    message(FATAL_ERROR "plain_weave_gen CLI:${_fail}")
endif()
