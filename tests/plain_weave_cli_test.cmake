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
# Одеяло из патчей: файлы .fib, .layout.txt и .patches.txt.
_run(0 --preset=twaron_ct717 --patch=34x34 --tile=3x3 --out=${_dir}/quilt.fib)
foreach(f quilt.fib quilt.fib.layout.txt quilt.fib.patches.txt)
    if(NOT EXISTS "${_dir}/${f}")
        set(_fail "${_fail}\nне создан файл ${f}")
    endif()
endforeach()
# Нечётный патч при tile > 1: ошибка построения (1); с --allow-odd-patch -- успех.
_run(1 --preset=test --patch=3x3 --tile=2x2 --out=${_dir}/odd.fib)
_run(0 --preset=test --patch=3x3 --tile=2x2 --allow-odd-patch --out=${_dir}/odd_ok.fib)
# Конфликт с --nx/--ny, неполная пара и разбор -- ошибка использования (2).
_run(2 --preset=test --patch=4x4 --tile=2x2 --nx=8 --out=${_dir}/c1.fib)
_run(2 --preset=test --patch=4x4 --out=${_dir}/c2.fib)
_run(2 --preset=test --patch=4 --tile=2x2 --out=${_dir}/c3.fib)
_run(2 --preset=test --allow-odd-patch --out=${_dir}/c4.fib)

# Случайная раскройка: файлы .fib, .layout.txt и .cut.json.
_run(0 --preset=twaron_ct717 --random-cut=9 --size=0.05x0.05 --seed=3 --out=${_dir}/rc.fib)
foreach(f rc.fib rc.fib.layout.txt rc.fib.cut.json)
    if(NOT EXISTS "${_dir}/${f}")
        set(_fail "${_fail}\nне создан файл ${f}")
    endif()
endforeach()
# Конфликты --random-cut с --nx/--patch/--tile -- ошибка использования (2).
_run(2 --preset=twaron_ct717 --random-cut=9 --nx=10 --out=${_dir}/rc_e1.fib)
_run(2 --preset=twaron_ct717 --random-cut=4 --patch=4x4 --tile=2x2 --out=${_dir}/rc_e2.fib)
# Многослойный пакет: две раскройки по два слоя.
_run(0 --preset=twaron_ct717 --random-cut=4 --size=0.02x0.02 --layers=4 --layers-per-cut=2 --out=${_dir}/rc4.fib)
_run(2 --preset=twaron_ct717 --random-cut=4 --layers=4 --layers-per-cut=5 --out=${_dir}/rc_e3.fib)
# Ключи раскройки без --random-cut и некорректный --size -- ошибка использования (2).
_run(2 --preset=test --layers=2 --out=${_dir}/rc_e4.fib)
_run(2 --preset=twaron_ct717 --random-cut=4 --size=abc --out=${_dir}/rc_e5.fib)
# gap < d без --allow-penetration -- ошибка построения (1).
_run(1 --preset=twaron_ct717 --random-cut=4 --size=0.02x0.02 --gap=1e-4 --out=${_dir}/rc_e6.fib)

file(REMOVE_RECURSE "${_dir}")
if(NOT "${_fail}" STREQUAL "")
    message(FATAL_ERROR "plain_weave_gen CLI:${_fail}")
endif()
