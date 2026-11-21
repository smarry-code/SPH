# bc2_plot.cmake
# Обёртка цели bc2_plot: четыре полных e2e-прогона (по одному на
# contact_force_type) с выгрузкой кривых, затем сборка сводного SVG из выгрузок.
# Судьба цели решается по ФАКТУ записи графика, а не по коду возврата тестов.
#
# Зачем так. Назначение цели — перестроить картинку для отчёта. Прогоны при этом
# могут закончиться неудачей по существу расчёта, но график в этом случае нужен
# тем более — именно по нему расхождение и разбирают. Поэтому ненулевые коды не
# гасят цель, а печатаются предупреждением. Сейчас все четыре прогона проходят;
# содержательное расхождение с приборным эталоном (занижённый пик) осталось, но
# грубые ворота сверки по форме оно проходит.
#
# Почему сборка отдельным шагом. Один прогон отвечает за один тип контакта
# (так устроены и тесты), значит четыре кривые живут в четырёх процессах.
# Выгрузка в файл — единственный способ свести их на один график, не заводя
# пятый прогон ради картинки.
#
# Ожидаемые переменные (-D): SPH_PLOT_EXE, SPH_PLOT_OUT, SPH_PLOT_TIME_LIMIT.

if(NOT SPH_PLOT_EXE OR NOT SPH_PLOT_OUT)
    message(FATAL_ERROR "bc2_plot.cmake: нужны -DSPH_PLOT_EXE и -DSPH_PLOT_OUT")
endif()

# Отметка времени «до»: по ней отличается свежая запись от старого файла,
# оставшегося с прошлого раза (иначе упавший на старте прогон выглядел бы успехом).
set(_stamp_before "")
if(EXISTS "${SPH_PLOT_OUT}")
    file(TIMESTAMP "${SPH_PLOT_OUT}" _stamp_before "%Y-%m-%d %H:%M:%S" UTC)
endif()

# Каталог обмена кривыми. Лежит в каталоге сборки, а не в дереве исходников:
# это промежуточные данные цели, а не результат.
set(_curves "${CMAKE_CURRENT_BINARY_DIR}/bc2_curves")
file(REMOVE_RECURSE "${_curves}")
file(MAKE_DIRECTORY "${_curves}")

set(_bad_types "")
foreach(_type RANGE 0 3)
    message(STATUS "bc2_plot: прогон contact_force_type = ${_type}")
    execute_process(
        COMMAND "${SPH_PLOT_EXE}" --full "--contact-type=${_type}"
                "--time-limit=${SPH_PLOT_TIME_LIMIT}"
                "--emit-curve=${_curves}/bc2_curve_${_type}.txt"
        RESULT_VARIABLE _code)
    if(NOT _code EQUAL 0)
        list(APPEND _bad_types "${_type}(код ${_code})")
    endif()
endforeach()

message(STATUS "bc2_plot: сборка сводного графика")
execute_process(
    COMMAND "${SPH_PLOT_EXE}" --full "--plot-from=${_curves}" "--emit-plot=${SPH_PLOT_OUT}"
    RESULT_VARIABLE _plot_code)

if(NOT EXISTS "${SPH_PLOT_OUT}")
    message(FATAL_ERROR "bc2_plot: график не создан (${SPH_PLOT_OUT}), код сборки ${_plot_code}")
endif()

file(TIMESTAMP "${SPH_PLOT_OUT}" _stamp_after "%Y-%m-%d %H:%M:%S" UTC)
file(SIZE "${SPH_PLOT_OUT}" _size)

if(_stamp_before STREQUAL _stamp_after AND NOT _plot_code EQUAL 0)
    message(FATAL_ERROR
        "bc2_plot: сборка графика завершилась с кодом ${_plot_code}, файл не обновлён "
        "(остался от прошлого запуска): ${SPH_PLOT_OUT}")
endif()

if(_bad_types)
    string(REPLACE ";" ", " _bad_list "${_bad_types}")
    message(WARNING
        "bc2_plot: график перестроен (${_size} байт), НО прогоны типов ${_bad_list} "
        "завершились неудачей — см. вывод выше. Для отчёта график годен; "
        "провалы разбираются отдельно.")
else()
    message(STATUS "bc2_plot: график перестроен (${_size} байт), все прогоны пройдены")
endif()
