# bc2_plot.cmake
# Обёртка цели bc2_plot: запускает полный e2e-прогон с --emit-plot и решает
# судьбу цели по ФАКТУ записи графика, а не по коду возврата теста.
#
# Зачем. Назначение цели — перестроить картинку для отчёта. Тест при этом может
# закончиться неудачей по существу расчёта (сейчас так и есть: сверка с
# приборным эталоном не проходит по регрессионным допускам), но график в этом
# случае нужен тем более — именно по нему расхождение и разбирают. Поэтому
# ненулевой код теста не гасит цель, а печатается как предупреждение.
#
# Ожидаемые переменные (-D): SPH_PLOT_EXE, SPH_PLOT_OUT, SPH_PLOT_TIME_LIMIT.

if(NOT SPH_PLOT_EXE OR NOT SPH_PLOT_OUT)
    message(FATAL_ERROR "bc2_plot.cmake: нужны -DSPH_PLOT_EXE и -DSPH_PLOT_OUT")
endif()

# Отметка времени «до»: по ней отличается свежая запись от старого файла,
# оставшегося с прошлого раза (иначе упавший на старте тест выглядел бы успехом).
set(_stamp_before "")
if(EXISTS "${SPH_PLOT_OUT}")
    file(TIMESTAMP "${SPH_PLOT_OUT}" _stamp_before "%Y-%m-%d %H:%M:%S" UTC)
endif()

execute_process(
    COMMAND "${SPH_PLOT_EXE}" --full "--time-limit=${SPH_PLOT_TIME_LIMIT}" "--emit-plot=${SPH_PLOT_OUT}"
    RESULT_VARIABLE _code)

if(NOT EXISTS "${SPH_PLOT_OUT}")
    message(FATAL_ERROR "bc2_plot: график не создан (${SPH_PLOT_OUT}), код теста ${_code}")
endif()

file(TIMESTAMP "${SPH_PLOT_OUT}" _stamp_after "%Y-%m-%d %H:%M:%S" UTC)
file(SIZE "${SPH_PLOT_OUT}" _size)

if(_stamp_before STREQUAL _stamp_after AND NOT _code EQUAL 0)
    message(FATAL_ERROR
        "bc2_plot: тест завершился с кодом ${_code}, файл графика не обновлён "
        "(остался от прошлого запуска): ${SPH_PLOT_OUT}")
endif()

if(NOT _code EQUAL 0)
    message(WARNING
        "bc2_plot: график перестроен (${_size} байт), НО тест завершился с кодом ${_code} — "
        "см. вывод выше. Для отчёта график годен; провал теста разбирается отдельно.")
else()
    message(STATUS "bc2_plot: график перестроен (${_size} байт), тест пройден")
endif()
