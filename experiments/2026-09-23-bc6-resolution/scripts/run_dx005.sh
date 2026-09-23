#!/usr/bin/env bash
# Шестая ступень лестницы: Δx = 0,05 мм, 706 406 частиц.
# ПОДГОТОВЛЕН, НЕ ЗАПУСКАЛСЯ: на машине, где считались остальные ступени
# (15,2 ГБ ОЗУ, 4,6 ГБ свободно), памяти впритык. Запуск — на более мощной.
#
# Оценки (не измерения), масштабированием от прогона 0,10 мм:
#   память  ≈ 5 ГБ: буфер пар 140·N·40 Б ≈ 4,0 ГБ выделяется заново на каждом
#           шаге (tests/test_e2e_bc6.cpp:186, src/single_step.cpp:114), при
#           0,10 мм рабочий набор был 1,2 ГБ; нужно ≥ 8 ГБ свободной памяти;
#   шагов   ≈ 6500 (вдвое больше, чем 3232 при 0,10 мм);
#   время   ≈ ×8 от 0,10 мм (4726 с на 6 потоках) — порядка 10 ч на 6 потоках.
#
# Запуск из любого каталога:
#   bash experiments/2026-09-23-bc6-resolution/scripts/run_dx005.sh
# Потоки — через OMP_NUM_THREADS (по умолчанию все логические ядра).
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -z "${OMP_NUM_THREADS:-}" ]; then
    OMP_NUM_THREADS=$(nproc 2>/dev/null || echo 6)
fi
export OMP_NUM_THREADS
# Сторож: сутки. Штатные 6 ч из one.sh этому прогону заведомо малы.
export SPH_TIME_LIMIT="${SPH_TIME_LIMIT:-86400}"
if [ -r /proc/meminfo ]; then
    avail_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
    if [ "${avail_kb:-0}" -lt 8000000 ]; then
        echo "ВНИМАНИЕ: свободно $((avail_kb/1024)) МБ, оценка потребности ≈ 5 ГБ, рекомендуется ≥ 8 ГБ" >&2
    fi
fi
echo "Δx = 0,05 мм, 706 406 частиц; OMP_NUM_THREADS=$OMP_NUM_THREADS; бюджет ${SPH_TIME_LIMIT} с"
bash "$HERE/one.sh" 0.5e-4 dx005
