#!/usr/bin/env bash
# ======================================================================
#  НЕ ЗАПУСКАТЬ на машине серии (15,2 ГБ ОЗУ). Решение владельца от
#  2026-09-23: ступень Δx = 0,05 мм считается позже, на более мощной
#  системе. В run_all.sh и data/todo.txt она намеренно не входит.
# ======================================================================
# ε = 0,05, Δx = 0,05 мм, 706 406 частиц. Оценки (не измерения):
#   память ≈ 5 ГБ (буфер пар 140·N·40 Б ≈ 4 ГБ заново на каждом шаге),
#          нужно ≥ 8 ГБ свободной; время порядка 10 ч на 6 потоках.
# Пара с ε = 0,25 при том же Δx — experiments/2026-09-23-bc6-resolution,
# scripts/run_dx005.sh; без неё этот прогон сравнить не с чем.
#
# Скрипт откажется работать без явного разрешения:
#   SPH_RUN_HEAVY=1 bash scripts/run_e005_dx005.sh
set -u
if [ "${SPH_RUN_HEAVY:-0}" != "1" ]; then
    echo "Отказ: ступень 0,05 мм не запускается на этой машине (см. шапку скрипта)." >&2
    echo "На мощной системе: SPH_RUN_HEAVY=1 bash $0" >&2
    exit 2
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -z "${OMP_NUM_THREADS:-}" ]; then OMP_NUM_THREADS=$(nproc 2>/dev/null || echo 6); fi
export OMP_NUM_THREADS
export SPH_TIME_LIMIT="${SPH_TIME_LIMIT:-86400}"
if [ -r /proc/meminfo ]; then
    avail_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
    if [ "${avail_kb:-0}" -lt 8000000 ]; then
        echo "ВНИМАНИЕ: свободно $((avail_kb/1024)) МБ, рекомендуется ≥ 8 ГБ" >&2
    fi
fi
bash "$HERE/one.sh" 0.05 0.5e-4 e005_dx005
