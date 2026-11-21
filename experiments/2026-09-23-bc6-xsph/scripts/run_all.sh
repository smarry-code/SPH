#!/usr/bin/env bash
# Все прогоны серии, КРОМЕ Δx = 0,05 мм: четыре строки data/todo.txt, по два
# одновременно, самый долгий (0,10 мм) первым. Около 80 мин по стенным часам
# на 12-ядерной машине (оценка по стоимости ступеней 2026-09-23-bc6-resolution).
# Требует наложенной заплаты 01 и пересобранного test_e2e_bc6 — см. REPRO.md.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXP="$(dirname "$HERE")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-6}"
xargs -P 2 -a "$EXP/data/todo.txt" -n 3 bash "$HERE/one.sh"
