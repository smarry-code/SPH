#!/usr/bin/env bash
# Один прогон: one.sh <contact_force_type> <delta_x, м>
# Кривая пишется в data/bc6_curve_<тип>.txt — под этим именем её ищет
# штатный режим --plot-from.
# Полный кейс BC6, штатная конфигурация; меняется только шаг решётки
# (штатный ключ --delta-x) и тип контакта (штатный --contact-type).
#
# Переменные окружения:
#   OMP_NUM_THREADS  потоки на прогон (по умолчанию 6)
#   SPH_TIME_LIMIT   бюджет сторожа, с (по умолчанию 21600 = 6 ч)
#   SPH_E2E_EXE      путь к бинарю (по умолчанию build/test_e2e_bc6[.exe])
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXP="$(dirname "$HERE")"
ROOT="$(cd "$EXP/../.." && pwd)"
if [ -n "${SPH_E2E_EXE:-}" ]; then E="$SPH_E2E_EXE"
elif [ -x "$ROOT/build/test_e2e_bc6.exe" ]; then E="$ROOT/build/test_e2e_bc6.exe"
else E="$ROOT/build/test_e2e_bc6"; fi

ty=$1; dx=$2; tag="bc6_curve_$ty"
limit="${SPH_TIME_LIMIT:-21600}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-6}"
LOG="$EXP/data/run_$tag.log"

mkdir -p "$EXP/data"
cd "$EXP/data"
t0=$(date +%s)
"$E" --full --contact-type="$ty" --delta-x="$dx" --time-limit="$limit" \
    --emit-curve="$EXP/data/$tag.txt" > "$LOG" 2>&1
rc=$?
t1=$(date +%s)
echo "wall_s=$((t1-t0)) rc=$rc omp=$OMP_NUM_THREADS type=$ty dx=$dx host=$(hostname)" >> "$LOG"
echo "готово: тип $ty, dx=$dx ($tag), код $rc, $((t1-t0)) с"
