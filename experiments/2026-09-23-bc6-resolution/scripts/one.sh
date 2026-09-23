#!/usr/bin/env bash
# Один прогон лестницы по разрешению: one.sh <delta_x, м> <тег>
# Тип контакта 0, полный кейс (t_total из tests/data/bc6_initial_conditions.txt),
# геометрия по штатному правилу раскладки — меняется только шаг решётки.
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

dx=$1; tag=$2
limit="${SPH_TIME_LIMIT:-21600}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-6}"

mkdir -p "$EXP/data"
cd "$EXP/data"
t0=$(date +%s)
"$E" --full --contact-type=0 --delta-x="$dx" --time-limit="$limit" \
    --emit-curve="$EXP/data/curve_$tag.txt" > "$EXP/data/run_$tag.log" 2>&1
rc=$?
t1=$(date +%s)
echo "wall_s=$((t1-t0)) rc=$rc omp=$OMP_NUM_THREADS" >> "$EXP/data/run_$tag.log"
echo "готово: dx=$dx ($tag), код $rc, $((t1-t0)) с"
