#!/usr/bin/env bash
# Один прогон серии XSPH: one.sh <epsilon> <delta_x, м> <тег>
# Тип контакта 0, полный кейс; меняются только ε (через SPH_BC6_EPSILON —
# требует заплаты 01 из repro-diff.zip) и шаг решётки (штатный --delta-x).
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

eps=$1; dx=$2; tag=$3
limit="${SPH_TIME_LIMIT:-21600}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-6}"
LOG="$EXP/data/run_$tag.log"

mkdir -p "$EXP/data"
cd "$EXP/data"
t0=$(date +%s)
SPH_BC6_EPSILON="$eps" "$E" --full --contact-type=0 --delta-x="$dx" --time-limit="$limit" \
    --emit-curve="$EXP/data/curve_$tag.txt" > "$LOG" 2>&1
rc=$?
t1=$(date +%s)
echo "wall_s=$((t1-t0)) rc=$rc omp=$OMP_NUM_THREADS eps=$eps dx=$dx" >> "$LOG"

# Без заплаты 01 бинарь молча игнорирует SPH_BC6_EPSILON и считает с ε из
# конфига (0.25) — такой прогон выглядел бы рабочим, но был бы не тем.
# Значение сравнивается как число: тест печатает его через %g ("0.1", а не "0.10").
got=$(sed -n 's/.*SPH_BC6_EPSILON = \([0-9.eE+-]*\).*/\1/p' "$LOG" | head -1)
if [ -z "$got" ] || ! awk -v a="$got" -v b="$eps" 'BEGIN { exit !(a == b + 0) }'; then
    echo "ОШИБКА: $tag — в логе нет строки о переопределении ε; заплата 01 не наложена?" >&2
    mv "$EXP/data/curve_$tag.txt" "$EXP/data/INVALID_curve_$tag.txt" 2>/dev/null
    exit 3
fi
echo "готово: eps=$eps dx=$dx ($tag), код $rc, $((t1-t0)) с"
