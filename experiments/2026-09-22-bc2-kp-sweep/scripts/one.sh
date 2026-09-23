#!/usr/bin/env bash
# Один прогон развёртки: one.sh <contact-type> <Kp>
# Требует собранного build/test_e2e_bc2.exe с временной заплатой SPH_BC2_KP
# (см. REPRO.md, раздел «Заплата»).
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXP="$(dirname "$HERE")"
ROOT="$(cd "$EXP/../.." && pwd)"
E="${SPH_E2E_EXE:-$ROOT/build/test_e2e_bc2.exe}"

ty=$1; kp=$2
tag="t${ty}_kp$(echo "$kp" | tr -d '.')"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

mkdir -p "$EXP/data"
cd "$EXP/data"
SPH_BC2_KP=$kp "$E" --full --contact-type="$ty" --time-limit=21600 \
    --emit-curve="$EXP/data/curve_$tag.txt" > "$EXP/data/run_$tag.log" 2>&1
echo "готово: тип $ty Kp=$kp (код $?)"
