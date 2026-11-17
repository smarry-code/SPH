#!/usr/bin/env bash
# Последовательный запуск metrics -> frames (по 2 процесса) -> plots на claude-compute.
# Использование: bash run_all.sh [ROOT=~/proj_analysis] [OUT=ROOT/../proj_analysis_out]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${1:-$HOME/proj_analysis}"
OUT="${2:-$HOME/proj_analysis_out}"
RUNS="C1 C2 C3 C4 C5 C1_k1 C1_k2 C1_k3"
mkdir -p "$OUT/frames" "$OUT/png"
cd "$HERE"
python3 metrics.py "$ROOT" --out "$OUT" > "$OUT/metrics.log"
# frames.py: не больше 2 процессов одновременно (сервер: 2 ядра)
for r in $RUNS; do echo "$r"; done | xargs -P2 -I{} sh -c \
  "python3 frames.py '$ROOT/{}' --out '$OUT/frames/{}.csv' 2>/dev/null | tee '$OUT/frames/{}.summary'"
python3 plots.py "$ROOT" "$OUT/frames" --out "$OUT/png"
