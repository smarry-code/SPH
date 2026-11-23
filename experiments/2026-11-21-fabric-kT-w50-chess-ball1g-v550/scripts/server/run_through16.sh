#!/usr/bin/env bash
# 16T шахматной ткани со сдвигом раскройки 2 и стяжками между слоями (варианты 1 + 2, решение владельца 2026-11-22):
# сначала flip, затем same, по одному выстрелу на 8 потоках. Постановка -- run_chess_k.sh (как у 16T без них).
# Вызов: run_through16.sh ROOT БИНАРНИК; готовый выстрел (есть result.txt) не пересчитывается.
set -u
ROOT=${1:?ROOT}
BIN=${2:?БИНАРНИК}
D=$(cd "$(dirname "$0")" && pwd)
for MODE in flip same; do
  W=$ROOT/${MODE}_16T_S2T
  [ -s "$W/result.txt" ] && continue
  "$D/run_chess_k.sh" "$MODE" 16 "$W" "$BIN" "--chess-shift=2 --chess-through"
done
