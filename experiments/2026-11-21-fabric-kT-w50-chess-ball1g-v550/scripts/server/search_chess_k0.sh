#!/usr/bin/env bash
# Поиск числа слоёв k шахматной ткани от стартового K0 (решение владельца 2026-11-22), шар 1 г, 550 м/с:
#   выстрел K0T; П -> K0+4, П -> +4 ... до первого НП (не выше KMAX);
#   НП найдено -> деление пополам на [последнее П; первое НП] до соседних k;
#   K0T НП -> K0-4, НП -> -4 ... до первого П (не ниже 2), затем деление пополам (мой выбор: владелец задал только ветку П);
#   «не опр.» к 400 мкс -> продолжение до 1000 мкс (run_chess_k.sh); всё ещё «не опр.» -> вердикт по monitor.csv
#   (shot_summary.py); не определён и там или выстрел упал -> стоп, решает владелец.
# Выстрелы по одному, каталог ROOT/<MODE>_<k>T; готовый (есть result.txt) не пересчитывается -- после перезагрузки
# задача продолжает с того же места. Итог -- ROOT/search_<MODE>/search_result.txt, ход -- search.log, таблица -- shots.txt.
# Вызов: search_chess_k0.sh MODE ROOT БИНАРНИК "ДОП_КЛЮЧИ" [K0=16] [KMAX=40]
set -u
MODE=${1:?MODE}
ROOT=${2:?ROOT}
BIN=${3:?БИНАРНИК}
EXTRA=${4:-}
K0=${5:-16}
KMAX=${6:-40}
D=$(cd "$(dirname "$0")" && pwd)
S=$ROOT/search_$MODE
mkdir -p "$S"
LOG=$S/search.log
say() { echo "$(date '+%F %T') $*" >> "$LOG"; }

shot() {
  local k=$1 w=$ROOT/${MODE}_${1}T
  if [ ! -s "$w/result.txt" ]; then
    say "старт ${k}T ($w), ключи: $EXTRA"
    "$D/run_chess_k.sh" "$MODE" "$k" "$w" "$BIN" "$EXTRA"
    [ -s "$w/result.txt" ] || { say "${k}T: нет result.txt (выстрел упал)"; echo X; return; }
  fi
  local line; line=$(python3 "$D/shot_summary.py" "$w" "$k")
  say "итог ${k}T: $(head -1 "$w/result.txt") | $line"
  grep -q "^$k " "$S/shots.txt" 2>/dev/null || echo "$line" >> "$S/shots.txt"
  echo "$line" | awk '{print $2}'
}

finish() { say "КОНЕЦ: $*"; echo "$*" > "$S/search_result.txt"; exit 0; }

say "поиск $MODE от ${K0}T, KMAX=$KMAX"
k=$K0
v=$(shot $k)
case $v in
  П)
    lo=$k
    while :; do
      k=$((k + 4))
      [ $k -le "$KMAX" ] || finish "до ${lo}T все пробиты, НП не найдено (KMAX=$KMAX) -- решает владелец"
      v=$(shot $k)
      case $v in
        П)  lo=$k ;;
        НП) hi=$k; break ;;
        *)  finish "${k}T не определён или упал; k >= $lo -- решает владелец" ;;
      esac
    done ;;
  НП)
    hi=$k
    while :; do
      k=$((k - 4)); [ $k -ge 2 ] || k=2
      v=$(shot $k)
      case $v in
        П)  lo=$k; break ;;
        НП) hi=$k; [ $k -gt 2 ] || finish "$MODE: даже 2T не пробит -- решает владелец" ;;
        *)  finish "${k}T не определён или упал; k <= $hi -- решает владелец" ;;
      esac
    done ;;
  *) finish "${k}T не определён или упал -- решает владелец" ;;
esac
while [ $((hi - lo)) -gt 1 ]; do
  k=$(( (lo + hi) / 2 ))
  v=$(shot $k)
  case $v in
    П)  lo=$k ;;
    НП) hi=$k ;;
    *)  finish "${k}T не определён или упал; отрезок [$lo; $hi] -- решает владелец" ;;
  esac
done
finish "$MODE: k = $lo -- ${lo}T пробит, ${hi}T не пробит"
