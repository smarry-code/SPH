#!/usr/bin/env bash
# Поиск числа слоёв k шахматной ткани (решение владельца 2026-11-22), шар 1 г, 550 м/с:
#   известно 16T = П (2026-11-21); выстрел 20T;
#   20T П -> 24T; П -> 28T и далее шагом +4 до первого НП (не выше KMAX);
#   найдено НП -> деление пополам на отрезке [последнее П; первое НП] до соседних k (П при k, НП при k+1);
#   20T НП -> отрезок [16; 20];
#   «не опр.» к 400 мкс -> продолжение из чекпоинта до 1000 мкс (run_chess_k.sh); если и тогда «не опр.» -- вердикт по
#   monitor.csv (shot_summary.py); не определён и там или выстрел упал -> стоп, решает владелец.
# Выстрелы строго по одному, каждый в своём каталоге ROOT/<MODE>_<k>T; готовый (есть result.txt) не пересчитывается,
# после перезагрузки задача продолжает с того же места. Итог -- ROOT/search_<MODE>/search_result.txt, ход -- search.log,
# таблица выстрелов -- shots.txt.
# Вызов: search_chess.sh MODE ROOT БИНАРНИК [KMAX=40]
set -u
MODE=${1:?MODE}
ROOT=${2:?ROOT}
BIN=${3:?БИНАРНИК}
KMAX=${4:-40}
D=$(cd "$(dirname "$0")" && pwd)
S=$ROOT/search_$MODE
mkdir -p "$S"
LOG=$S/search.log
say() { echo "$(date '+%F %T') $*" >> "$LOG"; }

# shot K -> печатает П / НП / X
shot() {
  local k=$1 w=$ROOT/${MODE}_${1}T
  if [ ! -s "$w/result.txt" ]; then
    say "старт ${k}T ($w)"
    "$D/run_chess_k.sh" "$MODE" "$k" "$w" "$BIN"
    [ -s "$w/result.txt" ] || { say "${k}T: нет result.txt (выстрел упал)"; echo X; return; }
  fi
  local line; line=$(python3 "$D/shot_summary.py" "$w" "$k")
  say "итог ${k}T: $(head -1 "$w/result.txt") | $line"
  grep -q "^$k " "$S/shots.txt" 2>/dev/null || echo "$line" >> "$S/shots.txt"
  echo "$line" | awk '{print $2}'
}

finish() { say "КОНЕЦ: $*"; echo "$*" > "$S/search_result.txt"; exit 0; }

lo=16
say "поиск $MODE: известно 16T = П; KMAX=$KMAX"
k=20
while :; do
  v=$(shot $k)
  case $v in
    П)  lo=$k; k=$((k + 4))
        [ $k -le "$KMAX" ] || finish "до ${lo}T все пробиты, НП не найдено (KMAX=$KMAX) -- решает владелец" ;;
    НП) hi=$k; break ;;
    *)  finish "${k}T не определён или упал; k >= $lo -- решает владелец" ;;
  esac
done
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
