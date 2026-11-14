#!/usr/bin/env bash
# Задача автозапуска fabric_w50_search (claude-compute-8t), решение владельца 2026-11-14:
# после 16T = П -- выстрел 20T; если П -- 24T; затем поиск делением на отрезке [последнее П; первое НП].
#   20T НП           -> отрезок [16; 20] -> 18T -> 17T или 19T;
#   20T П, 24T НП    -> отрезок [20; 24] -> 22T -> 21T или 23T;
#   20T П, 24T П     -> НП не найдено до 24T, стоп (решает владелец);
#   «не опр.» или падение выстрела -> стоп (решает владелец).
# Выстрелы идут строго по одному, каждый -- run_w50_k.sh в своём каталоге ~/fabric_kT/w50_<k>T.
# Готовый выстрел (есть result.txt) не пересчитывается: после перезагрузки задача продолжает с того же места.
# Итог серии -- ~/fabric_kT/w50_search/search_result.txt, ход -- search.log.
set -u
S=$HOME/fabric_kT/w50_search
LOG=$S/search.log
say() { echo "$(date '+%F %T') $*" >> "$LOG"; }

# shot K -> печатает П / НП / X (не определён или упал)
shot() {
  local k=$1 w=$HOME/fabric_kT/w50_${1}T
  if [ ! -s "$w/result.txt" ]; then
    say "старт ${k}T ($w)"
    "$S/run_w50_k.sh" "$k" "$w"
    [ -s "$w/result.txt" ] || { say "${k}T: нет result.txt (выстрел упал)"; echo X; return; }
  fi
  local line; line=$(head -1 "$w/result.txt")
  say "итог ${k}T: $line"
  case "$line" in
    *": П ("*)  echo П ;;
    *": НП ("*) echo НП ;;
    *)          echo X ;;
  esac
}

finish() { say "КОНЕЦ: $*"; echo "$*" > "$S/search_result.txt"; exit 0; }

lo=16   # 16T = П (2026-11-14, ~/fabric_kT/w50_16T)
hi=""
say "поиск: известно 16T = П"

v=$(shot 20)
case $v in
  НП) hi=20 ;;
  П)  lo=20
      v=$(shot 24)
      case $v in
        НП) hi=24 ;;
        П)  finish "16T, 20T, 24T пробиты; НП до 24T не найдено, k >= 24 -- решает владелец" ;;
        *)  finish "24T не определён или упал; k >= 20 -- решает владелец" ;;
      esac ;;
  *)  finish "20T не определён или упал; k >= 16 -- решает владелец" ;;
esac

while [ $((hi - lo)) -gt 1 ]; do
  k=$(( (lo + hi) / 2 ))
  v=$(shot $k)
  case $v in
    П)  lo=$k ;;
    НП) hi=$k ;;
    *)  finish "${k}T не определён или упал; отрезок [$lo; $hi) -- решает владелец" ;;
  esac
done
finish "k = $lo: ${lo}T пробит, ${hi}T не пробит"
