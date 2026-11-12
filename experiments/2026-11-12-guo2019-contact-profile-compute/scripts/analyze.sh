#!/usr/bin/env bash
# Разбор результатов prof_compute.sh: таблицы времени, профили по потокам и окнам, вызывающие memmove/memset
# (dwarf), сводка счётчиков. Всё строго последовательно, после окончания прогонов.
# Запуск: bash analyze.sh <каталог результатов> <каталог скриптов> <каталог вывода>
set -u
R=${1:?результаты}; S=${2:?скрипты}; O=${3:?вывод}; mkdir -p "$O"
log() { echo "$(date "+%F %T") $*" >> "$O/analyze.log"; }
log "старт"
python3 "$S/bench_tables.py" "$R" > "$O/bench_tables.txt" 2>&1
python3 "$S/bench_tables.py" "$R" --imp 1700-2733 > "$O/bench_tables_imp25.txt" 2>&1
python3 "$S/bench_tables.py" "$R" --imp 2733-999999 > "$O/bench_tables_imp25_40.txt" 2>&1
log "таблицы готовы"
# окна perf (доли записи) из колонки «окна perf»
wins() { awk -v r="$2" '$1 == r { for (i = 1; i <= NF; i++) if ($i ~ /^(pre|imp):/) printf "--win %s ", $i }' "$1"; }
for r in H.prof1 P.prof1 H.prof2 P.prof2; do
  [ -f "$R/$r/perf.data" ] || continue
  w=$(wins "$O/bench_tables_imp25.txt" "$r")
  late=$(awk -v r="$r" '$1 == r { for (i = 1; i <= NF; i++) if ($i ~ /^imp:/) { sub(/^imp:/, "late:", $i); printf "--win %s ", $i } }' "$O/bench_tables_imp25_40.txt")
  [ "$r" = H.prof1 ] && late=""
  python3 "$S/perf_threads.py" "$R/$r/perf.data" $w $late --top 45 > "$O/$r.threads.txt" 2>&1
  log "perf_threads $r ($w $late) rc=$?"
done
# dwarf: вызывающие функций памяти libc
if [ -f "$R/P.dwarf2/perf.data" ]; then
  P=$R/P.dwarf2/perf.data
  perf report -i "$P" --no-children --stdio --sort sym -g none --percent-limit 0.3 > "$O/P.dwarf2.top.txt" 2>/dev/null
  log "dwarf top rc=$?"
  for s in $(grep -oE "__(memmove|memset|memcpy)[A-Za-z0-9_]*" "$O/P.dwarf2.top.txt" | sort -u); do
    perf report -i "$P" --no-children --stdio --sort sym --symbols="$s" -g caller,0.5,callee,function,percent \
      --max-stack 12 > "$O/P.dwarf2.callers.$s.txt" 2>/dev/null
    log "dwarf callers $s rc=$?"
  done
  perf report -i "$P" --no-children --stdio --sort tid,sym -g none --percent-limit 0.3 > "$O/P.dwarf2.by_tid.txt" 2>/dev/null
  python3 "$S/perf_threads.py" "$P" --top 45 > "$O/P.dwarf2.threads.txt" 2>&1
  log "perf_threads dwarf rc=$?"
fi
python3 "$S/counts_summary.py" "$R"/P.cnt/counts.csv "$R"/H.cnt/counts.csv > "$O/counts_summary.txt" 2>&1
log "счётчики rc=$?"
cp "$R/driver.log" "$O/driver.log"
log "DONE"
