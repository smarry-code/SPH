#!/usr/bin/env bash
# Замеры и профили контакта нитей на claude-compute (KVM, 2 ядра, без аппаратных счётчиков): Guo 1T base,
# V0 = 79 м/с, до t = 40 мкс; ключи как в 2026-11-11/scripts/bench_8t.sh. Варианты: H -- закон 0 (Герц);
# P -- закон 3 + проекция (K = 8, тёплый старт). Все прогоны строго по одному.
#   t1, t2   -- время, 1 и 2 потока (1 поток -- ядро 1, 2 потока -- ядра 0-1, OMP_PROC_BIND=close);
#   cnt      -- 2 потока, сборка с диагностическим патчем счётчиков (SPH_PROF_COUNTS -> counts.csv);
#   prof1    -- perf record -e cpu-clock --call-graph fp, 999 Гц, 1 поток (H до 25 мкс, P до 40 мкс);
#   prof2    -- то же, 2 потока, 499 Гц;
#   dwarf2   -- 2 потока, --call-graph dwarf,8192, 100 Гц, запись с начала удара (шаг 1700, задержка по P.t2).
# Запуск: bash prof_compute.sh <src NEW> <src CNT (с патчем)> <fib-файл> <каталог результатов>
# Необязательно: RUNS -- список прогонов (по умолчанию все, в порядке ALL).
set -u
NEW=${1:?исходники NEW}; CNT=${2:?исходники CNT}; FIB=${3:?fib}; R=${4:?каталог результатов}
[ -f "$FIB" ] || { echo "нет файла $FIB" >&2; exit 2; }
mkdir -p "$R" && R=$(cd "$R" && pwd)
[ "$R" != "/" ] && [ "$R" != "$HOME" ] || { echo "каталог результатов $R недопустим" >&2; exit 2; }
log() { echo "$(date "+%F %T") $*" >> "$R/driver.log"; }

build() { # src dir flags...
  local src=$1 dir=$2; shift 2
  cmake -S "$src" -B "$src/$dir" -DCMAKE_BUILD_TYPE=Release "$@" > "$R/cmake_$(basename "$src")_$dir.log" 2>&1 &&
    cmake --build "$src/$dir" -j2 --target guo2019_shot > "$R/build_$(basename "$src")_$dir.log" 2>&1
  log "сборка $src/$dir rc=$? предупреждений: $(grep -c "warning:" "$R/build_$(basename "$src")_$dir.log")"
}
build "$NEW" build
build "$NEW" build-prof "-DCMAKE_CXX_FLAGS=-g -fno-omit-frame-pointer"
build "$CNT" build
B_NEW=$NEW/build/guo2019_shot; B_PROF=$NEW/build-prof/guo2019_shot; B_CNT=$CNT/build/guo2019_shot
for b in "$B_NEW" "$B_PROF" "$B_CNT"; do [ -x "$b" ] || { log "нет $b -- стоп"; exit 1; }; log "  $b md5=$(md5sum "$b" | cut -c1-32)"; done

COMMON="--v0=79 --config=1T --level=base --mu=0.27 --frames=off --save-t=1e-4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --fiber-inelastic=$FIB --eps-break=0.05 --fiber-E=8e10 --dt-scale=0.5 --impulse-push-tau=0"
declare -A V
V[H]="--ff-contact=0 --fs-contact=0"
V[P]="--ff-contact=3 --fs-contact=3 --ff-contact-predictor=1 --fs-contact-predictor=1 --contact-proj-iters=8 --contact-proj-warm=1"
MON=out/1T_base_v79.0_r0/monitor.csv

sampler() { while kill -0 "$2" 2>/dev/null; do
  [ -f "$1" ] && echo "$(date +%s.%N) $(tail -1 "$1" | cut -d, -f1,2)" >> "$3"; sleep 10; done; }

# секунды от старта прогона до шага $2 по wall_samples.txt прогона $1 (линейная интерполяция)
time_to_step() {
  awk -v st="$2" -v t0="$(cat "$1/start.txt")" '{ split($2, a, ","); s = a[1] + 0; t = $1 + 0;
    if (NR == 1 && s >= st) { print t - t0; exit }  # первый же сэмпл уже за целевым шагом
    if (s >= st && ps < st && s > ps) { print t0 == "" ? 0 : (pt + (t - pt) * (st - ps) / (s - ps)) - t0; exit }
    pt = t; ps = s }' "$1/wall_samples.txt"
}

# run ИМЯ ВАРИАНТ ПОТОКИ T_MAX РЕЖИМ(plain|cnt|fp:частота|dwarf:частота:задержка_мс)
run() {
  local name=$1 v=$2 nt=$3 tmax=$4 mode=$5 bin=$B_NEW cpus pre=() env=()
  if [ "$nt" = 1 ]; then cpus=1; else cpus="0-$((nt - 1))"; fi
  local d=$R/$name
  rm -rf -- "$d"; mkdir -p "$d"; cd "$d" || exit 1
  case $mode in
    cnt) bin=$B_CNT; env=(SPH_PROF_COUNTS="$d/counts.csv") ;;
    fp:*) bin=$B_PROF; pre=(perf record -e cpu-clock -F "${mode#fp:}" --call-graph fp -o "$d/perf.data" --) ;;
    dwarf:*) bin=$B_PROF; local f=${mode#dwarf:}; local dl=${f#*:}; f=${f%%:*}
      pre=(perf record -e cpu-clock -F "$f" --call-graph dwarf,8192 -D "$dl" -o "$d/perf.data" --) ;;
  esac
  log "run $name: $v, потоков $nt, t_max $tmax, режим $mode, cpus $cpus"
  local s; s=$(date +%s.%N); echo "$s" > "$d/start.txt"
  env OMP_NUM_THREADS="$nt" OMP_PROC_BIND=close "${env[@]}" "${pre[@]}" \
    taskset -c "$cpus" "$bin" $COMMON ${V[$v]} --threads="$nt" --t-max="$tmax" --out="$d/out" > "$d/stdout.log" 2>&1 &
  local pid=$!
  sampler "$d/$MON" $pid "$d/wall_samples.txt" &
  wait $pid; local rc=$?
  log "  rc=$rc wall=$(echo "$(date +%s.%N) - $s" | bc) $(grep -o 'шагов=[0-9]*, wall=[0-9.]* с ([0-9.]* с/шаг)' "$d/stdout.log" | head -1)"
}

ALL="H.t1 P.t1 H.t2 P.t2 P.cnt H.cnt H.prof1 P.prof1 H.prof2 P.prof2 P.dwarf2"
RUNS=${RUNS:-$ALL}
want() { case " $RUNS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
for name in $ALL; do
  want "$name" || continue
  v=${name%%.*}
  case ${name#*.} in
    t1) run "$name" "$v" 1 4e-5 plain ;;
    t2) run "$name" "$v" 2 4e-5 plain ;;
    cnt) run "$name" "$v" 2 4e-5 cnt ;;
    prof1) [ "$v" = H ] && tm=2.5e-5 || tm=4e-5; run "$name" "$v" 1 "$tm" fp:999 ;;
    prof2) run "$name" "$v" 2 4e-5 fp:499 ;;
    dwarf2)
      [ -s "$R/P.t2/wall_samples.txt" ] && [ -f "$R/P.t2/start.txt" ] ||
        { log "ОШИБКА: нет P.t2 для задержки dwarf2"; continue; }
      dsec=$(time_to_step "$R/P.t2" 1700 2>/dev/null)
      [ -n "$dsec" ] || { log "ОШИБКА: нет P.t2 для задержки dwarf2 (шаг 1700 не достигнут)"; continue; }
      dms=$(printf "%.0f" "$(echo "$dsec * 1000" | bc)")
      log "  задержка записи dwarf2: $dms мс (шаг 1700 по P.t2)"
      run "$name" "$v" 2 4e-5 "dwarf:100:$dms" ;;
  esac
done
# сверки monitor.csv
cmpl() { local a=$R/$1/$MON b=$R/$2/$MON r; [ -f "$a" ] && [ -f "$b" ] || return 0
  cmp -s "$a" "$b" && r=СОВПАДАЕТ || r="ОТЛИЧАЕТСЯ ($(cmp "$a" "$b" | grep -o "line [0-9]*"))"; log "сверка $1 / $2 monitor.csv: $r"; }
for v in H P; do cmpl "$v.t1" "$v.t2"; cmpl "$v.t2" "$v.cnt"; cmpl "$v.t1" "$v.prof1"; cmpl "$v.t2" "$v.prof2"; done
cmpl P.t2 P.dwarf2
if [ -n "${REF8T:-}" ]; then
  for v in H P; do a=$REF8T/$v.t1.monitor.csv b=$R/$v.t1/$MON
    [ -f "$a" ] && [ -f "$b" ] || continue
    cmp -s "$a" "$b" && r=СОВПАДАЕТ || r="ОТЛИЧАЕТСЯ ($(cmp "$a" "$b" | grep -o "line [0-9]*"))"; log "сверка $v.t1 против 8t monitor.csv: $r"
  done
fi
log "DONE"
