#!/usr/bin/env bash
# Бенчмарк и профиль контакта нитей на claude-compute-8t (8 ядер): Guo 1T base, V0 = 79 м/с, до t = 40 мкс.
# Варианты: H -- закон 0 (Герц); P -- закон 3 + проекция (K = 8, тёплый старт); L3 -- закон 3 без проекции;
# Pb -- P на базовой сборке (коммит до правок по ревью Codex), цена правок.
# Не больше 8 потоков одновременно: все прогоны строго по одному.
#   1) один поток: H, P, L3, Pb, затем повторы H, P (шум); процесс привязан к ядру 2;
#   2) масштабирование: H и P на 2, 4, 8 потоках (taskset 0..N-1, OMP_PROC_BIND=close);
#   3) профиль (сборка -g -fno-omit-frame-pointer): H и P в 1 потоке (до 25 мкс), H и P на 8 потоках (до 40 мкс).
# Запуск: bash bench_8t.sh <src NEW> <src BASE|-> <fib-файл ткани> <каталог результатов>
# Сборки делает сам: NEW/build (Release), NEW/build-prof, BASE/build (только guo2019_shot; BASE = "-" -- без неё).
# Необязательные переменные окружения:
#   RUNS    -- список прогонов через пробел (по умолчанию все: H.t1 P.t1 L3.t1 Pb.t1 H.t1b P.t1b H.t2 P.t2 ... P.prof8);
#   OLD_RES -- каталог результатов прежнего прогона: monitor.csv одноимённых прогонов сравниваются cmp (побитовая
#              регрессия на полной сцене).
set -u
NEW=${1:?исходники NEW}; BASE=${2:?исходники BASE}; FIB=${3:?fib}; R=${4:?каталог результатов}
[ -f "$FIB" ] || { echo "нет файла $FIB" >&2; exit 2; }
mkdir -p "$R" && R=$(cd "$R" && pwd)
[ "$R" != "/" ] && [ "$R" != "$HOME" ] || { echo "каталог результатов $R недопустим" >&2; exit 2; }
log() { echo "$(date "+%F %T") $*" >> "$R/driver.log"; }

# --- сборки ---------------------------------------------------------------------------------------
build() { # src dir flags...
  local src=$1 dir=$2; shift 2
  cmake -S "$src" -B "$src/$dir" -DCMAKE_BUILD_TYPE=Release "$@" > "$R/cmake_$(basename "$src")_$dir.log" 2>&1 &&
    cmake --build "$src/$dir" -j8 --target guo2019_shot > "$R/build_$(basename "$src")_$dir.log" 2>&1
  log "сборка $src/$dir rc=$?"
}
build "$NEW" build
build "$NEW" build-prof "-DCMAKE_CXX_FLAGS=-g -fno-omit-frame-pointer"
[ "$BASE" != - ] && build "$BASE" build
B_NEW=$NEW/build/guo2019_shot; B_PROF=$NEW/build-prof/guo2019_shot; B_BASE=$BASE/build/guo2019_shot
BINS=("$B_NEW" "$B_PROF"); [ "$BASE" != - ] && BINS+=("$B_BASE")
for b in "${BINS[@]}"; do [ -x "$b" ] || { log "нет $b -- стоп"; exit 1; }; log "  $b md5=$(md5sum "$b" | cut -c1-32)"; done

COMMON="--v0=79 --config=1T --level=base --mu=0.27 --frames=off --save-t=1e-4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --fiber-inelastic=$FIB --eps-break=0.05 --fiber-E=8e10 --dt-scale=0.5 --impulse-push-tau=0"
declare -A V
V[H]="--ff-contact=0 --fs-contact=0"
V[P]="--ff-contact=3 --fs-contact=3 --ff-contact-predictor=1 --fs-contact-predictor=1 --contact-proj-iters=8 --contact-proj-warm=1"
V[Pb]="${V[P]}"
V[L3]="--ff-contact=3 --fs-contact=3"

# раз в 10 с: время и последняя строка monitor.csv (шаг, t) -- с/шаг по окнам
sampler() { while kill -0 "$2" 2>/dev/null; do
  [ -f "$1" ] && echo "$(date +%s.%N) $(tail -1 "$1" | cut -d, -f1,2)" >> "$3"; sleep 10; done; }

# run ИМЯ ВАРИАНТ ПОТОКИ T_MAX PERF(0/1/частота)
run() {
  local name=$1 v=$2 nt=$3 tmax=$4 perf=$5 bin=$B_NEW cpus
  [ "$v" = Pb ] && bin=$B_BASE
  [ "$perf" != 0 ] && bin=$B_PROF
  if [ "$nt" = 1 ]; then cpus=2; else cpus="0-$((nt - 1))"; fi
  local d=$R/$name
  rm -rf -- "$d"; mkdir -p "$d"; cd "$d" || exit 1
  log "run $name: $v, потоков $nt, t_max $tmax, perf $perf, cpus $cpus"
  local s; s=$(date +%s.%N)
  if [ "$perf" != 0 ]; then
    OMP_NUM_THREADS=$nt OMP_PROC_BIND=close perf record -F "$perf" --call-graph fp -o "$d/perf.data" -- \
      taskset -c "$cpus" "$bin" $COMMON ${V[$v]} --threads="$nt" --t-max="$tmax" --out="$d/out" > "$d/stdout.log" 2>&1 &
  else
    OMP_NUM_THREADS=$nt OMP_PROC_BIND=close \
      taskset -c "$cpus" "$bin" $COMMON ${V[$v]} --threads="$nt" --t-max="$tmax" --out="$d/out" > "$d/stdout.log" 2>&1 &
  fi
  local pid=$!
  sampler "$d/out/1T_base_v79.0_r0/monitor.csv" $pid "$d/wall_samples.txt" &
  wait $pid; local rc=$?
  log "  rc=$rc wall=$(echo "$(date +%s.%N) - $s" | bc) $(grep -o 'шагов=[0-9]*, wall=[0-9.]* с ([0-9.]* с/шаг)' "$d/stdout.log" | head -1)"
}

ALL="H.t1 P.t1 L3.t1 Pb.t1 H.t1b P.t1b H.t2 P.t2 H.t4 P.t4 H.t8 P.t8 H.prof1 P.prof1 H.prof8 P.prof8"
RUNS=${RUNS:-$ALL}
want() { case " $RUNS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
# имя прогона -> вариант, потоки, t_max, perf
spec() {
  local name=$1 v=${1%%.*} rest=${1#*.}
  case $rest in
    t1|t1b) echo "$v 1 4e-5 0" ;;
    t2|t4|t8) echo "$v ${rest#t} 4e-5 0" ;;
    prof1) echo "$v 1 2.5e-5 999" ;;
    prof8) echo "$v 8 4e-5 499" ;;
  esac
}
# 1) один поток, 2) масштабирование, 3) профиль -- в порядке ALL
for name in $ALL; do
  want "$name" || continue
  [ "${name%%.*}" = Pb ] && [ "$BASE" = - ] && { log "пропуск $name: нет BASE"; continue; }
  run "$name" $(spec "$name")
done
# детерминизм: monitor.csv на 1 и 8 потоках
for v in H P; do
  a=$R/$v.t1/out/1T_base_v79.0_r0/monitor.csv; b=$R/$v.t8/out/1T_base_v79.0_r0/monitor.csv
  [ -f "$a" ] && [ -f "$b" ] || continue
  cmp -s "$a" "$b" && r=СОВПАДАЕТ || r=ОТЛИЧАЕТСЯ; log "детерминизм $v monitor.csv 1 против 8 потоков: $r"
done
# побитовая регрессия против прежнего прогона (OLD_RES)
if [ -n "${OLD_RES:-}" ]; then
  for name in $RUNS; do
    a=$OLD_RES/$name/out/1T_base_v79.0_r0/monitor.csv; b=$R/$name/out/1T_base_v79.0_r0/monitor.csv
    [ -f "$a" ] && [ -f "$b" ] || continue
    cmp -s "$a" "$b" && r=СОВПАДАЕТ || r="ОТЛИЧАЕТСЯ ($(cmp "$a" "$b" | grep -o "line [0-9]*"))"
    log "регрессия $name monitor.csv против $OLD_RES: $r"
  done
fi
log "DONE"
