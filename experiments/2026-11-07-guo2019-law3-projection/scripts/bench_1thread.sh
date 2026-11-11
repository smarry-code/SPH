#!/usr/bin/env bash
# Однопоточный бенчмарк контакта на Guo 1T base, V0 = 79 м/с, до t = 40 мкс (§7.1 отчёта).
# Варианты: H — закон 0 (Герц), P — закон 3 с проекцией (K = 8, тёплый старт, как C1),
# L3 — закон 3 без проекции (как C2). Прогоны строго по одному, процесс привязан к ядру 1.
# Проход t1 — только время; проход perf — perf record -F 499 (разбивка: bench_perf_categories.py).
# Запуск: nohup bash bench_1thread.sh <SHOT_BIN> <fib-файл ткани> <каталог> > nohup.out 2>&1 &
set -u
B=${1:?guo2019_shot}; FIB=${2:?twaron_s25_dmg0.00415.fib}; R=${3:?каталог результатов}
[ -x "$B" ] || { echo "нет исполняемого $B" >&2; exit 2; }
[ -f "$FIB" ] || { echo "нет файла $FIB" >&2; exit 2; }
mkdir -p "$R" && R=$(cd "$R" && pwd)   # абсолютный путь: rm -rf ниже -- только внутри него
[ "$R" != "/" ] && [ "$R" != "$HOME" ] || { echo "каталог результатов $R недопустим" >&2; exit 2; }
COMMON="--v0=79 --config=1T --level=base --mu=0.27 --frames=off --save-t=1e-4 --threads=1 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=4e-5 --fiber-inelastic=$FIB --eps-break=0.05 --fiber-E=8e10 --dt-scale=0.5 --impulse-push-tau=0"
declare -A V
V[H]="--ff-contact=0 --fs-contact=0"
V[P]="--ff-contact=3 --fs-contact=3 --ff-contact-predictor=1 --fs-contact-predictor=1 --contact-proj-iters=8 --contact-proj-warm=1"
V[L3]="--ff-contact=3 --fs-contact=3"
export OMP_NUM_THREADS=1 OMP_PROC_BIND=true
mkdir -p "$R"
log(){ echo "$(date "+%F %T") $*" >> "$R/driver.log"; }
# раз в 10 с: время и последняя строка monitor.csv (шаг, t) — для с/шаг по окнам
sampler(){ while kill -0 "$2" 2>/dev/null; do
  [ -f "$1" ] && echo "$(date +%s.%N) $(tail -1 "$1" | cut -d, -f1,2)" >> "$3"; sleep 10; done; }
log "start bin md5=$(md5sum "$B" | cut -c1-32)"
for job in H:t1 P:t1 H:perf P:perf L3:perf; do
  v=${job%%:*}; pass=${job##*:}; d=$R/$v.$pass
  rm -rf -- "$d"; mkdir -p "$d"; cd "$d" || exit 1
  log "run $v $pass"; s=$(date +%s.%N)
  if [ "$pass" = perf ]; then
    perf record -F 499 -o "$d/perf.data" -- taskset -c 1 "$B" $COMMON ${V[$v]} --out="$d/out" > "$d/stdout.log" 2>&1 &
  else
    taskset -c 1 "$B" $COMMON ${V[$v]} --out="$d/out" > "$d/stdout.log" 2>&1 &
  fi
  pid=$!; sampler "$d/out/1T_base_v79.0_r0/monitor.csv" $pid "$d/wall_samples.txt" &
  wait $pid; rc=$?
  log "  rc=$rc wall=$(echo "$(date +%s.%N) - $s" | bc)"
done
log "DONE"
