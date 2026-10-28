#!/usr/bin/env bash
# Ручной запуск mu=0.3 V0=70 (закон 4) на свободных 4 потоках + перезапуск драйвера.
set -eu
W=$HOME/vbl1T_damped
BIN=$HOME/SPH_damp/build-elog/guo2019_shot
FIB=$W/twaron_plastic_s25.fib
cd "$W"
# 1. Остановить прежний драйвер (только python; идущий guo2019_shot mu=0.1 V0=76 остаётся).
kill 7133 7132 2>/dev/null || true
sleep 2
pgrep -f "window_pinball.py" && { echo "driver still alive"; exit 1; }
pgrep -f -- "--out=$W/mu0.1_v76" >/dev/null || { echo "v76 died!"; exit 1; }
# 2. Выстрел mu=0.3 V0=70 -- те же ключи, что у драйвера.
S=$W/mu0.3_v70
[ -e "$S" ] && { echo "$S exists"; exit 1; }
mkdir -p "$S"
echo "$(date '+%F %T') старт (вручную) mu=0.3 V0=70 (4 потоков, t_max=0.0015)" >> queue.log
OMP_NUM_THREADS=4 setsid nohup "$BIN" --config=1T --level=base --v0=70 --mu=0.3 --out="$S" \
  --frames=off --save-t=2e-05 --threads=4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket \
  --t-max=0.0015 --fiber-inelastic="$FIB" --eps-break=0.05 --fiber-E=9.7e+10 \
  --fs-contact=4 --ff-contact=4 --energy-log=100 > "$S.log" 2>&1 < /dev/null &
sleep 3
# 3. Новый драйвер (ждёт идущие выстрелы, готовые пропускает).
setsid nohup python3 window_pinball.py --work="$W" --bin="$BIN" --fib="$FIB" --contact=4 \
  > driver.out 2>&1 < /dev/null &
sleep 5
ps -eo pid,etime,args | grep -E "guo2019_shot|window_pinball" | grep -v grep
tail -5 queue.log
