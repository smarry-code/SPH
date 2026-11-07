#!/usr/bin/env bash
# Подбор окна 1T base: закон 4, E = 78 ГПа, RATE JC (s25_rate), mu = 0.1 и 0.3 -- два драйвера параллельно,
# по 2 потока на выстрел (пара = 4 потока на драйвер, всего 8). V0 в [71; 81], первая пара 71 / 76, tol = 2.
set -eu
B=$HOME/vbl1T_law4_E78_rate
BIN=$HOME/SPH_damp/build-elog/guo2019_shot
mkdir -p "$B"
mv ~/twaron_plastic_s25_rate.fib "$B/"
cp $HOME/vbl1T_law4_E78/window_bisect.py $HOME/vbl1T_law4_E78/window_pinball.py "$B/"
if pgrep -f guo2019_shot >/dev/null; then echo "busy"; exit 1; fi
for MU in 0.1 0.3; do
  W=$B/mu$MU
  mkdir -p "$W"
  cd "$B"
  setsid nohup python3 window_bisect.py --work="$W" --bin="$BIN" --fib="$B/twaron_plastic_s25_rate.fib" \
    --fiber-E=78e9 --contact=4 --mu=$MU --vmin=71 --vmax=81 --v-first=76 --tol=2 \
    --threads=2 --threads-single=4 > "$W/driver.out" 2>&1 < /dev/null &
done
sleep 30
ps -eo pid,etime,args | grep -E "guo2019_shot|window_bisect" | grep -v grep | cut -c1-150
for MU in 0.1 0.3; do cat $B/mu$MU/queue.log; grep -a -m1 "fiber_E" $B/mu$MU/mu${MU}_v71.log | cut -c1-220; done
