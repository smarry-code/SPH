#!/usr/bin/env bash
# Перезапуск двух драйверов серии vbl1T_law4_E78_rate с --no-long --max-rounds=3 --adaptive.
# Идущие guo2019_shot не трогаются: новый драйвер их дожидается (shot_running) и берёт результат.
set -eu
B=$HOME/vbl1T_law4_E78_rate
BIN=$HOME/SPH_damp/build-elog/guo2019_shot
sed -i "s/\r$//" ~/window_bisect.py
python3 -m py_compile ~/window_bisect.py
mv ~/window_bisect.py "$B/window_bisect.py"
for pid in $(pgrep -f "window_bisect.py --work=$B/"); do kill "$pid"; done
sleep 3
pgrep -f "window_bisect.py --work=$B/" && { echo "drivers still alive"; exit 1; }
n=$(pgrep -fc -- "--out=$B/") || true
echo "shots still running: $n"
[ "$n" = 4 ] || { echo "expected 4 running shots"; exit 1; }
cd "$B"
for MU in 0.1 0.3; do
  W=$B/mu$MU
  setsid nohup python3 window_bisect.py --work="$W" --bin="$BIN" --fib="$B/twaron_plastic_s25_rate.fib" \
    --fiber-E=78e9 --contact=4 --mu=$MU --vmin=71 --vmax=81 --v-first=76 --tol=2 \
    --threads=2 --threads-single=4 --no-long --max-rounds=3 --adaptive --cores=8 >> "$W/driver.out" 2>&1 < /dev/null &
done
sleep 10
ps -eo pid,etime,args | grep -E "guo2019_shot|window_bisect" | grep -v grep | cut -c1-120
for MU in 0.1 0.3; do tail -3 $B/mu$MU/queue.log; done
