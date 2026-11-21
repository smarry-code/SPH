#!/usr/bin/env bash
# Подбор окна 1T base: закон 4, E = 78 ГПа, mu = 0.1, V0 в [71; 81], первая пара 71 / 76.
set -eu
W=$HOME/vbl1T_law4_E78
BIN=$HOME/SPH_damp/build-elog/guo2019_shot
mkdir -p "$W"
cd "$W"
cp $HOME/vbl1T_damped/twaron_plastic_s25.fib "$W/"
mv ~/window_bisect.py ~/window_pinball.py "$W/"
sed -i "s/\r$//" "$W/window_bisect.py" "$W/window_pinball.py"
python3 -m py_compile window_bisect.py window_pinball.py
if pgrep -f guo2019_shot >/dev/null; then echo "busy"; exit 1; fi
setsid nohup python3 window_bisect.py --work="$W" --bin="$BIN" --fib="$W/twaron_plastic_s25.fib" \
  --fiber-E=78e9 --contact=4 --mu=0.1 --vmin=71 --vmax=81 --v-first=76 --tol=2 > driver.out 2>&1 < /dev/null &
sleep 25
ps -eo pid,etime,args | grep -E "guo2019_shot|window_bisect" | grep -v grep | cut -c1-160
cat queue.log
grep -a -m1 "fiber_E" mu0.1_v71.log | cut -c1-200
