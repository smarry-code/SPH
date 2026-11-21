#!/usr/bin/env bash
# tau = 0,5 мкс (решение владельца 2026-10-20): этап A для него готов; B и C -- сразу.
# Выстрел A с tau = 2 мкс не трогаем: он досчитывается сам (дочерний процесс, переживает драйвер).
W=$HOME/vbl1T_law3
exec >> $W/switch_tau.log 2>&1
pkill -f switch_tau.sh
pkill -f "law3_program.py --work=$W"
sleep 3
pkill -f "guo2019_shot.*--out=$W/B_"
pkill -f "guo2019_shot.*--out=$W/C_"
[ -f $W/program.md ] && mv $W/program.md $W/program_run1_tau_auto.md
rm -rf $W/B_dt0.5 $W/C_window
mv $W/law3_program.py.new $W/law3_program.py
cd $W
echo "$(date '+%F %T') перезапуск с tau = 5e-7 (B на 4 потоках, рядом досчитывается A tau = 2e-6)"
setsid nohup python3 law3_program.py --work=$W --bin=$HOME/SPH_law3/build-elog/guo2019_shot --fib=$W/twaron_plastic_s25.fib \
  --fiber-E=8e10 --mu=0.1 --v-check=73 --taus=5e-7 --vmin=70 --vmax=81 --v-first=73 --tol=2 --max-rounds=4 \
  --threads-single=4 > $W/program2.out 2>&1 < /dev/null &
sleep 20
cat $W/program.md
pgrep -af guo2019_shot | grep -o -- "--out=[^ ]*\|--threads=[0-9]*\|--dt-scale=[^ ]*\|--impulse-push-tau=[^ ]*"
