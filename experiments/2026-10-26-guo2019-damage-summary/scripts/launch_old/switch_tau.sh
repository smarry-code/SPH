#!/usr/bin/env bash
# Решение владельца 2026-10-20: после этапа A остановить программу и перезапустить B и C с tau = 0,5 мкс.
W=$HOME/vbl1T_law3
exec >> $W/switch_tau.log 2>&1
echo "$(date '+%F %T') жду конца этапа A"
until [ -f $W/A_tau2e-06/mu0.1_v73/shots.csv ] && [ -f $W/A_tau5e-07/mu0.1_v73/shots.csv ]; do sleep 10; done
echo "$(date '+%F %T') этап A готов, останавливаю программу"
pkill -f "law3_program.py --work=$W"
sleep 2
pkill -f "guo2019_shot.*--out=$W/B_"
pkill -f "guo2019_shot.*--out=$W/C_"
sleep 5
pgrep -af "guo2019_shot|law3_program" || echo "процессов нет"
[ -f $W/program.md ] && mv $W/program.md $W/program_run1_tau_auto.md
rm -rf $W/B_dt0.5 $W/C_window
mv $W/law3_program.py.new $W/law3_program.py
cd $W
echo "$(date '+%F %T') перезапуск с tau = 5e-7"
setsid nohup python3 law3_program.py --work=$W --bin=$HOME/SPH_law3/build-elog/guo2019_shot --fib=$W/twaron_plastic_s25.fib \
  --fiber-E=8e10 --mu=0.1 --v-check=73 --taus=5e-7 --vmin=70 --vmax=81 --v-first=73 --tol=2 --max-rounds=4 \
  > $W/program2.out 2>&1 < /dev/null &
sleep 30
cat $W/program.md
pgrep -af guo2019_shot | grep -o -- "--out=[^ ]*"
