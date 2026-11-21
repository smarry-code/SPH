#!/usr/bin/env bash
# Решение владельца 2026-10-22: 1T 80 м/с (закон 3 без выталкивания, шаг x0,5, E = 80 ГПа, mu = 0,1)
# на 4 потоках, когда закончится проверочный 1T 71 с выталкиванием (~/vbl1T_law3_vmax20).
# Параллельно в той же папке уже идут 68 и 71 м/с (по 2 потока, отдельный экземпляр драйвера).
W=$HOME/vbl1T_law3_dt05_nopush
touch $HOME/server_busy.lock
exec >> $W/wait_then_80.log 2>&1
echo "$(date '+%F %T') жду конца проверочного 71 с выталкиванием"
while pgrep -f -- "--out=$HOME/vbl1T_law3_vmax20/" > /dev/null; do sleep 60; done
echo "$(date '+%F %T') старт 1T 80 без выталкивания (4 потока)"
cd $W && python3 law3_fixed.py --work=$W --bin=$HOME/SPH_law3/build-elog/guo2019_shot --fib=$W/twaron_plastic_s25.fib \
  --pair= --single=80 --threads-single=4 --extra="--dt-scale=0.5"
echo "$(date '+%F %T') 1T 80 готов, rc=$?"
rm -f $HOME/server_busy.lock
