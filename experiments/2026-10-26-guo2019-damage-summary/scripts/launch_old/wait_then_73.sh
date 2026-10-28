#!/usr/bin/env bash
# Решение владельца 2026-10-21: 1T 73 м/с (закон 3 + выталкивание, шаг x0,5) после 2T 142 на сервере.
# ~/server_busy.lock держит сервер «занятым» для координатора 2T, пока 73 не досчитается.
W=$HOME/vbl1T_law3_dt05
touch $HOME/server_busy.lock
exec >> $W/wait_then_73.log 2>&1
echo "$(date '+%F %T') жду конца 2T 142"
while pgrep -x guo2019_shot > /dev/null; do sleep 60; done
echo "$(date '+%F %T') старт 1T 73"
cd $W && python3 law3_fixed.py --work=$W --bin=$HOME/SPH_law3/build-elog/guo2019_shot --fib=$W/twaron_plastic_s25.fib \
  --pair=71,76 --single=78,73 --threads-single=8
echo "$(date '+%F %T') 1T 73 готов, rc=$?"
rm -f $HOME/server_busy.lock
