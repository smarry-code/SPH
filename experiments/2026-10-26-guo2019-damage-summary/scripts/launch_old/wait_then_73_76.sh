#!/usr/bin/env bash
# Решение владельца 2026-10-22: в очередь 1T 73 и 76 м/с, закон 3 БЕЗ выталкивания, шаг x0,5,
# E = 80 ГПа, mu = 0,1, нить s25; кадры 20 мкс, лог энергии. Старт, когда свободны все 8 ядер
# (по 4 потока на выстрел: на 2 потоках 1T идёт ~в 1,8 раза медленнее).
W=$HOME/vbl1T_law3_dt05_nopush
exec >> $W/wait_then_73_76.log 2>&1
used() {
  local s=0
  for a in $(ps -eo args | grep -E '^/home/claude/SPH_law3/build-elog/guo2019_shot ' | grep -o -- '--threads=[0-9]*' | cut -d= -f2); do
    s=$((s + a))
  done
  echo $s
}
echo "$(date '+%F %T') жду, пока освободятся все 8 ядер (занято $(used))"
while [ $(used) -gt 0 ]; do sleep 60; done
free=$((8 - $(used))); th=$((free / 2)); [ $th -lt 2 ] && th=2
echo "$(date '+%F %T') старт 73 и 76 по $th потока (свободно $free)"
cd $W && python3 law3_fixed.py --work=$W --bin=$HOME/SPH_law3/build-elog/guo2019_shot --fib=$W/twaron_plastic_s25.fib \
  --pair=73,76 --single= --threads=$th --extra="--dt-scale=0.5"
echo "$(date '+%F %T') 73 и 76 готовы, rc=$?"
