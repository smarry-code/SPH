#!/usr/bin/env bash
# Владелец 2026-10-26: после 1T (71/73/76/85, ~/vbl1T_law3push_mu03_71_73) -- 2T base, закон 3 С ВЫТАЛКИВАНИЕМ
# (tau 0,5 мкс), шаг x0,5, mu = 0,3, нить s25 + повреждённость (полный разрыв при 4 %), E = 80 ГПа; V0 = 140 и 142,
# по 4 потока. Каждый выстрел стартует, когда 76 и 85 уже запущены и свободны 4 ядра.
# t_max = 2 мс сразу: при 144 к 1 мс исход был «не определён» и понадобился повтор до 2 мс.
W=$HOME/vbl2T_law3push_mu03_140_142
W1=$HOME/vbl1T_law3push_mu03_71_73
B=$HOME/SPH_dmg/build-elog/guo2019_shot
cd $W
echo "$(date '+%F %T') жду старта 1T 76 и 85" >> queue.log
until grep -q "старт V0=85" $W1/queue.log 2>/dev/null; do sleep 60; done
for v in 140 142; do
  while [ "$(pgrep -c -x guo2019_shot)" -ge 2 ]; do sleep 60; done
  echo "$(date '+%F %T') старт V0=$v" >> queue.log
  ( OMP_NUM_THREADS=4 $B --config=2T --level=base --v0=$v --mu=0.3 --fs-contact=3 --ff-contact=3 --dt-scale=0.5 \
      --impulse-push-tau=5e-7 --fiber-inelastic=$W/twaron_s25_dmg0.00415.fib --eps-break=0.05 --fiber-E=8e+10 \
      --frames=off --save-t=2e-05 --threads=4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.002 \
      --energy-log=100 --out=$W/mu0.3_v$v > $W/mu0.3_v$v.log 2>&1
    echo "$(date '+%F %T') готово V0=$v rc=$?" >> queue.log ) &
  sleep 30
done
wait
echo "$(date '+%F %T') ALLDONE" >> queue.log
