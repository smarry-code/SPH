#!/usr/bin/env bash
# Владелец 2026-10-26: добавить к серии 1T 71/73/76 выстрел V0 = 85 (те же условия: закон 3 + выталкивание tau 0,5 мкс,
# шаг x0,5, mu 0,3, повреждённость 4 %, E 80 ГПа). Стартует после старта 76, как только освободятся 4 ядра
# (раньше 2T 140/142, которые ждут старта 85).
W=$HOME/vbl1T_law3push_mu03_71_73
B=$HOME/SPH_dmg/build-elog/guo2019_shot
cd $W
echo "$(date '+%F %T') V0=85 ждёт старта 76 и свободные 4 ядра" >> queue.log
until grep -q "старт V0=76" queue.log; do sleep 60; done
sleep 30
while [ "$(pgrep -c -x guo2019_shot)" -ge 2 ]; do sleep 60; done
echo "$(date '+%F %T') старт V0=85" >> queue.log
OMP_NUM_THREADS=4 $B --config=1T --level=base --v0=85 --mu=0.3 --fs-contact=3 --ff-contact=3 --dt-scale=0.5 \
  --impulse-push-tau=5e-7 --fiber-inelastic=$W/twaron_s25_dmg0.00415.fib --eps-break=0.05 --fiber-E=8e+10 \
  --frames=off --save-t=2e-05 --threads=4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.0015 \
  --energy-log=100 --out=$W/mu0.3_v85 > $W/mu0.3_v85.log 2>&1
echo "$(date '+%F %T') готово V0=85 rc=$?" >> queue.log
