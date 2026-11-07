#!/usr/bin/env bash
# Решение владельца 2026-10-27: подбор окна пробития / непробития Guo 2019 2T base, всё как в серии 1T mu = 0,27:
# закон 3 (импульсный) с выталкиванием tau 0,5 мкс на нить-нить и нить-шар, шаг x0,5, нить s25 + повреждённость
# (полный разрыв при 4 %, DEPS_F=0.00415), E = 80 ГПа, страховочный eps_break 5 %. Вердикт по шару (--verdict-ball:
# |v_z| < 5 % V0 за плоскостью ткани к t_max -> «не пробит (оценка)»). Половинное деление в [120; 150] (первая пара
# 120 и 150, как в 2T mu 0,3, где порог в области 135-144). t_max 1,5 мс, «не определён» -> 3 мс. Кадры 20 мкс,
# лог энергии. Стартует после окончания серии 1T mu 0,27 (иначе ядра заняты вдвойне).
W=$HOME/vbl2T_law3push_mu027
B=$HOME/SPH_dmg/build-elog/guo2019_shot
exec >> $W/queue_wait.log 2>&1
echo "$(date '+%F %T') жду окончания серии 1T mu 0,27"
until grep -q "серия завершена" $HOME/vbl1T_law3push_mu027/queue_wait.log; do
  sleep 120
done
while pgrep -x guo2019_shot > /dev/null; do
  sleep 120
done
echo "$(date '+%F %T') старт серии"
cd $W && python3 window_bisect.py --work=$W --bin=$B --fib=$W/twaron_s25_dmg0.00415.fib --fiber-E=8e10 --eps-break=0.05 \
  --config=2T --verdict-ball --contact=3 --mu=0.27 --vmin=120 --v-first=150 --vmax=150 --tol=1.5 --halves --max-rounds=6 \
  --extra="--dt-scale=0.5 --impulse-push-tau=5e-7" --threads=4 --threads-single=8 --t-max=1.5e-3 --t-max-long=3e-3 \
  --save-t=2e-5 --energy-log=100 > $W/driver.out 2>&1
echo "$(date '+%F %T') серия завершена, rc=$?"
