#!/usr/bin/env bash
# Решение владельца 2026-10-25: подбор окна пробития / непробития Guo 2019 2T base на законе 3 (импульсный) С
# ВЫТАЛКИВАНИЕМ (tau 0,5 мкс, по прямой просьбе владельца) для нить-нить и нить-шар, mu = 0,3, шаг x0,5,
# половинное деление в [120; 150] м/с (при mu 0,25 без выталкивания 140 уже пробит, V_R 35,6)
# (первая пара 120 и 150). Нить s25 + повреждённость: полный разрыв при 4 % (DEPS_F=0.00415), E = 80 ГПа,
# страховочный eps_break 5 %. Кадры 20 мкс, лог энергии. t_max = 1 мс как в прежних 2T, «не определён» -> 2 мс.
# Вердикт 2T -- по шару (not_perforated и V_R = 0), разрывы нитей -- отдельная метрика (владелец, 2026-10-10).
# Стартует, когда на сервере нет guo2019_shot.
W=$HOME/vbl2T_law3push_mu03
B=$HOME/SPH_dmg/build-elog/guo2019_shot
exec >> $W/queue_wait.log 2>&1
echo "$(date '+%F %T') жду свободного сервера"
while pgrep -x guo2019_shot > /dev/null; do
  sleep 120
done
echo "$(date '+%F %T') старт серии"
cd $W && python3 window_bisect.py --work=$W --bin=$B --fib=$W/twaron_s25_dmg0.00415.fib --fiber-E=8e10 --eps-break=0.05 \
  --config=2T --contact=3 --mu=0.3 --vmin=120 --v-first=150 --vmax=150 --tol=1.5 --halves --max-rounds=6 \
  --extra="--dt-scale=0.5 --impulse-push-tau=5e-7" --threads=4 --threads-single=8 --t-max=1e-3 --t-max-long=2e-3 --save-t=2e-5 \
  --energy-log=100 > $W/driver.out 2>&1
echo "$(date '+%F %T') серия завершена, rc=$?"
