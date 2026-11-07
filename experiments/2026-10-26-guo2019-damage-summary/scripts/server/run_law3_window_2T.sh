#!/usr/bin/env bash
# Решение владельца 2026-10-24: подбор окна пробития / непробития Guo 2019 2T base на законе 3 (импульсный) БЕЗ
# выталкивания для нить-нить и нить-шар, mu = 0,25 (оба контакта), шаг x0,5, половинное деление в [140; 150] м/с
# (первая пара 140 и 150). Нить s25 + повреждённость: полный разрыв при 4 % (DEPS_F=0.00415), E = 80 ГПа,
# страховочный eps_break 5 %. Кадры 20 мкс, лог энергии. t_max = 1 мс как в прежних 2T, «не определён» -> 2 мс.
# Вердикт 2T -- по шару (not_perforated и V_R = 0), разрывы нитей -- отдельная метрика (владелец, 2026-10-10).
# Стартует после окончания серии 1T mu = 0,25 (~/vbl1T_law3_mu025) и когда на сервере нет guo2019_shot.
W=$HOME/vbl2T_law3_mu025
B=$HOME/SPH_dmg/build-elog/guo2019_shot
exec >> $W/queue_wait.log 2>&1
echo "$(date '+%F %T') жду окончания серии 1T mu = 0,25 и свободного сервера"
while pgrep -f "^python3 window_bisect.py --work=$HOME/vbl1T_law3_mu025" > /dev/null || pgrep -x guo2019_shot > /dev/null; do
  sleep 120
done
echo "$(date '+%F %T') старт серии"
cd $W && python3 window_bisect.py --work=$W --bin=$B --fib=$W/twaron_s25_dmg0.00415.fib --fiber-E=8e10 --eps-break=0.05 \
  --config=2T --contact=3 --mu=0.25 --vmin=140 --v-first=150 --vmax=150 --tol=1.5 --halves --max-rounds=5 \
  --extra=--dt-scale=0.5 --threads=4 --threads-single=8 --t-max=1e-3 --t-max-long=2e-3 --save-t=2e-5 \
  --energy-log=100 > $W/driver.out 2>&1
echo "$(date '+%F %T') серия завершена, rc=$?"
