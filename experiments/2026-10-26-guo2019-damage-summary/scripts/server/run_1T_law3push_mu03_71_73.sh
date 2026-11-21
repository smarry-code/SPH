#!/usr/bin/env bash
# Владелец 2026-10-26: Guo 1T base, закон 3 С ВЫТАЛКИВАНИЕМ (tau 0,5 мкс, по прямой просьбе владельца), шаг x0,5,
# mu = 0,3, нить s25 + повреждённость (полный разрыв при 4 %), E = 80 ГПа; V0 = 71 и 73 параллельно по 4 потока;
# кадры 20 мкс, лог энергии. Перенесено с ноутбука: там оба выстрела оборвались на ~378 мкс (ноутбук уснул в 08:38).
W=$HOME/vbl1T_law3push_mu03_71_73
B=$HOME/SPH_dmg/build-elog/guo2019_shot
cd $W
for v in 71 73; do
  ( echo "$(date '+%F %T') старт V0=$v" >> queue.log
    OMP_NUM_THREADS=4 $B --config=1T --level=base --v0=$v --mu=0.3 --fs-contact=3 --ff-contact=3 --dt-scale=0.5 \
      --impulse-push-tau=5e-7 --fiber-inelastic=$W/twaron_s25_dmg0.00415.fib --eps-break=0.05 --fiber-E=8e+10 \
      --frames=off --save-t=2e-05 --threads=4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.0015 \
      --energy-log=100 --out=$W/mu0.3_v$v > $W/mu0.3_v$v.log 2>&1
    echo "$(date '+%F %T') готово V0=$v rc=$?" >> queue.log ) &
done
wait
echo "$(date '+%F %T') ALLDONE" >> queue.log
