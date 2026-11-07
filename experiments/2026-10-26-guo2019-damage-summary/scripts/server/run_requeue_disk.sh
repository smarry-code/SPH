#!/usr/bin/env bash
# 2026-10-26: повтор выстрелов, упавших из-за заполненного диска (15:08): 1T 76, 85 и 2T 140, 142 -- закон 3 +
# выталкивание tau 0,5 мкс, шаг x0,5, mu 0,3, повреждённость 4 %, E 80 ГПа, кадры 20 мкс, лог энергии.
# Очередь на 2 места по 4 потока: сначала 1T 76 и 85, затем 2T 140 и 142 (t_max 2 мс).
B=$HOME/SPH_dmg/build-elog/guo2019_shot
W1=$HOME/vbl1T_law3push_mu03_71_73
W2=$HOME/vbl2T_law3push_mu03_140_142
COMMON="--level=base --mu=0.3 --fs-contact=3 --ff-contact=3 --dt-scale=0.5 --impulse-push-tau=5e-7 --eps-break=0.05 \
 --fiber-E=8e+10 --frames=off --save-t=2e-05 --threads=4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --energy-log=100"
jobs_list=("1T 76 $W1 0.0015" "1T 85 $W1 0.0015" "2T 140 $W2 0.002" "2T 142 $W2 0.002")
for j in "${jobs_list[@]}"; do
  set -- $j; cfg=$1; v=$2; W=$3; tmax=$4
  while [ "$(pgrep -c -x guo2019_shot)" -ge 2 ]; do sleep 60; done
  echo "$(date '+%F %T') старт V0=$v (повтор после заполнения диска)" >> $W/queue.log
  ( cd $W && OMP_NUM_THREADS=4 $B --config=$cfg --v0=$v $COMMON --t-max=$tmax \
      --fiber-inelastic=$W/twaron_s25_dmg0.00415.fib --out=$W/mu0.3_v$v > $W/mu0.3_v$v.log 2>&1
    echo "$(date '+%F %T') готово V0=$v rc=$?" >> $W/queue.log ) &
  sleep 30
done
wait
echo "$(date '+%F %T') ALLDONE (повтор)" >> $W1/queue.log
