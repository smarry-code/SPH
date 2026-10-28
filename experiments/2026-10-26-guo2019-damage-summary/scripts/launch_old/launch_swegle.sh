#!/usr/bin/env bash
# Закон 2 (Swegle) на обоих контактах нитей, mu = 0.1, V0 = 70 и 76 -- для сравнения с pinball (серия vbl1T_pinball).
# Ключи совпадают с window_pinball.py (кадры 20 мкс, лог энергии, нить s25, E = 97 ГПа, t_max = 1.5 мс).
set -eu
W=$HOME/vbl1T_swegle
BIN=$HOME/SPH_damp/build-elog/guo2019_shot
mkdir -p "$W"
cp -n $HOME/vbl1T_damped/twaron_plastic_s25.fib "$W/"
FIB=$W/twaron_plastic_s25.fib
cd "$W"
if pgrep -f guo2019_shot >/dev/null; then echo "busy: guo2019_shot running"; exit 1; fi
for V in 70 76; do
  S=$W/mu0.1_v$V
  [ -e "$S" ] && { echo "$S exists"; exit 1; }
  mkdir -p "$S"
  echo "$(date '+%F %T') старт mu=0.1 V0=$V закон 2 (4 потоков, t_max=0.0015)" >> queue.log
  OMP_NUM_THREADS=4 setsid nohup "$BIN" --config=1T --level=base --v0=$V --mu=0.1 --out="$S" \
    --frames=off --save-t=2e-05 --threads=4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket \
    --t-max=0.0015 --fiber-inelastic="$FIB" --eps-break=0.05 --fiber-E=9.7e+10 \
    --fs-contact=2 --ff-contact=2 --energy-log=100 > "$S.log" 2>&1 < /dev/null &
done
sleep 20
ps -eo pid,etime,args | grep guo2019_shot | grep -v grep | cut -c1-140
grep -a -m3 -E "fs_contact|contact|закон" "$W/mu0.1_v70.log" | head -3
tail -2 "$W/mu0.1_v70.log"
