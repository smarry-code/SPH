#!/usr/bin/env bash
# Серверная очередь после перезагрузки 17.10 ~01:47: k0_hertz уже готов,
# k1_pinball и k2_swegle_hold считаются локально.
#   lane A: k2_swegle
#   lane B: k1_pinball_hold
set -u
BIN=$HOME/SPH_ctypes/build-cs/guo2019_shot
FIB=$HOME/ct1T/twaron_plastic_s23.fib
ROOT=$HOME/ct1T
COMMON="--config=1T --level=base --v0=73 --mu=0.1 --t-max=1.5e-3 --n-mon=10 --n-diag=100 --frames=off --fiber-inelastic=$FIB --eps-break=0.05 --contact-stats=10 --threads=4"

run() { # name k extra
  local name=$1 k=$2 extra=${3:-}
  local d=$ROOT/$name
  if [ -f "$d/shots.csv" ]; then echo "$(date +%T) $name: уже готов"; return 0; fi
  rm -rf "$d"; mkdir -p "$d"
  echo "$(date +%T) $name: старт"
  ( cd "$d" && OMP_NUM_THREADS=4 "$BIN" $COMMON --fs-contact=$k --ff-contact=$k $extra --out="$d" > "$d.log" 2>&1 )
  echo "$(date +%T) $name: rc=$?"
}

run k2_swegle 2 > $ROOT/laneA2.log 2>&1 &
run k1_pinball_hold 1 --hold-n > $ROOT/laneB2.log 2>&1 &
wait
echo "$(date +%T) DONE" >> $ROOT/laneA2.log
