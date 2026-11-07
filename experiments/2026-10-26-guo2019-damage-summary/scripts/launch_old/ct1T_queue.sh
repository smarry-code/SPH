#!/usr/bin/env bash
# Серия 1T base, μ = 0.1, V0 = 73, пластичная нить s23: две цепочки по 4 потока.
#   lane A: k0 (Герц) -> k2 (Swegle)
#   lane B: k1 (pinball) -> k1 hold-n
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

( run k0_hertz 0; run k2_swegle 2 ) > $ROOT/laneA.log 2>&1 &
( run k1_pinball 1; run k1_pinball_hold 1 --hold-n ) > $ROOT/laneB.log 2>&1 &
wait
echo "$(date +%T) DONE" >> $ROOT/laneA.log
