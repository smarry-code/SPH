#!/usr/bin/env bash
# Swegle (закон 2), mu = 0.1, V0 = 70/76: скачать, GIF ткани, анимация v_z(t), энергия.
set -u
export PYTHONIOENCODING=utf-8
SC=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_scripts
V=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_swegle
OUT=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_swegle
TP="Swegle (закон 2), пластичная нить, σ_u = 2,5 ГПа, E = 97 ГПа"
mkdir -p "$V" && cd "$V"
for run in mu0.1_v70 mu0.1_v76; do
  echo "== $run download"
  ok=0
  for try in 1 2 3 4 5 6; do
    rm -rf "./$run"
    if ssh -o ConnectTimeout=20 -o ServerAliveInterval=15 claude-compute "cd ~/vbl1T_swegle && tar czf - $run" | tar xzf - && [ -f "$(ls -d $run/*/ | head -1)energy_log.csv" ]; then ok=1; break; fi
    echo "retry $try $run"; sleep 20
  done
  [ $ok = 1 ] || { echo "FAIL download $run"; continue; }
  v=${run#*_v}
  R=$(ls -d $V/$run/*/ | head -1); R=${R%/}
  O=$OUT/$run
  mkdir -p "$O"
  echo "== $run vz"
  python "$SC/anim_vz.py" "$R" "$O/vz_anim_swegle_$run.gif" --law "Swegle (закон 2)"
  echo "== $run energy"
  python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "1T base, Swegle (закон 2), mu = 0.1, V0 = $v м/с"
  echo "== $run gif"
  python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_1t_swegle_$run.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "$TP"
done
echo ALLDONE
