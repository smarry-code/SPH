#!/usr/bin/env bash
# Закон 4, серия 1T base: скачать готовые выстрелы, GIF (деформация нити) и энергия по каждому.
set -u
export PYTHONIOENCODING=utf-8
SC=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_scripts
V=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_damped
OUT=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_law4
TP="закон 4 (Герц + демпфер), пластичная нить, σ_u = 2,5 ГПа, E = 97 ГПа"
cd "$V"
for run in mu0.1_v76 mu0.3_v70 mu0.3_v76 mu0.5_v70 mu0.5_v76 mu0.7_v76; do
  echo "== $run download"
  ok=0
  for try in 1 2 3 4 5 6; do
    rm -rf "$run"
    if ssh -o ConnectTimeout=20 -o ServerAliveInterval=15 claude-compute "cd ~/vbl1T_damped && tar czf - $run" | tar xzf - && [ -f "$(ls -d $run/*/ | head -1)energy_log.csv" ]; then ok=1; break; fi
    echo "retry $try $run"; sleep 20
  done
  [ $ok = 1 ] || { echo "FAIL download $run"; continue; }
  mu=${run#mu}; mu=${mu%_v*}; v=${run#*_v}
  R=$(ls -d $V/$run/*/ | head -1); R=${R%/}
  O=$OUT/$run
  mkdir -p "$O"
  echo "== $run energy"
  python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "1T base, закон 4 (Герц + демпфер), mu = $mu, V0 = $v м/с"
  echo "== $run gif"
  python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_1t_law4_$run.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "$TP"
done
echo ALLDONE
