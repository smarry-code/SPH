#!/usr/bin/env bash
# Локальные 1T (E = 78 ГПа): Swegle mu 0.3 V0 70; импульсный mu 0.1 V0 73 -- GIF ткани, v_z(t), энергия.
set -u
export PYTHONIOENCODING=utf-8
SC=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_scripts
L=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/local1T
OUT=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_local_E78
for spec in "swegle_mu0.3_v70_E78|Swegle (закон 2)|0.3|70" "impulse_mu0.1_v73_E78|импульсный (закон 3)|0.1|73"; do
  IFS='|' read -r run law mu v <<< "$spec"
  R=$(ls -d $L/$run/*/ | head -1); R=${R%/}
  O=$OUT/$run
  mkdir -p "$O"
  echo "== $run vz"
  python "$SC/anim_vz.py" "$R" "$O/vz_anim_$run.gif" --law "$law, E = 78 ГПа"
  echo "== $run energy"
  python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "1T base, $law, E = 78 ГПа, mu = $mu, V0 = $v м/с"
  echo "== $run gif"
  python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_1t_$run.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "$law, пластичная нить, σ_u = 2,5 ГПа, E = 78 ГПа"
done
echo ALLDONE
