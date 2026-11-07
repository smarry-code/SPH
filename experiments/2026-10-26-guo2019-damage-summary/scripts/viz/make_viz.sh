#!/usr/bin/env bash
# Визуализация выстрелов mu = 0.1 (70 и 76 м/с) серии pinball 1T base.
set -u
SC=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_scripts
V=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz
OUT=/c/dip/sph/SPH/.claude/worktrees/fiber-contact-types/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim
TP="pinball, пластичная нить, σ_u = 2,5 ГПа, E = 97 ГПа"
cd "$V"
for v in 70 76; do
  R=$V/mu0.1_v$v/1T_base_v$v.0_r0
  O=$OUT/mu0.1_v$v
  mkdir -p "$O"
  echo "== v$v vel"
  python "$SC/animate_vel.py" --run-dir "$R" --out "$O/vel_mu0.1_v$v.gif" --max-frames 40 --snapshots auto
  echo "== v$v 1t strain"
  python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_1t_mu0.1_v$v.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "$TP"
  echo "== v$v 1t eps_p"
  python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_1t_epsp_mu0.1_v$v.gif" --max-frames 40 --snapshots auto --color-by eps_p --vmax 0.005 --title-prefix "$TP"
  echo "== v$v weave"
  python "$SC/animate_weave.py" --run-dir "$R" --out "$O/weave_mu0.1_v$v.gif" --max-frames 40 --snapshots auto
  echo "== v$v energy"
  python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "1T base, pinball, mu = 0.1, V0 = $v м/с"
done
echo "== compare"
python "$SC/animate_vel.py" --compare --run-dir "$V/mu0.1_v70/1T_base_v70.0_r0" --run-dir "$V/mu0.1_v76/1T_base_v76.0_r0" --out "$OUT/compare_mu0.1_v70_v76.gif" --max-frames 40 --snapshots auto
python "$SC/animate_vel.py" --monitor-dirs "$V/mu0.1_v70/1T_base_v70.0_r0" "$V/mu0.1_v76/1T_base_v76.0_r0" --static-out "$OUT/vz_compare_mu0.1.svg"
echo ALLDONE
