#!/usr/bin/env bash
# 2T base, закон 3 без выталкивания, шаг x0,5, E = 80 ГПа, mu = 0,1 (ноутбук, local2T/nopush_dt05_v*):
# GIF двухслойной ткани (animate_ply.py), v_z(t), энергия. 150 м/с -- после окончания его расчёта.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
O=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-18-guo2019-2t-law4/anim_law3_nopush
law="закон 3 без выталкивания, шаг ×0,5"
for v in "$@"; do
  d=$SP/local2T/nopush_dt05_v$v
  until ls $d/*/shot.txt > /dev/null 2>&1; do sleep 120; done
  R=$(dirname "$(ls $d/*/shot.txt | head -1)")
  mkdir -p $O/v$v
  echo "== v$v vz";     python $SC/anim_vz.py "$R" "$O/v$v/vz_anim_2T_v$v.gif" --law "$law, E = 80 ГПа"
  echo "== v$v energy"; python $SC/plot_energy_all.py "$R/energy_log.csv" "$O/v$v" --title "2T base, $law, E = 80 ГПа, mu = 0.1, V0 = $v м/с"
  echo "== v$v gif";    python $SC/animate_ply.py --run-dir "$R" --out "$O/v$v/anim_2T_v$v.gif" --color strain --vmax 0.65 --max-frames 41 --snapshots auto
  echo "ALLDONE v$v"
done
