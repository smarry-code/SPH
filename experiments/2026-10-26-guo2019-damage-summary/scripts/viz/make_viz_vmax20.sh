#!/usr/bin/env bash
# Проверка предела скорости выталкивания (--impulse-push-vmax=20, сервер ~/vbl1T_law3_vmax20):
# скачать 1T 73 и 71 м/с, построить GIF ткани, v_z(t), энергию и сравнение с выстрелами без предела.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
O=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_law3_vmax20
law="закон 3 + выталкивание τ = 0,5 мкс, предел 20 м/с, шаг ×0,5"
mkdir -p "$O"
for v in 73 71; do
  bash $SP/fetch_run.sh /home/claude/vbl1T_law3_vmax20/mu0.1_v$v $SP/law3_vmax20/v$v || { echo "FETCH FAIL v$v"; continue; }
  R=$(dirname "$(ls $SP/law3_vmax20/v$v/*/shot.txt | head -1)")
  mkdir -p "$O/v$v"
  echo "== v$v vz";     python "$SC/anim_vz.py" "$R" "$O/v$v/vz_anim_1T_v$v.gif" --law "$law, E = 80 ГПа"
  echo "== v$v energy"; python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O/v$v" --title "1T base, $law, E = 80 ГПа, mu = 0.1, V0 = $v м/с"
  echo "== v$v gif";    python "$SC/animate_1t.py" --run-dir "$R" --out "$O/v$v/anim_1T_v$v.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "1T, $law, σ_u = 2,5 ГПа, E = 80 ГПа"
  echo "ALLDONE v$v"
done
W=$(cygpath -m $SP)
python "$SC/summary_dt05.py" "$(cygpath -m $O)/cmp_vmax20_vs_nolimit.png" \
  "71 без предела=$W/law3_dt05_srv/v71" "71 предел 20 м/с=$W/law3_vmax20/v71" \
  "73 без предела (4 потока)=$W/law3_srv/B_tau5e-07_dt0.5" "73 без предела (8 потоков)=$W/law3_dt05_srv/v73" \
  "73 предел 20 м/с=$W/law3_vmax20/v73"
python "$SC/breaks_law3.py" "71_vmax20=$W/law3_vmax20/v71" "73_vmax20=$W/law3_vmax20/v73" > "$SP/breaks_vmax20.txt" 2>&1
echo ALLDONE_ALL
