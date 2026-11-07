#!/usr/bin/env bash
# make_viz_nopush.sh <V0> [<V0> ...]: 1T, закон 3 БЕЗ выталкивания, шаг x0,5, E = 80 ГПа, mu = 0,1
# (сервер ~/vbl1T_law3_dt05_nopush). Скачать, GIF ткани, v_z(t), энергия; в конце -- сводный график.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
O=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_law3_nopush
law="закон 3 без выталкивания, шаг ×0,5"
mkdir -p "$O"
for v in "$@"; do
  ls $SP/law3_nopush/v$v/*/shot.txt > /dev/null 2>&1 || bash $SP/fetch_run.sh /home/claude/vbl1T_law3_dt05_nopush/mu0.1_v$v $SP/law3_nopush/v$v || { echo "FETCH FAIL v$v"; continue; }
  R=$(dirname "$(ls $SP/law3_nopush/v$v/*/shot.txt | head -1)")
  mkdir -p "$O/v$v"
  echo "== v$v vz";     python "$SC/anim_vz.py" "$R" "$O/v$v/vz_anim_1T_v$v.gif" --law "$law, E = 80 ГПа"
  echo "== v$v energy"; python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O/v$v" --title "1T base, $law, E = 80 ГПа, mu = 0.1, V0 = $v м/с"
  echo "== v$v gif";    python "$SC/animate_1t.py" --run-dir "$R" --out "$O/v$v/anim_1T_v$v.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "1T, $law, σ_u = 2,5 ГПа, E = 80 ГПа"
  echo "ALLDONE v$v"
done
W=$(cygpath -m $SP)
args=()
for d in $SP/law3_nopush/v*; do v=${d##*/v}; args+=("$v=$(cygpath -m $d)"); done
python "$SC/summary_dt05.py" "$(cygpath -m $O)/summary_1T_law3_nopush.png" "${args[@]}"
echo ALLDONE_ALL
