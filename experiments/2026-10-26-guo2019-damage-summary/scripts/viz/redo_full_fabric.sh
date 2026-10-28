#!/usr/bin/env bash
# redo_full_fabric.sh: перестроить GIF ткани (anim2_*) по ВСЕМ расчётам с повреждённостью -- вся ткань без обрезки,
# закрепление без выделения (полосы убраны по просьбе владельца), яркие нити, D > 0 розовым (владелец, 2026-10-27).
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
E=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments
O1=$E/2026-10-17-guo2019-1t-pinball-vbl/anim_damage
O2=$E/2026-10-18-guo2019-2t-law4/anim_damage
L0="закон 0 (Герц)"; NP="закон 3 без выталкивания, шаг ×0,5"; PU="закон 3 + выталкивание τ = 0,5 мкс, шаг ×0,5"
W=$SP/local1T/win1T_law3push_dmg_mu03
# имя | каталог прогона | конфигурация | закон | μ | папка
runs=(
  "law0_mu0.2_v73|$SP/dmg_viz/law0_mu0.2_v73|1T|$L0|0.2|$O1"
  "law0_mu0.2_v79|$SP/dmg_viz/law0_mu0.2_v79|1T|$L0|0.2|$O1"
  "law0_mu0.3_v73|$SP/dmg_viz/law0_mu0.3_v73|1T|$L0|0.3|$O1"
  "law0_mu0.3_v79|$SP/dmg_viz/law0_mu0.3_v79|1T|$L0|0.3|$O1"
  "law3_mu0.1_v73|$SP/dmg_viz/law3_mu0.1_v73|1T|$NP|0.1|$O1"
  "law3_mu0.25_v71|$SP/dmg_viz/law3_mu0.25_v71|1T|$NP|0.25|$O1"
  "law3_mu0.25_v72|$SP/dmg_viz/law3_mu0.25_v72|1T|$NP|0.25|$O1"
  "law3_mu0.25_v73|$SP/dmg_viz/law3_mu0.25_v73|1T|$NP|0.25|$O1"
  "law3_mu0.25_v75|$SP/dmg_viz/law3_mu0.25_v75|1T|$NP|0.25|$O1"
  "law3_mu0.25_v79|$SP/dmg_viz/law3_mu0.25_v79|1T|$NP|0.25|$O1"
  "law3_mu0.3_v71|$SP/local1T/law3_dmg_mu03_v71|1T|$NP|0.3|$O1"
  "law3_mu0.3_v76|$SP/local1T/law3_dmg_mu03_v76|1T|$NP|0.3|$O1"
  "law3_mu0.3_v79|$SP/local1T/law3_dmg_mu03_v79|1T|$NP|0.3|$O1"
  "law3push_mu0.25_v71|$SP/local1T/law3push_dmg_mu025_v71|1T|$PU|0.25|$O1"
  "law3push_mu0.25_v73|$SP/local1T/law3push_dmg_mu025_v73|1T|$PU|0.25|$O1"
  "law3push_mu0.25_v79|$SP/local1T/law3push_dmg_mu025_v79|1T|$PU|0.25|$O1"
  "law3push_mu0.3_v65|$W/mu0.3_v65|1T|$PU|0.3|$O1"
  "law3push_mu0.3_v71|$SP/dmg_viz/law3push_mu0.3_v71|1T|$PU|0.3|$O1"
  "law3push_mu0.3_v73|$SP/dmg_viz/law3push_mu0.3_v73|1T|$PU|0.3|$O1"
  "law3push_mu0.3_v79|$W/mu0.3_v79|1T|$PU|0.3|$O1"
  "law3push_mu0.3_v80.5|$W/mu0.3_v80.5|1T|$PU|0.3|$O1"
  "law3push_mu0.3_v82|$W/mu0.3_v82|1T|$PU|0.3|$O1"
  "2T_law3_mu0.25_v140|$SP/dmg_viz/2T_law3_mu0.25_v140|2T|$NP|0.25|$O2"
  "2T_law3_mu0.25_v150|$SP/dmg_viz/2T_law3_mu0.25_v150|2T|$NP|0.25|$O2"
  "2T_law3push_mu0.3_v120|$SP/dmg_viz/2T_law3push_mu0.3_v120|2T|$PU|0.3|$O2"
  "2T_law3push_mu0.3_v135|$SP/dmg_viz/2T_law3push_mu0.3_v135|2T|$PU|0.3|$O2"
  "2T_law3push_mu0.3_v142.5|$SP/dmg_viz/2T_law3push_mu0.3_v142.5|2T|$PU|0.3|$O2"
  "2T_law3push_mu0.3_v144|$SP/dmg_viz/2T_law3push_mu0.3_v144|2T|$PU|0.3|$O2"
  "2T_law3push_mu0.3_v144_long|$SP/dmg_viz/2T_law3push_mu0.3_v144_long|2T|$PU|0.3|$O2"
  "2T_law3push_mu0.3_v146|$SP/dmg_viz/2T_law3push_mu0.3_v146|2T|$PU|0.3|$O2"
  "2T_law3push_mu0.3_v150|$SP/dmg_viz/2T_law3push_mu0.3_v150|2T|$PU|0.3|$O2"
)
for r in "${runs[@]}"; do
  IFS='|' read -r name dir cfg law mu O <<< "$r"
  v=${name##*_v}; v=${v%_long}
  R=$(cygpath -m "$(dirname "$(ls $dir/*/shot.txt | head -1)")") || { echo "NO DATA $name"; continue; }
  D=$O/$name; mkdir -p "$D"; rm -f "$D"/anim2_*; Dm=$(cygpath -m "$D")
  full="$cfg base, $law, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = $mu, V0 = $v м/с"
  [ "${name%_long}" != "$name" ] && full="$full (повтор до 2 мс)"
  echo "== $name"
  python "$SC/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --title "$full" || echo "FAIL $name"
  echo "ALLDONE $name"
done
echo ALLDONE_ALL
