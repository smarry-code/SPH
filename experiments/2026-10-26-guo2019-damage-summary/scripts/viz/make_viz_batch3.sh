#!/usr/bin/env bash
# make_viz_batch3.sh: 2T mu 0,25 без выталкивания (140, 150), 2T mu 0,3 с выталкиванием (120, 150), 1T mu 0,25 с
# выталкиванием (71, 73, 79; ноутбук). Все -- закон 3, шаг x0,5, повреждённость 4 %. GIF ткани с двух сторон (крупно),
# v_z(t), энергия.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
E=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments
O1=$E/2026-10-17-guo2019-1t-pinball-vbl/anim_damage
O2=$E/2026-10-18-guo2019-2t-law4/anim_damage
NP="закон 3 без выталкивания, шаг ×0,5"; PU="закон 3 + выталкивание τ = 0,5 мкс, шаг ×0,5"
runs=(
  "2T|2T_law3_mu0.25_v140|dmg_viz/2T_law3_mu0.25_v140|/home/claude/vbl2T_law3_mu025/mu0.25_v140|$NP|0.25|$O2"
  "2T|2T_law3_mu0.25_v150|dmg_viz/2T_law3_mu0.25_v150|/home/claude/vbl2T_law3_mu025/mu0.25_v150|$NP|0.25|$O2"
  "2T|2T_law3push_mu0.3_v120|dmg_viz/2T_law3push_mu0.3_v120|/home/claude/vbl2T_law3push_mu03/mu0.3_v120|$PU|0.3|$O2"
  "2T|2T_law3push_mu0.3_v150|dmg_viz/2T_law3push_mu0.3_v150|/home/claude/vbl2T_law3push_mu03/mu0.3_v150|$PU|0.3|$O2"
  "1T|law3push_mu0.25_v71|local1T/law3push_dmg_mu025_v71|-|$PU|0.25|$O1"
  "1T|law3push_mu0.25_v73|local1T/law3push_dmg_mu025_v73|-|$PU|0.25|$O1"
  "1T|law3push_mu0.25_v79|local1T/law3push_dmg_mu025_v79|-|$PU|0.25|$O1"
)
for r in "${runs[@]}"; do
  IFS='|' read -r cfg name dir srv law mu O <<< "$r"
  v=${name##*_v}
  if [ "$srv" != "-" ]; then
    ls $SP/$dir/*/shot.txt > /dev/null 2>&1 || bash $SP/fetch_run.sh $srv $SP/$dir || { echo "FETCH FAIL $name"; continue; }
  fi
  R=$(cygpath -m "$(dirname "$(ls $SP/$dir/*/shot.txt | head -1)")")
  D=$O/$name; mkdir -p "$D"; Dm=$(cygpath -m "$D")
  full="$cfg base, $law, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = $mu, V0 = $v м/с"
  echo "== $name"
  python "$SC/anim_vz.py" "$R" "$Dm/vz_anim_$name.gif" --law "$law, повреждённость 4 %, E = 80 ГПа"
  python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$Dm" --title "$full"
  python "$SC/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --title "$full"
  echo "ALLDONE $name"
done
echo ALLDONE_ALL
