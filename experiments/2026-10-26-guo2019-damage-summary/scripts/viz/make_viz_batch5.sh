#!/usr/bin/env bash
# make_viz_batch5.sh: повторы после заполнения диска (mu = 0,3, закон 3 + выталкивание tau 0,5 мкс, шаг x0,5,
# повреждённость 4 %): 1T 76, 85; 2T 140, 142. GIF ткани с двух сторон (вся ткань), v_z(t), энергия; сводки 1T и 2T.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
E=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments
O1=$E/2026-10-17-guo2019-1t-pinball-vbl/anim_damage
O2=$E/2026-10-18-guo2019-2t-law4/anim_damage
L1=$SP/local1T/win1T_law3push_dmg_mu03
PU="закон 3 + выталкивание τ = 0,5 мкс, шаг ×0,5"
runs=(
  "1T|law3push_mu0.3_v76|$SP/dmg_viz/law3push_mu0.3_v76|/home/claude/vbl1T_law3push_mu03_71_73/mu0.3_v76|$O1"
  "1T|law3push_mu0.3_v85|$SP/dmg_viz/law3push_mu0.3_v85|/home/claude/vbl1T_law3push_mu03_71_73/mu0.3_v85|$O1"
  "2T|2T_law3push_mu0.3_v140|$SP/dmg_viz/2T_law3push_mu0.3_v140|/home/claude/vbl2T_law3push_mu03_140_142/mu0.3_v140|$O2"
  "2T|2T_law3push_mu0.3_v142|$SP/dmg_viz/2T_law3push_mu0.3_v142|/home/claude/vbl2T_law3push_mu03_140_142/mu0.3_v142|$O2"
)
for r in "${runs[@]}"; do
  IFS='|' read -r cfg name dir srv O <<< "$r"
  v=${name##*_v}
  ls $dir/*/shot.txt > /dev/null 2>&1 || bash $SP/fetch_run.sh $srv $dir || { echo "FETCH FAIL $name"; continue; }
  [ -f "$O/$name/anim2_$name.gif" ] && { echo "SKIP $name (уже построен)"; continue; }
  R=$(cygpath -m "$(dirname "$(ls $dir/*/shot.txt | head -1)")")
  D=$O/$name; mkdir -p "$D"; Dm=$(cygpath -m "$D")
  full="$cfg base, $PU, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = 0.3, V0 = $v м/с"
  echo "== $name"
  python "$SC/anim_vz.py" "$R" "$Dm/vz_anim_$name.gif" --law "$PU, повреждённость 4 %, E = 80 ГПа"
  python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$Dm" --title "$full"
  python "$SC/ball_defl.py" "$name=$R" 2>&1 | tail -3
  python "$SC/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --title "$full" || echo "FAIL $name"
  echo "ALLDONE $name"
done
m() { cygpath -m "$1"; }
SUMMARY_TITLE="1T base, $PU, повреждённость 4 %, E = 80 ГПа, μ = 0,3" python "$SC/summary_dt05.py" "$(m $O1)/summary_1T_law3push_dmg_mu03.png" \
  65=$(m $L1/mu0.3_v65) 71=$(m $SP/dmg_viz/law3push_mu0.3_v71) 73=$(m $SP/dmg_viz/law3push_mu0.3_v73) \
  76=$(m $SP/dmg_viz/law3push_mu0.3_v76) 79=$(m $L1/mu0.3_v79) 82=$(m $L1/mu0.3_v82) 85=$(m $SP/dmg_viz/law3push_mu0.3_v85)
SUMMARY_CFG=2T SUMMARY_TMAX=2000 SUMMARY_BREAK_XLIM=150,450 SUMMARY_TITLE="2T base, $PU, повреждённость 4 %, E = 80 ГПа, μ = 0,3" \
  python "$SC/summary_dt05.py" "$(m $O2)/summary_2T_law3push_dmg_mu03.png" \
  120=$(m $SP/dmg_viz/2T_law3push_mu0.3_v120) 135=$(m $SP/dmg_viz/2T_law3push_mu0.3_v135) 140=$(m $SP/dmg_viz/2T_law3push_mu0.3_v140) \
  142=$(m $SP/dmg_viz/2T_law3push_mu0.3_v142) 142.5=$(m $SP/dmg_viz/2T_law3push_mu0.3_v142.5) \
  "144 (2 мс)=$(m $SP/dmg_viz/2T_law3push_mu0.3_v144_long)" 146=$(m $SP/dmg_viz/2T_law3push_mu0.3_v146)
echo ALLDONE_ALL
