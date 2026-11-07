#!/usr/bin/env bash
# make_viz_batch4.sh: готовые выстрелы с повреждённостью 4 % при mu = 0,3, закон 3 + выталкивание tau 0,5 мкс, шаг x0,5:
# 1T 65, 79, 80.5, 82 (ноутбук, серия окна), 71, 73 (сервер); 2T 135, 142.5, 144 (1 мс), 144_long (2 мс), 146 (сервер;
# 120 и 150 построены раньше). GIF ткани с двух сторон, v_z(t), энергия; сводки по 1T и 2T.
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
  "1T|law3push_mu0.3_v65|$L1/mu0.3_v65|-|$O1"
  "1T|law3push_mu0.3_v71|$SP/dmg_viz/law3push_mu0.3_v71|/home/claude/vbl1T_law3push_mu03_71_73/mu0.3_v71|$O1"
  "1T|law3push_mu0.3_v73|$SP/dmg_viz/law3push_mu0.3_v73|/home/claude/vbl1T_law3push_mu03_71_73/mu0.3_v73|$O1"
  "1T|law3push_mu0.3_v79|$L1/mu0.3_v79|-|$O1"
  "1T|law3push_mu0.3_v80.5|$L1/mu0.3_v80.5|-|$O1"
  "1T|law3push_mu0.3_v82|$L1/mu0.3_v82|-|$O1"
  "2T|2T_law3push_mu0.3_v135|$SP/dmg_viz/2T_law3push_mu0.3_v135|/home/claude/vbl2T_law3push_mu03/mu0.3_v135|$O2"
  "2T|2T_law3push_mu0.3_v142.5|$SP/dmg_viz/2T_law3push_mu0.3_v142.5|/home/claude/vbl2T_law3push_mu03/mu0.3_v142.5|$O2"
  "2T|2T_law3push_mu0.3_v144|$SP/dmg_viz/2T_law3push_mu0.3_v144|/home/claude/vbl2T_law3push_mu03/mu0.3_v144|$O2"
  "2T|2T_law3push_mu0.3_v144_long|$SP/dmg_viz/2T_law3push_mu0.3_v144_long|/home/claude/vbl2T_law3push_mu03/mu0.3_v144_long|$O2"
  "2T|2T_law3push_mu0.3_v146|$SP/dmg_viz/2T_law3push_mu0.3_v146|/home/claude/vbl2T_law3push_mu03/mu0.3_v146|$O2"
)
for r in "${runs[@]}"; do
  IFS='|' read -r cfg name dir srv O <<< "$r"
  v=${name##*_v}; v=${v%_long}
  if [ "$srv" != "-" ]; then
    ls $dir/*/shot.txt > /dev/null 2>&1 || bash $SP/fetch_run.sh $srv $dir || { echo "FETCH FAIL $name"; continue; }
  fi
  [ -f "$O/$name/anim2_$name.gif" ] && { echo "SKIP $name (уже построен)"; continue; }
  R=$(cygpath -m "$(dirname "$(ls $dir/*/shot.txt | head -1)")")
  D=$O/$name; mkdir -p "$D"; Dm=$(cygpath -m "$D")
  full="$cfg base, $PU, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = 0.3, V0 = $v м/с"
  [ "${name%_long}" != "$name" ] && full="$full (повтор до 2 мс)"
  echo "== $name"
  python "$SC/anim_vz.py" "$R" "$Dm/vz_anim_$name.gif" --law "$PU, повреждённость 4 %, E = 80 ГПа"
  python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$Dm" --title "$full"
  python "$SC/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --title "$full"
  echo "ALLDONE $name"
done
m() { cygpath -m "$1"; }
SUMMARY_TITLE="1T base, $PU, повреждённость 4 %, E = 80 ГПа, μ = 0,3" python "$SC/summary_dt05.py" "$(m $O1)/summary_1T_law3push_dmg_mu03.png" \
  65=$(m $L1/mu0.3_v65) 71=$(m $SP/dmg_viz/law3push_mu0.3_v71) 73=$(m $SP/dmg_viz/law3push_mu0.3_v73) \
  79=$(m $L1/mu0.3_v79) 80.5=$(m $L1/mu0.3_v80.5) 82=$(m $L1/mu0.3_v82)
SUMMARY_CFG=2T SUMMARY_TMAX=2000 SUMMARY_BREAK_XLIM=150,450 SUMMARY_TITLE="2T base, $PU, повреждённость 4 %, E = 80 ГПа, μ = 0,3" \
  python "$SC/summary_dt05.py" "$(m $O2)/summary_2T_law3push_dmg_mu03.png" \
  120=$(m $SP/dmg_viz/2T_law3push_mu0.3_v120) 135=$(m $SP/dmg_viz/2T_law3push_mu0.3_v135) 142.5=$(m $SP/dmg_viz/2T_law3push_mu0.3_v142.5) \
  "144 (1 мс)=$(m $SP/dmg_viz/2T_law3push_mu0.3_v144)" "144 (2 мс)=$(m $SP/dmg_viz/2T_law3push_mu0.3_v144_long)" \
  146=$(m $SP/dmg_viz/2T_law3push_mu0.3_v146) 150=$(m $SP/dmg_viz/2T_law3push_mu0.3_v150)
echo ALLDONE_ALL
