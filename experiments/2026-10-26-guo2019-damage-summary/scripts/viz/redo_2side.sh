#!/usr/bin/env bash
# redo_2side.sh: перестроить GIF ткани с двух сторон (крупно, ярко, D > 0 розовым) по всем готовым 1T с повреждённостью.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
O=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_damage
L3="закон 3 без выталкивания, шаг ×0,5"; L0="закон 0 (Герц)"
runs=(
  "law3_mu0.25_v71|dmg_viz/law3_mu0.25_v71|$L3|0.25" "law3_mu0.25_v72|dmg_viz/law3_mu0.25_v72|$L3|0.25"
  "law3_mu0.25_v73|dmg_viz/law3_mu0.25_v73|$L3|0.25" "law3_mu0.25_v75|dmg_viz/law3_mu0.25_v75|$L3|0.25"
  "law3_mu0.25_v79|dmg_viz/law3_mu0.25_v79|$L3|0.25"
  "law3_mu0.1_v73|dmg_viz/law3_mu0.1_v73|$L3|0.1"
  "law3_mu0.3_v71|local1T/law3_dmg_mu03_v71|$L3|0.3" "law3_mu0.3_v76|local1T/law3_dmg_mu03_v76|$L3|0.3"
  "law3_mu0.3_v79|local1T/law3_dmg_mu03_v79|$L3|0.3"
  "law0_mu0.2_v73|dmg_viz/law0_mu0.2_v73|$L0|0.2" "law0_mu0.2_v79|dmg_viz/law0_mu0.2_v79|$L0|0.2"
  "law0_mu0.3_v73|dmg_viz/law0_mu0.3_v73|$L0|0.3" "law0_mu0.3_v79|dmg_viz/law0_mu0.3_v79|$L0|0.3"
)
for r in "${runs[@]}"; do
  IFS='|' read -r name dir law mu <<< "$r"
  v=${name##*_v}
  [ "$name" = law3_mu0.25_v72 ] && { ls $SP/$dir/*/shot.txt > /dev/null 2>&1 || bash $SP/fetch_run.sh /home/claude/vbl1T_law3_mu025/mu0.25_v72 $SP/$dir; }
  R=$(cygpath -m "$(dirname "$(ls $SP/$dir/*/shot.txt | head -1)")")
  D=$O/$name; mkdir -p "$D"; rm -f "$D"/anim2_*; Dm=$(cygpath -m "$D")
  full="1T base, $law, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = $mu, V0 = $v м/с"
  if [ ! -f "$D/vz_anim_$name.gif" ]; then
    python "$SP/viz_scripts/anim_vz.py" "$R" "$Dm/vz_anim_$name.gif" --law "$law, повреждённость 4 %, E = 80 ГПа"
    python "$SP/viz_scripts/plot_energy_all.py" "$R/energy_log.csv" "$Dm" --title "$full"
  fi
  python "$SP/viz_scripts/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --title "$full"
  echo "ALLDONE $name"
done
echo ALLDONE_ALL
