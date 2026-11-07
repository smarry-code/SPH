#!/usr/bin/env bash
# make_viz_damage.sh: все готовые расчёты 1T с повреждённостью (полный разрыв при 4 %): GIF ткани с лицевой и
# тыльной стороны (anim_2side.py), бегущий v_z(t) (anim_vz.py), графики энергии; в конце -- сводный график.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
O=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_damage
until grep -q ALLFETCHED $SP/dmg_viz/fetch.log; do sleep 30; done
mkdir -p "$O"
# имя | каталог с прогоном | подпись закона | μ
runs=(
  "law3_mu0.1_v73|$SP/dmg_viz/law3_mu0.1_v73|закон 3 без выталкивания, шаг ×0,5|0.1"
  "law3_mu0.3_v71|$SP/local1T/law3_dmg_mu03_v71|закон 3 без выталкивания, шаг ×0,5|0.3"
  "law3_mu0.3_v76|$SP/local1T/law3_dmg_mu03_v76|закон 3 без выталкивания, шаг ×0,5|0.3"
  "law0_mu0.2_v73|$SP/dmg_viz/law0_mu0.2_v73|закон 0 (Герц)|0.2"
  "law0_mu0.2_v79|$SP/dmg_viz/law0_mu0.2_v79|закон 0 (Герц)|0.2"
  "law0_mu0.3_v73|$SP/dmg_viz/law0_mu0.3_v73|закон 0 (Герц)|0.3"
  "law0_mu0.3_v79|$SP/dmg_viz/law0_mu0.3_v79|закон 0 (Герц)|0.3"
)
args=()
for r in "${runs[@]}"; do
  IFS='|' read -r name dir law mu <<< "$r"
  R=$(cygpath -m "$(dirname "$(ls $dir/*/shot.txt | head -1)")")
  v=${name##*_v}
  D=$O/$name; mkdir -p "$D"; Dm=$(cygpath -m "$D")
  full="1T base, $law, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = $mu, V0 = $v м/с"
  echo "== $name vz";     python "$SC/anim_vz.py" "$R" "$Dm/vz_anim_$name.gif" --law "$law, повреждённость 4 %, E = 80 ГПа"
  echo "== $name energy"; python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$Dm" --title "$full"
  echo "== $name 2side";  python "$SC/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --elev 18 --azim 45 --title "$full"
  echo "ALLDONE $name"
  args+=("$name=$(cygpath -m "$dir")")
done
SUMMARY_TITLE="1T base, повреждённость (полный разрыв при 4 %), E = 80 ГПа" python "$SC/summary_dt05.py" "$(cygpath -m $O)/summary_1T_damage.png" "${args[@]}"
echo ALLDONE_ALL
