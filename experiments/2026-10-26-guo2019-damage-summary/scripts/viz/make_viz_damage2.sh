#!/usr/bin/env bash
# make_viz_damage2.sh: серия закона 3 без выталкивания + повреждённость 4 %, mu = 0,25 (сервер ~/vbl1T_law3_mu025:
# 71, 75, 79, затем 73 -- скачивается, когда готов) и локальный mu = 0,3, 79 м/с. GIF ткани с двух сторон, v_z(t),
# энергия; сводки по mu = 0,25 и по локальной mu = 0,3 (71, 76, 79).
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
O=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_damage
law="закон 3 без выталкивания, шаг ×0,5"
one() {  # имя каталог μ
  local name=$1 dir=$2 mu=$3 v=${1##*_v}
  local R; R=$(cygpath -m "$(dirname "$(ls $dir/*/shot.txt | head -1)")")
  local D=$O/$name; mkdir -p "$D"; local Dm; Dm=$(cygpath -m "$D")
  local full="1T base, $law, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = $mu, V0 = $v м/с"
  echo "== $name vz";     python "$SC/anim_vz.py" "$R" "$Dm/vz_anim_$name.gif" --law "$law, повреждённость 4 %, E = 80 ГПа"
  echo "== $name energy"; python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$Dm" --title "$full"
  echo "== $name 2side";  python "$SC/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --elev 18 --azim 45 --title "$full"
  echo "ALLDONE $name"
}
one law3_mu0.3_v79 $SP/local1T/law3_dmg_mu03_v79 0.3
for v in 71 75 79 73; do
  ls $SP/dmg_viz/law3_mu0.25_v$v/*/shot.txt > /dev/null 2>&1 || bash $SP/fetch_run.sh /home/claude/vbl1T_law3_mu025/mu0.25_v$v $SP/dmg_viz/law3_mu0.25_v$v
  one law3_mu0.25_v$v $SP/dmg_viz/law3_mu0.25_v$v 0.25
done
m() { cygpath -m "$1"; }
SUMMARY_TITLE="1T base, $law, повреждённость 4 %, E = 80 ГПа, μ = 0,25" python "$SC/summary_dt05.py" "$(m $O)/summary_1T_law3_dmg_mu025.png" \
  71=$(m $SP/dmg_viz/law3_mu0.25_v71) 73=$(m $SP/dmg_viz/law3_mu0.25_v73) 75=$(m $SP/dmg_viz/law3_mu0.25_v75) 79=$(m $SP/dmg_viz/law3_mu0.25_v79)
SUMMARY_TITLE="1T base, $law, повреждённость 4 %, E = 80 ГПа, μ = 0,3 (ноутбук)" python "$SC/summary_dt05.py" "$(m $O)/summary_1T_law3_dmg_mu03.png" \
  71=$(m $SP/local1T/law3_dmg_mu03_v71) 76=$(m $SP/local1T/law3_dmg_mu03_v76) 79=$(m $SP/local1T/law3_dmg_mu03_v79)
echo ALLDONE_ALL
