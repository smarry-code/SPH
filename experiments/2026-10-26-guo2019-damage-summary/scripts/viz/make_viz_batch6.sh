#!/usr/bin/env bash
# make_viz_batch6.sh: серия окна 1T при mu = 0,27 (закон 3 + выталкивание tau 0,5 мкс, шаг x0,5, повреждённость 4 %,
# вердикт по шару): готовые выстрелы 65 (повтор до 3 мс), 79, 82 и, если готов, 85 (повтор до 3 мс).
# GIF ткани с двух сторон (вся ткань), v_z(t), энергия, прогиб под шаром; сводка серии.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
E=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments
O1=$E/2026-10-17-guo2019-1t-pinball-vbl/anim_damage
SRV=/home/claude/vbl1T_law3push_mu027
PU="закон 3 + выталкивание τ = 0,5 мкс, шаг ×0,5"
names=(${@:-65_long 79 82})
sumargs=()
for v in "${names[@]}"; do
  name=law3push_mu0.27_v$v; dir=$SP/dmg_viz/$name
  ls $dir/*/shot.txt > /dev/null 2>&1 || bash $SP/fetch_run.sh $SRV/mu0.27_v$v $dir || { echo "FETCH FAIL $name"; continue; }
  R=$(cygpath -m "$(dirname "$(ls $dir/*/shot.txt | head -1)")")
  vv=${v%_long}; lab=$vv; [ "$vv" != "$v" ] && lab="$vv (до 3 мс)"
  sumargs+=("$lab=$(cygpath -m $dir)")
  [ -f "$O1/$name/anim2_$name.gif" ] && { echo "SKIP $name (уже построен)"; continue; }
  D=$O1/$name; mkdir -p "$D"; Dm=$(cygpath -m "$D")
  full="1T base, $PU, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = 0.27, V0 = $vv м/с"
  [ "$vv" != "$v" ] && full="$full (повтор до 3 мс)"
  echo "== $name"
  python "$SC/anim_vz.py" "$R" "$Dm/vz_anim_$name.gif" --law "$PU, повреждённость 4 %, E = 80 ГПа"
  python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$Dm" --title "$full"
  python "$SC/ball_defl.py" "$name=$R" 2>&1 | tail -3
  python "$SC/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --title "$full" || echo "FAIL $name"
  echo "ALLDONE $name"
done
SUMMARY_TITLE="1T base, $PU, повреждённость 4 %, E = 80 ГПа, μ = 0,27" SUMMARY_TMAX=3000 \
  python "$SC/summary_dt05.py" "$(cygpath -m $O1)/summary_1T_law3push_dmg_mu027.png" "${sumargs[@]}"
echo ALLDONE_ALL
