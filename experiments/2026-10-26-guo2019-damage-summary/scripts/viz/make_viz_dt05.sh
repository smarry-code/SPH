#!/usr/bin/env bash
# make_viz_dt05.sh <run dir (с подкаталогом r0)> <out dir> <name> <V0> <config 1T|2T>
# Закон 3 + выталкивание tau 0,5 мкс, шаг x0,5, E = 80 ГПа, mu = 0,1: GIF ткани, v_z(t), энергия.
set -u
export PYTHONIOENCODING=utf-8
SC=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_scripts
R=$(dirname "$(ls $1/*/shot.txt | head -1)")
O=$2; name=$3; V=$4; CFG=$5
law="закон 3 + выталкивание τ = 0,5 мкс, шаг ×0,5"
mkdir -p "$O"
echo "== $name vz";     python "$SC/anim_vz.py" "$R" "$O/vz_anim_$name.gif" --law "$law, E = 80 ГПа"
echo "== $name energy"; python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "$CFG base, $law, E = 80 ГПа, mu = 0.1, V0 = $V м/с"
echo "== $name gif";    python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_${CFG}_$name.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "$CFG, $law, σ_u = 2,5 ГПа, E = 80 ГПа"
echo "ALLDONE $name"
