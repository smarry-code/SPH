#!/usr/bin/env bash
# make_viz_law3_one.sh <run dir with r0 subdir> <out dir> <name> <law title>; закон 3 + выталкивание, E = 80, mu 0.1, 73 м/с.
set -u
export PYTHONIOENCODING=utf-8
SC=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_scripts
R=$(ls -d $1/*/ | head -1); R=${R%/}
O=$2; name=$3; law=$4
mkdir -p "$O"
echo "== $name vz";     python "$SC/anim_vz.py" "$R" "$O/vz_anim_$name.gif" --law "$law, E = 80 ГПа"
echo "== $name energy"; python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "1T base, $law, E = 80 ГПа, mu = 0.1, V0 = 73 м/с"
echo "== $name gif";    python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_1t_$name.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "$law, пластичная нить, σ_u = 2,5 ГПа, E = 80 ГПа"
echo ALLDONE
