#!/usr/bin/env bash
# Ждёт окончания локального 2T (mu=0.1, V0=144, закон 4) и строит: анимацию v_z(t), GIF ткани (2 слоя), энергию.
set -u
export PYTHONIOENCODING=utf-8
SC=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_scripts
R=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/local2T/mu0.1_v144/2T_base_v144.0_r0
O=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-18-guo2019-2t-law4/anim/mu0.1_v144
until [ -f "$R/shot.txt" ] && ! tasklist //FI "IMAGENAME eq guo2019_shot.exe" | grep -q guo2019_shot; do sleep 30; done
echo "run finished: $(date)"
tr ' ' '\n' < "$R/shot.txt" | grep -E '^(outcome|V_R|t_end|n_broken_v|n_broken_h|wall_s)='
mkdir -p "$O"
echo "== vz anim"
python "$SC/anim_vz.py" "$R" "$O/vz_anim_law4_2T_mu0.1_v144.gif"
echo "== energy"
python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "2T base, закон 4 (Герц + демпфер), mu = 0.1, V0 = 144 м/с"
echo "== gif"
python "$SC/animate_ply.py" --run-dir "$R" --out "$O/anim_2t_law4_mu0.1_v144.gif" --color strain --vmax 0.65 --max-frames 41 --snapshots auto
echo ALLDONE
