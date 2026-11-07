#!/usr/bin/env bash
# Ждёт окончания локального 1T (pinball + RATE, mu = 0.1, V0 = 70, E = 78) и строит GIF ткани, v_z(t), энергию.
set -u
export PYTHONIOENCODING=utf-8
SC=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/viz_scripts
R=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad/local1T/pinball_rate_mu0.1_v70_E78/1T_base_v70.0_r0
O=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl/anim_pinball_rate/mu0.1_v70_E78
until [ -f "$R/shot.txt" ] && ! tasklist //FI "IMAGENAME eq guo2019_shot.exe" | grep -q guo2019_shot; do sleep 30; done
echo "run finished: $(date)"
tr ' ' '\n' < "$R/shot.txt" | grep -E '^(outcome|V_R|t_end|n_broken_v|n_broken_h|t_first_break|w_f_at_first_break|max_delta_ff/d|max_delta_fs/d|n_tunnel|W_p|wall_s)='
mkdir -p "$O"
echo "== vz"
python "$SC/anim_vz.py" "$R" "$O/vz_anim_pinball_rate_mu0.1_v70_E78.gif" --law "pinball (закон 1) + RATE, E = 78 ГПа"
echo "== energy"
python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "1T base, pinball + RATE (JC C=0.005), E = 78 ГПа, mu = 0.1, V0 = 70 м/с"
echo "== gif"
python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_1t_pinball_rate_mu0.1_v70_E78.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "pinball + RATE, пластичная нить, σ_u = 2,5 ГПа, E = 78 ГПа"
echo ALLDONE
