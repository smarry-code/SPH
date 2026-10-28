#!/usr/bin/env bash
# Готовые расчёты 2026-10-20: закон 3 + выталкивание tau = 0,5 мкс (сервер, E = 80, 73 м/с);
# Swegle mu 0.1 V0 70 E 78 (локально). GIF ткани, v_z(t) с бегущим маркером, энергия.
set -u
export PYTHONIOENCODING=utf-8
SP=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
SC=$SP/viz_scripts
EXP=/c/dip/sph/SPH/.claude/worktrees/contact-damped/experiments/2026-10-17-guo2019-1t-pinball-vbl
run() { # dir out_dir name law title_law mu v E
  local R=$1 O=$2 name=$3 law=$4 mu=$5 v=$6 E=$7
  mkdir -p "$O"
  echo "== $name vz";     python "$SC/anim_vz.py" "$R" "$O/vz_anim_$name.gif" --law "$law, E = $E ГПа"
  echo "== $name energy"; python "$SC/plot_energy_all.py" "$R/energy_log.csv" "$O" --title "1T base, $law, E = $E ГПа, mu = $mu, V0 = $v м/с"
  echo "== $name gif";    python "$SC/animate_1t.py" --run-dir "$R" --out "$O/anim_1t_$name.gif" --max-frames 40 --snapshots auto --vmax 0.65 --title-prefix "$law, пластичная нить, σ_u = 2,5 ГПа, E = $E ГПа"
}
run "$SP/law3_srv/A_tau5e-07/1T_base_v73.0_r0" "$EXP/anim_law3_push/tau0.5us_mu0.1_v73_E80" impulse_push_tau0.5us_v73_E80 \
    "импульсный (закон 3) + выталкивание τ = 0,5 мкс" 0.1 73 80
R=$(ls -d $SP/local1T/swegle_mu0.1_v70_E78/*/ | head -1); R=${R%/}
run "$R" "$EXP/anim_local_E78/swegle_mu0.1_v70_E78" swegle_mu0.1_v70_E78 "Swegle (закон 2)" 0.1 70 78
echo ALLDONE
