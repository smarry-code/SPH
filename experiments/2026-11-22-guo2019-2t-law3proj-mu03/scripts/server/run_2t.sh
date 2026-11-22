#!/usr/bin/env bash
# Задача автозапуска guo2t_proj_mu03 (claude-compute): Guo 2019 2T base, закон 3 на нить-нить и нить-шар
# В РЕЖИМЕ ПРОЕКЦИИ (--ff/fs-contact-predictor=1, K = 8, тёплый старт), без выталкивания, mu = 0,3,
# пара V0 = 144 / 150 м/с (2 выстрела по 4 потока). Остальное -- как в серии 2026-11-22-guo2019-1t-law3proj-mu-sweep:
# нить twaron_s25_dmg0.00415.fib, E = 80 ГПа, eps_break = 0,05, dt-scale 0,5, вердикт по шару,
# t_max 1,5 мс (+ продолжение до 3 мс при «не определён»), кадры раз в 100 мкс.
# Аргументы: W MU VLO VHI (по умолчанию ~/guo2t_proj_mu03 0.3 144 150).
set -u
W=${1:-$HOME/guo2t_proj_mu03}; MU=${2:-0.3}; VLO=${3:-144}; VHI=${4:-150}
BIN=$HOME/guo1t_proj/src/build-elog/guo2019_shot
cd "$W" || exit 1
log() { echo "$(date '+%F %T') $*" >> "$W/sweep.log"; }
log "старт/продолжение: 2T, mu = $MU, V0 = $VLO / $VHI"
env -u OMP_NUM_THREADS python3 "$W/window_bisect.py" --work="$W" --bin="$BIN" \
  --fib="$W/twaron_s25_dmg0.00415.fib" --fiber-E=8e10 --eps-break=0.05 --config=2T --level=base \
  --contact=3 --mu="$MU" --verdict-ball --vmin="$VLO" --v-first="$VHI" --vmax="$VHI" --tol=100 --max-rounds=1 \
  --extra="--dt-scale=0.5 --ff-contact-predictor=1 --fs-contact-predictor=1" --threads=4 --threads-single=4 \
  --t-max=1.5e-3 --t-max-long=3e-3 --save-t=1e-4 --energy-log=100 --checkpoint-every=20000 \
  --slot-script=/nonexistent --proj-status=/nonexistent \
  --min-uptime=600 --slot-poll=300 --slot-verify-s=120 \
  >> "$W/driver.out" 2>&1
rc=$?
log "серия завершена, rc=$rc; $(grep -h 'mu = ' "$W/window.md" 2>/dev/null | tail -1)"
exit $rc
