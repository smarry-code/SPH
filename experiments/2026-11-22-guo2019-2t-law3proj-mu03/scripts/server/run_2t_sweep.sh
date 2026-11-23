#!/usr/bin/env bash
# Задача автозапуска guo2t_proj_mu (claude-compute): подбор mu для Guo 2019 2T base, закон 3 с проекцией,
# без выталкивания, пара V0 = 144 / 150 м/с на каждое mu, mu по очереди. Каждое mu -- run_2t.sh в подкаталоге W/muX
# (условия -- в шапке run_2t.sh). Аргументы: W MU... (по умолчанию ~/guo2t_proj_mu 0.2 0.27 0.25).
# STOP_ON_WINDOW=1 -- остановиться на первом mu, где 144 не пробит, а 150 пробит (окно).
set -u
W=${1:-$HOME/guo2t_proj_mu}; shift || true
MUS=("$@"); [ ${#MUS[@]} -eq 0 ] && MUS=(0.2 0.27 0.25)
cd "$W" || exit 1
log() { echo "$(date '+%F %T') $*" >> "$W/sweep.log"; }
log "старт/продолжение: 2T, mu = ${MUS[*]}, V0 = 144 / 150"
rc_all=0
for MU in "${MUS[@]}"; do
  D="$W/mu$MU"; mkdir -p "$D"
  for f in window_bisect.py window_pinball.py twaron_s25_dmg0.00415.fib; do [ -e "$D/$f" ] || cp "$W/$f" "$D/"; done
  log "mu = $MU: пара 144 / 150 ($D)"
  bash "$W/run_2t.sh" "$D" "$MU" 144 150
  rc=$?
  log "mu = $MU: rc=$rc; $(grep -h 'mu = ' "$D/window.md" 2>/dev/null | tail -1)"
  [ $rc -ne 0 ] && rc_all=$rc
  if [ "${STOP_ON_WINDOW:-0}" = 1 ] && grep -q "(window)" "$D/window.md" 2>/dev/null; then
    log "mu = $MU: окно найдено, остальные mu не считаю"; break
  fi
done
log "серия завершена, rc=$rc_all"
exit $rc_all
