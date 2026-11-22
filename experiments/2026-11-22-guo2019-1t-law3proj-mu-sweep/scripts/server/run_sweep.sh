#!/usr/bin/env bash
# Задача автозапуска guo1t_proj_mu (claude-compute): подбор коэффициента трения Guo 2019 1T base,
# закон 3 на нить-нить и нить-шар В РЕЖИМЕ ПРОЕКЦИИ (--ff/fs-contact-predictor=1, K = 8, тёплый старт --
# умолчания), без выталкивания (с проекцией оно запрещено). Остальное -- как в серии
# 2026-11-08-guo2019-1t-law3push-mu027: нить twaron_s25_dmg0.00415.fib, E = 80 ГПа, eps_break = 0,05,
# dt-scale 0,5, вердикт по шару, t_max 1,5 мс (+ продолжение до 3 мс при «не определён»), кадры раз в 100 мкс (как в серии C; по слову владельца -- меньше кадров).
# Для каждого mu -- пара V0 = 73 (ожидается непробитие) и 79 (ожидается пробитие), 2 выстрела по 4 потока.
#
# Аргументы: W MU... (по умолчанию ~/guo1t_proj_mu и 0.2 0.3 0.4 0.5). Каждый mu -- свой подкаталог W/muX.
# Готовые выстрелы не пересчитываются, оборванные продолжаются из checkpoint.ckpt (драйвер window_bisect.py).
set -u
W=${1:-$HOME/guo1t_proj_mu}; shift || true
MUS=("$@"); [ ${#MUS[@]} -eq 0 ] && MUS=(0.2 0.3 0.4 0.5)
BIN=$HOME/guo1t_proj/src/build-elog/guo2019_shot
cd "$W" || exit 1
log() { echo "$(date '+%F %T') $*" >> "$W/sweep.log"; }
log "старт/продолжение: mu = ${MUS[*]}, V0 = 73 / 79"
rc_all=0
for MU in "${MUS[@]}"; do
  D="$W/mu$MU"; mkdir -p "$D"
  log "mu = $MU: пара 73 / 79 ($D)"
  env -u OMP_NUM_THREADS python3 "$W/window_bisect.py" --work="$D" --bin="$BIN" \
    --fib="$W/twaron_s25_dmg0.00415.fib" --fiber-E=8e10 --eps-break=0.05 --config=1T --level=base \
    --contact=3 --mu="$MU" --verdict-ball --vmin=73 --v-first=79 --vmax=79 --tol=10 --max-rounds=1 \
    --extra="--dt-scale=0.5 --ff-contact-predictor=1 --fs-contact-predictor=1" --threads=4 --threads-single=4 \
    --t-max=1.5e-3 --t-max-long=3e-3 --save-t=1e-4 --energy-log=100 --checkpoint-every=20000 \
    --slot-script=/nonexistent --proj-status=/nonexistent \
    --min-uptime=600 --slot-poll=300 --slot-verify-s=120 \
    >> "$D/driver.out" 2>&1
  rc=$?
  log "mu = $MU: rc=$rc; $(grep -h 'mu = ' "$D/window.md" 2>/dev/null | tail -1)"
  [ $rc -ne 0 ] && rc_all=$rc
done
log "серия завершена, rc=$rc_all"
exit $rc_all
