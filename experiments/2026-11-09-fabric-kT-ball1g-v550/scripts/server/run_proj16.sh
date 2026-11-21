#!/usr/bin/env bash
# Задача автозапуска fabric_proj16 (claude-compute-8t), решение владельца 2026-11-13: 16T (после 24T -- не пробит),
# закон 3 как проекция скорости (--ff/fs-contact-predictor=1, K = 8, тёплый старт), без выталкивания (с проекцией
# запрещено), шар 1 г, V0 = 550 м/с, mu = 0,27, 8 потоков. Остальное -- как у 24T с проекцией (run_proj24.sh).
# Поиск k на отрезке [12; 24): 12T пробит (закон 3 + выталкивание), 24T не пробит (проекция).
# Чекпоинт раз в 500 шагов; после перезагрузки продолжается из чекпоинта.
set -u
W=${1:-$HOME/fabric_kT/proj16}
cd "$W" || exit 1
env -u OMP_NUM_THREADS python3 "$W/one_shot.py" --work="$W" --bin="$HOME/fabric_kT/src_p/build-elog/guo2019_shot" \
  --fib="$W/twaron_s25_dmg0.00415.fib" --config=16T --contact=3 --v0=550 --mu=0.27 --threads=8 --t-max=4e-4 \
  --fiber-E=8e10 --eps-break=0.05 --level=base --save-t=7.5e-6 --energy-log=100 --checkpoint-every=500 \
  --extra="--dt-scale=0.5 --ball-d=6.253e-3 --ball-mass=1e-3 --ball-dx=3e-4 --ff-contact-predictor=1 --fs-contact-predictor=1" \
  --slot-script="$HOME/proj_runs/threads_free.sh" --proj-status="$HOME/proj_runs/proj_status.sh" \
  --min-uptime=600 --slot-poll=60 --slot-verify-s=120 >> "$W/driver.out" 2>&1
