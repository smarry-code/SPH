#!/usr/bin/env bash
# Один выстрел 16T шахматной ткани (решение владельца 2026-11-21): flip -- на claude-compute-8t, same -- на claude-compute.
# Вызов: run_chess16.sh MODE КАТАЛОГ БИНАРНИК
#   MODE      flip | same (--chess=MODE);
#   КАТАЛОГ   рабочий каталог выстрела; рядом со скриптом лежат one_shot.py, window_pinball.py, twaron_s25_dmg0.00415.fib;
#   БИНАРНИК  guo2019_shot из сборки 83ba942 (Release, SPH_ENERGY_LOG=ON).
# Постановка -- как у 16T серии w50 (2026-11-14-fabric-kT-w50-ball1g-v550/scripts/server/run_w50_k.sh), кроме ткани:
#   полотно 43 x 43 нити, окно 50 мм, зажим |y| >= 25 мм, выступ 5,22 мм (50,7 x 60 мм); шахматная ткань --chess=MODE,
#   кусочки 7 x 7 (ядро 5 x 5), подслои через 0,859 мм, стяжки (материал стёжки по умолчанию);
#   шар 1 г, D = 6,253 мм, dx = 0,3 мм; V0 = 550 м/с; закон 3 нить-нить и нить-шар как проекция скорости
#   (--ff/fs-contact-predictor=1, K = 8, тёплый старт), без выталкивания; mu = 0,27; dt-scale 0,5; нить Twaron s25 +
#   повреждённость, E = 80 ГПа, eps_break 0,05; t_max 400 мкс, «не опр.» -> продолжение до 1000 мкс; кадры раз в 20 мкс,
#   energy_log раз в 100 шагов, чекпоинт раз в 500 шагов (после перезагрузки -- продолжение); 8 потоков.
set -u
MODE=${1:?MODE}
W=${2:?КАТАЛОГ}
BIN=${3:?БИНАРНИК}
S=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$W" && cd "$W" || exit 1
env -u OMP_NUM_THREADS python3 "$S/one_shot.py" --work="$W" --bin="$BIN" \
  --fib="$S/twaron_s25_dmg0.00415.fib" --config=16T --contact=3 --v0=550 --mu=0.27 --threads=8 --t-max=4e-4 --t-max-long=1e-3 \
  --fiber-E=8e10 --eps-break=0.05 --level=base --save-t=2e-5 --energy-log=100 --checkpoint-every=500 \
  --extra="--dt-scale=0.5 --ball-d=6.253e-3 --ball-mass=1e-3 --ball-dx=3e-4 --ff-contact-predictor=1 --fs-contact-predictor=1 --fabric-n=43 --fabric-clamp-edge=25e-3 --fabric-overhang-y=5.22e-3 --chess=$MODE" \
  --slot-script="$HOME/proj_runs/threads_free.sh" --proj-status="$HOME/proj_runs/proj_status.sh" \
  --min-uptime=600 --slot-poll=60 --slot-verify-s=120 >> "$W/driver.out" 2>&1
