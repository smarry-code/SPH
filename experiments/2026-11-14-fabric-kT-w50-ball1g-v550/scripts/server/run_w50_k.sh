#!/usr/bin/env bash
# Один выстрел kT на малом полотне (claude-compute-8t), те же параметры, что у run_w50_16T.sh, кроме числа слоёв.
# Вызов: run_w50_k.sh K [каталог]; каталог по умолчанию ~/fabric_kT/w50_${K}T.
# 43 x 43 нити, окно 50 мм, зажим |y| >= 25 мм, выступ 5,22 мм (полотно 50,7 x 60 мм), 8 потоков; закон 3 как
# проекция скорости (--ff/fs-contact-predictor=1, K = 8, тёплый старт), без выталкивания, шар 1 г, V0 = 550 м/с,
# mu = 0,27. Сборка ~/fabric_kT/src_w (86485d6). Чекпоинт раз в 500 шагов; после перезагрузки -- продолжение.
set -u
K=${1:?K}
W=${2:-$HOME/fabric_kT/w50_${K}T}
mkdir -p "$W" && cd "$W" || exit 1
S=$HOME/fabric_kT/w50_search
env -u OMP_NUM_THREADS python3 "$S/one_shot.py" --work="$W" --bin="$HOME/fabric_kT/src_w/build-elog/guo2019_shot" \
  --fib="$S/twaron_s25_dmg0.00415.fib" --config=${K}T --contact=3 --v0=550 --mu=0.27 --threads=8 --t-max=4e-4 \
  --fiber-E=8e10 --eps-break=0.05 --level=base --save-t=7.5e-6 --energy-log=100 --checkpoint-every=500 \
  --extra="--dt-scale=0.5 --ball-d=6.253e-3 --ball-mass=1e-3 --ball-dx=3e-4 --ff-contact-predictor=1 --fs-contact-predictor=1 --fabric-n=43 --fabric-clamp-edge=25e-3 --fabric-overhang-y=5.22e-3" \
  --slot-script="$HOME/proj_runs/threads_free.sh" --proj-status="$HOME/proj_runs/proj_status.sh" \
  --min-uptime=600 --slot-poll=60 --slot-verify-s=120 >> "$W/driver.out" 2>&1
