#!/usr/bin/env bash
# Пересчёт артефактов, которые не хранятся в каталоге (experiments/README.md §3):
#   GIF и PNG выстрелов 69 и 70 м/с (anim2_*.gif ткани с двух сторон, PNG-снимки кадров, vz_anim_*.gif,
#   energy_*.png, summary_1T_law3push_mu027_v69_70.png) -- то, что лежало в anim/ этого каталога.
#
#   scripts/rebuild_artifacts.sh <корень репозитория> <рабочий каталог вне репозитория> [all|run|anim|plots]
#
#   run   -- сборка guo2019_shot (Release, SPH_ENERGY_LOG=ON) в <рабочий>/build и пара 69 / 70 тем же драйвером
#            и с теми же ключами, что задача mu027_v69_70 на claude-compute-8t: 2 выстрела по 4 потока параллельно,
#            до 1,5 мс, затем «не определён» -> продолжение до 3 мс из финального чекпоинта. Linux, 8 ядер:
#            около 5 ч (3,4 ч база + 1,2-1,5 ч продолжение); кадры около 0,8 ГБ на выстрел (база + _long,
#            кадры базы в _long -- жёсткие ссылки). Выход -- <рабочий>/run (queue.log, window.md, mu0.27_v*).
#   anim  -- анимации по <рабочий>/run (Python 3, numpy, matplotlib, Pillow; около 2 мин) -> <рабочий>/anim
#   plots -- SVG и data/metrics.txt по хранимым data/ (секунды): scripts/ и data/ копируются в <рабочий>/plots,
#            графики пишутся туда же (в каталоге эксперимента лежат такие же)
#   all   -- все шаги (по умолчанию)
#
# Код: коммит 9437d7d (ветка guo1t_law3push_mu027; код расчёта = fiber_new ec2439f). Счёт детерминирован при
# любом числе потоков (fiber_new 6f1b63f): monitor.csv, energy_log.csv и кадры должны совпасть с прогоном
# 2026-11-08 побайтово при том же компиляторе (REPRO.md §6). Анимации строились по каждому второму кадру
# (53-54 из 106-107); по полному набору anim_2side.py берёт 40 кадров из всех, поэтому кадры GIF могут
# отличаться выбором, но не содержанием.
set -euo pipefail

REPO=$(cd "${1:?укажите корень репозитория}" && pwd)
WORK=${2:?укажите рабочий каталог вне репозитория}
STAGE=${3:-all}
HERE=$(cd "$(dirname "$0")" && pwd)
EXP=$(dirname "$HERE")
PY=${PYTHON:-}
if [ -z "$PY" ]; then   # python3 в Windows может быть заглушкой магазина, проверяем запуском
    for c in python3 python; do "$c" -c "" 2>/dev/null && { PY=$c; break; }; done
fi
[ -n "$PY" ] || { echo "не найден Python 3 (задайте PYTHON=...)" >&2; exit 2; }
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8

mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
case "$WORK/" in "$REPO"/*) echo "рабочий каталог должен быть вне репозитория" >&2; exit 2 ;; esac
want() { [ "$STAGE" = all ] || [ "$STAGE" = "$1" ]; }

if want run; then
    cmake -S "$REPO" -B "$WORK/build" -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPH_ENERGY_LOG=ON
    cmake --build "$WORK/build" --target guo2019_shot -j8
    R=$WORK/run
    mkdir -p "$R"
    cp "$EXP/scripts/window_bisect.py" "$EXP/scripts/window_pinball.py" "$EXP/twaron_s25_dmg0.00415.fib" "$R/"
    # Как scripts/server/run_series.sh с аргументами "W 69 70 70". Без ~/proj_runs/*.sh драйвер считает свободные
    # потоки сам (8 минус --threads идущих guo2019_shot) и не ждёт чужих серий. OMP_NUM_THREADS драйвер задаёт
    # каждому выстрелу сам (4).
    env -u OMP_NUM_THREADS "$PY" "$R/window_bisect.py" --work="$R" --bin="$WORK/build/guo2019_shot" \
        --fib="$R/twaron_s25_dmg0.00415.fib" --fiber-E=8e10 --eps-break=0.05 --config=1T --level=base \
        --contact=3 --mu=0.27 --verdict-ball --vmin=69 --v-first=70 --vmax=70 --tol=1.5 \
        --extra="--dt-scale=0.5 --impulse-push-tau=5e-7" --threads=4 --threads-single=4 \
        --t-max=1.5e-3 --t-max-long=3e-3 --save-t=2e-5 --energy-log=100 --checkpoint-every=20000 \
        --slot-script="$WORK/no-threads-free.sh" --proj-status="$WORK/no-proj-status.sh" \
        --min-uptime=0 --slot-poll=60 --slot-verify-s=120 2>&1 | tee "$WORK/run.out"
fi

if want anim; then
    PY=$PY DATA="$WORK/run" OUT="$WORK/anim" bash "$EXP/scripts/viz/make_viz_69_70.sh" 69 70
fi

if want plots; then
    P=$WORK/plots
    mkdir -p "$P"
    cp -r "$EXP/scripts" "$EXP/data" "$P/"
    "$PY" "$P/scripts/plot_series.py" "$P"
fi
echo "готово: $WORK"
