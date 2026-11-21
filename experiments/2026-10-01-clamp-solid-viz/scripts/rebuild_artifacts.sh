#!/usr/bin/env bash
# Пересчёт артефактов, которые не хранятся в каталоге (experiments/README.md §3):
#   GIF (24 кадра из 80, dpi 55) и PNG-снимок кадра 16 прогонов ES1 (PC, pull, Эйлер), ES2, ES3
#   (σ_xx и v_z), EF3 (Б, В(а), В(б));
#   графики es1_reaction.svg, es1_variants.svg, es2_cantilever.svg, es3_pinned.svg, ef3_reaction.svg.
#
#   scripts/rebuild_artifacts.sh <корень репозитория> <рабочий каталог вне репозитория> [all|run|anim|plots]
#
#   run   — сборка test_e2e_clamp и запуски с кадрами (80 на прогон), OMP_NUM_THREADS=1, параллельно;
#           на 8 ядрах около 2 мин (ES3 — 101 с), кадры 137 МБ → <рабочий>/frames, логи → <рабочий>/logs
#   anim  — GIF и PNG по кадрам из <рабочий>/frames (Python 3, numpy, matplotlib, Pillow;
#           несколько секунд на прогон) → <рабочий>/anim
#   plots — SVG по хранимым data/*/reaction.csv и disp.csv; scripts/ и data/ копируются
#           в <рабочий>/plots, туда же кладутся логи шага run (plot_reaction.py читает из
#           data/logs/ES1.log полупериод эталона), графики и metrics.txt пишутся туда же
#   all   — все шаги (по умолчанию)
#
# Код: коммит 8ae109b (ветка clamp; входит в историю fiber_new). Расчёт детерминирован
# при OMP_NUM_THREADS=1: reaction.csv должны совпасть с data/ побайтово (REPRO.md §6).
# Флаг -include cmath нужен g++ 15; переопределяется переменной CLAMP_CXX_FLAGS.
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
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8   # Windows: иначе файлы пишутся в cp1251
CXXF=${CLAMP_CXX_FLAGS--include cmath}

mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
case "$WORK/" in "$REPO"/*) echo "рабочий каталог должен быть вне репозитория" >&2; exit 2 ;; esac
want() { [ "$STAGE" = all ] || [ "$STAGE" = "$1" ]; }

if want run; then
    cmake -S "$REPO" -B "$REPO/build" -DCMAKE_BUILD_TYPE=Release "-DCMAKE_CXX_FLAGS=$CXXF"
    cmake --build "$REPO/build" --target test_e2e_clamp -j8
    B="$REPO/build/test_e2e_clamp"
    F="$WORK/frames"
    mkdir -p "$WORK/logs"
    cd "$WORK"
    export OMP_NUM_THREADS=1
    "$B" --case=ES1                --frames-dir="$F" --frames=80 > logs/ES1.log 2>&1 &
    "$B" --case=ES1 --variant=pull --frames-dir="$F" --frames=80 > logs/ES1-pull.log 2>&1 &
    "$B" --case=ES1 --solver=euler --frames-dir="$F" --frames=80 > logs/ES1-euler.log 2>&1 &
    "$B" --case=ES2                --frames-dir="$F" --frames=80 > logs/ES2.log 2>&1 &
    "$B" --case=ES3                --frames-dir="$F" --frames=80 > logs/ES3.log 2>&1 &
    "$B" --case=EF3                --frames-dir="$F" --frames=80 > logs/EF3.log 2>&1 &   # EF3, EF3_2, EF3_3
    "$B" --case=EK1 > logs/EK1.log 2>&1 &                                                 # только лог
    "$B" --case=ES4 > logs/ES4.log 2>&1 &                                                 # отложен, только WARN
    fail=0
    for job in $(jobs -p); do wait "$job" || fail=1; done
    [ "$fail" = 0 ] || { echo "запуск теста завершился с ошибкой, см. $WORK/logs" >&2; exit 1; }
fi

if want anim; then
    F="$WORK/frames"
    A="$WORK/anim"
    mkdir -p "$A"
    for r in ES1_pc_ES1 ES1_pc_pull_ES1 ES1_euler_ES1 ES2_pc_ES2 ES3_pc_ES3 EF3_pc_EF3 EF3_pc_EF3_2 EF3_pc_EF3_3; do
        "$PY" "$REPO/tools/clamp_animate.py" --max-frames 24 --dpi 55 --snapshot 16 \
            --run-dir "$F/$r" --out "$A/$r.gif"
    done
    for r in ES2_pc_ES2 ES3_pc_ES3; do
        "$PY" "$REPO/tools/clamp_animate.py" --max-frames 24 --dpi 55 --snapshot 16 --color-by vz \
            --run-dir "$F/$r" --out "$A/${r}_vz.gif"
    done
fi

if want plots; then
    [ -f "$WORK/logs/ES1.log" ] || { echo "нет $WORK/logs/ES1.log: сначала шаг run" >&2; exit 1; }
    rm -rf "$WORK/plots"
    mkdir -p "$WORK/plots"
    cp -r "$EXP/scripts" "$EXP/data" "$WORK/plots/"
    mkdir -p "$WORK/plots/data/logs"
    cp "$WORK/logs/"*.log "$WORK/plots/data/logs/"
    "$PY" "$WORK/plots/scripts/plot_reaction.py"
fi
