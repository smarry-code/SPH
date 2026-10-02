#!/usr/bin/env bash
# Пересчёт артефактов, которые не хранятся в каталоге (experiments/README.md §3):
#   GIF и PNG-снимки (кадры 15 и 40) шести прогонов EF1, EF1-b, EF2, EF4 (mu = 0 и 0,5), EF5;
#   графики ef1_reaction.svg, ef2_reaction.svg.
#
#   scripts/rebuild_artifacts.sh <корень репозитория> <рабочий каталог вне репозитория> [all|run|anim|plots]
#
#   run   — сборка test_e2e_clamp и пять запусков с кадрами (80 на прогон), OMP_NUM_THREADS=1;
#           0–1 с на запуск → <рабочий>/fr, логи → <рабочий>/logs
#   anim  — GIF (60 кадров из 80, dpi 80) и PNG по кадрам из <рабочий>/fr
#           (Python 3, numpy, matplotlib, Pillow; около 1 мин) → <рабочий>/anim
#   plots — SVG по хранимым data/*/reaction.csv (секунды); scripts/ и data/ копируются
#           в <рабочий>/plots, графики и metrics.txt пишутся туда же
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
    cmake -S "$REPO" -B "$REPO/build_clamp" -DCMAKE_BUILD_TYPE=Release "-DCMAKE_CXX_FLAGS=$CXXF"
    cmake --build "$REPO/build_clamp" --target test_e2e_clamp -j8
    B="$REPO/build_clamp/test_e2e_clamp"
    mkdir -p "$WORK/logs"
    cd "$WORK"    # тест запускается не из корня репозитория
    export OMP_NUM_THREADS=1
    "$B" --case=EF1             --solver=pc --frames-dir=fr --frames=80 | tee logs/EF1.log
    "$B" --case=EF1 --variant=b --solver=pc --frames-dir=fr --frames=80 | tee logs/EF1-b.log
    "$B" --case=EF2             --solver=pc --frames-dir=fr --frames=80 | tee logs/EF2.log
    "$B" --case=EF4             --solver=pc --frames-dir=fr --frames=80 | tee logs/EF4.log   # mu = 0 и mu > 0
    "$B" --case=EF5             --solver=pc --frames-dir=fr --frames=80 | tee logs/EF5.log
fi

if want anim; then
    "$PY" "$REPO/tools/clamp_animate.py" --all "$WORK/fr" --out-dir "$WORK/anim" \
        --max-frames 60 --dpi 80 --snapshot 15,40
fi

if want plots; then
    rm -rf "$WORK/plots"
    mkdir -p "$WORK/plots"
    cp -r "$EXP/scripts" "$EXP/data" "$WORK/plots/"
    "$PY" "$WORK/plots/scripts/plot_reaction.py"
fi
