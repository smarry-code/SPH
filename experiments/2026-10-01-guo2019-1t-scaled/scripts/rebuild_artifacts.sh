#!/usr/bin/env bash
# Пересчёт артефактов, которые не хранятся в каталоге (experiments/README.md §3):
#   anim_1t.gif — анимация выстрела 1T, уровень scaled, V0 = 76 м/с (по кадрам прогона 2).
# PNG-снимки anim_1t_frame{0,21,28,33}.png хранятся в каталоге; скрипт строит их заново рядом с GIF.
#
#   scripts/rebuild_artifacts.sh <корень репозитория> <рабочий каталог вне репозитория> [all|run|anim]
#
#   run  — сборка guo2019_shot и прогон с кадрами раз в 2,5 мкс
#          (Linux, 8 потоков: ~23 мин, 13 058 шагов; кадры 68 + 68, 72 МБ) → <рабочий>/run
#   anim — анимация и снимки по кадрам из <рабочий>/run (Python 3, numpy, matplotlib, Pillow)
#          → <рабочий>/artifacts
#   all  — оба шага (по умолчанию)
#
# Код: коммит 841300c ветки fiber_new или e6e2d21 + repro-diff.zip (см. REPRO.md §3).
# Расчёт с OpenMP побитово не воспроизводится: GIF совпадёт с прежним качественно,
# числа — в пределах REPRO.md §6. Флаг -include cmath нужен g++ 15 (REPRO.md §4);
# переопределяется переменной GUO_CXX_FLAGS.
set -euo pipefail

REPO=$(cd "${1:?укажите корень репозитория}" && pwd)
WORK=${2:?укажите рабочий каталог вне репозитория}
STAGE=${3:-all}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-}
if [ -z "$PY" ]; then   # python3 в Windows может быть заглушкой магазина, проверяем запуском
    for c in python3 python; do "$c" -c "" 2>/dev/null && { PY=$c; break; }; done
fi
[ -n "$PY" ] || { echo "не найден Python 3 (задайте PYTHON=...)" >&2; exit 2; }
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8   # Windows: иначе файлы пишутся в cp1251
CXXF=${GUO_CXX_FLAGS--include cmath}

mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
case "$WORK/" in "$REPO"/*) echo "рабочий каталог должен быть вне репозитория" >&2; exit 2 ;; esac

if [ "$STAGE" = all ] || [ "$STAGE" = run ]; then
    cmake -S "$REPO" -B "$REPO/build-guo" -G Ninja -DCMAKE_BUILD_TYPE=Release "-DCMAKE_CXX_FLAGS=$CXXF"
    cmake --build "$REPO/build-guo" -j8 --target guo2019_shot
    "$REPO/build-guo/guo2019_shot" --config=1T --level=scaled --v0=76 --save-t=2.5e-6 --out="$WORK/run"
fi

if [ "$STAGE" = all ] || [ "$STAGE" = anim ]; then
    mkdir -p "$WORK/artifacts"
    "$PY" "$HERE/animate_1t.py" \
        --run-dir "$WORK/run/1T_scaled_v76.0_r0" --out "$WORK/artifacts/anim_1t.gif" \
        --fps 6 --dpi 55 --max-frames 34 --zoom 9e-3 --zlim=-12,5 --snapshots auto
fi
