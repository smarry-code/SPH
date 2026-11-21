#!/usr/bin/env bash
# Пересчёт артефактов, которые не хранятся в каталоге (experiments/README.md §3):
#   anim_1t_base.gif — анимация выстрела 1T, уровень base, V0 = 76 м/с.
# PNG-снимки anim_1t_base_frame{0,41,61}.png и full_field_frame{81,122}.png хранятся в каталоге;
# скрипт строит их заново рядом с GIF.
#
#   scripts/rebuild_artifacts.sh <корень репозитория> <рабочий каталог вне репозитория> [all|run|anim]
#
#   run  — сборка guo2019_shot и прогон с кадрами раз в 5 мкс
#          (Linux, 8 потоков: ~2,4 ч, 33 302 шага; кадры 123 + 123, 749 МБ) → <рабочий>/run
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
    "$REPO/build-guo/guo2019_shot" --config=1T --level=base --v0=76 --save-t=5e-6 --out="$WORK/run"
fi

if [ "$STAGE" = all ] || [ "$STAGE" = anim ]; then
    RUN="$WORK/run/1T_base_v76.0_r0"
    mkdir -p "$WORK/artifacts"
    "$PY" "$HERE/animate_1t.py" --run-dir "$RUN" --out "$WORK/artifacts/anim_1t_base.gif" \
        --zoom 30e-3 --yarn-stride 2 --dpi 55 --zlim=-45,12 --vmax 0.5 --elev 30 \
        --max-frames 62 --fps 6 --snapshots auto
    # только снимки всего поля: full_field_frame81.png, full_field_frame122.png
    "$PY" "$HERE/animate_1t.py" --run-dir "$RUN" --out "$WORK/artifacts/full_field.gif" --no-gif \
        --yarn-stride 4 --dpi 60 --zlim=-45,12 --vmax 0.5 --snapshots 81,122
fi
