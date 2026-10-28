#!/usr/bin/env bash
# Визуализация в ParaView по кадрам прогона 1T scaled, V0 = 76 м/с (см. PARAVIEW.md).
#   scripts/paraview_scene.sh [run|convert|scene|anim|all]   (по умолчанию: convert scene)
#   PVPYTHON=<путь к pvpython> — если pvpython не в PATH
# Всё лежит в data/run (кадры), data/vtk (VTK, сцена), data/csv (таблицы); эти каталоги не в git.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); EXP=$(dirname "$HERE"); REPO=$(cd "$EXP/../.." && pwd)
STAGE=${1:-prep}
PY=${PYTHON:-python}; PV=${PVPYTHON:-pvpython}
RUN=$EXP/data/run/1T_scaled_v76.0_r0; VTK=$EXP/data/vtk; CSV=$EXP/data/csv
export PYTHONUTF8=1
if [ "$STAGE" = run ] || [ "$STAGE" = all ]; then   # ~25 мин, 8 потоков, 91 МБ кадров
    "${GUO_SHOT:?укажите GUO_SHOT=<путь к guo2019_shot>}" --config=1T --level=scaled --v0=76 \
        --save-t=2.5e-6 --threads=8 --out="$EXP/data/run" > "$EXP/data/run.log"
fi
if [ "$STAGE" != anim ] && [ "$STAGE" != scene ]; then
    "$PY" "$REPO/tools/paraview/sph2vtk.py" --run-dir "$RUN" --out-dir "$VTK" --jobs 6
    mkdir -p "$CSV"; cp "$VTK"/{frames,breaks}.csv "$VTK/columns.txt" "$RUN/monitor.csv" "$CSV/"
fi
if [ "$STAGE" != convert ] && [ "$STAGE" != run ]; then
    "$PV" "$REPO/tools/paraview/fiber_scene.py" --vtk-dir "$VTK" --zoom 0.02 --state "$VTK/scene.pvsm"
    "$PV" "$REPO/tools/paraview/fiber_scene.py" --vtk-dir "$VTK" --zoom 0.02 --snapshots 0,21,33,-1 \
        --out "$VTK/snap.png" --size 1600x800
fi
if [ "$STAGE" = anim ] || [ "$STAGE" = all ]; then
    "$PV" "$REPO/tools/paraview/fiber_scene.py" --vtk-dir "$VTK" --zoom 0.02 --out "$VTK/anim.mp4" --fps 15
fi
