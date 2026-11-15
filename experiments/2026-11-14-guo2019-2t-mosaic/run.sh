#!/usr/bin/env bash
# Guo 2019 2T, мозаичная ткань (решение владельца 2026-11-13): случайная раскройка N = 30, одна раскройка на оба слоя
# (k = 2), зажим |y| >= 60 мм для всех бусин (как в генераторе мозаики), стежки через швы и сквозные (шаг p = 1,18 мм,
# материал стёжки «одеяла»: E = 8,5e10 Па, d = 0,2 мм, eps_b = 3 %). Закон 3 как проекция скорости (K = 8, тёплый
# старт), без выталкивания, mu = 0,27. Нить twaron_s25_dmg0.00415.fib, E = 8e10, eps_break = 5 %, dt-scale = 0,5.
# Локально (серверы заняты): 6 потоков, лог энергии, чекпоинт раз в 500 шагов.
# Использование: [LPC=k] run.sh <V0, м/с> [--restart]   (из корня репозитория); LPC -- слоёв на раскройку, по умолчанию 2
# (одна общая раскройка); LPC=1 -- у каждого слоя своя раскройка (одна последовательность ГСЧ: раскройка слоя 1 та же, что при k=2)
set -u
V0=${1:?нужна скорость V0}
shift
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$ROOT/experiments/2026-11-14-guo2019-2t-mosaic/runs"
mkdir -p "$OUT"
export OMP_NUM_THREADS=6
exec "$ROOT/build-elog/guo2019_shot.exe" --config=2T --level=base --v0="$V0" --out="$OUT" \
  --frames=off --save-t=7.5e-6 --threads=6 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=8e-4 \
  --random-cut=30 --cut-seed=1 --layers-per-cut="${LPC:-2}" \
  --mu=0.27 --ff-contact=3 --fs-contact=3 --ff-contact-predictor=1 --fs-contact-predictor=1 --dt-scale=0.5 \
  --fiber-inelastic="$ROOT/experiments/2026-11-09-fabric-kT-ball1g-v550/twaron_s25_dmg0.00415.fib" \
  --fiber-E=8e10 --eps-break=0.05 --energy-log=100 --checkpoint-every=500 "$@"
