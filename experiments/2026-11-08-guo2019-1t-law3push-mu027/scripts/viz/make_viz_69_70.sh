#!/usr/bin/env bash
# make_viz_69_70.sh: анимации выстрелов 69 и 70 м/с (окно [69; 70], 1T base, закон 3 + выталкивание tau 0,5 мкс,
# шаг x0,5, повреждённость 4 %, E = 80 ГПа, mu = 0,27, вердикт по шару). По образцу
# experiments/2026-10-26-guo2019-damage-summary/scripts/viz/make_viz_batch6.sh (скрипты рисования -- оттуда).
# Каталоги _long -- полная история 0..~2 мс (кадры базового прогона + продолжение из чекпоинта).
#
# Вход ($DATA), любой из двух видов на каждую скорость:
#   $DATA/law3push_mu0.27_v<V>_long/<прогон>/  -- скачанные логи и каждый второй кадр (так строилось 2026-11-08);
#   $DATA/mu0.27_v<V>_long/<прогон>/           -- рабочий каталог драйвера как есть (scripts/rebuild_artifacts.sh).
# Выход: $OUT (по умолчанию <рабочий каталог>/anim; в каталог эксперимента не пишется):
#   law3push_mu0.27_v<V>_long/{anim2_*.gif, anim2_*_frame*.png, vz_anim_*.gif, energy_*.png},
#   summary_1T_law3push_mu027_v69_70.png
# Запуск: PY=<python с matplotlib> DATA=<данные> OUT=<куда> bash make_viz_69_70.sh [69 70]
set -u
export PYTHONIOENCODING=utf-8 PYTHONUTF8=1
PY=${PY:-python3}
DATA=${DATA:?каталог со скачанными выстрелами}
OUT=${OUT:?каталог для анимаций (вне репозитория)}
HERE=$(cd "$(dirname "$0")" && pwd)
EXP=$(cd "$HERE/../.." && pwd)
SC=$(cd "$EXP/../2026-10-26-guo2019-damage-summary/scripts/viz" && pwd)
win() { if command -v cygpath > /dev/null; then cygpath -m "$1"; else echo "$1"; fi; }   # Windows-путь для python
PU="закон 3 + выталкивание τ = 0,5 мкс, шаг ×0,5"
names=(${@:-69 70})
sumargs=()
mkdir -p "$OUT"
for vv in "${names[@]}"; do
  name=law3push_mu0.27_v${vv}_long; dir=$DATA/$name
  [ -d "$dir" ] || dir=$DATA/mu0.27_v${vv}_long
  R=$(win "$(dirname "$(ls $dir/*/shot.txt | head -1)")")
  sumargs+=("$vv=$(win $dir)")
  D=$OUT/$name; mkdir -p "$D"; Dm=$(win "$D")
  full="1T base, $PU, повреждённость (разрыв при 4 %), E = 80 ГПа, μ = 0.27, V0 = $vv м/с (до 3 мс)"
  echo "== $name"
  "$PY" "$SC/anim_vz.py" "$R" "$Dm/vz_anim_$name.gif" --law "$PU, повреждённость 4 %, E = 80 ГПа"
  "$PY" "$SC/plot_energy_all.py" "$R/energy_log.csv" "$Dm" --title "$full"
  "$PY" "$SC/ball_defl.py" "$name=$R" 2>&1 | tail -3
  "$PY" "$SC/anim_2side.py" --run-dir "$R" --out "$Dm/anim2_$name.gif" --max-frames 40 --title "$full" || echo "FAIL $name"
  echo "ALLDONE $name"
done
SUMMARY_TITLE="1T base, $PU, повреждённость 4 %, E = 80 ГПа, μ = 0,27: окно [69; 70]" SUMMARY_TMAX=2200 \
  "$PY" "$SC/summary_dt05.py" "$(win $OUT)/summary_1T_law3push_mu027_v69_70.png" "${sumargs[@]}"
echo ALLDONE_ALL
