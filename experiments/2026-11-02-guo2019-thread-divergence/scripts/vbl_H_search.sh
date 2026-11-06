#!/usr/bin/env bash
# Подбор порога пробития V_BL: 1T, контакт 0 (Герц) нить-нить и нить-шар, mu=0.1, нить twaron_s25_dmg0.00415.
# Раунд = 2 выстрела параллельно по 4 потока. Границы: lo = наибольшая V0 «не пробит», hi = наименьшая V0 «пробит»
# (66 м/с пробит по серии late_H). Пока lo неизвестна, раунд = hi-3, hi-6; затем отрезок [lo; hi] делится на три.
# Останов: hi - lo <= TOL, противоречие (lo > hi), исход undetermined/ошибка или MAXR раундов.
# Решения выводятся только из завершённых выстрелов, поэтому после перезагрузки скрипт повторяет тот же раунд,
# а незавершённые выстрелы продолжаются из чекпоинта (shot_resume.sh). Запуск -- под guarded_run.sh (задача vbl_H).
set -u
D=$HOME/vbl_H
HI0=66
TOL=${TOL:-1}
MAXR=${MAXR:-6}
BIN=$HOME/SPH_det/build-elog/guo2019_shot
KEYS=(--config=1T --level=base --mu=0.1 --frames=off --save-t=1e-4 --threads=4 --n-mon=10 --n-diag=100 --repeat=1
      --no-bracket --t-max=3e-3 --fiber-inelastic=$HOME/late/twaron_s25_dmg0.00415.fib --eps-break=0.05
      --fiber-E=8e10 --fs-contact=0 --ff-contact=0 --energy-log=100 --dt-scale=1 --impulse-push-tau=0
      --checkpoint-every=20000)
LOG=$D/search.log
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

shoot() {
  local v=$1 out=$D/v$1
  ls "$out"/*/shot.txt > /dev/null 2>&1 && return 0
  env OMP_NUM_THREADS=4 SHOT_BIN=$BIN "$HOME/autostart/shot_resume.sh" "$out" --v0="$v" "${KEYS[@]}" > "$D/v$v.log" 2>&1
}

# Сводка завершённых выстрелов: "V0 outcome V_R t_end wall_s"
table() {
  for f in "$D"/v*/shots.csv; do
    [ -e "$f" ] || continue
    tail -n 1 "$f" | awk -F, '{printf "%s %s %s %s %s\n", $3, $5, $6, $8, $10}'
  done | sort -g
}

for r in $(seq 1 "$MAXR"); do
  T=$(table)
  read -r lo hi bad < <(echo "$T" | awk -v hi="$HI0" '
    BEGIN { lo = "nan" }
    $2 == "perforated"     { if ($1 + 0 < hi + 0) hi = $1 }
    $2 == "not_perforated" { if (lo == "nan" || $1 + 0 > lo + 0) lo = $1 }
    NF && $2 != "perforated" && $2 != "not_perforated" { bad = bad $1 ":" $2 ";" }
    END { print lo, hi, (bad == "" ? "-" : bad) }')
  say "раунд $r: lo=$lo hi=$hi"
  if [ "$bad" != "-" ]; then say "останов: неопределённый исход $bad"; break; fi
  if [ "$lo" != nan ] && awk -v a="$lo" -v b="$hi" 'BEGIN { exit !(a > b) }'; then
    say "останов: противоречие, не пробит при $lo > пробит при $hi (полоса неустойчивости исхода)"; break
  fi
  if [ "$lo" != nan ] && awk -v a="$lo" -v b="$hi" -v t="$TOL" 'BEGIN { exit !(b - a <= t + 1e-9) }'; then
    say "готово: V_BL в [$lo; $hi]"; break
  fi
  if [ "$lo" = nan ]; then
    read -r v1 v2 < <(awk -v h="$hi" 'BEGIN { printf "%.1f %.1f\n", h - 3, h - 6 }')
  else
    read -r v1 v2 < <(awk -v a="$lo" -v b="$hi" 'BEGIN { printf "%.1f %.1f\n", a + (b - a) / 3, a + 2 * (b - a) / 3 }')
  fi
  say "раунд $r: выстрелы $v1 и $v2"
  shoot "$v1" & p1=$!
  shoot "$v2" & p2=$!
  wait $p1; rc1=$?; wait $p2; rc2=$?
  if [ $rc1 -ne 0 ] || [ $rc2 -ne 0 ]; then say "останов: выстрел не завершён (rc $rc1 $rc2)"; exit 1; fi
  table | while read -r a b c d e; do say "  V0=$a $b V_R=$c t_end=$d wall=${e}s"; done
done
table > "$D/summary.txt"
exit 0
