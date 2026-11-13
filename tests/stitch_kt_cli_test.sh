#!/usr/bin/env bash
# tests/stitch_kt_cli_test.sh <путь к guo2019_shot>
# Стёжка в пакете kT (--config=3T, Config::TK) через guo2019_shot:
#   1) --stitch=grid --stitch-tile=4 на smoke (12 x 12 нитей, кусочки 3 x 3 нити): n_stitch в shot.txt =
#      (3 + 3 линии) x 12 проколов x (3 - 1) пары слоёв = 144; столбцы стежков в monitor.csv -- последние;
#   2) раздельный прогон (--checkpoint-every=4 --max-steps=10, затем --restart --max-steps=15) против непрерывного
#      (--max-steps=15): финальный кадр нитей побайтно и последняя строка monitor.csv (столбцы стежков, с n_stitch
#      и slip_mean, которые после рестарта берутся из сцены, а не из чекпоинта) совпадают;
#   3) --restart с другой стёжкой (--stitch-tile=6) отклоняется по scene_id.txt;
#   4) --stitch=grid при 1 слое -- ошибка разбора (код 2).
set -u
SHOT="${1:?нужен путь к guo2019_shot}"
export OMP_NUM_THREADS=1
T="$(mktemp -d "${TMPDIR:-/tmp}/sph_stitch_kt.XXXXXX")"
trap 'rm -rf "$T"' EXIT
fail=0
ok()  { echo "  ok   : $1"; }
bad() { echo "  FAIL : $1"; fail=1; }

A=(--config=3T --level=smoke --v0=100 --frames=off --threads=1 --no-bracket --stitch=grid --stitch-tile=4)

"$SHOT" --config=1T --level=smoke --v0=100 --frames=off --threads=1 --no-bracket --stitch=grid --out="$T/e1" \
    --max-steps=2 >"$T/e1.log" 2>&1
[ $? -eq 2 ] && ok "--stitch=grid при 1T -- код 2" || bad "--stitch=grid при 1T: ожидался код 2"

"$SHOT" "${A[@]}" --out="$T/c" --max-steps=15 >"$T/c.log" 2>&1 || bad "непрерывный прогон упал"
"$SHOT" "${A[@]}" --out="$T/s" --checkpoint-every=4 --max-steps=10 >"$T/s1.log" 2>&1 || bad "часть 1 упала"
"$SHOT" "${A[@]}" --out="$T/s" --restart --max-steps=15 >"$T/s2.log" 2>&1 || bad "часть 2 (рестарт) упала"

CD="$(ls -d "$T"/c/3T_smoke_* 2>/dev/null | head -1)"; SD="$(ls -d "$T"/s/3T_smoke_* 2>/dev/null | head -1)"
if [ -z "$CD" ] || [ -z "$SD" ]; then bad "нет каталогов выстрелов 3T"; echo "FAIL"; exit 1; fi

ns="$(tr ' ' '\n' < "$CD/shot.txt" | sed -n 's/^n_stitch=//p')"
[ "$ns" = "144" ] && ok "n_stitch = 144" || bad "n_stitch = '$ns', ожидалось 144"
hdr="$(head -1 "$CD/monitor.csv")"
case "$hdr" in
    *",n_stitch,n_stitch_broken,F_stitch_max,F_stitch_sum,slip_mean,gap_mean") ok "столбцы стежков -- последние" ;;
    *) bad "столбцы стежков не в конце заголовка monitor.csv" ;;
esac
grep -q "stitch=grid tile=4" "$SD/scene_id.txt" && ok "стёжка в scene_id.txt" || bad "нет стёжки в scene_id.txt"

if cmp -s "$CD/fiber_output_99999.txt" "$SD/fiber_output_99999.txt"; then ok "финальный кадр нитей совпадает"
else bad "финальный кадр нитей отличается"; fi
lastc="$(tail -1 "$CD/monitor.csv" | awk -F, '{print $(NF-5)","$(NF-4)","$(NF-3)","$(NF-2)","$(NF-1)","$NF}')"
lasts="$(tail -1 "$SD/monitor.csv" | awk -F, '{print $(NF-5)","$(NF-4)","$(NF-3)","$(NF-2)","$(NF-1)","$NF}')"
if [ -n "$lastc" ] && [ "$lastc" = "$lasts" ]; then ok "столбцы стежков последней строки совпадают ($lastc)"
else bad "столбцы стежков: непрерывный '$lastc', с рестартом '$lasts'"; fi

B=("${A[@]/--stitch-tile=4/--stitch-tile=6}")
"$SHOT" "${B[@]}" --out="$T/s" --restart --max-steps=20 >"$T/s3.log" 2>&1
if grep -q "сцена не совпадает" "$T/s3.log"; then ok "рестарт с другой стёжкой отклонён"
else bad "рестарт с другой стёжкой не отклонён"; fi

[ "$fail" -eq 0 ] && echo "PASS" || echo "FAIL"
exit "$fail"
