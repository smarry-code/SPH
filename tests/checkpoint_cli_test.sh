#!/usr/bin/env bash
# tests/checkpoint_cli_test.sh <путь к guo2019_shot>
# CLI чекпоинтов/рестарта guo2019_shot: коды ошибок и сравнение раздельного прогона
# («--checkpoint-every=4 --max-steps=10», затем «--restart --max-steps=15») с непрерывным
# («--max-steps=15») по детерминированному результату: финальные кадры sph_output_99999.txt и
# fiber_output_99999.txt (побайтно) и t_end в shots.csv. Остальные столбцы shots.csv содержат
# время стенда (wall_s, ...) и счётчики монитора, который в чекпоинт не входит, -- не сравниваются.
set -u
SHOT="${1:?нужен путь к guo2019_shot}"
export OMP_NUM_THREADS=1
T="$(mktemp -d "${TMPDIR:-/tmp}/sph_ckpt_cli.XXXXXX")"
trap 'rm -rf "$T"' EXIT
fail=0
ok()  { echo "  ok   : $1"; }
bad() { echo "  FAIL : $1"; fail=1; }

A=(--config=1T --level=smoke --v0=100 --frames=off --threads=1 --no-bracket)

expect_code() { # ожидаемый_код описание -- команда...
    local want="$1" what="$2"; shift 3
    "$@" >"$T/out.log" 2>&1
    local got=$?
    if [ "$got" -eq "$want" ]; then ok "$what (код $got)"; else bad "$what: код $got, ожидался $want"; fi
}

expect_code 2 "--force без --restart" -- "$SHOT" "${A[@]}" --out="$T/e1" --max-steps=3 --force
expect_code 2 "--checkpoint-every=0 -- ошибка разбора" -- "$SHOT" "${A[@]}" --out="$T/e2" --checkpoint-every=0
# --restart без чекпоинта: ошибка выстрела (ОШИБКА + «нет чекпоинта»). Код возврата утилиты при
# исключении в выстреле 0 (так устроен guo2019_shot: серия продолжается, итог -- "выстрелов с
# исключением"), поэтому проверяется текст и отсутствие чекпоинта/кадров, а не код.
"$SHOT" "${A[@]}" --out="$T/e3" --restart --max-steps=3 >"$T/out.log" 2>&1
if grep -q "ОШИБКА" "$T/out.log" && grep -q "нет чекпоинта" "$T/out.log" && grep -q "выстрелов с исключением: 1" "$T/out.log"
then ok "--restart без файла -- ошибка выстрела"; else bad "--restart без файла: нет сообщения об ошибке"; fi

"$SHOT" "${A[@]}" --out="$T/c" --max-steps=15 >"$T/c.log" 2>&1 || bad "непрерывный прогон упал"
"$SHOT" "${A[@]}" --out="$T/s" --checkpoint-every=4 --max-steps=10 >"$T/s1.log" 2>&1 || bad "часть 1 упала"
"$SHOT" "${A[@]}" --out="$T/s" --restart --max-steps=15 >"$T/s2.log" 2>&1 || bad "часть 2 (рестарт) упала"

CD="$(ls -d "$T"/c/1T_smoke_*)"; SD="$(ls -d "$T"/s/1T_smoke_*)"
[ -f "$SD/checkpoint.ckpt" ] && ok "checkpoint.ckpt записан" || bad "checkpoint.ckpt нет"
[ -f "$SD/checkpoint.ckpt.prev" ] && ok "checkpoint.ckpt.prev есть" || bad "checkpoint.ckpt.prev нет"
for f in sph_output_99999.txt fiber_output_99999.txt; do
    if cmp -s "$CD/$f" "$SD/$f"; then ok "финальный кадр $f совпадает"; else bad "финальный кадр $f отличается"; fi
done
tend() { awk -F, 'NR==2 {print $8}' "$1/shots.csv"; }
if [ "$(tend "$T/c")" = "$(tend "$T/s")" ] && [ -n "$(tend "$T/c")" ]; then ok "t_end в shots.csv совпадает"; else bad "t_end отличается"; fi

[ "$fail" -eq 0 ] && echo "PASS" || echo "FAIL"
exit "$fail"
