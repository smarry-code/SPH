#!/usr/bin/env bash
# tests/ceramic_shot_cli_test.sh <путь к bc2_shot|bc6_shot> <кейс: bc2|bc6>
# Проверка экспериментальных утилит керамики (cases/common/ceramic_shot.hpp) на грубой сетке
# (Δx = 0,8 мм, t_total = 0,3 мкс, импульсный контакт):
#   * коды выхода: справка 0, ошибка ключей 2, --restart без каталога 1, повторный запуск в
#     занятый каталог 1, --overwrite 0;
#   * выходные файлы прогона и сводка runs.csv; outcome=completed в shot.txt;
#   * серия --contact-list=0,3 -- два каталога;
#   * рестарт: «--checkpoint-every=8 --max-steps=20» (outcome=stopped), затем «--restart» --
#     кривая curve.txt и финальный кадр совпадают с непрерывным прогоном побайтно.
# Физика не проверяется: это делают test_e2e_bc2 / test_e2e_bc6.
set -u
SHOT="${1:?нужен путь к утилите}"
CASE="${2:?нужен кейс bc2|bc6}"
export OMP_NUM_THREADS=1
T="$(mktemp -d "${TMPDIR:-/tmp}/sph_${CASE}_shot_cli.XXXXXX")"
trap 'rm -rf "$T"' EXIT
fail=0
ok()  { echo "  ok   : $1"; }
bad() { echo "  FAIL : $1"; fail=1; }

A=(--delta-x=8e-4 --t-total=3e-7 --contact-type=3 --threads=1 --quiet)
RUN="${CASE}_dx800_c3"

expect_code() { # ожидаемый_код описание -- команда...
    local want="$1" what="$2"; shift 3
    "$@" >"$T/out.log" 2>&1
    local got=$?
    if [ "$got" -eq "$want" ]; then ok "$what (код $got)"; else bad "$what: код $got, ожидался $want"; tail -5 "$T/out.log"; fi
}
data() { grep -v '^#' "$1"; }   # строки кривой без заголовка
kv() { sed -n "s/^$2=//p" "$1/shot.txt"; }

expect_code 0 "--help" -- "$SHOT" --help
expect_code 2 "неизвестный ключ" -- "$SHOT" --no-such-key
expect_code 2 "--contact-type=4 вне 0..3" -- "$SHOT" --contact-type=4
expect_code 2 "--force без --restart" -- "$SHOT" "${A[@]}" --out="$T/e" --force
expect_code 2 "--restart вместе с --overwrite" -- "$SHOT" "${A[@]}" --out="$T/e" --restart --overwrite
expect_code 1 "--restart без каталога прогона" -- "$SHOT" "${A[@]}" --out="$T/e" --restart

# Непрерывный прогон и его файлы.
expect_code 0 "непрерывный прогон" -- "$SHOT" "${A[@]}" --out="$T/c"
for f in log.txt curve.txt shot.txt compare.svg initial_conditions.used.txt; do
    [ -s "$T/c/$RUN/$f" ] && ok "файл $f записан" || bad "файла $f нет"
done
ls "$T/c/$RUN"/sph_output_*.txt >/dev/null 2>&1 && ok "кадры sph_output_*.txt записаны" || bad "кадров нет"
[ "$(kv "$T/c/$RUN" outcome)" = completed ] && ok "outcome=completed" || bad "outcome=$(kv "$T/c/$RUN" outcome)"
[ "$(kv "$T/c/$RUN" contact_force_type)" = 3 ] && ok "contact_force_type=3 в shot.txt" || bad "contact_force_type в shot.txt не 3"
steps="$(kv "$T/c/$RUN" steps)"
[ "${steps:-0}" -gt 20 ] && ok "шагов $steps > 20" || bad "шагов ${steps:-0}, нужно > 20 для проверки рестарта"
[ "$(wc -l < "$T/c/runs.csv")" -eq 2 ] && ok "runs.csv: заголовок + строка" || bad "runs.csv: $(wc -l < "$T/c/runs.csv") строк"

expect_code 1 "повторный запуск в занятый каталог" -- "$SHOT" "${A[@]}" --out="$T/c"
expect_code 0 "--overwrite" -- "$SHOT" "${A[@]}" --out="$T/c" --overwrite
if [ "$(data "$T/c/$RUN/curve.txt" | wc -l)" = "$steps" ]; then ok "после --overwrite кривая не дописана к старой"; else bad "после --overwrite в кривой не $steps точек"; fi

# Серия по типам контакта.
expect_code 0 "серия --contact-list=0,3" -- "$SHOT" "${A[@]}" --contact-list=0,3 --out="$T/l"
[ -f "$T/l/${CASE}_dx800_c0/shot.txt" ] && [ -f "$T/l/${CASE}_dx800_c3/shot.txt" ] && ok "серия: два каталога" || bad "серия: нет каталогов c0/c3"

# Раздельный прогон с рестартом против непрерывного.
expect_code 0 "часть 1: --max-steps=20" -- "$SHOT" "${A[@]}" --out="$T/s" --checkpoint-every=8 --max-steps=20
[ "$(kv "$T/s/$RUN" outcome)" = stopped ] && ok "часть 1: outcome=stopped" || bad "часть 1: outcome=$(kv "$T/s/$RUN" outcome)"
[ -f "$T/s/$RUN/checkpoint.ckpt" ] && ok "checkpoint.ckpt записан" || bad "checkpoint.ckpt нет"
[ -f "$T/s/$RUN/scene_id.txt" ] && ok "scene_id.txt записан" || bad "scene_id.txt нет"
expect_code 0 "часть 2: --restart" -- "$SHOT" "${A[@]}" --out="$T/s" --checkpoint-every=8 --restart
if cmp -s <(data "$T/c/$RUN/curve.txt") <(data "$T/s/$RUN/curve.txt"); then ok "curve.txt совпадает с непрерывным"; else bad "curve.txt отличается от непрерывного"; fi
last="$(ls "$T/c/$RUN"/sph_output_*.txt | sort | tail -n 1)"
if cmp -s "$last" "$T/s/$RUN/$(basename "$last")"; then ok "финальный кадр $(basename "$last") совпадает"; else bad "финальный кадр $(basename "$last") отличается"; fi
[ "$(wc -l < "$T/s/runs.csv")" -eq 3 ] && ok "runs.csv: строки обеих частей" || bad "runs.csv: $(wc -l < "$T/s/runs.csv") строк"

[ "$fail" -eq 0 ] && echo "PASS" || echo "FAIL"
exit "$fail"
