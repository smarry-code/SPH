#!/usr/bin/env bash
# tests/mosaic_kt_cli_test.sh <путь к guo2019_shot>
# Мозаичная ткань (--random-cut) в пакете 2T через guo2019_shot (smoke):
#   1) --random-cut=4 --layers-per-cut=2: одна раскройка на оба слоя; в shot.txt n_stitch = n_stitch_in +
#      n_stitch_through, обе части > 0; рядом с выстрелом guo2019_fabric.cut.json; столбцы стежков в monitor.csv --
#      последние;
#   2) раздельный прогон (--checkpoint-every=4 --max-steps=10, затем --restart --max-steps=15) против непрерывного
#      (--max-steps=15): финальный кадр нитей побайтно и столбцы стежков последней строки monitor.csv совпадают;
#   3) --restart с другой раскройкой (--cut-seed=2) отклоняется по scene_id.txt;
#   4) ошибки разбора (код 2): --random-cut вместе с --stitch=grid; --cut-seed без --random-cut; --random-cut при A.
set -u
SHOT="${1:?нужен путь к guo2019_shot}"
export OMP_NUM_THREADS=1
T="$(mktemp -d "${TMPDIR:-/tmp}/sph_mosaic_kt.XXXXXX")"
trap 'rm -rf "$T"' EXIT
fail=0
ok()  { echo "  ok   : $1"; }
bad() { echo "  FAIL : $1"; fail=1; }

C=(--level=smoke --v0=100 --frames=off --threads=1 --no-bracket)
A=(--config=2T "${C[@]}" --random-cut=4 --layers-per-cut=2)

"$SHOT" "${A[@]}" --stitch=grid --out="$T/e1" --max-steps=2 >"$T/e1.log" 2>&1
[ $? -eq 2 ] && ok "--random-cut с --stitch=grid -- код 2" || bad "--random-cut с --stitch=grid: ожидался код 2"
"$SHOT" --config=2T "${C[@]}" --cut-seed=3 --out="$T/e2" --max-steps=2 >"$T/e2.log" 2>&1
[ $? -eq 2 ] && ok "--cut-seed без --random-cut -- код 2" || bad "--cut-seed без --random-cut: ожидался код 2"
"$SHOT" --config=A "${C[@]}" --random-cut=4 --out="$T/e3" --max-steps=2 >"$T/e3.log" 2>&1
[ $? -eq 2 ] && ok "--random-cut при A -- код 2" || bad "--random-cut при A: ожидался код 2"

"$SHOT" "${A[@]}" --out="$T/c" --max-steps=15 >"$T/c.log" 2>&1 || bad "непрерывный прогон упал"
"$SHOT" "${A[@]}" --out="$T/s" --checkpoint-every=4 --max-steps=10 >"$T/s1.log" 2>&1 || bad "часть 1 упала"
"$SHOT" "${A[@]}" --out="$T/s" --restart --max-steps=15 >"$T/s2.log" 2>&1 || bad "часть 2 (рестарт) упала"

CD="$(ls -d "$T"/c/2T_smoke_rc4s1k2_* 2>/dev/null | head -1)"; SD="$(ls -d "$T"/s/2T_smoke_rc4s1k2_* 2>/dev/null | head -1)"
if [ -z "$CD" ] || [ -z "$SD" ]; then bad "нет каталогов выстрелов 2T_smoke_rc4s1k2_*"; echo "FAIL"; exit 1; fi

kv() { tr ' ' '\n' < "$CD/shot.txt" | sed -n "s/^$1=//p"; }
ns="$(kv n_stitch)"; ni="$(kv n_stitch_in)"; nt="$(kv n_stitch_through)"; cc="$(kv cut_count)"
if [ -n "$ns" ] && [ -n "$ni" ] && [ -n "$nt" ] && [ "$ni" -gt 0 ] && [ "$nt" -gt 0 ] && [ "$ns" -eq $((ni + nt)) ]; then
    ok "n_stitch = $ns = $ni в слоях + $nt сквозных"
else bad "стежки: n_stitch='$ns' n_stitch_in='$ni' n_stitch_through='$nt'"; fi
[ "$cc" = "1" ] && ok "одна раскройка на 2 слоя" || bad "cut_count = '$cc', ожидалось 1"
[ -s "$CD/guo2019_fabric.cut.json" ] && ok "guo2019_fabric.cut.json записан" || bad "нет guo2019_fabric.cut.json"
hdr="$(head -1 "$CD/monitor.csv")"
case "$hdr" in
    *",n_stitch,n_stitch_broken,F_stitch_max,F_stitch_sum,slip_mean,gap_mean") ok "столбцы стежков -- последние" ;;
    *) bad "столбцы стежков не в конце заголовка monitor.csv" ;;
esac
grep -q "random_cut=4 seed=1 layers_per_cut=2" "$SD/scene_id.txt" && ok "мозаика в scene_id.txt" || bad "нет мозаики в scene_id.txt"

if cmp -s "$CD/fiber_output_99999.txt" "$SD/fiber_output_99999.txt"; then ok "финальный кадр нитей совпадает"
else bad "финальный кадр нитей отличается"; fi
lastc="$(tail -1 "$CD/monitor.csv" | awk -F, '{print $(NF-5)","$(NF-4)","$(NF-3)","$(NF-2)","$(NF-1)","$NF}')"
lasts="$(tail -1 "$SD/monitor.csv" | awk -F, '{print $(NF-5)","$(NF-4)","$(NF-3)","$(NF-2)","$(NF-1)","$NF}')"
if [ -n "$lastc" ] && [ "$lastc" = "$lasts" ]; then ok "столбцы стежков последней строки совпадают ($lastc)"
else bad "столбцы стежков: непрерывный '$lastc', с рестартом '$lasts'"; fi

# другая раскройка: seed входит в имя каталога, поэтому каталог переименован под seed = 2 -- проверяется
# scene_id.txt, а не имя каталога
mv "$SD" "${SD/rc4s1k2/rc4s2k2}"
"$SHOT" "${A[@]}" --cut-seed=2 --out="$T/s" --restart --max-steps=20 >"$T/s3.log" 2>&1
if grep -q "сцена не совпадает" "$T/s3.log"; then ok "рестарт с другой раскройкой отклонён"
else bad "рестарт с другой раскройкой не отклонён"; fi

[ "$fail" -eq 0 ] && echo "PASS" || echo "FAIL"
exit "$fail"
