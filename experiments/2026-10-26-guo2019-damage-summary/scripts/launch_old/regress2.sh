#!/usr/bin/env bash
# Побитовое сравнение двух сборок: regress2.sh <bin_base> <bin_new> <out>
set -u
B=$1; N=$2; OUT=$3
rm -rf "$OUT"; mkdir -p "$OUT"
export OMP_NUM_THREADS=1
filt() { grep -v -E '/(shot\.txt|shots\.csv|run\.log)$'; }
run_set() { # $1 bin, $2 dir
  local bin=$1 d=$2; mkdir -p "$d"
  for cfg in A 1T_A; do
    mkdir -p "$d/e2e_$cfg"; ( cd "$d/e2e_$cfg" && "$bin/test_e2e_guo2019.exe" --config=$cfg --keep-workdir --time-limit=600 > run.log 2>&1; echo "e2e_$cfg rc=$?" )
    wd=$(grep -o -E '[A-Za-z]:[\\/][^ ]*sph_e2e_guo2019_[A-Za-z0-9_]*' "$d/e2e_$cfg/run.log" | head -1)
    cp -r "$(cygpath -u "$wd")" "$d/e2e_$cfg/workdir"
  done
  ( cd "$d" && "$bin/guo2019_shot.exe" --config=1T_A --level=smoke --v0=150 --t-max=2e-5 --frames=off --threads=1 --out="$d/shot" > shot_run.log 2>&1; echo "shot rc=$?" )
}
run_set "$B" "$OUT/base"
run_set "$N" "$OUT/new"
st=0
for sub in e2e_A/workdir e2e_1T_A/workdir shot; do
  (cd "$OUT/base/$sub" && find . -type f | sort | filt | xargs sha256sum) > "$OUT/b.sha"
  (cd "$OUT/new/$sub" && find . -type f | sort | filt | xargs sha256sum) > "$OUT/n.sha"
  if diff -q "$OUT/b.sha" "$OUT/n.sha" > /dev/null; then echo "$sub: BITWISE OK ($(wc -l < "$OUT/n.sha") files)"; else echo "$sub: DIFF"; diff "$OUT/b.sha" "$OUT/n.sha" | head; st=1; fi
done
exit $st
