#!/usr/bin/env bash
# Побитовая регрессия Guo против эталона scratchpad/baseline. Аргумент: метка прогона.
set -u
S=/c/Users/Maria/AppData/Local/Temp/claude/c--dip-sph-SPH/9c93d159-f8d5-4f88-9269-5c9cdc59e2d2/scratchpad
W=/c/dip/sph/SPH/.claude/worktrees/fiber-contact-types
TAG=${1:-new}
OUT=$S/regress_$TAG
rm -rf "$OUT"; mkdir -p "$OUT"
cd "$W" && cmake --build build --target guo2019_shot test_e2e_guo2019 > "$OUT/build.log" 2>&1 || { echo "BUILD FAILED"; exit 1; }
export OMP_NUM_THREADS=1
# Файлы с wall-временем исключаются (см. отчёт базовой линии).
filt() { grep -v -E '/(shot\.txt|shots\.csv|run\.log|rc\.txt|workdir_src\.txt|sha256\.txt)$'; }
status=0
for cfg in A 1T_A; do
  d=$OUT/e2e_$cfg; mkdir -p "$d"; cd "$d"
  "$W/build/test_e2e_guo2019.exe" --config=$cfg --keep-workdir --time-limit=600 > run.log 2>&1
  echo "e2e_$cfg rc=$?"
  wd=$(grep -o -E '[A-Za-z]:[\\/][^ ]*sph_e2e_guo2019_[A-Za-z0-9_]*' run.log | head -1)
  wdu=$(cygpath -u "$wd" 2>/dev/null || echo "$wd")
  cp -r "$wdu" "$d/workdir"
  (cd "$d/workdir" && find . -type f | sort | filt | xargs sha256sum) > "$d/new.sha"
  (cd "$S/baseline/e2e_$cfg/workdir" && find . -type f | sort | filt | xargs sha256sum) > "$d/base.sha"
  if diff -q "$d/base.sha" "$d/new.sha" > /dev/null; then echo "e2e_$cfg: BITWISE OK ($(wc -l < "$d/new.sha") files)"; else echo "e2e_$cfg: DIFF"; diff "$d/base.sha" "$d/new.sha" | head; status=1; fi
done
d=$OUT/shot; mkdir -p "$d"; cd "$d"
"$W/build/guo2019_shot.exe" --config=1T_A --level=smoke --v0=150 --t-max=2e-5 --frames=off --threads=1 --out="$d/out" > run.log 2>&1
echo "shot rc=$?"
(cd "$d/out" && find . -type f | sort | filt | xargs sha256sum) > "$d/new.sha"
(cd "$S/baseline/shot/out" && find . -type f | sort | filt | xargs sha256sum) > "$d/base.sha"
if diff -q "$d/base.sha" "$d/new.sha" > /dev/null; then echo "shot: BITWISE OK ($(wc -l < "$d/new.sha") files)"; else echo "shot: DIFF"; diff "$d/base.sha" "$d/new.sha" | head; status=1; fi
exit $status
