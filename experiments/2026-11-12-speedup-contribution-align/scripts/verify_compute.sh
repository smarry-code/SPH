#!/usr/bin/env bash
# Проверка пункта 2 аудита (alignas(64) у ContributionBuffer) на claude-compute: сборки, тесты, побитово против f1298aa, замер 2 потока база/новая попеременно.
set -u
R=$HOME/c3fs_tests; mkdir -p $R
log() { echo "$(date "+%F %T") $*" | tee -a "$R/verify.log"; }
NEW=$HOME/SPH_c3opt2; REF=$HOME/SPH_c3opt
# база: бинарь 0b9d3b0 (A–D) до правки
[ -x $R/base_shot ] || cp $NEW/build/guo2019_shot $R/base_shot
cp $R/contributions.hpp $NEW/include/rhs/contributions.hpp && echo "0b9d3b0 + alignas(64) ContributionBuffer (аудит п.2), незакоммичено" > $NEW/SOURCE_COMMIT
log "=== старт: NEW=$(cat $NEW/SOURCE_COMMIT); REF=$(cat $REF/SOURCE_COMMIT)"
# 1. сборки
for b in build build-elog; do
  ex=""; [ $b = build-elog ] && ex="-DSPH_ENERGY_LOG=ON"
  cmake -S $NEW -B $NEW/$b -DCMAKE_BUILD_TYPE=Release $ex > $R/cmake_$b.log 2>&1
  cmake --build $NEW/$b -j2 -- -k > $R/build_$b.log 2>&1; log "сборка NEW $b rc=$? предупреждений: $(grep -c "warning:" $R/build_$b.log)"
  grep "warning:" $R/build_$b.log | sed "s#.*/SPH_c3opt2/##" | sort -u > $R/warn_$b.txt
done
cmake -S $REF -B $REF/build-elog -DCMAKE_BUILD_TYPE=Release -DSPH_ENERGY_LOG=ON > $R/cmake_ref.log 2>&1
cmake --build $REF/build-elog -j2 --target guo2019_shot > $R/build_ref.log 2>&1; log "сборка REF build-elog (guo2019_shot) rc=$? предупреждений: $(grep -c "warning:" $R/build_ref.log)"
log "предупреждения в изменённых файлах NEW:"
grep -hE "spatial_index|cpu_executor|contributions.hpp|fiber_contact.hpp|fiber_solid_contact.hpp|fiber_contact_law|fiber_friction|fiber_normal_projection|ode_solver|single_step|solid_sph|diagnostics.hpp|time_integration|test_spatial_index|test_executor|test_fiber_contact|test_fiber_friction|test_fiber_solid|test_fiber_contact_projection" $R/warn_*.txt | sort -u | tee -a $R/verify.log
# 2. изменённые юнит-тесты -- с выводом
for t in test_spatial_index test_executor test_fiber_contact test_fiber_friction test_fiber_solid test_fiber_contact_projection test_fiber_time_integration_equivalence; do
  for b in build build-elog; do
    if [ ! -x $NEW/$b/$t ]; then log "$t ($b): нет бинаря"; continue; fi
    $NEW/$b/$t > $R/${t}_$b.log 2>&1; log "$t ($b) rc=$?: $(tail -1 $R/${t}_$b.log)"
  done
  grep -E "FAIL" $R/${t}_build.log 2>/dev/null | head -20 | tee -a $R/verify.log
done
# 3. ctest целиком, по одному тесту
for b in build build-elog; do
  ctest --test-dir $NEW/$b -j1 --output-on-failure > $R/ctest_$b.log 2>&1
  log "ctest $b rc=$?: $(grep -E "tests passed" $R/ctest_$b.log)"; grep -E "\*\*\*|Failed  |Timeout" $R/ctest_$b.log | tee -a $R/verify.log
done
( cd $R && mkdir -p ckpt_test && cd ckpt_test && $NEW/build/test_checkpoint > $R/ckpt_test.log 2>&1 ); log "test_checkpoint (вывод) rc=$?"
grep -E "^PASS|^FAIL" $R/ckpt_test.log | tee -a $R/verify.log
# 4. побитово против REF, один поток: законы 0/3/4 (проекция выключена)
export OMP_NUM_THREADS=1
SKIP="^(wall|cost|monitor_share_pct|scene_build_s|frame_|monitor_wall)"
shotcmp() { diff <(tr " " "\n" < $1 | grep -vE "$SKIP") <(tr " " "\n" < $2 | grep -vE "$SKIP") > $3 && echo СОВПАДАЕТ || echo ОТЛИЧАЕТСЯ; }
P0=(--config=1T --level=smoke --v0=79 --mu=0.27 --dt-scale=0.5 --threads=1 --t-max=1e-4 --energy-log=100)
for L in 0 3 4; do
  for s in ref new; do
    bin=$REF/build-elog/guo2019_shot; [ $s = new ] && bin=$NEW/build-elog/guo2019_shot
    rm -rf $R/cmp/law$L/$s; $bin "${P0[@]}" --ff-contact=$L --fs-contact=$L --out=$R/cmp/law$L/$s > $R/cmp_law${L}_$s.log 2>&1
  done
  A=$(ls -d $R/cmp/law$L/ref/1T_smoke_*); B=$(ls -d $R/cmp/law$L/new/1T_smoke_*)
  for f in monitor.csv energy_log.csv; do cmp -s $A/$f $B/$f && r=СОВПАДАЕТ || r=ОТЛИЧАЕТСЯ; log "  закон $L $f: $r"; done
  log "  закон $L shot.txt (без времени): $(shotcmp $A/shot.txt $B/shot.txt $R/cmp_law$L.shot.diff)"
done
# 5. проекция: NEW против REF -- теперь тоже должно совпасть побитово
PP=(--config=1T --level=smoke --v0=79 --mu=0.27 --dt-scale=0.5 --energy-log=100 --max-steps=2000 --no-bracket
    --ff-contact=3 --fs-contact=3 --ff-contact-predictor=1 --fs-contact-predictor=1)
rm -rf $R/proj
$REF/build-elog/guo2019_shot "${PP[@]}" --threads=1 --out=$R/proj/ref > $R/proj_ref.log 2>&1; log "проекция REF rc=$?"
$NEW/build-elog/guo2019_shot "${PP[@]}" --threads=1 --out=$R/proj/new > $R/proj_new.log 2>&1; log "проекция NEW rc=$?"
OMP_NUM_THREADS=2 $NEW/build-elog/guo2019_shot "${PP[@]}" --threads=2 --out=$R/proj/new2 > $R/proj_new2.log 2>&1; log "проекция NEW 2 потока rc=$?"
A=$(ls -d $R/proj/ref/1T_smoke_*); B=$(ls -d $R/proj/new/1T_smoke_*); B2=$(ls -d $R/proj/new2/1T_smoke_*)
for f in monitor.csv energy_log.csv; do
  cmp -s $A/$f $B/$f && r=СОВПАДАЕТ || r="ОТЛИЧАЕТСЯ ($(cmp $A/$f $B/$f | grep -o "line [0-9]*"))"; log "  проекция NEW/REF $f: $r"
  cmp -s $B/$f $B2/$f && r=СОВПАДАЕТ || r="ОТЛИЧАЕТСЯ ($(cmp $B/$f $B2/$f | grep -o "line [0-9]*"))"; log "  проекция NEW 1/2 потока $f: $r"
done
log "  проекция NEW/REF shot.txt (без времени): $(shotcmp $A/shot.txt $B/shot.txt $R/proj.shot.diff)"
log "  NEW shot.txt проекция: $(tr " " "\n" < $B/shot.txt | grep -E "^(proj_slots_max|proj_warm_scaled_steps|W_proj_pos_max)" | tr "\n" " ")"
# 6. рестарт NEW с чекпоинта против непрерывного
C=$R/ckpt; rm -rf $C
$NEW/build-elog/guo2019_shot "${PP[@]}" --threads=1 --out=$C/s --checkpoint-every=1200 --max-steps=1300 > $C.s1.log 2>&1; log "чекпоинт ч.1 rc=$?"
$NEW/build-elog/guo2019_shot "${PP[@]}" --threads=1 --out=$C/s --restart > $C.s2.log 2>&1; log "чекпоинт ч.2 rc=$?: $(grep "монитор:" $C.s2.log)"
S=$(ls -d $C/s/1T_smoke_*)
for f in monitor.csv energy_log.csv; do cmp -s $B/$f $S/$f && r=СОВПАДАЕТ || r=ОТЛИЧАЕТСЯ; log "  рестарт $f: $r"; done
SKIP="^(wall|cost|monitor_share_pct|scene_build_s|frame_|monitor_state|monitor_wall)"
log "  рестарт shot.txt: $(shotcmp $B/shot.txt $S/shot.txt $R/ckpt.shot.diff)"
# 7. полная сцена Guo 1T на 2 потоках: база (0b9d3b0) и новая попеременно, по 2 повтора; monitor.csv против ~/c3prof/res (f1298aa)
FIB=$HOME/c3prof/twaron_s25_dmg0.00415.fib
COMMON="--v0=79 --config=1T --level=base --mu=0.27 --frames=off --save-t=1e-4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --fiber-inelastic=$FIB --eps-break=0.05 --fiber-E=8e10 --dt-scale=0.5 --impulse-push-tau=0"
VH="--ff-contact=0 --fs-contact=0"
VP="--ff-contact=3 --fs-contact=3 --ff-contact-predictor=1 --fs-contact-predictor=1 --contact-proj-iters=8 --contact-proj-warm=1"
for rep in 1 2; do
 for v in H P; do
  for s in base new; do
    bin=$R/base_shot; [ $s = new ] && bin=$NEW/build/guo2019_shot
    d=$R/full/$v.$s.$rep; rm -rf "$d"; mkdir -p "$d"
    V=$VH; [ $v = P ] && V=$VP
    OMP_NUM_THREADS=2 OMP_PROC_BIND=close taskset -c 0-1 $bin $COMMON $V --threads=2 --t-max=4e-5 --out=$d/out > $d/stdout.log 2>&1; rc=$?
    log "полная сцена $v.$s.$rep rc=$rc: $(grep -o 'шагов=[0-9]*, wall=[0-9.]* с ([0-9.]* с/шаг)' $d/stdout.log)"
    M=$d/out/1T_base_v79.0_r0/monitor.csv; MR=$HOME/c3prof/res/$v.t2/out/1T_base_v79.0_r0/monitor.csv
    cmp -s "$M" "$MR" && r=СОВПАДАЕТ || r="ОТЛИЧАЕТСЯ ($(cmp "$M" "$MR" 2>&1 | grep -o "line [0-9]*"))"
    log "  $v.$s.$rep monitor.csv против f1298aa: $r"
  done
 done
done
log "готово"
