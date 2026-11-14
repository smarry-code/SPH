#!/usr/bin/env bash
# Цепочка этапов в одной сессии claude-compute (2 потока, полная сцена Guo 1T, t_max 4e-5, 4372 шага):
# 0b9d3b0 (base_shot) -> alignas(64) (bin_align, 5995594) -> H3 (bin_h3, 10fbcb1) -> H1+H4+правки Codex (bin_fix, c8973f7).
# Бинарники попеременно, порядок внутри повтора сдвигается (снять дрейф ВМ); monitor.csv каждого прогона против f1298aa.
set -u
T=$HOME/c3sp_tests; R=$T/chain; mkdir -p $R
log() { echo "$(date '+%F %T') $*" | tee -a $R/chain.log; }
FIB=$HOME/c3prof/twaron_s25_dmg0.00415.fib
COMMON="--v0=79 --config=1T --level=base --mu=0.27 --frames=off --save-t=1e-4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --fiber-inelastic=$FIB --eps-break=0.05 --fiber-E=8e10 --dt-scale=0.5 --impulse-push-tau=0"
VH="--ff-contact=0 --fs-contact=0"
VP="--ff-contact=3 --fs-contact=3 --ff-contact-predictor=1 --fs-contact-predictor=1 --contact-proj-iters=8 --contact-proj-warm=1"
declare -A BIN=([base]=$HOME/c3fs_tests/base_shot [align]=$T/bin_align [h3]=$T/bin_h3 [fix]=$T/bin_fix)
log "=== старт цепочки; $(for k in base align h3 fix; do printf '%s=%s ' $k $(md5sum ${BIN[$k]} | cut -c1-8); done)"
ORD1="base align h3 fix"; ORD2="fix h3 align base"
for rep in 1 2; do
 ORD=$ORD1; [ $rep = 2 ] && ORD=$ORD2
 for v in H P; do
  for s in $ORD; do
    d=$R/$v.$s.$rep; rm -rf "$d"; mkdir -p "$d"
    V=$VH; [ $v = P ] && V=$VP
    OMP_NUM_THREADS=2 OMP_PROC_BIND=close taskset -c 0-1 ${BIN[$s]} $COMMON $V --threads=2 --t-max=4e-5 --out=$d/out > $d/stdout.log 2>&1; rc=$?
    M=$d/out/1T_base_v79.0_r0/monitor.csv; MR=$HOME/c3prof/res/$v.t2/out/1T_base_v79.0_r0/monitor.csv
    cmp -s "$M" "$MR" && r=СОВПАДАЕТ || r=ОТЛИЧАЕТСЯ
    log "$v.$s.$rep rc=$rc: $(grep -o 'шагов=[0-9]*, wall=[0-9.]* с ([0-9.]* с/шаг)' $d/stdout.log); monitor против f1298aa: $r"
  done
 done
done
log "готово"
