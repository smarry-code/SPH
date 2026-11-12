#!/usr/bin/env bash
# Очередь на claude-compute-8t: ждать, пока сервер свободен (threads_free.sh: FREE=8, нет guo2019_shot), затем
# замер правок A–D (отчёт 2026-11-12, коммит 5da639c) и побитовая сверка с f1298aa/1658756.
# Эталон OLD_RES -- каталог ссылок ~/c3opt2_bench/oldres: прогон из ~/c3opt_bench/res (f1298aa), иначе из
# ~/c3fix_bench/res (1658756; результат побитово тот же -- отчёт 2026-11-11, раздел «Внедрение»).
# Регистрация: ~/autostart/register_job.sh c3opt2_bench bash ~/c3opt2_bench/queue.sh
set -u
R=$HOME/c3opt2_bench; mkdir -p "$R/res" "$R/oldres"
log() { echo "$(date "+%F %T") $*" >> "$R/queue.log"; }
RUNS="H.t1 P.t1 H.t1b P.t1b H.t2 P.t2 H.t4 P.t4 H.t8 P.t8 H.prof1 P.prof1 H.prof8 P.prof8"
for n in $RUNS; do
  if [ -d "$HOME/c3opt_bench/res/$n" ]; then ln -sfn "$HOME/c3opt_bench/res/$n" "$R/oldres/$n"
  elif [ -d "$HOME/c3fix_bench/res/$n" ]; then ln -sfn "$HOME/c3fix_bench/res/$n" "$R/oldres/$n"; fi
done
log "в очереди; эталон: $(cd "$R/oldres" && for n in *; do printf "%s->%s " "$n" "$(readlink "$n" | sed "s#$HOME/##")"; done)"
while :; do
  free=$($HOME/proj_runs/threads_free.sh 2>/dev/null | grep -o "FREE=[0-9]*" | cut -d= -f2)
  shots=$(pgrep -c -x guo2019_shot || true)
  [ "${free:-0}" = 8 ] && [ "${shots:-0}" = 0 ] && break
  sleep 300
done
log "сервер свободен (FREE=8, guo2019_shot нет) -- старт"
RUNS="$RUNS" OLD_RES=$R/oldres bash "$R/bench_8t.sh" "$HOME/SPH_c3opt2" - "$HOME/late/twaron_s25_dmg0.00415.fib" "$R/res"
rc=$?
log "bench_8t.sh rc=$rc"
exit $rc
