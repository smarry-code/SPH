#!/usr/bin/env bash
# Задача автозапуска mu027_window (claude-compute-8t): подбор окна непробитие / пробитие Guo 2019 1T base,
# закон 3 (импульсный) на нить-нить и нить-шар + выталкивание tau = 0,5 мкс, mu = 0,27, вердикт по шару,
# V0 в [71; 82] м/с, деление на трети, 2 выстрела x 4 потока (решения владельца 2026-11-08).
# Код: fiber_new ec2439f + ветка guo1t_law3push_mu027 (~/mu027/src/SOURCE_COMMIT).
#
# Регистрация: ~/autostart/register_job.sh mu027_window $W/run_series.sh
# Запуск:      setsid nohup ~/autostart/jobs/mu027_window.job > ~/autostart/mu027_window.out 2>&1 < /dev/null &
# После перезагрузки задачу поднимает @reboot resume_all.sh (метки .done нет до rc = 0). Драйвер восстанавливается
# сам: готовые выстрелы (shots.csv) не пересчитываются, оборванные продолжаются из checkpoint.ckpt.
set -u
W=$HOME/vbl1T_law3push_mu027_det
BIN=$HOME/mu027/src/build-elog/guo2019_shot
SELF=mu027_window
cd "$W" || exit 1
log() { echo "$(date '+%F %T') $*" >> "$W/queue_wait.log"; }

# --- ожидание свободного сервера (правила 8t): proj_* завершены (rc 0), нет guo2019_shot, у всех чужих задач
#     автозапуска есть .done, load average за 1 мин < 1. При rc 2 серии proj_* сервер не занимаем.
if [ -f "$W/series.started" ] && [ -f "$W/window.json" ]; then
  log "рестарт задачи: серия уже шла -- жду освобождения сервера от чужих задач и продолжаю"
fi
last=""
while true; do
  "$HOME/proj_runs/proj_status.sh" > /dev/null 2>&1
  prc=$?
  [ -x "$HOME/proj_runs/proj_status.sh" ] || prc=0        # серия proj_* убрана с сервера -- не ждём её
  nshot=$(pgrep -c -x guo2019_shot)
  pend=""
  for j in "$HOME"/autostart/jobs/*.job; do
    [ -e "$j" ] || continue
    n=$(basename "$j" .job)
    [ "$n" = "$SELF" ] && continue
    [ -e "$HOME/autostart/$n.done" ] || pend="$pend $n"
  done
  load=$(cut -d' ' -f1 /proc/loadavg)
  lowload=$(awk -v l="$load" 'BEGIN { print (l < 1.0) ? 1 : 0 }')
  if [ "$prc" = 0 ] && [ "$nshot" = 0 ] && [ -z "$pend" ] && [ "$lowload" = 1 ]; then
    break
  fi
  state="proj rc=$prc, guo2019_shot=$nshot, load=$load, без .done:${pend:- нет}"
  [ "$prc" = 2 ] && state="$state -- proj: ОСТАНОВЛЕНО С ОШИБКОЙ, жду решения владельца"
  # в лог -- при смене состояния и раз в час
  if [ "$state" != "$last" ] || [ $(( $(date +%s) / 3600 )) != "${hour:-}" ]; then
    log "жду: $state"
    last=$state
    hour=$(( $(date +%s) / 3600 ))
  fi
  sleep 300
done

log "сервер свободен -- старт/продолжение серии"
touch "$W/series.started"
OMP_NUM_THREADS=4 python3 "$W/window_bisect.py" --work="$W" --bin="$BIN" \
  --fib="$W/twaron_s25_dmg0.00415.fib" --fiber-E=8e10 --eps-break=0.05 --config=1T --level=base \
  --contact=3 --mu=0.27 --verdict-ball --vmin=71 --v-first=82 --vmax=82 --tol=1.5 \
  --extra="--dt-scale=0.5 --impulse-push-tau=5e-7" --threads=4 --threads-single=4 \
  --t-max=1.5e-3 --t-max-long=3e-3 --save-t=2e-5 --energy-log=100 --checkpoint-every=20000 \
  >> "$W/driver.out" 2>&1
rc=$?
log "серия завершена, rc=$rc (итог: $W/window.md)"
exit $rc
