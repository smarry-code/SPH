#!/usr/bin/env bash
# Задача автозапуска mu027_window (claude-compute-8t): подбор окна непробитие / пробитие Guo 2019 1T base,
# закон 3 (импульсный) на нить-нить и нить-шар + выталкивание tau = 0,5 мкс, mu = 0,27, вердикт по шару,
# V0 в [71; 82] м/с, деление на трети, выстрелы по 4 потока, не больше двух сразу (решения владельца 2026-11-08).
# Код: fiber_new ec2439f + ветка guo1t_law3push_mu027 (~/mu027/src/SOURCE_COMMIT).
#
# Регистрация: ~/autostart/register_job.sh mu027_window $W/run_series.sh
# Запуск:      setsid nohup ~/autostart/jobs/mu027_window.job > ~/autostart/mu027_window.out 2>&1 < /dev/null &
# После перезагрузки задачу поднимает @reboot resume_all.sh (метки .done нет до rc = 0). Драйвер восстанавливается
# сам: готовые выстрелы (shots.csv) не пересчитываются, оборванные продолжаются из checkpoint.ckpt.
#
# Правила 8t (инструкция владельца 2026-11-08 04:05): каждый выстрел -- 4 потока OpenMP; ставится, только если
# uptime > 600 с, ~/proj_runs/proj_status.sh != 2 и ~/proj_runs/threads_free.sh даёт FREE >= 4; через 2 мин
# проверяется, что выстрел занял 4 потока. Это делает драйвер перед каждым запуском guo2019_shot
# (window_pinball.start_gated). OMP_NUM_THREADS=4 драйвер задаёт каждому выстрелу сам, а у python его нет
# (env -u): threads_free.sh считает потоки по OMP_NUM_THREADS в окружении и принял бы драйвер за 4 потока.
#
# Аргументы (необязательные): W VMIN VFIRST VMAX -- рабочий каталог и скорости (по умолчанию первая серия [71; 82]).
#   Пример (2026-11-08, пара 69 / 70): register_job.sh mu027_v69_70 ~/vbl1T_law3push_mu027_v69_70/run_series.sh #            ~/vbl1T_law3push_mu027_v69_70 69 70 70
set -u
W=${1:-$HOME/vbl1T_law3push_mu027_det}
VMIN=${2:-71}; VFIRST=${3:-82}; VMAX=${4:-82}
BIN=$HOME/mu027/src/build-elog/guo2019_shot
cd "$W" || exit 1
log() { echo "$(date '+%F %T') $*" >> "$W/queue_wait.log"; }

log "старт/продолжение серии; выстрелы ставит драйвер по свободным потокам (queue.log)"
touch "$W/series.started"
env -u OMP_NUM_THREADS python3 "$W/window_bisect.py" --work="$W" --bin="$BIN" \
  --fib="$W/twaron_s25_dmg0.00415.fib" --fiber-E=8e10 --eps-break=0.05 --config=1T --level=base \
  --contact=3 --mu=0.27 --verdict-ball --vmin="$VMIN" --v-first="$VFIRST" --vmax="$VMAX" --tol=1.5 \
  --extra="--dt-scale=0.5 --impulse-push-tau=5e-7" --threads=4 --threads-single=4 \
  --t-max=1.5e-3 --t-max-long=3e-3 --save-t=2e-5 --energy-log=100 --checkpoint-every=20000 \
  --slot-script="$HOME/proj_runs/threads_free.sh" --proj-status="$HOME/proj_runs/proj_status.sh" \
  --min-uptime=600 --slot-poll=300 --slot-verify-s=120 \
  >> "$W/driver.out" 2>&1
rc=$?
log "серия завершена, rc=$rc (итог: $W/window.md)"
exit $rc
