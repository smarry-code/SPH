#!/usr/bin/env bash
# Задача автозапуска fabric_kT (claude-compute-8t): подбор числа слоёв ткани k (конфигурация kT), при котором шар 1 г
# (d = 6,253 мм, dx = 0,3 мм) на V0 = 550 м/с ещё пробивает ткань; Guo 2019 base, закон 3 + выталкивание tau = 0,5 мкс,
# mu = 0,27, вердикт по шару. Алгоритм -- ply_search.py (пара 4T / 8T, затем деление на трети), выстрелы по 4 потока.
#
# Регистрация: ~/autostart/register_job.sh fabric_kT $W/run_series.sh
# Запуск:      setsid nohup ~/autostart/jobs/fabric_kT.job > ~/autostart/fabric_kT.out 2>&1 < /dev/null &
# После перезагрузки задачу поднимает @reboot resume_all.sh. Драйвер восстанавливается сам: готовые выстрелы
# (shots.csv) не пересчитываются, оборванные продолжаются из checkpoint.ckpt.
#
# Правила 8t: каждый выстрел -- 4 потока; ставится при uptime > 600 с, proj_status.sh != 2 и FREE >= 4 в threads_free.sh,
# через 2 мин проверяется занятость потоков (window_pinball.start_gated). OMP_NUM_THREADS=4 драйвер задаёт выстрелу сам,
# у python его нет (env -u).
#
# Аргументы (необязательные): W K_FIRST -- рабочий каталог и первая пара k (по умолчанию 4,8).
set -u
W=${1:-$HOME/fabric_kT/run}
K_FIRST=${2:-4,8}
cd "$W" || exit 1
log() { echo "$(date '+%F %T') $*" >> "$W/queue_wait.log"; }

log "старт/продолжение серии; выстрелы ставит драйвер по свободным потокам (queue.log)"
touch "$W/series.started"
env -u OMP_NUM_THREADS python3 "$W/ply_search.py" --work="$W" \
  --bin="$HOME/fabric_kT/src/build-elog/guo2019_shot" \
  --fib="$W/twaron_s25_dmg0.00415.fib" --fiber-E=8e10 --eps-break=0.05 --contact=3 --mu=0.27 --verdict-ball --stop-vz-frac=0 \
  --level=base --v0=550 --k-first="$K_FIRST" --k-step=2 \
  --extra="--dt-scale=0.5 --impulse-push-tau=5e-7 --ball-d=6.253e-3 --ball-mass=1e-3 --ball-dx=3e-4" \
  --threads=4 --t-max=4e-4 --t-max-long=1e-3 --save-t=7.5e-6 --energy-log=100 --checkpoint-every=2000 \
  --slot-script="$HOME/proj_runs/threads_free.sh" --proj-status="$HOME/proj_runs/proj_status.sh" \
  --min-uptime=600 --slot-poll=300 --slot-verify-s=120 \
  >> "$W/driver.out" 2>&1
rc=$?
log "серия завершена, rc=$rc (итог: $W/search.md)"
exit $rc
