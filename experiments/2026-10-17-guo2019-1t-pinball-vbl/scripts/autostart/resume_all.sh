#!/usr/bin/env bash
# Вызывается при загрузке сервера (crontab @reboot) и вручную. Для каждой зарегистрированной задачи
# ~/autostart/jobs/ИМЯ.job (исполняемый скрипт, который сам вызывает guarded_run.sh ИМЯ ...) без ИМЯ.done
# запускает её в фоне. Скрипт идемпотентен: занятый замок guarded_run.sh отсекает дубли.
A=$HOME/autostart
mkdir -p "$A/jobs"
echo "$(date '+%F %T') resume_all: проверка задач" >> "$A/autostart.log"
sleep "${RESUME_DELAY:-20}"   # дать системе подняться (сеть, диск)
for job in "$A"/jobs/*.job; do
  [ -e "$job" ] || continue
  name=$(basename "$job" .job)
  if [ -e "$A/$name.done" ]; then continue; fi
  echo "$(date '+%F %T') resume_all: поднимаю $name" >> "$A/autostart.log"
  setsid nohup "$job" > "$A/$name.out" 2>&1 < /dev/null &
done
