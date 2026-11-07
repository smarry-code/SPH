#!/usr/bin/env bash
# Запуск задачи под замком: guarded_run.sh ИМЯ КОМАНДА [АРГУМЕНТЫ...]
# - flock ~/autostart/ИМЯ.lock: второй экземпляр той же задачи не стартует (защита от дублей, которые пишут в один
#   каталог выстрела; guo2019_shot такого не допускает);
# - при rc = 0 создаётся ~/autostart/ИМЯ.done, и resume_all.sh больше эту задачу не поднимает.
set -u
A=$HOME/autostart
name=$1; shift
mkdir -p "$A"
exec 9>"$A/$name.lock"
flock -n 9 || { echo "$(date '+%F %T') $name: уже запущена, выхожу" >> "$A/autostart.log"; exit 0; }
echo "$(date '+%F %T') $name: старт ($*)" >> "$A/autostart.log"
"$@"
rc=$?
echo "$(date '+%F %T') $name: завершено rc=$rc" >> "$A/autostart.log"
[ $rc -eq 0 ] && touch "$A/$name.done"
exit $rc
