#!/usr/bin/env bash
# shot_resume.sh ВЫХОДНОЙ_КАТАЛОГ КЛЮЧИ... : запуск guo2019_shot; если в каталоге выстрела есть чекпоинт --
# добавляет --restart. Использовать под guarded_run.sh (регистрация через register_job.sh).
# guo2019_shot возвращает 0 даже при ошибке выстрела, поэтому успех = появился shot.txt и нет чекпоинта «в работе»:
# без shot.txt код выхода 1, и guarded_run.sh не ставит метку .done (после перезагрузки задача поднимется снова).
OUT=$1; shift
CK=$(ls "$OUT"/*/checkpoint.ckpt "$OUT"/*/checkpoint.ckpt.prev 2>/dev/null | head -1)
EXTRA=()
if [ -n "$CK" ]; then
  EXTRA=(--restart)
  echo "$(date '+%F %T') $OUT: рестарт из $CK" >> "$HOME/autostart/autostart.log"
fi
"$HOME/SPH_ckpt/build-elog/guo2019_shot" --out="$OUT" "$@" "${EXTRA[@]}"
ls "$OUT"/*/shot.txt > /dev/null 2>&1 || { echo "$(date '+%F %T') $OUT: нет shot.txt -- выстрел не завершён" >> "$HOME/autostart/autostart.log"; exit 1; }
exit 0
