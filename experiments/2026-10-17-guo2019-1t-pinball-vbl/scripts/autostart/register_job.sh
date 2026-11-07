#!/usr/bin/env bash
# register_job.sh ИМЯ КОМАНДА... -- записать задачу в ~/autostart/jobs/ИМЯ.job (снять метку .done, если была).
# Сама задачу не запускает: после регистрации выполните  ~/autostart/jobs/ИМЯ.job  (или resume_all.sh).
set -eu
A=$HOME/autostart
name=$1; shift
mkdir -p "$A/jobs"
rm -f "$A/$name.done"
{
  echo '#!/usr/bin/env bash'
  printf 'exec %s/guarded_run.sh %q' "$A" "$name"
  for a in "$@"; do printf ' %q' "$a"; done
  echo
} > "$A/jobs/$name.job"
chmod +x "$A/jobs/$name.job"
echo "зарегистрировано: $A/jobs/$name.job"
