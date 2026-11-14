#!/bin/bash
# Ждёт, пока текущий экземпляр fabric_w50_search отпустит замок (20T «не опр.» к 400 мкс -> поиск стоп),
# затем снимает его итог и метку .done и поднимает задачу заново -- уже с продолжением до 1000 мкс.
A=$HOME/autostart; S=$HOME/fabric_kT/w50_search
exec 8>"$A/fabric_w50_search.lock"; flock 8; flock -u 8; exec 8>&-
echo "$(date "+%F %T") перезапуск поиска (продолжение «не опр.» до 1000 мкс; было: $(cat $S/search_result.txt 2>/dev/null))" >> $S/search.log
mv -f $S/search_result.txt $S/search_result_stop400.txt 2>/dev/null
rm -f "$A/fabric_w50_search.done"
exec "$A/jobs/fabric_w50_search.job" >> "$A/fabric_w50_search.out" 2>&1
