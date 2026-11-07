#!/usr/bin/env bash
# fetch_run.sh <server shot dir (с shots.csv)> <local dir>: ждёт shots.csv, скачивает логи и каждый второй кадр.
set -u
SRC=$1; DST=$2
until timeout 40 ssh claude-compute "test -f $SRC/shots.csv" < /dev/null; do sleep 300; done
timeout 900 ssh claude-compute "cd $SRC && r=\$(ls -d */ | head -1); r=\${r%/}; cd \$r && (ls fiber_output_* | awk 'NR%2==1'; ls sph_output_* | awk 'NR%2==1'; ls fiber_output_* | tail -1; ls sph_output_* | tail -1) | sort -u | sed \"s|^|\$r/|\" > /tmp/sel_fetch.txt; cd .. && tar -czf /tmp/fetch.tgz shots.csv vbl.txt \$r/energy_log.csv \$r/monitor.csv \$r/shot.txt \$r/log.txt \$r/guo2019_fabric_layout.txt -T /tmp/sel_fetch.txt" < /dev/null
mkdir -p "$DST"
for i in 1 2 3 4 5; do timeout 1200 scp claude-compute:/tmp/fetch.tgz "$DST/fetch.tgz" && break; sleep 20; done
tar -xzf "$DST/fetch.tgz" -C "$DST" && rm "$DST/fetch.tgz"
echo "FETCHED $DST $(ls $DST/*/ | wc -l) files"
