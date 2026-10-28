#!/usr/bin/env bash
# fetch_partial.sh <server shot dir> <local dir> <note>: скачать НЕЗАКОНЧЕННЫЙ выстрел (нет shots.csv/shot.txt),
# каждый второй кадр + логи, и написать локальную заглушку shot.txt с outcome=interrupted по последней строке
# monitor.csv (только для подписей визуализации; это НЕ результат расчёта).
set -u
SRC=$1; DST=$2; NOTE=$3
timeout 900 ssh -o ConnectTimeout=30 claude-compute "cd $SRC && r=\$(ls -d */ | head -1); r=\${r%/}; cd \$r && (ls fiber_output_* | awk 'NR%2==1'; ls sph_output_* | awk 'NR%2==1'; ls fiber_output_* | tail -1; ls sph_output_* | tail -1) | sort -u | sed \"s|^|\$r/|\" > /tmp/sel_partial.txt; cd .. && tar -czf /tmp/partial.tgz \$r/energy_log.csv \$r/monitor.csv \$r/log.txt \$r/guo2019_fabric_layout.txt -T /tmp/sel_partial.txt" < /dev/null
rm -rf "$DST"; mkdir -p "$DST"
for i in 1 2 3; do timeout 1200 scp -o ConnectTimeout=30 claude-compute:/tmp/partial.tgz "$DST/p.tgz" && break; sleep 20; done
tar -xzf "$DST/p.tgz" -C "$DST" && rm "$DST/p.tgz"
R=$(cygpath -m "$(ls -d $DST/*/ | head -1)")
python - "$R" "$NOTE" <<'EOF'
import csv, os, sys
rd, note = sys.argv[1], sys.argv[2]
m = [q for q in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if q.get("w_p")]
last = m[-1]
name = os.path.basename(os.path.normpath(rd))           # 1T_base_v68.0_r0
v0 = name.split("_v")[1].split("_")[0]
tfb = next((q["t"] for q in m if int(float(q["n_broken_v"])) + int(float(q["n_broken_h"])) > 0), "-1")
open(os.path.join(rd, "shot.txt"), "w", encoding="utf-8").write(
    "config=1T level=base V0=%s mu=0.1 outcome=interrupted V_R=nan t_end=%s n_broken_v=%s n_broken_h=%s "
    "t_first_break=%s w_f_at_first_break=nan note=%s stub=1\n" % (
        v0, last["t"], last["n_broken_v"], last["n_broken_h"], tfb, note.replace(" ", "_")))
print("заглушка shot.txt:", rd, "t_end", last["t"])
EOF
echo "FETCHED_PARTIAL $DST"
