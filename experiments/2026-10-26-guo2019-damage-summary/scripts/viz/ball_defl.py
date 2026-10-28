"""Прогиб ткани под шаром = ход низа шара за плоскость ткани: z_fab0 - (z_c - r) по monitor.csv; z_fab0 и r -- по
кадру 0 (sph_output / fiber_output). Для остановленного шара -- максимум и его момент.
Usage: python ball_defl.py label=run_dir ..."""
import csv, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from animate_1t import discover_frame_ids, parse_sph_frame, parse_fiber_frame
for a in sys.argv[1:]:
    lab, rd = a.split("=", 1)
    i0 = discover_frame_ids(rd)[0]
    b = parse_sph_frame(os.path.join(rd, "sph_output_%s.txt" % i0))
    f = parse_fiber_frame(os.path.join(rd, "fiber_output_%s.txt" % i0))
    r = 0.5 * np.ptp(b[:, 2]) + 0.5 * 1.2e-3 / 2  # радиус по частицам + полчастицы (dx 0,6 мм)
    zf = np.median(f.pos[:, 2])
    m = [q for q in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if q.get("w_p")]
    t = np.array([float(q["t"]) for q in m]); zc = np.array([float(q["z_c"]) for q in m])
    vz = np.array([float(q["v_z"]) for q in m])
    d = zf - (zc - r)
    k = int(np.argmax(d))
    stop = np.where(vz >= 0)[0]
    print("%s gap0=%.2f мм  max=%.1f мм при t=%.0f мкс  v_z=0 при t=%s" % (lab, -d[0] * 1e3, d[k] * 1e3, t[k] * 1e6,
          "%.0f мкс, прогиб %.1f мм" % (t[stop[0]] * 1e6, d[stop[0]] * 1e3) if len(stop) else "—"))
