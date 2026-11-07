"""Где и в каком кадре рвутся нити: по кадрам fiber_output_* -- для каждой порванной нити (пропала хотя бы одна
связь) номер нити, направление, смещение нити от центра, место разрыва вдоль нити (от центра удара), кадр
разрыва; для повреждённости -- связи с 0 < D < 1 в каждом кадре (сколько кадров живёт разупрочнение).
Usage: python break_where.py <run_dir>
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from animate_1t import discover_frame_ids, parse_fiber_frame  # noqa: E402

rd = sys.argv[1]
ids = discover_frame_ids(rd)
lay = []
for ln in open(os.path.join(rd, "guo2019_fabric_layout.txt")):
    p = ln.split()
    if len(p) >= 7 and not p[0].startswith("#"):
        lay.append((p[1], int(p[2]), float(p[3]), int(p[5]), int(p[6])))
nb = max(l[4] for l in lay) + 1
yarn_of = np.empty(nb, dtype=int)
for k, l in enumerate(lay):
    yarn_of[l[3]:l[4] + 1] = k


def bonds_with_d(path):
    f = parse_fiber_frame(path)
    with open(path) as fh:
        for ln in fh:
            if ln.strip() == "BONDS":
                break
        rows = [next(fh).split() for _ in range(len(f.ij))]
    d = np.array([float(r[10]) if len(r) > 10 else 0.0 for r in rows])
    return f, d


f0, _ = bonds_with_d(os.path.join(rd, "fiber_output_%s.txt" % ids[0]))
ctr = 0.5 * (f0.pos.min(axis=0) + f0.pos.max(axis=0))
alive_prev = {tuple(b) for b in f0.ij}
broken = {}  # yarn -> (frame t, x along, bond)
dmg_life = {}
for fid in ids[1:]:
    f, d = bonds_with_d(os.path.join(rd, "fiber_output_%s.txt" % fid))
    alive = {tuple(b) for b in f.ij}
    for b in alive_prev - alive:
        y = yarn_of[b[0]]
        if y not in broken:
            p = f0.pos[b[0]]
            along = (p[1] - ctr[1]) if lay[y][0] == "y" else (p[0] - ctr[0])
            broken[y] = (f.time * 1e6, along * 1e3)
    for k in np.nonzero((d > 0) & (d < 1))[0]:
        key = tuple(f.ij[k])
        dmg_life.setdefault(key, []).append((round(f.time * 1e6), round(float(d[k]), 3)))
    alive_prev = alive
print("кадров %d, шаг кадра %.1f мкс" % (len(ids), (f.time - f0.time) * 1e6 / max(len(ids) - 1, 1)))
print("порванные нити (направл., смещение нити от центра, мм; место разрыва вдоль нити от центра, мм; кадр, мкс):")
for y, (t, al) in sorted(broken.items(), key=lambda kv: kv[1][0]):
    print("  %s  смещение %+6.2f   разрыв на %+6.2f   t <= %4.0f" % (lay[y][0], (lay[y][2]) * 1e3, al, t))
print("связи с 0 < D < 1 (связь: [(t мкс, D), ...]):")
for k, v in sorted(dmg_life.items(), key=lambda kv: kv[1][0][0]):
    print("  ", k, v)
