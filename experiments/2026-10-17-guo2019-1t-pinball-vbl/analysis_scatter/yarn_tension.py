"""Распределение напряжения вдоль центральных вертикальных нитей в последнем кадре перед первым разрывом:
доля связей с sigma >= 0.9*max и >= 0.95*max, длина зоны с sigma >= 0.95*max (мм)."""
import glob, os, sys

for label, root in zip(sys.argv[1::2], sys.argv[2::2]):
    layout = [l.split() for l in open(os.path.join(root, "guo2019_fabric_layout.txt"), encoding="utf-8") if not l.startswith("#") and l.strip()]
    yarns = {int(s[0]): (s[1], float(s[3]), int(s[5]), int(s[6])) for s in layout}
    frames = sorted(glob.glob(os.path.join(root, "fiber_output_*.txt")))
    nb0 = None
    chosen = None
    for ff in frames:
        with open(ff, encoding="utf-8", errors="replace") as f:
            h = f.readline().split()
        nb = int(h[h.index("N_BONDS") + 1])
        if nb0 is None: nb0 = nb
        if nb < nb0: break
        chosen = ff
    t = float(open(chosen, encoding="utf-8", errors="replace").readline().split()[1])
    sig = {}
    mode = None
    for line in open(chosen, encoding="utf-8", errors="replace"):
        if line.startswith("BONDS"): mode = "b"; continue
        if mode == "b" and line.strip() and line[0] != "#":
            s = line.split()
            sig[int(s[0])] = (float(s[7]), float(s[2]))
    print("== %s, кадр t = %.0f мкс (последний до первого разрыва)" % (label, t * 1e6))
    for y in (150, 151, 152, 153, 154, 155):
        d, off, a, b = yarns[y]
        vals = [sig[k] for k in range(a, b) if k in sig]
        smax = max(v[0] for v in vals)
        n90 = sum(v[0] >= 0.9 * smax for v in vals)
        z95 = sum(v[1] for v in vals if v[0] >= 0.95 * smax) * 1e3
        z90 = sum(v[1] for v in vals if v[0] >= 0.90 * smax) * 1e3
        print("  нить %d (%s, смещение %+.2f мм): max sigma %.2f ГПа; связей >= 0.9 max: %d из %d; длина зоны >=0.95 max %.1f мм, >=0.9 max %.1f мм" % (
            y, d, off * 1e3, smax / 1e9, n90, len(vals), z95, z90))
