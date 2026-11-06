#!/usr/bin/env python3
"""diff_frames.py DIR_A DIR_B [--every N] [--top K]
Сравнение кадров нитей двух прогонов по бусинам (кадры с SPH_FRAME_PRECISION=17).
Для каждого общего кадра: max и RMS |dx|, max |dv|, где максимум (нить, бусина, расстояние
от оси шара), max |d sigma| по связям, доля смен знака dz вдоль нити с максимумом
(чередование через бусину -> около 1), число бусин с |dx| > 1e-3 max.
"""
import glob, math, os, re, sys

def frames(d):
    out = {}
    for f in glob.glob(os.path.join(d, "fiber_output_*.txt")):
        out[int(re.search(r"_(\d+)\.txt$", f).group(1))] = f
    return out

def read_frame(path):
    with open(path) as fh:
        head = fh.readline().split()
        t = float(head[1]); n = int(head[5]); nb = int(head[7])
        fh.readline(); fh.readline()
        X = [None] * n; V = [None] * n; obj = [0] * n
        for i in range(n):
            p = fh.readline().split()
            obj[i] = int(p[1])
            X[i] = (float(p[3]), float(p[4]), float(p[5]))
            V[i] = (float(p[6]), float(p[7]), float(p[8]))
        fh.readline(); fh.readline()
        B = [None] * nb
        for k in range(nb):
            p = fh.readline().split()
            B[k] = (int(p[0]), int(p[1]), float(p[7]) if len(p) > 7 else float(p[5]) * 8e10)
    return t, X, V, obj, B

def layout(d):
    lay = {}
    f = glob.glob(os.path.join(d, "guo2019_fabric_layout.txt"))
    if not f: return lay
    for line in open(f[0]):
        if line.startswith("#"): continue
        p = line.split()
        lay[int(p[0])] = (p[1], int(p[2]), int(p[5]), int(p[6]))
    return lay

def ball_axis(d):
    f = sorted(glob.glob(os.path.join(d, "sph_output_*.txt")))
    if not f: return (0.0, 0.0)
    sx = sy = 0.0; n = 0
    with open(f[0]) as fh:
        fh.readline()
        for line in fh:
            p = line.split()
            if len(p) < 6: continue
            sx += float(p[3]); sy += float(p[4]); n += 1
    return (sx / n, sy / n) if n else (0.0, 0.0)

def main():
    a, b = sys.argv[1], sys.argv[2]
    every = int(sys.argv[sys.argv.index("--every") + 1]) if "--every" in sys.argv else 1
    fa, fb = frames(a), frames(b)
    steps = sorted(set(fa) & set(fb))[::every]
    lay = layout(a); cx, cy = ball_axis(a)
    print("# step t max_dx rms_dx max_dv fiber(dir,idx) bead_in_fiber r_axis_mm max_dsig bond(i-j) sig_ref alt_frac n_spread")
    for s in steps:
        t, Xa, Va, obj, Ba = read_frame(fa[s])
        _, Xb, Vb, _, Bb = read_frame(fb[s])
        n = len(Xa)
        dx = [math.dist(Xa[i], Xb[i]) for i in range(n)]
        dv = max(math.dist(Va[i], Vb[i]) for i in range(n))
        im = max(range(n), key=dx.__getitem__); mx = dx[im]
        rms = math.sqrt(sum(x * x for x in dx) / n)
        ds = [abs(Ba[k][2] - Bb[k][2]) for k in range(len(Ba))]
        km = max(range(len(ds)), key=ds.__getitem__)
        o = obj[im]; dirn, idx, f0, f1 = lay.get(o, ("?", o, 0, n - 1))
        dz = [Xa[i][2] - Xb[i][2] for i in range(f0, f1 + 1)]
        nz = [z for z in dz if z != 0.0]
        alt = sum(1 for u, w in zip(nz, nz[1:]) if u * w < 0) / max(1, len(nz) - 1)
        spread = sum(1 for x in dx if x > 1e-3 * mx) if mx > 0 else 0
        r = math.hypot(Xa[im][0] - cx, Xa[im][1] - cy) * 1e3
        print(f"{s} {t:.6e} {mx:.3e} {rms:.3e} {dv:.3e} {dirn}{idx} {im - f0} {r:.2f} "
              f"{ds[km]:.3e} {Ba[km][0]}-{Ba[km][1]} {Ba[km][2]:.4e} {alt:.2f} {spread}", flush=True)

if __name__ == "__main__":
    main()
