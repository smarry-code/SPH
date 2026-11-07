"""Где оторванные бусины (без связей): расстояние до центра шара и скорость относительно шара.
Центр шара -- средняя координата частиц sph_output того же шага; скорость шара -- средняя скорость частиц."""
import glob, math, os, sys


def read_fiber(path):
    pos, vel, deg = [], [], {}
    mode = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("PARTICLES"):
                mode = "p"; continue
            if line.startswith("BONDS"):
                mode = "b"; continue
            if not line.strip() or line[0] in "#T":
                continue
            s = line.split()
            if mode == "p":
                pos.append((float(s[3]), float(s[4]), float(s[5])))
                vel.append((float(s[6]), float(s[7]), float(s[8])))
            elif mode == "b":
                for k in (int(s[0]), int(s[1])):
                    deg[k] = deg.get(k, 0) + 1
    return pos, vel, deg


def read_ball(path):
    xs, vs = [], []
    with open(path, encoding="utf-8", errors="replace") as f:
        hdr = None
        for line in f:
            if hdr is None and " vx " in line:
                hdr = line.lstrip("#").split(); continue
            s = line.split()
            if hdr and len(s) == len(hdr):
                try:
                    d = dict(zip(hdr, map(float, s)))
                except ValueError:
                    continue
                xs.append((d["x"], d["y"], d["z"])); vs.append((d["vx"], d["vy"], d["vz"]))
    n = len(xs)
    c = tuple(sum(p[i] for p in xs) / n for i in range(3))
    v = tuple(sum(p[i] for p in vs) / n for i in range(3))
    return c, v, n


for label, root in zip(sys.argv[1::2], sys.argv[2::2]):
    ff = sorted(glob.glob(os.path.join(root, "fiber_output_*.txt")))[-1]
    sf = ff.replace("fiber_output_", "sph_output_")
    pos, vel, deg = read_fiber(ff)
    c, vb, nb = read_ball(sf)
    free = [k for k in range(len(pos)) if deg.get(k, 0) == 0]
    inside = near = 0
    rel = []
    for k in free:
        r = math.dist(pos[k], c)
        if r < 6.0e-3: inside += 1
        elif r < 8.0e-3: near += 1
        rel.append(math.dist(vel[k], vb))
    rel.sort()
    med = rel[len(rel) // 2] if rel else float("nan")
    print("%-8s кадр %s: шар n=%d центр z=%.1f мм v_z=%.1f; бусин без связей %d: внутри шара (r<6 мм) %d, у поверхности (6-8 мм) %d, дальше %d; |v - v_шара| медиана %.0f м/с, макс %.0f" % (
        label, os.path.basename(ff), nb, c[2] * 1e3, vb[2], len(free), inside, near, len(free) - inside - near, med,
        rel[-1] if rel else float("nan")))
