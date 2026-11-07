"""Анализ разрывов серии с законом 4 (Guo 1T base): где, когда и в каком состоянии рвутся нити.

Для каждого выстрела с кадрами: связи, исчезнувшие между кадрами -> время разрыва (точность = шаг кадров),
нить (layout), положение разрыва вдоль нити относительно оси удара, смещение нити от центра,
напряжение/eps_p связей в кадре перед первым разрывом. Плюс энергия шара к первому разрыву (monitor).
"""
import csv, glob, json, math, os, sys

V = sys.argv[1]
OUT = sys.argv[2]
M_BALL = 0.007


def read_frame(path, want_particles):
    t = None
    pos = {}
    bonds = {}
    mode = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("TIME"):
                t = float(line.split()[1])
                continue
            if line.startswith("PARTICLES"):
                mode = "p"
                continue
            if line.startswith("BONDS"):
                mode = "b"
                continue
            if not line or line[0] == "#":
                continue
            s = line.split()
            if mode == "p" and want_particles:
                pos[int(s[0])] = (float(s[3]), float(s[4]), float(s[5]))
            elif mode == "b":
                bonds[(int(s[0]), int(s[1]))] = (float(s[5]), float(s[6]), float(s[7]))  # elongation, eps_p, sigma
    return t, pos, bonds


def layout(run):
    yarns = []
    for line in open(os.path.join(run, "guo2019_fabric_layout.txt"), encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        s = line.split()
        yarns.append((s[1], float(s[3]), int(s[5]), int(s[6])))
    return yarns


def yarn_of(yarns, i):
    for k, (dr, off, a, b) in enumerate(yarns):
        if a <= i <= b:
            return k
    return -1


res = {}
for run in sorted(glob.glob(os.path.join(V, "mu*_v*", "*", "shot.txt"))):
    rd = os.path.dirname(run)
    name = os.path.basename(os.path.dirname(rd))
    kv = dict(t.split("=", 1) for t in open(run, encoding="utf-8").read().split() if "=" in t)
    mon = list(csv.DictReader(open(os.path.join(rd, "monitor.csv"))))
    t_fb = float(kv["t_first_break"])
    r = {"V0": float(kv["V0"]), "mu": float(kv["mu"]), "outcome": kv["outcome"], "V_R": float(kv["V_R"]),
         "nv": int(float(kv["n_broken_v"])), "nh": int(float(kv["n_broken_h"])), "t_fb_us": t_fb * 1e6,
         "w_fb_mm": float(kv["w_f_at_first_break"]) * 1e3, "t_end_us": float(kv["t_end"]) * 1e6,
         "W_p": float(kv["W_p"]), "eps_p_max": float(kv["eps_p_max"]), "max_dff": float(kv["max_delta_ff/d"]),
         "max_dfs": float(kv["max_delta_fs/d"])}
    E0 = 0.5 * M_BALL * r["V0"] ** 2
    row = min(mon, key=lambda q: abs(float(q["t"]) - t_fb))
    vfb = float(row["v_z"])
    r["v_fb"] = vfb
    r["E_ball_fb_frac"] = 0.5 * M_BALL * vfb ** 2 / E0
    # хронология разрывов нитей (monitor)
    tl = []
    prev = 0
    for q in mon:
        n = int(float(q["n_broken_v"])) + int(float(q["n_broken_h"]))
        if n > prev:
            tl.append((round(float(q["t"]) * 1e6, 1), n))
            prev = n
    r["timeline"] = tl
    frames = sorted(glob.glob(os.path.join(rd, "fiber_output_*.txt")))
    if frames:
        yarns = layout(rd)
        t0, pos0, b0 = read_frame(frames[0], True)
        prev_b = b0
        prev_t = t0
        before = None
        breaks = []
        for fp in frames[1:]:
            t, _, bb = read_frame(fp, False)
            gone = [k for k in prev_b if k not in bb]
            if gone and before is None:
                before = (prev_t, prev_b)
            for (i, j) in gone:
                y = yarn_of(yarns, i)
                dr, off = yarns[y][0], yarns[y][1]
                x, yy, _ = pos0[i]
                along = yy if dr == "y" else x
                breaks.append({"t_us": round(t * 1e6, 1), "yarn": y, "dir": dr, "offset_mm": round(off * 1e3, 2),
                               "along_mm": round(along * 1e3, 2), "eps_prev": prev_b[(i, j)][0],
                               "sig_prev": prev_b[(i, j)][2]})
            prev_b, prev_t = bb, t
        r["breaks"] = breaks
        if before:
            tb, bb = before
            sig = sorted((v[2] for v in bb.values()), reverse=True)
            r["pre_break_t_us"] = tb * 1e6
            r["pre_break_sig_top"] = sig[:5]
            r["pre_break_n_plastic"] = sum(1 for v in bb.values() if v[1] > 0)
            r["pre_break_n_sig_gt_2GPa"] = sum(1 for v in bb.values() if v[2] > 2.0e9)
    res[name] = r
    print(name, json.dumps({k: v for k, v in r.items() if k not in ("breaks",)}, ensure_ascii=False))
    for b in r.get("breaks", [])[:40]:
        print("   ", b)
json.dump(res, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
