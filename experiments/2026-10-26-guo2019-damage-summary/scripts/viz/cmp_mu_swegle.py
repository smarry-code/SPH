"""Swegle, 1T, E = 78 ГПа, V0 = 70: mu = 0.1 против mu = 0.3 -- исход, энергия, v_z(t), разлёт бусин."""
import csv, glob, math, os, sys

L = sys.argv[1]
runs = [("mu=0.1", "swegle_mu0.1_v70_E78"), ("mu=0.3", "swegle_mu0.3_v70_E78")]
keys = ["outcome", "V_R", "t_end", "n_broken_v", "n_broken_h", "t_first_break", "w_f_at_first_break",
        "max_delta_ff/d", "max_delta_fs/d", "n_tunnel", "W_p", "wall_s"]
for name, d in runs:
    rd = glob.glob(os.path.join(L, d, "*"))[0]
    kv = dict(t.split("=", 1) for t in open(os.path.join(rd, "shot.txt"), encoding="utf-8").read().split() if "=" in t)
    print("==", name, " ".join("%s=%s" % (k, kv.get(k)) for k in keys))
    e = list(csv.DictReader(open(os.path.join(rd, "energy_log.csv"))))
    E0 = float(e[0]["E_mech"])
    tfb = float(kv["t_first_break"])
    a = min(e, key=lambda q: abs(float(q["time"]) - tfb))
    b = e[-1]
    ch = ("KE_solid", "KE_fiber", "U_bond", "D_friction", "D_contact_fiber", "W_p_fiber", "residual")
    print("   к первому разрыву, %% E0: " + " ".join("%s=%.1f" % (k, 100 * float(a[k]) / E0) for k in ch))
    print("   в конце,          %% E0: " + " ".join("%s=%.1f" % (k, 100 * float(b[k]) / E0) for k in ch))
    m = [q for q in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if q.get("w_p")]
    vz = {}
    for q in m:
        t = float(q["t"]) * 1e6
        for tt in (200, 300, 400, 450, 500, 550, 600, 700, 800, 1000):
            if abs(t - tt) < 1.0 and tt not in vz:
                vz[tt] = float(q["v_z"])
    print("   v_z(t): " + ", ".join("%d мкс: %.1f" % (t, v) for t, v in sorted(vz.items())))
    # разлёт: последний кадр
    ff = sorted(glob.glob(os.path.join(rd, "fiber_output_*.txt")))[-1]
    speeds, deg, mode = [], {}, None
    for line in open(ff, encoding="utf-8", errors="replace"):
        if line.startswith("PARTICLES"): mode = "p"; continue
        if line.startswith("BONDS"): mode = "b"; continue
        if not line.strip() or line[0] in "#T": continue
        s = line.split()
        if mode == "p":
            speeds.append(math.sqrt(sum(float(x) ** 2 for x in s[6:9])))
        elif mode == "b":
            for k in (int(s[0]), int(s[1])): deg[k] = deg.get(k, 0) + 1
    free = sum(1 for k in range(len(speeds)) if deg.get(k, 0) == 0)
    print("   последний кадр: бусин без связей %d, |v|>100 м/с: %d, max|v| %.0f м/с" % (
        free, sum(v > 100 for v in speeds), max(speeds)))
