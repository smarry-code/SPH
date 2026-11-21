"""Разбор одного выстрела: исход, каналы энергии, v_z(t), разрывы по времени, свободные и быстрые бусины по кадрам."""
import csv, glob, math, os, sys

rd = sys.argv[1]
kv = dict(t.split("=", 1) for t in open(os.path.join(rd, "shot.txt"), encoding="utf-8").read().split() if "=" in t)
keys = ["outcome", "V_R", "t_end", "n_broken_v", "n_broken_h", "t_first_break", "w_f_at_first_break", "w_f_max",
        "max_delta_ff/d", "max_delta_fs/d", "n_tunnel", "contacts_per_crossing", "eps_p_max", "W_p",
        "n_bonds_plastic", "dt_mean", "wall_s", "fiber_E", "mu_ff", "mu_fs", "pl_sy", "pl_H", "fail_sigma_u",
        "fiber_eps_break", "n_beads", "n_bonds", "n_clamped_f", "F_clamp_f_z_max"]
for k in keys:
    print("%-24s %s" % (k, kv.get(k)))

e = list(csv.DictReader(open(os.path.join(rd, "energy_log.csv"))))
print("energy columns:", list(e[0].keys()))
E0 = float(e[0]["E_mech"])
print("E0 = %.2f J" % E0)
tfb = float(kv["t_first_break"])
ch = [c for c in ("KE_solid", "KE_fiber", "U_bond", "U_hertz_ff", "U_hertz_fs", "D_friction", "D_contact_fiber",
                  "W_p_fiber", "W_break", "residual") if c in e[0]]
for tt in (100e-6, 200e-6, tfb, 600e-6, 800e-6, 1000e-6, float(e[-1]["time"])):
    a = min(e, key=lambda q: abs(float(q["time"]) - tt))
    print("t=%5.0f мкс  " % (float(a["time"]) * 1e6) + " ".join("%s=%.1f%%" % (k, 100 * float(a[k]) / E0) for k in ch))

m = [q for q in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if q.get("v_z")]
print("monitor columns:", list(m[0].keys()))
want = list(range(50, 1351, 50))
got = {}
for q in m:
    t = float(q["t"]) * 1e6
    for tt in want:
        if abs(t - tt) < 1.5 and tt not in got:
            got[tt] = q
print("v_z / w_p / разрывы по времени:")
for tt in sorted(got):
    q = got[tt]
    print("  %5d мкс  v_z=%7.2f  " % (tt, float(q["v_z"])) +
          " ".join("%s=%s" % (k, q[k]) for k in q if k in ("w_p", "w_f", "n_broken", "n_broken_v", "n_broken_h", "z")))

frames = sorted(glob.glob(os.path.join(rd, "fiber_output_*.txt")))
print("кадров:", len(frames))
for ff in frames[::6] + [frames[-1]]:
    speeds, deg, mode, t = [], {}, None, None
    for line in open(ff, encoding="utf-8", errors="replace"):
        if line.startswith("PARTICLES"): mode = "p"; continue
        if line.startswith("BONDS"): mode = "b"; continue
        if line.startswith("T") or line.startswith("#"):
            if t is None:
                for tok in line.replace("=", " ").split():
                    try: t = float(tok); break
                    except ValueError: pass
            continue
        if not line.strip(): continue
        s = line.split()
        if mode == "p":
            speeds.append(math.sqrt(sum(float(x) ** 2 for x in s[6:9])))
        elif mode == "b":
            for k in (int(s[0]), int(s[1])): deg[k] = deg.get(k, 0) + 1
    free = [k for k in range(len(speeds)) if deg.get(k, 0) == 0]
    print("  %s: t=%s свободных %d, |v|>100: %d, |v|>300: %d, max|v| %.0f" % (
        os.path.basename(ff), t, len(free), sum(v > 100 for v in speeds), sum(v > 300 for v in speeds), max(speeds)))
