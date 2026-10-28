"""Каналы энергии (% E0) в момент первого разрыва и в конце: pinball / Swegle / закон 4, 1T, mu = 0.1."""
import csv, glob, os, sys

for label, root in zip(sys.argv[1::2], sys.argv[2::2]):
    for f in sorted(glob.glob(os.path.join(root, "mu0.1_v7[06]", "*", "energy_log.csv"))):
        d = os.path.basename(os.path.dirname(os.path.dirname(f)))
        kv = dict(t.split("=", 1) for t in open(os.path.join(os.path.dirname(f), "shot.txt"), encoding="utf-8").read().split()
                  if "=" in t)
        r = list(csv.DictReader(open(f)))
        E0 = float(r[0]["E_mech"])
        tfb = float(kv["t_first_break"])
        a = min(r, key=lambda q: abs(float(q["time"]) - tfb))
        b = r[-1]
        g = lambda q, k: ("%5.1f" % (100 * float(q[k]) / E0)) if k in q and q[k] != "" else "  —  "
        keys = ["KE_solid", "KE_fiber", "U_bond", "U_hertz_ff", "D_friction", "D_contact_fiber", "W_p_fiber", "residual"]
        print("%-8s %s E0=%.2f | fb: %s | end: %s" % (label, d, E0, " ".join(g(a, k) for k in keys),
                                                     " ".join(g(b, k) for k in keys)))
print("keys:", "KE_solid KE_fiber U_bond U_hertz_ff D_friction D_contact_fiber W_p_fiber residual")
