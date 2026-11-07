"""Распределение энергии (доли E0, %) в момент первого разрыва и в конце, серия с законом 4."""
import csv, glob, os, sys

V = sys.argv[1]
print("run        E0  | первый разрыв: KEs  KEf  Ubond Uh  fric damp Wp  | конец: KEs KEf Ubond fric damp Wp")
for f in sorted(glob.glob(os.path.join(V, "mu*_v*", "*", "energy_log.csv"))):
    d = os.path.basename(os.path.dirname(os.path.dirname(f)))
    kv = dict(t.split("=", 1) for t in open(os.path.join(os.path.dirname(f), "shot.txt"), encoding="utf-8").read().split()
              if "=" in t)
    tfb = float(kv["t_first_break"])
    r = list(csv.DictReader(open(f)))
    E0 = float(r[0]["E_mech"])
    a = min(r, key=lambda q: abs(float(q["time"]) - tfb))
    b = r[-1]
    p = lambda q, k: 100 * float(q[k]) / E0
    uh = lambda q: p(q, "U_hertz_ff") + p(q, "U_hertz_fs")
    print("%-10s %.1f | %4.0f %4.0f %5.1f %4.1f %4.0f %4.0f %4.1f | %4.0f %4.0f %5.1f %4.0f %4.0f %4.1f" % (
        d, E0, p(a, "KE_solid"), p(a, "KE_fiber"), p(a, "U_bond"), uh(a), p(a, "D_friction"), p(a, "D_contact_fiber"),
        p(a, "W_p_fiber"), p(b, "KE_solid"), p(b, "KE_fiber"), p(b, "U_bond"), p(b, "D_friction"),
        p(b, "D_contact_fiber"), p(b, "W_p_fiber")))
