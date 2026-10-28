"""1T base, закон 4, mu = 0.1: E = 97 (70/76 м/с) против E = 78 ГПа (71/76 м/с) -- v_z(t) и каналы энергии."""
import csv, glob, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E97, E78, OUT = sys.argv[1:4]
runs = [("E = 97, 70 м/с", E97, "mu0.1_v70", "#e04a3f", "-"), ("E = 97, 76 м/с", E97, "mu0.1_v76", "#e04a3f", "--"),
        ("E = 78, 71 м/с", E78, "mu0.1_v71", "#2a78d6", "-"), ("E = 78, 76 м/с", E78, "mu0.1_v76", "#2a78d6", "--")]
fig, ax = plt.subplots(1, 3, figsize=(18, 5.6), facecolor="#fcfcfb")
for name, root, d, c, ls in runs:
    sh = glob.glob(os.path.join(root, d, "*", "shot.txt"))[0]
    rd = os.path.dirname(sh)
    kv = dict(t.split("=", 1) for t in open(sh, encoding="utf-8").read().split() if "=" in t)
    m = list(csv.DictReader(open(os.path.join(rd, "monitor.csv"))))
    ax[0].plot([float(q["t"]) * 1e6 for q in m], [float(q["v_z"]) for q in m], color=c, ls=ls, lw=2,
               label="%s: V_R = %.1f, разрывы %s/%s, 1-й разрыв %.0f мкс" % (
                   name, float(kv["V_R"]), kv["n_broken_v"].split(".")[0], kv["n_broken_h"].split(".")[0],
                   float(kv["t_first_break"]) * 1e6))
    e = list(csv.DictReader(open(os.path.join(rd, "energy_log.csv"))))
    E0 = float(e[0]["E_mech"])
    te = [float(q["time"]) * 1e6 for q in e]
    ax[1].plot(te, [100 * float(q["D_friction"]) / E0 for q in e], color=c, ls=ls, lw=2,
               label="%s: %.0f %%" % (name, 100 * float(e[-1]["D_friction"]) / E0))
    ax[2].plot(te, [100 * float(q["D_contact_fiber"]) / E0 for q in e], color=c, ls=ls, lw=2,
               label="%s: %.0f %%" % (name, 100 * float(e[-1]["D_contact_fiber"]) / E0))
ax[0].axhline(0, color="#52514e", lw=1)
ax[0].set_ylabel("v_z шара, м/с")
ax[0].set_title("скорость шара (все пробиты)", fontsize=11)
ax[1].set_title("трение, % E0", fontsize=11)
ax[2].set_title("демпфер контакта (D_contact_fiber), % E0", fontsize=11)
for a in ax:
    a.set_xlabel("время, мкс")
    a.grid(color="#e4e3df")
    a.legend(frameon=False, fontsize=8)
fig.suptitle("Guo 1T base, закон 4 (Герц + демпфер), μ = 0,1, σ_u = 2,5 ГПа: модуль Юнга 97 против 78 ГПа", fontsize=10.5)
fig.tight_layout()
fig.savefig(OUT, dpi=120)
print("ok", OUT)
