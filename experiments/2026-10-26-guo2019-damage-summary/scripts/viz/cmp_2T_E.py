"""2T base, закон 4, mu = 0.1, V0 = 144: E = 97 против E = 78 ГПа -- v_z(t), разрывы по слоям, каналы энергии."""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1]   # scratchpad/local2T
OUT = sys.argv[2]
runs = [("E = 97 ГПа", "mu0.1_v144", "#e04a3f"), ("E = 78 ГПа", "mu0.1_v144_E78", "#2a78d6")]
fig, ax = plt.subplots(1, 3, figsize=(18, 5.6), facecolor="#fcfcfb")
summary = []
for name, d, c in runs:
    rd = os.path.join(D, d, "2T_base_v144.0_r0")
    kv = dict(t.split("=", 1) for t in open(os.path.join(rd, "shot.txt"), encoding="utf-8").read().split() if "=" in t)
    m = [q for q in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if q.get("w_f_p2")]
    t = [float(q["t"]) * 1e6 for q in m]
    ax[0].plot(t, [float(q["v_z"]) for q in m], color=c, lw=2,
               label="%s: первый разрыв %.0f мкс, прогиб %.1f мм; в конце v_z = %.1f м/с" % (
                   name, float(kv["t_first_break"]) * 1e6, float(kv["w_f_at_first_break"]) * 1e3, float(m[-1]["v_z"])))
    ax[1].plot(t, [int(q["n_broken_v_p1"]) for q in m], color=c, lw=2, label="%s, слой 1" % name)
    ax[1].plot(t, [int(q["n_broken_v_p2"]) for q in m], color=c, lw=2, ls="--", label="%s, слой 2" % name)
    e = list(csv.DictReader(open(os.path.join(rd, "energy_log.csv"))))
    E0 = float(e[0]["E_mech"])
    te = [float(q["time"]) * 1e6 for q in e]
    for k, lab, ls in (("KE_solid", "KE шара", "-"), ("D_friction", "трение", "--"),
                       ("D_contact_fiber", "демпфер", ":"), ("KE_fiber", "KE ткани", "-.")):
        ax[2].plot(te, [100 * float(q[k]) / E0 for q in e], color=c, ls=ls, lw=1.8, label="%s: %s" % (name, lab))
ax[0].set_title("скорость шара", fontsize=11)
ax[0].set_ylabel("v_z, м/с")
ax[0].axhline(0, color="#52514e", lw=1)
ax[1].set_title("порвано вертикальных нитей (горизонтальных 0)", fontsize=11)
ax[1].set_ylabel("нитей")
ax[2].set_title("каналы энергии, % E0 (E0 = 72,6 Дж)", fontsize=11)
ax[2].set_ylabel("% E0")
for a in ax:
    a.set_xlabel("время, мкс")
    a.grid(color="#e4e3df")
    a.legend(frameon=False, fontsize=8)
fig.suptitle("Guo 2T base, закон 4 (Герц + демпфер), μ = 0,1, V0 = 144 м/с (статья: V_BL = 147, не пробит): "
             "влияние модуля Юнга нити (σ_u = 2,5 ГПа)", fontsize=10.5)
fig.tight_layout()
fig.savefig(OUT, dpi=120)
print("ok", OUT)
