"""v_z(t) шара по всем готовым выстрелам серии с законом 4: две панели (70 и 76 м/с), цвет -- mu."""
import csv, glob, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V = sys.argv[1]      # каталог с mu*_v*/<run>/
OUT = sys.argv[2]
COLORS = {0.1: "#2a78d6", 0.3: "#008300", 0.5: "#e04a3f", 0.7: "#8a4fd8", 0.9: "#c98a00"}
fig, axes = plt.subplots(1, 2, figsize=(14, 5.8), facecolor="#fcfcfb", sharey=True)
for ax, v0 in zip(axes, (70, 76)):
    for mu in sorted(COLORS):
        runs = glob.glob(os.path.join(V, "mu%g_v%d" % (mu, v0), "*", "shot.txt"))
        if not runs:
            continue
        d = os.path.dirname(runs[0])
        kv = dict(tok.split("=", 1) for tok in open(runs[0], encoding="utf-8").read().split() if "=" in tok)
        rows = list(csv.DictReader(open(os.path.join(d, "monitor.csv"))))
        t = [float(r["t"]) * 1e6 for r in rows]
        vz = [float(r["v_z"]) for r in rows]
        out = {"perforated": "пробит", "not_perforated": "шар остановлен"}.get(kv.get("outcome"), kv.get("outcome"))
        nv, nh = kv["n_broken_v"].split(".")[0], kv["n_broken_h"].split(".")[0]
        lab = "μ = %g: %s при %.0f мкс, V_R = %.1f м/с, разрывы %s/%s" % (
            mu, out, float(kv["t_end"]) * 1e6, float(kv["V_R"]), nv, nh)
        c = COLORS[mu]
        ax.plot(t, vz, color=c, linewidth=2, label=lab)
        ax.plot([t[-1]], [vz[-1]], "o", color=c)
    ax.axhline(0.0, color="#52514e", linewidth=1)
    ax.set_title("V0 = %d м/с" % v0, fontsize=11)
    ax.set_xlabel("время, мкс")
    ax.grid(color="#e4e3df")
    ax.legend(loc="lower right", frameon=False, fontsize=8.5)
axes[0].set_ylabel("v_z шара, м/с")
fig.suptitle("Guo 1T base, закон 4 (Герц + демпфер), пластичная нить σ_u = 2,5 ГПа, E = 97 ГПа: скорость шара"
             " (не пробит = V_R = 0 и разрывов ≤ 8/0)", fontsize=10)
fig.tight_layout()
fig.savefig(OUT, dpi=130)
print("ok", OUT)
