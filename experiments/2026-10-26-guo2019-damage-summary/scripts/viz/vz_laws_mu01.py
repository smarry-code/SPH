"""v_z(t) шара: pinball / Swegle / закон 4 при mu = 0.1, V0 = 70 и 76 (две панели). Нить s25, E = 97 ГПа."""
import csv, glob, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PIN, SWE, L4, OUT = sys.argv[1:5]  # каталоги с mu0.1_v70/<run>/, mu0.1_v76/<run>/
laws = [("pinball (1)", PIN, "#e04a3f"), ("Swegle (2)", SWE, "#008300"), ("Герц + демпфер (4)", L4, "#2a78d6")]
fig, axes = plt.subplots(1, 2, figsize=(14, 5.8), facecolor="#fcfcfb", sharey=True)
for ax, v0 in zip(axes, (70, 76)):
    for name, root, c in laws:
        sh = glob.glob(os.path.join(root, "mu0.1_v%d" % v0, "*", "shot.txt"))
        if not sh:
            continue
        d = os.path.dirname(sh[0])
        kv = dict(t.split("=", 1) for t in open(sh[0], encoding="utf-8").read().split() if "=" in t)
        rows = list(csv.DictReader(open(os.path.join(d, "monitor.csv"))))
        t = [float(r["t"]) * 1e6 for r in rows]
        vz = [float(r["v_z"]) for r in rows]
        lab = "%s: V_R = %.1f м/с, разрывы %s/%s, δ/d нить–нить %.2f" % (
            name, float(kv["V_R"]), kv["n_broken_v"].split(".")[0], kv["n_broken_h"].split(".")[0],
            float(kv["max_delta_ff/d"]))
        ax.plot(t, vz, color=c, linewidth=2, label=lab)
        ax.plot([t[-1]], [vz[-1]], "o", color=c)
    ax.axhline(0.0, color="#52514e", linewidth=1)
    ax.set_title("V0 = %d м/с (все пробиты)" % v0, fontsize=11)
    ax.set_xlabel("время, мкс")
    ax.grid(color="#e4e3df")
    ax.legend(loc="lower right", frameon=False, fontsize=8.5)
axes[0].set_ylabel("v_z шара, м/с")
fig.suptitle("Guo 1T base, μ = 0,1, пластичная нить σ_u = 2,5 ГПа, E = 97 ГПа: скорость шара по законам контакта нитей",
             fontsize=10)
fig.tight_layout()
fig.savefig(OUT, dpi=130)
print("ok", OUT)
