"""Сводка серии закона 3 + выталкивание (tau 0,5 мкс, шаг x0,5, E 80, mu 0,1): v_z(t), число разрывов(t),
энергия шара и импульс в закрепление; отметки статьи. Аргументы: out.png label=dir ..."""
import csv, glob, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = sys.argv[1]
runs = [a.split("=", 1) for a in sys.argv[2:]]
colors = ["#2a78d6", "#d6452a", "#2aa876", "#8a4fd6", "#d69a2a", "#555555", "#c2398a"]
fig, ax = plt.subplots(1, 3, figsize=(18, 5))
for (lab, d), c in zip(runs, colors):
    rd = os.path.dirname(glob.glob(os.path.join(d, "**", "shot.txt"), recursive=True)[0])
    kv = dict(t.split("=", 1) for t in open(os.path.join(rd, "shot.txt"), encoding="utf-8").read().split() if "=" in t)
    m = [q for q in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if q.get("w_p")]
    t = [float(q["t"]) * 1e6 for q in m]
    if kv["outcome"] == "not_perforated":
        res = "НП"
    elif kv["outcome"] == "interrupted":
        res = "оборван на %.0f мкс" % (float(kv["t_end"]) * 1e6)
    elif kv["outcome"] == "undetermined":
        res = "не определён к %.0f мкс" % (float(kv["t_end"]) * 1e6)
    else:
        res = "П, V_R %.1f" % float(kv["V_R"])
    full = "%s: %s, %s/%s" % (lab, res, kv["n_broken_v"], kv["n_broken_h"])
    ax[0].plot(t, [float(q["v_z"]) for q in m], color=c, label=full)
    ax[1].step(t, [int(float(q["n_broken_v"])) + int(float(q["n_broken_h"])) for q in m], where="post", color=c,
               label=lab)
    J, js = 0.0, [0.0]
    for a, b in zip(m, m[1:]):
        J += 0.5 * (float(a["Ff_z"]) + float(b["Ff_z"])) * (float(b["t"]) - float(a["t"]))
        js.append(J)
    ax[2].plot(t, js, color=c, label=lab)
ax[0].axhline(0, color="k", lw=0.6)
TMAX = float(os.environ.get("SUMMARY_TMAX", "1500"))
BX = [float(x) for x in os.environ.get("SUMMARY_BREAK_XLIM", "300,700").split(",")]
TWO_T = os.environ.get("SUMMARY_CFG", "1T") == "2T"
ax[0].set(title="v_z шара", xlabel="t, мкс", ylabel="м/с", xlim=(0, TMAX))
if not TWO_T:
    ax[1].axhline(3, color="k", ls="--", lw=0.8, label="статья: 3 / 0 (73 м/с, не пробит)")
ax[1].set(title="Число порванных нитей", xlabel="t, мкс", ylabel="нитей", xlim=tuple(BX))
ax[2].set(title="Импульс, переданный тканью в закрепление", xlabel="t, мкс", ylabel="Н·с", xlim=(0, TMAX))
for a in ax:
    a.grid(alpha=0.3)
    a.legend(fontsize=7.5)
fig.suptitle(os.environ.get("SUMMARY_TITLE", "1T base, закон 3 + выталкивание τ = 0,5 мкс, шаг ×0,5, E = 80 ГПа, μ = 0,1")
             + (". Статья: V_BL(2T) = 147 м/с, прогиб к первому разрыву 25–27 мм" if TWO_T else
                ". Статья: не пробит до 73 м/с, V_BL = 76 (+3/−3), 3 / 0 разрывов, прогиб к первому разрыву 25–27 мм"))
fig.tight_layout()
fig.savefig(out, dpi=110)
print(out)
