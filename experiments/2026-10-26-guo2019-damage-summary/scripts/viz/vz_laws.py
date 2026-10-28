"""v_z(t) шара для Герца, pinball и Swegle: Guo 1T base, mu = 0.1, V0 = 73 м/с, нить s23."""
import csv, glob, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1]
OUT = sys.argv[2]
runs = [("Герц", "k0_hertz", "#2a78d6"), ("pinball", "k1_pinball", "#e04a3f"), ("Swegle", "k2_swegle", "#008300")]
fig, ax = plt.subplots(figsize=(10, 5.6), facecolor="#fcfcfb")
for name, d, color in runs:
    mon = glob.glob(os.path.join(R, d, "*", "monitor.csv"))[0]
    shot = glob.glob(os.path.join(R, d, "*", "shot.txt"))[0]
    kv = dict(tok.split("=", 1) for tok in open(shot, encoding="utf-8").read().split() if "=" in tok)
    rows = list(csv.DictReader(open(mon)))
    t = [float(r["t"]) * 1e6 for r in rows]
    vz = [float(r["v_z"]) for r in rows]
    lab = "%s: %s при %.0f мкс, V_R = %.1f м/с, разрывы %s/%s" % (
        name, "пробит" if kv.get("outcome") == "perforated" else kv.get("outcome"), float(kv["t_end"]) * 1e6,
        float(kv["V_R"]), kv["n_broken_v"].split(".")[0], kv["n_broken_h"].split(".")[0])
    ax.plot(t, vz, color=color, linewidth=2, label=lab)
    ax.plot([t[-1]], [vz[-1]], "o", color=color)
ax.axhline(0.0, color="#52514e", linewidth=1)
ax.set_xlabel("время, мкс")
ax.set_ylabel("v_z шара, м/с")
ax.set_title("Guo 1T base, μ = 0,1, V0 = 73 м/с, пластичная нить σ_u = 2,3 ГПа: скорость шара по законам контакта нитей",
             fontsize=10)
ax.grid(color="#e4e3df")
ax.legend(loc="lower right", frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(OUT, dpi=130)
print("ok", OUT)
