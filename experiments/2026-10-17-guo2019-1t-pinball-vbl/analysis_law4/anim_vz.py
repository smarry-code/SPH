"""GIF: только график v_z(t) шара (monitor.csv) с бегущим маркером текущего времени; итог выстрела в заголовке.

Usage: python anim_vz.py <run-dir с monitor.csv и shot.txt> <out.gif> [--frames 60] [--law "закон 4 (Герц + демпфер)"]
"""
import argparse, csv, io, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("out")
ap.add_argument("--frames", type=int, default=60)
ap.add_argument("--law", default="закон 4 (Герц + демпфер)")
ap.add_argument("--fps", type=float, default=8.0)
a = ap.parse_args()

kv = dict(t.split("=", 1) for t in open(os.path.join(a.run, "shot.txt"), encoding="utf-8").read().split() if "=" in t)
rows = list(csv.DictReader(open(os.path.join(a.run, "monitor.csv"))))
t = [float(r["t"]) * 1e6 for r in rows]
vz = [float(r["v_z"]) for r in rows]
nb = [int(float(r["n_broken_v"])) + int(float(r["n_broken_h"])) for r in rows]
v0, mu = float(kv["V0"]), float(kv["mu"])
nv, nh = kv["n_broken_v"].split(".")[0], kv["n_broken_h"].split(".")[0]
vr = float(kv["V_R"])
if kv.get("outcome") == "perforated":
    verdict, col = "ПРОБИТ: V_R = %.1f м/с, разрывы %s/%s" % (vr, nv, nh), "#c0392b"
elif kv.get("outcome") == "not_perforated":
    ok = vr == 0.0 and int(nv) <= 8 and int(nh) == 0
    verdict = "шар остановлен, разрывы %s/%s -> %s" % (nv, nh, "НЕ ПРОБИТ" if ok else "пробит по критерию (> 8/0)")
    col = "#1e8449" if ok else "#c0392b"
else:
    verdict, col = "исход не определён", "#52514e"
tb = float(kv["t_first_break"]) * 1e6 if float(kv.get("t_first_break", "-1")) > 0 else None

idx = sorted({round(i * (len(t) - 1) / (a.frames - 1)) for i in range(a.frames)})
ylo = min(min(vz), 0.0) - 5.0
yhi = max(max(vz), 0.0) + 5.0
frames = []
for k in idx:
    fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=90, facecolor="#fcfcfb")
    ax.plot(t, vz, color="#c9c8c3", linewidth=1.5)
    ax.plot(t[:k + 1], vz[:k + 1], color="#2a78d6", linewidth=2.2)
    ax.plot([t[k]], [vz[k]], "o", color="#2a78d6", markersize=8)
    ax.axvline(t[k], color="#2a78d6", linewidth=0.8, alpha=0.5)
    ax.axhline(0.0, color="#52514e", linewidth=1)
    if tb is not None:
        ax.axvline(tb, color="#c98a00", linewidth=1, linestyle="--")
        ax.text(tb, yhi - 2, " 1-й разрыв %.0f мкс" % tb, color="#c98a00", fontsize=8, va="top")
    ax.text(0.98, 0.06, "t = %.0f мкс\nv_z = %.1f м/с\nразрывов: %d" % (t[k], vz[k], nb[k]), transform=ax.transAxes,
            ha="right", va="bottom", fontsize=10, family="monospace")
    ax.set_xlim(0, t[-1] * 1.02)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel("время, мкс")
    ax.set_ylabel("v_z шара, м/с")
    ax.grid(color="#e4e3df")
    fig.suptitle("Guo 1T base, %s, μ = %g, V0 = %g м/с" % (a.law, mu, v0), fontsize=10)
    ax.set_title(verdict, color=col, fontsize=11, fontweight="bold")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    buf.seek(0)
    frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))
dur = [int(1000 / a.fps)] * len(frames)
dur[-1] = 2500
frames[0].save(a.out, save_all=True, append_images=frames[1:], duration=dur, loop=0)
print("GIF", a.out, "%.2f МБ" % (os.path.getsize(a.out) / 1048576), len(frames), "кадров")
