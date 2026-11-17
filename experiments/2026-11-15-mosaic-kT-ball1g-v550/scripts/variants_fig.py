#!/usr/bin/env python3
"""Картинки вариантов сшивки A и B (16T) по реальной сцене: cut.json + начальный кадр нитей.

Для каждого варианта -- один PNG:
  слева  -- вид сверху: кусочки слоя 1 (заливка), швы слоёв 1..3 (цветные линии), точки сквозных стежков 1-2;
  справа -- сечение x-z полосой |y - y0| < p/2 около точки удара: бусины по кусочкам, стежки в слое (зелёные),
            сквозные стежки (чёрные отрезки между соседними слоями).
Запуск: python scripts/variants_fig.py   (из каталога эксперимента) -> fig/variant_A.png, fig/variant_B.png
"""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "fig")
SEAM_COL = ["#D55E00", "#0072B2", "#009E73"]


def load(var):
    d = glob.glob(os.path.join(HERE, "runs", "16T_%s" % var, "16T_*"))[0]
    j = json.load(open(os.path.join(d, "guo2019_fabric.cut.json"), encoding="utf-8"))
    n = j["n_beads"]
    pos = np.zeros((n, 3))
    with open(os.path.join(d, "fiber_output_00001.txt")) as f:
        for ln in f:
            if ln.strip() == "PARTICLES":
                break
        for i in range(n):
            s = ln = f.readline().split()
            pos[int(s[0])] = (float(s[3]), float(s[4]), float(s[5]))
    frag = -np.ones(n, int)
    ply = -np.ones(n, int)
    for r in j["ranges"]:
        frag[r["bead_first"]:r["bead_first"] + r["bead_count"]] = r["frag"]
        ply[r["bead_first"]:r["bead_first"] + r["bead_count"]] = r["ply"]
    cols = j["stitch_columns"]
    st = np.array(j["stitches"], dtype=np.int64).reshape(-1, len(cols))
    return j, pos, frag, ply, st, cols


def palette(n):
    base = [plt.get_cmap("tab10")(i) for i in range(10)]
    return [base[i % 10] for i in range(n)]


def figure(var, title):
    j, pos, frag, ply, st, cols = load(var)
    a, b = st[:, cols.index("a")], st[:, cols.index("b")]
    thr = st[:, cols.index("through")].astype(bool)
    pa, pb = st[:, cols.index("ply_a")], st[:, cols.index("ply_b")]
    pal = palette(j["n_frag"])
    p = j["p"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(15, 6.6), gridspec_kw={"width_ratios": [1, 1.35]})
    # --- вид сверху
    cut1 = j["cuts"][j["plies"][0]["cut"]]
    for fr in cut1["fragments"]:
        ax.add_patch(Polygon(np.array(fr["polygon"]) * 1e3, closed=True, fc=pal[fr["id"]], alpha=0.30, ec="none"))
        c = np.array(fr["center"]) * 1e3
        ax.text(c[0], c[1], str(fr["id"] + 1), ha="center", va="center", fontsize=11, color="0.25")
    ncut = len(j["cuts"])
    for k in range(min(3, ncut)):
        cut = j["cuts"][j["plies"][k]["cut"]]
        for s in cut["seams"]:
            q = np.array(s["points"]) * 1e3
            ax.plot(q[:, 0], q[:, 1], "-", color=SEAM_COL[k], lw=2.2 - 0.4 * k,
                    label=("швы слоя %d" % (k + 1)) if s is cut["seams"][0] else None)
    m = thr & (pa == 0)
    ax.plot(pos[a[m], 0] * 1e3, pos[a[m], 1] * 1e3, "k.", ms=2.5, label="сквозные стежки слой 1 – слой 2")
    ax.axhspan(25, 35, color="0.6", alpha=0.35)
    ax.axhspan(-35, -25, color="0.6", alpha=0.35)
    ax.text(0, 30, "зажим", ha="center", va="center", fontsize=10)
    ax.text(0, -30, "зажим", ha="center", va="center", fontsize=10)
    ax.plot([0], [0], marker="+", color="k", ms=16, mew=2)
    ax.add_patch(plt.Circle((0, 0), 6.253 / 2, fill=False, ls="--", color="k"))
    y0 = 0.0
    ax.axhline(y0, color="m", lw=0.8, ls=":")
    ax.set_xlim(-25, 25); ax.set_ylim(-35, 35); ax.set_aspect("equal")
    ax.set_xlabel("x, мм"); ax.set_ylabel("y, мм")
    ax.set_title("Вид сверху: кусочки слоя 1%s" % ("" if ncut == 1 else "; швы слоёв 1–3 не совпадают"), fontsize=10)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    # --- сечение x-z
    X0, X1 = -12.0, 12.0
    sel = (np.abs(pos[:, 1] - y0) < 0.5 * p) & (pos[:, 0] * 1e3 > X0) & (pos[:, 0] * 1e3 < X1)
    zl = np.array([pl["cz"] for pl in j["plies"]])
    for k in range(len(zl)):
        mk = sel & (ply == k)
        bx.scatter(pos[mk, 0] * 1e3, pos[mk, 2] * 1e3, s=4, c=[pal[f] for f in frag[mk]], lw=0)
    win = lambda i: (np.abs(pos[i, 1] - y0) < 1.5 * p) & (pos[i, 0] * 1e3 > X0) & (pos[i, 0] * 1e3 < X1)
    ms = win(a) & win(b)
    for i in np.where(ms & thr)[0]:
        bx.plot(pos[[a[i], b[i]], 0] * 1e3, pos[[a[i], b[i]], 2] * 1e3, "k-", lw=0.9)
    for i in np.where(ms & ~thr)[0]:
        bx.plot(pos[[a[i], b[i]], 0] * 1e3, pos[[a[i], b[i]], 2] * 1e3, "-", color="#009E73", lw=2.0)
    bx.plot([], [], "k-", label="сквозной стежок (соседние слои)")
    bx.plot([], [], "-", color="#009E73", lw=2, label="стежок в слое (через шов)")
    bx.add_patch(plt.Circle((0, zl[0] * 1e3 + 0.2 + 3.13), 3.13, fill=False, ls="--", color="k"))
    bx.set_xlim(X0, X1); bx.set_ylim(zl[-1] * 1e3 - 1.0, zl[0] * 1e3 + 7.0); bx.set_aspect("equal")
    bx.set_xlabel("x, мм"); bx.set_ylabel("z, мм")
    bx.set_title("Сечение y = 0 (|Δy| < %.1f мм): цвет — номер кусочка" % (1.5 * p * 1e3), fontsize=10)
    bx.legend(loc="upper right", fontsize=8, framealpha=0.9)
    n_in, n_th = int((~thr).sum()), int(thr.sum())
    fig.suptitle("%s   (16T: стежков в слоях %d, сквозных %d)" % (title, n_in, n_th), fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, "variant_%s.png" % var)
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    figure("A", "Вариант A: одна раскройка на весь пакет, кусочек — стопка из 16 слоёв")
    figure("B", "Вариант B: у каждого слоя своя раскройка; сквозные стежки по швам всех раскроек")
