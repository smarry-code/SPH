#!/usr/bin/env python3
"""Фрагмент сгенерированной «шахматной» ткани по файлам сцены guo2019_shot (--chess=flip|same) -- PDF.

Читает guo2019_fabric.fib (FIBER -- бусины нитей, BOND ... 1 -- стяжки) и guo2019_fabric_layout.txt (слой, направление,
индекс, первая и последняя бусина нити). Подслой нити определяется по её среднему z: s = round((c_z0 - z) / D).

Страницы на каждый режим:
  1. вид сверху на фрагмент: нити по подслоям (основа и уток), стяжки -- точки;
  2. сечение y = 0 (основа j = 21 и утки в пересечениях с ней) и сечение по крайней нити ядра;
  3. 3D-вид фрагмента у угла кусочков.

Запуск (claude-compute-8t: ~/fabric_kT/vizenv/bin/python):
  chess_fib_view.py --scene flip=<каталог выстрела> --scene same=<каталог> --out chess_generated.pdf [--png pg]
"""
import argparse, math, os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d.art3d import Line3DCollection

D_ = 0.39e-3
A_ = 0.6 * D_
DPLY = 2 * A_ + D_ + 1e-6
CZ0 = -(A_ + D_ / 2)
P_ = 1.18e-3
N = 43
SUBCOL = ["#8e44ad", "#2e86c1", "#d35400", "#16a085", "#c0392b", "#7f8c8d"]
NAMES = {"flip": "режим flip: у слоя i + 1 чётность доски обратная",
         "same": "режим same: у всех слоёв одна доска"}


def read_scene(d):
    fib = os.path.join(d, "guo2019_fabric.fib")
    pts, bonds = [], []
    with open(fib, encoding="utf-8") as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        t = lines[i].split()
        if t and t[0] == "FIBER":
            n = int(t[2])
            for k in range(n):
                x, y, z = map(float, lines[i + 1 + k].split())
                pts.append((x, y, z))
            i += n + 1
            continue
        if t and t[0] == "BOND" and len(t) >= 4 and t[3] == "1":
            bonds.append((int(t[1]), int(t[2])))
        i += 1
    X = np.array(pts)
    yarns = []
    with open(os.path.join(d, "guo2019_fabric_layout.txt"), encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            ply, oid, dirn, idx = int(t[0]), int(t[1]), t[2], int(t[3])
            a, b = int(t[6]), int(t[7])
            seg = X[a:b + 1]
            s = int(round((CZ0 - seg[:, 2].mean()) / DPLY))
            yarns.append(dict(ply=ply, dir=dirn, idx=idx, a=a, b=b, sub=s))
    return X, yarns, np.array(bonds, dtype=int).reshape(-1, 2)


def coord(i):
    return (i - (N - 1) / 2) * P_


def plan_page(pdf, X, yarns, bonds, mode, box):
    subs = sorted({y["sub"] for y in yarns})
    fig, axs = plt.subplots(1, len(subs), figsize=(11.7, 4.2))
    for ax, s in zip(axs, subs):
        for y in yarns:
            if y["sub"] != s:
                continue
            seg = X[y["a"]:y["b"] + 1] * 1e3
            ax.plot(seg[:, 0], seg[:, 1], lw=0.9, color=SUBCOL[s % len(SUBCOL)] if y["dir"] == "x" else "#f39c12")
        # стяжки, верхний конец которых на подслое s (нижняя нить верхнего кусочка: z = c_z,s - A)
        sel = bonds[(np.floor((CZ0 + A_ + 1e-9 - X[bonds[:, 0], 2]) / DPLY).astype(int) == s)]
        ax.plot(X[sel[:, 0], 0] * 1e3, X[sel[:, 0], 1] * 1e3, "o", ms=2.2, color="k", zorder=5)
        ax.set_xlim(-box, box)
        ax.set_ylim(-box, box)
        ax.set_aspect("equal")
        ax.set_title("подслой %d (слой %d, %s)" % (s + 1, s // 2 + 1, "i − h" if s % 2 == 0 else "i + h"), fontsize=9)
        ax.set_xlabel("x, мм")
        ax.set_ylabel("y, мм")
    fig.suptitle("Сгенерированная ткань 2T, %s.\nВид сверху на фрагмент |x|, |y| ≤ %g мм: основа — цвет подслоя, "
                 "уток — оранжевый, точки — стяжки (верхний конец)" % (NAMES[mode], box), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    pdf.savefig(fig)
    plt.close(fig)


def tube(x, z, r):
    P = np.stack([x, z], 1)
    t = np.gradient(P, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    nrm = np.stack([-t[:, 1], t[:, 0]], 1)
    poly = np.vstack([P + r * nrm, (P - r * nrm)[::-1]])
    return poly[:, 0], poly[:, 1]


def section(ax, X, yarns, bonds, j, xlim):
    yj = coord(j)
    for y in yarns:
        seg = X[y["a"]:y["b"] + 1]
        if y["dir"] == "x" and y["idx"] == j:
            ax.fill(*tube(seg[:, 0] * 1e3, seg[:, 2] * 1e3, D_ / 2 * 1e3), fc=SUBCOL[y["sub"] % len(SUBCOL)],
                    ec="k", lw=0.3, zorder=2)
        elif y["dir"] == "y":
            k = np.argmin(np.abs(seg[:, 1] - yj))
            if abs(seg[k, 1] - yj) < 1e-7:
                ax.add_patch(Circle((seg[k, 0] * 1e3, seg[k, 2] * 1e3), D_ / 2 * 1e3, fc="#f5b041", ec="k", lw=0.3,
                                    zorder=4))
    for a, b in bonds:
        if abs(X[a, 1] - yj) < 1e-7:
            ax.plot([X[a, 0] * 1e3] * 2, [X[a, 2] * 1e3, X[b, 2] * 1e3], color="#e74c3c", lw=1.8, zorder=5)
    ax.set_xlim(*xlim)
    ax.set_aspect("equal")
    ax.set_xlabel("x, мм")
    ax.set_ylabel("z, мм")


def section_page(pdf, X, yarns, bonds, mode):
    fig, axs = plt.subplots(2, 1, figsize=(11.7, 7.5))
    for ax, j, lab in ((axs[0], 21, "y = 0 (основа j = 21)"), (axs[1], 23, "y = %.2f мм (основа j = 23, крайняя нить ядра)"
                                                                       % (coord(23) * 1e3))):
        section(ax, X, yarns, bonds, j, (-9.5, 9.5))
        ax.set_title("Сечение %s" % lab, fontsize=10)
    fig.suptitle("Сгенерированная ткань 2T, %s.\nОснова — цвет подслоя (1 фиолетовый, 2 синий, 3 оранжевый, "
                 "4 зелёный), утки — круги, стяжки — красные (по файлу .fib)" % NAMES[mode], fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def view3d_page(pdf, X, yarns, bonds, mode, lo, hi):
    fig = plt.figure(figsize=(11.7, 8.0))
    for n, (elev, azim) in enumerate(((28, -60), (12, -20))):
        ax = fig.add_subplot(1, 2, n + 1, projection="3d")
        for y in yarns:
            seg = X[y["a"]:y["b"] + 1] * 1e3
            m = (seg[:, 0] >= lo[0]) & (seg[:, 0] <= hi[0]) & (seg[:, 1] >= lo[1]) & (seg[:, 1] <= hi[1])
            if m.sum() < 2:
                continue
            idx = np.where(m)[0]
            seg = seg[idx.min():idx.max() + 1]
            col = SUBCOL[y["sub"] % len(SUBCOL)]
            ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color=col, lw=2.2 if y["dir"] == "x" else 1.4,
                    alpha=1.0 if y["dir"] == "x" else 0.75)
        segs = []
        for a, b in bonds:
            p, q = X[a] * 1e3, X[b] * 1e3
            if lo[0] <= p[0] <= hi[0] and lo[1] <= p[1] <= hi[1]:
                segs.append([p, q])
        if segs:
            ax.add_collection3d(Line3DCollection(segs, colors="#e74c3c", linewidths=2.0))
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(-3.6, 0.3)
        ax.set_box_aspect((hi[0] - lo[0], hi[1] - lo[1], 3.9 * 1.6))
        ax.view_init(elev=elev, azim=azim)
        ax.set_xlabel("x, мм")
        ax.set_ylabel("y, мм")
        ax.set_zlabel("z, мм")
    fig.suptitle("Сгенерированная ткань 2T, %s: фрагмент у угла четырёх кусочков x ∈ [%g; %g], y ∈ [%g; %g] мм.\n"
                 "Цвет — подслой (1 фиолетовый, 2 синий, 3 оранжевый, 4 зелёный); толстые — основа, тонкие — уток; "
                 "красные — стяжки" % (NAMES[mode], lo[0], hi[0], lo[1], hi[1]), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    pdf.savefig(fig)
    plt.close(fig)


class PdfPng:
    def __init__(self, pdf, png):
        self.pdf, self.png, self.n = pdf, png, 0

    def savefig(self, fig):
        self.pdf.savefig(fig)
        if self.png:
            self.n += 1
            fig.savefig("%s_%d.png" % (self.png, self.n), dpi=110)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene", action="append", required=True, help="режим=каталог выстрела")
    ap.add_argument("--out", default="chess_generated.pdf")
    ap.add_argument("--box", type=float, default=9.0, help="полуширина фрагмента на виде сверху, мм")
    ap.add_argument("--png", default="")
    a = ap.parse_args()
    with PdfPages(a.out) as pdf0:
        pdf = PdfPng(pdf0, a.png)
        for item in a.scene:
            mode, d = item.split("=", 1)
            X, yarns, bonds = read_scene(d)
            subs = sorted({y["sub"] for y in yarns})
            print(mode, "бусин", len(X), "нитей", len(yarns), "стяжек", len(bonds), "подслои", subs)
            plan_page(pdf, X, yarns, bonds, mode, a.box)
            section_page(pdf, X, yarns, bonds, mode)
            view3d_page(pdf, X, yarns, bonds, mode, (0.6, 0.6), (6.5, 6.5))
    print("записан", a.out)


if __name__ == "__main__":
    main()
