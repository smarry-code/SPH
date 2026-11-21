#!/usr/bin/env python3
"""Сечение «шахматной» ткани kT с перекрытием кусочков (малое полотно 43 нити, окно 50 мм) -- PDF.

Ткань (решение владельца 2026-11-21):
  * слой i разбит на два подслоя: i - h (к лицу) и i + h (к тылу);
  * кусочки -- шахматная доска: ядро P x P нитей (P = 5) по линиям шва centered_patch_seams, кусочек 7 x 7 --
    ядро плюс по одной нити с каждой стороны, заходящей на соседний кусочек (перекрытие);
  * кусочки цвета c = +1 слоя i лежат на подслое i - h, цвета -1 -- на i + h; два варианта: у слоя i + 1
    чётность доски обратная (c = +1 -- на (i+1) + h, --variants flip) или та же (same);
  * расстояние между подслоями i + h и (i+1) - h -- шаг слоёв исходного эксперимента D = 2A + d + 1 мкм.

Что выбрано здесь (геометрия, а не решение владельца):
  * 2h = D: в перекрытии пересечение края одного кусочка лежит над пересечением ядра соседа; меньше D нельзя --
    нити войдут друг в друга (нижняя нить верхнего кусочка на z_u - A, верхняя нить нижнего на z_l + A,
    зазор 2h - 2A >= d). Тогда все подслои идут с шагом D, слой занимает 2D, пакет kT вдвое толще сплошного;
  * стяжка -- в каждом пересечении перекрытия: от нижней нити верхнего кусочка к верхней нити нижнего
    (длина 2h - 2A = d + 1 мкм);
  * углы: где сходятся четыре кусочка, квадрат 2 x 2 пересечений входит в расширения обоих диагональных соседей
    одного подслоя -- их нити совпали бы. Пересечение такого квадрата принадлежит только кусочку, в ядро
    которого оно входит; полосы перекрытия поэтому на 2 пересечения короче ядра (3 из 5 у внутренних).

Геометрия нитей -- как у генератора (include/fiber/plain_weave.hpp): p = 1,18 мм, d = 0,39 мм, A = 0,6 d,
seg = 2, профиль cos; концы нитей кусочка -- на полушаге за крайним пересечением.

Запуск (claude-compute-8t: ~/fabric_kT/vizenv/bin/python):
  chess_section.py --out ../fig/chess_section.pdf [--k 20] [--patch 5] [--variants flip,same] [--png pg]
"""
import argparse, math

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Rectangle

# ---- параметры сцены (мм)
P_ = 1.18
D_ = 0.39
A_ = 0.6 * D_
EPS_PLY = 1e-3
DPLY = 2 * A_ + D_ + EPS_PLY  # 0,859 -- шаг слоёв исходного эксперимента
H2 = DPLY  # 2h -- расстояние между подслоями одного слоя
N = 43
CE = 25.0
BALL_D = 6.253
BALL_GAP = 0.5
R_EFF = BALL_D / 2


def seams_idx(n, patch):
    """Линии шва в индексах нитей (полуцелые) -- как centered_patch_seams."""
    base2 = (n - 1) + (patch if (n - patch) % 2 == 0 else 0)
    j0 = 0
    while base2 + 2 * (j0 - 1) * patch > 0:
        j0 -= 1
    out, j = [], j0
    while True:
        s2 = base2 + 2 * j * patch
        j += 1
        if s2 <= 0:
            continue
        if s2 >= 2 * (n - 1):
            break
        out.append(0.5 * s2)
    return np.array(out)


def coord(u):
    return (np.asarray(u, dtype=float) - (N - 1) / 2) * P_


class Board:
    """Кусочки по одной оси: ядра [lo, hi] (индексы нитей), расширения на одну нить, нити у швов."""

    def __init__(self, seams):
        self.seams = seams
        edges = [-1] + [int(s - 0.5) for s in seams] + [N - 1]
        self.core = [(edges[m] + 1, edges[m + 1]) for m in range(len(edges) - 1)]
        self.M = len(self.core)
        self.mc = self.patch_of((N - 1) // 2)
        self.adj = set()  # нити по обе стороны шва (из них -- угловые квадраты 2 x 2)
        for s in seams:
            self.adj.update((int(s - 0.5), int(s + 0.5)))

    def patch_of(self, i):
        for m, (lo, hi) in enumerate(self.core):
            if lo <= i <= hi:
                return m
        raise ValueError(i)

    def ext(self, m):
        lo, hi = self.core[m]
        return (lo - 1 if m > 0 else lo, hi + 1 if m < self.M - 1 else hi)

    def color(self, mx, my):
        return 1 if (mx + my - 2 * self.mc) % 2 == 0 else -1


def owns(B, mx, my, ix, iy):
    """Входит ли пересечение (ix, iy) в кусочек (mx, my)."""
    cx, cy = B.core[mx], B.core[my]
    ex, ey = B.ext(mx), B.ext(my)
    in_cx, in_cy = cx[0] <= ix <= cx[1], cy[0] <= iy <= cy[1]
    in_ex, in_ey = ex[0] <= ix <= ex[1], ey[0] <= iy <= ey[1]
    if in_cx and in_cy:
        return True
    if not (in_ex and in_ey):
        return False
    if ix in B.adj and iy in B.adj:  # угловой квадрат: только кусочку, в ядро которого входит
        return False
    return (in_cx and in_ey) or (in_ex and in_cy)


FLIP = True  # True -- у слоя i + 1 чётность доски обратная; False -- у всех слоёв одна доска
VLABEL = ""


def sublayer(B, k, mx, my):
    """Номер подслоя (0 -- лицевой) кусочка (mx, my) слоя k (с нуля)."""
    c = B.color(mx, my) * (1 if (k % 2 == 0 or not FLIP) else -1)
    return 2 * k + (0 if c > 0 else 1)


def zc(s):
    return -(A_ + D_ / 2) - s * H2


def crimp(u, sig):
    z = sig * A_ * np.cos(math.pi * u)
    return np.where(np.abs(np.abs(u - np.floor(u)) - 0.5) < 1e-9, 0.0, z)


def piece_yarns(B, k, mx, my):
    """Нити кусочка: список (dir, idx, x, y, z). Основа -- вдоль x на y_idx, уток -- вдоль y на x_idx."""
    out = []
    s = sublayer(B, k, mx, my)
    ex, ey = B.ext(mx), B.ext(my)
    for dirn, rng_cross, rng_along in (("x", ey, ex), ("y", ex, ey)):
        for j in range(rng_cross[0], rng_cross[1] + 1):
            cr = [i for i in range(rng_along[0], rng_along[1] + 1)
                  if (owns(B, mx, my, i, j) if dirn == "x" else owns(B, mx, my, j, i))]
            if not cr:
                continue
            a, b = min(cr), max(cr)
            u = np.arange(2 * a - 1, 2 * b + 2) / 2.0
            sig = (1.0 if j % 2 == 0 else -1.0) * (1.0 if dirn == "x" else -1.0)
            z = zc(s) + crimp(u, sig)
            if dirn == "x":
                out.append((dirn, j, coord(u), np.full_like(u, coord(j)), z))
            else:
                out.append((dirn, j, np.full_like(u, coord(j)), coord(u), z))
    return out


def crossing_z(B, k, mx, my, ix, iy):
    """z основы и утка кусочка в пересечении (ix, iy)."""
    s = sublayer(B, k, mx, my)
    zw = zc(s) + (1.0 if iy % 2 == 0 else -1.0) * A_ * math.cos(math.pi * ix)
    zf = zc(s) - (1.0 if ix % 2 == 0 else -1.0) * A_ * math.cos(math.pi * iy)
    return zw, zf


def owners(B, ix, iy):
    mx0, my0 = B.patch_of(ix), B.patch_of(iy)
    res = []
    for mx in range(max(0, mx0 - 1), min(B.M, mx0 + 2)):
        for my in range(max(0, my0 - 1), min(B.M, my0 + 2)):
            if owns(B, mx, my, ix, iy):
                res.append((mx, my))
    return res


def stitches(B, k):
    """Стяжки слоя k: (ix, iy, z_верх, z_низ) -- от нижней нити верхнего кусочка к верхней нити нижнего."""
    out = []
    for ix in range(N):
        for iy in range(N):
            ow = owners(B, ix, iy)
            if len(ow) < 2:
                continue
            if len(ow) > 2 or len({sublayer(B, k, *m) for m in ow}) != 2:
                raise RuntimeError("пересечение (%d, %d): кусочки %s" % (ix, iy, ow))
            ow.sort(key=lambda m: sublayer(B, k, *m))
            zu = min(crossing_z(B, k, *ow[0], ix, iy))
            zl = max(crossing_z(B, k, *ow[1], ix, iy))
            out.append((ix, iy, zu, zl))
    return out


# ------------------------------------------------------------------ проверка проникания
def sample(xs, ys, zs, m=8):
    t = np.linspace(0, 1, m, endpoint=False)[:, None]
    P = np.stack([xs, ys, zs], 1)
    a, b = P[:-1], P[1:]
    pts = (a[None] * (1 - t[..., None]) + b[None] * t[..., None]).reshape(-1, 3)
    return np.vstack([pts, P[-1:]])


def min_dist(Pa, Pb, chunk=500):
    best = np.inf
    for s in range(0, len(Pa), chunk):
        A = Pa[s:s + chunk]
        d2 = ((A[:, None, :] - Pb[None, :, :]) ** 2).sum(-1)
        best = min(best, d2.min())
    return math.sqrt(best)


def check(B, half=9.0, layers=2):
    """Наименьшее расстояние осей: основа-уток одного кусочка; нити разных кусочков (любых слоёв)."""
    groups = []
    for k in range(layers):
        for mx in range(B.M):
            for my in range(B.M):
                for dirn in "xy":
                    arr = []
                    for d_, j, x, y, z in piece_yarns(B, k, mx, my):
                        if d_ != dirn:
                            continue
                        m = (np.abs(x) <= half + P_) & (np.abs(y) <= half + P_)
                        if m.sum() >= 2:
                            arr.append(sample(x[m], y[m], z[m]))
                    if arr:
                        groups.append(((k, mx, my), dirn, np.vstack(arr)))
    same, diff = np.inf, np.inf
    for a in range(len(groups)):
        for b in range(a + 1, len(groups)):
            pa, da, Pa = groups[a]
            pb, db, Pb = groups[b]
            if pa == pb and da == db:
                continue
            lo = np.maximum(Pa.min(0), Pb.min(0)) - 1.0
            hi = np.minimum(Pa.max(0), Pb.max(0)) + 1.0
            if np.any(lo > hi):
                continue
            d = min_dist(Pa, Pb)
            if pa == pb:
                same = min(same, d)
            else:
                diff = min(diff, d)
    return {"в кусочке": same, "между кусочками": diff}


# ------------------------------------------------------------------ рисунки
COL = {1: "#c0392b", -1: "#2e86c1"}


def plan_page(pdf, B, patch):
    fig, axs = plt.subplots(1, 2, figsize=(11.7, 6.8))
    for ax, k in zip(axs, (0, 1)):
        for mx in range(B.M):
            for my in range(B.M):
                c = 1 if sublayer(B, k, mx, my) % 2 == 0 else -1
                for ix in range(N):
                    for iy in range(N):
                        if owns(B, mx, my, ix, iy):
                            ax.add_patch(Rectangle((coord(ix) - P_ / 2, coord(iy) - P_ / 2), P_, P_, fc=COL[c],
                                                   alpha=0.28, lw=0))
                cx, cy = B.core[mx], B.core[my]
                ax.add_patch(Rectangle((coord(cx[0]) - P_ / 2, coord(cy[0]) - P_ / 2),
                                       (cx[1] - cx[0] + 1) * P_, (cy[1] - cy[0] + 1) * P_, fc="none", ec="k", lw=0.5))
        st = stitches(B, k)
        ax.plot([coord(a) for a, b, _, _ in st], [coord(b) for a, b, _, _ in st], ".", color="k", ms=2.0)
        for yc in (-CE, CE):
            ax.axhline(yc, color="k", lw=1.0, ls="--")
        ax.add_patch(Circle((0, 0), BALL_D / 2, fc="none", ec="k", lw=1.3))
        for yy, lab in ((0.0, "A"), (coord((N - 1) // 2 + 2), "B")):
            ax.plot([-26, 26], [yy, yy], color="0.2", lw=0.8)
            ax.text(26.3, yy, lab, va="center", fontsize=10, weight="bold")
        ax.set_aspect("equal")
        ax.set_xlim(-27, 28)
        ax.set_ylim(-27, 27)
        ax.set_xlabel("x, мм")
        ax.set_ylabel("y, мм")
        ax.set_title(("Слои 1, 3, 5, … " if k == 0 else "Слои 2, 4, 6, … (чётность доски обратная)") if FLIP else
                     ("Слой 1" if k == 0 else "Слой 2 (и все остальные): доска та же"), fontsize=10)
    fig.suptitle("Кусочки 7×7 (ядро %d×%d в рамке), шахматный порядок: красный — подслой i − h (к лицу), "
                 "синий — i + h (к тылу).\nТёмные полосы — перекрытие на одну нить (два кусочка друг над другом), "
                 "точки — стяжки;\nугловые квадраты 2×2 не перекрываются. Пунктир — зажим |y| ≥ 25 мм; "
                 "A, B — сечения" % (patch, patch), fontsize=9)
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


def section_axes(ax, B, j, k_list, xlim, show_ball=True, label=True):
    """Сечение y = y_j: основа j -- трубки диаметра d (цвет подслоя), утки -- круги, стяжки -- красные."""
    for k in k_list:
        for mx in range(B.M):
            for my in range(B.M):
                s = sublayer(B, k, mx, my)
                fc = "#a569bd" if s % 2 == 0 else "#5dade2"
                for d_, jj, x, y, z in piece_yarns(B, k, mx, my):
                    if d_ == "x" and jj == j:
                        ax.fill(*tube(x, z, D_ / 2), fc=fc, ec="k", lw=0.3, zorder=2)
                for ix in range(N):
                    if owns(B, mx, my, ix, j):
                        zw, zf = crossing_z(B, k, mx, my, ix, j)
                        ax.add_patch(Circle((float(coord(ix)), zf), D_ / 2, fc="#f5b041", ec="k", lw=0.3, zorder=4))
        for ix, iy, zu, zl in stitches(B, k):
            if iy == j:
                ax.plot([coord(ix)] * 2, [zu, zl], color="#e74c3c", lw=1.8, zorder=5, solid_capstyle="butt")
        if label and xlim[1] - xlim[0] < 30:
            ax.text(xlim[0] + 0.1, zc(2 * k) + 0.3, "слой %d" % (k + 1), fontsize=7, va="bottom",
                    bbox=dict(fc="w", ec="none", alpha=0.7, pad=0.5), zorder=6)
    for sv in B.seams:
        ax.axvline(coord(sv), color="0.4", lw=0.5, ls=":", zorder=1)
    if show_ball:
        ax.add_patch(Circle((0, R_EFF + BALL_GAP), BALL_D / 2, fc="0.85", ec="k", lw=0.8, zorder=1))
    ax.set_xlim(*xlim)
    ax.set_aspect("equal")
    ax.set_xlabel("x, мм")
    ax.set_ylabel("z, мм")


def zb(k_total):
    return zc(2 * k_total - 1) - A_ - D_ / 2, -(k_total * (2 * A_ + D_) + (k_total - 1) * EPS_PLY)


def footer(chk, k_total, n_st):
    zb_new, zb_old = zb(k_total)
    return ("Подслои через D = %.3f мм (2h = D), слой занимает 2D; тыл %dT z_b = %.2f мм (сплошная ткань %.2f мм). "
            "Стяжек на слой %d, длина %.3f мм.\nНаименьшее расстояние осей нитей (|x|,|y| ≤ 9 мм, слои 1–2), "
            "d = %.2f мм: основа–уток одного кусочка %.3f мм; нити разных кусочков %.3f мм."
            % (H2, k_total, zb_new, zb_old, n_st, H2 - 2 * A_, D_, chk["в кусочке"], chk["между кусочками"]))


def section_pages(pdf, B, patch, k_total, chk, n_st):
    j0 = (N - 1) // 2
    jB = j0 + 2  # крайняя нить ядра центрального кусочка (у шва по y): угловые квадраты
    w = max(8.0, 1.6 * patch * P_)
    fig, axs = plt.subplots(2, 1, figsize=(11.7, 8.3))
    for ax, j, lab in ((axs[0], j0, "A (y = 0, середина ряда кусочков)"),
                       (axs[1], jB, "B (y = %.2f мм, крайняя нить ядра: у швов угловые квадраты)" % coord(jB))):
        section_axes(ax, B, j, range(2), (-w, w))
        ax.set_ylim(zc(3) - A_ - D_, BALL_D + 1.0)
        ax.set_title("Сечение %s, слои 1–2" % lab, fontsize=10)
    fig.suptitle("Кусочки 7×7, ядро %d×%d, натуральный масштаб. Основа: фиолетовая — подслой i − h, голубая — "
                 "i + h;\nутки оранжевые; стяжки красные; пунктир — линии шва ядер" % (patch, patch), fontsize=10)
    fig.text(0.01, 0.005, footer(chk, k_total, n_st), fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)
    # один шов крупно
    sv = float(coord(B.seams[B.seams > (N - 1) / 2][0]))
    fig, ax = plt.subplots(figsize=(11.7, 7.0))
    section_axes(ax, B, j0, range(2), (sv - 3.2, sv + 3.2), show_ball=False, label=False)
    ax.set_ylim(zc(3) - A_ - D_ - 0.2, 0.5)
    for s in range(4):
        ax.text(sv + 3.25, zc(s), " подслой %d: слой %d, %s" % (s + 1, s // 2 + 1, "i − h" if s % 2 == 0 else "i + h"),
                fontsize=8, va="center")
    ax.axvline(sv, color="k", lw=0.8, ls=":")
    ax.set_title("Шов крупно (x = %.2f мм), сечение A, слои 1–2: перекрытие на одну нить и стяжки" % sv, fontsize=10)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)
    # весь пакет
    fig, ax = plt.subplots(figsize=(11.7, 8.0))
    section_axes(ax, B, j0, range(k_total), (-26, 26))
    zb_new, zb_old = zb(k_total)
    ax.axhline(0, color="k", lw=0.5, ls="--")
    for zz, lab in ((zb_old, "z_b сплошной %.2f" % zb_old), (zb_new, "z_b шахматной %.2f" % zb_new)):
        ax.axhline(zz, color="k", lw=0.7, ls="--")
        ax.text(-25.8, zz - 0.15, lab, fontsize=7, va="top", ha="left", zorder=7,
                bbox=dict(fc="w", ec="none", alpha=0.8, pad=0.5))
    ax.set_ylim(zb_new - 1.5, BALL_D + 1.5)
    ax.set_title("Сечение A (y = 0): %dT, кусочки 7×7 (ядро %d×%d), весь пролёт" % (k_total, patch, patch),
                 fontsize=10)
    pdf.savefig(fig)
    plt.close(fig)


class PdfPng:
    """PdfPages + копия каждой страницы в <png>_<n>.png (для просмотра без растеризатора PDF)."""

    def __init__(self, pdf, png):
        self.pdf, self.png, self.n = pdf, png, 0

    def savefig(self, fig):
        if VLABEL:
            fig.text(0.995, 0.005, VLABEL, fontsize=10, weight="bold", va="bottom", ha="right", color="#922b21")
        self.pdf.savefig(fig)
        if self.png:
            self.n += 1
            fig.savefig("%s_%d.png" % (self.png, self.n), dpi=110)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="chess_section.pdf")
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--patch", type=int, default=5)
    ap.add_argument("--png", default="", help="префикс PNG-копий страниц")
    ap.add_argument("--variants", default="flip,same", help="flip -- со сменой чётности, same -- без")
    a = ap.parse_args()
    global FLIP, VLABEL
    B = Board(seams_idx(N, a.patch))
    names = {"flip": "Вариант 1: у слоя i + 1 чётность доски обратная",
             "same": "Вариант 2: у всех слоёв одна доска (без смены чётности)"}
    with PdfPages(a.out) as pdf0:
        pdf = PdfPng(pdf0, a.png)
        for v in a.variants.split(","):
            FLIP, VLABEL = v == "flip", names[v]
            n_st = len(stitches(B, 0))
            chk = check(B)
            print(v, "ядра:", B.core, "; проверка:", chk, "; стяжек на слой:", n_st)
            plan_page(pdf, B, a.patch)
            section_pages(pdf, B, a.patch, a.k, chk, n_st)
    print("записан", a.out)


if __name__ == "__main__":
    main()
