#!/usr/bin/env python3
"""GIF прогона Guo 1T: ткань в косом ракурсе с ЛИЦЕВОЙ (со стороны удара, +z) и ТЫЛЬНОЙ (−z) сторон рядом,
как на слайде лекции. Крупно: окно ±--zoom мм вокруг центра удара, толстые нити яркой шкалой eps/eps_break,
связи с начавшимся повреждением (D > 0) -- розовым поверх. Шарик -- зелёные частицы. Чтение кадров -- из
animate_1t.py (тот же каталог). Первая версия (вся ткань, тонкие линии) -- anim_2side_v1.py.

Usage: python anim_2side.py --run-dir <прогон> --out anim.gif [--max-frames 40] [--zoom 35] [--title "..."]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d.art3d import Line3DCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from animate_1t import discover_frame_ids, parse_sph_frame, parse_fiber_frame, read_monitor, ratio_of  # noqa: E402

# яркая шкала: фиолетовый (не нагружена) -> синий -> бирюзовый -> жёлтый -> красный (у предела)
CMAP = LinearSegmentedColormap.from_list("fabric", ["#8a5cf5", "#2d7dff", "#13c2c2", "#ffd21f", "#ff2a2a"])
DMG_COLOR = "#ff00c8"  # связи с D > 0 (разупрочнение) -- поверх, толще
CLAMP_COLOR = "#2b2b2b"  # закреплённые бусины (концы вертикальных нитей)


def clamped_beads(run_dir: str, n: int) -> np.ndarray:
    """Маска закреплённых бусин по guo2019_fabric_layout.txt (n_clamp_lo / n_clamp_hi у концов нити; у 2T первый
    столбец -- номер слоя)."""
    m = np.zeros(n, dtype=bool)
    path = os.path.join(run_dir, "guo2019_fabric_layout.txt")
    if not os.path.exists(path):
        return m
    off = 0
    for ln in open(path):
        p = ln.split()
        if ln.startswith("#"):
            off = 1 if "ply" in p else 0
            continue
        if len(p) < 9 + off:
            continue
        a0, b0, lo, hi = int(p[5 + off]), int(p[6 + off]), int(p[7 + off]), int(p[8 + off])
        if lo > 0:
            m[a0:a0 + lo] = True
        if hi > 0:
            m[b0 - hi + 1:b0 + 1] = True
    return m


def damage_of(path: str, nb: int) -> np.ndarray:
    """Колонка D (11-я) BONDS, если повреждённость включена; иначе нули."""
    with open(path) as f:
        for ln in f:
            if ln.strip() == "BONDS":
                break
        rows = [next(f).split() for _ in range(nb)]
    return np.array([float(r[10]) if len(r) > 10 else 0.0 for r in rows])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-frames", type=int, default=40)
    ap.add_argument("--fps", type=float, default=6.0)
    ap.add_argument("--dpi", type=int, default=90)
    ap.add_argument("--vmax", type=float, default=0.65)
    ap.add_argument("--elev", type=float, default=18.0)
    ap.add_argument("--azim", type=float, default=45.0)
    ap.add_argument("--zoom", type=float, default=0.0, help="полуокно вокруг центра удара по x, y, мм (0 = вся ткань)")
    ap.add_argument("--lw", type=float, default=0.7, help="толщина нитей")
    ap.add_argument("--title", default="")
    ap.add_argument("--snapshots", default="auto", help="'auto' (первый разрыв, макс. прогиб, последний) или ''")
    a = ap.parse_args()

    ids = discover_frame_ids(a.run_dir)
    if not ids:
        sys.exit(f"{a.run_dir}: нет кадров")
    if len(ids) > a.max_frames:
        sel = np.linspace(0, len(ids) - 1, a.max_frames).round().astype(int)
        ids = [ids[i] for i in sorted(set(sel.tolist()))]
    nf = len(ids)
    t0 = time.time()
    fpaths = [os.path.join(a.run_dir, f"fiber_output_{i}.txt") for i in ids]
    ball = [parse_sph_frame(os.path.join(a.run_dir, f"sph_output_{i}.txt")) for i in ids]
    fib = [parse_fiber_frame(p) for p in fpaths]
    dmg = [damage_of(p, len(f.ij)) for p, f in zip(fpaths, fib)]
    print(f"{nf} кадров, чтение {time.time() - t0:.0f} с")

    mon = read_monitor(a.run_dir)
    ms = None
    if mon is not None:
        steps = mon["step"].astype(int)
        rows = [int(np.argmin(np.abs(steps - f.step))) for f in fib]
        ms = {k: mon[k][rows] for k in ("t", "v_z", "w_f", "n_broken_v", "n_broken_h") if k in mon}
    times = ms["t"] if ms else np.array([f.time for f in fib])

    p0 = fib[0].pos
    lo, hi = p0.min(axis=0), p0.max(axis=0)
    ctr = 0.5 * (lo + hi)
    half = a.zoom * 1e-3 if a.zoom > 0 else 0.5 * (hi[:2] - lo[:2]).max()
    # бусины окна -- по начальному положению (набор не меняется от кадра к кадру)
    win = (np.abs(p0[:, 0] - ctr[0]) <= half) & (np.abs(p0[:, 1] - ctr[1]) <= half)
    keep = [win[f.ij[:, 0]] & win[f.ij[:, 1]] for f in fib]
    clamp = clamped_beads(a.run_dir, len(p0))
    ckeep = [clamp[f.ij[:, 0]] & clamp[f.ij[:, 1]] & kk for f, kk in zip(fib, keep)]
    r_ball = 0.5 * float(np.ptp(ball[0][:, 2])) if len(ball[0]) else 3e-3
    fz = np.concatenate([f.pos[win, 2] for f in fib])
    zlo = np.percentile(fz, 0.1) - 2.5 * r_ball
    zhi = max(np.percentile(fz, 99.9), ball[0][:, 2].max() if len(ball[0]) else 0.0) + 0.5 * r_ball

    fig = plt.figure(figsize=(16, 8.6), facecolor="white")
    lcs, scs, lds, lks = [], [], [], []
    # на экране: нормаль ткани z -> горизонталь X, x ткани -> Y, y ткани (вдоль вертикальных нитей) -> вверх
    for k, (azim, lab) in enumerate(((a.azim, "Лицевая сторона (со стороны удара)"),
                                     (180.0 - a.azim, "Тыльная сторона"))):
        ax = fig.add_axes((-0.02 + 0.47 * k, -0.04, 0.53, 0.94), projection="3d", computed_zorder=False)
        ax.set_xlim(zlo, zhi); ax.set_ylim(ctr[0] - half, ctr[0] + half); ax.set_zlim(ctr[1] - half, ctr[1] + half)
        ax.set_box_aspect(((zhi - zlo) / (2 * half), 1, 1))
        ax.view_init(elev=a.elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(lab, fontsize=13, fontweight="bold", y=1.02)
        lc = Line3DCollection(np.zeros((1, 2, 3)), cmap=CMAP, clim=(0.0, a.vmax), linewidths=a.lw, zorder=1)
        ax.add_collection3d(lc, autolim=False)
        ld = Line3DCollection(np.zeros((1, 2, 3)), colors=DMG_COLOR, linewidths=3.0 * a.lw, zorder=1.5)
        ax.add_collection3d(ld, autolim=False)
        lk = Line3DCollection(np.zeros((1, 2, 3)), colors=CLAMP_COLOR, linewidths=5.0, zorder=1.6)
        ax.add_collection3d(lk, autolim=False)
        sc = ax.scatter(*ball[0][:, [2, 0, 1]].T, c="#1f9a50", s=6, alpha=0.9, depthshade=True, zorder=2,
                        edgecolors="none")
        lcs.append(lc); scs.append(sc); lds.append(ld); lks.append(lk)
    cax = fig.add_axes((0.935, 0.22, 0.011, 0.5))
    cb = fig.colorbar(lcs[0], cax=cax)
    cb.set_label(r"$\varepsilon/\varepsilon_{break}$ связи;  розовым -- связи с D > 0", fontsize=9)
    head = fig.text(0.02, 0.975, "", fontsize=11, va="top", ha="left")
    zoom_txt = ("окно ±%.0f мм вокруг центра удара" % (half * 1e3) + ("" if clamp[win].any() else ", закрепление вне окна")
                if a.zoom > 0 else "вся ткань, концы вертикальных нитей закреплены")

    def update(k: int):
        f, b, kk = fib[k], ball[k], keep[k]
        seg = np.stack([f.pos[f.ij[:, 0]], f.pos[f.ij[:, 1]]], axis=1)[:, :, [2, 0, 1]]
        col = np.minimum(ratio_of(f), a.vmax)
        dsel = kk & (dmg[k] > 0)
        bvis = b[b[:, 2] >= zlo] if len(b) else b
        # шарик спереди или сзади ткани (по самой прогнутой точке под шариком)
        bc = b.mean(axis=0) if len(b) else ctr
        near = np.hypot(f.pos[:, 0] - bc[0], f.pos[:, 1] - bc[1]) < r_ball
        z_fab = f.pos[near, 2].min() if near.any() else ctr[2]
        in_front = bc[2] > z_fab
        for j, (lc, sc, ld, lk) in enumerate(zip(lcs, scs, lds, lks)):
            lk.set_segments(np.zeros((0, 2, 3)))  # закрепление не выделяется (владелец, 2026-10-27)
            lc.set_segments(seg[kk])
            lc.set_array(col[kk])
            ld.set_segments(seg[dsel] if dsel.any() else np.zeros((0, 2, 3)))
            sc._offsets3d = tuple(bvis[:, [2, 0, 1]].T)
            sc.set_zorder(2 if (j == 0) == in_front else 0)
        s = (a.title + "\n" if a.title else "") + f"t = {times[k] * 1e6:.0f} мкс"
        if ms:
            s += (f"   v_z шара = {ms['v_z'][k]:.1f} м/с   прогиб = {ms['w_f'][k] * 1e3:.1f} мм"
                  f"   порвано нитей верт./гориз. = {int(ms['n_broken_v'][k])}/{int(ms['n_broken_h'][k])}")
        s += f"   связей с D > 0: {int((dmg[k] > 0).sum())}   ({zoom_txt})"
        head.set_text(s)

    base = os.path.splitext(a.out)[0]
    if a.snapshots == "auto" and ms:
        nb = ms["n_broken_v"] + ms["n_broken_h"]
        hit = np.where(nb > 0)[0]
        snaps = sorted({int(hit[0]) if len(hit) else 0, int(np.argmax(ms["w_f"])), nf - 1})
        for si in snaps:
            update(si)
            fig.savefig(f"{base}_frame{si}.png", dpi=int(a.dpi * 1.4))
    anim = FuncAnimation(fig, update, frames=nf, blit=False)
    t0 = time.time()
    anim.save(a.out, writer=PillowWriter(fps=a.fps), dpi=a.dpi)
    print(f"GIF {a.out} {os.path.getsize(a.out) / 1048576:.2f} МБ, {nf} кадров, рендер {time.time() - t0:.0f} с")


if __name__ == "__main__":
    main()
