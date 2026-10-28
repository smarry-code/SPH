#!/usr/bin/env python3
"""Анимация выстрела Guo 2019 (1T) с акцентом на нити и их переплетение. СКОРОСТЬ НЕ ОТОБРАЖАЕТСЯ.

Нити раскрашены по направлению: основа (вдоль y, зажатая) -- тёмно-синяя, уток (вдоль x, свободная) --
оранжевая. Порванные связи в кадре отсутствуют; окно рисования ограничено (--weave-window по x,y и
пределы по z), поэтому оторванные концы вне окна не засоряют кадр.

Три панели: (а) 3D в окне вокруг места удара; (б) вид сверху (x-y), сегменты рисуются в порядке по высоте z
(выше -- сверху), яркость слегка зависит от высоты бусины над локальной средней по нити (так читаются
верх/низ переплетения), контур шарика -- пунктир; (в) боковое сечение y-z (полоса около центра).
Заголовок: V0, mu (из ../shots.csv), время, прогиб w_f, число порванных вертикальных/горизонтальных нитей,
итог ('НЕ ПРОБИТ' зелёный / 'ПРОБИТ, V_R = ... м/с' красный).

Режим --compare: несколько --run-dir, в строку панели "вид сверху" на одном моменте времени (ближайший кадр).

Usage:
    python animate_weave.py --run-dir mu18/v76.0/1T_base_v76.0_r0 --out weave_out/v76_weave.gif --max-frames 40 --snapshots auto
    python animate_weave.py --compare --run-dir A --run-dir B --out weave_out/compare_weave.gif --max-frames 40
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Line3DCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from animate_ply import (discover_frame_ids, parse_sph_frame, parse_fiber_frame,  # noqa: E402
                         read_monitor, read_layout)

C_WARP = "#16348c"   # основа (вдоль y, зажата)
C_WEFT = "#f07a10"   # уток (вдоль x, свободна)
BG_TOP = "#eceef2"
GREEN, RED = "#1a8a2e", "#c4161c"


# ---------------------------------------------------------------------------

def read_shot(run_dir: str) -> dict:
    """V0, mu, outcome, V_R из shots.csv рядом с каталогом выстрела (../shots.csv)."""
    for p in (os.path.join(run_dir, "..", "shots.csv"), os.path.join(run_dir, "shots.csv")):
        if os.path.exists(p):
            with open(p, newline="") as f:
                rows = list(csv.DictReader(f))
            if rows:
                return rows[0]
    return {}


def yarn_ranges(run_dir: str, n: int):
    """Для каждой бусины первая/последняя бусина её нити (по раскладке)."""
    first = np.arange(n)
    last = np.arange(n)
    p = os.path.join(run_dir, "guo2019_fabric_layout.txt")
    if os.path.exists(p):
        for ln in open(p):
            t = ln.split()
            if t and not t[0].startswith("#") and len(t) >= 7:
                a, b = int(t[5]), int(t[6])
                first[a:b + 1] = a
                last[a:b + 1] = b
    return first, last


class Run:
    pass


def load_run(run_dir: str, max_frames: int) -> Run:
    r = Run()
    r.dir = run_dir
    ids = discover_frame_ids(run_dir)
    if not ids:
        sys.exit(f"{run_dir}: нет пар sph_output/fiber_output")
    if max_frames and len(ids) > max_frames:
        sel = np.linspace(0, len(ids) - 1, max_frames).round().astype(int)
        ids = [ids[i] for i in sorted(set(sel.tolist()))]
    t0 = time.time()
    r.ball = [parse_sph_frame(os.path.join(run_dir, f"sph_output_{i}.txt")) for i in ids]
    r.fib = [parse_fiber_frame(os.path.join(run_dir, f"fiber_output_{i}.txt")) for i in ids]
    print(f"  {run_dir}: {len(ids)} кадров, чтение {time.time() - t0:.1f} с")
    n = len(r.fib[0].pos)
    r.dirs, _, _, r.pitch = read_layout(run_dir, n)
    r.pitch = r.pitch or 1.18e-3
    r.first, r.last = yarn_ranges(run_dir, n)
    r.idx = np.arange(n)
    mon = read_monitor(run_dir)
    r.shot = read_shot(run_dir)
    r.times = np.array([f.time for f in r.fib])
    r.wf = np.zeros(len(ids)); r.nbv = np.zeros(len(ids), int); r.nbh = np.zeros(len(ids), int)
    if mon is not None:
        steps = mon["step"].astype(int)
        for k, f in enumerate(r.fib):
            j = int(np.argmin(np.abs(steps - f.step)))
            r.times[k] = mon["t"][j]
            r.wf[k] = mon["w_f"][j]
            r.nbv[k] = int(mon["n_broken_v"][j]); r.nbh[k] = int(mon["n_broken_h"][j])
    c0 = r.fib[0].pos
    r.ctr = 0.5 * (c0.min(axis=0) + c0.max(axis=0))[:2]
    r.r_ball = 0.5 * float(np.ptp(r.ball[0][:, 2])) if len(r.ball[0]) else 3e-3
    return r


def local_z(pos: np.ndarray, run: Run) -> np.ndarray:
    """Высота бусины над скользящей средней вдоль нити (ядро [.5 1 1 1 .5]/4 гасит период-4 переплетения)."""
    z = pos[:, 2]
    w = {-2: 0.5, -1: 1.0, 0: 1.0, 1: 1.0, 2: 0.5}
    m = np.zeros_like(z)
    for k, wk in w.items():
        m += wk * z[np.clip(run.idx + k, run.first, run.last)]
    return z - m / 4.0


def shade(base: str, s: np.ndarray, amt: float = 0.40) -> np.ndarray:
    """Осветление (s>0) / затемнение (s<0) базового цвета; s в [-1,1]."""
    b = np.array(to_rgb(base))
    s = np.clip(s, -1, 1)[:, None]
    up = b + (1 - b) * (s * amt)
    dn = b * (1 + s * amt)
    return np.where(s >= 0, up, dn)


def window_mask(pos, run, W):
    return (np.abs(pos[:, 0] - run.ctr[0]) <= W) & (np.abs(pos[:, 1] - run.ctr[1]) <= W)


def frame_geometry(run: Run, k: int, W: float, zlo: float, zhi: float, shade_mm: float):
    """Сегменты связей в окне, цвета (по направлению + оттенок по z), порядок по z."""
    f = run.fib[k]
    ok = window_mask(f.pos, run, W) & (f.pos[:, 2] >= zlo) & (f.pos[:, 2] <= zhi)
    keep = ok[f.ij[:, 0]] & ok[f.ij[:, 1]]
    ij = f.ij[keep]
    seg = np.stack([f.pos[ij[:, 0]], f.pos[ij[:, 1]]], axis=1)
    warp = run.dirs[ij[:, 0]] == 1
    lz = local_z(f.pos, run)
    s = 0.5 * (lz[ij[:, 0]] + lz[ij[:, 1]]) / (shade_mm * 1e-3)
    col = np.empty((len(ij), 3))
    col[warp] = shade(C_WARP, s[warp], 0.55)
    col[~warp] = shade(C_WEFT, s[~warp], 0.45)
    zm = seg[:, :, 2].mean(axis=1)
    return seg, col, warp, np.argsort(zm, kind="stable"), s


def sphere_mesh(c, r, n=14):
    u = np.linspace(0, 2 * np.pi, 2 * n)
    v = np.linspace(0, np.pi, n)
    return (c[0] + r * np.outer(np.cos(u), np.sin(v)), c[1] + r * np.outer(np.sin(u), np.sin(v)),
            c[2] + r * np.outer(np.ones_like(u), np.cos(v)))


def ball_center(run, k):
    b = run.ball[k]
    return np.array([b[:, 0].mean(), b[:, 1].mean(), 0.5 * (b[:, 2].min() + b[:, 2].max())])


def title_lines(run: Run, k: int, short=False):
    sh = run.shot
    v0 = float(sh["V0"]) if "V0" in sh else float("nan")
    mu = sh.get("mu", "?")
    l1 = f"V0 = {v0:g} м/с   mu = {mu}   t = {run.times[k] * 1e6:.0f} мкс   прогиб w_f = {run.wf[k] * 1e3:.1f} мм"
    l2 = f"порвано нитей: вертикальных (основа) {run.nbv[k]}, горизонтальных (уток) {run.nbh[k]}"
    oc = sh.get("outcome", "")
    if oc == "not_perforated":
        res, colr = "НЕ ПРОБИТ", GREEN
    elif oc:
        res, colr = f"ПРОБИТ, V_R = {float(sh.get('V_R', 0)):.1f} м/с", RED
    else:
        res, colr = "", "k"
    return l1, l2, res, colr


def auto_zlim(runs, W):
    lo, hi = 1e9, -1e9
    for r in runs:
        for f in r.fib:
            m = window_mask(f.pos, r, W)
            zz = f.pos[m, 2]
            zz = zz[(zz > np.percentile(zz, 1.0)) & (zz < np.percentile(zz, 99.8))] if len(zz) > 20 else zz
            lo, hi = min(lo, zz.min()), max(hi, zz.max())
        hi = max(hi, max(b[:, 2].max() for b in r.ball[:1]))
        rb = r.r_ball
        lo = min(lo, min(b[:, 2].min() for b in r.ball)) if False else lo
    return lo - 0.08 * (hi - lo), hi + 0.06 * (hi - lo)


# ---------------------------------------------------------------------------
# Панели

class TopPanel:
    def __init__(self, ax, run, W, zlo, zhi, args, label=True):
        self.ax, self.run, self.W, self.zlo, self.zhi, self.args = ax, run, W, zlo, zhi, args
        ax.set_facecolor(BG_TOP)
        ax.set_xlim(run.ctr[0] - W, run.ctr[0] + W)
        ax.set_ylim(run.ctr[1] - W, run.ctr[1] + W)
        ax.set_aspect("equal", adjustable="box")
        ax.xaxis.set_major_locator(MaxNLocator(5)); ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.tick_params(labelsize=7)
        ax.set_xlabel("x, м", fontsize=8); ax.set_ylabel("y, м", fontsize=8)
        self.lc = LineCollection([], capstyle="butt")
        ax.add_collection(self.lc)
        self.circ = Circle((0, 0), run.r_ball, fill=False, ec="k", lw=1.2, ls=(0, (4, 3)), zorder=5)
        ax.add_patch(self.circ)
        self.fixed_lw = None

    def lw_pts(self):
        # толщина в пунктах ~ 0.78 шага нитей (чтобы между нитями оставались просветы)
        self.ax.figure.canvas.draw_idle()
        bb = self.ax.get_window_extent()
        pts_per_m = min(bb.width, bb.height) * 72.0 / self.ax.figure.dpi / (2 * self.W)
        return max(self.args.top_lw * self.run.pitch * pts_per_m, 1.0)

    def update(self, k):
        r = self.run
        seg, col, warp, order, _ = frame_geometry(r, k, self.W, self.zlo, self.zhi, self.args.shade_mm)
        if self.fixed_lw is None:
            self.fixed_lw = self.lw_pts()
        self.lc.set_segments(seg[order][:, :, :2])
        self.lc.set_color(col[order])
        self.lc.set_linewidths(self.fixed_lw)
        c = ball_center(r, k)
        self.circ.center = (c[0], c[1])


class SectionPanel:
    def __init__(self, ax, run, W, zlo, zhi, args):
        self.ax, self.run, self.W, self.zlo, self.zhi = ax, run, W, zlo, zhi
        self.band = args.band or 1.5 * run.pitch
        ax.set_xlim(run.ctr[1] - W, run.ctr[1] + W)
        ax.set_ylim(zlo, zhi)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=7)
        ax.set_xlabel("y, м", fontsize=8); ax.set_ylabel("z, м", fontsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(5)); ax.yaxis.set_major_locator(MaxNLocator(5))
        self.lc = LineCollection([], linewidths=1.3)
        ax.add_collection(self.lc)
        self.dots = ax.scatter([], [], s=7, zorder=3)
        self.ballp, = ax.plot([], [], ".", ms=1.5, color="0.55", zorder=1)
        self.circ = Circle((0, 0), run.r_ball, fill=False, ec="0.2", lw=1.0, ls="--", zorder=4)
        ax.add_patch(self.circ)

    def update(self, k):
        r = self.run
        f = r.fib[k]
        inb = (np.abs(f.pos[:, 0] - r.ctr[0]) < self.band) & (np.abs(f.pos[:, 1] - r.ctr[1]) <= self.W)
        bi = inb[f.ij[:, 0]] & inb[f.ij[:, 1]] & (r.dirs[f.ij[:, 0]] == 1)
        seg = np.stack([f.pos[f.ij[bi, 0]], f.pos[f.ij[bi, 1]]], axis=1)[:, :, [1, 2]]
        self.lc.set_segments(seg)
        self.lc.set_color(C_WARP)
        weft = inb & (r.dirs == 0)
        self.dots.set_offsets(np.column_stack([f.pos[weft, 1], f.pos[weft, 2]]))
        self.dots.set_facecolor(C_WEFT); self.dots.set_edgecolor("none")
        b = r.ball[k]
        sel = np.abs(b[:, 0] - r.ctr[0]) < self.band
        self.ballp.set_data(b[sel, 1], b[sel, 2])
        c = ball_center(r, k)
        self.circ.center = (c[1], c[2])


class Panel3D:
    def __init__(self, ax, run, W, zlo, zhi, args):
        self.ax, self.run, self.W, self.zlo, self.zhi, self.args = ax, run, W, zlo, zhi, args
        ax.set_xlim(run.ctr[0] - W, run.ctr[0] + W)
        ax.set_ylim(run.ctr[1] - W, run.ctr[1] + W)
        ax.set_zlim(zlo, zhi)
        ax.set_box_aspect((1, 1, max((zhi - zlo) / (2 * W), 0.35)))
        ax.view_init(elev=args.elev, azim=args.azim)
        ax.set_xlabel("x, м", fontsize=8); ax.set_ylabel("y, м", fontsize=8); ax.set_zlabel("z, м", fontsize=8)
        ax.tick_params(labelsize=6)
        for a in (ax.xaxis, ax.yaxis, ax.zaxis):
            a.set_major_locator(MaxNLocator(4))
        self.lc = Line3DCollection(np.zeros((0, 2, 3)))
        ax.add_collection3d(self.lc, autolim=False)
        self.surf = None
        ax.legend(handles=[Line2D([], [], color=C_WARP, lw=3, label="основа (вдоль y, зажата)"),
                           Line2D([], [], color=C_WEFT, lw=3, label="уток (вдоль x, свободна)"),
                           Line2D([], [], color="0.55", marker="o", ls="", ms=6, alpha=0.6, label="шарик")],
                  loc="upper left", fontsize=7, framealpha=0.85)

    def update(self, k):
        r = self.run
        seg, col, warp, order, _ = frame_geometry(r, k, self.W, self.zlo, self.zhi, self.args.shade_mm)
        self.lc.set_segments(seg[order])
        self.lc.set_color(col[order])
        self.lc.set_linewidths(self.args.lw * np.where(warp[order], 1.5, 1.4))
        if self.surf is not None:
            self.surf.remove()
        c = ball_center(r, k)
        X, Y, Z = sphere_mesh(c, r.r_ball)
        self.surf = self.ax.plot_surface(X, Y, Z, color="0.45", alpha=self.args.ball_alpha, linewidth=0.2,
                                         edgecolor="0.3", shade=True, rcount=14, ccount=28, zorder=1)


# ---------------------------------------------------------------------------

def pick_snaps(run: Run, spec: str, nf: int):
    if spec == "auto":
        s = [0]
        hit = np.where(run.nbv + run.nbh > 0)[0]
        if len(hit):
            s.append(int(hit[0]))
        s.append(int(np.argmax(run.wf)))
        s.append(nf - 1)
        return sorted(set(s))
    return [int(x) for x in spec.split(",") if x.strip()]


def main():
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(encoding="utf-8")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", action="append", required=True, help="каталог прогона (для --compare можно несколько)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--compare", action="store_true", help="сравнение: вид сверху для нескольких выстрелов в строку")
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--fps", type=float, default=6.0)
    ap.add_argument("--dpi", type=int, default=80)
    ap.add_argument("--snap-dpi-factor", type=float, default=2.0)
    ap.add_argument("--snapshots", default=None, help="индексы или 'auto'")
    ap.add_argument("--no-gif", action="store_true")
    ap.add_argument("--zoom", default=None, help="синоним --weave-window, мм->м не конвертируется: значение в метрах")
    ap.add_argument("--weave-window", type=float, default=15e-3, help="полуразмер окна, м (по умолчанию 0.015)")
    ap.add_argument("--zlim", default="auto", help="'lo,hi' в мм или auto")
    ap.add_argument("--elev", type=float, default=32.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--lw", type=float, default=1.0, help="множитель толщины линий 3D (база 1.5/1.4 pt)")
    ap.add_argument("--band", type=float, default=0.0, help="полуширина полосы сечения, м (0 = 1.5 шага)")
    ap.add_argument("--shade-mm", type=float, default=0.25, help="превышение над средней по нити (мм), дающее полный оттенок")
    ap.add_argument("--top-lw", type=float, default=0.5, help="толщина нитей на виде сверху, доля шага нитей")
    ap.add_argument("--ball-alpha", type=float, default=0.28)
    args = ap.parse_args()
    W = float(args.zoom) if args.zoom else args.weave_window

    runs = [load_run(d, args.max_frames) for d in (args.run_dir if args.compare else args.run_dir[:1])]
    if args.zlim != "auto":
        zlo, zhi = (float(v) * 1e-3 for v in args.zlim.split(","))
    else:
        zlo, zhi = auto_zlim(runs, W)
    print(f"Окно ±{W * 1e3:.1f} мм, z [{zlo * 1e3:.2f}, {zhi * 1e3:.2f}] мм")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    base = os.path.splitext(args.out)[0]

    if not args.compare:
        run = runs[0]
        nf = len(run.fib)
        fig = plt.figure(figsize=(17, 6.6))
        ax3 = fig.add_axes((0.0, 0.02, 0.40, 0.84), projection="3d")
        axT = fig.add_axes((0.405, 0.08, 0.31, 0.76))
        axS = fig.add_axes((0.735, 0.08, 0.255, 0.76))
        p3, pT, pS = Panel3D(ax3, run, W, zlo, zhi, args), TopPanel(axT, run, W, zlo, zhi, args), \
            SectionPanel(axS, run, W, zlo, zhi, args)
        ax3.set_title("(а) 3D", fontsize=9, y=0.98)
        axT.set_title("(б) вид сверху: нити сверху ярче, снизу темнее; пунктир -- шарик", fontsize=8)
        axS.set_title(f"(в) сечение y-z, |x| < {pS.band * 1e3:.1f} мм", fontsize=9)
        t1 = fig.text(0.015, 0.975, "", fontsize=12, va="top")
        t2 = fig.text(0.015, 0.930, "", fontsize=10.5, va="top")
        t3 = fig.text(0.015, 0.888, "", fontsize=13, va="top", fontweight="bold")

        def update(k):
            p3.update(k); pT.update(k); pS.update(k)
            l1, l2, res, colr = title_lines(run, k)
            t1.set_text(l1); t2.set_text(l2); t3.set_text(res); t3.set_color(colr)

        snaps = pick_snaps(run, args.snapshots, nf) if args.snapshots else []
        total = nf
    else:
        n = len(runs)
        fig = plt.figure(figsize=(5.6 * n, 6.9))
        panels, tt = [], []
        lm, gap = 0.04, 0.02
        wpan = (1 - lm - 0.01 - gap * (n - 1)) / n
        for i, r in enumerate(runs):
            ax = fig.add_axes((lm + i * (wpan + gap), 0.10, wpan, 0.70))
            panels.append(TopPanel(ax, r, W, zlo, zhi, args))
            tt.append((fig.text(lm + i * (wpan + gap), 0.965, "", fontsize=11, va="top"),
                       fig.text(lm + i * (wpan + gap), 0.915, "", fontsize=9.5, va="top"),
                       fig.text(lm + i * (wpan + gap), 0.865, "", fontsize=13, va="top", fontweight="bold")))
        ref = max(runs, key=lambda r: len(r.fib))
        tgrid = ref.times
        idx_map = [[int(np.argmin(np.abs(r.times - t))) for t in tgrid] for r in runs]
        sup = fig.text(0.5, 0.008, "", fontsize=11, ha="center")

        def update(k):
            for i, (p, r) in enumerate(zip(panels, runs)):
                kk = idx_map[i][k]
                p.update(kk)
                l1, l2, res, colr = title_lines(r, kk)
                tt[i][0].set_text(f"V0 = {float(r.shot.get('V0', 'nan')):g} м/с   mu = {r.shot.get('mu', '?')}   "
                                  f"w_f = {r.wf[kk] * 1e3:.1f} мм")
                tt[i][1].set_text(l2.replace("порвано нитей: ", "порвано: ").replace("вертикальных (основа)", "верт.")
                                  .replace("горизонтальных (уток)", "гориз."))
                tt[i][2].set_text(res); tt[i][2].set_color(colr)
            sup.set_text(f"вид сверху, t = {tgrid[k] * 1e6:.0f} мкс (ближайший кадр каждого выстрела; "
                         f"основа -- синяя, уток -- оранжевая)")

        total = len(tgrid)
        if args.snapshots == "auto":
            snaps = sorted({0, total // 3, total // 2, 2 * total // 3, total - 1})
        else:
            snaps = [int(x) for x in args.snapshots.split(",") if x.strip()] if args.snapshots else []

    for si in snaps:
        if 0 <= si < total:
            t0 = time.time()
            update(si)
            p = f"{base}_frame{si}.png"
            fig.savefig(p, dpi=int(args.dpi * args.snap_dpi_factor))
            print(f"снимок {p} ({time.time() - t0:.1f} с)")
    if not args.no_gif:
        t0 = time.time()
        anim = FuncAnimation(fig, update, frames=total, interval=1000 / args.fps, blit=False)
        anim.save(args.out, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
        dt = time.time() - t0
        print(f"GIF {args.out}: {dt:.1f} с ({dt / total:.2f} с/кадр), {os.path.getsize(args.out) / 1e6:.2f} МБ, {total} кадров")


if __name__ == "__main__":
    main()
