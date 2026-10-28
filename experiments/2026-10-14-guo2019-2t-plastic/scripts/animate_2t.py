#!/usr/bin/env python3
"""Визуализация выстрелов 2T (две несклеенные ткани Twaron, Guo 2019) с пластично-вязкой нитью.

Версия для двух слоёв (animate_1t.py / animate_vel.py -- для 1T, не изменены). Использует animate_ply.read_layout
(колонка ply в guo2019_fabric_layout.txt). Каталог прогона: fiber_output_*.txt, sph_output_*.txt (шар),
monitor.csv, shot.txt, guo2019_fabric_layout.txt.

Подкоманды:
  vel      GIF + снимки: 3D (слой 1 и слой 2 разными цветами), шар, сечение yz с цветом по v_z, график v_z(t) с маркером
  strain   GIF + снимки: слой 1 и слой 2 рядом (вид сверху, со стороны шара), цвет eps/eps_break (--what eps)
           или eps_p (--what epsp); порванные связи -- красные крестики
  summary  PNG: v_z(t), прогиб w_f_p1/w_f_p2 (t), число порванных нитей по слоям (t)
  static   PNG+SVG: v_z(t) шара для набора каталогов с monitor.csv
  compare  GIF: два выстрела рядом, синхронно по времени

Примеры:
  python animate_2t.py vel --run-dir R --out anim/mu0.1_v137/vel_mu0.1_v137.gif --snapshots auto
  python animate_2t.py strain --run-dir R --out anim/.../anim_2t_mu0.1_v137.gif --what eps --vmax 0.75
  python animate_2t.py strain --run-dir R --out anim/.../anim_2t_epsp_mu0.1_v137.gif --what epsp --vmax 0.008
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Line3DCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from animate_ply import read_layout, read_monitor_file  # noqa: E402

PLY_COL = ["#D55E00", "#0072B2"]   # слой 1 (лицевой) -- вермильон, слой 2 (тыльный) -- синий
PLY_NAME = ["слой 1 (лицевой)", "слой 2 (тыльный)"]
CACHE_DIR = os.environ.get("ANIM2T_CACHE", "")


# ---------------------------------------------------------------------------
# Чтение и кэш
# ---------------------------------------------------------------------------

def parse_shot(run_dir):
    d = {}
    p = os.path.join(run_dir, "shot.txt")
    if os.path.exists(p):
        for tok in open(p).read().split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                d[k] = v
    return d


def fl(d, k, default=float("nan")):
    try:
        return float(d[k])
    except Exception:
        return default


def _parse_fiber(path):
    import pandas as pd
    with open(path) as f:
        h = f.readline().split()
    hd = {h[i]: h[i + 1] for i in range(0, len(h) - 1, 2)}
    n, nb = int(hd["N_PARTICLES"]), int(hd["N_BONDS"])
    pr = pd.read_csv(path, sep=r"\s+", header=None, skiprows=3, nrows=n, engine="c").values
    pos = np.zeros((n, 3), np.float32)
    vel = np.zeros((n, 3), np.float32)
    ids = pr[:, 0].astype(int)
    pos[ids], vel[ids] = pr[:, 3:6], pr[:, 6:9]
    br = pd.read_csv(path, sep=r"\s+", header=None, skiprows=3 + n + 2, nrows=nb, engine="c").values
    ij = br[:, :2].astype(np.int32)
    l0, lb, el = br[:, 2], br[:, 4], br[:, 5]
    eb = np.maximum(lb / l0 - 1.0, 1e-12)
    ratio = np.clip(el / eb, 0, None).astype(np.float32)
    epsp = br[:, 6].astype(np.float32) if br.shape[1] > 6 else np.zeros(nb, np.float32)
    return dict(time=float(hd["TIME"]), step=int(float(hd["STEP"])), pos=pos, vel=vel, ij=ij, ratio=ratio, epsp=epsp)


def _parse_ball(path):
    with open(path) as f:
        cols = f.readline().lstrip("# ").split()
    import pandas as pd
    d = pd.read_csv(path, sep=r"\s+", header=None, skiprows=1, engine="c").values
    ci = {c: i for i, c in enumerate(cols)}
    return d[:, [ci["x"], ci["y"], ci["z"]]].astype(np.float32)


def _load_one(args):
    fp, sp = args
    d = _parse_fiber(fp)
    d["ball"] = _parse_ball(sp)
    return d


class Run:
    pass


def load_run(run_dir, label=None):
    r = Run()
    r.dir = run_dir
    r.shot = parse_shot(run_dir)
    r.V0, r.mu = fl(r.shot, "V0"), fl(r.shot, "mu")
    r.outcome, r.VR = r.shot.get("outcome", ""), fl(r.shot, "V_R")
    r.t_break = fl(r.shot, "t_first_break")
    r.name = label or f"mu{r.mu:g}_v{r.V0:g}"
    r.mon = read_monitor_file(os.path.join(run_dir, "monitor.csv"))
    ids = sorted({re.search(r"_(\d+)\.txt$", f).group(1) for f in glob.glob(os.path.join(run_dir, "fiber_output_*.txt"))} &
                 {re.search(r"_(\d+)\.txt$", f).group(1) for f in glob.glob(os.path.join(run_dir, "sph_output_*.txt"))}, key=int)
    cache = os.path.join(CACHE_DIR, f"{os.path.basename(os.path.abspath(run_dir))}_{r.name}.npz") if CACHE_DIR else ""
    jobs = [(os.path.join(run_dir, f"fiber_output_{i}.txt"), os.path.join(run_dir, f"sph_output_{i}.txt")) for i in ids]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=min(8, len(jobs))) as ex:
        r.fr = list(ex.map(_load_one, jobs))
    print(f"{run_dir}: {len(r.fr)} кадров, {time.time() - t0:.0f} с")
    steps = r.mon["step"].astype(int)
    for f in r.fr:
        j = int(np.argmin(np.abs(steps - f["step"])))
        f["t"] = float(r.mon["t"][j]) if steps[j] == f["step"] else f["time"]
    r.times = np.array([f["t"] for f in r.fr])
    n = len(r.fr[0]["pos"])
    r.dirs, r.yidx, r.ply, r.pitch = read_layout(run_dir, n)
    r.ply = np.clip(r.ply, 0, 1)
    # исходные связи (кадр 0): ключи для поиска порванных
    ij0 = r.fr[0]["ij"].astype(np.int64)
    r.key0 = ij0[:, 0] * n + ij0[:, 1]
    r.ij0 = ij0
    r.n = n
    for f in r.fr:
        ij = f["ij"].astype(np.int64)
        f["alive"] = np.isin(r.key0, ij[:, 0] * n + ij[:, 1])
    ball0 = r.fr[0]["ball"]
    r.rb = 0.5 * float(np.ptp(ball0[:, 2]))
    c0 = r.fr[0]["pos"]
    r.ctr = 0.5 * (c0.min(axis=0) + c0.max(axis=0))[:2]
    r.half = 0.5 * float((c0.max(axis=0) - c0.min(axis=0))[:2].max())
    return r


def mon_at(r, t, key):
    if key not in r.mon:
        return float("nan")
    return float(np.interp(t, r.mon["t"], r.mon[key]))


def zoom_of(r, arg):
    if str(arg) != "auto":
        return float(arg)
    c0 = r.fr[0]["pos"]
    dz = np.max([np.abs(f["pos"][:, 2] - c0[:, 2]) for f in r.fr], axis=0)
    moved = dz > 0.05 * dz.max()
    rr = float(np.max(np.hypot(c0[moved, 0] - r.ctr[0], c0[moved, 1] - r.ctr[1])))
    z = min(max(1.25 * rr, 4 * r.rb), r.half)
    return 0.0 if z > 0.9 * r.half else z


def key_frames(r):
    """Индексы: до контакта, первый разрыв, максимальный прогиб, выход шара (шар целиком ниже ткани), конец."""
    nf = len(r.fr)
    out = {"до контакта": 0}
    c0z = r.fr[0]["pos"][:, 2].mean()
    pre = [k for k in range(nf) if mon_at(r, r.times[k], "z_c") - r.rb > c0z + 1e-4]
    if pre:
        out["до контакта"] = pre[-1]
    nb = r.mon["n_broken_v"] + r.mon["n_broken_h"]
    hit = np.where(nb > 0)[0]
    if len(hit):
        out["первый разрыв"] = int(np.argmin(np.abs(r.times - r.mon["t"][hit[0]])))
    out["макс. прогиб"] = int(np.argmin(np.abs(r.times - r.mon["t"][int(np.argmax(r.mon["w_f"]))])))
    ex = []
    for k in range(nf):
        f = r.fr[k]
        s_ = np.hypot(f["pos"][:, 0] - r.ctr[0], f["pos"][:, 1] - r.ctr[1]) < 2 * r.rb
        if f["ball"][:, 2].max() < np.median(f["pos"][s_, 2]):
            ex.append(k)
            break
    if ex:
        out["шар прошёл ткань"] = ex[0]   # верх шара ниже медианной z бусин ткани в круге 2*R_шара
    out["конец"] = nf - 1
    return out


def save_all(fig, update, nf, a, kf):
    base = os.path.splitext(a.out)[0]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    if a.snapshots == "auto":
        snaps = sorted(set(kf.values()))
    elif a.snapshots:
        snaps = [int(s) for s in a.snapshots.split(",")]
    else:
        snaps = []
    for si in snaps:
        update(si)
        p = f"{base}_frame{si}.png"
        fig.savefig(p, dpi=int(a.dpi * a.snap_dpi_factor))
        print("снимок:", p)
    if not a.no_gif:
        t0 = time.time()
        FuncAnimation(fig, update, frames=nf, interval=1000.0 / a.fps, blit=False).save(
            a.out, writer=PillowWriter(fps=a.fps), dpi=a.dpi)
        print(f"GIF: {a.out} ({os.path.getsize(a.out) / 1e6:.2f} МБ, {nf} кадров, {time.time() - t0:.0f} с)")


def outcome_str(r):
    if r.outcome == "perforated":
        return f"ПРОБИТ, V_R = {r.VR:.1f} м/с"
    if r.outcome == "not_perforated":
        return "НЕ ПРОБИТ"
    return ""


def head(r, k, extra=True):
    t = r.times[k]
    s = (f"2T, пластично-вязкая нить   V0 = {r.V0:g} м/с   μ = {r.mu:g}   t = {t * 1e6:.0f} мкс   "
         f"v_z шара = {mon_at(r, t, 'v_z'):.1f} м/с")
    if extra:
        s += (f"   порвано (верт./гориз.): сл.1 {int(round(mon_at(r, t, 'n_broken_v_p1')))}/{int(round(mon_at(r, t, 'n_broken_h_p1')))}"
              f", сл.2 {int(round(mon_at(r, t, 'n_broken_v_p2')))}/{int(round(mon_at(r, t, 'n_broken_h_p2')))}")
    return s


class VzPlot:
    def __init__(self, ax, r):
        self.r = r
        ax.axhline(0, color="0.3", lw=1)
        ax.plot(r.mon["t"] * 1e6, r.mon["v_z"], color="#1f77b4", lw=1.6)
        if np.isfinite(r.t_break):
            ax.axvline(r.t_break * 1e6, color="0.6", ls=":", lw=1)
            ax.text(r.t_break * 1e6, 0.98, " первый разрыв", transform=ax.get_xaxis_transform(), fontsize=7,
                    color="0.4", va="top", rotation=90)
        ax.set_xlim(0, float(r.mon["t"].max()) * 1e6)
        ax.set_xlabel("t, мкс")
        ax.set_ylabel("v_z шара, м/с")
        ax.grid(alpha=0.3)
        self.vline = ax.axvline(0, color="#d62728", lw=1.2)
        self.dot, = ax.plot([], [], "o", color="#d62728", ms=6, zorder=5)
        self.txt = ax.text(0.97, 0.06, "", transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
                           fontweight="bold", bbox=dict(fc="white", ec="0.7", alpha=0.9, boxstyle="round,pad=0.25"))

    def set(self, t):
        v = mon_at(self.r, t, "v_z")
        self.vline.set_xdata([t * 1e6] * 2)
        self.dot.set_data([t * 1e6], [v])
        self.txt.set_text(f"v_z = {v:.1f} м/с")


class Section:
    """Сечение yz (полоса по x около центра): слой 1 -- кружки, слой 2 -- ромбы, цвет по v_z; контур шара."""

    def __init__(self, ax, r, lo, hi, vzmax, band):
        self.r, self.band = r, band
        ax.set_xlim(lo[1], hi[1])
        ax.set_ylim(lo[2], hi[2])
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("y, м")
        ax.set_ylabel("z, м")
        ax.grid(alpha=0.3)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        self.sc = [ax.scatter([], [], s=s, c=[], cmap="RdBu_r", vmin=-vzmax, vmax=vzmax, marker=m, zorder=3 + i)
                   for i, (s, m) in enumerate([(11, "o"), (14, "D")])]
        self.ballp, = ax.plot([], [], ".", ms=2, color="0.6")
        self.circ = Circle((0, 0), r.rb, fill=False, ec="0.15", lw=1.2, ls="--")
        ax.add_patch(self.circ)
        ax.scatter([], [], s=11, c="0.5", marker="o", label="слой 1")
        ax.scatter([], [], s=14, c="0.5", marker="D", label="слой 2")
        ax.legend(fontsize=7, loc="upper right")

    def set(self, f):
        p = f["pos"]
        inb = np.abs(p[:, 0] - self.r.ctr[0]) < self.band
        for i in range(2):
            m = inb & (self.r.ply == i)
            self.sc[i].set_offsets(np.column_stack([p[m, 1], p[m, 2]]))
            self.sc[i].set_array(f["vel"][m, 2])
        bb = f["ball"]
        sel = np.abs(bb[:, 0] - self.r.ctr[0]) < self.band
        self.ballp.set_data(bb[sel, 1], bb[sel, 2])
        self.circ.center = (float(bb[:, 1].mean()), float(0.5 * (bb[:, 2].min() + bb[:, 2].max())))


def box(r, zoom, zlim=None):
    lo, hi = np.zeros(3), np.zeros(3)
    if zoom > 0:
        lo[:2], hi[:2] = r.ctr - zoom, r.ctr + zoom
    else:
        lo[:2], hi[:2] = r.ctr - r.half, r.ctr + r.half
    zb = np.concatenate([f["ball"][:, 2] for f in r.fr])
    zf = np.concatenate([f["pos"][::7, 2] for f in r.fr])   # прореженно; отдельные улетевшие бусины не растягивают ось
    lo[2], hi[2] = min(zb.min(), np.percentile(zf, 0.3)), max(zb.max(), zf.max())
    pad = 0.06 * (hi - lo)
    return lo - pad, hi + pad


# ---------------------------------------------------------------------------
# vel
# ---------------------------------------------------------------------------

def cmd_vel(a):
    r = load_run(a.run_dir)
    nf = len(r.fr)
    zoom = zoom_of(r, a.zoom)
    lo, hi = box(r, zoom)
    span = hi - lo
    aspect = np.maximum(span / span[:2].max(), 0.25)
    vzmax = a.vzmax or r.V0
    band = 1.5 * r.pitch if r.pitch > 0 else 1.5e-3
    print(f"zoom={zoom * 1e3:.1f} мм, полоса ±{band * 1e3:.2f} мм")
    fig = plt.figure(figsize=(17, 5.8))
    ax = fig.add_axes((0.0, 0.0, 0.34, 0.86), projection="3d")
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(tuple(aspect))
    ax.view_init(elev=a.elev, azim=a.azim)
    ax.set_xlabel("x, м"); ax.set_ylabel("y, м"); ax.set_zlabel("z, м")
    ax.tick_params(labelsize=7)
    for asp, axis in zip(aspect, (ax.xaxis, ax.yaxis, ax.zaxis)):
        axis.set_major_locator(MaxNLocator(nbins=5 if asp > 0.5 else 3))
    ax.set_title("нити: слой 1 (оранжевый), слой 2 (синий); шар -- серый", fontsize=9, y=1.0)
    lcs = []
    for i in range(2):
        lc = Line3DCollection(np.zeros((0, 2, 3)), colors=PLY_COL[i], linewidths=0.9 * a.lw, alpha=0.6 if i == 0 else 0.6)
        ax.add_collection3d(lc, autolim=False)
        lcs.append(lc)
    sc3 = ax.scatter(*r.fr[0]["ball"].T, c="0.35", s=3, alpha=0.35, depthshade=True)
    ax.plot([], [], color=PLY_COL[0], label=PLY_NAME[0]); ax.plot([], [], color=PLY_COL[1], label=PLY_NAME[1])
    ax.legend(loc="lower left", fontsize=8)
    axS = fig.add_axes((0.40, 0.11, 0.26, 0.70))
    sec = Section(axS, r, lo, hi, vzmax, band)
    axS.set_title(f"сечение yz, полоса ±{band * 1e3:.1f} мм: цвет = v_z", fontsize=10)
    fig.colorbar(sec.sc[0], ax=axS, shrink=0.8, pad=0.02).set_label("v_z бусины, м/с")
    axV = fig.add_axes((0.755, 0.13, 0.225, 0.68))
    vzp = VzPlot(axV, r)
    axV.set_title("v_z шара(t)", fontsize=10)
    t1 = fig.text(0.01, 0.985, "", fontsize=10.5, va="top")
    t2 = fig.text(0.01, 0.91, "", fontsize=15, va="top", fontweight="bold")
    ij0 = r.ij0
    stride_ok = np.ones(len(ij0), bool)
    ya, yb = r.yidx[ij0[:, 0]], r.yidx[ij0[:, 1]]
    stride_ok = ((ya % a.yarn_stride == 0) | (ya < 0)) & ((yb % a.yarn_stride == 0) | (yb < 0))
    plyb = r.ply[ij0[:, 0]]
    col = r.outcome == "perforated"

    def update(k):
        f = r.fr[k]
        pos = f["pos"]
        if zoom > 0:
            w = (np.abs(pos[:, 0] - r.ctr[0]) <= zoom) & (np.abs(pos[:, 1] - r.ctr[1]) <= zoom)
        else:
            w = np.ones(len(pos), bool)
        keep = f["alive"] & stride_ok & w[ij0[:, 0]] & w[ij0[:, 1]]
        for i in range(2):
            m = keep & (plyb == i)
            lcs[i].set_segments(np.stack([pos[ij0[m, 0]], pos[ij0[m, 1]]], axis=1))
        sc3._offsets3d = tuple(f["ball"].T)
        sec.set(f)
        vzp.set(r.times[k])
        t1.set_text(head(r, k))
        late = k == nf - 1 or (np.isfinite(r.t_break) and r.times[k] >= r.t_break)
        t2.set_text("Итог: " + outcome_str(r))
        t2.set_color(("#c62828" if col else "#1a7f37") if late else "0.25")

    save_all(fig, update, nf, a, key_frames(r))


# ---------------------------------------------------------------------------
# strain
# ---------------------------------------------------------------------------

def cmd_strain(a):
    r = load_run(a.run_dir)
    nf = len(r.fr)
    zoom = zoom_of(r, a.zoom)
    half = zoom if zoom > 0 else r.half
    ij0 = r.ij0
    plyb = r.ply[ij0[:, 0]]
    # индекс связей кадра 0 -> положение в списке связей кадра k (для живых)
    n = r.n
    fig, axs = plt.subplots(1, 2, figsize=(13.2, 6.9))
    fig.subplots_adjust(left=0.05, right=0.88, top=0.84, bottom=0.07, wspace=0.12)
    cm = plt.get_cmap("turbo")
    lcs, crs, balls = [], [], []
    for i, ax in enumerate(axs):
        ax.set_xlim(r.ctr[0] - half, r.ctr[0] + half)
        ax.set_ylim(r.ctr[1] - half, r.ctr[1] + half)
        ax.set_aspect("equal")
        ax.set_xlabel("x, м"); ax.set_ylabel("y, м")
        ax.set_title(PLY_NAME[i] + ": нити вдоль y -- зажатые, вдоль x -- свободные", fontsize=10, color=PLY_COL[i])
        ax.grid(alpha=0.25)
        lc = LineCollection([], cmap=cm, clim=(0, a.vmax), linewidths=0.9)
        ax.add_collection(lc)
        lcs.append(lc)
        crs.append(ax.scatter([], [], marker="x", c="k", s=36, linewidths=1.6, zorder=6))
        balls.append(Circle((0, 0), r.rb, fill=False, ec="0.2", lw=1.0, ls="--", zorder=5))
        ax.add_patch(balls[-1])
    cax = fig.add_axes((0.90, 0.2, 0.015, 0.55))
    if a.what == "eps":
        lab = r"$\varepsilon/\varepsilon_{break}$ связи ($\varepsilon_{break}$ = 0.06), шкала 0..%g" % a.vmax
    else:
        lab = r"$\varepsilon_p$ связи (пластическая деформация), шкала 0..%g" % a.vmax
    fig.colorbar(lcs[0], cax=cax).set_label(lab)
    fig.text(0.89, 0.12, "x -- порванные связи", fontsize=9, color="k")
    t1 = fig.text(0.01, 0.985, "", fontsize=10.5, va="top")
    t2 = fig.text(0.01, 0.925, "", fontsize=14, va="top", fontweight="bold")
    fig.text(0.01, 0.88, "вид сверху (со стороны шара); пунктир -- проекция шара", fontsize=9, color="0.35", va="top")
    # соответствие: для кадра k значения живых связей по порядку ij (ключи)
    keys_all = ij0[:, 0] * n + ij0[:, 1]
    sorter = np.argsort(keys_all)
    sk = keys_all[sorter]

    def update(k):
        f = r.fr[k]
        pos = f["pos"][:, :2]
        ijk = f["ij"].astype(np.int64)
        val = f["ratio"] if a.what == "eps" else f["epsp"]
        val = np.minimum(val, a.vmax)
        pk = r.ply[ijk[:, 0]]
        dead = ~f["alive"]
        for i in range(2):
            m = pk == i
            seg = np.stack([pos[ijk[m, 0]], pos[ijk[m, 1]]], axis=2).transpose(0, 2, 1)
            lcs[i].set_segments(seg)
            lcs[i].set_array(val[m])
            md = dead & (plyb == i)
            if md.any():
                mid = 0.5 * (pos[ij0[md, 0]] + pos[ij0[md, 1]])
                crs[i].set_offsets(mid)
            else:
                crs[i].set_offsets(np.zeros((0, 2)))
            bb = f["ball"]
            balls[i].center = (float(bb[:, 0].mean()), float(bb[:, 1].mean()))
        t = r.times[k]
        t1.set_text(head(r, k, extra=False) + f"   прогиб w_f = {mon_at(r, t, 'w_f') * 1e3:.1f} мм")
        late = k == nf - 1 or (np.isfinite(r.t_break) and t >= r.t_break)
        s = (f"порвано нитей верт./гориз.: слой 1 -- {int(round(mon_at(r, t, 'n_broken_v_p1')))}/{int(round(mon_at(r, t, 'n_broken_h_p1')))}"
             f", слой 2 -- {int(round(mon_at(r, t, 'n_broken_v_p2')))}/{int(round(mon_at(r, t, 'n_broken_h_p2')))}")
        if a.what == "epsp":
            s += f"   max eps_p кадра = {float(f['epsp'].max()):.4f}"
        t2.set_text(s)
        t2.set_color(("#c62828" if r.outcome == "perforated" else "#1a7f37") if late else "0.25")

    save_all(fig, update, nf, a, key_frames(r))


# ---------------------------------------------------------------------------
# summary / static / compare
# ---------------------------------------------------------------------------

def cmd_summary(a):
    mon = read_monitor_file(os.path.join(a.run_dir, "monitor.csv"))
    sh = parse_shot(a.run_dir)
    t = mon["t"] * 1e6
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.8))
    ax = axs[0]
    ax.plot(t, mon["v_z"], color="#1f77b4", lw=1.8)
    ax.axhline(0, color="0.3", lw=1)
    ax.set_title("скорость шара v_z(t)"); ax.set_ylabel("v_z, м/с")
    ax = axs[1]
    for i, c in enumerate(PLY_COL):
        ax.plot(t, mon[f"w_f_p{i + 1}"] * 1e3, color=c, lw=1.8, label=PLY_NAME[i])
    ax.set_title("прогиб слоя w_f_p1, w_f_p2 (t)"); ax.set_ylabel("прогиб, мм"); ax.legend(fontsize=9)
    ax = axs[2]
    for i, c in enumerate(PLY_COL):
        ax.step(t, mon[f"n_broken_v_p{i + 1}"], where="post", color=c, lw=1.8, label=f"{PLY_NAME[i]}: вертикальные")
        ax.step(t, mon[f"n_broken_h_p{i + 1}"], where="post", color=c, lw=1.8, ls="--", label=f"{PLY_NAME[i]}: горизонтальные")
    ax.set_title("число порванных нитей (накопленное)"); ax.set_ylabel("нитей"); ax.legend(fontsize=8)
    tb = fl(sh, "t_first_break") * 1e6
    for ax in axs:
        ax.set_xlabel("t, мкс"); ax.grid(alpha=0.3); ax.set_xlim(0, t.max())
        if np.isfinite(tb):
            ax.axvline(tb, color="0.6", ls=":", lw=1)
    fig.suptitle(f"2T base, μ = {fl(sh, 'mu'):g}, V0 = {fl(sh, 'V0'):g} м/с, пластично-вязкая нить; первый разрыв {tb:.0f} мкс (пунктир-точки)",
                 fontsize=11, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=130)
    print("записано:", a.out)


def cmd_static(a):
    """--spec 'каталог|подпись|цвет|стиль|t_обрыва_мкс(опц.)' ..."""
    fig, ax = plt.subplots(figsize=(10.5, 6))
    for s in a.spec:
        p = s.split("|")
        mon = read_monitor_file(os.path.join(p[0], "monitor.csv"))
        ax.plot(mon["t"] * 1e6, mon["v_z"], lw=1.8, color=p[2], ls=p[3], label=p[1])
        if p[3] != "-":
            ax.plot([mon["t"][-1] * 1e6], [mon["v_z"][-1]], "x", color=p[2], ms=9, mew=2)
    ax.axhline(0, color="0.3", lw=1)
    ax.set_xlabel("t, мкс"); ax.set_ylabel("v_z шара, м/с"); ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="lower right")
    ax.set_title("2T base, пластично-вязкая нить: скорость шара v_z(t). Пунктир -- прерванные расчёты (крест -- точка обрыва)", fontsize=10)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out)
    fig.savefig(os.path.splitext(a.out)[0] + ".png", dpi=150)
    print("записано:", a.out)


def cmd_compare(a):
    runs = [load_run(d) for d in a.run_dir]
    tmax = max(r.times.max() for r in runs)
    t_ref = np.arange(0, tmax + 1e-9, 20e-6)
    zoom = max(zoom_of(r, "auto") or r.half for r in runs)
    lo, hi = np.min([box(r, zoom)[0] for r in runs], axis=0), np.max([box(r, zoom)[1] for r in runs], axis=0)
    n = len(runs)
    fig = plt.figure(figsize=(6.2 * n + 0.8, 9.0))
    gs = fig.add_gridspec(2, n, left=0.06, right=0.90, top=0.82, bottom=0.07, hspace=0.28, wspace=0.28, height_ratios=[1.35, 1.0])
    secs, vzs, heads, outs = [], [], [], []
    vzmax = max(r.V0 for r in runs)
    for j, r in enumerate(runs):
        axS = fig.add_subplot(gs[0, j])
        secs.append(Section(axS, r, lo, hi, vzmax, 1.5 * r.pitch))
        axS.set_title(f"μ = {r.mu:g}, V0 = {r.V0:g} м/с", fontsize=10)
        vzs.append(VzPlot(fig.add_subplot(gs[1, j]), r))
        pos = axS.get_position()
        heads.append(fig.text(pos.x0, 0.955, "", fontsize=9.5, va="top"))
        outs.append(fig.text(pos.x0, 0.865, "", fontsize=14, va="top", fontweight="bold"))
    cax = fig.add_axes((0.92, 0.45, 0.012, 0.35))
    fig.colorbar(secs[0].sc[0], cax=cax).set_label("v_z бусины, м/с")
    ttl = fig.text(0.5, 0.995, "", ha="center", va="top", fontsize=11, color="0.3")

    def update(m):
        t = t_ref[m]
        ttl.set_text(f"t = {t * 1e6:.0f} мкс (сечение yz: слой 1 -- кружки, слой 2 -- ромбы)")
        for r, s, v, h, o in zip(runs, secs, vzs, heads, outs):
            k = int(np.argmin(np.abs(r.times - t)))
            s.set(r.fr[k])
            tt = min(t, r.times[-1])
            v.set(tt)
            h.set_text(head(r, k).replace("   t = ", "\nt = ", 1).replace("   порвано", "\nпорвано", 1).replace("   v_z", "\nv_z", 1))
            late = t >= r.t_break
            o.set_text(outcome_str(r))
            o.set_color(("#c62828" if r.outcome == "perforated" else "#1a7f37") if late else "0.25")

    class A: pass
    b = A(); b.out, b.dpi, b.fps, b.snap_dpi_factor, b.no_gif, b.snapshots = a.out, a.dpi, a.fps, a.snap_dpi_factor, False, None
    save_all(fig, update, len(t_ref), b, {})


def main():
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(encoding="utf-8")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["vel", "strain", "summary", "static", "compare"])
    ap.add_argument("--run-dir", action="append", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--spec", nargs="+", default=[])
    ap.add_argument("--what", choices=["eps", "epsp"], default="eps")
    ap.add_argument("--vmax", type=float, default=0.75)
    ap.add_argument("--vzmax", type=float, default=0.0)
    ap.add_argument("--fps", type=float, default=6.0)
    ap.add_argument("--dpi", type=int, default=65)
    ap.add_argument("--snap-dpi-factor", type=float, default=1.6)
    ap.add_argument("--snapshots", default=None)
    ap.add_argument("--zoom", default="auto")
    ap.add_argument("--lw", type=float, default=1.0)
    ap.add_argument("--yarn-stride", type=int, default=2)
    ap.add_argument("--elev", type=float, default=24.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--no-gif", action="store_true")
    a = ap.parse_args()
    if a.cmd in ("vel", "strain", "summary"):
        a.run_dir = a.run_dir[0]
    {"vel": cmd_vel, "strain": cmd_strain, "summary": cmd_summary, "static": cmd_static, "compare": cmd_compare}[a.cmd](a)


if __name__ == "__main__":
    main()
