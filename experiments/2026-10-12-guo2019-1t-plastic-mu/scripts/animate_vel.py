#!/usr/bin/env python3
"""Анимация выстрела Guo 2019 со ШКАЛОЙ СКОРОСТИ (видно пробитие / непробитие).

Одиночный режим -- кадр из трёх панелей:
  (1) 3D нити, цвет линии = |v| бусины (фиксированная шкала 0..vmax, по умолчанию V0), шарик -- серые точки;
  (2) боковое сечение (yz/xz, полоса около центра), бусины по v_z (diverging, симметричная шкала +-vzmax), контур шарика;
  (3) график v_z шарика от t (monitor.csv), маркер текущего времени, линия v_z = 0, число текущей скорости.
Заголовок: V0, mu, t, v_z шарика, прогиб w_f, итог (shots.csv рядом с каталогом прогона: <run-dir>/../shots.csv).
Итог в заголовке всегда; цветом (зелёный -- не пробит, красный -- пробит) выделяется после первого разрыва и на последнем кадре.

--compare: несколько --run-dir (2-3); один GIF, колонка на выстрел: сечение по v_z + график v_z(t);
кадры выбираются по общему времени (по времени первого выстрела, у остальных -- ближайший кадр).

--monitor-dirs d1 d2 ... --static-out vz_compare.svg : статическая фигура v_z(t) шара для нескольких выстрелов
(PNG рядом), исход в легенде. Работает только по monitor.csv (кадры не нужны).

Usage:
  python animate_vel.py --run-dir mu18/v79.0/1T_base_v79.0_r0 --out vel_out/v79.gif --max-frames 40 --snapshots auto
  python animate_vel.py --compare --run-dir A --run-dir B --run-dir C --out vel_out/cmp.gif --max-frames 40
  python animate_vel.py --monitor-dirs A B C --static-out vel_out/vz_compare.svg
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time

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
from animate_ply import (discover_frame_ids, parse_sph_frame, read_monitor_file,  # noqa: E402
                         read_layout)

FRAME_NUM_RE = re.compile(r"_output_(\d+)\.txt$")
GREEN, RED, NEUTRAL = "#1a7f37", "#c62828", "0.25"


# ---------------------------------------------------------------------------
# Чтение
# ---------------------------------------------------------------------------

def parse_fiber_vel(path: str):
    """(time, step, pos(N,3), vel(N,3), ij(M,2)) из fiber_output."""
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]
    h = lines[0].split()
    hd = {h[i]: h[i + 1] for i in range(0, len(h) - 1, 2)}
    n, nb = int(hd["N_PARTICLES"]), int(hd["N_BONDS"])
    idx = 1
    while idx < len(lines) and lines[idx].lstrip().startswith("#"):
        idx += 1
    if lines[idx].strip() == "PARTICLES":
        idx += 1
    pr = np.array([ln.split()[:9] for ln in lines[idx:idx + n]], dtype=float)
    pos, vel = np.zeros((n, 3)), np.zeros((n, 3))
    ids = pr[:, 0].astype(int)
    pos[ids], vel[ids] = pr[:, 3:6], pr[:, 6:9]
    idx += n
    while idx < len(lines) and lines[idx].lstrip().startswith("#"):
        idx += 1
    if idx < len(lines) and lines[idx].strip() == "BONDS":
        idx += 1
    rows = np.array([ln.split()[:2] for ln in lines[idx:idx + nb]], dtype=float).reshape(-1, 2)
    return float(hd["TIME"]), int(float(hd["STEP"])), pos, vel, rows.astype(int)


def read_shot(run_dir: str) -> dict | None:
    """Строка shots.csv из <run-dir>/../shots.csv, лучше всего подходящая по V0 из имени каталога."""
    path = os.path.join(os.path.dirname(os.path.abspath(run_dir)), "shots.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    m = re.search(r"_v([0-9.]+?)(?:_rep\w*)?_r\d", os.path.basename(os.path.abspath(run_dir)))
    if m:
        v = float(m.group(1))
        rows.sort(key=lambda r: abs(float(r.get("V0", "nan")) - v))
    return rows[0]


def fnum(x, default=float("nan")):
    try:
        return float(x)
    except Exception:
        return default


def outcome_text(shot: dict | None) -> tuple[str, str | None]:
    """(текст, цвет) итога."""
    if shot is None:
        return "", None
    o = shot.get("outcome", "")
    if o == "perforated":
        return f"ПРОБИТ, V_R = {fnum(shot.get('V_R')):.1f} м/с", RED
    if o == "not_perforated":
        return "НЕ ПРОБИТ", GREEN
    return "исход не определён", None


class Run:
    pass


def load_run(run_dir: str, max_frames: int, want_times: np.ndarray | None = None) -> Run:
    r = Run()
    r.dir = run_dir
    r.name = os.path.basename(os.path.abspath(run_dir))
    r.shot = read_shot(run_dir)
    r.V0 = fnum(r.shot["V0"]) if r.shot else float("nan")
    r.mu = fnum(r.shot["mu"]) if r.shot else float("nan")
    r.otext, r.ocol = outcome_text(r.shot)
    r.t_break = fnum(r.shot.get("t_first_break")) if r.shot else float("nan")
    mp = os.path.join(run_dir, "monitor.csv")
    r.mon = read_monitor_file(mp) if os.path.exists(mp) else None
    if r.mon is not None and np.isnan(r.V0):
        r.V0 = -float(r.mon["v_z"][0])
    ids = discover_frame_ids(run_dir)
    if not ids:
        raise SystemExit(f"{run_dir}: нет пар sph_output/fiber_output")
    r.all_ids = ids
    r.ids = ids
    if want_times is None and max_frames and len(ids) > max_frames:
        sel = np.linspace(0, len(ids) - 1, max_frames).round().astype(int)
        r.ids = [ids[i] for i in sorted(set(sel.tolist()))]
    return r


def frame_time(r: Run, step: int, hdr_t: float) -> float:
    if r.mon is None:
        return hdr_t
    steps = r.mon["step"].astype(int)
    j = int(np.argmin(np.abs(steps - step)))
    return float(r.mon["t"][j]) if steps[j] == step else hdr_t


def read_frames(r: Run, ids: list[str]):
    r.fr = []
    for i in ids:
        t, st, pos, vel, ij = parse_fiber_vel(os.path.join(r.dir, f"fiber_output_{i}.txt"))
        ball = parse_sph_frame(os.path.join(r.dir, f"sph_output_{i}.txt"))
        r.fr.append(dict(t=frame_time(r, st, t), step=st, pos=pos, vel=vel, ij=ij, ball=ball))
    r.times = np.array([f["t"] for f in r.fr])


def mon_at(r: Run, t: float, key: str):
    if r.mon is None or key not in r.mon:
        return float("nan")
    return float(np.interp(t, r.mon["t"], r.mon[key]))


def auto_window(r: Run, zoom_arg: str, zlim_arg: str):
    c0 = r.fr[0]["pos"]
    ctr = 0.5 * (c0.min(axis=0) + c0.max(axis=0))[:2]
    ball0 = r.fr[0]["ball"]
    rb = 0.5 * float(np.ptp(ball0[:, 2])) if len(ball0) else 3e-3
    half = 0.5 * float((c0.max(axis=0) - c0.min(axis=0))[:2].max())
    if str(zoom_arg).lower() == "auto":
        dz = np.max([np.abs(f["pos"][:, 2] - c0[:, 2]) for f in r.fr], axis=0)
        moved = dz > 0.05 * max(dz.max(), 1e-12)
        rr = float(np.max(np.hypot(c0[moved, 0] - ctr[0], c0[moved, 1] - ctr[1]))) if moved.any() else 0.0
        zoom = float(min(max(1.25 * rr, 4.0 * rb), half))
        if zoom > 0.9 * half:
            zoom = 0.0
    else:
        zoom = float(zoom_arg)
    lo, hi = np.zeros(3), np.zeros(3)
    if zoom > 0:
        lo[:2], hi[:2] = ctr - zoom, ctr + zoom
    else:
        allxy = np.concatenate([f["pos"] for f in r.fr])
        lo[:2], hi[:2] = allxy.min(axis=0)[:2], allxy.max(axis=0)[:2]

    def inw(p):
        if zoom <= 0:
            return np.ones(len(p), bool)
        return (np.abs(p[:, 0] - ctr[0]) <= zoom) & (np.abs(p[:, 1] - ctr[1]) <= zoom)

    if zlim_arg != "auto":
        lo[2], hi[2] = (float(v) * 1e-3 for v in zlim_arg.split(","))
    else:
        zs = np.concatenate([f["pos"][inw(f["pos"]), 2] for f in r.fr] +
                            [f["ball"][inw(f["ball"]), 2] for f in r.fr if len(f["ball"])])
        lo[2], hi[2] = zs.min(), zs.max()
    r.ctr, r.rb, r.zoom, r.inw = ctr, rb, zoom, inw
    return lo, hi


def pad_box(lo, hi):
    pad = 0.06 * (hi - lo)
    return lo - pad, hi + pad


# ---------------------------------------------------------------------------
# Общие элементы рисования
# ---------------------------------------------------------------------------

class VzPlot:
    """График v_z(t) шарика с маркером текущего времени."""

    def __init__(self, ax, r: Run):
        self.r, self.ax = r, ax
        ax.axhline(0.0, color="0.3", lw=1.0)
        if r.mon is not None:
            ax.plot(r.mon["t"] * 1e6, r.mon["v_z"], color="#1f77b4", lw=1.6)
            if np.isfinite(r.t_break):
                ax.axvline(r.t_break * 1e6, color="0.6", ls=":", lw=1.0)
                ax.text(r.t_break * 1e6, 0.98, " первый разрыв", transform=ax.get_xaxis_transform(),
                        fontsize=7, color="0.4", va="top", rotation=90)
            ax.set_xlim(0, float(r.mon["t"].max()) * 1e6)
        ax.set_xlabel("t, мкс")
        ax.set_ylabel("v_z шарика, м/с")
        ax.grid(alpha=0.3)
        self.vline = ax.axvline(0, color="#d62728", lw=1.2)
        self.dot, = ax.plot([], [], "o", color="#d62728", ms=6, zorder=5)
        self.txt = ax.text(0.97, 0.06, "", transform=ax.transAxes, ha="right", va="bottom",
                           fontsize=11, fontweight="bold",
                           bbox=dict(fc="white", ec="0.7", alpha=0.9, boxstyle="round,pad=0.25"))

    def set(self, t: float):
        v = mon_at(self.r, t, "v_z")
        self.vline.set_xdata([t * 1e6, t * 1e6])
        self.dot.set_data([t * 1e6], [v])
        self.txt.set_text(f"v_z = {v:.1f} м/с")


class Section:
    """Сечение с бусинами по v_z и контуром шарика."""

    def __init__(self, ax, r: Run, lo, hi, vzmax, section, band, cmap):
        self.r, self.ax = r, ax
        self.ia = 1 if section == "yz" else 0
        self.ic = 0 if section == "yz" else 1
        self.band = band
        ax.set_xlim(lo[self.ia], hi[self.ia])
        ax.set_ylim(lo[2], hi[2])
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"{'xyz'[self.ia]}, м")
        ax.set_ylabel("z, м")
        ax.grid(alpha=0.3)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        self.sc = ax.scatter([], [], s=11, c=[], cmap=cmap, vmin=-vzmax, vmax=vzmax, zorder=3)
        self.ballp, = ax.plot([], [], ".", ms=2, color="0.6")
        self.circ = Circle((0, 0), r.rb, fill=False, ec="0.15", lw=1.2, ls="--")
        ax.add_patch(self.circ)

    def set(self, f):
        inb = np.abs(f["pos"][:, self.ic] - self.r.ctr[self.ic]) < self.band
        self.sc.set_offsets(np.column_stack([f["pos"][inb, self.ia], f["pos"][inb, 2]]))
        self.sc.set_array(f["vel"][inb, 2])
        bb = f["ball"]
        if len(bb):
            sel = np.abs(bb[:, self.ic] - self.r.ctr[self.ic]) < self.band
            self.ballp.set_data(bb[sel, self.ia], bb[sel, 2])
            self.circ.center = (float(bb[:, self.ia].mean()), float(0.5 * (bb[:, 2].min() + bb[:, 2].max())))


def get_band(r: Run, a) -> float:
    if a.band:
        return a.band
    _, _, _, pitch = read_layout(r.dir, len(r.fr[0]["pos"]))
    return 1.5 * pitch if pitch > 0 else 1.5e-3


def status(r: Run, k: int) -> tuple[str, str]:
    """Итог для кадра k и цвет: нейтральный до первого разрыва, цветной после и в конце."""
    t = r.times[k]
    late = (k == len(r.fr) - 1) or (np.isfinite(r.t_break) and t >= r.t_break)
    return r.otext, (r.ocol if (late and r.ocol) else NEUTRAL)


def head(r: Run, k: int) -> str:
    t = r.times[k]
    s = f"пластично-вязкая нить, σ_u = 3.0 ГПа, V0 = {r.V0:g} м/с   mu = {r.mu:g}   t = {t * 1e6:.1f} мкс   v_z шарика = {mon_at(r, t, 'v_z'):.1f} м/с"
    w = mon_at(r, t, "w_f")
    if np.isfinite(w):
        s += f"   w_f = {w * 1e3:.2f} мм"
    return s


# ---------------------------------------------------------------------------
# Режимы
# ---------------------------------------------------------------------------

def save_outputs(fig, update, nf, a, snaps_default_fn):
    base = os.path.splitext(a.out)[0]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    snaps = []
    if a.snapshots == "auto":
        snaps = snaps_default_fn()
    elif a.snapshots:
        snaps = [int(s) for s in a.snapshots.split(",") if s.strip()]
    for si in snaps:
        if 0 <= si < nf:
            update(si)
            p = f"{base}_frame{si}.png"
            fig.savefig(p, dpi=int(a.dpi * a.snap_dpi_factor))
            print("снимок:", p)
    if not a.no_gif:
        t0 = time.time()
        FuncAnimation(fig, update, frames=nf, interval=1000.0 / a.fps, blit=False).save(
            a.out, writer=PillowWriter(fps=a.fps), dpi=a.dpi)
        print(f"GIF: {a.out} ({os.path.getsize(a.out) / 1e6:.2f} МБ, {nf} кадров, {time.time() - t0:.1f} с)")


def run_single(a):
    r = load_run(a.run_dir[0], a.max_frames)
    print(f"Загрузка {len(r.ids)} кадров из {r.dir} ...")
    read_frames(r, r.ids)
    nf = len(r.fr)
    lo, hi = pad_box(*auto_window(r, a.zoom, a.zlim))
    span = hi - lo
    aspect = np.maximum(span / span[:2].max(), 0.25)
    vmax = a.vmax or (r.V0 if np.isfinite(r.V0) else 100.0)
    vzmax = a.vzmax or vmax
    band = get_band(r, a)
    print(f"V0={r.V0} mu={r.mu} итог: {r.otext or '-'}; шкала |v| 0..{vmax:g}, v_z +-{vzmax:g}; кадров {nf}")

    fig = plt.figure(figsize=(17, 5.8))
    ax = fig.add_axes((0.0, 0.0, 0.34, 0.86), projection="3d")
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(tuple(aspect))
    ax.view_init(elev=a.elev, azim=a.azim)
    ax.set_xlabel("x, м"); ax.set_ylabel("y, м"); ax.set_zlabel("z, м")
    ax.tick_params(labelsize=7)
    for asp, axis in zip(aspect, (ax.xaxis, ax.yaxis, ax.zaxis)):
        axis.set_major_locator(MaxNLocator(nbins=5 if asp > 0.5 else 3))
    ax.set_title("нити: цвет = |v|", fontsize=10, y=0.98)
    lc3 = Line3DCollection(np.zeros((0, 2, 3)), cmap=a.cmap, clim=(0.0, vmax), linewidths=1.3 * a.lw)
    ax.add_collection3d(lc3, autolim=False)
    sc3 = ax.scatter(*r.fr[0]["ball"].T, c="0.35", s=3, alpha=a.ball_alpha, depthshade=True)
    cb1 = fig.colorbar(lc3, ax=ax, shrink=0.55, pad=0.0, location="left")
    cb1.set_label("|v| бусины, м/с")

    axS = fig.add_axes((0.40, 0.11, 0.26, 0.70))
    sec = Section(axS, r, lo, hi, vzmax, a.section, band, a.vz_cmap)
    axS.set_title(f"сечение {a.section}, полоса ±{band * 1e3:.1f} мм: цвет = v_z", fontsize=10)
    cb2 = fig.colorbar(sec.sc, ax=axS, shrink=0.8, pad=0.02)
    cb2.set_label("v_z бусины, м/с")

    axV = fig.add_axes((0.755, 0.13, 0.225, 0.68))
    vzp = VzPlot(axV, r)
    axV.set_title("v_z шарика(t)", fontsize=10)
    t1 = fig.text(0.01, 0.985, "", fontsize=11, va="top", ha="left")
    t2 = fig.text(0.01, 0.91, "", fontsize=15, va="top", ha="left", fontweight="bold")

    def update(k):
        f = r.fr[k]
        pos, ij = f["pos"], f["ij"]
        seg = np.stack([pos[ij[:, 0]], pos[ij[:, 1]]], axis=1)
        sp = np.linalg.norm(f["vel"], axis=1)
        keep = np.ones(len(ij), bool)
        if r.zoom > 0 and not a.no_clip:
            w = r.inw(pos)
            keep = w[ij[:, 0]] & w[ij[:, 1]]
        lc3.set_segments(seg[keep])
        lc3.set_array(0.5 * (sp[ij[keep, 0]] + sp[ij[keep, 1]]))
        sc3._offsets3d = tuple(f["ball"].T)
        sec.set(f)
        vzp.set(r.times[k])
        t1.set_text(head(r, k))
        s, c = status(r, k)
        t2.set_text("Итог: " + s if s else "")
        t2.set_color(c)

    def defaults():
        s = [0]
        if r.mon is not None and "n_broken_v" in r.mon:
            hit = np.where(r.mon["n_broken_v"] + r.mon["n_broken_h"] > 0)[0]
            if len(hit):
                s.append(int(np.argmin(np.abs(r.times - r.mon["t"][hit[0]]))))
        s.append(nf // 2)
        s.append(nf - 1)
        return sorted(set(s))

    save_outputs(fig, update, nf, a, defaults)


def run_compare(a):
    runs = [load_run(d, a.max_frames, want_times=np.zeros(1)) for d in a.run_dir]
    # опорные кадры: из первого выстрела (с учётом --max-frames), времена по монитору
    r0 = runs[0]
    ids0 = r0.all_ids
    if a.max_frames and len(ids0) > a.max_frames:
        sel = np.linspace(0, len(ids0) - 1, a.max_frames).round().astype(int)
        ids0 = [ids0[i] for i in sorted(set(sel.tolist()))]
    # времена всех кадров каждого выстрела (дёшево: только заголовок fiber_output)
    def hdr_times(r):
        out = []
        for i in r.all_ids:
            with open(os.path.join(r.dir, f"fiber_output_{i}.txt")) as f:
                h = f.readline().split()
            hd = {h[j]: h[j + 1] for j in range(0, len(h) - 1, 2)}
            out.append(frame_time(r, int(float(hd["STEP"])), float(hd["TIME"])))
        return np.array(out)
    for r in runs:
        r.all_t = hdr_times(r)
    t_ref = r0.all_t[[r0.all_ids.index(i) for i in ids0]]
    for r in runs:
        pick = [int(np.argmin(np.abs(r.all_t - t))) for t in t_ref]
        read_frames(r, [r.all_ids[p] for p in pick])
        print(f"{r.name}: {len(set(pick))} различных кадров из {len(r.all_ids)}")
    nf = len(t_ref)
    wins = [auto_window(r, a.zoom, a.zlim) for r in runs]
    lo = np.min([w[0] for w in wins], axis=0)
    hi = np.max([w[1] for w in wins], axis=0)
    lo, hi = pad_box(lo, hi)
    n = len(runs)
    vzmax = a.vzmax or max(r.V0 for r in runs if np.isfinite(r.V0))
    fig = plt.figure(figsize=(5.6 * n + 0.8, 8.4))
    gs = fig.add_gridspec(2, n, left=0.06, right=0.90, top=0.84, bottom=0.07, hspace=0.28, wspace=0.28,
                          height_ratios=[1.35, 1.0])
    secs, vzps, heads, outs = [], [], [], []
    for j, r in enumerate(runs):
        axS = fig.add_subplot(gs[0, j])
        secs.append(Section(axS, r, lo, hi, vzmax, a.section, get_band(r, a), a.vz_cmap))
        axS.set_title(r.name, fontsize=9)
        axV = fig.add_subplot(gs[1, j])
        vzps.append(VzPlot(axV, r))
        pos = axS.get_position()
        heads.append(fig.text(pos.x0, 0.955, "", fontsize=10, va="top"))
        outs.append(fig.text(pos.x0, 0.905, "", fontsize=14, va="top", fontweight="bold"))
    cax = fig.add_axes((0.92, 0.45, 0.012, 0.35))
    cb = fig.colorbar(secs[0].sc, cax=cax)
    cb.set_label("v_z бусины, м/с")
    ttl = fig.text(0.5, 0.995, "", ha="center", va="top", fontsize=11, color="0.3")

    def update(k):
        ttl.set_text(f"t = {t_ref[k] * 1e6:.1f} мкс (кадры ближайшие по времени)")
        for r, s, v, h, o in zip(runs, secs, vzps, heads, outs):
            f = r.fr[k]
            s.set(f)
            v.set(r.times[k])
            h.set_text(head(r, k).replace("   t = ", "\nt = ", 1))
            txt, c = status(r, k)
            o.set_text(txt)
            o.set_color(c)

    save_outputs(fig, update, nf, a, lambda: [0, nf // 2, nf - 1])


def run_static(a):
    fig, ax = plt.subplots(figsize=(9, 5.2))
    cols = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b"]
    for j, d in enumerate(a.monitor_dirs):
        mp = os.path.join(d, "monitor.csv")
        if not os.path.exists(mp):
            print("нет monitor.csv:", d, file=sys.stderr)
            continue
        mon = read_monitor_file(mp)
        shot = read_shot(d)
        txt, _ = outcome_text(shot)
        v0 = fnum(shot["V0"]) if shot else -float(mon["v_z"][0])
        mu = fnum(shot["mu"]) if shot else float("nan")
        lab = f"пласт. нить, V0 = {v0:g} м/с" + (f", mu = {mu:g}" if np.isfinite(mu) else "") + (f": {txt}" if txt else "")
        if "_rep" in os.path.basename(os.path.abspath(d)) or "_rep" in os.path.basename(os.path.dirname(os.path.abspath(d))):
            lab += " (повтор)"
        ax.plot(mon["t"] * 1e6, mon["v_z"], lw=1.8, color=cols[j % len(cols)], label=lab,
                ls="--" if "(повтор)" in lab else "-")
    ax.axhline(0, color="0.3", lw=1)
    ax.set_xlabel("t, мкс")
    ax.set_ylabel("v_z шарика, м/с")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="best")
    ax.set_title("Скорость шарика v_z(t), Guo 2019")
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.static_out)), exist_ok=True)
    fig.savefig(a.static_out)
    png = os.path.splitext(a.static_out)[0] + ".png"
    fig.savefig(png, dpi=150)
    print("записано:", a.static_out, png)


def main():
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(encoding="utf-8")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", action="append", default=[], help="каталог прогона (повторять 2-3 раза с --compare)")
    ap.add_argument("--out", help="выходной .gif; снимки рядом: <имя>_frame<k>.png")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--monitor-dirs", nargs="+", help="каталоги с monitor.csv для статической фигуры")
    ap.add_argument("--static-out", default=None, help="vz_compare.svg (PNG рядом)")
    ap.add_argument("--fps", type=float, default=6.0)
    ap.add_argument("--dpi", type=int, default=65)
    ap.add_argument("--snap-dpi-factor", type=float, default=1.6)
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--snapshots", default=None, help="индексы кадров или 'auto'")
    ap.add_argument("--zoom", default="auto")
    ap.add_argument("--section", choices=["yz", "xz"], default="yz")
    ap.add_argument("--band", type=float, default=0.0)
    ap.add_argument("--zlim", default="auto")
    ap.add_argument("--no-clip", action="store_true")
    ap.add_argument("--lw", type=float, default=1.0)
    ap.add_argument("--vmax", type=float, default=0.0, help="верх шкалы |v|, м/с (0 = V0)")
    ap.add_argument("--vzmax", type=float, default=0.0, help="симметричная шкала v_z +-, м/с (0 = vmax / V0)")
    ap.add_argument("--cmap", default="turbo")
    ap.add_argument("--vz-cmap", default="RdBu_r")
    ap.add_argument("--ball-alpha", type=float, default=0.35)
    ap.add_argument("--elev", type=float, default=24.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--no-gif", action="store_true")
    a = ap.parse_args()
    if a.monitor_dirs or a.static_out:
        if not (a.monitor_dirs and a.static_out):
            ap.error("--monitor-dirs и --static-out нужны вместе")
        run_static(a)
        if not a.run_dir:
            return
    if not a.run_dir or not a.out:
        ap.error("нужны --run-dir и --out")
    if a.compare:
        run_compare(a)
    else:
        run_single(a)


if __name__ == "__main__":
    main()
