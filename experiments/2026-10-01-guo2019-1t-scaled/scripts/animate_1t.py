#!/usr/bin/env python3
"""Анимация и снимки прогона Guo 2019, конфигурация 1T (level scaled), в стиле
experiments/2026-09-27-fiber-layers-demo/scripts/animate.py и tools/clamp_animate.py.

Каталог прогона (--run-dir) содержит
    sph_output_<NNNNN>.txt     -- частицы Solid (здесь только шарик),
    fiber_output_<NNNNN>.txt   -- бусины и связи ткани (формат src/fiber_output.cpp),
    monitor.csv                -- построчный монитор (step,t,dt,z_c,v_z,...,w_f,...),
    guo2019_fabric_layout.txt  -- раскладка нитей (направление, диапазон бусин); необязательно.

Рисуется (только реальные кадры, без интерполяции):
  * слева -- 3D: нити ломаными по связям, цвет -- eps/eps_break (фиксированная шкала
    0..1); порванные связи в кадре отсутствуют. Вертикальные нити (вдоль y, зажаты) -- толще,
    горизонтальные (вдоль x, свободны) -- тоньше. Частицы шарика -- серые точки;
  * справа -- боковое сечение (проекция y-z или x-z бусин из полосы около центральной
    плоскости) с контуром шарика: виден прогиб ткани.
Время на кадре -- время монитора для STEP кадра (кадр шага k пишется после advance, но до
time += dt, поэтому TIME в заголовке fiber_output может быть сдвинут на шаг; расхождение
печатается при загрузке). Оси фиксированы по всем кадрам.

Usage:
    python animate_1t.py --run-dir out/1T_scaled_v76.0_r0 --out anim_1t.gif
    python animate_1t.py --run-dir out/... --out anim.gif --snapshots auto --zoom 8e-3
    python animate_1t.py --run-dir out/... --out anim.gif --snapshots 0,10,20 --max-frames 21
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (регистрирует проекцию 3d)
from mpl_toolkits.mplot3d.art3d import Line3DCollection

FRAME_RE = re.compile(r"(?:sph|fiber)_output_(\d+)\.txt$")


# ---------------------------------------------------------------------------
# Чтение
# ---------------------------------------------------------------------------

def discover_frame_ids(run_dir: str) -> list[str]:
    """Номера кадров, для которых есть ОБА файла (sph и fiber)."""
    sph = {FRAME_RE.search(os.path.basename(f)).group(1)
           for f in glob.glob(os.path.join(run_dir, "sph_output_*.txt"))}
    fib = {FRAME_RE.search(os.path.basename(f)).group(1)
           for f in glob.glob(os.path.join(run_dir, "fiber_output_*.txt"))}
    return sorted(sph & fib, key=int)


def parse_sph_frame(path: str) -> np.ndarray:
    """Позиции частиц (N,3). Колонки берём по заголовку."""
    with open(path) as f:
        cols = f.readline().lstrip("# ").split()
        rows = [ln.split() for ln in f if ln.strip()]
    if not rows:
        return np.zeros((0, 3))
    ci = {n: i for i, n in enumerate(cols)}
    d = np.array(rows, dtype=float)
    return d[:, [ci["x"], ci["y"], ci["z"]]]


class FiberFrame:
    __slots__ = ("time", "step", "pos", "ij", "l0", "lbreak", "elong")


def parse_fiber_frame(path: str) -> FiberFrame:
    """Формат src/fiber_output.cpp: TIME/STEP/N_PARTICLES/N_BONDS, PARTICLES, BONDS."""
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
    pos = np.zeros((n, 3))
    for k in range(n):
        p = lines[idx + k].split()
        pos[int(p[0])] = (float(p[3]), float(p[4]), float(p[5]))
    idx += n
    while idx < len(lines) and lines[idx].lstrip().startswith("#"):
        idx += 1
    if idx < len(lines) and lines[idx].strip() == "BONDS":
        idx += 1
    rows = np.array([ln.split() for ln in lines[idx:idx + nb]], dtype=float).reshape(-1, 6)
    ff = FiberFrame()
    ff.time, ff.step = float(hd["TIME"]), int(float(hd["STEP"]))
    ff.pos = pos
    ff.ij = rows[:, :2].astype(int)
    ff.l0, ff.lbreak, ff.elong = rows[:, 2], rows[:, 4], rows[:, 5]
    return ff


def read_monitor_file(path: str) -> dict[str, np.ndarray]:
    """monitor.csv; недописанная последняя строка (прогон ещё идёт) отбрасывается."""
    with open(path) as f:
        names = f.readline().strip().split(",")
        rows = [ln.split(",") for ln in f]
    rows = [r for r in rows if len(r) == len(names)]
    d = np.array(rows, dtype=float).reshape(-1, len(names))
    return {n: d[:, i] for i, n in enumerate(names)}


def read_monitor(run_dir: str) -> dict[str, np.ndarray] | None:
    path = os.path.join(run_dir, "monitor.csv")
    return read_monitor_file(path) if os.path.exists(path) else None


def read_layout(run_dir: str, n_beads: int) -> tuple[np.ndarray, float]:
    """Направление нити для каждой бусины: 0 -- x (горизонтальная), 1 -- y (вертикальная, зажатая), -1 -- неизвестно.
    Второй результат -- шаг нитей (pitch), м, либо 0."""
    dirs = -np.ones(n_beads, dtype=int)
    path = os.path.join(run_dir, "guo2019_fabric_layout.txt")
    offs = []
    if os.path.exists(path):
        with open(path) as f:
            for ln in f:
                p = ln.split()
                if len(p) < 7 or p[0].startswith("#"):
                    continue
                a, b = int(p[5]), int(p[6])
                dirs[a:b + 1] = 0 if p[1] == "x" else 1
                if p[1] == "x":
                    offs.append(float(p[3]))
    pitch = float(np.median(np.diff(sorted(offs)))) if len(offs) > 1 else 0.0
    return dirs, pitch


def ratio_of(ff: FiberFrame) -> np.ndarray:
    """eps / eps_break в [0,1]; сжатие (elong<0) даёт 0."""
    eb = np.maximum(ff.lbreak / ff.l0 - 1.0, 1e-12)
    return np.clip(ff.elong / eb, 0.0, 1.0)


# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="каталог прогона (sph_output_*, fiber_output_*, monitor.csv)")
    ap.add_argument("--out", required=True, help="выходной .gif; снимки -- рядом: <имя>_frame<k>.png")
    ap.add_argument("--fps", type=float, default=6.0)
    ap.add_argument("--dpi", type=int, default=75)
    ap.add_argument("--max-frames", type=int, default=0,
                    help="если кадров больше -- равномерно проредить до стольких (0 = все)")
    ap.add_argument("--snapshots", default=None,
                    help="индексы кадров через запятую (в выбранной последовательности) или 'auto' "
                         "(первый, первый разрыв, максимальный прогиб, последний)")
    ap.add_argument("--zoom", type=float, default=0.0,
                    help="полуразмер окна вокруг места удара по x,y, м (0 = вся ткань)")
    ap.add_argument("--section", choices=["yz", "xz"], default="yz",
                    help="проекция бокового сечения (по умолчанию yz: вертикальные нити вдоль y)")
    ap.add_argument("--band", type=float, default=0.0,
                    help="полуширина полосы сечения по поперечной оси, м (0 = 1.5 шага нитей)")
    ap.add_argument("--zlim", default=None,
                    help="пределы по z в мм, 'lo,hi' (по умолчанию -- по всем кадрам, включая начальное положение шарика)")
    ap.add_argument("--elev", type=float, default=24.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--no-gif", action="store_true", help="только снимки")
    args = ap.parse_args()

    ids = discover_frame_ids(args.run_dir)
    if not ids:
        print(f"{args.run_dir}: нет пар sph_output/fiber_output", file=sys.stderr)
        sys.exit(1)
    if args.max_frames and len(ids) > args.max_frames:
        sel = np.linspace(0, len(ids) - 1, args.max_frames).round().astype(int)
        ids = [ids[i] for i in sorted(set(sel.tolist()))]
    nf = len(ids)
    print(f"Загрузка {nf} кадров из {args.run_dir} ...")
    ball = [parse_sph_frame(os.path.join(args.run_dir, f"sph_output_{i}.txt")) for i in ids]
    fib = [parse_fiber_frame(os.path.join(args.run_dir, f"fiber_output_{i}.txt")) for i in ids]
    mon = read_monitor(args.run_dir)
    dirs, pitch = read_layout(args.run_dir, len(fib[0].pos))
    band = args.band or (1.5 * pitch if pitch > 0 else 1.5e-3)

    # --- время и величины монитора для каждого кадра (по STEP) ---
    hdr_t = np.array([f.time for f in fib])
    times = hdr_t.copy()
    ms = None
    if mon is not None:
        ms = {}
        steps = mon["step"].astype(int)
        order = {s: r for r, s in enumerate(steps)}
        rows, exact = [], []
        for f in fib:
            r = order.get(f.step)
            exact.append(r is not None)
            rows.append(r if r is not None else int(np.argmin(np.abs(steps - f.step))))
        exact = np.array(exact)
        for key in ("t", "v_z", "w_f", "n_broken_v", "n_broken_h", "z_c"):
            if key in mon:
                ms[key] = mon[key][rows]
        if "t" in ms and exact.any():
            times = np.where(exact, ms["t"], hdr_t)
            diff = np.abs(ms["t"][exact] - hdr_t[exact]).max()
            print(f"Время: монитор для STEP кадра; max|t_монитор - TIME_заголовка| = {diff * 1e6:.4f} мкс "
                  f"(кадров без строки монитора: {int((~exact).sum())}, для них берётся TIME заголовка)")
    else:
        print("monitor.csv нет: время берётся из заголовка fiber_output", file=sys.stderr)

    # --- границы осей по всем кадрам (робастно: процентили) ---
    allp = np.concatenate([f.pos for f in fib] + ball, axis=0)
    lo = np.percentile(allp, 0.2, axis=0)
    hi = np.percentile(allp, 99.8, axis=0)
    ctr_xy = 0.5 * (fib[0].pos.min(axis=0) + fib[0].pos.max(axis=0))[:2]
    if args.zoom > 0:
        lo[:2], hi[:2] = ctr_xy - args.zoom, ctr_xy + args.zoom
    if args.zlim:
        lo[2], hi[2] = (float(v) * 1e-3 for v in args.zlim.split(","))
    pad = 0.05 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    span = hi - lo
    aspect = np.maximum(span / span[:2].max(), 0.25)

    r_ball = 0.5 * float(np.ptp(ball[0][:, 2])) if len(ball[0]) else 3e-3
    print(f"Пределы (мм): x[{lo[0]*1e3:.2f},{hi[0]*1e3:.2f}] y[{lo[1]*1e3:.2f},{hi[1]*1e3:.2f}] "
          f"z[{lo[2]*1e3:.2f},{hi[2]*1e3:.2f}]; шаг нитей {pitch*1e3:.3f} мм; полоса сечения ±{band*1e3:.2f} мм")

    # --- фигура ---
    fig = plt.figure(figsize=(13, 6.6))
    ax = fig.add_axes((0.0, 0.03, 0.55, 0.86), projection="3d")
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(tuple(aspect))
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_xlabel("x, м"); ax.set_ylabel("y, м"); ax.set_zlabel("z, м")
    ax.tick_params(labelsize=7)
    for a, axis in zip(aspect, (ax.xaxis, ax.yaxis, ax.zaxis)):
        axis.set_major_locator(MaxNLocator(nbins=5 if a > 0.5 else 3))

    c0 = fib[0]
    seg0 = np.stack([c0.pos[c0.ij[:, 0]], c0.pos[c0.ij[:, 1]]], axis=1)
    lc3 = Line3DCollection(seg0, cmap="turbo", clim=(0.0, 1.0), linewidths=1.2)
    lc3.set_array(np.zeros(len(seg0)))
    ax.add_collection3d(lc3, autolim=False)
    sc3 = ax.scatter(*ball[0].T, c="0.45", s=4, depthshade=True, label="шарик")
    ax.legend(loc="upper left", fontsize=8)
    cb = fig.colorbar(lc3, ax=ax, shrink=0.5, pad=0.0, location="left")
    cb.set_label(r"$\varepsilon / \varepsilon_{разр}$ связи")

    # боковое сечение
    ia = 1 if args.section == "yz" else 0     # горизонтальная ось сечения
    ic = 0 if args.section == "yz" else 1     # поперечная ось полосы
    axS = fig.add_axes((0.62, 0.10, 0.36, 0.78))
    axS.set_xlim(lo[ia], hi[ia]); axS.set_ylim(lo[2], hi[2])
    axS.set_aspect("equal", adjustable="box")
    axS.set_xlabel(f"{'xyz'[ia]}, м"); axS.set_ylabel("z, м")
    axS.set_title(f"Сечение {args.section}: бусины с |{'xyz'[ic]}| < {band*1e3:.2f} мм", fontsize=10)
    axS.grid(alpha=0.3)
    lcS = LineCollection([], cmap="turbo", clim=(0.0, 1.0), linewidths=1.0)
    axS.add_collection(lcS)
    dotsS, = axS.plot([], [], ".", ms=2.5, color="k")
    ballS, = axS.plot([], [], ".", ms=2, color="0.6")
    circ = Circle((0, 0), r_ball, fill=False, ec="tab:red", lw=1.2, ls="--")
    axS.add_patch(circ)
    title = fig.text(0.02, 0.965, "", fontsize=11, va="top", ha="left")

    def update(k: int):
        f = fib[k]
        seg = np.stack([f.pos[f.ij[:, 0]], f.pos[f.ij[:, 1]]], axis=1)
        r = ratio_of(f)
        lc3.set_segments(seg)
        lc3.set_array(r)
        if len(dirs) == len(f.pos) and (dirs >= 0).any():
            lc3.set_linewidths(np.where(dirs[f.ij[:, 0]] == 1, 1.6, 0.9))  # вертикальные (y) толще горизонтальных (x)
        sc3._offsets3d = tuple(ball[k].T)
        inb = np.abs(f.pos[:, ic] - ctr_xy[ic]) < band
        bi = inb[f.ij[:, 0]] & inb[f.ij[:, 1]]
        lcS.set_segments(seg[bi][:, :, [ia, 2]])
        lcS.set_array(r[bi])
        dotsS.set_data(f.pos[inb, ia], f.pos[inb, 2])
        bb = ball[k]
        if len(bb):
            bsel = np.abs(bb[:, ic] - ctr_xy[ic]) < band
            ballS.set_data(bb[bsel, ia], bb[bsel, 2])
            circ.center = (float(bb[:, ia].mean()), float(0.5 * (bb[:, 2].min() + bb[:, 2].max())))
        s = f"t = {times[k] * 1e6:.2f} мкс   шаг = {f.step}"
        if ms is not None:
            if "v_z" in ms:
                s += f"   v_z шарика = {ms['v_z'][k]:.1f} м/с"
            if "w_f" in ms:
                s += f"   прогиб w_f = {ms['w_f'][k] * 1e3:.2f} мм"
            if "n_broken_v" in ms:
                s += f"   порвано нитей верт./гориз. = {int(ms['n_broken_v'][k])}/{int(ms['n_broken_h'][k])}"
        title.set_text(s)

    # --- снимки ---
    base = os.path.splitext(args.out)[0]
    snaps: list[int] = []
    if args.snapshots == "auto":
        snaps = [0]
        if ms is not None and "n_broken_v" in ms:
            hit = np.where(ms["n_broken_v"] + ms["n_broken_h"] > 0)[0]
            if len(hit):
                snaps.append(int(hit[0]))
        if ms is not None and "w_f" in ms:
            snaps.append(int(np.argmax(ms["w_f"])))
        snaps.append(nf - 1)
        snaps = sorted(set(snaps))
    elif args.snapshots:
        snaps = [int(s) for s in args.snapshots.split(",") if s.strip()]
    for si in snaps:
        if 0 <= si < nf:
            update(si)
            p = f"{base}_frame{si}.png"
            fig.savefig(p, dpi=int(args.dpi * 1.5))
            print(f"снимок: {p}")
        else:
            warnings.warn(f"--snapshots {si} вне [0, {nf - 1}]")

    if not args.no_gif:
        anim = FuncAnimation(fig, update, frames=nf, interval=1000.0 / args.fps, blit=False)
        anim.save(args.out, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
        print(f"GIF: {args.out} ({os.path.getsize(args.out) / 1e6:.2f} МБ, {nf} кадров)")


if __name__ == "__main__":
    main()
