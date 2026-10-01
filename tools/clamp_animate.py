#!/usr/bin/env python3
"""Анимация e2e-прогонов закрепления (tests/test_e2e_clamp.cpp) в стиле демо нитей
(experiments/2026-09-27-fiber-layers-demo/scripts/animate.py).

Кадры пишет сам тест:
    test_e2e_clamp --case=ES1 --frames-dir=frames [--frames=80]
В frames/<кейс>_<решатель>[_<вариант>]_<тег прогона>/ лежат
    sph_output_<NNNNN>.txt, fiber_output_<NNNNN>.txt  -- кадры боевого цикла,
    clamp_mask.txt                                    -- "solid|fiber <индекс> <маска>",
    reaction.csv                                      -- step,t,Rs_x..Rf_z по шагам.

Рисуется:
  * слева -- 3D-вид: частицы Solid (цвет -- --color-by), бонды нитей (цвет --
    eps/eps_break, как в демо), бусины. Зажатые частицы -- чёрные квадраты
    (маска xyz) или оранжевые треугольники (частичная маска). Смещения от
    первого кадра умножаются на --exaggerate (auto: максимальное смещение =
    8 % размера сцены), иначе деформации ~1e-3 не видны;
  * справа сверху -- реакция опоры R(t) = R_solid + R_fiber по осям и курсор
    текущего кадра;
  * справа снизу -- профиль вдоль длинной оси сцены: поле --color-by у Solid,
    eps/eps_break связей и поперечное смещение у нитей.

Usage:
    python tools/clamp_animate.py --run-dir frames/ES1_pc_es1 --out es1.gif
    python tools/clamp_animate.py --run-dir frames/EF2_pc_ef2 --out ef2.gif --snapshot 0,20,40
    python tools/clamp_animate.py --all frames --out-dir anim   # все прогоны каталога
"""
from __future__ import annotations

import argparse
import csv
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
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)
from mpl_toolkits.mplot3d.art3d import Line3DCollection


# ---------------------------------------------------------------------------
# Чтение кадров
# ---------------------------------------------------------------------------

FRAME_RE = re.compile(r"(?:sph|fiber)_output_(\d{5})\.txt$")
SPH_FIELDS = ["speed", "p", "sxx", "vx", "vy", "vz", "rho"]


def discover_frame_ids(run_dir: str) -> list[str]:
    ids = set()
    for f in glob.glob(os.path.join(run_dir, "*_output_*.txt")):
        m = FRAME_RE.search(os.path.basename(f))
        if m:
            ids.add(m.group(1))
    return sorted(ids)


class SphFrame:
    __slots__ = ("pos", "fields")

    def __init__(self, pos: np.ndarray, fields: dict[str, np.ndarray]):
        self.pos = pos
        self.fields = fields


def parse_sph_frame(path: str) -> SphFrame | None:
    if not os.path.exists(path):
        return None
    with open(path, "r") as f:
        cols = f.readline().split()
        rows = [line.split() for line in f if line.strip()]
    col = {name: i for i, name in enumerate(cols)}
    if not rows:
        return SphFrame(np.zeros((0, 3)), {k: np.zeros(0) for k in SPH_FIELDS})
    data = np.array(rows, dtype=float)
    pos = np.stack([data[:, col["x"]], data[:, col["y"]], data[:, col["z"]]], axis=1)
    v = np.stack([data[:, col["vx"]], data[:, col["vy"]], data[:, col["vz"]]], axis=1)
    fields = {
        "speed": np.linalg.norm(v, axis=1),
        "vx": v[:, 0], "vy": v[:, 1], "vz": v[:, 2],
        "p": data[:, col["p"]],
        "sxx": data[:, col["sxx"]] - data[:, col["p"]],  # полное напряжение sigma_xx = s_xx - p
        "rho": data[:, col["rho"]],
    }
    return SphFrame(pos, fields)


class FiberFrame:
    __slots__ = ("time", "step", "bead_pos", "bonds_ij", "bonds_l0", "bonds_lbreak", "bonds_elong")

    def __init__(self, time, step, bead_pos, bonds_ij, bonds_l0, bonds_lbreak, bonds_elong):
        self.time = time
        self.step = step
        self.bead_pos = bead_pos
        self.bonds_ij = bonds_ij
        self.bonds_l0 = bonds_l0
        self.bonds_lbreak = bonds_lbreak
        self.bonds_elong = bonds_elong


def parse_fiber_frame(path: str) -> FiberFrame | None:
    """Формат -- src/fiber_output.cpp (тот же разбор, что в демо)."""
    if not os.path.exists(path):
        return None
    with open(path, "r") as f:
        lines = [line.rstrip("\n") for line in f]
    header = lines[0].split()
    h = {header[i]: header[i + 1] for i in range(0, len(header) - 1, 2)}
    time, step = float(h["TIME"]), int(float(h["STEP"]))
    n, nb = int(h["N_PARTICLES"]), int(h["N_BONDS"])
    idx = 1
    while idx < len(lines) and lines[idx].strip().startswith("#"):
        idx += 1
    if idx < len(lines) and lines[idx].strip() == "PARTICLES":
        idx += 1
    pos = np.zeros((n, 3))
    for k in range(n):
        p = lines[idx + k].split()
        pos[int(p[0])] = [float(p[3]), float(p[4]), float(p[5])]
    idx += n
    while idx < len(lines) and lines[idx].strip().startswith("#"):
        idx += 1
    if idx < len(lines) and lines[idx].strip() == "BONDS":
        idx += 1
    ij, l0, lb, el = [], np.zeros(nb), np.zeros(nb), np.zeros(nb)
    for k in range(nb):
        p = lines[idx + k].split()
        ij.append((int(p[0]), int(p[1])))
        l0[k], lb[k], el[k] = float(p[2]), float(p[4]), float(p[5])
    return FiberFrame(time, step, pos, ij, l0, lb, el)


def read_masks(run_dir: str) -> tuple[dict[int, int], dict[int, int]]:
    solid, fib = {}, {}
    path = os.path.join(run_dir, "clamp_mask.txt")
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                p = line.split()
                if len(p) == 3 and not p[0].startswith("#"):
                    (solid if p[0] == "solid" else fib)[int(p[1])] = int(p[2])
    return solid, fib


def read_reaction(run_dir: str) -> dict[str, np.ndarray] | None:
    path = os.path.join(run_dir, "reaction.csv")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


# ---------------------------------------------------------------------------
# Вспомогательное
# ---------------------------------------------------------------------------

def format_time(t: float, t_max: float) -> str:
    if t_max < 1e-3:
        return f"{t * 1e6:.2f} мкс"
    if t_max < 1.0:
        return f"{t * 1e3:.3f} мс"
    return f"{t:.4g} с"


def mask_label(mk: int) -> str:
    return "".join(c for b, c in enumerate("xyz") if mk & (1 << b))


def render_run(run_dir: str, out: str, args) -> None:
    ids = discover_frame_ids(run_dir)
    if not ids:
        print(f"{run_dir}: кадров нет", file=sys.stderr)
        return
    if args.max_frames and len(ids) > args.max_frames:
        sel = np.linspace(0, len(ids) - 1, args.max_frames).round().astype(int)
        ids = [ids[i] for i in sorted(set(sel.tolist()))]

    sph = [parse_sph_frame(os.path.join(run_dir, f"sph_output_{i}.txt")) for i in ids]
    fib = [parse_fiber_frame(os.path.join(run_dir, f"fiber_output_{i}.txt")) for i in ids]
    has_sph = any(s is not None and len(s.pos) > 0 for s in sph)
    has_fib = any(f is not None and len(f.bead_pos) > 0 for f in fib)
    smask, fmask = read_masks(run_dir)
    react = read_reaction(run_dir)

    # Время кадра: из заголовка fiber_output, иначе по reaction.csv. Кадр с
    # номером шага k пишется ДО прибавления dt, т. е. в момент t после шага k-1.
    def frame_time(k: int) -> float:
        if fib[k] is not None:
            return fib[k].time
        step = int(ids[k])
        if react is not None and step >= 2:
            j = min(step - 2, len(react["t"]) - 1)
            return float(react["t"][j])
        return 0.0

    times = [frame_time(k) for k in range(len(ids))]
    t_max = max(max(times), float(react["t"][-1]) if react is not None else 0.0, 1e-30)
    name = os.path.basename(os.path.normpath(run_dir))

    # --- опорная конфигурация и усиление смещений ---
    def positions(k):
        parts = []
        if has_sph and sph[k] is not None:
            parts.append(sph[k].pos)
        if has_fib and fib[k] is not None:
            parts.append(fib[k].bead_pos)
        return np.concatenate(parts, axis=0) if parts else np.zeros((0, 3))

    ref = positions(0)
    lo0, hi0 = ref.min(axis=0), ref.max(axis=0)
    extent = max(float((hi0 - lo0).max()), 1e-9)
    max_disp = 0.0
    for k in range(len(ids)):
        P = positions(k)
        if P.shape == ref.shape:
            max_disp = max(max_disp, float(np.abs(P - ref).max()))
    if args.exaggerate == "auto":
        K = 1.0 if max_disp <= 0 else float(np.clip(0.08 * extent / max_disp, 1.0, 1e5))
        K = float(10 ** np.floor(np.log10(K))) if K >= 10 else (1.0 if K < 2 else float(np.floor(K)))
    else:
        K = float(args.exaggerate)

    def warp(P, R):
        return R + K * (P - R) if P.shape == R.shape else P

    sph_ref = sph[0].pos if has_sph else None
    fib_ref = fib[0].bead_pos if has_fib else None

    all_pts = []
    for k in range(len(ids)):
        if has_sph:
            all_pts.append(warp(sph[k].pos, sph_ref))
        if has_fib:
            all_pts.append(warp(fib[k].bead_pos, fib_ref))
    all_pts = np.concatenate(all_pts, axis=0)
    lo = np.percentile(all_pts, 0.5, axis=0)
    hi = np.percentile(all_pts, 99.5, axis=0)
    span = np.maximum(hi - lo, 1e-9)
    pad = 0.08 * span.max()
    lo, hi = lo - pad, hi + pad
    aspect = np.maximum((hi - lo) / (hi - lo).max(), 0.18)

    # Профиль нитей строится по бусинам, зажатым НЕ полностью (полностью
    # зажатая нить -- опора, её профиль тождественно ноль и только мешает,
    # как нить B в EF4). Длинная ось профиля -- по этим бусинам, если Solid
    # нет, иначе по всей сцене.
    prof_beads = None
    if has_fib:
        nb0 = len(fib[0].bead_pos)
        prof_beads = np.array([fmask.get(k, 0) != 7 for k in range(nb0)])
        if not prof_beads.any():
            prof_beads[:] = True
    if has_fib and not has_sph:
        P = fib[0].bead_pos[prof_beads]
        lo0, hi0 = P.min(axis=0), P.max(axis=0)
    long_axis = int(np.argmax(hi0 - lo0))
    axis_name = "xyz"[long_axis]

    # --- шкалы цвета ---
    field = args.color_by
    if has_sph:
        vals = np.concatenate([s.fields[field] for s in sph if s is not None and len(s.pos)])
        vmin, vmax = np.percentile(vals, [1, 99])
        if field in ("p", "sxx", "vx", "vy", "vz") and vmin < 0 < vmax:
            a = max(-vmin, vmax)
            vmin, vmax, cmap = -a, a, "coolwarm"
        else:
            cmap = "viridis"
        if vmin == vmax:
            vmin, vmax = vmin - 1.0, vmax + 1.0
    field_label = {"speed": "|v|, м/с", "p": "p, Па", "sxx": "σ_xx = s_xx − p, Па", "vx": "v_x, м/с",
                   "vy": "v_y, м/с", "vz": "v_z, м/с", "rho": "ρ, кг/м³"}[field]

    # --- фигура ---
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_axes((0.0, 0.03, 0.56, 0.86), projection="3d")
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(tuple(aspect))
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_xlabel("x, м"); ax.set_ylabel("y, м"); ax.set_zlabel("z, м")
    from matplotlib.ticker import MaxNLocator
    for a, axis in zip(aspect, (ax.xaxis, ax.yaxis, ax.zaxis)):
        axis.set_major_locator(MaxNLocator(nbins=6 if a > 0.5 else 2))
    ax.tick_params(labelsize=7)
    axR = fig.add_axes((0.64, 0.56, 0.33, 0.33))
    axP = fig.add_axes((0.64, 0.09, 0.33, 0.33))
    title = fig.text(0.02, 0.965, "", fontsize=11, va="top", ha="left")

    s_free = s_full = s_part = None
    if has_sph:
        n = len(sph[0].pos)
        clamped = np.array([i in smask for i in range(n)])
        full = np.array([smask.get(i, 0) == 7 for i in range(n)])
        free_idx = np.where(~clamped)[0]
        full_idx = np.where(full)[0]
        part_idx = np.where(clamped & ~full)[0]
        P0 = sph[0].pos
        s_free = ax.scatter(P0[free_idx, 0], P0[free_idx, 1], P0[free_idx, 2], c=sph[0].fields[field][free_idx],
                            cmap=cmap, vmin=vmin, vmax=vmax, s=args.size, depthshade=True, label="Solid")
        s_full = ax.scatter(P0[full_idx, 0], P0[full_idx, 1], P0[full_idx, 2], c=sph[0].fields[field][full_idx],
                            cmap=cmap, vmin=vmin, vmax=vmax, s=args.size * 1.4, marker="s", edgecolors="k",
                            linewidths=0.6, depthshade=False, label="Solid, зажат xyz")
        if len(part_idx):
            lbl = sorted({mask_label(smask[i]) for i in part_idx})
            s_part = ax.scatter(P0[part_idx, 0], P0[part_idx, 1], P0[part_idx, 2], c=sph[0].fields[field][part_idx],
                                cmap=cmap, vmin=vmin, vmax=vmax, s=args.size * 1.4, marker="^",
                                edgecolors="darkorange", linewidths=0.6, depthshade=False,
                                label="Solid, зажат " + "/".join(lbl))
        cb = fig.colorbar(s_free, ax=ax, shrink=0.5, pad=0.0, location="left")
        cb.set_label(field_label)

    bond_lc = b_free = b_full = b_part = None
    if has_fib:
        nb = len(fib[0].bead_pos)
        fclamped = np.array([k in fmask for k in range(nb)])
        ffull = np.array([fmask.get(k, 0) == 7 for k in range(nb)])
        bf_idx = np.where(~fclamped)[0]
        bfull_idx = np.where(ffull)[0]
        bpart_idx = np.where(fclamped & ~ffull)[0]
        B0 = fib[0].bead_pos
        # Начальный сегмент-заглушка: add_collection3d не принимает пустую коллекцию.
        bond_lc = Line3DCollection([np.stack([B0[0], B0[0]])], linewidths=2.0, cmap="RdYlGn_r", clim=(0.0, 1.0))
        bond_lc.set_array(np.zeros(1))
        ax.add_collection3d(bond_lc, autolim=False)
        b_free = ax.scatter(B0[bf_idx, 0], B0[bf_idx, 1], B0[bf_idx, 2], c="tab:blue", s=6, depthshade=False,
                            label="бусины")
        b_full = ax.scatter(B0[bfull_idx, 0], B0[bfull_idx, 1], B0[bfull_idx, 2], c="k", marker="s", s=28,
                            depthshade=False, label="бусины, зажаты xyz")
        if len(bpart_idx):
            lbl = sorted({mask_label(fmask[k]) for k in bpart_idx})
            b_part = ax.scatter(B0[bpart_idx, 0], B0[bpart_idx, 1], B0[bpart_idx, 2], c="darkorange", marker="^",
                                s=36, depthshade=False, label="бусины, зажаты " + "/".join(lbl))
        cb2 = fig.colorbar(bond_lc, ax=ax, shrink=0.5, pad=0.06)
        cb2.set_label(r"$\varepsilon / \varepsilon_{break}$ связи")
    ax.legend(loc="upper right", fontsize=8)

    # --- реакция ---
    cursor = None
    if react is not None:
        tR = react["t"]
        for d, c in zip("xyz", ("tab:red", "tab:green", "tab:blue")):
            R = react[f"Rs_{d}"] + react[f"Rf_{d}"]
            if np.any(R != 0.0):
                axR.plot(tR, R, color=c, lw=1.0, label=f"R_{d}")
        axR.axhline(0.0, color="0.6", lw=0.6)
        cursor = axR.axvline(0.0, color="k", lw=1.0, ls="--")
        axR.set_xlim(0, t_max)
        axR.set_xlabel("t, с"); axR.set_ylabel("R = R_solid + R_fiber, Н")
        axR.set_title("Реакция опоры", fontsize=10)
        if axR.get_legend_handles_labels()[0]:
            axR.legend(fontsize=8, loc="best")
        axR.ticklabel_format(axis="both", style="sci", scilimits=(-3, 3))
        axR.grid(alpha=0.3)
    else:
        axR.text(0.5, 0.5, "reaction.csv нет", ha="center", va="center", transform=axR.transAxes)

    # --- профиль: границы по всем кадрам ---
    prof_lines = {}
    if has_sph:
        allv = np.concatenate([s.fields[field] for s in sph])
        prof_lines["sph"], = axP.plot([], [], ".", ms=2.5, color="tab:purple", label=f"Solid: {field_label}")
        if len(smask):
            prof_lines["sph_c"], = axP.plot([], [], "s", ms=3, mfc="none", mec="k", label="зажатые")
        axP.set_ylim(*(np.percentile(allv, [0.5, 99.5]) + np.array([-1, 1]) * 0.05 * (np.ptp(allv) + 1e-30)))
        axP.set_ylabel(field_label)
    axP2 = None
    if has_fib:
        axP2 = axP.twinx() if has_sph else axP
        prof_lines["eps"], = axP2.plot([], [], "-", color="tab:red", lw=1.2, label="ε/ε_break связей")
        trans = [a for a in range(3) if a != long_axis]
        prof_lines["lat"], = axP2.plot([], [], "-", color="tab:blue", lw=1.0)
        ratios = []
        for f in fib:
            if f is not None and len(f.bonds_l0):
                eb = np.maximum(f.bonds_lbreak / f.bonds_l0 - 1.0, 1e-12)
                ratios.append(f.bonds_elong / eb)
        rr = np.concatenate(ratios) if ratios else np.zeros(1)
        # Поперечное смещение -- по той поперечной оси, где оно больше всего за прогон.
        lat_by_axis = {a: max(float(np.abs((f.bead_pos - fib_ref)[prof_beads][:, a]).max()) for f in fib)
                       for a in trans}
        lat_axis = max(lat_by_axis, key=lat_by_axis.get)
        lat_max = lat_by_axis[lat_axis]
        prof_lines["lat"].set_label(f"смещение по {'xyz'[lat_axis]} / max ({lat_max:.2e} м)")
        # Снизу -- не ниже -1.5: провисшие связи (сжатие, сила 0) уходят в большие
        # отрицательные eps и сжимают полезную часть графика.
        axP2.set_ylim(max(-1.5, min(-1.05, rr.min() * 1.05)), max(1.05, rr.max() * 1.05))
        axP2.axhline(0, color="0.6", lw=0.6)
        axP2.set_ylabel("ε/ε_break, смещение/max")
    axP.set_xlim(lo0[long_axis] - 0.03 * extent, hi0[long_axis] + 0.03 * extent)
    axP.set_xlabel(f"{axis_name} (начальная конфигурация), м")
    axP.set_title(f"Профиль вдоль {axis_name}", fontsize=10)
    h1, l1 = axP.get_legend_handles_labels()
    h2, l2 = (axP2.get_legend_handles_labels() if axP2 is not None and axP2 is not axP else ([], []))
    axP.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    axP.grid(alpha=0.3)

    first_bonds = set(fib[0].bonds_ij) if has_fib else set()

    def update(k: int):
        artists = []
        if has_sph:
            S = sph[k]
            W = warp(S.pos, sph_ref)
            vals = S.fields[field]
            for sc, idx in ((s_free, free_idx), (s_full, full_idx), (s_part, part_idx)):
                if sc is None:
                    continue
                sc._offsets3d = (W[idx, 0], W[idx, 1], W[idx, 2])
                sc.set_array(vals[idx])
                artists.append(sc)
            prof_lines["sph"].set_data(sph_ref[:, long_axis], vals)
            if "sph_c" in prof_lines:
                ci = np.where(clamped)[0]
                prof_lines["sph_c"].set_data(sph_ref[ci, long_axis], vals[ci])
        n_broken = 0
        if has_fib:
            F = fib[k]
            W = warp(F.bead_pos, fib_ref)
            if F.bonds_ij:
                ii = np.array([i for i, _ in F.bonds_ij]); jj = np.array([j for _, j in F.bonds_ij])
                segs = np.stack([W[ii], W[jj]], axis=1)
                eb = np.maximum(F.bonds_lbreak / F.bonds_l0 - 1.0, 1e-12)
                ratio = F.bonds_elong / eb
                bond_lc.set_segments(segs)
                bond_lc.set_array(np.clip(ratio, 0.0, 1.0))
                sel = prof_beads[ii] & prof_beads[jj]
                mid = 0.5 * (fib_ref[ii, long_axis] + fib_ref[jj, long_axis])[sel]
                order = np.argsort(mid)
                prof_lines["eps"].set_data(mid[order], ratio[sel][order])
            else:
                bond_lc.set_segments([])
                prof_lines["eps"].set_data([], [])
            for sc, idx in ((b_free, bf_idx), (b_full, bfull_idx), (b_part, bpart_idx)):
                if sc is None:
                    continue
                sc._offsets3d = (W[idx, 0], W[idx, 1], W[idx, 2])
                artists.append(sc)
            lat = F.bead_pos[:, lat_axis] - fib_ref[:, lat_axis]
            pb = np.where(prof_beads)[0]
            o = pb[np.argsort(fib_ref[pb, long_axis])]
            prof_lines["lat"].set_data(fib_ref[o, long_axis], lat[o] / (lat_max if lat_max > 0 else 1.0))
            n_broken = len(first_bonds - set(F.bonds_ij))
            artists.append(bond_lc)
        if cursor is not None:
            cursor.set_xdata([times[k], times[k]])
        extra = f"   разорвано связей = {n_broken}" if has_fib else ""
        title.set_text(f"{name}   t = {format_time(times[k], t_max)}   шаг = {int(ids[k])}   "
                       f"смещения ×{K:g}{extra}")
        return artists

    if args.snapshot:
        base, _ = os.path.splitext(out)
        for si in (int(s) for s in args.snapshot.split(",") if s.strip()):
            if 0 <= si < len(ids):
                update(si)
                fig.savefig(f"{base}_frame{si}.png", dpi=args.dpi)
                print(f"снимок: {base}_frame{si}.png")
            else:
                warnings.warn(f"--snapshot {si} вне [0, {len(ids) - 1}]")

    anim = FuncAnimation(fig, update, frames=len(ids), interval=1000.0 / args.fps, blit=False)
    print(f"{name}: {len(ids)} кадров -> {out}")
    anim.save(out, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--run-dir", help="каталог одного прогона (frames/<прогон>)")
    src.add_argument("--all", help="каталог --frames-dir теста: анимировать все прогоны в нём")
    ap.add_argument("--out", help="файл .gif (для --run-dir)")
    ap.add_argument("--out-dir", help="каталог для .gif (для --all)")
    ap.add_argument("--color-by", choices=SPH_FIELDS, default="sxx",
                    help="поле цвета частиц Solid (по умолчанию sxx -- видна волна напряжений)")
    ap.add_argument("--exaggerate", default="auto", help="множитель смещений от первого кадра или auto")
    ap.add_argument("--fps", type=float, default=12.0)
    ap.add_argument("--elev", type=float, default=22.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--size", type=float, default=14.0, help="размер точек Solid")
    ap.add_argument("--max-frames", type=int, default=80)
    ap.add_argument("--dpi", type=int, default=90)
    ap.add_argument("--snapshot", type=str, default=None, help="номера кадров (через запятую) для PNG рядом с .gif")
    args = ap.parse_args()

    if args.run_dir:
        out = args.out or os.path.join(args.run_dir, "anim.gif")
        render_run(args.run_dir, out, args)
    else:
        out_dir = args.out_dir or args.all
        os.makedirs(out_dir, exist_ok=True)
        for d in sorted(glob.glob(os.path.join(args.all, "*"))):
            if os.path.isdir(d) and discover_frame_ids(d):
                render_run(d, os.path.join(out_dir, os.path.basename(d) + ".gif"), args)


if __name__ == "__main__":
    main()
