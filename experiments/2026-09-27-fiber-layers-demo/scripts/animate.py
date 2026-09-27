#!/usr/bin/env python3
"""Render a 3D animation of "projectile through fiber net" frames produced by
the SPH solver (sph_output_<NNNNN>.txt + fiber_output_<NNNNN>.txt pairs).

See fiber/DEMO_VISUALIZATION_TASK.md sections 5-6 for the exact file formats
this reads, and experiments/2026-09-26-fiber-net-demo/scripts/make_synthetic_frames.py
for a synthetic generator that produces frames in the same format for testing.

Usage:
    python animate.py --run-dir data/synthetic --out anim.gif
    python animate.py --run-dir data/run --out anim.gif --color-by p \
        --snapshot 0,34,59
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
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)
from mpl_toolkits.mplot3d.art3d import Line3DCollection


# ---------------------------------------------------------------------------
# Frame discovery / parsing
# ---------------------------------------------------------------------------

FRAME_RE = re.compile(r"sph_output_(\d{5})\.txt$")


def discover_frame_ids(run_dir: str) -> list[str]:
    sph_files = glob.glob(os.path.join(run_dir, "sph_output_*.txt"))
    ids = []
    for f in sph_files:
        m = FRAME_RE.search(os.path.basename(f))
        if m:
            ids.append(m.group(1))
    ids.sort()
    return ids


class SphFrame:
    __slots__ = ("pos", "p", "speed")

    def __init__(self, pos: np.ndarray, p: np.ndarray, speed: np.ndarray):
        self.pos = pos
        self.p = p
        self.speed = speed


def parse_sph_frame(path: str) -> SphFrame | None:
    with open(path, "r") as f:
        header_line = f.readline()
        cols = header_line.split()
        try:
            rows = [line.split() for line in f if line.strip()]
        except Exception:
            return None
    if not rows:
        return None
    col_idx = {name: i for i, name in enumerate(cols)}
    required = ["x", "y", "z", "vx", "vy", "vz", "p"]
    if not all(c in col_idx for c in required):
        warnings.warn(f"{path}: header missing expected columns, got {cols}")
        return None

    data = np.array(rows, dtype=float)
    x = data[:, col_idx["x"]]
    y = data[:, col_idx["y"]]
    z = data[:, col_idx["z"]]
    vx = data[:, col_idx["vx"]]
    vy = data[:, col_idx["vy"]]
    vz = data[:, col_idx["vz"]]
    p = data[:, col_idx["p"]]
    pos = np.stack([x, y, z], axis=1)
    speed = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    return SphFrame(pos, p, speed)


class FiberFrame:
    __slots__ = ("time", "step", "bead_pos", "bead_obj", "bonds_ij", "bonds_l0",
                 "bonds_l", "bonds_lbreak", "bonds_elong", "n_bonds")

    def __init__(self, time, step, bead_pos, bead_obj, bonds_ij, bonds_l0, bonds_l,
                 bonds_lbreak, bonds_elong):
        self.time = time
        self.step = step
        self.bead_pos = bead_pos
        self.bead_obj = bead_obj
        self.bonds_ij = bonds_ij
        self.bonds_l0 = bonds_l0
        self.bonds_l = bonds_l
        self.bonds_lbreak = bonds_lbreak
        self.bonds_elong = bonds_elong
        self.n_bonds = len(bonds_ij)


def parse_fiber_frame(path: str) -> FiberFrame | None:
    with open(path, "r") as f:
        lines = [line.rstrip("\n") for line in f]
    if not lines:
        return None

    header = lines[0].split()
    # TIME <t> STEP <step> N_PARTICLES <n> N_BONDS <nb>
    try:
        h = {header[i]: header[i + 1] for i in range(0, len(header) - 1, 2)}
        time = float(h["TIME"])
        step = int(float(h["STEP"]))
        n_particles = int(h["N_PARTICLES"])
        n_bonds = int(h["N_BONDS"])
    except (KeyError, IndexError, ValueError):
        warnings.warn(f"{path}: malformed header line: {lines[0]!r}")
        return None

    idx = 1
    # skip "# id object_id ..." comment then "PARTICLES"
    while idx < len(lines) and lines[idx].strip().startswith("#"):
        idx += 1
    if idx < len(lines) and lines[idx].strip() == "PARTICLES":
        idx += 1

    bead_pos = np.zeros((n_particles, 3))
    bead_obj = np.zeros(n_particles, dtype=int)
    for k in range(n_particles):
        parts = lines[idx + k].split()
        bid = int(parts[0])
        bead_obj[bid] = int(parts[1])
        bead_pos[bid, 0] = float(parts[3])
        bead_pos[bid, 1] = float(parts[4])
        bead_pos[bid, 2] = float(parts[5])
    idx += n_particles

    while idx < len(lines) and lines[idx].strip().startswith("#"):
        idx += 1
    if idx < len(lines) and lines[idx].strip() == "BONDS":
        idx += 1

    bonds_ij = []
    bonds_l0 = np.zeros(n_bonds)
    bonds_l = np.zeros(n_bonds)
    bonds_lbreak = np.zeros(n_bonds)
    bonds_elong = np.zeros(n_bonds)
    for k in range(n_bonds):
        parts = lines[idx + k].split()
        i, j = int(parts[0]), int(parts[1])
        bonds_ij.append((i, j))
        bonds_l0[k] = float(parts[2])
        bonds_l[k] = float(parts[3])
        bonds_lbreak[k] = float(parts[4])
        bonds_elong[k] = float(parts[5])

    return FiberFrame(time, step, bead_pos, bead_obj, bonds_ij, bonds_l0, bonds_l,
                       bonds_lbreak, bonds_elong)


# ---------------------------------------------------------------------------
# Time formatting
# ---------------------------------------------------------------------------

def format_time(t: float, t_max: float) -> str:
    if t_max < 1e-3:
        return f"{t * 1e6:.1f} us"
    if t_max < 1.0:
        return f"{t * 1e3:.3f} ms"
    return f"{t:.4g} s"


# ---------------------------------------------------------------------------
# Main rendering
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="Directory with sph_output_*.txt / fiber_output_*.txt")
    ap.add_argument("--out", required=True, help="Output path (.gif always works; .mp4 needs ffmpeg)")
    ap.add_argument("--color-by", choices=["p", "speed"], default="speed",
                    help="SPH scatter color source (default: speed)")
    ap.add_argument("--fps", type=float, default=12.0)
    ap.add_argument("--elev", type=float, default=20.0)
    ap.add_argument("--azim", type=float, default=-60.0)
    ap.add_argument("--max-frames", type=int, default=None,
                    help="Subsample to at most this many frames")
    ap.add_argument("--dpi", type=int, default=110)
    ap.add_argument("--snapshot", type=str, default=None,
                    help="Comma-separated frame indices (into the selected sequence, 0-based) "
                         "to also save as standalone PNGs next to --out")
    ap.add_argument("--fibers-per-layer", type=int, default=5,
                    help="Beads colored by LAYER (frame-0 object_id // this), not by object_id -- "
                         "with 20 fibers the tab10 colormap keyed on object_id would be unreadable "
                         "and it drifts after breaks reassign object_id. Bead's frame-0 layer is "
                         "stable across the whole run (bead id itself never changes).")
    ap.add_argument("--n-layers", type=int, default=4)
    args = ap.parse_args()

    frame_ids = discover_frame_ids(args.run_dir)
    if not frame_ids:
        print(f"No sph_output_*.txt frames found in {args.run_dir}", file=sys.stderr)
        sys.exit(1)

    pairs = []  # (frame_id, sph_path, fiber_path)
    for fid in frame_ids:
        sph_path = os.path.join(args.run_dir, f"sph_output_{fid}.txt")
        fiber_path = os.path.join(args.run_dir, f"fiber_output_{fid}.txt")
        if not os.path.exists(fiber_path):
            warnings.warn(f"frame {fid}: sph_output present but fiber_output missing, skipping")
            continue
        pairs.append((fid, sph_path, fiber_path))

    if not pairs:
        print("No complete sph+fiber frame pairs found.", file=sys.stderr)
        sys.exit(1)

    if args.max_frames and len(pairs) > args.max_frames:
        idx = np.linspace(0, len(pairs) - 1, args.max_frames).round().astype(int)
        idx = sorted(set(idx.tolist()))
        pairs = [pairs[i] for i in idx]

    print(f"Loading {len(pairs)} frame pairs from {args.run_dir} ...")
    sph_frames: list[SphFrame] = []
    fiber_frames: list[FiberFrame] = []
    kept_pairs = []
    for fid, sph_path, fiber_path in pairs:
        sf = parse_sph_frame(sph_path)
        ff = parse_fiber_frame(fiber_path)
        if sf is None or ff is None:
            warnings.warn(f"frame {fid}: failed to parse, skipping")
            continue
        sph_frames.append(sf)
        fiber_frames.append(ff)
        kept_pairs.append(fid)

    if not sph_frames:
        print("No frames parsed successfully.", file=sys.stderr)
        sys.exit(1)

    n_frames = len(sph_frames)
    print(f"{n_frames} frames parsed OK.")

    # --- global bounding box (SPH + fiber beads), equal aspect, small margin ---
    # ROBUST percentiles (0.5/99.5), not raw min/max: a handful of broken-bond
    # fragments get flung out at speeds far above the projectile's own (see
    # README.md "странности физики" -- a freed single bead has ~1e-6 kg mass,
    # so even a brief residual contact force gives it a huge acceleration).
    # Frame-0 (nothing has happened yet) or any early frame would otherwise
    # be squeezed into a tiny blob by a bounding box sized for that one late
    # outlier bead -- raw min/max was tried first and produced exactly that.
    all_pts = [f.pos for f in sph_frames] + [f.bead_pos for f in fiber_frames]
    all_pts = np.concatenate(all_pts, axis=0)
    lo = np.percentile(all_pts, 0.5, axis=0)
    hi = np.percentile(all_pts, 99.5, axis=0)
    center = (lo + hi) / 2.0
    half_range = max((hi - lo).max() / 2.0, 1e-6) * 1.15
    xlim = (center[0] - half_range, center[0] + half_range)
    ylim = (center[1] - half_range, center[1] + half_range)
    zlim = (center[2] - half_range, center[2] + half_range)

    # --- global color scales ---
    t_max = max(f.time for f in fiber_frames)

    if args.color_by == "p":
        all_p = np.concatenate([f.p for f in sph_frames])
        lo, hi = np.percentile(all_p, [1, 99])
        if lo == hi:
            lo, hi = lo - 1.0, hi + 1.0
        sph_vmin, sph_vmax = lo, hi
        sph_label = "Pressure p, Pa"
        sph_cmap = "coolwarm" if lo < 0 < hi else "viridis"
    else:
        all_speed = np.concatenate([f.speed for f in sph_frames])
        lo, hi = np.percentile(all_speed, [1, 99])
        if lo == hi:
            lo, hi = 0.0, hi + 1.0
        sph_vmin, sph_vmax = max(0.0, lo), hi
        sph_label = "Speed |v|, m/s"
        sph_cmap = "viridis"

    elong_cmap = "RdYlGn_r"  # green (0) -> red (1, at/over break)

    # bond-break bookkeeping relative to first frame
    first_bonds = set(fiber_frames[0].bonds_ij)

    # Layer index per bead, from FRAME 0's object_id (= FIBER-block index at
    # load time, see demo_fiber_layers.cpp: fibers are written layer-by-layer,
    # so object_id_at_t0 // fibers_per_layer == layer). Bead id (array index)
    # never changes across a break -- only object_id does (broken-off pieces
    # get a fresh, larger object_id) -- so this per-bead layer mapping stays
    # correct and STABLE for the whole run, unlike coloring by the live
    # object_id (which would (a) run out of tab10's 10 colors with 20 fibers,
    # and (b) reassign colors mid-animation whenever a bond snaps).
    bead_layer0 = (fiber_frames[0].bead_obj // args.fibers_per_layer).astype(float)
    layer_vmax = max(args.n_layers - 1, 1)

    def ranked_obj(ff: "FiberFrame") -> np.ndarray:  # noqa: ARG001 (kept for call-site symmetry)
        return bead_layer0

    # --- figure setup (created once; reused each frame) ---
    # Explicit axes rect (not add_subplot(111)) -- the original layout left a
    # large blank margin above a small plot in the default 10x8 figure (3D
    # axes shrink to keep box_aspect square, and the figure wasn't resized to
    # compensate); this fills more of the canvas so beads/net are legible in
    # the still snapshots this scene relies on (20 fibers vs. the original's 2).
    fig = plt.figure(figsize=(11, 8.5))
    ax = fig.add_axes((0.02, 0.02, 0.62, 0.88), projection="3d")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_xlabel("x, m")
    ax.set_ylabel("y, m")
    ax.set_zlabel("z, m")

    sf0 = sph_frames[0]
    c0 = sf0.p if args.color_by == "p" else sf0.speed
    sph_scatter = ax.scatter(sf0.pos[:, 0], sf0.pos[:, 1], sf0.pos[:, 2],
                              c=c0, cmap=sph_cmap, vmin=sph_vmin, vmax=sph_vmax,
                              s=14, depthshade=True, label="projectile (SPH)")

    ff0 = fiber_frames[0]
    bead_scatter = ax.scatter(ff0.bead_pos[:, 0], ff0.bead_pos[:, 1], ff0.bead_pos[:, 2],
                               c=ranked_obj(ff0), cmap="tab10", vmin=0, vmax=layer_vmax,
                               s=8, depthshade=False,
                               label="fiber beads (color = layer, 0=top)")

    def bond_segments_and_colors(ff: FiberFrame):
        if ff.n_bonds == 0:
            return np.zeros((0, 2, 3)), np.zeros(0)
        segs = np.stack([ff.bead_pos[[i for i, j in ff.bonds_ij]],
                          ff.bead_pos[[j for i, j in ff.bonds_ij]]], axis=1)
        eps_break = ff.bonds_lbreak / ff.bonds_l0 - 1.0
        eps_break = np.where(eps_break > 1e-12, eps_break, 1e-12)
        ratio = ff.bonds_elong / eps_break
        ratio = np.clip(ratio, 0.0, 1.0)  # compression (elong<0) also clamps to 0
        return segs, ratio

    segs0, ratio0 = bond_segments_and_colors(ff0)
    bond_lc = Line3DCollection(segs0, linewidths=2.0, cmap=elong_cmap, clim=(0.0, 1.0))
    bond_lc.set_array(ratio0)
    ax.add_collection3d(bond_lc)

    # fig.text (not ax.set_title): mplot3d places the axes title relative to
    # the actual rendered 3D box, which floats around inside the nominal axes
    # rect depending on view_init/box_aspect -- with the enlarged axes rect
    # above that left a large, unpredictable gap. A fixed figure-relative
    # text sits exactly where placed, every frame.
    title = fig.text(0.02, 0.965, "", fontsize=11, va="top", ha="left")

    cbar_sph = fig.colorbar(sph_scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar_sph.set_label(sph_label)
    cbar_bond = fig.colorbar(bond_lc, ax=ax, shrink=0.6, pad=0.10)
    cbar_bond.set_label(r"$\varepsilon / \varepsilon_{break}$ (bond elongation)")

    def update(frame_idx: int):
        sf = sph_frames[frame_idx]
        ff = fiber_frames[frame_idx]

        sph_scatter._offsets3d = (sf.pos[:, 0], sf.pos[:, 1], sf.pos[:, 2])
        sph_scatter.set_array(sf.p if args.color_by == "p" else sf.speed)

        bead_scatter._offsets3d = (ff.bead_pos[:, 0], ff.bead_pos[:, 1], ff.bead_pos[:, 2])
        bead_scatter.set_array(ranked_obj(ff))

        segs, ratio = bond_segments_and_colors(ff)
        bond_lc.set_segments(segs)
        bond_lc.set_array(ratio)

        broken = first_bonds - set(ff.bonds_ij)
        broken_per_layer = [0] * args.n_layers
        for (i, _j) in broken:
            broken_per_layer[int(bead_layer0[i])] += 1
        broken_str = "/".join(str(c) for c in broken_per_layer)
        t_str = format_time(ff.time, t_max)
        z_min = sf.pos[:, 2].min()
        title.set_text(
            f"t = {t_str}   step = {ff.step}   N_BONDS = {ff.n_bonds}   "
            f"broken/layer(0..{args.n_layers - 1}) = {broken_str}   z_min(proj) = {z_min * 1e3:.2f} mm"
        )
        return sph_scatter, bead_scatter, bond_lc, title

    # --- snapshots (rendered directly, without going through the animation) ---
    if args.snapshot:
        snap_indices = [int(s) for s in args.snapshot.split(",") if s.strip() != ""]
        out_base, _ = os.path.splitext(args.out)
        for si in snap_indices:
            if si < 0 or si >= n_frames:
                warnings.warn(f"--snapshot index {si} out of range [0,{n_frames - 1}], skipping")
                continue
            update(si)
            fig.canvas.draw()
            snap_path = f"{out_base}_frame{si}.png"
            fig.savefig(snap_path, dpi=args.dpi)
            print(f"Saved snapshot: {snap_path}")

    # --- full animation ---
    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000.0 / args.fps, blit=False)

    ext = os.path.splitext(args.out)[1].lower()
    use_mp4 = False
    if ext == ".mp4":
        import shutil
        if shutil.which("ffmpeg") is not None:
            use_mp4 = True
        else:
            print("ffmpeg not found on PATH; falling back to GIF next to --out", file=sys.stderr)
            args.out = os.path.splitext(args.out)[0] + ".gif"

    print(f"Rendering animation ({n_frames} frames) to {args.out} ...")
    if use_mp4:
        from matplotlib.animation import FFMpegWriter
        try:
            anim.save(args.out, writer=FFMpegWriter(fps=args.fps), dpi=args.dpi)
        except Exception as e:
            print(f"ffmpeg encoding failed ({e}); falling back to GIF next to --out", file=sys.stderr)
            args.out = os.path.splitext(args.out)[0] + ".gif"
            anim.save(args.out, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    else:
        anim.save(args.out, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    print("Done.")


if __name__ == "__main__":
    main()
