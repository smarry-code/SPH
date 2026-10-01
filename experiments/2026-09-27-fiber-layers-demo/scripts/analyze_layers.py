#!/usr/bin/env python3
"""Numeric check of pierce/stop criteria for the 4-layer fiber net demo.

Reads sph_output_*.txt / fiber_output_*.txt frame pairs from --run-dir
(produced by demo_fiber_layers) and prints, per frame or as a summary:
  - projectile mean/min z and mean vz
  - number of surviving bonds, broken bonds per layer (layer = frame-0
    object_id // fibers_per_layer, since bead id never changes even after
    a break reassigns object_id -- see demo_fiber_layers.cpp header comment)
  - NaN/Inf check across all parsed arrays

Usage:
    python analyze_layers.py --run-dir data/pierce --fibers-per-layer 5
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np

FRAME_RE = re.compile(r"sph_output_(\d+)\.txt$")


def discover_frame_ids(run_dir: str) -> list[str]:
    ids = []
    for f in glob.glob(os.path.join(run_dir, "sph_output_*.txt")):
        m = FRAME_RE.search(os.path.basename(f))
        if m:
            ids.append(m.group(1))
    ids.sort(key=int)
    return ids


def parse_sph_frame(path: str):
    with open(path) as f:
        cols = f.readline().split()
        rows = [line.split() for line in f if line.strip()]
    col_idx = {name: i for i, name in enumerate(cols)}
    data = np.array(rows, dtype=float)
    x = data[:, col_idx["x"]]
    y = data[:, col_idx["y"]]
    z = data[:, col_idx["z"]]
    vx = data[:, col_idx["vx"]]
    vy = data[:, col_idx["vy"]]
    vz = data[:, col_idx["vz"]]
    return dict(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)


def parse_fiber_frame(path: str):
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    header = lines[0].split()
    h = {header[i]: header[i + 1] for i in range(0, len(header) - 1, 2)}
    time = float(h["TIME"])
    step = int(float(h["STEP"]))
    n_particles = int(h["N_PARTICLES"])
    n_bonds = int(h["N_BONDS"])

    idx = 1
    while lines[idx].strip().startswith("#"):
        idx += 1
    if lines[idx].strip() == "PARTICLES":
        idx += 1
    bead_obj = np.zeros(n_particles, dtype=int)
    bead_pos = np.zeros((n_particles, 3))
    for k in range(n_particles):
        parts = lines[idx + k].split()
        bid = int(parts[0])
        bead_obj[bid] = int(parts[1])
        bead_pos[bid] = [float(parts[3]), float(parts[4]), float(parts[5])]
    idx += n_particles

    while lines[idx].strip().startswith("#"):
        idx += 1
    if lines[idx].strip() == "BONDS":
        idx += 1
    bonds_ij = []
    for k in range(n_bonds):
        parts = lines[idx + k].split()
        bonds_ij.append((int(parts[0]), int(parts[1])))

    return dict(time=time, step=step, bead_obj=bead_obj, bead_pos=bead_pos, bonds_ij=bonds_ij, n_bonds=n_bonds)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--fibers-per-layer", type=int, default=5)
    ap.add_argument("--n-layers", type=int, default=4)
    ap.add_argument("--every", type=int, default=1, help="print every Nth frame (default: all)")
    args = ap.parse_args()

    frame_ids = discover_frame_ids(args.run_dir)
    if not frame_ids:
        print(f"No frames in {args.run_dir}", file=sys.stderr)
        sys.exit(1)

    sph0 = None
    fiber0 = None
    any_nan = False
    print(f"{'frame':>5} {'time_us':>10} {'z_mean':>10} {'z_min':>10} {'vz_mean':>10} "
          f"{'n_bonds':>8} {'broken_L0':>9} {'broken_L1':>9} {'broken_L2':>9} {'broken_L3':>9}")
    for k, fid in enumerate(frame_ids):
        sph_path = os.path.join(args.run_dir, f"sph_output_{fid}.txt")
        fib_path = os.path.join(args.run_dir, f"fiber_output_{fid}.txt")
        if not os.path.exists(fib_path):
            continue
        sf = parse_sph_frame(sph_path)
        ff = parse_fiber_frame(fib_path)
        if sph0 is None:
            sph0, fiber0 = sf, ff

        for key in ("x", "y", "z", "vx", "vy", "vz"):
            if not np.all(np.isfinite(sf[key])):
                any_nan = True
        if not np.all(np.isfinite(ff["bead_pos"])):
            any_nan = True

        if k % args.every != 0 and k != len(frame_ids) - 1:
            continue

        layer0 = fiber0["bead_obj"] // args.fibers_per_layer  # id -> layer, stable across breaks
        first_bonds = set(fiber0["bonds_ij"])
        cur_bonds = set(ff["bonds_ij"])
        broken = first_bonds - cur_bonds
        broken_per_layer = [0] * args.n_layers
        for (i, j) in broken:
            broken_per_layer[layer0[i]] += 1

        z_mean = sf["z"].mean()
        z_min = sf["z"].min()
        vz_mean = sf["vz"].mean()
        print(f"{k:5d} {ff['time']*1e6:10.2f} {z_mean*1e3:10.4f} {z_min*1e3:10.4f} {vz_mean:10.4f} "
              f"{ff['n_bonds']:8d} " + " ".join(f"{c:9d}" for c in broken_per_layer))

    print()
    print("NaN/Inf detected:", any_nan)
    print(f"Total bonds at frame 0: {len(fiber0['bonds_ij'])}  (bonds/layer nominal = "
          f"{len(fiber0['bonds_ij']) // args.n_layers})")


if __name__ == "__main__":
    main()
