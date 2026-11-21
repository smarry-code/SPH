#!/usr/bin/env python3
"""Generate synthetic sph_output_*.txt / fiber_output_*.txt frame pairs in the
EXACT text format produced by the real solver, for debugging animate.py before
real simulation data is available.

Formats reproduced (see src/output.cpp, include/particle_system.hpp::vuvod_index,
src/fiber_output.cpp):

sph_output_<NNNNN>.txt (dim == 3):
    line 1 (header, space separated):
        N Nmat M x y z vx vy vz rho p u c sigma_i e_pl D hsml sxx syy szz sxy syz sxz
    then one line per particle:
        <1-based id> <itype> <mass> <x> <y> <z> <vx> <vy> <vz>
        <rho> <p> <u> <c> <sigma_i> <e_pl> <D> <hsml>
        <sxx> <syy> <szz> <sxy> <syz> <sxz>
    (23 columns total, matches vuvod_index: " itype mass x y z vx vy vz rho p u c
    sigma_i e_pl D hsml sxx syy szz sxy syz sxz").

fiber_output_<NNNNN>.txt:
    TIME <t> STEP <step> N_PARTICLES <n> N_BONDS <nb>
    # id object_id material_id x y z vx vy vz
    PARTICLES
    <n lines: id(0-based) object_id material_id x y z vx vy vz>
    # i j l0 l l_break elongation
    BONDS
    <nb lines: i j l0 l l_break elongation>   (i, j 0-based bead ids)

Scene (simple prescribed kinematics, NOT a physical simulation):
  - A cube-shaped SPH projectile (~500 particles) falls along -z towards a
    cross-shaped net made of two straight fiber chains (~20 beads each) lying
    in the z ~ 0 plane, crossing at the origin.
  - The beads sag towards -z as the projectile passes, with the sag amplitude
    peaking when the projectile crosses z = 0 (a Gaussian bump in time and in
    planar distance from the impact axis).
  - One bond near the crossing point is scripted to exceed l_break at a fixed
    frame; from the next frame onward it is dropped from BONDS, and the beads
    on one side of the break are reassigned a new object_id (and given extra
    downward velocity), so the net visibly splits into pieces.
"""
from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# Scene parameters (SI units; only relative magnitudes matter for the demo)
# ---------------------------------------------------------------------------

N_FRAMES = 60
DT = 2.0e-5          # s, time between saved frames
T0 = 0.0             # s, time of frame 0

# Projectile: cube grid of SPH particles
PROJ_NX = 8
PROJ_NY = 8
PROJ_NZ = 8
PROJ_DX = 0.0025          # m, lattice spacing -> cube side 0.02 m
PROJ_RHO = 2700.0         # kg/m3, aluminum-like
PROJ_ITYPE = 4            # body/solid material id, matches build_plate() convention
PROJ_HSML_FACTOR = 1.5
PROJ_Z0 = 0.035           # m, initial center height above net (net at z=0)
PROJ_VZ = -60.0           # m/s, downward velocity (tens of m/s, per task doc)
PROJ_C0 = 5200.0          # m/s, nominal sound speed (aluminum-like), for display only

# Net: two crossing chains ("+"), each with BEADS_PER_CHAIN beads
BEADS_PER_CHAIN = 21          # odd -> one bead sits exactly at the crossing arm center
CHAIN_HALF_LEN = 0.05         # m, half length of each chain arm
FIBER_MATERIAL_ID = 1
BEAD_MASS = 3.0e-6            # kg, nominal bead mass
ANCHOR_MASS_FACTOR = 50.0     # pseudo-anchor: end beads made much heavier
ANCHOR_BEADS = 2              # number of beads at each end that are heavy

SAG_AMPLITUDE = 0.03           # m, peak sag at the crossing point
SAG_SIGMA_R = 0.02             # m, planar Gaussian width of the sag bump
SAG_SIGMA_Z = 0.018            # m, Gaussian width in projectile height (time-like)

EPS_BREAK = 0.12               # elongation threshold used for l_break = l0*(1+EPS_BREAK)
BREAK_FRAME = 34               # frame index at which the critical bond's l exceeds l_break
BREAK_CHAIN = "A"              # which chain loses a bond
BREAK_BOND_LOCAL_I = 9         # local bond index (connects bead 9 and bead 10) in that chain

FINITE_DIFF_H = 1.0e-7         # s, small step for numerical velocity of fiber beads


@dataclass
class Chain:
    name: str
    object_id: int
    axis: int              # 0 -> beads spread along x, 1 -> along y
    beads_xy0: np.ndarray = field(default_factory=lambda: np.zeros(0))  # (n,2) rest x,y (z=0)
    l0: float = 0.0


def build_chains() -> list[Chain]:
    n = BEADS_PER_CHAIN
    s = np.linspace(-CHAIN_HALF_LEN, CHAIN_HALF_LEN, n)
    chain_a = Chain(name="A", object_id=0, axis=0)
    chain_a.beads_xy0 = np.stack([s, np.zeros(n)], axis=1)
    chain_a.l0 = s[1] - s[0]

    chain_b = Chain(name="B", object_id=1, axis=1)
    chain_b.beads_xy0 = np.stack([np.zeros(n), s], axis=1)
    chain_b.l0 = s[1] - s[0]
    return [chain_a, chain_b]


def sag_amplitude_at(proj_z: float) -> float:
    """Gaussian-in-height bump: sag peaks when the projectile is near z=0."""
    return SAG_AMPLITUDE * math.exp(-(proj_z / SAG_SIGMA_Z) ** 2)


def bead_positions(chain: Chain, t: float, broken: bool, split_local_idx: int | None,
                    frame_idx: int) -> np.ndarray:
    """Return (n,3) bead positions at time t for a chain, given prescribed kinematics."""
    proj_z = PROJ_Z0 + PROJ_VZ * t
    amp = sag_amplitude_at(proj_z)
    xy = chain.beads_xy0
    r = np.sqrt(xy[:, 0] ** 2 + xy[:, 1] ** 2)
    shape = np.exp(-(r / SAG_SIGMA_R) ** 2)
    z = -amp * shape

    if broken and split_local_idx is not None:
        # Beads beyond the break point on the "downstream" side keep falling
        # faster/further once detached, to make the split visually obvious.
        frames_since_break = max(0, frame_idx - BREAK_FRAME)
        extra = 0.004 * frames_since_break
        idx = np.arange(len(xy))
        mask = idx > split_local_idx
        z = z.copy()
        z[mask] -= extra

    pos = np.stack([xy[:, 0], xy[:, 1], z], axis=1)
    return pos


def compute_object_ids(chain: Chain, frame_idx: int) -> np.ndarray:
    n = BEADS_PER_CHAIN
    obj = np.full(n, chain.object_id, dtype=int)
    if chain.name == BREAK_CHAIN and frame_idx > BREAK_FRAME:
        # Beads downstream of the broken bond become a new, separate object.
        new_obj_id = 100 + chain.object_id
        idx = np.arange(n)
        obj[idx > BREAK_BOND_LOCAL_I] = new_obj_id
    return obj


def chain_bonds(chain: Chain, frame_idx: int) -> list[tuple[int, int]]:
    """Local (i, i+1) bond list for a chain at a given frame (0-based local ids)."""
    n = BEADS_PER_CHAIN
    bonds = [(i, i + 1) for i in range(n - 1)]
    if chain.name == BREAK_CHAIN and frame_idx > BREAK_FRAME:
        bonds = [b for b in bonds if b != (BREAK_BOND_LOCAL_I, BREAK_BOND_LOCAL_I + 1)]
    return bonds


def build_projectile_positions(t: float) -> np.ndarray:
    nx, ny, nz = PROJ_NX, PROJ_NY, PROJ_NZ
    idx = np.array([(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)])
    cx = (nx - 1) / 2.0
    cy = (ny - 1) / 2.0
    cz = (nz - 1) / 2.0
    x = (idx[:, 0] - cx) * PROJ_DX
    y = (idx[:, 1] - cy) * PROJ_DX
    z0_local = (idx[:, 2] - cz) * PROJ_DX
    z_center = PROJ_Z0 + PROJ_VZ * t
    z = z0_local + z_center
    return np.stack([x, y, z], axis=1)


def write_sph_frame(path: str, t: float) -> np.ndarray:
    pos = build_projectile_positions(t)
    n = pos.shape[0]
    vx = np.zeros(n)
    vy = np.zeros(n)
    vz = np.full(n, PROJ_VZ)

    proj_z_center = PROJ_Z0 + PROJ_VZ * t
    # Synthetic pressure: builds up as the projectile approaches/crosses the
    # net plane, with a bit of radial variation so the color map is not flat.
    r_local = np.sqrt((pos[:, 0]) ** 2 + (pos[:, 1]) ** 2)
    impact_factor = math.exp(-(proj_z_center / SAG_SIGMA_Z) ** 2)
    p = 5.0e7 * impact_factor * (1.0 - 0.5 * (r_local / (PROJ_NX * PROJ_DX)))
    p = np.clip(p, 0.0, None)

    rho = np.full(n, PROJ_RHO)
    u = p / (rho * 500.0)
    c = np.full(n, PROJ_C0)
    sigma_i = p * 0.3
    e_pl = np.clip(impact_factor - 0.3, 0.0, None) * 0.02
    e_pl = np.full(n, e_pl)
    D = np.zeros(n)
    hsml = np.full(n, PROJ_DX * PROJ_HSML_FACTOR)
    sxx = -p * 0.3
    syy = -p * 0.3
    szz = -p * 0.6
    sxy = np.zeros(n)
    syz = np.zeros(n)
    sxz = np.zeros(n)

    header = ("N Nmat M x y z vx vy vz rho p u c sigma_i e_pl D hsml "
              "sxx syy szz sxy syz sxz\n")
    with open(path, "w") as f:
        f.write(header)
        for i in range(n):
            f.write(
                f"{i + 1} {PROJ_ITYPE} {BEAD_MASS_PROJ:.8g} "
                f"{pos[i, 0]:.8g} {pos[i, 1]:.8g} {pos[i, 2]:.8g} "
                f"{vx[i]:.8g} {vy[i]:.8g} {vz[i]:.8g} "
                f"{rho[i]:.8g} {p[i]:.8g} {u[i]:.8g} {c[i]:.8g} {sigma_i[i]:.8g} "
                f"{e_pl[i]:.8g} {D[i]:.8g} {hsml[i]:.8g} "
                f"{sxx[i]:.8g} {syy[i]:.8g} {szz[i]:.8g} "
                f"{sxy[i]:.8g} {syz[i]:.8g} {sxz[i]:.8g}\n"
            )
    return pos


BEAD_MASS_PROJ = PROJ_RHO * PROJ_DX ** 3


def write_fiber_frame(path: str, t: float, step: int, frame_idx: int,
                       chains: list[Chain]) -> None:
    all_pos = []
    all_vel = []
    all_obj = []
    all_mat = []
    offsets = {}
    next_id = 0
    for chain in chains:
        offsets[chain.name] = next_id
        broken = chain.name == BREAK_CHAIN
        split_idx = BREAK_BOND_LOCAL_I if broken else None
        pos = bead_positions(chain, t, broken and frame_idx > BREAK_FRAME, split_idx, frame_idx)
        pos_h = bead_positions(chain, t + FINITE_DIFF_H, broken and frame_idx > BREAK_FRAME,
                                split_idx, frame_idx)
        vel = (pos_h - pos) / FINITE_DIFF_H
        obj = compute_object_ids(chain, frame_idx)

        all_pos.append(pos)
        all_vel.append(vel)
        all_obj.append(obj)
        all_mat.append(np.full(len(pos), FIBER_MATERIAL_ID, dtype=int))
        next_id += len(pos)

    pos = np.concatenate(all_pos, axis=0)
    vel = np.concatenate(all_vel, axis=0)
    obj = np.concatenate(all_obj, axis=0)
    mat = np.concatenate(all_mat, axis=0)
    n = pos.shape[0]

    bonds = []  # (i, j, l0, l_break)
    for chain in chains:
        base = offsets[chain.name]
        local_bonds = chain_bonds(chain, frame_idx)
        l_break = chain.l0 * (1.0 + EPS_BREAK)
        for (li, lj) in local_bonds:
            bonds.append((base + li, base + lj, chain.l0, l_break))
    nb = len(bonds)

    with open(path, "w") as f:
        f.write(f"TIME {t:.8g} STEP {step} N_PARTICLES {n} N_BONDS {nb}\n")
        f.write("# id object_id material_id x y z vx vy vz\n")
        f.write("PARTICLES\n")
        for i in range(n):
            f.write(
                f"{i} {obj[i]} {mat[i]} "
                f"{pos[i, 0]:.8g} {pos[i, 1]:.8g} {pos[i, 2]:.8g} "
                f"{vel[i, 0]:.8g} {vel[i, 1]:.8g} {vel[i, 2]:.8g}\n"
            )
        f.write("# i j l0 l l_break elongation\n")
        f.write("BONDS\n")
        for (i, j, l0, l_break) in bonds:
            l = float(np.linalg.norm(pos[j] - pos[i]))
            elongation = (l - l0) / l0 if l0 > 0 else 0.0
            f.write(f"{i} {j} {l0:.8g} {l:.8g} {l_break:.8g} {elongation:.8g}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir",
                    default="experiments/2026-09-26-fiber-net-demo/data/synthetic",
                    help="Directory to write sph_output_*.txt / fiber_output_*.txt into")
    ap.add_argument("--n-frames", type=int, default=N_FRAMES)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    chains = build_chains()

    for frame_idx in range(args.n_frames):
        t = T0 + frame_idx * DT
        step = frame_idx * 10  # arbitrary but monotonic "solver step" number

        sph_path = os.path.join(args.out_dir, f"sph_output_{frame_idx:05d}.txt")
        fiber_path = os.path.join(args.out_dir, f"fiber_output_{frame_idx:05d}.txt")

        write_sph_frame(sph_path, t)
        write_fiber_frame(fiber_path, t, step, frame_idx, chains)

    print(f"Wrote {args.n_frames} synthetic frame pairs to {args.out_dir}")


if __name__ == "__main__":
    main()
