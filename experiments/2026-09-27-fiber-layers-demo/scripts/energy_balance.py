#!/usr/bin/env python3
"""Energy-balance checker for the fiber demo runs (SPH solid + discrete fiber
beads/bonds).

Reads paired sph_output_NNNNN.txt / fiber_output_NNNNN.txt frames from a run
directory (produced by demo_fiber_net.cpp or demo_fiber_layers.cpp) and
tracks, per frame:

    KE_solid   = sum 0.5*m_solid*|v|^2   (solid SPH particles)
    KE_fiber   = sum 0.5*m_bead*|v|^2    (fiber beads)
    dU_int     = sum m_solid*u  - (same at frame 0)   (SPH internal energy,
                 elastic+plastic+artificial-viscosity heating, taken relative
                 to frame 0 since only its *change* is a like-for-like energy
                 term -- absolute u already includes whatever reference the
                 EOS uses)
    U_bond     = sum over CURRENT bonds of 0.5*E*A/l0*max(l-l0,0)^2
                 (tension-only bond law, matches
                 include/interactions/fiber_bond.hpp: F = E*A*(l/l0-1) for
                 l>l0, else 0 -- so U = 0.5*(E*A/l0)*(l-l0)^2 is the spring
                 potential whose derivative reproduces that force)
    E_tot      = KE_solid + KE_fiber + dU_int + U_bond

Hertz contact potential energy between solid and fiber (or solid-solid) is
NOT tracked anywhere in the output files, so transient dips in E_tot during
contact are expected and not a bug. What must NOT happen is E_tot drifting
UP beyond frame 0 over the whole run -- that signals energy being created
(see fiber/fiber-fiber-contact-energy-creation.md in the project memory).
Friction heat and the energy stored in a bond at the instant it breaks are
real, physical losses (they leave U_bond/KE without appearing anywhere else
in this tally) -- they only ever show up as *drops* in E_tot, never rises.

Bead masses are NOT in the fiber_output files. They are reconstructed from
the FIRST frame's bond topology (before any breaks), matching how the C++
demo drivers do it (fiber/fiber_bond.hpp for the axial stiffness law,
demo_fiber_net.cpp / demo_fiber_layers.cpp for the anchor-mass hack):

    A = pi*d^2/4
    each bond (i,j,l0) contributes rho_f*A*l0/2 to EACH of beads i and j
    (interior bead of a chain ends up with rho_f*A*l0 from its two bonds,
    a chain-end bead only gets the half from its single bond)
    then beads with degree 1 in the frame-0 bond graph ("anchors" -- the
    two ends of every straight fiber in these demos) have their mass
    multiplied by --anchor-factor (the demos use 50 to pin the fiber ends
    in place without adding explicit boundary conditions)

Bead ids are stable across frames (fiber_output.cpp writes "id" as the
particle's index inside FiberState, which never renumbers on a bond break --
only topology.bonds shrinks), so this frame-0 mass/layer assignment stays
valid for the whole run.

Usage examples:

  # experiments/2026-09-26-fiber-net-demo/data/run (see demo_fiber_net.cpp):
  python energy_balance.py --run-dir ../../2026-09-26-fiber-net-demo/data/run \\
      --solid-mass-dx 1.2e-3 --solid-rho 2700 \\
      --fiber-E 2e7 --fiber-d 1e-3 --fiber-rho 1200 --anchor-factor 50 \\
      --csv /tmp/run_energy.csv --plot /tmp/run_energy.png --json /tmp/run_energy.json

  # experiments/2026-09-27-fiber-layers-demo/data/pierce (see
  # demo_fiber_layers.cpp):
  python energy_balance.py --run-dir ../data/pierce \\
      --solid-mass-dx 1.3e-3 --solid-rho 2700 \\
      --fiber-E 5e7 --fiber-d 1.5e-3 --fiber-rho 1200 --anchor-factor 50 \\
      --layer-of-object 5 \\
      --csv /tmp/pierce_energy.csv --plot /tmp/pierce_energy.png --json /tmp/pierce_energy.json

--solid-mass (a plain per-particle mass in kg) can be used instead of
--solid-mass-dx/--solid-rho. If the sph_output header has an "M" column
(observed in both demos: it already equals rho_s*dx^3 exactly), that
per-particle value is read directly from the file and takes precedence over
both -- --solid-mass/--solid-mass-dx are then only a fallback/cross-check.
"""
from __future__ import annotations

import argparse
import csv as csv_module
import glob
import json
import math
import os
import re
import sys

import numpy as np

FRAME_RE = re.compile(r"sph_output_(\d{5})\.txt$")


# ---------------------------------------------------------------------
# Frame discovery / parsing
# ---------------------------------------------------------------------
def discover_frame_ids(run_dir: str) -> list[str]:
    ids = []
    for f in glob.glob(os.path.join(run_dir, "sph_output_*.txt")):
        m = FRAME_RE.search(os.path.basename(f))
        if m:
            fid = m.group(1)
            if os.path.exists(os.path.join(run_dir, f"fiber_output_{fid}.txt")):
                ids.append(fid)
    ids.sort()
    return ids


def parse_sph_frame(path: str):
    """Returns dict with x,y,z,vx,vy,vz,u arrays and mass (per-particle
    array, from the 'M' header column if present, else None)."""
    with open(path) as f:
        cols = f.readline().split()
        rows = [line.split() for line in f if line.strip()]
    col_idx = {name: i for i, name in enumerate(cols)}
    if not rows:
        raise ValueError(f"{path}: no particle rows")
    data = np.array(rows, dtype=float)

    def col(name):
        if name not in col_idx:
            return None
        return data[:, col_idx[name]]

    mass = col("M")
    u = col("u")
    if u is None:
        raise ValueError(f"{path}: no 'u' column in header {cols!r}")
    return dict(
        x=col("x"), y=col("y"), z=col("z"),
        vx=col("vx"), vy=col("vy"), vz=col("vz"),
        u=u, mass=mass,
    )


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
    assert lines[idx].strip() == "PARTICLES", f"{path}: expected PARTICLES at line {idx+1}"
    idx += 1
    bead_obj = np.zeros(n_particles, dtype=int)
    pos = np.zeros((n_particles, 3))
    vel = np.zeros((n_particles, 3))
    for k in range(n_particles):
        parts = lines[idx + k].split()
        bid = int(parts[0])
        bead_obj[bid] = int(parts[1])
        pos[bid] = [float(parts[3]), float(parts[4]), float(parts[5])]
        vel[bid] = [float(parts[6]), float(parts[7]), float(parts[8])]
    idx += n_particles

    while lines[idx].strip().startswith("#"):
        idx += 1
    assert lines[idx].strip() == "BONDS", f"{path}: expected BONDS at line {idx+1}"
    idx += 1
    bond_i = np.zeros(n_bonds, dtype=int)
    bond_j = np.zeros(n_bonds, dtype=int)
    bond_l0 = np.zeros(n_bonds)
    bond_l = np.zeros(n_bonds)
    for k in range(n_bonds):
        parts = lines[idx + k].split()
        bond_i[k] = int(parts[0])
        bond_j[k] = int(parts[1])
        bond_l0[k] = float(parts[2])
        bond_l[k] = float(parts[3])
        # parts[4] = l_break, parts[5] = elongation -- both derivable, not needed here

    return dict(
        time=time, step=step, n_particles=n_particles, n_bonds=n_bonds,
        bead_obj=bead_obj, pos=pos, vel=vel,
        bond_i=bond_i, bond_j=bond_j, bond_l0=bond_l0, bond_l=bond_l,
    )


# ---------------------------------------------------------------------
# Bead mass reconstruction from frame-0 topology
# ---------------------------------------------------------------------
def compute_bead_masses(fiber0: dict, fiber_area: float, fiber_rho: float, anchor_factor: float):
    n = fiber0["n_particles"]
    mass = np.zeros(n)
    degree = np.zeros(n, dtype=int)
    for i, j, l0 in zip(fiber0["bond_i"], fiber0["bond_j"], fiber0["bond_l0"]):
        m_half = fiber_rho * fiber_area * l0 / 2.0
        mass[i] += m_half
        mass[j] += m_half
        degree[i] += 1
        degree[j] += 1
    is_anchor = degree == 1
    mass[is_anchor] *= anchor_factor
    n_anchors = int(is_anchor.sum())
    n_isolated = int((degree == 0).sum())
    if n_isolated:
        print(f"WARNING: {n_isolated} fiber bead(s) with no bonds in frame 0 "
              f"(mass=0) -- excluded from any per-bead diagnostics that divide by mass",
              file=sys.stderr)
    return mass, is_anchor, n_anchors


def bond_energy(fiber_frame: dict, E: float, A: float) -> float:
    l0 = fiber_frame["bond_l0"]
    l = fiber_frame["bond_l"]
    stretch = np.maximum(l - l0, 0.0)
    # avoid div-by-zero for a degenerate l0==0 bond (shouldn't occur in practice)
    safe_l0 = np.where(l0 > 0.0, l0, 1.0)
    u = 0.5 * E * A / safe_l0 * stretch * stretch
    return float(np.sum(np.where(l0 > 0.0, u, 0.0)))


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="directory with sph_output_*.txt / fiber_output_*.txt pairs")

    ap.add_argument("--solid-mass", type=float, default=None,
                     help="per-particle solid mass in kg (fallback if sph_output has no 'M' column)")
    ap.add_argument("--solid-mass-dx", type=float, default=None,
                     help="solid particle spacing dx (m); combined with --solid-rho as rho*dx^3 (fallback)")
    ap.add_argument("--solid-rho", type=float, default=None, help="solid material density kg/m^3 (fallback)")

    ap.add_argument("--fiber-E", type=float, required=True, help="fiber Young's modulus, Pa")
    ap.add_argument("--fiber-d", type=float, required=True, help="fiber diameter, m")
    ap.add_argument("--fiber-rho", type=float, required=True, help="fiber material density, kg/m^3")
    ap.add_argument("--anchor-factor", type=float, default=50.0,
                     help="mass multiplier for frame-0 degree-1 (chain-end) beads (default: 50, matches the demos)")

    ap.add_argument("--layer-of-object", type=int, default=None, metavar="K",
                     help="if set, layer = object_id_in_frame0 // K; reports per-layer KE_fiber and bonds remaining "
                          "(K=fibers_per_layer=5 for the layers demo)")

    ap.add_argument("--tol", type=float, default=0.02, help="FAIL if max E_tot/E_tot(0) > 1+tol (default 0.02)")
    ap.add_argument("--every", type=int, default=5, help="print every Nth row to stdout (default 5)")

    ap.add_argument("--csv", default=None, help="write per-frame CSV here")
    ap.add_argument("--plot", default=None, help="write PNG (E_tot + components vs time) here")
    ap.add_argument("--json", default=None, help="write summary JSON here")
    args = ap.parse_args()

    if args.solid_mass is None and (args.solid_mass_dx is not None and args.solid_rho is not None):
        args.solid_mass = args.solid_rho * args.solid_mass_dx ** 3

    frame_ids = discover_frame_ids(args.run_dir)
    if not frame_ids:
        print(f"No paired sph_output/fiber_output frames found in {args.run_dir}", file=sys.stderr)
        sys.exit(1)

    fiber_area = math.pi * args.fiber_d ** 2 / 4.0

    # --- frame 0: bead masses from topology, layer assignment ---
    fiber0_path = os.path.join(args.run_dir, f"fiber_output_{frame_ids[0]}.txt")
    fiber0 = parse_fiber_frame(fiber0_path)
    bead_mass, is_anchor, n_anchors = compute_bead_masses(fiber0, fiber_area, args.fiber_rho, args.anchor_factor)
    bead_layer0 = fiber0["bead_obj"] // args.layer_of_object if args.layer_of_object else None
    n_layers = int(bead_layer0.max()) + 1 if bead_layer0 is not None and len(bead_layer0) else 0

    print(f"run_dir={args.run_dir}  frames={len(frame_ids)}  "
          f"fiber beads={fiber0['n_particles']} (anchors={n_anchors})  "
          f"bonds@0={fiber0['n_bonds']}")
    if args.solid_mass is not None:
        print(f"solid mass (CLI fallback) = {args.solid_mass:.6e} kg/particle")

    rows = []
    solid_mass_source_reported = False
    for k, fid in enumerate(frame_ids):
        sph_path = os.path.join(args.run_dir, f"sph_output_{fid}.txt")
        fib_path = os.path.join(args.run_dir, f"fiber_output_{fid}.txt")
        sf = parse_sph_frame(sph_path)
        ff = parse_fiber_frame(fib_path)

        if sf["mass"] is not None:
            solid_mass = sf["mass"]
            if not solid_mass_source_reported:
                uniq = np.unique(solid_mass)
                print(f"solid mass: read from sph_output 'M' column "
                      f"({'uniform ' + format(uniq[0], '.6e') if len(uniq) == 1 else f'{len(uniq)} distinct values'})")
                if args.solid_mass is not None and len(uniq) == 1 and not math.isclose(uniq[0], args.solid_mass, rel_tol=1e-6):
                    print(f"  NOTE: differs from --solid-mass/--solid-mass-dx*rho fallback "
                          f"({args.solid_mass:.6e}) -- using the file's 'M' column.", file=sys.stderr)
                solid_mass_source_reported = True
        elif args.solid_mass is not None:
            solid_mass = args.solid_mass
            if not solid_mass_source_reported:
                print("solid mass: no 'M' column in sph_output header -- using CLI fallback for all particles")
                solid_mass_source_reported = True
        else:
            print(f"{sph_path}: no 'M' column and no --solid-mass/--solid-mass-dx+--solid-rho given", file=sys.stderr)
            sys.exit(1)

        v2_solid = sf["vx"] ** 2 + sf["vy"] ** 2 + sf["vz"] ** 2
        KE_solid = float(np.sum(0.5 * solid_mass * v2_solid))
        U_int_solid = float(np.sum(solid_mass * sf["u"]))
        max_solid_speed = float(np.sqrt(v2_solid.max())) if len(v2_solid) else 0.0
        P_solid = np.array([
            np.sum(solid_mass * sf["vx"]),
            np.sum(solid_mass * sf["vy"]),
            np.sum(solid_mass * sf["vz"]),
        ])

        v2_fiber = np.sum(ff["vel"] ** 2, axis=1)
        KE_fiber = float(np.sum(0.5 * bead_mass * v2_fiber))
        max_bead_speed = float(np.sqrt(v2_fiber.max())) if len(v2_fiber) else 0.0
        P_fiber = (bead_mass[:, None] * ff["vel"]).sum(axis=0)

        U_bond = bond_energy(ff, args.fiber_E, fiber_area)

        P = P_solid + P_fiber

        row = dict(
            frame=k, file_id=fid, time=ff["time"],
            KE_solid=KE_solid, KE_fiber=KE_fiber, U_int_solid=U_int_solid, U_bond=U_bond,
            n_bonds=ff["n_bonds"], max_solid_speed=max_solid_speed, max_bead_speed=max_bead_speed,
            Px=P[0], Py=P[1], Pz=P[2],
        )

        if bead_layer0 is not None:
            ke_layer = np.zeros(n_layers)
            bonds_layer = np.zeros(n_layers, dtype=int)
            np.add.at(ke_layer, bead_layer0, 0.5 * bead_mass * v2_fiber)
            layer_of_bond = bead_layer0[ff["bond_i"]]  # both endpoints of a bond share a layer in these demos
            for lay in layer_of_bond:
                bonds_layer[lay] += 1
            row["KE_fiber_layers"] = ke_layer.tolist()
            row["bonds_layers"] = bonds_layer.tolist()

        rows.append(row)

    # --- derived quantities relative to frame 0 ---
    dU_int = np.array([r["U_int_solid"] for r in rows]) - rows[0]["U_int_solid"]
    KE_solid_arr = np.array([r["KE_solid"] for r in rows])
    KE_fiber_arr = np.array([r["KE_fiber"] for r in rows])
    U_bond_arr = np.array([r["U_bond"] for r in rows])
    E_tot = KE_solid_arr + KE_fiber_arr + dU_int + U_bond_arr
    KE_tot = KE_solid_arr + KE_fiber_arr

    P0 = np.array([rows[0]["Px"], rows[0]["Py"], rows[0]["Pz"]])
    P0_norm = np.linalg.norm(P0)
    P_drift = []
    for r in rows:
        P = np.array([r["Px"], r["Py"], r["Pz"]])
        d = np.linalg.norm(P - P0)
        P_drift.append(d / P0_norm if P0_norm > 0 else d)
    P_drift = np.array(P_drift)

    E_tot0 = E_tot[0]
    E_ratio = E_tot / E_tot0 if E_tot0 != 0 else E_tot
    KE_tot0 = KE_tot[0]
    KE_ratio = KE_tot / KE_tot0 if KE_tot0 != 0 else KE_tot

    imax = int(np.argmax(E_ratio))
    max_E_ratio = float(E_ratio[imax])
    final_E_ratio = float(E_ratio[-1])
    final_KE_ratio = float(KE_ratio[-1]) if KE_tot0 != 0 else float("nan")
    max_bead_speed_final = rows[-1]["max_bead_speed"]
    max_bead_speed_overall = max(r["max_bead_speed"] for r in rows)
    max_solid_speed_final = rows[-1]["max_solid_speed"]
    max_solid_speed_overall = max(r["max_solid_speed"] for r in rows)
    max_P_drift = float(P_drift.max())

    passed = max_E_ratio <= 1.0 + args.tol

    # --- stdout table (subsampled + frame of max always included) ---
    header = (f"{'frame':>5} {'time':>12} {'KE_solid':>12} {'KE_fiber':>12} {'dU_int':>12} "
              f"{'U_bond':>12} {'E_tot':>12} {'E/E0':>8} {'n_bonds':>8} {'maxbead_v':>10} {'|dP|/P0':>10}")
    print()
    print(header)
    printed = set()
    for k, r in enumerate(rows):
        if k % args.every == 0 or k == len(rows) - 1 or k == imax:
            if k in printed:
                continue
            printed.add(k)
            print(f"{k:5d} {r['time']:12.5e} {r['KE_solid']:12.5e} {r['KE_fiber']:12.5e} {dU_int[k]:12.5e} "
                  f"{r['U_bond']:12.5e} {E_tot[k]:12.5e} {E_ratio[k]:8.4f} {r['n_bonds']:8d} "
                  f"{r['max_bead_speed']:10.4f} {P_drift[k]:10.3e}")

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"E_tot(0)                    = {E_tot0:.6e} J")
    print(f"max E_tot/E_tot(0)          = {max_E_ratio:.6f}  (frame {imax}, t={rows[imax]['time']:.6e})")
    print(f"final E_tot/E_tot(0)        = {final_E_ratio:.6f}")
    print(f"final KE_total/KE_total(0)  = {final_KE_ratio:.6f}")
    print(f"max bead speed              = {max_bead_speed_overall:.4f} m/s (final: {max_bead_speed_final:.4f})")
    print(f"max solid speed             = {max_solid_speed_overall:.4f} m/s (final: {max_solid_speed_final:.4f})")
    print(f"max |P-P0|/|P0|             = {max_P_drift:.6e}")
    print(f"tol                         = {args.tol}")
    print(f"{'PASS' if passed else 'FAIL'}: max E_tot/E_tot(0) {'<=' if passed else '>'} 1+tol "
          f"({max_E_ratio:.6f} vs {1.0 + args.tol:.6f})")

    if bead_layer0 is not None:
        print()
        print(f"Per-layer breakdown (layer = object_id@frame0 // {args.layer_of_object}, {n_layers} layers):")
        ke_layers_final = np.array(rows[-1]["KE_fiber_layers"])
        ke_layers_all = np.array([r["KE_fiber_layers"] for r in rows])
        ke_layers_max = ke_layers_all.max(axis=0)
        bonds_layers_0 = rows[0]["bonds_layers"]
        bonds_layers_final = rows[-1]["bonds_layers"]
        for lay in range(n_layers):
            print(f"  layer {lay}: KE_fiber final={ke_layers_final[lay]:.5e} J  max={ke_layers_max[lay]:.5e} J  "
                  f"bonds {bonds_layers_0[lay]} -> {bonds_layers_final[lay]}")

    # --- CSV ---
    if args.csv:
        fieldnames = ["frame", "file_id", "time", "KE_solid", "KE_fiber", "dU_int_solid", "U_bond", "E_tot",
                      "E_tot_ratio", "KE_total", "KE_total_ratio", "n_bonds", "max_solid_speed", "max_bead_speed",
                      "Px", "Py", "Pz", "P_drift_rel"]
        if bead_layer0 is not None:
            for lay in range(n_layers):
                fieldnames += [f"KE_fiber_layer{lay}", f"bonds_layer{lay}"]
        with open(args.csv, "w", newline="") as f:
            w = csv_module.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for k, r in enumerate(rows):
                out = dict(
                    frame=r["frame"], file_id=r["file_id"], time=r["time"],
                    KE_solid=r["KE_solid"], KE_fiber=r["KE_fiber"], dU_int_solid=dU_int[k], U_bond=r["U_bond"],
                    E_tot=E_tot[k], E_tot_ratio=E_ratio[k], KE_total=KE_tot[k], KE_total_ratio=KE_ratio[k],
                    n_bonds=r["n_bonds"], max_solid_speed=r["max_solid_speed"], max_bead_speed=r["max_bead_speed"],
                    Px=r["Px"], Py=r["Py"], Pz=r["Pz"], P_drift_rel=P_drift[k],
                )
                if bead_layer0 is not None:
                    for lay in range(n_layers):
                        out[f"KE_fiber_layer{lay}"] = r["KE_fiber_layers"][lay]
                        out[f"bonds_layer{lay}"] = r["bonds_layers"][lay]
                w.writerow(out)
        print(f"\nCSV written: {args.csv}")

    # --- PNG plot ---
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        times = np.array([r["time"] for r in rows])
        fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
        ax = axes[0]
        ax.plot(times, KE_solid_arr, label="KE_solid")
        ax.plot(times, KE_fiber_arr, label="KE_fiber")
        ax.plot(times, dU_int, label="dU_int_solid")
        ax.plot(times, U_bond_arr, label="U_bond")
        ax.plot(times, E_tot, label="E_tot", color="black", linewidth=2)
        ax.axhline(E_tot0, color="gray", linestyle="--", linewidth=1, label="E_tot(0)")
        ax.set_ylabel("Energy, J")
        ax.set_title(f"Energy balance: {args.run_dir}")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        ax2 = axes[1]
        ax2.plot(times, E_ratio, label="E_tot / E_tot(0)", color="black")
        ax2.plot(times, KE_ratio, label="KE_total / KE_total(0)", color="tab:blue")
        ax2.axhline(1.0, color="gray", linestyle="--", linewidth=1)
        ax2.axhline(1.0 + args.tol, color="red", linestyle=":", linewidth=1, label=f"1+tol ({args.tol})")
        ax2.set_xlabel("time, s")
        ax2.set_ylabel("ratio")
        ax2.legend(fontsize=8)
        ax2.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(args.plot, dpi=130)
        print(f"Plot written: {args.plot}")

    # --- JSON summary ---
    if args.json:
        summary = dict(
            run_dir=args.run_dir,
            n_frames=len(rows),
            E_tot0=E_tot0,
            max_E_tot_ratio=max_E_ratio,
            max_E_tot_ratio_frame=imax,
            max_E_tot_ratio_time=rows[imax]["time"],
            final_E_tot_ratio=final_E_ratio,
            final_KE_total_ratio=final_KE_ratio,
            max_bead_speed=max_bead_speed_overall,
            final_bead_speed=max_bead_speed_final,
            max_solid_speed=max_solid_speed_overall,
            final_solid_speed=max_solid_speed_final,
            max_momentum_drift_rel=max_P_drift,
            tol=args.tol,
            passed=bool(passed),
            n_bonds_initial=int(rows[0]["n_bonds"]),
            n_bonds_final=int(rows[-1]["n_bonds"]),
            n_fiber_beads=int(fiber0["n_particles"]),
            n_anchors=n_anchors,
        )
        if bead_layer0 is not None:
            summary["layers"] = [
                dict(layer=lay, KE_fiber_final=float(ke_layers_final[lay]), KE_fiber_max=float(ke_layers_max[lay]),
                     bonds_initial=int(bonds_layers_0[lay]), bonds_final=int(bonds_layers_final[lay]))
                for lay in range(n_layers)
            ]
        with open(args.json, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"JSON written: {args.json}")


if __name__ == "__main__":
    main()
