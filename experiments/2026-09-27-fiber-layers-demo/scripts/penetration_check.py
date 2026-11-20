#!/usr/bin/env python3
"""Contact-penetration checker for fiber demo runs (SPH solid + discrete
fiber beads/bonds), companion to energy_balance.py in the same directory
(same --run-dir frame discovery, same sph_output/fiber_output parsing
approach -- discover_frame_ids/parse_fiber_frame are imported from it
directly; see the note on parse_sph_frame_geom() below for why the
solid-particle parser is a separate, self-contained function instead).

Reads paired sph_output_NNNNN.txt / fiber_output_NNNNN.txt frames from a run
directory and reports, per frame:

  1. Fiber-Fiber-Contact penetration. Primitive = pair of BOND SEGMENTS
     belonging to different CURRENT fibers: object_id is re-read from
     EVERY frame, never frozen at frame 0 -- when a bond breaks, the tail of
     the fiber gets a NEW object_id (include/interactions/postprocessing.hpp,
     around "topology.next[i] = -1; topology.prev[j] = -1; ... object_id[...
     ] = new_id"), so two fragments of what was originally one fiber count
     as "different fibers" for contact purposes exactly as the running code
     does. Closest points between the two segments follow Ericson,
     "Real-Time Collision Detection" S5.1.9 (ClosestPtSegmentSegment),
     vectorised over ALL cross-fiber bond pairs at once with numpy (a
     python-level loop only over frames, not over pairs -- at most a few
     hundred bonds per frame here, brute-force O(N^2) in one shot is cheap
     and keeps CPU use modest). Contact distance matches
     include/interactions/fiber_contact.hpp::compute_fiber_contact_geometry
     for a single fiber material: d_eff = r_a+r_b = (d_a+d_b)/2 = d
     (--fiber-d), contact iff D < d_eff, penetration delta = d_eff - D.

     CAVEAT vs. the real kernel: this is a standalone reimplementation for
     post-processing, not a bit-for-bit port. It skips (a) the spatial-index
     search radius (moot here -- we just brute-force every cross-fiber bond
     pair), (b) snapping s/t to the 0/1 boundary within eta=1e-9, and (c) the
     "ownership" rule that accepts exactly one primitive per physical
     intersection (notes/fiber-contact/fiber-contact-segment-pairs.md SS3). Near a node
     shared by two bonds of the SAME fiber, both of that fiber's segments
     can therefore report a contact against the same segment of the other
     fiber, where the real kernel accepts only one. This inflates the raw
     "pairs in contact" COUNT in that situation but does not change
     max(delta)/d (the true worst overlap is still found -- generally both
     neighbours report almost the same delta there).

  2. Pass-through / tunnelling detection. For a fiber-fiber pair that is in
     contact OR within --near-mult * d (default 2*d, per the task spec) in
     frame k, AND whose two bonds (identified by their stable bead-id pairs,
     not by array position -- bond arrays shrink as bonds break) both still
     exist (unbroken) in frame k+1: freeze the frame-k contact normal
     n_k = (qA-qB)/|qA-qB| and the frame-k material parameters s,t, then
     re-evaluate the SAME material points at frame k+1's bead positions
     (linear interpolation with the frozen s,t along the SAME two bonds).
     If (qA'-qB').n_k < 0, point A ended up on the other side of B along the
     frame-k separation direction -- i.e. the two fibres tunnelled through
     each other between the two sampled frames without the contact force
     ever separating them in the (unsampled) intervening steps. This is the
     geometric signature behind fiber/bc2-contact-type2-mechanism.md.

  3. Fiber-Solid-Contact penetration. Same bond segments vs. every real SPH
     solid-particle position (point-to-segment, Ericson S5.1.2,
     ClosestPtPointSegment, also vectorised). Contact radius replicates
     include/interactions/fiber_solid_contact.hpp /
     include/interactions/solid_contact.hpp EXACTLY:
       r_f = d/2                                  (--fiber-d)
       r_s = solid_contact_radius(sys,k)
           = 0.5 * (mass_k / max(rho_k, rho_0))^(1/3)   (dim=3; fibers only
             exist at dim==3 in this codebase, see fiber_contact.hpp's
             "if constexpr (dim != 3)" guard, so hardcoding the cube root
             here is safe)
     using the frame's own 'M' and 'rho' sph_output columns for mass_k/rho_k,
     and --solid-rho0 for the reference density rho_0 (default 2700 kg/m^3
     = aluminium, matching sys.rho_0[k]=rho0=2700 set in
     build_projectile_cube() in both demo_fiber_net.cpp and
     demo_fiber_layers.cpp -- confirmed by reading those files).
     d_eff_fs = r_f + r_s, delta = d_eff_fs - D, contact iff delta > 0.
     CAVEAT: omits point-contacts from fully-isolated fiber beads with NO
     remaining bond at all (both neighbours broken) -- a rare edge case not
     covered by the "bond segment" primitive used here (the real kernel
     handles it as a degenerate point segment, fiber_solid_contact.hpp
     decision 2).

  4. Optional --layer-of-object K: layer = object_id_at_FRAME_0 // K, looked
     up per bead id (bead ids are stable across frames -- exactly
     energy_balance.py's convention, see its compute_bead_masses/
     bead_layer0). Fiber-fiber pairs are broken down by the sorted pair of
     layers of their two segments' anchor (first) bead.

Validated on (see final report for numbers):
  experiments/2026-09-26-fiber-net-demo/data/run       (d=1e-3, old contact code, complete: t=4.96e-4 of 5e-4)
  experiments/2026-09-27-fiber-layers-demo/data/pierce  (d=1.5e-3, --layer-of-object 5, complete: t=6.00e-4 of 6e-4)

Usage:
  python penetration_check.py --run-dir ../../2026-09-26-fiber-net-demo/data/run \\
      --fiber-d 1e-3 --solid-rho0 2700 \\
      --csv /tmp/run_penetration.csv --json /tmp/run_penetration.json

  python penetration_check.py --run-dir ../data/pierce \\
      --fiber-d 1.5e-3 --solid-rho0 2700 --layer-of-object 5 \\
      --csv /tmp/pierce_penetration.csv --json /tmp/pierce_penetration.json
"""
from __future__ import annotations

import argparse
import csv as csv_module
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from energy_balance import discover_frame_ids, parse_fiber_frame  # noqa: E402


# ---------------------------------------------------------------------
# Solid-frame parsing (self-contained: needs 'rho' in addition to the
# columns energy_balance.parse_sph_frame extracts, so re-reading via that
# function and re-opening the file for 'rho' would double the I/O -- this
# mirrors its header-driven column-index approach instead, see module
# docstring).
# ---------------------------------------------------------------------
def parse_sph_frame_geom(path: str) -> dict:
    """Returns dict with x,y,z (each (N,)), mass (from 'M' column), rho."""
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
    rho = col("rho")
    if mass is None or rho is None:
        raise ValueError(f"{path}: sph_output header {cols!r} missing 'M' and/or 'rho' "
                          f"(needed to reconstruct solid_contact_radius per interactions/solid_contact.hpp)")
    return dict(x=col("x"), y=col("y"), z=col("z"), mass=mass, rho=rho)


def solid_contact_radius(mass: np.ndarray, rho: np.ndarray, rho0: float) -> np.ndarray:
    """0.5*(mass/max(rho,rho0))^(1/3) -- interactions/solid_contact.hpp::
    solid_contact_radius, dim=3 hardcoded (fibers only exist at dim==3)."""
    rho_eff = np.maximum(rho, rho0)
    vol = np.divide(mass, rho_eff, out=np.zeros_like(mass), where=rho_eff > 0.0)
    r = 0.5 * np.cbrt(vol)
    return np.where(rho > 0.0, r, 0.0)


# ---------------------------------------------------------------------
# Vectorised geometry primitives (Ericson, Real-Time Collision Detection)
# ---------------------------------------------------------------------
def _dot_rows(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    return np.sum(u * v, axis=-1)


def closest_seg_seg(A0: np.ndarray, A1: np.ndarray, B0: np.ndarray, B1: np.ndarray, eps: float = 1e-20):
    """Closest points between every segment in A=(A0,A1) and every segment
    in B=(B0,B1) (SS5.1.9, ClosestPtSegmentSegment). A0,A1: (M,3); B0,B1: (N,3).
    Returns s (M,N), t (M,N) in [0,1], qA (M,N,3), qB (M,N,3), D (M,N).
    See module docstring for what this omits vs. the real kernel (snapping,
    ownership rule)."""
    d1 = A1 - A0                          # (M,3)
    d2 = B1 - B0                          # (N,3)
    r = A0[:, None, :] - B0[None, :, :]   # (M,N,3)

    a = _dot_rows(d1, d1)                 # (M,)
    e = _dot_rows(d2, d2)                 # (N,)
    b = d1 @ d2.T                         # (M,N)
    c = np.einsum('ik,ijk->ij', d1, r)    # (M,N)
    f = np.einsum('jk,ijk->ij', d2, r)    # (M,N)

    a_col = a[:, None]                    # (M,1)
    e_row = e[None, :]                    # (1,N)
    a_safe = np.where(a_col > eps, a_col, 1.0)
    e_safe = np.where(e_row > eps, e_row, 1.0)

    denom = a_col * e_row - b * b
    denom_safe = np.where(denom > 0.0, denom, 1.0)
    s_gen = np.clip((b * f - c * e_row) / denom_safe, 0.0, 1.0)
    s_gen = np.where(denom > 0.0, s_gen, 0.0)  # parallel segments -> s=0 (Ericson)

    t_gen = (b * s_gen + f) / e_safe

    s_tlt0 = np.clip(-c / a_safe, 0.0, 1.0)
    s_tgt1 = np.clip((b - c) / a_safe, 0.0, 1.0)

    t = np.where(t_gen < 0.0, 0.0, np.where(t_gen > 1.0, 1.0, t_gen))
    s = np.where(t_gen < 0.0, s_tlt0, np.where(t_gen > 1.0, s_tgt1, s_gen))

    a_deg = a_col <= eps   # (M,1), broadcasts against (M,N)
    e_deg = e_row <= eps   # (1,N)
    t_adeg = np.clip(f / e_safe, 0.0, 1.0)      # segment A degenerate -> project onto B
    s_edeg = np.clip(-c / a_safe, 0.0, 1.0)     # segment B degenerate -> project onto A

    both_deg = a_deg & e_deg
    s = np.where(both_deg, 0.0, np.where(a_deg, 0.0, np.where(e_deg, s_edeg, s)))
    t = np.where(both_deg, 0.0, np.where(a_deg, t_adeg, np.where(e_deg, 0.0, t)))

    qA = A0[:, None, :] + s[..., None] * d1[:, None, :]
    qB = B0[None, :, :] + t[..., None] * d2[None, :, :]
    D = np.sqrt(np.maximum(_dot_rows(qA - qB, qA - qB), 0.0))
    return s, t, qA, qB, D


def closest_seg_point(A0: np.ndarray, A1: np.ndarray, P: np.ndarray, eps: float = 1e-20):
    """Closest point on every segment in A=(A0,A1) to every point in P
    (SS5.1.2, ClosestPtPointSegment). A0,A1: (M,3); P: (N,3).
    Returns t (M,N) in [0,1], Q (M,N,3), D (M,N)."""
    d1 = A1 - A0                                    # (M,3)
    a = _dot_rows(d1, d1)                            # (M,)
    a_col = a[:, None]
    a_safe = np.where(a_col > eps, a_col, 1.0)
    diff = P[None, :, :] - A0[:, None, :]            # (M,N,3)
    numer = np.einsum('ik,ijk->ij', d1, diff)        # (M,N)
    t = np.clip(numer / a_safe, 0.0, 1.0)
    t = np.where(a_col <= eps, 0.0, t)               # degenerate segment -> point A0
    Q = A0[:, None, :] + t[..., None] * d1[:, None, :]
    D = np.sqrt(np.maximum(_dot_rows(P[None, :, :] - Q, P[None, :, :] - Q), 0.0))
    return t, Q, D


# ---------------------------------------------------------------------
# Per-frame analysis
# ---------------------------------------------------------------------
def build_bond_segments(ff: dict):
    """Bond segments for one fiber frame: endpoints A0,A1 (M,3), the bead-id
    pairs (bi,bj) (stable across frames -- used as bond identity), and the
    CURRENT per-segment object_id (bead_obj of the first bead -- both ends
    of a bond always share the same fiber/fragment object_id)."""
    bi = ff["bond_i"]
    bj = ff["bond_j"]
    A0 = ff["pos"][bi]
    A1 = ff["pos"][bj]
    obj = ff["bead_obj"][bi]
    return A0, A1, bi, bj, obj


def fiber_fiber_frame(ff: dict, fiber_d: float, near_mult: float, bead_layer0):
    """Returns (summary dict, list of 'near' pair records for pass-through
    tracking). summary contains max_delta_ratio, n_contact, n_near, and (if
    bead_layer0 is not None) a dict layer_pair -> (max_delta_ratio, n_contact)."""
    A0, A1, bi, bj, obj = build_bond_segments(ff)
    M = len(bi)
    summary = dict(max_delta_ratio=-np.inf, n_contact=0, n_near=0, layers={})
    records = []
    if M < 2:
        return summary, records

    s, t, qA, qB, D = closest_seg_seg(A0, A1, A0, A1)
    iu, ju = np.triu_indices(M, k=1)
    Diu = D[iu, ju]
    diff_obj = obj[iu] != obj[ju]
    share_bead = (bi[iu] == bi[ju]) | (bi[iu] == bj[ju]) | (bj[iu] == bi[ju]) | (bj[iu] == bj[ju])
    valid = diff_obj & ~share_bead
    if not np.any(valid):
        return summary, records

    delta = fiber_d - Diu
    delta_v = delta[valid]
    ratio_v = delta_v / fiber_d
    summary["max_delta_ratio"] = float(ratio_v.max())
    summary["n_contact"] = int(np.sum(delta_v > 0.0))

    near = valid & (Diu < near_mult * fiber_d)
    summary["n_near"] = int(np.sum(near))

    if bead_layer0 is not None:
        la = bead_layer0[bi[iu]]
        lb = bead_layer0[bi[ju]]
        lay_lo = np.minimum(la, lb)
        lay_hi = np.maximum(la, lb)
        for lo, hi, d_ratio, is_c in zip(lay_lo[valid], lay_hi[valid], ratio_v, delta_v > 0.0):
            key = f"L{int(lo)}-L{int(hi)}"
            cur = summary["layers"].setdefault(key, [-np.inf, 0])
            if d_ratio > cur[0]:
                cur[0] = float(d_ratio)
            if is_c:
                cur[1] += 1

    # Pass-through candidate records: near pairs, with frozen s,t,n,D and
    # bead-id identity of both bonds.
    near_idx = np.nonzero(near)[0]
    for kk in near_idx:
        i, j = iu[kk], ju[kk]
        Dij = D[i, j]
        if Dij < 1e-15:
            continue  # normal undefined (coincident points) -- skip, can't track direction
        n = (qA[i, j] - qB[i, j]) / Dij
        records.append(dict(
            keyA=(int(bi[i]), int(bj[i])), keyB=(int(bi[j]), int(bj[j])),
            s=float(s[i, j]), t=float(t[i, j]), n=n, D=float(Dij),
            layer_a=(int(bead_layer0[bi[i]]) if bead_layer0 is not None else None),
            layer_b=(int(bead_layer0[bi[j]]) if bead_layer0 is not None else None),
        ))

    return summary, records


def fiber_solid_frame(ff: dict, sf: dict, fiber_d: float, rho0: float):
    A0, A1, bi, bj, obj = build_bond_segments(ff)
    M = len(bi)
    summary = dict(max_delta_ratio=-np.inf, n_contact=0)
    if M == 0 or sf["x"] is None or len(sf["x"]) == 0:
        return summary

    P = np.stack([sf["x"], sf["y"], sf["z"]], axis=1)
    r_s = solid_contact_radius(sf["mass"], sf["rho"], rho0)   # (N,)

    _, _, D = closest_seg_point(A0, A1, P)   # (M,N)
    r_f = 0.5 * fiber_d
    d_eff = r_f + r_s[None, :]
    delta = d_eff - D
    ratio = delta / fiber_d
    summary["max_delta_ratio"] = float(ratio.max())
    summary["n_contact"] = int(np.sum(delta > 0.0))
    return summary


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="directory with sph_output_*.txt / fiber_output_*.txt pairs")
    ap.add_argument("--fiber-d", type=float, required=True, help="fiber diameter, m (uniform material, matches Meyer d_eff=d)")
    ap.add_argument("--solid-rho0", type=float, default=2700.0,
                     help="solid reference/initial density rho_0, kg/m^3, used in "
                          "solid_contact_radius=0.5*(mass/max(rho,rho0))^(1/3) (default 2700 = aluminium, "
                          "matches build_projectile_cube() in both demo drivers)")
    ap.add_argument("--near-mult", type=float, default=2.0,
                     help="pass-through tracking considers a fiber-fiber pair if D < near_mult*d (default 2.0, per spec)")
    ap.add_argument("--layer-of-object", type=int, default=None, metavar="K",
                     help="if set, layer = object_id_at_frame0 // K; breaks fiber-fiber contact down by layer pair "
                          "(K=fibers_per_layer=5 for the layers demo)")
    ap.add_argument("--every", type=int, default=5, help="print every Nth row to stdout (default 5)")
    ap.add_argument("--max-events", type=int, default=20, help="max pass-through events kept/printed (default 20)")
    ap.add_argument("--csv", default=None, help="write per-frame CSV here")
    ap.add_argument("--json", default=None, help="write summary JSON here")
    args = ap.parse_args()

    frame_ids = discover_frame_ids(args.run_dir)
    if not frame_ids:
        print(f"No paired sph_output/fiber_output frames found in {args.run_dir}", file=sys.stderr)
        sys.exit(1)

    bead_layer0 = None
    n_layers = 0
    if args.layer_of_object:
        fiber0 = parse_fiber_frame(os.path.join(args.run_dir, f"fiber_output_{frame_ids[0]}.txt"))
        bead_layer0 = fiber0["bead_obj"] // args.layer_of_object
        n_layers = int(bead_layer0.max()) + 1 if len(bead_layer0) else 0

    print(f"run_dir={args.run_dir}  frames={len(frame_ids)}  fiber_d={args.fiber_d:.4e}  "
          f"solid_rho0={args.solid_rho0:.1f}  near_mult={args.near_mult}")
    if bead_layer0 is not None:
        print(f"layer_of_object=K={args.layer_of_object}  n_layers={n_layers}")

    rows = []
    all_layer_keys = set()
    passthrough_events = []
    prev_records = None
    prev_time = None
    prev_frame_idx = None

    for k, fid in enumerate(frame_ids):
        sph_path = os.path.join(args.run_dir, f"sph_output_{fid}.txt")
        fib_path = os.path.join(args.run_dir, f"fiber_output_{fid}.txt")
        sf = parse_sph_frame_geom(sph_path)
        ff = parse_fiber_frame(fib_path)

        ff_summary, ff_records = fiber_fiber_frame(ff, args.fiber_d, args.near_mult, bead_layer0)
        fs_summary = fiber_solid_frame(ff, sf, args.fiber_d, args.solid_rho0)
        all_layer_keys.update(ff_summary["layers"].keys())

        bonds_set = set(zip(ff["bond_i"].tolist(), ff["bond_j"].tolist()))

        if prev_records is not None:
            for rec in prev_records:
                if rec["keyA"] in bonds_set and rec["keyB"] in bonds_set:
                    a0i, a1i = rec["keyA"]
                    b0i, b1i = rec["keyB"]
                    qA2 = (1.0 - rec["s"]) * ff["pos"][a0i] + rec["s"] * ff["pos"][a1i]
                    qB2 = (1.0 - rec["t"]) * ff["pos"][b0i] + rec["t"] * ff["pos"][b1i]
                    signed = float(np.dot(qA2 - qB2, rec["n"]))
                    if signed < 0.0:
                        passthrough_events.append(dict(
                            frame_k=prev_frame_idx, frame_k1=k,
                            time_k=prev_time, time_k1=ff["time"],
                            bond_a=rec["keyA"], bond_b=rec["keyB"],
                            layer_a=rec["layer_a"], layer_b=rec["layer_b"],
                            D_k=rec["D"], signed_separation_k1=signed,
                        ))

        row = dict(
            frame=k, file_id=fid, time=ff["time"],
            n_bonds=ff["n_bonds"],
            ff_max_delta_ratio=ff_summary["max_delta_ratio"],
            ff_n_contact=ff_summary["n_contact"],
            ff_n_near=ff_summary["n_near"],
            fs_max_delta_ratio=fs_summary["max_delta_ratio"],
            fs_n_contact=fs_summary["n_contact"],
            ff_layers=ff_summary["layers"],
        )
        rows.append(row)

        prev_records = ff_records
        prev_time = ff["time"]
        prev_frame_idx = k

    # --- stdout table ---
    header = (f"{'frame':>5} {'time':>12} {'n_bonds':>8} {'ff_max_d/d':>11} {'ff_ncnt':>8} "
              f"{'ff_nnear':>9} {'fs_max_d/d':>11} {'fs_ncnt':>8}")
    print()
    print(header)
    imax_ff = max(range(len(rows)), key=lambda i: rows[i]["ff_max_delta_ratio"])
    imax_fs = max(range(len(rows)), key=lambda i: rows[i]["fs_max_delta_ratio"])
    printed = set()
    for k, r in enumerate(rows):
        if k % args.every == 0 or k == len(rows) - 1 or k in (imax_ff, imax_fs):
            if k in printed:
                continue
            printed.add(k)
            print(f"{k:5d} {r['time']:12.5e} {r['n_bonds']:8d} {r['ff_max_delta_ratio']:11.4f} "
                  f"{r['ff_n_contact']:8d} {r['ff_n_near']:9d} {r['fs_max_delta_ratio']:11.4f} "
                  f"{r['fs_n_contact']:8d}")

    max_ff_ratio = rows[imax_ff]["ff_max_delta_ratio"]
    max_fs_ratio = rows[imax_fs]["fs_max_delta_ratio"]

    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"max fiber-fiber delta/d   = {max_ff_ratio:.6f}  (frame {imax_ff}, t={rows[imax_ff]['time']:.6e})")
    print(f"max fiber-solid  delta/d  = {max_fs_ratio:.6f}  (frame {imax_fs}, t={rows[imax_fs]['time']:.6e})")
    print(f"pass-through events       = {len(passthrough_events)}")
    for ev in passthrough_events[:5]:
        print(f"  frame {ev['frame_k']}->{ev['frame_k1']} (t={ev['time_k']:.6e}->{ev['time_k1']:.6e}): "
              f"bond {ev['bond_a']} vs {ev['bond_b']}  layers "
              f"{ev['layer_a']}/{ev['layer_b']}  D_k={ev['D_k']:.4e}  signed(k+1)={ev['signed_separation_k1']:.4e}")
    if len(passthrough_events) > 5:
        print(f"  ... and {len(passthrough_events) - 5} more")

    if bead_layer0 is not None and all_layer_keys:
        print()
        print(f"Per-layer-pair breakdown (layer = object_id@frame0 // {args.layer_of_object}, {n_layers} layers), "
              f"max over all frames:")
        agg = {}
        for r in rows:
            for key, (ratio, ncontact) in r["ff_layers"].items():
                cur = agg.setdefault(key, [-np.inf, 0])
                if ratio > cur[0]:
                    cur[0] = ratio
                cur[1] = max(cur[1], ncontact)
        for key in sorted(agg.keys()):
            ratio, ncontact = agg[key]
            print(f"  {key}: max delta/d = {ratio:.4f}  max n_contact (single frame) = {ncontact}")

    # --- CSV ---
    if args.csv:
        layer_keys_sorted = sorted(all_layer_keys)
        fieldnames = ["frame", "file_id", "time", "n_bonds",
                      "ff_max_delta_ratio", "ff_n_contact", "ff_n_near",
                      "fs_max_delta_ratio", "fs_n_contact"]
        for key in layer_keys_sorted:
            fieldnames += [f"ff_delta_ratio_{key}", f"ff_ncontact_{key}"]
        with open(args.csv, "w", newline="") as f:
            w = csv_module.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in rows:
                out = {k2: r[k2] for k2 in fieldnames if k2 in r}
                out["frame"] = r["frame"]; out["file_id"] = r["file_id"]; out["time"] = r["time"]
                out["n_bonds"] = r["n_bonds"]
                out["ff_max_delta_ratio"] = r["ff_max_delta_ratio"]
                out["ff_n_contact"] = r["ff_n_contact"]
                out["ff_n_near"] = r["ff_n_near"]
                out["fs_max_delta_ratio"] = r["fs_max_delta_ratio"]
                out["fs_n_contact"] = r["fs_n_contact"]
                for key in layer_keys_sorted:
                    ratio, ncontact = r["ff_layers"].get(key, (float("nan"), 0))
                    out[f"ff_delta_ratio_{key}"] = ratio
                    out[f"ff_ncontact_{key}"] = ncontact
                w.writerow(out)
        print(f"\nCSV written: {args.csv}")

    # --- JSON ---
    if args.json:
        summary = dict(
            run_dir=args.run_dir,
            n_frames=len(rows),
            fiber_d=args.fiber_d,
            solid_rho0=args.solid_rho0,
            near_mult=args.near_mult,
            max_ff_delta_ratio=max_ff_ratio,
            max_ff_delta_ratio_frame=imax_ff,
            max_ff_delta_ratio_time=rows[imax_ff]["time"],
            max_fs_delta_ratio=max_fs_ratio,
            max_fs_delta_ratio_frame=imax_fs,
            max_fs_delta_ratio_time=rows[imax_fs]["time"],
            n_passthrough_events=len(passthrough_events),
            passthrough_events=passthrough_events[:args.max_events],
        )
        if bead_layer0 is not None and all_layer_keys:
            summary["layer_pairs"] = {
                key: dict(max_delta_ratio=agg[key][0], max_n_contact_single_frame=agg[key][1])
                for key in sorted(agg.keys())
            }
        with open(args.json, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"JSON written: {args.json}")


if __name__ == "__main__":
    main()
