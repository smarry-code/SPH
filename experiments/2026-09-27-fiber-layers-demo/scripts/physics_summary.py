#!/usr/bin/env python3
"""Per-run physics summary for the fiber demos (net / layers): one markdown
column + JSON block per run.

Reuses (does not duplicate) energy_balance.py and penetration_check.py:
  * parsers / bead-mass reconstruction / frame discovery are imported;
  * max E/E0, final E/E0, max |P-P0|/|P0| come from energy_balance.main()
    and max delta/d (fiber-fiber, fiber-solid) + pass-through count from
    penetration_check.main(). Both are run in-process with patched argv and
    captured stdout; their --json output is read back (their pass-through
    logic lives inside main(), so this is the minimal wrapper).
NEW metrics computed here (single pass over frames, numpy only):
  projectile trajectory / arrested / pierced-vs-carried, anchors, angular
  momentum, broken bonds (per fiber-group and per chain position), fragment
  speeds.

Assumptions: every particle of sph_output is the projectile (true for both
demos); beads are fiber beads with masses rebuilt from frame-0 topology
(anchors = frame-0 degree-1 beads, mass x --anchor-factor).

Usage:
  python physics_summary.py \\
     --run ../data/v3_mu0/stop:layers:stop_mu0 --run ../data/v3_mu0/pierce:layers:pierce_mu0 \\
     --out-md table.md --out-json summary.json
  RUN is DIR:KIND:LABEL (KIND = net | layers; 'C:/...' drive colons are fine:
  the LAST two ':' split KIND/LABEL).
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import energy_balance as eb  # noqa: E402
import penetration_check as pc  # noqa: E402

KIND_DEFAULTS = {
    "net": dict(fiber_E=2e7, fiber_d=1e-3, fiber_rho=1200.0, layer_of_object=None),
    "layers": dict(fiber_E=5e7, fiber_d=1.5e-3, fiber_rho=1200.0, layer_of_object=5),
}


def _run_quiet(fn, argv):
    old = sys.argv
    sys.argv = ["x"] + argv
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn()
    finally:
        sys.argv = old


def chain_positions(fiber0):
    """Walk frame-0 chains. Returns dict bond(i,j sorted)->(chain_id, seg_index),
    and list of chain lengths in beads. Chain starts at the lower-id end."""
    adj = {}
    for i, j in zip(fiber0["bond_i"].tolist(), fiber0["bond_j"].tolist()):
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
    seen, pos, lengths = set(), {}, []
    ends = sorted(b for b, nb in adj.items() if len(nb) == 1)
    for s in ends:
        if s in seen:
            continue
        chain, prev, cur = [s], None, s
        while True:
            nxt = [n for n in adj[cur] if n != prev]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            chain.append(cur)
        seen.update(chain)
        cid = len(lengths)
        lengths.append(len(chain))
        for k in range(len(chain) - 1):
            pos[tuple(sorted((chain[k], chain[k + 1])))] = (cid, k)
    return pos, lengths


def analyse_run(run_dir, kind, args):
    P = dict(KIND_DEFAULTS[kind])
    for k in ("fiber_E", "fiber_d", "fiber_rho", "layer_of_object"):
        v = getattr(args, k)
        if v is not None:
            P[k] = v
    K = P["layer_of_object"]
    fids = eb.discover_frame_ids(run_dir)
    if not fids:
        raise SystemExit(f"no frames in {run_dir}")

    res = {"run_dir": run_dir, "kind": kind, "n_frames": len(fids)}

    # ---- reuse energy_balance / penetration_check ----
    with tempfile.TemporaryDirectory() as td:
        ej, pj = os.path.join(td, "e.json"), os.path.join(td, "p.json")
        eargs = ["--run-dir", run_dir, "--fiber-E", str(P["fiber_E"]), "--fiber-d", str(P["fiber_d"]),
                 "--fiber-rho", str(P["fiber_rho"]), "--anchor-factor", str(args.anchor_factor),
                 "--json", ej, "--tol", "1e9"]
        pargs = ["--run-dir", run_dir, "--fiber-d", str(P["fiber_d"]), "--solid-rho0", str(args.solid_rho0),
                 "--json", pj]
        if K:
            eargs += ["--layer-of-object", str(K)]
            pargs += ["--layer-of-object", str(K)]
        _run_quiet(eb.main, eargs)
        _run_quiet(pc.main, pargs)
        E = json.load(open(ej))
        Pn = json.load(open(pj))
    res["energy"] = dict(
        E0_J=E["E_tot0"], max_E_ratio=E["max_E_tot_ratio"], max_E_ratio_time=E["max_E_tot_ratio_time"],
        final_E_ratio=E["final_E_tot_ratio"], max_momentum_drift_rel=E["max_momentum_drift_rel"],
        n_bonds_initial=E["n_bonds_initial"], n_bonds_final=E["n_bonds_final"])
    res["penetration"] = dict(
        max_ff_delta_over_d=Pn["max_ff_delta_ratio"], max_fs_delta_over_d=Pn["max_fs_delta_ratio"],
        n_passthrough=Pn["n_passthrough_events"])

    # ---- own pass ----
    fa = math.pi * P["fiber_d"] ** 2 / 4.0
    f0 = eb.parse_fiber_frame(os.path.join(run_dir, f"fiber_output_{fids[0]}.txt"))
    bmass, is_anchor, _ = eb.compute_bead_masses(f0, fa, P["fiber_rho"], args.anchor_factor)
    grp0 = f0["bead_obj"] // K if K else f0["bead_obj"]
    bpos, chain_len = chain_positions(f0)
    bonds0 = {tuple(sorted(b)): k for k, b in enumerate(zip(f0["bond_i"].tolist(), f0["bond_j"].tolist()))}
    anchors = np.nonzero(is_anchor)[0]
    x0_anch = f0["pos"][anchors].copy()

    nF = len(fids)
    t = np.zeros(nF); zmin = np.zeros(nF); zmean = np.zeros(nF); vz = np.zeros(nF)
    bead_vmax = np.zeros(nF); n_fast = np.zeros(nF, dtype=int)
    ever_fast = np.zeros(len(bmass), dtype=bool)
    anch_disp = np.zeros(nF); anch_v = np.zeros(nF); anch_vz = np.zeros(nF)
    dL0 = np.zeros(nF); dLc = np.zeros(nF)
    foot_gap = np.full(nF, np.nan)   # min over footprint beads of (z_min_cube - z_bead)
    foot_n = np.zeros(nF, dtype=int); foot_above = np.zeros(nF, dtype=int)
    break_frame = {}   # bond -> first frame index where absent

    for k, fid in enumerate(fids):
        sf = eb.parse_sph_frame(os.path.join(run_dir, f"sph_output_{fid}.txt"))
        ff = f0 if k == 0 else eb.parse_fiber_frame(os.path.join(run_dir, f"fiber_output_{fid}.txt"))
        m = sf["mass"]
        X = np.stack([sf["x"], sf["y"], sf["z"]], 1)
        V = np.stack([sf["vx"], sf["vy"], sf["vz"]], 1)
        if k == 0:
            mp = float(m.sum())
            v0vec = (m[:, None] * V).sum(0) / mp
            v0 = float(np.linalg.norm(v0vec))
            sgn = -1.0 if v0vec[2] < 0 else 1.0
            rho0 = pc.parse_sph_frame_geom(os.path.join(run_dir, f"sph_output_{fid}.txt"))["rho"]
            dx = float(np.cbrt(m.mean() / rho0.mean()))
            ext = X.max(0) - X.min(0) + dx
            L_char = args.cube_side if args.cube_side else float(ext[0])
            Mtot0 = mp + bmass.sum()
            com0 = ((m[:, None] * X).sum(0) + (bmass[:, None] * ff["pos"]).sum(0)) / Mtot0
            Lscale = mp * v0 * L_char
            fast_thr = args.fast_factor * v0
        t[k] = ff["time"]
        zmin[k], zmean[k] = X[:, 2].min(), X[:, 2].mean()
        vz[k] = (m * V[:, 2]).sum() / mp
        sp = np.linalg.norm(ff["vel"], axis=1)
        bead_vmax[k] = sp.max()
        fast = sp > fast_thr
        n_fast[k] = int(fast.sum()); ever_fast |= fast
        d = ff["pos"][anchors] - x0_anch
        anch_disp[k] = np.linalg.norm(d, axis=1).max()
        anch_v[k] = np.linalg.norm(ff["vel"][anchors], axis=1).max()
        anch_vz[k] = ff["vel"][anchors, 2].mean()
        # angular momentum (solid + beads)
        Ptot = (m[:, None] * V).sum(0) + (bmass[:, None] * ff["vel"]).sum(0)
        Lo = (m[:, None] * np.cross(X - com0, V)).sum(0) + (bmass[:, None] * np.cross(ff["pos"] - com0, ff["vel"])).sum(0)
        Rcm = ((m[:, None] * X).sum(0) + (bmass[:, None] * ff["pos"]).sum(0)) / Mtot0
        Lc = Lo - np.cross(Rcm - com0, Ptot)   # = sum m (r-Rcm) x v
        if k == 0:
            Lo0, Lc0 = Lo.copy(), Lc.copy()
        dL0[k] = np.linalg.norm(Lo - Lo0); dLc[k] = np.linalg.norm(Lc - Lc0)
        # bonds
        cur = {tuple(sorted(b)) for b in zip(ff["bond_i"].tolist(), ff["bond_j"].tolist())}
        for b in bonds0:
            if b not in cur and b not in break_frame:
                break_frame[b] = k
        # pierced vs carried: intact beads (>=1 bond) under the cube footprint
        has_bond = np.zeros(len(bmass), dtype=bool)
        for i, j in cur:
            has_bond[i] = has_bond[j] = True
        xlo, xhi = X[:, 0].min(), X[:, 0].max(); ylo, yhi = X[:, 1].min(), X[:, 1].max()
        bp = ff["pos"]
        inside = has_bond & (bp[:, 0] >= xlo) & (bp[:, 0] <= xhi) & (bp[:, 1] >= ylo) & (bp[:, 1] <= yhi)
        foot_n[k] = int(inside.sum())
        if inside.any():
            foot_gap[k] = float(zmin[k] - bp[inside, 2].max())
        foot_above[k] = int((inside & (bp[:, 2] > zmin[k] - args.margin)).sum())

    # trajectory
    s_vz = sgn * vz   # >0 while moving in initial direction
    cross = np.nonzero(s_vz < 0)[0]
    rebound_t = float(t[cross[0]]) if len(cross) else None
    imin = int(np.argmin(np.abs(vz)))
    arrested = bool(len(cross) or abs(vz[-1]) < 0.05 * v0)
    pier = np.nonzero(foot_above > 0)[0]
    plane_z0 = np.array([f0["pos"][grp0 == g, 2].mean() for g in range(int(grp0.max()) + 1)])
    res["projectile"] = dict(
        groups_plane_z0=plane_z0.tolist(),
        groups_planes_below_final_zmin=int((plane_z0 > zmin[-1]).sum()),
        footprint_min_gap_final=None if np.isnan(foot_gap[-1]) else float(foot_gap[-1]),
        v0=v0, v0_vec=v0vec.tolist(), m_proj=mp, cube_side=L_char,
        z_min_initial=float(zmin[0]), z_min_final=float(zmin[-1]), z_min_lowest=float(zmin.min()),
        z_mean_initial=float(zmean[0]), z_mean_final=float(zmean[-1]),
        vz_final=float(vz[-1]),
        rebound_time=rebound_t, min_abs_vz=float(abs(vz[imin])), min_abs_vz_time=float(t[imin]),
        final_over_v0=float(vz[-1] / v0) if v0 else None,
        arrested=arrested,
        arrested_reason=("sign change" if len(cross) else ("|vz_end|<0.05|v0|" if arrested else "none")),
        pierced_ever=bool(len(pier)), pierced_first_time=float(t[pier[0]]) if len(pier) else None,
        footprint_beads_final=int(foot_n[-1]), footprint_beads_above_zmin_final=int(foot_above[-1]),
        max_footprint_beads_above_zmin=int(foot_above.max()),
        verdict_pierced_or_carried="pierced" if foot_above[-1] > 0 else "carried")
    res["anchors"] = dict(n=int(len(anchors)), max_disp=float(anch_disp.max()), final_disp_max=float(anch_disp[-1]),
                          max_speed=float(anch_v.max()), final_mean_vz=float(anch_vz[-1]))
    res["angular_momentum"] = dict(
        scale_m_v0_L=Lscale, max_dL_about_COM0_rel=float(dL0.max() / Lscale),
        max_dL_about_current_COM_rel=float(dLc.max() / Lscale),
        L0_about_COM0_norm_rel=float(np.linalg.norm(Lo0) / Lscale),
        final_dL_about_COM0_rel=float(dL0[-1] / Lscale), final_dL_about_current_COM_rel=float(dLc[-1] / Lscale))
    # broken bonds
    grp_counts, pos_counts = {}, {}
    times_b = []
    for b, kf in break_frame.items():
        g = int(grp0[b[0]]); grp_counts[g] = grp_counts.get(g, 0) + 1
        cid, sidx = bpos.get(b, (-1, -1))
        pos_counts[sidx] = pos_counts.get(sidx, 0) + 1
        times_b.append(float(t[kf]))
    n_groups = int(grp0.max()) + 1
    res["broken_bonds"] = dict(
        n_broken=len(break_frame), n_total0=len(bonds0),
        group_def=(f"object_id@0 // {K}" if K else "object_id@0"),
        per_group=[grp_counts.get(g, 0) for g in range(n_groups)],
        chain_beads=sorted(set(chain_len)),
        position_hist={str(k2): v for k2, v in sorted(pos_counts.items())},
        first_break_time=min(times_b) if times_b else None,
        broken_list=[dict(bond=list(b), group=int(grp0[b[0]]), seg=bpos.get(b, (-1, -1))[1],
                          frame=kf, time=float(t[kf])) for b, kf in sorted(break_frame.items(), key=lambda x: x[1])])
    res["fragments"] = dict(
        max_bead_speed=float(bead_vmax.max()), max_bead_speed_over_v0=float(bead_vmax.max() / v0),
        max_bead_speed_time=float(t[int(np.argmax(bead_vmax))]),
        thr_factor=args.fast_factor, n_beads_over_thr_max_single_frame=int(n_fast.max()),
        n_beads_ever_over_thr=int(ever_fast.sum()))
    res["series"] = dict(time=t.tolist(), z_min=zmin.tolist(), z_mean=zmean.tolist(), vz_mean=vz.tolist(),
                         max_bead_speed=bead_vmax.tolist(), anchor_max_disp=anch_disp.tolist())
    return res


def fmt(x, f="{:.4g}"):
    return "-" if x is None else f.format(x)


def rows_for(r):
    p, a, L, b, fr, e, pe = (r["projectile"], r["anchors"], r["angular_momentum"], r["broken_bonds"],
                             r["fragments"], r["energy"], r["penetration"])
    return [
        ("v0, m/s (frame 0)", fmt(p["v0"])),
        ("z_min initial / final / lowest, mm", f"{p['z_min_initial']*1e3:.3f} / {p['z_min_final']*1e3:.3f} / {p['z_min_lowest']*1e3:.3f}"),
        ("z_mean initial -> final, mm", f"{p['z_mean_initial']*1e3:.3f} -> {p['z_mean_final']*1e3:.3f}"),
        ("final vz_mean, m/s", fmt(p["vz_final"], "{:.3f}")),
        ("final vz / v0", fmt(p["final_over_v0"], "{:.3f}")),
        ("rebound time (vz crosses 0), s", fmt(p["rebound_time"], "{:.4e}")),
        ("min abs(vz_mean), m/s @ t", f"{p['min_abs_vz']:.3f} @ {p['min_abs_vz_time']:.3e}"),
        ("arrested (sign change or <5% v0)", f"{p['arrested']} ({p['arrested_reason']})"),
        ("footprint beads above cube z_min (final / max)", f"{p['footprint_beads_above_zmin_final']} / {p['max_footprint_beads_above_zmin']}"),
        ("pierced ever (first t)", f"{p['pierced_ever']} ({fmt(p['pierced_first_time'], '{:.3e}')})"),
        ("fiber groups whose frame-0 plane is above final z_min", f"{p['groups_planes_below_final_zmin']} of {len(p['groups_plane_z0'])}"),
        ("footprint min gap z_min - z_bead (final), mm", fmt(None if p["footprint_min_gap_final"] is None else p["footprint_min_gap_final"]*1e3, "{:.3f}")),
        ("verdict at end: pierced / carried", p["verdict_pierced_or_carried"]),
        ("anchors n", str(a["n"])),
        ("anchor max displacement, mm", f"{a['max_disp']*1e3:.3f}"),
        ("anchor max speed, m/s", fmt(a["max_speed"], "{:.3f}")),
        ("anchor final mean vz, m/s", fmt(a["final_mean_vz"], "{:.3f}")),
        ("max abs(dL)/(m v0 L) about COM0", fmt(L["max_dL_about_COM0_rel"], "{:.3e}")),
        ("max abs(dL)/(m v0 L) about current COM", fmt(L["max_dL_about_current_COM_rel"], "{:.3e}")),
        ("abs(L0)/(m v0 L) about COM0", fmt(L["L0_about_COM0_norm_rel"], "{:.3e}")),
        ("bonds broken / total", f"{b['n_broken']} / {b['n_total0']}"),
        (f"broken per group ({b['group_def']})", ",".join(map(str, b["per_group"]))),
        (f"broken by segment idx (chain beads {b['chain_beads']})", ", ".join(f"{k}:{v}" for k, v in b["position_hist"].items()) or "-"),
        ("first break time, s", fmt(b["first_break_time"], "{:.4e}")),
        ("max bead speed, m/s (t)", f"{fr['max_bead_speed']:.3f} ({fr['max_bead_speed_time']:.3e})"),
        ("max bead speed / v0", fmt(fr["max_bead_speed_over_v0"], "{:.3f}")),
        (f"beads > {fr['thr_factor']} v0: max in frame / ever", f"{fr['n_beads_over_thr_max_single_frame']} / {fr['n_beads_ever_over_thr']}"),
        ("max E/E0 (t)", f"{e['max_E_ratio']:.4f} ({e['max_E_ratio_time']:.3e})"),
        ("final E/E0", fmt(e["final_E_ratio"], "{:.4f}")),
        ("max abs(P-P0)/abs(P0)", fmt(e["max_momentum_drift_rel"], "{:.3e}")),
        ("max delta/d fiber-fiber", fmt(pe["max_ff_delta_over_d"], "{:.4f}")),
        ("max delta/d fiber-solid", fmt(pe["max_fs_delta_over_d"], "{:.4f}")),
        ("pass-through events", str(pe["n_passthrough"])),
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="DIR:KIND:LABEL")
    ap.add_argument("--out-md"); ap.add_argument("--out-json")
    ap.add_argument("--fiber-E", dest="fiber_E", type=float); ap.add_argument("--fiber-d", dest="fiber_d", type=float)
    ap.add_argument("--fiber-rho", dest="fiber_rho", type=float)
    ap.add_argument("--layer-of-object", dest="layer_of_object", type=int)
    ap.add_argument("--anchor-factor", type=float, default=50.0)
    ap.add_argument("--solid-rho0", type=float, default=2700.0)
    ap.add_argument("--cube-side", type=float, default=None, help="L_char override (default: frame-0 x-extent + dx)")
    ap.add_argument("--fast-factor", type=float, default=1.5)
    ap.add_argument("--margin", type=float, default=0.0, help="pierced if footprint bead z > z_min - margin")
    args = ap.parse_args()

    results, labels = [], []
    for spec in args.run:
        parts = spec.rsplit(":", 2)
        if len(parts) != 3 or parts[1] not in KIND_DEFAULTS:
            raise SystemExit(f"bad --run {spec!r}; want DIR:KIND:LABEL with KIND in {list(KIND_DEFAULTS)}")
        d, kind, label = parts
        r = analyse_run(d, kind, args)
        r["label"] = label
        results.append(r); labels.append(label)
    tabs = [rows_for(r) for r in results]
    md = ["| metric | " + " | ".join(labels) + " |", "|---|" + "---|" * len(labels)]
    for i, (name, _) in enumerate(tabs[0]):
        md.append(f"| {name} | " + " | ".join(t[i][1] for t in tabs) + " |")
    text = "\n".join(md)
    print(text)
    if args.out_md:
        open(args.out_md, "w", encoding="utf-8").write(text + "\n")
    if args.out_json:
        json.dump({r["label"]: r for r in results}, open(args.out_json, "w"), indent=1)


if __name__ == "__main__":
    main()
