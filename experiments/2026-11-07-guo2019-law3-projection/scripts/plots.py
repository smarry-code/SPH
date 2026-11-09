#!/usr/bin/env python3
"""plots.py <root> <frames_csv_dir> --out DIR
PNG-графики серии. Кадры (frames_csv_dir/<RUN>.csv из frames.py): δ_ff/d, δ_fs/d, n_tunnel, проходы — C1..C5.
monitor.csv: proj_resid, W_proj_pos_max — C1, C4, C5. energy_log.csv: E_mech/E0, D_contact_fiber/E0, residual/E0 — C1, C2, C3."""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import read_csv_dedup

LABEL = {"C1": "C1 закон 3 + проекция K=8", "C2": "C2 закон 3 без проекции", "C3": "C3 закон 4",
         "C4": "C4 проекция K=2", "C5": "C5 проекция K=8, без тёплого старта"}
COL = {"C1": "#1f77b4", "C2": "#d62728", "C3": "#2ca02c", "C4": "#ff7f0e", "C5": "#9467bd"}
MS = 1e3


def setup(ax, ylabel, title, logy=False):
    ax.set_xlabel("t, мс"); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.grid(alpha=0.3)
    if logy:
        ax.set_yscale("log")
    ax.legend(fontsize=8)


def save(fig, out, name):
    fig.tight_layout(); fig.savefig(os.path.join(out, name), dpi=130); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root"); ap.add_argument("frames_dir"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    # --- кадры ---
    fr = {}
    for r in ("C1", "C2", "C3", "C4", "C5"):
        p = os.path.join(a.frames_dir, r + ".csv")
        if os.path.isfile(p):
            fr[r] = np.genfromtxt(p, delimiter=",", names=True)
    items = [("dff_max", "δ_ff / d (максимум по кадру)", "Проникание нить–нить", "frames_delta_ff.png", 0.3),
             ("dfs_max", "δ_fs / d (максимум по кадру)", "Проникание нить–шар", "frames_delta_fs.png", 0.3),
             ("n_tunnel", "n_tunnel, бусин", "Туннелирование в шар (по кадрам)", "frames_n_tunnel.png", 5),
             ("n_cross_flip", "число перекрестий со сменой знака", "Проходы (а): смена знака z_уток − z_основа в плане", "frames_pass_cross.png", None),
             ("dff_n_gt09", "число пар сегментов", "Проходы (б): пары с δ_ff/d > 0,9", "frames_pass_ff09.png", None)]
    for col, yl, ti, name, thr in items:
        fig, ax = plt.subplots(figsize=(7, 4.2))
        for r, A in fr.items():
            ax.plot(A["t"] * MS, A[col], "o-", ms=3, color=COL[r], label=LABEL[r])
        if thr is not None:
            ax.axhline(thr, color="k", ls="--", lw=0.8, label=f"порог ТЗ {thr}")
        setup(ax, yl, ti)
        save(fig, a.out, name)
    # --- монитор (проекция) ---
    mo = {}
    for r in ("C1", "C4", "C5"):
        if os.path.isfile(os.path.join(a.root, r, "monitor.csv")):
            mo[r] = read_csv_dedup(os.path.join(a.root, r, "monitor.csv"))[1]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for r, m in mo.items():
        act = (np.nan_to_num(m["proj_n_ff"]) + np.nan_to_num(m["proj_n_fs"]) > 0) & (m["proj_resid"] > 0)
        ax.plot(m["t"][act] * MS, m["proj_resid"][act], ".", ms=2, color=COL[r], label=LABEL[r])
    setup(ax, "proj_resid (макс. в окне)", "Невязка проекции", logy=True)
    save(fig, a.out, "monitor_proj_resid.png")
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for r, m in mo.items():
        ax.plot(m["t"] * MS, np.nan_to_num(m["W_proj_pos_max"]), ".-", ms=3, lw=0.6, color=COL[r], label=LABEL[r])
    setup(ax, "W_proj_pos_max (колонка monitor, окно между строками)", "Положительная работа проекции")
    ax.set_yscale("symlog", linthresh=1e-8)
    ax.legend(fontsize=8)
    save(fig, a.out, "monitor_W_proj_pos_max.png")
    # --- энергия ---
    en = {}
    for r in ("C1", "C2", "C3"):
        if os.path.isfile(os.path.join(a.root, r, "energy_log.csv")):
            en[r] = read_csv_dedup(os.path.join(a.root, r, "energy_log.csv"))[1]
    for col, yl, name in (("E_mech", "E_mech / E0", "energy_Emech.png"),
                          ("D_contact_fiber", "D_contact_fiber / E0", "energy_Dcf.png"),
                          ("residual", "residual / E0", "energy_residual.png")):
        fig, ax = plt.subplots(figsize=(7, 4.2))
        for r, e in en.items():
            ax.plot(e["time"] * MS, e[col] / e["E_mech"][0], color=COL[r], lw=1.2, label=LABEL[r])
        setup(ax, yl, yl + " (энергетический баланс)")
        save(fig, a.out, name)
    print("PNG ->", a.out, sorted(os.listdir(a.out)))


if __name__ == "__main__":
    main()
