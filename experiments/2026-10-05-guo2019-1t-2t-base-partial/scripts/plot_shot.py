#!/usr/bin/env python3
"""Графики одного выстрела по monitor.csv (колонки n_broken_v_p<k>, n_broken_h_p<k>, w_f_p<k>).

Выход: breaks_by_ply.svg (вертикальные и горизонтальные нити по слоям, вертикали -- первый разрыв слоя),
       deflection_by_ply.svg (прогиб слоёв, мм, полоса 25-27 мм статьи).
Usage:
    python plot_shot.py --shot-dir server_series/4T/4T_scaled_v268.7_r0 [--out DIR] [--png]
"""
from __future__ import annotations
import argparse, os, re
import numpy as np
from viz_common import *


def n_plies(m):
    k = 1
    while f"n_broken_v_p{k}" in m:
        k += 1
    return k - 1


def title(shot_dir):
    p = os.path.join(shot_dir, "shot.txt")
    if os.path.isfile(p):
        t = open(p, encoding="utf-8").read().split("\n")[0]
        return " ".join(w for w in t.split() if re.match(r"(config|level|V0|mu|outcome|V_R)=", w))
    return os.path.basename(os.path.normpath(shot_dir))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shot-dir", required=True)
    ap.add_argument("--out")
    ap.add_argument("--png", action="store_true")
    a = ap.parse_args()
    out = a.out or os.path.join(a.shot_dir, "plots")
    m = read_monitor(os.path.join(a.shot_dir, "monitor.csv"))
    t = m["t"] * 1e6
    n = n_plies(m)
    if n == 0:
        raise SystemExit("в monitor.csv нет колонок n_broken_v_p<k> (однослойный выстрел?)")
    tt = title(a.shot_dir)
    firsts = {}
    for k in range(1, n + 1):
        hit = np.where(m[f"n_broken_v_p{k}"] + m[f"n_broken_h_p{k}"] > 0)[0]
        firsts[k] = hit[0] if len(hit) else None

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True)
    for ax, (pref, lab) in zip(axs, (("n_broken_v_p", "вертикальных (зажатых)"),
                                     ("n_broken_h_p", "горизонтальных (свободных)"))):
        for k in range(1, n + 1):
            c = LAYER_COLORS[(k - 1) % len(LAYER_COLORS)]
            ax.step(t, m[f"{pref}{k}"], where="post", color=c, label=f"слой {k}")
            hit = np.where(m[f"{pref}{k}"] > 0)[0]
            if len(hit):
                ax.axvline(t[hit[0]], color=c, ls=":", lw=1)
        ax.set_xlabel("время t, мкс")
        ax.set_ylabel("число порванных нитей, шт.")
        ax.set_title(f"порванные нити {lab}")
        ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(tt + "   (пунктир: первый разрыв нитей этого типа в слое; слой 1 — лицевой)", fontsize=9)
    save(fig, out, "breaks_by_ply", a.png)

    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.axhspan(25, 27, color="k", alpha=0.12, label="статья: 25–27 мм перед разрушением")
    for k in range(1, n + 1):
        c = LAYER_COLORS[(k - 1) % len(LAYER_COLORS)]
        y = m[f"w_f_p{k}"] * 1e3
        ax.plot(t, y, color=c, label=f"слой {k}")
        if firsts[k] is not None:
            ax.plot(t[firsts[k]], y[firsts[k]], "o", color=c, ms=6)
    ax.plot([], [], "ko", label="первый разрыв слоя (любой нити)")
    ax.set_xlabel("время t, мкс")
    ax.set_ylabel("прогиб слоя w_f, мм")
    ax.set_title(tt, fontsize=9)
    ax.legend(fontsize=8, loc="upper left")
    save(fig, out, "deflection_by_ply", a.png)

    for k in range(1, n + 1):
        y = m[f"w_f_p{k}"] * 1e3
        j = firsts[k]
        print(f"слой {k}: max w_f = {y.max():.2f} мм; первый разрыв: "
              + (f"{t[j]:.2f} мкс, w_f = {y[j]:.2f} мм" if j is not None else "нет")
              + f"; порвано v/h в конце: {int(m[f'n_broken_v_p{k}'][-1])}/{int(m[f'n_broken_h_p{k}'][-1])}")


if __name__ == "__main__":
    main()
