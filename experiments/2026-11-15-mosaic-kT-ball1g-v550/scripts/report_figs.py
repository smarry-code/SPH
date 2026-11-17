#!/usr/bin/env python3
"""Рисунки отчёта серии «мозаика kT, шар 1 г, 550 м/с» -> fig/*.png.

  vr_vs_n.png       -- V_R по числу слоёв n (вариант B) + 16T A;
  vz_t.png          -- v_z(t) и z_c(t) шара для всех выстрелов;
  energy.png        -- куда ушла энергия шара E0 (на конец выстрела);
  damage.png        -- разрывы нитей и стежков по n;
  sections.png      -- сечения y = 0 пакета в момент t ~ 60 мкс: 14T B, 15T B, 16T B, 16T A.
Запуск: python scripts/report_figs.py (из каталога эксперимента). Читает runs/*/shots.csv, monitor.csv, energy_log.csv,
fiber_output_*.txt; незаконченный выстрел (нет shots.csv) рисуется по monitor.csv как «идёт».
"""
import csv
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(HERE, "fig")
V0, M = 550.0, 1e-3
E0 = 0.5 * M * V0 ** 2
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE = {8: "#9ec5f4", 12: "#6da7ec", 14: "#3987e5", 15: "#256abf", 16: "#184f95", 17: "#0d366b"}
ORANGE = "#eb6834"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb"})


def shots():
    out = []
    for o in sorted(glob.glob(os.path.join(HERE, "runs", "*T_[AB]"))):
        n, v = re.match(r"(\d+)T_([AB])", os.path.basename(o)).groups()
        d = glob.glob(os.path.join(o, "*T_base_*"))[0]
        row = None
        if os.path.exists(os.path.join(o, "shots.csv")):
            row = list(csv.DictReader(open(os.path.join(o, "shots.csv"), encoding="utf-8")))[-1]
        mon = [r for r in csv.DictReader(open(os.path.join(d, "monitor.csv"))) if r.get("gap_mean")]
        en = list(csv.DictReader(open(os.path.join(d, "energy_log.csv"))))
        out.append(dict(n=int(n), v=v, dir=d, row=row, mon=mon, en=en))
    return sorted(out, key=lambda s: (s["v"] != "B", s["n"]))


def label(s):
    return "%dT %s%s" % (s["n"], s["v"], "" if s["row"] else " (идёт)")


def color(s):
    return BLUE[s["n"]] if s["v"] == "B" else ORANGE


def fig_vr(S):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    B = [s for s in S if s["v"] == "B" and s["row"]]
    n = [s["n"] for s in B]
    vr = [float(s["row"]["V_R"]) for s in B]
    ax.plot(n, vr, "-", color=BLUE[14], zorder=2)
    for s, x, y in zip(B, n, vr):
        perf = s["row"]["outcome"] == "perforated"
        ax.plot(x, y, "o", ms=9, mfc=BLUE[14] if perf else "#fcfcfb", mec=BLUE[14], mew=2, zorder=3)
        ax.annotate("%.0f" % y if perf else "стоп", (x, y), textcoords="offset points", xytext=(0, 10), ha="center",
                    color=INK2, fontsize=9)
    A = [s for s in S if s["v"] == "A" and s["row"]]
    for s in A:
        y = float(s["row"]["V_R"])
        ax.plot(s["n"], y, "s", ms=9, color=ORANGE, zorder=3)
        ax.annotate("16T A: %.0f" % y, (s["n"], y), textcoords="offset points", xytext=(0, -18), ha="center",
                    color=INK2, fontsize=9)
    ax.axvspan(14, 15, color="#f0efec", zorder=0)
    ax.text(14.5, 520, "граница\n14 | 15", ha="center", va="top", fontsize=9, color=INK2)
    ax.plot([], [], "o", color=BLUE[14], ms=8, label="вариант B, пробит")
    ax.plot([], [], "o", mfc="#fcfcfb", mec=BLUE[14], mew=2, ms=8, label="вариант B, не пробит (V_R = 0)")
    ax.plot([], [], "s", color=ORANGE, ms=8, label="вариант A")
    ax.set_xlabel("число слоёв n")
    ax.set_ylabel("скорость шара за мишенью V_R, м/с")
    ax.set_ylim(-20, 560)
    ax.set_xticks(sorted(set(n + [s["n"] for s in S])))
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    ax.set_title("Остаточная скорость шара (1 г, V0 = 550 м/с)", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "vr_vs_n.png"), dpi=140)


def fig_vz(S):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.6))
    for s in S:
        t = np.array([float(r["t"]) for r in s["mon"]]) * 1e6
        vz = np.array([float(r["v_z"]) for r in s["mon"]])
        zc = np.array([float(r["z_c"]) for r in s["mon"]]) * 1e3
        ls = "--" if s["v"] == "A" else "-"
        a.plot(t, vz, ls, color=color(s), label=label(s))
        b.plot(t, zc, ls, color=color(s), label=label(s))
    a.axhline(0, color=INK2, lw=0.8)
    a.set_xlabel("t, мкс"); a.set_ylabel("v_z шара, м/с")
    a.set_title("Скорость шара", loc="left", fontsize=11)
    b.set_xlabel("t, мкс"); b.set_ylabel("z_c центра шара, мм")
    b.set_title("Положение шара (лицевая сторона пакета — z = 0)", loc="left", fontsize=11)
    b.legend(fontsize=8, frameon=False, loc="lower left")
    a.set_xlim(0, None); b.set_xlim(0, None)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "vz_t.png"), dpi=140)


def fig_energy(S):
    S = [s for s in S if s["row"]]
    parts = [("энергия шара", lambda e: float(e["KE_solid"]), BLUE[14]),
             ("кинетика и упругость ткани", lambda e: float(e["E_mech"]) - float(e["KE_solid"]), "#86b6ef"),
             ("трение", lambda e: float(e["D_friction"]), ORANGE),
             ("контакт нитей", lambda e: float(e["D_contact_fiber"]), "#eda100"),
             ("прочее (XSPH, разрывы, пластика)", lambda e: sum(float(e[k]) for k in
              ("D_xsph", "D_break", "Q_av", "W_plastic", "W_p_fiber", "W_v_fiber", "S_u_clamp")), "#c3c2b7")]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    names = [label(s) for s in S]
    left = np.zeros(len(S))
    for nm, f, c in parts:
        w = np.array([f(s["en"][-1]) for s in S]) / E0 * 100
        ax.barh(names, w, left=left, color=c, edgecolor="#fcfcfb", linewidth=2, label=nm, height=0.62)
        for i, (l, x) in enumerate(zip(left, w)):
            if x > 7:
                ax.text(l + x / 2, i, "%.0f" % x, ha="center", va="center", fontsize=8, color="#0b0b0b")
        left += w
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("доля начальной энергии шара E0 = 151 Дж, %")
    ax.grid(axis="y", visible=False)
    ax.legend(fontsize=8, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    ax.set_title("Куда ушла энергия шара (на конец выстрела)", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "energy.png"), dpi=140)


def fig_damage(S):
    S = [s for s in S if s["row"]]
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.0))
    names = [label(s) for s in S]
    y = np.arange(len(S))
    bv = np.array([int(s["row"]["n_broken_v"]) for s in S])
    bh = np.array([int(s["row"]["n_broken_h"]) for s in S])
    a.barh(y - 0.18, bv, 0.34, color=BLUE[14], label="основа")
    a.barh(y + 0.18, bh, 0.34, color=ORANGE, label="уток")
    a.set_yticks(y, names); a.invert_yaxis(); a.grid(axis="y", visible=False)
    a.set_xlabel("порвано нитей"); a.legend(frameon=False, fontsize=9)
    a.set_title("Разрывы нитей", loc="left", fontsize=11)
    ns = np.array([int(s["row"]["n_stitch"]) for s in S])
    nb = np.array([int(s["row"]["n_stitch_broken_end"]) for s in S])
    b.barh(y, nb, 0.6, color=[color(s) for s in S])
    for i in range(len(S)):
        b.text(nb[i] + 15, i, "%d из %d (%.1f %%)" % (nb[i], ns[i], 100 * nb[i] / ns[i]), va="center", fontsize=8,
               color=INK2)
    b.set_yticks(y, names); b.invert_yaxis(); b.grid(axis="y", visible=False)
    b.set_xlim(0, nb.max() * 1.6)
    b.set_xlabel("порвано стежков")
    b.set_title("Разрывы стежков", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "damage.png"), dpi=140)


def read_frame(path, n):
    pos = np.zeros((n, 3))
    with open(path) as f:
        head = f.readline().split()
        t = float(head[1])
        for ln in f:
            if ln.strip() == "PARTICLES":
                break
        arr = np.loadtxt(f, max_rows=n, usecols=(0, 3, 4, 5))
    pos[arr[:, 0].astype(int)] = arr[:, 1:]
    return t, pos


def fig_sections(S, t_want=60e-6):
    pick = [s for s in S if (s["n"], s["v"]) in ((14, "B"), (15, "B"), (16, "B"), (16, "A"))]
    fig, axs = plt.subplots(2, 2, figsize=(13, 7.4), sharex=True, sharey=True)
    axs = axs.ravel()
    for ax, s in zip(axs, pick):
        import json
        j = json.load(open(os.path.join(s["dir"], "guo2019_fabric.cut.json"), encoding="utf-8"))
        n = j["n_beads"]
        ply = np.zeros(n, int)
        for r in j["ranges"]:
            ply[r["bead_first"]:r["bead_first"] + r["bead_count"]] = r["ply"]
        fr = sorted(glob.glob(os.path.join(s["dir"], "fiber_output_*.txt")))
        times = []
        for p in fr:
            with open(p) as f:
                times.append(float(f.readline().split()[1]))
        k = int(np.argmin(np.abs(np.array(times) - t_want)))
        t, pos = read_frame(fr[k], n)
        sel = (np.abs(pos[:, 1]) < 0.6e-3) & (np.abs(pos[:, 0]) < 20e-3)
        cmap = plt.get_cmap("Blues")
        nl = ply.max() + 1
        ax.scatter(pos[sel, 0] * 1e3, pos[sel, 2] * 1e3, s=2.5, c=[cmap(0.35 + 0.6 * p / max(nl - 1, 1))
                                                                    for p in ply[sel]], lw=0)
        mon = s["mon"]
        r = min(mon, key=lambda r: abs(float(r["t"]) - t))
        zc = float(r["z_c"]) * 1e3
        ax.add_patch(plt.Circle((0, zc), 6.253 / 2, fc="#eb683440", ec=ORANGE, lw=1.5))
        ax.set_aspect("equal"); ax.grid(False)
        ax.set_ylim(-36, 6)
        ax.set_ylabel("z, мм")
        res = s["row"]["outcome"] if s["row"] else "идёт"
        res = {"perforated": "пробит", "not_perforated": "не пробит"}.get(res, res)
        ax.set_title("%s — %s;  t = %.0f мкс, v_z = %.0f м/с" % (label(s), res, t * 1e6, float(r["v_z"])),
                     loc="left", fontsize=10)
    for ax in axs[2:]:
        ax.set_xlabel("x, мм  (сечение |y| < 0,6 мм; светлее — верхние слои)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "sections.png"), dpi=130)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    S = shots()
    fig_vr(S); fig_vz(S); fig_energy(S); fig_damage(S); fig_sections(S)
    print("ok", [label(s) for s in S])
