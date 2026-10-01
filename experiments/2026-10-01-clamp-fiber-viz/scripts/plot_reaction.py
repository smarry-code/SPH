#!/usr/bin/env python3
"""R(t) опоры против аналитики для EF1, EF1-b, EF2 (e2e закрепления нитей).

Читает data/<прогон>/reaction.csv (колонки step,t,Rs_x..Rf_z; R = R_solid + R_fiber),
константы берёт из tests/clamp_case.hpp (секция fib), считает аналитику заново и
пишет ef1_reaction.svg, ef2_reaction.svg и data/metrics.txt (печатается и в stdout).
Запуск из любого места: python scripts/plot_reaction.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "clamp-fiber-viz"
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")

# --- clamp_case::fib -------------------------------------------------------
E, d, rho, eps_break, l0 = 1.0e7, 1.0e-3, 1000.0, 0.05, 0.0078125
A = np.pi * d * d / 4
c = np.sqrt(E / rho)
m = rho * A * l0
k = E * A / l0
sqrt_km = np.sqrt(k * m)            # = A sqrt(E rho)
n_ef1, n_c = 50, 3
L = (n_ef1 - n_c + 1) * l0
n_ef2, eps0, amp_over_l = 20, 0.3 * eps_break, 0.01

C = dict(R_x="#1f77b4", R_y="#d95f02", ref="#444444")


def load(run):
    a = np.genfromtxt(os.path.join(DATA, run, "reaction.csv"), delimiter=",", names=True)
    R = {ax: a["Rs_" + ax] + a["Rf_" + ax] for ax in "xyz"}
    t = a["t"]
    tm = t - 0.5 * np.diff(np.r_[0.0, t])   # R шага n действует на [t_{n-1}, t_n]: как step_mid() в тесте
    return tm, R


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.5)


out = []

# ---------------- EF1 / EF1-b ----------------
t1, R1 = load("EF1_pc_ef1")
v0 = 0.5 * eps_break * c
plateau = -sqrt_km * v0
tb, Rb = load("EF1_pc_b_ef1b")
v0b = 1.5 * eps_break * c
plateau_b = -sqrt_km * v0b
mask = (t1 >= 0.2 * L / c) & (t1 < 1.8 * L / c)
mean_pl = R1["x"][mask].mean()
below = np.where((t1 >= 1.5 * L / c) & (np.abs(R1["x"]) < 0.1 * abs(plateau)))[0]
t_drop = t1[below[0]]
out.append(f"EF1: -sqrt(km) v0 = {plateau:.6e} Н; среднее R_x на [0.2,1.8]L/c = {mean_pl:.6e} Н "
           f"({(mean_pl/plateau-1)*100:+.3e} %); 2L/c = {2*L/c:.4e} с; R_x < 10 % плато с t_mid = {t_drop:.4e} с "
           f"({(t_drop/(2*L/c)-1)*100:+.3f} %); max|R_x| после 2.2L/c = {np.abs(R1['x'][t1>2.2*L/c]).max():.1e} Н")
out.append(f"EF1-b: -sqrt(km) v0b = {plateau_b:.6e} Н; R_x на 1-м шаге = {Rb['x'][0]:.6e}, на 2-м = {Rb['x'][1]:.6e}, "
           f"ненулевых шагов R_x: {(Rb['x']!=0).sum()} из {len(tb)}; max|R_x| после 2-го шага = {np.abs(Rb['x'][2:]).max():.1e}")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.9), constrained_layout=True)
a1.plot(t1 * 1e3, R1["x"], color=C["R_x"], lw=1.6, label="R_x, прогон EF1 (PC)")
a1.axhline(plateau, color=C["ref"], ls="--", lw=1.1, label=f"-√(km)·v₀ = {plateau:.4f} Н")
a1.axvline(2 * L / c * 1e3, color=C["ref"], ls=":", lw=1.1, label=f"2L/c = {2*L/c*1e3:.2f} мс")
a1.axvline(2.2 * L / c * 1e3, color="#999999", ls=":", lw=0.8)
a1.text(2.2 * L / c * 1e3, plateau * 0.5, "  2.2 L/c: R ≡ 0", fontsize=8, color="#666666")
a1.set(xlabel="t, мс", ylabel="R_x, Н", title="EF1: рывок нити от зажима, v₀ = 0.5 ε_break c")
a1.legend(fontsize=8, loc="center left", bbox_to_anchor=(0.02, 0.62), frameon=False)
style(a1)
a2.plot(np.arange(1, len(tb) + 1)[:12], Rb["x"][:12], "o-", color=C["R_x"], ms=4, lw=1.2, label="R_x, прогон EF1-b (PC)")
a2.axhline(plateau_b, color=C["ref"], ls="--", lw=1.1, label=f"-√(km)·v₀ (v₀ = 1.5 ε_break c) = {plateau_b:.3f} Н")
a2.set(xlabel="номер шага", ylabel="R_x, Н", title="EF1-b: разрыв связи (2,3) на шаге 2, первые 12 шагов")
a2.legend(fontsize=8, loc="center right", frameon=False)
style(a2)
fig.savefig(os.path.join(ROOT, "ef1_reaction.svg"), format="svg")
plt.close(fig)

# ---------------- EF2 ----------------
t2, R2 = load("EF2_pc_ef2")
N = E * A * eps0
l = l0 * (1 + eps0)
w1 = 2 * np.sqrt(N / (m * l)) * np.sin(np.pi / (2 * n_ef2))
T = 2 * np.pi / w1
Ry_amp = 2 * N * (amp_over_l * l) * np.sin(np.pi / n_ef2) / l
# амплитуда и фаза R_y: подгонка a sin(w1 t) + b cos(w1 t)
M = np.c_[np.sin(w1 * t2), np.cos(w1 * t2)]
(sa, sb), *_ = np.linalg.lstsq(M, R2["y"], rcond=None)
fit_amp = np.hypot(sa, sb)
res = R2["y"] - M @ np.array([sa, sb])
out.append(f"EF2: N = {N:.6e} Н; T = 2π/ω1 = {T:.6e} с; ω1 = {w1:.6e} 1/с; амплитуда R_y аналит. = {Ry_amp:.6e} Н; "
           f"подгонка на ω1: {fit_amp:.6e} Н ({(fit_amp/Ry_amp-1)*100:+.3e} %), sa = {sa:.4e}, sb = {sb:.4e}; "
           f"max|невязка| = {np.abs(res).max():.2e} Н ({np.abs(res).max()/Ry_amp*100:.3f} % амплитуды); "
           f"max|R_x| суммарный = {np.abs(R2['x']).max():.2e} Н; max|R_z| = {np.abs(R2['z']).max():.2e} Н; t_end = {t2[-1]:.6e} с")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.9), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.6, 1]))
a1.plot(t2 * 1e3, R2["y"] * 1e4, color=C["R_y"], lw=1.6, label="R_y из прогона EF2 (PC)")
tt = np.linspace(0, t2[-1], 800)
a1.plot(tt * 1e3, (sa * np.sin(w1 * tt) + sb * np.cos(w1 * tt)) * 1e4, color=C["ref"], ls="--", lw=1.0,
        label=f"синус на ω₁ дискретной цепочки, T = {T*1e3:.3f} мс")
for j in range(1, 5):
    a1.axvline(j * T / 2 * 1e3, color="#999999", ls=":", lw=0.7)
a1.set(xlabel="t, мс", ylabel="R_y, 10⁻⁴ Н", title=f"EF2: R_y суммы двух зажимов; аналитическая амплитуда {Ry_amp*1e4:.3f}·10⁻⁴ Н")
a1.set_ylim(-4.2, 5.6)
a1.legend(fontsize=8, loc="upper center", frameon=False)
style(a1)
a2.plot(t2 * 1e3, R2["x"], color=C["R_x"], lw=1.4, label="R_x суммарная (лев. + прав.)")
a2.axhline(+N, color=C["ref"], ls="--", lw=1.0)
a2.axhline(-N, color=C["ref"], ls="--", lw=1.0, label=f"±N = ±{N:.4f} Н (аналитика по зажимам)")
a2.set(xlabel="t, мс", ylabel="R_x, Н", title="EF2: R_x — отчёт суммирует оба зажима",
       ylim=(-1.2 * N, 2.0 * N))
a2.legend(fontsize=8, loc="upper center", frameon=False)
style(a2)
fig.savefig(os.path.join(ROOT, "ef2_reaction.svg"), format="svg")
plt.close(fig)

txt = "\n".join(out) + "\n"
open(os.path.join(DATA, "metrics.txt"), "w").write(txt)
print(txt)
