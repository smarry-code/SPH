#!/usr/bin/env python3
"""R(t) опоры и смещения против аналитики для ES1, ES2, ES3, EF3 (e2e закрепления Solid).

Читает data/<прогон>/reaction.csv (step,t,Rs_x..Rf_z; R = R_solid + R_fiber) и disp.csv,
константы берёт из tests/clamp_case.hpp (секция solid) и materials.hpp (алюминий, itype 4),
аналитику считает заново. Пишет es1_reaction.svg, es1_variants.svg, es2_cantilever.svg,
es3_pinned.svg, ef3_reaction.svg и data/metrics.txt (печатается и в stdout).
Запуск из любого места: python scripts/plot_reaction.py
Время на графиках -- середина шага t - h/2 (как в тесте): R шага n действует на [t_{n-1}, t_n].
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "clamp-solid-viz"
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")

# --- clamp_case::solid и Materials (алюминий, itype 4) ------------------------
rho = 2700.0
K = 1.97e10 * 4.2           # aa[1] * zz[1] (itype 4 - n_Me = 1)
G = 26.0e9
E = 9 * K * G / (3 * K + G)
c_L = np.sqrt((K + 4 * G / 3) / rho)
dx = 1.0e-3
strain = 1.0e-3
v0 = strain * c_L
# ES1 (smoke): L = 30 dx свободная часть, сечение 4x4
L1, B1 = 30 * dx, 4 * dx
A1 = B1 * B1
Rp = rho * c_L * v0 * A1
M_free = rho * A1 * L1
# ES2: H = 4 dx, L = 10 H
H2 = 4 * dx
L2 = 10 * H2
T_EB = 1.0 / ((1.875104 ** 2 / (2 * np.pi)) * (H2 / L2 ** 2) * np.sqrt(E / (12 * rho)))
amp2 = 0.01 * L2
# ES3: H = 5 dx, L = 10 H, L_eff = (nx-1) dx
H3 = 5 * dx
nx3 = 10 * 5
Leff3 = (nx3 - 1) * dx
f3 = (np.pi / 2) * (H3 / Leff3 ** 2) * np.sqrt(E / (12 * rho))
T3 = 1 / f3
amp3 = 0.01 * Leff3

C = dict(R="#1f77b4", ref="#d95f02", ana="#444444", c2="#2a9d55", c3="#8a4fbf")


def load(run, col="Rs_x"):
    a = np.genfromtxt(os.path.join(DATA, run, "reaction.csv"), delimiter=",", names=True)
    t = a["t"]
    tm = t - 0.5 * np.diff(np.r_[0.0, t])
    R = {ax: a["Rs_" + ax] + a["Rf_" + ax] for ax in "xyz"}
    Rs = {ax: a["Rs_" + ax] for ax in "xyz"}
    return t, tm, R, Rs


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.5)


def first_cross(t, y, sign0, thr):
    """Первая смена знака y (с гистерезисом thr), линейная интерполяция."""
    state = 0
    for i in range(len(y)):
        if state == 0 and sign0 * y[i] > thr:
            state = 1
        elif state == 1 and sign0 * y[i] < -thr:
            j = i - 1
            while j > 0 and np.sign(y[j]) == np.sign(y[i]):
                j -= 1
            return t[j] - y[j] * (t[j + 1] - t[j]) / (y[j + 1] - y[j])
    return np.nan


def crossings(t, y, thr):
    out, state = [], 0
    for i in range(len(y)):
        s = 1 if y[i] > thr else (-1 if y[i] < -thr else 0)
        if s != 0 and s != state:
            if state != 0:
                j = i - 1
                while j > 0 and np.sign(y[j]) == np.sign(y[i]):
                    j -= 1
                out.append(t[j] - y[j] * (t[j + 1] - t[j]) / (y[j + 1] - y[j]))
            state = s
    return out


out = []
out.append(f"Константы: rho = {rho:g}, K = {K:.4e}, G = {G:.4e}, E = {E:.4e} Па, c_L = {c_L:.4f} м/с, v0 = {v0:.4f} м/с")

# ---------------- ES1 ----------------
t1, tm1, R1, _ = load("ES1_pc_ES1")
# Полупериод зеркального эталона (тот же SPH без заделки) тест считает сам и печатает в лог
# («эталон: первая смена знака на ... с»); кривая эталона в reaction.csv не пишется.
import re
def ref_halfperiod(log):
    with open(os.path.join(DATA, "logs", log), encoding="utf-8") as f:
        return float(re.search(r"эталон: первая смена знака на ([0-9.e+-]+) с", f.read()).group(1))
Tc = ref_halfperiod("ES1.log")
c_eff = 2 * L1 / Tc
Tc_cl = first_cross(tm1, R1["x"], +1, 0.1 * Rp)
Rp_eff = rho * c_eff * v0 * A1
out.append(f"ES1: A = {A1:.4e} м2, M_free = {M_free:.4e} кг, rho c_L v0 A = {Rp:.4f} Н, 2L/c_L = {2*L1/c_L:.4e} с")
out.append(f"ES1: полупериод эталона (из лога теста) = {Tc:.4e} с -> c_eff = {c_eff:.2f} м/с = {c_eff/c_L:.4f} c_L; "
           f"с заделкой {Tc_cl:.4e} с (+{(Tc_cl/Tc-1)*100:.2f} % от эталона); rho c_eff v0 A = {Rp_eff:.4f} Н (оценка)")
ES1_pl = {}
for name, tt, RR, e in (("ES1", tm1, R1["x"], 1), ):
    for k in (0, 1):
        m = (tt >= (k + 0.2) * Tc) & (tt <= (k + 0.8) * Tc)
        pl = np.average(RR[m], weights=np.gradient(tt)[m])
        ES1_pl[k] = pl
        out.append(f"ES1: плато полупериода {k+1} (окно 0.2..0.8 по c_eff): {pl:.4f} Н = {pl/(e*(-1)**k*Rp)*100-100:+.2f} % от ±rho c_L v0 A, "
                   f"{pl/(e*(-1)**k*Rp_eff)*100-100:+.2f} % от ±rho c_eff v0 A")
for k in (0, 1):
    m = (tm1 >= k * Tc) & (tm1 <= (k + 1) * Tc)
    imp = np.sum(np.diff(np.r_[0.0, t1])[m] * R1["x"][m])
    out.append(f"ES1: ∫R_x dt на [{k}, {k+1}]·T_c = {imp:.4e} Н·с против ±2 M v0 = {2*M_free*v0:.4e} ({imp/((-1)**k*2*M_free*v0)*100-100:+.2f} %)")
out.append(f"ES1: max|R_y| = {np.abs(R1['y']).max():.2e}, max|R_z| = {np.abs(R1['z']).max():.2e} Н; max|R_x| = {np.abs(R1['x']).max():.2f} Н")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.1), constrained_layout=True)
for ax, pl_lv, T_half, lab in ((a1, Rp, 2 * L1 / c_L, "c_L"), (a2, Rp_eff, Tc, "c_eff")):
    ax.plot(tm1 * 1e6, R1["x"], color=C["R"], lw=1.6, label="R_x, заделка (PC)")
    ax.axhline(pl_lv, color=C["ana"], ls=":", lw=1.1)
    ax.axhline(-pl_lv, color=C["ana"], ls=":", lw=1.1, label=f"±ρ₀·{lab}·v₀·A = ±{pl_lv:.0f} Н")
    for k in range(1, 4):
        ax.axvline(k * T_half * 1e6, color=C["ana"], ls="-.", lw=0.8, alpha=0.8,
                   label=f"2kL/{lab}, k = 1,2,3: {T_half*1e6:.2f}, {2*T_half*1e6:.2f}, {3*T_half*1e6:.2f} мкс" if k == 1 else None)
    if lab == "c_eff":
        ax.axvline(Tc * 1e6, color=C["ref"], ls="--", lw=1.0, label="полупериод эталона из лога теста (= 2L/c_eff)")
    ax.set_xlim(0, tm1.max() * 1e6)
    ax.set(xlabel="t, мкс", ylabel="R_x, Н")
    style(ax)
    ax.legend(fontsize=7.5, loc="lower left", frameon=False)
a1.set_title("ES1: плато и времена по c_L (аналитика постановки)", fontsize=10)
a2.set_title("ES1: те же по c_eff = 2L/T_c, T_c — полупериод эталона", fontsize=10)
fig.savefig(os.path.join(ROOT, "es1_reaction.svg"))
plt.close(fig)

# варианты: сжатие, растяжение, Эйлер
fig, ax = plt.subplots(figsize=(7.6, 3.9), constrained_layout=True)
for run, lab, col, ls in (("ES1_pc_ES1", "PC, сжатие первым", C["R"], "-"),
                          ("ES1_pc_pull_ES1", "PC, --variant=pull (растяжение первым)", C["c2"], "-"),
                          ("ES1_euler_ES1", "Эйлер, сжатие первым", C["c3"], "--")):
    t_, tm_, R_, _ = load(run)
    ax.plot(tm_ * 1e6, R_["x"] / Rp, color=col, lw=1.3, ls=ls, label=lab)
    tc_ = first_cross(tm_, R_["x"], 1 if "pull" not in run else -1, 0.1 * Rp)
    out.append(f"ES1 [{run}]: 1-я смена знака R_x при t = {tc_:.4e} с ({tc_/(2*L1/c_L)*100-100:+.1f} % от 2L/c_L); "
               f"max|R_x|/Rp = {np.abs(R_['x']).max()/Rp:.3f}")
ax.axhline(1, color=C["ana"], ls=":", lw=1); ax.axhline(-1, color=C["ana"], ls=":", lw=1, label="±ρ₀ c_L v₀ A")
for k in (1, 2):
    ax.axvline(k * 2 * L1 / c_L * 1e6, color=C["ana"], ls="-.", lw=0.8)
ax.set(xlabel="t, мкс", ylabel="R_x / (ρ₀ c_L v₀ A)", title="ES1: три прогона, R_x в долях аналитического плато по c_L", ylim=(-1.15, 1.15))
style(ax)
ax.legend(fontsize=8, loc="lower left", frameon=False)
fig.savefig(os.path.join(ROOT, "es1_variants.svg"))
plt.close(fig)

# ---------------- ES2 ----------------
t2, tm2, R2, _ = load("ES2_pc_ES2")
d2 = np.genfromtxt(os.path.join(DATA, "ES2_pc_ES2", "disp.csv"), delimiter=",", names=True)
tz2 = crossings(d2["t"], d2["z_mean"], 0.02 * amp2)
out.append(f"ES2: T_EB = {T_EB:.4e} с, амплитуда конца {amp2:.2e} м; нули смещения конца по 80 кадрам: {', '.join(f'{x:.4e}' for x in tz2)}")
if len(tz2) >= 2:
    T2m = 2 * (tz2[-1] - tz2[0]) / (len(tz2) - 1)
    out.append(f"ES2: период по кадрам = {T2m:.4e} с ({(T2m/T_EB-1)*100:+.2f} % от T_EB); в логе теста по шагам 4.5150e-04 с (-6.76 %)")
out.append(f"ES2: макс. |z| конца по кадрам: {np.abs(d2['z_mean']).max():.4e} м против A = {amp2:.4e} м (формально V/omega_EB; {np.abs(d2['z_mean']).max()/amp2*100-100:+.1f} %, оценка по 80 точкам)")
out.append(f"ES2: max|R_z| = {np.abs(R2['z']).max():.4f} Н; max|R_x| = {np.abs(R2['x']).max():.4f}, max|R_y| = {np.abs(R2['y']).max():.2e} Н")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 3.9), constrained_layout=True)
tt = np.linspace(0, d2["t"].max(), 600)
a1.plot(d2["t"] * 1e6, d2["z_mean"] * 1e3, "o", color=C["R"], ms=3, label="SPH: средний z концевого слоя (80 кадров)")
a1.plot(tt * 1e6, amp2 * np.sin(2 * np.pi * tt / T_EB) * 1e3, color=C["ana"], ls="--", lw=1.2, label=f"A·sin(2πt/T_EB), T_EB = {T_EB*1e6:.1f} мкс")
for k in range(1, 3):
    a1.axvline(k * T_EB / 2 * 1e6, color=C["ana"], ls=":", lw=0.8)
a1.set(xlabel="t, мкс", ylabel="смещение конца z, мм", title="ES2: консоль, смещение конца против Эйлера–Бернулли", ylim=(-0.5, 0.5) if False else None)
style(a1); a1.legend(fontsize=8, frameon=False, loc="lower left")
a2.plot(tm2 * 1e6, R2["z"], color=C["R"], lw=1.2, label="R_z (по шагам)")
for k in range(0, 3):
    a2.axvline((k + 0.25) * T_EB * 1e6 if False else (k * T_EB / 2) * 1e6, color=C["ana"], ls=":", lw=0.8,
               label=f"k·T_EB/2, T_EB = {T_EB*1e6:.1f} мкс" if k == 1 else None)
a2.axhline(0, color="#999999", lw=0.6)
a2.set(xlabel="t, мкс", ylabel="R_z, Н", title="ES2: реакция заделки")
style(a2); a2.legend(fontsize=8, frameon=False, loc="upper right")
fig.savefig(os.path.join(ROOT, "es2_cantilever.svg"))
plt.close(fig)

# ---------------- ES3 ----------------
t3, tm3, R3, _ = load("ES3_pc_ES3")
d3 = np.genfromtxt(os.path.join(DATA, "ES3_pc_ES3", "disp.csv"), delimiter=",", names=True)
tz3 = crossings(d3["t"], d3["z_mean"], 0.02 * amp3)
out.append(f"ES3: T_an = {T3:.4e} с (L_eff = {Leff3:.3e} м), амплитуда середины {amp3:.2e} м; нули смещения по 80 кадрам: {', '.join(f'{x:.4e}' for x in tz3)}")
if len(tz3) >= 2:
    T3m = 2 * (tz3[-1] - tz3[0]) / (len(tz3) - 1)
    out.append(f"ES3: период по кадрам = {T3m:.4e} с ({(T3m/T3-1)*100:+.2f} % от T_an); в логе теста по шагам 2.5097e-04 с (+21.2 %)")
out.append(f"ES3: макс. |z| середины по кадрам: {np.abs(d3['z_mean']).max():.4e} м против A = {amp3:.4e} м ({np.abs(d3['z_mean']).max()/amp3*100-100:+.1f} %, оценка по 80 точкам)")
out.append(f"ES3: max|R_z| = {np.abs(R3['z']).max():.2f} Н; max|R_x| = {np.abs(R3['x']).max():.1e}, max|R_y| = {np.abs(R3['y']).max():.1e} Н")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 3.9), constrained_layout=True)
tt = np.linspace(0, d3["t"].max(), 600)
a1.plot(d3["t"] * 1e6, d3["z_mean"] * 1e3, "o", color=C["R"], ms=3, label="SPH: средний z слоёв 24–25 (80 кадров)")
a1.plot(tt * 1e6, amp3 * np.sin(2 * np.pi * tt / T3) * 1e3, color=C["ana"], ls="--", lw=1.2,
        label=f"A·sin(2πt/T_an), T_an = {T3*1e6:.1f} мкс (шарнирная балка)")
for k in range(1, 3):
    a1.axvline(k * T3 / 2 * 1e6, color=C["ana"], ls=":", lw=0.8)
a1.set(xlabel="t, мкс", ylabel="смещение середины z, мм", title="ES3: балка на шарнирах, смещение середины против аналитики")
style(a1); a1.legend(fontsize=8, frameon=False, loc="lower left")
a2.plot(tm3 * 1e6, R3["z"], color=C["R"], lw=1.2, label="R_z опор (по шагам)")
a2.plot(tm3 * 1e6, R3["x"], color=C["ref"], lw=1.0, label="R_x ≡ 0")
for k in range(0, 3):
    a2.axvline(k * T3 / 2 * 1e6, color=C["ana"], ls=":", lw=0.8, label=f"k·T_an/2, T_an = {T3*1e6:.1f} мкс" if k == 1 else None)
a2.set(xlabel="t, мкс", ylabel="R, Н", title="ES3: реакция опор")
style(a2); a2.legend(fontsize=8, frameon=False, loc="upper right")
fig.savefig(os.path.join(ROOT, "es3_pinned.svg"))
plt.close(fig)

# ---------------- EF3 ----------------
fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.7), constrained_layout=True)
for ax, run, ttl in zip(axs, ("EF3_pc_EF3", "EF3_pc_EF3_2", "EF3_pc_EF3_3"),
                        ("Б: падение без трения", "В(а): mu = 0, v_t0 = 2 м/с", "В(б): mu = 0.2, v_t0 = 2 м/с")):
    t_, tm_, R_, Rs_ = load(run)
    dts = np.diff(np.r_[0.0, t_])
    Jz, Jx = np.sum(dts * R_["z"]), np.sum(dts * R_["x"])
    out.append(f"EF3 [{ttl}]: ∫R_z dt = {Jz:.4e} Н·с, ∫R_x dt = {Jx:.4e} Н·с; max R_z = {R_['z'].max():.4e} Н, шагов {len(t_)}")
    ax.plot(tm_ * 1e3, R_["z"], color=C["R"], lw=1.3, label=f"R_z, ∫ = {Jz:.4e} Н·с")
    ax.plot(tm_ * 1e3, R_["x"], color=C["ref"], lw=1.1, label=f"R_x, ∫ = {Jx:.4e} Н·с")
    ax.set(xlabel="t, мс", ylabel="R, Н", title=ttl, xlim=(0, tm_.max() * 1e3))
    style(ax); ax.legend(fontsize=8, frameon=False, loc="upper right")
fig.savefig(os.path.join(ROOT, "ef3_reaction.svg"))
plt.close(fig)

txt = "\n".join(out) + "\n"
print(txt)
with open(os.path.join(DATA, "metrics.txt"), "w") as f:
    f.write(txt)
