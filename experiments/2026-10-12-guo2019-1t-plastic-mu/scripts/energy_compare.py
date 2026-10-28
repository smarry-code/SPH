"""Сравнение энергобаланса: упругая 73 (Guo2T), пластичные 66 mu0.3 и mu0.5; запасы пластичной 73 по monitor.csv.
Запуск: python energy_compare.py   (пути заданы ниже). Пишет energy/table.md, energy/friction_compare.png,
energy/mu0.3_v73_plastic_monitor/{stocks_v73_compare,etot_v73_compare}.png
"""
import csv
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EN = os.path.join(ROOT, "energy")
EL = r"C:\dip\sph\SPH\.claude\worktrees\Guo2T\experiments\2026-10-12-energy-log-1t-v73-full"
RUNS = [("упругая 73, μ=0.3", os.path.join(EL, "energy_log.csv"), os.path.join(EL, "shot.txt"), "#2a78d6"),
        ("пластичная 66, μ=0.3", os.path.join(EN, "mu0.3_v66", "energy_log.csv"),
         os.path.join(EN, "mu0.3_v66", "shot.txt"), "#eb6834"),
        ("пластичная 66, μ=0.5", os.path.join(EN, "mu0.5_v66", "energy_log.csv"),
         os.path.join(EN, "mu0.5_v66", "shot.txt"), "#008300")]
INK2, GRID, SURF = "#52514e", "#e4e3df", "#fcfcfb"


def load(p):
    rows = list(csv.DictReader(open(p, newline="")))
    c = {k: [float(r[k]) for r in rows] for k in rows[0]}
    for k in ("W_p_fiber", "W_v_fiber"):
        c.setdefault(k, [0.0] * len(rows))
    return c


def shot(p):
    return dict(re.findall(r"(\w+)=(\S+)", open(p).read()))


def at(c, k, t):
    T = c["time"]
    if t >= T[-1]:
        return c[k][-1]
    for i in range(1, len(T)):
        if T[i] >= t:
            f = (t - T[i - 1]) / (T[i] - T[i - 1])
            return c[k][i - 1] + f * (c[k][i] - c[k][i - 1])


KEYS = ("time", "KE_solid", "KE_fiber", "U_bond", "U_hertz_ff", "U_hertz_fs", "U_int", "D_friction", "D_break", "D_xsph",
        "D_contact_ss", "W_p_fiber", "W_v_fiber", "residual", "S_u_clamp", "E_mech")
data = [(n, load(pc), shot(ps), col) for n, pc, ps, col in RUNS]
rows_def = [("t, мкс", None), ("KE шара", "KE_solid"), ("поглощено E0 - KE шара", "ABS"), ("D_friction (трение)", "D_friction"),
            ("D_break", "D_break"), ("W_p_fiber", "W_p_fiber"), ("W_v_fiber", "W_v_fiber"), ("D_xsph", "D_xsph"),
            ("D_contact_ss", "D_contact_ss"), ("U_bond", "U_bond"), ("KE нитей", "KE_fiber"), ("невязка", "residual")]
fmt = lambda v: "%.4g" % v
out = []
for label, which in (("КОНЕЦ РАСЧЁТА", "end"), ("ПЕРВЫЙ РАЗРЫВ СВЯЗЕЙ", "break")):
    out.append("### %s\n" % label)
    out.append("| величина (Дж и % от E0) | " + " | ".join(n for n, *_ in data) + " |")
    out.append("|---|" + "---|" * len(data))
    sts = []
    for n, c, sh, col in data:
        t = c["time"][-1] if which == "end" else float(sh["t_first_break"])
        sts.append({k: at(c, k, t) for k in KEYS})
    e0s = [c["E_mech"][0] for _, c, _, _ in data]
    out.append("| E0, Дж | " + " | ".join(fmt(e) for e in e0s) + " |")
    for rl, k in rows_def:
        cells = []
        for s, e0 in zip(sts, e0s):
            if k is None:
                cells.append("%.1f" % (s["time"] * 1e6))
                continue
            v = (e0 - s["KE_solid"]) if k == "ABS" else s[k]
            cells.append("%s (%.2f %%)" % (fmt(v), 100 * v / e0))
        out.append("| %s | " % rl + " | ".join(cells) + " |")
    out.append("| доля трения D_friction / (E0 - KE шара) | " + " | ".join(
        "%.1f %%" % (100 * s["D_friction"] / (e0 - s["KE_solid"])) for s, e0 in zip(sts, e0s)) + " |")
    out.append("")
out.append("Из shot.txt:\n")
for n, c, sh, _ in data:
    out.append("- %s: outcome=%s V_R=%s e=%s t_end=%s t_first_break=%s n_broken_v=%s n_broken_h=%s" % (
        n, sh["outcome"], sh["V_R"], sh["e"], sh["t_end"], sh["t_first_break"], sh["n_broken_v"], sh["n_broken_h"]))
open(os.path.join(EN, "table.md"), "w", encoding="utf-8").write("\n".join(out))
print("\n".join(out))

# ---- friction_compare
plt.rcParams.update({"font.size": 10})
fig, ax = plt.subplots(2, 2, figsize=(13, 9), facecolor=SURF)
for a in ax.flat:
    a.set_facecolor(SURF)
    a.grid(True, color=GRID)
    a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
for n, c, sh, col in data:
    t = [x * 1e6 for x in c["time"]]
    e0 = c["E_mech"][0]
    tb = float(sh["t_first_break"]) * 1e6
    ax[0, 0].plot(t, c["D_friction"], color=col, lw=2, label=n)
    ax[0, 1].plot(t, [100 * x / e0 for x in c["D_friction"]], color=col, lw=2, label=n)
    lab = n if max(c["W_p_fiber"]) > 0 else n + " (упругая нить: W_p_fiber = 0)"
    ax[1, 0].plot(t, [x * 1e3 for x in c["W_p_fiber"]], color=col, lw=2, label=lab)
    ax[1, 1].plot([x - tb for x in t], [100 * x / e0 for x in c["D_friction"]], color=col, lw=2, label=n)
    for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
        a.axvline(tb, color=col, ls=":", lw=1.3)
ax[0, 0].set_title("D_friction(t), Дж", loc="left")
ax[0, 0].set_ylabel("Дж")
ax[0, 1].set_title("D_friction(t) / E0, %  (E0 = 18.65 / 15.25 / 15.25 Дж)", loc="left")
ax[0, 1].set_ylabel("% от E0")
ax[1, 0].set_title("W_p_fiber(t): пластическая работа нити, мДж", loc="left")
ax[1, 0].set_ylabel("мДж")
ax[1, 1].set_title("D_friction / E0 от момента первого разрыва", loc="left")
ax[1, 1].set_ylabel("% от E0")
ax[1, 1].axvline(0, color=INK2, ls=":")
for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
    a.set_xlabel("время, мкс (пунктир: первый разрыв связей)")
ax[1, 1].set_xlabel("t - t первого разрыва, мкс")
ax[0, 0].legend(frameon=False, fontsize=9, loc="upper left")
ax[0, 1].legend(frameon=False, fontsize=9, loc="lower right")
ax[1, 0].legend(frameon=False, fontsize=9, loc="center left")
ax[1, 1].legend(frameon=False, fontsize=9, loc="lower right")
fig.suptitle("Трение: упругая нить 73 м/с (код Guo2T) и пластичная 66 м/с. Условия не совпадают (V0, нить, t_end), см. README",
             x=0.01, ha="left", color=INK2)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(EN, "friction_compare.png"), dpi=150, facecolor=SURF)
plt.close(fig)

# ---- запасы 73: пластичная (monitor) и упругая (лог)
m = {k: [] for k in ("t", "KE_s", "KE_p", "KE_f", "U_b", "U_int", "E_tot", "n_broken_v")}
for r in csv.DictReader(open(os.path.join(EN, "mu0.3_v73_plastic_monitor", "monitor.csv"))):
    m["t"].append(float(r["t"]) * 1e6)
    for k in m:
        if k != "t":
            m[k].append(float(r[k]))
tbm = next(t for t, n in zip(m["t"], m["n_broken_v"]) if n > 0)
el = data[0][1]
te = [x * 1e6 for x in el["time"]]
tbe = float(data[0][2]["t_first_break"]) * 1e6
fig, ax = plt.subplots(1, 3, figsize=(16, 5.6), facecolor=SURF)
pairs = [("KE шара (+пластины)", m["KE_s"], el["KE_solid"]), ("KE нитей (бусин)", m["KE_f"], el["KE_fiber"]),
         ("упругая энергия связей U_bond", m["U_b"], el["U_bond"])]
for a, (nm, yp, ye) in zip(ax, pairs):
    a.set_facecolor(SURF)
    a.plot(te, ye, color="#2a78d6", lw=2, label="упругая 73 (лог энергии, Guo2T)")
    a.plot(m["t"], yp, color="#e04a3f", lw=2, label="пластичная 73 (monitor.csv)")
    a.axvline(tbe, color="#2a78d6", ls=":")
    a.axvline(tbm, color="#e04a3f", ls=":")
    a.set_title(nm, loc="left")
    a.set_ylabel("Дж")
    a.set_xlabel("время, мкс (пунктир: первый разрыв)")
    a.grid(True, color=GRID)
    a.legend(frameon=False, fontsize=8.5)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
fig.suptitle("Запасы, V0 = 73 м/с, μ = 0.3: пластичная нить (только monitor) и упругая. "
             "Каналов трения для пластичной 73 нет (лога энергии не было)", x=0.01, ha="left", color=INK2, fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(EN, "mu0.3_v73_plastic_monitor", "stocks_v73_compare.png"), dpi=150, facecolor=SURF)
plt.close(fig)
fig, ax = plt.subplots(figsize=(9, 5.4), facecolor=SURF)
ax.set_facecolor(SURF)
ax.plot(te, el["E_mech"], color="#2a78d6", lw=2, label="упругая 73: E_mech (лог)")
ax.plot(m["t"], m["E_tot"], color="#e04a3f", lw=2, label="пластичная 73: E_tot (monitor)")
ax.axvline(tbe, color="#2a78d6", ls=":")
ax.axvline(tbm, color="#e04a3f", ls=":")
ax.set_xlabel("время, мкс (пунктир: первый разрыв)")
ax.set_ylabel("Дж")
ax.grid(True, color=GRID)
ax.set_title("Полный запас, 73 м/с, μ = 0.3. E_tot монитора и E_mech лога определены по-разному\n"
             "(E_mech включает Герц); каналов трения для пластичной нет", loc="left", fontsize=10)
ax.legend(frameon=False)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(EN, "mu0.3_v73_plastic_monitor", "etot_v73_compare.png"), dpi=150, facecolor=SURF)
print("monitor first break t=%.1f us; elastic %.1f; plastic end t=%.1f; elastic end %.1f" % (tbm, tbe, m["t"][-1], te[-1]))
j = min(range(len(m["t"])), key=lambda j: abs(m["t"][j] - tbm))
print("plastic73 at first break KE_s,KE_f,U_b,E_tot", m["KE_s"][j], m["KE_f"][j], m["U_b"][j], m["E_tot"][j])
print("plastic73 end KE_s,KE_f,U_b,E_tot", m["KE_s"][-1], m["KE_f"][-1], m["U_b"][-1], m["E_tot"][-1], m["n_broken_v"][-1])
