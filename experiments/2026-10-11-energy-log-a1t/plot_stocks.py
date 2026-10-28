"""Запасы механической энергии во времени (без каналов потерь): smoke-выстрел Guo 2019 A_1T.

Данные: energy_log.csv (include/energy_log.hpp), см. plot_energy.py.
Запуск: python plot_stocks.py [energy_log.csv] [energy_stocks.png]

Слева -- Solid (шар и пластина) и полная механическая энергия E_mech, Дж.
Справа -- ткань, мДж: её запасы на два-три порядка меньше, в общей шкале их не видно.
"""
import csv
import sys

import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "energy_log.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "energy_stocks.png"

with open(src, newline="") as f:
    rows = list(csv.DictReader(f))
col = {k: [float(r[k]) for r in rows] for k in rows[0]}
t = [x * 1e6 for x in col["time"]]  # мкс

# Эталонная палитра dataviz, слоты в фиксированном порядке.
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2), facecolor=SURF)
for ax in (a1, a2):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("время, мкс")
    ax.set_xlim(t[0], t[-1])


# --- Solid и полная ---
solid = [
    ("E_mech (полная механическая)", col["E_mech"], INK, "-"),
    ("кинетическая шара и пластины", col["KE_solid"], BLUE, "-"),
    ("внутренняя Solid, Σm·u", col["U_int"], ORANGE, "-"),
    ("в т.ч. упругая сдвига", col["U_shear_el"], ORANGE, ":"),
]
for name, y, c, ls in solid:
    a1.plot(t, y, color=c, linewidth=2, linestyle=ls, label=f"{name}: {y[-1]:.3f}")
a1.set_ylim(0, max(col["E_mech"]) * 1.08)
a1.set_ylabel("энергия, Дж")
a1.set_title("Шар и пластина (Solid)", loc="left", fontsize=12, color=INK)
a1.legend(loc="center left", frameon=False, fontsize=9)

# --- ткань, мДж ---
mj = lambda v: [x * 1e3 for x in v]
fabric = [
    ("кинетическая бусин", mj(col["KE_fiber"]), AQUA),
    ("упругая связей", mj(col["U_bond"]), YELLOW),
    ("Герц нить–нить", mj(col["U_hertz_ff"]), MAGENTA),
    ("Герц нить–Solid", mj(col["U_hertz_fs"]), GREEN),
]
for name, y, c in fabric:
    a2.plot(t, y, color=c, linewidth=2, label=f"{name}: {y[-1]:.2f}")
a2.set_ylim(bottom=0)
a2.set_ylabel("энергия, мДж")
a2.set_title("Ткань", loc="left", fontsize=12, color=INK)
a2.legend(loc="upper left", frameon=False, fontsize=9)

fig.suptitle("Guo 2019 smoke A_1T, V0 = 100 м/с: запасы механической энергии (без потерь); "
             "в подписях -- значения в конце", x=0.01, ha="left", fontsize=11, color=INK2)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(out, dpi=150, facecolor=SURF)
print("ok", out)
