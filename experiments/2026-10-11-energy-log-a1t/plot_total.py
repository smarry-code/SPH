"""Суммарная механическая энергия E_mech во времени, разбитая по носителям: smoke-выстрел Guo 2019 A_1T.

Данные: energy_log.csv (include/energy_log.hpp), см. plot_energy.py.
Запуск: python plot_total.py [energy_log.csv] [energy_total.png]

Слева -- E_mech в Дж: высота стопки равна E_mech, слои -- носители.
Справа -- те же слои в долях E_mech (100 %), чтобы видеть перераспределение
и тогда, когда сама E_mech падает.
Носители: кинетическая Solid; внутренняя Solid за вычетом упругой энергии
сдвига; упругая энергия сдвига (часть Σm·u); ткань целиком (кинетическая
бусин + связи + Герц нить-нить и нить-Solid -- по отдельности она тонка).
"""
import csv
import sys

import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "energy_log.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "energy_total.png"

with open(src, newline="") as f:
    rows = list(csv.DictReader(f))
col = {k: [float(r[k]) for r in rows] for k in rows[0]}
t = [x * 1e6 for x in col["time"]]  # мкс
n = len(t)

# Эталонная палитра dataviz, слоты 1-4 в фиксированном порядке (соседние пары проверены).
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

fabric = [col["KE_fiber"][i] + col["U_bond"][i] + col["U_hertz_ff"][i] + col["U_hertz_fs"][i] for i in range(n)]
u_rest = [col["U_int"][i] - col["U_shear_el"][i] for i in range(n)]
layers = [
    ("кинетическая шара и пластины", col["KE_solid"], BLUE),
    ("внутренняя Solid без упругой сдвига", u_rest, ORANGE),
    ("упругая энергия сдвига Solid", col["U_shear_el"], YELLOW),
    ("ткань: кинетическая + связи + Герц", fabric, AQUA),
]
total = col["E_mech"]
share = [[100.0 * y[i] / total[i] for i in range(n)] for _, y, _ in layers]

plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.6), facecolor=SURF)
for ax in (a1, a2):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("время, мкс")
    ax.set_xlim(t[0], t[-1])

a1.stackplot(t, [y for _, y, _ in layers], colors=[c for _, _, c in layers], edgecolor=SURF, linewidth=1.0,
             labels=[f"{name}: {y[-1]:.3f} Дж ({100 * y[-1] / total[-1]:.0f} %)" for name, y, _ in layers])
a1.plot(t, total, color=INK, linewidth=1.5)
a1.annotate(f"E_mech = {total[-1]:.3f} Дж", (t[-1], total[-1]), xytext=(-20, 45), textcoords="offset points",
            ha="right", fontsize=9, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
a1.annotate(f"E_mech = {total[0]:.3f} Дж", (t[0], total[0]), xytext=(4, 8), textcoords="offset points",
            fontsize=9, color=INK)
a1.set_ylim(0, total[0] * 1.12)
a1.set_ylabel("энергия, Дж")
a1.set_title("Суммарная механическая энергия по носителям", loc="left", fontsize=12, color=INK)

a2.stackplot(t, share, colors=[c for _, _, c in layers], edgecolor=SURF, linewidth=1.0)
a2.set_ylim(0, 100)
a2.set_ylabel("доля E_mech, %")
a2.set_title("Те же носители в долях", loc="left", fontsize=12, color=INK)

handles, labels = a1.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.06, 0.0), ncol=2, frameon=False, fontsize=9,
           title="носитель: значение в конце (доля E_mech)", title_fontsize=9)
fig.suptitle("Guo 2019 smoke A_1T, V0 = 100 м/с: распределение механической энергии во времени",
             x=0.01, ha="left", fontsize=11, color=INK2)
fig.tight_layout(rect=(0, 0.13, 1, 0.95))
fig.savefig(out, dpi=150, facecolor=SURF)
print("ok", out)
