"""Энергия по каналам во времени: smoke-выстрел Guo 2019 A_1T, V0 = 100 м/с.

Данные: energy_log.csv (include/energy_log.hpp), получен командой
    guo2019_shot --config=A_1T --level=smoke --v0=100 --t-max=2e-5 --energy-log=10
из сборки с -DSPH_ENERGY_LOG=ON.

Запуск: python plot_energy.py [energy_log.csv] [energy_channels.png]

Слева -- бюджет: куда ушла начальная энергия E0 (сложенные запасы и
накопленные потери; разрыв между верхом стопки и E0 -- невязка + S_u_clamp).
Справа -- накопленные каналы потерь по отдельности, в той же шкале Дж.
"""
import csv
import sys

import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "energy_log.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "energy_channels.png"

with open(src, newline="") as f:
    rows = list(csv.DictReader(f))
col = {k: [float(r[k]) for r in rows] for k in rows[0]}
t = [x * 1e6 for x in col["time"]]  # мкс
e0 = col["E_mech"][0]

# Эталонная палитра dataviz, слоты 1-7 в фиксированном порядке; цвет следует за сущностью в обоих графиках.
C = {
    "ke_s": "#2a78d6",   # 1 синий
    "u_int": "#eb6834",  # 2 оранжевый
    "fabric": "#1baf7a", # 3 бирюзовый
    "css": "#eda100",    # 4 жёлтый
    "xsph": "#e87ba4",   # 5 пурпурный
    "fric": "#008300",   # 6 зелёный
    "ucl": "#6250d6",    # 7 фиолетовый
}
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

fabric = [a + b + c + d for a, b, c, d in zip(col["KE_fiber"], col["U_bond"], col["U_hertz_ff"], col["U_hertz_fs"])]
fric_brk = [a + b for a, b in zip(col["D_friction"], col["D_break"])]
stack = [
    ("KE шара и пластины", col["KE_solid"], C["ke_s"]),
    ("Внутренняя Solid, Σm·u", col["U_int"], C["u_int"]),
    ("Ткань: KE + связи + Герц", fabric, C["fabric"]),
    ("Потери: контакт", col["D_contact_ss"], C["css"]),
    ("Потери: XSPH", col["D_xsph"], C["xsph"]),
    ("Потери: трение + разрыв", fric_brk, C["fric"]),
]

plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.8), facecolor=SURF)
for ax in (a1, a2):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("время, мкс")
    ax.set_xlim(t[0], t[-1])

# --- бюджет энергии ---
polys = a1.stackplot(t, [s[1] for s in stack], colors=[s[2] for s in stack],
                     labels=[f"{s[0]}: {s[1][-1]:.3f}" for s in stack],
                     edgecolor=SURF, linewidth=1.0)
a1.axhline(e0, color=INK2, linewidth=1.2, linestyle="--")
a1.text(t[0] + 0.3, e0 * 1.012, f"E0 = {e0:.3f} Дж", color=INK2, fontsize=9, va="bottom")
a1.set_ylim(0, e0 * 1.15)
a1.set_ylabel("энергия, Дж")
a1.set_title("Куда уходит энергия удара", loc="left", fontsize=12, color=INK)
top = sum(v[-1] for _, v, _ in stack)
a1.annotate("", xy=(t[-1], top), xytext=(t[-1], e0), arrowprops=dict(arrowstyle="<->", color=INK2, lw=1))
a1.text(t[-1], top + 0.015, f"невязка + S_u_clamp = {top - e0:.3f} Дж", ha="right", va="bottom",
        fontsize=8.5, color=INK2)
a1.legend(loc="upper left", bbox_to_anchor=(0.0, -0.13), ncol=2, frameon=False, fontsize=9)

# --- накопленные каналы ---
lines = [
    ("контакт шар–пластина", col["D_contact_ss"], C["css"], "-"),
    ("XSPH", col["D_xsph"], C["xsph"], "-"),
    ("трение", col["D_friction"], C["fric"], "-"),
    ("обрезка u ≥ 0 (источник)", col["S_u_clamp"], C["ucl"], "-"),
    ("невязка баланса", col["residual"], INK2, "--"),
]
for name, y, c, ls in lines:
    a2.plot(t, y, color=c, linewidth=2, linestyle=ls, label=f"{name}: {y[-1]:.3f}")
# прямые подписи у концов линий: только для крупных каналов, мелкие -- в легенде
for name, y, c, ls in lines[:2] + lines[4:]:
    a2.annotate(f"{y[-1]:.3f}", (t[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                va="center", fontsize=9, color=INK2, annotation_clip=False)
a2.set_ylim(bottom=min(0.0, min(col["residual"])) - 0.01)
a2.set_ylabel("накопленная энергия, Дж")
a2.set_title("Накопленные каналы потерь", loc="left", fontsize=12, color=INK)
a2.legend(loc="upper left", frameon=False, fontsize=9)  # до удара (t < 6 мкс) панель пуста

d_brk = col["D_break"][-1]
fig.suptitle("Guo 2019 smoke A_1T, V0 = 100 м/с: механическая энергия по каналам "
             f"(разрывов нет: D_break = {d_brk:g})", x=0.01, ha="left", fontsize=11, color=INK2)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(out, dpi=150, facecolor=SURF)
print("записан", out)
