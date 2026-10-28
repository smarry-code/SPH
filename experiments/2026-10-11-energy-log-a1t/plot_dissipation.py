"""Суммарная диссипация во времени, разбитая по каналам: smoke-выстрел Guo 2019 A_1T.

Данные: energy_log.csv (include/energy_log.hpp), см. plot_energy.py.
Запуск: python plot_dissipation.py [energy_log.csv] [energy_dissipation.png]

Слева -- накопленные потери в Дж: высота стопки -- сумма каналов D_*, чёрная
линия -- фактическая убыль механической энергии E0 - E_mech. Разрыв между ними --
невязка баланса + S_u_clamp (ошибка схемы и источник от обрезки u >= 0).
Справа -- те же каналы в долях их суммы (100 %); до контакта потерь нет.
Цвета каналов -- как в plot_energy.py.
"""
import csv
import math
import sys

import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "energy_log.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "energy_dissipation.png"

with open(src, newline="") as f:
    rows = list(csv.DictReader(f))
col = {k: [float(r[k]) for r in rows] for k in rows[0]}
t = [x * 1e6 for x in col["time"]]  # мкс
n = len(t)

# Эталонная палитра dataviz, слоты 4-6 (соседние пары проверены); те же цвета, что в plot_energy.py.
YELLOW, MAGENTA, GREEN = "#eda100", "#e87ba4", "#008300"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

fric_brk = [col["D_friction"][i] + col["D_break"][i] for i in range(n)]
layers = [
    ("контакт шар–пластина (pinball)", col["D_contact_ss"], YELLOW),
    ("XSPH (сглаживание скоростей)", col["D_xsph"], MAGENTA),
    ("трение + разрыв связей", fric_brk, GREEN),
]
total = [sum(y[i] for _, y, _ in layers) for i in range(n)]
e0 = col["E_mech"][0]
actual = [e0 - e for e in col["E_mech"]]
# Доли -- только когда потери уже заметны (до контакта сумма -- шум округления).
eps = 1e-6 * e0
share = [[100.0 * y[i] / total[i] if total[i] > eps else math.nan for i in range(n)] for _, y, _ in layers]

def pct(f):
    return "<1 %" if 0 < f < 0.005 else f"{100 * f:.0f} %"


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
             labels=[f"{name}: {y[-1]:.3f} Дж ({pct(y[-1] / total[-1])})" for name, y, _ in layers])
a1.plot(t, actual, color=INK, linewidth=1.5, linestyle="--", label=f"фактическая убыль E0 − E_mech: {actual[-1]:.3f} Дж")
# Подписи итогов -- в пустой области до контакта, со стрелкой-выноской к концу кривой.
a1.annotate(f"сумма каналов = {total[-1]:.3f} Дж", (t[-1], total[-1]), xytext=(0.5, 0.92 * total[-1]),
            textcoords="data", fontsize=9, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
a1.text(0.5, 0.80 * total[-1], f"фактическая убыль = {actual[-1]:.3f} Дж\n"
        f"разница (невязка + S_u_clamp) = {total[-1] - actual[-1]:.3f} Дж", fontsize=9, color=INK, va="top")
a1.set_ylim(0, total[-1] * 1.15)
a1.set_ylabel("накопленная энергия, Дж")
a1.set_title("Накопленная диссипация по каналам", loc="left", fontsize=12, color=INK)

a2.stackplot(t, share, colors=[c for _, _, c in layers], edgecolor=SURF, linewidth=1.0)
a2.set_ylim(0, 100)
a2.set_ylabel("доля суммы каналов, %")
a2.set_title("Те же каналы в долях", loc="left", fontsize=12, color=INK)
t_on = next(t[i] for i in range(n) if total[i] > eps)
a2.text(t_on / 2, 50, "потерь нет\n(до контакта)", ha="center", va="center", fontsize=9, color=INK2)

handles, labels = a1.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.06, 0.0), ncol=2, frameon=False, fontsize=9,
           title="канал: значение в конце (доля суммы каналов)", title_fontsize=9)
fig.suptitle(f"Guo 2019 smoke A_1T, V0 = 100 м/с: диссипация во времени (D_break = {col['D_break'][-1]:g}: "
             "разрывов нет)", x=0.01, ha="left", fontsize=11, color=INK2)
fig.tight_layout(rect=(0, 0.13, 1, 0.95))
fig.savefig(out, dpi=150, facecolor=SURF)
print("ok", out)
