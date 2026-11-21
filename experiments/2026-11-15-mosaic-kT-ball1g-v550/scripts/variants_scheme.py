#!/usr/bin/env python3
"""Схема (не расчёт) сшивки вариантов A и B на 4 слоях, разрез поперёк швов -> fig/variants_scheme.png.

A: швы всех слоёв в одних и тех же местах (одна раскройка); в точке шва -- стежок в слое через шов и цепочка
   сквозных стежков 1-2, 2-3, 3-4.
B: у каждого слоя свои швы; в точке шва ЛЮБОГО слоя ставится та же цепочка сквозных стежков через весь пакет,
   в т.ч. в местах, где остальные слои целые.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fig", "variants_scheme.png")
NL, L = 4, 20.0
COL = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
GAP = 0.35


def draw(ax, seams, title):
    """seams[k] -- координаты швов слоя k (слой 0 -- верхний)."""
    allpts = sorted({x for s in seams for x in s})
    for k in range(NL):
        z = -k
        edges = [0.0] + sorted(seams[k]) + [L]
        for i in range(len(edges) - 1):
            x0 = edges[i] + (GAP if i > 0 else 0)
            x1 = edges[i + 1] - (GAP if i < len(edges) - 2 else 0)
            c = COL[(i + 2 * k) % len(COL)] if len(set(map(tuple, seams))) > 1 else COL[i % len(COL)]
            ax.fill_between([x0, x1], z - 0.17, z + 0.17, color=c, alpha=0.75, lw=0)
        for x in seams[k]:  # стежок в слое через шов
            ax.plot([x - GAP - 0.25, x + GAP + 0.25], [z, z], color="#009E73", lw=3.2, solid_capstyle="round")
            ax.plot([x - GAP - 0.25, x + GAP + 0.25], [z, z], color="white", lw=1.0, ls=(0, (1, 1.5)))
        ax.text(-0.6, z, "слой %d" % (k + 1), ha="right", va="center", fontsize=10)
    for x in allpts:  # сквозные стежки: цепочка между соседними слоями
        for k in range(NL - 1):
            ax.plot([x, x], [-k - 0.17, -k - 1 + 0.17], color="k", lw=1.6)
            ax.plot([x], [-k - 0.5], marker="o", ms=3, color="k")
    ax.set_xlim(-3.5, L + 0.5)
    ax.set_ylim(-NL + 0.3, 0.9)
    ax.axis("off")
    ax.set_title(title, fontsize=11, loc="left")


fig, (a, b) = plt.subplots(2, 1, figsize=(11, 7.2))
draw(a, [[6.0, 13.0]] * NL,
     "A — одна раскройка на пакет: швы всех слоёв друг над другом, сквозная цепочка только на этом шве")
draw(b, [[6.0, 13.0], [3.5, 10.0, 16.5], [8.0, 15.0], [5.0, 11.5, 18.0]],
     "B — у каждого слоя своя раскройка: сквозная цепочка ставится на шве ЛЮБОГО слоя, через весь пакет")
a.plot([], [], color="#009E73", lw=3, label="стежок в слое: соединяет два соседних кусочка через шов")
a.plot([], [], color="k", lw=1.6, marker="o", ms=3, label="сквозной стежок: соединяет соседние слои (1–2, 2–3, …)")
a.legend(loc="upper center", bbox_to_anchor=(0.5, 1.32), ncol=1, fontsize=9, frameon=False)
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=130)
print(OUT)
