"""Графики лога механической энергии (include/energy_log.hpp, сборка -DSPH_ENERGY_LOG=ON).

Запуск: python plot_energy_all.py <energy_log.csv> [каталог_вывода] [--title "заголовок"]
Пять файлов в каталоге вывода: energy_channels.png, energy_dissipation.png, energy_stocks.png, energy_total.png,
energy_transfer.png.

Каналы пластичной нити W_p_fiber (пластическая работа связей) и W_v_fiber (вязкая диссипация Максвелл + Кельвин) -- два
последних столбца CSV. В старых CSV без них (до переноса лога в ветку plastic) столбцы считаются нулевыми, а нулевые на
всём интервале каналы в легенды не выводятся.
Невязка: residual = (E_mech - E_mech0) + D_friction + D_xsph + D_contact_ss + D_break - S_u_clamp + W_p_fiber + W_v_fiber.

Образец оформления -- plot_energy/plot_dissipation/plot_stocks/plot_total (Guo2T, experiments/2026-10-11-energy-log-a1t);
исходный plot_energy_all.py в Guo2T не найден, скрипт собран из них заново.
"""
import argparse
import csv
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
                                                      "#6250d6")
RED, BROWN, BREAK = "#e04a3f", "#8a5a2b", "#8a8a8a"  # BREAK -- разрыв связей (серый)
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def load(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    col = {k: [float(r[k]) for r in rows] for k in rows[0]}
    n = len(rows)
    for k in ("W_p_fiber", "W_v_fiber"):  # старые CSV: каналов нет
        col.setdefault(k, [0.0] * n)
    return col


def style():
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
                         "ytick.color": INK2, "text.color": INK})


def axes(ncols=2, h=5.6, t=None):
    fig, ax = plt.subplots(1, ncols, figsize=(13, h), facecolor=SURF)
    for a in ax:
        a.set_facecolor(SURF)
        a.grid(True, color=GRID, linewidth=0.8)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.set_xlabel("время, мкс")
        if t:
            a.set_xlim(t[0], t[-1])
    return fig, ax


def nonzero(y):
    return any(abs(v) > 0.0 for v in y)


def fmt(v):
    return "%.3g" % v if abs(v) < 0.01 else "%.3f" % v


def plot_channels(col, t, out, title):
    n = len(t)
    e0 = col["E_mech"][0]
    fabric = [col["KE_fiber"][i] + col["U_bond"][i] + col["U_hertz_ff"][i] + col["U_hertz_fs"][i] for i in range(n)]
    stack = [("KE шара и пластины", col["KE_solid"], BLUE), ("Внутренняя Solid, Σm·u", col["U_int"], ORANGE),
             ("Ткань: KE + связи + Герц", fabric, AQUA), ("Потери: контакт Solid-Solid", col["D_contact_ss"], YELLOW),
             ("Потери: XSPH", col["D_xsph"], MAGENTA), ("Потери: трение", col["D_friction"], GREEN),
             ("Потери: разрыв связей", col["D_break"], BREAK),
             ("Потери: пластика нити W_p_fiber", col["W_p_fiber"], RED),
             ("Потери: вязкость нити W_v_fiber", col["W_v_fiber"], BROWN)]
    stack = [s for s in stack if nonzero(s[1])]
    fig, (a1, a2) = axes(2, 5.8, t)
    a1.stackplot(t, [s[1] for s in stack], colors=[s[2] for s in stack],
                 labels=["%s: %s" % (s[0], fmt(s[1][-1])) for s in stack], edgecolor=SURF, linewidth=1.0)
    a1.axhline(e0, color=INK2, linewidth=1.2, linestyle="--")
    a1.text(t[0] + 0.3, e0 * 1.012, "E0 = %s Дж" % fmt(e0), color=INK2, fontsize=9, va="bottom")
    a1.set_ylim(0, e0 * 1.15)
    a1.set_ylabel("энергия, Дж")
    a1.set_title("Куда уходит энергия удара", loc="left", fontsize=12, color=INK)
    top = sum(s[1][-1] for s in stack)
    a1.text(t[-1], e0 * 1.02, "невязка + S_u_clamp = %s Дж" % fmt(top - e0), ha="right", va="bottom", fontsize=8.5,
            color=INK2)
    a1.legend(loc="upper left", bbox_to_anchor=(0.0, -0.13), ncol=2, frameon=False, fontsize=9)
    lines = [("контакт Solid-Solid", col["D_contact_ss"], YELLOW, "-"), ("XSPH", col["D_xsph"], MAGENTA, "-"),
             ("трение", col["D_friction"], GREEN, "-"), ("разрыв связей", col["D_break"], RED, ":"),
             ("пластика нити W_p_fiber", col["W_p_fiber"], RED, "-"),
             ("вязкость нити W_v_fiber", col["W_v_fiber"], BROWN, "-"),
             ("обрезка u >= 0 (источник)", col["S_u_clamp"], VIOLET, "-"),
             ("невязка баланса", col["residual"], INK2, "--")]
    for name, y, c, ls in lines:
        if nonzero(y):
            a2.plot(t, y, color=c, linewidth=2, linestyle=ls, label="%s: %s" % (name, fmt(y[-1])))
    a2.set_ylabel("накопленная энергия, Дж")
    a2.set_title("Накопленные каналы потерь", loc="left", fontsize=12, color=INK)
    a2.legend(loc="upper left", frameon=False, fontsize=9)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)


def plot_dissipation(col, t, out, title):
    n = len(t)
    e0 = col["E_mech"][0]
    layers = [("контакт Solid-Solid", col["D_contact_ss"], YELLOW), ("XSPH", col["D_xsph"], MAGENTA),
              ("трение", col["D_friction"], GREEN), ("разрыв связей", col["D_break"], BREAK),
              ("пластика нити W_p_fiber", col["W_p_fiber"], RED),
              ("вязкость нити W_v_fiber", col["W_v_fiber"], BROWN)]
    layers = [s for s in layers if nonzero(s[1])]
    total = [sum(y[i] for _, y, _ in layers) for i in range(n)]
    actual = [e0 - e for e in col["E_mech"]]
    eps = 1e-6 * e0
    share = [[100.0 * y[i] / total[i] if total[i] > eps else math.nan for i in range(n)] for _, y, _ in layers]
    fig, (a1, a2) = axes(2, 5.6, t)
    a1.stackplot(t, [y for _, y, _ in layers], colors=[c for _, _, c in layers], edgecolor=SURF, linewidth=1.0,
                 labels=["%s: %s Дж (%.0f %%)" % (nm, fmt(y[-1]), 100 * y[-1] / total[-1] if total[-1] > 0 else 0)
                         for nm, y, _ in layers])
    a1.plot(t, actual, color=INK, linewidth=1.5, linestyle="--", label="убыль E0 - E_mech: %s Дж" % fmt(actual[-1]))
    a1.set_ylim(0, max(total[-1], actual[-1], eps) * 1.15)
    a1.set_ylabel("накопленная энергия, Дж")
    a1.set_title("Накопленная диссипация по каналам", loc="left", fontsize=12, color=INK)
    a2.stackplot(t, share, colors=[c for _, _, c in layers], edgecolor=SURF, linewidth=1.0)
    a2.set_ylim(0, 100)
    a2.set_ylabel("доля суммы каналов, %")
    a2.set_title("Те же каналы в долях", loc="left", fontsize=12, color=INK)
    h, lb = a1.get_legend_handles_labels()
    fig.legend(h, lb, loc="lower left", bbox_to_anchor=(0.06, 0.0), ncol=2, frameon=False, fontsize=9,
               title="канал: значение в конце (доля суммы каналов)", title_fontsize=9)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=INK2)
    fig.tight_layout(rect=(0, 0.14, 1, 0.95))
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)


def plot_stocks(col, t, out, title):
    mj = lambda v: [x * 1e3 for x in v]
    fig, (a1, a2) = axes(2, 5.2, t)
    for nm, y, c, ls in (("E_mech (полная механическая)", col["E_mech"], INK, "-"),
                         ("кинетическая шара и пластины", col["KE_solid"], BLUE, "-"),
                         ("внутренняя Solid, Σm·u", col["U_int"], ORANGE, "-"),
                         ("в т.ч. упругая сдвига", col["U_shear_el"], ORANGE, ":")):
        a1.plot(t, y, color=c, linewidth=2, linestyle=ls, label="%s: %s" % (nm, fmt(y[-1])))
    a1.set_ylim(0, max(col["E_mech"]) * 1.08)
    a1.set_ylabel("энергия, Дж")
    a1.set_title("Шар и пластина (Solid)", loc="left", fontsize=12, color=INK)
    a1.legend(loc="center left", frameon=False, fontsize=9)
    for nm, y, c in (("кинетическая бусин", mj(col["KE_fiber"]), AQUA), ("упругая связей", mj(col["U_bond"]), YELLOW),
                     ("Герц нить-нить", mj(col["U_hertz_ff"]), MAGENTA),
                     ("Герц нить-Solid", mj(col["U_hertz_fs"]), GREEN)):
        a2.plot(t, y, color=c, linewidth=2, label="%s: %s" % (nm, fmt(y[-1])))
    for nm, y, c in (("пластика связей W_p_fiber (накоплено)", mj(col["W_p_fiber"]), RED),
                     ("вязкость связей W_v_fiber (накоплено)", mj(col["W_v_fiber"]), BROWN)):
        if nonzero(y):
            a2.plot(t, y, color=c, linewidth=2, linestyle="--", label="%s: %s" % (nm, fmt(y[-1])))
    a2.set_ylim(bottom=0)
    a2.set_ylabel("энергия, мДж")
    a2.set_title("Ткань (пунктир -- накопленные потери нити)", loc="left", fontsize=12, color=INK)
    a2.legend(loc="upper left", frameon=False, fontsize=9)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)


def plot_total(col, t, out, title):
    n = len(t)
    fabric = [col["KE_fiber"][i] + col["U_bond"][i] + col["U_hertz_ff"][i] + col["U_hertz_fs"][i] for i in range(n)]
    u_rest = [col["U_int"][i] - col["U_shear_el"][i] for i in range(n)]
    layers = [("кинетическая шара и пластины", col["KE_solid"], BLUE),
              ("внутренняя Solid без упругой сдвига", u_rest, ORANGE),
              ("упругая энергия сдвига Solid", col["U_shear_el"], YELLOW),
              ("ткань: кинетическая + связи + Герц", fabric, AQUA)]
    total = col["E_mech"]
    share = [[100.0 * y[i] / total[i] if total[i] > 0 else math.nan for i in range(n)] for _, y, _ in layers]
    fig, (a1, a2) = axes(2, 5.6, t)
    a1.stackplot(t, [y for _, y, _ in layers], colors=[c for _, _, c in layers], edgecolor=SURF, linewidth=1.0,
                 labels=["%s: %s Дж (%.0f %%)" % (nm, fmt(y[-1]), 100 * y[-1] / total[-1]) for nm, y, _ in layers])
    a1.plot(t, total, color=INK, linewidth=1.5)
    a1.set_ylim(0, total[0] * 1.12)
    a1.set_ylabel("энергия, Дж")
    a1.set_title("Суммарная механическая энергия по носителям", loc="left", fontsize=12, color=INK)
    a2.stackplot(t, share, colors=[c for _, _, c in layers], edgecolor=SURF, linewidth=1.0)
    a2.set_ylim(0, 100)
    a2.set_ylabel("доля E_mech, %")
    a2.set_title("Те же носители в долях", loc="left", fontsize=12, color=INK)
    h, lb = a1.get_legend_handles_labels()
    fig.legend(h, lb, loc="lower left", bbox_to_anchor=(0.06, 0.0), ncol=2, frameon=False, fontsize=9,
               title="носитель: значение в конце (доля E_mech)", title_fontsize=9)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=INK2)
    fig.tight_layout(rect=(0, 0.13, 1, 0.95))
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)


def plot_transfer(col, t, out, title):
    """Куда ушла кинетическая энергия Solid: E0 - KE_solid по получателям."""
    n = len(t)
    ke0 = col["KE_solid"][0]
    bond_h = [col["U_bond"][i] + col["U_hertz_ff"][i] + col["U_hertz_fs"][i] for i in range(n)]
    rec = [("кинетическая ткани", col["KE_fiber"], AQUA), ("упругая ткани: связи + Герц", bond_h, VIOLET),
           ("внутренняя Solid, Σm·u", col["U_int"], ORANGE), ("потери: трение", col["D_friction"], GREEN),
           ("потери: разрыв связей", col["D_break"], BREAK),
           ("потери: XSPH", col["D_xsph"], MAGENTA), ("потери: контакт Solid-Solid", col["D_contact_ss"], YELLOW),
           ("потери: пластика нити W_p_fiber", col["W_p_fiber"], RED),
           ("потери: вязкость нити W_v_fiber", col["W_v_fiber"], BROWN)]
    zero = [r[0] for r in rec if not nonzero(r[1])]
    rec = [r for r in rec if nonzero(r[1])]
    lost = [ke0 - k for k in col["KE_solid"]]
    tot = [sum(y[i] for _, y, _ in rec) for i in range(n)]
    eps = 1e-9 * max(ke0, 1e-300)
    share = [[100.0 * y[i] / tot[i] if tot[i] > eps else math.nan for i in range(n)] for _, y, _ in rec]
    fig, (a1, a2) = axes(2, 5.8, t)
    a1.stackplot(t, [y for _, y, _ in rec], colors=[c for _, _, c in rec], edgecolor=SURF, linewidth=1.0,
                 labels=["%s: %s Дж (%.0f %%)" % (nm, fmt(y[-1]), 100 * y[-1] / tot[-1] if tot[-1] > 0 else 0)
                         for nm, y, _ in rec])
    a1.plot(t, lost, color=INK, linewidth=1.5, linestyle="--",
            label="убыль KE Solid, KE0 - KE_solid: %s Дж" % fmt(lost[-1]))
    a1.text(0.02, 0.97, "разница стопки и пунктира (невязка + S_u_clamp) = %s Дж" % fmt(tot[-1] - lost[-1]),
            transform=a1.transAxes, fontsize=9, va="top")
    a1.set_ylim(0, max(tot[-1], lost[-1], eps) * 1.15)
    a1.set_ylabel("энергия, Дж")
    a1.set_title("Куда ушла кинетическая энергия Solid", loc="left", fontsize=12, color=INK)
    a2.stackplot(t, share, colors=[c for _, _, c in rec], edgecolor=SURF, linewidth=1.0)
    a2.set_ylim(0, 100)
    a2.set_ylabel("доля, %")
    a2.set_title("Те же получатели в долях", loc="left", fontsize=12, color=INK)
    h, lb = a1.get_legend_handles_labels()
    fig.legend(h, lb, loc="lower left", bbox_to_anchor=(0.06, 0.0), ncol=2, frameon=False, fontsize=9,
               title="получатель: значение в конце (доля)", title_fontsize=9)
    sub = ("\nнулевые на всём интервале (не показаны): " +
           ", ".join(z.replace("потери: ", "") for z in zero)) if zero else ""
    fig.suptitle(title + "\nубыль кинетической энергии Solid по получателям" + sub, x=0.01, ha="left", fontsize=11,
                 color=INK2)
    fig.tight_layout(rect=(0, 0.16, 1, 0.9))
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("outdir", nargs="?", default=".")
    ap.add_argument("--title", default=None)
    a = ap.parse_args()
    col = load(a.csv)
    t = [x * 1e6 for x in col["time"]]
    os.makedirs(a.outdir, exist_ok=True)
    title = a.title or "Лог энергии: " + os.path.basename(os.path.dirname(os.path.abspath(a.csv)))
    style()
    for fn, f in (("energy_channels.png", plot_channels), ("energy_dissipation.png", plot_dissipation),
                  ("energy_stocks.png", plot_stocks), ("energy_total.png", plot_total),
                  ("energy_transfer.png", plot_transfer)):
        p = os.path.join(a.outdir, fn)
        f(col, t, p, title)
        print("ok", p)


if __name__ == "__main__":
    main()
