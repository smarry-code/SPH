# -*- coding: utf-8 -*-
# Метрики и фигура серии BC6 по коэффициенту XSPH (ε) и шагу решётки.
# Запуск:  python scripts/fig.py   (пути считаются от положения скрипта)
# Читает data/curve_e<ε>_dx<Δx>.txt и эталон tests/data/bc6_reference.txt;
# печатает таблицы, пишет data/summary.txt и bc6_xsph.svg.
#
# ПОЧЕМУ МЕТРИКИ НЕ ТЕ, ЧТО ЗАДУМЫВАЛИСЬ. План (README, «Постановка») брал
# определения из 2026-09-23-bc6-resolution: пик, уровень плато 1 и все
# производные от них считаются по самой кривой. При ε < 0,25 решение после
# выхода на плато разваливается (всплески до 375 м/с, обрыв по ρ < 0), и эти
# метрики берут опорные уровни из развалившейся части. Поэтому здесь:
#   * фронт 10-90 % строится по ФИКСИРОВАННОМУ уровню — плато 1 кривой с
#     ε = 0,25 на той же сетке; фронт при всех ε лежит раньше развала;
#   * кривые совмещаются по t50 этого уровня и сравниваются с ε = 0,25:
#     момент отхода больше 10 и 20 м/с после выхода на плато (t50 + 0,15 мкс).
#     При ε = 0,25 реплики одной сетки до 1,7 мкс совпадают с точностью
#     печати, так что отход раньше этого — действие ε, а не шум;
#   * метрики stage A (длительность плато, спуск в провал, момент провала)
#     печатаются только для ε = 0,25.
import os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
EXP  = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(EXP, "..", ".."))
D    = os.path.join(EXP, "data")
REF  = os.path.join(ROOT, "tests", "data", "bc6_reference.txt")
OUT  = os.path.join(EXP, "bc6_xsph.svg")

EPS = [0.25, 0.10, 0.05]
DXS = [0.20, 0.14, 0.10, 0.05]
# e005_dx005 подготовлен, но НЕ ЗАПУСКАЕТСЯ на машине серии
# (scripts/run_e005_dx005.sh); его пара с ε = 0,25 — в соседнем эксперименте.
EXTERNAL = {(0.25, 0.05): os.path.join(ROOT, "experiments", "2026-09-23-bc6-resolution",
                                       "data", "curve_dx005.txt")}

def tag(eps, dx):
    return "e%03d_dx%03d" % (round(eps * 100), round(dx * 100))

def load(p):
    T = []; V = []; steps = None; aborted = None
    for L in open(p, encoding="utf-8", errors="replace"):
        L = L.strip().lstrip("﻿")
        if L.startswith("#"):
            if "steps" in L:   steps = int(L.split("=")[1])
            if "aborted" in L: aborted = int(L.split("=")[1])
            continue
        if not L or L[0] == "*": continue
        a = L.split()
        try: T.append(float(a[0]) * 1e6); V.append(float(a[1]))
        except (ValueError, IndexError): continue
    a = np.array(T); b = np.array(V)
    o = np.argsort(a, kind="stable")
    return a[o], b[o], steps, aborted

def cross_up(t, v, lvl, i0=0):
    for i in range(max(i0, 1), len(v)):
        if v[i - 1] < lvl <= v[i]:
            return t[i - 1] + (lvl - v[i - 1]) / (v[i] - v[i - 1]) * (t[i] - t[i - 1])
    return math.nan

def cross_down(t, v, lvl, i0, i1):
    for i in range(i1, max(i0, 1) - 1, -1):
        if v[i - 1] >= lvl > v[i]:
            return t[i - 1] + (v[i - 1] - lvl) / (v[i - 1] - v[i]) * (t[i] - t[i - 1])
    return math.nan

def med(t, v, a, b):
    m = (t >= a) & (t <= b)
    return float(np.median(v[m])) if m.any() else math.nan

def stage_a_features(t, v):
    """Определения stage A — только для ε = 0,25 и для эталона."""
    nan = math.nan; f = {}
    pk = v.max(); t50 = cross_up(t, v, 0.5 * pk)
    pl1 = med(t, v, t50 + 0.15, t50 + 0.60)
    f.update(peak=pk, t50=t50, pl1=pl1,
             rise=(cross_up(t, v, 0.9 * pl1) - cross_up(t, v, 0.1 * pl1)) * 1e3)
    i0 = np.searchsorted(t, t50 + 0.1); j = np.where(v[i0:] < 0.9 * pl1)[0]
    f["hold"] = t[i0 + j[0]] - t50 if len(j) else nan
    f.update(dip1_t=nan, fall=nan)
    j = np.where((t > t50 + 0.6) & (t < t50 + 1.2))[0]
    if len(j) >= 5:
        jd = j[np.argmin(v[j])]
        if jd - j[0] >= 2 and j[-1] - jd >= 2:
            vd = v[jd]; ip = np.searchsorted(t, t50 + 0.15)
            f.update(dip1_t=t[jd] - t50,
                     fall=(cross_down(t, v, pl1 - 0.9 * (pl1 - vd), ip, jd) -
                           cross_down(t, v, pl1 - 0.1 * (pl1 - vd), ip, jd)) * 1e3)
    return f

rt, rv, _, _ = load(REF)
R = stage_a_features(rt, rv)

curves = {}; pending = []
for dx in DXS:
    for eps in EPS:
        p = EXTERNAL.get((eps, dx)) or os.path.join(D, "curve_%s.txt" % tag(eps, dx))
        if os.path.exists(p):
            curves[(eps, dx)] = load(p)
        elif (eps, dx) in EXTERNAL or eps == 0.05 or (eps == 0.10 and dx == 0.20) or eps == 0.25:
            if not (eps == 0.10 and dx != 0.20):
                pending.append("%s%s" % (tag(eps, dx), " (resolution)" if (eps, dx) in EXTERNAL else ""))

# Уровень и t50 кривой ε = 0,25 на каждой сетке — опора для всех ε.
base = {}
for dx in DXS:
    if (0.25, dx) not in curves: continue
    t, v, _, _ = curves[(0.25, dx)]
    f = stage_a_features(t, v)
    L = f["pl1"]
    base[dx] = dict(L=L, t50=cross_up(t, v, 0.5 * L), t=t, v=v, sa=f)

rows = []
for (eps, dx), (t, v, st, ab) in sorted(curves.items(), key=lambda k: (k[0][1] * -1, -k[0][0])):
    if dx not in base: continue
    b = base[dx]; L = b["L"]
    t50 = cross_up(t, v, 0.5 * L)
    rise = (cross_up(t, v, 0.9 * L) - cross_up(t, v, 0.1 * L)) * 1e3
    s = t - t50
    d = v - np.interp(s + b["t50"], b["t"], b["v"])
    after = s > 0.15
    def first(thr):
        j = np.where(after & (np.abs(d) > thr))[0]
        return s[j[0]] if len(j) else math.nan
    w = (s > 0.15) & (s < 0.5)
    rows.append(dict(eps=eps, dx=dx, steps=st, aborted=ab, t_end=t[-1], s_end=s[-1],
                     rise=rise, shift=(t50 - b["t50"]) * 1e3, dplat=np.abs(d[w]).max(),
                     div10=first(10.0), div20=first(20.0), vmax=v.max(), t=t, v=v, t50=t50))

def fmt(x, w=7, d=1):
    return ("%*.*f" % (w, d, x)) if not (x is None or (isinstance(x, float) and math.isnan(x))) else " " * (w - 1) + "—"

lines = ["1. Фронт по фиксированному уровню (плато 1 кривой ε = 0,25 на той же сетке) и отход от ε = 0,25",
         "   ε    Δx,мм  шагов обрыв t_конец,мкс  уровень  фронт10-90,нс  сдвиг t50,нс  |Δv| плато  отход>10  отход>20  обрыв от t50  v_max"]
for f in rows:
    lines.append("%5.2f %6.2f %6s %5s %10.3f %8.1f %13s %13s %11s %9s %9s %12s %7.1f" % (
        f["eps"], f["dx"], f["steps"], f["aborted"], f["t_end"], base[f["dx"]]["L"],
        fmt(f["rise"], 13), fmt(f["shift"], 13), fmt(f["dplat"], 11),
        fmt(f["div10"], 9, 3), fmt(f["div20"], 9, 3),
        fmt(f["s_end"], 12, 3) if f["aborted"] else "           —", f["vmax"]))
lines.append("   сдвиг t50 — относительно ε = 0,25 на той же сетке; |Δv| плато — максимум отклонения от ε = 0,25 "
             "на t50 + 0,15…0,5 мкс; отход — первый момент после t50 + 0,15 мкс, когда |Δv| больше порога, "
             "мкс от t50; кривые совмещены по t50")

lines.append("")
fits = {}
for eps in EPS:
    pts = sorted([(f["dx"], f["rise"]) for f in rows if f["eps"] == eps and not math.isnan(f["rise"])])
    if len(pts) < 3:
        lines.append("ε = %.2f: точек фронта %d — показатель не считается" % (eps, len(pts))); continue
    p = np.polyfit(np.log([a for a, _ in pts]), np.log([b for _, b in pts]), 1)[0]
    fits[eps] = p
    lines.append("ε = %.2f: фронт ∝ Δx^p, p = %.2f по Δx = %s мм" %
                 (eps, p, " / ".join("%.2f" % a for a, _ in sorted(pts, reverse=True))))
for dx in DXS:
    a = next((f for f in rows if f["eps"] == 0.25 and f["dx"] == dx), None)
    c = next((f for f in rows if f["eps"] == 0.05 and f["dx"] == dx), None)
    if a and c:
        lines.append("Δx = %.2f мм: фронт ε 0,05 / ε 0,25 = %.3f" % (dx, c["rise"] / a["rise"]))

lines.append("")
lines.append("2. Метрики stage A (только ε = 0,25; при меньших ε решение разваливается после выхода на плато)")
lines.append("   Δx,мм   пик  плато1  длит1  фронт10-90  спуск10-90  t_пр1")
lines.append("эталон  " + fmt(R["peak"]) + fmt(R["pl1"]) + fmt(R["hold"], 7, 3) + fmt(R["rise"], 12) +
             fmt(R["fall"], 12) + fmt(R["dip1_t"], 7, 3))
for dx in DXS:
    if dx in base:
        f = base[dx]["sa"]
        lines.append("%6.2f  " % dx + fmt(f["peak"]) + fmt(f["pl1"]) + fmt(f["hold"], 7, 3) +
                     fmt(f["rise"], 12) + fmt(f["fall"], 12) + fmt(f["dip1_t"], 7, 3))
if pending:
    lines.append("")
    lines.append("не посчитаны (кривой нет): " + ", ".join(pending))

txt = "\n".join(lines)
print(txt)
with open(os.path.join(D, "summary.txt"), "w", encoding="utf-8") as fo:
    fo.write(txt + "\n")

# ---------------------------------------------------------------- фигура
plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})
fig, ax = plt.subplots(1, 3, figsize=(18, 6.0), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
col = {0.25: "#1f77b4", 0.10: "#ff7f0e", 0.05: "#d62728"}
for a, dx in zip(ax[:2], (0.20, 0.10)):
    a.plot(rt - R["t50"], rv, color="black", lw=2.2, label="эталон (прибор)", zorder=5)
    for f in rows:
        if f["dx"] != dx: continue
        a.plot(f["t"] - f["t50"], f["v"], color=col[f["eps"]], lw=1.4,
               label="ε = %.2f%s" % (f["eps"], ", обрыв при %.2f мкс" % f["t_end"] if f["aborted"] else ""))
        if f["aborted"]:
            a.plot(f["t"][-1] - f["t50"], f["v"][-1], "x", color=col[f["eps"]], ms=11, mew=2.5)
    a.axvline(1.7 - base[dx]["t50"], color="grey", ls=":", lw=1)
    a.text(1.7 - base[dx]["t50"] + 0.03, 385, "1,7 мкс: дальше при ε = 0,25\nреплики расходятся",
           color="grey", fontsize=8, va="top")
    a.set_xlim(-0.3, 2.2); a.set_ylim(-5, 400)
    a.set_xlabel("время от прихода фронта (50 %), мкс"); a.set_ylabel("осевая скорость монитора, м/с")
    a.set_title("Δx = %.2f мм: кривые при разных ε" % dx)
    a.grid(alpha=0.3); a.legend(loc="upper right", fontsize=8.5, bbox_to_anchor=(1.0, 0.86))
b = ax[2]
for eps in EPS:
    pts = sorted([(f["dx"], f["rise"]) for f in rows if f["eps"] == eps])
    if not pts: continue
    x, y = zip(*pts)
    b.loglog(x, y, "o-", color=col[eps], lw=1.5, ms=6,
             label="ε = %.2f%s" % (eps, "" if eps not in fits else ", p = %.2f" % fits[eps]))
b.axhline(R["rise"], color="black", ls=":", lw=1)
b.text(0.093, R["rise"] * 1.1, "эталон %.0f нс" % R["rise"], fontsize=8)
b.set_xlim(0.09, 0.22); b.set_ylim(10, 300)
b.set_xticks([0.1, 0.14, 0.2]); b.set_xticklabels(["0,10", "0,14", "0,20"])
b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
b.set_xlabel("шаг решётки Δx, мм"); b.set_ylabel("фронт 10–90 %, нс")
b.set_title("Фронт по фиксированному уровню")
b.grid(alpha=0.3, which="both"); b.legend(fontsize=8.5, loc="lower right")
fig.tight_layout()
fig.savefig(OUT, facecolor="white")
print("\nфигура:", OUT)
