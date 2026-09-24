# -*- coding: utf-8 -*-
# Метрики и фигура серии BC6 по линейному коэффициенту искусственной вязкости α.
# Запуск:  python scripts/fig.py   (пути считаются от положения скрипта)
# Читает data/curve_a<α>_dx020.txt и эталон tests/data/bc6_reference.txt;
# печатает таблицы, пишет data/summary.txt и bc6_art_visc.svg.
#
# Определения — как в experiments/2026-09-23-bc6-xsph/scripts/fig.py:
#   * фронт 10-90 % по ФИКСИРОВАННОМУ уровню — плато 1 кривой α = 1,0;
#   * кривые совмещаются по t50 этого уровня и сравниваются с α = 1,0;
#   * метрики stage A (длительность плато, спуск, провал) — по своему плато
#     каждой кривой; осмысленны, только если отход от α = 1,0 наступает позже.
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
OUT  = os.path.join(EXP, "bc6_art_visc.svg")

ALPHAS = [1.0, 0.5, 0.25]
# Прогнозы, записанные до прогонов (README, «Постановка»), нс
PRED = {"остаток ∝ α, квадратура": {0.5: 150.5, 0.25: 129.5},
        "остаток ∝ √α, сложение":  {0.5: 140.2, 0.25: 108.1},
        "вся ширина ∝ √α":         {0.5: 131.2, 0.25: 92.8}}

def tag(a): return "a%03d_dx020" % round(a * 100)

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
    nan = math.nan; f = {}
    pk = v.max(); t50 = cross_up(t, v, 0.5 * pk)
    pl1 = med(t, v, t50 + 0.15, t50 + 0.60)
    f.update(peak=pk, t50=t50, pl1=pl1,
             rise=(cross_up(t, v, 0.9 * pl1) - cross_up(t, v, 0.1 * pl1)) * 1e3)
    i0 = np.searchsorted(t, t50 + 0.1); j = np.where(v[i0:] < 0.9 * pl1)[0]
    f["hold"] = t[i0 + j[0]] - t50 if len(j) else nan
    f.update(dip1_t=nan, dip1_v=nan, fall=nan)
    j = np.where((t > t50 + 0.6) & (t < t50 + 1.2))[0]
    if len(j) >= 5:
        jd = j[np.argmin(v[j])]
        if jd - j[0] >= 2 and j[-1] - jd >= 2:
            vd = v[jd]; ip = np.searchsorted(t, t50 + 0.15)
            f.update(dip1_t=t[jd] - t50, dip1_v=vd,
                     fall=(cross_down(t, v, pl1 - 0.9 * (pl1 - vd), ip, jd) -
                           cross_down(t, v, pl1 - 0.1 * (pl1 - vd), ip, jd)) * 1e3)
    return f

rt, rv, _, _ = load(REF)
R = stage_a_features(rt, rv)

curves = {}; pending = []
for a in ALPHAS:
    p = os.path.join(D, "curve_%s.txt" % tag(a))
    if os.path.exists(p): curves[a] = load(p)
    else: pending.append(tag(a))
if 1.0 not in curves:
    sys.exit("нет опорной кривой α = 1,0 (data/curve_a100_dx020.txt)")

tb, vb, _, _ = curves[1.0]
sa_b = stage_a_features(tb, vb)
L = sa_b["pl1"]
t50b = cross_up(tb, vb, 0.5 * L)

rows = []
for a in ALPHAS:
    if a not in curves: continue
    t, v, st, ab = curves[a]
    t50 = cross_up(t, v, 0.5 * L)
    rise = (cross_up(t, v, 0.9 * L) - cross_up(t, v, 0.1 * L)) * 1e3
    s = t - t50
    d = v - np.interp(s + t50b, tb, vb)
    after = s > 0.15
    def first(thr):
        j = np.where(after & (np.abs(d) > thr))[0]
        return s[j[0]] if len(j) else math.nan
    w = (s > 0) & (s < 0.3); ip = np.argmax(v[w]); pk = v[w][ip]; tpk = s[w][ip]
    vmin = v[(s > tpk) & (s < tpk + 0.2)].min()
    rows.append(dict(a=a, steps=st, aborted=ab, t_end=t[-1], rise=rise,
                     shift=(t50 - t50b) * 1e3, peak=pk, swing=pk - vmin,
                     plat=med(s, v, 0.15, 0.30), div10=first(10.0), div20=first(20.0),
                     sa=stage_a_features(t, v), t=t, v=v, t50=t50))

def fmt(x, w=7, d=1):
    return ("%*.*f" % (w, d, x)) if not (x is None or (isinstance(x, float) and math.isnan(x))) else " " * (w - 1) + "—"

lines = ["1. Фронт по фиксированному уровню %.1f м/с (плато 1 кривой α = 1,0), Δx = 0,20 мм, ε = 0,25" % L,
         "    α   шагов обрыв t_конец  фронт10-90,нс  доля от α=1  сдвиг t50,нс  выброс  размах  плато(t50+.15….30)  отход>10  отход>20"]
for f in rows:
    lines.append("%5.2f %7s %5s %8.3f %14s %12.3f %13s %7.1f %7.1f %19s %9s %9s" % (
        f["a"], f["steps"], f["aborted"], f["t_end"], fmt(f["rise"], 14), f["rise"] / rows[0]["rise"],
        fmt(f["shift"], 13), f["peak"], f["swing"], fmt(f["plat"], 19),
        fmt(f["div10"], 9, 3), fmt(f["div20"], 9, 3)))
lines.append("   выброс — максимум за фронтом (t50…t50+0,3); размах — выброс минус первый минимум после него;"
             " отход — мкс от t50, когда |v − v(α=1)| впервые больше порога после t50 + 0,15")
lines.append("")
ok = [(f["a"], f["rise"]) for f in rows if not math.isnan(f["rise"])]
if len(ok) >= 3:
    q = np.polyfit(np.log([a for a, _ in ok]), np.log([r for _, r in ok]), 1)[0]
    lines.append("фронт ∝ α^q по α = %s: q = %.2f  (√α дало бы 0,5)" %
                 (" / ".join("%.2f" % a for a, _ in ok), q))
else:
    q = math.nan
lines.append("прогнозы, записанные до прогонов, нс: " + "; ".join(
    "%s — %s" % (k, " / ".join("α %.2f: %.0f" % (a, v) for a, v in d.items())) for k, d in PRED.items()))
lines.append("")
lines.append("2. Метрики stage A по собственному плато каждой кривой")
lines.append("    α     пик  плато1  длит1  фронт10-90  спуск10-90  провал1: t / v")
lines.append("эталон " + fmt(R["peak"]) + fmt(R["pl1"]) + fmt(R["hold"], 7, 3) + fmt(R["rise"], 12) +
             fmt(R["fall"], 12) + "   %s / %s" % (fmt(R["dip1_t"], 5, 3), fmt(R["dip1_v"], 5)))
for f in rows:
    s_ = f["sa"]
    lines.append("%5.2f  " % f["a"] + fmt(s_["peak"]) + fmt(s_["pl1"]) + fmt(s_["hold"], 7, 3) +
                 fmt(s_["rise"], 12) + fmt(s_["fall"], 12) +
                 "   %s / %s" % (fmt(s_["dip1_t"], 5, 3), fmt(s_["dip1_v"], 5)))
if pending:
    lines.append(""); lines.append("не посчитаны (кривой нет): " + ", ".join(pending))

txt = "\n".join(lines)
print(txt)
with open(os.path.join(D, "summary.txt"), "w", encoding="utf-8") as fo:
    fo.write(txt + "\n")

# ---------------------------------------------------------------- фигура
plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})
fig, ax = plt.subplots(1, 3, figsize=(18, 6.0), gridspec_kw={"width_ratios": [1.4, 1, 0.9]})
col = {1.0: "#1f77b4", 0.5: "#ff7f0e", 0.25: "#d62728"}
for k, (a_, xlim, ylim, title) in enumerate((
        (ax[0], (-0.3, 2.2), (-5, 230), "Δx = 0,20 мм, ε = 0,25: кривые при разных α"),
        (ax[1], (-0.25, 0.45), (-5, 230), "Фронт крупно"))):
    a_.plot(rt - R["t50"], rv, color="black", lw=2.2, label="эталон (прибор)", zorder=5)
    for f in rows:
        a_.plot(f["t"] - f["t50"], f["v"], color=col[f["a"]], lw=1.5,
                label="α = %.2f, фронт %.0f нс%s" % (f["a"], f["rise"],
                      ", обрыв при %.2f мкс" % f["t_end"] if f["aborted"] else ""))
        if f["aborted"] and f["t"][-1] - f["t50"] < xlim[1]:
            a_.plot(f["t"][-1] - f["t50"], f["v"][-1], "x", color=col[f["a"]], ms=11, mew=2.5)
    if k == 0:
        a_.axvline(1.7 - t50b, color="grey", ls=":", lw=1)
        a_.text(1.7 - t50b + 0.03, 222, "1,7 мкс: дальше реплики\nрасходятся", color="grey",
                fontsize=8, va="top")
    a_.set_xlim(*xlim); a_.set_ylim(*ylim)
    a_.set_xlabel("время от прихода фронта (50 %), мкс"); a_.set_ylabel("осевая скорость монитора, м/с")
    a_.set_title(title); a_.grid(alpha=0.3); a_.legend(loc="lower right", fontsize=8.5)
b = ax[2]
x = np.array([f["a"] for f in rows]); y = np.array([f["rise"] for f in rows])
b.loglog(x, y, "o-", color="black", lw=1.6, ms=7, zorder=5,
         label="расчёт%s" % ("" if math.isnan(q) else ", q = %.2f" % q))
xx = np.array([0.2, 1.1])
b.loglog(xx, rows[0]["rise"] * np.sqrt(xx), color="grey", ls="--", lw=1, label="∝ √α от α = 1")
for (name, d), mk in zip(PRED.items(), ("s", "^", "v")):
    b.loglog(list(d), list(d.values()), mk, mfc="none", color="grey", ms=7, label="прогноз: " + name)
b.axhline(R["rise"], color="black", ls=":", lw=1)
b.text(0.21, R["rise"] * 1.1, "эталон %.0f нс" % R["rise"], fontsize=8)
b.set_xlim(0.2, 1.1); b.set_ylim(10, 300)
b.set_xticks([0.25, 0.5, 1.0]); b.set_xticklabels(["0,25", "0,5", "1"])
b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
b.set_xlabel("линейный коэффициент вязкости α"); b.set_ylabel("фронт 10–90 %, нс")
b.set_title("Фронт против α"); b.grid(alpha=0.3, which="both"); b.legend(fontsize=7.5, loc="lower right")
fig.tight_layout()
fig.savefig(OUT, facecolor="white")
print("\nфигура:", OUT)
