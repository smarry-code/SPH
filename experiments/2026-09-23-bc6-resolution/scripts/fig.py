# -*- coding: utf-8 -*-
# Метрики и фигура лестницы BC6 по разрешению.
# Запуск:  python scripts/fig.py   (пути считаются от положения скрипта)
# Читает data/curve_dx*.txt и эталон tests/data/bc6_reference.txt;
# печатает сводную таблицу, пишет data/summary.txt и bc6_resolution.svg.
import os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
EXP  = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(EXP, "..", ".."))
D    = os.path.join(EXP, "data")
REF  = os.path.join(ROOT, "tests", "data", "bc6_reference.txt")
OUT  = os.path.join(EXP, "bc6_resolution.svg")

# (тег, Δx в мм)
RUNS = [("dx040", 0.40), ("dx028", 0.28), ("dx020", 0.20), ("dx014", 0.14), ("dx010", 0.10),
        ("dx005", 0.05)]   # 0,05 мм подготовлен (scripts/run_dx005.sh); пока кривой нет — пропускается
# Реплики при Δx = 0.2 мм: та же конфигурация, другой запуск. Расчёт
# недетерминирован (порядок суммирования OpenMP), и реплики дают оценку шума
# метрик. В подгонку по Δx и в левую панель фигуры не входят.
REPS = [("dx020_rep", 0.20, "повтор, 6 потоков"),
        ("dx020_bc6plot", 0.20, "цель bc6_plot, потоки по умолчанию"),
        ("dx020_keep", 0.20, "повтор с --keep-workdir, 6 потоков")]
T_IMPACTOR = 2.004   # мм, номинальная толщина ударника (tests/bc6_case.hpp)

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
    """Первое пересечение уровня снизу вверх начиная с i0, с интерполяцией."""
    for i in range(max(i0, 1), len(v)):
        if v[i - 1] < lvl <= v[i]:
            return t[i - 1] + (lvl - v[i - 1]) / (v[i] - v[i - 1]) * (t[i] - t[i - 1])
    return math.nan

def cross_down(t, v, lvl, i0, i1):
    """Последнее пересечение уровня сверху вниз на отрезке [i0, i1]."""
    for i in range(i1, max(i0, 1) - 1, -1):
        if v[i - 1] >= lvl > v[i]:
            return t[i - 1] + (v[i - 1] - lvl) / (v[i - 1] - v[i]) * (t[i] - t[i - 1])
    return math.nan

def med(t, v, a, b):
    m = (t >= a) & (t <= b)
    return float(np.median(v[m])) if m.any() else math.nan

def features(t, v):
    f = {}
    pk = v.max()
    t50 = cross_up(t, v, 0.5 * pk)
    pl1 = med(t, v, t50 + 0.15, t50 + 0.60)
    t10 = cross_up(t, v, 0.1 * pl1); t90 = cross_up(t, v, 0.9 * pl1)
    f.update(peak=pk, t50=t50, pl1=pl1, rise=(t90 - t10) * 1e3)
    # длительность плато 1: от фронта (t50) до первого падения ниже 0.9·плато
    # после выхода на плато (t50 + 0.1 мкс)
    i0 = np.searchsorted(t, t50 + 0.1)
    j = np.where(v[i0:] < 0.9 * pl1)[0]
    f["hold"] = t[i0 + j[0]] - t50 if len(j) else math.nan
    # Провал — локальный минимум ВНУТРИ окна поиска. Минимум на краю окна
    # означает монотонный спад, а не провал: тогда провала нет, и зависящие
    # от него величины не определены. Выраженность — высота меньшего из двух
    # соседних максимумов (0.4 мкс до, 0.3 мкс после) над дном.
    def dip(a, b):
        j = np.where((t > a) & (t < b))[0]
        if len(j) < 5: return None
        jd = j[np.argmin(v[j])]
        if jd - j[0] < 2 or j[-1] - jd < 2: return None
        before = v[(t >= t[jd] - 0.4) & (t <= t[jd])].max()
        after  = v[(t >= t[jd]) & (t <= t[jd] + 0.3)].max()
        return jd, min(before, after) - v[jd], after - v[jd]
    nan = math.nan
    f.update(dip1_t=nan, dip1_v=nan, fall=nan, prom1=nan, pullback=nan,
             dip2_t=nan, dip2_v=nan, prom2=nan)
    d1 = dip(t50 + 0.6, t50 + 1.2)
    if d1:
        jd, prom, pb = d1
        vd1 = v[jd]
        i_pl = np.searchsorted(t, t50 + 0.15)
        hi = cross_down(t, v, pl1 - 0.1 * (pl1 - vd1), i_pl, jd)
        lo = cross_down(t, v, pl1 - 0.9 * (pl1 - vd1), i_pl, jd)
        f.update(dip1_t=t[jd] - t50, dip1_v=vd1, fall=(lo - hi) * 1e3, prom1=prom, pullback=pb)
        ta = t[jd]
    else:
        ta = t50 + 0.9          # опорный момент для окон плато 2 и провала 2
    pl2 = med(t, v, ta + 0.25, ta + 0.65)
    f.update(pl2=pl2, ratio=pl2 / pl1)
    d2 = dip(ta + 0.6, ta + 1.1)
    if d2:
        jd, prom, pb = d2
        f.update(dip2_t=t[jd] - t50, dip2_v=v[jd], prom2=prom)
    return f

def load_particles():
    """data/particles.txt — вывод scripts/particles.cpp (bc6::plan_geometry)."""
    out = {}
    p = os.path.join(D, "particles.txt")
    if not os.path.exists(p): return out
    for L in open(p, encoding="utf-8"):
        if L.startswith("#") or not L.strip(): continue
        a = L.split("|")
        dx, tot = a[0].split()
        bodies = {}
        for part in a[1].split(";"):
            part = part.strip()
            if not part: continue
            name, spec = part.split()
            nrnz, cnt = spec.split("=")
            nr, nz = nrnz.split("×")
            bodies[name] = (int(nr), int(nz), int(cnt))
        out[round(float(dx), 2)] = (int(tot), bodies)
    return out

PART = load_particles()

def impactor_n(dx):
    return math.floor(T_IMPACTOR / dx + 1e-9) + 1

rt, rv, _, _ = load(REF)
R = features(rt, rv)
t_end_ref = rt[-1] - R["t50"]

def rms_vs_ref(t, v, f):
    """RMS отклонения от эталона после совмещения по t50, на всей эталонной записи."""
    tt = rt - R["t50"]
    m = (tt >= 0) & (tt <= t_end_ref)
    vi = np.interp(tt[m] + f["t50"], t, v)
    return float(np.sqrt(np.mean((vi - rv[m]) ** 2)))

rows = []
pending = []
for tag, dx in RUNS:
    p = os.path.join(D, "curve_%s.txt" % tag)
    if not os.path.exists(p):
        pending.append(tag); continue
    t, v, st, ab = load(p)
    f = features(t, v); f["rms"] = rms_vs_ref(t, v, f)
    n = impactor_n(dx)
    f.update(tag=tag, dx=dx, steps=st, aborted=ab, t=t, v=v,
             imp=n * dx, period_norm=f["dip1_t"] * T_IMPACTOR / (n * dx))
    rows.append(f)

reps = []
for tag, dx, note in REPS:
    p = os.path.join(D, "curve_%s.txt" % tag)
    if not os.path.exists(p): continue
    t, v, st, ab = load(p)
    f = features(t, v); f["rms"] = rms_vs_ref(t, v, f)
    n = impactor_n(dx)
    f.update(tag=tag, dx=dx, steps=st, aborted=ab, t=t, v=v, note=note,
             imp=n * dx, period_norm=f["dip1_t"] * T_IMPACTOR / (n * dx))
    reps.append(f)

def fmt(x, w=7, d=1):
    return ("%*.*f" % (w, d, x)) if not (x is None or (isinstance(x, float) and math.isnan(x))) else " " * (w - 1) + "—"

hdr = ("  Δx,мм  частиц n_z(уд/обр) шагов прерв  пик   плато1 длит1 фронт10-90 спуск10-90 провал1  t_пр1  выр1 возврат "
       "плато2 отн1→2 выр2  ударник t_пр1·2.004/h  RMS")
lines = [hdr]
lines.append("эталон        —         —         —     —  " + fmt(R["peak"]) + fmt(R["pl1"], 8) + fmt(R["hold"], 6, 3) + fmt(R["rise"], 10) +
             fmt(R["fall"], 11) + fmt(R["dip1_v"], 8) + fmt(R["dip1_t"], 7, 3) + fmt(R["prom1"], 6) +
             fmt(R["pullback"], 7) + fmt(R["pl2"], 8) + fmt(R["ratio"], 7, 3) + fmt(R["prom2"], 6) +
             "   2.004 " + fmt(R["dip1_t"], 9, 3) + "      —")
def part_cols(dx):
    if round(dx, 2) not in PART: return "%6.2f %8s %13s " % (dx, "—", "—")
    tot, b = PART[round(dx, 2)]
    return "%6.2f %8d %7s " % (dx, tot, "%d/%d" % (b["ударник"][1], b["образец"][1]))

for f in rows:
    lines.append(part_cols(f["dx"]) + "%6d %4s " % (f["steps"], f["aborted"]) + fmt(f["peak"]) + fmt(f["pl1"], 8) + fmt(f["hold"], 6, 3) +
                 fmt(f["rise"], 10) + fmt(f["fall"], 11) + fmt(f["dip1_v"], 8) + fmt(f["dip1_t"], 7, 3) +
                 fmt(f["prom1"], 6) + fmt(f["pullback"], 7) + fmt(f["pl2"], 8) + fmt(f["ratio"], 7, 3) +
                 fmt(f["prom2"], 6) + "   " + fmt(f["imp"], 5, 3) + fmt(f["period_norm"], 9, 3) +
                 fmt(f["rms"], 9, 2))

for f in reps:
    lines.append(part_cols(f["dx"]) + "%6d %4s " % (f["steps"], f["aborted"]) + fmt(f["peak"]) + fmt(f["pl1"], 8) + fmt(f["hold"], 6, 3) +
                 fmt(f["rise"], 10) + fmt(f["fall"], 11) + fmt(f["dip1_v"], 8) + fmt(f["dip1_t"], 7, 3) +
                 fmt(f["prom1"], 6) + fmt(f["pullback"], 7) + fmt(f["pl2"], 8) + fmt(f["ratio"], 7, 3) +
                 fmt(f["prom2"], 6) + "   " + fmt(f["imp"], 5, 3) + fmt(f["period_norm"], 9, 3) +
                 fmt(f["rms"], 9, 2) + "   ← реплика: " + f["note"])

# Где реплики расходятся с основным прогоном 0.2 мм
base = next((f for f in rows if f["tag"] == "dx020"), None)
if base is not None and reps:
    lines.append("")
    for f in reps:
        n = min(len(base["v"]), len(f["v"]))
        d = np.abs(base["v"][:n] - f["v"][:n])
        tb = base["t"][:n]
        first = tb[np.argmax(d > 0)] if (d > 0).any() else math.nan
        w = lambda a, b: np.sqrt(np.mean(d[(tb >= a) & (tb < b)] ** 2))
        lines.append("реплика %-14s против dx020: первое расхождение t = %.3f мкс; RMS разности "
                     "0-1 / 1-2 / 2-3 / 3-4.3 мкс = %.2f / %.2f / %.2f / %.2f м/с; максимум %.1f м/с"
                     % (f["tag"], first, w(0, 1), w(1, 2), w(2, 3), w(3, 4.31), d.max()))

def fit(key):
    ok = [(f["dx"], f[key]) for f in rows if f[key] > 0 and not math.isnan(f[key])]
    if len(ok) < 3: return math.nan, math.nan
    x = np.log([a for a, _ in ok]); y = np.log([b for _, b in ok])
    p, c = np.polyfit(x, y, 1)
    return p, c

fits = {}
lines.append("")
for key, name in (("rise", "фронт 10-90"), ("fall", "спуск в провал 1, 10-90")):
    p, c = fit(key)
    fits[key] = (p, c)
    if not math.isnan(p):
        need = math.exp((math.log(R[key]) - c) / p) if p > 0.05 else math.nan
        lines.append("%s: показатель p = %.2f  (время ∝ Δx^p); эталон %.1f нс достигается при Δx = %s мм (экстраполяция)"
                     % (name, p, R[key], ("%.3g" % need) if not math.isnan(need) else "— (p ≈ 0)"))
lines.append("провал — локальный минимум внутри окна поиска (фронт + 0.6…1.2 мкс для первого, "
             "провал 1 + 0.6…1.1 мкс для второго); выраженность — высота меньшего из соседних "
             "максимумов (0.4 мкс до, 0.3 мкс после) над дном; «—» — минимума внутри окна нет")
lines.append("частиц — всего, по bc6::plan_geometry (data/particles.txt); n_z(уд/обр) — слоёв "
             "по толщине ударника и образца")
lines.append("длит1 — время от фронта (50 %) до первого падения ниже 0.9 от уровня плато 1, мкс")
lines.append("t_пр1·2.004/h — момент провала 1, пересчитанный на номинальную толщину ударника "
             "(оценка: считает период целиком пропорциональным толщине ударника)")

if pending:
    lines.append("не посчитаны (кривой нет): " + ", ".join(pending))

txt = "\n".join(lines)
print(txt)
with open(os.path.join(D, "summary.txt"), "w", encoding="utf-8") as fo:
    fo.write(txt + "\n")

# ---------------------------------------------------------------- фигура
plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})
fig, ax = plt.subplots(1, 2, figsize=(15.5, 6.2), gridspec_kw={"width_ratios": [1.7, 1]})
cmap = plt.get_cmap("viridis")
a = ax[0]
a.plot(rt - R["t50"], rv, color="black", lw=2.4,
       label="эталон (прибор), плато %.2f мкс, спуск %.0f нс" % (R["hold"], R["fall"]), zorder=5)
for n, f in enumerate(rows):
    a.plot(f["t"] - f["t50"], f["v"], color=cmap(n / max(1, len(rows) - 1) * 0.9), lw=1.4,
           label="Δx = %.2f мм, плато %.2f мкс, %s" % (f["dx"], f["hold"],
               "провала 1 нет" if math.isnan(f["fall"]) else "спуск %.0f нс" % f["fall"]))
a.set_xlim(-0.3, t_end_ref + 0.1); a.set_ylim(-5, 215)
a.set_xlabel("время от прихода фронта (50 %), мкс"); a.set_ylabel("осевая скорость монитора, м/с")
a.set_title("BC6, тип контакта 0: кривые при пяти шагах решётки")
a.grid(alpha=0.3); a.legend(loc="lower left", fontsize=9)

b = ax[1]
dxs = np.array([f["dx"] for f in rows])
for key, mk, name in (("rise", "o", "фронт 10–90 %"), ("fall", "s", "спуск в провал 1, 10–90 %")):
    y = np.array([f[key] for f in rows])
    l, = b.loglog(dxs, y, mk + "-", lw=1.6, ms=6, label="%s, p = %.2f" % (name, fits[key][0]))
    b.axhline(R[key], color=l.get_color(), ls=":", lw=1.2)
    b.text(dxs.min() * 0.92, R[key] * 1.08, "эталон %.0f нс" % R[key], color=l.get_color(),
           fontsize=8, ha="right")
xx = np.array([0.04, 0.45])
y0 = rows[[f["dx"] for f in rows].index(0.20)]["fall"] if 0.20 in dxs else 200
for p, ls in ((1.0, "--"), (0.5, "-.")):
    b.loglog(xx, y0 * (xx / 0.2) ** p, color="grey", ls=ls, lw=0.9, label="наклон %.1f" % p)
b.set_xlim(0.035 if dxs.min() < 0.07 else 0.06, 0.5); b.set_ylim(5, 400)
b.set_xlabel("шаг решётки Δx, мм"); b.set_ylabel("время перехода, нс")
b.set_title("Крутизна переходов против Δx")
b.grid(alpha=0.3, which="both"); b.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(OUT, facecolor="white")
print("\nфигура:", OUT)
