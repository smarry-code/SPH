# -*- coding: utf-8 -*-
# Метрики кривых монитора BC6 по четырём типам контакта при Δx = 0,1 мм
# и сравнение с теми же метриками при штатных 0,2 мм (build/bc6_curves —
# выгрузки цели bc6_plot; при их отсутствии сравнение пропускается).
# Запуск:  python scripts/metrics.py   → печать и data/summary.txt
import os, sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); EXP = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(EXP, "..", ".."))
SETS = [("0,1 мм", os.path.join(EXP, "data")), ("0,2 мм", os.path.join(ROOT, "build", "bc6_curves"))]

def load(p):
    T = []; V = []; st = ab = None
    for L in open(p, encoding="utf-8", errors="replace"):
        L = L.strip()
        if L.startswith("#"):
            if "steps" in L: st = int(L.split("=")[1])
            if "aborted" in L: ab = int(L.split("=")[1])
            continue
        if not L or L[0] == "*": continue
        a = L.split(); T.append(float(a[0]) * 1e6); V.append(float(a[1]))
    t = np.array(T); v = np.array(V); o = np.argsort(t, kind="stable")
    return t[o], v[o], st, ab

def cross_up(t, v, L):
    i = np.argmax(v >= L)
    return t[i - 1] + (L - v[i - 1]) / (v[i] - v[i - 1]) * (t[i] - t[i - 1]) if i > 0 else math.nan

def ripple(t, v, a, b):
    """СКО отклонения от скользящего среднего по 21 точке и автокорреляция лага 1.
    Среднее считается только там, где окно из 21 точки целиком внутри записи:
    режим 'same' дополняет края нулями и давал бы ложную рябь в конце записи."""
    k = 21; h = k // 2
    sm = np.full_like(v, np.nan)
    sm[h:len(v) - h] = np.convolve(v, np.ones(k) / k, "valid")
    m = (t >= a) & (t < b) & ~np.isnan(sm)
    if m.sum() < 30: return math.nan, math.nan
    r = (v - sm)[m]
    return r.std(), np.corrcoef(r[:-1], r[1:])[0, 1]

lines = []
for name, d in SETS:
    if not os.path.isdir(d): lines.append("%s: нет каталога %s" % (name, d)); continue
    lines.append("Δx = %s  (%s)" % (name, os.path.relpath(d, ROOT)))
    lines.append(" тип  шагов обрыв t_конец   пик   плато1  фронт10-90,нс   рябь 1,5–2,0 / 2,0–2,5 / 2,5–3,0 / 3,0–3,3 мкс, м/с")
    for ty in range(4):
        p = os.path.join(d, "bc6_curve_%d.txt" % ty)
        if not os.path.exists(p): continue
        t, v, st, ab = load(p)
        pk = v.max(); t50 = cross_up(t, v, 0.5 * pk)
        m = (t > t50 + 0.15) & (t < t50 + 0.6); pl = float(np.median(v[m]))
        rise = (cross_up(t, v, 0.9 * pl) - cross_up(t, v, 0.1 * pl)) * 1e3
        # Окна выбраны так, чтобы отделить устойчивую пилу (типа 2 — с 1,5 мкс)
        # от поздней неустойчивости, с которой при 0,1 мм зашумлены все типы
        # (у типа 3 — с 3,0 мкс, у типа 0 — с 3,3 мкс).
        rs = [ripple(t, v, a, b) for a, b in ((1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 3.3))]
        lines.append("  %d  %6s %4s %7.3f %7.2f %7.1f %11.1f        %s" % (
            ty, st, ab, t[-1], pk, pl, rise, " / ".join("%.3f" % s for s, _ in rs)))
    lines.append("")
lines.append("рябь — СКО отклонения от скользящего среднего по 21 точке (только там, где окно "
             "целиком внутри записи)")
txt = "\n".join(lines); print(txt)
open(os.path.join(EXP, "data", "summary.txt"), "w", encoding="utf-8").write(txt + "\n")
