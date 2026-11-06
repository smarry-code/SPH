#!/usr/bin/env python3
"""growth.py REF.csv PERT.csv [COL ...]
По energy_log двух прогонов: для каждой колонки (по умолчанию KE_fiber, U_bond, KE_solid)
шаг первого отличия, показатель роста лямбда на фазе 2 (наклон log10|d| в диапазоне
1e-13..1e-8, порядков на 100 шагов и e-кратное время в мкс), скачок (наибольший прирост
log10|d| за 50 шагов, шаг и высота) и уровень насыщения (медиана log10|d| за последние 1000 шагов).
"""
import csv, math, statistics, sys

def load(p):
    return {int(r["step"]): r for r in csv.DictReader(open(p))}

def main():
    ref, per = load(sys.argv[1]), load(sys.argv[2])
    cols = sys.argv[3:] or ["KE_fiber", "U_bond", "KE_solid"]
    steps = sorted(set(ref) & set(per))
    dt = statistics.median(float(ref[s]["dt"]) for s in steps[len(steps) // 2:])
    print(f"# dt ~ {dt:.3e} s, шагов {steps[-1]}")
    print("# col first_diff lambda[ord/100st] efold[us] fit_pts jump_step jump_h[ord] sat_log10")
    for c in cols:
        L = []
        for s in steps:
            d = abs(float(ref[s][c]) - float(per[s][c]))
            L.append((s, math.log10(d) if d > 0 else None))
        first = next((s for s, l in L if l is not None), None)
        best = (0, None)
        ld = dict((s, l) for s, l in L if l is not None)
        for s in ld:
            if s + 50 in ld and ld[s + 50] - ld[s] > best[0]:
                best = (ld[s + 50] - ld[s], s)
        # фаза 2: от плато фазы 1 (медиана 500 шагов после первого отличия) + 0,7 порядка
        # до плато + 5 порядков (не выше 1e-6), только до скачка
        plat = statistics.median([l for s, l in L if l is not None and first <= s < first + 500] or [0])
        jump = best[1] if best[1] is not None else steps[-1]
        pts = [(s, l) for s, l in L if l is not None and s < jump and plat + 0.7 <= l <= min(plat + 5, -6)]
        lam = ef = None
        if len(pts) >= 5:
            xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
            mx, my = statistics.mean(xs), statistics.mean(ys)
            sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
            lam = sl * 100
            ef = (1 / (sl * math.log(10))) * dt * 1e6 if sl > 0 else float("inf")
        tail = [l for s, l in L if l is not None and s > steps[-1] - 1000]
        sat = statistics.median(tail) if tail else None
        f = lambda v, fmt: (fmt.format(v) if v is not None else "-")
        print(f"{c} {first} {f(lam,'{:.2f}')} {f(ef,'{:.2f}')} {len(pts)} {best[1]} {best[0]:.1f} {f(sat,'{:.1f}')}")

if __name__ == "__main__":
    main()
