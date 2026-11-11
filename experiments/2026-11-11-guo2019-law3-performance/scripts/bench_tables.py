#!/usr/bin/env python3
"""Таблицы бенчмарка bench_8t.sh: с/шаг по окнам шагов (по wall_samples.txt) и ускорение по числу потоков.

wall_samples.txt -- строки "<unix-время> <шаг>,<t>" раз в 10 с. Окна по шагам: до контакта и удар
(границы -- как в §7.1 отчёта закона 3). Для окон perf выводятся доли длительности прогона, на которых
окно начинается и кончается (аргумент --win для perf_threads.py).

Запуск: python3 bench_tables.py <каталог результатов> [--pre 130-1300] [--imp 1700-999999]
"""
import argparse
import os
import re


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 2 or "," not in p[1]:
                continue
            try:
                step = int(p[1].split(",")[0])
            except ValueError:
                continue
            rows.append((float(p[0]), step))
    return rows


def interp_time(rows, step):
    for (ta, sa), (tb, sb) in zip(rows, rows[1:]):
        if sa <= step <= sb and sb > sa:
            return ta + (tb - ta) * (step - sa) / (sb - sa)
    return None


def window(rows, lo, hi):
    """с/шаг в окне шагов [lo, hi]; hi обрезается последним отмеченным шагом."""
    if len(rows) < 2:
        return None, None
    hi = min(hi, rows[-1][1])
    ta, tb = interp_time(rows, lo), interp_time(rows, hi)
    if ta is None or tb is None or hi <= lo:
        return None, None
    return (tb - ta) / (hi - lo), (ta, tb)


def total(run):
    try:
        s = open(os.path.join(run, "stdout.log"), errors="replace").read()
    except OSError:
        return None, None
    m = re.search(r"шагов=(\d+), wall=([\d.]+) с", s)
    return (int(m.group(1)), float(m.group(2))) if m else (None, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--pre", default="130-1300")
    ap.add_argument("--imp", default="1700-999999")
    a = ap.parse_args()
    pre = tuple(map(int, a.pre.split("-")))
    imp = tuple(map(int, a.imp.split("-")))
    runs = sorted(d for d in os.listdir(a.root) if os.path.isfile(os.path.join(a.root, d, "wall_samples.txt")))
    res = {}
    print(f"{'прогон':10s} {'шагов':>6s} {'wall,с':>8s} {'с/шаг':>8s} {'до конт.':>9s} {'удар':>9s}  окна perf (доли)")
    for r in runs:
        d = os.path.join(a.root, r)
        rows = load(os.path.join(d, "wall_samples.txt"))
        n, w = total(d)
        p, pw = window(rows, *pre)
        i, iw = window(rows, *imp)
        res[r] = (n, w, p, i)
        span = ""
        if rows and pw and iw:
            t0, t1 = rows[0][0], rows[-1][0]
            f = lambda x: (x - t0) / (t1 - t0)
            span = f"pre:{f(pw[0]):.2f}-{f(pw[1]):.2f} imp:{f(iw[0]):.2f}-{f(iw[1]):.2f}"
        fmt = lambda x, k: f"{x:{k}.4f}" if x is not None else f"{'-':>{k}s}"
        print(f"{r:10s} {n if n else '-':>6} {w if w else '-':>8} {fmt(w / n if n and w else None, '8')} "
              f"{fmt(p, '9')} {fmt(i, '9')}  {span}")
    print("\nмасштабирование (ускорение S = T1/TN, эффективность S/N, f_s по Амдалу из N = 8):")
    for v in ("H", "P"):
        base = res.get(f"{v}.t1")
        if not base:
            continue
        for k, idx in (("до контакта", 2), ("удар", 3), ("весь прогон", None)):
            def val(x):
                if x is None:
                    return None
                if idx is None:
                    return x[1] / x[0] if x[0] and x[1] else None
                return x[idx]
            t1 = val(base)
            line = f"  {v} {k:12s} T1={t1:.4f}" if t1 else f"  {v} {k}"
            for nt in (2, 4, 8):
                x = val(res.get(f"{v}.t{nt}"))
                if t1 and x:
                    s = t1 / x
                    line += f"  N={nt}: S={s:.2f} E={s / nt:.2f}"
                    if nt == 8:
                        line += f" f_s={(1 / s - 1 / nt) / (1 - 1 / nt):.3f}"
            print(line)


if __name__ == "__main__":
    main()
