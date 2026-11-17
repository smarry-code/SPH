#!/usr/bin/env python3
"""Общие читатели для metrics.py / frames.py / plots.py (серия «закон 3 как проекция скорости»).
Читатели кадров и раскладки — по образцу experiments/2026-11-02-guo2019-thread-divergence/scripts/diff_frames.py
(read_frame / layout), но на numpy."""
import csv, glob, os, re
import numpy as np

RUNS = ["C1", "C2", "C3", "C4", "C5", "C1_k1", "C1_k2", "C1_k3"]


def find_runs(root):
    """Каталоги прогонов в root (плоско: root/<RUN>/shot.txt)."""
    return [r for r in RUNS if os.path.isfile(os.path.join(root, r, "shot.txt"))]


def read_shot(path):
    """shot.txt: одна строка key=value через пробел."""
    d = {}
    with open(path) as fh:
        for tok in fh.read().split():
            if "=" not in tok:
                continue
            k, v = tok.split("=", 1)
            try:
                d[k] = float(v)
            except ValueError:
                d[k] = v
    return d


def read_csv_dedup(path):
    """CSV с заголовком -> (имена, dict колонка->np.array). Повторы шагов (рестарты): остаётся
    последнее вхождение; пустые/недостающие поля -> nan. Строки отсортированы по шагу."""
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        head = next(rd)
        rows = {}
        for r in rd:
            if not r:
                continue
            r = r + [""] * (len(head) - len(r))
            vals = []
            for x in r[:len(head)]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(np.nan)
            rows[int(vals[0])] = vals
    steps = sorted(rows)
    A = np.array([rows[s] for s in steps]) if steps else np.zeros((0, len(head)))
    return head, {h: A[:, i] for i, h in enumerate(head)}


def read_fiber_frame(path):
    """Кадр нитей -> dict(t, step, X[n,3], V[n,3], obj[n], mat[n], bi[nb], bj[nb])."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    h = lines[0].split()
    t, step, n, nb = float(h[1]), int(h[3]), int(h[5]), int(h[7])
    ip = lines.index("PARTICLES")
    P = np.array([l.split() for l in lines[ip + 1: ip + 1 + n]], dtype=float)
    ib = lines.index("BONDS")
    B = np.array([l.split()[:2] for l in lines[ib + 1: ib + 1 + nb]], dtype=np.int64).reshape(-1, 2)
    return dict(t=t, step=step, X=P[:, 3:6], V=P[:, 6:9], obj=P[:, 1].astype(int), mat=P[:, 2].astype(int),
                bi=B[:, 0], bj=B[:, 1])


def read_sph_frame(path):
    """Кадр шара: колонки N Nmat M x y z vx vy vz rho p u c sigma_i e_pl D hsml ... -> dict(m, X, V, rho)."""
    A = np.loadtxt(path, skiprows=1, usecols=(2, 3, 4, 5, 6, 7, 8, 9))
    return dict(m=A[:, 0], X=A[:, 1:4], V=A[:, 4:7], rho=A[:, 7])


def read_layout(path):
    """guo2019_fabric_layout.txt: # object_id direction index offset clamped first_bead last_bead ..."""
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        out.append(dict(obj=int(p[0]), dir=p[1], idx=int(p[2]), off=float(p[3]), clamped=int(p[4]),
                        first=int(p[5]), last=int(p[6])))
    return out


def frame_steps(rundir):
    d = {}
    for f in glob.glob(os.path.join(rundir, "fiber_output_*.txt")):
        d[int(re.search(r"_(\d+)\.txt$", f).group(1))] = f
    return dict(sorted(d.items()))


def fiber_diameter(rundir, default=3.9e-4):
    """Диаметр нити из шапки guo2019_fabric.fib: '# d=...' (MATERIAL-строка: id E nu d rho ...)."""
    p = os.path.join(rundir, "guo2019_fabric.fib")
    try:
        for line in open(p):
            if line.startswith("MATERIAL"):
                return float(line.split()[4])
    except OSError:
        pass
    return default
