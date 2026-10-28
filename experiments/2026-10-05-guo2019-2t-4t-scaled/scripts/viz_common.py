"""Общие функции для plot_series.py и plot_shot.py (Guo 2019, серии 1T/2T/4T)."""
from __future__ import annotations
import csv, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 10
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.3

DEFAULT_CASE_HPP = "C:/dip/sph/SPH/.claude/worktrees/Guo2T/tests/guo2019_case.hpp"
CFG = {"1T": "T1", "2T": "T2", "4T": "T4"}
# e при V_BL статьи, Дж: ТЗ 2T/4T §3.1
PAPER_E_BL = {"1T": 20.2, "2T": 75.6, "4T": 233.0}
PAPER_VTILDE = {2: 1.9, 4: 3.4}   # Fig. 12a
PAPER_ETILDE = {2: 3.7, 4: 11.5}  # Fig. 12b
LAYER_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b"]


def save(fig, out_dir, name, png=False):
    os.makedirs(out_dir, exist_ok=True)
    p = os.path.join(out_dir, name + ".svg")
    fig.savefig(p, bbox_inches="tight")
    if png:
        fig.savefig(os.path.join(out_dir, name + ".png"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("записан", p)


def paper_data(case_hpp, series):
    """(shots[(V0, perforated, V_R)], (V_BL, +, -)) из tests/guo2019_case.hpp."""
    txt = open(case_hpp, encoding="utf-8").read()
    c = CFG[series]
    m = re.search(r"case Config::%s:\s*return \{(.*?)\};" % c, txt, re.S)
    shots = [(float(a), b == "true", float(r))
             for a, b, r in re.findall(r"\{\s*([\d.]+)\s*,\s*(true|false)\s*,\s*([\d.]+)\s*\}", m.group(1))]
    m = re.search(r"case Config::%s:\s*return \{([\d.]+),\s*([\d.]+),\s*([\d.]+)\}" % c, txt.split("paper_vbl")[1])
    return shots, tuple(float(x) for x in m.groups())


def read_shots(path):
    """shots.csv -> список словарей; None, если файла нет или он пуст."""
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r.get("V0")]
    return rows or None


def read_vbl(path):
    """vbl.txt -> dict(vbl, ebar); None, если не определены."""
    out = {"vbl": None, "ebar": None}
    if not os.path.isfile(path):
        return out
    t = open(path, encoding="utf-8").read()
    m = re.search(r"V_BL\s*=\s*([\d.]+)\s*м/с", t)
    if m:
        out["vbl"] = float(m.group(1))
        m = re.search(r"ē\s*=\s*([\d.eE+-]+)\s*Дж", t)
        if m:
            out["ebar"] = float(m.group(1))
    return out


def load_series(series_dir):
    """{'1T': (rows|None, vbl)}; неполные данные допустимы."""
    res = {}
    for s in ("1T", "2T", "4T"):
        d = os.path.join(series_dir, s)
        res[s] = (read_shots(os.path.join(d, "shots.csv")), read_vbl(os.path.join(d, "vbl.txt")))
    return res


def read_monitor(path):
    """monitor.csv; недописанная последняя строка отбрасывается."""
    with open(path) as f:
        names = f.readline().strip().split(",")
        rows = [ln.split(",") for ln in f]
    rows = [r for r in rows if len(r) == len(names)]
    d = np.array(rows, dtype=float).reshape(-1, len(names))
    return {n: d[:, i] for i, n in enumerate(names)}
