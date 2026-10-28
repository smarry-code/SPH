#!/usr/bin/env python3
"""Графики серий Guo 2019 1T/2T/4T (ТЗ 2T/4T §11).

--series-dir содержит подкаталоги 1T, 2T, 4T с shots.csv и vbl.txt (любой может отсутствовать).
Выход (SVG, при --png ещё и PNG):
    vr_v0.svg        -- V_R(V_0) модели поверх точек статьи (как Fig. 10), три панели;
    vbl_energy_n.svg -- Ṽ_BL(n) и ẽ(n), логарифмические оси, прямые статьи n^0.887, n^1.773;
    --table FILE     -- сводная таблица (.md или .csv по расширению).
Usage:
    python plot_series.py --series-dir server_series --out plots_out [--png] [--table t.md]
"""
from __future__ import annotations
import argparse, csv, os
import numpy as np
from viz_common import *


def f(r, k):
    try:
        return float(r[k])
    except (ValueError, TypeError, KeyError):
        return float("nan")


def plot_vr(data, case_hpp, out, png):
    fig, axs = plt.subplots(1, 3, figsize=(15, 6.2))
    for ax, s in zip(axs, ("1T", "2T", "4T")):
        shots, (vp, pp, pm) = paper_data(case_hpp, s)
        rows, vbl = data[s]
        ax.axvspan(vp - pm, vp + pp, color="k", alpha=0.08)
        ax.axvline(vp, color="k", ls="--", lw=1, label=f"V_BL статьи = {vp:g} м/с")
        first = True
        for v0, perf, vr in shots:
            if perf:
                ax.plot(v0, vr, "o", mfc="k", mec="k", ms=7, label="статья: пробит (V_R)" if first else None)
                first = False
            else:
                ax.plot(v0, 0, "o", mfc="none", mec="k", ms=7)
        ax.plot([], [], "o", mfc="none", mec="k", ms=7, label="статья: не пробит (V_R = 0)")
        if rows:
            ok = [r for r in rows if r["outcome"] == "perforated"]
            no = [r for r in rows if r["outcome"] == "not_perforated"]
            un = [r for r in rows if r["outcome"] not in ("perforated", "not_perforated")]
            if ok:
                ax.plot([f(r, "V0") for r in ok], [f(r, "V_R") for r in ok], "s", color="#d62728", ms=6,
                        label=f"модель: пробит ({len(ok)})")
            if no:
                ax.plot([f(r, "V0") for r in no], [0] * len(no), "s", mfc="none", mec="#d62728", ms=6,
                        label=f"модель: не пробит ({len(no)})")
            if un:
                ax.plot([f(r, "V0") for r in un], [0] * len(un), "D", color="#7f7f7f", ms=6,
                        label=f"модель: исход не определён, V_R = 0 условно ({len(un)})")
            if vbl["vbl"] is not None:
                ax.axvline(vbl["vbl"], color="#d62728", ls=":", lw=1.5, label=f"V_BL модели = {vbl['vbl']:.1f} м/с")
            else:
                ax.text(0.98, 0.5, "V_BL модели не определён\n(нет «не пробит»)",
                        transform=ax.transAxes, ha="right", fontsize=8, color="#d62728")
        else:
            ax.text(0.5, 0.5, "данных модели пока нет", transform=ax.transAxes, ha="center", color="#7f7f7f")
        ax.set_title(s)
        ax.set_xlabel("начальная скорость V_0, м/с")
        ax.set_ylabel("остаточная скорость V_R, м/с")
        ax.set_ylim(-8, None)
        ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    fig.suptitle("V_R(V_0): модель (красный) и опыт Guo 2019, Fig. 10 (чёрный)")
    save(fig, out, "vr_v0", png)


def fit_b(ns, vals):
    ln, lv = np.log(ns), np.log(vals)
    return float(np.sum(ln * lv) / np.sum(ln * ln))


def plot_norm(data, out, png):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
    xs = np.linspace(1, 4.4, 50)
    items = [("Ṽ_BL = V_BL(nT)/V_BL(1T)", 0.887, PAPER_VTILDE, "vbl", "V_BL"),
             ("ẽ = e(nT)/e(1T)", 1.773, PAPER_ETILDE, "ebar", "ē")]
    base = data["1T"][1]
    for ax, (lab, b, pp, key, nm) in zip(axs, items):
        ax.plot(xs, xs ** b, "k-", lw=1, label=f"статья: n^{b}")
        ax.plot(list(pp), list(pp.values()), "ko", ms=7, label="статья: точки")
        pts = []
        if base[key] is not None:
            for n, s in ((2, "2T"), (4, "4T")):
                if data[s][1][key] is not None:
                    pts.append((n, data[s][1][key] / base[key]))
        if pts:
            ns, vs = zip(*pts)
            ax.plot(ns, vs, "s", color="#d62728", ms=7, label="модель")
            bm = fit_b(np.array(ns), np.array(vs))
            ax.plot(xs, xs ** bm, "-", color="#d62728", lw=1,
                    label=f"модель: n^{bm:.3f} (по {len(pts)} точкам)")
        miss = [f"{s}: {nm} не определён или серии нет" for s in ("1T", "2T", "4T") if data[s][1][key] is None]
        if miss:
            ax.text(0.98, 0.04, "Не показано:\n" + "\n".join(miss), transform=ax.transAxes,
                    ha="right", fontsize=8, color="#d62728")
        ax.set_xscale("log"); ax.set_yscale("log")
        from matplotlib.ticker import ScalarFormatter, NullFormatter
        ax.yaxis.set_major_formatter(ScalarFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xticks([1, 2, 4]); ax.set_xticklabels(["1", "2", "4"])
        ax.set_xlabel("число слоёв n"); ax.set_ylabel(lab + ", безразм.")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(True, which="both", alpha=0.3)
    fig.suptitle("Нормировка на 1T: баллистический предел и энергия при V_BL")
    save(fig, out, "vbl_energy_n", png)


def table(data, path):
    hdr = ["серия", "V_0, м/с", "повтор", "исход", "V_R, м/с", "e, Дж", "ē серии, Дж", "V_BL модели, м/с"]
    names = {"perforated": "пробит", "not_perforated": "не пробит", "undetermined": "не определён"}
    out = []
    for s in ("1T", "2T", "4T"):
        rows, vbl = data[s]
        if not rows:
            out.append([s, "—", "—", "данных нет", "—", "—", "—", "—"])
            continue
        for r in sorted(rows, key=lambda r: (f(r, "V0"), f(r, "rep"))):
            vr, e = f(r, "V_R"), f(r, "e")
            out.append([s, f"{f(r, 'V0'):g}", f"{int(f(r, 'rep'))}", names.get(r["outcome"], r["outcome"]),
                        "—" if vr != vr else f"{vr:.1f}", "—" if e != e else f"{e:.2f}",
                        "—" if vbl["ebar"] is None else f"{vbl['ebar']:.2f}",
                        "не определён" if vbl["vbl"] is None else f"{vbl['vbl']:.1f}"])
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if path.endswith(".csv"):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh); w.writerow(hdr); w.writerows(out)
    else:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr) + "\n")
            for r in out:
                fh.write("| " + " | ".join(r) + " |\n")
    print("записан", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series-dir", required=True)
    ap.add_argument("--out", default="plots_out")
    ap.add_argument("--table")
    ap.add_argument("--png", action="store_true")
    ap.add_argument("--case-hpp", default=DEFAULT_CASE_HPP)
    a = ap.parse_args()
    data = load_series(a.series_dir)
    for s, (rows, v) in data.items():
        print(s, "выстрелов:", 0 if not rows else len(rows), "V_BL:", v["vbl"], "ē:", v["ebar"])
    plot_vr(data, a.case_hpp, a.out, a.png)
    plot_norm(data, a.out, a.png)
    if a.table:
        table(data, a.table)


if __name__ == "__main__":
    main()
