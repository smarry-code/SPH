#!/usr/bin/env python3
"""Графики и метрики по monitor.csv прогона Guo 2019 (1T, level scaled).

Колонки: step,t,dt,z_c,v_z,KE_s,KE_p,KE_f,U_b,U_int,E_tot,Fs_x..Ff_z,n_bonds,
n_broken_v,n_broken_h,w_f,w_p.

Выход (SVG, подписи по-русски):
    v_z.svg                -- v_z(t) шарика, м/с, от t, мкс;
    energy.svg             -- KE_s, KE_f, U_b, E_tot, делённые на E_tot(0);
    deflection_breaks.svg  -- прогиб w_f(t), мм, и число порванных нитей (вертикальные/горизонтальные);
    reaction.svg           -- реакция заделки Ff_z(t), Н.
Метрики печатаются в stdout и пишутся в --metrics.
--thin=N пишет прореженную копию monitor.csv (каждая N-я строка + последняя) для data/.

Usage:
    python plot_monitor.py --run-dir out/1T_scaled_v76.0_r0 --out-dir figs \\
        --metrics figs/metrics.txt --thin 20 --thin-out data/monitor_thin.csv
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"  # текст остаётся текстом, файлы компактнее


def read_monitor(path: str) -> dict[str, np.ndarray]:
    """monitor.csv; недописанная последняя строка (прогон ещё идёт) отбрасывается."""
    with open(path) as f:
        names = f.readline().strip().split(",")
        rows = [ln.split(",") for ln in f]
    rows = [r for r in rows if len(r) == len(names)]
    d = np.array(rows, dtype=float).reshape(-1, len(names))
    return {n: d[:, i] for i, n in enumerate(names)}


def thin_csv(src: str, dst: str, n: int) -> None:
    """Каждая n-я строка данных + последняя; заголовок сохраняется как есть."""
    with open(src) as f:
        lines = f.read().splitlines()
    head, rows = lines[0], [ln for ln in lines[1:] if ln.strip()]
    keep = rows[::max(n, 1)]
    if rows and keep[-1] is not rows[-1]:
        keep.append(rows[-1])
    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
    with open(dst, "w") as f:
        f.write(head + "\n" + "\n".join(keep) + "\n")
    print(f"thin: {len(keep)} из {len(rows)} строк -> {dst} ({os.path.getsize(dst) / 1e3:.0f} КБ)")


def metrics(m: dict[str, np.ndarray]) -> list[str]:
    t_us = m["t"] * 1e6
    nb = m["n_broken_v"] + m["n_broken_h"]
    E0, E = m["E_tot"][0], m["E_tot"]
    out = [
        f"строк монитора            : {len(t_us)}",
        f"t_конец, мкс              : {t_us[-1]:.3f}",
        f"max w_f, мм               : {m['w_f'].max() * 1e3:.4f} (при t = {t_us[np.argmax(m['w_f'])]:.3f} мкс)",
    ]
    hit = np.where(nb > 0)[0]
    out.append("первый разрыв нити, мкс   : " + (f"{t_us[hit[0]]:.3f} (шаг {int(m['step'][hit[0]])})" if len(hit)
                                                  else "разрывов нет"))
    out += [
        f"порвано нитей верт./гориз. (кон.): {int(m['n_broken_v'][-1])}/{int(m['n_broken_h'][-1])}",
        f"v_z в конце, м/с          : {m['v_z'][-1]:.3f}",
        f"min v_z, м/с              : {m['v_z'].min():.3f}",
        f"max v_z, м/с              : {m['v_z'].max():.3f}",
        f"E_tot(0), Дж              : {E0:.6g}",
        f"E_tot_max / E_tot_0       : {E.max() / E0:.6f}",
        f"E_tot_end / E_tot_0       : {E[-1] / E0:.6f}",
        f"max |Ff_z|, Н             : {np.abs(m['Ff_z']).max():.4g}",
    ]
    return out


def style(ax):
    ax.grid(alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", help="каталог прогона с monitor.csv")
    ap.add_argument("--monitor", help="путь к monitor.csv (вместо --run-dir), в т.ч. прореженному")
    ap.add_argument("--out-dir", default=".", help="каталог для SVG")
    ap.add_argument("--metrics", default=None, help="файл для текста метрик")
    ap.add_argument("--thin", type=int, default=0, help="N: записать прореженную копию monitor.csv")
    ap.add_argument("--thin-out", default=None, help="путь прореженной копии (по умолчанию <out-dir>/monitor_thin.csv)")
    args = ap.parse_args()

    src = args.monitor or (os.path.join(args.run_dir, "monitor.csv") if args.run_dir else None)
    if not src or not os.path.exists(src):
        print("нужен --run-dir или --monitor с существующим monitor.csv", file=sys.stderr)
        sys.exit(1)
    os.makedirs(args.out_dir, exist_ok=True)
    m = read_monitor(src)
    t = m["t"] * 1e6

    # (a) v_z
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, m["v_z"], color="tab:blue", lw=1.3)
    ax.set_xlabel("t, мкс"); ax.set_ylabel("v_z шарика, м/с")
    ax.set_title("Скорость шарика"); style(ax)
    fig.tight_layout(); fig.savefig(os.path.join(args.out_dir, "v_z.svg")); plt.close(fig)

    # (b) энергии
    E0 = m["E_tot"][0]
    fig, ax = plt.subplots(figsize=(7, 4))
    for key, lab, c in (("KE_s", "KE шарика", "tab:blue"), ("KE_f", "KE ткани", "tab:green"),
                        ("U_b", "U связей", "tab:orange"), ("E_tot", "E полная", "k")):
        ax.plot(t, m[key] / E0, label=lab, color=c, lw=1.6 if key == "E_tot" else 1.2)
    ax.set_xlabel("t, мкс"); ax.set_ylabel("E / E_tot(0)")
    ax.set_title("Баланс энергии"); ax.legend(fontsize=8); style(ax)
    fig.tight_layout(); fig.savefig(os.path.join(args.out_dir, "energy.svg")); plt.close(fig)

    # (c) прогиб и разрывы
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True)
    a1.plot(t, m["w_f"] * 1e3, color="tab:red", lw=1.3)
    a1.set_ylabel("прогиб w_f, мм"); a1.set_title("Прогиб ткани и разрывы нитей"); style(a1)
    a2.plot(t, m["n_broken_v"], label="вертикальные (вдоль y, зажаты)", color="tab:purple", lw=1.3)
    a2.plot(t, m["n_broken_h"], label="горизонтальные (вдоль x)", color="tab:brown", lw=1.3, ls="--")
    a2.set_xlabel("t, мкс"); a2.set_ylabel("разорванные нити, шт."); a2.legend(fontsize=8); style(a2)
    fig.tight_layout(); fig.savefig(os.path.join(args.out_dir, "deflection_breaks.svg")); plt.close(fig)

    # (d) реакция заделки
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, m["Ff_z"], color="tab:blue", lw=1.0)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("t, мкс"); ax.set_ylabel("Ff_z, Н")
    ax.set_title("Сила на заделке (нити), z-компонента"); style(ax)
    fig.tight_layout(); fig.savefig(os.path.join(args.out_dir, "reaction.svg")); plt.close(fig)

    lines = metrics(m)
    print("\n".join(lines))
    if args.metrics:
        os.makedirs(os.path.dirname(os.path.abspath(args.metrics)), exist_ok=True)
        with open(args.metrics, "w", encoding="utf-8") as f:
            f.write(f"# метрики из {src}\n" + "\n".join(lines) + "\n")
        print(f"метрики -> {args.metrics}")
    if args.thin > 0:
        thin_csv(src, args.thin_out or os.path.join(args.out_dir, "monitor_thin.csv"), args.thin)


if __name__ == "__main__":
    main()
