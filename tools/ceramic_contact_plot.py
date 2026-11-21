#!/usr/bin/env python3
"""Сводный график BC2 / BC6: кривые монитора при contact_force_type 0-3 и приборный эталон.

Читает каталоги прогонов bc2_shot / bc6_shot (<runs>/<кейс>_dx<мкм>_c<тип>/curve.txt, shot.txt)
и эталон cases/<кейс>/experiment/reference.txt. Формат вывода -- по расширению (.pdf, .svg, .png).

Совмещение эталона -- как у сводных SVG целей bc2_plot / bc6_plot:
  * BC6 -- эталон сдвинут на ref_shift_t10 из shot.txt базового типа 0 (совмещение по t10 фронта);
  * BC2 -- эталон на своей шкале времени (начало отсчёта не привязано к удару), сдвиг не вносится.
  --align=none|front перекрывает правило кейса.

Пример (нужен matplotlib):
  python3 tools/ceramic_contact_plot.py --case=bc6 --runs=runs --dx-um=200 \\
      --reference=cases/bc6/experiment/reference.txt --out=bc6_contact_types.pdf
"""
import argparse
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TYPE_NAME = {0: "выключен", 1: "пинбол", 2: "Свигл", 3: "импульс"}
# Цвет + штрих + толщина: кривые типов совпадают местами до 4-й значащей цифры.
TYPE_STYLE = {
    0: dict(color="#1f77b4", lw=2.2, ls="-"),
    1: dict(color="#ff7f0e", lw=1.8, ls=(0, (6, 4))),
    2: dict(color="#2ca02c", lw=1.4, ls="-"),
    3: dict(color="#d62728", lw=1.6, ls=(0, (1.5, 2.5))),
}
TITLE = {
    "bc2": "BC2: осевая скорость частицы-монитора (граница «керамика / окно»)",
    "bc6": "BC6: осевая скорость частицы-монитора (граница «керамика / окно»)",
}


def read_curve(path):
    t, v = [], []
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                t.append(float(parts[0]))
                v.append(float(parts[1]))
            except ValueError:
                continue
    return t, v


def read_kv(path):
    kv = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if "=" in line:
                k, val = line.rstrip("\n").split("=", 1)
                kv[k] = val
    return kv


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, choices=["bc2", "bc6"])
    ap.add_argument("--runs", required=True, help="корень --out утилиты")
    ap.add_argument("--dx-um", type=int, default=200, help="Δx прогонов, мкм (имя каталога)")
    ap.add_argument("--types", default="0,1,2,3")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--align", choices=["case", "none", "front"], default="case")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    runs = pathlib.Path(a.runs)
    types = [int(x) for x in a.types.split(",")]
    series, kv0, ntotal = [], None, None
    for c in types:
        d = runs / f"{a.case}_dx{a.dx_um}_c{c}"
        if not (d / "curve.txt").exists():
            print(f"нет прогона: {d}", file=sys.stderr)
            continue
        kv = read_kv(d / "shot.txt")
        if c == 0:
            kv0 = kv
        ntotal = kv.get("ntotal", ntotal)
        t, v = read_curve(d / "curve.txt")
        series.append((c, t, v, kv.get("outcome", "?")))
    if not series:
        print("ни одного прогона не найдено", file=sys.stderr)
        return 1

    rt, rv = read_curve(a.reference)
    align = a.align if a.align != "case" else ("front" if a.case == "bc6" else "none")
    shift = 0.0
    if align == "front":
        key = "ref_shift_t10" if a.case == "bc6" else "shift_t50"
        if kv0 is None or key not in kv0:
            print(f"совмещение: нет {key} в shot.txt типа 0, эталон без сдвига", file=sys.stderr)
        else:
            shift = float(kv0[key]) if a.case == "bc6" else -float(kv0[key])

    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.plot([(x + shift) * 1e6 for x in rt], rv, color="black", lw=2.6, alpha=0.55,
            label="эталон (прибор)" + (f", сдвиг {shift * 1e6:+.2f} мкс".replace(".", ",").replace("-", "−") if shift else ""))
    for c, t, v, outcome in series:
        lab = f"тип {c} — {TYPE_NAME.get(c, '?')}"
        if outcome != "completed":
            lab += f" ({outcome})"
        ax.plot([x * 1e6 for x in t], v, label=lab, **TYPE_STYLE.get(c, {}))
        if outcome != "completed" and t:
            ax.plot([t[-1] * 1e6], [v[-1]], "x", color=TYPE_STYLE.get(c, {}).get("color", "k"), ms=9, mew=2)

    ax.set_xlabel("время, мкс")
    ax.set_ylabel("осевая скорость, м/с")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    sub = f"Δx = {a.dx_um / 1000:g} мм".replace(".", ",")
    if ntotal:
        n = int(ntotal)
        word = "частица" if n % 10 == 1 and n % 100 != 11 else (
            "частицы" if 2 <= n % 10 <= 4 and not 12 <= n % 100 <= 14 else "частиц")
        sub += ", " + f"{n:,}".replace(",", "\u00a0") + " " + word
    if a.case == "bc2" and not shift:
        sub += "; эталон на своей шкале времени (не привязан к удару)"
    elif shift:
        sub += "; эталон совмещён по фронту"
    ax.set_title(TITLE[a.case] + "\n" + sub, fontsize=11)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(a.out)
    print(f"записан {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
