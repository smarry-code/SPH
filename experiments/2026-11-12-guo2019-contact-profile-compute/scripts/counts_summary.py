#!/usr/bin/env python3
"""Сводка диагностических счётчиков (counts.csv из сборки с патчем SPH_PROF_COUNTS) по окнам шагов.

Колонки: step,time,ff_pairs,ff_prims,fs_prims,normal_ff,normal_fs,fric_ff_active,fric_fs_active (по каждому шагу).
normal_* -- активные нормальные слоты проекции (только P); fric_*_active -- активные слоты трения (контакт и mu > 0).
Окна по шагам -- как в bench_tables.py: до контакта 130-1300, удар 1700-2733 (до 25 мкс), поздний удар 2733-конец.

Запуск: python3 counts_summary.py <counts.csv> [<counts.csv> ...]
"""
import csv
import sys

WINDOWS = (("до контакта", 130, 1300), ("удар до 25 мкс", 1700, 2733), ("удар 25-40 мкс", 2733, 10 ** 9),
           ("весь прогон", 1, 10 ** 9))
COLS = ("ff_pairs", "ff_prims", "fs_prims", "normal_ff", "normal_fs", "fric_ff_active", "fric_fs_active")


def main():
    for path in sys.argv[1:]:
        rows = list(csv.DictReader(open(path)))
        print(f"== {path}: шагов {len(rows)}")
        print(f"{'окно':16s} {'шагов':>6s} " + " ".join(f"{c:>15s}" for c in COLS) + "  active/prims ff, fs")
        for name, lo, hi in WINDOWS:
            w = [r for r in rows if lo <= int(r["step"]) <= hi]
            if not w:
                continue
            mean = {c: sum(int(r[c]) for r in w) / len(w) for c in COLS}
            mx = {c: max(int(r[c]) for r in w) for c in COLS}
            act_ff = max(mean["normal_ff"], mean["fric_ff_active"])
            act_fs = max(mean["normal_fs"], mean["fric_fs_active"])
            rff = act_ff / mean["ff_prims"] if mean["ff_prims"] else 0.0
            rfs = act_fs / mean["fs_prims"] if mean["fs_prims"] else 0.0
            print(f"{name:16s} {len(w):6d} " + " ".join(f"{mean[c]:9.1f}/{mx[c]:<5d}" for c in COLS)
                  + f"  {rff:.4f}, {rfs:.4f}")
        print("(среднее/максимум; active = max(нормальные слоты, слоты трения))\n")


if __name__ == "__main__":
    main()
