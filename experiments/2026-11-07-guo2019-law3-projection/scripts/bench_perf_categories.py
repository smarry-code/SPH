#!/usr/bin/env python3
"""Разбивка профиля perf бенчмарка (§7.1) по категориям, в окнах до контакта и во время удара.

Доли каждой категории делятся на долю поиска соседей нитей: эта работа одинакова во всех режимах,
поэтому отношение не зависит от шума виртуальной машины (разброс абсолютного времени до 25 %).
Окна заданы в процентах длительности записи perf: 3–30 % — до контакта (шаги ~130–1300),
40–100 % — удар (шаги ~1700–4372). Запуск: python3 bench_perf_categories.py <каталог> [H L3 P]
"""
import re
import subprocess
import sys

CATS = [
    ("search", r"find_pairs|build_fiber_contact_search|SpatialIndex<3ul>::build\b"),
    ("proj", r"run_normal_projection|NormalWarmStart|compute_normal_counts|fill_friction_slots_from_normal|detail_normal"),
    ("contact", r"fiber_segment_closest_params|compute_fiber_contact_geometry|compute_fiber_contact_contributions"),
    ("solid", r"build_solid_sph_search|compute_solid_sph|SolidDelta|kernel_cubic"),
]
LINE = re.compile(r"\s*([\d.]+)%\s+\[[.k]\]\s+(.*)")


def main():
    root = sys.argv[1]
    variants = sys.argv[2:] or ["H", "L3", "P"]
    for v in variants:
        for w in ["3%-30%", "40%-100%"]:
            out = subprocess.run(
                ["perf", "report", "-i", f"{root}/{v}.perf/perf.data", "--time", w,
                 "--stdio", "--sort", "symbol", "--percent-limit", "0"],
                capture_output=True, text=True).stdout
            acc = {c: 0.0 for c, _ in CATS}
            acc["other"] = 0.0
            for line in out.splitlines():
                m = LINE.match(line)
                if not m:
                    continue
                p, sym = float(m.group(1)), m.group(2)
                for c, rx in CATS:
                    if re.search(rx, sym):
                        acc[c] += p
                        break
                else:
                    acc["other"] += p
            s = acc["search"]
            print(f"{v:3s} {w:9s} " + " ".join(f"{k}={x:5.1f}%" for k, x in acc.items())
                  + "  | /search: " + " ".join(f"{k}={x / s:.3f}" for k, x in acc.items() if k != "search")
                  + f"  step/search={sum(acc.values()) / s:.3f}")


if __name__ == "__main__":
    main()
