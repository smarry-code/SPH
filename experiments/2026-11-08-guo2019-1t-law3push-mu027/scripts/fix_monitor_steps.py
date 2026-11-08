#!/usr/bin/env python3
"""Восстановить сквозной номер шага в monitor.csv выстрелов, склеенных до коммита с fix_monitor_steps.

После --restart монитор считает шаги с 1, а кадры fiber_output_<шаг> -- сквозные, и скрипты анимаций
(anim_2side.py, animate_1t.py) сопоставляли кадры после рестарта с последней строкой до него (t = 1500 мкс).
Исходный файл сохраняется как monitor.steps_raw.csv (если его ещё нет). Повторный запуск безопасен.

Usage: python fix_monitor_steps.py <каталог прогона с monitor.csv> ...
"""
import os, shutil, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import window_pinball as wp  # noqa: E402

for d in sys.argv[1:]:
    m = os.path.join(d, "monitor.csv")
    header, rows = wp._read_csv(m)
    n = wp.fix_monitor_steps(header, rows)
    if n:
        raw = os.path.join(d, "monitor.steps_raw.csv")
        if not os.path.exists(raw):
            shutil.copy2(m, raw)
        wp._write_csv(m, header, rows)
    print("%s: стыков %d, последний шаг %s" % (m, n, rows[-1][header.index("step")] if rows else "-"))
