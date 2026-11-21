"""Сводная таблица расчётов 1T с повреждённостью: исход, V_R, разрывы, первый разрыв, прогиб (к первому разрыву и
максимальный), сквозные проходы. Usage: python dmg_table.py label=shot.txt ..."""
import sys
print("| Расчёт | Исход | V_R, м/с | Разрывы в/г | 1-й разрыв, мкс | Прогиб к 1-му разрыву, мм | Макс. прогиб, мм | t_end, мкс | n_tunnel |")
print("|---|---|---|---|---|---|---|---|---|")
for a in sys.argv[1:]:
    lab, p = a.split("=", 1)
    kv = dict(t.split("=", 1) for t in open(p, encoding="utf-8").read().split() if "=" in t)
    f = lambda k, s=1.0, fmt="%.1f": (fmt % (float(kv[k]) * s)) if kv.get(k) not in (None, "nan") else "—"
    o = {"perforated": "пробит", "not_perforated": "шар остановлен", "undetermined": "не определён"}.get(kv["outcome"], kv["outcome"])
    print("| %s | %s | %s | %s/%s | %s | %s | %s | %s | %s |" % (
        lab, o, f("V_R"), kv["n_broken_v"], kv["n_broken_h"], f("t_first_break", 1e6, "%.0f"),
        f("w_f_at_first_break", 1e3), f("w_f_max", 1e3), f("t_end", 1e6, "%.0f"), kv.get("n_tunnel", "—")))
