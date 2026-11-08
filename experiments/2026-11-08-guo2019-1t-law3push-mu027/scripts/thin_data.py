#!/usr/bin/env python3
"""Сводные данные каталога эксперимента (data/) из рабочих каталогов серии на claude-compute-8t.

Вход: корень со скачанными рабочими каталогами серий (как на сервере):
  <src>/vbl1T_law3push_mu027_det/       -- первая серия, V0 = 71 / 82 (задача mu027_window)
  <src>/vbl1T_law3push_mu027_v69_70/    -- вторая серия, V0 = 69 / 70 (задача mu027_v69_70)
в каждом: queue.log, window.md, mu0.27_v<V>.log, mu0.27_v<V>[_long]/shots.csv, .../1T_base_v<V>.0_r0/shot.txt,
в _long также monitor.csv, energy_log.csv (полная история 0..t_end: база + продолжение из чекпоинта).

Выход (data/):
  shot_v<V>_base.txt, shot_v<V>_long.txt -- сводки выстрела (long -- после склейки монитора драйвером)
  monitor_v<V>.csv  -- монитор, каждая --mon-every строка (строка монитора -- каждые 10 шагов) + последняя;
                       только столбцы MON_COLS
  energy_v<V>.csv   -- energy_log.csv, каждая --en-every строка (строка лога -- каждые 100 шагов) + последняя
  log_v<V>.txt, log_v<V>_long.txt -- лог guo2019_shot; queue_<серия>.log, window_<серия>.md -- журнал и итог драйвера
  guo2019_fabric_layout.txt -- раскладка нитей (одинакова во всех выстрелах)

Usage: python thin_data.py <src> <data> [--mon-every 50] [--en-every 10]
"""
import argparse, csv, glob, os, shutil

SERIES = {"vbl1T_law3push_mu027_det": "71_82", "vbl1T_law3push_mu027_v69_70": "69_70"}
MON_COLS = ["step", "t", "z_c", "v_z", "KE_s", "KE_f", "U_b", "E_tot", "Ff_z", "n_bonds", "n_broken_v", "n_broken_h",
            "w_f"]


def thin(src, dst, every, cols=None):
    with open(src, newline="") as f:
        rows = list(csv.reader(f))
    head, body = rows[0], [r for r in rows[1:] if len(r) == len(rows[0])]
    keep = body[::every]
    if body and keep[-1] is not body[-1]:
        keep.append(body[-1])
    idx = [head.index(c) for c in cols] if cols else list(range(len(head)))
    with open(dst, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow([head[i] for i in idx])
        for r in keep:
            w.writerow([r[i] for i in idx])
    return len(body), len(keep)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("data")
    ap.add_argument("--mon-every", type=int, default=50)
    ap.add_argument("--en-every", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(a.data, exist_ok=True)
    for ser, tag in SERIES.items():
        w = os.path.join(a.src, ser)
        shutil.copy(os.path.join(w, "queue.log"), os.path.join(a.data, "queue_%s.log" % tag))
        shutil.copy(os.path.join(w, "window.md"), os.path.join(a.data, "window_%s.md" % tag))
        for d in sorted(glob.glob(os.path.join(w, "mu0.27_v*"))):
            if not os.path.isdir(d):
                continue
            name = os.path.basename(d)              # mu0.27_v69 / mu0.27_v69_long
            v = name.split("_v")[1]
            kind = "long" if v.endswith("_long") else "base"
            v = v.replace("_long", "")
            run = glob.glob(os.path.join(d, "1T_*", ""))[0]
            shutil.copy(os.path.join(run, "shot.txt"), os.path.join(a.data, "shot_v%s_%s.txt" % (v, kind)))
            shutil.copy(d + ".log", os.path.join(a.data, "log_v%s%s.txt" % (v, "_long" if kind == "long" else "")))
            lay = os.path.join(run, "guo2019_fabric_layout.txt")
            if os.path.exists(lay) and not os.path.exists(os.path.join(a.data, "guo2019_fabric_layout.txt")):
                shutil.copy(lay, a.data)
            if kind == "long":
                n, k = thin(os.path.join(run, "monitor.csv"), os.path.join(a.data, "monitor_v%s.csv" % v),
                            a.mon_every, MON_COLS)
                m, j = thin(os.path.join(run, "energy_log.csv"), os.path.join(a.data, "energy_v%s.csv" % v),
                            a.en_every)
                print("V0=%s: monitor %d -> %d строк, energy_log %d -> %d" % (v, n, k, m, j))


if __name__ == "__main__":
    main()
