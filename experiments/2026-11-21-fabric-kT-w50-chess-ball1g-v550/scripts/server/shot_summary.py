#!/usr/bin/env python3
"""Итог выстрела kT шахматной ткани по каталогу драйвера one_shot.py: печатает одну строку
   K MODE ВЕРДИКТ v_end z_min t_end nv nh n_stitch_broken ИСТОЧНИК
ВЕРДИКТ: П / НП / X. Берётся из result.txt; если там «не опр.» после продолжения до 1000 мкс (ошибка продолжения
серии w50: продолжение не ставит вердикт), -- по monitor.csv каталога _long тем же критерием, что у монитора:
  П  -- z_c < z_b - R, v_z < 0 и |Δv_z| < 0,1 м/с за 20 мкс;  НП -- v_z >= 0 20 мкс подряд.
v_end -- v_z шара в конце (у П -- минус V_R, у НП -- скорость отскока, >= 0), z_min -- наибольшее погружение центра.
Вызов: shot_summary.py КАТАЛОГ K
"""
import csv, glob, os, sys

D, A, R = 0.39e-3, 0.6 * 0.39e-3, 6.253e-3 / 2


def main():
    w, k = sys.argv[1], int(sys.argv[2])
    res = open(os.path.join(w, "result.txt"), encoding="utf-8").read().strip() if os.path.exists(
        os.path.join(w, "result.txt")) else ""
    dirs = sorted(glob.glob(os.path.join(w, "k%dT_*" % k, "*", "monitor.csv")))
    longs = [d for d in dirs if "_long" in d]
    mon = (longs or dirs or [None])[-1]
    rows = []
    if mon:
        with open(mon, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                rows.append(r)
    v = "X"
    src = "result.txt"
    if ": П (" in res:
        v = "П"
    elif ": НП (" in res:
        v = "НП"
    if not rows:
        print(k, "X нет monitor.csv")
        return
    t = [float(r["t"]) for r in rows]
    z = [float(r["z_c"]) for r in rows]
    vz = [float(r["v_z"]) for r in rows]
    if v == "X" and res:
        zb = -(2 * k * (2 * A + D) + (2 * k - 1) * 1e-6)
        src = "monitor.csv"
        for i in range(len(t)):
            j = i
            while j > 0 and t[i] - t[j - 1] <= 2e-5:
                j -= 1
            if t[i] - t[j] < 2e-5 * 0.99:
                continue
            if z[i] < zb - R and vz[i] < 0 and abs(vz[i] - vz[j]) < 0.1:
                v = "П"
                break
            if all(x >= 0 for x in vz[j:i + 1]):
                v = "НП"
                break
        else:
            i = len(t) - 1
        t_end = t[i]
        v_end = vz[i]
    else:
        t_end, v_end = t[-1], vz[-1]
    last = rows[-1]
    print("%d %s v_end=%.2f z_min=%.3e t_end=%.4e nv=%s nh=%s stitch_broken=%s src=%s" % (
        k, v, v_end, min(z), t_end, last["n_broken_v"], last["n_broken_h"], last.get("n_stitch_broken", "-"), src))


if __name__ == "__main__":
    main()
