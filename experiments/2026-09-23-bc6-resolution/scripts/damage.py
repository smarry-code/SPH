# -*- coding: utf-8 -*-
# Разбор повреждения JH-2 в образце B4C по кадрам прогона.
# Запуск:  python scripts/damage.py <каталог прогона с sph_output_*.txt>
# Каталог получается ключом --keep-workdir (см. REPRO.md). Кадры в
# эксперимент не входят (752 МБ); вывод скрипта сохранён в data/damage_dx020.txt.
import sys, glob, os
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
W = sys.argv[1]
# Основной монитор при Δx = 0.2 мм. В логе теста он назван #11440 — это индекс
# с нуля; колонка N кадра нумерует с единицы, поэтому в кадре это N = 11441.
# Соседняя по оси частица на слой глубже — N = 11440; она выводится тоже.
MONITOR = 11441
INNER   = 11440
SAVE_T = 0.05            # мкс, save_t из tests/data/bc6_initial_conditions.txt
fs = sorted(glob.glob(os.path.join(W, "sph_output_*.txt")), key=lambda p: int(p[-9:-4]))
print("кадров:", len(fs), "; время кадра k принимается равным k·save_t (с точностью до шага)")
print("\n1. Счёт разрушенных частиц B4C по времени")
print("   «у монитора» — r < 2 мм и не далее 2 мм от тыльной грани образца")
for k, p in enumerate(fs):
    if k % 2: continue
    a = np.loadtxt(p, skiprows=1, usecols=(0, 1, 3, 4, 8, 13))
    N, mat, r, z, pr, D = a.T
    s = mat == 10; zb = z[s].max()
    near = s & (r < 2e-3) & (z > zb - 2e-3)
    Dm = D[N == MONITOR][0]; Di = D[N == INNER][0]
    print(f"   t≈{k*SAVE_T:4.2f} мкс  всего B4C {int(s.sum())}: D>0 {int((D[s]>0).sum()):5d}, "
          f"D≥0.99 {int((D[s]>=0.99).sum()):5d} | у монитора D>0 {int((D[near]>0).sum()):3d} из {int(near.sum())} | "
          f"D монитора {Dm:.2f}, на слой глубже {Di:.2f}")
print("\n2. Доля D ≥ 0.99 в B4C: строки — глубина от плоскости удара, столбцы — радиус")
for k in (4, 20, 22, 26, 30, 40, 50, len(fs) - 1):
    a = np.loadtxt(fs[k], skiprows=1, usecols=(1, 3, 4, 13))
    mat, r, z, D = a.T; s = mat == 10
    zs, rs, Ds = z[s], r[s], D[s]; z0 = zs.min()
    zb = np.linspace(z0, zs.max() + 1e-9, 6); rb = [0, 2e-3, 10e-3, 20e-3, 30e-3, 40e-3]
    print(f"   t≈{k*SAVE_T:.2f} мкс   r, мм:     0-2  2-10 10-20 20-30 30-38")
    for i in range(5):
        row = []
        for j in range(5):
            m = (zs >= zb[i]) & (zs < zb[i+1]) & (rs >= rb[j]) & (rs < rb[j+1])
            row.append(f"{(Ds[m] >= 0.99).mean():5.2f}" if m.any() else "   — ")
        print(f"     z {(zb[i]-z0)*1e3:5.1f}–{(zb[i+1]-z0)*1e3:5.1f} мм: " + " ".join(row))
