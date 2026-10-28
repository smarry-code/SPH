"""A1: геометрия ткани по guo2019_fabric.fib (1T base), без расчёта."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

SCR = os.path.dirname(DATA)
fib = os.path.join(DATA, '1T_v65', 'guo2019_fabric.fib')
lay_f = os.path.join(DATA, '1T_v65', 'guo2019_fabric_layout.txt')
if not os.path.exists(fib): fib = os.path.join(SCR, 'fib1T', 'f.fib'); lay_f = os.path.join(SCR, 'fib1T', 'l.txt')
Y = read_fib(fib)
dirs = [l.split()[1] for l in open(lay_f) if not l.startswith('#')]
p = 1.18e-3; A = 0.234e-3; d = D
WIN = 0.0600   # окно = зажимы, |x|,|y| <= 60 мм
out = {}
for dr in 'xy':
    ax = 0 if dr == 'x' else 1
    sumL = sumS = Lwin = Swin = 0.0; zs = []; n = 0
    for y, dd in zip(Y, dirs):
        if dd != dr: continue
        seg = np.linalg.norm(np.diff(y, axis=0), axis=1)
        mid = 0.5 * (y[1:, ax] + y[:-1, ax])
        inw = np.abs(mid) <= WIN
        Lwin += seg[inw].sum(); Swin += np.abs(np.diff(y[:, ax]))[inw].sum()
        sumL += seg.sum(); sumS += abs(y[-1, ax] - y[0, ax]); zs.append(y[:, 2]); n += 1
    zs = np.concatenate(zs)
    out[dr] = dict(n=n, crimp_win=Lwin / Swin - 1, crimp_all=sumL / sumS - 1, Lwin=Lwin, Swin=Swin,
                   zmin=zs.min(), zmax=zs.max())
mass_win = mass_all = 0.0
for y in Y:
    seg = np.linalg.norm(np.diff(y, axis=0), axis=1)
    mid = 0.5 * (y[1:, :2] + y[:-1, :2])
    inw = (np.abs(mid) <= WIN).all(axis=1)
    mass_win += RHO * AREA * seg[inw].sum(); mass_all += RHO * AREA * seg.sum()
out['mass_win_kg'] = mass_win; out['mass_all_kg'] = mass_all
out['rho_s_win'] = mass_win / (2 * WIN) ** 2
out['rho_s_all_over_win_area'] = mass_all / (2 * WIN) ** 2
allz = np.concatenate([y[:, 2] for y in Y])
out['thick_total'] = float(allz.max() - allz.min() + d)
zx = np.concatenate([y[:, 2] for y, dd in zip(Y, dirs) if dd == 'x'])
zy = np.concatenate([y[:, 2] for y, dd in zip(Y, dirs) if dd == 'y'])
out['z_range_x'] = (zx.min(), zx.max()); out['z_range_y'] = (zy.min(), zy.max())
out['thick_x_swing_plus_d'] = float(zx.max() - zx.min() + d)


def zig(Aa, p2=p / 2): return np.hypot(p2, Aa) / p2 - 1


xx = np.linspace(0, p, 200001)


def cosarc(Aa):
    z = Aa * np.cos(np.pi * xx / p); return np.sum(np.hypot(np.diff(xx), np.diff(z))) / p - 1


out['analytic'] = {nm: dict(zigzag=zig(a), cosine=cosarc(a)) for nm, a in
                   [('A=0.6d (0.234 mm)', 0.6 * d), ('A=0.5d (0.195 mm)', 0.5 * d), ('A=0.4d', 0.4 * d), ('A=0.25d', 0.25 * d)]}
tex = RHO * AREA * 1e6
rho_flat = 2 * (1 / p) * RHO * AREA
cm = out['x']['crimp_win']
out.update(tex=tex, yarns_per_cm=1 / (p * 100), rho_s_noCrimp=rho_flat, rho_s_model_formula=rho_flat * (1 + cm),
           crimp_for_280=0.280 / rho_flat - 1,
           pitch_for_280_mm=2 * RHO * AREA * (1 + cm) / 0.280 * 1e3,
           pitch_for_280_flat_mm=2 * RHO * AREA / 0.280 * 1e3,
           tex_for_280=0.280 * p / 2 / (1 + cm) * 1e6, tex_for_280_nocrimp=0.280 * p / 2 * 1e6,
           min_thick_round=2 * d, A_min_round=d / 2, gap_cross_model=2 * A - d,
           flat_width_for_t02=AREA / (np.pi / 4 * 0.2e-3))
json.dump(out, open(os.path.join(OUT, 'data', 'a1_geometry.json'), 'w'), indent=1, default=float)
for k, v in out.items(): print(k, v)

fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
y = Y[50]; m = np.abs(y[:, 0]) < 4e-3
ax[0].plot(y[m, 0] * 1e3, y[m, 2] * 1e3, 'o-', label='x-нить (по x)')
yy = [y2 for y2, dd in zip(Y, dirs) if dd == 'y'][50]
m2 = np.abs(yy[:, 1]) < 4e-3
ax[0].plot(yy[m2, 1] * 1e3, yy[m2, 2] * 1e3, 's--', label='y-нить (по y)')
ax[0].set_xlabel('координата вдоль нити, мм'); ax[0].set_ylabel('z, мм'); ax[0].set_title('Узлы нитей 1T base (k=2)')
ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
Aa = np.linspace(0.1, 0.8, 100) * d
ax[1].plot(Aa / d, [zig(a) * 100 for a in Aa], label='ломаная k=2')
ax[1].plot(Aa / d, [cosarc(a) * 100 for a in Aa], label='косинус')
ax[1].axvline(0.5, color='r', ls=':', label='A=d/2 (мин. для круглой нити)')
ax[1].axvline(0.6, color='k', ls=':', label='A=0.6d (модель)')
ax[1].axhline(3.0, color='g', ls='--', label='eps_break 3 %')
ax[1].set_xlabel('A/d'); ax[1].set_ylabel('излишек длины, %'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'a1_geometry.png'), dpi=130)
