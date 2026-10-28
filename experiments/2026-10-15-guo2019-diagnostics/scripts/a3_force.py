"""A3 (часть 1): сила на шаре F(t), F(x), e(x) по monitor.csv."""
import os, sys, json
import numpy as np
from scipy.signal import savgol_filter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

R_BALL = 6e-3
cases = [c for c in ['1T_v65', '1T_v71', '2T_m01', '2T_m09'] if os.path.exists(os.path.join(DATA, c, 'monitor.csv'))]
out = {}
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for c in cases:
    m = monitor(c); sh = shot(c)
    t = m.t.values; vz = m.v_z.values; zc = m.z_c.values
    v0 = -vz[0]; E0 = 0.5 * MBALL * v0**2
    dt_s = np.median(np.diff(t))
    win = 41 if len(t) > 100 else 5
    vs = savgol_filter(vz, win, 3)
    F = MBALL * np.gradient(vs, t)           # >0: сила вверх на шар (тормозит)
    x = np.maximum(0.0, R_BALL - zc)          # перемещение шара после касания верха ткани (z=0)
    nb = (m.n_broken_v + m.n_broken_h).values
    ifb = int(np.argmax(nb > 0)) if (nb > 0).any() else len(t) - 1
    KE = 0.5 * MBALL * vz**2
    e_abs = E0 - KE                           # энергия, отданная шаром
    # e(x) интегралом F dx (трапеции по x монотонному)
    xm = np.maximum.accumulate(x)
    dx = np.diff(xm)
    e_int = np.concatenate([[0], np.cumsum(0.5 * (F[1:] + F[:-1]) * dx)])
    wf = m.w_f.values
    def at_w(w):
        i = int(np.argmax(wf >= w)) if (wf >= w).any() else None
        return None if i is None else dict(t_us=t[i] * 1e6, x_mm=x[i] * 1e3, e_abs=float(e_abs[i]), e_int=float(e_int[i]), share=float(e_abs[i] / E0), v=float(vz[i]))
    o = dict(V0=v0, E0=E0, n_samples=len(t), t_first_break_us=t[ifb] * 1e6 if nb[ifb] > 0 else None,
             x_fb_mm=x[ifb] * 1e3, wf_fb_mm=wf[ifb] * 1e3, e_fb=float(e_abs[ifb]), share_fb=float(e_abs[ifb] / E0),
             v_fb=float(vz[ifb]), e_fb_int=float(e_int[ifb]),
             F_max=float(F[:ifb + 1].max()), F_at_fb=float(F[ifb]), x_at_Fmax_mm=float(x[int(np.argmax(F[:ifb + 1]))] * 1e3),
             w25=at_w(25e-3), w27=at_w(27e-3), w_f_max_mm=float(wf.max() * 1e3),
             x_max_mm=float(x.max() * 1e3), e_end=float(e_abs[-1]), v_end=float(vz[-1]),
             ball_stop_t_us=float(t[np.argmax(vz > -0.02 * v0)] * 1e6) if (vz > -0.02 * v0).any() else None,
             n_broken_end=int(nb[-1]), outcome=sh.get('outcome'), e_shot=sh.get('e'), V_R=sh.get('V_R'),
             t_contact_us=float(t[np.argmax(zc <= R_BALL)] * 1e6))
    # жёсткость: F/x на участке 5-15 мм (линейный участок) и секущая при первом разрыве
    mk = (x > 5e-3) & (x < 15e-3) & (np.arange(len(t)) <= ifb)
    if mk.sum() > 3: o['k_5_15_N_per_mm'] = float(np.polyfit(x[mk] * 1e3, F[mk], 1)[0])
    out[c] = o
    lab = CASES[c]
    ax[0].plot(t * 1e6, F / 1e3, label=lab); ax[1].plot(x * 1e3, F / 1e3, label=lab); ax[2].plot(x * 1e3, e_abs, label=lab)
    ax[1].plot([x[ifb] * 1e3], [F[ifb] / 1e3], 'kx'); ax[2].plot([x[ifb] * 1e3], [e_abs[ifb]], 'kx')
ax[0].set_xlabel('t, мкс'); ax[0].set_ylabel('F, кН'); ax[0].set_title('Сила на шаре m dv/dt (сглаж. 41 отсчёт)')
ax[1].set_xlabel('перемещение шара x, мм'); ax[1].set_ylabel('F, кН'); ax[1].set_title('F(x), крестик = 1-й разрыв')
ax[2].set_xlabel('x, мм'); ax[2].set_ylabel('E0 - KE шара, Дж'); ax[2].set_title('e(x) (отданная шаром энергия)')
ax[2].axhline(20.2, color='gray', ls='--'); ax[2].text(1, 20.6, 'статья 1T: 20.2 Дж', fontsize=8)
ax[2].axvspan(25, 27, color='y', alpha=.25)
for a in ax: a.grid(alpha=.3); a.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'a3_force_energy.png'), dpi=130)
json.dump(out, open(os.path.join(OUT, 'data', 'a3_force.json'), 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
