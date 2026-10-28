"""A3 (часть 2): поперечная волна, профили eps, извитость, места разрыва, втягивание. Вход: data/<case>_frames.pkl, a3_force.json"""
import os, sys, json, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

THR = 0.2e-3
C0 = np.sqrt(E_Y / RHO)
force = json.load(open(os.path.join(OUT, 'data', 'a3_force.json')))


def smith(V, c=C0):
    """Смит: деформация eps в нити при ударе скоростью V и лабораторная скорость поперечной волны."""
    from scipy.optimize import brentq
    V = max(V, 0.5)
    f = lambda e: c * np.sqrt(2 * e * np.sqrt(e * (1 + e)) - e**2) - V
    e = brentq(f, 1e-9, 0.5)
    U = c * (np.sqrt(e * (1 + e)) - e)
    return e, U


def front(s, dzs, sign, thr=THR, gap=8):
    """Положение фронта прогиба: идём от центра наружу; фронт -- последняя точка выше порога перед gap подряд ниже порога."""
    m = (s * sign > 0)
    ss = s[m] * sign; v = np.abs(dzs[m])
    o = np.argsort(ss); ss = ss[o]; v = v[o]
    below = 0; last = 0.0
    for sv, vv in zip(ss, v):
        if vv > thr: last = sv; below = 0
        else:
            below += 1
            if below >= gap and last > 0: break
    return last


def crimp_removed_and_share(p, p0):
    return None


summary = {}
cases = [c for c in ['1T_v65', '1T_v71', '2T_m01', '2T_m09'] if os.path.exists(os.path.join(OUT, 'data', c + '_frames.pkl'))]
figF, axF = plt.subplots(1, 2, figsize=(13, 4.6))
figP, axP = plt.subplots(1, 2, figsize=(13, 4.2))
for c in cases:
    R = pickle.load(open(os.path.join(OUT, 'data', c + '_frames.pkl'), 'rb'))
    fr = R['frames']; nply = max(l['ply'] for l in R['lay'])
    t = np.array([f['t'] for f in fr]) * 1e6
    fx = np.zeros(len(fr)); fy = np.zeros(len(fr))
    for i, f in enumerate(fr):
        px = f['prof'][(1, 'x', 50)]; py = f['prof'][(1, 'y', 50)]
        fx[i] = 0.5 * (front(px['s'], px['dzs'], 1) + front(px['s'], px['dzs'], -1))
        fy[i] = 0.5 * (front(py['s'], py['dzs'], 1) + front(py['s'], py['dzs'], -1))
    # скорость шара по monitor
    m = monitor(c); V0 = -m.v_z.values[0]
    vball = np.interp(t * 1e-6, m.t.values, -m.v_z.values)
    # фронт по Смиту (интеграл U(V(t)) от касания), нить и ткань c/sqrt2
    tc = force[c]['t_contact_us']
    ts = np.linspace(tc, t.max(), 400)
    vb = np.interp(ts * 1e-6, m.t.values, -m.v_z.values)
    U_y = np.array([smith(v)[1] for v in vb]); U_f = np.array([smith(v, C0 / np.sqrt(2))[1] for v in vb])
    S_y = np.concatenate([[0], np.cumsum(0.5 * (U_y[1:] + U_y[:-1]) * np.diff(ts) * 1e-6)])
    S_f = np.concatenate([[0], np.cumsum(0.5 * (U_f[1:] + U_f[:-1]) * np.diff(ts) * 1e-6)])
    # оценка "щелевой" волны: прогиб w = int V dt, front = w / sqrt(2c+c^2), c = 7.58 %
    cr = 0.0757786
    w_ball = np.concatenate([[0], np.cumsum(0.5 * (vb[1:] + vb[:-1]) * np.diff(ts) * 1e-6)])
    S_slack = w_ball / np.sqrt(2 * cr + cr**2)
    # линейные подгонки скорости фронта
    def fit(tt, ff, lo, hi):
        k = (tt >= lo) & (tt <= hi) & (ff < 0.055)
        return np.polyfit(tt[k] * 1e-6, ff[k], 1)[0] if k.sum() >= 3 else np.nan
    sm = dict(V0=V0, c0=C0)
    for V in (V0, 0.9 * V0):
        e, U = smith(V); sm['smith_yarn_V%.0f' % V] = dict(eps=e, U=U)
    e, U = smith(V0, C0 / np.sqrt(2)); sm['smith_fabric_eps'] = e; sm['smith_fabric_U'] = U
    sm['U_slack_est'] = V0 / np.sqrt(2 * cr + cr**2)
    sm['front_speed_x_40_120'] = fit(t, fx, 40, 120); sm['front_speed_y_40_120'] = fit(t, fy, 40, 120)
    sm['front_speed_x_120_220'] = fit(t, fx, 120, 220); sm['front_speed_y_120_220'] = fit(t, fy, 120, 220)
    tt_clamp = None
    for arr, nm in ((fy, 'y'), (fx, 'x')):
        k = np.where(arr > 0.055)[0]
        sm['t_front_reaches_clamp_%s_us' % nm] = float(t[k[0]]) if len(k) else None
    # Смит: время достижения 60 мм
    k = np.where(S_y > 0.06)[0]; sm['smith_yarn_t60_us'] = float(ts[k[0]] - tc) if len(k) else None
    k = np.where(S_f > 0.06)[0]; sm['smith_fabric_t60_us'] = float(ts[k[0]] - tc) if len(k) else None
    sm['fx_us'] = fx.tolist(); sm['fy_us'] = fy.tolist(); sm['t_us'] = t.tolist()
    # ---- втягивание
    pm = np.array([[f['pull'][:, 2:].mean(), f['pull'][:, 2:].max(), f['pull'][np.abs(f['pull'][:, 1] - 50.5) < 12][:, 2:].mean()] for f in fr])
    sm['pull_mean_end_mm'] = (pm[-1, 0] * 1e3); sm['pull_max_end_mm'] = pm[-1, 1] * 1e3
    tfb = force[c]['t_first_break_us']
    ifb = int(np.max(np.where(t < tfb)[0])) if tfb else len(t) - 1       # последний кадр до первого разрыва
    sm['frame_before_break_us'] = float(t[ifb]); sm['pull_mean_at_fb_mm'] = pm[ifb, 0] * 1e3; sm['pull_max_at_fb_mm'] = pm[ifb, 1] * 1e3
    sm['pull_center_at_fb_mm'] = pm[ifb, 2] * 1e3
    lab = CASES[c]
    l1, = axF[0].plot(t, fy * 1e3, 'o-', label=lab + ' y')
    axF[0].plot(t, fx * 1e3, 's--', color=l1.get_color(), alpha=.6, label=lab + ' x')
    if c == '1T_v65':
        axF[0].plot(ts, S_y * 1e3, 'k-', lw=1, label='Смит, нить c=sqrt(E/rho)')
        axF[0].plot(ts, S_f * 1e3, 'k--', lw=1, label='Смит, ткань c/sqrt2')
        axF[0].plot(ts, S_slack * 1e3, 'k:', lw=1.2, label='геометрия: w/sqrt(2c+c^2), c=7.6 %')
    axP[0].plot(t, pm[:, 0] * 1e3, label=lab + ' среднее'); axP[0].plot(t, pm[:, 1] * 1e3, '--', color=axP[0].lines[-1].get_color(), label=lab + ' макс')
    # ---- профили eps к кадру до разрыва
    f = fr[ifb]; f0 = fr[0]
    prof_stats = {}
    for ply in range(1, nply + 1):
        for dr in 'xy':
            for idx in (50, 51):
                k = (ply, dr, idx); p = f['prof'][k]; p0 = f0['prof'][k]
                e = p['eps']; ok = ~np.isnan(e)
                Da = p['La'] - p['La0']; Ds = p['Ls'] - p['Ls0']
                prof_stats['%d%s%d' % (ply, dr, idx)] = dict(
                    eps_mean=float(np.nanmean(e)), eps_max=float(np.nanmax(e)), s_of_max_mm=float(p['s'][:-1][np.nanargmax(e)] * 1e3),
                    r_of_max_mm=float(np.hypot(p['s'][:-1][np.nanargmax(e)], p['other'][:-1][np.nanargmax(e)]) * 1e3),
                    epsp_max=float(np.nanmax(p['epp'])), frac_len_eps_gt_1pct=float(np.mean(e[ok] > 0.01)),
                    dLa_mm=Da * 1e3, dLs_mm=Ds * 1e3, share_straightening=float((Ds - Da) / Ds) if Ds > 1e-9 else None,
                    La0_mm=p['La0'] * 1e3, crimp_t=float(p['La'] / p['Ls'] - 1), crimp_0=float(p0['La'] / p0['Ls'] - 1),
                    nbonds=p['nbonds'], nbeads=p['nbeads'])
    sm['profiles_fb'] = prof_stats
    # ---- где рвётся: связи, исчезнувшие в первых кадрах с потерями
    gone_frames = [(f_['t'] * 1e6, f_['gone']) for f_ in fr if f_.get('gone')]
    gl = []
    for tg, g in gone_frames[:2]:
        for b in g:
            gl.append(dict(t_frame_us=tg, ply=b['ply'], dir=b['dir'], idx=b['idx'], pos_in_yarn=b['pos_in_yarn'],
                           x_mm=b['xyz_prev'][0] * 1e3, y_mm=b['xyz_prev'][1] * 1e3, z_mm=b['xyz_prev'][2] * 1e3,
                           r_mm=float(np.hypot(*b['xyz_prev'][:2]) * 1e3), eps_prev=b['eps_prev'], epsp_prev=b['epsp_prev'], t_prev_us=b['t_prev'] * 1e6))
    sm['first_breaks'] = gl
    # график eps(s) для центральных y- и x-нитей
    for j, dr in enumerate('xy'):
        pr = f['prof'][(1, dr, 50)]
        s = pr['s'][:-1] + 0
        axP[1].plot(s * 1e3 if dr == 'y' else np.nan * s, pr['eps'] * 100, color=axP[1].lines[-1].get_color() if False else None, label=(lab + ' y-нить 50, t=%.0f мкс' % t[ifb]) if dr == 'y' else None)
    summary[c] = sm
    # индивидуальные графики по кейсу: eps, dz, локальная извитость
    fg, ag = plt.subplots(1, 3, figsize=(16, 4.2))
    for dr, st in (('y', '-'), ('x', '--')):
        for ti in sorted(set([max(1, ifb - 6), max(1, ifb - 3), ifb])):
            pr = fr[ti]['prof'][(1, dr, 50)]
            ag[0].plot(pr['s'][:-1] * 1e3, pr['eps'] * 100, st, label='%s-нить, t=%.0f мкс' % (dr, fr[ti]['t'] * 1e6))
            ag[1].plot(pr['s'] * 1e3, pr['dzs'] * 1e3, st)
            lc = pr['loc']
            if len(lc): ag[2].plot(lc[:, 0] * 1e3, lc[:, 2] * 100, st)
    ag[0].set_title(lab + ': eps связей вдоль центральной нити'); ag[0].set_xlabel('координата вдоль нити, мм'); ag[0].set_ylabel('eps, %'); ag[0].legend(fontsize=7)
    ag[1].set_title('прогиб нити (сглаж. по периоду)'); ag[1].set_xlabel('координата, мм'); ag[1].set_ylabel('dz, мм')
    ag[2].set_title('локальная извитость (длина/хорда-1 на периоде)'); ag[2].set_xlabel('координата, мм'); ag[2].set_ylabel('%')
    ag[2].axhline(7.58, color='gray', ls=':')
    for a in ag: a.grid(alpha=.3)
    plt.tight_layout(); fg.savefig(os.path.join(FIG, 'a3_profiles_%s.png' % c), dpi=120); plt.close(fg)
axF[0].set_xlabel('t, мкс'); axF[0].set_ylabel('фронт прогиба (|dz|>0.2 мм), мм'); axF[0].set_title('Фронт поперечной волны (вдоль осей через центр)')
axF[0].set_ylim(0, 65); axF[0].grid(alpha=.3); axF[0].legend(fontsize=6)
axP[0].set_xlabel('t, мкс'); axP[0].set_ylabel('втягивание свободных концов x-нитей, мм'); axP[0].grid(alpha=.3); axP[0].legend(fontsize=7)
axP[1].axis('off')
axF[1].axis('off')
figF.tight_layout(); figF.savefig(os.path.join(FIG, 'a3_front.png'), dpi=130)
figP.tight_layout(); figP.savefig(os.path.join(FIG, 'a3_pullin.png'), dpi=130)
json.dump(summary, open(os.path.join(OUT, 'data', 'a3_wave.json'), 'w'), indent=1, default=float)
for c, s in summary.items():
    print('=====', c)
    for k, v in s.items():
        if k in ('fx_us', 'fy_us', 't_us'): continue
        print(k, v)
