"""B5: смятие / самопересечение нитей. Вход: data/<case>_frames.pkl, energy_log.csv, monitor.csv. Также баланс энергии на 1-м разрыве."""
import os, sys, json, pickle
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

force = json.load(open(os.path.join(OUT, 'data', 'a3_force.json')))
cases = [c for c in ['1T_v65', '1T_v71', '2T_m01', '2T_m09'] if os.path.exists(os.path.join(OUT, 'data', c + '_frames.pkl'))]
summ = {}
fig, ax = plt.subplots(2, 2, figsize=(14, 8))
for c in cases:
    R = pickle.load(open(os.path.join(OUT, 'data', c + '_frames.pkl'), 'rb'))
    fr = R['frames']; t = np.array([f['t'] for f in fr]) * 1e6
    e = energy(c); m = monitor(c); E0 = e.E_mech.iloc[0]
    S = {}
    # --- кадры
    tab = []
    for f in fr:
        sd = f['short_d']; sh = f['short_half']
        tab.append(dict(t_us=f['t'] * 1e6, nb=f['nb'], n_lt_d=sd['n'], n_lt_half=sh['n'], lt_d_ball=sd['ball'], lt_d_edge=sd['edge'],
                        lt_d_brokenyarn=sd['brk'], lt_d_intact=sd['intact'], min_l_over_l0=f['min_l_ratio'],
                        hz_same_n=f.get('hz_same_n', 0), hz_same_U=f.get('hz_same_U', 0.0), hz_same_n_diffoid=f.get('hz_same_n_diff_oid', 0),
                        hz_same_Udiffoid=f.get('hz_same_Udiff_oid', 0.0), hz_same_maxdelta_d=f.get('hz_same_maxdelta', 0.0),
                        hz_diff_n=f.get('hz_diff_n', 0), hz_diff_U=f.get('hz_diff_U', 0.0), hz_diff_maxdelta_d=f.get('hz_diff_maxdelta', 0.0)))
    T = pd.DataFrame(tab)
    T.to_csv(os.path.join(OUT, 'data', 'b5_frames_%s.csv' % c), index=False)
    S['frames_table'] = T.to_dict('list')
    S['max_n_lt_d'] = int(T.n_lt_d.max()); S['max_n_lt_half'] = int(T.n_lt_half.max())
    S['end_n_lt_d'] = int(T.n_lt_d.iloc[-1]); S['end_n_lt_half'] = int(T.n_lt_half.iloc[-1])
    S['end_hz_same_U'] = float(T.hz_same_U.iloc[-1]); S['max_hz_same_U'] = float(T.hz_same_U.max())
    S['max_hz_same_Udiffoid'] = float(T.hz_same_Udiffoid.max())
    S['E0'] = E0; S['max_hz_same_U_over_E0'] = S['max_hz_same_U'] / E0
    S['first_frame_lt_d_us'] = float(T.t_us[T.n_lt_d > 0].iloc[0]) if (T.n_lt_d > 0).any() else None
    S['first_frame_lt_half_us'] = float(T.t_us[T.n_lt_half > 0].iloc[0]) if (T.n_lt_half > 0).any() else None
    # --- energy_log: D_break_contact
    e['dBC'] = e.D_break_contact.diff()
    mm = m.set_index('step')
    def at(step, col):
        i = mm.index.get_indexer([step], method='nearest')[0]; return mm[col].iloc[i]
    e['n_bonds'] = [at(s, 'n_bonds') for s in e.step]; e['nbrk'] = [at(s, 'n_broken_v') + at(s, 'n_broken_h') for s in e.step]
    e['dnb'] = e.nbrk.diff(); e['dDbreak'] = e.D_break.diff()
    top = e.loc[e.dBC.abs().nlargest(6).index.sort_values(), ['step', 'time', 'D_break', 'D_break_contact', 'dBC', 'dnb', 'dDbreak', 'n_bonds', 'residual']]
    S['top_dBC'] = top.to_dict('records')
    S['D_break_contact_end'] = float(e.D_break_contact.iloc[-1]); S['D_break_end'] = float(e.D_break.iloc[-1])
    S['sum_pos_dBC'] = float(e.dBC[e.dBC > 0].sum()); S['sum_neg_dBC'] = float(e.dBC[e.dBC < 0].sum())
    S['max_abs_resid'] = float(e.residual.abs().max()); S['max_abs_resid_over_E0'] = S['max_abs_resid'] / E0
    S['clip_loss_end'] = float(e.clip_loss.iloc[-1]) if 'clip_loss' in e else None
    S['n_bonds_end'] = int(m.n_bonds.iloc[-1]); S['n_bonds_start'] = int(m.n_bonds.iloc[0])
    S['max_dnb_per_row'] = float(e.dnb.abs().max())
    # соотношение сразу после первого разрыва
    # --- баланс энергии на моменте первого разрыва
    tfb = force[c]['t_first_break_us'] * 1e-6
    i = int(np.argmin(np.abs(e.time.values - tfb)))
    r = e.iloc[i]
    S['partition_fb'] = dict(t_us=r.time * 1e6, KE_ball=float(r.KE_solid), KE_fiber=float(r.KE_fiber), U_bond=float(r.U_bond), U_hertz_ff=float(r.U_hertz_ff),
                             U_hertz_fs=float(r.U_hertz_fs), D_friction=float(r.D_friction), D_break=float(r.D_break), W_p_fiber=float(r.W_p_fiber),
                             W_v_fiber=float(r.W_v_fiber), D_xsph=float(r.D_xsph), D_contact_ss=float(r.D_contact_ss), E_mech=float(r.E_mech), residual=float(r.residual), E0=float(E0))
    r = e.iloc[-1]
    S['partition_end'] = dict(t_us=r.time * 1e6, KE_ball=float(r.KE_solid), KE_fiber=float(r.KE_fiber), U_bond=float(r.U_bond), U_hertz_ff=float(r.U_hertz_ff),
                              U_hertz_fs=float(r.U_hertz_fs), D_friction=float(r.D_friction), D_break=float(r.D_break), W_p_fiber=float(r.W_p_fiber),
                              E_mech=float(r.E_mech), residual=float(r.residual))
    summ[c] = S
    lab = CASES[c]
    ax[0, 0].plot(T.t_us, T.n_lt_d, 'o-', label=lab + ' l<d'); ax[0, 0].plot(T.t_us, T.n_lt_half, 'x--', color=ax[0, 0].lines[-1].get_color(), label=lab + ' l<0.5 l0')
    ax[0, 1].semilogy(T.t_us, T.hz_same_U + 1e-9, 'o-', label=lab)
    ax[1, 0].plot(e.time * 1e6, e.D_break_contact, label=lab + ' D_break_contact'); ax[1, 0].plot(e.time * 1e6, e.D_break, '--', color=ax[1, 0].lines[-1].get_color(), alpha=.6)
    ax[1, 1].plot(m.t * 1e6, m.n_bonds.iloc[0] - m.n_bonds, label=lab)
ax[0, 0].set_title('Короткие связи в кадрах (сплошные: l<d, пунктир: l<0.5 l0)'); ax[0, 0].set_xlabel('t, мкс'); ax[0, 0].set_ylabel('число связей')
ax[0, 1].set_title('Виртуальная энергия Герца перекрытий бусин одной нити (индекс >=2), Дж'); ax[0, 1].set_xlabel('t, мкс'); ax[0, 1].set_ylabel('Дж')
ax[1, 0].set_title('D_break_contact (сплошн.) и D_break (пунктир), Дж'); ax[1, 0].set_xlabel('t, мкс')
ax[1, 1].set_title('Число разорванных связей (n_bonds0 - n_bonds)'); ax[1, 1].set_xlabel('t, мкс')
for a in ax.flat: a.grid(alpha=.3); a.legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'b5_crush.png'), dpi=130)
json.dump(summ, open(os.path.join(OUT, 'data', 'b5_crush.json'), 'w'), indent=1, default=float)
for c, S in summ.items():
    print('=====', c)
    for k, v in S.items():
        if k != 'frames_table': print(k, v)
    print(pd.DataFrame(S['frames_table'])[['t_us', 'n_lt_d', 'n_lt_half', 'lt_d_ball', 'lt_d_edge', 'lt_d_brokenyarn', 'min_l_over_l0', 'hz_same_n', 'hz_same_U', 'hz_same_n_diffoid', 'hz_same_maxdelta_d', 'hz_diff_n', 'hz_diff_U']].to_string())
