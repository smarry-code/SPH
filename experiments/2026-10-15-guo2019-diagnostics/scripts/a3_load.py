"""A3 (часть 3): какие нити несут нагрузку к кадру перед 1-м разрывом: средняя сила T = sigma*A по нити."""
import os, sys, json, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

force = json.load(open(os.path.join(OUT, 'data', 'a3_force.json')))
cases = [c for c in ['1T_v65', '1T_v71', '2T_m01', '2T_m09'] if os.path.exists(os.path.join(OUT, 'data', c + '_frames.pkl'))]
res = {}
fig, ax = plt.subplots(1, len(cases), figsize=(4.6 * len(cases), 4), squeeze=False)
for k, c in enumerate(cases):
    R = pickle.load(open(os.path.join(OUT, 'data', c + '_frames.pkl'), 'rb'))
    fr = R['frames']; t = np.array([f['t'] for f in fr]) * 1e6
    ifb = int(np.max(np.where(t < force[c]['t_first_break_us'])[0]))
    f = fr[ifb]
    isx = f['yarn_dir_x']; ply = f['yarn_ply']; T = f['yarn_T_mean']; em = f['yarn_eps_mean']
    o = {}
    for p in sorted(set(ply)):
        for nm, mk in (('x', isx), ('y', ~isx)):
            s = (ply == p) & mk
            o['ply%d_%s' % (p, nm)] = dict(n_T_gt_20N=int((T[s] > 20).sum()), n_T_gt_100N=int((T[s] > 100).sum()), T_max=float(T[s].max()), T_sum=float(T[s].sum()),
                                           eps_mean_max=float(em[s].max()), n_eps_gt_1pct=int((em[s] > 0.01).sum()))
            ax[0, k].plot(f['yarn_xy'][s] * 1e3, T[s], 'o-' if nm == 'y' else 's--', ms=3, label='ply%d %s' % (p, nm))
    o['t_us'] = float(t[ifb]); res[c] = o
    ax[0, k].set_title('%s, t=%.0f мкс' % (CASES[c], t[ifb])); ax[0, k].set_xlabel('поперечная координата нити, мм'); ax[0, k].set_ylabel('средняя сила в нити, Н')
    ax[0, k].grid(alpha=.3); ax[0, k].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'a3_load_distribution.png'), dpi=120)
json.dump(res, open(os.path.join(OUT, 'data', 'a3_load.json'), 'w'), indent=1, default=float)
print(json.dumps(res, indent=1))
