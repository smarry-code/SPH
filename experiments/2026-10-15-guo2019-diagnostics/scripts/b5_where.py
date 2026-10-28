"""B5: где находятся пары перекрытия одной нити (по кадрам)."""
import os, sys, pickle, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
force = json.load(open(os.path.join(OUT, 'data', 'a3_force.json')))
res = {}
for c in CASES:
    R = pickle.load(open(os.path.join(OUT, 'data', c + '_frames.pkl'), 'rb'))
    tfb = force[c]['t_first_break_us']
    for f in R['frames']:
        t = f['t'] * 1e6
        if 'hz_same_pos' not in f or len(f['hz_same_pos']) == 0: continue
        p = f['hz_same_pos'] * 1e3; isx = f['hz_same_dir']
        r = np.hypot(p[:, 0], p[:, 1]); corner = (np.abs(p[:, 0]) > 40) & (np.abs(p[:, 1]) > 40)
        res.setdefault(c, []).append(dict(t=round(t), n=len(p), frac_x_yarns=float(isx.mean()), frac_corner_gt40=float(corner.mean()),
                                          frac_r_lt15=float((r < 15).mean()), median_r=float(np.median(r)), after_break=bool(t > tfb)))
for c, L in res.items():
    print(c)
    for d in L: print('  ', d)
json.dump(res, open(os.path.join(OUT, 'data', 'b5_where.json'), 'w'), indent=1)
