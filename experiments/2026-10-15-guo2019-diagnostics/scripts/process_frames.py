"""Один проход по кадрам fiber_output_* выстрела -> data/<case>_frames.pkl
Для A3 (фронт прогиба, профили eps, извитость, втягивание, места разрыва) и B5 (смятие нитей).
Запуск: python process_frames.py 1T_v65 [...]"""
import os, sys, pickle
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

KRES = 4   # узлов на период извитости (2 сегмента на полуволну)
K_HERTZ = 4 / 3 * (E_Y / (2 * (1 - 0.3**2))) * np.sqrt(D / 4)   # как в предыдущем разборе (точечная модель)


def process(case):
    fl = frames(case)
    lay = layout(case)
    ny = len(lay)
    X0 = None
    nbead = lay[-1]['b'] + 1
    yarn_of = np.zeros(nbead, int)
    for k, l in enumerate(lay): yarn_of[l['a']:l['b'] + 1] = k
    isx = np.array([l['dir'] == 'x' for l in lay])
    key = {(l['ply'], l['dir'], l['idx']): k for k, l in enumerate(lay)}
    nply = max(l['ply'] for l in lay)
    mon = monitor(case)
    res = dict(case=case, frames=[], lay=lay)
    prev_bonds = None; prevF = None
    for fn in fl:
        F = load_frame(fn)
        X = F['X']; t = F['t']
        if X0 is None: X0 = X.copy(); oid0 = F['oid'].copy(); ids0 = None
        dz = X[:, 2] - X0[:, 2]
        B = F['B']; l = F['l']; l0 = F['l0']
        r = dict(t=t, step=F['step'], nb=len(B), wmax=float(-dz.min()))
        # ---- центральные нити: профили
        prof = {}
        for ply in range(1, nply + 1):
            for dr in 'xy':
                for idx in (50, 51):
                    k = key[(ply, dr, idx)]
                    a, b = lay[k]['a'], lay[k]['b']
                    Xi = X[a:b + 1]; Xi0 = X0[a:b + 1]
                    ax = 0 if dr == 'x' else 1
                    # связи этой нити
                    m = (B[:, 0] >= a) & (B[:, 0] < b)
                    bi = B[m, 0] - a
                    eps_b = np.full(b - a, np.nan); eps_b[bi] = F['eps'][m]
                    epp_b = np.full(b - a, np.nan); epp_b[bi] = F['eps_p'][m]
                    # сглаженный профиль dz (скользящее среднее по периоду)
                    kern = np.ones(KRES) / KRES
                    dzs = np.convolve(dz[a:b + 1], kern, mode='same')
                    # локальная извитость на периоде: длина по узлам / хорда
                    loc = []
                    for s in range(0, len(Xi) - KRES, KRES):
                        seg = np.linalg.norm(np.diff(Xi[s:s + KRES + 1], axis=0), axis=1).sum()
                        ch = np.linalg.norm(Xi[s + KRES] - Xi[s])
                        loc.append((0.5 * (Xi[s, ax] + Xi[s + KRES, ax]), 0.5 * (Xi[s, 1 - ax] + Xi[s + KRES, 1 - ax]), seg / ch - 1))
                    # полная длина, сглаженная длина (через каждый KRES-й узел), упругое удлинение
                    idxs = list(range(0, len(Xi), KRES))
                    if idxs[-1] != len(Xi) - 1: idxs.append(len(Xi) - 1)
                    Ls = np.linalg.norm(np.diff(Xi[idxs], axis=0), axis=1).sum()
                    Ls0 = np.linalg.norm(np.diff(Xi0[idxs], axis=0), axis=1).sum()
                    La = np.linalg.norm(np.diff(Xi, axis=0), axis=1).sum()
                    La0 = np.linalg.norm(np.diff(Xi0, axis=0), axis=1).sum()
                    nint = int(m.sum())
                    prof[(ply, dr, idx)] = dict(s=Xi[:, ax], s0=Xi0[:, ax], other=Xi[:, 1 - ax], z=Xi[:, 2], dz=dz[a:b + 1], dzs=dzs,
                                                eps=eps_b, epp=epp_b, loc=np.array(loc), Ls=Ls, Ls0=Ls0, La=La, La0=La0,
                                                nbonds=nint, nbeads=b - a + 1, vx=F['V'][a:b + 1, 0], vz=F['V'][a:b + 1, 2])
        r['prof'] = prof
        # ---- втягивание свободных концов x-нитей
        pull = []
        for k, ll in enumerate(lay):
            if ll['dir'] != 'x': continue
            a, b = ll['a'], ll['b']
            pl = (abs(X0[a, 0]) - abs(X[a, 0])) , (abs(X0[b, 0]) - abs(X[b, 0]))
            pull.append((ll['ply'], ll['idx'], pl[0], pl[1]))
        r['pull'] = np.array(pull)
        # ---- целостность: связи по нитям
        cnt = np.bincount(yarn_of[B[:, 0]], minlength=ny)
        full = np.array([l_['b'] - l_['a'] for l_ in lay])
        r['yarn_broken'] = np.where(cnt < full)[0]
        yb = yarn_of[B[:, 0]]
        r['yarn_eps_mean'] = np.bincount(yb, weights=F['eps'], minlength=ny) / np.maximum(cnt, 1)
        r['yarn_T_mean'] = np.bincount(yb, weights=F['sigma'] * AREA, minlength=ny) / np.maximum(cnt, 1)
        r['yarn_xy'] = np.array([X[l_['a']:l_['b'] + 1, 1 if l_['dir'] == 'x' else 0].mean() for l_ in lay])
        r['yarn_dir_x'] = isx.copy(); r['yarn_ply'] = np.array([l_['ply'] for l_ in lay])
        # ---- B5: короткие связи
        ratio_d = l / D; ratio_l0 = l / l0
        sh = np.where(l < D)[0]
        shh = np.where(l < 0.5 * l0)[0]
        mid = 0.5 * (X[B[:, 0]] + X[B[:, 1]])
        broken_y = np.zeros(ny, bool); broken_y[r['yarn_broken']] = True
        def classify(ii):
            if len(ii) == 0: return dict(n=0, ball=0, edge=0, other=0, brk=0, intact=0, minratio=np.nan)
            mm = mid[ii]; rr = np.hypot(mm[:, 0], mm[:, 1])
            edge = (np.abs(mm[:, 0]) > 0.054) | (np.abs(mm[:, 1]) > 0.054)
            ball = rr < 0.010
            by = broken_y[yarn_of[B[ii, 0]]]
            return dict(n=len(ii), ball=int(ball.sum()), edge=int(edge.sum()), other=int((~ball & ~edge).sum()),
                        brk=int(by.sum()), intact=int((~by).sum()), minratio=float(ratio_d[ii].min()),
                        brk_ball=int((by & ball).sum()), brk_edge=int((by & edge).sum()))
        r['short_d'] = classify(sh); r['short_half'] = classify(shh)
        r['short_pos'] = mid[sh][:, :2] if len(sh) < 5000 else mid[sh][::max(1, len(sh) // 5000), :2]
        r['min_l_ratio'] = float(ratio_l0.min())
        # ---- B5: виртуальная энергия Герца для бусин одной нити (индекс-разница >= 2)
        tree = cKDTree(X)
        pr = tree.query_pairs(D, output_type='ndarray')
        if len(pr):
            same = yarn_of[pr[:, 0]] == yarn_of[pr[:, 1]]
            far = np.abs(pr[:, 0] - pr[:, 1]) >= 2
            sel = pr[same & far]
            Dp = np.linalg.norm(X[sel[:, 0]] - X[sel[:, 1]], axis=1)
            dl = D - Dp
            oid_diff = F['oid'][sel[:, 0]] != F['oid'][sel[:, 1]]
            r['hz_same_n'] = len(sel)
            r['hz_same_U'] = float((0.4 * K_HERTZ * dl**2.5).sum())
            r['hz_same_Udiff_oid'] = float((0.4 * K_HERTZ * dl[oid_diff]**2.5).sum())
            r['hz_same_n_diff_oid'] = int(oid_diff.sum())
            r['hz_same_maxdelta'] = float(dl.max() / D) if len(dl) else 0.0
            # где: центр пар
            r['hz_same_pos'] = 0.5 * (X[sel[:, 0], :2] + X[sel[:, 1], :2])
            r['hz_same_dir'] = isx[yarn_of[sel[:, 0]]]
            # пары разных нитей (реальный контакт; для справки)
            diff = ~same
            Dd = np.linalg.norm(X[pr[diff, 0]] - X[pr[diff, 1]], axis=1)
            r['hz_diff_n'] = int(diff.sum())
            r['hz_diff_U'] = float((0.4 * K_HERTZ * (D - Dd)**2.5).sum())
            r['hz_diff_maxdelta'] = float((D - Dd).max() / D) if diff.any() else 0.0
        # ---- разорванные в этом кадре (относительно предыдущего)
        if prevF is not None:
            sp = set(map(tuple, prevF['B'])); sc = set(map(tuple, B))
            gone = sorted(sp - sc)
            gl = []
            pbond = {tuple(b): k for k, b in enumerate(prevF['B'])}
            for (i, j) in gone:
                k = pbond[(i, j)]
                yk = yarn_of[i]
                gl.append(dict(i=i, j=j, yarn=yk, ply=lay[yk]['ply'], dir=lay[yk]['dir'], idx=lay[yk]['idx'], pos_in_yarn=i - lay[yk]['a'],
                               xyz_prev=0.5 * (prevF['X'][i] + prevF['X'][j]), xyz_cur=0.5 * (X[i] + X[j]),
                               eps_prev=float(prevF['eps'][k]), epsp_prev=float(prevF['eps_p'][k]), t_prev=prevF['t']))
            r['gone'] = gl
            # профиль eps всех связей предыдущего кадра для разорванных нитей: максимум eps по нити
        # распределение eps по всем связям (гистограмма) и максимум
        r['eps_max'] = float(F['eps'].max()); r['eps_p_max'] = float(F['eps_p'].max())
        r['eps_hist'] = np.histogram(F['eps'], bins=np.linspace(-0.2, 0.08, 57))[0]
        # поле dz (грубо) для карты: центры связей на сетке не храним -- только осевые срезы выше
        res['frames'].append(r)
        prevF = dict(B=B, X=X.copy(), eps=F['eps'].copy(), eps_p=F['eps_p'].copy(), t=t)
        print(case, 'frame', F['step'], 't=%.1f us' % (t * 1e6), 'nb', len(B), 'wmax=%.2f mm' % (r['wmax'] * 1e3),
              'short<d', len(sh), 'short<.5l0', len(shh), flush=True)
    pickle.dump(res, open(os.path.join(OUT, 'data', case + '_frames.pkl'), 'wb'))


if __name__ == '__main__':
    for c in sys.argv[1:]: process(c)
