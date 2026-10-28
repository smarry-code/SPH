"""Общие читатели: кадры fiber_output_*.txt, .fib, layout, monitor/energy_log."""
import os, glob
import numpy as np

DATA = os.environ.get('GUO_DATA', r'C:\Users\Maria\AppData\Local\Temp\claude\c--dip-sph-SPH--claude-worktrees-Guo2T\f944ccdc-8e09-4e7b-b37b-73a448cdb648\scratchpad\dd')
OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
FIG = os.path.join(OUT, 'fig')
CASES = {'1T_v65': '1T mu0.1 V65', '1T_v71': '1T mu0.1 V71',
         '2T_m01': '2T mu0.1 V137', '2T_m09': '2T mu0.9 V137'}
D = 3.9e-4; RHO = 1440.0; E_Y = 85e9; MBALL = 7e-3
AREA = np.pi * D**2 / 4


def load_frame(fn):
    """-> dict(t, step, oid, X(n,3), V(n,3), B(nb,2), l0, l, eps, eps_p)"""
    with open(fn) as f:
        head = f.readline().split()
        t = float(head[1]); step = int(head[3]); n = int(head[5]); nb = int(head[7])
        f.readline(); f.readline()          # '# id ...' , 'PARTICLES'
        P = np.loadtxt((f.readline() for _ in range(n)), ndmin=2)
        f.readline(); f.readline()          # '# i j ...', 'BONDS'
        Bd = np.loadtxt((f.readline() for _ in range(nb)), ndmin=2)
    return dict(t=t, step=step, oid=P[:, 1].astype(int), X=P[:, 3:6], V=P[:, 6:9],
                B=Bd[:, :2].astype(int), l0=Bd[:, 2], l=Bd[:, 3], eps=Bd[:, 5],
                eps_p=Bd[:, 6] if Bd.shape[1] > 6 else np.zeros(nb),
                sigma=Bd[:, 7] if Bd.shape[1] > 7 else E_Y * np.maximum(Bd[:, 5], 0))


def frames(case):
    return sorted(glob.glob(os.path.join(DATA, case, 'fiber_output_*.txt')))


def layout(case):
    rows = [l.split() for l in open(os.path.join(DATA, case, 'guo2019_fabric_layout.txt')) if not l.startswith('#')]
    L = []
    for r in rows:
        ply = 1
        if len(r) == 10:  # 2T: ply obj dir idx off clamped first last nlo nhi
            ply = int(r[0]); r = r[1:]
        L.append(dict(ply=ply, obj=int(r[0]), dir=r[1], idx=int(r[2]), off=float(r[3]), cl=int(r[4]),
                      a=int(r[5]), b=int(r[6])))
    return L


def read_fib(path):
    """-> list of arrays (n,3) по нитям в порядке файла"""
    yarns = []; cur = None; need = 0
    with open(path) as f:
        for ln in f:
            if ln.startswith('FIBER'):
                need = int(ln.split()[2]); cur = []; yarns.append(cur); continue
            if need > 0 and cur is not None and ln[0] in '-0123456789':
                cur.append([float(v) for v in ln.split()[:3]]); need -= 1
    return [np.array(y) for y in yarns]


def monitor(case):
    import pandas as pd
    return pd.read_csv(os.path.join(DATA, case, 'monitor.csv'))


def energy(case):
    import pandas as pd
    return pd.read_csv(os.path.join(DATA, case, 'energy_log.csv'))


def shot(case):
    s = open(os.path.join(DATA, case, 'shot.txt')).read().split()
    return dict(p.split('=', 1) for p in s if '=' in p)
