"""B2: обработка data/b2/<tag>/{drum.csv,run/energy_log.csv} -> data/b2_table.md/.json, fig/b2_*.png.
Метрики: доля начальной кинетической энергии, ушедшая в D_friction к моментам 5, 20, 50, 100, 250, 450 мкс; период и огибающая
колебаний центральной бусины (если колебания есть); максимальное относительное смещение пар нитей на перекрёстках (признак макроскольжения)."""
import os, csv, json, glob, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
B2 = os.path.join(DATA, 'b2')
FIG = os.path.join(HERE, '..', 'fig')
D = 0.39e-3


def kv(path):
    d = {}
    for ln in open(path, encoding='utf-8'):
        if '=' in ln:
            k, v = ln.strip().split('=', 1)
            d[k] = v
    return d


def col(rows, k):
    return np.array([float(r[k]) for r in rows])


def main():
    tags = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(B2, '*', 'summary.txt')))
    res, curves = {}, {}
    for t in tags:
        S = kv(os.path.join(B2, t, 'summary.txt'))
        r = list(csv.DictReader(open(os.path.join(B2, t, 'drum.csv'))))
        e = list(csv.DictReader(open(os.path.join(B2, t, 'run', 'energy_log.csv'))))
        tt, KE, zc = col(r, 't'), col(r, 'KE'), col(r, 'zc')
        slip_max, slip_rms, frac = col(r, 'slip_max'), col(r, 'slip_rms'), col(r, 'frac_slip')
        te, KEf, Em, Df, Ub = col(e, 'time'), col(e, 'KE_fiber'), col(e, 'E_mech'), col(e, 'D_friction'), col(e, 'U_bond')
        KE0 = float(S['KE0'])
        out = dict(tag=t, mu=float(S['mu']), stretch=float(S['stretch']), v0=float(S['v0']), KE0=KE0, wall_s=float(S['wall_s']), steps=int(S['steps']))
        for tus in (2, 5, 20, 50, 100, 250, 450):
            j = np.searchsorted(te, tus * 1e-6)
            if j < len(te):
                out['D_frac_%dus' % tus] = float(Df[j] / KE0)
                out['KE_frac_%dus' % tus] = float(KEf[j] / KE0)
        # E_mech: изменение запаса (без D) от начала: E_mech(t) + D(t) - E_mech(0) должно быть ~ 0 (баланс), невязка относительно KE0
        out['balance_resid_over_KE0'] = float(np.max(np.abs(Em + Df - Em[0])) / KE0)
        out['E_mech_end_frac'] = float((Em[-1] - (Em[0] - KE0)) / KE0)  # доля KE0, оставшаяся в колебательной энергии
        out['slip_max_um'] = float(slip_max.max() * 1e6)
        out['slip_max_over_d'] = float(slip_max.max() / D)
        out['slip_rms_um_end'] = float(slip_rms[-1] * 1e6)
        out['frac_slip_max'] = float(frac.max())
        out['zc_amp_um'] = float(np.abs(zc).max() * 1e6)
        # колебания центральной бусины: пики |zc|
        pk, _ = find_peaks(zc, prominence=0.1 * max(np.abs(zc).max(), 1e-12))
        if len(pk) >= 3:
            out['period_us'] = float(np.mean(np.diff(tt[pk])) * 1e6)
            out['n_peaks'] = int(len(pk))
            a = zc[pk]
            out['amp_ratio_first_last'] = float(a[-1] / a[0])
            out['log_decrement'] = float(np.mean(np.log(a[:-1] / a[1:]))) if np.all(a > 0) else None
        else:
            out['period_us'] = None
            out['n_peaks'] = int(len(pk))
            out['log_decrement'] = None
        T_ref = 38.2e-6  # период колебаний этого лоскута при mu=0 (из прогонов mu=0)
        d20 = out.get('D_frac_20us')
        out['delta_eq_20us'] = float(0.5 * math.log(1.0 / max(1.0 - d20, 1e-6)) * T_ref / 20e-6) if d20 is not None and d20 > 0 else None
        res[t] = out
        curves[t] = (tt, KE, zc, te, KEf, Df, Em, slip_max)
    json.dump(res, open(os.path.join(DATA, 'b2_table.json'), 'w'), indent=1, ensure_ascii=False)
    L = ['| прогон | mu | стяжка | v0, м/с | KE0, Дж | D_fric/KE0 к 5 мкс | к 20 мкс | к 100 мкс | к 450 мкс | осталось в колебаниях, % KE0 | макс. смещение пар, мкм (доля d) | макс. доля контактов со смещением > 0.05d | период, мкс | эквив. лог. декремент (по 20 мкс) |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for t, o in res.items():
        L.append('| %s | %g | %g | %g | %.3g | %.3f | %.3f | %.3f | %.3f | %.1f | %.1f (%.3f) | %.3f | %s | %s |' % (
            t, o['mu'], o['stretch'], o['v0'], o['KE0'], o.get('D_frac_5us', float('nan')), o.get('D_frac_20us', float('nan')),
            o.get('D_frac_100us', float('nan')), o.get('D_frac_450us', float('nan')), 100 * o['E_mech_end_frac'], o['slip_max_um'],
            o['slip_max_over_d'], o['frac_slip_max'], '-' if o['period_us'] is None else '%.1f' % o['period_us'],
            '-' if o['delta_eq_20us'] is None else '%.2f' % o['delta_eq_20us']))
    open(os.path.join(DATA, 'b2_table.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    print('\n'.join(L))

    for stretch, v0 in ((0.02, 3.0), (0.01, 3.0), (0.02, 0.5)):
        sel = [t for t, o in res.items() if abs(o['stretch'] - stretch) < 1e-9 and abs(o['v0'] - v0) < 1e-9]
        if not sel:
            continue
        fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
        for t in sorted(sel):
            tt, KE, zc, te, KEf, Df, Em, sm = curves[t]
            mu = res[t]['mu']
            KE0 = res[t]['KE0']
            ax[0].plot(te * 1e6, KEf / KE0, label='mu=%g' % mu)
            ax[1].plot(te * 1e6, Df / KE0, label='mu=%g' % mu)
            ax[2].plot(tt * 1e6, zc * 1e6, label='mu=%g' % mu, lw=0.8)
        ax[0].set_ylabel('KE_fiber / KE0')
        ax[0].set_yscale('symlog', linthresh=1e-3)
        ax[1].set_ylabel('D_friction / KE0')
        ax[2].set_ylabel('z центральной бусины, мкм')
        for a_ in ax:
            a_.set_xlabel('время, мкс')
            a_.grid(alpha=0.3)
            a_.legend()
        fig.suptitle('лоскут 18x18, стяжка %g, v0=%g м/с' % (stretch, v0))
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, 'b2_s%g_v%g.png' % (stretch, v0)), dpi=130)
        plt.close(fig)


if __name__ == '__main__':
    main()
