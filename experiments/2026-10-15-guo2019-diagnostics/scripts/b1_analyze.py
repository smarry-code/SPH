"""B1: обработка data/b1/<tag>/capstan.csv -> таблица data/b1_table.md/.json, графики fig/b1_*.png.
Отношение T2/T1 берётся двумя способами: пик сглаженного отношения (возможность трения по всей дуге в момент выхода на скольжение)
и среднее за последние 25 % прогона (устоявшееся скольжение, T растёт по мере вытяжки нити)."""
import os, csv, json, glob, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
B1 = os.path.join(DATA, 'b1')
FIG = os.path.join(HERE, '..', 'fig')


def kv(path):
    d = {}
    for ln in open(path, encoding='utf-8'):
        if '=' in ln:
            k, v = ln.strip().split('=', 1)
            d[k] = v
    return d


def smooth(y, w=40):
    k = np.ones(w) / w
    return np.convolve(y, k, mode='same')


def main():
    tags = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(B1, '*', 'summary.txt')))
    res, curves = {}, {}
    for t in tags:
        S = kv(os.path.join(B1, t, 'summary.txt'))
        r = list(csv.DictReader(open(os.path.join(B1, t, 'capstan.csv'))))
        tt = np.array([float(x['t']) for x in r])
        T1 = np.array([float(x['T1']) for x in r])
        T2 = np.array([float(x['T2']) for x in r])
        ratio = np.where(T1 > 0, T2 / np.maximum(T1, 1e-12), np.nan)
        rs = smooth(np.nan_to_num(ratio, nan=1.0))
        n = len(rs)
        i0, i1 = n // 10, n - 20
        ipk = i0 + int(np.argmax(rs[i0:i1]))
        theta = float(S['theta_deg'])
        mu = float(S['mu'])
        th = math.exp(mu * math.radians(theta))
        end = float(np.nanmean(ratio[3 * n // 4:]))
        res[t] = dict(tag=t, drum=S['drum'], theta=theta, mu=mu, ds_mm=float(S['chord']) * 1e3, theory=th, discrete=float(S['exp_discrete']),
                      peak=float(rs[ipk]), peak_dev=100 * (rs[ipk] / th - 1), t_peak_us=float(tt[ipk]) * 1e6,
                      end=end, end_dev=100 * (end / th - 1), T1_init=float(S['T1_init']), T1_end=float(np.nanmean(T1[3 * n // 4:])),
                      T2_end=float(np.nanmean(T2[3 * n // 4:])), wall_s=float(S['wall_s']), steps=int(S['steps']), dt=float(S['dt_first']))
        curves[t] = (tt, ratio, rs)
    json.dump(res, open(os.path.join(DATA, 'b1_table.json'), 'w'), indent=1, ensure_ascii=False)
    L = ['| цилиндр | theta, град | mu | шаг бусин, мм | exp(mu*theta) | T2/T1 пик | откл., % | T2/T1 конец | откл., % | T1 нач./кон., Н | T2 кон., Н | wall, с |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for t, r in sorted(res.items(), key=lambda kv_: (kv_[1]['drum'], kv_[1]['ds_mm'], kv_[1]['theta'], kv_[1]['mu'])):
        L.append('| %s | %g | %g | %.3f | %.3f | %.3f | %+.1f | %.3f | %+.1f | %.2f / %.2f | %.1f | %.0f |' % (
            r['drum'], r['theta'], r['mu'], r['ds_mm'], r['theory'], r['peak'], r['peak_dev'], r['end'], r['end_dev'], r['T1_init'], r['T1_end'], r['T2_end'], r['wall_s']))
    open(os.path.join(DATA, 'b1_table.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    print('\n'.join(L))

    # график 1: отношение против exp(mu*theta), основная серия
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    xs = np.logspace(0, math.log10(80), 50)
    ax[0].plot(xs, xs, 'k-', lw=1, label='теория Эйлера')
    mk = {90: 'o', 180: 's', 270: '^'}
    for t, r in res.items():
        if r['drum'] == 'solid' and abs(r['ds_mm'] - 0.1) < 0.02:
            ax[0].plot(r['theory'], r['peak'], mk[int(r['theta'])], color='C0', mfc='none')
            ax[0].plot(r['theory'], r['end'], mk[int(r['theta'])], color='C1')
        if r['drum'] == 'fiber':
            ax[0].plot(r['theory'], r['peak'], mk[int(r['theta'])], color='C2', mfc='none')
            ax[0].plot(r['theory'], r['end'], mk[int(r['theta'])], color='C3')
    ax[0].plot([], [], 'o', color='C0', mfc='none', label='Solid, пик')
    ax[0].plot([], [], 'o', color='C1', label='Solid, конец')
    ax[0].plot([], [], 'o', color='C2', mfc='none', label='нить-нить, пик')
    ax[0].plot([], [], 'o', color='C3', label='нить-нить, конец')
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('exp(mu*theta)')
    ax[0].set_ylabel('T2/T1 в модели')
    ax[0].set_title('маркеры: o 90, квадрат 180, треугольник 270 град (шаг бусин 0.1 мм)')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3, which='both')
    # график 2: зависимость от шага бусин, theta = 180
    for mu, c in ((0.1, 'C0'), (0.3, 'C1'), (0.9, 'C2')):
        pts = [(r['ds_mm'], r['peak'] / r['theory']) for r in res.values() if r['drum'] == 'solid' and r['theta'] == 180 and abs(r['mu'] - mu) < 1e-9]
        pts.sort()
        if pts:
            ax[1].plot([p[0] for p in pts], [p[1] for p in pts], 'o-', color=c, label='mu=%g (пик / теория)' % mu)
    ax[1].axhline(1.0, color='k', lw=1)
    ax[1].set_xscale('log')
    ax[1].set_xlabel('шаг бусин нити (хорда), мм')
    ax[1].set_ylabel('(T2/T1 модели) / exp(mu*theta)')
    ax[1].set_title('theta = 180 град, цилиндр Solid')
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'b1_capstan.png'), dpi=130)
    plt.close(fig)

    # график 3: кривые отношения во времени (theta = 180, ds = 0.1)
    fig, ax = plt.subplots(figsize=(8, 5))
    for t, (tt, ratio, rs) in curves.items():
        r = res[t]
        if r['drum'] == 'solid' and r['theta'] == 180 and abs(r['ds_mm'] - 0.1) < 0.02 and r['mu'] > 0:
            ax.plot(tt * 1e6, rs, label='mu=%g' % r['mu'])
            ax.axhline(r['theory'], color='0.7', lw=0.7, ls=':')
    ax.set_yscale('log')
    ax.set_xlabel('время, мкс')
    ax.set_ylabel('T2/T1')
    ax.set_title('theta = 180 град: отношение натяжений во времени (пунктир -- exp(mu*theta))')
    ax.legend()
    ax.grid(alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'b1_time.png'), dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    main()
