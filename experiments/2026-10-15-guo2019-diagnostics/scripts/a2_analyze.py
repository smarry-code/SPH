"""A2: обработка data/a2/<tag>/strip.csv -> toe-участок, таблица data/a2_toe.json/.md, графики fig/a2_*.png."""
import os, csv, json, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
A2 = os.path.join(DATA, 'a2')
FIG = os.path.join(HERE, '..', 'fig')


def kv(path):
    d = {}
    for ln in open(path, encoding='utf-8'):
        if '=' in ln:
            k, v = ln.strip().split('=', 1)
            d[k] = v
    return d


def load(tag):
    r = list(csv.DictReader(open(os.path.join(A2, tag, 'strip.csv'))))
    eps = np.array([float(x['eps']) for x in r])
    f = np.array([float(x['F']) for x in r])
    nb = np.array([int(x['bonds']) for x in r])
    return eps, f, nb, kv(os.path.join(A2, tag, 'summary.txt'))


def smooth_grid(eps, f, nEA, h=0.0005):
    """Среднее по бинам eps шириной h -> (e, F/nEA)."""
    e_max = eps.max()
    edges = np.arange(0.0, e_max + h, h)
    idx = np.digitize(eps, edges) - 1
    e, y = [], []
    for k in range(len(edges) - 1):
        m = idx == k
        if m.any():
            e.append(eps[m].mean())
            y.append(f[m].mean() / nEA)
    return np.array(e), np.array(y)


def tangent(e, y, w=0.0075):
    """Локальная линейная регрессия y(e) в окне +-w -> dy/de (безразмерно: K_t/nEA)."""
    K = np.full(len(e), np.nan)
    for i in range(len(e)):
        m = np.abs(e - e[i]) <= w
        if m.sum() >= 4:
            K[i] = np.polyfit(e[m], y[m], 1)[0]
    return K


def analyze(tag):
    eps, f, nb, S = load(tag)
    nEA = float(S['nEA'])
    e, y = smooth_grid(eps, f, nEA)
    K = tangent(e, y)
    ipk = int(np.nanargmax(y))  # область до пика силы
    res = dict(tag=tag, kind=S['kind'], amp_d=float(S['amp_d']), material=S['material'], rate=float(S['rate']), mu=float(S['mu']),
               L0_mm=float(S['L0']) * 1e3, nEA=nEA, F_max=float(S['F_max']), eps_Fmax=float(S['eps_at_Fmax']),
               eps_first_break=float(S['eps_first_break']), wall_s=float(S['wall_s']), steps=int(S['steps']))
    for th in (0.1, 0.5):
        j = next((i for i in range(ipk + 1) if not np.isnan(K[i]) and K[i] >= th), None)
        res['eps_K%02d' % int(th * 100)] = float(e[j]) if j is not None else None
    jm = int(np.nanargmax(K[:ipk + 1])) if ipk > 0 else 0
    res['K_max_over_nEA'] = float(K[jm])
    res['eps_Kmax'] = float(e[jm])
    res['eps_intercept'] = float(e[jm] - y[jm] / K[jm]) if K[jm] > 0 else None
    return res, (e, y, K)


def fmt(x, p=4):
    return '-' if x is None else ('%.*f' % (p, x))


def main():
    tags = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(A2, '*', 'summary.txt')))
    allres, curves = {}, {}
    for t in tags:
        try:
            allres[t], curves[t] = analyze(t)
        except Exception as ex:
            print('ошибка', t, ex)
    json.dump(allres, open(os.path.join(DATA, 'a2_toe.json'), 'w'), indent=1, ensure_ascii=False)
    lines = ['| вариант | rate, 1/с | eps при K_t=10 % nEA | eps при K_t=50 % nEA | eps пересечения касательной | K_max/nEA | F_max/nEA | eps при F_max | eps 1-го разрыва | wall, с |',
             '|---|---|---|---|---|---|---|---|---|---|']
    for t, r in allres.items():
        lines.append('| %s | %g | %s | %s | %s | %.3f | %.4f | %.4f | %s | %.0f |' % (
            t, r['rate'], fmt(r['eps_K10']), fmt(r['eps_K50']), fmt(r['eps_intercept']), r['K_max_over_nEA'], r['F_max'] / r['nEA'],
            r['eps_Fmax'], fmt(r['eps_first_break']) if r['eps_first_break'] >= 0 else '-', r['wall_s']))
    open(os.path.join(DATA, 'a2_toe.md'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

    def has(t):
        return t in curves

    colors = {'w06_el': 'C0', 'w06_pl': 'C2', 'w05k4_el': 'C1', 'w05k4_pl': 'C3', 'w06k4_el': 'C4', 's_el': 'k', 's_pl': 'gray'}
    labels = {'w06_el': 'ткань A=0.6d k=2, упругая', 'w06_pl': 'ткань A=0.6d k=2, s30', 'w06k4_el': 'ткань A=0.6d k=4, упругая', 'w05k4_el': 'ткань A=0.5d k=4, упругая', 'w05k4_pl': 'ткань A=0.5d k=4, s30',
              's_el': 'прямые нити, упругие', 's_pl': 'прямые нити, s30'}
    rates = [300, 30, 10]
    fig, axs = plt.subplots(1, 3, figsize=(17, 5), sharey=True)
    for ax, rt in zip(axs, rates):
        for k, c in colors.items():
            t = '%s_r%d' % (k, rt)
            if has(t):
                e, y, K = curves[t]
                ax.plot(e * 100, y, color=c, lw=1.4, ls='--' if k.startswith('s_') else '-', label=labels[k])
        ax.set_title('скорость деформации %d 1/с' % rt)
        ax.set_xlabel('деформация eps, %')
        ax.set_xlim(0, 14)
        ax.set_ylim(0, 0.14)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel('F / (n E A)  (прямая упругая нить: = eps)')
    axs[0].legend(fontsize=8, loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'a2_force_strain.png'), dpi=130)
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), sharey=True)
    for ax, rt in zip(axs, rates):
        for k, c in colors.items():
            t = '%s_r%d' % (k, rt)
            if has(t) and not k.startswith('s_'):
                e, y, K = curves[t]
                ax.plot(e * 100, K, color=c, lw=1.4, label=labels[k])
        ax.axhline(0.1, color='0.5', lw=0.7, ls=':')
        ax.axhline(0.5, color='0.5', lw=0.7, ls=':')
        ax.set_title('K_t / (nEA), %d 1/с' % rt)
        ax.set_xlabel('деформация eps, %')
        ax.set_xlim(0, 14)
        ax.set_ylim(-0.1, 1.3)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel('касательная жёсткость / (n E A)')
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'a2_tangent.png'), dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.6))
    for t, c in (('w06_el_r300_mu0', 'C4'), ('w06_el_r300', 'C0'), ('w06_el_r300_mu09', 'C5')):
        if has(t):
            e, y, K = curves[t]
            ax.plot(e * 100, y, color=c, label=t)
    ax.set_xlabel('деформация eps, %')
    ax.set_ylabel('F / (n E A)')
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 0.14)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'a2_mu.png'), dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    main()
