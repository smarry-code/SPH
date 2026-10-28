"""B1: серия капстана (тест Эйлера). Запуск: python b1_run.py [--jobs 4]. Результаты: data/b1/<tag>/{capstan.csv,summary.txt}.
Готовые прогоны (есть summary.txt) пропускаются."""
import argparse, math, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
FC = os.environ.get('FC', os.path.join(ROOT, 'build', 'fabric_checks.exe'))
OUT = os.path.normpath(os.path.join(HERE, '..', 'data', 'b1'))


def tend(theta, mu):
    e = math.exp(mu * math.radians(theta))
    return '1e-4' if e < 1.001 else ('3e-4' if e <= 3 else ('6e-4' if e <= 12 else '1.5e-3'))


def cases():
    L = []
    for th in (90, 180, 270):  # основная серия: цилиндр Solid, шаг бусин 0.1 мм
        for mu in (0, 0.1, 0.3, 0.9):
            L.append(('solid_th%d_mu%g_ds0.1' % (th, mu), ['--drum=solid', '--theta-deg=%d' % th, '--mu=%g' % mu, '--ds=0.1e-3', '--t-end=' + tend(th, mu)]))
    for mu in (0.1, 0.3, 0.9):  # зависимость от шага бусин нити, theta = 180
        L.append(('solid_th180_mu%g_ds0.2' % mu, ['--drum=solid', '--theta-deg=180', '--mu=%g' % mu, '--ds=0.2e-3', '--t-end=' + tend(180, mu)]))
    L.append(('solid_th180_mu0.3_ds0.05', ['--drum=solid', '--theta-deg=180', '--mu=0.3', '--ds=0.05e-3', '--t-end=3e-4']))
    L.append(('solid_th180_mu0.3_ds0.6', ['--drum=solid', '--theta-deg=180', '--mu=0.3', '--ds=0.6e-3', '--t-end=3e-4']))
    for th in (90, 180):  # нить-нить: барабан из зажатых тонких нитей
        for mu in (0.3, 0.9):
            L.append(('fiber_th%d_mu%g_ds0.1' % (th, mu), ['--drum=fiber', '--theta-deg=%d' % th, '--mu=%g' % mu, '--ds=0.1e-3', '--t-end=' + tend(th, mu)]))
    return L


def run(c):
    tag, args = c
    if os.path.exists(os.path.join(OUT, tag, 'summary.txt')):
        return tag, 'skip'
    with open(os.path.join(OUT, tag + '.log'), 'w') as lg:
        rc = subprocess.run([FC, 'capstan', '--out=' + OUT, '--tag=' + tag, '--threads=3'] + args, stdout=lg, stderr=subprocess.STDOUT).returncode
    return tag, rc


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=4)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for tag, rc in ex.map(run, cases()):
            print(tag, rc, flush=True)
