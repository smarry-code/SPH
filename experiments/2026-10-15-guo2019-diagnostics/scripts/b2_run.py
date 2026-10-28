"""B2: серия колебаний лоскута. Запуск: python b2_run.py [--jobs 3]. Нужна сборка с -DSPH_ENERGY_LOG=ON (build-elog).
Результаты: data/b2/<tag>/{drum.csv,summary.txt,run/energy_log.csv}. Готовые прогоны пропускаются."""
import argparse, os, subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
FC = os.environ.get('FC', os.path.join(ROOT, 'build-elog', 'fabric_checks.exe'))
OUT = os.path.normpath(os.path.join(HERE, '..', 'data', 'b2'))


def cases():
    L = []
    for stretch in (0.02, 0.01):
        for mu in (0, 0.1, 0.9):
            L.append(('s%g_mu%g_v3' % (stretch, mu), ['--mu=%g' % mu, '--stretch=%g' % stretch, '--v0=3']))
    for mu in (0.003, 0.01, 0.03):  # где кончается режим прилипания: убывает ли затухание с mu
        L.append(('s0.02_mu%g_v3' % mu, ['--mu=%g' % mu, '--stretch=0.02', '--v0=3']))
    for mu in (0, 0.1, 0.9):
        L.append(('s0.02_mu%g_v0.5' % mu, ['--mu=%g' % mu, '--stretch=0.02', '--v0=0.5']))
    return L


def run(c):
    tag, args = c
    if os.path.exists(os.path.join(OUT, tag, 'summary.txt')):
        return tag, 'skip'
    with open(os.path.join(OUT, tag + '.log'), 'w') as lg:
        rc = subprocess.run([FC, 'drum', '--out=' + OUT, '--tag=' + tag, '--threads=2', '--t-end=5e-4', '--energy-log=10'] + args,
                            stdout=lg, stderr=subprocess.STDOUT).returncode
    return tag, rc


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=3)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for tag, rc in ex.map(run, cases()):
            print(tag, rc, flush=True)
