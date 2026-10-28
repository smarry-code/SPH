"""A2: серия растяжения полоски. Запуск: python a2_run.py [--jobs 6]. Результаты: data/a2/<tag>/{strip.csv,summary.txt}.
Готовые прогоны (есть summary.txt) пропускаются. Тяжёлая группа подвижного конца: M_total = 1e5 кг (по умолчанию в fabric_checks)."""
import argparse, os, subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
FC = os.environ.get('FC', os.path.join(ROOT, 'build', 'fabric_checks.exe'))
OUT = os.path.normpath(os.path.join(HERE, '..', 'data', 'a2'))
FIB = os.path.join(ROOT, 'experiments', '2026-10-12-guo2019-1t-plastic-mu', 'twaron_plastic_visco_rate_s30.fib')


def cases():
    L = []

    def add(tag, *args):
        L.append((tag, list(args)))

    for rate in (300, 30, 10):
        for amp, a_tag in ((0.6, '06'),):
            add('w%s_el_r%d' % (a_tag, rate), '--kind=weave', '--amp-d=%g' % amp, '--material=elastic', '--rate=%d' % rate, '--mu=0.1', '--eps-max=0.16')
            add('w%s_pl_r%d' % (a_tag, rate), '--kind=weave', '--amp-d=%g' % amp, '--material=plastic', '--rate=%d' % rate, '--mu=0.1', '--eps-max=0.16', '--fib-inelastic=' + FIB)
    for rate in (300, 30):
        add('s_el_r%d' % rate, '--kind=straight', '--material=elastic', '--rate=%d' % rate, '--mu=0.1', '--eps-max=0.05')
        add('s_pl_r%d' % rate, '--kind=straight', '--material=plastic', '--rate=%d' % rate, '--mu=0.1', '--eps-max=0.12', '--fib-inelastic=' + FIB)
    # минимум амплитуды: ломаная k=2 при A < ~0.57 d перекрывается (см. README), поэтому A = 0.5 d -- при k = 4 (гладкий профиль)
    for rate in (300, 30):
        add('w05k4_el_r%d' % rate, '--kind=weave', '--amp-d=0.5', '--seg=4', '--material=elastic', '--rate=%d' % rate, '--mu=0.1', '--eps-max=0.16')
        add('w05k4_pl_r%d' % rate, '--kind=weave', '--amp-d=0.5', '--seg=4', '--material=plastic', '--rate=%d' % rate, '--mu=0.1', '--eps-max=0.16', '--fib-inelastic=' + FIB)
        add('w06k4_el_r%d' % rate, '--kind=weave', '--amp-d=0.6', '--seg=4', '--material=elastic', '--rate=%d' % rate, '--mu=0.1', '--eps-max=0.16')
    add('w057_el_r300', '--kind=weave', '--amp-d=0.57', '--material=elastic', '--rate=300', '--mu=0.1', '--eps-max=0.16')
    add('w05_el_r300', '--kind=weave', '--amp-d=0.5', '--material=elastic', '--rate=300', '--mu=0.1', '--eps-max=0.02')  # ожидается срыв (перекрытие)
    for mu in (0, 0.9):
        add('w06_el_r300_mu%g' % mu, '--kind=weave', '--amp-d=0.6', '--material=elastic', '--rate=300', '--mu=%g' % mu, '--eps-max=0.16')
    return L


def run(c):
    tag, args = c
    if os.path.exists(os.path.join(OUT, tag, 'summary.txt')):
        return tag, 'skip'
    with open(os.path.join(OUT, tag + '.log'), 'w') as lg:
        rc = subprocess.run([FC, 'strip', '--out=' + OUT, '--tag=' + tag, '--threads=1'] + args, stdout=lg, stderr=subprocess.STDOUT).returncode
    return tag, rc


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=6)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        for tag, rc in ex.map(run, cases()):
            print(tag, rc, flush=True)
