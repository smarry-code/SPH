# Воспроизведение

Код: `contact3_fix` @ `2f21706` (код A–D — `cf4a264`, инструменты — `5da639c`, `2f21706`). Инструментирования нет, патчей к коду нет.

## Проверка на claude-compute

```bash
# дерево ~/SPH_c3opt2 = HEAD (md5 изменённых файлов), эталон ~/SPH_c3opt = f1298aa, ~/c3prof/res = полная сцена f1298aa на 2 потоках
setsid bash scripts/verify_compute.sh > ~/c3opt2_tests/run.out 2>&1 < /dev/null &
# итог: ~/c3opt2_tests/verify.log -> results/verify_log.txt
```

## Бенчмарк на claude-compute-8t (только через очередь)

```bash
# ~/c3opt2_bench/oldres -- ссылки на эталонные прогоны (c3opt_bench/res, иначе c3fix_bench/res), их пишет queue_8t.sh
~/autostart/register_job.sh c3opt2_bench bash ~/c3opt2_bench/queue.sh   # scripts/queue_8t.sh
# queue.sh ждёт FREE=8 и отсутствия guo2019_shot, затем:
#   RUNS="H.t1 P.t1 H.t1b P.t1b H.t2 P.t2 H.t4 P.t4 H.t8 P.t8 H.prof1 P.prof1 H.prof8 P.prof8" OLD_RES=~/c3opt2_bench/oldres \
#   bash bench_8t.sh ~/SPH_c3opt2 - ~/late/twaron_s25_dmg0.00415.fib
```

## Разбор

```bash
S=../2026-11-11-guo2019-law3-performance/scripts
python3 $S/bench_tables.py ~/c3opt2_bench/res > results/bench_tables.txt
for v in H P; do
  python3 $S/perf_threads.py ~/c3opt2_bench/res/$v.prof1/perf.data --freq 999 --threads 1 --win pre:0.02-0.45 --win imp:0.60-1.0 --top 40 > results/$v.prof1.threads.txt
  python3 $S/perf_threads.py ~/c3opt2_bench/res/$v.prof8/perf.data --freq 499 --threads 8 --win pre:0.03-0.28 --win imp:0.40-1.0 --top 30 > results/$v.prof8.threads.txt
done
python3 $S/perf_threads.py ~/c3opt_bench/res/P.prof8/perf.data --freq 499 --threads 8 --win pre:0.01-0.25 --win imp:0.35-1.0 --top 30 > results/old_P.prof8.threads.txt
python3 $S/perf_threads.py ~/c3fix_bench/res/H.prof8/perf.data --freq 499 --threads 8 --win pre:0.00-0.28 --win imp:0.37-1.0 --top 30 > results/old_H.prof8.threads.txt
```

Окна 8-поточных профилей взяты по долям шагов (130–1300 и 1700–4372 из 4372), а не по `wall_samples.txt`: прогон ~110 с, выборка раз в 10 с.
