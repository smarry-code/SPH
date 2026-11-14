# Воспроизведение (2026-11-13)

## claude-compute: проверка этапа и парный замер

Дерево `~/SPH_c3sp` — копия `~/SPH_c3opt2` без каталогов сборки, поверх tar изменённых файлов этапа (сверка md5).
Эталон `~/SPH_c3opt` (код равен `f1298aa`); эталонные `monitor.csv` полной сцены — `~/c3prof/res`.

```bash
# на claude-compute, ~/c3sp_tests/run.sh = scripts/verify_compute.sh
STAGE=h3  PREV=bin_align bash run.sh   # bin_align -- сборка 5995594
STAGE=h1  PREV=bin_h3    bash run.sh
STAGE=fix PREV=bin_h1    bash run.sh
FROM=6 STAGE=... PREV=... bash run.sh  # только рестарт и замер (сборки и сверки 1-5 уже есть)
```

Бинарь этапа сохраняется в `~/c3sp_tests/bin_$STAGE`; журнал — `~/c3sp_tests/$STAGE/verify.log`.

## claude-compute: цепочка этапов в одной серии

```bash
bash ~/c3sp_tests/chain.sh   # = scripts/chain_compute.sh; журнал ~/c3sp_tests/chain/chain.log
```

Бинарники: `~/c3fs_tests/base_shot` (`0b9d3b0`), `bin_align` (`5995594`), `bin_h3` (`10fbcb1`), `bin_fix` (`c8973f7`).

## claude-compute-8t: бенчмарк (только через очередь)

```bash
# дерево ~/SPH_c3sp = rsync ~/SPH_c3opt2 (5da639c, код == 0b9d3b0) без build* + tar файлов 5da639c..c8973f7;
# сверка md5 всех git ls-files include src tests CMakeLists.txt
~/autostart/register_job.sh c3sp_bench bash ~/c3sp_bench/queue.sh   # queue.sh = scripts/queue_8t.sh
setsid nohup ~/autostart/jobs/c3sp_bench.job > /dev/null 2>&1 < /dev/null &
```

`bench_8t.sh`, `bench_tables.py`, `perf_threads.py` — те же, что в отчёте 2026-11-12 (`2f21706`, `0b9d3b0`).

Разбор (на 8t, после прогона):

```bash
cd ~/c3sp_bench
python3 bench_tables.py res > ana/bench_tables.txt
for v in H P; do
  python3 perf_threads.py res/$v.prof8/perf.data --freq 499 --threads 8 --win pre:0.03-0.28 --win imp:0.40-1.0 --top 30 > ana/$v.prof8.threads.txt
done
```

## Ревью Codex

`mcp__codex__codex`, `sandbox: read-only`, дифф `git diff 5995594..7147083`; ответ — `results/codex_review.txt`.
