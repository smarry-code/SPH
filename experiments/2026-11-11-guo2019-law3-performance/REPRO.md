# Воспроизведение: производительность закона 3 с проекцией (2026-11-11)

## 1. Окружение

- **Машина:** `claude-compute-8t`, Intel Xeon (Cascadelake), 2 сокета × 4 ядра, без SMT, 1 узел NUMA.
- **ПО:** Linux 7.0.0-38-generic, g++ 15.2.0 (Ubuntu 15.2.0-16ubuntu1), CMake 4.2.3, Python 3.14.4, perf 7.0.14.
- **Потоки OpenMP:** 1, 2, 4, 8 — строго по одному прогону, `OMP_PROC_BIND=close`, `taskset` (один поток — ядро 2, N потоков — ядра 0..N−1).
- **Сборки:**
  - Release `-O2` — `build`;
  - `-O2 -g -fno-omit-frame-pointer` — `build-prof`, для perf;
  - цель только `guo2019_shot`, без `SPH_ENERGY_LOG`.

## 2. Коммит

- NEW: `16587564845ee03c7a348888a6b9e32cc218c1eb` «Проекция нормального контакта: правки по ревью Codex».
- BASE (вариант Pb): `6edcdc54201ee43309f6840b873600d33e0905e9` «Проекция нормального контакта: правки по ревью (метрики остатка, предупреждения, тесты Т4/Т6)».

## 3. Diff

Нет: прогоны шли на закоммиченных исходниках без изменений, временной инструментации не было. Архив `repro-diff.zip` не нужен.

## 4. Команды

```bash
# исходники NEW и BASE -- копии деревьев указанных коммитов в ~/SPH_c3fix и ~/SPH_c3fix_base
bash scripts/bench_8t.sh ~/SPH_c3fix ~/SPH_c3fix_base ~/late/twaron_s25_dmg0.00415.fib ~/c3fix_bench/res
python3 scripts/bench_tables.py ~/c3fix_bench/res > results/bench_tables.txt
for v in H P; do
  python3 scripts/perf_threads.py ~/c3fix_bench/res/$v.prof8/perf.data --win pre:0.01-0.26 --win imp:0.35-1.0 --top 30 > results/$v.prof8.threads.txt
  python3 scripts/perf_threads.py ~/c3fix_bench/res/$v.prof1/perf.data --win pre:0.03-0.45 --win imp:0.60-0.98 --top 40 > results/$v.prof1.threads.txt
done
```

- Ключи расчёта — в `scripts/bench_8t.sh` (`COMMON` и `V[...]`).
- Окна perf — доли длительности записи. Они пересчитаны из окон по шагам (130–1300 и 1700–конец) в колонке «окна perf» таблицы `bench_tables.py`.
- Проверка правок до бенчмарка — `~/c3fix_tests/run.sh` на том же дереве, лог в `results/verify.log`.

## 5. Стоимость

- Бенчмарк: 4 ч 47 мин стенного времени (16:04–20:51), на 8t больше ничего не работало.
- Разбор профилей: ~10 мин на файл.
- Место: 669 МБ в `~/c3fix_bench/res`, из них основное — `perf.data`. В git — только текстовые `results/`.

## 6. Проверка

Повтор должен дать (шум 1,5–2 %, абсолютное время зависит от машины):

- H.t1 и P.t1: 4372 шага, 0,37 и 0,44 с/шаг;
- P/H = 1,20 в одном потоке и 1,29 на 8 потоках;
- ускорение на 8 потоках: H ≈ 5,2, P ≈ 4,8;
- `monitor.csv` H и P на 1 и 8 потоках побайтно равны (строки «детерминизм» в `results/driver.log`);
- в профиле на 8 потоках последовательный код главного потока занимает ≈ 26 % его времени, а ожидание на барьерах ≈ 0,1 %.
