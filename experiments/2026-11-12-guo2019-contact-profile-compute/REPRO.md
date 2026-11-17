# Воспроизведение: профиль контакта нитей после ускорения на claude-compute (2026-11-12)

## 1. Окружение

- **Машина:** `claude-compute`, KVM, Intel Xeon (Icelake), 2 ядра, 7 ГБ, без аппаратных счётчиков производительности (`perf stat -e cycles` — «not supported»), `perf_event_paranoid` = 1.
- **ПО:** Ubuntu, g++ 15.2.0 (Ubuntu 15.2.0-16ubuntu1), perf 7.0.14, Python 3.14.4.
- **Потоки OpenMP:** 1 (`taskset` ядро 1) и 2 (ядра 0–1, `OMP_PROC_BIND=close`), строго по одному прогону.
- **Сборки:**
  - Release `-O2` — `build`;
  - `-O2 -g -fno-omit-frame-pointer` — `build-prof`;
  - цель `guo2019_shot`.

  md5 бинарников `build` и `build-prof` совпадают со сборками на 8t в отчёте 2026-11-11 (`6677c73c…`, `96794909…`).

## 2. Коммит

`f1298aa12b4e2a9db0ddc747c75ed5cc822dab49` «Сжатые нормальные слоты: постоянная память сбора по отрезкам, тест Т8 (ревью Codex)», ветка `contact3_fix`. Копия дерева — `~/SPH_c3opt`.

## 3. Diff

`repro-diff.zip`:

| Файл | Что |
|---|---|
| `01-prof-counts.patch` | временные счётчики пар, примитивов и слотов (`SPH_PROF_COUNTS`), только для P.cnt и H.cnt |
| `MANIFEST.md` | состав, наложение и снятие |
| `HEAD.txt` | коммит базы |

```bash
cp -a ~/SPH_c3opt ~/SPH_c3cnt   # без каталогов build*
cd ~/SPH_c3cnt && git apply 01-prof-counts.patch   # или patch -p1 < 01-prof-counts.patch
```

## 4. Команды

```bash
# ткань -- копия ~/late/twaron_s25_dmg0.00415.fib с 8t; эталон 8t -- monitor.csv H.t1 и P.t1 из ~/c3opt_bench/res
REF8T=~/c3prof/ref8t bash scripts/prof_compute.sh ~/SPH_c3opt ~/SPH_c3cnt ~/c3prof/twaron_s25_dmg0.00415.fib ~/c3prof/res
# bench_tables.py и perf_threads.py берутся из ../2026-11-11-guo2019-law3-performance/scripts (или из $COMMON_SCRIPTS);
# counts_summary.py и libc_callers.py -- из каталога analyze.sh. Вызывающие libc в dwarf-профиле (главный tid --
# поток с наибольшим числом образцов) analyze.sh выполняет сам: <вывод>/P.dwarf2.libc_callers.txt
bash scripts/analyze.sh ~/c3prof/res ~/c3prof/an
# то же вручную:
T=$(perf script -i ~/c3prof/res/P.dwarf2/perf.data -F tid | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
perf script -i ~/c3prof/res/P.dwarf2/perf.data -F tid,ip,sym,dso | python3 scripts/libc_callers.py "$T" > results/P.dwarf2.libc_callers.txt
# повторный разбор профиля P.prof8 с 8t (отчёт 2026-11-11) исправленным perf_threads.py
python3 ../2026-11-11-guo2019-law3-performance/scripts/perf_threads.py P.prof8/perf.data --win pre:0.01-0.25 --win imp:0.35-1.0 --top 30 > results/P.prof8_8t.threads.txt
```

- `bench_tables.py` и `perf_threads.py` — из `../2026-11-11-guo2019-law3-performance/scripts/`.
- `perf_threads.py` исправлен в этом эксперименте: образцы с листовой функцией libc внутри региона OpenMP теперь относятся к работе региона. `P.prof2`, `H.prof2` и `P.prof8_8t` разобраны исправленной версией. В однопоточных `H.prof1` и `P.prof1` регионов нет, исправление на них не влияет.
- Профиль P.prof8 скопирован с 8t через локальную машину: прямой связи claude-compute → 8t нет. Бинарник по тому же пути и с тем же md5 есть на claude-compute.

## 5. Стоимость

- Прогоны: 1 ч 7 мин (08:39–09:46), машина больше ничем не занята.
- Разбор: ~10 мин, плюс повторный разбор ~5 мин.
- Место: ~0,5 ГБ в `~/c3prof/res` (`perf.data`; dwarf — 299 МБ). В git — только текстовые `results/`.

## 6. Проверка

- В `results/driver.log` все сверки дают «СОВПАДАЕТ», включая «H.t1 / P.t1 против 8t». Единственное исключение — H.t1 / H.prof1 из-за длины: первые 275 строк совпадают.
- Время на claude-compute: H.t1 ≈ 0,123 и P.t1 ≈ 0,111 с/шаг; H.t2 ≈ 0,069 и P.t2 ≈ 0,064 с/шаг.
- `counts_summary.txt`: ff-примитивов 429 260 до контакта (ровно, каждый шаг), активных нормальных слотов ff в позднем ударе в среднем ≈ 637, максимум 1188.
- H.prof1, до контакта: геометрия под `fiber_contact_pair_dt_limit` ≈ 15,8 % (10,1 + 5,7).
