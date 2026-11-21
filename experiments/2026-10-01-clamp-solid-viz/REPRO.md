# Воспроизведение: закрепление Solid, ES1 / ES2 / ES3 / EF3 / EK1 / ES4

## 1. Окружение

| | |
|---|---|
| Расчёт (удалённая машина `claude-compute`) | Ubuntu 26.04.1 LTS, Linux 7.0.0-34-generic x86_64, 8 ядер |
| Компилятор | g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0 |
| CMake | 4.2.3 |
| Потоки | `OMP_NUM_THREADS=1` на каждый процесс (до 8 процессов параллельно) |
| Анимации, графики (локально, macOS Darwin 24.6.0) | Python 3.14.7, numpy 2.5.3, matplotlib 3.11.2, Pillow 12.3.0 (venv; в системном Python этих пакетов нет) |

Сборка на g++ 15 идёт с `-DCMAKE_CXX_FLAGS="-include cmath"` — обход давнего дефекта `param.hpp`; исходники не правились.

## 2. Коммит

```
8ae109beaa7d34380efe163597c1f65f170da91b
clamp e2e: сцены Solid и Fiber (ES1-ES3, EF1-EF6, EK1), правки по ревью Codex
ветка clamp
```

## 3. Отличия рабочего дерева от этого коммита

Патчей нет: прогоны на чистом дереве HEAD. [`repro-diff.zip`](repro-diff.zip) — `MANIFEST.md` («патчей нет») и `HEAD.txt` (хэш, автор, дата, заголовок).

| Файл | Содержание |
|---|---|
| `HEAD.txt` | полный хэш, автор, дата, заголовок коммита |
| `MANIFEST.md` | состав архива, пояснение «патчей нет», проверка |

Проверено на момент прогонов: `git diff --quiet HEAD -- CMakeLists.txt src include tests tools` — различий нет (незакоммиченными были только документы `Clamp/`, `fiber/*.md`, каталоги `experiments/…` и `build_clamp/`; на сборку не влияют). Архив распакован `python -m zipfile -e` в отдельный каталог и сверен с исходными файлами (`diff -r`), наложение не требуется. На удалённую машину дерево скопировано `rsync -a --delete --exclude build_clamp --exclude .git --exclude 'build*' --exclude '*.gif' --exclude '*.png' ./ claude-compute:~/SPH_e2e_all/` (вместе с незакоммиченными документами, в сборке они не участвуют).

## 4. Команды

```bash
# 1. сборка (на удалённой машине, ~/SPH_e2e_all)
cd ~/SPH_e2e_all
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-include cmath"
cmake --build build --target test_e2e_clamp -j8

# 2. прогоны с кадрами; каждый — свой процесс, параллельно
cd build
export OMP_NUM_THREADS=1
F=$HOME/solid_frames
./test_e2e_clamp --case=ES1                   --frames-dir=$F --frames=80 > ES1.log
./test_e2e_clamp --case=ES1 --variant=pull    --frames-dir=$F --frames=80 > ES1-pull.log
./test_e2e_clamp --case=ES1 --solver=euler    --frames-dir=$F --frames=80 > ES1-euler.log
./test_e2e_clamp --case=ES2                   --frames-dir=$F --frames=80 > ES2.log
./test_e2e_clamp --case=ES3                   --frames-dir=$F --frames=80 > ES3.log
./test_e2e_clamp --case=EF3                   --frames-dir=$F --frames=80 > EF3.log   # три сцены: EF3, EF3_2, EF3_3
./test_e2e_clamp --case=EK1 > EK1.log                                                  # только лог
./test_e2e_clamp --case=ES4 > ES4.log                                                  # отложен, только WARN

# 3. на локальной машине: скопировать кадры в каталог ВНЕ репозитория (scratchpad)
rsync -a claude-compute:~/solid_frames/ <scratch>/solid_frames/
# оставить только основные прогоны (без *_free и *ref) ссылками в <scratch>/main/
# анимации и снимки (venv с numpy/matplotlib/pillow), из корня репозитория:
python tools/clamp_animate.py --max-frames 24 --dpi 55 --snapshot 16 --run-dir <scratch>/main/ES1_pc_ES1 --out ES1_pc_ES1.gif
#   так же для ES1_pc_pull_ES1, ES1_euler_ES1, ES2_pc_ES2, ES3_pc_ES3, EF3_pc_EF3, EF3_pc_EF3_2, EF3_pc_EF3_3;
#   для ES2/ES3 дополнительно --color-by vz --out ES2_pc_ES2_vz.gif / ES3_pc_ES3_vz.gif

# 4. данные и графики (из корня репозитория)
python experiments/2026-10-01-clamp-solid-viz/scripts/extract_displacement.py <scratch>/main   # data/ES2…/disp.csv, data/ES3…/disp.csv
python experiments/2026-10-01-clamp-solid-viz/scripts/plot_reaction.py                        # SVG и data/metrics.txt
```

В `data/<прогон>/` лежат `reaction.csv` и `clamp_mask.txt` из `solid_frames/<прогон>/`; логи — в `data/logs/`. Кадры `sph_output_*.txt`, `fiber_output_*.txt` и «свободные»/`ref`-каталоги в репозиторий не копировались.

Долгие прогоны запускались через `nohup sh -c "... > log 2>&1; date +%s > log.done" &` и опрашивались короткими `ssh` (каждый вызов `ssh` короче 2 минут).

## 5. Стоимость

Одновременно шла работа 8 процессов на 8 ядрах: ES1 37 с, ES1 pull 36 с, ES1 Euler 36 с, ES2 91 с, ES3 101 с, EF3 2 с, EK1 и ES4 < 1 с (стенное время, один поток на процесс). Кадры на удалённой машине — 137 МБ (не хранятся). Анимация одного прогона (24 кадра, dpi 55) — несколько секунд. Каталог — 1,52 МБ (данные 1,5 МБ, скрипты, описание, архив 1 КБ), в пределах 2,5 МБ (`experiments/README.md` §3). GIF (4,4 МБ), PNG (0,7 МБ), SVG (0,2 МБ) и логи в каталоге не хранятся. `bash scripts/rebuild_artifacts.sh <корень репозитория> <рабочий каталог> [all|run|anim|plots]` повторяет шаги 1–4 из §4 на одной машине: run — сборка и восемь запусков параллельно (кадры в `frames/`, логи в `logs/`), anim — GIF и PNG в `anim/`, plots — SVG и `metrics.txt` в `plots/` по хранимым `data/` и логам шага run; всё пишется в рабочий каталог, не в каталог эксперимента. Шаг 4 `extract_displacement.py` скрипт не повторяет: `data/*/disp.csv` хранятся.

## 6. Проверка повтора

Проверено 2026-10-01: ES1 и EF3 (все три сцены) повторены на той же сборке в отдельном каталоге кадров; `reaction.csv` ES1 и трёх сцен EF3 совпали с `data/*/reaction.csv` побайтово (`cmp`), логи отличаются только стенным временем (ES1: 1,07·10¹ с против 1,11·10¹ с в строке «сторож не оборвал расчёт»). Расчёт детерминирован при `OMP_NUM_THREADS=1`.

Должно совпасть:

- ES1: c_eff/c_L = 0,7757 (T_c = 1,1731·10⁻⁵ с); плато 1,4762·10³ / −1,4875·10³ Н (−21 % от 1,8785·10³); смена знака №1 1,2164·10⁻⁵ с; ∫R_x dt −4,46 % / −11,7 % от ±2 M v0; рост энергии 0,1654 %; 1221 шаг. ES1 pull: c_eff/c_L = 0,7695, ∫: −4,20 % / −11,0 %. ES1 Euler: c_eff/c_L = 0,7733, ∫: −3,43 % / −9,22 %.
- ES2: T = 4,5150·10⁻⁴ с (−6,76 % от T_EB = 4,8422·10⁻⁴ с), 7568 шагов, R_z амплитуда 54,619 Н, среднее −5,9364·10⁻² Н.
- ES3: T = 2,5097·10⁻⁴ с (+21,2 % от 2,0709·10⁻⁴ с), 3833 шага, R_x ≡ R_y ≡ 0,0, max \|R_z\| = 793,25 Н (из `data/metrics.txt`).
- EF3: vz после отскока 0,49690 / 0,49542 / 0,48783 м/с, ∫R_z dt = 9,3955·10⁻⁶ / 9,3816·10⁻⁶ / 9,3101·10⁻⁶ Н·с, шагов 1142 / 571 / 571; в В(б) Δv_x = 0,24768 м/с.
- EK1: 11 проверок ok, solid_max_h/dx = 2,0 на старте, 1,3 после 100 шагов. ES4: одна строка WARN «ES4 отложен».
- Все запуски заканчиваются `PASS (провалов 0)`. Побитовое совпадение логов не обещается на другом компиляторе или числе потоков.
