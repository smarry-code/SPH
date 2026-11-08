# Воспроизведение: Guo 2019, 1T `base`, окно при законе 3 + выталкивание, μ = 0,27 (V0 = 69, 70, 71, 82 м/с)

## 1. Окружение

| | |
|---|---|
| Расчёт: машина | claude-compute-8t: Linux Ubuntu 26.04.1 LTS, 8 ядер |
| Компилятор | g++ 15.2.0 (Ubuntu 15.2.0-16ubuntu1) |
| Сборка | CMake 4.2.3, Ninja, `-DCMAKE_BUILD_TYPE=Release -DSPH_ENERGY_LOG=ON`, `CMAKE_CXX_FLAGS` пусто |
| Потоки | по 4 потока OpenMP на выстрел (`OMP_NUM_THREADS=4`, `--threads=4`), не больше двух выстрелов сразу (≤ 8 потоков на сервере); рядом с V0 = 82 и 71 до 11:17 шла чужая серия proj_* (4 × 1 поток) |
| Драйвер | Python 3.14.4 (сервер), `scripts/window_bisect.py` + `scripts/window_pinball.py` |
| Графики, анимации | Windows 11, Python 3.12.10, matplotlib 3.11.2, numpy 2.5.4, Pillow 13.0.0 |

## 2. Коммит

```
9437d7d627db2621f05fbcb15a7f46a1bc0ab26a
guo2019 1T mu 0,27: данные для анимаций во всех выстрелах, кадры _long жёсткими ссылками
ветка guo1t_law3push_mu027 (от fiber_new ec2439fabfd263cc64c3ca11c33bd9acdf2a5dfc)
```

Бинарник `guo2019_shot` собран из `9437d7d` (`~/mu027/src/SOURCE_COMMIT` на сервере). `git diff ec2439f 9437d7d` вне `experiments/` пуст: код расчёта — `fiber_new` `ec2439f` (детерминированное суммирование OpenMP — `6f1b63f`, слияние `fiber_divergence_diag` — `cd25e5f`).

Скрипты драйвера в рабочих каталогах сервера (изменения только в `experiments/2026-11-08-guo2019-1t-law3push-mu027/scripts/`):

| Серия | Скрипты | Отличие от `9437d7d` |
|---|---|---|
| 71 / 82 (`mu027_window`) | `e51ed6c68a59f8a2e0254abf0bf8f6d4cac4a9fc` | допуск выстрела по свободным потокам (`start_gated`), `run_series.sh` без ожидания конца proj_* |
| 69 / 70 (`mu027_v69_70`) | `3ae6e18d22b2672757c77dd9bd47c834e19eca69` | то же + `run_series.sh` с каталогом и скоростями в аргументах |
| постобработка `_long` (все 4) | `42ea87aac21f6be490c621894a4654f556039312` | `fix_monitor_steps.py`: сквозной номер шага в склеенном `monitor.csv` |

Аргументы `guo2019_shot`, которые строит драйвер, во всех трёх версиях одинаковы.

## 3. Отличия рабочего дерева от этого коммита

Патчей нет: бинарник собран из архива чистого коммита `9437d7d`, скрипты — закоммиченные версии из §2. [`repro-diff.zip`](repro-diff.zip) содержит `HEAD.txt` (хэш, автор, дата, заголовок) и `MANIFEST.md` с этим утверждением.

| Файл архива | Содержание |
|---|---|
| `HEAD.txt` | полный хэш, автор, дата и заголовок `9437d7d` |
| `MANIFEST.md` | «патчей нет», версии скриптов по сериям, проверка состояния |

Наложение не требуется: `git checkout 9437d7d627db2621f05fbcb15a7f46a1bc0ab26a` (или `git worktree add --detach <путь> <хэш>`); для повтора с теми же скриптами — `e51ed6c` / `3ae6e18`.

## 4. Команды

```bash
# 1. исходники на сервер: архив коммита без тяжёлых каталогов experiments/ (как разворачивался ~/mu027/src);
#    core.autocrlf=false обязателен на Windows: иначе в архив попадают CRLF и ломается checkpoint_cli_test
git -c core.autocrlf=false archive --format=tar.gz -o src.tgz 9437d7d -- <пути дерева, см. ниже>
scp src.tgz claude-compute-8t:mu027/ && ssh claude-compute-8t 'cd ~/mu027 && mkdir src && tar -xzf src.tgz -C src'
#    Пути: всё дерево, кроме experiments/, плюс файлы experiments/, на которые ссылается CMakeLists.txt
#    (сборка их требует), и каталог этого эксперимента. Полный архив (342 МБ) тоже годится.

# 2. сборка (сервер; во время чужой серии -- nice -n 19, -j1)
cmake -S ~/mu027/src -B ~/mu027/src/build-elog -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPH_ENERGY_LOG=ON
cmake --build ~/mu027/src/build-elog --target guo2019_shot -j8

# 3. рабочие каталоги и задачи автозапуска (скрипты -- из scripts/ этого каталога в версии §2)
W1=~/vbl1T_law3push_mu027_det; W2=~/vbl1T_law3push_mu027_v69_70
for W in $W1 $W2; do mkdir -p $W; cp window_bisect.py window_pinball.py server/run_series.sh ../twaron_s25_dmg0.00415.fib $W/; done
~/autostart/register_job.sh mu027_window $W1/run_series.sh                 # [71; 82]: пара 71 / 82
~/autostart/register_job.sh mu027_v69_70 $W2/run_series.sh $W2 69 70 70     # пара 69 / 70
setsid nohup ~/autostart/jobs/mu027_window.job > ~/autostart/mu027_window.out 2>&1 < /dev/null &
setsid nohup ~/autostart/jobs/mu027_v69_70.job > ~/autostart/mu027_v69_70.out 2>&1 < /dev/null &

# run_series.sh W VMIN VFIRST VMAX запускает (без OMP_NUM_THREADS у python; выстрелам драйвер ставит 4):
env -u OMP_NUM_THREADS python3 $W/window_bisect.py --work=$W --bin=~/mu027/src/build-elog/guo2019_shot \
  --fib=$W/twaron_s25_dmg0.00415.fib --fiber-E=8e10 --eps-break=0.05 --config=1T --level=base \
  --contact=3 --mu=0.27 --verdict-ball --vmin=$VMIN --v-first=$VFIRST --vmax=$VMAX --tol=1.5 \
  --extra="--dt-scale=0.5 --impulse-push-tau=5e-7" --threads=4 --threads-single=4 \
  --t-max=1.5e-3 --t-max-long=3e-3 --save-t=2e-5 --energy-log=100 --checkpoint-every=20000 \
  --slot-script=~/proj_runs/threads_free.sh --proj-status=~/proj_runs/proj_status.sh \
  --min-uptime=600 --slot-poll=300 --slot-verify-s=120

# Каждый выстрел драйвер запускает так (V -- скорость, D -- каталог выстрела):
OMP_NUM_THREADS=4 guo2019_shot --config=1T --level=base --v0=V --mu=0.27 --out=D --frames=off --save-t=2e-05 \
  --threads=4 --n-mon=10 --n-diag=100 --repeat=1 --no-bracket --t-max=0.0015 \
  --fiber-inelastic=twaron_s25_dmg0.00415.fib --eps-break=0.05 --fiber-E=8e+10 --fs-contact=3 --ff-contact=3 \
  --energy-log=100 --dt-scale=0.5 --impulse-push-tau=5e-7 --checkpoint-every=20000
# «не определён» -> каталог D копируется в D_long (кадры -- жёсткими ссылками), запуск с --restart --t-max=0.003;
# monitor.csv базы откладывается в monitor.part1.csv и после выстрела склеивается с новым (window_pinball.splice_shot).

# 4. постобработка (2026-11-08, после серий): сквозной шаг монитора в _long (в версии 42ea87a склейка делает это сама)
python3 fix_monitor_steps.py $W1/mu0.27_v*_long/1T_*/ ; python3 fix_monitor_steps.py $W2/mu0.27_v*_long/1T_*/

# 5. данные каталога: скачать queue.log, window.md, mu0.27_v*.log, shots.csv, shot.txt, monitor.csv и energy_log.csv
#    из _long в <src>/<имя рабочего каталога>/... (структура как на сервере), затем
python scripts/thin_data.py <src> data          # monitor -- каждая 50-я строка, energy_log -- каждая 10-я
python scripts/plot_series.py .                  # v_z.svg, breaks.svg, energy.svg, data/metrics.txt

# 6. анимации 69 и 70: скачать из $W2/mu0.27_v<V>_long/1T_base_v<V>.0_r0 логи и каждый второй кадр
#    (fiber_output_* и sph_output_* с нечётным порядковым номером + последний) в <DATA>/law3push_mu0.27_v<V>_long/,
#    затем
PY=<python с matplotlib> DATA=<DATA> OUT=<куда> bash scripts/viz/make_viz_69_70.sh 69 70
```

## 5. Стоимость

| V0, м/с | Прогон | Время (сервер, MSK) | Стенное время, с | Шагов | с/шаг | Кадров (нити) | Кадры, ГБ |
|---|---|---|---|---|---|---|---|
| 82 | база | 04:23:58–07:53:38 | 12 579,8 | 162 433 | 0,0774 | 75 | 0,54 |
| 82 | продолжение | 07:53:39–08:58:43 | 3 904,1 | 44 564 | 0,0876 | 96 (75 — ссылки) | 0,69 |
| 71 | база | 06:06:23–09:37:04 | 12 640,5 | 162 238 | 0,0779 | 75 | 0,54 |
| 71 | продолжение | 09:37:04–11:17:12 | 6 007,5 | 56 224 | 0,1069 | 102 (75 — ссылки) | 0,74 |
| 70 | база | 14:44:13–18:05:46 | 12 092,1 | 162 222 | 0,0745 | 75 | 0,54 |
| 70 | продолжение | 18:05:47–19:16:58 | 4 270,6 | 53 819 | 0,0794 | 101 (75 — ссылки) | 0,73 |
| 69 | база | 14:46:15–18:08:49 | 12 154,3 | 162 216 | 0,0749 | 75 | 0,54 |
| 69 | продолжение | 18:08:50–19:35:52 | 5 221,8 | 64 486 | 0,0810 | 106 (75 — ссылки) | 0,76 |

Серия 71 / 82 — 6,9 ч (04:23–11:17; V0 = 71 ждал свободных потоков 1 ч 42 мин), серия 69 / 70 — 4,9 ч (14:44–19:35). Доля монитора 2,4–2,8 %. Рабочие каталоги на сервере: `~/vbl1T_law3push_mu027_det` 1,5 ГБ, `~/vbl1T_law3push_mu027_v69_70` 1,6 ГБ (кадры; в репозитории не хранятся). Чекпоинты (≈ 11 МБ) удалены драйвером после определённого исхода.

Каталог эксперимента в git 0,90 МБ (данные 0,67 МБ, SVG 0,17 МБ, скрипты, README/REPRO, архив) — в пределах 2,5 МБ. Анимации и PNG `anim/` (37 МБ) не хранятся: `bash scripts/rebuild_artifacts.sh <корень репозитория> <рабочий каталог> [all|run|anim|plots]` — сборка, пара 69 / 70 заново тем же драйвером (≈ 5 ч на 8 ядрах, кадры ≈ 1,6 ГБ), анимации (≈ 2 мин), SVG по `data/` (секунды); всё пишется в рабочий каталог (`build/`, `run/`, `anim/`, `plots/`), не в каталог эксперимента. Стадии `anim` (по скачанным кадрам) и `plots` проверены 2026-11-08: SVG и `data/metrics.txt` совпали с хранимыми побайтово. Стадия `run` целиком не прогонялась (те же ключи, что у серии; без `~/proj_runs/*.sh` драйвер считает свободные потоки по `--threads=` идущих `guo2019_shot`).

## 6. Проверка повтора

Должно получиться (`data/shot_v*_long.txt`):

| V0, м/с | outcome | V_R, м/с | t_end, с | n_broken_v / h | t_first_break (база), с |
|---|---|---|---|---|---|
| 69 | not_perforated | 0 | 0,00210080845 | 10 / 0 | 0,00045405781 |
| 70 | perforated | 0,324422314 | 0,00200137725 | 10 / 0 | 0,00044837415 |
| 71 | perforated | 0,23334177 | 0,00202375377 | 10 / 0 | 0,000441369045 |
| 82 | perforated | 6,97951923 | 0,0019145669 | 11 / 1 | 0,000378203143 |

Счёт детерминирован при любом числе потоков (`6f1b63f`): с тем же компилятором и флагами повтор должен дать эти значения точно, а `energy_log.csv` и кадры — побайтово (на дымовом прогоне с `kill -9` и рестартом из чекпоинта `energy_log.csv` и все кадры совпали с прогоном без обрыва побайтово). Независимого повтора выстрелов этой серии не было. С другим компилятором или флагами побитового совпадения ждать нельзя; у границы (V_R 0,2–0,3 м/с при пороге 0,1 м/с на 20 мкс) исход 70 и 71 м/с может измениться.

## 7. Драйвер: перезагрузка сервера и рестарт выстрела

- Задача автозапуска: `~/autostart/jobs/<имя>.job` → `guarded_run.sh` (flock, метка `.done` при rc 0) → `run_series.sh`; `@reboot ~/autostart/resume_all.sh` поднимает задачи без `.done`.
- Перед каждым запуском `guo2019_shot` драйвер ждёт: uptime > 600 с, `proj_status.sh` даёт 0 или 1, `threads_free.sh` даёт FREE ≥ 4 (опрос раз в 5 мин); через 2 мин проверяет, что выстрел занял 4 потока, иначе снимает его и останавливает серию (rc 1). `threads_free.sh` считает потоки по `OMP_NUM_THREADS` в окружении, поэтому python-драйвер запускается без этой переменной.
- Ход деления восстанавливается по готовым выстрелам на диске; оборванный выстрел продолжается с `--restart` из `checkpoint.ckpt` (`.prev` — запасной), повторное падение — с нуля.
- Монитор в чекпоинт не входит (`contact3_fix` не влит): перед `--restart` `monitor.csv` откладывается в `monitor.partK.csv`; после выстрела `monitor.csv` и `energy_log.csv` склеиваются по времени, в `shot.txt` пересчитываются `t_first_break`, `w_f_at_first_break`, `w_f_max`, `w_p_max`, добавляются `restarts`, `monitor_gaps` (исходный — `shot.raw.txt`). `n_tunnel`, max δ/d, `E_tot_max`, реакции зажимов, `steps`, `wall_s` в `shot_v*_long.txt` — только по продолжению. `shots.csv` в `_long` драйвер не правит (в нём `t_first_break` = время рестарта).
- В этих сериях рестарт был только плановый (продолжение «не определён» из финального чекпоинта): `restarts=1`, `monitor_gaps=0` во всех четырёх `_long`.
