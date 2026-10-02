# Воспроизведение: Guo 2019, 1T, уровень `scaled`, V0 = 76 м/с

## 1. Окружение

| | |
|---|---|
| Расчёт: ОС | Linux Ubuntu 26.04.1 LTS, ядро 7.0.0, 8 ядер, 15 ГБ ОЗУ |
| Компилятор | g++ 15.2.0 |
| Сборка | CMake 4.2.3, Ninja, `-DCMAKE_BUILD_TYPE=Release` |
| Потоки | OpenMP по умолчанию, 8 потоков (`OMP_NUM_THREADS` не задана) |
| Графики и анимация | Windows 11, Python 3.11, matplotlib 3.10.7, numpy 2.3.3 |

## 2. Коммит

```
e6e2d218ca2b837dfbbe6f5c64b05e99c6b286ab
Merge pull request #2 from smarry-code/weave
ветка fiber_new
```

## 3. Отличия рабочего дерева от этого коммита

Весь код Guo 2019 не закоммичен (владелец коммитит по своей команде), поэтому он целиком лежит в [`repro-diff.zip`](repro-diff.zip) (84 165 байт, sha256 `35d8073abc6742b21692c54daa8f05c974ce89e4ea7c3706097949792bb154b9`). Патчи сняты настоящим `git diff`; порядок наложения — по числовому префиксу.

| Файл архива | Что это |
|---|---|
| `01-sph-same-body-only.patch` | `include/simulation_config.hpp`, `src/input.cpp`, `include/interactions/solid_sph.hpp`, `tests/test_sph_body_isolation.cpp` (новый): опция «SPH только внутри одного тела» + тест изоляции |
| `02-fiber-frames-reader.patch` | `tests/e2e_fiber_frames.hpp`, `tests/test_e2e_fiber_frames.cpp` (новые): чтение кадров нитей |
| `03-guo2019-scene-monitor-smoke.patch` | `tests/guo2019_case.hpp`, `guo2019_scene.hpp`, `guo2019_monitor.hpp`, `test_e2e_guo2019.cpp`, `tests/data/guo2019_initial_conditions.txt` (новые): сцена, монитор, smoke-тест, уровни, включая `scaled` |
| `04-guo2019-shot-tool.patch` | `tools/guo2019/guo2019_shot.cpp` (новый): инструмент одного выстрела |
| `05-cmake-targets.patch` | `CMakeLists.txt`: цели сборки |
| `06-docs-guo2019-e2e.patch` | `fiber/guo2019-e2e.md` (новый): документация, §8.2 с отступлениями О1–О14 |
| `HEAD.txt`, `MANIFEST.md` | хэш, автор, дата, заголовок коммита; состав и команды наложения/снятия |

Наложение:

```bash
git checkout --detach e6e2d218ca2b837dfbbe6f5c64b05e99c6b286ab
mkdir /tmp/rd && python -m zipfile -e repro-diff.zip /tmp/rd && cd /tmp/rd   # либо unzip
for p in 01 02 03 04 05 06; do git -C <repo> apply --check $p-*.patch && git -C <repo> apply $p-*.patch; done
```

(пути к патчам — абсолютные или команды выполнять из каталога с патчами с `-C <repo>`). Снятие — в обратном порядке, `git apply -R` (см. `MANIFEST.md`).

Прогон 1 выполнен бинарником, собранным до появления опции `--save-t` в `guo2019_shot.cpp`; опция лишь задаёт `cfg.save_t`, расчёт не меняет.

## 4. Команды

```bash
# 1. сборка (в корне репозитория с наложенными патчами)
# -include cmath обязателен для g++ 15: include/param.hpp:241 использует sqrt без <cmath> (прежний дефект)
cmake -S . -B build-guo -G Ninja -DCMAKE_BUILD_TYPE=Release "-DCMAKE_CXX_FLAGS=-include cmath"
cmake --build build-guo -j8 --target guo2019_shot

# 2. прогон 1: кадры раз в 25 мкс (--frames=camera)
build-guo/guo2019_shot --config=1T --level=scaled --v0=76 --frames=camera --out=<dir1>

# 3. прогон 2: кадры раз в 2,5 мкс (для анимации; заодно повтор)
build-guo/guo2019_shot --config=1T --level=scaled --v0=76 --save-t=2.5e-6 --out=<dir2>

# 4. анимация и снимки (по кадрам прогона 2) и графики (по monitor.csv прогона N)
PYTHONIOENCODING=utf-8 python scripts/animate_1t.py --run-dir <dir2>/1T_scaled_v76.0_r0 --out anim_1t.gif --fps 6 --dpi 55 --max-frames 34 --zoom 9e-3 --zlim=-12,5 --snapshots auto
python scripts/plot_monitor.py --run-dir <dirN>/1T_scaled_v76.0_r0 --out-dir . --metrics data/metrics_runN.txt --thin 40 --thin-out data/monitor_runN_thin.csv
```

Графики `v_z.svg`, `energy.svg`, `deflection_breaks.svg`, `reaction.svg` построены по прогону 2 (`--run-dir <dir2>/1T_scaled_v76.0_r0 --out-dir .`); для прогона 1 тот же скрипт дал `data/metrics_run1.txt` и `data/monitor_run1_thin.csv`. Строки `shot.txt` обоих прогонов — `data/shot_run1.txt`, `data/shot_run2.txt`.

## 5. Стоимость

| Прогон | Начало–конец | Стенное время | Шагов | Время на шаг |
|---|---|---|---|---|
| 1 | 22:11:12–22:33:27 | 1334,6 с | 13 088 | 0,102 с |
| 2 | 22:36:16–22:58:58 | 1361,7 с | 13 058 | 0,104 с |

Кадры: прогон 1 — 7 + 7 пар, 7,2 МБ; прогон 2 — 68 + 68 пар, 72,2 МБ; в репозитории не хранятся. Каталог эксперимента — 2,81 МБ (GIF 1,52 МБ, PNG 0,84 МБ, SVG 0,12 МБ, данные 0,22 МБ, скрипты 0,02 МБ, архив 0,08 МБ) при ориентире ~1 МБ; исключение владельца от 2026-10-01 (до 2,9 МБ) записано в `experiments/README.md` §3.

## 6. Проверка повтора

Должно получиться: исход «пробит»; V_R = 59,1 ± 0,05 м/с (59,111 и 59,131 м/с); порвано 4 вертикальные / 0 горизонтальных нитей; первый разрыв на 103,61 мкс; E_tot,конец/E_tot,0 ≈ 0,66 (0,6607 и 0,6599). Расчёт с OpenMP (8 потоков) побитово не воспроизводится: два прогона различаются в V_R на 0,02 м/с, в максимальном прогибе на 0,9 %, в t_end на 0,4 мкс; ожидать совпадения до 5 знаков нельзя. Шаг первого разрыва в обоих прогонах одинаков (8081).
