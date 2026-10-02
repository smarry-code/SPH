# Воспроизведение: Guo 2019, 1T, уровень `base`, V0 = 76 м/с

## 1. Окружение

| | |
|---|---|
| Расчёт: ОС | Linux Ubuntu 26.04.1 LTS, 8 ядер |
| Компилятор | g++ 15.2.0 |
| Сборка | CMake 4.2.3, Ninja, `-DCMAKE_BUILD_TYPE=Release` |
| Потоки | OpenMP по умолчанию, 8 потоков (`OMP_NUM_THREADS` не задана) |
| Графики и анимация | Windows 11, Python 3.11, matplotlib 3.10.7 |

## 2. Коммит

```
e6e2d218ca2b837dfbbe6f5c64b05e99c6b286ab
Merge pull request #2 from smarry-code/weave
ветка fiber_new
```

## 3. Отличия рабочего дерева от этого коммита

Код тот же, что у `2026-10-01-guo2019-1t-scaled`: [`repro-diff.zip`](repro-diff.zip) (84 165 байт, sha256 `35d8073abc6742b21692c54daa8f05c974ce89e4ea7c3706097949792bb154b9`, совпадает с архивом `scaled`; между двумя прогонами код не менялся). Патчи сняты настоящим `git diff`; порядок наложения — по числовому префиксу.

| Файл архива | Что это |
|---|---|
| `01-sph-same-body-only.patch` | `include/simulation_config.hpp`, `src/input.cpp`, `include/interactions/solid_sph.hpp`, `tests/test_sph_body_isolation.cpp` (новый): опция «SPH только внутри одного тела» + тест изоляции |
| `02-fiber-frames-reader.patch` | `tests/e2e_fiber_frames.hpp`, `tests/test_e2e_fiber_frames.cpp` (новые): чтение кадров нитей |
| `03-guo2019-scene-monitor-smoke.patch` | `tests/guo2019_case.hpp`, `guo2019_scene.hpp`, `guo2019_monitor.hpp`, `test_e2e_guo2019.cpp`, `tests/data/guo2019_initial_conditions.txt` (новые): сцена, монитор, smoke-тест, уровни |
| `04-guo2019-shot-tool.patch` | `tools/guo2019/guo2019_shot.cpp` (новый): инструмент одного выстрела |
| `05-cmake-targets.patch` | `CMakeLists.txt`: цели сборки |
| `06-docs-guo2019-e2e.patch` | `fiber/guo2019-e2e.md` (новый): документация |
| `HEAD.txt`, `MANIFEST.md` | хэш, автор, дата, заголовок коммита; состав и команды наложения/снятия |

Наложение:

```bash
git checkout --detach e6e2d218ca2b837dfbbe6f5c64b05e99c6b286ab
mkdir /tmp/rd && python -m zipfile -e repro-diff.zip /tmp/rd && cd /tmp/rd   # либо unzip
for p in 01 02 03 04 05 06; do git -C <repo> apply --check $p-*.patch && git -C <repo> apply $p-*.patch; done
```

(пути к патчам — абсолютные или команды выполнять из каталога с патчами с `-C <repo>`). Снятие — в обратном порядке, `git apply -R` (см. `MANIFEST.md`).

Скрипт `scripts/animate_1t.py` в этом каталоге отличается от версии в `scaled`: добавлены опции `--yarn-stride`, `--vmax`, `--lw`, `--cmap`, `--ball-alpha`, `--snap-dpi-factor`, обрезка по `--zoom`, векторизованный разбор кадров и исправление маски прореживания нитей на каждом кадре (число связей падает после разрывов). Расчёт скрипт не затрагивает.

## 4. Команды

```bash
# 1. сборка (в корне репозитория с наложенными патчами)
# -include cmath обязателен для g++ 15: include/param.hpp:241 использует sqrt без <cmath> (прежний дефект)
cmake -S . -B build-guo -G Ninja -DCMAKE_BUILD_TYPE=Release "-DCMAKE_CXX_FLAGS=-include cmath"
cmake --build build-guo -j8 --target guo2019_shot

# 2. прогон: кадры раз в 5 мкс (для анимации)
build-guo/guo2019_shot --config=1T --level=base --v0=76 --save-t=5e-6 --out=<dir>

# 3. анимация и снимки (по кадрам прогона)
PYTHONIOENCODING=utf-8 python scripts/animate_1t.py --run-dir <dir>/1T_base_v76.0_r0 --out anim_1t_base.gif --zoom 30e-3 --yarn-stride 2 --dpi 55 --zlim=-45,12 --vmax 0.5 --elev 30 --max-frames 62 --fps 6 --snapshots auto
python scripts/animate_1t.py --run-dir <dir>/1T_base_v76.0_r0 --out full_field.gif --no-gif --yarn-stride 4 --dpi 60 --zlim=-45,12 --vmax 0.5 --snapshots 81,122

# 4. графики (по monitor.csv прогона) и прореженная выборка
python scripts/plot_monitor.py --run-dir <dir>/1T_base_v76.0_r0 --out-dir . --metrics data/metrics.txt --thin 100 --thin-out data/monitor_thin.csv
```

Номера снимков: `--snapshots auto` для `--max-frames 62` дал кадры анимации 0, 41, 61 (полные кадры 0, 82, 122); `--snapshots 81,122` — полные кадры 81 (405,04 мкс) и 122 (610,03 мкс). Файлы `full_field.gif` при `--no-gif` не создаются; снимки сохранены как `full_field_frame81.png`, `full_field_frame122.png`. Сводка `shot.txt` и лог `shot.log` скопированы из каталога вывода `guo2019_shot`; `guo2019_fabric_layout.txt` — оттуда же.

## 5. Стоимость

| Прогон | Начало–конец | Стенное время | Шагов | Время на шаг |
|---|---|---|---|---|
| 1 | 23:24:51–01:46:31 | 8499,9 с (2,36 ч) | 33 302 | 0,255 с |

Монитор — 956 с (11,25 % времени). Кадры: 123 + 123 пар, 749 МБ; в репозитории не хранятся. Каталог эксперимента — 3,09 МБ (GIF 1,84 МБ, PNG 0,85 МБ, SVG 0,13 МБ, данные 0,12 МБ, скрипты 0,03 МБ, архив 0,08 МБ) при ориентире ~1 МБ; исключение владельца от 2026-10-02 (до 3,2 МБ) записано в `experiments/README.md` §3.

## 6. Проверка повтора

Должно получиться: исход «пробит»; V_R = 37,3 ± ~0,1 м/с (37,288); первый разрыв около 404 мкс (403,95); порвано 11 вертикальных / 1 горизонтальная нить; w_f при первом разрыве около 26,2 мм; E_tot,конец/E_tot,0 ≈ 0,32 (0,319). Расчёт с OpenMP (8 потоков) побитово не воспроизводится; разброс на `base` не измерен (один прогон), на `scaled` V_R между двумя прогонами расходился на 0,03 %, поэтому допуск ±0,1 м/с — ориентир, а не измерение.
