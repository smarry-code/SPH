# SPH — Симулятор гидродинамики сглаженных частиц

Вычислительный код на C++20/OpenMP для моделирования ударно-волновых процессов,
высокоскоростного соударения и поведения керамических/металлических материалов
методом сглаженных частиц (SPH).

Репозиторий: [github.com/smarry-code/SPH](https://github.com/smarry-code/SPH).

## Возможности

- **Сплошная среда (Solid).** Прочность по Мизесу, JH-2 для керамики, Джонсон–Кук
  для металлов; размерность 1, 2 или 3 выбирается при сборке (`SPH_DIM`), в 2D
  доступна осесимметричная постановка.
- **Нити и ткани (Fiber).** Нити из бусин со связями и пластичностью, генератор
  полотняного переплетения, пакеты слоёв, стёганое одеяло и мозаичная ткань со
  стежками; только `SPH_DIM=3`.
- **Контакты.** Частица–частица (пинбол, Свигл, импульсный закон), нить–тело и
  нить–нить с трением (законы 1–3 нитей — также как проекция скорости); задаются
  ключами `Initial_Conditions.txt`.
- **Прочее.** Закрепления (Clamp), импорт LS-DYNA `.k`, чекпоинты и рестарт,
  лог энергии, конвертация кадров в VTK для ParaView.

Полный список — в [`docs/index.md`](docs/index.md).

## Быстрый старт

Требуется Clang (C++20), CMake 3.21+, Ninja, OpenMP (libomp).

```bash
# Настройка и сборка
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPH_DIM=1
cmake --build build

# Запуск: аргументов нет, читается Initial_Conditions.txt из текущего каталога
./build/SPH

# Тесты
ctest --test-dir build --output-on-failure
```

Готовые конфигурации (`clang-windows`, `clang-windows-debug`,
`clang-windows-native`) заданы в `CMakePresets.json`:
`cmake --preset clang-windows`. Подробности — в
[Быстром старте](docs/01-tutorials/01-quick-start.md), параметры запуска — в
[`docs/03-reference/cli.md`](docs/03-reference/cli.md), тесты — в
[Запуске тестов](docs/02-how-to/01-run-tests.md).

### Пример: эксперимент Guo 2019

Серийные выстрелы стального шара по ткани и пластине Guo 2019 запускает утилита
`guo2019_shot`. Она собирается из экспериментальной копии постановки
[`cases/guo2019/experiment/`](cases/guo2019/experiment/) и всегда в 3D, при любом `SPH_DIM`.

```bash
# Сборка утилиты
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target guo2019_shot

# Один пробный выстрел: 1 слой ткани, грубая сцена, 147 м/с, 4 потока
OMP_NUM_THREADS=4 ./build/guo2019_shot --config=1T --level=smoke --v0=147 \
    --out=runs/guo-1T-smoke --threads=4

# Серия по скоростям статьи с поиском V_BL и чекпоинтами
OMP_NUM_THREADS=4 ./build/guo2019_shot --config=1T --level=base --paper-shots \
    --out=runs/guo-1T-base --threads=4 --checkpoint-every=5000 --frames=camera
```

Каждый выстрел пишет результаты в свой каталог `<out>/<cfg>_<level>_v<V0>_r<повтор>/`,
например `runs/guo-1T-smoke/1T_smoke_v147.0_r0/`:

- `shot.txt` — исход и итоговые величины;
- `monitor.csv` — монитор по шагам;
- кадры `sph_output_*.txt` и `fiber_output_*.txt`;
- `log.txt` — журнал расчёта.

Сводки серии `shots.csv` и `vbl.txt` (баллистический предел $V_{BL}$) лежат в `<out>`.

Керамику BC2 и BC6 так же считают утилиты `bc2_shot` и `bc6_shot`
(`--contact-type=3 --out=runs/bc2 --threads=4`). Все опыты с запуском и ключами — раздел
[«Эксперименты»](docs/05-experiments/index.md). Справочник по утилите — [`docs/03-reference/guo2019-shot.md`](docs/03-reference/guo2019-shot.md):
там ключи мишени (`A`, `1T_A`, `<k>T`), трения, контакта, стёжки и мозаики. Если расчёт
оборвался, его можно продолжить с того же места, добавив к тем же ключам `--restart`, см.
[Рестарт расчёта из чекпоинта](docs/02-how-to/05-restart-checkpoints.md).

Чтобы поменять постановку под расчёт, правьте
[`cases/guo2019/experiment/`](cases/guo2019/experiment/), а не `regression/`.
Замороженную копию в `regression/` берёт `ctest`, правила —
[`cases/README.md`](cases/README.md#роли-постановок). Журналы прошлых серий собраны в
[`experiments/`](experiments/README.md).

### Пример: геометрия из LS-DYNA `.k`

Произвольную геометрию можно взять из keyword-файла `.k`, например экспорта
LS-PrePost. Из `.k` читаются узлы, SPH-частицы (`*ELEMENT_SPH`), нити
(цепочки `*ELEMENT_BEAM`) и закрепления (`*BOUNDARY_SPC_*`); каждая деталь
(PID) — отдельное тело. Материал и начальную скорость детали задаёт таблица
деталей, контакт — ключи `Initial_Conditions.txt`; `*MAT_*`, `*CONTACT_*`,
`*INITIAL_VELOCITY*` пропускаются.

Пример: стальной ударник (PID 1) летит на нить (PID 3) над алюминиевой
пластиной (PID 2), модель в миллиметрах. В каталоге запуска лежат три файла.

`model.k`:

```text
*KEYWORD
*PART
Ударник
1
*PART
Пластина
2
*PART
Нить
3
*NODE
1, 0.0, 0.0, 3.0
2, 0.0, 0.0, 2.0
11, -1.0, 0.0, 0.0
12, 0.0, 0.0, 0.0
13, 1.0, 0.0, 0.0
21, -1.0, 0.0, 1.0
22, 0.0, 0.0, 1.0
23, 1.0, 0.0, 1.0
*ELEMENT_SPH
1, 1, -1.0, 2
11, 2, -1.0, 13
*ELEMENT_BEAM
1, 3, 21, 22
2, 3, 22, 23
*BOUNDARY_SPC_NODE
11, 0, 1, 1, 1, 0, 0, 0
13, 0, 1, 1, 1, 0, 0, 0
*END
```

`model.parts` — физика деталей (`itype` — номер материала проекта):

```text
SOLID_PART 1 12 7830 v=0 0 -200   # ударник: сталь 4340, 200 м/с против z
SOLID_PART 2 4  2700              # пластина: алюминий
# нить: E, nu, d, rho, eps_break, mu_ff, mu_fs
MATERIAL 0 1.0e11 0.3 2.0e-4 1440 0.03 0.1 0.2
FIBER_PART 3 0
```

`Initial_Conditions.txt` — обычные параметры решателя (за основу можно взять
[`cases/guo2019/experiment/initial_conditions.txt`](cases/guo2019/experiment/initial_conditions.txt))
и ключи импорта, контакта и времени:

```text
k_file = model.k
k_mode = solid_fiber
k_parts_file = model.parts
k_length_scale = 1.0e-3
contact_force_type = 3
fiber_solid_contact_type = 3
t_total = 2.0e-6
save_t = 1.0e-6
```

```bash
# Нити требуют 3D-сборки
cmake -S . -B build3 -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPH_DIM=3
cmake --build build3 --target SPH
cd <каталог с тремя файлами> && /путь/к/build3/SPH
```

Перед расчётом печатается сводка импорта (число частиц и нитей, габарит в
метрах), результаты — кадры `sph_output_*.txt`, `fiber_output_*.txt` и
`log.txt` в том же каталоге. Что где задаётся (тела, нити, контакт) —
[«Тела, нити и контакт»](docs/03-reference/file-formats/k-file.md#тела-нити-и-контакт),
пошагово — [Из LS-PrePost в расчёт](docs/02-how-to/03-ls-prepost-to-calculation.md).

## Документация

Сайт документации (MkDocs Material) — каталог [`docs/`](docs/):

```bash
pip install -r requirements-docs.txt
mkdocs serve        # http://127.0.0.1:8000
```

Математические модели — [`docs/04-explanation/models/`](docs/04-explanation/models/index.md),
ключи конфигурации — [`docs/03-reference/config/initial-conditions.md`](docs/03-reference/config/initial-conditions.md),
конвейер расчёта — [`docs/04-explanation/architecture/pipeline-overview.md`](docs/04-explanation/architecture/pipeline-overview.md).

## Структура репозитория

```
SPH/
├── src/, include/   # Код C++20 (include/: rhs, interactions, fiber, ode, strength, k_input)
├── tests/           # Модульные, интеграционные и e2e-тесты (постановки опытов — cases/*/regression/)
├── tools/           # Утилиты: plain_weave, paraview, fabric_checks, скрипты (guo2019_shot, bc2_shot, bc6_shot — cases/<кейс>/experiment/)
├── cmake/           # Скрипты CMake (графики, регрессия по эталону)
├── docs/            # Документация (Markdown, MkDocs Material)
├── mkdocs_hooks/    # Хуки сайта документации
├── notes/           # Рабочие заметки, ТЗ, планы и отчёты по темам (оглавление — notes/README.md)
├── experiments/     # Журналы экспериментов
├── cases/           # Опытные кейсы: исходные данные, постановки для регрессии и экспериментов (cases/README.md)
├── CMakeLists.txt, CMakePresets.json
├── Doxyfile         # Конфигурация автогенерации API
└── mkdocs.yml       # Конфигурация сайта документации
```
