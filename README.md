# SPH — Симулятор гидродинамики сглаженных частиц

Вычислительный код на C++20/OpenMP для моделирования ударно-волновых процессов,
высокоскоростного соударения и поведения керамических/металлических материалов
методом сглаженных частиц (SPH).

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

Справочник по утилите — [`docs/03-reference/guo2019-shot.md`](docs/03-reference/guo2019-shot.md):
там ключи мишени (`A`, `1T_A`, `<k>T`), трения, контакта, стёжки и мозаики. Если расчёт
оборвался, его можно продолжить с того же места, добавив к тем же ключам `--restart`, см.
[Рестарт расчёта из чекпоинта](docs/02-how-to/05-restart-checkpoints.md).

Чтобы поменять постановку под расчёт, правьте
[`cases/guo2019/experiment/`](cases/guo2019/experiment/), а не `regression/`.
Замороженную копию в `regression/` берёт `ctest`, правила —
[`cases/README.md`](cases/README.md#роли-постановок). Журналы прошлых серий собраны в
[`experiments/`](experiments/README.md).

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
├── tools/           # Утилиты: plain_weave, paraview, fabric_checks, скрипты (guo2019_shot — cases/guo2019/experiment/)
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
