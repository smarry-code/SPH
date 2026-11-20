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

Долгие расчёты (`guo2019_shot`) можно сохранять в чекпоинты и продолжать после обрыва:
[Рестарт расчёта из чекпоинта](docs/02-how-to/05-restart-checkpoints.md).

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
├── tests/           # Модульные, интеграционные и e2e-тесты, эталонные данные
├── tools/           # Утилиты: guo2019, plain_weave, paraview, fabric_checks, скрипты
├── cmake/           # Скрипты CMake (графики, регрессия по эталону)
├── docs/            # Документация (Markdown, MkDocs Material)
├── mkdocs_hooks/    # Хуки сайта документации
├── notes/           # Рабочие заметки, ТЗ, планы и отчёты по темам (оглавление — notes/README.md)
├── experiments/     # Журналы экспериментов
├── SPH_Lee/         # Исходная Fortran-реализация (справочно)
├── Керамика/, Соударение двух стальных цилиндров/    # Исходные материалы и данные
├── CMakeLists.txt, CMakePresets.json
├── Doxyfile         # Конфигурация автогенерации API
└── mkdocs.yml       # Конфигурация сайта документации
```
