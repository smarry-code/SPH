# Керамика: BC2 — удар по карбиду бора

Карточка кейса — [`cases/bc2/README.md`](../../cases/bc2/README.md), постановка —
[`notes/contact/06-e2e-case.md`](../../notes/contact/06-e2e-case.md).

## Опыт

Пластина из пенополиуретана с летящим слоем B₄C ударяет по покоящемуся образцу B₄C с окном
LiF. Интерферометр регистрирует скорость границы «керамика / окно» на оси. Источник эталона —
Grady & Moody, SAND96-0551 (1996), опыт BC2 (в статье — тест CE6). Расчёт осесимметричный
(2D, `axis_sim`), сборка слоёв встык, без зазоров.

| Слой | Материал (`itype`) | Толщина, мм | Радиус, мм | ρ₀, кг/м³ | Скорость, м/с |
|---|---|---|---|---|---|
| Ударник | полиуретан (8) | 6,000 | 43,75 | 640 | 2210 |
| Керамика летящая | B₄C, JH-2 (10) | 3,917 | 34,6 | 2517 | 2210 |
| Керамика покоящаяся | B₄C, JH-2 (10) | 9,033 | 34,6 | 2517 | 0 |
| Окно | LiF (5) | 25,4 | 25,4 | 2640 | 0 |

Монитор — частица на оси, на тыльном слое покоящейся керамики (граница с окном). Числа —
[`cases/bc2/experiment/bc2_geometry.hpp`](../../cases/bc2/experiment/bc2_geometry.hpp).

## Запуск

Расчёты ведёт утилита `bc2_shot` ([`cases/bc2/experiment/bc2_shot.cpp`](../../cases/bc2/experiment/bc2_shot.cpp)),
по образцу `guo2019_shot`. Она берёт экспериментальную копию постановки
[`cases/bc2/experiment/`](../../cases/bc2/experiment/) и собирается с `SPH_DIM=2` независимо от общей
настройки. Общий драйвер — [`cases/common/ceramic_shot.hpp`](../../cases/common/ceramic_shot.hpp).

```bash
cmake --build build --target bc2_shot

# Полный кейс: Δx = 0,2 мм, до t_total = 3,5 мкс, импульсный контакт, чекпоинт каждые 500 шагов
OMP_NUM_THREADS=4 ./build/bc2_shot --contact-type=3 --threads=4 --out=runs/bc2 --checkpoint-every=500

# Оборвался — продолжить с теми же ключами
OMP_NUM_THREADS=4 ./build/bc2_shot --contact-type=3 --threads=4 --out=runs/bc2 --checkpoint-every=500 --restart

# Серия: все типы контакта; сходимость по Δx
OMP_NUM_THREADS=4 ./build/bc2_shot --contact-list=0,1,2,3 --threads=4 --out=runs/bc2
OMP_NUM_THREADS=4 ./build/bc2_shot --contact-type=3 --dx-list=4e-4,2e-4,1e-4 --threads=4 --out=runs/bc2-dx

# Быстрая проверка, что всё собирается и считается (Δx = 0,8 мм, 0,3 мкс)
./build/bc2_shot --delta-x=8e-4 --t-total=3e-7 --out=runs/bc2-smoke
```

Полный прогон на 4 потоках (Release) — около 7 минут.

### Выходные файлы

Каталог прогона — `<out>/bc2_dx<Δx, мкм>_c<тип>[_<tag>]/`, например `runs/bc2/bc2_dx200_c3/`:

- `log.txt` — монитор конвейера, строка на шаг (время, осевая скорость, `itype`);
- `sph_output_*.txt` — кадры раз в `save_t` (в ParaView — через `tools/paraview`);
- `curve.txt` — кривая монитора без повторов после рестарта;
- `compare.svg` — расчёт против приборного эталона;
- `shot.txt` — итог в виде `ключ=значение`: параметры, число частиц, исход
  (`completed`, `stopped`, `guard`, `error`), `t_reached`, `steps`, `wall_s` и метрики сверки
  (признаки фронта и пика расчёта и эталона, сдвиг по t50, отношения амплитуды, фронта 10–90 % и спада, признак попадания в ворота теста);
- `initial_conditions.used.txt` — копия применённого файла параметров решателя;
- `scene_id.txt`, `checkpoint.ckpt[.prev]` — при `--checkpoint-every`.

Сводка серии — `<out>/runs.csv`, строка на прогон.

Время эталона не привязано к моменту удара, поэтому сверка идёт по форме: эталон совмещается с
расчётом по t50 переднего фронта, сдвиг — справочная величина (`shift_t50`).

### Ключи `bc2_shot`

| Ключ | По умолчанию | Описание |
|---|---|---|
| `--out=<каталог>` | `runs/bc2` | Корень результатов |
| `--contact-type=N` | из конфига (0) | `contact_force_type`: 0 — выкл., 1 — пинбол, 2 — Свигл, 3 — импульсный ([ADR-007](../04-explanation/decisions/adr-007-contact-force-types.md)) |
| `--contact-list=N,M,…` | — | Серия по типам контакта |
| `--delta-x=МЕТРЫ` | 2·10⁻⁴ (полный кейс) | Шаг решётки |
| `--dx-list=М,М,…` | — | Серия по шагу решётки (сходимость) |
| `--t-total=С`, `--save-t=С` | из конфига | Время расчёта, интервал кадров |
| `--ic=<файл>` | `cases/bc2/experiment/initial_conditions.txt` | Другой файл параметров решателя |
| `--threads=N` | `OMP_NUM_THREADS` | Число потоков OpenMP |
| `--time-limit=СЕК` | 0 (без лимита) | Бюджет стенного времени на прогон |
| `--max-steps=N` | — | Останов после N шагов (штатный, исход `stopped`) |
| `--checkpoint-every=N` | выкл. | Чекпоинт `checkpoint.ckpt` каждые N шагов |
| `--restart` [`--force`] | — | Продолжить прогон из чекпоинта с теми же ключами сцены ([рестарт](../02-how-to/05-restart-checkpoints.md)) |
| `--overwrite` | — | Стереть прежний каталог прогона и считать заново (без ключа повторный запуск в занятый каталог — ошибка) |
| `--tag=<строка>` | — | Суффикс имени каталога прогона |
| `--quiet` | — | Не печатать вывод конвейера |
| `-h`, `--help` | — | Справка |

Коды выхода: 0 — все прогоны дошли до `t_total` (или до `--max-steps`), 1 — хотя бы один оборван
(сторож времени или застоя, исключение), 2 — ошибка ключей.

## Регрессия: сквозной тест `test_e2e_bc2`

Тест [`tests/test_e2e_bc2.cpp`](../../tests/test_e2e_bc2.cpp) считает замороженную постановку
[`cases/bc2/regression/`](../../cases/bc2/regression/) и проверяет её воротами; его меняют только по
[правилам регрессии](../../cases/README.md#роли-постановок).

Программа — сквозной тест [`tests/test_e2e_bc2.cpp`](../../tests/test_e2e_bc2.cpp), сборка
с `SPH_DIM=2` независимо от общей настройки.

```bash
cmake --build build --target test_e2e_bc2

# Быстрый прогон (smoke): Δx = 0,8 мм, до 0,3 мкс — проверка, что всё собирается и считается
./build/test_e2e_bc2

# Полный кейс: Δx = 0,2 мм, до t_total = 3,5 мкс, без контактных сил, сверка с эталоном
OMP_NUM_THREADS=4 ./build/test_e2e_bc2 --full --contact-type=0 --time-limit=0 --keep-workdir

# Тот же кейс с импульсным контактом и выгрузкой кривой монитора
OMP_NUM_THREADS=4 ./build/test_e2e_bc2 --full --contact-type=3 --time-limit=0 \
    --emit-curve=bc2_curve_3.txt
```

Прогон идёт во временном каталоге `$TMPDIR/sph_e2e_bc2_<режим>_<метка>/`: там `log.txt`
(монитор на каждом шаге) и кадры `sph_output_*.txt`. По умолчанию каталог удаляется после
прогона; `--keep-workdir` его сохраняет и печатает путь.

Сводный график «эталон против типов контакта 0–3» строит отдельная цель: четыре полных
прогона, затем сборка SVG в
[`notes/contact/e2e/bc2_monitor_velocity.svg`](../../notes/contact/e2e/bc2_monitor_velocity.svg)
(путь — переменная кэша `SPH_BC2_PLOT`):

```bash
cmake --build build --target bc2_plot
```

Через `ctest`: smoke-прогоны `e2e_bc2_smoke_contact0…3` идут всегда, полные
`e2e_bc2_full_contact0…3` — только в сборке с `-DSPH_E2E_SLOW_TESTS=ON`
(`ctest --test-dir build -L slow`).

### Ключи `test_e2e_bc2`

| Ключ | По умолчанию | Описание |
|---|---|---|
| `--full` | smoke | Полный кейс: Δx = 0,2 мм, `t_total` из конфигурации (3,5 мкс). Без ключа — smoke: Δx = 0,8 мм, `t_total` = 0,3 мкс |
| `--contact-type=N` | 0 | `contact_force_type` прогона, N = 0…3: 0 — выкл., 1 — пинбол, 2 — Свигл, 3 — импульсный ([ADR-007](../04-explanation/decisions/adr-007-contact-force-types.md)). Сверка с эталоном — только при N = 0, для остальных проверяется устойчивость |
| `--time-limit=СЕК` | 600 | Бюджет стенного времени на прогон; 0 — без ограничения. При исчерпании расчёт обрывается с разбором причины |
| `--keep-workdir` | выкл. | Не удалять временный каталог прогона (`log.txt`, кадры) |
| `--verbose` | выкл. | Не глушить вывод конвейера |
| `--emit-curve=<файл>` | — | Выгрузить кривую монитора (время, скорость, `itype`) для сводного графика |
| `--plot-from=<каталог>` | — | Не считать, а собрать график из `bc2_curve_0…3.txt` в каталоге; требует `--emit-plot` |
| `--emit-plot=<файл.svg>` | — | Записать SVG «расчёт против эталона» |
| `--emit-reference` | выкл. | Напечатать в stdout эталон текущего режима, снятый расчётом |
| `-h`, `--help` | — | Справка |

Шага решётки ключом не задать: Δx определяется режимом (`--full` или smoke).

## Параметры решателя

Файл [`cases/bc2/experiment/initial_conditions.txt`](../../cases/bc2/experiment/initial_conditions.txt)
(у теста — копия в `regression/`) в формате `Initial_Conditions.txt` ([справочник ключей](../03-reference/config/initial-conditions.md)).
Его разбирает тот же `parse_sph_config()`, что и основной `SPH`.

| Ключ | Значение | Смысл |
|---|---|---|
| `pa_sph` | 2 | Аппроксимация давления |
| `nnps` | 2 | Поиск соседей: связные списки |
| `sle` | 1 | Обновление длины сглаживания |
| `skf` | 1 | Ядро: кубический сплайн |
| `CFL` | 0,1 | Число Куранта |
| `epsilon` | 0,25 | XSPH (Monaghan 1992) |
| `alpha`, `beta`, `etq` | 1,0; 1,0; 0,01 | Псевдовязкость |
| `t_total` | 3,5·10⁻⁶ с | Время полного режима; в архиве опыта 3·10⁻⁶, продлено под окно эталона |
| `save_t` | 2,5·10⁻⁸ с | Интервал кадров |
| `print_step` | 50 | Шаг вывода прогресса |
| `contact_force_type` | 0 | Базовый сценарий; перекрывается `--contact-type` |
| `contact_penalty_scale` | 2,5 | Kp для типов 1–2 |
| `contact_impulse_gamma` | закомментирован (1,0) | γ для типа 3 |

Геометрия, Δx и монитор задаются кодом постановки, а не этим файлом.

## BC2 в 3D

Уменьшенная модель BC2 в полном 3D
([`bc2_mini3d_geometry.hpp`](../../cases/bc2/regression/bc2_mini3d_geometry.hpp)): материалы,
толщины и скорости BC2, радиусы уменьшены до R керамики = 4 мм. С эталоном она не сверяется;
тест [`tests/test_e2e_3d.cpp`](../../tests/test_e2e_3d.cpp) проверяет законы сохранения,
симметрию решётки и приход фронта.

```bash
cmake --build build --target test_e2e_3d
OMP_NUM_THREADS=4 ./build/test_e2e_3d --full --contact-type=0 --time-limit=0
```

| Ключ | По умолчанию | Описание |
|---|---|---|
| `--full` | smoke (до 1 мкс) | До `t_total` из конфигурации BC2 (3,5 мкс) |
| `--contact-type=N` | 0 | `contact_force_type`, N = 0…3 |
| `--no-xsph` | выкл. | `epsilon = 0`: проверка сохранения полной энергии |
| `--time-limit=СЕК` | 600 | Бюджет времени; 0 — без ограничения |
| `--keep-workdir`, `--verbose` | выкл. | Как у `test_e2e_bc2` |

## Журналы расчётов

- [`2026-09-19-bc2-discrepancy`](../../experiments/2026-09-19-bc2-discrepancy/README.md)
- [`2026-09-20-bc2-abort-contact-type2`](../../experiments/2026-09-20-bc2-abort-contact-type2/README.md)
- [`2026-09-22-bc2-contact-types`](../../experiments/2026-09-22-bc2-contact-types/README.md)
- [`2026-09-22-bc2-kp-sweep`](../../experiments/2026-09-22-bc2-kp-sweep/README.md)
- [`2026-09-22-bc2-resolution-scheme`](../../experiments/2026-09-22-bc2-resolution-scheme/README.md)
