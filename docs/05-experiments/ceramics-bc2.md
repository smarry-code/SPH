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
[`cases/bc2/regression/bc2_geometry.hpp`](../../cases/bc2/regression/bc2_geometry.hpp).

## Запуск

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

## Ключи командной строки

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

Файл [`cases/bc2/regression/initial_conditions.txt`](../../cases/bc2/regression/initial_conditions.txt)
в формате `Initial_Conditions.txt` ([справочник ключей](../03-reference/config/initial-conditions.md)).
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
