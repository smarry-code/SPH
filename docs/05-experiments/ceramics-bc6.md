# Керамика: BC6 — откольный опыт на карбиде бора

Карточка кейса — [`cases/bc6/README.md`](../../cases/bc6/README.md), постановка —
[`notes/contact/06-e2e-bc6.md`](../../notes/contact/06-e2e-bc6.md).

## Опыт

Ударник из PMMA на подложке из пены ударяет по образцу B₄C с окном из PMMA. Сжатие упругое
(около 3,6 ГПа, ниже HEL). Растяжение создаёт низкий импеданс окна, поэтому опыт проверяет
отклик JH-2 на растяжение. VISAR регистрирует скорость границы «керамика / окно»: убывающую
лесенку из трёх ступеней. Эталон — оцифровка опубликованного профиля опыта BC6. Расчёт
осесимметричный (2D, `axis_sim`), сборка встык, удар при t = 0.

| Тело | Материал (`itype`) | Толщина, мм | Радиус, мм | ρ₀, кг/м³ | Скорость, м/с |
|---|---|---|---|---|---|
| Подложка | пена, полиуретан 160 (6) | 8,000 | 43,75 | 160 | 913 |
| Ударник | PMMA (9) | 2,004 | 43,75 | 1186 | 913 |
| Образец | B₄C, JH-2 (10) | 10,526 | 38,05 | 2506 | 0 |
| Окно | PMMA (9) | 24,200 | 38,05 | 1186 | 0 |

Мониторы — на тыльном слое образца: основной на оси, контрольный на r = 12,5 мм. Числа опыта
и допуски приёмки (приход фронта ±0,05 мкс, плато 1 — 170–190 м/с, отношение ступеней
0,82 ± 0,02, период 0,88 ± 0,15 мкс) —
[`cases/bc6/regression/bc6_case.hpp`](../../cases/bc6/regression/bc6_case.hpp).

## Запуск

Расчёты ведёт утилита `bc6_shot` ([`cases/bc6/experiment/bc6_shot.cpp`](../../cases/bc6/experiment/bc6_shot.cpp)),
по образцу `guo2019_shot`. Она берёт экспериментальную копию постановки
[`cases/bc6/experiment/`](../../cases/bc6/experiment/) и собирается с `SPH_DIM=2` независимо от общей
настройки. Общий драйвер — [`cases/common/ceramic_shot.hpp`](../../cases/common/ceramic_shot.hpp).

```bash
cmake --build build --target bc6_shot

# Полный кейс: Δx = 0,2 мм, до t_total = 4,3 мкс, импульсный контакт, чекпоинт каждые 500 шагов
OMP_NUM_THREADS=4 ./build/bc6_shot --contact-type=3 --threads=4 --out=runs/bc6 --checkpoint-every=500

# Оборвался — продолжить с теми же ключами
OMP_NUM_THREADS=4 ./build/bc6_shot --contact-type=3 --threads=4 --out=runs/bc6 --checkpoint-every=500 --restart

# Сходимость по Δx
OMP_NUM_THREADS=4 ./build/bc6_shot --contact-type=3 --dx-list=4e-4,2e-4,1e-4 --threads=4 --out=runs/bc6-dx

# Быстрая проверка (Δx = 0,8 мм, 1 мкс)
./build/bc6_shot --delta-x=8e-4 --t-total=1e-6 --out=runs/bc6-smoke
```

Полный прогон на 4 потоках (Release) — около 6–8 минут.

### Выходные файлы

Каталог прогона — `<out>/bc6_dx<Δx, мкм>_c<тип>[_<tag>]/`, например `runs/bc6/bc6_dx200_c3/`:

- `log.txt` — монитор конвейера, строка на шаг (время, осевая скорость, `itype`);
- `sph_output_*.txt` — кадры раз в `save_t` (в ParaView — через `tools/paraview`);
- `curve.txt` — кривая монитора без повторов после рестарта;
- `compare.svg` — расчёт против приборного эталона;
- `shot.txt` — итог в виде `ключ=значение`: параметры, число частиц, исход
  (`completed`, `stopped`, `guard`, `error`), `t_reached`, `steps`, `wall_s` и метрики сверки
  (лесенка расчёта и эталона — фронт по 10 % и 50 %, ступени, провалы, отношения ступеней, периоды; признаки допусков `front_ok`, `plateau1_ok`, `ratio12_ok`, `period1_ok` и др.; отношение амплитуд к эталону);
- `initial_conditions.used.txt` — копия применённого файла параметров решателя;
- `scene_id.txt`, `checkpoint.ckpt[.prev]` — при `--checkpoint-every`.

Сводка серии — `<out>/runs.csv`, строка на прогон.

Допуски — те же, что у теста (`bc6_case.hpp`), но утилита ими не проваливается: признак `*_ok`
только сообщает, попало ли значение в допуск. Эталон на графике совмещён по t10 фронта
(`ref_shift_t10`).

### Ключи `bc6_shot`

| Ключ | По умолчанию | Описание |
|---|---|---|
| `--out=<каталог>` | `runs/bc6` | Корень результатов |
| `--contact-type=N` | из конфига (0) | `contact_force_type`: 0 — выкл., 1 — пинбол, 2 — Свигл, 3 — импульсный ([ADR-007](../04-explanation/decisions/adr-007-contact-force-types.md)) |
| `--contact-list=N,M,…` | — | Серия по типам контакта |
| `--delta-x=МЕТРЫ` | 2·10⁻⁴ (полный кейс) | Шаг решётки |
| `--dx-list=М,М,…` | — | Серия по шагу решётки (сходимость) |
| `--t-total=С`, `--save-t=С` | из конфига | Время расчёта, интервал кадров |
| `--ic=<файл>` | `cases/bc6/experiment/initial_conditions.txt` | Другой файл параметров решателя |
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

## Регрессия: сквозной тест `test_e2e_bc6`

Тест [`tests/test_e2e_bc6.cpp`](../../tests/test_e2e_bc6.cpp) считает замороженную постановку
[`cases/bc6/regression/`](../../cases/bc6/regression/) и проверяет её воротами; его меняют только по
[правилам регрессии](../../cases/README.md#роли-постановок).

Программа — сквозной тест [`tests/test_e2e_bc6.cpp`](../../tests/test_e2e_bc6.cpp), сборка
с `SPH_DIM=2` независимо от общей настройки.

```bash
cmake --build build --target test_e2e_bc6

# Быстрый прогон (smoke): Δx = 0,8 мм, до 1 мкс; физика не проверяется
./build/test_e2e_bc6

# Полный кейс: Δx = 0,2 мм, до t_total = 4,3 мкс, сверка лесенки с эталоном
OMP_NUM_THREADS=4 ./build/test_e2e_bc6 --full --contact-type=0 --time-limit=0 --keep-workdir

# Исследование сходимости: свой шаг решётки (контрольные числа и физические проверки пропускаются)
OMP_NUM_THREADS=4 ./build/test_e2e_bc6 --full --delta-x=1.0e-4 --time-limit=0 \
    --emit-curve=bc6_dx010.txt
```

Прогон идёт во временном каталоге (`log.txt` — основной монитор на каждом шаге, кадры
`sph_output_*.txt` раз в 50 нс). По умолчанию каталог удаляется; `--keep-workdir` его
сохраняет и печатает путь.

Сводный график «эталон против типов контакта 0–3» — четыре полных прогона и сборка SVG в
[`notes/contact/e2e/bc6_monitor_velocity.svg`](../../notes/contact/e2e/bc6_monitor_velocity.svg)
(путь — переменная кэша `SPH_BC6_PLOT`, бюджет прогона — `SPH_BC6_RUN_TIME_LIMIT`, 1800 с):

```bash
cmake --build build --target bc6_plot
```

Через `ctest`: smoke-прогоны `e2e_bc6_smoke_contact0…3` идут всегда, полные
`e2e_bc6_full_contact0…3` — только с `-DSPH_E2E_SLOW_TESTS=ON`.

Эталон [`regression/reference.txt`](../../cases/bc6/regression/reference.txt) этой программой
не генерируется: его строит из приборной записи
[`source/reference-curve.txt`](../../cases/bc6/source/reference-curve.txt) скрипт
[`tools/bc6_make_reference.py`](../../tools/bc6_make_reference.py), запускаемый вручную.

### Ключи `test_e2e_bc6`

| Ключ | По умолчанию | Описание |
|---|---|---|
| `--full` | smoke | Полный кейс: Δx = 0,2 мм, `t_total` из конфигурации (4,3 мкс). Без ключа — smoke: Δx = 0,8 мм, `t_total` = 1 мкс |
| `--contact-type=N` | 0 | `contact_force_type` прогона, N = 0…3 ([ADR-007](../04-explanation/decisions/adr-007-contact-force-types.md)). Сверка формы кривой с эталоном — только при N = 0 |
| `--delta-x=МЕТРЫ` | по режиму | Свой шаг решётки (сходимость). Контрольные числа и физические проверки при этом пропускаются |
| `--time-limit=СЕК` | 600 | Бюджет стенного времени; 0 — без ограничения |
| `--keep-workdir` | выкл. | Не удалять временный каталог прогона |
| `--verbose` | выкл. | Не глушить вывод конвейера |
| `--emit-curve=<файл>` | — | Выгрузить кривую монитора (время, скорость, `itype`) |
| `--plot-from=<каталог>` | — | Не считать, а собрать график из `bc6_curve_0…3.txt`; требует `--emit-plot` |
| `--emit-plot=<файл.svg>` | — | Записать SVG «расчёт против эталона» |
| `-h`, `--help` | — | Справка |

## Параметры решателя

Файл [`cases/bc6/experiment/initial_conditions.txt`](../../cases/bc6/experiment/initial_conditions.txt)
(у теста — копия в `regression/`)
([справочник ключей](../03-reference/config/initial-conditions.md)).

| Ключ | Значение | Смысл |
|---|---|---|
| `pa_sph`, `nnps`, `sle`, `skf` | 2, 2, 1, 1 | Как у BC2 |
| `CFL` | 0,1 | Число Куранта |
| `epsilon` | 0,25 | XSPH |
| `alpha`, `beta`, `etq` | 1,0; 1,0; 0,01 | Псевдовязкость |
| `t_total` | 4,3·10⁻⁶ с | Время полного режима; в архиве 3·10⁻⁶, продлено, чтобы захватить обе ступени реверберации |
| `save_t` | 5·10⁻⁸ с | Интервал кадров (86 кадров) |
| `print_step` | 50 | Шаг вывода прогресса |
| `contact_force_type` | 0 | Базовый сценарий; перекрывается `--contact-type` |
| `contact_penalty_scale` | закомментирован (1,0) | Kp для типов 1–2; на профиль BC6 не влияет |
| `contact_impulse_gamma` | 0,5 | γ для типа 3; при 1,0 кривая идёт пилой с периодом 2dt |

## Журналы расчётов

- [`2026-09-23-bc6-art-visc`](../../experiments/2026-09-23-bc6-art-visc/README.md)
- [`2026-09-23-bc6-butt-joint`](../../experiments/2026-09-23-bc6-butt-joint/README.md)
- [`2026-09-23-bc6-resolution`](../../experiments/2026-09-23-bc6-resolution/README.md)
- [`2026-09-23-bc6-xsph`](../../experiments/2026-09-23-bc6-xsph/README.md)
- [`2026-09-24-bc6-contact-dx010`](../../experiments/2026-09-24-bc6-contact-dx010/README.md)
- [`2026-11-21-bc6-shot-contact3`](../../experiments/2026-11-21-bc6-shot-contact3/README.md)
