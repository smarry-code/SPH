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

## Ключи командной строки

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

Файл [`cases/bc6/regression/initial_conditions.txt`](../../cases/bc6/regression/initial_conditions.txt)
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
