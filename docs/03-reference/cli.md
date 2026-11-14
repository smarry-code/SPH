# Параметры запуска

Исполняемый файл `SPH` (на Windows — `SPH.exe`) не принимает аргументов командной
строки (`main()` без `argc/argv`, [`main.cpp`](../../src/main.cpp:42)). Вся
конфигурация берётся из файлов и compile-time констант. Исключения конфигурации
и импорта перехватываются в `main()`: печатается `Ошибка: <текст>`, код
возврата 1.

!!! note "`guo2019_shot` и рестарт"
    Утилита [`guo2019_shot`](../../tools/guo2019/guo2019_shot.cpp) принимает аргументы,
    в том числе ключи рестарта из чекпоинта: `--checkpoint`, `--checkpoint-every=<шаги>`,
    `--restart`, `--force`. Правила, риски `--force` и рекомендации —
    в [«Рестарт расчёта из чекпоинта»](../02-how-to/05-restart-checkpoints.md).
    Остальные утилиты с аргументами перечислены в разделе
    [«Другие исполняемые файлы»](#другие-исполняемые-файлы).

## Каталог запуска

Программа использует **текущий рабочий каталог** для размещения входных и
выходных файлов:

| Файл | Роль | Создаётся/читается |
|---|---|---|
| `Initial_Conditions.txt` | Входная конфигурация | Читается в `input()`; если файла нет, расчёт молча идёт на значениях по умолчанию |
| `sph_output_NNNNN.txt` | Кадры состояния | Создаётся в `output()` |
| `sph_gp.csv` | Накопительный дамп | Дозаписывается в `output()` |
| `log.txt` | Журнал | Создаётся в `main()`, дозаписывается |
| `../data/particles.dat` | Загрузка конфигурации | Читается при `config_input=true` |

## Консольный вывод

В процессе расчёта выводится прогресс (шаг, время, шаг по времени `dt`), а при
сохранении кадра — имя выходного файла. Опция `-DSPH_PROFILE=ON` определяет макрос
`SPH_PROFILE`, но блоков `#ifdef SPH_PROFILE` в текущем коде нет, поэтому
`force_timing.csv` программа `SPH` не пишет (это имя лишь исключается из сравнения
каталогов в тестах).

## Сценарии расчёта

Сценарий задаётся compile-time флагами в [`param.hpp`](../../include/param.hpp):

| Флаг | Описание |
|---|---|
| [`shocktube`](../../include/param.hpp:283) | Тест «ударная труба» (по умолчанию `true`) |
| [`shearcavity`](../../include/param.hpp:284) | Тест «сдвиговая каверна» |
| [`dim`](../../include/param.hpp:16) | Размерность задачи (1/2/3), из макроса `SPH_DIM` (опция CMake `SPH_DIM`) |
| [`config_input`](../../include/param.hpp:259) | Загрузка из `../data/particles.dat` |

## Другие исполняемые файлы

Цели с собственным разбором аргументов (`--ключ=значение`). Полный список ключей
печатает `-h` / `--help`.

| Цель | Исходник | Назначение и ключи |
|---|---|---|
| `guo2019_shot` | [`tools/guo2019/guo2019_shot.cpp`](../../tools/guo2019/guo2019_shot.cpp) | Серийные выстрелы Guo 2019. Обязательны `--config=A\|1T\|1T_A\|A_1T\|2T\|4T`, `--level=base\|fine\|smoke\|scaled`, `--out=<каталог>` и ровно одно из `--v0=`, `--paper-shots`, `--v0-list=a,b,c`. Необязательные: `--frames=off\|camera`, `--save-t`, `--t-max`, `--time-limit`, `--max-steps`, `--no-body-isolation`, `--mu` / `--mu-ff` / `--mu-fs`, `--repeat`, `--bracket` / `--no-bracket`, `--n-mon`, `--n-diag`, `--threads`, `--probe-frame` / `--no-probe-frame`, `--plate-nz`, `--fs-contact` / `--ff-contact` (0..4), `--fiber-inelastic`, `--eps-break`, `--fiber-E`, `--impulse-gamma`, `--impulse-push-tau` / `-depth` / `-vmax`, `--dt-scale`, `--checkpoint`, `--checkpoint-every`, `--restart`, `--force`. Только в сборках с диагностикой: `--energy-log=N` (`SPH_ENERGY_LOG`), `--contact-stats=N`, `--hold-n` (`SPH_CONTACT_STATS`); без неё ключи дают ошибку. Неизвестный аргумент или ошибка значения — код 2 |
| `fabric_checks` | [`tools/fabric_checks/fabric_checks.cpp`](../../tools/fabric_checks/fabric_checks.cpp) | Проверки лоскута ткани. Первый аргумент — подкоманда `strip`, `capstan` или `drum`, далее `--ключ=значение` (`--out`, `--tag`, `--mu`, `--threads`, `--time-limit` и ключи подкоманды) |
| `plain_weave_gen` | [`tools/plain_weave/plain_weave_gen.cpp`](../../tools/plain_weave/plain_weave_gen.cpp) | Генератор полотняного переплетения `.fib`; обязателен `--out=<файл.fib>`. [Руководство](../02-how-to/02-plain-weave-generator.md) |
| `test_e2e_*`, `test_regression*`, `test_energy_log`, `fiber_t0_dump` | `tests/` | Свои ключи (`--full`, `--contact-type=N`, `--time-limit=<с>`, `--keep-workdir`, `--verbose`, `--case=`, `--emit-*` и др.), см. [«Запуск тестов»](../02-how-to/01-run-tests.md) |

`guo2019_shot` и `fabric_checks` на время расчёта переходят в рабочий подкаталог
(`chdir`), поэтому файлы конвейера (`log.txt`, кадры) пишутся туда, а не в каталог
запуска.
