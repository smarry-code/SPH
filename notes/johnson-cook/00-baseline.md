# J0. Базовые результаты до внедрения Джонсона-Кука

- Дата: 2026-11-03
- Коммит: `6f1b63f` (ветка `fiber_johnson_cook` до первых JC-коммитов)
- Хост: `claude-compute` (2 ядра, 7.7 ГиБ), g++ 15.2.0, CMake 4.2.3, Ninja
- Каталог на сервере: `~/SPH_jc`, логи `ctest_buildN.log`, `build_buildN.log`

## Сборки

Все сборки Release, сборка `-j2`, `SPH_E2E_SLOW_TESTS` выключен.

| Сборка | Флаги | Время сборки, с | Тестов | Прошло | Время ctest, с |
|---|---|---|---|---|---|
| build1 | `-DSPH_DIM=1` | 966 | 100 | 100 | 348.4 |
| build2 | `-DSPH_DIM=2` | 968 | 100 | 100 | 352.6 |
| build3 | `-DSPH_DIM=3` | 971 | 100 | 100 | 349.5 |
| build3e | `-DSPH_DIM=3 -DSPH_ENERGY_LOG=ON` | 986 | 100 | 100 | 348.9 |

Упавших тестов нет.

## Замечание о размерностях

Опции `SPH_DIM` и `SPH_ENERGY_LOG` влияют только на цели `SPH`, `guo2019_shot`
и `fabric_checks`. Каждая тестовая цель задаёт `SPH_DIM` (и при необходимости
`SPH_ENERGY_LOG`) сама. Поэтому набор тестов во всех четырёх сборках одинаков.
Для проверки изменений достаточно полного `ctest` в одной сборке. В остальных
сборках достаточно скомпилировать `SPH`.

## Самые долгие тесты (build1)

| Тест | с |
|---|---|
| e2e_clamp_ES3_pc_smoke | 48.2 |
| e2e_clamp_ES2_pc_smoke | 43.4 |
| energy_log_test | 36.4 |
| e2e_clamp_ES1_pc_smoke | 16.4 |
| e2e_clamp_ES1_euler_smoke | 16.3 |
| e2e_clamp_ES3_euler_smoke | 16.2 |
| e2e_clamp_ES2_euler_smoke | 14.1 |
| e2e_clamp_EF2_euler_smoke | 13.4 |
| e2e_guo2019_smoke_4T | 12.3 |
| e2e_clamp_ES1_pc_smoke_t4 | 9.7 |

## Команды

```
cmake -S ~/SPH_jc -B buildN -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPH_DIM=N [-DSPH_ENERGY_LOG=ON]
cmake --build buildN -j2
ctest --test-dir buildN --output-on-failure -j1
```
