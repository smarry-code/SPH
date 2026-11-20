# Запуск тестов и бенчмарков

В проекте настроены тесты через CTest (регистрация — после [`enable_testing()`](../../CMakeLists.txt:72)).

## Запуск всех тестов

```bash
ctest --test-dir build --output-on-failure
```

## Состав тестов

| Тест | Исполняемый файл | Что проверяет |
|---|---|---|
| `executor_test` | `test_executor` | `CpuExecutor` ([`include/rhs/cpu_executor.hpp`](../../include/rhs/cpu_executor.hpp)): `parallel_for`, `parallel_reduce`, `parallel_scan`, `parallel_sort` против последовательных эталонов при 1, 2 и 7 потоках |
| `spatial_index_test` | `test_spatial_index` | `SpatialIndex<Dim>`: найденные пары совпадают с полным перебором («свои» и «чужие» типы, разная плотность) |
| `contributions_test` | `test_contributions` | `ContributionWriter`, `reduce_contributions` и `RhsWorkspace` (`include/rhs/`) |
| `axisym_test` | `test_axisym` | Осесимметричная ветка (`axis_sim`, dim = 2): соглашение об индексах девиатора, бесследовость, равенство S_rr и S_theta, нулевой радиальный источник при одноосном сжатии |
| `particle_system_test` | `test_particle_system` | Тождество моста AoS↔SoA (точность 1e-12) |
| `contact_test` | `test_contact` | Контактные силы «частица-частица» (ADR-007): аналитика типов 1-3, баланс импульса, ограничитель шага |
| `e2e_bc2_smoke_contact0..3` | `test_e2e_bc2 --contact-type=N` | Сквозной прогон кейса BC2 на огрублённой сетке: геометрия, живучесть конвейера. По тесту на `contact_force_type`; сверка с эталоном — только в `contact0` |
| `e2e_bc2_full_contact0..3` | `test_e2e_bc2 --full --contact-type=N` | Тот же кейс на рабочей сетке. Метка `slow`, по умолчанию **не регистрируются** |
| `e2e_bc6_smoke_contact0..3` | `test_e2e_bc6 --contact-type=N` | Сквозной прогон кейса BC6 (откол на B4C, ударник и окно PMMA) на огрублённой сетке: геометрия, живучесть конвейера, самопроверки компаратора и извлекателя лесенки. Физические проверки там пропускаются — зазор удара 1,0 мм лежит внутри носителя ядра |
| `e2e_bc6_full_contact0..3` | `test_e2e_bc6 --full --contact-type=N` | Тот же кейс на рабочей сетке: приход фронта, уровни ступеней, отношения, период лесенки, `D = 0` до разгрузки, контакт на всех трёх границах. Метка `slow`, по умолчанию **не регистрируются** |
| `regression_selfcheck` | `test_regression --selfcheck` | Самопроверка компаратора снимков без расчёта: ловит возмущение выше допуска и называет частицу, пропускает ниже допуска, ловит NaN, расхождение `itype`, `ntotal`, числа шагов; круг записи и чтения эталона |
| `regression_bc{2,6}_contact0..3` | `test_regression --case=bcN --contact-type=N` | Регрессия: все поля всех частиц после 1 и 100 шагов сверяются с эталоном, снятым этим же кодом (`cases/<кейс>/regression/snapshots/`). Метка `regression` |
| `kinematics_3d_test` | `test_kinematics_3d` | Ветка `dim = 3` силового расчёта на решётке 11³: простой сдвиг по трём парам осей, объёмное и одноосное сжатие, жёсткое вращение (производная Яуманна: dS/dt = R·S − S·R, след и S:S неизменны), покой при однородном и линейном поле напряжений |
| `e2e_3d_smoke_contact0..3` | `test_e2e_3d --contact-type=N` | Сквозной прогон уменьшенного BC2 в полном 3D до 1 мкс: здоровье, сохранение массы и импульса, энергия не растёт, симметрия решётки x → −x, y → −y, x ↔ y, приход фронта к монитору |
| `e2e_3d_energy_noxsph` | `test_e2e_3d --no-xsph` | Тот же прогон при `epsilon = 0`: полная энергия сохраняется в пределах 1 % |
| `e2e_3d_full_contact0..3` | `test_e2e_3d --full --contact-type=N` | Тот же кейс до 3,5 мкс. Метка `slow`, по умолчанию **не регистрируются** |
| `regression_selfcheck_3d` | `test_regression_3d --selfcheck` | Самопроверка компаратора в сборке `dim = 3` (19 полей) |
| `regression_bc2mini3d_contact0..3` | `test_regression_3d --case=bc2mini3d --contact-type=N` | Регрессия уменьшенного BC2 в 3D: все поля всех частиц после 1 и 50 шагов. Метки `regression` и `3d` |
| `checkpoint_test` | `test_checkpoint` | Рестарт из чекпоинта ([руководство](05-restart-checkpoints.md)): побитовое совпадение «2N шагов» и «N + рестарт + N» для предиктора–корректора и Эйлера, в том числе с разрывом связей до чекпоинта и со стежком; JC-сцена с нагревом; порча, обрезание, неверная сигнатура и чужая версия формата; политика изменения настроек и `--force`; период и ротация `.prev`; таблица тёплого старта проекции и блок `USER`; рестарт сцены с контактом нитей при заполненной таблице и независимость от числа потоков (1 и 4). Метка `3d` |
| `checkpoint_energy_test` | `test_checkpoint_energy` | То же с `SPH_ENERGY_LOG`: CSV лога энергии после рестарта совпадает с непрерывным (один заголовок, без дублей строк) |
| `checkpoint_cli_test` | `tests/checkpoint_cli_test.sh` | Ключи `guo2019_shot`: `--force` без `--restart` и `--checkpoint-every=0` дают код 2, `--restart` без файла — ошибку, разбитый на два запуска выстрел даёт те же итоговые кадры. Требует собранный `guo2019_shot` |
| `stitch_kt_cli_test` | `tests/stitch_kt_cli_test.sh` | Стёганое одеяло в `guo2019_shot`: 3T smoke со стёжкой `grid` (tile 4) — 144 стежка, рестарт побайтно равен непрерывному, рестарт с другой стёжкой отклонён, `grid` на 1T — код 2, `--stitch-patch=3` ≡ `--stitch-tile=4`, 11 нитей: `P = 4` — код 2, `P = 5` — 88 стежков |
| `mosaic_kt_cli_test` | `tests/mosaic_kt_cli_test.sh` | Мозаичная ткань в `guo2019_shot` (2T smoke, `--random-cut=4 --layers-per-cut=2`): стежки в слое и сквозные, `cut.json`, столбцы стежков в `monitor.csv`, рестарт побайтно, рестарт с другим `--cut-seed` отклонён, ошибки ключей — код 2 |

!!! note "Почему по тесту на тип контакта"
    Прежде один тест делал базовый прогон и следом три сценария §6. Обрыв любого
    из них прекращал разбор остальных, а в отчёте CTest всё сводилось к одной
    строке «провален», по которой нельзя сказать, какой именно тип сломался.
    Отдельные тесты падают независимо, идут параллельно под `ctest -j` и
    называют виновника именем.

!!! note "Где лежат исходники тестов"
    Все тесты — в [`tests/`](../../tests/), эталонные данные сквозного теста — в
    [`tests/data/`](../../tests/data/). Путь к данным передаётся в тест
    макросом `SPH_TEST_DATA_DIR` из CMake и не зависит от рабочего каталога.

### Тесты нитей, импорта `.k`, закреплений и Джонсона–Кука

Эти тесты зарегистрированы в [`CMakeLists.txt`](../../CMakeLists.txt) наравне с перечисленными
выше (метка `3d` у большинства; без метки — `fiber_plasticity_test`, `fiber_dim_guard_test`, `k_reader_test`, `k_solid_{1d,2d}_test`, `jc_strength_*`, `jc_eos_test`, `jc_urs_consistency_test`, `jc_integration_*`, `jh2_test`, `jh2_element_test`; `art_stress_*` — `jc`; `e2e_jc_*` — `e2e;jc`; `e2e_clamp_*` — `e2e;3d`; `fiber_plasticity_guo_smoke`, `e2e_fiber_frames_selfcheck`, `e2e_guo2019_smoke_*`, `energy_log_test`, `fabric_checks_smoke_*` — `3d;e2e`). Список сверен с `add_test` на дату правки; полный
перечень всегда выдаёт `ctest --test-dir build -N`.

| Тест (имя в ctest) | Исполняемый файл | Что проверяет |
|---|---|---|
| `fiber_model_test` | `test_fiber_model` | Модель данных нитей: разбор `.fib`, топология, `l0`/`l_break`, масса бусин, `fiber_output` |
| `fiber_bond_test`, `fiber_break_test` | `test_fiber_bond`, `test_fiber_break` | Сила и предел шага связи Fiber-Bond; разрыв ровно при `l > l_break` |
| `fiber_contact_projection_test` | `test_fiber_contact_projection` | Нормальный контакт как проекция скорости (законы 1–3 при `*_contact_predictor = 1`): формула одного контакта, тёплый старт, слоты трения (U1–U6); статика под нагрузкой, лобовой удар, связанные контакты, энергия и импульс стадии (Т1–Т7) |
| `fiber_stitch_test` | `test_fiber_stitch` | Директива `BOND` и материал связи `mat`: разбор и ошибки, сила стежка $E_1A_1(l-l_0)/l_0$ и ноль при сжатии, разрыв стежка без изменения `prev`/`next`/`object_id` |
| `fiber_contact_test`, `fiber_contact_types_test`, `fiber_friction_test` | `test_fiber_contact`, `test_fiber_contact_types`, `test_fiber_friction` | Контакт нить–нить (пары сегментов), законы нормальной силы (Герц, пинбол, Свигл, импульсный), трение проекцией скоростей |
| `fiber_solid_test`, `fiber_solid_types_test` | `test_fiber_solid`, `test_fiber_solid_types` | Контакт нить–Solid и его законы нормальной силы, включая закон 4 «Герц + демпфер» (формула $F_H + c\cdot\max(\dot p,0)$, разбор `fiber_*_contact_type = 4`; `fiber_solid_types_test`) |
| `fiber_plasticity_test`, `fiber_plasticity_integration_test`, `fiber_plasticity_validation_test`, `fiber_input_inelastic_test`, `fiber_damage_test` | `test_fiber_plasticity`, `..._integration`, `..._validation`, `test_fiber_input_inelastic`, `test_fiber_damage` | Пластичность, вязкость и повреждённость связи, директивы `PLASTIC`/`RATE`/`VISCO`/`FAILURE`/`DAMAGE` в `.fib` |
| `fiber_plasticity_guo_smoke` | `test_fiber_plasticity_guo` | Сцена Guo 2019 с упругой и упруго-пластичной нитью (метки `3d;e2e`) |
| `fiber_dim_guard_test` | `test_fiber_dim_guard` | Нити при `SPH_DIM != 3` дают исключение (цель собрана без `SPH_DIM`) |
| `fiber_time_integration_equivalence_test`, `ode_solver_test` | `test_fiber_time_integration_equivalence`, `test_ode_solver` | Эквивалентность Solid-цикла двух перегрузок `time_integration`; Euler и предиктор–корректор |
| `sph_body_isolation_test` | `test_sph_body_isolation` | Опция `sph_same_body_only`: SPH-пары Solid только внутри одного тела |
| `clamp_test_t1/t2/t4` | `test_clamp` | Кинематическое закрепление: список, проекция, реакция, `CLAMP_BOX`; при 1, 2 и 4 потоках |
| `e2e_clamp_<ES1..ES3,EF1..EF6,EK1>_...smoke` | `test_e2e_clamp --case=...` | Закрепления против аналитики (реакция опоры, отражение волны, частота); дополнительно `e2e_clamp_EF1b_pc_smoke` (разрыв у кромки зажима) и `e2e_clamp_{ES1,EF1}_pc_smoke_t4` (4 потока). Полные прогоны `e2e_clamp_<ES1,ES2,ES3,EF1,EF2,EF3>_pc_full` — метка `slow` |
| `e2e_fiber_contact_test`, `e2e_fiber_solid_test`, `e2e_fiber_frames_selfcheck` | `test_e2e_fiber_contact`, `test_e2e_fiber_solid`, `test_e2e_fiber_frames` | Две нити отскакивают без проникания; нить не проникает в пластину; читатель кадров нитей |
| `e2e_guo2019_smoke_{A,1T,1T_A,A_1T,2T,4T}` | `test_e2e_guo2019 --config=...` | Сцены Guo 2019 (шар, ткань, Al-пластина) на уровне smoke, 2 потока (метки `3d;e2e`) |
| `energy_log_test` | `test_energy_log` | Лог механической энергии (формат CSV, замыкание баланса); цель всегда собрана с `SPH_ENERGY_LOG` |
| `fabric_checks_smoke_strip`, `fabric_checks_smoke_capstan` | `fabric_checks` (`tools/fabric_checks/`) | Короткие прогоны проверочных сцен ткани: растяжение полоски, капстан |
| `plain_weave_test`, `plain_weave_gen_cli_test` | `test_plain_weave`, `plain_weave_gen` | Геометрия генератора плетения ([руководство](02-plain-weave-generator.md)) и его CLI |
| `plain_weave_quilt_test` | `test_plain_weave_quilt` | Стёганое одеяло: полотно из патчей побитово равно сплошному, метки патчей, линии шва, стежки `grid` (число, вертикальность, соседние слои) и `seam` ([ткани](../04-explanation/models/07-fabrics.md)) |
| `plain_weave_random_cut_test` | `test_plain_weave_random_cut` | Мозаика (случайная раскройка): детерминизм по `seed`, площади ячеек, зазор между фрагментами, хорды, стежки в слое и сквозные, пакеты с `layers_per_cut`, зажим, `.cut.json` |
| `k_reader_test`, `k_parts_test` | `test_k_reader`, `test_k_parts` | Чтение LS-DYNA `.k`; таблица деталей |
| `k_solid_{1d,2d,3d}_test` | `test_k_solid_{1d,2d,3d}` (один исходник `test_k_solid.cpp`) | Сборка Solid из `.k` для `SPH_DIM` = 1, 2, 3 |
| `k_legacy_equiv_test`, `k_fiber_test`, `k_clamp_test`, `k_solid_elements_test`, `k_e2e_test` | одноимённые `test_k_*` | Эквивалентность прежнему разбору; нити из `*ELEMENT_BEAM`; закрепления из `*BOUNDARY_SPC_*`; `*ELEMENT_SOLID` в частицы; сквозной импорт `.k` через `input()` и `time_integration` |
| `jc_strength_{1,2,3}d_test`, `jc_eos_test`, `jc_urs_consistency_test` | `test_jc_strength_Nd`, `test_jc_eos`, `test_jc_urs_consistency` | Модель Джонсона–Кука: прочность и разрушение (U1–U14), УРС Linear/Shock (U15–U20), совпадение `urs()` и `solid_eos()` ([руководство](06-johnson-cook-material.md)) |
| `jc_integration_{1d,3d,3d_elog}_test` | `test_jc_integration_*` | Интеграционные проверки I1–I9 (в том числе нагрев Тейлора–Куинни) через `apply_material_strength`; вариант `3d_elog` собран с `SPH_ENERGY_LOG` |
| `e2e_jc_hel`, `e2e_jc_spall_nodamage`, `e2e_jc_spall` | `test_e2e_jc_plate --case=...` | Плоский удар стали 4340 (`SPH_DIM=1`): упругий предвестник, откол без разрушения и с ним (метки `e2e;jc`) |
| `art_stress_{1,2,3}d_test` | `test_art_stress_{1,2,3}d` | Искусственное напряжение Монагана: чистые функции (A1–A5) и член в конвейере — при `art_stress_eps = 0` результат побитово прежний, сохранение импульса и энергии, нулевой вклад при сжатии (P1–P5). Метка `jc` |
| `jh2_test` | `test_jh2` | Керамика JH-2 (`SPH_DIM=1`, B₄C): приращение пластической деформации, `sigma_i` в плоском случае, разрыхление и повреждение, флаги `jh2_clamp_low_rate` и `jh2_efmin` (T1–T7, [ADR-010](../04-explanation/decisions/adr-010-jh2-fix.md)) |
| `jh2_element_test` | `test_jh2_element` | Одноэлементный приёмочный расчёт JH-2: одномерное сжатие B₄C со скоростью $10^5$ с⁻¹ до $\mu = 0{,}1$ и разгрузка |
| `e2e_jc_taylor_smoke` | `test_e2e_jc_taylor --case=smoke` | Короткий прогон теста Тейлора (`SPH_DIM=2`). Полные `e2e_jc_taylor_e4` (цель `test_e2e_jc_taylor_elog`) и `e2e_jc_taylor_e5` — метка `slow` |

Медленные тесты (`slow`) регистрируются только при `-DSPH_E2E_SLOW_TESTS=ON`; сюда входят и
`e2e_clamp_*_pc_full`, и `e2e_jc_taylor_e4/e5`, и 3D-прогоны `e2e_3d_full_contact*`.
Цели без `add_test`: `guo2019_shot`, `demo_fiber_net`, `demo_fiber_layers` и `fiber_t0_dump` (последняя
с `EXCLUDE_FROM_ALL`).

### Полный прогон BC2

Регистрируется только при явном включении, потому что считает минуты:

```bash
cmake -S . -B build -DSPH_E2E_SLOW_TESTS=ON
cmake --build build -j8
ctest --test-dir build -L slow --output-on-failure
```

Каждый тест — один прогон, по тесту на `contact_force_type`. Все четыре можно
пустить параллельно:

```bash
ctest --test-dir build -L slow -j4 --output-on-failure
```

Прогон ограничен по времени сторожем. Бюджет задаётся при конфигурации:

```bash
cmake -S . -B build -DSPH_E2E_SLOW_TESTS=ON -DSPH_BC2_RUN_TIME_LIMIT=1800
```

По умолчанию 900 с на прогон. Это **не** оценка длительности, а backstop на
случай зависания, поэтому он взят с запасом: на машине разработки полный прогон
идёт около 4,3 шага/с, то есть 1531 шагу нужно примерно 355 с. На более
медленном железе бюджет стоит поднять — иначе сторож оборвёт здоровый расчёт.
Коллапс шага ловится независимо от этого бюджета, вторым механизмом сторожа
(детектор застоя), поэтому запас по времени защиту не ослабляет.

При исчерпании бюджета расчёт обрывается и сообщает, что именно произошло — не
хватило времени при здоровом шаге или расчёт перестал продвигаться.

### Полный прогон BC6

Устроен так же, но кейс тяжелее: 44 638 частиц против 33 944 и `t_total`
4,3 мкс против 3,5 мкс. Своя переменная бюджета, по умолчанию 1800 с:

```bash
cmake -S . -B build -DSPH_E2E_SLOW_TESTS=ON -DSPH_BC6_RUN_TIME_LIMIT=1800
cmake --build build -j8
ctest --test-dir build -R e2e_bc6_full --output-on-failure -j1
```

Здесь уместно `-j1`, а не `-j4`: каждый прогон держит около 300 МБ временных
буферов на шаг, пишет во временный каталог порядка 1 ГБ кадров и запускает
семь потоков OpenMP.

Исследование сходимости (§6 постановки, «по требованию») в ctest не
регистрируется и запускается вручную:

```bash
./build/test_e2e_bc6 --full --contact-type=0 --delta-x=1e-4 --time-limit=7200
```

### Графики для отчётов

```bash
cmake --build build --target bc2_plot   # notes/contact/e2e/bc2_monitor_velocity.svg
cmake --build build --target bc6_plot   # notes/contact/e2e/bc6_monitor_velocity.svg
```

Каждая цель делает четыре полных прогона (по одному на `contact_force_type`)
и собирает из их выгрузок сводный SVG. В обычном прогоне тестов графики не
трогаются: запись в дерево исходников — побочный эффект, которого тест не
должен иметь.

## Трёхмерный расчёт

```bash
ctest --test-dir build -R "kinematics_3d|e2e_3d|_3d$|bc2mini3d" --output-on-failure -j8
```

Метка `3d` теперь шире: её несут и тесты нитей, импорта `.k`, закреплений, чекпоинтов
и сцен Guo 2019, поэтому `ctest -L 3d` запускает почти весь набор. Регулярное выражение
выше выбирает только тесты ветки `dim = 3` силового расчёта Solid.

Сквозные тесты BC2 и BC6 и 2D-регрессия считают только осесимметричную
постановку. Ветка `dim = 3` проходит через конвейер другими путями: шесть
компонент девиатора, три поворота, трёхмерная сетка поиска соседей, контакт
без зеркальных частиц. Нити (notes/fiber-model/fiber-plan.md) живут только в 3D, поэтому эта
ветка проверена отдельно, тремя программами со своими объектниками конвейера
(`SPH_DIM=3`):

* `test_kinematics_3d` — алгебра 3D-ветки на решётке, по образцу
  `test_axisym`. Проверки опираются на симметрию и тензорные тождества, а не
  на точность SPH-градиента;
* `test_e2e_3d` — уменьшенный BC2
  ([`cases/bc2/regression/bc2_mini3d_geometry.hpp`](../../cases/bc2/regression/bc2_mini3d_geometry.hpp)):
  материалы, толщины и скорости BC2, радиус керамики 4 мм, 3648 частиц. Полный
  BC2 в 3D даже на сетке smoke — сотни тысяч частиц;
* `test_regression_3d` — та же программа, что `test_regression`, в сборке
  `dim = 3`, кейс `bc2mini3d`, K = 50.

!!! note "Почему энергия проверяется без XSPH"
    Поправка XSPH `av` прибавляется к **скорости**, а не только к
    перемещению (`time_integration`). Она сглаживает поле скоростей без
    нагрева, и на уменьшенном BC2 за 1 мкс уносит около четверти полной
    энергии, при этом импульс сохраняется точно. Поэтому обычный прогон
    проверяет лишь, что энергия не растёт, а сохранение полной энергии —
    отдельный прогон с `epsilon = 0` (+0,4 % за 1 мкс).

!!! note "Что нашли 3D-тесты"
    Юнит-тест вращения нашёл неверный знак поворотного члена
    `syz·ryz` в `dszz` (`src/int_force.cpp`) и в его копии в JH-2
    (`src/ppochnost.cpp`). При вращении вокруг оси x девиатор переставал быть
    бесследовым, а сквозной 3D-прогон терял симметрию x ↔ y на уровне 6e-3.
    Исправлено; в 2D `ryz = 0`, и 2D-регрессия подтвердила, что осесимметричная
    ветка не изменилась.

## Регрессия

```bash
ctest --test-dir build -L regression --output-on-failure -j8
```

Набор 2D (9 тестов) идёт 10–15 с, 3D (5 тестов) — около минуты. Эталоны —
24 файла в `cases/bc2/regression/snapshots/` (BC2 и уменьшенный BC2 в 3D) и
`cases/bc6/regression/snapshots/`, около 20 МБ. Постановки, которые берёт ctest, —
замороженные копии в `cases/<кейс>/regression/` ([правила](../../cases/README.md#роли-постановок)).

### Чем отличается от e2e

Сквозные тесты сверяют расчёт с **приборной** кривой в грубых воротах:
амплитуда может отличаться до двух раз. Они отвечают на вопрос «похоже ли на
физику» и ловят только кратную поломку. Регрессия отвечает на другой вопрос:
«изменился ли расчёт». Ошибка, которая сдвигает, например, осевой член
импульса на 5 %, e2e проходит, а регрессию — нет.

### Как устроена

* Кейсы BC2 и BC6 на сетке smoke (Δx = 0,8 мм, 2152 и 2930 частиц), по тесту
  на `contact_force_type` 0..3. В сборке `dim = 3` (`test_regression_3d`) —
  уменьшенный BC2 в 3D (3648 частиц) с контрольными точками N = 1 и N = 50.
* На кейс — два независимых прогона с нуля, **N = 1** и **N = 100** шагов.
  Останов делает предел шагов сторожа (`RunGuard::set_step_limit`), а не
  `t_total`: так сравниваются состояния после одинакового числа шагов, даже
  если `dt` чуть разошёлся. N = 1 проверяет одно вычисление правой части и
  один шаг интегратора, N = 100 — накопление по времени.
  Контакт проверяется только на N = 100: в начальной сборке тела не
  перекрываются, и на первом шаге все четыре типа контакта дают побитово
  одинаковый результат. По той же причине N = 1 не видит членов, зависящих
  от девиатора (он в начальном состоянии нулевой).
* Снимок берётся прямо из памяти, а не из кадров: `sph_output_*` и `log.txt`
  пишутся с 6 значащими цифрами.
* Сравниваются 15 полей (`x`, `v`, `rho`, `p`, `u`, `c`, `hsml`, `e_pl`, `D`
  и девиатор) по нормам Linf и L2 относительно масштаба поля. `ntotal`,
  число шагов и `itype` каждой частицы должны совпасть точно. Масса и
  нулевые при dim = 2 компоненты `syz`, `sxz` проверяются отдельно.

### Как читать отчёт

Отчёт печатается и при успехе. По каждому полю видны Linf и L2 против
допуска и **худшая частица**: номер, `itype`, координаты, значения в эталоне
и в расчёте. По ней сразу видно, где разошлось: на контактной границе, на
оси, в одном материале. Если снимки несравнимы (разошлись `itype` или
`ntotal`), отчёт называет первую такую частицу.

### Допуски

Допуск — это разброс результата при разном числе потоков OpenMP, умноженный
на 100, но не меньше 1e-11. Поиск соседей раздаёт номера пар через atomic,
поэтому при нескольких потоках порядок суммирования меняется, и результат
плавает на уровне округления. Ровно такой же шум даёт любая правка, которая
только переставляет слагаемые (перенос Solid на owner-computes, этап Ф4
плана нитей). Поэтому допуск её пропускает, а изменение формулы — нет.

Числа калибровки и дата лежат в комментарии у таблицы допусков в
[`tests/test_regression.cpp`](../../tests/test_regression.cpp).

### Когда и как перегенерировать эталон

Только при **осознанном** изменении физики или численной схемы, и с причиной
в сообщении коммита. Если регрессия упала после правки, которая по замыслу
ничего не меняет, эталон не трогают: ищут ошибку.

```bash
cmake --build build --target regression_golden
```

Цель ставит `OMP_NUM_THREADS=1`, пересчитывает все 24 файла в
`cases/<кейс>/regression/snapshots/` (2D и 3D) и вписывает в их шапку commit и дату. Ни один
зарегистрированный тест в дерево исходников не пишет.

Проверка, что регрессия вообще что-то ловит: прогон с другим конфигом должен
провалиться.

```bash
cp cases/bc2/regression/initial_conditions.txt /tmp/bc2_alpha.txt
sed -i 's/^alpha = .*/alpha = 1.001/' /tmp/bc2_alpha.txt
./build/test_regression --case=bc2 --contact-type=0 --config=/tmp/bc2_alpha.txt
```

## Бенчмарк

Отдельного бенчмарка в сборке нет; замеры ускорения лежат в `notes/speedup/*.md`.
