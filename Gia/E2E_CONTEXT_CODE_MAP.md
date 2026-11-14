# Карта «E2E_EXPERIMENT_CONTEXT.md → код»

> **Примечание (2026-11-14):** карта — снимок на `2b67fd6` (2026-09-30), сверка с тех пор не повторялась: после него в код вошло более 140 изменений (контакт нитей по парам сегментов, трение, закрепление, пластичность и повреждение нитей, лог энергии, Johnson-Cook, k-файлы, удаление старого пути RHS в `a4ce0b9`). Номера строк `file:line` ниже, прежде всего в `src/time_integration.cpp`, `fiber_contact.hpp`, `fiber_solid_contact.hpp`, `CMakeLists.txt`, недействительны; пункты §B.3, §6 («Дефекты») и §7 («Пробелы») описывают состояние до исправлений. Пути и символы, а также утверждения про `solid_contact.hpp`, `body[]`, `RunGuard`, формат кадров — по смыслу в основном сохранились, но по номерам не проверялись. Актуальное описание — [Конвейер](../docs/04-explanation/architecture/pipeline-overview.md), [Контакт нитей](../docs/04-explanation/models/05-fiber-contact.md).

Сверено по `smarry-code/SPH`, ветка `fiber_new`, HEAD `2b67fd6` (совпадает с документом).
Проверка статическая: код прочитан, ничего не собиралось и не запускалось. Числа из прогонов
(скорости, ×2.54, ×18.3, 165/141 связей, времена) не перепроверялись.

Легенда: ✅ подтверждено кодом · ⚠️ верно с оговоркой / неточность · ❌ расходится с кодом · ⏸ не проверялось (прогон или локальная машина).

## Главное

- Документ в целом **точен**: конвейер шага, формулы нитей, форматы файлов, флаги тестов, дефекты §6 и отсутствие их исправлений подтверждаются кодом.
- Настоящих расхождений нашёл четыре, все мелкие (§A ниже). Для ТЗ важнее три **дополнения**, которых в документе нет (§B).

### A. Неточности документа

| § | Утверждение | Что в коде |
|---|---|---|
| 4.4 | «Без `sys.body[k]` не работают фильтры тела (Solid-Solid-Contact, **Fiber-Solid**)» | ❌ `body[]` фильтрует только Solid-Solid (`include/interactions/solid_contact.hpp:105,196`). Fiber-Solid отбирает Solid-частицы по `elastic_by_itype()` (`include/interactions/fiber_solid_contact.hpp:213`), `body` там не встречается. |
| 6.1 п.1 | «в `test_e2e_cross_friction` — 7 primitives, тест допускает до 20» | ⚠️ Проверка `<= 20` живёт в `tests/test_fiber_contact.cpp:462`, а не в сценарии кросс-трения. `test_e2e_cross_friction` — не отдельный файл, а функция внутри `tests/test_e2e_fiber_contact.cpp:391`. Цифра «7» взята из `fiber/friction_chattering_analysis.md:85`. |
| 4.3 | Имя кадра `fiber_output_<ts>.txt` | ⚠️ Номер дополняется нулями до 5 цифр: `fiber_output_00001.txt` (`src/fiber_output.cpp:38`), так же `sph_output_00001.txt` (`src/output.cpp:22`). Читалке нужно это учесть. |
| 3.4 | «14 тестов регрессии: BC2, BC6 × 4, bc2mini3d × 4» | ⚠️ Кейсов 12; до 14 счёт доводят два `regression_selfcheck*` (`CMakeLists.txt:612,639-647`). Эталоны: `step100` для BC2/BC6, `step50` для bc2mini3d. |

### B. Чего в документе нет, а для ТЗ важно

1. **Детектор застоя в демо выключен.** Демо создают `RunGuard guard(900.0, 0.0, 1)` (`demo_fiber_layers.cpp:300`); при `max_projected_steps = 0` ветка Stalled не срабатывает (`include/run_guard.hpp:80`). Новый e2e с нитями должен задавать ненулевой порог явно, как делают BC-тесты.
2. **Время в `log.txt` и `TIME` в кадре отличаются на один `dt`.** Кадр пишется до `time += dt` (`src/time_integration.cpp:238-244`), а строка `log.txt` после (`:247-253`). При совмещении кривой из лога с кадрами это сдвиг на шаг.
3. **Радиус поиска Fiber-Fiber = `d_max + max(l0)`** (`include/interactions/fiber_contact.hpp:179-187`), Fiber-Solid = `r_f_max + r_s_max + l_max` по `l0` (`fiber_solid_contact.hpp:137-151`). Это код за дефектом §6.1 п.4. Там же видно, что primitive — это **пара бусин** из spatial search (`fiber_contact.hpp:222-251`), что и есть корень кратного контакта §6.1 п.1.

## Карта по разделам

### §1. Проект

| Утверждение | Статус | Код |
|---|---|---|
| Реестр `itype` 3–11, JH-2 для 10/11 | ✅ | `include/materials.hpp:12-23,59-68,102-127` |
| `rho_0[i] = rho[i]` при построении | ✅ | `src/input.cpp:487,496,504` |
| `h = 1.3·(m/ρ)^(1/dim)` | ✅ | `src/h_upgrade.cpp:19` |
| `dim` = `constexpr` из `SPH_DIM`, по умолчанию 1 | ✅ | `include/param.hpp:13-16` |
| Каждая цель компилирует свои `src/*.cpp` | ✅ | `CMakeLists.txt:54` (`SPH_PIPELINE_SOURCES`), напр. `:304` |
| Нити при `dim != 3` запрещены | ✅ | `include/fiber/fiber_system.hpp:42`; вызовы `src/fiber_input.cpp:85`, `src/fiber_output.cpp:35`, `src/input.cpp:405` |
| Единый цикл `detail::run_time_integration<Solver>`, две перегрузки-обёртки | ✅ | `src/time_integration.cpp:70-262` (цикл), `:270-304` (обёртки) |
| Порядок шага (predict → single_step → update_h → Bond → FF → FS → dt → advance → strength → разрыв → вывод → log → guard) | ✅ | `src/time_integration.cpp:138, 141, 152, 165-182, 192-205, 209-211, 216, 233, 236, 238-242, 244-253, 259` |
| Кадр на шаге 1 при `time = 0`, общий номер у пары кадров | ✅ | `:111-117, 238-242` |
| `log.txt`: `время v_axis[moni_particle] itype`, только Solid | ✅ | `:247-253` (ось `dim-1`) |
| Файлы пишутся в CWD жёсткими именами | ✅ | `:247`; `src/output.cpp:22,88` (`sph_gp.csv`); `src/fiber_output.cpp:38` |
| `build/` устарел, macOS, python без numpy | ⏸ | Про локальную машину, в репозитории `build/` нет |
| Медленные BC только при `SPH_E2E_SLOW_TESTS=ON`, метка `slow` | ✅ | `CMakeLists.txt:365, 414-435, 524-530` |

### §2. BC2 / BC6

| Утверждение | Статус | Код |
|---|---|---|
| BC6: 4 тела, ρ, `v = 913`, толщины, радиусы | ✅ | `tests/bc6_case.hpp:43-72` |
| BC6: Δx 0.2 / 0.8 мм, smoke `t_total = 1e-6`, full 4.3 мкс | ✅ | `bc6_case.hpp:79-80,124`; `tests/data/bc6_initial_conditions.txt:45` |
| BC6: контрольный монитор `r = 12.5 мм` | ✅ | `bc6_case.hpp:132` |
| BC6: ворота (фронт ±0.05 мкс, плато [170;190], 0.82±0.02, период 0.88±0.15) | ✅ | `bc6_case.hpp:146-152` |
| BC6: эталонные ступени/провалы (182.85 / 171.52 / 141.19 / 115.84; 114.93 и 105.79 м/с) | ✅ | `bc6_case.hpp:179-186` |
| Слой делится на целое число ячеек, узлы по центрам | ✅ | `tests/bc2_geometry.hpp:69-82` |
| Монитор ищется геометрически | ✅ | `bc2_geometry.hpp:124-128,204` (`pick_moni_particle`) |
| BC2: `t_total = 3.5e-6`, smoke `3e-7` | ✅ | `tests/data/bc2_initial_conditions.txt:43`; `tests/test_e2e_bc2.cpp:2062` |
| Конфиг `pa_sph=2, nnps=2, sle=1, skf=1, CFL=0.1, eps=0.25, α=β=1, etq=0.01` | ✅ | `tests/data/bc6_initial_conditions.txt:22-37` (BC2 аналогично) |
| B4C↔B4C не виден: критерий `body` различен, а `body = itype` | ✅ | `solid_contact.hpp:105`; `bc2_geometry.hpp:187` (`body = itype`); `test_e2e_bc2.cpp:25` |
| Текущие результаты (0.659, плато 182.8, 0.703 и др.) | ⏸ | Из `contact/06-e2e-case.md:240`, `contact/06-e2e-bc6.md`, не перепрогонялось |

### §3. Шаблон e2e

| Утверждение | Статус | Код |
|---|---|---|
| Файлы BC6 (`bc6_case.hpp`, `bc6_geometry.hpp`, `test_e2e_bc6.cpp`, данные, `tools/bc6_make_reference.py`, `cmake/bc6_plot.cmake`) | ✅ | все существуют |
| `tests/e2e_*.hpp`: harness / curve / shape / staircase / frames / plot | ✅ | `check/note/warn`, `failures` — `tests/e2e_harness.hpp:38-65`; `ScopedWorkDir` `:78`; `ScopedSilence` `:115` |
| Флаги `--full --contact-type --keep-workdir --emit-curve --plot-from --emit-plot --time-limit` | ✅ | `tests/test_e2e_bc6.cpp:1127-1143` |
| Прогон во временном каталоге | ✅ | `test_e2e_bc6.cpp:206,228` |
| `RunGuard`: WallClock / Stalled / StepLimit, `set_step_limit` | ✅ | `include/run_guard.hpp:13-56,60-100` |
| `e2e_frames.hpp` читает только `sph_output_*` | ✅ | `tests/e2e_frames.hpp:3,110-113` |
| Регрессия: `test_regression`, цель `regression_golden` | ✅ (⚠️ счёт — §A) | `CMakeLists.txt:596-665`, `cmake/regression_golden.cmake` |
| `test_e2e_3d`: `--no-xsph`, законы сохранения | ✅ | `CMakeLists.txt:555-586`; вызов цикла `tests/test_e2e_3d.cpp:306` |
| Продовый цикл с непустыми нитями — только демо | ✅ | Все вызовы `time_integration` в тестах идут с пустыми нитями: `test_e2e_fiber_solid.cpp:130`, `test_fiber_time_integration_equivalence.cpp:182-187` |
| Локальные интеграторы `FiberBondContactIntegrator`, `FiberOnPlateIntegrator` | ✅ | `tests/test_e2e_fiber_contact.cpp` (упоминание `test_e2e_fiber_solid.cpp:27`), `tests/test_e2e_fiber_solid.cpp:199` |

### §4. Нити

| Утверждение | Статус | Код |
|---|---|---|
| Масса бусины `ρA(l_prev + l_next)/2`, у конца одно слагаемое | ✅ | `src/fiber_input.cpp:227,254-257` |
| Bond: только растяжение, `F = EA(l/l0 − 1)` | ✅ | `include/interactions/fiber_bond.hpp:93-96` |
| `dt_bond = 0.8/√(k_eff/m_red)` | ✅ | `fiber_bond.hpp:131-159` |
| Разрыв `l > l_break`, `next/prev = −1`, новый `object_id` у j-стороны, `version++` | ✅ | `include/interactions/postprocessing.hpp:183, 258-278, 295` |
| FF: Герц `(4/3)E*√R*·δ^1.5`, 4 случая, пары только разных `object_id` | ✅ | `fiber_contact.hpp:160-165, 230, 248-251, 439-441` |
| `dt` контакта FF: `k_eff = 1.5F/δ` | ✅ | `fiber_contact.hpp:678-711` |
| Трение `mu_ff`: явный Кулон, смешение `√(μ_a μ_b)` | ✅ | `fiber_contact.hpp:94-104, 485-496, 602` |
| FS: `q = a + λ(b−a)`, `K_fs = (4/3)E*√R*`, `F = K δ^1.5`, вклад в обе стороны | ✅ | `fiber_solid_contact.hpp:316-348, 356-369, 435-436`; сложение `time_integration.cpp:196-201` |
| `r_s^contact = ½(m/max(ρ,ρ₀))^(1/dim)` | ✅ | `solid_contact.hpp:50-57` |
| FS выходит при `ntotal == 0` | ✅ | `fiber_solid_contact.hpp:198` |
| `compute_cfl_timestep` при `ntotal = 0` → `1e6` | ✅ | `time_integration.cpp:23-26` |
| Трения FS нет, `mu_fs` нет | ✅ | `FiberMaterial` без `mu_fs`: `include/fiber/fiber_system.hpp:145-153` |
| Формат файла нитей, 7-е поле `mu_ff`, `object_id` = номер блока | ✅ | `src/fiber_input.cpp:8-43, 108-154, 157-208, 258, 271` |
| Формат кадра нитей | ✅ (⚠️ имя — §A) | `src/fiber_output.cpp:4-11, 38-75` |
| Заголовок кадра Solid в 3D | ✅ | `src/output.cpp:51-75` |
| Схема сборки сцены (куб, анкеры ×50, `maxn = 1.4·ntotal`, `urs`, chdir, `RunGuard(900,0,1)`) | ✅ | `demo_fiber_layers.cpp:130-153, 197, 268-300` |
| `body[]` нужен для Fiber-Solid | ❌ | см. §A |

### §5. Демо

| Утверждение | Статус | Код |
|---|---|---|
| Крест: `E=2e7, ν=0.3, d=1 мм, ρ=1200, ε=0.05`, 21 бусина, `z=±0.5 мм`, куб 8³, dx 1.2 мм, `v0=30`, зазор 4 мм, `t_total=5e-4`, `save_t=8e-6` | ✅ | `experiments/2026-09-26-fiber-net-demo/demo_fiber_net.cpp:162-186, 212-213`; `data/run/net.fib:3-4` |
| 4 слоя: `E=5e7, d=1.5 мм, ε=0.15`, 17 бусин, 5 нитей/слой, шаг 3 мм, зазор слоёв `d`, куб 6³, dx 1.3 мм, зазор 4 мм | ✅ | `experiments/2026-09-27-fiber-layers-demo/demo_fiber_layers.cpp:187-212` |
| pierce `35 м/с, 6e-4, 8e-6`; stop `10 м/с, 1.8e-3, 2.2e-5` | ✅ | `demo_fiber_layers.cpp:219-236`; переопределяются из argv `:241-243` |
| Цели CMake без `add_test` | ✅ | `CMakeLists.txt:735-765` |
| `animate.py`: цвет по слою `object_id // fibers_per_layer`, перцентили 0.5/99.5 | ✅ | `scripts/animate.py:251-261, 295` |
| `analyze_layers.py`: `z_min`, `vz_mean`, NaN | ✅ | `scripts/analyze_layers.py:104-141` |
| Результаты прогонов (×2.54, ×18.3, 165/141, времена) | ⏸ | Требуют прогона |

### §6. Дефекты

| Утверждение | Статус | Код |
|---|---|---|
| Primitive = пара бусин, одно пересечение → несколько primitives | ✅ | `fiber_contact.hpp:222-251`; решение Р8 — `fiber/PLAN.md:20,91` |
| Дребезг явного Кулона | ✅ | Риск описан в самом коде `fiber_contact.hpp:97-102`; диагностика в `test_e2e_fiber_contact.cpp:414-545` |
| FS: пропуск в выпуклом изломе (`t ≥ 1` не последний / `t < 0` не первый) | ✅ | `fiber_solid_contact.hpp:334-342` |
| Радиус поиска по `l0`, без растяжения | ✅ | `fiber_contact.hpp:179-187`, `fiber_solid_contact.hpp:137-151` |
| Инструкции 1–4 лежат в `fiber/`, в коде не внедрены | ✅ | Файлы есть; признаки невнедрения: primitive — пара бусин, нет `mu_fs`, нет стадии трения после `dt` в цикле |

### §7. Пробелы

| Утверждение | Статус | Код |
|---|---|---|
| Монитора нитей нет, хука на шаг нет (`should_stop` получает только `dt, time, t_total`) | ✅ | `time_integration.cpp:253,259`; `run_guard.hpp:60` |
| Нет C++-читалки `fiber_output_*` | ✅ | Общей читалки нет; есть только разовая проверка одного файла в `tests/test_fiber_model.cpp:162-172` |
| Нет энергобаланса для нитей | ✅ | в `tests/test_fiber_*` / `test_e2e_fiber_*` нет |
| `ReduceOptions::deterministic = false`, вызовы по умолчанию | ✅ | `include/rhs/contributions.hpp:126,217`; вызовы `time_integration.cpp:173`, `single_step.cpp:137`, `fiber_bond.hpp:127`, `fiber_contact.hpp:673`, `fiber_solid_contact.hpp:465-466`; побитность проверяет `tests/test_contributions.cpp:182-189` |
| Метки `3d`, `e2e`, `slow`, `regression`; e2e нитей помечены только `3d` | ✅ | `CMakeLists.txt:279,320` (у `e2e_fiber_*` метки `e2e` нет, поэтому `ctest -L e2e` их не запускает) |

### §7.5. Устаревшие документы

| Утверждение | Статус | Где |
|---|---|---|
| `DEMO_VISUALIZATION_TASK.md` §4: «`fiber_output` не вызывается» устарело | ✅ | `fiber/DEMO_VISUALIZATION_TASK.md:79` против `src/time_integration.cpp:241` |
| `IMPLEMENTATION_REPORT.md` §4.5 «ни разу не запускается с непустыми нитями» устарело | ✅ | `fiber/IMPLEMENTATION_REPORT.md:291` (демо запускают) |
| BC2: 985.0 против 868.6 м/с | ✅ (сам факт расхождения) | `IMPLEMENTATION_REPORT.md:17`, `contact/06-e2e-case.md:240`; какое число актуально — ⏸ |
| Зеркальные частицы живут в `solid_sph.hpp` | ✅ | `include/interactions/solid_sph.hpp:22-33` |
