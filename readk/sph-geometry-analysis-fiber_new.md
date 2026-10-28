# Работа с геометрией в smarry-code/SPH, ветка fiber_new

Ветка `origin/fiber_new`, коммит e930bc3. Сравнивается с master (отчёт [sph-geometry-analysis.md](sph-geometry-analysis.md)). Ссылки `файл:строка` относятся к fiber_new.

## Главное

- **Появилась вторая геометрия: нити (бусины на ломаных).** Это отдельная модель `fiber::Model` с собственным файловым форматом `.fib`. Солид-часть (`ParticleSystem`) осталась прежней, добавлено поле `body[]`.
- **Файловый `.k` не подключён и здесь:** `parse_k_file` теперь объявлен в `include/input.hpp:21`, но вызывается только из теста; в `input()` он по-прежнему в комментарии (`src/input.cpp:857`).
- **Рабочая геометрия** в основном строится не в `src/`, а в генераторах `tests/*.hpp` и `include/fiber/plain_weave.hpp` (e2e-сцены BC2, BC6, Guo 2019). Production-`main` (`src/main.cpp`) по-прежнему вызывает только `input()` → `shock_tube`/`shear_cavity`; нити через `main` не загружаются.

## 1. Загрузка из внешних файлов

### 1.1 `.k` — `parse_k_file`, `src/input.cpp:121-260`
Поведение то же, что в master (блоки `*NODE`, `*ELEMENT_SPH`, free-format через запятую, три координаты, слияние по NID). Отличия:
- добавлено `sys.body[idx] = e.pid` (`input.cpp:246`): номер тела = PID;
- в `trim_spaces` исправлено UB `std::isspace` для кириллицы (`input.cpp:64-75`);
- функция объявлена в `include/input.hpp:21` (для теста `tests/test_fiber_model.cpp:408-440`, `test_body_from_k_file`).
Прежние замечания остаются: запись `x[1]`, `x[2]` при `dim=1` (UB), нет фиксированных колонок, жёстко зашитое условие мониторной частицы, вызов в `input()` закомментирован.

### 1.2 `.fib` — файл нитей, `fiber::load_fiber_file`, `src/fiber_input.cpp:102-278` (объявление `include/fiber/fiber_input.hpp`)
**Новый формат.** Текстовый, поля через пробелы, комментарии `#`, формат описан в шапке `src/fiber_input.cpp:1-60`. Директивы:
- `MATERIAL <id> <E> <nu> <d> <rho> <eps_break> [mu_ff [mu_fs]]` — материал нити;
- `FIBER <material_id> <point_count>` + ровно `point_count` строк `x y z` — ломаная (≥ 2 точек), вершины в недеформированной конфигурации;
- `CLAMP_BOX xmin ymin zmin xmax ymax zmax [маска xyz]` — закрепление бусин в боксе (необязательно).
Любая ошибка формата бросает `std::runtime_error`; работает только при `SPH_DIM=3` (`check_fiber_requires_dim3`, `include/fiber/fiber_system.hpp:42`).
Что вычисляется при чтении (`fiber_input.cpp:330-400`): бусина на каждую вершину; `l0` связи — расстояние между точками, `l_break = (1+eps_break)·l0`; масса бусины `ρ·A·(l_prev+l_next)/2`, `A = πd²/4`; `object_id` = порядковый номер блока FIBER; `prev/next` по цепочке; `CLAMP_BOX` → `model.clamps.fiber_mask` через `interactions::mark_clamp_box`.
Вызывается только из тестов/сцен: `tests/guo2019_scene.hpp:422`, `tests/test_fiber_model.cpp`, `tests/test_clamp.cpp:554`. Ключ конфигурации `fiber_file` разбирается (`src/input.cpp:421-427`), но в `cfg.fiber_file` его больше никто не читает.

### 1.3 `particles.dat` — как в master (`src/input.cpp:771`, `config_input=false`); после чтения `body[i] = itype[i]` (`input.cpp:843`).

### 1.4 Конфигурации
`Initial_Conditions.txt` и `tests/data/*_initial_conditions.txt` разбирает `parse_sph_config`; добавлены ключи `contact_force_type`, `contact_penalty_scale`, `contact_impulse_gamma`, `ode_solver`, `sph_same_body_only`, `fiber_file` (`input.cpp:348-427`). Геометрии в них нет.

### 1.5 Эталонные данные
`tests/data/bc2_reference.txt`, `bc6_reference.txt`, `regression/*.txt` — кривые и снимки для сравнения, не входные данные расчёта. Файл `1.k` в корне (тестовый набор SPH-элементов) на master и fiber_new тот же.

## 2. Генерация геометрии силами кода

### 2.1 Во встроенном пути `src/` (как в master + `body`)
- `shock_tube`, `src/input.cpp:432`; `body[i] = itype[i]` в конце (`input.cpp:551-558`);
- `shear_cavity`, `src/input.cpp:687`; `body[i] = itype[i]` (`input.cpp:727`);
- диспетчер `input()` — `src/input.cpp:732`, выбор по `shocktube`/`shearcavity` (`input.cpp:926-933`).

### 2.2 Генератор тканого полотна (нити) — `include/fiber/plain_weave.hpp`
Header-only, не зависит от `dim`. `build_plain_weave(PlainWeaveParams)` (`plain_weave.hpp:417`) строит полотняное переплетение: основа и уток как ломаные с профилем `cos` или `arcs`, бусины в каждом пересечении, между ними `seg` хорд равной длины, выступы за окно, боксы закрепления (`clamp_boxes`). Пресеты — `apply_plain_weave_preset` (`:106`). Запись в `.fib` — `write_plain_weave_fib` (`:738`), раскладка — `write_plain_weave_layout` (`:787`), сводка — `plain_weave_summary` (`:690`).
Консольная утилита: `tools/plain_weave/plain_weave_gen.cpp` (флаги `--out`, `--preset`, `--d`, `--pitch`, `--amp`, `--nx/--ny`, `--seg`, `--clamp` и др.; сборка — `CMakeLists.txt:466`). Документация: `docs/02-how-to/02-plain-weave-generator.md`.

### 2.3 Генераторы тел для тестовых сцен (Solid, `tests/`)
Каждый заполняет `ParticleSystem` напрямую, ставит `body[k] = itype`:
- `tests/bc2_geometry.hpp:134` `build_bc2_geometry` — осесимметричная 2D-сетка, четыре слоя `bc2_layers()` (`:86`: ударник itype 8, керамика 10 ×2, окно 5; радиусы 0.04375/0.0346/0.0254 м, Z1…Z4), число ячеек `cell_count` (`:79`), `axis[k]` на оси, `mass = ρ·dr·dz`;
- `tests/bc2_mini3d_geometry.hpp:71` `build_geometry` — 3D-вариант BC2: круглые сечения на решётке `(ix+½)dx`, масштаб радиуса `kRadiusScale`, симметрийные индексы `mirror_x/mirror_y/swap_xy`;
- `tests/bc6_geometry.hpp:130` `build_geometry` — BC6: четыре тела `kBodies` (`bc6_case.hpp:66`, пена/PMMA/B4C/PMMA), 2D-осесимметрия;
- `tests/guo2019_scene.hpp:205` `build_guo_scene` — удар шара по мишени: шар решёткой с отбором по радиусу (`for_each_ball_node`, `:176`), пластина `nx·ny·nz` (`:362-392`), ткань из `build_plain_weave` (`:272`) через запись и повторное чтение `.fib` (`:411-422`);
- `tests/clamp_case.hpp`, `clamp_e2e_*.hpp`, `clamp_integration_*.hpp` — малые сцены для проверки закрепления; `tools/guo2019/guo2019_shot.cpp` — серия выстрелов Guo 2019.
Регрессионный набор подбирает сцену по `Case` (`tests/regression_cases.hpp:77`).

## 3. Внутренний формат геометрии (цель конвертера)

Сцена теперь — `fiber::Model` (`include/fiber/model.hpp:19`):

| Часть | Структура | Поля |
|---|---|---|
| Solid | `ParticleSystem` (`include/particle_system.hpp`) | как в master + `body[]` (`:63`, номер тела; из `.k` это PID, у генераторов = itype) |
| Fiber | `FiberState` (`include/fiber/fiber_system.hpp:58`) | `x[dim]`, `vx[dim]`, `mass`, `object_id`, `material_id` |
| Топология нитей | `FiberTopology` (`:101`) | `prev[]`, `next[]`, `bonds` (`FiberBond{i, j, l0, l_break}`), `version` |
| Материалы нитей | `FiberMaterial` (`:157`) | `E, nu, d, rho, eps_break, mu_ff, mu_fs` |
| Закрепление | `model.clamps` (`include/interactions/clamp.hpp`) | маска по бусинам и частицам |

Для конвертера `.k` → внутренний формат из `.k` выходят: Solid-частицы (`x`, `mass`, `itype`, `body`); нити в LS-DYNA (`*ELEMENT_BEAM`/`*ELEMENT_DISCRETE`) в парсере не поддерживаются, их надо либо переводить в `.fib` (ломаные + `MATERIAL`), либо расширять формат. Закрепления (`*BOUNDARY_SPC_*`) соответствуют `CLAMP_BOX` (бокс + маска xyz).

## 4. Замечания

1. `fiber_file` из конфигурации никуда не передаётся: в `main` нити не загружаются, работает только через тесты/сцены.
2. `parse_k_file` по-прежнему не вызывается в `input()`.
3. `body[]` из `.k` — это PID, а у генераторов — itype; смысл поля разный в разных путях (комментарии в `particle_system.hpp:55-63`, `input.cpp:242-246`).
4. `.fib` требует `dim==3`, а `.k`-парсер пишет три координаты безусловно.
5. Единицы — СИ, значения материалов нити задаются прямо в `.fib`.
