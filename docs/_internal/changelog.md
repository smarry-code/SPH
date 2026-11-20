# Журнал изменений документации

<!-- Этот документ живёт внутри docs/; ссылки на репозиторий ведут за пределы
     docs_dir и проверяются линкером lychee в CI (см. mkdocs.yml: validation). -->

Хронология изменений системы документации. Записи ведутся в формате
*Keep a Changelog*.

## [Неопубликовано]

### Добавлено

- Каркас каталога `docs/` по модели Diátaxis.
- Главная страница сайта документации ([`index.md`](../index.md)).
- Корневой `README.md` (корень репозитория) — точка входа для отображения на GitHub.
- Глоссарий терминов (RU/EN).
- Конфигурация MkDocs Material (`mkdocs.yml`) и `requirements-docs.txt`.
- Конфигурация Doxygen (`Doxyfile`) для автогенерации API.
- Быстрый старт — пошаговое руководство сборки и запуска.
- Спецификации форматов файлов: `Initial_Conditions.txt`, `*.k`, выходные `.dat`/`.txt`.
- Раздел архитектуры: обзор расчётного конвейера, SoA, поиск соседей, параллелизм.
- Записи архитектурных решений ADR-001…ADR-005 и шаблон ADR (ADR-006 удалён вместе с раскраской графа, нумерация продолжается с ADR-007).
- Контактные силы «частица-частица» (типы 1-3: пинбол Белицко, Свигл,
  импульсный) — `include/contact_force.hpp`, оркестратор в `ext_force.cpp`,
  ограничитель шага и диагностика в `WorkBuffers`, ключи конфигурации
  `contact_force_type`/`contact_penalty_scale`, тесты `tests/test_contact.cpp`
  (`contact_test`), [ADR-007](../04-explanation/decisions/adr-007-contact-force-types.md)
  и [описание архитектуры](../04-explanation/architecture/contact-particle-particle.md).
  По умолчанию контакт выключен (`contact_force_type = 0`), поведение
  существующих кейсов не изменилось.
- Рестарт расчёта из чекпоинта: модуль `include/checkpoint.hpp` / `src/checkpoint.cpp`,
  хуки в `time_integration`, ключи `--checkpoint`, `--checkpoint-every`, `--restart`,
  `--force` в `guo2019_shot`, поля `checkpoint_file`/`checkpoint_every_steps`/`restart_file`/
  `restart_force` в `SimulationConfig`, тесты `checkpoint_test`, `checkpoint_energy_test`,
  `checkpoint_cli_test`. Руководство с правилами и рекомендациями —
  [`05-restart-checkpoints.md`](../02-how-to/05-restart-checkpoints.md).
- Модель Джонсона-Кука (itype 12–15: STEEL 4340, STEEL 1006, AL 7039, AL5083H116) —
  прочность, температура, разрушение, УРС Linear и Shock, `p_min`; опция нагрева
  Тейлора–Куинни (`jc_heating`, `jc_chi`); поля частицы `T`, `edot`, `r0` (в чекпоинте); в 2D осесимметрии вес энергии `2π·r0·m` в логе энергии; тесты
  U1–U20, I1–I9, E1–E5; ссылки на [ADR-009](../04-explanation/decisions/adr-009-johnson-cook.md)
  и [how-to 06](../02-how-to/06-johnson-cook-material.md).

- Ключ `xsph_mode = velocity | position` (по умолчанию `velocity`, поведение прежнее): `position` —
  классический XSPH Монагана, диагностический режим; без стабилизатора расчёты Тейлора и BC2/BC6
  обрываются по $\rho<0$ (см. `notes/johnson-cook/05-xsph-axisym-energy.md`).
- Лог энергии (`SPH_ENERGY_LOG`): диагностические каналы мощности `P_pair_*`, `P_hoop_v`, `P_av_*`,
  `P_dev_u`, `P_pdv_u`, `P_szz_u` и `closure` (только чтение состояния).
- Тесты `e2e_clamp_*_pc_full`: бюджет `--time-limit` 12000 с, TIMEOUT 14400 с (ES2 — 10 568 с).
- Раздел [«Математические модели»](../04-explanation/models/index.md) (`docs/04-explanation/models/`):
  ядро SPH, определяющие соотношения (УРС, Мизес, JH-2, Джонсон–Кук), контакт SPH-тел и зажимы,
  механика нитей, контакт нитей (законы 0–4, трение), интегрирование по времени и журнал энергии;
  сводка замечаний к коду, найденных при сверке формул.
- Ключи `art_stress_eps`, `art_stress_n` — искусственное напряжение Монагана (Monaghan 2000; Gray,
  Monaghan, Swift 2001), в осесимметрии с окружной компонентой; по умолчанию выключено, поведение
  прежнее. Неустойчивость у оси в тесте Тейлора не устраняет (`notes/johnson-cook/06-artificial-stress.md`).
  Тест `art_stress_{1,2,3}d_test`; ключи e2e Тейлора `--art-eps`, `--save-t`, `--no-damage`, `--pmin`.
- Глава [«Ткани: одеяло и мозаика»](../04-explanation/models/07-fabrics.md) с пятью схемами
  (`docs/assets/img/fabrics/`): полотняное переплетение 1/1, пакет kT, стёганое одеяло (патчи, линии шва,
  проколы, стежки `grid` и `seam`, формулы числа стежков), мозаичная ткань (раскройка Вороного, стежки в слое
  и сквозные, варианты A/B), стежок в ядре (`FiberBond::mat`).
- Справочник утилиты [`guo2019_shot`](../03-reference/guo2019-shot.md): все ключи (`<k>T`, шар, габариты ткани,
  проекция контакта, стёжка, мозаика, чекпоинт), имя каталога, `shot.txt`, `monitor.csv`, `scene_id.txt`.
- Формат [`.cut.json`](../03-reference/file-formats/cut-json.md) (`random_cut_v1`); директива `BOND` в генераторах.
- Нормальный контакт нитей как проекция скорости (законы 1–3, ключи `fiber_*_contact_predictor`,
  `fiber_contact_projection_iters`, `fiber_contact_projection_warm`; тёплый старт, метки 6, 7, 14) —
  [глава 5](../04-explanation/models/05-fiber-contact.md), справочник конфигурации.
- [ADR-010](../04-explanation/decisions/adr-010-jh2-fix.md) «Исправление JH-2»; раздел «Искусственное
  напряжение Монагана» в [главе 1](../04-explanation/models/01-sph-core.md).
- Тесты в how-to: `fiber_contact_projection_test`, `fiber_stitch_test`, `plain_weave_quilt_test`,
  `plain_weave_random_cut_test`, `stitch_kt_cli_test`, `mosaic_kt_cli_test`, `art_stress_{1,2,3}d_test`,
  `jh2_test`, `jh2_element_test`.

### Изменено

- Рабочие заметки (`Gia/`, `contact/`, `fiber/`, `johnson-cook/`, `jh2/`, `speedup/`, `part_0.md`,
  `Report_SPH_contact.md`) перенесены в [`notes/`](../../notes/README.md) по темам (2026-11-20), имена — kebab-case;
  ссылки в коде, `docs/`, `experiments/` и самих заметках обновлены; графики e2e — `notes/contact/e2e/`.
- Документация после слияния веток `fabric_penetration`, `fiber_fabric`, `fabric_mosaic`,
  `fiber_johnson_cook`, `fiber_johnson_holmquist` (2026-11-18):
  - чекпоинт: одна версия формата (6), старые не читаются; блоки `PROJ` (тёплый старт) и `USER`
    (состояние монитора выстрела восстанавливается после `--restart`), материал связи в `TOPO`, `scene_id.txt`;
  - параллелизм: новые примитивы `CpuExecutor`, параллельные циклы ODE, `SpatialIndex::build`,
    `build_solid_extended`, `Diagnostics::validate` (результат побитово прежний); ускорение контакта нитей
    (`cull_apart`, активные слоты трения);
  - `monitor.csv`: столбцы реакций зажимов, проекции и стежков;
  - руководство по генератору ткани: стёжка kT, `--stitch-patch`, мозаика в `guo2019_shot`;
  - глава 4 (связи вне цепочки, `mat`), глоссарий (раздел «Ткани», проекция контакта, искусственное
    напряжение), номера строк в ссылках на код;
  - комментарии и справка `guo2019_shot` (без «формата v3»), комментарий `artificial_stress.hpp`
    (окружная компонента в осесимметрии).

- STEEL 4340: `p_min` = −5,5 ГПа (было подобранное −4 ГПа) — откольная прочность 4340 RC30 по
  Bless 1981 (AFWAL-TR-81-4040); пометка `TODO(verify)` снята. Константы 4340 сверены с библиотекой
  Autodyn (`notes/sources/steel-4340-autodyn-constants.md`).
- Аудит документации на соответствие коду `fiber_new` (2026-11-14):
  - конвейер, поиск соседей и параллелизм переписаны под `rhs::SpatialIndex`, `CpuExecutor`
    и reduce по шардам; главы о `link_list` и `FlatGrid` помечены как исторические
    (код удалён в a4ce0b9);
  - справочник конфигурации дополнен ключами `ode_solver`, `dt_scale`, `sph_same_body_only`,
    контактов нитей, `push_*`, `jc_*`, `energy_log_*`; помечены ключи, которые читаются, но
    не влияют на расчёт (`nnps`, `coloring_niac_threshold`), и поля, не читаемые из файла
    (`num_threads`, `xi`, `epsilon_h`); `sle = 2` и `SPH_PROFILE` описаны как нерабочие;
  - список тестов в how-to приведён к `add_test` в `CMakeLists.txt`; удалены несуществующие
    `flat_grid_test`, `pair_graph_test`, `bench_color`;
  - ADR-001…009: исправлены ссылки на код, добавлены разделы «Актуальное состояние (2026-11)»;
  - глоссарий дополнен терминами RHS, контактов нитей и моделей; на главной — раздел «Возможности»;
  - обновлены номера строк в ссылках на код во всех главах.
- `mkdocs.yml`: юникодные якоря заголовков (`pymdownx.slugs`), раздел «Математические модели» в nav.

- `README.md` переписан (возможности, сборка и запуск, тесты, ссылки на сайт документации,
  структура репозитория). Рабочие заметки (`contact/`, `fiber/`, `Gia/`, `johnson-cook/`,
  `speedup/`, `Report_SPH_contact.md`; с 2026-11-20 — в `notes/`): планы и ТЗ получили врезку «Статус (2026-11-14)»
  по факту реализации, отчёты с устаревшим описанием кода — «Примечание»; индекс
  `contact/README.md` со статусами; исправлены битые относительные ссылки.

### Удалено

- Каталог `SPH_Lee/` (исходная Fortran-реализация, справочно) и ссылки на него в `README.md`, `Doxyfile`,
  ADR-001 (2026-11-20); код доступен в истории git.
- Главы «Связный список (link_list)», «Сетка FlatGrid» и ADR-006 (раскраска графа пар):
  описанный в них код удалён в a4ce0b9 (этап Ф4-B2), тексты доступны в истории git.
- Исторические врезки о старом конвейере в главах архитектуры, раздел о `bench_color`,
  строка «Graph coloring» глоссария.
- Рабочие заметки: `Gia/E2E_EXPERIMENT_CONTEXT.md` (принципы e2e-тестов и чек-лист автора ТЗ
  перенесены приложениями А и Б в `Gia/TZ_GUO2019_EXPERIMENTS.md`, ныне `notes/guo2019/guo2019-experiments-tz.md`), `Gia/E2E_CONTEXT_CODE_MAP.md`,
  `fiber/DEMO_VISUALIZATION_TASK.md`, `contact/03-review-brief.md`, `contact/10-debug-handoff.md`;
  ссылки на них заменены текстом с пометкой «(удалён)».

### Исправлено

- Хук `mkdocs_hooks/source_links.py` публикует исходники рекурсивно (`glob` с `recursive=True`): файлы из подкаталогов
  `include/` (`rhs`, `interactions`, `fiber`, `ode`, `strength`, `k_input`) теперь попадают в сайт, ссылки на них не дают 404.
- JH-2 (ветка `fiber_johnson_holmquist`, ADR-010): приращение пластической деформации радиальным возвратом,
  разрыхление по потере упругой энергии, `sigma_i` в плоском случае; результаты керамики изменились.
  **Эталоны регрессии BC2/BC6 (2D) после слияния нужно переснять**: взяты из ветки JH-2 и не учитывают
  член $-(p/\rho)\,v_r/r$ осесимметрии.

- Реестр материалов: УРС `itype` 3–9 (Ми–Грюнайзен — 3, 4; УА Тэта — 5–9) были указаны неверно;
  константы B₄C приведены к `materials.hpp` (2517 кг/м³, 199,4 ГПа, K₁ = 230,5 ГПа).
- Формат кадра `.dat`: добавлены колонки `T`, `edot`; описание `log.txt` приведено к коду.
- `.fib`: описана директива `DAMAGE` и колонка `D` в выводе связей.
- Ссылки на код из `02-how-to/` и `03-reference/cli.md` уходили выше корня репозитория (`../../../`).
- Якорные ссылки с кириллицей не работали (стандартный slugify выбрасывал буквы).

- В 2D осесимметрии добавлен член −(p/ρ)·v_r/r в уравнение энергии (работа давления на
  окружной деформации); регрессионные снимки BC2/BC6 (N = 100) обновлены после него (затем заменены эталонами ветки JH-2, см. выше).
- Совместимость чекпоинта с Windows: запись через `_commit`, ротация файла `.prev`.
