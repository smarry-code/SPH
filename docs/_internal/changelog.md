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
- Записи архитектурных решений ADR-001…ADR-006 и шаблон ADR.
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
  Тейлора–Куинни (`jc_heating`, `jc_chi`); поля частицы `T`, `edot`, `r0`; чекпоинт
  формата v3; в 2D осесимметрии вес энергии `2π·r0·m` в логе энергии; тесты
  U1–U20, I1–I9, E1–E5; ссылки на [ADR-009](../04-explanation/decisions/adr-009-johnson-cook.md)
  и [how-to 06](../02-how-to/06-johnson-cook-material.md).

### Исправлено

- В 2D осесимметрии добавлен член −(p/ρ)·v_r/r в уравнение энергии (работа давления на
  окружной деформации); регрессионные снимки BC2/BC6 (N = 100) обновлены после него.
- Совместимость чекпоинта с Windows: запись через `_commit`, ротация файла `.prev`.
