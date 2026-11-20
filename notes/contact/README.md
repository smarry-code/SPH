# Контактные силы «частица–частица»: рабочие материалы внедрения

Каталог — архив планов, ТЗ и отчётов по внедрению контакта (ветка `contact`, сентябрь 2026).
Внедрение завершено и влито; **актуальное описание кода** — в документации `docs/`:

- [Контакт Solid–Solid (модель)](../docs/04-explanation/models/03-solid-contact.md);
- [Контакт нитей (модель)](../docs/04-explanation/models/05-fiber-contact.md);
- [Контакт «частица–частица» (архитектура)](../docs/04-explanation/architecture/contact-particle-particle.md);
- [ADR-007](../docs/04-explanation/decisions/adr-007-contact-force-types.md);
- ключи — [initial-conditions.md](../docs/03-reference/config/initial-conditions.md); тесты — [01-run-tests.md](../docs/02-how-to/01-run-tests.md).

Где в коде сейчас: формулы типов 1–3 — `include/contact_force.hpp`; продовый контакт Solid–Solid —
`include/interactions/solid_contact.hpp` (критерий пары `body[i] != body[j]`);
`src/ext_force.cpp` — только стенд для `tests/test_contact.cpp`; контакт нитей —
`include/interactions/fiber_contact*.hpp`, `fiber_solid_contact.hpp`.
Старый путь RHS (`direct_find`, `FlatGrid`, `PairGraph`, `virt_part.cpp`, `density.cpp` и др.)
удалён в a4ce0b9 (этап Ф4-B2), поэтому описания «как было» в файлах ниже — исторические.

## Документы

Статус на 2026-11-14: **план** — план/ТЗ, у файла врезка «Статус»; **отчёт** — исторический снимок,
тело не правится, при расхождении с кодом есть врезка «Примечание»; **актуален** — совпадает с кодом.

| Файл | Назначение | Статус |
|---|---|---|
| [00-audit.md](00-audit.md) | Проверка `Report_SPH_contact.md` против базы `a4147b9`: что неверно, чего не хватает | отчёт, снимок 2026-09-18 |
| [01-plan.md](01-plan.md) | План внедрения: шаги 0–6, код, тесты, чек-лист | план: реализовано с отличиями (оркестратор перенесён в `solid_contact.hpp`) |
| [02-decisions.md](02-decisions.md) | Шесть решений D1–D6 с рекомендованными вариантами | план: реализовано, D2 заменено на критерий `body[]` |
| [04-e2e-test-brief.md](04-e2e-test-brief.md) | ТЗ на e2e-тест BC2 | план: реализовано частично (`tests/test_e2e_bc2.cpp` есть, `bc2_reference_smoke.txt` нет) |
| [05-e2e-math.md](05-e2e-math.md) | Численная методика расчёта BC2 (уравнения, УРС, прочность, интегрирование) | отчёт; конвейер устарел (до Ф4-B2) |
| [06-e2e-case.md](06-e2e-case.md) | Постановка задачи BC2: геометрия, условия, результат с графиками | отчёт; постановка соответствует `tests/bc2_geometry.hpp` |
| [06-e2e-bc6.md](06-e2e-bc6.md) | Постановка задачи BC6 (откол на B4C): геометрия, метрики лесенки, результат | отчёт; постановка соответствует `tests/bc6_case.hpp` |
| [07-changes.md](07-changes.md) | Разбор изменений ветки `contact` относительно `master` (27 коммитов) | отчёт, снимок на момент слияния |
| [08-phase-shift.md](08-phase-shift.md) | Разбор расхождения BC2 с приборным эталоном (фаза, амплитуда) | отчёт, старый путь RHS |
| [09-article-conformance.md](09-article-conformance.md) | Запрос на сверку реализации со статьёй (BC2) | отчёт |
| [11-contact-report.md](11-contact-report.md) | Полный отчёт о внедрении: файлы, путь данных, отступления от статьи, параметры, результаты BC2/BC6 | отчёт; расположение кода изменилось (см. врезку) |
| [12-fiber-contact-types.md](12-fiber-contact-types.md) | Законы pinball / Swegle / импульсный для контактов нитей (Fiber–Solid, Fiber–Fiber) | отчёт; диапазон типов теперь `0..4` (см. врезку) |
| [e2e/](e2e/) | Данные и схемы BC2/BC6: профили скорости, SVG, `mkgeom.py` | данные |
| [tests/data/](../tests/data/) | Эталоны e2e-тестов BC2/BC6 (`bc2_reference.txt`, `bc6_reference.txt`) | актуален |

## Что нужно знать

- `Report_SPH_contact.md` (корень репозитория) описывал другое состояние проекта и опираться на него
  нельзя — см. [00-audit.md](00-audit.md).
- Типы контакта Solid–Solid: `contact_force_type` 0 = выключен, 1 = pinball, 2 = Свигл, 3 = импульсный
  (`contact::ForceType`, `include/contact_force.hpp`); ключи `contact_penalty_scale`,
  `contact_impulse_gamma`, `sph_same_body_only` разбираются в `src/input.cpp`.
- Контакт нитей: `fiber_solid_contact_type` и `fiber_fiber_contact_type`, диапазон `0..4`,
  **0 = Герц** (не «выключено»).
- Стартовый промпт и порядок передачи в «следующий запуск» (из прежней версии README) устарели:
  план выполнен.
