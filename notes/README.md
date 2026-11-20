# Рабочие заметки

ТЗ, планы, инструкции и отчёты, написанные по ходу работы над кодом. Это **не** документация: актуальное описание
модели и кода — в [`docs/`](../docs/index.md), журналы расчётов — в [`experiments/`](../experiments/README.md).
Заметка фиксирует состояние на дату написания; где код с тех пор изменился, верен код и `docs/`.

## Правила

- Файл кладётся в тематический каталог; имя — латиница, нижний регистр, через дефис (`law3-projection-tz.md`);
  номер в начале (`01-…`) — только если важен порядок чтения.
- Вид документа виден из имени или заголовка: ТЗ (`-tz`), план (`plan-…`, `-plan`), отчёт (`-report`), справка.
- Выполненный план или ТЗ не переносится и не удаляется: в начале файла ставится врезка
  `> **Статус (дата):** реализовано / частично / заменено …` со ссылкой на код и `docs/`, и статус в таблице ниже.
- Внешние источники и их сводки — в `sources/`.

## Каталоги

| Каталог | Тема |
|---|---|
| [`sources/`](sources/) | Внешние источники: статья о контакте в SPH, сводка статьи Guo 2019, константы Autodyn |
| [`contact/`](contact/README.md) | Контакт частица–частица (Solid): план, e2e-тесты BC2/BC6, отчёты; `e2e/` — эталон и графики |
| [`fiber-model/`](fiber-model/) | Модель нити: постановка, план этапов Ф0–Ф9, отчёты, пластичность, повреждённость |
| [`fiber-contact/`](fiber-contact/) | Контакт и трение нитей: законы 1–4, проекция трения, закон 3 как проекция скорости |
| [`fabrics/`](fabrics/) | Ткани: генератор полотна, закрепление, одеяло, мозаика, серия kT |
| [`guo2019/`](guo2019/) | Эксперименты Guo 2019 (1T, 1T/A, A/1T, 2T/4T): ТЗ, постановка, отчёты |
| [`johnson-cook/`](johnson-cook/) | Модель Джонсона–Кука: анализ, план, отчёты J0–J7/E1–E5, искусственное напряжение |
| [`jh2/`](jh2/) | JH-2: обзор кода и план исправления |
| [`speedup/`](speedup/) | Ускорение многопоточного расчёта |

## Файлы

| Файл | Вид | Статус |
|---|---|---|
| **sources/** | | |
| [sph-contact-article.md](sources/sph-contact-article.md) | источник | текст статьи об алгоритме контакта в SPH (основа ADR-007) |
| [guo2019-article-summary.md](sources/guo2019-article-summary.md) | сводка | справка |
| [steel-4340-autodyn-constants.md](sources/steel-4340-autodyn-constants.md) | справка | действует: константы STEEL 4340 приоритетнее литературы |
| **contact/** | | статусы — в [contact/README.md](contact/README.md) |
| [report-sph-contact.md](contact/report-sph-contact.md) | план | реализовано → [ADR-007](../docs/04-explanation/decisions/adr-007-contact-force-types.md) |
| `00-audit.md` … `12-fiber-contact-types.md` | ТЗ, постановки, отчёты | см. [contact/README.md](contact/README.md) |
| **fiber-model/** | | |
| [fiber-problem-statement.md](fiber-model/fiber-problem-statement.md) | постановка | реализовано частично (врезка в файле) |
| [fiber-plan.md](fiber-model/fiber-plan.md) | план | реализовано, этапы Ф0–Ф9 → [ADR-008](../docs/04-explanation/decisions/adr-008-rhs-architecture-fiber.md) |
| [fiber-implementation-report.md](fiber-model/fiber-implementation-report.md) | отчёт | исторический (2026-10-01) |
| [fiber-model-report.md](fiber-model/fiber-model-report.md) | отчёт | исторический; актуально — [глава 4](../docs/04-explanation/models/04-fiber-mechanics.md) |
| [fiber-plasticity-design.md](fiber-model/fiber-plasticity-design.md) | дизайн | реализовано |
| [fiber-damage-plan.md](fiber-model/fiber-damage-plan.md) | план | реализовано частично |
| **fiber-contact/** | | актуально — [глава 5](../docs/04-explanation/models/05-fiber-contact.md) |
| [friction-report.md](fiber-contact/friction-report.md) | отчёт | исторический (задачи 1–4 трения) |
| [fiber-friction-implementation.md](fiber-contact/fiber-friction-implementation.md) | инструкция | заменено проекцией трения |
| [fiber-friction-projection.md](fiber-contact/fiber-friction-projection.md) | инструкция | реализовано |
| [fiber-solid-friction-implementation.md](fiber-contact/fiber-solid-friction-implementation.md) | инструкция | реализовано |
| [friction-chattering-analysis.md](fiber-contact/friction-chattering-analysis.md) | анализ | исторический |
| [fiber-contact-segment-pairs.md](fiber-contact/fiber-contact-segment-pairs.md) | инструкция | реализовано |
| [fiber-solid-contact-fix.md](fiber-contact/fiber-solid-contact-fix.md) | инструкция | реализовано |
| [plan-fiber-contact-types.md](fiber-contact/plan-fiber-contact-types.md) | план | реализовано, с расширением |
| [plan-fiber-contact-laws12-fix.md](fiber-contact/plan-fiber-contact-laws12-fix.md) | план | не реализовано (заменено) |
| [law3-projection-tz.md](fiber-contact/law3-projection-tz.md) | ТЗ | реализовано (проекция нормального контакта) |
| [law3-projection-report.md](fiber-contact/law3-projection-report.md) | отчёт | исторический |
| **fabrics/** | | актуально — [глава 7](../docs/04-explanation/models/07-fabrics.md) |
| [plain-weave-generator-tz.md](fabrics/plain-weave-generator-tz.md) | ТЗ | реализовано |
| [clamp-solid-fiber-tz.md](fabrics/clamp-solid-fiber-tz.md) | ТЗ | реализовано |
| [clamp-tests.md](fabrics/clamp-tests.md) | план тестов | реализовано частично |
| [plan-quilt-fabric.md](fabrics/plan-quilt-fabric.md) | план | реализовано (стёганое одеяло) |
| [tz-random-cut-quilt.md](fabrics/tz-random-cut-quilt.md) | ТЗ | реализовано (мозаичная ткань) |
| [fabric-kt-penetration-tz.md](fabrics/fabric-kt-penetration-tz.md) | ТЗ серии расчётов | расчёты — `experiments/2026-11-09-fabric-kT-ball1g-v550/` и далее |
| **guo2019/** | | |
| [guo2019-experiments-tz.md](guo2019/guo2019-experiments-tz.md) | ТЗ | реализовано (и расширено) |
| [guo2019-2t-4t-tz.md](guo2019/guo2019-2t-4t-tz.md) | ТЗ | реализовано (2T, 4T, затем `<k>T`) |
| [guo2019-e2e.md](guo2019/guo2019-e2e.md) | постановка | действует |
| [exp-1t-single-ply.md](guo2019/exp-1t-single-ply.md), [exp-1t-a-fabric-front.md](guo2019/exp-1t-a-fabric-front.md), [exp-a-1t-metal-front.md](guo2019/exp-a-1t-metal-front.md) | описания опытов | справка |
| [guo2019-plastic-1t-2t-report.md](guo2019/guo2019-plastic-1t-2t-report.md) | отчёт | исторический (2026-10-16) |
| **johnson-cook/** | | актуально — [ADR-009](../docs/04-explanation/decisions/adr-009-johnson-cook.md) |
| [jc-analysis.md](johnson-cook/jc-analysis.md), [jc-implementation-plan.md](johnson-cook/jc-implementation-plan.md) | анализ, план | реализовано |
| `00-baseline.md` … `06-artificial-stress.md` | отчёты | исторические |
| [handoff.md](johnson-cook/handoff.md) | передача контекста | снимок на 2026-11-09 |
| **jh2/** | | |
| [jh2-code-review.md](jh2/jh2-code-review.md), [jh2-fix-plan.md](jh2/jh2-fix-plan.md) | обзор, план | реализовано → [ADR-010](../docs/04-explanation/decisions/adr-010-jh2-fix.md) |
| **speedup/** | | |
| `base.md`, `memory.md`, `reduce-targets.md`, `sph-contrib.md`, `prims-parallel.md`, `final.md`, `contact.md` | отчёты | исторические; актуально — [параллелизм](../docs/04-explanation/architecture/parallelism.md) |
