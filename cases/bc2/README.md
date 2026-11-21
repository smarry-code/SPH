# Кейс BC2

1. **Опыт.** Удар пластины по карбиду бора (B₄C) с регистрацией скорости свободной поверхности
   лазерным интерферометром; в расчёте — осевая скорость частицы-монитора.
2. **Источник.** Эталон тестов — приборный профиль
   [`notes/contact/e2e/gold_standart.txt`](../../notes/contact/e2e/gold_standart.txt) (Grady & Moody,
   SAND96-0551; начало отсчёта времени произвольно). Происхождение `source/reference-curve.txt` не
   установлено, в тестах он не используется.
3. **Статус.** В коде и тестах.
4. **Файлы.**
   - [`source/initial-conditions.txt`](source/initial-conditions.txt) — параметры решателя старого кода (формат `Initial_Conditions.txt`): `t_total = 3e-6`, `save_t = 2.5e-8`; геометрии, материалов и скорости удара в нём нет.
   - [`source/reference-curve.txt`](source/reference-curve.txt) — эталонная кривая: 340 точек, две колонки (время, с; скорость, м/с), время 1,88–3,49 мкс.
5. **Код и тесты.**
   - **Регрессия** (замороженная копия для `ctest`, [правила](../README.md#роли-постановок)) —
     [`regression/`](regression/): геометрия [`bc2_geometry.hpp`](regression/bc2_geometry.hpp) (перенесена
     из старого `src/input.cpp`), уменьшенный 3D-вариант
     [`bc2_mini3d_geometry.hpp`](regression/bc2_mini3d_geometry.hpp), параметры решателя
     [`initial_conditions.txt`](regression/initial_conditions.txt) (значения из `source/initial-conditions.txt`),
     эталон [`reference.txt`](regression/reference.txt), снимки регрессии [`snapshots/`](regression/snapshots/)
     (`bc2_*`, `bc2mini3d_*`).
   - **Эксперимент** (живая копия, ответвлена от `regression/` в 265f959) — [`experiment/`](experiment/):
     геометрия [`bc2_geometry.hpp`](experiment/bc2_geometry.hpp), [`initial_conditions.txt`](experiment/initial_conditions.txt),
     [`reference.txt`](experiment/reference.txt) и утилита [`bc2_shot.cpp`](experiment/bc2_shot.cpp)
     (цель `bc2_shot`, драйвер [`cases/common/ceramic_shot.hpp`](../common/ceramic_shot.hpp)); запуск и ключи —
     [docs/05-experiments/ceramics-bc2.md](../../docs/05-experiments/ceramics-bc2.md).
   - Тесты — [`tests/test_e2e_bc2.cpp`](../../tests/test_e2e_bc2.cpp), `test_e2e_3d`, `test_regression`; цели `ctest`:
     `e2e_bc2_smoke_contact0…3`, `e2e_bc2_full_contact0…3`, `e2e_3d_smoke_contact0…3`,
     `regression_bc2_contact0…3`, `regression_bc2mini3d_contact0…3`; конфиг BC2 берут также
     `test_e2e_fiber_solid`, `test_fiber_time_integration_equivalence`, `test_sph_body_isolation`.
6. **Эксперименты.**
   - [`2026-09-19-bc2-discrepancy`](../../experiments/2026-09-19-bc2-discrepancy/README.md)
   - [`2026-09-20-bc2-abort-contact-type2`](../../experiments/2026-09-20-bc2-abort-contact-type2/README.md)
   - [`2026-09-22-bc2-contact-types`](../../experiments/2026-09-22-bc2-contact-types/README.md)
   - [`2026-09-22-bc2-kp-sweep`](../../experiments/2026-09-22-bc2-kp-sweep/README.md)
   - [`2026-09-22-bc2-resolution-scheme`](../../experiments/2026-09-22-bc2-resolution-scheme/README.md)
   - [`2026-11-21-bc2-shot-contact3`](../../experiments/2026-11-21-bc2-shot-contact3/README.md)
   - [`2026-11-21-bc2-shot-contact-types`](../../experiments/2026-11-21-bc2-shot-contact-types/README.md)
7. **Постановка.** [`04-e2e-test-brief.md`](../../notes/contact/04-e2e-test-brief.md),
   [`05-e2e-math.md`](../../notes/contact/05-e2e-math.md), [`06-e2e-case.md`](../../notes/contact/06-e2e-case.md),
   [`08-phase-shift.md`](../../notes/contact/08-phase-shift.md).
