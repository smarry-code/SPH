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
   - Геометрия — [`tests/bc2_geometry.hpp`](../../tests/bc2_geometry.hpp) (перенесена из старого
     `src/input.cpp`), уменьшенный 3D-вариант — [`tests/bc2_mini3d_geometry.hpp`](../../tests/bc2_mini3d_geometry.hpp).
   - Параметры решателя — [`tests/data/bc2_initial_conditions.txt`](../../tests/data/bc2_initial_conditions.txt)
     (значения из `source/initial-conditions.txt`), эталон — [`tests/data/bc2_reference.txt`](../../tests/data/bc2_reference.txt).
   - Тест — [`tests/test_e2e_bc2.cpp`](../../tests/test_e2e_bc2.cpp); цели `ctest`:
     `e2e_bc2_smoke_contact0…3`, `e2e_bc2_full_contact0…3`, `e2e_3d_smoke_contact0…3`,
     `regression_bc2_contact0…3`, `regression_bc2mini3d_contact0…3`.
6. **Эксперименты.**
   - [`2026-09-19-bc2-discrepancy`](../../experiments/2026-09-19-bc2-discrepancy/README.md)
   - [`2026-09-20-bc2-abort-contact-type2`](../../experiments/2026-09-20-bc2-abort-contact-type2/README.md)
   - [`2026-09-22-bc2-contact-types`](../../experiments/2026-09-22-bc2-contact-types/README.md)
   - [`2026-09-22-bc2-kp-sweep`](../../experiments/2026-09-22-bc2-kp-sweep/README.md)
   - [`2026-09-22-bc2-resolution-scheme`](../../experiments/2026-09-22-bc2-resolution-scheme/README.md)
7. **Постановка.** [`04-e2e-test-brief.md`](../../notes/contact/04-e2e-test-brief.md),
   [`05-e2e-math.md`](../../notes/contact/05-e2e-math.md), [`06-e2e-case.md`](../../notes/contact/06-e2e-case.md),
   [`08-phase-shift.md`](../../notes/contact/08-phase-shift.md).
