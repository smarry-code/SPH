# Кейс BC6

1. **Опыт.** Откольный опыт на карбиде бора (B₄C): скорость границы «керамика / окно» регистрирует
   интерферометр (VISAR); в расчёте — осевая скорость частицы-монитора на тыльном слое образца.
2. **Источник.** `source/reference-curve.txt` — оцифровка опубликованного профиля опыта BC6; точка в
   точку совпадает с таблицей реперных точек §5.2 [`06-e2e-bc6.md`](../../notes/contact/06-e2e-bc6.md).
   Начало отсчёта времени в первоисточнике произвольно; передний фронт оцифрован облаком точек.
3. **Статус.** В коде и тестах.
4. **Файлы.**
   - [`source/initial-conditions.txt`](source/initial-conditions.txt) — параметры решателя старого кода (формат `Initial_Conditions.txt`): `t_total = 3.0e-6`, `save_t = 5e-8`; геометрии, материалов и скорости удара в нём нет.
   - [`source/reference-curve.txt`](source/reference-curve.txt) — эталонная кривая: 236 точек, две колонки (время, с; скорость, м/с), время 1,08–3,22 мкс.
   - [`source/reference-curve.jpg`](source/reference-curve.jpg) — график к кривой.
5. **Код и тесты.**
   - Числа опыта — [`tests/bc6_case.hpp`](../../tests/bc6_case.hpp), геометрия —
     [`tests/bc6_geometry.hpp`](../../tests/bc6_geometry.hpp).
   - Параметры решателя — [`tests/data/bc6_initial_conditions.txt`](../../tests/data/bc6_initial_conditions.txt),
     эталон — [`tests/data/bc6_reference.txt`](../../tests/data/bc6_reference.txt); его строит из
     `source/reference-curve.txt` скрипт [`tools/bc6_make_reference.py`](../../tools/bc6_make_reference.py)
     (запускается вручную).
   - Тест — [`tests/test_e2e_bc6.cpp`](../../tests/test_e2e_bc6.cpp); цели `ctest`:
     `e2e_bc6_smoke_contact0…3`, `e2e_bc6_full_contact0…3`, `regression_bc6_contact0…3`.
6. **Эксперименты.**
   - [`2026-09-23-bc6-art-visc`](../../experiments/2026-09-23-bc6-art-visc/README.md)
   - [`2026-09-23-bc6-butt-joint`](../../experiments/2026-09-23-bc6-butt-joint/README.md)
   - [`2026-09-23-bc6-resolution`](../../experiments/2026-09-23-bc6-resolution/README.md)
   - [`2026-09-23-bc6-xsph`](../../experiments/2026-09-23-bc6-xsph/README.md)
   - [`2026-09-24-bc6-contact-dx010`](../../experiments/2026-09-24-bc6-contact-dx010/README.md)
7. **Постановка.** [`06-e2e-bc6.md`](../../notes/contact/06-e2e-bc6.md).
