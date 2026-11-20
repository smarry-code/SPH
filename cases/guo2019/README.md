# Кейс Guo 2019

1. **Опыт.** Удар стального шара (12 мм, 7 г) по мишеням из ткани Twaron CT 717 и пластины
   Al 6061-T6: A, 1T, 1T/A, A/1T, 2T, 4T. В коде есть ещё `<k>T`, пакет из k слоёв ткани (до 64).
2. **Источник.** Guo et al., International Journal of Impact Engineering 127 (2019) 1–16. Статья в
   репозитории не хранится, сводка —
   [`notes/sources/guo2019-article-summary.md`](../../notes/sources/guo2019-article-summary.md).
   Данные статьи (скорости выстрелов, V_BL, поглощённая энергия, число порванных нитей) записаны
   константами с пометкой источника в `guo2019_case.hpp`.
3. **Статус.** В коде и тестах.
4. **Файлы.** Каталога `source/` нет: кроме статьи, исходных данных нет.
5. **Код и тесты.**
   - **Регрессия** (замороженная копия для `ctest`, [правила](../README.md#роли-постановок)) —
     [`regression/`](regression/): числа опыта и допуски [`guo2019_case.hpp`](regression/guo2019_case.hpp),
     сцена [`guo2019_scene.hpp`](regression/guo2019_scene.hpp), монитор
     [`guo2019_monitor.hpp`](regression/guo2019_monitor.hpp), параметры решателя
     [`initial_conditions.txt`](regression/initial_conditions.txt). Её берут
     [`tests/test_e2e_guo2019.cpp`](../../tests/test_e2e_guo2019.cpp) (цели `e2e_guo2019_smoke_A`, `_1T`,
     `_1T_A`, `_A_1T`, `_2T`, `_4T`), `test_energy_log`, `test_fiber_plasticity_guo`, `test_plain_weave_quilt`.
   - **Эксперимент** — [`experiment/`](experiment/): те же файлы (ответвлены от `regression/` в `a9c3964`) и
     драйверы: серийные выстрелы [`guo2019_shot.cpp`](experiment/guo2019_shot.cpp)
     ([справочник](../../docs/03-reference/guo2019-shot.md); ткань, раскладку и параметры каждого выстрела
     пишет в его каталог), дамп нитей [`fiber_t0_dump.cpp`](experiment/fiber_t0_dump.cpp). Её же подключает
     `tools/fabric_checks`. Проверки этих утилит в `ctest` (`checkpoint_cli_test`, `stitch_kt_cli_test`,
     `mosaic_kt_cli_test`, `fabric_checks_smoke_*`) идут на экспериментальной копии.
   - Отличия копий: `diff -r cases/guo2019/regression cases/guo2019/experiment`.
6. **Эксперименты.** Серии `fabric-kT` и `mosaic-kT` (шар 1 г, 550 м/с) идут на той же сцене через
   `guo2019_shot --config=<k>T`.
   - [`2026-10-01-guo2019-1t-scaled`](../../experiments/2026-10-01-guo2019-1t-scaled/README.md)
   - [`2026-10-02-guo2019-1t-base`](../../experiments/2026-10-02-guo2019-1t-base/README.md)
   - [`2026-10-12-fiber-plasticity-guo`](../../experiments/2026-10-12-fiber-plasticity-guo/README.md)
   - [`2026-10-16-guo2019-1t-contact-types-friction`](../../experiments/2026-10-16-guo2019-1t-contact-types-friction/README.md)
   - [`2026-10-16-guo2019-fiber-contact-types`](../../experiments/2026-10-16-guo2019-fiber-contact-types/README.md)
   - [`2026-10-17-guo2019-1t-pinball-vbl`](../../experiments/2026-10-17-guo2019-1t-pinball-vbl/)
   - [`2026-10-18-guo2019-2t-law4`](../../experiments/2026-10-18-guo2019-2t-law4/)
   - [`2026-10-18-guo2019-law4-summary`](../../experiments/2026-10-18-guo2019-law4-summary/README.md)
   - [`2026-10-22-guo2019-law3-summary`](../../experiments/2026-10-22-guo2019-law3-summary/README.md)
   - [`2026-10-26-guo2019-damage-summary`](../../experiments/2026-10-26-guo2019-damage-summary/README.md)
   - [`2026-11-02-guo2019-thread-divergence`](../../experiments/2026-11-02-guo2019-thread-divergence/)
   - [`2026-11-07-guo2019-law3-projection`](../../experiments/2026-11-07-guo2019-law3-projection/README.md)
   - [`2026-11-08-guo2019-1t-law3push-mu027`](../../experiments/2026-11-08-guo2019-1t-law3push-mu027/README.md)
   - [`2026-11-09-fabric-kT-ball1g-v550`](../../experiments/2026-11-09-fabric-kT-ball1g-v550/)
   - [`2026-11-11-guo2019-law3-performance`](../../experiments/2026-11-11-guo2019-law3-performance/README.md)
   - [`2026-11-12-guo2019-contact-ad-bench`](../../experiments/2026-11-12-guo2019-contact-ad-bench/README.md)
   - [`2026-11-12-guo2019-contact-profile-compute`](../../experiments/2026-11-12-guo2019-contact-profile-compute/README.md)
   - [`2026-11-14-fabric-kT-w50-ball1g-v550`](../../experiments/2026-11-14-fabric-kT-w50-ball1g-v550/README.md)
   - [`2026-11-14-guo2019-2t-mosaic`](../../experiments/2026-11-14-guo2019-2t-mosaic/README.md)
   - [`2026-11-15-mosaic-kT-ball1g-v550`](../../experiments/2026-11-15-mosaic-kT-ball1g-v550/README.md)
7. **Постановка.** ТЗ [`guo2019-experiments-tz.md`](../../notes/guo2019/guo2019-experiments-tz.md)
   и [`guo2019-2t-4t-tz.md`](../../notes/guo2019/guo2019-2t-4t-tz.md), постановка
   [`guo2019-e2e.md`](../../notes/guo2019/guo2019-e2e.md), описания опытов
   [`exp-1t-single-ply.md`](../../notes/guo2019/exp-1t-single-ply.md),
   [`exp-1t-a-fabric-front.md`](../../notes/guo2019/exp-1t-a-fabric-front.md),
   [`exp-a-1t-metal-front.md`](../../notes/guo2019/exp-a-1t-metal-front.md).
