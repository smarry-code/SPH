# Контекст для продолжения работы над моделью Джонсона-Кука

Дата: 2026-11-09. Ветка `fiber_johnson_cook` (на `origin`, GitHub `smarry-code/SPH`).
Файл для новой сессии Claude Code: прочитать целиком перед началом работы.

> **Примечание (2026-11-14):** документ — снимок на 2026-11-09. Хеши коммитов ниже — из ветки `fiber_johnson_cook`; в `fiber_new`/`documentation` те же изменения лежат под другими хешами (напр. `c2733a0` + `51f9371` вместо `6be23bb` + `e17ea4d`; `887d6b8` — отчёт 05, ADR-009, HANDOFF). Ключ `axisym_pair_mass` откатан и в коде отсутствует; ключи `xsph_mode`, `jc_heating`, `jc_chi` и каналы мощности лога энергии (`include/energy_log.hpp`) в коде есть (`src/input.cpp`, `include/simulation_config.hpp`). Старый путь RHS (`density.cpp`, `int_force()` и т. п.) удалён в a4ce0b9. Сведения о сервере (2 ядра, правила про `claude-compute-8t`) верны на 2026-11-09 и могли измениться. Актуальное описание — [ADR-009](../docs/04-explanation/decisions/adr-009-johnson-cook.md), [how-to](../docs/02-how-to/06-johnson-cook-material.md), [запуск тестов](../docs/02-how-to/01-run-tests.md).

## 1. Рабочая машина и сборка

- Windows 11, репозиторий `C:\dip\SPH\SPH`.
- Сборка через MSYS2 clang64: `C:/msys64/clang64/bin/clang++.exe`. Пресет `clang-windows` из `CMakePresets.json` не работает без Windows SDK.
- Рабочая команда (Git Bash; каталоги `build-*` в `.gitignore`):
  ```
  export PATH=/c/msys64/clang64/bin:$PATH
  cmake -G Ninja -B build-jc2e -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER=C:/msys64/clang64/bin/clang.exe \
        -DCMAKE_CXX_COMPILER=C:/msys64/clang64/bin/clang++.exe \
        -DSPH_DIM=2 -DSPH_ENERGY_LOG=ON
  cmake --build build-jc2e -j 12
  ```
  Полная сборка ~3.5 мин на 14 ядрах. Локально — только сборка и чистые unit-тесты (`jc_*`, `axisym`, `particle_system`);
  `test_checkpoint` гоняет шаги SPH, его — на сервере.
- git: `user.name` не настроен. Коммиты так:
  `git -c user.name=smarry-code -c user.email=<адрес как в прежних коммитах> commit ...`
- Обязательно прочитать:
  1. `fiber/jc-implementation-plan.md` (ТЗ), `fiber/jc-analysis.md`;
  2. `johnson-cook/00-baseline.md` ... `04-pdv-energy-j7-bc6.md` по порядку;
  3. `docs/04-explanation/decisions/adr-009-johnson-cook.md`.

## 2. Правила владельца (обязательны)

- **Не запускать расчёты локально.** Любой прогон SPH, включая короткие e2e и переборы, идёт на `ssh claude-compute`. Локально — только компиляция и чистые unit-тесты.
- **`claude-compute-8t` не трогать**: он занят другой задачей (`guo2019_shot`).
- Код на сервер доставляется через `git bundle` + `scp`. На сервере нет доступа к GitHub.
- Коммиты делать по ходу работы. `git push` и PR — только по просьбе.
- В `experiments/` не коммитить картинки (gif/svg/png).
- Константы материалов не подгонять под ворота тестов. Непроверенные значения помечать `TODO(verify)`. Заменять их только литературой со ссылкой на источник.
- **Данных Autodyn не будет** (решение 2026-11-04).
- Подзадачи можно отдавать агентам Sonnet/Haiku.
- Язык общения и документов — русский.

## 3. Что сделано

| Этап | Содержание | Коммит |
|---|---|---|
| J1+J2a | `include/strength/johnson_cook.hpp`, УРС Linear/Shock, U1–U17 | a64075f |
| J2b+J3 | материалы itype 12–15, диспетчеризация УРС, поля `T`, `edot`, чекпоинт v2 | a2060b9 |
| J4+J5 | JC-прочность и разрушение в конвейере, I1–I8 | 32f2032, 7733304 |
| J0, E1 | baseline, регрессия BC2/BC6 | 75e4b16 |
| E2–E5 | HEL, откол (1D), Тейлор (2D осесимметрия), сходимость | b7801b8, 2645142, a3a8ddb |
| p_min | `JcMaterial::p_min`, −4 ГПа для 4340, U20 | 0de8845, 4a85ec1 |
| разбор энергии E4 | `03-pmin-energy.md` | 1fff62e |
| ТЗ в репозитории | `fiber/jc-*.md` | 28e19c2 |
| HANDOFF (прежний) | | 9cfbf86 |
| pdV | член −(p/ρ)·v_r/r в `du`, тест G | a686ee9 |
| вес энергии | `energy_mass = 2π·r0·m`, поле `r0`, чекпоинт v3 | 9cc2266 |
| эталоны | снимки BC2/BC6 N = 100 (8 файлов `*_step100`) | 88a3b2f |
| чекпоинт Windows | `_commit`, ротация `.prev` | 54b61a9 |
| диагностика | `--cfl` в тесте Тейлора | 8639e59 |
| J7a | Тейлор-Куинни, I9, предупреждение Tm_f ≤ 0 | 89555c2 |
| J7b | ADR-009, how-to, реестр, changelog | 04dbe6f |
| рестарт | JC 4340 + TQ в `test_checkpoint` | e638823 |
| G4 | `energy_mass` в `test_axisym` | cf9133e |
| ADR | XSPH и зависимость от dt | 78bb0ab |
| BC6 | устойчивый извлекатель плато 1 | ed5ac88 |
| BC6 отчёт | оговорка об извлекателе; бюджет clamp full 12000/14400 с | fd76598 |
| XSPH | ключ `xsph_mode = velocity \| position` | 30f76b2 |
| G5 | диагностика баланса пар в `test_axisym` | 4fe6e0d |
| ring | `axisym_pair_mass = ring` и откат (Тейлор хуже); в коде ключа нет | 6be23bb, e17ea4d |
| лог энергии | каналы мощности по механизмам, `closure` | f99b11d |

Подробности и числа — `04-pdv-energy-j7-bc6.md`, `05-xsph-axisym-energy.md`.

Состояние (сервер, `cf9133e`): быстрые тесты 112/112 в `build1/2/3`; `ctest -L regression` 14/14; E4 — энергия в режиме `--energy-warn-only`; `e2e_clamp_ES2_pc_full` PASS при `--time-limit=20000` (10 568 с). Полные прогоны BC6 после `ed5ac88` — раздел 7 в `04-...md`. На `f99b11d`: выборочные тесты 35/35.

## 4. Открытые вопросы

1. **Tensile instability / искусственное напряжение Монагана** (следующая задача). Без «скоростного» XSPH невязка энергии растёт до ~9 % к 10 мкс и обрыв по ρ<0 у оси (не зависит от dt и ε); вероятная причина — неустойчивость при растяжении, не доказана. Шаги: напряжение по ключу, повтор матрицы (ε=0 / position / CFL), при успехе — XSPH по умолчанию `position` или меньше ε, перегенерация эталонов. Литература: Monaghan 2000; Gray, Monaghan, Swift 2001 (не сверено). См. `05-xsph-axisym-energy.md`.
2. **XSPH зависит от dt**: ключ `position` есть, но без стабилизатора неустойчив (BC2/BC6 обрываются).
3. **E4 n_r = 24: 2.25 % > 2 %**, остаётся `--energy-warn-only`.
4. **Ворота D = 0 в BC6** — краевой артефакт (частица #11388), жёсткие по решению владельца.
5. **Литература**: откольная прочность 4340 вместо p_min = −4 ГПа; эталон Тейлора (`tests/data/jc/README.md`); значения `TODO(verify)`. Только с источником.
6. **Push и PR** — только по просьбе владельца.

## 5. Сервер `claude-compute`

- 2 ядра, 7.7 ГБ, `OMP_NUM_THREADS=2`. Связь иногда рвётся: всё запускать через `setsid nohup` и очереди со `STATUS`. Образец — `~/jc_runs/queue_cf9133e.sh`.
- `~/SPH_jc` — на HEAD. Сборки `build1/2/3` (SPH_DIM = 1/2/3) и `build3e` (SPH_ENERGY_LOG=ON), Ninja, Release.
- Каталоги прогонов `~/jc_runs/`: `q88`, `q_cf9`, `xsph_dt`, `bc6_bisect`, `q_bc6`, `pmin_pdv` (прошлые сессии), `q_diag`, `q_ring`, `q_hoop` (XSPH, ring, разложение мощности).

Обновить код на сервере:
```
git bundle create /tmp/jc.bundle <последний_коммит_на_сервере>..fiber_johnson_cook
scp /tmp/jc.bundle claude-compute:~/
ssh claude-compute 'cd ~/SPH_jc && git pull ~/jc.bundle fiber_johnson_cook'
```

Полезные запуски:
```
ctest --test-dir build1 -L jc --output-on-failure
cmake -B build2 -DSPH_E2E_SLOW_TESTS=ON && ctest --test-dir build2 -R e2e_jc_taylor_e
```

## 6. Технические заметки

- В 2D у решателя всегда осесимметрия (`axis_sim` compile-time), масса «плоская»: m = ρ0·dx². Физическая масса кольца 2π·r0·m; лог энергии теперь её учитывает.
- XSPH меняет саму скорость, поэтому реально рассеивает кинетическую энергию (канал `D_xsph`), и его вклад зависит от dt.
- Режим нагрева JC: `jc_heating = energy` (по умолчанию) или `taylor_quinney` (`jc_chi`, 0.9).
- Извлекатель ступеней BC6 — `tests/e2e_staircase.hpp`.
- У `e2e::ScopedWorkDir` в имени каталога случайный суффикс.
- Тест Тейлора печатает энергию в трёх вариантах (без веса, вес r0, вес r), см. `tests/jc_taylor.hpp`.
