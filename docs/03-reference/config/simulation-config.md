# Справочник `SimulationConfig`

Объект [`SimulationConfig`](../../../include/simulation_config.hpp:31)
инкапсулирует runtime-параметры расчёта. Значения по умолчанию задаются в
[`apply_defaults()`](../../../include/simulation_config.hpp:189) (и дублируются
инициализаторами полей) и могут быть переопределены в
[`Initial_Conditions.txt`](initial-conditions.md). Колонка «Ключ» ниже показывает,
читается ли поле из файла: «да» — одноимённый ключ, «нет» — поле задаётся только
в коде (тесты, утилиты) или вычисляется.

## Поля

### SPH-метод и алгоритмы

| Поле | Тип | По умолчанию | Ключ | Описание |
|---|---|---|---|---|
| `pa_sph` | int | `2` | да | Метод аппроксимации давления (1/2) |
| `nnps` | int | `1` | да | Метод поиска соседей (1=direct, 2=link_list, 3=tree); **в конвейере не используется** |
| `sle` | int | `1` | да | Алгоритм обновления длины сглаживания (0 — не менять, 1 — $1.3\,(m/\varrho)^{1/\dim}$; 2 — ошибка выполнения) |
| `skf` | int | `1` | да | Тип ядра сглаживания |
| `density_method` | int | `1` | да | Метод расчёта плотности (1/2/3) |
| `coloring_niac_threshold` | int | `0` | да | Порог graph coloring (0 = всегда atomic); **в конвейере не используется** |
| `ode_solver` | int | `0` | да | Решатель: 0 — предиктор-корректор, 1 — явный Эйлер (эталон для тестов порядка, не для рабочих расчётов) |
| `sph_same_body_only` | bool | `false` | да | Отвергать SPH-пары частиц разных тел; тела взаимодействуют только через контакт (требует `contact_force_type != 0`) |

### Контакт «частица-частица» ([ADR-007](../../04-explanation/decisions/adr-007-contact-force-types.md))

| Поле | Тип | По умолчанию | Ключ | Описание |
|---|---|---|---|---|
| `contact_force_type` | int | `0` | да | 0 — выключен, 1 — пинбол Белицко, 2 — Свигл, 3 — импульсный |
| `contact_penalty_scale` | double | `1.0` | да | Масштабный коэффициент $K_p$ штрафной силы (типы 1-2) |
| `contact_impulse_gamma` | double | `1.0` | да | Доля нормальной скорости сближения, снимаемая за шаг (тип 3), рабочий диапазон $(0;\,1.6]$ |

### Контакты нитей (закон нормальной силы)

`0` — **Герц**, а не «выключено» (в отличие от `contact_force_type`); 1 — пинбол,
2 — Свигл, 3 — импульсный, 4 — Герц + демпфер. Все поля читаются одноимёнными ключами.

| Поле | Тип | По умолчанию | Описание |
|---|---|---|---|
| `fiber_solid_contact_type` | int | `0` | Закон для Fiber-Solid |
| `fiber_fiber_contact_type` | int | `0` | Закон для Fiber-Fiber |
| `fiber_solid_contact_penalty_scale` | double | `1.0` | $K_p$ (типы 1, 2) для Fiber-Solid |
| `fiber_fiber_contact_penalty_scale` | double | `1.0` | $K_p$ (типы 1, 2) для Fiber-Fiber |
| `fiber_solid_contact_impulse_gamma` | double | `1.0` | $\gamma$ типа 3 для Fiber-Solid, $(0;\,1.6]$ |
| `fiber_fiber_contact_impulse_gamma` | double | `1.0` | $\gamma$ типа 3 для Fiber-Fiber, $(0;\,1.6]$ |
| `fiber_solid_contact_impulse_push_tau` | double | `0.0` | Время выталкивания $\tau$ типа 3 (с); `0` — выключено |
| `fiber_fiber_contact_impulse_push_tau` | double | `0.0` | То же для Fiber-Fiber |
| `fiber_contact_impulse_push_depth` | double | `0.1` | Глубина выталкивания — доля суммы радиусов пары |
| `fiber_contact_impulse_push_vmax` | double | `0.0` | Предел скорости выталкивания, м/с; `0` — без предела |
| `fiber_fiber_contact_predictor` | int | `0` | `1` — нормальный контакт Fiber-Fiber законов 1–3 как проекция скорости (стадия между выбором $\Delta t$ и трением); `0` — прежняя явная сила, побитово. Требует закона 1–3 |
| `fiber_solid_contact_predictor` | int | `0` | То же для Fiber-Solid |
| `fiber_contact_projection_iters` | int | `8` | Число итераций Якоби $K$, `1..100` |
| `fiber_contact_projection_warm` | int | `1` | Тёплый старт импульсов (`1`) или холодный (`0`, предупреждение) |

Правила проверки ([`validate_fiber_contact_projection`](../../../include/simulation_config.hpp:224),
исключение `std::runtime_error`; при обоих `*_predictor = 0` не выполняется; вызывается из конца
`parse_sph_config`, из `guo2019_shot` и в начале `run_time_integration`, поэтому предупреждения
печатаются один раз за процесс): `*_predictor` — только 0/1;
`predictor = 1` требует закона контакта 1, 2 или 3; при `predictor = 1` значение
`*_contact_impulse_push_tau > 0` — ошибка; `*_contact_impulse_gamma != 1` — предупреждение, $\gamma$ игнорируется
(проекция требует $\gamma = 1$); `iters` вне 1..100 и `warm` не 0/1 — ошибка; `warm = 0` — предупреждение.
Все четыре поля сверяются при рестарте (менять только с `restart_force`). Описание метода —
[Контакт нитей, §7](../../04-explanation/models/05-fiber-contact.md#7-нормальный-контакт-как-проекция-скорости).

### Шаг по времени и вывод

| Поле | Тип | По умолчанию | Ключ | Описание |
|---|---|---|---|---|
| `CFL` | double | `0.1` | да | Число Куранта |
| `dt_scale` | double | `1.0` | да | Множитель общего шага, $(0;\,1]$ (проверка сходимости по шагу) |
| `t_total` | double | `3e-6` | да | Полное время расчёта, с |
| `save_t` | double | `2.5e-8` | да | Интервал сохранения кадров, с |
| `print_step` | size_t | `50` | да | Шаг вывода прогресса |

### Вязкость / диффузия / усреднение

| Поле | Тип | По умолчанию | Ключ | Описание |
|---|---|---|---|---|
| `alpha` | double | `1.0` | да | Линейный коэффициент псевдовязкости |
| `beta` | double | `1.0` | да | Квадратичный коэффициент псевдовязкости |
| `etq` | double | `0.01` | да | Параметр отсечки псевдовязкости |
| `epsilon` | double | `0.3` | да | Усреднение скоростей (XSPH) |
| `xsph_mode` | int | `0` | да (`velocity`/`position`) | Куда входит поправка XSPH: 0 — в скорость (прежнее), 1 — только в перенос позиций (классический XSPH); только предиктор-корректор |
| `art_stress_eps` | double | `0.0` | да | Коэффициент $\varepsilon$ искусственного напряжения Монагана (только твёрдые, `itype` > 2); `0` — выключено (побитово прежнее), типично `0.3`; $\ge 0$. Для каждого растягивающего главного напряжения $\sigma'_k>0$ добавка $R'_k=-\varepsilon\sigma'_k$. В осесимметрии учитывается окружная компонента |
| `art_stress_n` | double | `4.0` | да | Показатель $n$ в $f^n$, $f = W(r)/W(\Delta p)$, $\Delta p=(m/\rho_0)^{1/\dim}$; $> 0$ |
| `xi` | double | `0.001` | нет | Интенсивность диффузии плотности |
| `epsilon_h` | double | `0.01` | нет | Регуляризатор знаменателя `psi_ij` |

### Нагрев Джонсона–Кука

| Поле | Тип | По умолчанию | Ключ | Описание |
|---|---|---|---|---|
| `jc_heating` | int | `0` | да (`energy`/`taylor_quinney`) | 0 — $T = T_r + u/c_p$, 1 — Тейлор–Куинни (температура интегрируется отдельно) |
| `jc_chi` | double | `0.9` | да | Доля пластической работы $\chi$, идущая в тепло, $[0;\,1]$ |

Применяются к `cfg.mats.jc[*].th` в конце `parse_sph_config`. Подробнее —
[материал Джонсона–Кука](../../02-how-to/06-johnson-cook-material.md).

### Импорт геометрии и нити

Все поля читаются одноимёнными ключами (см. [формат `Initial_Conditions.txt`](initial-conditions.md)).

| Поле | Тип | По умолчанию | Описание |
|---|---|---|---|
| `fiber_file` | string | `""` | Файл нитей `.fib`; пусто — нитей нет; требует `SPH_DIM=3` |
| `k_file` | string | `""` | Файл LS-DYNA `.k`; пусто — импорта нет |
| `k_mode` | string | `"solid"` | `"solid"` или `"solid_fiber"` |
| `k_parts_file` | string | `""` | Таблица деталей |
| `k_length_scale` | double | `1.0` | Множитель координат, $>0$ |
| `k_mass_scale` | double | `1.0` | Множитель масс SPH-элементов, $>0$ |
| `k_monitor_nid` | long long | `-1` | NID мониторной частицы; $<0$ — не задан |
| `k_solid_elements` | string | `"skip"` | `"skip"` или `"particles"` (`*ELEMENT_SOLID`) |

### Параллелизм

| Поле | Тип | По умолчанию | Ключ | Описание |
|---|---|---|---|---|
| `num_threads` | int | `7` | нет | Число потоков OpenMP (`omp_set_num_threads` в `main`; утилиты задают своё, например `--threads`) |

### Чекпоинты и рестарт

Настройки запуска; в сам чекпоинт не записываются и при рестарте не сравниваются.
Порядок работы — в [«Рестарт расчёта из чекпоинта»](../../02-how-to/05-restart-checkpoints.md).
Из `Initial_Conditions.txt` не читаются.

| Поле | Тип | По умолчанию | Описание |
|---|---|---|---|
| `checkpoint_file` | string | `""` | Путь к файлу чекпоинта; пусто — не писать |
| `checkpoint_every_steps` | size_t | `0` | Период записи в шагах; `0` — не писать |
| `restart_file` | string | `""` | Файл, с которого продолжить расчёт; пусто — старт с нуля |
| `restart_force` | bool | `false` | Разрешить менять физические настройки при рестарте (`--force`) |

Без изменения этих полей поведение цикла прежнее. Остальные поля `SimulationConfig` перечислены в
`visit_cfg` и сохраняются в блоке `CONF` чекпоинта (единственный формат — версия 6; старые версии не читаются;
из `mats` туда входят только `G`, `sigy`, `aa`, `zz`, `gamm`, `const_JH_2`, константы Джонсона–Кука — нет); при рестарте свободно
меняются только `t_total`, `save_t`, `print_step`, `num_threads`, `energy_log_*` и `contact_stats_*`,
остальные — только с `restart_force`. Среди сверяемых — четыре поля проекции нормального контакта
(`fiber_*_contact_predictor`, `fiber_contact_projection_iters`, `fiber_contact_projection_warm`), `jc_heating`,
`jc_chi`, `xsph_mode`, `art_stress_eps` и `art_stress_n`
(`visit_cfg` в [checkpoint.cpp](../../../src/checkpoint.cpp)). Число полей записывается в чекпоинт и
сверяется при чтении («другой состав настроек» — ошибка).

### Диагностика (только в особых сборках)

Поля существуют, только если сборка выполнена с соответствующей опцией CMake
([опции](cmake-options.md)). `SPH_CONTACT_STATS` определяется только для `guo2019_shot`, поэтому поля `contact_stats_*` и `contact_hold_n` относятся к нему, а не к исполняемому `SPH`.

| Поле | Тип | По умолчанию | Условие | Описание |
|---|---|---|---|---|
| `energy_log_every` | size_t | `0` | `SPH_ENERGY_LOG` | Строка лога энергии каждые k шагов; `0` — выключен. Ключ файла — да |
| `energy_log_file` | string | `"energy_log.csv"` | `SPH_ENERGY_LOG` | Файл лога. Ключ файла — да |
| `contact_stats_every` | size_t | `0` | `SPH_CONTACT_STATS` | Строка CSV статистики контактов нитей каждые k шагов. Ключа нет (ставит `guo2019_shot --contact-stats`) |
| `contact_stats_file` | string | `"contact_stats.csv"` | `SPH_CONTACT_STATS` | Файл статистики |
| `contact_stats_summary` | string | `"contact_stats_summary.txt"` | `SPH_CONTACT_STATS` | Файл итоговой сводки |
| `contact_hold_n` | bool | `false` | `SPH_CONTACT_STATS` | Контрфакт «удержание N» (меняет физику, только диагностика); ставит `guo2019_shot --hold-n` |

### Runtime-размеры (вычисляются в `main`)

| Поле | Тип | Назначение |
|---|---|---|
| `maxn` | size_t | Ёмкость под частицы (`ntotal * 1.4`) |
| `max_interaction` | size_t | Ёмкость под пары (`100 * maxn`) |
| `moni_particle` | size_t | Индекс мониторируемой частицы (пишется в `input`, читается в `time_integration`) |
| `mats` | `Materials` | Материальные константы (см. [реестр](../materials/registry.md)) |
