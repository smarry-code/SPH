# Интегрирование по времени, шаг и энергетический баланс

Страница описывает, как проект собирает правую часть (RHS) из вкладов отдельных
взаимодействий и продвигает состояние во времени: явные решатели
(`PredictorCorrector` и `Euler`), порядок операций внутри одного шага,
критерий выбора шага $\Delta t$ (CFL и ограничители от контакта, нитей),
постобработку после шага, журнал энергетического баланса (`SPH_ENERGY_LOG`) и
сторож расчёта (`RunGuard`). Физику самих вкладов описывают соседние страницы:
[01-sph-core](01-sph-core.md) (SPH-суммы, EOS, вязкость, XSPH),
[02-constitutive](02-constitutive.md) (прочность Мизеса, JH-2, Джонсон-Кук),
[03-solid-contact](03-solid-contact.md) (контакт Solid-Solid),
[04-fiber-mechanics](04-fiber-mechanics.md) (связи нитей),
[05-fiber-contact](05-fiber-contact.md) (контакт и трение нитей). Здесь они
рассматриваются только как источники ускорений и ограничителей шага.

Источник истины — код; ссылки на файлы даны от каталога `docs/04-explanation/models/`.

## 1. Интегрируемые переменные

Состояние делится на две compute-структуры.

**Solid** (частицы, `ParticleSystem`, реальные индексы $i<n_{\rm total}$):

| Переменная | Поле | Производная (буфер `SolidRhsBuffers`) | Интегрируется |
|---|---|---|---|
| положение $\mathbf x$ | `x[d]` | $\dot{\mathbf x}=\mathbf v$ (+ XSPH-поправка) | да |
| скорость $\mathbf v$ | `vx[d]` | `dvx` | да |
| удельная внутренняя энергия $e$ (в коде `u`) | `u` | `du` | да, с отсечкой $u\ge 0$ |
| плотность $\rho$ | `rho` | `drho` | да (при `summation_density = false`, это `constexpr` в `param.hpp`) |
| девиатор напряжений $s^{\alpha\beta}$ | `sxx,syy,szz,sxy,syz,sxz` | `dsxx..dsxz` | да |
| давление $p$, скорость звука $c$ | `p`, `c` | — | нет: вычисляются из EOS на каждом RHS |
| длина сглаживания $h$ | `hsml` | — | нет: пересчитывается алгебраически (раздел 5) |
| температура $T$ (Джонсон-Кук), накопленная пластическая деформация, повреждённость $D$ | поля `jc_*` | — | нет: обновляются в постобработке (раздел 5) |

Для $d=1$ в предикторе и корректоре к $du$ и $d\rho$ добавляются члены симметрии
$-n_{\rm sym}\,p\,v/(x\rho)$ и $-n_{\rm sym}\,\rho\,v/x$; в сборке `nsym = 0`
(`param.hpp`), поэтому они нулевые. В осесимметрии ($d=2$, `axis_sim`)
осевые члены уже включены в `finalize_solid` (см. [01-sph-core](01-sph-core.md)).

**Fiber** (бусины нитей, `fiber::FiberState`, только 3D): интегрируются только
$\mathbf x$ и $\mathbf v$; у бусин нет $e,\rho,s$, поправки XSPH нет. Внутренние
переменные связей (пластическая деформация $\varepsilon_p$, ветвь Максвелла $q$,
$w_p$, $w_v$, повреждённость) — не часть решателя: они обновляются в постобработке
(раздел 5), а топология связей меняется при разрыве.

Буферы производных прошлого шага $f^{n-1}$ (`solid_rhs`, `fiber_rhs`) принадлежат
решателю и живут между итерациями цикла; они сохраняются в чекпоинте
(`checkpoint::Payload::hist`, см. [рестарт](../../02-how-to/05-restart-checkpoints.md)).

## 2. Решатели по времени

Выбор — ключом `ode_solver` (0 — `PredictorCorrector`, по умолчанию; 1 — `Euler`;
иное значение — предупреждение и откат к 0, `src/input.cpp`). Решатели не
полиморфны: цикл `detail::run_time_integration<Solver>` инстанцируется нужным типом
один раз перед циклом `while (time < t_total)` в
[time_integration.cpp](../../../src/time_integration.cpp). Оба решателя имеют
единый интерфейс `predict(...)` / `advance(...)`
([ode_solver.hpp](../../../include/ode/ode_solver.hpp)).

### 2.1. PredictorCorrector (схема проекта)

Обозначим $\mathbf y=(\mathbf v,e,\rho,s^{\alpha\beta})$, $f(\cdot)$ — производные,
$\Delta t_{n-1}$ — шаг предыдущего цикла, $\Delta t_n$ — шаг текущего (вычисляется
после RHS, раздел 4). Для шагов $n\ge2$:

1. Сохранить $\mathbf y^n$ (`v_min_`, `u_min_`, `rho_min_`, `s*_min_`) и сделать предиктор
   на полшага **по производной прошлого шага**:
$$
\mathbf y^{*}=\mathbf y^{n}+\tfrac12\,\Delta t_{n-1}\,f^{\,n-1}.
$$
   Положение $\mathbf x$ в предикторе не меняется; $e$ после предиктора отсекается
   $e^*\!\leftarrow\max(e^*,0)$ (канал энергии эту отсечку не считает — $e$
   перезаписывается из $e^n$ в корректоре).
2. Вычислить $f^{\,n}=F(\mathbf x^{n},\mathbf y^{*})$ — один раз на шаг (раздел 3);
   по предсказанному состоянию вычисляется и $\Delta t_n$.
3. Корректор — от сохранённого $\mathbf y^n$ на полный шаг:
$$
e^{n+1}=\max\!\big(e^n+\Delta t_n f_e^{\,n},0\big),\quad
\rho^{n+1}=\rho^n+\Delta t_n f_\rho^{\,n},\quad
s^{n+1}=s^n+\Delta t_n f_s^{\,n},
$$
$$
\mathbf v^{n+1}=\mathbf v^{n}+\Delta t_n\,\mathbf a^{\,n}+\tilde{\mathbf v},\qquad
\mathbf x^{n+1}=\mathbf x^{n}+\Delta t_n\,\mathbf v^{n+1},
$$
   где $\mathbf a^n$ — полное ускорение (`dvx`) и $\tilde{\mathbf v}$ — XSPH-поправка
   (`av`, см. ниже).

**Первый шаг** (`first_step() == true`): предиктора нет, $f^0=F(\mathbf x^0,\mathbf y^0)$,
все переменные продвигаются на **полшага**
$\mathbf y^{1}=\mathbf y^{0}+\tfrac12\Delta t_0 f^0$ (коэффициент $\Delta t/2$
вместо $\Delta t$, база — текущее состояние), а координата — полным шагом по уже
полушаговой скорости: $\mathbf x^1=\mathbf x^0+\Delta t_0\mathbf v^{1}$. Признак
сбрасывается в `mark_step_done()` в конце шага; восстанавливается из чекпоинта.

Интерпретация схемы (вывод из кода, в ТЗ не заявлено): после первого шага
сохранённые переменные находятся на уровне «полуцелого» времени, и каждый
следующий шаг — «кик–дрейф» (leapfrog): $\mathbf v^{n+1}=\mathbf v^{n}+\Delta t\,f(\mathbf x^n,\mathbf y^*)$,
$\mathbf x^{n+1}=\mathbf x^n+\Delta t\,\mathbf v^{n+1}$, а предиктор $\mathbf y^*$
нужен, чтобы вычислять скоростно-зависимые члены (вязкость, девиатор, давление
через $\rho$, контакт, трение) на «синхронизированном» состоянии,
экстраполированном по $f^{n-1}$. Поэтому в выводе и в энергетическом балансе
скорость и положение сдвинуты на полшага.

**XSPH** (`epsilon`, `xsph_mode`). Поправка $\tilde{\mathbf v}_i=\varepsilon\,\mathrm{av}_i$,
где $\mathrm{av}$ — сырая SPH-сумма из `SolidDelta::xsph`, умножаемая на
$\varepsilon=$`cfg.epsilon` при записи в буфер (`single_step.cpp`). Два режима
(только для PC, у Эйлера XSPH нет):

| `xsph_mode` | Скорость | Положение | Функция |
|---|---|---|---|
| `velocity` (0, по умолчанию) | $\mathbf v=\text{base}+c\,\mathbf a+\tilde{\mathbf v}$ | $\mathbf x\mathrel{+}=\Delta t\,\mathbf v$ | `pc_finalize_velocity_position_with_av` |
| `position` (1) | $\mathbf v=\text{base}+c\,\mathbf a$ | $\mathbf x\mathrel{+}=\Delta t(\mathbf v+\tilde{\mathbf v})$ | `pc_finalize_velocity_position_xsph_position` |

Здесь $c=\Delta t/2$ на первом шаге (base — текущая `vx`) и $c=\Delta t$ позже
(base — сохранённая $\mathbf v^n$). В режиме `position` кинетическая энергия не
получает вклада от XSPH (канал $D_{\rm xsph}=0$), `SolidVelocityUpdate::av = nullptr`.
Для нитей поправки нет ни в одном режиме (`pc_finalize_velocity_position`).

### 2.2. Euler (эталон для тестов сходимости)

$$
\mathbf y^{n+1}=\mathbf y^{n}+\Delta t\,f(\mathbf y^{n}),\qquad
\mathbf x^{n+1}=\mathbf x^{n}+\Delta t\,\mathbf v^{n}
$$
— положение по **старой** скорости (`euler_step_position_velocity`: сначала $\mathbf x$,
затем $\mathbf v$), XSPH не применяется (`av` считается, но не используется),
отсечка $e\ge0$ сохраняется. `predict()` — пустая операция, `first_step()` всегда
`false`. Код помечает схему как неустойчивую на незатухающих упругих связях;
для рабочих расчётов не предназначена (комментарий в `simulation_config.hpp`).

## 3. Порядок операций одного шага и сборка RHS

Один проход тела цикла в
[time_integration.cpp](../../../src/time_integration.cpp) (`detail::run_time_integration`):

| № | Операция | Где |
|---|---|---|
| 0 | вывод прогресса каждые `print_step` шагов | цикл |
| 1 | `solver.predict` (предиктор Solid + Fiber; на 1-м шаге и у Euler — пусто) | `ode_solver.hpp` |
| 2 | RHS Solid: `single_step` — EOS, ghost-частицы, поиск пар, ядра SPH и контакта Solid-Solid, reduce, `finalize_solid` | `single_step.cpp` |
| 3 | `update_smoothing_length` (`sle == 1`) — пересчёт $h$ **между RHS и $\Delta t$** (см. 5.2) | `postprocessing.hpp` |
| 4 | если есть нити: RHS Fiber-Bond, Fiber-Fiber-Contact, Fiber-Solid-Contact; пределы $\Delta t$ от них | `time_integration.cpp` |
| 5 | [1] проекция закрепления (clamp) на $\mathbf a$ и $\mathrm{av}$, реакция опоры | `clamp.hpp` |
| 6 | $\Delta t=\min(\dots)$ (раздел 4) | `time_integration.cpp` |
| 6а | при `fiber_*_contact_predictor = 1` (законы 1–3): **проекция нормального контакта** — тёплый старт, $K$ итераций Якоби, $\Delta v/h$ в ускорения, [3] повторная проекция зажимов, слоты трения из нормальных | `fiber_normal_projection.hpp` |
| 7 | стадия трения нить-нить и Fiber-Solid ([2] $c=0$ у зажатых; [3] повторная проекция) | `fiber_friction.hpp` |
| 8 | (при `SPH_ENERGY_LOG`) накопление каналов потерь за шаг | `energy_log.hpp` |
| 9 | `solver.advance` — корректор (или шаг Эйлера) Solid и Fiber | `ode_solver.hpp` |
| 10 | постобработка: `apply_material_strength`, затем `run_fiber_topology_postprocessing` | `postprocessing.hpp` |
| 11 | вывод кадра, если `time >= next_output_time` (до `time += dt`) | `output.hpp` |
| 12 | `time += dt`, запись монитора в `log.txt` (осевая компонента скорости частицы `moni_particle`) | цикл |
| 13 | строка лога энергии (каждые `energy_log_every` шагов), `RunGuard::should_stop`, чекпоинт | цикл |

### 3.1. Состав RHS

Вклады складываются в ускорения по схеме «ядро пишет вклады → reduce по шардам»
([contributions.hpp](../../../include/rhs/contributions.hpp),
[cpu_executor.hpp](../../../include/rhs/cpu_executor.hpp)): ядро формирует
`Contribution` через `ContributionWriter`, а редукция владельцем-вычислителем
(owner-computes, диапазон индексов режется на шарды) суммирует вклады каждой
частицы без мьютексов и атомиков. Порядок суммирования канонический и не зависит
от числа потоков (`parallel_for_segments`, `kDetSegments`), поэтому результат
детерминирован. Параллельны и поэлементные циклы решателя (`predict`/`correct`/`advance`, шаг
Эйлера, режим `xsph_mode = position`), `build_solid_extended`, `SpatialIndex::build` и
`Diagnostics::validate`; результат побитово прежний (последовательным остаётся накопление
`u_clamp_acc` при `SPH_ENERGY_LOG`) — см. [Параллелизм](../architecture/parallelism.md).

**Solid** (`single_step`, один `reduce_contributions_into`):

1. `compute_solid_derived` — $p,c$ из EOS (результат копируется в `sys.p`, `sys.c`: «шим совместимости»,
   нужен CFL и прочности);
2. `build_ghost_set`, `build_solid_extended`, `build_art_stress` (искусственное напряжение, см.
   [01-sph-core.md](01-sph-core.md)) — зеркальные частицы живут только в
   `SolidExtended`, в `sys` не пишутся;
3. `build_solid_sph_search` — пары (общие для SPH и контакта);
4. `compute_solid_sph_contributions` — давление, девиатор, вязкость, XSPH, $d\rho$, $de$, скорости деформации;
5. `compute_solid_contact_contributions` — контакт Solid-Solid (только $\Delta\mathbf v$);
   порядок склейки: сначала SPH, затем контакт;
6. `finalize_solid` — $\dot\rho,\dot e$ (в т.ч. осесимметричные добавки), $\dot s$ (закон Гука + производная Яуманна), запись
   $p,c,\dot\varepsilon,\omega$ в `sys` для постобработки.

Затем в `single_step` $\dot x=\mathbf v$ и `av = xsph·ε` записываются в буферы решателя, а
`Diagnostics::validate` проверяет конечность $p,\rho,c,e,h$ (только сообщение в
`stderr`; ни остановки, ни уменьшения $\Delta t$ — комментарий в `diagnostics.hpp`
о «адаптивном dt» остаётся задачей на будущее).

**Fiber** (только при `nfiber > 0`), в порядке вызова:
`compute_fiber_bond` → `build_fiber_contact_search` + `compute_fiber_contact_contributions`
(нормальная сила; заполняет слоты трения) → `compute_fiber_solid_contact`. Ускорения бусин
складываются в `fiber_rhs.dv`:
$\mathbf a^{\rm fib}=\mathbf a^{\rm bond}+\mathbf a^{\rm ff}+\mathbf a^{\rm fs}$; вклад
Fiber-Solid в частицу Solid прибавляется к `solid_rhs.dvx`. Контактные законы
типа 3 (импульсный) строятся по $\Delta t_{n-1}$ (`fiber_*_params.dt_prev = dt`
на этот момент — шаг прошлого цикла): новый шаг в момент RHS ещё неизвестен;
то же верно для Solid-Solid-контакта (`single_step` получает `dt` предыдущего цикла).

Для видов контакта с `predictor = 1` ядро силы не даёт, а заполняет нормальные слоты; сила
(точнее, импульс $P$, $N=P/h$) определяется отдельной стадией 6а **после** выбора $\Delta t$ и
до трения (она читает $h$ решателя и добавляет $\Delta v/h$ в `fiber_rhs.dv` и `solid_rhs.dvx`).
Параметр $\Delta t_{n-1}$ для такого вида не используется (см. [05-fiber-contact](05-fiber-contact.md#7-нормальный-контакт-как-проекция-скорости)).

После стадии трения (пункт 7) в `dvx` и `dv` лежит **полное ускорение шага**; именно оно
становится $f^{n-1}$ для предиктора следующего шага. Стадия трения и её
проекция на уровне скоростей описаны в [05-fiber-contact](05-fiber-contact.md);
здесь важен только порядок: она выполняется **после** вычисления $\Delta t$ и
использует формулу обновления скорости самого решателя (`fiber_velocity_update`,
`solid_velocity_update`: база, коэффициент, XSPH).

**Закрепление (clamp).** Три точки вставки (`clamp.hpp`): [1] после полной сборки
RHS и до $\Delta t$ — проекция $a\leftarrow0$ в зажатых компонентах (у Solid и
`av`), с накоплением реакции; [2] в стадии трения $c_k=c_s=0$ у полностью
зажатых (бесконечная масса); [3] повторная проекция после прибавления трения.
При пустом списке зажимов цикл побитово прежний.

## 4. Выбор шага по времени

После RHS и закрепления
$$
\Delta t=\min\!\Big(\Delta t_{\rm CFL},\;\Delta t_{\rm cont},\;\Delta t_{\rm fib}\Big),\qquad
\text{и если }s_{\Delta t}<1:\ \Delta t\leftarrow s_{\Delta t}\,\Delta t.
$$

**CFL** (`compute_cfl_timestep`):
$$
\Delta t_{\rm CFL}=\mathrm{CFL}\;\min_{i<n_{\rm total}}\frac{h_i}{c_i+|\mathbf v_i|},
$$
минимум по реальным частицам; начальное значение минимума $10^6$. Скорость звука
$c_i$ — из EOS этого RHS (предсказанное состояние), $\mathbf v_i$ — предсказанная
скорость, $h_i$ — **уже обновлённая** на шаге 3 длина сглаживания.
Член искусственной вязкости в критерии отсутствует. Если $c_i+|\mathbf v_i|=0$,
получается $+\infty$ (частица шаг не ограничивает); если $h_i=0$ и знаменатель
положителен — шаг нулевой (коллапс, ловится сторожем).

**Контакт Solid-Solid** $\Delta t_{\rm cont}$ — `wb.min_dt_contact`, минимум
`dt_limit` по межтельным парам с проникновением ($p>0$), считается каждый шаг
в `single_step` при `contact_force_type != 0` (иначе $10^{30}$, не ограничивает).
Предел зависит от типа ([contact_force.hpp](../../../include/contact_force.hpp)),
`kOscSafety = 0.8`:

- герцевская ветвь (пинбол): $0.8/\sqrt{k_{\rm eff}/m_{\rm red}}$, $k_{\rm eff}=1.5F/p$;
- вязкая ветвь: $0.8\cdot2\,\Delta t_{\rm ref}/k_p$;
- импульсный тип 3: $0.8\cdot2\,\Delta t/\gamma$ с $\Delta t$ — шагом прошлого цикла
  (при $\gamma\le1.6$ не срабатывает);
- Свиглов тип: линейная ветвь с $k_{\rm eff}=F/p$ (детали — [03-solid-contact](03-solid-contact.md)).

Здесь $\Delta t_{\rm ref}$ — локальный CFL-шаг пары: минимум по двум частицам
$\mathrm{CFL}\,h_k/(c_k+|\mathbf v_k|)$ (`solid_contact_pair_dt_ref`; $10^{30}$ при нулевом знаменателе).

**Нити** $\Delta t_{\rm fib}$ — минимум из (нити в CFL **не** входят):

- связи, `compute_fiber_bond_dt_limit`: для каждой связи
  $0.8/\omega$, $\omega=\sqrt{k_{\rm eff}/m_{\rm red}}$, $k_{\rm eff}=E_{\rm inst}A/l_0$,
  $m_{\rm red}=m_im_j/(m_i+m_j)$ ($E_{\rm inst}=E+E_1$ для неупругого материала); при
  демпфере Кельвина дополнительно
  $0.8\,h_{\rm cr}(\sqrt2\,\xi)/(\sqrt2\,\omega)$, $h_{\rm cr}(x)=2/(\sqrt{1+4x^2}+2x)$
  (старшая мода цепочки, предел Белычко). Считается по всем связям независимо от знака деформации;
  без связей — `numeric_limits::max()`;
- контакт нить-нить, `compute_fiber_contact_dt_limit`: Герц —
  $0.8/\sqrt{k_{\rm eff}/m_{\rm red}}$, $k_{\rm eff}=1.5F/\delta$ по эффективным массам точек контакта,
  только при $\delta>0$; законы 1–4 — `res.dt_limit` из слота пары (см. [05-fiber-contact](05-fiber-contact.md));
  **в режиме проекции** (`predictor = 1`) слот пары получает `kNoLimit`: нормальный контакт $\Delta t$ не
  ограничивает (жёсткости у контакта нет), шаг определяют связи, CFL и контакт другого вида;
- контакт Fiber-Solid, `compute_fiber_solid_contact_dt_limit`: та же формула, $m_{\rm red}=m_am_s/(m_a+m_s)$.

Известное ограничение (комментарий в коде): предел контакта считается по одному
контакту, а жёсткости нескольких контактов одной бусины складываются — при высокой
кратности контактов оценка может быть недостаточно консервативной.

**Масштаб** `dt_scale` $s_{\Delta t}\in(0;1]$ применяется последним и только при
$s<1$; вне диапазона — сообщение и откат к 1. Предел импульсного закона нитей
делится на $s$, чтобы замкнутая цепочка «шаг → сила → шаг» не схлопывала $\Delta t$
(`impulse_unscaled_limit`).

**Начальный шаг.** `dt` задаётся в `main.cpp`: 0.005 (shocktube) / $5\cdot10^{-5}$
(shearcavity); для импорта `.k` — по `compute_cfl_timestep` реальной геометрии, если
значение конечно и положительно. Для предиктора первого шага оно не нужно
(предиктора нет), но служит «шагом прошлого цикла» $\Delta t_{\rm prev}$ для
контактных законов типа 3 в первом RHS.

**Диагностика на экране**: при `contact_force_type != 0` каждые `print_step` шагов
выводится, чем ограничен шаг («КОНТАКТОМ» если $\Delta t_{\rm cont}\le\Delta t_{\rm CFL}$,
иначе «CFL»), число контактных пар и максимальное проникание.

**Удалённые/разрушенные частицы.** Частицы Solid из расчёта не удаляются:
ни один критерий разрушения не сокращает массивы, число $n_{\rm total}$ постоянно,
и CFL/контактный предел учитывают все частицы. Для Джонсон-Кука разрушенная
частица ($D\ge1$) остаётся в расчёте с нулевым модулем сдвига
(`effective_shear_modulus`, [02-constitutive](02-constitutive.md)). Удаляются
**связи нитей** (раздел 5.3), бусины остаются; ограничитель
$\Delta t_{\rm fib}$ строится по текущему списку `topology.bonds`, а для бусины без связей
контактный закон использует виртуальную связь длины $d$ (`compute_fiber_bead_dt_ref`).

## 5. Постобработка шага

Все дискретные пороговые правила выполняются в
[postprocessing.hpp](../../../include/interactions/postprocessing.hpp).

### 5.1. Прочность Solid — `apply_material_strength`

Вызывается **после** корректора, читает финальные $s^{\alpha\beta}$ и по типу материала:

- `itype` 3–9 — возврат на поверхность текучести Мизеса (`prochnost_mises`);
- `itype` 10–11 — JH-2 (`prochnost_cem_JH2`, получает `dt` и `drho`; читает $\dot\varepsilon,\omega$ из `sys`);
- `itype` 12–15 (`is_jc_itype`) — Джонсон-Кук (`apply_jc_strength`: пластика, температура, повреждение).

Физика — в [02-constitutive](02-constitutive.md). Под `SPH_ENERGY_LOG` Мизес
накапливает $W_{\rm plastic}=\sum w_i(\,s{:}s|_{\rm до}-s{:}s|_{\rm после})/(4G\rho)$
(информационно, уже внутри $U_{\rm int}$), для JC — пластическая работа $Y\,d\varepsilon_p$; JH-2 в канал не входит.

### 5.2. Длина сглаживания — `update_smoothing_length` (ключ `sle`)

| `sle` | Действие |
|---|---|
| 0 | $h$ постоянна |
| 1 (по умолчанию) | $h_i=1{,}3\,(m_i/\rho_i)^{1/d}$ (`h_upgrade`, коэффициент зашит в код); после — проверка на NaN/Inf, при обнаружении — исключение |
| 2 | в рабочем пути — `std::runtime_error` (нужна сумма по парам $\nabla\!\cdot\mathbf v$, RHS-конвейер её не строит); сама формула $\dot h=(h/d)\,\nabla\!\cdot\mathbf v$ в `h_upgrade.cpp` осталась |

Особенность, оговоренная в шапке файла: вызов стоит **сразу после `single_step`
и до CFL/корректора**, а не после корректора. В этот момент `sys.rho` —
предсказанная плотность, и именно от неё считается $h$ — так же, как это делал
исходный `h_upgrade` внутри `single_step`. Следствие: RHS шага $n$ использует
$h$ предыдущего шага (диагностика `Diagnostics::validate` тоже), а CFL этого же шага — свежий.
Константа `rho_floor` в `h_upgrade.cpp` объявлена, но не применяется:
при $\rho\le0$ получится NaN/Inf, и сработает исключение.

### 5.3. Нити

`run_fiber_topology_postprocessing` (при `nfiber > 0`) выполняет по порядку:

1. `apply_fiber_bond_constitutive` — обновление внутренних переменных
   неупругих связей ($\varepsilon_p,q,\varepsilon_{\rm prev},w_p,w_v,D$) по
   $\varepsilon=l/l_0-1$ от $\mathbf x^{n+1}$ с шагом $\Delta t_n$; параллельно по связям (связь $k$ пишет
   только `bonds[k].st`); для упругих материалов — пустая операция, при
   отсутствии неупругих материалов выход сразу;
2. `find_fiber_break_events` — разрыв, если выполняется первый сработавший критерий:
   $l>l_{\rm break}$ (`Strain`); $\varepsilon_p\ge\varepsilon_p^{\rm fail}$ (`PlasticStrain`);
   при включённой повреждённости — $D\ge1$ (`Damage`, критерий $\sigma_u$ тогда не действует);
   иначе $E(\varepsilon-\varepsilon_p)+q\ge\sigma_u$ (`Stress`; член Кельвина не входит).
   Поиск — только чтение, события сортируются по индексу связи;
3. `apply_fiber_break_events` — последовательно: связь удаляется из `topology.bonds`,
   у $j$-стороны (и всей цепочки `next[]`) новый `object_id` = (максимум по бусинам)+1,
   `topology.version` увеличивается один раз на проход. Кинематические структуры
   контакта и связей строятся заново на каждом RHS, отдельной инвалидации не нужно.

Порядок «constitutive → поиск → применение» зафиксирован внутри обёртки.
`init_fiber_bond_state` (один раз до цикла, не при рестарте) задаёт
$\varepsilon_{\rm prev}=l/l_0-1$, чтобы предварительно натянутая связь не давала фиктивного скачка $\Delta\varepsilon$.

Начальная температура JC-частиц ($T=T_r$) и опорный радиус $r_0=x^0$ (вес энергии в
осесимметрии) задаются перед циклом, если расчёт не рестартуется.

## 6. Журнал энергетического баланса

Только в сборке `-DSPH_ENERGY_LOG=ON`, включается ключом `energy_log_every = k > 0`
(строка CSV каждые $k$ шагов после шага интегратора; файл — `energy_log_file`,
по умолчанию `energy_log.csv`) —
[energy_log.hpp](../../../include/energy_log.hpp). Все врезки только читают состояние:
траектория побитово не меняется. Без макроса заголовок пуст.

### 6.1. Масса и вес энергии

Масса Solid-частицы в суммах энергии — `energy_mass(sys, i)`
([particle_system.hpp](../../../include/particle_system.hpp)):
$$
w_i=\begin{cases}2\pi\,r_{0,i}\,m_i & \text{2D осесимметрия (axis\_sim, dim = 2)},\\ m_i & \text{иначе}.\end{cases}
$$
$r_{0,i}=x^0_i$ — начальный радиус (`sys.r0`, задаётся `init_reference_radius`,
хранится в чекпоинте). «Плоская» масса $m=\rho_0\,\Delta x^2$ не даёт сохраняющейся
энергии в осесимметрии: $2\pi r_0 m$ — масса кольца, постоянная во времени, и столбцы Solid
дают энергию всего тела (Дж). Для нитей — $m_a$ без изменений.

### 6.2. Запасы (мгновенные)

`compute_energy_stocks`:
$$
\begin{aligned}
KE_{\rm solid}&=\sum_i \tfrac12 w_i|\mathbf v_i|^2,\qquad
KE_{\rm fiber}=\sum_a \tfrac12 m_a|\mathbf v_a|^2,\\
U_{\rm int}&=\sum_i w_i\,e_i,\qquad
U_{\rm shear}=\sum_i w_i\,\frac{s^{\alpha\beta}s^{\alpha\beta}}{4G\rho_i}\ \ (\text{информационно, уже в }U_{\rm int}),
\end{aligned}
$$
где $s^{\alpha\beta}s^{\alpha\beta}=s_{xx}^2+s_{yy}^2+s_{zz}^2+2(s_{xy}^2+s_{yz}^2+s_{xz}^2)$,
а $G$ — модуль сдвига материала (при $G\le0$ или $\rho\le0$ член нулевой).
Нити: $U_{\rm bond}=\sum_{\rm bonds}$ `fiber_bond_stored_energy`:

- упругая связь: $\tfrac12\frac{EA}{l_0}(l-l_0)^2$ при $l>l_0$ (сжатие запаса не даёт);
- неупругая: $\tfrac12E(\varepsilon-\varepsilon_p)^2Al_0$ при $\varepsilon>\varepsilon_p$ плюс
  $\tfrac12\frac{q^2}{E_1}Al_0$ при $E_1>0$ (пружина Максвелла), умноженное на $(1-D)$ при $D>0$.
  Накопленные $w_p,w_v$ в запас не входят.

Потенциал Герца $U_{\rm hertz}=\sum\tfrac25F\delta$ ($F=k\delta^{3/2}$) отдельно для
Fiber-Fiber и Fiber-Solid, считается собственным поиском контактов **только на
записываемых строках** и только для законов 0 и 4 (`fiber_contact_law_has_hertz_potential`).
Полная механическая энергия
$$
E_{\rm mech}=KE_{\rm solid}+KE_{\rm fiber}+U_{\rm int}+U_{\rm bond}+U_{\rm hertz}^{ff}+U_{\rm hertz}^{fs}.
$$

### 6.3. Каналы потерь и источников (накапливаются на каждом шаге)

Принцип: работа каждого механизма берётся **на средней скорости шага**
(`work_mid_solid`, `work_mid_fiber`). Если скорость обновляется как
$\mathbf v_{\rm new}=\mathbf v_{\rm base}+h\,\mathbf a+\tilde{\mathbf v}$ ($h=\Delta t/2$ на первом шаге, иначе $\Delta t$)
и вклад $X$ даёт приращение $\Delta\mathbf v_X$, то
$$
W_X=\sum_i w_i\,\Delta\mathbf v_{X,i}\cdot\mathbf v^{\rm mid}_i,\qquad
\mathbf v^{\rm mid}=\mathbf v_{\rm base}+\tfrac12\big(h\,\mathbf a+\tilde{\mathbf v}\big),\qquad
\sum_XW_X=KE(\mathbf v_{\rm new})-KE(\mathbf v_{\rm base})\ \text{точно}.
$$
Такое разбиение не зависит от порядка вкладов (бусина, удерживаемая трением против
натянутой связи, не получает ложной «потери» $\tfrac12mh^2a^2$). Потеря канала $= -W_X$.

| Канал (CSV) | Что накапливает | Знак |
|---|---|---|
| `D_friction` | $-W$ ускорения трения (приращение $\mathbf a$ за стадию трения, Fiber и Solid) | потеря $>0$ |
| `D_xsph` | $-W$ поправки $\tilde{\mathbf v}$ (у Euler и при `xsph_mode = position` — 0) | потеря |
| `D_contact_ss` | $-W$ ускорения контакта Solid-Solid (после проекции зажимов) | потеря |
| `D_contact_fiber` | $-W$ нормальной силы контактов нитей законов 1–3; в режиме проекции (`predictor = 1`) — $-W$ ускорения стадии проекции (приращение после проекции зажимов, т. е. работа нормальных импульсов); у закона 4 — только демпфера (отдельный проход ядра `damper_only`); у Герца — 0. Положительная работа стадии за шаг (теоретически 0) видна в мониторе Guo как `W_proj_pos_max` | потеря |
| `D_break` | $\sum$ `fiber_bond_stored_energy` порвавшихся связей по состоянию после шага | $>0$ |
| `D_break_contact` | $U_{\rm hertz}^{\rm до}-U_{\rm hertz}^{\rm после}$ удаления связей (исчезающий/появляющийся контактный потенциал) | знаковая |
| `S_u_clamp` | $\sum w_i\,(-e_i)$ при отсечке $e<0\to0$ в корректоре/Эйлере | источник $>0$ |
| `W_p_fiber` | $\sum\Delta(w_p+w_d)\,A\,l_0$ пластическая работа и работа разупрочнения связей | потеря |
| `W_v_fiber` | $\sum\Delta w_v\,A\,l_0$ вязкая диссипация (Максвелл + Кельвин) | потеря $\ge0$ |
| `clip_loss` | $\sum\sigma_m\,\Delta\varepsilon\,A\,l_0$ на отсечке $\sigma\ge0$ (быстрая разгрузка) | знаковая |
| `Q_av` | $\sum w_i\,h\,\dot e^{\rm av}_i$ — тепло искусственной вязкости | информационно, уже в $U_{\rm int}$ |
| `W_plastic` | пластическая работа Solid (Мизес, JC) | информационно, уже в $U_{\rm int}$ |

Невязка баланса, первая записанная строка задаёт $E_{\rm mech,0}$:
$$
\text{residual}=(E_{\rm mech}-E_{\rm mech,0})+D_{\rm friction}+D_{\rm xsph}+D_{\rm contact\_ss}+D_{\rm break}
-S_{u,{\rm clamp}}+W_{p}^{\rm fib}+W_{v}^{\rm fib}+D_{\rm break\_contact}-\text{clip\_loss}+D_{\rm contact\_fiber}.
$$
Для замкнутого баланса $\text{residual}\approx0$. Для нитей с вязкостью и пластикой
остаётся ошибка схемы $O(\Delta t)$ (силы по состоянию на начало шага, Кельвин по
скорости предиктора, запас связи по $\mathbf x^{n+1}$); в комментарии заголовка
на сцене E5 $\max|\text{residual}|$ падает примерно пропорционально $\Delta t$
(2.0e-4 $E_0$ при шаге лога до 5.3e-6 при $\Delta t/64$). Каналы в
чекпоинте: 12 значений (`EnergyChannels::kCount`) пишутся и восстанавливаются
вместе с $E_{\rm mech,0}$.

Столбцы CSV (порядок): `step,time,dt,KE_solid,KE_fiber,U_int,U_shear_el,U_bond,U_hertz_ff,U_hertz_fs,E_mech,D_friction,D_xsph,D_contact_ss,D_break,S_u_clamp,Q_av,W_plastic,residual,W_p_fiber,W_v_fiber,D_break_contact,clip_loss,D_contact_fiber`, затем диагностический хвост.

### 6.4. Диагностические каналы мощности по механизмам

`DiagChannels` (`p[2][8]`) накапливают за шаг **мощность** каждого механизма SPH-Solid,
энергия, полученная системой (знак «+» — рождается). Работа ускорения $X$ —
$\sum w_i\,h\,\mathbf a_{X,i}\cdot\mathbf v^{\rm mid}_i$, вклад в $e$ —
$\sum w_i\,h\,\dot e_{X,i}$. Строка 0 использует вес $w_i=$ `energy_mass`, строка 1
(суффикс `_flat`) — плоскую массу $m_i$.

| Столбец | Механизм |
|---|---|
| `P_pair_v` | парная часть $d\mathbf v$: давление + девиатор в плоской постановке, без осевых источников и без ИВ (`dv_pair = dv − dv_hoop − dv_av`) |
| `P_pair_u` | парная часть $de$ без ИВ: $\tfrac12(\sum\dot e_{\rm pair}-\dot e^{\rm av})$ |
| `P_hoop_v` | осевые (hoop) источники в $d\mathbf v$: $(s_{xx}-s_{zz})/r$ в $r$ и $s_{xy}/r$ в $z$, сглаженные ядром |
| `P_av_v`, `P_av_u` | искусственная вязкость в $d\mathbf v$ и в $de$ ($=Q_{\rm av}$) |
| `P_dev_u` | девиаторная работа $(s_{xx}\dot\varepsilon_{xx}+s_{yy}\dot\varepsilon_{yy}+s_{zz}\dot\varepsilon_{zz}+2s_{xy}\dot\varepsilon_{xy}[+2s_{yz}\dot\varepsilon_{yz}+2s_{xz}\dot\varepsilon_{xz}])/\rho$ (в осесимметрии без компонент $yz,xz$) |
| `P_pdv_u` | $-(p/\rho)\,v_r/r$ — работа давления на окружной деформации (только осесимметрия) |
| `P_szz_u` | часть `P_dev_u`: $s_{zz}\dot\varepsilon_{zz}/\rho$; в сумму `closure` **не входит** (справочно) |
| `closure` | $\text{residual}-(P_{\rm pair\_v}+P_{\rm pair\_u}+P_{\rm hoop\_v}+P_{\rm av\_v}+P_{\rm av\_u}+P_{\rm dev\_u}+P_{\rm pdv\_u})$ при весе `energy_mass` |

По тождеству $E_{\rm mech}-E_0=\sum_{\rm мех}P-(\text{потери})+S_{u,{\rm clamp}}$ величина
`closure` $\approx$ шум; её отличие от нуля означает пропущенный механизм
(контакт, зажим, проекция опоры). Консервативность пар проверяется
(тест G5): `P_pair_v_flat + P_pair_u_flat ≈ 0`. Накопители `DiagChannels`
**не** сохраняются в чекпоинте: после рестарта они начинаются с нуля и `closure`
осмыслен только по приращениям. Диагностика выполняется лишь при
`wb.energy_log` и ненулевых буферах (`dv_pair`, `du_pair` и т.д.).

## 7. Сторож расчёта (`RunGuard`)

[run_guard.hpp](../../../include/run_guard.hpp). `should_stop(dt, time, t_total)`
вызывается один раз за шаг после `time += dt` (и после записи энергетической строки);
по умолчанию сторож выключен (нулевые лимиты, `main.cpp` его не передаёт). Порядок проверок:

1. вызывается наблюдатель шага (`set_step_observer`), если задан;
2. обновляются счётчики: число шагов, первый $\Delta t$, минимальный $\Delta t$ (учитываются только $\Delta t>0$);
3. **Stalled** («застой»): оценка остатка
   $N_{\rm proj}=(t_{\rm total}-t)/\Delta t$ (при $\Delta t\le0$ — $10^{300}$, при $t\ge t_{\rm total}$ — 0).
   Если `max_projected_steps > 0` и $N_{\rm proj}>$ потолка, счётчик превышения растёт;
   при достижении `patience` подряд — остановка (счётчик сбрасывается при первом «хорошем» шаге);
4. **WallClock**: `wall_limit_s > 0` и прошло больше секунд по `steady_clock` с момента создания сторожа;
5. **StepLimit** — последним: `max_steps > 0` и `steps >= max_steps`. Это штатный останов (не авария);
   проверяется после аварийных, чтобы при их совпадении причиной осталась авария.

`failed()` — истина для WallClock и Stalled. В цикле при остановке выводится
`message()`: для WallClock сообщение различает «шаг здоров» (худший $\Delta t$ не
просел более чем в 10 раз относительно первого, `step_healthy`) и «коллапс шага»;
для Stalled — «потолок превышен». Показатель просадки $\mathrm{dt\_drop}=\Delta t_{\rm first}/\Delta t_{\min}$
считается по худшему шагу, а не по последнему. Счётчики шагов и $\Delta t_{\rm first/min}$ сохраняются
в чекпоинте (`restore_progress`); часы WallClock и счётчик застоя — нет. Чекпоинт
пишется после `should_stop` (в том числе на останове по любой причине) и каждые
`checkpoint_every_steps` шагов; при нормальном завершении — финальный, если последний шаг не записан.

Дополнительно к сторожу: `Diagnostics::validate` (раздел 3.1) проверяет конечность
полей после RHS, но лишь печатает NaN/Inf в `stderr`; `h_upgrade` при NaN/Inf в $h$ бросает исключение;
`finalize_solid` бросает исключение при $r=0$ в осесимметрии.

## 8. Параметры и ключи конфигурации

| Символ | Ключ / поле | По умолчанию | Ед. | Смысл |
|---|---|---|---|---|
| — | `ode_solver` / `cfg.ode_solver` | 0 | — | 0 PredictorCorrector, 1 Euler |
| CFL | `CFL` / `cfg.CFL` | 0,1 | — | число Куранта |
| $s_{\Delta t}$ | `dt_scale` / `cfg.dt_scale` | 1 | — | множитель итогового шага, $(0;1]$ |
| $\varepsilon$ | `epsilon` / `cfg.epsilon` | 0,3 | — | коэффициент XSPH |
| — | `xsph_mode` / `cfg.xsph_mode` | `velocity` (0) | — | `velocity` / `position` |
| — | `sle` / `cfg.sle` | 1 | — | пересчёт $h$ (0, 1; 2 — ошибка) |
| $t_{\rm total}$ | `t_total` | $3\cdot10^{-6}$ | с | время расчёта |
| — | `save_t` | $2.5\cdot10^{-8}$ | с | интервал вывода кадров |
| — | `print_step` | 50 | шагов | печать прогресса |
| — | `checkpoint_file`, `checkpoint_every_steps` | пусто, 0 | — | чекпоинты (0 — выкл.) |
| — | `restart_file`, `restart_force` | пусто, false | — | рестарт |
| — | `energy_log_every`, `energy_log_file` | 0, `energy_log.csv` | шагов | лог энергии (только `SPH_ENERGY_LOG`) |
| — | `summation_density` | `false` (`constexpr`) | — | если `true`, $\rho$ не интегрируется |
| — | `nsym` | 0 (`constexpr`) | — | симметрия для $d=1$ |

Ключи контакта, нитей и закрепления — на соответствующих страницах. Полный перечень
ключей: [simulation-config](../../03-reference/config/simulation-config.md),
[initial-conditions](../../03-reference/config/initial-conditions.md). Рестарт
и ограничения на изменение `CFL`/`dt_scale`/`ode_solver` — [05-restart-checkpoints](../../02-how-to/05-restart-checkpoints.md).

## 9. Реализация в коде

| Модель | Файл | Функция / тип |
|---|---|---|
| Главный цикл, порядок шага | [time_integration.cpp](../../../src/time_integration.cpp) | `detail::run_time_integration<Solver>`, `time_integration` (две перегрузки) |
| CFL-шаг | [time_integration.cpp](../../../src/time_integration.cpp) | `compute_cfl_timestep` |
| Предиктор-корректор | [ode_solver.hpp](../../../include/ode/ode_solver.hpp) | `PredictorCorrectorSolver::predict_solid/correct_solid/correct_fiber`, `pc_finalize_velocity_position*` |
| Явный Эйлер | [ode_solver.hpp](../../../include/ode/ode_solver.hpp) | `EulerSolver::step_solid/step_fiber`, `euler_step_position_velocity` |
| Формула обновления скорости | [ode_solver.hpp](../../../include/ode/ode_solver.hpp) | `SolidVelocityUpdate`, `FiberVelocityUpdate` |
| RHS Solid | [single_step.cpp](../../../src/single_step.cpp) | `single_step` |
| Редукция вкладов | [contributions.hpp](../../../include/rhs/contributions.hpp), [cpu_executor.hpp](../../../include/rhs/cpu_executor.hpp) | `reduce_contributions_into`, `CpuExecutor::parallel_for/parallel_reduce` |
| Предел шага контакта Solid | [solid_contact.hpp](../../../include/interactions/solid_contact.hpp) | `compute_solid_contact_diagnostics`, `solid_contact_pair_dt_ref` |
| Предел шага связей | [fiber_bond.hpp](../../../include/interactions/fiber_bond.hpp) | `compute_fiber_bond_dt_limit`, `fiber_bond_pair_dt` |
| Проекция нормального контакта (стадия 6а) | [fiber_normal_projection.hpp](../../../include/interactions/fiber_normal_projection.hpp) | `NormalWarmStart::apply/rebuild`, `compute_normal_counts`, `run_normal_projection`, `normal_projection_post_residual`, `fill_friction_slots_from_normal_*` |
| Пределы шага контактов нитей | [fiber_contact.hpp](../../../include/interactions/fiber_contact.hpp), [fiber_solid_contact.hpp](../../../include/interactions/fiber_solid_contact.hpp) | `compute_fiber_contact_dt_limit`, `compute_fiber_solid_contact_dt_limit` |
| Законы контакта (`dt_limit`) | [contact_force.hpp](../../../include/contact_force.hpp) | `pinball_force`, `impulse_force` и др. |
| Пересчёт $h$ | [postprocessing.hpp](../../../include/interactions/postprocessing.hpp), [h_upgrade.cpp](../../../src/h_upgrade.cpp) | `update_smoothing_length`, `h_upgrade` |
| Прочность | [postprocessing.hpp](../../../include/interactions/postprocessing.hpp) | `apply_material_strength` |
| Нити: constitutive и разрыв | [postprocessing.hpp](../../../include/interactions/postprocessing.hpp) | `apply_fiber_bond_constitutive`, `find_fiber_break_events`, `apply_fiber_break_events`, `run_fiber_topology_postprocessing` |
| Закрепление | [clamp.hpp](../../../include/interactions/clamp.hpp) | `apply_clamp_projection`, `apply_clamp_friction_counts` |
| Энергия | [energy_log.hpp](../../../include/energy_log.hpp) | `compute_energy_stocks`, `work_mid_solid`, `work_mid_fiber`, `EnergyLogger::write_row`, `EnergyChannels`, `DiagChannels` |
| Масса в энергии | [particle_system.hpp](../../../include/particle_system.hpp) | `energy_mass`, `init_reference_radius` |
| Проверка конечности | [diagnostics.hpp](../../../include/diagnostics.hpp) | `Diagnostics::validate` |
| Сторож | [run_guard.hpp](../../../include/run_guard.hpp) | `RunGuard::should_stop`, `message` |
