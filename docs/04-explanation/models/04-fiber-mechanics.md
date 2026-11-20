# Механика нитей (дискретная модель волокон)

Раздел описывает, как в проекте представлена нить и как вычисляются силы внутри неё: дискретизацию нити цепочкой
узлов-«бусин», сосредоточенные массы, одноосный закон связи (упругость, пластичность, скоростная зависимость,
вязкоупругость, повреждённость), критерии обрыва, условие устойчивости явной схемы, геометрию полотняного
переплетения (генератор `.fib`) и физические параметры файла нитей. Контакты нить–нить, нить–тело и трение
описаны в [05-fiber-contact.md](05-fiber-contact.md); кинематические зажимы — в [03-solid-contact.md](03-solid-contact.md#2-закрепление-clamp).

Нити существуют **только при `SPH_DIM=3`** (`fiber::check_fiber_requires_dim3`, [fiber_system.hpp](../../../include/fiber/fiber_system.hpp)).

!!! note "Чего в модели нет"
    В коде нет ни изгибной, ни крутильной жёсткости: нить — цепочка осевых связей, работающих
    **только на растяжение**. Стёжка (шов) есть, но это не жёсткость нити, а отдельные связи между бусинами
    разных нитей — см. [07-fabrics.md](07-fabrics.md#6-стежок-в-ядре). Коэффициент Пуассона `nu` читается из файла и хранится в `FiberMaterial::nu`; в силах
    связей не используется (комментарий в [fiber_bond.hpp](../../../include/interactions/fiber_bond.hpp),
    `fiber_cross_section_area`), но входит в контакты нитей: модуль Герца $E^*$ и модуль сдвига
    $G_f=E/(2(1+\nu))$ ([05-fiber-contact.md](05-fiber-contact.md)). Изгибную «жёсткость» полотна создаёт только геометрия (волнистость
    переплетения) и контакт нитей.

Обозначения: $a, b$ — узлы нити (бусины), $k$ — связь между узлами $a,b$; прочие — как в общем разделе.

## 1. Дискретизация нити

### 1.1. Данные

Нить — ломаная $\mathbf x_0,\dots,\mathbf x_{n-1}$ в ненагруженной конфигурации. Каждая вершина становится
узлом (бусиной) с состоянием `FiberState`: положение $\mathbf x_a$, скорость $\mathbf v_a$, постоянная масса $m_a$,
номер нити `object_id`, индекс материала `material_id`. Каждый сегмент ломаной становится связью
`FiberBond`: пара $(a,b)$, $b=a+1$, длина $l_0=|\mathbf x_b-\mathbf x_a|$ в ненагруженной конфигурации и порог
геометрического обрыва

$$
l_\mathrm{break}=(1+\varepsilon_\mathrm{break})\,l_0 .
$$

Связь несёт индекс материала `mat`: $-1$ — материал бусины $a$ (звено цепочки, поведение прежнее), $\ge0$ —
собственный материал. Связь с собственным материалом, у которой нет пары `next[a]==b`, `prev[b]==a`
(`is_chain_bond` ложно), — **связь вне цепочки (стежок)** между бусинами разных нитей: `prev/next` она не меняет,
при разрыве только удаляется из списка. Подробнее — [07-fabrics.md](07-fabrics.md#6-стежок-в-ядре).

Топология хранится как `prev[]/next[]` (соседи по нити, $-1$ на конце) и список `bonds`; поле `version`
растёт при обрыве. Внутренние переменные закона связи (`BondState st`) лежат **внутри** связи: при обрыве
`topology.bonds` пересобирается копированием, и состояние оставшихся связей переезжает вместе с ними.

Сборка выполняется `fiber::build_model_from_polylines` ([fiber_build.hpp](../../../include/fiber/fiber_build.hpp),
определение в [fiber_input.cpp](../../../src/fiber_input.cpp)): нулевая начальная скорость, `object_id` =
номер ломаной, нулевая или отрицательная длина сегмента, менее двух точек или неизвестный материал — ошибка.

### 1.2. Сосредоточенные массы

Площадь сечения $A=\pi d^2/4$ (`fiber_cross_section_area`). Масса узла — половина массы двух прилежащих
сегментов:

$$
m_a=\rho\,A\,\frac{l_{a-1,a}+l_{a,a+1}}{2},
$$

для концевого узла берётся один прилежащий сегмент без деления пополам (т.е. $m=\rho A l/2$ для $l$ — длины
единственного сегмента). Сумма масс узлов равна $\rho A\sum l_0$. Масса постоянна и не пересчитывается при
деформации и обрыве. Связи-стежки (`mat >= 0`) в массу не входят: масса швейной нити не учитывается.

### 1.3. Закон Ньютона для узлов

Нить интегрируется (при `ode_solver = 0`; для Euler см. [06-time-integration.md](06-time-integration.md#22-euler-эталон-для-тестов-сходимости)) той же схемой предиктор–корректор, что и Solid (`ode::PredictorCorrectorSolver`,
[ode_solver.hpp](../../../include/ode/ode_solver.hpp)), без поправки XSPH: на корректоре
$\mathbf v^{n+1}=\mathbf v^n+\Delta t\,\dot{\mathbf v}+\ldots$, $\mathbf x^{n+1}=\mathbf x^n+\Delta t\,\mathbf v^{n+1}$
(первый шаг — с $\Delta t/2$). Ускорение узла — сумма вкладов связей и контактов, делённая на $m_a$.
Вклады связей — `FiberDelta` (только `dv[dim]`), суммируются детерминированно (`rhs::reduce_contributions`).

## 2. Сила осевой связи

Для связи $k=(a,b)$ с $\mathbf x_{ab}=\mathbf x_b-\mathbf x_a$, $l=|\mathbf x_{ab}|$, деформацией
$\varepsilon=l/l_0-1$ и скоростью деформации
$\dot\varepsilon=\big((\mathbf v_b-\mathbf v_a)\cdot\mathbf x_{ab}\big)/(l\,l_0)$ (по текущим скоростям
`FiberState::vx` на момент RHS, т.е. по скорости предиктора) сила на узел $a$

$$
\mathbf F_{a\leftarrow b}=A\,\sigma\,\frac{\mathbf x_{ab}}{l},\qquad
\mathbf F_{b\leftarrow a}=-\mathbf F_{a\leftarrow b},\qquad
\dot{\mathbf v}_a \mathrel{+}= \mathbf F_{a\leftarrow b}/m_a,\quad
\dot{\mathbf v}_b \mathrel{-}= \mathbf F_{a\leftarrow b}/m_b .
$$

Сила центральная и равна по модулю и противоположна на концах: сохраняются импульс и момент импульса
(третий закон Ньютона выполнен точно). Реализация — `interactions::compute_fiber_bond_contributions`
([fiber_bond.hpp](../../../include/interactions/fiber_bond.hpp)); ядро только читает состояние, пишет
вклады.

Два режима, переключаемые не ключом, а свойствами материала (`FiberMaterial::inelastic()`):

**Упругий (по умолчанию).** Только растяжение:

$$
\sigma=E\,\varepsilon\ \ (l>l_0),\qquad \sigma=0\ \ (l\le l_0).
$$

Сила линейна по $\varepsilon$ (жёсткость $k=EA/l_0$); при сжатии нить «провисает» и силы нет. Состояние
`st` в этом режиме нулевое и не обновляется.

**Неупругий** (включена пластика, ветвь Максвелла, демпфер Кельвина через `ETA`/`XI` или повреждённость) —
$\sigma=$ `fiber::bond_stress` (раздел 3); ветка выполняется при любом $l>0$.

## 3. Определяющий закон связи

Закон одноосный; реализован чистыми функциями в [fiber_plasticity.hpp](../../../include/fiber/fiber_plasticity.hpp)
(контекст: [fiber-plasticity-design.md](../../../notes/fiber-model/fiber-plasticity-design.md)).

### 3.1. Реологическая схема

$$
\sigma=\sigma_{ep}+q+\eta_v\dot\varepsilon,\qquad
\sigma_{ep}=E\,(\varepsilon-\varepsilon_p),\qquad
\dot q=E_1\dot\varepsilon-\frac{q}{\tau},
$$

где $\varepsilon_p\ge 0$ — накопленная пластическая деформация, $q$ — напряжение ветви Максвелла
($E_1$, время релаксации $\tau$), $\eta_v$ — вязкость демпфера Кельвина. Пружина $E$ стоит последовательно с
пластическим элементом; ветвь Максвелла и демпфер Кельвина параллельны им. Любой блок выключается
нулевыми параметрами.

Напряжение, которое возвращает `bond_stress(law, st, ε, ε̇)`:

$$
\sigma=\begin{cases}
0, & \varepsilon\le\varepsilon_p,\\[2pt]
\max\!\Big[(1-D)\big(E(\varepsilon-\varepsilon_p)+q+\eta_v\dot\varepsilon\big),\,0\Big], & \varepsilon>\varepsilon_p .
\end{cases}
$$

Особые случаи: при $\eta_v=0$ слагаемое Кельвина не вычисляется (чтобы не получить $0\cdot\infty$); при $D=0$
множитель $(1-D)$ не применяется (побитовая совместимость с упругим законом); отсечка $\sigma\ge0$ означает,
что связь никогда не толкает. Демпфер Кельвина на провисшей связи не работает.

### 3.2. Условие текучести и упрочнение

Пластика включена, если `PlasticParams::enabled()`: для `Linear`/`Power` — $\sigma_y>0$ (`sy`), для
`Table` — не менее двух точек таблицы. Статический предел текучести $\sigma_s(\varepsilon_p)$
(`HardeningType`):

| Тип | Формула | Параметры |
|---|---|---|
| `Linear` | $\sigma_s=\sigma_{y0}+H\varepsilon_p$ ($H=0$ — идеальная пластичность) | `SY`, `H` |
| `Power` | $\sigma_s=\sigma_{y0}+K\varepsilon_p^{\,n}$, $0<n\le1$ | `SY`, `K`, `N` |
| `Table` | кусочно-линейная по $(\varepsilon_p^{(i)},\sigma_s^{(i)})$; за последней точкой — константа | `PLASTIC_TABLE` |

В `Power` при $n<1$ производная в $\varepsilon_p=0$ бесконечна; код возвращает $+\infty$ и решатель переходит
к делению отрезка пополам.

Скоростное упрочнение (`RateModel`) даёт динамический предел $\sigma_y(\varepsilon_p,\dot\varepsilon_p)$, где
$\dot\varepsilon_p=\Delta\varepsilon_p/\Delta t$ (`plastic_yield_stress`):

| Модель | Формула | Параметры ($a$, $b$) |
|---|---|---|
| `CowperSymonds` | $\sigma_y=\sigma_s\left[1+(\dot\varepsilon_p/C)^{1/p}\right]$ | $a=C$ [1/с], $b=p$ |
| `JohnsonCook` | $\sigma_y=\sigma_s\left[1+c\ln\max(\dot\varepsilon_p/\dot\varepsilon_0,1)\right]$ | $a=c$, $b=\dot\varepsilon_0$ |
| `Perzyna` | $\sigma_y=\sigma_s+\eta_p\,\dot\varepsilon_p^{\,m}$ (аддитивно) | $a=\eta_p$, $b=m$ |

Скорость пластической деформации в формулах ограничивается снизу нулём. Модель скорости разрешена только
вместе с пластикой (проверка при чтении файла).

### 3.3. Возврат (радиальный возврат по одной компоненте)

`bond_constitutive_update(law, st, ε_new, dt)` — обновление за шаг (вызывается **после** `advance`, по
$\mathbf x^{n+1}$, с $\Delta t$ текущего шага). При $\Delta t\le0$ состояние не меняется.

1. Пробное напряжение $\sigma_{tr}=E(\varepsilon-\varepsilon_p)$. Если $\sigma_{tr}\le\sigma_s(\varepsilon_p)$ —
   шаг упругий.
2. Иначе приращение $\Delta=\Delta\varepsilon_p$ — единственный корень $g(\Delta)=\sigma_{tr}-E\Delta-\sigma_y(\varepsilon_p+\Delta,\Delta/\Delta t)=0$
   на $[0,\sigma_{tr}/E]$ ($g(0)>0$, $g$ строго убывает):
    - линейное упрочнение без скоростной зависимости — аналитически $\Delta=(\sigma_{tr}-\sigma_s)/(E+H)$;
    - иначе метод Ньютона с защитой: если шаг выходит из скобки или производная не конечна, берётся
      середина скобки; остановка при $|g|\le10^{-14}\sigma_{tr}$ или схлопывании скобки, не более 200 итераций.
3. $\varepsilon_p\mathrel{+}=\Delta$ (накопление монотонное, $\varepsilon_p\ge0$, пластика только при растяжении).
4. Пластическая работа на единицу объёма (трапеция) $w_p\mathrel{+}=\tfrac12(\sigma_\mathrm{start}+\sigma_\mathrm{end})\Delta$, где
   $\sigma_\mathrm{end}=\sigma_{tr}-E\Delta$, $\sigma_\mathrm{start}=E(\varepsilon_\mathrm{prev}-\varepsilon_p^\mathrm{old})$,
   зажатое в $[\sigma_s,\sigma_\mathrm{end}]$. Энергия связи — $w\cdot A\,l_0$.

### 3.4. Вязкоупругие ветви

Из-за односторонности связи ветви видят не полное $\Delta\varepsilon$, а приращение натянутой части

$$
\Delta\varepsilon_\mathrm{eff}=\max(\varepsilon^{n+1},\varepsilon_p^{n+1})-\max(\varepsilon^{n},\varepsilon_p^{n}).
$$

Если связь провисла на обоих концах шага, $\Delta\varepsilon_\mathrm{eff}=0$: $q$ только релаксирует, диссипация
Кельвина не копится.

**Максвелл** ($E_1>0$, $\tau>0$): точное решение при линейном $\varepsilon(t)$ на шаге, $r=\Delta t/\tau$,

$$
q^{n+1}=e^{-r}q^n+E_1\,\Delta\varepsilon_\mathrm{eff}\,\frac{1-e^{-r}}{r},
$$

при малом $r$ — ряды (без $0/0$). Диссипация ветви $\int q^2/(E_1\tau)\,dt$ вычисляется в замкнутой форме;
к ней добавляется диссипация Кельвина $\eta_v\Delta\varepsilon_\mathrm{eff}^2/\Delta t$; сумма — $w_v$.

**Кельвин:** $\eta_v$ задаётся напрямую (`ETA`) либо долей критического демпфирования `XI`
(приоритет у `XI`):

$$
c=2\xi\sqrt{k_\mathrm{inst}\,m_\mathrm{red}},\quad k_\mathrm{inst}=\frac{(E+E_1)A}{l_0},\quad
m_\mathrm{red}=\frac{m_am_b}{m_a+m_b},\quad \eta_v=\frac{c\,l_0}{A}
$$

(`fiber_bond_visco`, `kelvin_eta_from_xi`). В сумме с `E1`: мгновенный модуль $E_\mathrm{inst}=E+E_1$
(`fiber_instant_modulus`).

### 3.5. Повреждённость с разупрочнением

Включается директивой `DAMAGE` (требует `FAILURE SIGMA_U`). Эффективное напряжение
$\sigma_\mathrm{eff}=E(\varepsilon-\varepsilon_p)+q$. Начало повреждения — первый момент, когда
$\sigma_\mathrm{eff}\ge\sigma_u$; $\varepsilon_0$ — деформация, при которой $\sigma_\mathrm{eff}$ в точности равна
$\sigma_u$ (внутри шага, при $q$ конца шага, ограничена отрезком шага). После начала пластика замораживается.
Далее $\kappa=\max\varepsilon$ и линейная огибающая

$$
\sigma_\mathrm{env}(\kappa)=\sigma_u\max\!\Big(0,\,1-\frac{\kappa-\varepsilon_0}{\Delta\varepsilon_f}\Big),\qquad
D=\min\!\Big(\max\big(D^{\,\mathrm{old}},\,1-\frac{\sigma_\mathrm{env}}{\sigma_\mathrm{eff}(\kappa)}\big),\,1\Big)
$$

($D=1$, если $\sigma_\mathrm{eff}\le 0$). $D$ не убывает; напряжение связи — $(1-D)\sigma_\mathrm{eff}$ (на
огибающей при $\varepsilon=\kappa$; разгрузка и повторная нагрузка идут по секущей к $\varepsilon_p$).
Участок разупрочнения $\Delta\varepsilon_f$ задаётся напрямую (`DEPS_F`) или через энергию разрушения
$G_f$ [Дж/м²]:

$$
\Delta\varepsilon_f=\frac{2G_f}{\sigma_u\,l_0}
$$

(`fiber_bond_damage_deps_f`; для каждой связи свой $l_0$, поэтому энергия разрыва не зависит от сетки).
Работа разупрочнения $w_d\mathrel{+}=\tfrac14E(e_0^2+e_1^2)(D-D^\mathrm{old})$, $e=\kappa-\varepsilon_p$
(среднее $Y=\tfrac12E(\kappa-\varepsilon_p)^2$ по концам шага) — в лог энергии идёт в канал потерь $W_p$ нитей.

### 3.6. Параметры закона

| Символ | Директива .fib / поле кода | По умолчанию | Единицы |
|---|---|---|---|
| $E$ | `MATERIAL … E` / `FiberMaterial::E` | — (обязательно) | Па |
| $\sigma_{y0}$ | `PLASTIC SY` / `plastic.sy` | 0 (пластики нет) | Па |
| $H$ | `PLASTIC H` / `plastic.H` | 0 | Па |
| $K$, $n$ | `PLASTIC K, N` / `plastic.K, .n` | —, 1 | Па, — |
| таблица | `PLASTIC_TABLE <id> <N>` / `plastic.table_*` | пусто | —, Па |
| $C$, $p$ | `RATE CS C= P=` / `rate.a, .b` | выкл. | 1/с, — |
| $c$, $\dot\varepsilon_0$ | `RATE JC C= EDOT0=` | выкл. | —, 1/с |
| $\eta_p$, $m$ | `RATE PERZYNA ETA= M=` | выкл. | Па·с$^m$, — |
| $\eta_v$ | `VISCO ETA=` / `visco.eta` | 0 | Па·с |
| $\xi$ | `VISCO XI=` / `visco.xi` | 0 | — |
| $E_1$, $\tau$ | `VISCO E1= TAU=` | 0, 0 | Па, с |
| $\varepsilon_p^\mathrm{max}$ | `FAILURE EPS_P=` / `failure.eps_p` | 0 (выкл.) | — |
| $\sigma_u$ | `FAILURE SIGMA_U=` / `failure.sigma_u` | 0 (выкл.) | Па |
| $\Delta\varepsilon_f$ / $G_f$ | `DAMAGE DEPS_F=` или `GF=` | выкл. | —, Дж/м² |
| $\varepsilon_\mathrm{break}$ | `MATERIAL … eps_break` | — | — |

## 4. Разрушение (обрыв связи)

Критерии проверяются после конститутивного прохода, по состоянию на $\mathbf x^{n+1}$
(`find_fiber_break_events`, [postprocessing.hpp](../../../include/interactions/postprocessing.hpp)). Порядок
проверки в коде:

1. **Геометрический (всегда)**: $l>l_\mathrm{break}=(1+\varepsilon_\mathrm{break})l_0$ — причина `Strain`.
2. Если материал имеет `FAILURE` (`has_failure()`):
    - **по пластической деформации**: $\varepsilon_p\ge\varepsilon_p^\mathrm{max}$ (`EPS_P`) — `PlasticStrain`;
    - иначе, если включена повреждённость и задан $\sigma_u$: обрыв при $D\ge1$ — `Damage` (критерий по
      $\sigma_u$ при этом **не** применяется: связь не рвётся в момент $\sigma=\sigma_u$);
    - иначе, при $\sigma_u>0$: $E(\varepsilon-\varepsilon_p)+q\ge\sigma_u$ — `Stress`. Член Кельвина
      $\eta_v\dot\varepsilon$ в критерий **не** включается (он считается численным демпфированием).

`apply_fiber_break_events` выполняется строго последовательно: события сортируются по индексу связи, на связи
$(a,b)$ обнуляются `next[a]`, `prev[b]`, а **b-сторона** (хвост нити) получает новый `object_id`
(максимум $+1$; корректно и при каскаде); связь удаляется из `topology.bonds`; `topology.version`
увеличивается один раз за проход. Для связи вне цепочки (стежок, `is_chain_bond` ложно) `next`, `prev` и
`object_id` не меняются: связь только удаляется из `topology.bonds`. Тип материала при обрыве не меняется. Оборванный сегмент пропадает из
контактного поиска нитей, а фильтр «соседи одной нити» учитывает новый `object_id`.

Обрыв — потеря запасённой энергии связи. Для лога энергии она равна
$\tfrac12E(\varepsilon-\varepsilon_p)^2Al_0\,(1-D)+\tfrac12\,q^2/E_1\,Al_0\,(1-D)$ (упругая связь —
$\tfrac12(EA/l_0)(l-l_0)^2$), `fiber_bond_stored_energy` (под `SPH_ENERGY_LOG`).

## 5. Устойчивость явной схемы

Нить даёт ограничение на шаг по времени, которое берётся минимумом с CFL/контактными ограничениями Solid
(`time_integration.cpp`: `dt = min(dt, fiber_dt_limit)`). Для пары узлов с жёсткостью и приведённой массой

$$
k_\mathrm{eff}=\frac{E_\mathrm{inst}A}{l_0},\quad m_\mathrm{red}=\frac{m_am_b}{m_a+m_b},\quad
\omega=\sqrt{k_\mathrm{eff}/m_\mathrm{red}}
$$

упругий предел

$$
\Delta t_\mathrm{lim}=\frac{0.8}{\omega}
$$

(множитель 0.8 — запас; предел самой схемы для консервативной пары — $2/\omega$). Считается **для всех**
связей независимо от знака деформации (чтобы шаг не «прыгал» при переходе через $l_0$); результат — минимум
по связям (`compute_fiber_bond_dt_limit`, `fiber_bond_pair_dt`), без связей — максимум `double`.

**Демпфер Кельвина.** При $\xi>0$ или $\eta_v>0$ ($\xi=\eta_vA/(2l_0\sqrt{k_\mathrm{eff}m_\mathrm{red}})$,
если задана `ETA`) добавляется предел по старшей моде цепочки: $\omega_\mathrm{max}=\sqrt2\,\omega$,
$\xi_\mathrm{max}=\sqrt2\,\xi$ (бусина в двух связях, оба конца с $m_\mathrm{red}$):

$$
\Delta t_\mathrm{lim}=\min\!\Big(\frac{0.8}{\omega},\ \frac{0.8\,h_\mathrm{cr}(\sqrt2\xi)}{\sqrt2\,\omega}\Big),\qquad
h_\mathrm{cr}(\xi)=2\big(\sqrt{1+4\xi^2}-2\xi\big)=\frac{2}{\sqrt{1+4\xi^2}+2\xi}.
$$

$h_\mathrm{cr}$ — критическое $h=\omega\Delta t$ схемы проекта (предиктор скорости, RHS с демпфером по скорости
предиктора, корректор): сила Кельвина берётся по запаздывающей скорости, поэтому предел отличается от формулы
Белычко (в коде указано как «Белычко с удвоенным $\xi$»; подтверждено тестом Т9 пары и цепочки,
`tests/test_fiber_plasticity_integration.cpp`). Демпфер никогда не ослабляет упругий предел. Без `VISCO`
(в т.ч. при одной пластике) значение прежнее побитово. Приведённая масса при $E_1>0$ считается по
$E_\mathrm{inst}=E+E_1$.

Тот же член используется как `dt_ref` бусины в законах контакта нитей
(`compute_fiber_bead_dt_ref`, см. [05-fiber-contact.md](05-fiber-contact.md)).

## 6. Геометрия полотняного переплетения

!!! info "Дальше по теме"
    Пакет слоёв kT, стёганое одеяло (патчи, линии шва, стежки `grid`/`seam`) и мозаичная ткань описаны в
    [07-fabrics.md](07-fabrics.md); здесь — только базовое полотно.

Генератор [plain_weave.hpp](../../../include/fiber/plain_weave.hpp) (`build_plain_weave`,
`write_plain_weave_fib`) строит `.fib` ткани полотняного переплетения; не зависит от `dim` и не входит в
расчёт — он порождает входной файл. Плоскость ткани $xy$, толщина вдоль $z$, центр $(c_x,c_y,c_z)$,
$n_x$ нитей основы (вдоль $x$) и $n_y$ нитей утка (вдоль $y$), шаг $p$, амплитуда волнистости $A$,
положения осей

$$
y_j=c_y+\Big(j-\tfrac{n_x-1}{2}\Big)p,\qquad x_i=c_x+\Big(i-\tfrac{n_y-1}{2}\Big)p .
$$

Профиль волнистости (`PlainWeaveProfile`):

- `Cos`: основа $j$: $z=c_z+A(-1)^j\cos\!\big(\pi(x-x_0)/p\big)$; уток $i$: $z=c_z-A(-1)^i\cos\!\big(\pi(y-y_0)/p\big)$;
- `Arcs`: ось собрана из полуэллипсов с центрами в вершинах $x_0+mp$, полуоси $p/2$ (вдоль нити) и $A$ (по $z$);
  на стыке $z=c_z$.

В пересечении $(i,j)$ основа находится на $z=c_z+A(-1)^{i+j}$, уток — на $c_z-A(-1)^{i+j}$.

**Дискретизация.** Бусины стоят точно в каждом пересечении; между соседними пересечениями ровно $k$
(`seg`$\ge2$) сегментов **равной длины хорды**: на полуволне (параметр $T\in[0,1]$) длина хорды $L$
подбирается стрельбой (бисекция $L$ так, чтобы расстояние от последней точки до конечной вершины тоже было
$L$; внутри — бисекция по параметру кривой). Профиль полуволны считается один раз и переиспользуется со сдвигом
и знаком. Выступ $e$ за крайнюю поперечную нить делится на $n_e=\max(1,\mathrm{round}(e/(p/k)))$ равных хорд
($e=0$ — выступа нет); по умолчанию $e=p/2$. Следствие равных хорд: все $l_0$ нити одинаковы, поэтому
массы внутренних узлов равны, а шаг по времени (раздел 5) определяется одним значением.

**Проверки.** Ошибка (кроме режима `--allow-penetration`, тогда предупреждение), если $p<d$ (соседние
параллельные нити перекрываются) или $2A<d$ (нити проникают друг в друга в пересечении); $E,\rho,\varepsilon_\mathrm{break}>0$,
$\mu\ge0$; предохранители: число бусин $\le5\cdot10^7$, `seg`$\le10^4$, хорд выступа $\le10^7$.

**Зажим полотна** (опция `clamp = X|Y|XY`, граница $c_e$): для каждого зажимаемого направления пишутся два
`CLAMP_BOX` (бусины с $|{\rm координата}-c|\ge c_e$ закрепляются), допустимый диапазон
$c_e\in[\,y_l+d/2,\ y_l+e\,]$ (с $y_l$ — координата крайней нити). Генератор сверяет, что в бокс попали ровно
те бусины, которые должны (иначе ошибка), и предупреждает, если бусина лежит в пределах $10^{-9}p$ от границы.
Физика зажима — отдельный раздел.

**Сводка** (`plain_weave_summary`): толщина $2A+d$, поверхностная плотность $\sum m/(L_xL_y)$ с
$\sum m=\rho A\sum l_0$, $L_x=(n_y-1)p$, $L_y=(n_x-1)p$ (при $n=1$ — $p$).

## 7. Ввод `.fib`: какие физические параметры читаются

Разбор — `fiber::load_fiber_file` ([fiber_input.cpp](../../../src/fiber_input.cpp)); формат подробно:
[fib-file.md](../../03-reference/file-formats/fib-file.md). Здесь — только физика.

- `MATERIAL <id> <E> <nu> <d> <rho> <eps_break> [mu_ff [mu_fs]]`: модуль Юнга, коэффициент Пуассона (в связях не
  используется, входит в контактные формулы), диаметр, плотность, критическое удлинение обрыва; необязательно коэффициенты трения нить–нить
  $\mu_{ff}$ и нить–Solid $\mu_{fs}$ (по умолчанию 0, должны быть $\ge0$) — их использование описано в
  [05-fiber-contact.md](05-fiber-contact.md).
- `FIBER <material_id> <N>` + $N\ge2$ строк `x y z` — ломаная в ненагруженной конфигурации.
- `PLASTIC`, `PLASTIC_TABLE`, `RATE`, `VISCO`, `FAILURE`, `DAMAGE` — раздел 3.6 (формат `КЛЮЧ=значение`, ключи без
  учёта регистра, повтор директивы для материала — ошибка). Проверки диапазонов: `SY>0`; `H>=0` либо пара
  `K>=0`, $0<N\le1$ (взаимоисключающие); таблица: первая точка $\varepsilon_p=0$, $\varepsilon_p$ строго
  возрастает, $\sigma_s>0$ не убывает; `RATE CS`: $C>0$, $P\ge1$; `RATE JC`: $C>0$, $\dot\varepsilon_0>0$;
  `RATE PERZYNA`: $\eta_p,m>0$; `VISCO`: `ETA` и `XI` одновременно нельзя, `E1>0` требует `TAU>0` и наоборот;
  `FAILURE`: `EPS_P>0`, `SIGMA_U>0`; `DAMAGE`: ровно один из `DEPS_F`, `GF` (оба $>0$) и обязательный
  `FAILURE SIGMA_U`; `THETA_C` (угловой критерий) зарезервирован — ошибка; `RATE` без `PLASTIC` — ошибка.
  `FAILURE EPS_P` без пластики — только предупреждение (критерий не сработает: $\varepsilon_p\equiv0$).
- `CLAMP_BOX` — закрепление (отдельный раздел).

Те же `FiberMaterial`/сборка используются при импорте LS-DYNA `.k` (`build_model_from_polylines`,
`parse_material_directive`).

## 8. Реализация в коде

| Модель | Файл | Функция / структура |
|---|---|---|
| Узлы, связи, топология, материал | [fiber_system.hpp](../../../include/fiber/fiber_system.hpp) | `FiberState`, `FiberTopology`, `FiberBond`, `FiberMaterial` |
| Объединение данных расчёта | [model.hpp](../../../include/fiber/model.hpp) | `fiber::Model` |
| Массы узлов, $l_0$, $l_\mathrm{break}$ | [fiber_input.cpp](../../../src/fiber_input.cpp) | `build_model_from_polylines` |
| Разбор `MATERIAL` | [fiber_input.cpp](../../../src/fiber_input.cpp) | `parse_material_directive` |
| Разбор неупругих директив | [fiber_input.cpp](../../../src/fiber_input.cpp) | `parse_plastic_directive`, `parse_rate_directive`, `parse_visco_directive`, `parse_failure_directive`, `parse_damage_directive`, `apply_inelastic` |
| Сила связи, вклады в ускорение | [fiber_bond.hpp](../../../include/interactions/fiber_bond.hpp) | `compute_fiber_bond_contributions`, `compute_fiber_bond` |
| Напряжение связи | [fiber_plasticity.hpp](../../../include/fiber/fiber_plasticity.hpp) | `bond_stress` |
| Возврат, Максвелл, повреждённость | [fiber_plasticity.hpp](../../../include/fiber/fiber_plasticity.hpp) | `bond_constitutive_update`, `plastic_yield_stress`, `static_yield_stress` |
| Кельвин через $\xi$ | [fiber_plasticity.hpp](../../../include/fiber/fiber_plasticity.hpp), [fiber_bond.hpp](../../../include/interactions/fiber_bond.hpp) | `kelvin_eta_from_xi`, `fiber_bond_visco` |
| Проход обновления состояния связей | [postprocessing.hpp](../../../include/interactions/postprocessing.hpp) | `apply_fiber_bond_constitutive` |
| Критерии и применение обрыва | [postprocessing.hpp](../../../include/interactions/postprocessing.hpp) | `find_fiber_break_events`, `apply_fiber_break_events`, `run_fiber_topology_postprocessing` |
| Запасённая энергия связи | [postprocessing.hpp](../../../include/interactions/postprocessing.hpp) | `fiber_bond_stored_energy` |
| Критический шаг | [fiber_bond.hpp](../../../include/interactions/fiber_bond.hpp) | `compute_fiber_bond_dt_limit`, `fiber_bond_pair_dt`, `fiber_kelvin_h_crit` |
| Интегрирование нитей | [ode_solver.hpp](../../../include/ode/ode_solver.hpp) | `PredictorCorrectorSolver::predict_fiber`, `correct_fiber` |
| Генерация полотняного переплетения | [plain_weave.hpp](../../../include/fiber/plain_weave.hpp) | `build_plain_weave`, `write_plain_weave_fib` |
| Одеяло, стежки, мозаика, стежок в ядре | см. [07-fabrics.md](07-fabrics.md#7-реализация-в-коде) | `build_plain_weave_quilt`, `build_plain_weave_random_cut`, `bond_material_id`, `is_chain_bond` |
