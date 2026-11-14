# Ядро SPH-расчёта сплошной среды

Глава описывает ту часть расчёта тела («Solid»), которая не зависит от материала:
функции ядра сглаживания, парную SPH-аппроксимацию уравнений неразрывности,
движения и энергии, скорости деформации и поворота, искусственную вязкость,
поправку XSPH, обновление длины сглаживания и осесимметричную постановку.
Давление $p$, скорость звука $c$ и девиатор напряжений $s^{\alpha\beta}$
берутся из определяющих соотношений (уравнение состояния, закон Гука с
производной Яуманна, пластичность, повреждение); их формулы приведены в
[02-constitutive.md](02-constitutive.md). Здесь они рассматриваются как входные
поля частиц. Контакт между телами и нити описаны в других главах раздела.

Источник истины — код ветки `fiber_new`: [solid_sph.hpp](../../../include/interactions/solid_sph.hpp),
[kernel.hpp](../../../include/kernel.hpp), [h_upgrade.cpp](../../../src/h_upgrade.cpp),
[ode_solver.hpp](../../../include/ode/ode_solver.hpp). Файлы
[int_force.cpp](../../../src/int_force.cpp) и [int_force.hpp](../../../include/int_force.hpp)
после этапа Ф4-B2 содержат только диспетчер УРС `urs()`; прежние циклы
`body1/body2` удалены, а их математика перенесена в `solid_sph.hpp` построчно.

## 1. Общая схема стадий

Один вызов RHS для тела состоит из стадий (`single_step`, комментарий в шапке
[solid_sph.hpp](../../../include/interactions/solid_sph.hpp)):

1. `compute_solid_derived` — по $\rho_i$, $u_i$ и состоянию материала считаются
   $p_i$, $c_i$ (УРС, см. 02);
2. `build_ghost_set`, `build_solid_extended` — зеркальные частицы у оси (п. 8);
3. `build_solid_sph_search` — пары соседей и примитивы
   $(i,j,W_{ij},\nabla_iW_{ij})$;
4. `compute_solid_sph_contributions` — по каждому примитиву вклады
   (`SolidDelta`) в частицу $i$ и, если $j$ реальна, в частицу $j$;
5. `rhs::reduce_contributions` — детерминированная (не зависящая от числа потоков)
   свёртка вкладов по частицам;
6. `finalize_solid` — осевые члены, девиаторная работа, множитель $1/2$ в энергии,
   закон Гука с производной Яуманна.

Вклады формируются только из read-only состояния; прямой записи в общий RHS из
парного цикла нет. Пара $(i,j)$ обрабатывается один раз (канонический порядок
$i<j$), оба вклада вычисляются внутри одного примитива.

## 2. Функции ядра

Ядро $W(r,h)=h^{-d}f(q)$, $q=r/h$, где в качестве $h$ передаётся средняя длина
сглаживания пары $h_{ij}=(h_i+h_j)/2$ (`mhsml` в вызове). Градиент ядра по
координате $\mathbf x_i$:

$$
\nabla_i W_{ij}=\frac{1}{h}\frac{\mathrm dW}{\mathrm dq}\frac{\mathbf x_{ij}}{r_{ij}}.
$$

Градиент антисимметричен по паре ($\nabla_jW_{ji}=-\nabla_iW_{ij}$), поскольку
$h_{ij}$ симметрично, а $\mathbf x_{ij}$ меняет знак. Ядро выбирается ключом
`skf` (константа в `SimulationConfig`, по умолчанию 1). Нормировочные множители
зависят от `dim` (константа сборки `SPH_DIM`, 1, 2 или 3); при другом `dim`
кубическое и квинтическое ядра бросают исключение.

### 2.1. Кубический сплайн (`skf = 1`, `kernel_cubic`, Monaghan 1985)

Носитель $2h$:

$$
W=\sigma_d\begin{cases}\dfrac23-q^2+\dfrac{q^3}{2}, & 0\le q\le1,\\[2mm]
\dfrac16(2-q)^3, & 1<q\le2,\\[1mm] 0, & q>2,\end{cases}
\qquad
\sigma_d=\frac1h,\ \frac{15}{7\pi h^2},\ \frac{3}{2\pi h^3}\ \ (d=1,2,3).
$$

Градиент в коде: при $q\le1$ $\nabla W=\sigma_d(-2+\tfrac32 q)\,\mathbf x_{ij}/h^2$
(в нуле нет деления на $r$); при $1<q\le2$ $\nabla W=-\tfrac12\sigma_d(2-q)^2\,\mathbf x_{ij}/(r\,h)$.
Формулы согласованы с аналитической производной.

### 2.2. Гауссово ядро (`skf = 2`, `kernel_gauss`)

$$
W=\frac{e^{-q^2}}{\pi^{d/2}h^d},\qquad \nabla_iW=-\frac{2\,\mathbf x_{ij}}{h^2}\,W,
$$

обрезка при $q>3$ (при $q\le3$ значение не нулевое, разрыв на границе носителя
порядка $e^{-9}$). Нормировка записана одной формулой для любого `dim`.

### 2.3. Квинтическое ядро (`skf = 3`, `kernel_quintic`, Morris 1997)

Носитель $3h$:

$$
W=\sigma_d\Big[(3-q)^5-6(2-q)^5+15(1-q)^5\Big],\qquad
\sigma_d=\frac{1}{120h},\ \frac{7}{478\pi h^2},\ \frac{1}{120\pi h^3}.
$$

Слагаемое с $(2-q)^5$ действует при $q\le2$, с $(1-q)^5$ — при $q\le1$. Для
$1<q\le3$ градиент записан как $\sigma_d\,[-5(3-q)^4+30(2-q)^4]\,\mathbf x_{ij}/(rh)$
(при $q>2$ только первый член), с защитой $r>10^{-8}$.

При $0\le q\le1$ сумма трёх слагаемых даёт полином
$\mathrm dW/\mathrm dq=\sigma_d(-120q+120q^2-50q^3)$, поэтому код записывает градиент без
деления на $r$: $\nabla W=\sigma_d(-120+120q-50q^2)\,\mathbf x_{ij}/h^2$ (то же самое, что
$\frac{\mathrm dW}{\mathrm dq}\frac{\mathbf x_{ij}}{rh}$, но без особенности при $r\to0$).

### 2.4. Радиус поиска

Пара принимается, если $r_{ij}<k_s h_{ij}$ с $k_s=2$ при `skf = 1` и $k_s=3$ при
`skf = 2, 3` (`scale_k` в `build_solid_sph_search`); критерий применяется к
$h_{ij}$ каждой пары, поэтому корректен для неоднородного поля $h$. Ячейки
пространственного индекса строятся по радиусу $k_s\max_i h_i$. Для неизвестного
`skf` функция `kernel()` обнуляет $W$ и $\nabla W$.

## 3. SPH-аппроксимация

Значение и градиент поля $f$ в частице $i$:

$$
\langle f\rangle_i=\sum_j\frac{m_j}{\rho_j}f_jW_{ij},\qquad
\langle\nabla f\rangle_i=\sum_j\frac{m_j}{\rho_j}f_j\nabla_iW_{ij}.
$$

Все дифференциальные операторы в коде записаны в форме разностей
$(f_j-f_i)$ или симметричных пар, без нормализации ядра по Шепарду и без
коррекции градиента (Либерски/ядерной коррекции в коде нет). Следовательно,
градиент скорости у границ тела и при нерегулярном расположении частиц
аппроксимируется с ошибкой, не равной нулю даже для линейного поля.
Нормализация плотности (`nor_density`) и суммирование плотности
(`summation_density`) — константы `param.hpp`, обе равны `false`, соответствующего
кода в новом пути нет.

Обозначения для пары: $\mathbf v_{ij}=\mathbf v_i-\mathbf v_j$,
$\mathbf v_{ji}=-\mathbf v_{ij}$, $\nabla W\equiv\nabla_iW_{ij}$.

## 4. Уравнение неразрывности

Плотность не пересчитывается суммированием, а интегрируется по времени из
$\mathrm d\rho/\mathrm dt$ (`summation_density = false`, constexpr). Ключ
`density_method` выбирает форму парной суммы.

| `density_method` | $\dot\rho_i$ (сырая парная сумма) | Комментарий |
|---|---|---|
| 1 (по умолчанию) | $\displaystyle\sum_j m_j\,\mathbf v_{ij}\cdot\nabla_iW_{ij}$ | классическая форма |
| 2 | $\displaystyle\rho_i\sum_j\dfrac{m_j}{\rho_j}\,\mathbf v_{ij}\cdot\nabla_iW_{ij}$ | форма через объёмы $m_j/\rho_j$ |
| 3 (любое значение, отличное от 1 и 2) | $\displaystyle\sum_j m_j\,\mathbf v_{ij}\cdot\nabla_iW_{ij}-\rho_i D_i$ | с искусственной диффузией плотности |

Для частицы $j$ в той же паре вклад получается заменой $i\leftrightarrow j$;
скалярное произведение $\mathbf v_{ji}\cdot\nabla_jW_{ji}=\mathbf v_{ij}\cdot\nabla_iW_{ij}$
не меняет знака, поэтому в `dj.drho` стоит тот же `vcc`.

Диффузионный член метода 3 (форма Molteni–Colagrossi; в комментарии кода также
названа «Bonet–Lok/Molteni–Colagrossi»):

$$
D_i=\xi\,h_{ij}\,\bar c\,\psi_{ij}\,V_j,\quad
\psi_{ij}=\frac{2\,(\rho_j/\rho_i-1)}{r_{ij}^2+\varepsilon_hh_{ij}^2}\,\mathbf x_{ij}\cdot\nabla_iW_{ij},\quad
V_j=\frac{m_j}{\rho_j},\quad \bar c=\frac{c_i+c_j}{2}.
$$

Для $j$: $D_j=\xi h_{ij}\bar c\,\psi_{ji}V_i$, $\psi_{ji}=2(\rho_i/\rho_j-1)(-\mathbf x_{ij})\cdot\nabla_iW_{ij}/(\dots)$.

Осесимметричная добавка вносится в `finalize_solid` один раз на частицу
(п. 8): $\dot\rho_i\mathrel{-}=\rho_iv_{r,i}/r_i$.

Параметры:

| Символ | Ключ / поле | По умолчанию | Единицы |
|---|---|---|---|
| выбор формы | `density_method` (читается из `Initial_Conditions.txt`) | 1 | — |
| $\xi$ | `cfg.xi` | 0.001 | — |
| $\varepsilon_h$ | `cfg.epsilon_h` | 0.01 | — |

`xi` и `epsilon_h` есть в `SimulationConfig` и чекпоинте, но
[input.cpp](../../../src/input.cpp) их не читает: из файла условий они не
задаются, остаются умолчаниями.

## 5. Скорость деформации и вихрь

Градиент скорости в частице $i$ без коррекции (`compute_solid_sph_contributions`,
блок «скорости деформаций/поворота»):

$$
\dot\varepsilon_i^{\alpha\beta}=\frac12\sum_j\frac{m_j}{\rho_j}
\Big(v_{ji}^\alpha\,\partial^\beta W_{ij}+v_{ji}^\beta\,\partial^\alpha W_{ij}\Big),\qquad
\omega_i^{\alpha\beta}=\frac12\sum_j\frac{m_j}{\rho_j}
\Big(v_{ji}^\alpha\,\partial^\beta W_{ij}-v_{ji}^\beta\,\partial^\alpha W_{ij}\Big).
$$

Здесь $\partial^\beta W_{ij}$ — компонента $\nabla_iW_{ij}$. Тензор $\omega$ —
антисимметричная часть градиента $\partial^\beta v^\alpha$; например
$\omega^{xy}=\tfrac12(\partial_yv_x-\partial_xv_y)$. Для $j$ используются те же
$\mathbf v_{ji}$, $\nabla_iW_{ij}$ с заменой масс и плотностей ($m_i/\rho_i$), что
эквивалентно $\mathbf v_{ij}\cdot\nabla_jW_{ji}$ с двумя сменами знака. Хранятся
независимые компоненты: `dexx,deyy,dezz,dexy,dexz,deyz` и `drxy,drxz,dryz`;
число активных компонент определяется `dim` (1D — только `dexx`; 2D — `dexx`,
`deyy`, `dexy`, `drxy`; 3D — все). Компонента `dezz` при `dim = 2` из парной суммы
не приходит (п. 8). Скорость деформации и вихрь записываются в `sys.exx…ryz` как
«шим совместимости» для JH-2 ([ppochnost.cpp](../../../src/ppochnost.cpp)).

## 6. Уравнение движения и энергии

Уравнение движения в форме тензора напряжений $\sigma^{\alpha\beta}=-p\delta^{\alpha\beta}+s^{\alpha\beta}$:

$$
\frac{\mathrm dv^\alpha}{\mathrm dt}=\frac1\rho\,\partial_\beta\sigma^{\alpha\beta}.
$$

Давление и девиатор берутся из определяющих соотношений (02-constitutive.md);
в силовом блоке они используются только как данные. Ключ `pa_sph` выбирает форму
симметризации (оба варианта антисимметричны по паре, поэтому суммарный импульс
сохраняется до машинной точности; силы пары равны по модулю и противоположны,
но для девиаторных членов не направлены строго вдоль линии $\mathbf x_{ij}$, так что
сохранение момента импульса не гарантируется):

| `pa_sph` | ускорение частицы $i$ от пары |
|---|---|
| 2 (по умолчанию; также любое значение, отличное от 1) | $\displaystyle \frac{\mathrm dv_i^\alpha}{\mathrm dt}=\sum_j m_j\Big[-\Big(\frac{p_i}{\rho_i^2}+\frac{p_j}{\rho_j^2}\Big)\delta^{\alpha\beta}+\frac{s_i^{\alpha\beta}}{\rho_i^2}+\frac{s_j^{\alpha\beta}}{\rho_j^2}\Big]\partial^\beta W_{ij}$ |
| 1 | $\displaystyle \frac{\mathrm dv_i^\alpha}{\mathrm dt}=\sum_j\frac{m_j}{\rho_i\rho_j}\Big[-(p_i+p_j)\delta^{\alpha\beta}+s_i^{\alpha\beta}+s_j^{\alpha\beta}\Big]\partial^\beta W_{ij}$ |

Для частицы $j$ вклад равен $-m_i\,h^{\alpha}$ с тем же $h^\alpha$ (то есть
противоположен по знаку при той же «силе пары»), поэтому $\sum_i m_i\dot{\mathbf v}_i=0$
для парных членов. Это относится и к осесимметричному случаю, но не к hoop-источникам
(п. 8) и не к пустым вкладам ghost-частиц.

Энергия. Парная часть работы давления и вязкости накапливается в `du` сырой
суммой, затем в `finalize_solid` умножается на $1/2$ (так же, как делал прежний
`int_force.cpp`):

$$
\dot e_i=\frac12\sum_j m_j\Big[\Big(\frac{p_i}{\rho_i^2}+\frac{p_j}{\rho_j^2}\Big)\mathbf v_{ij}\cdot\nabla_iW_{ij}+\Pi_{ij}\,\mathbf v_{ij}\cdot\nabla_iW_{ij}\Big]
+\frac{s_i^{\alpha\beta}\dot\varepsilon_i^{\alpha\beta}}{\rho_i}.
$$

Для `pa_sph = 1` вместо скобки давлений стоит $(p_i+p_j)/(\rho_i\rho_j)$.
Парный вклад давления записывается как $\sum_d v_{ji}^d\,h_p^d$ с $h_p^d=-(\dots)\partial^dW$
(только давление, без напряжений). Для $j$ получается то же значение (знаки
$\mathbf v_{ji}$ и $\nabla_j W$ меняются одновременно). Девиаторная работа
добавляется отдельно (`finalize_solid`):
$s:\dot\varepsilon/\rho=(s^{xx}\dot\varepsilon^{xx}+s^{yy}\dot\varepsilon^{yy}+s^{zz}\dot\varepsilon^{zz}+2s^{xy}\dot\varepsilon^{xy}+2s^{yz}\dot\varepsilon^{yz}+2s^{xz}\dot\varepsilon^{xz})/\rho$
(в 2D-осесимметрии — только первые четыре слагаемых: $s^{xx}\dot\varepsilon^{xx}+s^{yy}\dot\varepsilon^{yy}+s^{zz}\dot\varepsilon^{zz}+2s^{xy}\dot\varepsilon^{xy}$).
Напряжения в парной форме работы не участвуют.

Физическая вязкость ($\mu\nabla^2\mathbf v$, константа `visc = false`) в новом пути не
реализована.

### 6.1. Эволюция девиатора (граница с 02)

В `finalize_solid` для $\text{itype}>2$ вычисляется
$\dot s^{\alpha\beta}=2G\big(\dot\varepsilon^{\alpha\beta}-\tfrac13\mathrm{tr}\,\dot\varepsilon\,\delta^{\alpha\beta}\big)+\text{(производная Яуманна через }\omega)$,
причём след скорости деформации берётся как $-\dot\rho/(3\rho)$ (в коде
$\dot\varepsilon^{\alpha\alpha}+\dot\rho/(3\rho)$, где $\dot\rho$ уже включает осевую добавку).
Здесь $G$ — эффективный модуль сдвига
(`effective_shear_modulus`; у разрушенной частицы JC $G=0$). Подробности и
пластическая коррекция — в 02-constitutive.md. Для $\text{itype}\le2$ (газ, вода)
$\dot s=0$. Для `dim = 1` Яуманна нет, $\dot s^{yy}=\dot s^{zz}=2G\dot\rho/(3\rho)$.
Поворотные члены в 2D содержат только $\omega^{xy}$.

## 7. Искусственная вязкость Монагана

Применяется только к сближающимся парам ($\mathbf v_{ij}\cdot\mathbf x_{ij}<0$):

$$
\mu_{ij}=\frac{h_{ij}\,\mathbf v_{ij}\cdot\mathbf x_{ij}}{r_{ij}^2+\eta^2h_{ij}^2},\qquad
\Pi_{ij}=\frac{\beta\mu_{ij}^2-\alpha\bar c\,\mu_{ij}}{\bar\rho},\qquad
\bar c=\frac{c_i+c_j}2,\ \ \bar\rho=\frac{\rho_i+\rho_j}2.
$$

Ускорение и тепло (по частице $i$; для $j$ симметрично):

$$
\frac{\mathrm d\mathbf v_i}{\mathrm dt}\mathrel{+}=-\sum_jm_j\Pi_{ij}\nabla_iW_{ij},\qquad
\dot e_i\mathrel{+}=\frac12\sum_jm_j\Pi_{ij}\,\mathbf v_{ij}\cdot\nabla_iW_{ij}
$$

(множитель $1/2$ даёт общий `0.5*du` в `finalize_solid`). При $\mathbf v_{ij}\cdot\mathbf x_{ij}\ge0$
вклад нулевой (нет порога по $c$ или лимитера Бальсары). Вклад антисимметричен
по паре, импульс сохраняется; вся работа $\Pi$ превращается в тепло в обеих
частицах, поэтому $\dot e_{AV}\ge0$ для сближения.

| Символ | Поле | По умолчанию | Единицы |
|---|---|---|---|
| $\alpha$ | `cfg.alpha` (`alpha`) | 1.0 | — |
| $\beta$ | `cfg.beta` (`beta`) | 1.0 | — |
| $\eta$ | `cfg.etq` (`etq`) | 0.01 | — (в коде $\eta^2h_{ij}^2$, не $\eta h$) |

Условный параметр $\varepsilon$ в формуле знаменателя Монагана обозначен в коде
`etq`. Значения читаются из `Initial_Conditions.txt`. Искусственное напряжение
(artificial stress) в коде отсутствует; вязкость включена константой
`visc_artificial = true` (на выбор кода в новом пути не влияет), искусственный нагрев
`heat_artificial = false` не реализован.

При сборке с `SPH_ENERGY_LOG` вклад вязкости в $\dot v$ и $\dot e$ накапливается
отдельно (`dv_av`, `du_av`) только для диагностики; на расчёт не влияет.

## 8. XSPH

Средняя скорость соседей:

$$
\tilde{\mathbf v}_i=\varepsilon\sum_j\frac{2m_j}{\rho_i+\rho_j}\,(\mathbf v_j-\mathbf v_i)\,W_{ij}.
$$

В парном цикле накапливается сырая сумма (`xsph`), а `single_step` домножает её
на `cfg.epsilon` (по умолчанию $\varepsilon=0.3$). Вклад получают только реальные
частицы; ghost-частицы поправку не набирают. Ключ `xsph_mode`
(`velocity` или `position`; другое значение — ошибка разбора) задаёт, куда
входит $\tilde{\mathbf v}$. Предиктор-корректор
([ode_solver.hpp](../../../include/ode/ode_solver.hpp): `pc_finalize_velocity_position_with_av`,
`pc_finalize_velocity_position_xsph_position`):

| `xsph_mode` | скорость | положение | энергетика |
|---|---|---|---|
| `velocity` (0, по умолчанию) | $\mathbf v^{n+1}=\mathbf v_{\rm base}+c_t\,\dot{\mathbf v}+\tilde{\mathbf v}$ | $\mathbf x^{n+1}=\mathbf x^n+\Delta t\,\mathbf v^{n+1}$ | поправка меняет кинетическую энергию (диссипация $D_{xsph}>0$ в логе энергии) |
| `position` (1) | $\mathbf v^{n+1}=\mathbf v_{\rm base}+c_t\,\dot{\mathbf v}$ | $\mathbf x^{n+1}=\mathbf x^n+\Delta t\,(\mathbf v^{n+1}+\tilde{\mathbf v})$ | скорость и кинетическая энергия не затрагиваются, $D_{xsph}=0$ |

Здесь $c_t=\Delta t/2$ на первом шаге (база — текущая скорость) и $c_t=\Delta t$ на
следующих (база $\mathbf v^n$). Режим `position` — только для PC-решателя; у
явного Эйлера XSPH нет вовсе. Для трения Solid–Fiber и лога энергии при
`position` поправка не видна (`av = nullptr`). При $\varepsilon=0$ режимы совпадают.

## 9. Длина сглаживания

`h_upgrade` ([h_upgrade.cpp](../../../src/h_upgrade.cpp)) вызывается через
`interactions::update_smoothing_length` сразу после `single_step` (до CFL):

| `sle` | действие |
|---|---|
| 0 | $h$ не меняется |
| 1 (по умолчанию) | $h_i=1.3\,(m_i/\rho_i)^{1/d}$ (коэффициент 1.3 зашит; $m_i/\rho_i$ — объём в 3D, «площадь» в плоской/осесимметричной 2D, длина в 1D) |
| 2 | в `h_upgrade` предусмотрено $\dot h_i=\tfrac{h_i}{d}\,\sum_j\tfrac{m_j}{\rho_j}\mathbf v_{ij}\cdot\nabla W$, но `update_smoothing_length` для `sle = 2` бросает исключение (суммы по парам в новом конвейере нет) |

После обновления все $h_i$ проверяются на NaN/Inf (ошибка при обнаружении).
Для ghost-частиц $h$ наследуется от источника на каждом шаге.

!!! warning "Наблюдения по `sle`"
    1. Ветка `sle = 2` в `h_upgrade.cpp` недостижима из конвейера; кроме того, сумма
       $\sum_j (m_j/\rho_j)\mathbf v_{ij}\cdot\nabla_iW_{ij}$ равна **$-\nabla\cdot\mathbf v$**,
       поэтому $\dot h=(h/d)\cdot(\text{сумма})$ имеет знак, противоположный
       физическому ($\dot h=(h/d)\nabla\cdot\mathbf v$ при $h\propto\rho^{-1/d}$);
       комментарий в коде называет множитель $h/\rho$, а код — $h/d$.
    2. Константа `rho_floor` в ветке `sle = 1` объявлена, но не используется.
    3. [initial-conditions.md](../../03-reference/config/initial-conditions.md)
       (строка `sle`) описывает `sle = 2` как рабочий режим; фактически он даёт ошибку.

## 10. Постановка в зависимости от размерности

`dim` — константа сборки (`SPH_DIM`). Выбор ветки определяют `dim` и константа
`axis_sim` ([param.hpp](../../../include/param.hpp), по умолчанию `true`).

* **3D** — полный тензор $\dot\varepsilon$ и $\omega$ (шесть и три компоненты),
  девиаторная работа с шестью компонентами, ядро нормировано на $\pi h^3$.
* **1D** — компонента `exx`; ядро $1/h$; сферическая/цилиндрическая симметрия
  (`nsym`) в коде 1D обрабатывается добавками в предикторе/корректоре
  (`temp_u`, `temp_rho`), не в `solid_sph.hpp`.
* **2D** — плоская постановка ($axis\_sim=false$) в коде сохранена как ветка
  `else` у девиаторной работы (учтены `syz`, `sxz`), но константа
  `axis_sim = true` делает **любой 2D-расчёт осесимметричным**; плоской 2D-сборки
  без правки `param.hpp` нет.

### 10.1. Осесимметричная постановка ($d=2$, `axis_sim`)

Координата 0 — радиус $r=x^0$, координата 1 — ось $z=x^1$. Соглашение о
хранении (закреплено тестом `test_axisym.cpp`):

$$
s^{xx}=S_{rr},\quad s^{yy}=S_{zz},\quad s^{zz}=S_{\theta\theta},\quad s^{xy}=S_{rz},\qquad
\dot\varepsilon^{xx}=\dot\varepsilon_{rr},\ \dot\varepsilon^{yy}=\dot\varepsilon_{zz},\ \dot\varepsilon^{zz}=\dot\varepsilon_{\theta\theta}.
$$

Осевая компонента лежит в `syy`, окружная — в `szz`.

**Масса.** Частицам приписывается «плоская» масса $m=\rho_0\Delta x^2$ (масса на
единицу угла по $\theta$ и без множителя $2\pi r$); ядро — двумерное; $\rho$ в
SPH-суммах — «плоская» плотность $m/(\Delta x^2)$. Вес $2\pi r$ в SPH-суммах
не входит; он используется лишь для учёта энергии: масса кольца
$m_{\rm ring}=2\pi r_0m$ с начальным радиусом $r_0$ (`ParticleSystem::r0`,
`init_reference_radius`, не меняется при перемещении; функция `energy_mass`) —
для лога энергии, чтобы $\sum m_{\rm ring}(v^2/2+u)$ было физической энергией тела.

**Окружные члены.** Добавки к парным суммам (в обеих ветках `pa_sph`, в
`compute_solid_sph_contributions`):

$$
\frac{\mathrm dv_{r,i}}{\mathrm dt}\mathrel{+}=\sum_j\frac{m_j}{\rho_i\rho_j}\frac{s^{xx}_j-s^{zz}_j}{r_j}W_{ij},\qquad
\frac{\mathrm dv_{z,i}}{\mathrm dt}\mathrel{+}=\sum_j\frac{m_j}{\rho_i\rho_j}\frac{s^{xy}_j}{r_j}W_{ij},
$$

то есть $(S_{rr}-S_{\theta\theta})/r$ и $S_{rz}/r$ из уравнения движения в
цилиндрических координатах. Давление в радиальный окружной источник не
входит ($-p-(-p)=0$). В `finalize_solid` добавляются:

$$
\dot\varepsilon^{zz}_i=\frac{v_{r,i}}{r_i}\ \ (\text{частичная формула, не парная сумма}),\qquad
\dot\rho_i\mathrel{-}=\rho_i\frac{v_{r,i}}{r_i},\qquad
\dot e_i\mathrel{-}=\frac{p_i}{\rho_i}\frac{v_{r,i}}{r_i}.
$$

Последний член — работа давления на окружной деформации
($\dot e=-(p/\rho)\,\mathrm{div}\,\mathbf v$, $\mathrm{div}\,\mathbf v=$ плоская часть $+\,v_r/r$);
плоская часть уже в парной сумме. Девиаторная работа берёт
$s^{zz}\dot\varepsilon^{zz}$ с $\dot\varepsilon^{zz}=v_r/r$. Если у реальной частицы
$r=0$ точно, `finalize_solid` бросает исключение (осевой член не определён);
в парном ядре деления на $r_j$ без защиты, поэтому закрыто только условием
в `finalize_solid` для реальных частиц.

**Ось симметрии (ghost-частицы).** Для реальных частиц с $x^0<2h_i$
(`build_ghost_set`) создаётся зеркальная копия: $r\to-r$, $v_r\to-v_r$,
$s^{xy}\to-s^{xy}$, остальные поля (масса, $\rho$, $p$, $c$, $h$, остальные
напряжения, номер тела) — те же; $p$ и $c$ копируются без повторного вызова УРС.
Ghost участвует в парах с реальными ($i$ реальна, $j$ — ghost), вклад получает
только $i$; пары ghost–ghost отбрасываются. Источники hoop для ghost берутся с
$r_j<0$, что даёт нечётное продолжение $(S_{rr}-S_{\theta\theta})/r$ и чётное — $S_{rz}/r$,
согласованное с зеркалированием. Порог $2h_i$ соответствует носителю
кубического сплайна; при `skf = 2, 3` (носитель $3h$) часть соседей,
лежащих за порогом, зеркальных копий не имеет.

### 10.2. Что проверяет `tests/test_axisym.cpp`

Сборка с `SPH_DIM=2`; регулярная решётка $14\times14$ с шагом $0.1$ мм, B4C
(`itype = 10`), кубический сплайн, пробная частица в глубине, конвейер —
тот же, что в `single_step`, без контакта. Проверки основаны на симметрии
и не зависят от точности SPH-аппроксимации:

| Блок | Состояние | Проверяется |
|---|---|---|
| A1–A4 | $v_r=0$, $v_z=-kz$ | $\dot s^{xx}=\dot s^{zz}$, сумма диагонали $=0$, $\dot s^{yy}=-2\dot s^{xx}$, $\dot s^{yy}<0$ (допуск $10^{-12}$) |
| B | $v_r=ar$, $v_z=0$ | то же для радиального расширения; ловит $\dot\varepsilon^{zz}=0$ (допуск $10\%$ и $10^{-6}$) |
| C1–C3 | однородное объёмное сжатие | все три $\dot s\approx0$ (допуск 5 %) |
| D | покой, однородный девиатор | $\mathrm dv_r/\mathrm dt=0$ (ловит $2s^{xx}+s^{zz}$ вместо $s^{xx}-s^{zz}$) |
| E | покой, однородный $S_{rz}$ | $\mathrm dv_z/\mathrm dt=S_{rz}/(\rho r)$, $\mathrm dv_r/\mathrm dt=0$ |
| F1–F6 | поиск пар | 1 и 7 потоков: примитивы и `real_pairs` совпадают с последовательной сборкой побайтно; `contact_pairs` пуст для одного тела |
| G1–G3 | одна частица, $p=P$ | $\dot e=-(P/\rho)v_r/r$ с допуском $10^{-14}$; $\dot\rho$ не зависит от $p$ |
| G4 | сдвиг частиц | `r0` не меняется; `energy_mass` $=2\pi r_0m$ |
| G5 | гомологическое расширение $v_r=kr$, $\lambda=1,\,1.2$ | диагностика парного баланса энергии (веса $1$, $2\pi r_0$, $2\pi r$); ворот нет, проверяется только однородность $p$ |

G5 печатает данные для анализа, но энергетического баланса парной суммы в
осесимметрии не гарантирует: с плоской массой и двумерным ядром сумма
$\sum m_i\mathbf v_i\cdot\dot{\mathbf v}_i$ с весом $2\pi r$ не нулевая.

## 11. Параметры сводно

| Символ | Ключ / поле | По умолчанию | Единицы / допустимые значения |
|---|---|---|---|
| тип ядра | `skf` | 1 | 1 — кубический, 2 — гауссово, 3 — квинтическое |
| форма давления | `pa_sph` | 2 | 1, 2 (иное трактуется как 2) |
| плотность | `density_method` | 1 | 1, 2, 3 (иное трактуется как 3) |
| длина сглаживания | `sle` | 1 | 0, 1 (2 — ошибка) |
| $\alpha,\beta,\eta$ | `alpha`, `beta`, `etq` | 1, 1, 0.01 | — |
| $\varepsilon$ (XSPH) | `epsilon` | 0.3 | — |
| режим XSPH | `xsph_mode` | `velocity` | `velocity` \| `position` |
| $\xi$, $\varepsilon_h$ | `xi`, `epsilon_h` (не читаются из входного файла) | 0.001, 0.01 | — |
| изоляция тел | `sph_same_body_only` | 0 | 0/1; требует `contact_force_type != 0` |
| осесимметрия, `dim` | `axis_sim` (constexpr), `SPH_DIM` | true, из сборки | — |

Ключ `sph_same_body_only = 1` отбрасывает SPH-примитивы между частицами разных
тел (такие пары остаются только в списке контакта).

## Реализация в коде

| Модель | Файл | Функция / структура |
|---|---|---|
| Кубическое, гауссово, квинтическое ядра | [kernel.hpp](../../../include/kernel.hpp) | `kernel_cubic`, `kernel_gauss`, `kernel_quintic`, `kernel` |
| Поиск пар, критерий $r<k_sh_{ij}$, сборка примитивов | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `build_solid_sph_search`, `SolidSphPrimitive` |
| Ghost-частицы у оси | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `build_ghost_set`, `build_solid_extended` |
| Неразрывность (методы 1–3) | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `compute_solid_sph_contributions`, блок «неразрывность» |
| $\dot\varepsilon$, $\omega$ | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `compute_solid_sph_contributions`, блок «скорости деформаций/поворота» |
| Импульс и энергия давления (`pa_sph`), hoop-члены | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `compute_solid_sph_contributions`, блок «импульс + энергия» |
| Искусственная вязкость | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `compute_solid_sph_contributions`, блок «искусственная вязкость» |
| XSPH: сумма | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | блок «XSPH»; масштабирование на `epsilon` — [single_step.cpp](../../../src/single_step.cpp) |
| XSPH: применение | [ode_solver.hpp](../../../include/ode/ode_solver.hpp) | `pc_finalize_velocity_position_with_av`, `pc_finalize_velocity_position_xsph_position` |
| Осевые члены, девиаторная работа, Яуманн | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `finalize_solid` |
| УРС (диспетчер) | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp), [int_force.cpp](../../../src/int_force.cpp) | `solid_eos`, `compute_solid_derived`, `urs` |
| Длина сглаживания | [h_upgrade.cpp](../../../src/h_upgrade.cpp), [postprocessing.hpp](../../../include/interactions/postprocessing.hpp) | `h_upgrade`, `update_smoothing_length` |
| Масса кольца для энергии | [particle_system.hpp](../../../include/particle_system.hpp) | `init_reference_radius`, `energy_mass` |
| Тест осесимметрии | [test_axisym.cpp](../../../tests/test_axisym.cpp) | `test_*` (A–G5) |
