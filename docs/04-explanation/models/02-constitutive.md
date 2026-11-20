# Определяющие соотношения: уравнения состояния, прочность, разрушение

Глава описывает, как код замыкает систему уравнений движения твёрдого тела: по
плотности $\rho_i$ и удельной внутренней энергии $e_i$ (в коде `sys.u`) уравнение
состояния (УРС) даёт давление $p_i$ и скорость звука $c_i$; закон Гука в скоростях
с объективной производной Яуманна задаёт эволюцию девиатора напряжений
$s^{\alpha\beta}_i$; модели прочности (Мизес, JH-2, Джонсон–Кук) ограничивают
девиатор и, для двух последних, накапливают повреждение. Тензор напряжений
$\sigma^{\alpha\beta}=-p\,\delta^{\alpha\beta}+s^{\alpha\beta}$. Как $p$ и $s$ входят
в уравнения движения и энергии, описано в [01-sph-core.md](01-sph-core.md);
пластичность нитей здесь не рассматривается.

Все материальные константы лежат в структуре
[`Materials`](../../../include/materials.hpp) (подобъект `SimulationConfig::mats`).
Через файл конфигурации меняются только `jc_heating` и `jc_chi`
([input.cpp](../../../src/input.cpp)); остальные константы зашиты в код (значения
по умолчанию в `Materials`/`make_jc_library`). В чекпоинте сверяются массивы
`G, sigy, aa, zz, gamm` и `const_JH_2` ([checkpoint.cpp](../../../src/checkpoint.cpp)); библиотека
`make_jc_library` и флаги JH-2 не сверяются, сверяются только ключи `jc_heating`, `jc_chi`.

## 1. Конвейер замыкания на шаге

На каждом вызове правой части:

1. `compute_solid_derived` → `solid_eos` ([solid_sph.hpp](../../../include/interactions/solid_sph.hpp))
   вычисляет $p_i, c_i$ для всех реальных частиц; зеркальные (осевые) частицы берут
   $p,c$ источника (УРС — чистая функция аргументов). Старый диспетчер `urs()`
   ([int_force.cpp](../../../src/int_force.cpp)) устроен идентично.
2. `finalize_solid` формирует скорости девиатора $\dot s^{\alpha\beta}$ (закон Гука +
   Яуманн, разд. 3) и вклад в $\dot e$.
3. После шага интегрирования (после корректора) `apply_material_strength`
   ([postprocessing.hpp](../../../include/interactions/postprocessing.hpp)) применяет
   прочность по `itype` частицы:

| `itype` | материал | УРС | прочность |
|---|---|---|---|
| 1 | идеальный газ | `p_gas` | нет (девиатор не считается, $\texttt{itype}>2$) |
| 2 | вода | `p_art_water` | нет |
| 3, 4 | тантал, алюминий | `eos_mi_grunaizen` | Мизес (`prochnost_mises`) |
| 5–9 | LiF, ПУ 160/320/626, PMMA | `eos_ua_teta` | Мизес |
| 10, 11 | B₄C, Al₂O₃ | `urs_cem_JH2` | JH-2 (`prochnost_cem_JH2`) |
| 12–15 | сталь 4340, сталь 1006, Al 7039, Al 5083-H116 | `eos_jc` (Linear/Shock) | Джонсон–Кук (`apply_jc_strength`) |

Неизвестный `itype` — исключение (`solid_itype_known`, `default:` в `solid_eos`).
Индексация массивов `G, sigy, aa, zz, gamm`: `itype - n_Me`, `n_Me = 3`; JH-2:
`itype - 10`; JC: `itype - n_JC`, `n_JC = 12`.

## 2. Уравнения состояния

### 2.1. Идеальный газ и «вода» (itype 1, 2)

`p_gas`: $p=(\gamma-1)\rho e$, $\gamma=1{,}4$; $c=\sqrt{\gamma p/\rho}$ при $p\ge0$,
иначе $c=100$ м/с. `p_art_water`: $p=B[(\rho/1000)^7-1]$, $B=1000$ Па (значение в коде),
$c=\sqrt{7B(\rho/1000)^6/1000}$. Для этих материалов прочность не применяется; в
расчётах с твёрдыми телами не используются.

### 2.2. УА Тэта (itype 5–9)

Функция `eos_ua_teta`. Обозначения: $w=\rho/\rho_0$, $a=$ `aa`, $z=$ `zz`
($\rho_0$ — `sys.rho_0[i]` частицы):

$$
p=a\,(w^{z}-1),\qquad c^{2}=\frac{a z}{\rho_0}\,w^{z-1}.
$$

Энергия в давление не входит. Если $c^2\le0$ (практически недостижимо), $c=100$ м/с.
При $w=1$: $\rho_0c^2=az=K$ — объёмный модуль; он же используется в
`elastic_by_itype` для получения $E,\nu$ (для контакта нитей с твёрдым телом).

### 2.3. Ми–Грюнайзен (itype 3, 4)

Функция `eos_mi_grunaizen`. Дополнительно параметр `gamm` $=\gamma$ (в код входит как
$\gamma-1$ перед энергией, т.е. параметр Грюнайзена $\Gamma=\gamma-1$):

$$
\begin{aligned}
p &= p_c(w)+(\gamma-1)\rho e,\\
p_c(w) &= a\left[\frac{z-\gamma}{z-1}\,w^{z}+\frac{(\gamma-1)z}{z-1}\,w-\gamma\right].
\end{aligned}
$$

Свойства (проверены подстановкой): $p_c(1)=0$, $p_c'(1)=az$ (тот же объёмный модуль,
что у Тэта). Скорость звука в коде:

$$
c^{2}=\frac{az}{\rho_0}\,w^{z-1}+\gamma\,\frac{p-p_x}{\rho},\qquad p_x=a\,(w^{z}-1),
$$

при $c^2\le0$ — $c=100$ м/с. Форма $c^2$ эвристическая: на $\rho=\rho_0$, $e=0$ она
даёт $c^2=az/\rho_0$; формального вывода из $\partial p/\partial\rho|_s$ в коде нет.

### 2.4. Линейный УРС (JC, Linear)

`eos_linear`: $\mu=\rho/\rho_0-1$,

$$
p=K\mu,\qquad c=\sqrt{K/\rho_0}.
$$

Энергия не входит.

### 2.5. УРС Shock (Ми–Грюнайзен по адиабате $U_s=c_0+S_1u_p$, JC)

`eos_shock`. Обозначения: $\mu=\rho/\rho_0-1$, $A=\rho_0c_0^2$, $\Gamma=\Gamma_0\rho_0$
(как в Autodyn: $\rho\Gamma$ принято постоянным $=\Gamma_0\rho_0$), $D=1-(S_1-1)\mu$.

При $\mu>0$ (сжатие):

$$
p_H=\frac{A\mu(1+\mu)}{D^{2}},\qquad e_H=\frac{c_0^{2}\mu^{2}}{2D^{2}};
$$

при $\mu\le0$: $p_H=A\mu$, $e_H=0$. Давление

$$
p=p_H+\Gamma\,(e-e_H).
$$

Скорость звука: $c^2=\dfrac{1}{\rho_0}\left(\dfrac{dp_H}{d\mu}-\Gamma\dfrac{de_H}{d\mu}\right)+\dfrac{p}{\rho^{2}}\Gamma$,
где $\dfrac{dp_H}{d\mu}=A\dfrac{(1+2\mu)D+2(S_1-1)\mu(1+\mu)}{D^{3}}$,
$\dfrac{de_H}{d\mu}=\dfrac{c_0^2\mu}{D^3}$ (при $\mu\le0$: $A$ и $0$). Если $c^2$ не
конечно или $\le0$ — $c=c_0$.
Сингулярность $D\to0$ ($\mu\to1/(S_1-1)$) устранена ограничением $D\ge k_{D\min}=10^{-3}$.

### 2.6. Единый УРС JC: `eos_jc`

Диспетчер по `m.eos` (`EosKind::Linear` / `Shock`), затем два ограничителя:

1. $p\ge p_{\min}$ (поле `JcMaterial::p_min`, Па; по умолчанию $-1{,}01\cdot10^{20}$ —
   фактически без ограничения; задан только для стали 4340: $-5{,}5\cdot10^{9}$ Па —
   откольная прочность 4340 RC30, критерий $\sigma_s=-55$ кбар по Bless 1981,
   AFWAL-TR-81-4040; ранее стояло подобранное $-4\cdot10^{9}$ Па);
2. если повреждение $D\ge1$ и $p<0$, то $p=0$ (разрушенная частица не воспринимает
   растяжение).

$\rho_0$ берётся из частицы (`sys.rho_0[i]`), а не из `m.rho0`. Параметр `rho0`
таблицы используется только в `elastic_by_itype` для $K=\rho_0c_0^2$ (Shock).

### 2.7. УРС керамики JH-2: `urs_cem_JH2`

Константы $K_1,K_2,K_3,T$ — столбцы 13–16 `const_JH_2`. $\mu=\rho/\rho_0-1$,
$D$ — повреждение, $\Delta p$ (`sys.dp`) — добавка давления от булкинга (разд. 5.4).

Целая (до разрушения, $D=0$):

$$
p=\begin{cases}K_1\mu+K_2\mu^2+K_3\mu^3,&\mu\ge0,\\ \max\,(K_1\mu,\,-T),&\mu<0,\end{cases}
$$

Повреждённая ($D\neq0$):

$$
p=\begin{cases}K_1\mu+K_2\mu^2+K_3\mu^3+\Delta p,&\mu\ge0,\\ \max\,(K_1\mu,\,-T(1-D)),&\mu<0.\end{cases}
$$

Скорость звука: $c=\sqrt{(K_1+2K_2\mu+3K_3\mu^2)/\rho_0}$ при $\mu\ge0$ и
$c=\sqrt{K_1/\rho_0}$ при $\mu<0$. Энергия в давление не входит.
При $D=1$ предел растяжения равен 0.

## 3. Упругая часть девиатора: закон Гука в скоростях и производная Яуманна

Реализация — конец `finalize_solid` ([solid_sph.hpp](../../../include/interactions/solid_sph.hpp)),
выполняется для частиц с $\texttt{itype}>2$.

### 3.1. Кинематика

Скорости деформации и поворота — «сырые» парные суммы по соседям $j$ (при
$\mathbf v_{ji}=\mathbf v_j-\mathbf v_i$, градиент ядра $\nabla_iW_{ij}$, без
поправки на несовершенство ядра):

$$
\dot\varepsilon^{\alpha\beta}_i=\frac12\sum_j\frac{m_j}{\rho_j}\left(v_{ji}^{\alpha}\,\partial^{\beta}_iW_{ij}+v_{ji}^{\beta}\,\partial^{\alpha}_iW_{ij}\right),\quad
\omega^{\alpha\beta}_i=\frac12\sum_j\frac{m_j}{\rho_j}\left(v_{ji}^{\alpha}\,\partial^{\beta}_iW_{ij}-v_{ji}^{\beta}\,\partial^{\alpha}_iW_{ij}\right).
$$

Диагональные компоненты $\dot\varepsilon^{\alpha\alpha}=\sum_j\frac{m_j}{\rho_j}v_{ji}^\alpha\partial^\alpha W$.
В коде имена: `exx, eyy, ezz, exy, exz, eyz` и `rxy, rxz, ryz`
(`rxy`$=\omega^{xy}$). Компоненты есть только для $\alpha,\beta<d$; в 2D
(плоскость) $\dot\varepsilon^{zz}=0$.

**Осесимметрия** ($d=2$, `axis_sim`): окружная скорость деформации — почастичная формула
$\dot\varepsilon^{zz}_i=v_r/r$ (а не парная сумма), в уравнение непрерывности
добавляется $-\rho v_r/r$, а в уравнение энергии — $-(p/\rho)\,v_r/r$. При $r=0$
`finalize_solid` бросает исключение.

### 3.2. Закон Гука

Скорость изменения девиатора:

$$
\dot s^{\alpha\beta}=2G\left(\dot\varepsilon^{\alpha\beta}-\tfrac13\delta^{\alpha\beta}\dot\varepsilon^{\gamma\gamma}\right)+s^{\alpha\gamma}\omega^{\beta\gamma}+s^{\gamma\beta}\omega^{\alpha\gamma}.
$$

В коде след $\dot\varepsilon^{\gamma\gamma}$ заменён величиной $-\dot\rho/\rho$ из
уравнения непрерывности (`drho` уже содержит осевой член и возможную диффузию
плотности из ядра), т.е. $\dot\varepsilon^{\alpha\alpha}-\frac13\dot\varepsilon^{\gamma\gamma}\to\dot\varepsilon^{\alpha\alpha}+\dot\rho/(3\rho)$.
Модуль сдвига $G$ — `effective_shear_modulus`: равен нулю для JC-частицы с $D\ge1$
(разрушенная не набирает девиатор), иначе `shear_modulus_by_itype`
(`G[itype-3]`, `const_JH_2[itype-10][1]` или `jc[itype-12].s.G`).

Компоненты в коде (3D и 2D, $d>1$):

$$
\begin{aligned}
\dot s^{xx}&=2G\left(\dot\varepsilon^{xx}+\tfrac{\dot\rho}{3\rho}\right)+2\left(s^{xy}\omega^{xy}+s^{xz}\omega^{xz}\right),\\
\dot s^{yy}&=2G\left(\dot\varepsilon^{yy}+\tfrac{\dot\rho}{3\rho}\right)+2\left(s^{yz}\omega^{yz}-s^{xy}\omega^{xy}\right),\\
\dot s^{zz}&=2G\left(\dot\varepsilon^{zz}+\tfrac{\dot\rho}{3\rho}\right)-2\left(s^{xz}\omega^{xz}+s^{yz}\omega^{yz}\right),\\
\dot s^{xy}&=2G\dot\varepsilon^{xy}+\omega^{xy}(s^{yy}-s^{xx})+s^{xz}\omega^{yz}+s^{yz}\omega^{xz},\\
\dot s^{yz}&=2G\dot\varepsilon^{yz}-s^{xy}\omega^{xz}-s^{xz}\omega^{xy}+\omega^{yz}(s^{zz}-s^{yy}),\\
\dot s^{xz}&=2G\dot\varepsilon^{xz}+\omega^{xz}(s^{zz}-s^{xx})-s^{xy}\omega^{yz}+s^{yz}\omega^{xy}.
\end{aligned}
$$

Эти выражения совпадают с общей формулой выше (проверено для $xx$ и $xy$). Знак при
$s^{yz}\omega^{yz}$ в $\dot s^{zz}$ исправлен (комментарий в коде и
`tests/test_kinematics_3d.cpp`): прежний знак нарушал бесследовость поворотной части.
Для $d=1$ вращения нет: $\dot s^{xx}=2G(\dot\varepsilon^{xx}+\dot\rho/3\rho)$,
$\dot s^{yy}=\dot s^{zz}=2G\,\dot\rho/(3\rho)$.
В осесимметрии поворотные члены $\omega^{xz},\omega^{yz}=0$, поэтому $\dot s^{zz}$
поворота не содержит; поправки полярной системы координат к производной Яуманна для
окружной компоненты не вводятся.

Интегрирование $s$ по времени (предиктор–корректор) — в
[ode_solver.hpp](../../../include/ode/ode_solver.hpp) и 01-sph-core; прочность
применяется уже к скорректированным $s$.

### 3.3. Энергия

Вклад девиаторной работы: $\dot e_i\mathrel{+}=\frac1{\rho_i}s^{\alpha\beta}\dot\varepsilon^{\alpha\beta}$ (в 2D плоскости и осесимметрии
$s^{xx}\dot\varepsilon^{xx}+s^{yy}\dot\varepsilon^{yy}+s^{zz}\dot\varepsilon^{zz}+2s^{xy}\dot\varepsilon^{xy}$,
в 3D — со всеми сдвигами). Давление и вязкость — в парной сумме, умноженной на 0,5.

## 4. Упругопластическая модель Мизеса (itype 3–9)

### 4.1. Условие текучести и возврат

Функции `prochnost_mises`, `sigma_i` ([ppochnost.cpp](../../../src/ppochnost.cpp)).
Квадрат второго инварианта девиатора (функция `sigma_i`, зависит от `dim` и `axis_sim`):

| геометрия | $J_2$ в коде |
|---|---|
| $d=1$ | $\tfrac34(s^{xx})^2$ |
| $d=2$, осесимметрия | $(s^{xx})^2+s^{xx}s^{yy}+(s^{yy})^2+(s^{xy})^2$ |
| $d=2$, плоскость | $(s^{xx})^2+s^{xx}s^{yy}+(s^{yy})^2+(s^{xy})^2$ (как в осесимметрии, $s^{zz}=-(s^{xx}+s^{yy})$) |
| $d=3$ | $(s^{xx})^2+s^{xx}s^{yy}+(s^{yy})^2+(s^{xy})^2+(s^{xz})^2+(s^{yz})^2$ |

(для бесследового девиатора $J_2=\tfrac12 s:s$, эквивалентное напряжение $\sigma_{eq}=\sqrt{3J_2}$).
Радиальный возврат по Мизесу (идеальная пластичность, без упрочнения):

$$
\text{если }\;J_2>\frac{\sigma_y^{2}}{3}:\quad s^{\alpha\beta}\leftarrow s^{\alpha\beta}\,\frac{\sigma_y/\sqrt3}{\sqrt{J_2}}\quad\text{(все 6 компонент).}
$$

Модуль сдвига в самой функции не используется. Повреждения нет, пластическая деформация
не накапливается; в `sys.e_pl` Мизес ничего не пишет. Под `SPH_ENERGY_LOG` упругая энергия
сдвига $\bigl((s{:}s)_{\rm до}-(s{:}s)_{\rm после}\bigr)/(4G\rho)$ учитывается как диссипация (справочная величина).

Замечание о плоском случае: в `sigma_i` плоская форма исправлена (этап A исправления JH-2,
`notes/jh2/jh2-fix-plan.md`). В `jc::j2_of_deviator` (`Geom::D2Plane`) осталась прежняя
$J_2=(s^{xx})^2+(s^{xy})^2$, верная лишь при $s^{yy}=-s^{xx}$; ветка недостижима, пока
`axis_sim` — `constexpr true`.

### 4.2. Параметры (Мизес, УА Тэта, Ми–Грюнайзен)

Значения идентичны [registry.md](../../03-reference/materials/registry.md) (сверено с
кодом). Индекс массива $=\texttt{itype}-3$.

| `itype` | материал | УРС | $G$, Па (`G`) | $\sigma_y$, Па (`sigy`) | $a$, Па (`aa`) | $z$ (`zz`) | $\gamma$ (`gamm`) |
|---|---|---|---|---|---|---|---|
| 3 | тантал | Ми–Грюнайзен | $6{,}9\cdot10^{10}$ | $4{,}5\cdot10^{8}$ | $1{,}4506\cdot10^{11}$ | 1,987 | 2,6 |
| 4 | алюминий | Ми–Грюнайзен | $2{,}6\cdot10^{10}$ | $2{,}4\cdot10^{8}$ | $1{,}97\cdot10^{10}$ | 4,2 | 1,97 |
| 5 | LiF | Тэт | $5{,}5\cdot10^{10}$ | $2{,}0\cdot10^{7}$ | $1{,}6659\cdot10^{10}$ | 4,272 | — |
| 6 | ПУ 160 | Тэт | $5{,}27\cdot10^{7}$ | $1{,}0\cdot10^{6}$ | $2{,}1953\cdot10^{8}$ | 1,245 | — |
| 7 | ПУ 320 | Тэт | $1{,}836\cdot10^{8}$ | $3{,}0\cdot10^{6}$ | $3{,}595\cdot10^{8}$ | 1,736 | — |
| 8 | ПУ 626 | Тэт | $4{,}0\cdot10^{8}$ | $4{,}0\cdot10^{6}$ | $3{,}0227\cdot10^{8}$ | 3,5 | — |
| 9 | PMMA | Тэт | $3{,}0\cdot10^{9}$ | $7{,}7\cdot10^{7}$ | $1{,}8797\cdot10^{9}$ | 4,668 | — |

Плотность $\rho_0$ для этих материалов задаётся при построении частиц (`sys.rho_0`),
в `Materials` не хранится.

## 5. Керамика JH-2 (itype 10, 11)

### 5.1. Физическая модель

Реализована схема Джонсона–Холмквиста (JH-2) для хрупкого материала с
накоплением повреждения $D\in[0,1]$: предел прочности целого материала
$\sigma_i^*$ и разрушенного $\sigma_f^*$, интерполяция по $D$, скоростной множитель,
булкинг (рост давления при разрушении) и предел растяжения $T(1-D)$.
Нормировка — на $\sigma_{HEL}$ и $P_{HEL}$.
Файл: [ppochnost.cpp](../../../src/ppochnost.cpp), функция `prochnost_cem_JH2`
(вызывается из `apply_material_strength` после корректора); УРС — разд. 2.7.

### 5.2. Предел текучести

$P_n=p/P_{HEL}$ (текущее $p$ частицы), $T_n=$ `const_JH_2[...][4]` (уже нормированное).

$$
\sigma_{0}=\begin{cases}\sigma_{HEL}A\,(P_n+T_n)^{N},&P_n>0\ \text{или}\ P_n+T_n>0\;(P_n\le0),\\0,&P_n\le0,\ P_n+T_n\le0,\end{cases}\qquad
\sigma_{f}=\begin{cases}\sigma_{HEL}B\,P_n^{M},&P_n>0,\\0,&P_n\le0,\end{cases}
$$

$\sigma_f\leftarrow\min(\sigma_f,\ s_{f\max}\sigma_{HEL})$. Текущий предел:

$$
\sigma_Y=\bigl[\sigma_0-D(\sigma_0-\sigma_f)\bigr]\bigl(1+C\ln\dot\varepsilon_{eq}\bigr)\quad(\text{множитель применяется только при }\dot\varepsilon_{eq}>0),
$$

где $\dot\varepsilon_{eq}=\frac{\sqrt2}{3}\sqrt{(\dot\varepsilon^{xx}-\dot\varepsilon^{yy})^2+(\dot\varepsilon^{yy}-\dot\varepsilon^{zz})^2+(\dot\varepsilon^{zz}-\dot\varepsilon^{xx})^2+6(\dot\varepsilon_{xy}^2+\dot\varepsilon_{xz}^2+\dot\varepsilon_{yz}^2)}$
из `sys.exx..eyz` (скорости, записанные `finalize_solid`; нормировка на $\dot\varepsilon_0=1$ с$^{-1}$ неявная).
При $\dot\varepsilon_{eq}<1$ множитель $<1$, как в оригинале J&H 1999; флаг `Materials::jh2_clamp_low_rate`
(по умолчанию `false`) включает отсечку $\dot\varepsilon^*=\max(\dot\varepsilon_{eq},1)$.
Предел считает функция `jh2_strength`.
Возврат: если $\sqrt{J_2}>\sigma_Y/\sqrt3$, все компоненты $s$ умножаются на
$\sigma_Y/(\sqrt3\sqrt{J_2})$ (флаг `flag_e = false`).

### 5.3. Накопление повреждения

Если возврат сработал и $D<1$, приращение эквивалентной пластической деформации берётся
из радиального возврата (Wilkins):

$$
\Delta\varepsilon^p=\frac{\sigma_{eq}^{tr}-\sigma_Y}{3G},\qquad \sigma_{eq}^{tr}=\sqrt{3J_2^{tr}},
$$

где $J_2^{tr}$ — пробный (до возврата) второй инвариант. Затем

$$
\Delta D=\frac{\Delta\varepsilon^p}{\varepsilon_f^p},\quad \varepsilon_f^p=\max\bigl(D_1\,(P_n+T_n)^{D_2},\ \varepsilon_{f,\min}\bigr),\quad e_{pl}\mathrel{+}=\Delta\varepsilon^p,\quad D\leftarrow\min(D+\Delta D,1),
$$

$\varepsilon_{f,\min}$ — `Materials::jh2_efmin` (по умолчанию 0, как в оригинале).
Если $P_n+T_n\le0$ (всестороннее растяжение при текучести) — $D=1$ сразу (в этой ветке
$e_{pl}$ не наращивается). Теперь $e_{pl}$ у JH-2, как и у JC, — накопленная эквивалентная
пластическая деформация $\sum\Delta\varepsilon^p$ (раньше — интеграл скорости с ошибкой в 4 раза).

### 5.4. Разрыхление (булкинг)

Добавка давления растёт только если повреждение выросло за шаг ($D_{\rm new}>D_{\rm old}$) и $\mu\ge0$.
Потеря упругой энергии сдвига при том же $P_n$ и $\dot\varepsilon_{eq}$ (J&H 1999, (11)–(14);
Cronin и др., 2003):

$$
\Delta U=\frac{\sigma_Y(D_{\rm old})^2-\sigma_Y(D_{\rm new})^2}{6G}\ \ge 0,\qquad
\Delta p_{\rm new}=-K_1\mu+\sqrt{(K_1\mu+\Delta p)^{2}+2\beta K_1\,\Delta U},\quad\beta=1.
$$

Результат записывается в `sys.dp` и входит в давление (разд. 2.7) при $D\ne0$. Ограничение
$\mu\ge0$ оставлено сознательно (в оригинале его нет).

> До исправления (ветка `fiber_johnson_holmquist`, этапы 0, A–E плана `notes/jh2/jh2-fix-plan.md`,
> см. [ADR-010](../decisions/adr-010-jh2-fix.md))
> упругая часть скорости деформации делилась на $G$ вместо $2G$, а энергии $u_0$, $u_1$
> в булкинге имели разную размерность, из-за чего $\Delta p$ практически оставалось нулевым.
> Константы Al₂O₃ (`const_JH_2[1]`) не проверены — см. [реестр материалов](../../03-reference/materials/registry.md).

### 5.5. Параметры JH-2

Порядок столбцов `const_JH_2`: `[0]` $\rho$, `[1]` $G$, `[2]` $P_{HEL}$, `[3]` $\sigma_{HEL}$,
`[4]` $T_n$, `[5]` $A$, `[6]` $N$, `[7]` $B$, `[8]` $M$, `[9]` $s_{f\max}$, `[10]` $C$,
`[11]` $D_1$, `[12]` $D_2$, `[13]` $K_1$, `[14]` $K_2$, `[15]` $K_3$, `[16]` $T$.

| символ | ключ | B₄C (itype 10) | Al₂O₃ (itype 11) | единицы |
|---|---|---|---|---|
| $\rho$ | [0] | 2517 | 3555 | кг/м³ |
| $G$ | [1] | $1{,}994\cdot10^{11}$ | $1{,}165\cdot10^{11}$ | Па |
| $P_{HEL}$ | [2] | $8{,}71\cdot10^{9}$ | $2{,}626\cdot10^{9}$ | Па |
| $\sigma_{HEL}$ | [3] | $1{,}544\cdot10^{10}$ | $3{,}561\cdot10^{9}$ | Па |
| $T_n$ | [4] | 0,03 | 0,137 | — |
| $A$ | [5] | 0,927 | 0,930 | — |
| $N$ | [6] | 0,67 | 0,600 | — |
| $B$ | [7] | 0,70 | 0,110 | — |
| $M$ | [8] | 0,85 | 0,600 | — |
| $s_{f\max}$ | [9] | 0,2 | 1,0 | — |
| $C$ | [10] | 0,005 | 0 | — |
| $D_1$ | [11] | 0,001 | 0,005 | — |
| $D_2$ | [12] | 0,5 | 1,0 | — |
| $K_1$ | [13] | $2{,}305\cdot10^{11}$ | $1{,}692\cdot10^{11}$ | Па |
| $K_2$ | [14] | $-5{,}93\cdot10^{11}$ | 0 | Па |
| $K_3$ | [15] | $2{,}8\cdot10^{12}$ | 0 | Па |
| $T$ | [16] | $2{,}6\cdot10^{8}$ | $2{,}0\cdot10^{8}$ | Па |
| $\beta$ | константа в коде | 1 | 1 | — |

Реестр совпадает с кодом (см. врезку про смешанный набор B₄C в
[registry.md](../../03-reference/materials/registry.md)). Набор намеренно смешанный: упругие константы
по Eagle Picher (Grady & Moody, SAND96-0551), прочностные прежние (комментарий в `materials.hpp`).
Столбец `[0]` (плотность) в вычислениях не читается (используется `sys.rho_0`). Для Al₂O₃
$T_n=0{,}137$ не равно $T/P_{HEL}=0{,}2\cdot10^9/2{,}626\cdot10^9\approx0{,}076$ (для B₄C равенство
выполняется: $0{,}03$); $T_n$ и $T$ используются в разных формулах (прочность/повреждение и УРС).

## 6. Модель Джонсона–Кука (itype 12–15)

### 6.1. Физический смысл

Вязкопластическая модель с изотропно-деформационным упрочнением, логарифмической
скоростной чувствительностью и термическим разупрочнением; повреждение накапливается
по деформации разрушения, зависящей от трёхосности, скорости и температуры.
Чистые функции — [johnson_cook.hpp](../../../include/strength/johnson_cook.hpp)
(пространство имён `jc`), применение к частице —
[jc_strength.hpp](../../../include/interactions/jc_strength.hpp)
(`apply_jc_strength`), нач. состояние — [jc_state.hpp](../../../include/jc_state.hpp)
(`init_jc_state`: $T=T_r$, $\dot\varepsilon=0$; при рестарте не вызывается).

### 6.2. Предел текучести

$$
Y=\bigl[A+B\,\varepsilon_p^{\,n}\bigr]\,R(\dot\varepsilon)\,\Theta(T),
$$

$\varepsilon_p=\max(e_{pl},0)$ (при $B=0$ упрочнение отсутствует). Скоростной множитель
(`rate_factor`, $x=\dot\varepsilon/\dot\varepsilon_0$):

$$
R=\begin{cases}1,&x\le1,\\1+C\ln x&\text{(форма `FirstOrder`, во всей библиотеке)},\\ x^{C}&\text{(форма `Exponential`, в библиотеке не используется)}.\end{cases}
$$

Температурный множитель (`thermal_factor`), $T^*=(T-T_r)/(T_m-T_r)$ (при $T_m\le T_r$ $T^*=0$):

$$
\Theta=\begin{cases}1,&T^*\le0,\\1-(T^*)^{m},&0<T^*<1,\\0,&T^*\ge1\ \text{(расплав)}.\end{cases}
$$

Запись «$1-T^{*m}$» со слайда лекции без вычитания $T_r$ не используется (комментарий в коде).

### 6.3. Возврат на поверхность текучести с упрочнением

`jc_radial_return`. Пробное эквивалентное напряжение
$\sigma_{tr}=\sqrt{3J_2}$ по текущему $s$ (формулы `j2_of_deviator`, как в `sigma_i`, с
выбором геометрии `jc_build_geom()` по `dim`/`axis_sim`).

1. Если $\sigma_{tr}\le Y(\varepsilon_p,\dot\varepsilon,T)$ — упругий шаг, $s$ не меняется.
2. Если $\Theta=0$ (расплав) — девиатор обнуляется, $\Delta\varepsilon_p=0$.
3. Иначе ищется $\Delta\varepsilon_p\in[0,\sigma_{tr}/3G]$ из условия
   $f(\Delta)=\sigma_{tr}-3G\Delta-Y(\varepsilon_p+\Delta)=0$ безопасным Ньютоном (при выходе
   из скобки или нефинитной производной — бисекция; допуск $10^{-12}Y_0$, до 50 итераций,
   затем откат на чистую бисекцию до относительной ширины $10^{-14}$; если не сошлось —
   одно предупреждение в stderr на запуск).
4. Новый предел $Y=\sigma_{tr}-3G\Delta$, масштаб $s\leftarrow s\cdot Y/\sigma_{tr}$ для
   всех шести компонент, $e_{pl}\mathrel{+}=\Delta\varepsilon_p$.

Здесь $G=$ `s.G` — табличный $G$ неповреждённого материала.

### 6.4. Скорость деформации

`equivalent_strain_rate`: $\dot\varepsilon_{eq}=\sqrt{\tfrac23 e'{:}e'}$ по девиатору скорости
полной деформации $e'=\dot\varepsilon-\frac13\mathrm{tr}\,\dot\varepsilon$ (сдвиги входят как $2\dot\varepsilon_{xy}^2$ и т.д.).
Геометрия: $d=1$ — $\mathrm{diag}(\dot\varepsilon^{xx},0,0)$; плоскость — $\dot\varepsilon^{zz}=0$, сдвиг
только $xy$; осесимметрия — $\dot\varepsilon^{xx},\dot\varepsilon^{yy},\dot\varepsilon^{zz}\,(=v_r/r),\dot\varepsilon^{xy}$; 3D — все шесть.
Источник (`RateSource`): `Total` — как выше из `sys.exx..eyz` (во всей библиотеке);
`Plastic` — $\Delta\varepsilon_p/\Delta t$ прошлого шага из `sys.edot` (в коде реализовано,
ключом конфигурации не включается). Результат пишется в `sys.edot`.

### 6.5. Температура и нагрев

Ключ `jc_heating` (применяется ко всем 4 материалам при разборе конфига):

| значение | режим | формула |
|---|---|---|
| `energy` (по умолч.) | `Heating::Energy` | $T=T_r+e/c_p$, пишется в `sys.T`; $e$ — полная внутренняя энергия частицы (включая работу давления и упругую девиаторную работу) |
| `taylor_quinney` | `Heating::TaylorQuinney` | $T$ интегрируется явно: берётся $T$ с начала шага, после возврата $T\mathrel{+}=\chi\,Y\,\Delta\varepsilon_p/(\rho c_p)$ |

Ключ `jc_chi` — доля $\chi\in[0;1]$ пластической работы, переходящая в тепло (по умолч.
0,9; значение вне $[0,1]$ — ошибка разбора). В режиме `energy` $\chi$ не используется.

### 6.6. Разрушение JC

Активно, если `enabled` и все $D_1..D_5<10^{20}$ (константа $\ge10^{20}$ — плейсхолдер
Autodyn «не задано», `jc_damage_active`). Накопление — только при $\Delta\varepsilon_p>0$ и
$Y>10^{-12}A$:

$$
\sigma^*=-\frac{p}{Y},\quad
\varepsilon_f=\max\Bigl\{\bigl[D_1+D_2e^{D_3\sigma^*}\bigr]\,R_f\,(1+D_5T^*_f),\ \varepsilon_{f,\min}\Bigr\},\quad
D\leftarrow\min\Bigl(D+\frac{\Delta\varepsilon_p}{\varepsilon_f},\,1\Bigr).
$$

$R_f=1+D_4\ln\max(\dot\varepsilon/\dot\varepsilon_{0f},1)$ для формы `JC`;
$R_f=(1+\dot\varepsilon/\dot\varepsilon_{0f})^{D_4}$ для формы `Borvik`.
$T_f^*=\max\bigl((T-T_r)/(T_{m,f}-T_r),0\bigr)$ ($T_{m,f}\le0$ означает «взять $T_m$ прочности», предупреждение
в stderr один раз), $\varepsilon_{f,\min}=10^{-6}$, $\dot\varepsilon_{0f}=1$ с$^{-1}$. При $D\ge1$
девиатор сразу обнуляется; на следующих шагах: $D\ge1$ ⇒ девиатор $=0$, $e_{pl}$ не растёт,
$G_{\rm eff}=0$ (разд. 3.2) и растяжение отрезается в `eos_jc` ($p\ge0$). Это единственный
способ «отрыва» в JC: частица остаётся в расчёте, теряя прочность и растягивающее давление;
дополнительного разрыва связей/удаления частиц для JC нет.
Порог $p_{\min}$ (разд. 2.6) задаёт откол по давлению независимо от $D$.

### 6.7. Порядок операций `apply_jc_strength`

1) температура (разд. 6.5); 2) $\dot\varepsilon$ (6.4); 3) если $D\ge1$ — $s=0$, выход;
4) возврат (6.3), $e_{pl}$, при Тейлоре–Куинни — нагрев; 5) повреждение (6.6);
6) под `SPH_ENERGY_LOG` — пластическая работа $m\,Y\,\Delta\varepsilon_p/\rho$ в лог энергии.

### 6.8. Параметры JC (значения из `make_jc_library`)

Единицы СИ. Общие для всех: $\dot\varepsilon_0=1$ с$^{-1}$, форма скорости `FirstOrder`,
источник `Total`, $\chi=0{,}9$ (перекрывается `jc_chi`), нагрев `Energy` (перекрывается
`jc_heating`).

| символ (поле) | STEEL 4340 (12) | STEEL 1006 (13) | AL 7039 (14) | AL5083H116 (15) | единицы |
|---|---|---|---|---|---|
| $\rho_0$ (`rho0`) | 7830 | 7896 | 2770 | 2700 | кг/м³ |
| УРС (`eos`) | Linear | Shock | Shock | Linear | — |
| $K$ (`lin.K`) | $1{,}59\cdot10^{11}$ | — | — | $5{,}833\cdot10^{10}$ | Па |
| $\Gamma_0$ (`shock.gamma0`) | — | 2,17 | 2,0 | — | — |
| $c_0$ (`shock.c0`) | — | 4569 | 5328 | — | м/с |
| $S_1$ (`shock.s1`) | — | 1,49 | 1,338 | — | — |
| $G$ (`s.G`) | $8{,}18\cdot10^{10}$ | $8{,}18\cdot10^{10}$ | $2{,}76\cdot10^{10}$ | $2{,}692\cdot10^{10}$ | Па |
| $A$ | $7{,}92\cdot10^{8}$ | $3{,}5\cdot10^{8}$ | $3{,}37\cdot10^{8}$ | $1{,}67\cdot10^{8}$ | Па |
| $B$ | $5{,}1\cdot10^{8}$ | $2{,}75\cdot10^{8}$ | $3{,}43\cdot10^{8}$ | $5{,}96\cdot10^{8}$ | Па |
| $n$ | 0,26 | 0,36 | 0,41 | 0,551 | — |
| $C$ | 0,014 | 0,022 | 0,010 | 0,001 | — |
| $m$ | 1,03 | 1,00 | 1,00 | 0,859 | — |
| $T_m$ | 1793 | 1811 | 877 | 893 | К |
| $T_r$ | 300 | 300 | 300 | 293 | К |
| $c_p$ | 477 | 452 | 875 | 910 | Дж/(кг·К) |
| разрушение | включено, форма `JC` | выкл. | выкл. | включено, форма `Borvik` | — |
| $D_1$ | 0,05 | $1{,}01\cdot10^{20}$ | $1{,}01\cdot10^{20}$ | 0,178 | — |
| $D_2$ | 3,44 | то же | то же | 0,389 | — |
| $D_3$ | −2,12 | то же | то же | −2,25 | — |
| $D_4$ | 0,002 | то же | то же | 0,147 | — |
| $D_5$ | 0,61 | то же | то же | 0 | — |
| $T_{m,f}$ | 1793 | — | — | 893 | К |
| $p_{\min}$ | $-5{,}5\cdot10^{9}$ | $-1{,}01\cdot10^{20}$ (нет) | $-1{,}01\cdot10^{20}$ (нет) | $-1{,}01\cdot10^{20}$ (нет) | Па |

Таблица совпадает с разделом «Константы Джонсона–Кука» в
[registry.md](../../03-reference/materials/registry.md) (сверено). Все константы
STEEL 4340, включая $D_1..D_5$, $T_{m,f}$ и $\dot\varepsilon_{0,f}$, сверены с библиотекой
Autodyn (`notes/sources/steel-4340-autodyn-constants.md`, скриншоты от 15.11.2026; эти значения
важнее литературных); $D_1..D_5$ совпадают с Johnson & Cook 1985. $p_{\min}$ стали 4340 —
Bless 1981 (AFWAL-TR-81-4040), пометки `TODO(verify)` для стали 4340 сняты. В `materials.hpp`
остаются `TODO(verify)`: $T_r$ и $c_p$ стали 1006, $\Gamma_0,c_0,S_1$ для 1006 и 7039,
плотность 7039.

## 7. Сводка режимов и ключей

| ключ / поле | значения | эффект |
|---|---|---|
| `jc_heating` | `energy` (по умолч.), `taylor_quinney` | способ расчёта $T$ в JC |
| `jc_chi` | $[0,1]$, 0,9 | доля Тейлора–Куинни |
| `dim`, `axis_sim` (компиляция) | 1, 2, 3; осесимметрия | форма $J_2$, $\dot\varepsilon_{eq}$, окружные члены |
| `JcMaterial::eos` | `Linear`, `Shock` | УРС JC |
| `JcMaterial::p_min` | Па | нижняя граница давления |
| `JcStrength::rate_form`, `rate_source` | `FirstOrder/Exponential`, `Total/Plastic` | только в коде (нет ключа) |
| `JcDamage::form` | `JC`, `Borvik` | форма скоростного члена $\varepsilon_f$ |

## Реализация в коде

| модель | файл | функция |
|---|---|---|
| диспетчер УРС | [int_force.cpp](../../../src/int_force.cpp), [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `urs`, `solid_eos`, `compute_solid_derived` |
| газ, вода | [eos.cpp](../../../src/eos.cpp) | `p_gas`, `p_art_water` |
| Ми–Грюнайзен, УА Тэта | [eos.cpp](../../../src/eos.cpp) | `eos_mi_grunaizen`, `eos_ua_teta` |
| Linear, Shock, JC | [eos.cpp](../../../src/eos.cpp), [eos_types.hpp](../../../include/eos_types.hpp) | `eos_linear`, `eos_shock`, `eos_jc` |
| УРС JH-2 | [eos.cpp](../../../src/eos.cpp) | `urs_cem_JH2` |
| константы материалов | [materials.hpp](../../../include/materials.hpp) | `Materials`, `make_jc_library`, `shear_modulus_by_itype`, `elastic_by_itype` |
| закон Гука + Яуманн | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `finalize_solid`, `effective_shear_modulus` |
| Мизес | [ppochnost.cpp](../../../src/ppochnost.cpp) | `prochnost_mises`, `sigma_i` |
| JH-2 прочность/повреждение/булкинг | [ppochnost.cpp](../../../src/ppochnost.cpp) | `prochnost_cem_JH2` |
| JC (чистые функции) | [johnson_cook.hpp](../../../include/strength/johnson_cook.hpp) | `jc_yield`, `jc_radial_return`, `jc_fracture_strain`, `equivalent_strain_rate`, `j2_of_deviator` |
| JC (частица) | [jc_strength.hpp](../../../include/interactions/jc_strength.hpp) | `apply_jc_strength` |
| JC (начальное состояние) | [jc_state.hpp](../../../include/jc_state.hpp) | `init_jc_state` |
| вызов прочности по `itype` | [postprocessing.hpp](../../../include/interactions/postprocessing.hpp) | `apply_material_strength` |
| ключи `jc_heating`, `jc_chi` | [input.cpp](../../../src/input.cpp), [simulation_config.hpp](../../../include/simulation_config.hpp) | `parse_sph_config` |
