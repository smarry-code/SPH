# Математическое описание расчёта в e2e-тесте BC2

> **Примечание (2026-11-14):** документ описывает конвейер до этапа Ф4-B2: `direct_find`, `density.cpp`, `art_visc.cpp`, `av_vel.cpp`, `virt_part.cpp` и `int_force()` (кроме `urs()`) удалены в a4ce0b9, ссылки на них убраны. Актуальные уравнения — в [docs/04-explanation/models](../docs/04-explanation/models/index.md), конвейер — в [pipeline-overview.md](../docs/04-explanation/architecture/pipeline-overview.md).

Что именно считает `test_e2e_bc2` — от постановки до дискретных уравнений, в том
виде, в каком они реализованы в коде. Каждый блок снабжён ссылкой на реализацию;
при расхождении формулы и кода истиной является код, а документ подлежит правке.

Конфигурация теста: `dim = 2`, `axis_sim = true`, `pa_sph = 2`, `skf = 1`, `sle = 1`,
`nnps = 2`, `CFL = 0.1`, `density_method = 1`, `visc = false`, `visc_artificial = true`,
`average_velocity = true`, `summation_density = false`.

---

## 1. Постановка

Осесимметричное соударение по оси $z$: полиуретановый ударник разгоняет пластину
карбида бора, та бьёт по неподвижной пластине B4C, за которой стоит окно LiF.
Регистрируется история скорости частицы-монитора на тыльной поверхности
неподвижной керамики.

Координаты: $x^0 \equiv r$ (радиус), $x^1 \equiv z$ (ось). Осевая симметрия означает
$\partial/\partial\theta = 0$ и $v_\theta = 0$; расчётная область — полуплоскость $r > 0$.

| Слой | Материал | `itype` | $\varrho_0$, кг/м³ | $v_z(0)$, м/с | $r$, м | $z$, м |
|---|---|---|---|---|---|---|
| 1 | полиуретан 626 | 8 | 626 | 2210 | $(0,\,0.04375]$ | $[0,\,0.006]$ |
| 2 | B4C | 10 | 2510 | 2210 | $(0,\,0.0346]$ | $(Z_1,\,Z_2]$ |
| 3 | B4C | 10 | 2510 | 0 | $(0,\,0.0346]$ | $[Z_2{+}4\Delta x,\ Z_3{+}4\Delta x]$ |
| 4 | LiF | 5 | 2638 | 0 | $(0,\,0.0254]$ | $[Z_3{+}5\Delta x,\ Z_4{+}4\Delta x]$ |

$Z_1 = 0.006$, $Z_2 = 0.009917$, $Z_3 = 0.018950$, $Z_4 = 0.044350$ м.

Начальная дискретизация — регулярная сетка с шагом $\Delta x$; частица с индексом
$(i_r, i_z)$ ставится в $r = (i_r + \tfrac12)\Delta x$, что удерживает все частицы
строго внутри $r > 0$. Начальные данные: $u = 0$, $p = 0$, $\varrho = \varrho_0$,
$v_r = 0$, $\boldsymbol{S} = 0$, $D = 0$.

$$h_i = 2\Delta x,\qquad m_i = \varrho_i\,\Delta x^2 .$$

Масса взята по **плоской** формуле $\varrho\,\Delta x^2$, а не по кольцевой
$2\pi r\,\varrho\,\Delta x^2$ (оба кольцевых варианта присутствуют в исходном коде
закомментированными, [input.cpp:516-517](../src/input.cpp#L516-L517)). Это согласовано
с тем, что осевая симметрия ниже входит только через добавочные члены $\propto 1/r$,
а не через веса частиц.

---

## 2. SPH-аппроксимация

Ядро — кубический сплайн Монагана ([kernel.cpp:21-62](../src/kernel.cpp#L21-L62)),
$q = r_{ij}/h_{ij}$, $h_{ij} = \tfrac12(h_i + h_j)$, при `dim = 2` нормировка
$\sigma = \dfrac{15}{7\pi h_{ij}^2}$:

$$W(q) = \sigma\begin{cases}
\tfrac23 - q^2 + \tfrac12 q^3, & 0 \le q \le 1,\\[2pt]
\tfrac16 (2-q)^3, & 1 < q \le 2,\\[2pt]
0, & q > 2.
\end{cases}$$

Носитель $2h_{ij}$ задаёт и радиус поиска соседей
(direct_find.cpp:34-37).

Длина сглаживания при `sle = 1` — не эволюционное уравнение, а алгебраическая связь
([h_upgrade.cpp:17-25](../src/h_upgrade.cpp#L17-L25)):

$$h_i = 1{,}3\left(\frac{m_i}{\varrho_i}\right)^{1/\mathrm{dim}} .$$

---

## 3. Уравнение неразрывности

Плоская часть при `density_method = 1`
(density.cpp:126-130):

$$\left.\frac{d\varrho_i}{dt}\right|_{\text{пар}} = \sum_j m_j\,(\mathbf v_i - \mathbf v_j)\cdot\nabla_i W_{ij},$$

что аппроксимирует $-\varrho\,(\partial_r v_r + \partial_z v_z)$.

### Осесимметричная поправка

В цилиндрических координатах

$$\nabla\cdot\mathbf v = \frac{\partial v_r}{\partial r} + \frac{v_r}{r} + \frac{\partial v_z}{\partial z},$$

поэтому парная сумма недобирает член $v_r/r$. Полное уравнение:

$$\frac{d\varrho_i}{dt} = \sum_j m_j\,(\mathbf v_i - \mathbf v_j)\cdot\nabla_i W_{ij}\;-\;\varrho_i\,\frac{v_{r,i}}{r_i}.$$

Второе слагаемое относится к **частице**, а не к паре, и вычисляется в отдельном
пер-частичном цикле после парного (density.cpp:235-247).

> **История.** До правки этот член стоял внутри цикла по парам и не имел множителя
> $\varrho_i$: на каждую пару добавлялось $-v_{r,i}/r_i$. Частица с $N_i$ соседями
> получала $-N_i v_{r,i}/r_i$ вместо $-\varrho_i v_{r,i}/r_i$, то есть при
> $N_i \approx 50$ и $\varrho = 2510$ вклад был занижен примерно в 50 раз, а его
> размерность (1/с вместо кг·м⁻³·с⁻¹) не отвечала левой части. Та же поправка в
> уравнении энергии всегда стояла правильно — пер-частично
> ([int_force.cpp:620-624](../src/int_force.cpp#L620-L624)).

---

## 4. Уравнение движения

Полное напряжение раскладывается на шаровую и девиаторную части:

$$\sigma^{\alpha\beta} = -p\,\delta^{\alpha\beta} + S^{\alpha\beta}.$$

При `pa_sph = 2` ([int_force.cpp:338](../src/int_force.cpp#L338), `rhoij` $= 1/(\varrho_i\varrho_j)$)
ускорение частицы:

$$\frac{dv_i^\alpha}{dt} = \sum_j \frac{m_j}{\varrho_i\varrho_j}\Bigl[-(p_i + p_j)\,\delta^{\alpha\beta} + (S_i^{\alpha\beta} + S_j^{\alpha\beta})\Bigr]\frac{\partial W_{ij}}{\partial x_i^\beta}\;+\;a^\alpha_{\text{сим}}\;+\;a^\alpha_{\text{иск.в}}\;+\;a^\alpha_{\text{внеш}} .$$

Осесимметричные добавки ([int_force.cpp:417-430](../src/int_force.cpp#L417-L430)) —
вклад окружного напряжения $S^{\theta\theta}$, в коде хранимого как `szz`:

$$a^r_{\text{сим},i} = \sum_j \frac{m_j}{\varrho_i\varrho_j}\,\frac{2S^{rr}_j + S^{\theta\theta}_j}{r_j}\,W_{ij},
\qquad
a^z_{\text{сим},i} = \sum_j \frac{m_j}{\varrho_i\varrho_j}\,\frac{S^{rz}_j}{r_j}\,W_{ij}.$$

### Искусственная вязкость

Монагановская форма (art_visc.cpp:42-62), включается
только при сближении $\mathbf v_{ij}\cdot\mathbf x_{ij} < 0$:

$$\mu_{ij} = \frac{h_{ij}\,\mathbf v_{ij}\cdot\mathbf x_{ij}}{|\mathbf x_{ij}|^2 + \eta_q^2 h_{ij}^2},
\qquad
\Pi_{ij} = \frac{\beta\,\mu_{ij}^2 - \alpha\,\bar c_{ij}\,\mu_{ij}}{\bar\varrho_{ij}},$$

$\bar c_{ij}$, $\bar\varrho_{ij}$ — средние арифметические по паре,
$\alpha = \beta = 1$, $\eta_q = 0.01$. Вклады:

$$\frac{dv_i^\alpha}{dt}\mathrel{+}= -\sum_j m_j \Pi_{ij}\frac{\partial W_{ij}}{\partial x_i^\alpha},
\qquad
\frac{de_i}{dt}\mathrel{+}= \tfrac12\sum_j m_j \Pi_{ij}\,\mathbf v_{ij}\cdot\nabla_i W_{ij}.$$

### Внешние силы

При `contact_force_type = 0` (базовый сценарий) контактных сил нет; действует только
пристеночное отталкивание и, при `self_gravity`, тяжесть — обе в текущей постановке
неактивны. Контактные сценарии §6 ТЗ включают штрафные силы по ADR-007.

---

## 5. Деформации, вращение и девиатор напряжений

Скорости деформации и вращения ([int_force.cpp:280-307](../src/int_force.cpp#L280-L307)):

$$\dot\varepsilon_i^{\alpha\beta} = \frac12\sum_j \frac{m_j}{\varrho_j}\left(v_{ji}^\alpha\frac{\partial W_{ij}}{\partial x_i^\beta} + v_{ji}^\beta\frac{\partial W_{ij}}{\partial x_i^\alpha}\right),
\qquad
\dot\omega_i^{\alpha\beta} = \frac12\sum_j \frac{m_j}{\varrho_j}\left(v_{ji}^\alpha\frac{\partial W_{ij}}{\partial x_i^\beta} - v_{ji}^\beta\frac{\partial W_{ij}}{\partial x_i^\alpha}\right).$$

Гипоупругое уравнение Джаумана для девиатора
([int_force.cpp:656-670](../src/int_force.cpp#L656-L670)); шаровая часть выражена через
$\dot\varrho/\varrho$, а не через след тензора деформаций:

$$\frac{dS^{\alpha\beta}}{dt} = 2G\left(\dot\varepsilon^{\alpha\beta} - \frac{\delta^{\alpha\beta}}{3}\,\frac{\dot\varrho}{\varrho}\cdot(-1)\right) + S^{\alpha\gamma}\dot\omega^{\beta\gamma} + S^{\gamma\beta}\dot\omega^{\alpha\gamma},$$

в коде — покомпонентно, например
$\dot S^{rr} = 2G\bigl(\dot\varepsilon^{rr} + \tfrac{\dot\varrho}{3\varrho}\bigr) + 2\bigl(S^{rz}\dot\omega^{rz} + \dots\bigr)$.

Модуль сдвига $G$ берётся из реестра материалов
([materials.hpp:47-60](../include/materials.hpp#L47-L60)): для `itype` 3–9 — из `G`,
для 10–11 — из `const_JH_2[·][1]`.

---

## 6. Уравнения состояния

Выбор по `itype` ([int_force.cpp:14-35](../src/int_force.cpp#L14-L35)).

**LiF и полиуретаны (`itype` 5–9), УРС «ударная адиабата Тэта»**
([eos.cpp:55-75](../src/eos.cpp#L55-L75)), $w = \varrho/\varrho_0$:

$$p = A\left(w^{Z} - 1\right),\qquad c^2 = \frac{AZ}{\varrho_0}\,w^{Z-1}.$$

**Тантал, алюминий (`itype` 3–4), Ми–Грюнайзен**
([eos.cpp:29-54](../src/eos.cpp#L29-L54)): в кейсе BC2 не используются.

**B4C (`itype` 10), полиномиальный УРС модели JH-2**
([eos.cpp:76-110](../src/eos.cpp#L76-L110)), $\mu = \varrho/\varrho_0 - 1$:

$$p = \begin{cases}
K_1\mu + K_2\mu^2 + K_3\mu^3 + \Delta p, & \mu \ge 0,\\[2pt]
\max(K_1\mu,\ -T), & \mu < 0,
\end{cases}
\qquad
c = \sqrt{\frac{K_1 + 2K_2\mu + 3K_3\mu^2}{\varrho_0}} .$$

$\Delta p$ — приращение давления от объёмного раздува при накоплении повреждения
(bulking, §7). $T$ — предел прочности на отрыв.

---

## 7. Модели прочности

### Мизес (металлы и полимеры)

([ppochnost.cpp:31-54](../src/ppochnost.cpp#L31-L54)) Радиальный возврат: если
$J_2 > \sigma_Y^2/3$, весь девиатор умножается на
$\sqrt{\sigma_Y^2/3}\,/\sqrt{J_2}$. В осесимметричном случае
([ppochnost.cpp:17-20](../src/ppochnost.cpp#L17-L20)):

$$J_2 = S_{rr}^2 + S_{rr}S_{zz} + S_{zz}^2 + S_{rz}^2 .$$

### JH-2 (керамика B4C)

([ppochnost.cpp:55-121](../src/ppochnost.cpp#L55-L121)) Нормированные давление и
прочность, $P^* = p/P_{HEL}$:

$$\sigma^*_{i} = A\,(P^* + T^*)^N \quad\text{(неповреждённая)},\qquad
\sigma^*_{f} = \min\bigl(B\,(P^*)^M,\ \sigma^*_{f,\max}\bigr) \quad\text{(разрушенная)},$$

текущая прочность с учётом повреждения $D \in [0,1]$ и скорости деформации:

$$\sigma^* = \bigl[\sigma^*_i - D\,(\sigma^*_i - \sigma^*_f)\bigr]\,\bigl(1 + C\ln\dot\varepsilon_{\text{экв}}\bigr),
\qquad
\dot\varepsilon_{\text{экв}} = \frac{\sqrt2}{3}\sqrt{(\dot\varepsilon_{rr}-\dot\varepsilon_{\theta\theta})^2 + \dots + 6\dot\varepsilon_{rz}^2}.$$

Предел текучести $\sigma_t = \sigma^*\,\sigma_{HEL}$, возврат девиатора — как у Мизеса.

**Bulking.** При $\mu \ge 0$ и падении упругой энергии девиатора $\Delta U > 0$
приращение давления находится из
$\Delta p_{\text{нов}} = -K_1\mu + \sqrt{(K_1\mu + \Delta p)^2 + 2\beta K_1 \Delta U}$
([ppochnost.cpp:122-138](../src/ppochnost.cpp#L122-L138)), $\beta = 1$.

**Накопление повреждения** идёт по приращению пластической деформации, вычисляемому
как разность полного и упругого приращения девиатора
([ppochnost.cpp:139-155](../src/ppochnost.cpp#L139-L155)).

---

## 8. Уравнение энергии

([int_force.cpp:619-630](../src/int_force.cpp#L619-L630)) В осесимметричном случае:

$$\frac{de_i}{dt} = \underbrace{\dot e_i^{\,p}}_{\text{давление}} + \frac{1}{\varrho_i}\Bigl[S^{rr}\bigl(\dot\varepsilon^{rr} - \tfrac{v_r}{r}\bigr) + S^{\theta\theta}\bigl(\dot\varepsilon^{\theta\theta} - \tfrac{v_r}{r}\bigr) + 2S^{rz}\dot\varepsilon^{rz}\Bigr],$$

где член $v_r/r$ входит **один раз на частицу** — эталонный образец того, как должна
была выглядеть и поправка в уравнении неразрывности (§3).

---

## 9. Сглаживание скоростей (XSPH)

При `average_velocity = true` (av_vel.cpp:24-36):

$$\delta\mathbf v_i = -\varepsilon\sum_j \frac{2m_j}{\varrho_i + \varrho_j}\,(\mathbf v_i - \mathbf v_j)\,W_{ij},\qquad \varepsilon = 0.25,$$

добавляется к скорости при переносе частиц.

---

## 10. Интегрирование по времени

Схема «чехарда» с предиктором на полушаге
([time_integration.cpp:100-205](../src/time_integration.cpp#L100-L205)). Первый шаг
особый (нет предыстории), далее для каждого шага:

$$\varrho^{n+1/2} = \varrho^{n-1/2} + \tfrac{\Delta t}{2}\dot\varrho^{\,n},\quad
\mathbf v^{n+1/2} = \mathbf v^{n-1/2} + \tfrac{\Delta t}{2}\dot{\mathbf v}^{\,n},\quad
\boldsymbol S^{n+1/2} = \boldsymbol S^{n-1/2} + \tfrac{\Delta t}{2}\dot{\boldsymbol S}^{\,n},$$

затем пересчёт правых частей и полный шаг с добавлением XSPH-поправки:

$$\mathbf v^{n+1} = \mathbf v^{n} + \Delta t\,\dot{\mathbf v}^{\,n+1/2} + \delta\mathbf v,
\qquad
\mathbf x^{n+1} = \mathbf x^{n} + \Delta t\,\mathbf v^{n+1}.$$

### Шаг по времени

([time_integration.cpp:13-27](../src/time_integration.cpp#L13-L27)) Условие Куранта
по всем частицам:

$$\Delta t = \min_i \frac{\mathrm{CFL}\cdot h_i}{c_i + \lVert\mathbf v_i\rVert},\qquad \mathrm{CFL} = 0.1,$$

и, если включён контакт, дополнительно ограничивается $\Delta t \le \Delta t_{\text{контакт}}$
из ADR-007.

---

## 11. Граничные условия

- **Ось $r = 0$.** Частицы стоят в $r = (i_r + \tfrac12)\Delta x$, поэтому $r > 0$ всегда
  и деления на нуль в членах $\propto 1/r$ не возникает. Зеркальные частицы, которые
  восполняли бы нехватку соседей у оси, **не генерируются**: вызов `virt_part` в
  [single_step.cpp:103-106](../src/single_step.cpp#L103-L106) закомментирован (и, как
  записан, не скомпилировался бы — четыре аргумента против трёх в объявлении). Осевой
  столбец остаётся с неполным носителем ядра.
- **Свободные поверхности** специально не обрабатываются: отсутствие соседей снаружи
  само даёт нулевое напряжение в среднем, без коррекции ядра.
- **Границы раздела слоёв** при `contact_force_type = 0` взаимодействуют только через
  общую SPH-сумму: отдельные тела для расчёта не различаются, соседство определяется
  исключительно расстоянием.

---

## 12. Что сверяет тест

Величина, попадающая в эталон, — вторая колонка `log.txt`
([time_integration.cpp:304-305](../src/time_integration.cpp#L304-L305)):

$$y(t) = v_r\bigl(t,\ i_{\text{мон}}\bigr),$$

то есть **радиальная** компонента скорости частицы-монитора (не осевая). Монитор —
частица на тыльной поверхности неподвижной керамики, $r \approx 12{,}5$ мм.

Кривые расчёта и эталона снимаются на разных временны́х сетках (шаг адаптивный),
поэтому обе линейно интерполируются на общую равномерную сетку из $N = 200$ узлов, и
сравниваются две метрики:

$$\text{rel}_{L_2} = \frac{\lVert y_{\text{расч}} - y_{\text{эт}}\rVert_2}{\lVert y_{\text{эт}}\rVert_2},
\qquad
\text{max\_abs\_dev} = \max_k \bigl|y_{\text{расч}}(t_k) - y_{\text{эт}}(t_k)\bigr|.$$

Допуски регрессионного режима — $\text{rel}_{L_2} < 10^{-3}$ и
$\text{max\_abs\_dev} < 1$ м/с.

> **Ограничение.** Эти допуски достижимы только при воспроизводимом расчёте. При
> включённом OpenMP порядок атомарных операций в парных суммах не фиксирован, из-за
> чего два прогона одного бинарника расходятся; на длинном интервале расхождение
> превышает допуск. До устранения недетерминизма полный режим пригоден как проверка
> устойчивости (отсутствие NaN, завершаемость, некатастрофический шаг), но не как
> побитовая регрессия.

---

## Обозначения

| Символ | Смысл | Поле в коде |
|---|---|---|
| $\varrho$, $\varrho_0$ | текущая и начальная плотность | `rho`, `rho_0` |
| $p$ | давление | `p` |
| $c$ | скорость звука | `c` |
| $e$ | удельная внутренняя энергия | `u` |
| $S^{\alpha\beta}$ | девиатор напряжений | `sxx`, `syy`, `szz`, `sxy`, `syz`, `sxz` |
| $\dot\varepsilon^{\alpha\beta}$ | скорость деформации | `exx`, … |
| $\dot\omega^{\alpha\beta}$ | скорость вращения | `rxy`, `rxz`, `ryz` |
| $D$ | повреждение JH-2, $[0,1]$ | `D` |
| $\Delta p$ | давление раздува (bulking) | `dp` |
| $h$ | длина сглаживания | `hsml` |
| $\theta\theta$ | окружная компонента | хранится в `szz` / `ezz` |

Индекс $z$ в коде — это `x[1]`; окружное направление явных координат не имеет и
входит только через члены $\propto 1/r$ и компоненту `szz`.
