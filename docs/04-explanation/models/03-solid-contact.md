# Контакт между SPH-телами, внешние силы и граничные условия

Раздел описывает три группы моделей, которые не входят в континуальные SPH-суммы
твёрдого тела: (1) штрафной контакт «частица–частица» между разными телами
(`contact_force_type` = 1, 2, 3); (2) кинематическое закрепление (Clamp) частиц
Solid и бусин нитей; (3) внешние силы и граничные условия, которые есть в коде
(и какие из них реально участвуют в расчёте). Контакты нитей (нить–нить,
нить–Solid) здесь не рассматриваются, см. [05-fiber-contact.md](05-fiber-contact.md).

Обозначения: $i,j$ — частицы разных тел, $\mathbf x_{ij}=\mathbf x_i-\mathbf x_j$,
$r_{ij}=|\mathbf x_{ij}|$, $\mathbf v_{ij}=\mathbf v_i-\mathbf v_j$,
$\mathbf n=\mathbf x_{ij}/r_{ij}$ — нормаль (от $j$ к $i$). Проникание ($\delta$ в
общих обозначениях) в коде называется $p$; $p>0$ при перекрытии (не путать с
давлением).

!!! note "Какой код продовый"
    Продовый путь контакта Solid–Solid — ядро
    [solid_contact.hpp](../../../include/interactions/solid_contact.hpp)
    (`compute_solid_contact_contributions`), вызываемое из
    [single_step.cpp](../../../src/single_step.cpp). Старый оркестратор
    [ext_force.cpp](../../../src/ext_force.cpp) (`ext_force`) в `src/` нигде не
    вызывается; он остался как тестовый стенд (`tests/test_contact.cpp`,
    `test_e2e_bc2.cpp`, `test_e2e_bc6.cpp`). Формулы силы у обоих путей общие
    ([contact_force.hpp](../../../include/contact_force.hpp)); различия путей
    указаны в п. 1.7.

---

## 1. Штрафной контакт между телами

### 1.1. Физический смысл

SPH-суммы по соседям связывают все частицы в носителе ядра, поэтому границу двух
тел SPH «не видит»: напряжение через неё передаётся обычными уравнениями
сохранения. Штрафной контакт добавляет к уравнению движения силу отталкивания
между дискретными частицами разных тел, как между телами конечного размера:
частица $k$ считается шаром (диском при `dim == 2`) радиуса $R_k$, и при
перекрытии шаров на пару действует сила вдоль линии центров. Сила только
отталкивающая; трения и тангенциальной компоненты нет; нормаль — по линии центров
пары (нормали к поверхности тела не строятся).

Источник постановки — `notes/sources/sph-contact-article.md`, разд. 2 («пинбол» Белицко, Свигл, импульсный
контакт); решения по отступлениям — [ADR-007](../decisions/adr-007-contact-force-types.md),
отчёт [notes/contact/11-contact-report.md](../../../notes/contact/11-contact-report.md).

### 1.2. Контактные пары

Пара $(i,j)$ — контактная, если выполнены все условия:

1. `contact_force_type != 0`;
2. $i<j$, обе частицы реальные (зеркальные/ghost-частицы осевой симметрии из
   контакта исключены: `contact_pairs` строится только из `real_pairs`, а в
   поиске пары с $i\ge n_\text{real}$ отбрасываются);
3. частицы принадлежат разным телам: `body[i] != body[j]`;
4. пара попала в SPH-поиск соседей (критерий $r<\text{scale}_k\,(h_i+h_j)/2$,
   `scale_k` = 2 для кубического сплайна, 3 для остальных ядер) — отдельного
   контактного поиска нет, список `contact_pairs` — подмножество `real_pairs`,
   построенное в том же проходе
   ([solid_sph.hpp](../../../include/interactions/solid_sph.hpp),
   `build_solid_sph_search`, структура `SolidContactPair`);
5. имеется перекрытие $p>0$ и $r>0$ (см. 1.3).

Идентификатор тела `body[]` для k-файла равен PID (`k_geometry.cpp`), для
встроенных генераторов и `particles.dat` — `itype` (`input.cpp`). Следовательно,
на встроенных кейсах «разные тела» = «разные материалы»; в k-файле два тела из
одного материала с разными PID контактируют, а частицы одного PID — нет.
Внутри `compute_solid_contact_contributions` стоит повторная защитная проверка
`body[i] == body[j]` (ранний выход).

Ключ `sph_same_body_only` (bool, по умолчанию `false`; принимает `0/1/true/false`)
отвергает SPH-пары реальных частиц разных тел: примитив SPH для такой пары не
строится, пара остаётся только в списке контакта. Тогда тела взаимодействуют
исключительно через штрафной контакт; при `sph_same_body_only = 1` и
`contact_force_type = 0` `build_solid_sph_search` бросает исключение. Пары
«реальная — ghost» остаются в SPH-суммах (условие отбраковки требует $j<n_\text{real}$).

### 1.3. Радиус, проникание, скорость сближения

Радиус частицы определяется из её объёма (при `dim == 2` — площади меридионального
сечения: масса в осесимметричном генераторе равна $\rho\,\Delta x^2$), а не из $h$:

$$
R_k=\frac12\left(\frac{m_k}{\max(\rho_k,\rho_{0,k})}\right)^{1/d},
\qquad R_k=0\ \text{при}\ \rho_k\le 0.
$$

Потолок по начальной плотности $\rho_0$ (`sys.rho_0`) означает, что сжатие радиус
уменьшает, а разрежение его не увеличивает (иначе расходящиеся частицы
«проникают», см. комментарий в `ext_force.cpp::particle_radius`). Проникание и
скорость сближения:

$$
p=(R_i+R_j)-r_{ij},\qquad
\dot p=-\frac{\mathbf x_{ij}\cdot\mathbf v_{ij}}{r_{ij}}\quad(\dot p>0\ \text{при сближении}).
$$

Для соприкасающихся в покое частиц решётки $R_i+R_j=\Delta x$ и $p=0$: контакт
включается только при реальном сближении. Дополнительно в пару передаются:

| Величина | Формула в коде |
|---|---|
| приведённая масса | $m_{ij}=m_im_j/(m_i+m_j)$ |
| скорость звука пары | $c_\text{pair}=\min(c_i,c_j)$ ($c$ — из `SolidDerived`, т.е. EOS текущего RHS) |
| эффективный модуль сдвига | $G_\text{eff}=G_iG_j/(G_i+G_j)$, если $G_i>0$ и $G_j>0$, иначе $0$ («G неизвестен») |
| эффективный радиус | $R_\text{eff}=R_iR_j/(R_i+R_j)$ |
| опорный шаг пары | $\Delta t_\text{ref}=\min_{k\in\{i,j\}}\dfrac{\mathrm{CFL}\,h_k}{c_k+\lVert\mathbf v_k\rVert}$ (при нулевом знаменателе $10^{30}$) |

$G_k$ берётся из `shear_modulus_by_itype(mats, itype)`. Жёсткость типов 1–2
строится от $\Delta t_\text{ref}$ (локальный CFL пары), а не от применяемого шага:
так $\Delta t_\text{applied}\le\Delta t_\text{ref}$, и петля «жёсткость → шаг → жёсткость»
не возникает. Поле `PairGeom::v_rel` (полная относительная скорость) заполняется,
но формулами не используется.

### 1.4. Типы силы

`contact::penalty_force(type, g, pr, kp, dt_applied, gamma)` — диспетчер;
результат `Result{fn, dt_limit}`: модуль нормальной силы $F_n\ge0$ и ограничитель
шага. Тип 0 и любое иное значение дают $F_n=0$ без ограничителя. Параметр
`kp` = `contact_penalty_scale` ($K_p$), `gamma` = `contact_impulse_gamma` ($\gamma$).

Сила прикладывается вдоль нормали:

$$
\mathbf F_i=+F_n\,\mathbf n,\qquad \mathbf F_j=-F_n\,\mathbf n,\qquad
\Delta\dot{\mathbf v}_i=\frac{\mathbf F_i}{m_i},\quad \Delta\dot{\mathbf v}_j=\frac{\mathbf F_j}{m_j}.
$$

Равные и противоположные силы на линии центров: импульс и момент импульса пары
сохраняются точно (с массами $m_k$, как они заведены в расчёт; в осесимметрии
это масса на единицу меридионального сечения, дополнительных радиальных/окружных
членов контакт не вводит).

#### Тип 1 — «пинбол» Белицко (`pinball_force`)

Условие: $p>0$, $m_{ij}>0$, $\Delta t_\text{ref}>0$, $K_p>0$. Две ветви:

$$
F_1=\begin{cases}m_{ij}\,\dot p/\Delta t_\text{ref},&\dot p>0\\0,&\dot p\le0\end{cases},\qquad
F_2=G_\text{eff}\sqrt{R_\text{eff}}\;p^{3/2}\ \ (\text{Герц}),
$$

$$
F_n=K_p\begin{cases}\min(F_1,F_2),&G_\text{eff}>0\ \text{и}\ R_\text{eff}>0\\F_1,&\text{иначе (откат: }F_2\text{ исключается из }\min).\end{cases}
$$

Множитель $p$ из буквальной формулы первоисточника опущен (иначе размерность
Н·м; решение D1 в `notes/contact/02-decisions.md`). $F_1$ — вязкая («гасящая») сила
(демпфер с $c_d=K_pm_{ij}/\Delta t_\text{ref}$), $F_2$ — герцевский потолок. Так как
$F_1=0$ при $\dot p\le0$, сила на разгрузке и в покое равна нулю — контакт активен
только при сближении.

Ограничитель шага. Если выбрана ветка Герца ($G_\text{eff}>0$ и $F_2<F_1$):
$k_\text{eff}=1.5F_n/p$, $\Delta t_\text{lim}=0.8/\sqrt{k_\text{eff}/m_{ij}}$.
Иначе (вязкая): $\Delta t_\text{lim}=0.8\cdot2\Delta t_\text{ref}/K_p$
(предел устойчивости явной схемы для демпфера, $\Delta t<2m_{ij}/c_d$).
Константа запаса `kOscSafety = 0.8`.

#### Тип 2 — Свигл (`swegle_force`)

Условие: $p>0$, $\dot p>0$, $m_{ij}>0$, $\Delta t_\text{ref}>0$, $r>0$, $K_p>0$.

$$
F_1=\frac{m_{ij}\,p}{\Delta t_\text{ref}^2},\qquad
F_2=\frac{2m_{ij}\,c_\text{pair}\,p}{\Delta t_\text{ref}\,r_{ij}},\qquad
F_\text{cap}=\frac{m_{ij}\,\dot p}{\Delta t_\text{ref}},
$$

$$
F_n=K_p\cdot\begin{cases}F_\text{cap},&\max(F_1,F_2)>F_\text{cap}\\\max(F_1,F_2),&\text{иначе.}\end{cases}
$$

Потолок $F_\text{cap}$ — добавка проекта (в формуле статьи его нет): формула
выведена для $p$, накопленного за шаг ($p=\dot p\,\Delta t_\text{ref}$), тогда
$F_1=F_\text{cap}$ и потолок неактивен; он срабатывает, когда это допущение
нарушено (на BC2 без него $p$ у расходящихся частиц давало силу, завышенную в
$\sim6\cdot10^3$ раз). Как следствие, при $\dot p\le0$ сила нулевая.
Ограничитель: на потолке — как у вязкой ветки типа 1,
$\Delta t_\text{lim}=0.8\cdot2\Delta t_\text{ref}/K_p$; иначе
$k_\text{eff}=F_n/p$ (сила линейна по $p$), $\Delta t_\text{lim}=0.8/\sqrt{k_\text{eff}/m_{ij}}$.
Во второй формуле статьи опечатка (помечена $F_1$, по смыслу $F_2$) — в коде
использована как $F_2$.

#### Тип 3 — импульсный (`impulse_force`)

Условие: $p>0$, $\dot p>0$, $m_{ij}>0$, $\Delta t_\text{applied}>0$, $\gamma>0$.

$$
F_n=\gamma\,\frac{m_{ij}\,\dot p}{\Delta t_\text{applied}},\qquad
\Delta t_\text{lim}=\frac{0.8\cdot2\,\Delta t_\text{applied}}{\gamma}.
$$

За шаг снимается доля $\gamma$ нормальной скорости сближения пары ($\gamma=1$ —
полное гашение за шаг). Отступления от буквальной (2.10): вместо модуля полной
относительной скорости берётся нормальная составляющая $\dot p$ (иначе
тангенциальное скольжение превращается в отталкивание, а сила разрывна в точке
включения), и введена доля гашения $\gamma$. $K_p$ тип 3 не использует.
Приращение скорости за шаг $\Delta v_i=\gamma\,(m_{ij}/m_i)\,\dot p$ от шага не
зависит, поэтому по шагу тип 3 не сходится ни при каком $\gamma$ (комментарий
в `contact_force.hpp`). Рабочий диапазон $\gamma\in(0;1.6]$: при $\gamma>1.6$
ограничитель становится строже применяемого шага ($1.6/\gamma<1$), а жёсткость
сама зависит от шага (петля); разбор конфигурации предупреждает, но значение
принимает.

Здесь $\Delta t_\text{applied}$ — значение `dt`, переданное в `single_step`
(в цикле `time_integration.cpp` это `dt` с конца предыдущей итерации, до
пересчёта после RHS), а не $\Delta t_\text{ref}$ пары.

### 1.5. Ограничитель шага от контакта

Диагностический проход `compute_solid_contact_diagnostics` (отдельный
`parallel_reduce` по `contact_pairs` с повторным вычислением геометрии,
выполняется каждый шаг при `contact_force_type != 0`) собирает:
`contact_pair_count` — число пар с $p>0$ и $r>0$, `contact_max_delta` — $\max p$,
`min_dt_contact` — $\min \Delta t_\text{lim}$ (поля `WorkBuffers`,
[buffers.hpp](../../../include/buffers.hpp)). Пары с $F_n=0$ вносят
$\Delta t_\text{lim}=10^{30}$ (`kNoLimit`), то есть не ограничивают, но считаются
в `contact_pair_count`. Итоговый шаг в
[time_integration.cpp](../../../src/time_integration.cpp):

$$
\Delta t=\min\bigl(\Delta t_\text{CFL},\ \Delta t_\text{contact},\ \Delta t_\text{fiber}\bigr)\cdot s,
$$

где $s$ = `dt_scale` (применяется только если `dt_scale < 1`). В блок печати
прогресса при `contact_force_type != 0` добавляется строка «Шаг ограничен
КОНТАКТОМ/CFL» с `contact_pair_count` и `contact_max_delta`. При $\gamma\le1.6$
ограничитель типа 3 не строже $\Delta t_\text{applied}$ и является страховкой.

### 1.6. Параметры

| Символ | Ключ / поле | По умолчанию | Единицы | Замечание |
|---|---|---|---|---|
| тип силы | `contact_force_type` (`cfg.contact_force_type`) | `0` | — | 0 выкл., 1 пинбол, 2 Свигл, 3 импульсный; вне 0..3 — откат к 0 с сообщением (`input.cpp`) |
| $K_p$ | `contact_penalty_scale` | `1.0` | — | типы 1–2; значение не проверяется, $K_p\le0$ даёт $F_n=0$ |
| $\gamma$ | `contact_impulse_gamma` | `1.0` | — | тип 3; $\le0$ — тип 3 фактически выключен, $>1.6$ — предупреждение |
| — | `sph_same_body_only` | `false` | — | изоляция SPH-сумм по телам (п. 1.2) |
| CFL | `cfg.CFL` | см. конфиг | — | входит в $\Delta t_\text{ref}$ |
| $G_k$, $\rho_{0}$ | материал (`Materials`), `sys.rho_0` | из материала | Па, кг/м³ | $G=0$ отключает ветку Герца типа 1 |

Типы 4 (Герц + демпфер) и «0 = Герц» относятся только к контактам нитей
(`fiber_*_contact_type`), к Solid–Solid не применяются.

### 1.7. Сводка свойств и различия двух путей

- Импульс и момент импульса пары сохраняются (силы равны и противоположны, лежат
  на линии центров). Энергия не сохраняется: потенциальной энергии контакта нет,
  вязкие ветви (тип 1 на $F_1$, тип 3, потолок типа 2) диссипативны; потери
  выводятся каналом `D_contact_ss` лога энергии (сборка с `SPH_ENERGY_LOG`).
- Условие активации: $p>0$ и $\dot p>0$ (для типа 1 при известном $G$ сила
  $\min(F_1,F_2)$ тоже равна нулю при $\dot p\le0$) — притяжения нет.
- Вклад контакта записывается как `dv[dim]` в `SolidDelta` и сводится общим
  `reduce_contributions_into` вместе с Solid-SPH; затем попадает в `dvx` через
  `finalize_solid`. Плотность, энергию и XSPH-поправку контакт напрямую не меняет.
- Детерминизм: разбиение `exec.det_segments()` не зависит от числа потоков;
  `interaction_type` вкладов = 11 (`kSolidContactInteractionType`).

| Аспект | Продовый путь (`solid_contact.hpp`) | Старый (`ext_force.cpp`) |
|---|---|---|
| критерий пары | `body[i] != body[j]` | `itype[i] != itype[j]`, оба `> 0` |
| ghost-частицы | исключены (список из `real_pairs`) | исключены условием $i,j<n_\text{real}$ |
| скорость звука | `derived.c` (EOS этого RHS) | `sys.c` |
| диагностика | отдельный проход `parallel_reduce` | накопление в том же цикле |

---

## 2. Закрепление (Clamp)

### 2.1. Постановка

Clamp — кинематическое условие Дирихле для частиц Solid и бусин нитей:

$$
x_{i}^{\alpha}(t)=x_{i}^{\alpha}(0),\qquad v_{i}^{\alpha}(t)=0\quad\text{для отмеченных компонент }\alpha .
$$

Это не тип взаимодействия (поиска и ядра нет), а стадия-ограничение над уже
сведённым RHS ([clamp.hpp](../../../include/interactions/clamp.hpp)). Условие
реализуется проекцией ускорений: у закреплённой компоненты $a_{i}^{\alpha}\leftarrow0$
(у Solid дополнительно обнуляется XSPH-поправка скорости `av`). Так как начальная
скорость закреплённых компонент обязана быть нулевой, оба решателя
(предиктор–корректор и Эйлер) дают $v=0$, $x\mathrel{+}=\Delta t\cdot0$ побитово, и
`advance()` не меняется.

### 2.2. Как задаётся закрепление

Маска — `std::uint8_t` на частицу Solid (размер $n_\text{total}$) или бусину
(размер — число бусин): бит $\alpha$ = 1 — компонента $\alpha$ закреплена;
`kClampAll` $=(1\ll d)-1$ — закреплена целиком; биты $\ge d$ запрещены. Пустая
маска — закреплений нет, пустой `ClampSet` — стадия не выполняется вообще
(результат побитово совпадает с расчётом без Clamp). Маски хранятся в
`fiber::Model::clamps` (`ClampSet{solid_mask, fiber_mask}`). Источники:

| Источник | Что закрепляет | Где |
|---|---|---|
| `CLAMP_BOX xmin … zmax [mask]` в `.fib` | бусины в замкнутом боксе $lo_\alpha\le x^\alpha\le hi_\alpha$; `mask` — буквы `x,y,z` без повторов, по умолчанию все; несколько боксов объединяются по «или» | `fiber_input.cpp`, `interactions::mark_clamp_box` |
| `*BOUNDARY_SPC_NODE`, `*BOUNDARY_SPC_SET` в k-файле | узел с DOF — SPH-частица Solid и/или бусина (по NID); компоненты $\ge d$ пропускаются с предупреждением; `CID != 0` — ошибка | `k_geometry.cpp`, `apply_spc_clamps` |
| код сцены | любая маска напрямую | `model.clamps` |

Проверки при старте (`validate_clamps`): размеры масок, биты $\ge d$ и нулевая
начальная скорость закреплённых компонент (иначе исключение с индексом
частицы). Печатается число закреплённых частиц/бусин и $\max h$ среди
закреплённых Solid (`ClampReport`).

### 2.3. Точки вставки в шаг и реакция

Все точки выполняются только при непустом `ClampSet`
([time_integration.cpp](../../../src/time_integration.cpp)):

1. После RHS (Solid-SPH, контакты, нити полностью сложены), до вычисления $\Delta t$ —
   `apply_clamp_projection` для Solid (`dvx` и `av`) и для нитей (`dv`);
   реакция опоры накапливается до обнуления: $R^{\alpha}=-\sum_i m_i\,a_i^{\alpha}$ по
   закреплённым компонентам. При сборке с `SPH_ENERGY_LOG` те же компоненты
   обнуляются и в отдельных буферах ускорений контакта (иначе в канале
   контакта появилась бы работа неподвижной опоры).
2. В стадии трения нитей, после подсчёта счётчиков mass splitting, —
   `apply_clamp_friction_counts`: $c_k=c_s=0$ у **полностью** закреплённых
   (бесконечная масса в проекции трения).
3. После сложения вкладов трения, до `advance()` — повторная проекция `dvx`, `dv`
   (без `av`); реакция прибавляется к накопленной в п. 1.

Итоговая реакция — в `ClampReport` (`solid_force[dim]`, `fiber_force[dim]`, а также сумма и максимум
модулей реакций `solid_force_abs`, `fiber_force_abs`, `solid_force_max`, `fiber_force_max`).
Реакция считается детерминированной двухуровневой суммой с блоками по
`kClampReduceBlock = 256` (не `parallel_reduce`), списки короче
`kClampSerialThreshold = 8192` обходятся последовательно; результат не зависит
от числа потоков. `clamp_project_item` — операция над одной частицей без
исполнителя и аллокаций (задел под GPU). В логе энергии работа закрепления
учитывается источником `S_u_clamp` (накопленный, $>0$).

Проекция выполняется после контакта, поэтому штрафной контакт на закреплённую
частицу входит в реакцию опоры, а само ускорение закреплённой компоненты
обнуляется.

---

## 3. Прочие внешние силы и граничные условия

| Элемент | Код | Статус |
|---|---|---|
| Гравитация («самогравитация») | `ext_force.cpp`: `dvxdt[dim-1] = -9.8` при `self_gravity` | `self_gravity` — `constexpr bool = false` ([param.hpp](../../../include/param.hpp)); `ext_force` не вызывается в продовом пути, гравитации в расчёте нет |
| Отталкивание от виртуальной границы (степени 12 и 4) | `ext_force.cpp`: пары «реальная ($\mathrm{itype}>0$) — виртуальная ($<0$)» | $\mathbf f = d_d\,\mathbf x_{ij}\,\big[(r_0/r)^{12}-(r_0/r)^{4}\big]/r^2$ при $r<r_0$, $d_d=10^{-2}$, $r_0=1.25\cdot10^{-5}$ м (константы зашиты в код); в продовом пути не вызывается |
| `virtual_symm_2d` (зеркальный образ на оси при `axis_sim`) | `ext_force.cpp` | определена, нигде не вызывается |
| Стенки, твёрдые поверхности | — | в коде нет |
| Начальные скорости | генераторы/начальные условия | отдельной модели нет; для закреплённых компонент обязаны быть нулевыми |
| Заданное движение границы | — | в коде нет; единственное кинематическое ГУ — Clamp (нулевая скорость) |

Осесимметричное замыкание ядра (ghost-частицы) относится к построению соседей
Solid-SPH, а не к внешним силам и здесь не описывается.

---

## 4. Реализация в коде

| Модель | Файл | Функция / сущность |
|---|---|---|
| Силы типов 1–3, диспетчер, ограничитель | [contact_force.hpp](../../../include/contact_force.hpp) | `pinball_force`, `swegle_force`, `impulse_force`, `penalty_force`, `contact::Result` |
| Радиус и $\Delta t_\text{ref}$ пары (продовый путь) | [solid_contact.hpp](../../../include/interactions/solid_contact.hpp) | `solid_contact_radius`, `solid_contact_pair_dt_ref` |
| Вклад контакта в RHS | [solid_contact.hpp](../../../include/interactions/solid_contact.hpp) | `compute_solid_contact_contributions` |
| Диагностика и $\Delta t_\text{contact}$ | [solid_contact.hpp](../../../include/interactions/solid_contact.hpp) | `compute_solid_contact_diagnostics` |
| Список контактных пар | [solid_sph.hpp](../../../include/interactions/solid_sph.hpp) | `build_solid_sph_search`, `SolidContactPair` |
| Вызов из конвейера RHS | [single_step.cpp](../../../src/single_step.cpp) | `single_step` |
| Применение $\Delta t_\text{contact}$ | [time_integration.cpp](../../../src/time_integration.cpp) | `dt = min(CFL, min_dt_contact, …)` |
| Разбор ключей | [input.cpp](../../../src/input.cpp) | `contact_force_type`, `contact_penalty_scale`, `contact_impulse_gamma`, `sph_same_body_only` |
| Старый оркестратор (тесты), внешние силы | [ext_force.cpp](../../../src/ext_force.cpp), [ext_force.hpp](../../../include/ext_force.hpp) | `ext_force`, `particle_radius`, `pair_dt_ref`, `is_contact_material_pair`, `virtual_symm_2d` |
| Clamp: проекция, реакция, маски | [clamp.hpp](../../../include/interactions/clamp.hpp) | `clamp_project_item`, `apply_clamp_projection`, `apply_clamp_friction_counts`, `mark_clamp_box`, `validate_clamps`, `prepare_clamp_workspace` |
| Источники масок | [k_geometry.cpp](../../../src/k_geometry.cpp), [fiber_input.cpp](../../../src/fiber_input.cpp) | `apply_spc_clamps`; директива `CLAMP_BOX` |

Контакты нитей — [05-fiber-contact.md](05-fiber-contact.md). Архитектурное
описание пути — [contact-particle-particle.md](../architecture/contact-particle-particle.md).
