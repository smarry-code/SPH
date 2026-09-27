# План: нити (Fiber) и вычисление RHS по схеме SPH.md

Документ описывает, как постановка из [SPH.md](SPH.md) реализуется в текущем коде проекта: что меняется в хранении частиц, как устроено вычисление правой части (RHS), в каком порядке идут этапы и чем проверяется каждый из них.

Уравнения всех взаимодействий и общие уравнения динамики — в SPH.md. Здесь только то, как их встроить в код.

---

## 1. Принятые решения

| № | Решение | Следствие |
|---|---|---|
| Р1 | Параллельность строго по SPH.md: `CpuExecutor`, поиск → primitives → interaction kernel → contributions → reduce по схеме owner computes | Solid-SPH тоже переносится на эту схему; `PairGraph`, раскраска и `AccAtomic`/`AccPlain` удаляются |
| Р2 | ODESolver — отдельный, последний этап | До него работает текущий цикл `time_integration`, в который встраиваются нити |
| Р3 | Реализуются два решателя: **явный Эйлер** и **PredictorCorrector в варианте проекта** | Эталоны BC2/BC6 сохраняются на PredictorCorrector |
| Р4 | Эйлер строго явный: $y^{n+1}=y^n+\Delta t\,f(y^n)$, координаты — по старой скорости | На незатухающих связях нити неустойчив; используется как эталонный решатель для тестов |
| Р5 | `include/particle_system_new.hpp` удаляется | Черновик противоречит SPH.md (полиморфный базовый класс, `float`, одна точка вместо массивов) и нигде не подключается |
| Р6 | Нить работает только на растяжение, изгиба нет | Сила связи равна нулю при $l\le l^0$ |
| Р7 | Разрыв связи — по расстоянию $l>l^{break}=(1+\varepsilon_b)l^0$ | Не по углу изгиба, как у Meyer |
| Р8 | `Fiber-Fiber-Contact` — буквально по SPH.md: primitive из пары $(i,j)$ и их топологических соседей, веса $w_{ij}$ Meyer | Одно пересечение учитывается несколькими парами — свойство модели Meyer, не ошибка. **Заменено (2026-09-27):** пара сегментов, один контакт на пересечение — `fiber/fiber_contact_segment_pairs.md`, SPH.md «Постановка проекта: пара сегментов» |
| Р9 | Нити только при `SPH_DIM=3` | Нужна новая 3D-цель сборки; это первый 3D-кейс Solid в проекте |
| Р10 | $E$ и $\nu$ Solid выводятся из $K$ и $G$ | Формулы — в SPH.md, раздел «Fiber-Solid-Contact» |

---

## 2. Исходное состояние кода

| Требование SPH.md | Что есть сейчас | Что делать |
|---|---|---|
| Compute-структура `Solid` | Весь текущий SPH-расчёт: `ParticleSystem` ([particle_system.hpp](../include/particle_system.hpp)), `int_force`, `con_density`, `art_visc`, `av_vel`, EOS, прочность | Перенести на новую схему RHS (этап Ф4) |
| `Solid-Solid-Contact` | Штрафной контакт в [ext_force.cpp](../src/ext_force.cpp), формулы в [contact_force.hpp](../include/contact_force.hpp), ADR-007 | Обернуть в primitive; формулы не меняются |
| `object_id` | Нет. Тело определяется по материалу `itype`; PID из k-файла затирается | Добавить `body[]` в Solid |
| Состояние только читается во время RHS | Нарушено: `virt_part` дописывает виртуальные частицы в `sys`, `int_force` пишет `p`, `c`, `exx…ryz`, `h_upgrade` пишет `hsml` | Вынести производные величины в `RhsWorkspace` (этап Ф4) |
| Параллельность | OpenMP + раскраска графа пар [pair_graph.hpp](../include/pair_graph.hpp); контакт последовательный | `CpuExecutor` поверх OpenMP (этап Ф1) |
| `SpatialIndex` | `FlatGrid::build` принимает только `ParticleSystem`, размер ячейки — из `hsml` | Обобщить: массивы координат + явный размер ячейки |
| ODESolver | Схема с полушагом зашита в [time_integration.cpp](../src/time_integration.cpp) | Выделить в решатель (этап Ф8) |
| Postprocessing | Прочность `prochnost_*` уже вызывается после шага в `time_integration` | Выделить в общий каркас (этап Ф7) |

---

## 3. Архитектура

### 3.1. Данные

Все поля — `double`. Хранение — Structure of Arrays, частицы всех объектов одного типа лежат в общих массивах.

**Solid** — текущий `ParticleSystem` с изменениями:

* добавляется `body[]` — номер тела (`object_id`). Для k-файла это PID; для встроенных генераторов (`shock_tube`) `body = itype`, чтобы критерий `Solid-Solid-Contact` на BC2/BC6 совпал с нынешним «разные материалы»;
* производные величины `p`, `c`, `exx…eyz`, `rxy…ryz` переезжают из состояния в `RhsWorkspace` (этап Ф4). Пока Ф4 не завершён, они остаются на месте.

**Fiber:**

```text
FiberState:
    x[dim], vx[dim]          // состояние
    mass                     // постоянна, из ненагруженной конфигурации
    object_id                // номер нити; меняется при разрыве
    material_id              // индекс в таблице FiberMaterial; при разрыве не меняется

FiberTopology:
    prev[], next[]           // соседи по нити, -1 на конце
    bonds: {i, j, l0, l_break}   // каждая связь один раз
    version                  // topology_version

FiberMaterial:
    E, nu, d, rho, eps_break
```

Частица — конец нити, если `prev == -1` или `next == -1`.

**Model** объединяет `SolidState`, `FiberState`, `FiberTopology` и таблицу материалов.

### 3.2. Типы взаимодействий и primitives

| Тип | Источник пар | Primitive | Вклады |
|---|---|---|---|
| `Solid-SPH` | `SpatialIndex<Solid>` + виртуальные частицы, `i<j` | `(i, j, W, ∇W)` | Solid `i`, `j` (только реальным) |
| `Solid-Solid-Contact` | `SpatialIndex<Solid>`, `i<j`, `body[i] != body[j]`, только реальные | `(i, j)` | Solid `i`, `j` |
| `Fiber-Bond` | `FiberTopology.bonds` | `(i, j, l0, l_break)` | Fiber `i`, `j` |
| `Fiber-Fiber-Contact` | `SpatialIndex<Fiber>`, `i<j`, `object_id[i] != object_id[j]` | `(i, i_prev, i_next, j, j_prev, j_next, case)` | до 4 частиц Fiber |
| `Fiber-Solid-Contact` | `Fiber → Solid`, только реальные Solid | `(a, b, s)` | Fiber `a`, `b`; Solid `s` |

Разные primitives — разные структуры и разные массивы; общего базового класса нет.

**Радиусы поиска:**

* `Fiber-Fiber-Contact`: $d+l_{max}$ — пара бусин должна находиться, даже если ближайшие точки осей лежат в середине сегментов. **Исправлено (переход на пары сегментов, `fiber/fiber_contact_segment_pairs.md`):** радиус теперь $d_\max+2L_\max$, где $L_\max$ — **текущая** (не номинальная $l_0$) максимальная длина связи `Fiber-Bond` — см. `fiber/SPH.md`, раздел «Fiber-Fiber-Contact», подраздел «Постановка проекта: пара сегментов»;
* `Fiber-Solid-Contact`: $r_f+r^{contact}_{s,max}+l_{max}$.

> **Заменено (2026-09-27).** Описанный ниже kernel (четыре случая, $\mathbf p_i$, веса $w_{ij}$) — историческая реализация Ф5. Действующая схема — пара сегментов с правилом владения: `fiber/fiber_contact_segment_pairs.md`, SPH.md «Постановка проекта: пара сегментов».

**`Fiber-Fiber-Contact` (Р8).** Вид контакта (`SurfaceSurface`, `SurfaceEnd`, `EndSurface`, `EndEnd`) определяется при построении primitive по `prev`/`next`, чтобы kernel не ветвился по топологии. Kernel:

1. строит локальные направления $\mathbf p_i$, $\mathbf p_j$ по соседям. Ориентация фиксируется так, чтобы $P_{ij}>0$ соответствовало стороне $i-1$ — этого требует формула весов с $d_p=\|\mathbf x_i-\mathbf x_{i-1}\|$;
2. вычисляет $D$, $\mathbf n$, $P_{ij}$, $P_{ji}$ по формулам своего вида контакта. Если $\|\mathbf p_i\times\mathbf p_j\|<\varepsilon$ (нити почти параллельны), система Surface-Surface вырождается, и используется формула Surface-End;
3. при $D<d$ вычисляет силу Герца $F$;
4. распределяет силу: частице $i$ — $w_{ij}F\mathbf n$, остаток $(1-w_{ij})F\mathbf n$ — соседу ($i-1$ при $P_{ij}>0$, $i+1$ при $P_{ij}\le0$); на нити $j$ — то же с $-\mathbf n$.

Каждый primitive сохраняет импульс точно: нить $i$ получает $F\mathbf n$, нить $j$ — $-F\mathbf n$.

**`Fiber-Solid-Contact`.** Пара `(f, s)` порождает только сегмент `(f, next(f))` — сегмент «принадлежит» первой бусине, поэтому дубликатов нет и сортировка не нужна. Одиночная бусина без соседей даёт точечный контакт. Чтобы узел нити не получал двойную силу, $\lambda$ принимается в полуинтервале $[0,1)$; у последнего сегмента нити — в $[0,1]$.

### 3.3. Delta, reduce, finalize

**SolidDelta** — то, что суммируется по primitives:

```text
dv[dim]      ускорение
drho         скорость изменения плотности
du           часть du/dt от давления и искусственной вязкости
eps[6]       накопители скорости деформаций
rot[3]       накопители скорости поворота
xsph[dim]    поправка XSPH
dh           только при sle = 2
```

**FiberDelta:** `dv[dim]`.

После reduce выполняется **finalize** — `parallel_for` по частицам своего типа:

* Solid: осевые члены (при $dim=2$), $\dot s$ по Яуманну, девиаторная работа $s:\dot\varepsilon/\rho$ в $\dot u$, $\dot{\mathbf x}=\mathbf v$;
* Fiber: $\dot{\mathbf x}=\mathbf v$.

Сложение `Delta` определено в своей compute-структуре; interaction type о нём не знает.

### 3.4. Граф одного вычисления RHS

```text
Solid derived (EOS: p, c; виртуальные частицы) ──┐
SpatialIndex<Solid> ─┬─ Solid-SPH ─────────────┐ │
                     ├─ Solid-Solid-Contact ───┼─ Solid reduce ─ Solid finalize ─┐
SpatialIndex<Fiber> ─┼─ Fiber-Solid ───────────┤                                 ├─ dt_limit
                     └─ Fiber-Fiber ───────────┼─ Fiber reduce ─ Fiber finalize ─┘
FiberTopology ────────── Fiber-Bond ───────────┘
```

Ветви разных interaction types независимы и запускаются как задачи `CpuExecutor`. `Fiber reduce` ждёт все ветви, дающие вклады в Fiber; `Solid reduce` — все ветви, дающие вклады в Solid. Сами reduce независимы.

### 3.5. CPU backend

**`CpuExecutor`** — единый объект на весь расчёт, реализация на OpenMP (он уже подключён; второй пул потоков, например TBB, конкурировал бы с OpenMP за ядра).

```text
parallel_for(n, grain, f)
parallel_reduce(n, init, map, combine)
parallel_scan(in, out, n)
parallel_sort(first, last, less)
task_group: run(f), wait()      // одновременное выполнение ветвей графа
num_workers(), worker_id()
```

Interaction kernels не содержат `#pragma omp`, mutex и atomic.

**Поиск.** Диапазон частиц или ячеек режется на небольшие блоки с динамическим распределением. Каждый поток пишет пары в свой буфер, затем «подсчёт → `parallel_scan` → копирование» в один непрерывный массив. Так же строятся массивы primitives.

**Вклады.** Kernel пишет через `ContributionWriter` своего потока: `emit<Fiber>(id, delta)`, `emit<Solid>(id, delta)`. Буферы — в `RhsWorkspace`, по одному на поток и compute-структуру.

**Reduce.** Диапазон ID каждой compute-структуры делится на шарды (число шардов заметно больше числа потоков). Вклады раскладываются по шардам: подсчёт по (поток, шард) → префиксная сумма → копирование. Затем каждый шард суммирует свой диапазон один, обычным `+=`.

**Детерминированный режим** (флаг конфигурации): внутри шарда вклады сортируются по `(interaction_type, primitive_id, slot)` перед суммированием. В обычном режиме порядок суммирования может меняться между запусками.

**`RhsWorkspace`** владеет пространственными индексами, буферами пар, primitives, contributions, шардов, производными величинами Solid и виртуальными частицами. Память переиспользуется между вызовами RHS и только растёт.

**Риск памяти.** SolidDelta — около 17 `double` в 3D (около 11 в 2D) на каждую сторону пары. Объём нужно замерить на BC6 до начала переноса Solid. Если он неприемлем, запасной вариант — суммирование на стороне владельца по списку соседей: договор SPH.md сохраняется, но каждая пара считается дважды.

### 3.6. Solid: что перестаёт писаться в состояние

| Сейчас | После Ф4 |
|---|---|
| `urs()` пишет `p`, `c` в `sys` внутри `int_force` | стадия Solid derived, результат в `RhsWorkspace` |
| `int_force` пишет `exx…eyz`, `rxy…ryz` в `sys` | накопители SolidDelta, затем finalize |
| `virt_part` дописывает виртуальные частицы в `sys` за `ntotal` | `GhostSet` в `RhsWorkspace`; вкладов к виртуальным частицам нет |
| `h_upgrade` пишет `hsml` | `sle=1` — Postprocessing; `sle=2` — переменная состояния |
| `prochnost_*` читают `sys.exx…` | Postprocessing получает последний RHS (только чтение) |

CFL-ограничение читает $c$ из производных величин последнего RHS.

### 3.7. Шаг по времени

`dt_limit` каждого вычисления RHS — минимум (через `parallel_reduce`) по:

* `Solid-Solid-Contact` — как сейчас (`contact::Result::dt_limit`);
* `Fiber-Bond` — $0.8\sqrt{m_{red}/k}$, $k=EA/l^0$; считается для всех связей независимо от знака деформации, чтобы шаг не прыгал;
* `Fiber-Fiber-Contact`, `Fiber-Solid-Contact` — $0.8/\sqrt{k_{eff}/m_{red}}$, $k_{eff}=1.5F/\delta$.

Итоговый шаг: $\Delta t=\min(\text{CFL},\ \texttt{dt\_limit})$.

### 3.8. Postprocessing

Выполняется после шага ODESolver в две фазы:

1. поиск событий — `parallel_for` только на чтение, события в буферы потоков;
2. применение — последовательно.

| Правило | Событие / действие |
|---|---|
| Мизес (itype 3–9) | возврат девиатора на поверхность текучести, по частице |
| JH-2 (itype 10–11) | прочность и повреждённость, по частице; использует $\dot\varepsilon$, $\dot\rho$ последнего RHS |
| `sle = 1` | пересчёт $h$ |
| Разрыв нити | `BreakFiberBond(k)`: связь удаляется, `next[i] = prev[j] = -1`, отделившейся части присваивается новый `object_id`, `topology_version++` |

При изменении `topology_version` списки связей, incidence-структуры `Fiber-Bond` и всё, что от них зависит, строятся заново перед следующим RHS.

### 3.9. ODESolver

Интерфейс:

```text
OdeSolver::step(Model&, RhsPipeline&, double& dt)
```

Операции над состоянием определены в каждой compute-структуре и выполняются `parallel_for`:

```text
State::save(y_n)                      // сохранить y^n
State::advance(y, y_base, a, f)       // y = y_base + a * f  (+ правила структуры)
```

Правила, специфичные для Solid, живут в его `advance`, а не в решателе: ограничение $u\ge0$, прибавление поправки XSPH к скорости, обновление координат по скорости.

**PredictorCorrector в варианте проекта** — перенос схемы из [time_integration.cpp](../src/time_integration.cpp) без изменений. Решатель хранит производные прошлого шага $f^{n-1}$.

Первый шаг ($n=1$):

1. $f^0=\mathrm{RHS}(y^0)$, $\Delta t$ из CFL и `dt_limit`;
2. $u,\rho,s\mathrel{+}=\tfrac{\Delta t}{2}f^0$; $\mathbf v\mathrel{+}=\tfrac{\Delta t}{2}\dot{\mathbf v}^0+\tilde{\mathbf v}$; $\mathbf x\mathrel{+}=\Delta t\,\mathbf v$.

Шаги $n\ge2$:

1. сохранить $y^n$ ($u,\rho,\mathbf v,s$);
2. предиктор: $u,\rho,\mathbf v,s\mathrel{+}=\tfrac{\Delta t_{n-1}}{2}f^{n-1}$;
3. $f^n=\mathrm{RHS}(y^{*})$, новый $\Delta t$;
4. корректор: $u,\rho,s=y^n+\Delta t\,f^n$; $\mathbf v=\mathbf v^n+\Delta t\,\dot{\mathbf v}^n+\tilde{\mathbf v}$; $\mathbf x\mathrel{+}=\Delta t\,\mathbf v$.

Нити интегрируются той же схемой (только $\mathbf x$, $\mathbf v$). Одно вычисление RHS на шаг.

**Явный Эйлер:**

1. $f=\mathrm{RHS}(y^n)$, $\Delta t$;
2. $y^{n+1}=y^n+\Delta t\,f$ для всех переменных, включая $\mathbf x^{n+1}=\mathbf x^n+\Delta t\,\mathbf v^n$.

На незатухающих упругих связях энергия растёт в $(1+\omega^2\Delta t^2)$ раз за шаг при любом $\Delta t$. Эйлер — эталон для тестов порядка сходимости, рабочий решатель — PredictorCorrector.

Порядок одного шага: ODESolver (RHS внутри) → Postprocessing → вывод кадров.

---

## 4. Этапы

У каждого этапа — содержание, тесты и критерий готовности. Тесты оформляются как остальные в проекте: отдельные исполняемые файлы в `tests/`, счётчик `failures`, регистрация в CTest.

### Ф0. Документы

* удалить `include/particle_system_new.hpp`;
* SPH.md: разделы «Solid -- Solid», «Общие уравнения динамики», постановка Fiber-Bond (только растяжение, разрыв по расстоянию), коэффициент $K_{fs}$ — **выполнено**;
* ADR-008 «Архитектура RHS и нити»: решения Р1–Р10, пометка, что ADR-006 (раскраска графа пар) заменяется после Ф4.

**Готово, когда** ADR-008 принят, черновик удалён, сборка и все тесты проходят как раньше.

### Ф1. Инфраструктура без физики

* `CpuExecutor`;
* `ContributionWriter`, буферы вкладов, reduce по шардам, детерминированный режим;
* `SpatialIndex` — обобщённая `FlatGrid` (массивы координат + размер ячейки), запросы «свой тип» и «чужой тип»;
* `RhsWorkspace`.

Тесты:

* `test_executor`: `parallel_for`/`reduce`/`scan`/`sort` против последовательных эталонов;
* `test_spatial_index`: пары совпадают с полным перебором, в том числе для двух разных типов;
* `test_contributions`: сумма вкладов совпадает с последовательной; в детерминированном режиме результат побитно одинаков при 1, 2 и 7 потоках.

**Готово, когда** тесты проходят; существующий код не изменён.

### Ф2. Модель данных

* `body[]` в Solid; заполнение из PID k-файла и в генераторах;
* `FiberState`, `FiberTopology`, `FiberMaterial`, `Model`;
* ввод нитей: отдельный файл ломаных (материал, точки каждой нити), ссылка из `Initial_Conditions.txt`; вычисление масс и $l^0$, $l^{break}$;
* вывод `fiber_output_*.txt` (`object_id`, координаты, скорости, удлинение связей);
* проверка конфигурации: нити заданы, а `dim != 3` — ошибка;
* CMake-цель со `SPH_DIM=3`;
* перегрузка `time_integration`: прежняя сигнатура продолжает работать для `test_e2e_bc2` и `test_e2e_bc6`.

Тесты: разбор файла нитей, массы, топология `prev`/`next`, отказ при `dim != 3`.

**Готово, когда** без нитей BC2/BC6 проходят без изменений: `ctest -L regression` проходит без перегенерации эталона.

### Ф3. Fiber-Bond

* primitive и kernel: сила только на растяжение;
* `dt_limit` связей;
* нити интегрируются в текущем цикле `time_integration` той же схемой, что Solid; шаг — общий минимум.

Тесты (`test_fiber_bond`):

* сила двух бусин против формулы; при сжатии сила равна нулю;
* частота колебаний в растянутой фазе;
* волна растяжения в прямой нити идёт со скоростью $\sqrt{E/\rho_f}$;
* сумма импульсов сохраняется.

**Готово, когда** тесты проходят; BC2/BC6 не изменились (`ctest -L regression` проходит без перегенерации эталона).

### Ф4. Перенос Solid на новую схему RHS

Самый крупный и рискованный этап.

* стадия Solid derived: EOS в `RhsWorkspace`, `GhostSet` вместо дописывания в `sys`;
* `Solid-SPH`: поиск, primitive `(i, j, W, ∇W)`, kernel со всеми парными членами SPH.md, reduce, finalize;
* `Solid-Solid-Contact` из `ext_force` в primitive; диагностика `contact_pair_count`, `contact_max_delta` — через `parallel_reduce`;
* до начала — замер памяти вкладов на BC6 (см. 3.5).

Проверка в два шага:

1. стенд `test_rhs_equivalence` (`SPH_DIM=2`): старый и новый RHS на одном состоянии, все производные совпадают с относительной точностью ~1e-12 (побитно — нет, порядок суммирования другой);
2. регрессия `ctest -L regression` проходит без перегенерации эталона (`tests/test_regression.cpp`: все поля всех частиц BC2/BC6 после N = 1 и N = K шагов, типы контакта 0–3; то же в 3D — `test_regression_3d`, уменьшенный BC2);
3. e2e BC2 и BC6 проходят по своим текущим критериям;
4. ветка `dim = 3`, в которой будут жить нити, не сломана: `ctest -L 3d` (`test_kinematics_3d` — алгебра девиатора и поворота, `test_e2e_3d` — сохранение импульса, симметрия решётки, баланс энергии без XSPH).

После этого удалить `direct_find`, `link_list`, `PairGraph`, `AccAtomic`/`AccPlain`, старый путь `single_step`; обновить документацию по параллельности и поиску соседей.

**Готово, когда** все четыре шага проверки пройдены и старый путь удалён.

### Ф5. Fiber-Fiber-Contact

* поиск, построение primitive с видом контакта, kernel по 3.2;
* `dt_limit` контакта.

Тесты (`test_fiber_contact`):

* перпендикулярные нити: $D$, $\mathbf n$, $P_{ij}$, $w_{ij}$ против ручного расчёта;
* почти параллельные нити: нет NaN, срабатывает переход на Surface-End;
* Surface-End и End-End;
* сумма импульсов по всем вкладам равна нулю;
* диагностика: число пар на одно пересечение.

**Готово, когда** тесты проходят; e2e: две нити, летящие навстречу крест-накрест, отскакивают без взаимного проникания.

### Ф6. Fiber-Solid-Contact

* `elastic_by_itype` в [materials.hpp](../include/materials.hpp): $E$, $\nu$ из $K$ и $G$; для itype 1–2 — «не определено», контакт не формируется, предупреждение при старте;
* поиск Fiber → Solid, primitive `(a, b, s)`, kernel, `dt_limit`;
* вклады в Solid идут в Solid reduce. Если Ф4 ещё не завершён — временно прибавляются к ускорению Solid в цикле сложения сил `single_step`.

Тесты (`test_fiber_solid`):

* $\lambda$ и clamp, контакт с концом нити;
* узел нити не получает двойную силу;
* $\mathbf F_a+\mathbf F_b+\mathbf F_s=0$;
* $E$, $\nu$ для материалов проекта против справочных значений.

e2e (`SPH_DIM=3`): сначала smoke-тест одного Solid в 3D, затем нить, падающая на пластину.

**Готово, когда** тесты проходят, e2e без проникания нити в пластину.

### Ф7. Postprocessing

* каркас «поиск событий → применение»;
* перенос прочности (Мизес, JH-2) и `sle = 1`;
* разрыв связей по $l>l^{break}$, новый `object_id`, `topology_version`, сброс устаревших структур.

Тесты (`test_fiber_break`):

* связь рвётся ровно при пороге;
* после разрыва две половины — разные `object_id` и взаимодействуют через `Fiber-Fiber-Contact`;
* после разрыва primitives `Fiber-Bond` построены заново.

**Готово, когда** тесты проходят; BC2/BC6 не изменились (`ctest -L regression` проходит без перегенерации эталона).

### Ф8. ODESolver

* интерфейс `OdeSolver`, операции `save`/`advance` для Solid и Fiber;
* `PredictorCorrector` (вариант проекта) и явный `Euler`;
* выбор решателя в `Initial_Conditions.txt`; по умолчанию — `PredictorCorrector`;
* `time_integration` сводится к циклу: `solver.step` → Postprocessing → вывод.

Тесты (`test_ode_solver`):

* Эйлер: первый порядок сходимости на задаче с аналитическим решением;
* PredictorCorrector: на BC2/BC6 результат совпадает с прежним циклом;
* нить без контактов: энергия ограничена на PredictorCorrector.

**Готово, когда** старый цикл удалён, e2e BC2/BC6 проходят, `ctest -L regression` проходит без перегенерации эталона.

### Зависимости этапов

```text
Ф0 → Ф1 → Ф2 → Ф3 ──────────┐
                └→ Ф4 ──────┼→ Ф5 → Ф6 → Ф7 → Ф8
```

Ф3 и Ф4 можно вести одновременно: нитям не нужен перенесённый Solid.

---

## 5. Риски

| Риск | Где | Что делаем |
|---|---|---|
| Память на вклады Solid-SPH | Ф4 | Замер на BC6 до переноса; запасной вариант — суммирование на стороне владельца |
| Расхождение BC2/BC6 после переноса Solid | Ф4 | Стенд эквивалентности RHS на одном состоянии до запуска e2e; регрессия `test_regression`: контрольная точка N = 1 — тот же RHS плюс один шаг, N = K — накопление |
| Трёхмерный Solid не проверялся | Ф6 | Smoke-тест одного Solid в 3D до контакта с нитью |
| Плотная сетка `FlatGrid` в 3D при длинных нитях | Ф1, Ф5 | Ограничение числа ячеек; при необходимости — хешированная сетка |
| ~~Кратный учёт пересечения в Fiber-Fiber (Meyer)~~ — **закрыто** | Ф5 | Исправлено переходом на схему «пара сегментов — одна сила», см. `fiber/fiber_contact_segment_pairs.md` и `fiber/SPH.md`, раздел Fiber-Fiber-Contact → «Постановка проекта: пара сегментов»: правило владения гарантирует ровно один контакт на пересечение |
| Неустойчивость явного Эйлера на нитях | Ф8 | Эйлер — только для тестов; рабочий решатель — PredictorCorrector |
| Шаг по времени от жёстких нитей | Ф3 | `dt_limit` связей, вывод причины ограничения шага, как сейчас для контакта |

**Заметка (возможная будущая оптимизация, не сделана).** Поиск `Fiber-Fiber-Contact` по-прежнему идёт по AABB **бусин** с радиусом $d_\max+2L_\max$ (см. таблицу радиусов поиска, п. 3.2) — `SpatialIndex` при переходе на пары сегментов не менялся. Поиск по AABB **сегментов** вместо бусин мог бы сузить радиус обратно к масштабу $d+l_\max$ и уменьшить число кандидатов широкой фазы; отдельная задача на будущее.

---

## 6. Размещение файлов

`CMakeLists.txt` собирает только `src/*.cpp` без подкаталогов, поэтому исходники лежат в `src/` с префиксами, а заголовки — в подкаталогах `include/`.

```text
include/rhs/cpu_executor.hpp
include/rhs/contributions.hpp
include/rhs/spatial_index.hpp
include/rhs/rhs_workspace.hpp
include/rhs/rhs_pipeline.hpp
include/fiber/fiber_system.hpp          FiberState, FiberTopology, FiberMaterial
include/interactions/solid_sph.hpp
include/interactions/solid_contact.hpp
include/interactions/fiber_bond.hpp
include/interactions/fiber_fiber.hpp
include/interactions/fiber_solid.hpp
include/interactions/hertz.hpp          общий hertz_force(E*, R*, delta)
include/ode/ode_solver.hpp              OdeSolver, Euler, PredictorCorrector
include/postprocessing.hpp

src/rhs_pipeline.cpp
src/fiber_input.cpp
src/fiber_output.cpp
src/postprocessing.cpp
src/ode_solver.cpp

tests/test_executor.cpp
tests/test_spatial_index.cpp
tests/test_contributions.cpp
tests/test_fiber_bond.cpp
tests/test_rhs_equivalence.cpp
tests/test_fiber_contact.cpp
tests/test_fiber_solid.cpp
tests/test_fiber_break.cpp
tests/test_ode_solver.cpp
```
