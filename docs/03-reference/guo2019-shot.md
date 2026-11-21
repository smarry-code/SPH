# Утилита `guo2019_shot`

Справочник ключей и выходных файлов утилиты
[`guo2019_shot`](../../cases/guo2019/experiment/guo2019_shot.cpp): серийные выстрелы стального
шара по мишени Guo 2019 (алюминиевая пластина и/или пакет слоёв ткани). Сама утилита
не входит в `ctest` как тест, но её CLI проверяют `checkpoint_cli_test`, `stitch_kt_cli_test`, `mosaic_kt_cli_test` (`ctest`). Сцену и монитор задают заголовки
[`guo2019_case.hpp`](../../cases/guo2019/experiment/guo2019_case.hpp),
[`guo2019_scene.hpp`](../../cases/guo2019/experiment/guo2019_scene.hpp),
[`guo2019_monitor.hpp`](../../cases/guo2019/experiment/guo2019_monitor.hpp).

Модель ткани — [«Ткани»](../04-explanation/models/07-fabrics.md), генератор
переплетения — [руководство](../02-how-to/02-plain-weave-generator.md), рестарт —
[«Рестарт из чекпоинта»](../02-how-to/05-restart-checkpoints.md).

## Назначение и запуск

Каждый выстрел — свежая сцена (`build_guo_scene`) в своём подкаталоге
`<out>/<имя выстрела>/`; на время выстрела утилита делает `chdir` в этот каталог, поэтому
`log.txt` и кадры пишутся туда же. Серия: основные скорости, затем (с `--bracket`) до 6
выстрелов бисекции границы «пробит / не пробит» (пока интервал больше 5 м/с), затем
повторы (`--repeat`) для скоростей в пределах 5 м/с от $V_{BL}$ модели, где $\bar e$ — средняя энергия пробивных
выстрелов не выше чем на 15 м/с над наибольшей скоростью «не пробит», $V_{BL}=\sqrt{2\bar e/m_s}$. Сводки серии — `shots.csv` и `vbl.txt` в `<out>`.

```bash
guo2019_shot --config=<имя> --level=base|fine|smoke|scaled \
    (--v0=<м/с> | --paper-shots | --v0-list=a,b,c) --out=<каталог> [ключи]
```

Обязательны `--config`, `--level`, `--out` и **ровно одно** из `--v0`, `--paper-shots`,
`--v0-list`. Пример (3 слоя ткани, стёжка сеткой патчей по 4 нити, чекпоинты):

```bash
guo2019_shot --config=3T --level=smoke --v0=147 --out=run1 \
    --stitch=grid --stitch-tile=4 --checkpoint-every=5000 --threads=4
```

`-h` / `--help` печатает справку и возвращает 0.

## Ключи

### Сцена и серия

| Ключ | Значение / диапазон | По умолчанию | Описание |
|---|---|---|---|
| `--config=` | `A`, `1T`, `1T_A` (`1T/A`), `A_1T` (`A/1T`), `<k>T`, k = 1..64 | обязателен | Мишень: `A` — пластина; `1T_A` — ткань спереди, Al сзади; `A_1T` — наоборот; `<k>T` — k слоёв ткани без пластины (1T, 2T, 4T — как в статье, прочие k — оценка без данных статьи). Предел k — `kMaxPlies` = 64 |
| `--level=` | `base`, `fine`, `smoke`, `scaled` | обязателен | Уровень разрешения/габаритов |
| `--v0=` | > 0, м/с | — | Один выстрел |
| `--v0-list=` | список > 0 через запятую | — | Несколько скоростей |
| `--paper-shots` | флаг | — | Скорости опытов статьи для конфигурации (дубликаты считаются один раз); включает `--bracket` и `--repeat=2`. Для `<k>T` с $k\notin\{1,2,4\}$ опытов статьи нет: список пуст, выстрелов не будет |
| `--out=` | путь | обязателен | Каталог серии |
| `--mu=` | $\ge 0$ | 0,3 (`kMuDefault`) | Трение: задаёт **оба** коэффициента, $\mu_{ff}$ (нить–нить) и $\mu_{fs}$ (нить–Solid) |
| `--mu-ff=`, `--mu-fs=` | $\ge 0$ | из `--mu` | Переопределяют только свой коэффициент, приоритет над `--mu` независимо от порядка ключей |
| `--t-max=` | > 0, с | по конфигурации и уровню (`A` 600 мкс, `1T` 1000, `1T_A`/`A_1T` 800, `2T` 800, `4T` и `<k>T` 600; на уровне `scaled` вдвое короче) | Потолок времени расчёта |
| `--no-body-isolation` | флаг | изоляция тел включена | Выключает `sph_same_body_only` |
| `--plate-nz=` | 1..64 | по уровню | Число слоёв пластины по толщине (сравнение стоимости); в имя добавляется `_nz` |
| `--repeat=` | 1..100 | 1; 2 с `--paper-shots` | Повторы около $V_{BL}$ |
| `--bracket` / `--no-bracket` | флаг | вкл. только с `--paper-shots` | Бисекция границы исхода |

### Шар

| Ключ | Значение | По умолчанию | Описание |
|---|---|---|---|
| `--ball-d=` | > 0, м | по уровню | Диаметр шара |
| `--ball-mass=` | > 0, кг | по уровню | Масса (плотность выводится) |
| `--ball-dx=` | > 0, м | по уровню | Шаг решётки частиц шара |

Любой из трёх меняет имя каталога (`_bD…m…dx…`) и добавляет в `shot.txt` ключи
`ball_D`, `ball_m`, `ball_dx`, `ball_rho`.

### Габариты ткани (только `<k>T`, без пластины)

| Ключ | Значение | По умолчанию | Описание |
|---|---|---|---|
| `--fabric-n=` | 2..1000 | по уровню | Число нитей $n_x = n_y$ |
| `--fabric-clamp-edge=` | > 0, м | по уровню | Граница зажима $c_e$ |
| `--fabric-overhang-y=` | > 0, м | по уровню | Выступ нитей $e_y$ |

Для конфигураций с пластиной — ошибка выстрела («габариты ткани … только для мишеней без
пластины»). Имя каталога: `_fN<n>ce<мм>ey<мм>`; в `shot.txt` — `fabric_n`,
`fabric_clamp_edge`, `fabric_overhang_y`.

### Материал нити

| Ключ | Значение | По умолчанию | Описание |
|---|---|---|---|
| `--eps-break=` | > 0 | из сцены | $\varepsilon_{break}$ материала 0; $l_{break} = (1+\varepsilon_{break})\,l_0$ пересчитывается для всех связей цепочек. Стежок (материал 1) не меняется |
| `--fiber-E=` | > 0, Па | из сцены | Модуль Юнга материала 0 |
| `--fiber-inelastic=` | файл | — | Файл директив `PLASTIC`/`PLASTIC_TABLE`/`RATE`/`VISCO`/`FAILURE`/`DAMAGE` (без `MATERIAL`/`FIBER`) для материала 0; разбирается парсером `.fib`. Файл должен существовать; директивы обязаны включить неупругость или критерий разрыва, иначе ошибка выстрела |

Ключи действуют только на материал 0; стежки остаются упругими.

### Контакт нитей

| Ключ | Значение | По умолчанию | Описание |
|---|---|---|---|
| `--fs-contact=`, `--ff-contact=` | 0..4 | как в `cases/guo2019/experiment/initial_conditions.txt` (ключи там закомментированы, то есть 0, Герц) | Закон нормальной силы нить–Solid и нить–нить (0 — Герц, 4 — Герц + демпфер) |
| `--ff-contact-predictor=`, `--fs-contact-predictor=` | 0, 1 | 0 | Законы 1–3 как проекция скорости нормального контакта |
| `--contact-proj-iters=` | 1..100 | 8 | Число итераций Якоби проекции |
| `--contact-proj-warm=` | 0, 1 | 1 | Тёплый старт проекции; 0 даёт предупреждение |
| `--impulse-gamma=` | $(0;1{,}6]$ | 1 | $\gamma$ законов 1–3 (в режиме проекции игнорируется с предупреждением) |
| `--impulse-push-tau=` | $\ge 0$, с | 0 | Выталкивание закона 3 (источник энергии, печатается предупреждение; в режиме проекции — ошибка) |
| `--impulse-push-depth=` | $[0;1)$ | 0,1 | Глубина выталкивания, доля радиуса |
| `--impulse-push-vmax=` | $\ge 0$, м/с | 0 | Ограничение скорости выталкивания |
| `--dt-scale=` | $(0;1]$ | 1 | Множитель шага |
| `--contact-stats=N`, `--hold-n` | $N \ge 1$ | — | Диагностика контактов; только сборка с `-DSPH_CONTACT_STATS=ON` (иначе код 2); в режиме проекции не поддерживается |

Правила режима проекции проверяет `validate_fiber_contact_projection`: предиктор 1
требует закона контакта 1–3; `*_impulse_push_tau > 0` — ошибка (исключение выстрела);
$\gamma \ne 1$ и `warm = 0` — предупреждения. В режиме проекции в `monitor.csv` и
`shot.txt` появляются ключи проекции (см. ниже), в `shot.txt` — `ff_contact_predictor`,
`fs_contact_predictor`, `contact_proj_iters`, `contact_proj_warm`.

### Стёжка слоёв

| Ключ | Значение | По умолчанию | Описание |
|---|---|---|---|
| `--stitch=` | `none`, `grid`, `seam` | `none` | `grid` — вертикальные стежки между соседними слоями на линиях шва (нужно $\ge 2$ слоёв); `seam` — шов-строчка в плоскости каждого слоя (любая конфигурация с тканью, в т. ч. 1T) |
| `--stitch-tile=` | 1..1000 | 3 | Сетка патчей $N \times N$; $n_x$, $n_y$ должны делиться на $N$ |
| `--stitch-patch=` | 1..1000 | — | Кусочки по $P$ нитей симметрично от центра, у краёв остаток; $P < n$; при нечётном $n$ — $P$ нечётно; несовместим с `--stitch-tile` |
| `--stitch-lines=` | `x`, `y`, `xy` | `xy` | Какие линии шва: `x` — линии $x=\text{const}$ (идут вдоль $y$, прокалывают нити основы), `y` — линии $y=\text{const}$ (вдоль $x$, прокалывают нити утка) |
| `--stitch-every=` | 1..100000 | 1 | Прокол на каждой $N$-й нити (индекс от края) |
| `--stitch-E=` | > 0, Па | $8{,}5\cdot10^{10}$ | Модуль стежка |
| `--stitch-d=` | > 0, м | $0{,}2\cdot10^{-3}$ | Диаметр стежка |
| `--stitch-eps-break=` | > 0 | 0,03 | Предельная деформация стежка |

Ключи `--stitch-*` без `--stitch=grid|seam` (и без `--random-cut`) — ошибка. Нити
проходят сквозь стыки патчей непрерывно; патчи и слои соединяют только стежки.
Подробности геометрии — в [главе о тканях](../04-explanation/models/07-fabrics.md).

### Мозаика (случайная раскройка)

| Ключ | Значение | По умолчанию | Описание |
|---|---|---|---|
| `--random-cut=` | 1..10000 | — | Число фрагментов раскройки Вороного; включает мозаичную ткань вместо сплошного полотна |
| `--cut-seed=` | $\ge 0$ | 1 | Зерно генератора |
| `--layers-per-cut=` | 1..64 | 0 (весь пакет — одна раскройка) | Слоёв на одну раскройку |
| `--cut-gap=` | > 0, м | $d$ нити | Зазор между фрагментами |
| `--cut-rot-max=` | 0..90, град | 25 | Максимальный поворот фрагмента |
| `--cut-stitch-step=` | $\ge 0$, м | шаг нитей $p$ | Шаг стежков вдоль шва; 0 — без стежков |
| `--cut-size=` | `WxH`, м, оба > 0 | габарит ткани уровня | Своё поле мозаики по центру (0, 0) |
| `--cut-clamp-edge=` | > 0, м | $c_e$ уровня | Граница зажима $\lvert y\rvert \ge c_e$ |

Материал стежков мозаики — `--stitch-E/-d/-eps-break`. `--cut-*` и `--layers-per-cut`
без `--random-cut` — ошибка.

### Вывод, время, потоки

| Ключ | Значение | По умолчанию | Описание |
|---|---|---|---|
| `--frames=` | `off`, `camera` | `off` | `camera` — кадр раз в 25 мкс (`kBaseSaveT`) |
| `--save-t=` | > 0, с | — | Интервал кадров (перекрывает `--frames`) |
| `--probe-frame` / `--no-probe-frame` | флаг | включён только при `--max-steps` | Кадр конечного состояния для замера размера кадра |
| `--n-mon=` | $\ge 1$ | 1 | Строка `monitor.csv` раз в $N$ шагов |
| `--n-diag=` | $\ge 1$ | 10 | Дорогая диагностика (туннелирование, $m_{ej}$ и др.) раз в $N$ шагов |
| `--energy-log=` | $N \ge 1$ | — | `energy_log.csv` каждые $N$ шагов; только сборка с `-DSPH_ENERGY_LOG=ON` (иначе код 2) |
| `--threads=` | 1..4096 | как в конфигурации | Число потоков OpenMP |
| `--time-limit=` | $\ge 0$, с | 0 (без лимита) | Лимит по стенному времени (`RunGuard`, останов `guard_reason=WallClock`). Контроль «застоя» включён всегда: оценка остатка более $2\cdot10^6$ шагов на 200 шагах подряд даёт `Stalled` |
| `--max-steps=` | $\ge 1$ | — | Лимит числа шагов; при рестарте — общий счётчик |

Первый кадр (время 0, номер `00001`: счётчик шагов растёт в начале шага) пишется всегда; при
`--frames=off`, без `--save-t` и без пробного кадра он после замера размера удаляется. Пробный
кадр конечного состояния имеет номер `99999`.

### Чекпоинт и рестарт

| Ключ | Значение | По умолчанию | Описание |
|---|---|---|---|
| `--checkpoint` | флаг | выкл. | Чекпоинт `checkpoint.ckpt` (+ `.prev`) в каталоге выстрела, финальный — при останове |
| `--checkpoint-every=` | $N \ge 1$ шагов | 10000 | Период; сам включает чекпоинты |
| `--restart` | флаг | — | Продолжить из `checkpoint.ckpt` (или `.prev`); нет обоих — ошибка выстрела |
| `--force` | флаг | — | Разрешить менять физические настройки при рестарте; только с `--restart` |

## Правила совместимости и коды выхода

- `--random-cut` несовместим с `--stitch=grid|seam` и с `--config=A`.
- `--stitch=grid` требует $\ge 2$ слоёв; `--stitch=seam` — конфигурации с тканью.
- `--stitch-patch` и `--stitch-tile` взаимоисключающие; проверка делимости использует
  $n_x$, $n_y$ уровня или `--fabric-n`.
- `--fabric-*` — только `<k>T`: иначе сцена бросает исключение (код выхода 0, см. ниже).
- `--force` без `--restart` — ошибка.

Код **2**: неизвестный аргумент, неверное значение, не заданы обязательные ключи, нарушены
правила выше (кроме `--fabric-*` не для `<k>T`), ключ диагностики в сборке без неё; перед выходом печатается справка. Код
**0**: нормальное завершение, `--help`, а также серия, где выстрел бросил исключение (в том
числе `--restart` без файла, `--restart` с другой сценой, `--fabric-*` не для `<k>T`, ошибки сцены) — ищите в выводе
`ОШИБКА` и `ВНИМАНИЕ: выстрелов с исключением`.

## Имя каталога выстрела

```text
<cfg>_<level>[_nz<N>][_bD<мм>m<г>dx<мм>][_fN<n>ce<мм>ey<мм>][_rc<N>s<S>k<k>[_f<W>x<H>c<c_e>]]_v<V0>_r<rep>
```

`V0` — с одним знаком после запятой (`v147.0`), `rep` — номер повтора. Для `--cut-size` /
`--cut-clamp-edge` по отдельности недостающая часть пишется как `L` (значение уровня);
размеры в мм. Примеры: `2T_base_rc30s1k2_v147.0_r0`, `…_rc7s1k0_f50x70c25_v550.0_r0`.

!!! warning "Тега стёжки в имени нет"
    Выстрелы с `--stitch=…` и без неё (и с разными `--stitch-*`) получают **одно и то же**
    имя каталога. Разные варианты стёжки запускайте в разные `--out`, иначе файлы
    перезапишутся.

## Выходные файлы

В `<out>`: `shots.csv` (по строке на выстрел; заголовок — ключи первого выстрела, порядок
как в `shot.txt`) и `vbl.txt` (интервал [н; п], $\bar e$, $V_{BL}$; для 2T/4T — значения
статьи). Обновляются после каждого выстрела.

В каталоге выстрела:

| Файл | Содержимое |
|---|---|
| `shot.txt` | Одна строка `ключ=значение` (см. ниже) |
| `monitor.csv` | Построчный монитор |
| `scene_id.txt` | Идентификатор сцены (только с `--checkpoint`) |
| `checkpoint.ckpt`, `.prev` | Чекпоинт (формат 6) |
| `sph_output_*.txt`, `fiber_output_*.txt` | Кадры Solid и нитей |
| `log.txt` | Журнал конвейера |
| `guo2019_fabric.fib`, `guo2019_fabric_layout.txt` | Сгенерированная ткань и раскладка нитей |
| `guo2019_fabric.cut.json` | Описание раскройки (только мозаика) |
| `fiber_inelastic_material.fib` | Временный файл разбора `--fiber-inelastic` (только с этим ключом) |
| `energy_log.csv`, `contact_stats*.csv` | Только в сборках с диагностикой |

### `shot.txt`

Первая часть пишется монитором (`mu` — значение `--mu`; `guard_reason` — `None`, `WallClock`,
`Stalled`, `StepLimit` или `Monitor` — штатный останов по исходу либо `--t-max`): `config level V0 mu outcome V_R e t_end steps wall_s
dt_first dt_min E_tot_max/E_tot_0 E_tot_end/E_tot_0 F_clamp_s_z_max F_clamp_f_z_max
n_clamped_s n_clamped_f n_broken_v n_broken_h w_f_max t_first_break w_f_at_first_break
w_p_max m_ej n_tunnel max_delta_ff/d max_delta_fs/d contacts_per_crossing guard_reason
eps_p_max W_p n_bonds_plastic`. Для пакета из $\ge 2$ слоёв далее: `n_plies` и по слою
`n_broken_v_pK n_broken_h_pK t_first_break_pK w_f_at_first_break_pK
pullin_at_first_break_pK pullin_end_pK`; со стежками (только в этом блоке) `n_stitch_broken_end
slip_mean_end`; затем `contacts_per_crossing_inter max_delta_ff_inter/d
ply_order_violations t_ply_order_first ply_gap_min ply_gap_max`. В режиме проекции в конце
монитор пишет `proj_resid_max proj_resid_rel_max W_proj_pos_max proj_slots_max
proj_resid_load_max proj_resid_post_max proj_resid_load_post_max W_warm_pos_max
proj_warm_scale_min proj_warm_scaled_steps`.

Утилита дописывает (в этом порядке): `rep t_max ntotal n_ball n_plate plate_nz n_beads
n_bonds n_mon n_diag scene_build_s wall_per_step_s dt_mean steps_forecast wall_forecast_h
monitor_wall_s monitor_share_pct`; `monitor_state` (только при `--restart`: `restored` —
вердикт и максимумы за весь выстрел, `partial` — только с шага рестарта);
`frames_sph_n frames_fiber_n frames_bytes_total frame_sph_bytes frame_fiber_bytes`;
`[fiber_inelastic pl_sy pl_H rate_model visco_eta fail_eps_p fail_sigma_u
(damage_deps_f damage_gf)] fiber_eps_break fiber_E impulse_push_tau impulse_push_depth
impulse_push_vmax impulse_gamma dt_scale` (первая группа — только при `--fiber-inelastic`,
`damage_*` — если у материала включена повреждённость) — вся часть только если заданы
`--eps-break`, `--fiber-E` или `--fiber-inelastic`; ключи проекции (см. выше); `mu_ff
mu_fs`; со стёжкой — `stitch stitch_patch|stitch_tile stitch_lines stitch_every stitch_E
stitch_d stitch_eps_break n_stitch`; с мозаикой — `random_cut cut_seed layers_per_cut
cut_count cut_size cut_clamp_edge cut_gap cut_rot_max cut_stitch_step cut_min_cross_dist
stitch_E stitch_d stitch_eps_break n_stitch n_stitch_in n_stitch_through`; с `--ball-*` —
`ball_D ball_m ball_dx ball_rho`; с `--fabric-*` — `fabric_n fabric_clamp_edge
fabric_overhang_y`.

### `monitor.csv`

Заголовок:

```text
step,t,dt,z_c,v_z,KE_s,KE_p,KE_f,U_b,U_int,E_tot,Fs_x,Fs_y,Fs_z,Ff_x,Ff_y,Ff_z,
n_bonds,n_broken_v,n_broken_h,w_f,w_p
```

Далее по порядку (номера прежних столбцов не сдвигаются):

1. Пакет из $\ge 2$ слоёв: `n_broken_v_p1..K`, `n_broken_h_p1..K`, `w_f_p1..K`.
2. Модули реакций зажимов: `Fs_norm,Fs_abs_sum,Fs_abs_max,Ff_norm,Ff_abs_sum,Ff_abs_max`
   — $\lvert\sum F\rvert$ векторной суммы, $\sum\lvert F_k\rvert$ и $\max\lvert F_k\rvert$
   по зажатым частицам (`Fs_*` — пластина/Solid, `Ff_*` — нити).
3. Режим проекции:
   `proj_n_ff,proj_n_fs,proj_max_ck,proj_max_cs,proj_resid,proj_resid_rel,proj_warm_frac,W_proj_pos_max,proj_resid_load,proj_resid_post,proj_warm_scale_min`.
   `W_proj_pos_max` — положительная работа стадии проекции (теоретически 0). После
   рестарта из состояния монитора старого формата `proj_resid_load` и/или
   `proj_resid_post,proj_warm_scale_min` не пишутся, чтобы шапка согласовалась со строками.
4. Со стёжкой или мозаикой (последними): `n_stitch,n_stitch_broken,F_stitch_max,F_stitch_sum,slip_mean,gap_mean`.
   Стежки исключены из `n_bonds`; у мозаики `slip_mean`/`gap_mean` усредняются только по
   сквозным стежкам.

### `scene_id.txt`

Пишется при `--checkpoint` (без `--restart`) одной строкой: `config`, `plies`, `ball_D`,
`ball_m`, `ball_dx`; при `--fabric-*` — `fabric_n fabric_ce fabric_ey`; при стёжке —
`stitch=grid|seam patch=|tile= lines= every= E= d= eps_break= n_stitch=`; при мозаике —
`random_cut= seed= layers_per_cut= gap= rot_max= step= E= d= eps_break= n_stitch=` (и
`field= clamp_edge=` при своём поле). При `--restart` и существующем файле строка должна
совпасть побайтно, иначе выстрел завершается ошибкой «сцена не совпадает с сохранённой»
(другое число слоёв, шар, ткань, стёжка, мозаика). Нет файла — проверка не делается.

## Рестарт

Состояние монитора (вердикт, максимумы, окна проекции, строки `monitor.csv`) хранится в
блоке `USER` чекпоинта (формат чекпоинта 6, старые версии не читаются; сам блок монитора имеет
сигнатуру `GUOMON03`, а `GUOMON01`/`GUOMON02` ещё читаются — недостающие поля проекции получают
значения по умолчанию, соответствующие столбцы `monitor.csv` не дописываются): после
`--restart` результат совпадает с непрерывным, `shot.txt` содержит `monitor_state=restored`.
`monitor_state=partial` — состояния в чекпоинте не было; монитор начинает с шага рестарта.
Состояние от другой сцены (частицы, бусины, слои, режим проекции) или с иным числом шагов
отвергается. Подробности — в [«Рестарт из чекпоинта»](../02-how-to/05-restart-checkpoints.md).
