# Выходные файлы расчёта

При каждом сохранении кадра функция [`output()`](../../../src/output.cpp:16)
формирует два файла: текстовый кадр состояния и накопительный CSV координат.
Также ведётся журнал расчёта. Если в модели есть нити, на тех же шагах пишется
`fiber_output_NNNNN.txt` (см. [`.fib`](fib-file.md#выходной-файл-нитей-fiber_output_nnnnntxt)).
Кадры выводятся, когда `time >= next_output_time` (период `save_t`).

## `sph_output_NNNNN.txt` — кадр состояния частиц

Создаётся на каждом шаге сохранения. Имя формируется как
`sph_output_` + 5-значный номер шага с ведущими нулями + `.txt`:

```
sph_output_00001.txt
sph_output_00050.txt
```

Номер — счётчик шагов цикла (`current_ts`), который увеличивается в начале шага, поэтому первый
кадр (время 0, состояние после первого шага) имеет номер `00001`, а не `00000`. После выхода из
цикла `main()` дописывает ещё один кадр с номером `maxtimestep` (3000). У `guo2019_shot` пробный
кадр конечного состояния (`--probe-frame`) имеет номер `99999`. Переменная окружения
`SPH_FRAME_PRECISION=<n>` задаёт число значащих цифр кадра (отладка расхождения траекторий;
действует и на `fiber_output_*.txt`).

!!! note "Расширение файлов"
    В коде ([`output.cpp`](../../../src/output.cpp:23)) создаётся `.txt`. Структура колонок
    ниже едина.

### Структура

**Строка 1** — заголовок с именами колонок (зависит от размерности `dim`).

При `dim == 1`:

```text
N Nmat M x vx rho p u c sigma_i e_pl D hsml sxx syy szz sxy syz sxz T edot
```

При `dim == 2` добавляются колонки `y`, `vy`; при `dim == 3` — `y z` и `vy vz`.

**Далее** — по одной строке на частицу (`ntotal` строк).

### Колонки

| Колонка | Обозначение | Единицы | Описание |
|---|---|---|---|
| `N` | — | — | Порядковый номер частицы (с 1) |
| `Nmat` | `itype` | — | Индекс материала |
| `M` | $m$ | кг | Масса |
| `x`, `y`, `z` | $x_i$ | м | Координаты (по `dim`) |
| `vx`, `vy`, `vz` | $v_i$ | м/с | Скорости (по `dim`) |
| `rho` | $\rho$ | кг/м³ | Плотность |
| `p` | $p$ | Па | Давление |
| `u` | $u$ | Дж/кг | Удельная внутренняя энергия |
| `c` | $c$ | м/с | Скорость звука |
| `sigma_i` | $\sigma_i$ | Па | Интенсивность напряжений (фон Мизес) |
| `e_pl` | $\varepsilon^{pl}$ | — | Накопленная эквивалентная пластическая деформация (JH-2: сумма приращений радиального возврата) |
| `D` | $D$ | — | Параметр повреждённости |
| `hsml` | $h$ | м | Длина сглаживания |
| `sxx`…`sxz` | $s_{ij}$ | Па | Компоненты девиатора напряжений (порядок `sxx syy szz sxy syz sxz`) |
| `T` | $T$ | К | Температура (модель Джонсона–Кука; у остальных материалов 0) |
| `edot` | $\dot{\varepsilon}$ | 1/с | Скорость деформации (Джонсон–Кук; у остальных 0) |

Интенсивность напряжений вычисляется как:

$$
\sigma_i = \frac{1}{\sqrt{2}} \sqrt{(s_{xx}-s_{yy})^2 + (s_{xx}-s_{zz})^2 + (s_{yy}-s_{zz})^2 + 6(s_{xy}^2 + s_{yz}^2 + s_{xz}^2)}
$$

(см. [`ParticleSystem::vuvod_index()`](../../../include/particle_system.hpp:139)).

## `sph_gp.csv` — накопительный дамп кинематики

Открывается в режиме дозаписи (`std::ios::app`) на каждом кадре. Каждая строка:

```text
<time> <x0> ... <x_{dim-1}> <vx0> ... <vx_{dim-1}>
```

Одна строка соответствует одной частице в данный момент времени. Поле `time`
позволяет реконструировать траектории отдельных частиц.

!!! note "В `.gitignore`"
    `sph_gp.csv` и `sph_output_*.txt` исключены из контроля версий как
    генерируемые артефакты.

## `log.txt` — журнал расчёта

Текстовый файл с UTF-8 BOM (создаётся в [`run_sph()`](../../../src/main.cpp:68)).
Содержит диагностические сообщения, в том числе:

- сообщения об ошибках инициализации;
- строку мониторинга на каждом шаге цикла (дозапись из
  [`time_integration.cpp`](../../../src/time_integration.cpp:983)):
  `<time> <v> <itype>` — время, осевая компонента скорости мониторной частицы
  (`vx[dim-1]`) и её материал.

Исключён из контроля версий ([`.gitignore`](../../../.gitignore:67)).

## Другие журналы

- `energy_log.csv` (имя — `energy_log_file`, период — `energy_log_every`) — CSV
  энергетического баланса, только в сборке `-DSPH_ENERGY_LOG=ON`
  ([`energy_log.hpp`](../../../include/energy_log.hpp:426)). Колонки:
  `step,time,dt,KE_solid,KE_fiber,U_int,U_shear_el,U_bond,U_hertz_ff,U_hertz_fs,E_mech,D_friction,D_xsph,D_contact_ss,D_break,S_u_clamp,Q_av,W_plastic,residual,W_p_fiber,W_v_fiber,D_break_contact,clip_loss,D_contact_fiber`,
  затем диагностические каналы мощности `P_*` (по механизмам; вариант с массой кольца
  и вариант `_flat`) и `closure`. При рестарте файл дописывается без заголовка.
- `monitor.csv` — построчный монитор сценария Guo 2019 (тестовый харнесс,
  [`guo2019_monitor.hpp`](../../../tests/guo2019_monitor.hpp:385), `open_csv_fresh`):
  `step,t,dt,z_c,v_z,KE_s,KE_p,KE_f,U_b,U_int,E_tot,Fs_x,Fs_y,Fs_z,Ff_x,Ff_y,Ff_z,n_bonds,n_broken_v,n_broken_h,w_f,w_p`
  (для многослойной ткани добавляются `n_broken_v_pK`, `n_broken_h_pK`, `w_f_pK`; затем
  `Fs_norm,Fs_abs_sum,Fs_abs_max,Ff_norm,Ff_abs_sum,Ff_abs_max` — модули реакций зажимов;
  в режиме проекции нормального контакта `proj_n_ff,…,W_proj_pos_max,proj_resid_load,proj_resid_post,proj_warm_scale_min`;
  со стёжкой или мозаикой последними `n_stitch,n_stitch_broken,F_stitch_max,F_stitch_sum,slip_mean,gap_mean`;
  полный список и `shot.txt`, `scene_id.txt` — в [справочнике `guo2019_shot`](../guo2019-shot.md));
  рядом пишется `guo2019_fabric_layout.txt` — раскладка нитей.
- `checkpoint.ckpt` и `checkpoint.ckpt.prev` — бинарный чекпоинт (сигнатура `SPHCKPT`,
  версия формата 6 — единственная, старые не читаются, [`checkpoint.cpp`](../../../src/checkpoint.cpp:35));
  см. [рестарт](../../02-how-to/05-restart-checkpoints.md).
- `contact_stats.csv`, `contact_stats_summary.txt` — диагностика контактов нитей, только
  при `-DSPH_CONTACT_STATS=ON` (ключи `contact_stats_*`).

## Ручная визуализация

Файлы `sph_output_*.txt` — текстовые таблицы с разделителем-пробелом, их можно
загрузить в ParaView, gnuplot, Python (`numpy.loadtxt`) или Excel. Например,
профиль плотности вдоль оси $x$ (одномерный случай):

```python
import numpy as np
data = np.loadtxt("sph_output_00050.txt", skiprows=1)
x, rho = data[:, 3], data[:, 5]
```
