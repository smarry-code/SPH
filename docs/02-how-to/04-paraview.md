# Визуализация в ParaView

Кадры расчёта (`sph_output_*.txt`, `fiber_output_*.txt`) переводятся в формат VTK XML
скриптом [`tools/paraview/sph2vtk.py`](../../tools/paraview/sph2vtk.py), а сцена в ParaView
(нити трубками, частицы Solid, сечение, подписи монитора) строится скриптом
[`tools/paraview/fiber_scene.py`](../../tools/paraview/fiber_scene.py). Проверено на ParaView 6.2.

## 1. Перевод кадров

```bash
python tools/paraview/sph2vtk.py --run-dir <каталог прогона> --out-dir <каталог vtk>
```

Нужен только `numpy` (`pandas` ускоряет чтение). Вид кадров определяется сам:

| Вид | Что в каталоге прогона | Нити |
|---|---|---|
| `pair` | `sph_output_<N>.txt` + `fiber_output_<N>.txt` (`guo2019_shot`, демо нитей) | связи из блока `BONDS` |
| `combined` | только `sph_output_<N>.txt` с колонками `... fiber i j eps eps_pl tension alive` | связь `i -> j` из строки бусины с `alive = 1` |
| `solid` | только `sph_output_<N>.txt` ([`src/output.cpp`](../../src/output.cpp)) | нет |

Если рядом лежат `monitor.csv` и `guo2019_fabric_layout.txt`, они тоже используются: время кадра
и величины монитора (`v_z`, `w_f`, `n_broken_v`, ...) берутся из строки с номером шага кадра,
а направление нити и зажатые бусины — из раскладки.

Результат:

| Файл | Содержимое |
|---|---|
| `solid_<N>.vtp`, `solid.pvd` | частицы Solid (вершины); PointData — все колонки кадра, `velocity`, `speed` |
| `fibers_<N>.vtp`, `fibers.pvd` | бусины (точки) и связи (отрезки). PointData бусин: `velocity`, `speed`, `object_id`/`fiber`, `direction`, `clamped`. CellData связей: `eps`, `eps_ratio` = ε/ε_разр, `l0`, `l`, `l_break` (pair) или `eps`, `tension`, `k`, `fiber` (combined) |

В FieldData каждого кадра лежат `TIME`, `STEP` и строка монитора. Порванная связь просто
отсутствует в кадре. Время в `.pvd` физическое; в кадрах `combined`/`solid` без `monitor.csv`
времени нет, и в `.pvd` пишется номер шага.

Ключи: `--every k` (каждый k-й кадр), `--jobs n` (процессов), `--double` (Float64 вместо
Float32), `--eps-break ε` (для `combined`: даёт `eps_ratio = eps/ε`). Файлы сжаты zlib: кадр
`combined` на 105 тыс. частиц — 1–4 МБ вместо 18 МБ текста.

### Журналы для графиков и анимаций

Вместе с `.vtp` конвертер пишет текстовые журналы:

| Файл | Содержимое |
|---|---|
| `frames.csv` | строка на кадр: `frame`, `step`, `time` (`time_us`); по нитям — `fiber_n_bonds`, `fiber_n_broken` (накопительно; `_x`/`_y` по раскладке), `fiber_w_max` (наибольший прогиб от первого кадра), `fiber_eps_max`, `fiber_eps_ratio_max`, `fiber_n_ratio_gt_0.9`, `fiber_tension_max`, `fiber_KE`; по каждому материалу Solid (`solid<Nmat>_…`) — центр масс, его скорость, KE, `z_min/z_max`, `p_min/p_max`, `sigma_i_max`; строка `monitor.csv` на шаге кадра — с приставкой `mon_` |
| `breaks.csv` | строка на разорванную связь: кадр и время (`time_before` … `time`), бусины `bead_i`, `bead_j`, исходный номер нити `fiber`, `direction`, `eps`/`eps_ratio`/`tension` в последнем кадре до разрыва, середина связи `x, y, z` |
| `breaks_<N>.vtp`, `breaks.pvd` | точки разрывов, накопленные к кадру N (`break_time`, `break_frame`, `fiber`, `direction`) |
| `columns.txt` | описание всех колонок обоих журналов |

Набор колонок зависит от того, что есть в кадрах: KE нитей — только если в кадре есть масса
(`combined`), `eps_ratio` — в `pair` или с `--eps-break`. Разрыв определяется как связь, которая
была в предыдущем кадре серии и пропала в текущем, поэтому его время известно с точностью до
интервала между кадрами. Номер нити берётся из раскладки: `object_id` в `fiber_output` после
разрыва получает оторванный кусок.

Как использовать в ParaView:

- **График по кадрам:** *File → Open* → `frames.csv` → *Apply*; затем *Line Chart View*, в
  свойствах снять *Use Index For X Axis*, выбрать *X Array Name* = `time_us` и отметить нужные
  ряды (например, `solid3_vz_c` и `mon_v_z`, `fiber_w_max`, `fiber_n_broken`).
- **Плотный график по шагам:** `monitor.csv` из каталога прогона открывается так же
  (*X Array Name* = `t`); в нём строка на каждый шаг расчёта.
- **Точки разрывов в 3D:** открыть `breaks.pvd`, представление *Point Gaussian*, цвет по
  `break_time_us` — при проигрывании анимации маркеры появляются в момент разрыва. Или
  `breaks.csv` → фильтр *Table To Points* (X/Y/Z = `x`, `y`, `z`) — все разрывы сразу.
- **Метка текущего кадра на графике:** *Line Chart View* сам рисует вертикальную черту
  времени анимации (отключается флажком *Hide Time Marker*). Черта стоит в единицах времени
  `.pvd`, поэтому для неё ось X должна быть `time` (секунды), а не `time_us`.

## 2. Сцена

```bash
pvpython tools/paraview/fiber_scene.py --vtk-dir <каталог vtk> --zoom 9e-3 --state scene.pvsm
pvpython tools/paraview/fiber_scene.py --vtk-dir <каталог vtk> --zoom 9e-3 --out anim.mp4
pvpython tools/paraview/fiber_scene.py --vtk-dir <каталог vtk> --zoom 9e-3 --snapshots 0,-1 --out snap.png
```

`pvpython` лежит в `bin/` комплекта ParaView (например `C:\Programs\ParaView\bin\pvpython.exe`).
Слева — 3D: нити трубками, цвет — `eps_ratio` (шкала 0..1), а если его нет — `tension`;
частицы Solid показаны сферами Point Gaussian. Справа — сечение `yz` (или `xz`) по полосе
`±--band` в параллельной проекции. Сверху — время и величины монитора.

| Ключ | Смысл |
|---|---|
| `--state f.pvsm` | сохранить состояние; открыть в GUI: *File → Load State* |
| `--out f.mp4/.avi/.ogv/.png` | анимация (`.png` — серия кадров) |
| `--snapshots 0,10,-1` | только снимки указанных кадров (`-1` — последний) |
| `--zoom h` | окно ±h м вокруг центра ткани |
| `--section yz/xz/none`, `--band b` | сечение и полуширина его полосы (по умолчанию 1,5 мм) |
| `--tube-radius r` | радиус трубок (по умолчанию 0,3 медианной длины связи) |
| `--solid-array p` | цвет частиц Solid по массиву (по умолчанию серые) |
| `--color-range lo,hi` | шкала `eps_ratio` |
| `--elev`, `--azim` | ракурс, как `view_init` в matplotlib |
| `--size WxH`, `--fps` | размер кадра и частота анимации |
| `--font f.ttf` | шрифт с кириллицей |

!!! note "Кириллица"
    Встроенный в ParaView шрифт Arial не содержит кириллицы: русские буквы в подписях
    пропадают. Скрипт подключает системный шрифт (`C:/Windows/Fonts/arial.ttf` или DejaVu Sans
    на Linux); другой можно задать `--font`.

## Ручная работа в GUI

Файлы `solid.pvd` и `fibers.pvd` открываются и без скрипта (*File → Open*). Полезные фильтры:
**Tube** — нити трубками, **Threshold** по `eps_ratio` — связи у порога разрыва,
**Python Annotation** с выражением вида `'%.2f мкс' % (TIME[0]*1e6)` — подпись из FieldData,
**Plot Data Over Time** — величины монитора по кадрам.
