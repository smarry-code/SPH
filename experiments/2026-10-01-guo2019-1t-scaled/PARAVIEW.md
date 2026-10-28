# Визуализация выстрела 1T scaled (V0 = 76 м/с) в ParaView

Кадры, VTK и таблицы лежат в каталоге эксперимента, но **не отслеживаются git**
(`experiments/.gitignore`). Всё получается заново скриптом `scripts/paraview_scene.sh`.

| Каталог | Содержимое |
|---|---|
| `data/run/1T_scaled_v76.0_r0/` | исходные кадры `sph_output_*.txt`, `fiber_output_*.txt` (67 кадров, 91 МБ), `monitor.csv`, раскладка нитей |
| `data/vtk/` | `solid.pvd`, `fibers.pvd`, `breaks.pvd` (+ `.vtp` кадры), `scene.pvsm`, снимки `snap_frame*.png` |
| `data/csv/` | `frames.csv` (строка на кадр, 55 колонок), `breaks.csv` (8 разрывов), `monitor.csv` (строка на шаг), `columns.txt` (описание колонок) |

## 1. Готовая сцена

Нужен ParaView 6.x (проверено на 6.2).

- **В GUI:** *File → Load State* → `data/vtk/scene.pvsm`. Путь к файлам VTK в состоянии
  сохранён как есть; если папку перенесли, выберите *Choose File Names* → укажите `data/vtk`.
  Кнопка Play проигрывает 67 кадров (0–165 мкс).
- **Снимки и ролик из командной строки:**
  `scripts/paraview_scene.sh scene` — пересохраняет `scene.pvsm` и снимки кадров 0, 21, 33 и последнего;
  `scripts/paraview_scene.sh anim` — дополнительно `data/vtk/anim.mp4`.
  Если `pvpython` не в PATH: `PVPYTHON="C:/Programs/ParaView/bin/pvpython.exe" scripts/paraview_scene.sh scene`.

Слева 3D: нити трубками, цвет — ε/ε_разр (0..1), частицы Solid (ударник) серые. Справа —
сечение `yz` по полосе x = 0 ± 1,5 мм. Сверху — время, v_z, прогиб, число порванных нитей.

## 2. Самостоятельная работа с данными

Откройте в ParaView (*File → Open*): `fibers.pvd`, `solid.pvd`, `breaks.pvd` из `data/vtk/`.

- **Нити:** `fibers.pvd` → фильтр *Tube* → раскраска по `eps_ratio` (CellData связей).
  *Threshold* по `eps_ratio` > 0,9 — связи у порога разрыва.
- **Ударник:** `solid.pvd` → представление *Point Gaussian*.
- **Разрывы в 3D:** `breaks.pvd` → *Point Gaussian*, цвет по `break_time_us`; маркеры появляются
  в момент разрыва. Либо `data/csv/breaks.csv` → *Table To Points* (X/Y/Z = `x`, `y`, `z`).
- **Графики по кадрам:** `data/csv/frames.csv` → *Apply* → *Line Chart View*, снять
  *Use Index For X Axis*, *X Array Name* = `time_us`; ряды: `solid3_vz_c`, `fiber_w_max`,
  `fiber_n_broken`, `fiber_eps_ratio_max`, `mon_*`.
- **Плотный график по шагам:** `data/csv/monitor.csv`, *X Array Name* = `t`.
- Описание всех колонок — `data/csv/columns.txt`.
- Подпись времени: *Python Annotation* с выражением `'%.2f мкс' % (TIME[0]*1e6)`.

## 3. Пересчёт с нуля

```bash
scripts/paraview_scene.sh convert   # кадры -> VTK + CSV (нужны python+numpy)
scripts/paraview_scene.sh scene     # сцена и снимки (нужен pvpython)
GUO_SHOT=<путь>/guo2019_shot scripts/paraview_scene.sh all   # + сам расчёт (~25 мин, 8 потоков) и mp4
```

Конвертер и сцена: `tools/paraview/sph2vtk.py`, `tools/paraview/fiber_scene.py`
(подробно — `docs/02-how-to/04-paraview.md`). Расчёт с OpenMP побитово не воспроизводится,
поэтому числа могут чуть отличаться от `README.md`.
