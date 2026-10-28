# 2T base (две несклеенные ткани Twaron, Guo 2019), пластично-вязкая нить

СТАТУС: НЕ ЗАВЕРШЕНО. Сервер claude-compute был недоступен (ssh: Connection timed out), поэтому готов только выстрел
μ=0.1, V0=137 м/с. Выстрел μ=0.08, V0=130 и сравнительные фигуры строятся после возврата сервера (команды ниже).

Постановка: нить `twaron_plastic_visco_rate_s30.fib` (eps_break 0.06, E = 8.5e10, σ_u = 3.0 ГПа, σ_y = 2.04 ГПа,
H = 1e11, rate_model 2), уровень base, 2 слоя, 86 496 бусин, 86 088 связей, шар 3544 частицы, кадры каждые 20 мкс (41 шт.).
Слой 1 -- лицевой (со стороны шара), слой 2 -- тыльный; вертикальные нити (вдоль y) зажаты, горизонтальные (вдоль x) свободны.

## μ=0.1, V0=137 м/с (`anim/mu0.1_v137/`)
Исход ПРОБИТ: V_R = 23.28 м/с, e = 63.79 Дж, t_end = 803 мкс. Порвано: слой 1 -- 9 верт. / 1 гориз., слой 2 -- 8 / 0;
первый разрыв 260 мкс (прогиб 27.4 мм); max eps_p = 0.72 %, W_p = 1.69 Дж.

| Файл | Что |
|---|---|
| `vel_mu0.1_v137.gif`, `vel_*_frame{0,3,13,18,40}.png` | 3D (слой 1 оранжевый, слой 2 синий, шар серый), сечение yz по v_z (слой 1 -- кружки, слой 2 -- ромбы), v_z(t) с маркером |
| `anim_2t_mu0.1_v137.gif`, `anim_2t_*_frame*.png` | слои 1 и 2 рядом, вид сверху, цвет eps/eps_break (шкала 0..0.75), красные/чёрные крестики -- порванные связи |
| `anim_2t_epsp_mu0.1_v137.gif`, `anim_2t_epsp_*_frame*.png` | то же по eps_p (шкала 0..0.008) |
| `summary_mu0.1_v137.png` | v_z(t); прогиб w_f_p1, w_f_p2 (t); число порванных нитей по слоям (t) |
| `energy_channels/dissipation/stocks/total/transfer.png` | баланс энергии по energy_log.csv |

Снимки: frame0 -- до контакта (t=0), frame3 -- шар в плоскости ткани, frame13 -- первый разрыв (260 мкс), frame18 -- шар прошёл ткань
(верх шара ниже медианной z бусин в круге 2R_шара), frame40 -- конец (800 мкс).

Энергия (E0 = 65.69 Дж, KE шара): к концу E_mech = 6.97 Дж (10.6 %: шар 1.90, нить 4.98); трение 49.36 Дж (75.1 %), вязкость нити
7.39 (11.2 %), пластика 1.70 (2.6 %), XSPH 0.28, разрыв связей 0.11, невязка -0.025 Дж, clip_loss 0.15 Дж.

## Скрипты (`scripts/`)
- `animate_2t.py` -- подкоманды vel / strain / summary / static / compare (для 2 слоёв; 1T-скрипты не тронуты)
- `animate_ply.py` -- копия из 1T-эксперимента (read_layout с колонкой ply), `plot_energy_all.py` -- копия
```
python scripts/animate_2t.py vel     --run-dir R --out anim/mu0.1_v137/vel_mu0.1_v137.gif --snapshots auto --dpi 70
python scripts/animate_2t.py strain  --run-dir R --out anim/mu0.1_v137/anim_2t_mu0.1_v137.gif --what eps --vmax 0.75 --snapshots auto
python scripts/animate_2t.py strain  --run-dir R --out anim/mu0.1_v137/anim_2t_epsp_mu0.1_v137.gif --what epsp --vmax 0.008 --snapshots auto
python scripts/animate_2t.py summary --run-dir R --out anim/mu0.1_v137/summary_mu0.1_v137.png
python scripts/plot_energy_all.py R/energy_log.csv anim/mu0.1_v137 --title "..."
```
