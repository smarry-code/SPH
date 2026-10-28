# Выстрел 1T base, mu = 0.1, V0 = 71 м/с, пластично-вязкая нить (σ_u = 3.0 ГПа)

Исход: ПРОБИТ, V_R = 9.15 м/с, e = 17.35 Дж, порвано 7 верт. / 0 гориз., первый разрыв 481 мкс (прогиб 27.9 мм),
max ε_p = 0.76 %, W_p = 0.73 Дж, t_end = 926 мкс. Нить: `../../twaron_plastic_visco_rate_s30.fib`.
Кадры (47 шт., каждые 20 мкс) -- из `~/shot1T_s30_mu0.1_v71/mu0.1_v71/1T_base_v71.0_r0/` на сервере (в репозиторий не входят).

| Файл | Что |
|---|---|
| `vel_mu0.1_v71.gif`, `vel_*_frame{0,23,24,46}.png` | 3D по \|v\|, сечение по v_z, график v_z(t) |
| `anim_1t_mu0.1_v71.gif`, `anim_1t_*_frame{0,25,46}.png` | 3D нити, цвет eps/eps_break (шкала 0..0.75), сечение |
| `anim_1t_epsp_mu0.1_v71.gif`, `anim_1t_epsp_*_frame*.png` | то же, цвет по eps_p (шкала 0..0.008) |
| `compare_65_71.gif` | сравнение 65 (не пробит, 31 кадр) и 71 (пробит); 31 различный кадр |
| `vz_compare_mu0.1.svg/.png` | v_z(t) шара, mu=0.1: 65 (не пробит), 67 (V_R=2.8), 71 (V_R=9.2), 73 (V_R=12.0) |
| `energy_channels/dissipation/stocks/total/transfer.png` | энергетический баланс по energy_log.csv |

Команды (из каталога эксперимента; `R` -- каталог прогона 71, `D65/D67/D73` -- каталоги прогонов с monitor.csv и shot.txt):

```
python scripts/animate_vel.py --run-dir $R --out anim_s30/mu0.1_v71/vel_mu0.1_v71.gif --max-frames 47 --snapshots auto
python scripts/animate_1t.py --run-dir $R --out anim_s30/mu0.1_v71/anim_1t_mu0.1_v71.gif --max-frames 47 --snapshots auto --vmax 0.75 --title-prefix "пластично-вязкая нить, σ_u = 3.0 ГПа"
python scripts/animate_1t.py --run-dir $R --out anim_s30/mu0.1_v71/anim_1t_epsp_mu0.1_v71.gif --max-frames 47 --snapshots auto --color-by eps_p --vmax 0.008 --title-prefix "пластично-вязкая нить, σ_u = 3.0 ГПа"
python scripts/animate_vel.py --compare --run-dir $D65 --run-dir $R --out anim_s30/mu0.1_v71/compare_65_71.gif --max-frames 40
python scripts/animate_vel.py --monitor-dirs $D65 $D67 $R $D73 --static-out anim_s30/mu0.1_v71/vz_compare_mu0.1.svg
python scripts/plot_energy_all.py $R/energy_log.csv anim_s30/mu0.1_v71 --title "1T base, mu=0.1, V0=71 м/с, пластично-вязкая нить"
```

## Что видно
- Шар идёт по v_z: -71 -> -50 м/с к ~390 мкс (все четыре выстрела почти совпадают), затем резкое торможение
  с разрывом 7 вертикальных нитей в центре около 481 мкс; после этого v_z ~ -12 м/с и дальше медленно до -9.2 м/с.
- На 460 мкс ткань прогнута «V»-образно (прогиб 27 мм), шар в вершине; на 920 мкс в ткани сквозное отверстие
  по центру, шар ниже ткани; нити разошлись по краям, ε_p сосредоточена в двух полосах вдоль линии разрыва (до 0.76 %).
- Энергия (E0 = 17.64 Дж): потеряно 15.6 Дж. Трение 12.93 (82 %), вязкость нити 1.97 (12 %), пластика 0.73 (5 %),
  XSPH 0.12, разрыв связей 0.04 Дж. Остаток к концу 2.04 Дж: ткань 1.72 (84 %), шар+пластины 0.29 Дж (14 %).
  Невязка баланса в energy_log: 0.145 Дж (~0.8 % E0), clip_loss 0.042 Дж.
