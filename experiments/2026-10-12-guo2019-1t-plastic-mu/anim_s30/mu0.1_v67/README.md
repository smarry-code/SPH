# Выстрел 1T base, mu = 0.1, V0 = 67 м/с, пластично-вязкая нить (σ_u = 3.0 ГПа)

Исход: ПРОБИТ, V_R = 2.80 м/с, e = 15.68 Дж, порвано 7 верт. / 0 гориз., первый разрыв 521 мкс (прогиб 27.9 мм),
max ε_p = 0.77 %, W_p = 0.70 Дж, t_end = 663 мкс. Нить: `../../twaron_plastic_visco_rate_s30.fib`
(PLASTIC SY=2.04e9 H=100e9; RATE JC C=0.005; VISCO E1=8.5e9 TAU=20e-6 XI=0.05; FAILURE SIGMA_U=3.0e9), eps_break = 0.06.
Кадры (34 шт., каждые 20 мкс) рендерились из копии `~/window1T_s30/mu0.1_v67/1T_base_v67.0_r0/` на сервере (в репозиторий не входят).

| Файл | Что |
|---|---|
| `vel_mu0.1_v67.gif`, `vel_*_frame{0,17,26,33}.png` | 3D по \|v\|, сечение по v_z, график v_z(t) |
| `anim_1t_mu0.1_v67.gif`, `anim_1t_*_frame*.png` | 3D нити, цвет eps/eps_break (шкала 0..0.75), сечение |
| `anim_1t_epsp_mu0.1_v67.gif`, `anim_1t_epsp_*_frame*.png` | то же, цвет по eps_p (шкала 0..0.008) |
| `vz_compare_s30.svg/.png` | v_z(t): mu=0.1 при V0=67 и 73, mu=0.3 при V0=73 |
| `energy_channels/dissipation/stocks/total/transfer.png` | энергетический баланс по energy_log.csv |

Команды (из каталога эксперимента; `R` -- каталог прогона, `D67/D73a/D73b` -- каталоги с monitor.csv):

```
python scripts/animate_vel.py --run-dir $R --out anim_s30/mu0.1_v67/vel_mu0.1_v67.gif --max-frames 40 --snapshots auto
python scripts/animate_vel.py --monitor-dirs $D67 $D73a $D73b --static-out anim_s30/mu0.1_v67/vz_compare_s30.svg
python scripts/animate_1t.py --run-dir $R --out anim_s30/mu0.1_v67/anim_1t_mu0.1_v67.gif --max-frames 40 --snapshots auto --vmax 0.75 --title-prefix "пластично-вязкая нить, σ_u = 3.0 ГПа"
python scripts/animate_1t.py --run-dir $R --out anim_s30/mu0.1_v67/anim_1t_epsp_mu0.1_v67.gif --max-frames 40 --snapshots auto --color-by eps_p --vmax 0.008 --title-prefix "пластично-вязкая нить, σ_u = 3.0 ГПа"
python scripts/plot_energy_all.py $R/energy_log.csv anim_s30/mu0.1_v67 --title "1T base, mu=0.1, V0=67 м/с, пластично-вязкая нить"
```

Правки скриптов: в `animate_vel.py` заголовок «пластично-вязкая нить, σ_u = 3.0 ГПа»; в `animate_1t.py` подпись шкалы eps_break = 0.06.
