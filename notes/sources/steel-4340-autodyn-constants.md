# STEEL 4340 — константы из библиотеки Autodyn

Источник: скриншоты окна Material Data Input в Autodyn (предоставила Mary, 27.10.2026 и 15.11.2026).
Значения переписаны с экрана; единицы — как в Autodyn (kPa, г/см³, К).
Эти значения имеют приоритет над литературными и над таблицей в `jc-implementation-plan.md` (§3.2), если они расходятся.

## Общие
| Параметр | Значение | Единицы Autodyn | СИ |
|---|---|---|---|
| Name | STEEL 4340 | | |
| Reference Density (ρ0) | 7.830 | г/см³ | 7830 кг/м³ |

## EOS: Linear
| Параметр | Значение | Единицы Autodyn | СИ |
|---|---|---|---|
| Bulk Modulus (K) | 1.590000e+08 | kPa | 159.0 ГПа |
| Reference Temperature (T_r) | 300.0 | К | 300 К |
| Specific Heat (c_p) | 477.000092 | Дж/(кг·К) | 477 Дж/(кг·К) |
| Thermal Conductivity | 0.0 | Дж/(м·К·с) | 0 |

Linear EOS: p = K·μ, μ = ρ/ρ0 − 1, скорость звука c = √(K/ρ0).

## Strength: Johnson Cook
| Параметр Autodyn | Обозначение | Значение | Единицы Autodyn | СИ |
|---|---|---|---|---|
| Shear Modulus | G | 8.180001e+07 | kPa | 81.8 ГПа |
| Yield Stress | A | 7.920000e+05 | kPa | 792 МПа |
| Hardening Constant | B | 5.100000e+05 | kPa | 510 МПа |
| Hardening Exponent | n | 0.26 | — | 0.26 |
| Strain Rate Constant | C | 0.014 | — | 0.014 |
| Thermal Softening Exponent | m | 1.03 | — | 1.03 |
| Melting Temperature | T_m | 1.793000e+03 | К | 1793 К |
| Ref. Strain Rate | ε̇_0 | 1.0 | 1/с | 1 с⁻¹ |
| Strain Rate Correction | | 1st Order | | |

Формула: Y = [A + B·ε_p^n] · [1 + C·ln(ε̇/ε̇_0)] · [1 − T*^m], T* = (T − T_r)/(T_m − T_r).
«1st Order» означает линейную по ln ε̇ форму 1 + C·ln ε̇*.

## Failure: Johnson Cook (включено, в отличие от прежних скриншотов)
| Параметр Autodyn | Значение | Единицы |
|---|---|---|
| Damage Constant D1 | 0.05 | — |
| Damage Constant D2 | 3.44 | — |
| Damage Constant D3 | −2.12 | — |
| Damage Constant D4 | 0.002 | — |
| Damage Constant D5 | 0.61 | — |
| Melting Temperature (Failure) | 1.793000e+03 | К |
| Ref. Strain Rate (Failure) | 1.0 | 1/с |

Формула: ε_f = [D1 + D2·exp(D3·σ*)] · [1 + D4·ln ε̇*] · [1 + D5·T*], σ* = −p/σ_eq (p > 0 при сжатии), D = ΣΔε_p/ε_f, разрушение при D = 1.

## Erosion
None.

## Что изменилось относительно плана
- D1–D5 для 4340 в Autodyn заданы и **совпадают** со значениями Johnson & Cook 1985, использованными в плане (0.05; 3.44; −2.12; 0.002; 0.61). Пометку «не сверено» для 4340 можно снять.
- Melting Temperature в разделе Failure = 1793 К (T_m из прочности), значит `Tm_f` не нужно брать по умолчанию из прочности — задавать явно 1793.
- Ref. Strain Rate для разрушения = 1 с⁻¹ (`edot0_f = 1`).
- Остальные параметры совпадают с таблицей в плане (G = 81.8 ГПа, A, B, n, C, m, T_m, c_p = 477, T_r = 300, ρ0 = 7830, K = 159 ГПа).

## Таблица в виде кода (СИ) для агента
```
itype 12 STEEL 4340
rho0    = 7830.0        // кг/м^3
K       = 159.0e9       // Па (Linear EOS)
Tr      = 300.0         // К
cp      = 477.0         // Дж/(кг·К)
G       = 81.8e9        // Па
A       = 792.0e6       // Па
B       = 510.0e6       // Па
n       = 0.26
C       = 0.014
m       = 1.03
Tm      = 1793.0        // К
edot0   = 1.0           // 1/с
rate_form = FirstOrder
D1..D5  = 0.05, 3.44, -2.12, 0.002, 0.61
Tm_f    = 1793.0        // К
edot0_f = 1.0           // 1/с
```
