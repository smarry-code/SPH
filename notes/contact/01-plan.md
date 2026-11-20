# План внедрения контактных сил «частица–частица» (типы 1–3 по `notes/sources/sph-contact-article.md`)

> **Статус (2026-11-14):** реализовано, с отличиями. Формулы типов 1–3 — `include/contact_force.hpp`, ключи — `src/input.cpp` (`contact_force_type`, `contact_penalty_scale`, `contact_impulse_gamma`). Оркестратор в `src/ext_force.cpp` (шаги 3–4 плана) больше не продовый: контакт Solid–Solid считает `include/interactions/solid_contact.hpp`, а `ext_force.cpp` линкуется только в `test_contact`; критерий пары «разные материалы» (D2=A) заменён на `body[i] != body[j]`. Номера строк и пути плана — на `a4147b9`. Актуальное описание — [Контакт Solid–Solid](../../docs/04-explanation/models/03-solid-contact.md), [ADR-007](../../docs/04-explanation/decisions/adr-007-contact-force-types.md).

**Адресат:** следующий запуск Claude Code в этом репозитории.
**Основание:** `notes/sources/sph-contact-article.md` (разд. 2.1, 2.3), проверено по первоисточнику.
**Не основание:** `notes/contact/report-sph-contact.md` — его раздел «что уже есть в кодовой базе»
не соответствует репозиторию, см. [00-audit.md](00-audit.md). Ниже — самодостаточный план.
**Перед началом:** прочитать [02-decisions.md](02-decisions.md). Если владелец не дал
ответов, применять рекомендованные варианты (D1=A, D2=A, D3=A, D4=A, D5=A, D6=A) и
фиксировать их в ADR-007 как принятые допущения.

Все ссылки на строки проверены на HEAD `a4147b9` (ветка `master`).

---

## 0. Что именно внедряется

Три штрафные контактные силы из `notes/sources/sph-contact-article.md` §2.3, обнаружение проникания по линии
центров (2.1). Вне объёма (фиксируется в ADR как отложенные фазы):

- (2.2) «средняя нормаль» — нужны нормали к свободной поверхности, которых в проекте нет;
- §2.4 «свободная поверхность + контакт» — нужны граничные виртуальные частицы;
- модель Wang, трение Кулона, контактное тепловыделение — в статье их нет;
- доля гашения $\gamma < 1$ для типа 3.

### Формулы (после решения D1=A)

Геометрия пары, (2.1):

$$p = \frac{h_i + h_j}{2} - r_{ij},\qquad
\hat{\mathbf n} = \frac{\mathbf x_i - \mathbf x_j}{r_{ij}},\qquad
\dot p = -\frac{\mathbf x_{ij}\cdot\mathbf v_{ij}}{r_{ij}},\qquad
R_k = \frac{h_k}{2}$$

$\dot p > 0$ означает сближение. Сила действует только при $p > 0$.
Приведённая масса $m_{ij} = \dfrac{m_im_j}{m_i+m_j}$ (замена $\varrho R^3 \to$ `sys.mass`,
обоснование — Report §2.3, проверено и принимается).

**Тип 1 — «пинбол» Белицко, (2.3)–(2.5):**

$$F = K_p\min(F_1,F_2),\qquad
F_1 = \begin{cases} m_{ij}\dot p/\Delta t_{\mathrm{ref}}, & \dot p > 0\\ 0,& \dot p \le 0\end{cases},\qquad
F_2 = \frac{G_iG_j}{G_i+G_j}\sqrt{\frac{R_iR_j}{R_i+R_j}}\,p^{3/2}$$

Откат: если $G_i$ или $G_j$ неизвестен (0), ветка $F_2$ **исключается из `min`**
(иначе $\min(F_1,0)=0$ обнулит контакт — ошибка в Report §2.8), остаётся $F = K_p F_1$.

**Тип 2 — Свигл, (2.6)–(2.8)** (опечатка первоисточника: вторая формула помечена $F_1$,
по смыслу (2.6) это $F_2$):

$$F = \begin{cases} K_p\max(F_1,F_2), & \dot p > 0\\ 0, & \dot p \le 0\end{cases},\qquad
F_1 = \frac{m_{ij}p}{\Delta t_{\mathrm{ref}}^2},\qquad
F_2 = \frac{2 m_{ij} c_{\mathrm{pair}} p}{\Delta t_{\mathrm{ref}}\, r_{ij}},\qquad
c_{\mathrm{pair}} = \min(c_i,c_j)$$

($\Delta x = |x_j - x_i| = r_{ij}$ по (2.8); обобщение $m\to m_{ij}$, $c\to\min(c_i,c_j)$.)

**Тип 3 — импульсный, (2.10):**

$$F = \begin{cases} \dfrac{\lVert \mathbf v_j - \mathbf v_i\rVert}{\Delta t_{\mathrm{applied}}}\,m_{ij}, & \dot p > 0\\ 0,& \dot p \le 0\end{cases}$$

Здесь **применяемый** шаг: импульс $F\Delta t = m_{ij}\lVert v_{ij}\rVert$ от шага не
зависит, обратной связи с ограничителем нет (ограничитель для типа 3 и не ставится).

### $\Delta t_{\mathrm{ref}}$ — локальный CFL-шаг пары (в проекте отсутствует, вводится здесь)

$$\Delta t_{\mathrm{ref}} = \min_{k\in\{i,j\}} \frac{\mathtt{cfg.CFL}\cdot h_k}{c_k + \lVert\mathbf v_k\rVert}$$

Это поформульно тот же критерий, что в [`compute_cfl_timestep`](../../src/time_integration.cpp#L13-L27),
суженный на пару. Отсюда **$\Delta t_{\mathrm{applied}} \le \Delta t_{\mathrm{ref}}$ всегда**
(глобальный шаг — минимум по всем частицам). Жёсткость строится от $\Delta t_{\mathrm{ref}}$,
а не от применяемого шага — иначе связка «жёсткость ← шаг ← ограничитель» даёт коллапс $\Delta t$.

### Ограничитель шага

Линеаризация по прониканию, $\omega = \sqrt{k_{\mathrm{eff}}/m_{ij}}$, запас
$\Delta t_{osc} = 0{,}8/\omega$:

| Тип | Активная ветка | $k_{\mathrm{eff}} = \partial F/\partial p$ | $\Delta t_{osc}$ |
|---|---|---|---|
| 1 | $F_2$ (Герц) | $1{,}5\,K_p F_2/p$ | $0{,}8/\omega$ |
| 1 | $F_1$ (вязкая, $F = c_d\dot p$, $c_d = K_p m_{ij}/\Delta t_{\mathrm{ref}}$) | позиционной жёсткости нет | $0{,}8\cdot 2\Delta t_{\mathrm{ref}}/K_p$ (предел явной схемы для демпфера $\Delta t < 2m/c_d$) |
| 2 | любая ($F$ линейна по $p$) | $F/p$ | $0{,}8/\omega$ |
| 3 | — | — | не ставится |

Оценка: тип 2, ветка $F_1$ → $\Delta t_{osc} = 0{,}8\,\Delta t_{\mathrm{ref}}/\sqrt{K_p}$;
ветка $F_2$ при $r\approx h$, `CFL = 0,1`, $K_p = 1$ → $\approx 0{,}18\,h/c$ против
CFL-шага $0{,}1\,h/c$ — не связывает. Ограничитель начинает резать шаг при $K_p \gtrsim 1{,}6$.

---

## Шаг 0. Починить baseline (решение D5=A)

**Цель:** проект собирается и `ctest` зелёный до начала работы над контактом.

1. [`include/param.hpp:12`](../../include/param.hpp#L12): `Particle::Particle();` → `Particle();`
   (лишняя квалификация; Clang — ошибка). Структура `Particle` — мёртвый код, первичное
   хранилище — `ParticleSystem`; удалять её в этом шаге **не надо** (вне объёма).
2. [`src/main.cpp:13`](../../src/main.cpp#L13): обернуть `#include <windows.h>` в `#ifdef _WIN32`
   — единственный незащищённый включатель, в остальных модулях он уже под охраной.
3. [`CMakeLists.txt`](../../CMakeLists.txt): ветка OpenMP для AppleClang. Сейчас
   [строки 49-52, 70-73, 104-107](../../CMakeLists.txt#L104-L107) выдают голый `-fopenmp`,
   который Apple Clang отвергает (`unsupported option '-fopenmp'`). Нужна ветка:
   `-Xpreprocessor -fopenmp`, `include`/`lib` из `brew --prefix libomp`, линковка `-lomp`.
   Windows- и GNU-пути **не менять**.

**Проверка (выполнялась, воспроизводима):**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8          # ожидается exit 0, цели SPH, test_*, bench_color
ctest --test-dir build --output-on-failure
```

Подтверждено: после этих трёх правок все цели, включая `SPH`, собираются (exit 0).
`libomp` в системе присутствует (`/opt/homebrew/opt/libomp`).

**Коммит:** отдельный, `fix(build): переносимость сборки на Apple Clang` — не смешивать
с контактом.

---

## Шаг 1. Конфигурация и доступ к модулю сдвига

**Цель:** ключи и материальные данные доступны коду силы.

### 1.1 `include/simulation_config.hpp`

Добавить в структуру (рядом с блоком «SPH-метод и алгоритмы», [строки 29-35](../../include/simulation_config.hpp#L29-L35)):

```cpp
// --- Контакт «частица-частица» (notes/sources/sph-contact-article.md, ADR-007) ---
// 0 - выключен (по умолчанию), 1 - пинбол Белицко, 2 - Свигл, 3 - импульсный.
int    contact_force_type   = 0;
double contact_penalty_scale = 1.0;   // масштабный коэффициент Kp (типы 1-2)
```

и **обязательно** продублировать в [`apply_defaults()`](../../include/simulation_config.hpp#L68-L77)
— по контракту, зафиксированному в комментарии к функции.

### 1.2 `src/input.cpp`

В `parse_sph_config` дописать два `else if` в конец цепочки, после ветки `print_step`
([input.cpp:327-330](../../src/input.cpp#L327-L330)), по образцу соседей:

```cpp
else if (key == "contact_force_type")   { cfg.contact_force_type   = std::stoi(value); }
else if (key == "contact_penalty_scale"){ cfg.contact_penalty_scale = std::stod(value); }
```

Валидация значения `contact_force_type` вне 0..3 — сообщение в `std::cerr` и откат к 0.

### 1.3 `Initial_Conditions.txt`

Дописать в конец (стиль комментариев — `*`, как в файле):

```text
* --- Контакт «частица-частица» (notes/sources/sph-contact-article.md, ADR-007) ---
* 0 - выключен, 1 - пинбол Белицко, 2 - Свигл, 3 - импульсный
contact_force_type = 0
* Масштабный коэффициент штрафной силы Kp (типы 1-2)
contact_penalty_scale = 1.0
```

### 1.4 `include/materials.hpp`

Добавить после структуры `Materials` ([materials.hpp:40](../../include/materials.hpp#L40)):

```cpp
// Модуль сдвига по номеру материала (поле itype частицы).
// Реестр разложен по двум массивам; диспетчер воспроизводит логику
// int_force.cpp:636-652. Возвращает 0.0 для материалов без G
// (газ/вода, itype 1-2, виртуальные частицы itype < 0) — вызывающая
// сторона обязана трактовать 0 как «G неизвестен», а не как нулевую жёсткость.
inline double shear_modulus_by_itype(const Materials& m, int itype)
{
    if (itype >= 3 && itype <= 9)
    {
        const size_t k = static_cast<size_t>(itype) - Materials::n_Me;
        return k < m.G.size() ? m.G[k] : 0.0;
    }
    if (itype == 10 || itype == 11)
    {
        const size_t k = static_cast<size_t>(itype - 10);
        return k < m.const_JH_2.size() ? m.const_JH_2[k][1] : 0.0;
    }
    return 0.0;
}
```

`src/int_force.cpp` **не трогать** — его паттерн сохраняется как есть (вне объёма).
Внимание: `itype 12` и поле `jc`, упомянутые в Report, в реестре отсутствуют.

**Проверка шага 1:** сборка зелёная; значения `shear_modulus_by_itype` для itype 4, 10, 11
совпадают с реестром (`26e9`, `197e9`, `1.165e11`); для itype 1, 2, 12, −3 возвращается 0.
Эти проверки войдут в `test_contact.cpp` на шаге 2.

---

## Шаг 2. Чистая математика контакта (новый заголовок) + её тесты

**Цель:** формулы типов 1–3 реализованы, dim-независимы и покрыты аналитически —
до всякой интеграции в конвейер.

### 2.1 Новый файл `include/contact_force.hpp`

Header-only, `inline`, без зависимости от `dim` и от `ParticleSystem` — по образцу
`flat_grid.hpp` / `pair_graph.hpp` (тестируются напрямую). Все величины скалярные:
направление силы применяет оркестратор.

```cpp
#pragma once
// contact_force.hpp
// Штрафные контактные силы «частица-частица» (notes/sources/sph-contact-article.md, разд. 2.3; ADR-007).
// Чистые функции: вход — скалярные характеристики пары, выход — модуль
// нормальной силы и ограничитель шага. Направление и раскладка по осям —
// в оркестраторе (ext_force.cpp), поэтому файл не зависит от dim.
#include <algorithm>
#include <cmath>

namespace contact
{

// Значение «ограничитель не установлен» (совместимо с min() в time_integration).
inline constexpr double kNoLimit = 1e30;
// Запас устойчивости явной схемы: dt*omega < 2, берём 0.8.
inline constexpr double kOscSafety = 0.8;

enum ForceType : int { Off = 0, Pinball = 1, Swegle = 2, Impulse = 3 };

// Геометрия пары по (2.1). Заполняется оркестратором.
struct PairGeom
{
    double r     = 0.0;  // |x_i - x_j|, > 0
    double p     = 0.0;  // проникание 0.5*(h_i+h_j) - r; контакт при p > 0
    double pdot  = 0.0;  // скорость проникания, > 0 при сближении
    double v_rel = 0.0;  // ||v_i - v_j||, полная относительная скорость (2.10)
};

// Свойства пары, не зависящие от геометрии.
struct PairProps
{
    double m_red  = 0.0; // приведённая масса m_i m_j/(m_i+m_j)
    double c_pair = 0.0; // min(c_i, c_j)
    double g_eff  = 0.0; // G_i G_j/(G_i+G_j); 0 => G неизвестен (ветка F2 отключается)
    double r_eff  = 0.0; // R_i R_j/(R_i+R_j), R = h/2
    double dt_ref = 0.0; // локальный CFL-шаг пары; НЕ применяемый шаг
};

struct Result
{
    double fn       = 0.0;      // модуль нормальной силы, Н (>= 0)
    double dt_limit = kNoLimit; // ограничитель шага; kNoLimit => не ограничивает
};

// Тип 1, (2.3)-(2.5). Множитель p из буквальной (2.3) опущен: F1 и F2 —
// силы, и p делал бы размерность Н*м (решение D1 в notes/contact/02-decisions.md).
inline Result pinball_force(const PairGeom& g, const PairProps& pr, double kp)
{
    Result out;
    if (g.p <= 0.0 || pr.m_red <= 0.0 || pr.dt_ref <= 0.0 || kp <= 0.0) return out;

    const double f1 = (g.pdot > 0.0) ? pr.m_red * g.pdot / pr.dt_ref : 0.0;
    const bool   hertz_known = (pr.g_eff > 0.0 && pr.r_eff > 0.0);
    const double f2 = hertz_known
                    ? pr.g_eff * std::sqrt(pr.r_eff) * g.p * std::sqrt(g.p)
                    : 0.0;

    // Откат при неизвестном G: F2 исключается из min, а не обнуляет силу.
    const bool hertz_active = hertz_known && (f2 < f1);
    const double f = hertz_known ? std::min(f1, f2) : f1;
    out.fn = kp * f;
    if (out.fn <= 0.0) return out;

    if (hertz_active)
    {
        // Герц: F ~ p^{3/2} => k_eff = dF/dp = 1.5 * F / p.
        const double k_eff = 1.5 * out.fn / g.p;
        out.dt_limit = kOscSafety / std::sqrt(k_eff / pr.m_red);
    }
    else
    {
        // Вязкая ветка: F = c_d * pdot, c_d = kp*m_red/dt_ref.
        // Предел явной схемы для демпфера: dt < 2*m_red/c_d = 2*dt_ref/kp.
        out.dt_limit = kOscSafety * 2.0 * pr.dt_ref / kp;
    }
    return out;
}

// Тип 2, (2.6)-(2.8). Во второй формуле первоисточника опечатка: она помечена
// F1, но по (2.6) это F2.
inline Result swegle_force(const PairGeom& g, const PairProps& pr, double kp)
{
    Result out;
    if (g.p <= 0.0 || g.pdot <= 0.0 || pr.m_red <= 0.0 ||
        pr.dt_ref <= 0.0 || g.r <= 0.0 || kp <= 0.0) return out;

    const double f1 = pr.m_red * g.p / (pr.dt_ref * pr.dt_ref);
    const double f2 = 2.0 * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r);
    out.fn = kp * std::max(f1, f2);
    if (out.fn <= 0.0) return out;

    const double k_eff = out.fn / g.p;   // F линейна по p
    out.dt_limit = kOscSafety / std::sqrt(k_eff / pr.m_red);
    return out;
}

// Тип 3, (2.10). От ПРИМЕНЯЕМОГО шага: импульс F*dt = m_red*||v_ij|| от шага
// не зависит, поэтому ограничитель не ставится и обратной связи нет.
inline Result impulse_force(const PairGeom& g, const PairProps& pr, double dt_applied)
{
    Result out;
    if (g.p <= 0.0 || g.pdot <= 0.0 || pr.m_red <= 0.0 || dt_applied <= 0.0) return out;
    out.fn = pr.m_red * g.v_rel / dt_applied;   // dt_limit остаётся kNoLimit
    return out;
}

inline Result penalty_force(int type, const PairGeom& g, const PairProps& pr,
                            double kp, double dt_applied)
{
    switch (type)
    {
        case Pinball: return pinball_force(g, pr, kp);
        case Swegle:  return swegle_force(g, pr, kp);
        case Impulse: return impulse_force(g, pr, dt_applied);
        default:      return Result{};
    }
}

} // namespace contact
```

### 2.2 Новый файл `test_contact.cpp` (математическая часть)

Стиль — как в [`test_particle_system.cpp`](../../tests/test_particle_system.cpp#L12-L20):
`int main()`, счётчик `failures`, `std::printf`, без фреймворка. Сравнение чисел —
относительное, допуск `1e-12`.

| № | Проверка | Ожидание |
|---|---|---|
| M1 | `shear_modulus_by_itype` для itype 4, 10, 11 | `26e9`, `197e9`, `1.165e11` |
| M2 | `shear_modulus_by_itype` для itype 1, 2, 12, 99, −3 | `0.0` (без выхода за границы массивов) |
| M3 | Тип 1, $p>0$, $\dot p = 0$ | `fn == 0` (ветка $F_1$ нулевая, `min` даёт 0) |
| M4 | Тип 1, малое $\dot p$ (доминирует Герц) | `fn == kp*f2` с аналитикой; `dt_limit == 0.8/sqrt(1.5*fn/p/m_red)` |
| M5 | Тип 1, большое $\dot p$ (доминирует $F_1$) | `fn == kp*f1`; `dt_limit == 1.6*dt_ref/kp` |
| M6 | Тип 1, `g_eff == 0` при $\dot p > 0$ | `fn == kp*f1` (**не 0** — проверка отката) |
| M7 | Тип 1/2/3, `p <= 0` | `fn == 0`, `dt_limit == kNoLimit` |
| M8 | Тип 2, $\dot p>0$, обе ветки по очереди доминируют | `fn == kp*max(f1,f2)` с аналитикой; `dt_limit == 0.8*sqrt(m_red*p/fn)` |
| M9 | Тип 2, $\dot p<0$ | `fn == 0` |
| M10 | Тип 3, $\dot p>0$ | `fn == m_red*v_rel/dt` точно; `dt_limit == kNoLimit` |
| M11 | Тип 3, $\dot p<0$ | `fn == 0` |
| M12 | Тип 0 / неизвестный тип | `fn == 0`, `dt_limit == kNoLimit` |

### 2.3 `CMakeLists.txt`

Добавить цель по образцу `test_particle_system` ([CMakeLists.txt:55-60](../../CMakeLists.txt#L55-L60)):

```cmake
# --- Тест контактных сил (notes/sources/sph-contact-article.md, ADR-007): аналитика типов 1-3,
#     баланс импульса, ограничитель шага. Линкуется с src/ext_force.cpp,
#     чтобы покрыть и оркестратор.
add_executable(test_contact test_contact.cpp src/ext_force.cpp)
target_include_directories(test_contact PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_test(NAME contact_test COMMAND test_contact)
```

Имя цели — `test_contact`, имя теста — `contact_test` (согласовано со стилем остальных).

**Проверка шага 2:** `ctest --test-dir build -R contact_test --output-on-failure` — зелёный.
На этом шаге `src/ext_force.cpp` ещё не содержит контакта, но уже линкуется — это нормально.

---

## Шаг 3. Оркестратор в `ext_force` + ограничитель шага

**Цель:** контакт считается в конвейере и влияет на шаг.

### 3.1 Диагностика и ограничитель в `include/buffers.hpp`

Добавить в `WorkBuffers` ([buffers.hpp:25](../../include/buffers.hpp#L25)):

```cpp
    // --- Контакт (ADR-007). Пишутся в ext_force, читаются в time_integration.
    // ext_force обходит пары последовательно (OpenMP там нет), гонок нет.
    double min_dt_contact    = 1e30;  // ограничитель шага; 1e30 => контакт не ограничивает
    size_t contact_pair_count = 0;    // число пар с p > 0 за шаг
    double contact_max_delta  = 0.0;  // max(p) за шаг, м

    // Сброс контактной диагностики перед накоплением на текущем шаге.
    void reset_contact()
    {
        min_dt_contact = 1e30;
        contact_pair_count = 0;
        contact_max_delta = 0.0;
    }
```

Значение `1e30` держать согласованным с `contact::kNoLimit`
(добавить `static_assert` или ссылку в комментарии).

### 3.2 Сигнатура `ext_force`

Текущая ([ext_force.hpp:16](../../include/ext_force.hpp#L16)) не пропускает ни конфиг,
ни шаг, ни буферы. Новая:

```cpp
void ext_force(const SimulationConfig& cfg, ParticleSystem& sys, double dt,
               size_t ntotal, size_t niac,
               const std::vector<size_t>& pair_i, const std::vector<size_t>& pair_j,
               std::vector<std::vector<double>>& dvxdt, WorkBuffers& wb);
```

Порядок аргументов — по образцу `int_force`/`art_visc` (`cfg, sys, ...` в начале,
`wb` в конце). В `ext_force.hpp` добавить `#include "simulation_config.hpp"`,
`#include "buffers.hpp"`, `#include "contact_force.hpp"` и описать контракт типов 0–3
в комментарии.

Точка вызова — [`single_step.cpp:205`](../../src/single_step.cpp#L205); `cfg`, `dt` и `wb`
там уже в области видимости (параметры `single_step`).

### 3.3 Тело контактной ветки в `src/ext_force.cpp`

Мёртвый блок `is_contact_pair_` ([ext_force.cpp:113-141](../../src/ext_force.cpp#L113-L141))
**заменить** контактной веткой. Пристеночное отталкивание и самогравитацию не трогать.

Вспомогательные чистые функции — в анонимном namespace файла:

```cpp
namespace {

// Признак контактной пары (решение D2=A): разные материалы, обе частицы реальные.
// Ограничение: два тела из одного материала контакта не увидят — в проекте нет
// идентификатора тела (part id из k-файла затирается номером материала,
// input.cpp:230 и 786-807). Смена критерия — правка только этой функции.
bool is_contact_pair(const ParticleSystem& sys, size_t i, size_t j)
{
    return sys.itype[i] > 0 && sys.itype[j] > 0 && sys.itype[i] != sys.itype[j];
}

// Локальный CFL-шаг пары: тот же критерий, что compute_cfl_timestep
// (time_integration.cpp:13-27), суженный на пару. Гарантирует
// dt_applied <= dt_ref, поэтому жёсткость от dt_ref не образует петли
// «жёсткость -> шаг -> жёсткость».
double pair_dt_ref(const SimulationConfig& cfg, const ParticleSystem& sys, size_t k)
{
    double v2 = 0.0;
    for (size_t d = 0; d < dim; ++d) v2 += sys.vx[d][k] * sys.vx[d][k];
    const double denom = sys.c[k] + std::sqrt(v2);
    return denom > 0.0 ? cfg.CFL * sys.hsml[k] / denom : 1e30;
}

} // namespace
```

Тело цикла по парам (внутри существующего цикла [ext_force.cpp:76](../../src/ext_force.cpp#L76),
после ветки `is_real_virtual_pair`):

```cpp
if (contact_type != contact::Off && is_contact_pair(sys, i, j))
{
    std::array<double, dim> dx{};      // x_i - x_j
    double r2 = 0.0;
    for (size_t d = 0; d < dim; ++d) { dx[d] = sys.x[d][i] - sys.x[d][j]; r2 += dx[d]*dx[d]; }
    const double r = std::sqrt(r2);
    const double p = 0.5 * (sys.hsml[i] + sys.hsml[j]) - r;    // (2.1)

    if (r > 0.0 && p > 0.0)
    {
        double dot = 0.0, v2 = 0.0;
        for (size_t d = 0; d < dim; ++d)
        {
            const double dv = sys.vx[d][i] - sys.vx[d][j];
            dot += dx[d] * dv;
            v2  += dv * dv;
        }
        contact::PairGeom g;
        g.r = r;  g.p = p;
        g.pdot  = -dot / r;                 // > 0 при сближении
        g.v_rel = std::sqrt(v2);

        const double Ri = 0.5 * sys.hsml[i], Rj = 0.5 * sys.hsml[j];
        const double Gi = shear_modulus_by_itype(cfg.mats, sys.itype[i]);
        const double Gj = shear_modulus_by_itype(cfg.mats, sys.itype[j]);

        contact::PairProps pr;
        pr.m_red  = sys.mass[i] * sys.mass[j] / (sys.mass[i] + sys.mass[j]);
        pr.c_pair = std::min(sys.c[i], sys.c[j]);
        pr.g_eff  = (Gi > 0.0 && Gj > 0.0) ? Gi * Gj / (Gi + Gj) : 0.0;
        pr.r_eff  = Ri * Rj / (Ri + Rj);
        pr.dt_ref = std::min(pair_dt_ref(cfg, sys, i), pair_dt_ref(cfg, sys, j));

        const contact::Result res =
            contact::penalty_force(contact_type, g, pr, cfg.contact_penalty_scale, dt);

        ++wb.contact_pair_count;
        if (p > wb.contact_max_delta) wb.contact_max_delta = p;
        if (res.dt_limit < wb.min_dt_contact) wb.min_dt_contact = res.dt_limit;

        if (res.fn > 0.0)
        {
            const double inv_r = 1.0 / r;
            for (size_t d = 0; d < dim; ++d)
            {
                const double f = res.fn * dx[d] * inv_r;    // вдоль n = (x_i - x_j)/r
                dvxdt[d][i] += f / sys.mass[i];             // dvxdt — ускорения
                dvxdt[d][j] -= f / sys.mass[j];             // точный баланс импульса
            }
        }
    }
}
```

`contact_type` читается из `cfg` **один раз до цикла** (`const int contact_type = cfg.contact_force_type;`);
`wb.reset_contact()` вызывается в начале `ext_force`, рядом с обнулением `dvxdt`
([ext_force.cpp:55-62](../../src/ext_force.cpp#L55-L62)).

Существующая проверка NaN/Inf в конце функции ([ext_force.cpp:144-165](../../src/ext_force.cpp#L144-L165))
покрывает и контактный вклад — оставить как есть.

### 3.4 Ограничение шага в `time_integration.cpp`

После [строки 144](../../src/time_integration.cpp#L144):

```cpp
dt = compute_cfl_timestep(cfg, sys, ntotal);
if (wb.min_dt_contact < dt) dt = wb.min_dt_contact;   // ADR-007: контактный предел
```

Диагностическую строку («шаг ограничен CFL / КОНТАКТОМ», `contact_pair_count`,
`contact_max_delta`) выводить в существующем блоке печати прогресса
([time_integration.cpp:87-88](../../src/time_integration.cpp#L87-L88)) — только при
`contact_force_type != 0`, чтобы не менять вывод существующих кейсов.

**Проверка шага 3:** сборка зелёная; `ctest` — все тесты проходят; запуск `SPH` с
`contact_force_type = 0` даёт тот же вывод, что до изменений.

---

## Шаг 4. Тесты оркестратора

Дописать в `test_contact.cpp` (цель уже создана на шаге 2):

| № | Проверка | Ожидание |
|---|---|---|
| O1 | Тип 3, две частицы навстречу, один вызов `ext_force` | `dvxdt[0][0]*m_0 + dvxdt[0][1]*m_1 == 0` (баланс импульса, допуск 1e-12) |
| O2 | Типы 1–3, знак силы | Ускорение частицы `i` направлено от `j` (отталкивание) |
| O3 | Фильтры: `itype[i] == itype[j]`; `itype[j] < 0` (виртуальная); `itype[j] == 0` | `contact_pair_count == 0`, вклад в `dvxdt` нулевой |
| O4 | `p <= 0` (частицы разведены) | `contact_pair_count == 0`, `min_dt_contact == 1e30` |
| O5 | `contact_force_type = 0` | `dvxdt` совпадает с расчётом без контакта бит-в-бит; `contact_pair_count == 0` |
| O6 | Ограничитель: типы 1–2 | `wb.min_dt_contact` равен аналитическому `dt_limit` из M4/M5/M8 |
| O7 | Ограничитель: тип 3 | `wb.min_dt_contact == 1e30` |
| O8 | Диагностика | `contact_max_delta` равен фактическому `p` пары |
| O9 | **Регрессия коллапса шага** | Цикл ×64: `ext_force(dt)` → `dt = min(dt_cfl_модельный, wb.min_dt_contact)`. Ожидание: `dt` выходит на плато (не убывает монотонно), ускорения конечны, NaN/Inf нет. Это тест на отсутствие петли «жёсткость ← применяемый шаг» |

Тесты строят `ParticleSystem` вручную (`resize` + присваивание полей), `SimulationConfig`
через `apply_defaults()` + нужные поля, `WorkBuffers` — локально. `niac`/`pair_i`/`pair_j`
задаются прямо, без вызова поиска соседей.

Замечание для O1: при `dim = 1` (текущее значение [param.hpp:8](../../include/param.hpp#L8))
тесты идут в 1D — это нормально, формулы dim-независимы. Ветку с `dim == 2` писать не надо.

**Проверка шага 4:** `ctest --test-dir build --output-on-failure` — все цели зелёные.

---

## Шаг 5. Документация

Параллелен шагу 4.

1. **`docs/04-explanation/decisions/adr-007-contact-force-types.md`** (новый; **007**,
   а не 012 — ADR доведены до 006). По [шаблону](../../docs/04-explanation/decisions/adr-template.md).
   Обязательно зафиксировать:
   - реестр типов и смысл `contact_force_type = 0` (решение D3);
   - $m_{ij}$ вместо $\varrho R^3$ и почему (размерность при `dim ≠ 3`, единственный
     источник истины по массе);
   - **отказ от множителя $p$ в (2.3)** и его обоснование (решение D1);
   - опечатку (2.7)/(2.8) первоисточника;
   - определение $\Delta t_{\mathrm{ref}}$ и аргумент об отсутствии петли с ограничителем;
   - отсутствие ограничителя у типа 3;
   - откат ветки $F_2$ при неизвестном $G$;
   - критерий контактной пары и его ограничение (решение D2);
   - отложенные фазы: средняя нормаль (2.2), свободная поверхность (§2.4), модель Wang,
     трение и тепловыделение, доля гашения $\gamma$ для типа 3, идентификатор тела.
2. **`docs/04-explanation/architecture/contact-particle-particle.md`** (новый): формулы,
   геометрия по линии центров, ограничения (сила только нормальная), диагностика,
   как выбирать $K_p$.
3. **`mkdocs.yml`**: две записи в `nav` — в раздел «Архитектура»
   ([строки 109-115](../../mkdocs.yml#L109-L115)) и в «Решения (ADR)»
   ([строки 116-123](../../mkdocs.yml#L116-L123)).
4. **`docs/03-reference/config/initial-conditions.md`**: две строки в таблицу ключей
   ([строки 21-38](../../docs/03-reference/config/initial-conditions.md#L21-L38)).
5. **`docs/03-reference/config/simulation-config.md`**: новые поля.
6. **`docs/02-how-to/01-run-tests.md`**: строка `contact_test` в таблицу состава тестов.
7. **`docs/_internal/changelog.md`**: запись в раздел «Добавлено».

**Проверка:** `mkdocs build` без ошибок (ссылки проверяются линкером в CI, см. `mkdocs.yml`).

---

## Шаг 6. Интеграционная проверка (выполняет владелец)

Автоматизации не подлежит: требует выбора физического кейса и суждения о результате.

1. Переключить [`param.hpp:8`](../../include/param.hpp#L8) на `dim = 2` или `3`, собрать,
   взять кейс соударения (например, из каталога «Соударение двух стальных цилиндров»
   или BC-серии по керамике).
2. Прогнать при `contact_force_type = 0, 1, 2, 3`, прочих равных.
3. Сравнить: `contact_max_delta` (цель — заметно меньше $\bar h$), долю шагов,
   ограниченных контактом, `vmax`, отсутствие NaN/Inf.
4. Калибровать $K_p$ в диапазоне 0,1–10 (в статье не калиброван).
5. **Отдельно проверить тип 3:** интегратор предиктор-корректор применяет ускорение
   половинами шага ([time_integration.cpp:125,178](../../src/time_integration.cpp#L178)),
   поэтому фактический импульс за шаг не строго равен $m_{ij}\lVert v_{ij}\rVert$;
   формула (2.10) выведена для одношагового приложения. Ожидаемое отклонение —
   неполное гашение относительной скорости за один шаг. Если это окажется существенным,
   опция «доля гашения $\gamma$» перестаёт быть отложенной.
6. Результаты и выбранные настройки — в changelog. Смена значения по умолчанию
   `contact_force_type` — отдельное решение, в этот план не входит.

---

## Порядок, зависимости, коммиты

```
Шаг 0 (сборка)  ──►  Шаг 1 (конфиг+материалы)  ──►  Шаг 2 (математика+её тесты)
                                                          │
                                                          ▼
                                            Шаг 3 (оркестратор+шаг по времени)
                                                          │
                                              ┌───────────┴───────────┐
                                              ▼                       ▼
                                      Шаг 4 (тесты)            Шаг 5 (документация)
                                              └───────────┬───────────┘
                                                          ▼
                                                 Шаг 6 (владелец)
```

Коммиты — по шагам, отдельно: `fix(build): ...`, `feat(contact): конфигурация и реестр G`,
`feat(contact): формулы типов 1-3`, `feat(contact): оркестратор и ограничитель шага`,
`test(contact): ...`, `docs(contact): ADR-007 и описание`.

`contact_force_type = 0` по умолчанию на каждом шаге сохраняет поведение существующих
кейсов бит-в-бит — это инвариант, который должен проверяться (тест O5) до конца работы.

---

## Чек-лист приёмки

- [ ] `cmake --build` без ошибок; `ctest` — все цели зелёные (включая `contact_test`)
- [ ] `contact_force_type = 0`: `dvxdt` бит-в-бит как до внедрения, `contact_pair_count == 0`
- [ ] Типы 1–3 совпадают с аналитикой в тестах M3–M11 (допуск 1e-12)
- [ ] Баланс импульса пары точный для всех типов (O1)
- [ ] Тип 3 не ставит ограничитель шага (O7); типы 1–2 ставят по формулам из §0 (O6)
- [ ] Регрессия коллапса шага пройдена (O9)
- [ ] Неизвестный $G$ не обнуляет контакт типа 1 (M6)
- [ ] Виртуальные частицы и пары одного материала отфильтрованы (O3)
- [ ] ADR-007 создан, включён в `mkdocs.yml`, все шесть решений из `02-decisions.md`
      зафиксированы с обоснованием
- [ ] `mkdocs build` без ошибок
