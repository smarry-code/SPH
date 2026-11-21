# ADR-005: Размыкание зависимости EOS (чистые функции)

**Статус:** Принято
**Дата:** Этап 6.3 «Архитектурная чистота»

## Контекст

Функции уравнений состояния (`eos_mi_grunaizen`, `eos_ua_teta`,
`urs_cem_JH2`) ранее читали материальные константы из неконстантных глобалов
`param.hpp` (`aa`, `zz`, `gamm`, `const_JH_2`, `n_Me`). Это создавало
циклическую зависимость `eos.hpp` ↔ `param.hpp` и делало функции
нетестируемыми изолированно.

## Решение

УРС выделены в отдельный модуль [`eos.hpp`](../../../include/eos.hpp) и
реализованы как **чистые функции**, принимающие материальные константы по
ссылке:

```cpp
void eos_mi_grunaizen(const Materials& mats, size_t nmat,
                      double rho, double rho_0, double u,
                      double& p, double& c);
```

Это размыкает зависимость от `param.hpp` и позволяет тестировать УРС с любым
набором констант.

## Последствия

**Положительные:**

- Изолируемое тестирование УРС (см. [`test_jc_eos`](../../../tests/test_jc_eos.cpp));
- Чёткое разделение: УРС — в `eos.hpp` (чистые функции), плотность — в `interactions/solid_sph.hpp`;
- Возможность подстановки моделей состояния без изменения остального конвейера.

**Отрицательные:**

- Увеличение числа параметров функций;
- Необходимость передавать `nmat` (индекс материала) явно.

## Связанные файлы

- [`eos.hpp`](../../../include/eos.hpp): объявления УРС.
- [`eos.cpp`](../../../src/eos.cpp): реализации.
- [`materials.hpp`](../../../include/materials.hpp:130): источник констант.

## Актуальное состояние (2026-11)

- Сигнатуры `eos_mi_grunaizen`, `eos_ua_teta`, `urs_cem_JH2` не менялись
  ([`eos.hpp`](../../../include/eos.hpp)). Добавлена `eos_jc(const JcMaterial&, ...)`,
  принимающая запись материала, а не `Materials` + индекс (учитывает `p_min` и `D`).
- Диспетчер по `itype` существует в двух копиях: `urs()`
  ([`int_force.cpp`](../../../src/int_force.cpp)) и `solid_eos()`
  ([`solid_sph.hpp`](../../../include/interactions/solid_sph.hpp)); их согласованность
  для JC проверяет `test_jc_urs_consistency`.
- Изолированно УРС проверяют [`test_jc_eos`](../../../tests/test_jc_eos.cpp) и
  `test_jc_integration`; `test_particle_system` проверяет только `ParticleSystem`.
