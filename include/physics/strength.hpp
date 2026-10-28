#pragma once
// physics/strength.hpp
// Прочность -- Мизес -- fiber/CUDA_PLAN.md §3.2, этап Г2д: построчный
// перенос src/ppochnost.cpp (sigma_i, prochnost_mises) в SPH_HD-функции над
// physics::MaterialTable вместо materials.hpp::Materials (std::vector --
// запрещён в SPH_HD-коде, план §3.2 "Правила кода"). Порядок операций не
// изменён -- требование побитной регрессии.
//
// JH-2 (prochnost_cem_JH2) -- в этом же файле по списку §3.2 ("strength.hpp
// -- Мизес, JH-2 с явными аргументами"), добавляется отдельным шагом Г2д
// (§4 Г2: "JH-2 на явных параметрах" -- перевод скоростей деформаций/
// поворота на явные аргументы делается ОДНИМ коммитом с удалением шима
// записи sys.exx..ryz в single_step, поэтому JH-2 переносится сюда следующим
// шагом, не вместе с Мизесом).
#include <cmath>
#include <cstddef>

#include "param.hpp"                  // dim, axis_sim
#include "physics/hd.hpp"             // SPH_HD
#include "physics/material_table.hpp" // MaterialTable

namespace physics
{

// src/ppochnost.cpp::sigma_i -- интенсивность девиатора напряжений.
// Параметры оригинала были by-ref (real_t&), хотя тело их не изменяет --
// здесь по значению (тот же читаемый результат, без изменения арифметики).
// Финальный return real_t(0.0) добавлен для корректной компиляции device-
// кода (SPH_HD-функция обязана возвращать значение на всех путях) -- ветка
// недостижима, т.к. dim всегда в {1,2,3} (param.hpp), поэтому на реальный
// результат не влияет.
SPH_HD inline real_t sigma_i(real_t sxx, real_t syy, real_t sxy, real_t syz, real_t sxz)
{
    if (dim == 1) // одномерная деформация
    {
        return real_t(0.75) * sxx * sxx;
    }
    else if (dim == 2)
    {
        if (axis_sim == true)
        {
            return (sxx * sxx + sxx * syy + syy * syy + sxy * sxy); // осесимметричный случай
        }
        else
        {
            return (sxx * sxx + sxy * sxy); // плоский случай
        }
    }
    else if (dim == 3) // трехмерный случай
    {
        return (sxx * sxx + sxx * syy + syy * syy + sxy * sxy + sxz * sxz + syz * syz);
    }
    return real_t(0.0); // недостижимо
}

// src/ppochnost.cpp::prochnost_mises -- побитная копия над MaterialTable
// вместо Materials; itype -- сама частица (nmat оригинала). gem/ge2m/flag_e
// вычисляются, но не используются -- как и в оригинале (сохранено ради
// побитного соответствия и минимального дифа, не удалено по своей воле).
SPH_HD inline void prochnost_mises(const MaterialTable& t, int itype, real_t& sxx, real_t& syy, real_t& szz,
                                    real_t& sxy, real_t& syz, real_t& sxz)
{
    const std::size_t k = static_cast<std::size_t>(itype) - MaterialTable::kNMe;
    real_t gem = t.G[k];
    real_t ge2m = real_t(2.0) * t.G[k];
    real_t smiz2m = std::pow(t.sigy[k], real_t(2)) / real_t(3.0);
    real_t smizm = sqrt(std::pow(t.sigy[k], real_t(2)) / real_t(3.0));
    real_t ripl = sigma_i(sxx, syy, sxy, syz, sxz);
    bool flag_e = true;
    if (ripl > smiz2m)
    {
        real_t corr = smizm / sqrt(ripl);
        sxx *= corr;
        syy *= corr;
        szz *= corr;
        sxy *= corr;
        syz *= corr;
        sxz *= corr;
        flag_e = false;
    }
    (void)gem; (void)ge2m; (void)flag_e; // не используются -- как и в оригинале
}

} // namespace physics
