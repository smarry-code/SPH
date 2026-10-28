#pragma once
// physics/eos.hpp
// Уравнения состояния (УРС) -- fiber/CUDA_PLAN.md §3.2, этап Г2д:
// построчный перенос src/eos.cpp (p_gas, p_art_water, eos_mi_grunaizen,
// eos_ua_teta, urs_cem_JH2) в SPH_HD-функции над physics::MaterialTable
// вместо materials.hpp::Materials (std::vector -- запрещён в SPH_HD-коде,
// план §3.2 "Правила кода"). Порядок операций и вызовы pow/sqrt (без
// std::, ровно как в оригинале -- namespace lookup находит их через ADL/
// глобальное <cmath>, тот же результат, что и раньше) не изменены -- этого
// требует побитная регрессия.
//
// solid_eos_dispatch -- SPH_HD-аналог диспетчера
// interactions::solid_sph.hpp::solid_eos (сам он теперь тонкая обёртка над
// этой функцией, см. там же) и src/int_force.cpp::urs() (тоже тонкая
// обёртка). Раньше оба диспетчера бросали исключение на неизвестном itype --
// SPH_HD-функциям исключения запрещены (план §3.2), поэтому здесь --
// SphErrorCode::kUnknownItype (план §3.7: "на GPU -- коды ошибок... исключение
// с тем же текстом бросает хост"); CPU-обёртки при kUnknownItype бросают то
// же самое исключение с тем же текстом, что и раньше.
#include <cmath>
#include <cstddef>

#include "physics/hd.hpp"
#include "physics/material_table.hpp"

namespace physics
{

// src/eos.cpp::p_gas -- адиабатическое УРС газа.
SPH_HD inline void p_gas(real_t rho, real_t u, real_t& p, real_t& c)
{
    const real_t gamma = real_t(1.4);
    p = (gamma - real_t(1.0)) * rho * u;
    if (p >= 0)
    {
        c = std::sqrt(gamma * p / rho);
    }
    else
    {
        c = 100;
    }
}

// src/eos.cpp::p_art_water -- упрощённое УРС воды.
SPH_HD inline void p_art_water(real_t rho, real_t& p, real_t& c)
{
    const real_t B = real_t(1000.0); // Модуль объёмной упругости
    p = B * (std::pow(rho / real_t(1000.0), real_t(7.0)) - real_t(1.0));
    c = std::sqrt(real_t(7.0) * B * std::pow(rho / real_t(1000.0), real_t(6.0)) / real_t(1000.0));
}

// src/eos.cpp::eos_mi_grunaizen (itype 3-4: тантал, алюминий) -- t.aa[k]/
// t.zz[k]/t.gamm[k] с k = itype - MaterialTable::kNMe -- побитная копия
// mats.aa[nmat-n_Me]/mats.zz[nmat-n_Me]/mats.gamm[nmat-n_Me] оригинала.
SPH_HD inline void eos_mi_grunaizen(const MaterialTable& t, int itype, real_t rho, real_t rho_0, real_t u, real_t& p, real_t& c)
{
    const std::size_t k = static_cast<std::size_t>(itype) - MaterialTable::kNMe;
    real_t w = rho / rho_0;
    real_t A1 = t.aa[k] * pow(w, t.zz[k]) * (t.zz[k] - t.gamm[k]) / (t.zz[k] - 1);
    real_t A2 = t.aa[k] * t.gamm[k];
    real_t A3 = w * t.aa[k] * t.zz[k] * (t.gamm[k] - 1) / (t.zz[k] - 1);
    p = A1 - A2 + A3 + u * rho * (t.gamm[k] - 1);
    real_t px = t.aa[k] * (pow(w, t.zz[k]) - 1);
    real_t C = (p - px) / real_t(rho);
    real_t zz1 = t.zz[k] - real_t(1.0);
    real_t aa1 = t.aa[k] * t.zz[k] / rho_0;
    real_t as2 = aa1 * pow(w, zz1) + C * t.gamm[k];
    if (as2 > 0)
    {
        c = sqrt(as2);
    }
    else
    {
        c = real_t(100.0);
    }
}

// src/eos.cpp::eos_ua_teta (itype 5-9: LiF, полиуретаны, PMMA).
SPH_HD inline void eos_ua_teta(const MaterialTable& t, int itype, real_t rho, real_t rho_0, real_t& p, real_t& c)
{
    const std::size_t k = static_cast<std::size_t>(itype) - MaterialTable::kNMe;
    real_t a1 = t.aa[k];
    real_t z1 = t.zz[k];
    real_t az1 = t.aa[k] * t.zz[k] / rho_0;
    real_t z11 = t.zz[k] - real_t(1.0);
    p = a1 * (pow((rho / rho_0), z1) - real_t(1.0));
    real_t as2 = az1 * pow((rho / rho_0), z11);
    if (as2 > 0)
    {
        c = sqrt(as2);
    }
    else
    {
        c = real_t(100.);
    }
}

// src/eos.cpp::urs_cem_JH2 (itype 10-11: B4C, Al2O3) -- t.jh2[mat][col] с
// mat = itype - 10, столбцы 13..16 (K1,K2,K3,T), ровно как
// const_JH_2[mat][13..16] оригинала. D передан по ссылке (как в оригинале)
// и НЕ изменяется -- см. комментарий у interactions::solid_eos.
SPH_HD inline void urs_cem_JH2(const MaterialTable& t, int itype, real_t qo, real_t qo0, real_t dp, real_t& as, real_t& pv, real_t& D)
{
    const std::size_t mat = static_cast<std::size_t>(itype) - 10;
    real_t K1 = t.jh2[mat][13];
    real_t K2 = t.jh2[mat][14];
    real_t K3 = t.jh2[mat][15];
    real_t T = t.jh2[mat][16];
    real_t mu = (qo / qo0) - 1;
    if (D == 0)
    {
        if (mu >= 0)
        {
            pv = K1 * mu + K2 * std::pow(mu, real_t(2)) + K3 * std::pow(mu, real_t(3));
            as = sqrt((K1 / qo0) + 2 * (K2 / qo0) * mu + 3 * (K3 / qo0) * std::pow(mu, real_t(2)));
        }
        else
        {
            pv = K1 * mu;
            if (pv < -T)
            {
                pv = -T;
            }
            as = sqrt(K1 / qo0);
        }
    }
    else
    {
        if (mu >= 0)
        {
            pv = K1 * mu + K2 * std::pow(mu, real_t(2)) + K3 * std::pow(mu, real_t(3)) + dp;
            as = sqrt((K1 / qo0) + 2 * (K2 / qo0) * mu + 3 * (K3 / qo0) * std::pow(mu, real_t(2)));
        }
        else
        {
            pv = K1 * mu;
            if (pv < (-T * (1 - D)))
            {
                pv = -T * (1 - D);
            }
            as = sqrt(K1 / qo0);
        }
    }
}

// Диспетчер по itype -- построчная копия interactions::solid_sph.hpp::
// solid_eos/src/int_force.cpp::urs() над MaterialTable вместо Materials.
// itype вне [1,11] -- SphErrorCode::kUnknownItype, p/c не трогает (как и
// раньше default: не писал в p/c, только бросал).
SPH_HD inline SphErrorCode solid_eos_dispatch(const MaterialTable& t, int itype, real_t rho, real_t rho_0, real_t u,
                                               real_t dp, real_t D, real_t& p, real_t& c)
{
    switch (itype)
    {
        case 1: // идеальный газ
            p_gas(rho, u, p, c);
            return SphErrorCode::kOk;
        case 2: // вода
            p_art_water(rho, p, c);
            return SphErrorCode::kOk;
        case 3: case 4: // тантал, алюминий
            eos_mi_grunaizen(t, itype, rho, rho_0, u, p, c);
            return SphErrorCode::kOk;
        case 5: case 6: case 7: case 8: case 9: // LiF, полиуретаны, PMMA
            eos_ua_teta(t, itype, rho, rho_0, p, c);
            return SphErrorCode::kOk;
        case 10: case 11: // B4C, Al2O3 (JH-2)
        {
            // urs_cem_JH2 принимает D по неконстантной ссылке, но не изменяет
            // её (см. src/eos.cpp) -- локальная копия достаточна.
            real_t D_local = D;
            urs_cem_JH2(t, itype, rho, rho_0, dp, c, p, D_local);
            return SphErrorCode::kOk;
        }
        default:
            return SphErrorCode::kUnknownItype;
    }
}

} // namespace physics
