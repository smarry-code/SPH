// int_force.cpp
// Определение urs() — диспетчера уравнения состояния (УРС). Объявление — в
// int_force.hpp (там же — почему файл сохранён после удаления int_force()
// на этапе Ф4-B2).
#include "int_force.hpp"
#include <stdexcept>
#include <string>

// Диспетчеризация уравнения состояния по типу материала (itype).
// цепочка if/else без else заменена на switch с default —
// управляемый отказ. Ранее при неизвестном itype УРС молчаливо не пересчитывался
// (p/c оставались прежними => некорректные силы, возможен NaN). Теперь — явное
// исключение с диагностикой. Материалы (mats) передаются по ссылке (ADR-4).
void urs(const Materials& mats, ParticleSystem& sys, size_t i)
{
    switch (sys.itype[i])
    {
        case 1: // идеальный газ
            p_gas(sys.rho[i], sys.u[i], sys.p[i], sys.c[i]);
            break;
        case 2: // вода
            p_art_water(sys.rho[i], sys.p[i], sys.c[i]);
            break;
        case 3: case 4: // тантал, алюминий
            eos_mi_grunaizen(mats, static_cast<size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.u[i], sys.p[i], sys.c[i]);
            break;
        case 5: case 6: case 7: case 8: case 9: // LiF, полиуретаны, PMMA
            eos_ua_teta(mats, static_cast<size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.p[i], sys.c[i]);
            break;
        case 10: case 11: // B4C, Al2O3
            urs_cem_JH2(mats, static_cast<size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.dp[i], sys.c[i], sys.p[i], sys.D[i]);
            break;
        default:
            throw std::runtime_error("urs(): неизвестный itype=" + std::to_string(sys.itype[i]));
    }
}
