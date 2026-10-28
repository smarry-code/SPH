// int_force.cpp
// Определение urs() — диспетчера уравнения состояния (УРС). Объявление — в
// int_force.hpp (там же — почему файл сохранён после удаления int_force()
// на этапе Ф4-B2).
//
// fiber/CUDA_PLAN.md §3.2, этап Г2д: сам диспетчер (switch по itype) —
// теперь тонкая обёртка над physics::solid_eos_dispatch (physics/eos.hpp,
// SPH_HD, используется и interactions::solid_sph.hpp::solid_eos — тот же
// диспетчер, тот же результат бит в бит для тех же itype). Неизвестный
// itype — SphErrorCode::kUnknownItype вместо default: в switch; здесь, в
// CPU-обёртке, по коду ошибки бросается ТО ЖЕ исключение с тем же текстом,
// что и раньше (план §3.7).
#include "int_force.hpp"
#include "physics/eos.hpp"
#include "physics/material_table.hpp"
#include <stdexcept>
#include <string>

void urs(const Materials& mats, ParticleSystem& sys, size_t i)
{
    const physics::MaterialTable t = physics::build_material_table(mats);
    const SphErrorCode err = physics::solid_eos_dispatch(
        t, sys.itype[i], sys.rho[i], sys.rho_0[i], sys.u[i], sys.dp[i], sys.D[i], sys.p[i], sys.c[i]);
    if (err == SphErrorCode::kUnknownItype)
    {
        throw std::runtime_error("urs(): неизвестный itype=" + std::to_string(sys.itype[i]));
    }
}
