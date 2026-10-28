#pragma once
// physics/solid_finalize.hpp
// Finalize одной частицы Solid -- fiber/CUDA_PLAN.md §3.2, этап Г2д:
// построчный перенос interactions/solid_sph.hpp::finalize_solid (тело
// parallel_for по частицам) в SPH_HD-функцию, общую для CPU- и (в будущей
// части Г4) GPU-путей. Порядок операций и типы промежуточных величин не
// изменены -- требование побитной регрессии (план §3.2).
//
// Осевые члены (axis_sim && dim==2), закон Гука с производной Яуманна
// (ИСПРАВЛЕННЫЙ знак у syz*ryz, см. tests/test_kinematics_3d.cpp) и
// девиаторная работа s:eps/rho в du -- та же математика, что в
// interactions::finalize_solid; сам вызов material_shear_modulus заменён на
// physics::material_shear_modulus (physics/material_table.hpp, этап Г2b) --
// побитная копия того же диспетчера над MaterialTable вместо Materials.
//
// SolidFinalizeInput/Output -- чистый вход/выход ОДНОЙ частицы (без
// std::vector/ParticleSystem): вызывающая сторона (CPU: interactions/
// solid_sph.hpp::finalize_solid из ParticleSystem + reduced[i]; GPU в
// будущей части Г4: DeviceSolidState) собирает вход и раскладывает выход.
#include <cstddef>

#include "param.hpp"                  // dim, axis_sim
#include "physics/hd.hpp"             // SPH_HD
#include "physics/material_table.hpp" // MaterialTable, material_shear_modulus

namespace physics
{

struct SolidFinalizeInput
{
    // "Сырые" накопители одной частицы -- interactions::SolidDelta reduced[i]
    // (сумма по всем primitives, до осевой поправки/домножения на 0.5).
    real_t dv[dim] = {};
    real_t drho_raw = 0;
    real_t du_raw = 0;
    real_t dexx = 0, deyy = 0, dezz = 0, dexy = 0, dexz = 0, deyz = 0;
    real_t drxy = 0, drxz = 0, dryz = 0;

    // Состояние частицы (readonly), нужное finalize.
    real_t x0 = 0;  // sys.x[0][i] -- читается ТОЛЬКО при axis_sim && dim==2
    real_t vx0 = 0; // sys.vx[0][i] -- читается ТОЛЬКО при axis_sim && dim==2
    real_t rho = 0; // sys.rho[i]
    real_t sxx = 0, syy = 0, szz = 0, sxy = 0, syz = 0, sxz = 0; // текущий (дособытийный) тензор напряжений
    int itype = 0;
};

struct SolidFinalizeOutput
{
    real_t drho = 0;
    real_t du = 0;
    real_t dv[dim] = {};
    real_t dsxx = 0, dsyy = 0, dszz = 0, dsxy = 0, dsyz = 0, dsxz = 0;
    real_t exx = 0, eyy = 0, ezz = 0, exy = 0, exz = 0, eyz = 0;
    real_t rxy = 0, rxz = 0, ryz = 0;
};

// Побитная копия тела parallel_for из interactions::finalize_solid. Защита
// от деления на 0 при r==0 (axis_sim && dim==2) остаётся ПОСЛЕДОВАТЕЛЬНОЙ
// CPU-проверкой ДО вызова этой функции (см. там же комментарий про UB-бросок
// из параллельного региона) -- сюда не перенесена.
SPH_HD inline void finalize_solid_particle(const MaterialTable& mtab, const SolidFinalizeInput& in, SolidFinalizeOutput& out)
{
    real_t exx = in.dexx, eyy = in.deyy, ezz = in.dezz;
    real_t exy = in.dexy, exz = in.dexz, eyz = in.deyz;
    const real_t rxy = in.drxy, rxz = in.drxz, ryz = in.dryz;

    real_t drho = in.drho_raw;
    if constexpr (axis_sim && dim == 2)
    {
        // Окружная скорость деформации -- ПОЧАСТИЧНАЯ формула, не парная
        // сумма (int_force.cpp, "Окружная скорость деформации").
        ezz = in.vx0 / in.x0;
        // Осевой член неразрывности (density.cpp, "Осесимметричная поправка
        // к уравнению непрерывности").
        drho -= in.rho * in.vx0 / in.x0;
    }
    out.drho = drho;

    for (std::size_t dd = 0; dd < dim; ++dd) out.dv[dd] = in.dv[dd];

    // Шим совместимости (fiber/PLAN.md, этап Ф4-B1; удалён из single_step на
    // этапе Г2д -- см. интерфейс prochnost_cem_JH2 в physics/strength.hpp) --
    // поле оставлено в Output, т.к. эти же величины нужны JH-2 явными
    // аргументами.
    out.exx = exx; out.eyy = eyy; out.ezz = ezz;
    out.exy = exy; out.exz = exz; out.eyz = eyz;
    out.rxy = rxy; out.rxz = rxz; out.ryz = ryz;

    // Энергия: 0.5*(сырая сумма давление+вязкость) + девиаторная работа.
    const real_t dedt_half = real_t(0.5) * in.du_raw;
    real_t deviatoric;
    if constexpr (axis_sim && dim == 2)
    {
        deviatoric = (in.sxx * exx + in.syy * eyy + in.szz * ezz + real_t(2.0) * in.sxy * exy) / in.rho;
    }
    else
    {
        deviatoric = (in.sxx * exx + in.syy * eyy + in.szz * ezz
                     + real_t(2.0) * in.sxy * exy + real_t(2.0) * in.syz * eyz + real_t(2.0) * in.sxz * exz) / in.rho;
    }
    out.du = dedt_half + deviatoric;

    // Закон Гука + производная Яуманна (int_force.cpp, финальный блок).
    if (in.itype > 2)
    {
        const real_t G_ = material_shear_modulus(mtab, in.itype);
        if constexpr (dim > 1)
        {
            out.dsxx = real_t(2.0) * G_ * (exx + drho / (real_t(3.0) * in.rho)) + real_t(2.0) * (in.sxy * rxy + in.sxz * rxz);
            out.dsyy = real_t(2.0) * G_ * (eyy + drho / (real_t(3.0) * in.rho)) + real_t(2.0) * (in.syz * ryz - in.sxy * rxy);
            // Исправленный знак у syz*ryz (см. tests/test_kinematics_3d.cpp,
            // проверки C2-C4).
            out.dszz = real_t(2.0) * G_ * (ezz + drho / (real_t(3.0) * in.rho)) - real_t(2.0) * (in.sxz * rxz + in.syz * ryz);
            out.dsxy = real_t(2.0) * G_ * exy + rxy * (in.syy - in.sxx) + in.sxz * ryz + in.syz * rxz;
            out.dsyz = real_t(2.0) * G_ * eyz - in.sxy * rxz - in.sxz * rxy + ryz * (in.szz - in.syy);
            out.dsxz = real_t(2.0) * G_ * exz + rxz * (in.szz - in.sxx) - in.sxy * ryz + in.syz * rxy;
        }
        else
        {
            out.dsxx = real_t(2.0) * G_ * (exx + drho / (real_t(3.0) * in.rho));
            out.dsyy = real_t(2.0) * G_ * (drho / (real_t(3.0) * in.rho));
            out.dszz = real_t(2.0) * G_ * (drho / (real_t(3.0) * in.rho));
        }
    }
}

} // namespace physics
