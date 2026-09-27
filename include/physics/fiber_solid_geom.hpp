#pragma once
// physics/fiber_solid_geom.hpp
// Геометрия одного Fiber-Solid-Contact primitive -- fiber/CUDA_PLAN.md §3.2,
// этап Г2 (часть "Fiber-Solid-Contact"): перевод математики
// interactions/fiber_solid_contact.hpp::compute_fiber_solid_contact_geometry
// (λ на сегменте/точечный контакт, D_fs, n_fs, δ_fs, сила Герца K_fs) в
// SPH_HD-функцию, общую для CPU- и (в будущей части Г8) GPU-путей.
//
// Построчная копия арифметики (порядок операций и типы промежуточных
// величин не изменены -- требование побитной регрессии, план §3.2). Разбор
// решений (правило "сегмент принадлежит первой бусине", полуинтервал
// λ∈[0,1), K_fs по материалу бусины a) -- в шапке include/interactions/
// fiber_solid_contact.hpp, здесь не повторяется.
//
// FiberSolidContactGeomInput -- чистый вход одной геометрии (без
// ParticleSystem/FiberState/Materials): r_s -- уже посчитанный вызывающей
// стороной physics::solid_contact_particle_radius (physics/
// solid_contact_pair.hpp, SPH_HD, часть Г2 "Solid-Contact"); E_s/nu_s/
// elastic_s_known -- уже посчитанный physics::elastic_table_lookup
// (physics/material_table.hpp, ElasticTable, ЭТА часть Г2) -- материал
// Solid-частицы s "непроизвольный": по решению 4 в шапке interactions/
// fiber_solid_contact.hpp K_fs использует именно материал бусины a, а
// материал s нужен только через (E_s, nu_s).
#include <cmath>
#include <cstddef>

#include "param.hpp"       // dim
#include "physics/hd.hpp"  // SPH_HD

namespace physics
{

// Защита от деления на ~0 при нормировке n (совпадающие точки).
constexpr real_t kFiberSolidContactMinD = real_t(1e-12);

// ----------------------------------------------------------------------
// Векторные помощники (3-компонентные, dim==3 -- fiber/SPH.md, Р9).
// Побитная копия локальных лямбд sub/dot/norm из interactions/
// fiber_solid_contact.hpp::compute_fiber_solid_contact_geometry. Имена с
// префиксом fiber_solid_ -- НЕ переиспользуют одноимённые функции
// physics/fiber_contact_geom.hpp (fiber_contact_vec3_*): небольшое
// сознательное дублирование трёх строчек ради независимости заголовков
// (тот же приём, что solid_contact_radius/solid_contact_pair_dt_ref в
// interactions/solid_contact.hpp) -- одинаковые имена в разных заголовках,
// оба включаемые в один TU (tests/cuda/physics_compile_probe.cu), дали бы
// переопределение.
// ----------------------------------------------------------------------
SPH_HD inline void fiber_solid_vec3_sub(const real_t a[3], const real_t b[3], real_t out[3])
{
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

SPH_HD inline real_t fiber_solid_vec3_dot(const real_t a[3], const real_t b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

SPH_HD inline real_t fiber_solid_vec3_norm(const real_t a[3])
{
    return std::sqrt(fiber_solid_vec3_dot(a, a));
}

struct FiberSolidContactGeomInput
{
    bool has_b = false;
    real_t x_a[dim] = {};
    real_t x_b[dim] = {}; // значим только при has_b
    real_t x_s[dim] = {};
    bool a_is_first = false; // значим всегда, даже при !has_b
    bool b_is_last = false;  // значим только при has_b

    real_t r_f = 0; // 0.5*d материала бусины a
    real_t r_s = 0; // physics::solid_contact_particle_radius(sys.mass[s],...)

    real_t E_a = 0, nu_a = 0; // fiber::FiberMaterial материала бусины a

    // Материал Solid-частицы s -- elastic_s_known отражает elastic_by_itype
    // (materials.hpp): false <=> материал s не определяет (E, nu) (газ/вода
    // itype 1-2 и т.п.) -- physics::elastic_table_lookup, ЭТА часть Г2.
    bool elastic_s_known = false;
    real_t E_s = 0, nu_s = 0;
};

struct FiberSolidContactGeomOutput
{
    bool contact = false;
    real_t n[dim] = {};
    real_t D = 0;
    real_t delta = 0;
    real_t F = 0;      // >= 0
    real_t lambda = 0; // значим только при has_b (входном)
};

// Побитная копия compute_fiber_solid_contact_geometry (interactions/
// fiber_solid_contact.hpp): λ на сегменте (или точечный контакт), D_fs,
// n_fs, δ_fs, сила Герца K_fs*delta^1.5. Порядок операций и типы
// промежуточных величин НЕ изменены.
SPH_HD inline void compute_fiber_solid_contact_geometry(const FiberSolidContactGeomInput& in,
                                                          FiberSolidContactGeomOutput& out)
{
    out = FiberSolidContactGeomOutput{};
    if constexpr (dim != 3)
    {
        return; // недостижимо при реальном использовании (нити только dim==3, Р9)
    }
    else
    {
        const real_t* xs = in.x_s;
        real_t q[3];
        real_t lambda = 0;

        if (!in.has_b)
        {
            q[0] = in.x_a[0]; q[1] = in.x_a[1]; q[2] = in.x_a[2];
        }
        else
        {
            const real_t* xa = in.x_a;
            const real_t* xb = in.x_b;
            real_t e[3];
            fiber_solid_vec3_sub(xb, xa, e);
            const real_t e2 = fiber_solid_vec3_dot(e, e);
            if (!(e2 > 0)) return; // вырожденный сегмент (совпадающие бусины) -- защита от деления на 0

            real_t xs_minus_xa[3];
            fiber_solid_vec3_sub(xs, xa, xs_minus_xa);
            // Сырой (до clamp) параметр t -- решение 3 в шапке interactions/
            // fiber_solid_contact.hpp: клэмп в границу разрешён ТОЛЬКО
            // когда дальше по нити клэмпить некому.
            const real_t t = fiber_solid_vec3_dot(xs_minus_xa, e) / e2;
            if (t >= real_t(1.0))
            {
                if (!in.b_is_last) return;
                lambda = real_t(1.0);
            }
            else if (t < 0)
            {
                if (!in.a_is_first) return;
                lambda = 0;
            }
            else
            {
                lambda = t;
            }
            q[0] = xa[0] + lambda * e[0];
            q[1] = xa[1] + lambda * e[1];
            q[2] = xa[2] + lambda * e[2];
        }

        real_t diff[3];
        fiber_solid_vec3_sub(xs, q, diff);
        const real_t D = fiber_solid_vec3_norm(diff);
        if (D < kFiberSolidContactMinD) return; // совпадающие точки -- пропуск, защита от NaN
        const real_t n[3] = {diff[0] / D, diff[1] / D, diff[2] / D};

        const real_t r_f = in.r_f;
        const real_t r_s = in.r_s;
        const real_t delta = r_f + r_s - D;
        if (!(delta > 0)) return; // нет проникновения

        if (!in.elastic_s_known) return; // защита (search уже отфильтровал такие s)

        const real_t inv_Estar = (in.E_a > 0 ? (real_t(1.0) - in.nu_a * in.nu_a) / in.E_a : real_t(0.0)) +
                                  (in.E_s > 0 ? (real_t(1.0) - in.nu_s * in.nu_s) / in.E_s : real_t(0.0));
        if (!(inv_Estar > 0)) return;
        const real_t E_star = real_t(1.0) / inv_Estar;
        const real_t R_star = (r_f > 0 && r_s > 0) ? (r_f * r_s) / (r_f + r_s) : real_t(0.0);
        const real_t K_fs = (real_t(4.0) / real_t(3.0)) * E_star * std::sqrt(R_star);
        const real_t F = K_fs * std::pow(delta, real_t(1.5));

        out.contact = true;
        out.n[0] = n[0]; out.n[1] = n[1]; out.n[2] = n[2];
        out.D = D;
        out.delta = delta;
        out.F = F;
        out.lambda = lambda;
    }
}

} // namespace physics
