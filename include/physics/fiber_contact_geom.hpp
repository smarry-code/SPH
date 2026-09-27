#pragma once
// physics/fiber_contact_geom.hpp
// Геометрия одного Fiber-Fiber-Contact primitive -- fiber/CUDA_PLAN.md §3.2,
// этап Г2 (часть "Fiber-Fiber-Contact"): перевод математики
// interactions/fiber_contact.hpp::compute_fiber_contact_geometry (4 случая
// Surface-Surface/Surface-End/End-Surface/End-End + вырожденный почти
// параллельный, веса Meyer "w_ij") в SPH_HD-функцию, общую для CPU- и (в
// будущей части Г7) GPU-путей.
//
// Построчная копия арифметики compute_fiber_contact_geometry (порядок
// операций и типы промежуточных величин не изменены -- требование побитной
// регрессии, план §3.2). Разбор знаков/ориентации (P_ij/P_ji, n_ij, почему
// именно эта формула Hertz для пары РАЗНЫХ материалов) НЕ повторяется здесь
// -- он в шапке include/interactions/fiber_contact.hpp, там же обоснование
// самого правила "w_ij" ниже.
//
// FiberContactGeomInput -- НЕ FiberView/FiberState из плана §3.2: чистый
// вход одной геометрии (без std::vector/topology/materials), собираемый
// вызывающей стороной (CPU: interactions/fiber_contact.hpp из FiberState +
// FiberTopology + FiberMaterial; GPU в будущей части Г7: DeviceFiberState +
// DeviceFiberTopology). x_i_prev/x_i_next (и j-сторона) читаются функцией
// ТОЛЬКО когда case_ отводит этой стороне Surface-роль (SurfaceEnd/
// EndSurface/SurfaceSurface, см. interactions/fiber_contact.hpp::
// build_fiber_contact_search -- Surface-роль назначается ТОЛЬКО частицам с
// prev!=-1 && next!=-1, поэтому вызывающая сторона гарантированно может
// заполнить эти поля валидными координатами до вызова).
//
// FiberContactGeomOutput -- дублирует поля interactions::FiberContactGeometry
// (contact, n, D, F, delta, P_ij, P_ji, w_i/w_j), кроме i_neighbor/j_neighbor
// (rhs::index_t -- индекс, не геометрия): вместо индекса -- флаг
// i_neighbor_is_prev/j_neighbor_is_prev (== "P>0" исходной формулы
// nb_out = P>0 ? prev : next), CPU-обёртка переводит его в rhs::index_t по
// prim.i_prev/i_next, зная prim (physics-слой topology не видит).
#include <cmath>
#include <cstddef>

#include "param.hpp"       // dim
#include "physics/hd.hpp"  // SPH_HD

namespace physics
{

// Случаи геометрии -- те же 4 варианта, что interactions::FiberContactCase
// (та же нумерация 0..3: SurfaceSurface/SurfaceEnd/EndSurface/EndEnd), но
// как отдельный физический тип (physics/* не знает про topology.prev/next
// вызывающей стороны) -- CPU-обёртка конвертирует один в другой через
// static_cast<int>.
enum class FiberContactCase : int
{
    SurfaceSurface = 0,
    SurfaceEnd = 1,
    EndSurface = 2,
    EndEnd = 3,
};

// Порог |p_i x p_j| < eps -- см. interactions/fiber_contact.hpp,
// kFiberContactParallelEps (перенесено сюда без изменения значения/смысла --
// это часть геометрии, которая теперь целиком живёт в physics/*). Известная
// неточность вырожденной ветки при смещении не вдоль нормали (D ~
// sqrt(y0^2+h^2) вместо истинного h) и малый запас порога 1e-3 во float --
// НЕ исправляются этой частью Г2c (fiber/CUDA_PLAN.md, задание координатора).
constexpr real_t kFiberContactParallelEps = real_t(1e-3);

// Защита от деления на ~0 при нормировке n (совпадающие точки/оси).
constexpr real_t kFiberContactMinD = real_t(1e-12);

// ----------------------------------------------------------------------
// Векторные помощники (3-компонентные, эта геометрия определена только при
// dim==3 -- fiber/SPH.md, Р9). Побитная копия локальных лямбд sub/dot/norm/
// cross из interactions/fiber_contact.hpp::compute_fiber_contact_geometry,
// вынесенных в отдельные SPH_HD-функции (безопасный перенос по правилам
// плана §3.2 -- "вынос в inline-функцию с теми же операциями"). Имена с
// префиксом fiber_contact_ -- physics/fiber_solid_geom.hpp объявляет свою
// копию с другим префиксом (сознательное небольшое дублирование, как
// solid_contact_radius/solid_contact_pair_dt_ref в interactions/
// solid_contact.hpp: не выносить в общий заголовок ради трёх строчек).
// ----------------------------------------------------------------------
SPH_HD inline void fiber_contact_vec3_sub(const real_t a[3], const real_t b[3], real_t out[3])
{
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

SPH_HD inline real_t fiber_contact_vec3_dot(const real_t a[3], const real_t b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

SPH_HD inline real_t fiber_contact_vec3_norm(const real_t a[3])
{
    return std::sqrt(fiber_contact_vec3_dot(a, a));
}

SPH_HD inline void fiber_contact_vec3_cross(const real_t a[3], const real_t b[3], real_t out[3])
{
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

// Побитная копия direction_backward: normalize(prev - next).
SPH_HD inline void fiber_contact_direction_backward(const real_t prev[3], const real_t next[3], real_t out[3])
{
    fiber_contact_vec3_sub(prev, next, out);
    const real_t len = fiber_contact_vec3_norm(out);
    if (len > 0)
    {
        out[0] /= len;
        out[1] /= len;
        out[2] /= len;
    }
}

struct FiberContactGeomInput
{
    int case_ = static_cast<int>(FiberContactCase::EndEnd);

    real_t x_i[dim] = {}, x_j[dim] = {};
    // "Назад" по нити (x[prev]/x[next]) соответствующей стороны -- значимы
    // ТОЛЬКО когда case_ отводит этой стороне Surface-роль (см. комментарий
    // у файла).
    real_t x_i_prev[dim] = {}, x_i_next[dim] = {};
    real_t x_j_prev[dim] = {}, x_j_next[dim] = {};

    real_t r_i = 0, r_j = 0;  // 0.5*d материала бусины (fiber::FiberMaterial::d)
    real_t E_i = 0, nu_i = 0; // fiber::FiberMaterial материала нити i
    real_t E_j = 0, nu_j = 0; // fiber::FiberMaterial материала нити j
};

struct FiberContactGeomOutput
{
    bool contact = false;
    real_t n[dim] = {};
    real_t D = 0;
    real_t F = 0;     // сила Герца, >=0
    real_t delta = 0; // проникновение d_eff - D_ij, >0 при контакте

    // P_ij/P_ji -- см. interactions::FiberContactGeometry (те же поля, тот
    // же смысл): заполняются ВСЕГДА, когда геометрия найдена (даже без
    // контакта), остаются 0 для End-роли соответствующей стороны.
    real_t P_ij = 0;
    real_t P_ji = 0;

    real_t w_i = 0;
    bool i_has_neighbor = false;
    bool i_neighbor_is_prev = false; // true -> сосед i_prev, false -> i_next

    real_t w_j = 0;
    bool j_has_neighbor = false;
    bool j_neighbor_is_prev = false; // true -> сосед j_prev, false -> j_next
};

// Побитная копия compute_fiber_contact_geometry (interactions/
// fiber_contact.hpp) -- 4 случая + вырожденный почти параллельный (SPH.md,
// «Fiber-Fiber-Contact»), сила Герца, веса "w_ij". Порядок операций и типы
// промежуточных величин НЕ изменены.
SPH_HD inline void compute_fiber_contact_geometry(const FiberContactGeomInput& in, FiberContactGeomOutput& out)
{
    out = FiberContactGeomOutput{};
    if constexpr (dim != 3)
    {
        return; // недостижимо при реальном использовании (нити только dim==3, Р9)
    }
    else
    {
        const real_t* xi = in.x_i;
        const real_t* xj = in.x_j;
        real_t r[3];
        fiber_contact_vec3_sub(xi, xj, r); // x_ij := x_i - x_j

        real_t P_ij = 0, P_ji = 0;
        bool i_surface = false;
        bool j_surface = false;
        real_t n[3] = {};
        real_t D = 0;

        const FiberContactCase case_v = static_cast<FiberContactCase>(in.case_);

        if (case_v == FiberContactCase::EndEnd)
        {
            real_t diff[3] = {r[0], r[1], r[2]};
            D = fiber_contact_vec3_norm(diff);
            if (D < kFiberContactMinD) return; // совпадающие бусины -- пропуск, защита от NaN
            n[0] = diff[0] / D; n[1] = diff[1] / D; n[2] = diff[2] / D;
        }
        else if (case_v == FiberContactCase::SurfaceEnd)
        {
            real_t p_i[3];
            fiber_contact_direction_backward(in.x_i_prev, in.x_i_next, p_i);
            real_t xj_minus_xi[3];
            fiber_contact_vec3_sub(xj, xi, xj_minus_xi);
            P_ij = fiber_contact_vec3_dot(p_i, xj_minus_xi);
            const real_t Q_i[3] = {xi[0] + P_ij * p_i[0], xi[1] + P_ij * p_i[1], xi[2] + P_ij * p_i[2]};
            real_t diff[3];
            fiber_contact_vec3_sub(Q_i, xj, diff);
            D = fiber_contact_vec3_norm(diff);
            if (D < kFiberContactMinD) return;
            n[0] = diff[0] / D; n[1] = diff[1] / D; n[2] = diff[2] / D;
            i_surface = true;
        }
        else if (case_v == FiberContactCase::EndSurface)
        {
            real_t p_j[3];
            fiber_contact_direction_backward(in.x_j_prev, in.x_j_next, p_j);
            real_t xi_minus_xj[3];
            fiber_contact_vec3_sub(xi, xj, xi_minus_xj);
            P_ji = fiber_contact_vec3_dot(p_j, xi_minus_xj);
            const real_t Q_j[3] = {xj[0] + P_ji * p_j[0], xj[1] + P_ji * p_j[1], xj[2] + P_ji * p_j[2]};
            real_t diff[3];
            fiber_contact_vec3_sub(xi, Q_j, diff);
            D = fiber_contact_vec3_norm(diff);
            if (D < kFiberContactMinD) return;
            n[0] = diff[0] / D; n[1] = diff[1] / D; n[2] = diff[2] / D;
            j_surface = true;
        }
        else // SurfaceSurface
        {
            real_t p_i[3], p_j[3];
            fiber_contact_direction_backward(in.x_i_prev, in.x_i_next, p_i);
            fiber_contact_direction_backward(in.x_j_prev, in.x_j_next, p_j);
            real_t cr[3];
            fiber_contact_vec3_cross(p_i, p_j, cr);
            const real_t cr_norm = fiber_contact_vec3_norm(cr);

            if (cr_norm < kFiberContactParallelEps)
            {
                // Почти параллельные нити -- переход на формулу Surface-End
                // (см. шапку interactions/fiber_contact.hpp, пункт про
                // вырожденную ветку): i остаётся "поверхностной" стороной, j
                // трактуется как точка.
                real_t xj_minus_xi[3];
                fiber_contact_vec3_sub(xj, xi, xj_minus_xi);
                P_ij = fiber_contact_vec3_dot(p_i, xj_minus_xi);
                const real_t Q_i[3] = {xi[0] + P_ij * p_i[0], xi[1] + P_ij * p_i[1], xi[2] + P_ij * p_i[2]};
                real_t diff[3];
                fiber_contact_vec3_sub(Q_i, xj, diff);
                D = fiber_contact_vec3_norm(diff);
                if (D < kFiberContactMinD) return;
                n[0] = diff[0] / D; n[1] = diff[1] / D; n[2] = diff[2] / D;
                i_surface = true;
            }
            else
            {
                // Ближайшие точки двух прямых (2x2 система -- см. шапку
                // interactions/fiber_contact.hpp).
                const real_t c = fiber_contact_vec3_dot(p_i, p_j);
                const real_t denom = real_t(1.0) - c * c; // == cr_norm^2 (>= eps^2 здесь)
                const real_t rpi = fiber_contact_vec3_dot(r, p_i);
                const real_t rpj = fiber_contact_vec3_dot(r, p_j);
                P_ji = (rpj - c * rpi) / denom;
                P_ij = c * P_ji - rpi;
                const real_t Q_i[3] = {xi[0] + P_ij * p_i[0], xi[1] + P_ij * p_i[1], xi[2] + P_ij * p_i[2]};
                const real_t Q_j[3] = {xj[0] + P_ji * p_j[0], xj[1] + P_ji * p_j[1], xj[2] + P_ji * p_j[2]};
                real_t diff[3];
                fiber_contact_vec3_sub(Q_i, Q_j, diff);
                D = fiber_contact_vec3_norm(diff);
                if (D < kFiberContactMinD)
                {
                    // Оси фактически пересекаются (D==0) -- нормаль из
                    // p_i x p_j, ориентированная "от j к i" (dot(n,r)>=0).
                    n[0] = cr[0] / cr_norm; n[1] = cr[1] / cr_norm; n[2] = cr[2] / cr_norm;
                    if (fiber_contact_vec3_dot(n, r) < 0) { n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2]; }
                    D = 0;
                }
                else
                {
                    n[0] = diff[0] / D; n[1] = diff[1] / D; n[2] = diff[2] / D;
                }
                i_surface = true;
                j_surface = true;
            }
        }

        // P_ij/P_ji заполняются ВСЕГДА, когда геометрия дошла до этой точки
        // (даже если ниже окажется, что контакта нет, D>=d_eff).
        out.P_ij = P_ij;
        out.P_ji = P_ji;

        // --- Hertz-сила (пара РАЗНЫХ материалов -- см. шапку interactions/
        // fiber_contact.hpp, пункт 2) ------------------------------------
        const real_t ri = in.r_i;
        const real_t rj = in.r_j;
        const real_t d_eff = ri + rj;
        if (!(D < d_eff)) return; // нет проникновения

        const real_t delta = d_eff - D;
        const real_t inv_Estar = (in.E_i > 0 ? (real_t(1.0) - in.nu_i * in.nu_i) / in.E_i : real_t(0.0)) +
                                  (in.E_j > 0 ? (real_t(1.0) - in.nu_j * in.nu_j) / in.E_j : real_t(0.0));
        if (!(inv_Estar > 0)) return; // некорректный материал (E<=0 с обеих сторон) -- защита
        const real_t E_star = real_t(1.0) / inv_Estar;
        const real_t R_star = (ri > 0 && rj > 0) ? (ri * rj) / (ri + rj) : real_t(0.0);
        const real_t F = (real_t(4.0) / real_t(3.0)) * E_star * std::sqrt(R_star) * std::pow(delta, real_t(1.5));

        out.contact = true;
        out.D = D;
        out.delta = delta;
        out.F = F;
        out.n[0] = n[0];
        out.n[1] = n[1];
        out.n[2] = n[2];

        // --- Веса и соседи (SPH.md, "w_ij") -- side_weight из interactions/
        // fiber_contact.hpp, побитная копия: surface_role==false (End) ->
        // w=1, соседа нет.
        if (!i_surface)
        {
            out.w_i = real_t(1.0);
            out.i_has_neighbor = false;
        }
        else
        {
            real_t d_p_vec[3], d_n_vec[3];
            fiber_contact_vec3_sub(xi, in.x_i_prev, d_p_vec);
            fiber_contact_vec3_sub(xi, in.x_i_next, d_n_vec);
            const real_t d_p = fiber_contact_vec3_norm(d_p_vec);
            const real_t d_n = fiber_contact_vec3_norm(d_n_vec);
            real_t w;
            if (P_ij > 0 && P_ij < d_p) w = (d_p - P_ij) / d_p;
            else if (P_ij <= 0 && P_ij > -d_n) w = (d_n + P_ij) / d_n;
            else w = 0; // контакт проецируется дальше ближайшего соседа
            out.w_i = w;
            out.i_has_neighbor = true;
            out.i_neighbor_is_prev = (P_ij > 0);
        }

        if (!j_surface)
        {
            out.w_j = real_t(1.0);
            out.j_has_neighbor = false;
        }
        else
        {
            real_t d_p_vec[3], d_n_vec[3];
            fiber_contact_vec3_sub(xj, in.x_j_prev, d_p_vec);
            fiber_contact_vec3_sub(xj, in.x_j_next, d_n_vec);
            const real_t d_p = fiber_contact_vec3_norm(d_p_vec);
            const real_t d_n = fiber_contact_vec3_norm(d_n_vec);
            real_t w;
            if (P_ji > 0 && P_ji < d_p) w = (d_p - P_ji) / d_p;
            else if (P_ji <= 0 && P_ji > -d_n) w = (d_n + P_ji) / d_n;
            else w = 0;
            out.w_j = w;
            out.j_has_neighbor = true;
            out.j_neighbor_is_prev = (P_ji > 0);
        }
    }
}

} // namespace physics
