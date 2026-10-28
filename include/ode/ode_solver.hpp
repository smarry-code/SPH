#pragma once
// ode/ode_solver.hpp
// OdeSolver (fiber/PLAN.md, этап Ф8; §3.9 «ODESolver»; fiber/SPH.md, разделы
// «Описание ODESolver» и «ODE solver»): два явных решателя -- PredictorCorrector
// (рабочий, схема проекта) и Euler (эталон для тестов порядка сходимости).
//
// Стиль -- по образцу include/rhs/cpu_executor.hpp/contributions.hpp:
// свободные шаблонные функции + простые struct/class без виртуальных методов,
// БЕЗ общего полиморфного интерфейса. Выбор решателя -- через enum SolverKind
// и switch/if на стороне вызывающего кода (src/time_integration.cpp), а не
// через runtime-полиморфизм: решателей всего два, они инстанцируются один раз
// перед циклом while (см. time_integration.cpp), и рантайм-диспетчеризация
// виртуальными методами не даёт здесь ничего, кроме лишней косвенности.
//
// Операции State::save/State::advance (fiber/PLAN.md §3.9) реализованы как
// методы PredictorCorrectorSolver::predict_*/correct_* и EulerSolver::step_*
// -- они и есть "save(y_n)"/"advance(y, y_base, a, f)" применительно к
// каждой из двух compute-структур (Solid, Fiber). Правила, специфичные для
// Solid (u >= 0, поправка XSPH к скорости, x += dt*v), живут ВНУТРИ этих
// методов (в *_solid), а не в общей части решателя -- ровно как требует
// fiber/PLAN.md §3.9.
//
// PredictorCorrectorSolver -- буквальный перенос схемы, ранее продублированной
// в обеих перегрузках time_integration (src/time_integration.cpp): арифметика
// НЕ меняется ни на бит (fiber/PLAN.md §3.9: «перенос схемы... без
// изменений»). Решатель хранит производные прошлого шага f^{n-1}
// (SolidRhsBuffers/FiberRhsBuffers -- то, что раньше было локальными
// std::vector внутри time_integration, живущими между итерациями while) и
// сохранённые y^n (u_min_/rho_min_/s*_min_/v_min_/v_min_fiber_).
#include <cstddef>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "physics/ode_update.hpp"

namespace ode
{

enum class SolverKind
{
    PredictorCorrector = 0, // по умолчанию (fiber/PLAN.md §Ф8)
    Euler = 1,
};

// ============================================================================
// Строительные блоки State::advance для PredictorCorrector (перенесены из
// src/time_integration.cpp без изменений арифметики). Общий шаблон по типу
// контейнера: у Solid внутри решателя x/vx -- std::vector<real_t>[dim]
// (ParticleSystem), у Fiber -- то же самое (fiber::FiberState); dv/av могут
// быть отдельным типом (std::vector<std::vector<real_t>>, как у
// SolidRhsBuffers) -- обоим достаточно operator[][].
// ============================================================================

// v_min = vx; vx += (dt_prev/2) * dv. Координаты не трогаются.
template <class N, class VxArr, class DvArr>
void pc_predict_velocity_half_step(N n, real_t dt_prev, VxArr& vx, const DvArr& dv, VxArr& v_min)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            v_min[d][i] = vx[d][i];
            vx[d][i] = physics::ode_velocity_half_step(vx[d][i], dt_prev, dv[d][i]);
        }
    }
}

// vx = base + coeff*dv + av; x += dt*vx (Solid -- с поправкой XSPH).
// use_base == false -> base = текущая vx (первый шаг, fiber/SPH.md «Первый
//   шаг»: v += dt/2*v̇^0 + ṽ, coeff = dt/2);
// use_base == true  -> base = v_base (сохранённая v^n, fiber/SPH.md
//   «Корректор»: v = v^n + dt*v̇^n + ṽ, coeff = dt).
template <class N, class XArr, class VxArr, class DvArr, class AvArr>
void pc_finalize_velocity_position_with_av(N n, real_t dt, real_t coeff, XArr& x, VxArr& vx, const DvArr& dv,
                                            const VxArr& v_base, bool use_base, const AvArr& av)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            const real_t base = use_base ? v_base[d][i] : vx[d][i];
            vx[d][i] = physics::ode_velocity_with_av(base, coeff, dv[d][i], av[d][i]);
            x[d][i] = physics::ode_position(x[d][i], dt, vx[d][i]);
        }
    }
}

// То же самое без поправки XSPH (Fiber -- fiber/SPH.md, «finalize»: у Fiber
// XSPH не определён, только ẋ = v).
template <class N, class XArr, class VxArr, class DvArr>
void pc_finalize_velocity_position(N n, real_t dt, real_t coeff, XArr& x, VxArr& vx, const DvArr& dv,
                                    const VxArr& v_base, bool use_base)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            const real_t base = use_base ? v_base[d][i] : vx[d][i];
            vx[d][i] = physics::ode_velocity(base, coeff, dv[d][i]);
            x[d][i] = physics::ode_position(x[d][i], dt, vx[d][i]);
        }
    }
}

// ============================================================================
// Явный Эйлер (fiber/PLAN.md §3.9, Р4): y^{n+1} = y^n + dt*f, координата --
// по СТАРОЙ скорости (x^{n+1} = x^n + dt*v^n), БЕЗ поправки XSPH (план не
// упоминает XSPH в формуле Эйлера -- она не добавляется). Поэтому x
// обновляется ДО перезаписи vx.
// ============================================================================
template <class N, class XArr, class VxArr, class DvArr>
void euler_step_position_velocity(N n, real_t dt, XArr& x, VxArr& vx, const DvArr& dv)
{
    for (std::size_t d = 0; d < dim; ++d)
        for (N i = 0; i < n; ++i) x[d][i] = physics::ode_position(x[d][i], dt, vx[d][i]); // x^{n+1} = x^n + dt*v^n (старая скорость)
    for (std::size_t d = 0; d < dim; ++d)
        for (N i = 0; i < n; ++i) vx[d][i] = physics::ode_velocity(vx[d][i], dt, dv[d][i]);
}

// Скалярная величина (Solid: u, rho, s-компоненты): y^{n+1} = y^n + dt*dy.
template <class N, class Arr, class DArr>
void euler_step_scalar(N n, real_t dt, Arr& y, const DArr& dy)
{
    for (N i = 0; i < n; ++i) y[i] += dt * dy[i];
}

// ============================================================================
// Буферы производных RHS одного шага (fiber/PLAN.md §3.9: «решатель хранит
// производные прошлого шага f^{n-1}») -- то, что раньше было локальными
// std::vector внутри time_integration (dvx, du, drho, ds*, av, а также
// dx/t/tdsdt/ds, нужные только сигнатуре single_step). Владеются решателем,
// живут между итерациями while (память переиспользуется, как и раньше).
// ============================================================================
struct SolidRhsBuffers
{
    std::vector<std::vector<real_t>> av, dx, dvx;
    std::vector<real_t> du, drho, ds, t, tdsdt;
    std::vector<real_t> dsxx, dsyy, dszz, dsxy, dsyz, dsxz;

    void resize(std::size_t nbuf)
    {
        av.assign(dim, std::vector<real_t>(nbuf, real_t(0.0)));
        dx.assign(dim, std::vector<real_t>(nbuf, real_t(0.0)));
        dvx.assign(dim, std::vector<real_t>(nbuf, real_t(0.0)));
        du.assign(nbuf, real_t(0.0));
        drho.assign(nbuf, real_t(0.0));
        ds.assign(nbuf, real_t(0.0));
        t.assign(nbuf, real_t(0.0));
        tdsdt.assign(nbuf, real_t(0.0));
        dsxx.assign(nbuf, real_t(0.0));
        dsyy.assign(nbuf, real_t(0.0));
        dszz.assign(nbuf, real_t(0.0));
        dsxy.assign(nbuf, real_t(0.0));
        dsyz.assign(nbuf, real_t(0.0));
        dsxz.assign(nbuf, real_t(0.0));
    }
};

struct FiberRhsBuffers
{
    std::vector<real_t> dv[dim];

    void resize(std::size_t nfiber)
    {
        for (std::size_t d = 0; d < dim; ++d) dv[d].assign(nfiber, real_t(0.0));
    }
};

// ============================================================================
// PredictorCorrectorSolver -- fiber/PLAN.md §3.9, «PredictorCorrector в
// варианте проекта»: буквальный перенос схемы time_integration.cpp. Первый
// шаг (first_step() == true) -- без предиктора, u/rho/s/v += (dt/2)*f^0,
// x += dt*v; шаги n>=2 -- сохранить y^n, предиктор на dt_prev/2 по f^{n-1},
// затем корректор y = y^n + dt*f^n, x += dt*v.
//
// Нити интегрируются той же схемой (fiber/SPH.md, «ODESolver»): только x, v,
// без u/rho/s, БЕЗ поправки XSPH (pc_finalize_velocity_position, не _with_av).
// ============================================================================
class PredictorCorrectorSolver
{
public:
    SolidRhsBuffers solid_rhs;
    FiberRhsBuffers fiber_rhs;

    // nbuf -- размер Solid RHS-буферов (>= ntotal, покрывает виртуальные
    // частицы, как раньше: max(cfg.maxn, ntotal)); ntotal/nfiber -- число
    // реальных частиц Solid/Fiber.
    void reset(std::size_t nbuf, std::size_t ntotal, std::size_t nfiber)
    {
        solid_rhs.resize(nbuf);
        fiber_rhs.resize(nfiber);
        u_min_.assign(ntotal, real_t(0.0));
        rho_min_.assign(ntotal, real_t(0.0));
        sxx_min_.assign(ntotal, real_t(0.0));
        syy_min_.assign(ntotal, real_t(0.0));
        szz_min_.assign(ntotal, real_t(0.0));
        sxy_min_.assign(ntotal, real_t(0.0));
        syz_min_.assign(ntotal, real_t(0.0));
        sxz_min_.assign(ntotal, real_t(0.0));
        for (std::size_t d = 0; d < dim; ++d) v_min_[d].assign(ntotal, real_t(0.0));
        for (std::size_t d = 0; d < dim; ++d) v_min_fiber_[d].assign(nfiber, real_t(0.0));
        first_step_ = true;
    }

    bool first_step() const { return first_step_; }

    // Предиктор Solid на полшага dt (по f^{n-1}, уже лежащей в solid_rhs с
    // прошлого вызова). Нет-оп на первом шаге -- fiber/SPH.md, «Первый шаг»:
    // предиктора нет вовсе.
    void predict_solid(ParticleSystem& sys, std::size_t ntotal, real_t dt)
    {
        if (first_step_) return;
        // half_dt вычисляется ОДИН раз (не на каждое поле, как раньше
        // текстуально) -- та же операция dt/real_t(2.0) даёт тот же результат
        // при каждом вычислении (детерминизм IEEE754), поэтому вынос за
        // пределы полей одной частицы не меняет ни одного бита.
        const real_t half_dt = dt / real_t(2.0);
        for (std::size_t i = 0; i < ntotal; ++i)
        {
            u_min_[i] = sys.u[i];
            const real_t temp_u = physics::ode_solid_energy_axis_term(sys.p[i], sys.vx[0][i], sys.x[0][i], sys.rho[i]);
            sys.u[i] = physics::ode_solid_update_nonneg(sys.u[i], half_dt, solid_rhs.du[i], temp_u);

            if (!summation_density)
            {
                rho_min_[i] = sys.rho[i];
                const real_t temp_rho = physics::ode_solid_density_axis_term(sys.rho[i], sys.vx[0][i], sys.x[0][i]);
                sys.rho[i] = physics::ode_solid_update(sys.rho[i], half_dt, solid_rhs.drho[i], temp_rho);
            }

            sxx_min_[i] = sys.sxx[i];
            sxy_min_[i] = sys.sxy[i];
            syy_min_[i] = sys.syy[i];
            syz_min_[i] = sys.syz[i];
            sxz_min_[i] = sys.sxz[i];
            szz_min_[i] = sys.szz[i];
            sys.sxx[i] = physics::ode_solid_advance(sys.sxx[i], half_dt, solid_rhs.dsxx[i]);
            sys.syy[i] = physics::ode_solid_advance(sys.syy[i], half_dt, solid_rhs.dsyy[i]);
            sys.szz[i] = physics::ode_solid_advance(sys.szz[i], half_dt, solid_rhs.dszz[i]);
            sys.sxy[i] = physics::ode_solid_advance(sys.sxy[i], half_dt, solid_rhs.dsxy[i]);
            sys.syz[i] = physics::ode_solid_advance(sys.syz[i], half_dt, solid_rhs.dsyz[i]);
            sys.sxz[i] = physics::ode_solid_advance(sys.sxz[i], half_dt, solid_rhs.dsxz[i]);
        }
        pc_predict_velocity_half_step(ntotal, dt, sys.vx, solid_rhs.dvx, v_min_);
    }

    // Предиктор Fiber -- только скорость (нет u/rho/s у Fiber).
    void predict_fiber(fiber::FiberState& fb, std::size_t nfiber, real_t dt)
    {
        if (first_step_) return;
        pc_predict_velocity_half_step(nfiber, dt, fb.vx, fiber_rhs.dv, v_min_fiber_);
    }

    void correct_solid(ParticleSystem& sys, std::size_t ntotal, real_t dt)
    {
        const bool fs = first_step_;
        const real_t half_dt = dt / real_t(2.0);
        const real_t coeff = fs ? half_dt : dt; // тот же выбор, что раньше выбирал if(fs)/else почленно
        for (std::size_t i = 0; i < ntotal; ++i)
        {
            const real_t temp_u = physics::ode_solid_energy_axis_term(sys.p[i], sys.vx[0][i], sys.x[0][i], sys.rho[i]);
            sys.u[i] = physics::ode_solid_update_nonneg(fs ? sys.u[i] : u_min_[i], coeff, solid_rhs.du[i], temp_u);

            if (!summation_density)
            {
                const real_t temp_rho = physics::ode_solid_density_axis_term(sys.rho[i], sys.vx[0][i], sys.x[0][i]);
                sys.rho[i] = physics::ode_solid_update(fs ? sys.rho[i] : rho_min_[i], coeff, solid_rhs.drho[i], temp_rho);
            }

            sys.sxx[i] = physics::ode_solid_advance(fs ? sys.sxx[i] : sxx_min_[i], coeff, solid_rhs.dsxx[i]);
            sys.syy[i] = physics::ode_solid_advance(fs ? sys.syy[i] : syy_min_[i], coeff, solid_rhs.dsyy[i]);
            sys.szz[i] = physics::ode_solid_advance(fs ? sys.szz[i] : szz_min_[i], coeff, solid_rhs.dszz[i]);
            sys.sxy[i] = physics::ode_solid_advance(fs ? sys.sxy[i] : sxy_min_[i], coeff, solid_rhs.dsxy[i]);
            sys.syz[i] = physics::ode_solid_advance(fs ? sys.syz[i] : syz_min_[i], coeff, solid_rhs.dsyz[i]);
            sys.sxz[i] = physics::ode_solid_advance(fs ? sys.sxz[i] : sxz_min_[i], coeff, solid_rhs.dsxz[i]);
        }
        // Скорость/координата -- через общий шаблон (axis-ветвление старого
        // кода арифметически совпадало в обеих ветках if/else, поэтому
        // безопасно заменено на безусловный цикл ещё в перегрузке
        // time_integration(fiber::Model&, ...), этап Ф3; здесь -- то же
        // самое, проверено tests/test_fiber_time_integration_equivalence.cpp).
        pc_finalize_velocity_position_with_av(ntotal, dt, coeff, sys.x, sys.vx, solid_rhs.dvx, v_min_,
                                               /*use_base=*/!fs, solid_rhs.av);
    }

    void correct_fiber(fiber::FiberState& fb, std::size_t nfiber, real_t dt)
    {
        const bool fs = first_step_;
        pc_finalize_velocity_position(nfiber, dt, fs ? dt / real_t(2.0) : dt, fb.x, fb.vx, fiber_rhs.dv, v_min_fiber_,
                                       /*use_base=*/!fs);
    }

    // Вызывается РОВНО один раз в конце шага (после correct_solid/correct_fiber).
    void mark_step_done() { first_step_ = false; }

    // ------------------------------------------------------------------
    // Единообразный интерфейс OdeSolver::step (fiber/PLAN.md §3.9), вызываемый
    // одинаково для ЛЮБОГО решателя из общего цикла (src/time_integration.cpp,
    // detail::run_time_integration<Solver>) -- без if/switch по типу решателя
    // в самом цикле: predict()/advance() у EulerSolver имеют ту же сигнатуру
    // (predict() -- нет-оп, у Эйлера нет предикторной фазы).
    // ------------------------------------------------------------------
    void predict(ParticleSystem& sys, std::size_t ntotal, fiber::FiberState& fb, std::size_t nfiber, real_t dt)
    {
        predict_solid(sys, ntotal, dt);
        predict_fiber(fb, nfiber, dt);
    }

    void advance(ParticleSystem& sys, std::size_t ntotal, fiber::FiberState& fb, std::size_t nfiber, real_t dt)
    {
        correct_solid(sys, ntotal, dt);
        correct_fiber(fb, nfiber, dt);
        mark_step_done();
    }

private:
    bool first_step_ = true;
    std::vector<real_t> u_min_, rho_min_;
    std::vector<real_t> sxx_min_, syy_min_, szz_min_, sxy_min_, syz_min_, sxz_min_;
    std::vector<real_t> v_min_[dim];
    std::vector<real_t> v_min_fiber_[dim];
};

// ============================================================================
// EulerSolver -- fiber/PLAN.md §3.9, «Явный Эйлер»: f = RHS(y^n) (одно
// вычисление), y^{n+1} = y^n + dt*f для ВСЕХ переменных (u, rho, v, s, x у
// Solid; x, v у Fiber), x^{n+1} = x^n + dt*v^n (СТАРАЯ скорость, БЕЗ поправки
// XSPH -- план не упоминает XSPH у Эйлера, Р4). Ограничение u >= 0 -- то же
// специфичное для Solid правило advance, что и у PredictorCorrector.
//
// Никакого сохранённого y^n не нужно (одностадийный метод), поэтому
// EulerSolver не хранит v_min/u_min/... -- только буферы RHS (f^n).
// ============================================================================
class EulerSolver
{
public:
    SolidRhsBuffers solid_rhs;
    FiberRhsBuffers fiber_rhs;

    // Сигнатура reset() намеренно совпадает с PredictorCorrectorSolver::reset
    // (nbuf, ntotal, nfiber) -- обе используются одним и тем же шаблонным
    // кодом time_integration.cpp (run_time_integration_impl<Solver>);
    // ntotal Эйлеру не нужен (нет сохранённого y^n), но параметр сохранён
    // для единообразия вызова.
    void reset(std::size_t nbuf, std::size_t /*ntotal*/, std::size_t nfiber)
    {
        solid_rhs.resize(nbuf);
        fiber_rhs.resize(nfiber);
    }

    void step_solid(ParticleSystem& sys, std::size_t ntotal, real_t dt)
    {
        for (std::size_t i = 0; i < ntotal; ++i)
        {
            const real_t temp_u = physics::ode_solid_energy_axis_term(sys.p[i], sys.vx[0][i], sys.x[0][i], sys.rho[i]);
            sys.u[i] = physics::ode_solid_update_nonneg(sys.u[i], dt, solid_rhs.du[i], temp_u);

            if (!summation_density)
            {
                const real_t temp_rho = physics::ode_solid_density_axis_term(sys.rho[i], sys.vx[0][i], sys.x[0][i]);
                sys.rho[i] = physics::ode_solid_update(sys.rho[i], dt, solid_rhs.drho[i], temp_rho);
            }

            sys.sxx[i] = physics::ode_solid_advance(sys.sxx[i], dt, solid_rhs.dsxx[i]);
            sys.syy[i] = physics::ode_solid_advance(sys.syy[i], dt, solid_rhs.dsyy[i]);
            sys.szz[i] = physics::ode_solid_advance(sys.szz[i], dt, solid_rhs.dszz[i]);
            sys.sxy[i] = physics::ode_solid_advance(sys.sxy[i], dt, solid_rhs.dsxy[i]);
            sys.syz[i] = physics::ode_solid_advance(sys.syz[i], dt, solid_rhs.dsyz[i]);
            sys.sxz[i] = physics::ode_solid_advance(sys.sxz[i], dt, solid_rhs.dsxz[i]);
        }
        euler_step_position_velocity(ntotal, dt, sys.x, sys.vx, solid_rhs.dvx);
    }

    void step_fiber(fiber::FiberState& fb, std::size_t nfiber, real_t dt)
    {
        euler_step_position_velocity(nfiber, dt, fb.x, fb.vx, fiber_rhs.dv);
    }

    // ------------------------------------------------------------------
    // Единообразный интерфейс OdeSolver::step (см. комментарий у
    // PredictorCorrectorSolver::predict/advance выше). У Эйлера предикторной
    // фазы нет -- RHS вычисляется прямо в y^n, поэтому predict() -- нет-оп.
    // ------------------------------------------------------------------
    void predict(ParticleSystem& /*sys*/, std::size_t /*ntotal*/, fiber::FiberState& /*fb*/, std::size_t /*nfiber*/,
                 real_t /*dt*/)
    {
        // нет-оп: явный Эйлер вычисляет f = RHS(y^n) без промежуточного
        // предиктора (fiber/PLAN.md §3.9, «Явный Эйлер»).
    }

    void advance(ParticleSystem& sys, std::size_t ntotal, fiber::FiberState& fb, std::size_t nfiber, real_t dt)
    {
        step_solid(sys, ntotal, dt);
        step_fiber(fb, nfiber, dt);
    }
};

} // namespace ode
