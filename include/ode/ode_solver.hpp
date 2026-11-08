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
// контейнера: у Solid внутри решателя x/vx -- std::vector<double>[dim]
// (ParticleSystem), у Fiber -- то же самое (fiber::FiberState); dv/av могут
// быть отдельным типом (std::vector<std::vector<double>>, как у
// SolidRhsBuffers) -- обоим достаточно operator[][].
// ============================================================================

// v_min = vx; vx += (dt_prev/2) * dv. Координаты не трогаются.
template <class N, class VxArr, class DvArr>
void pc_predict_velocity_half_step(N n, double dt_prev, VxArr& vx, const DvArr& dv, VxArr& v_min)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            v_min[d][i] = vx[d][i];
            vx[d][i] += (dt_prev / 2.0) * dv[d][i];
        }
    }
}

// vx = base + coeff*dv + av; x += dt*vx (Solid -- с поправкой XSPH).
// use_base == false -> base = текущая vx (первый шаг, fiber/SPH.md «Первый
//   шаг»: v += dt/2*v̇^0 + ṽ, coeff = dt/2);
// use_base == true  -> base = v_base (сохранённая v^n, fiber/SPH.md
//   «Корректор»: v = v^n + dt*v̇^n + ṽ, coeff = dt).
template <class N, class XArr, class VxArr, class DvArr, class AvArr>
void pc_finalize_velocity_position_with_av(N n, double dt, double coeff, XArr& x, VxArr& vx, const DvArr& dv,
                                            const VxArr& v_base, bool use_base, const AvArr& av)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            const double base = use_base ? v_base[d][i] : vx[d][i];
            vx[d][i] = base + coeff * dv[d][i] + av[d][i];
            x[d][i] += dt * vx[d][i];
        }
    }
}

// Режим xsph_mode = position (классический XSPH Монаганa): скорость БЕЗ av,
//   vx = base + coeff*dv;  x += dt*(vx + av)
// -- поправка av только переносит позиции, кинетическая энергия не меняется.
template <class N, class XArr, class VxArr, class DvArr, class AvArr>
void pc_finalize_velocity_position_xsph_position(N n, double dt, double coeff, XArr& x, VxArr& vx, const DvArr& dv,
                                                  const VxArr& v_base, bool use_base, const AvArr& av)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            const double base = use_base ? v_base[d][i] : vx[d][i];
            vx[d][i] = base + coeff * dv[d][i];
            x[d][i] += dt * (vx[d][i] + av[d][i]);
        }
    }
}

// То же самое без поправки XSPH (Fiber -- fiber/SPH.md, «finalize»: у Fiber
// XSPH не определён, только ẋ = v).
template <class N, class XArr, class VxArr, class DvArr>
void pc_finalize_velocity_position(N n, double dt, double coeff, XArr& x, VxArr& vx, const DvArr& dv,
                                    const VxArr& v_base, bool use_base)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (N i = 0; i < n; ++i)
        {
            const double base = use_base ? v_base[d][i] : vx[d][i];
            vx[d][i] = base + coeff * dv[d][i];
            x[d][i] += dt * vx[d][i];
        }
    }
}

// ============================================================================
// Явный Эйлер (fiber/PLAN.md §3.9, Р4): y^{n+1} = y^n + dt*f, координата --
// по СТАРОЙ скорости (x^{n+1} = x^n + dt*v^n), БЕЗ поправки XSPH (план не
// упоминает XSPH в формуле Эйлера -- она не добавляется). Поэтому x
// обновляется ДО перезаписи vx.
// ============================================================================
// dv_coeff -- коэффициент при dv в обновлении скорости (у Эйлера равен dt; для
// нити берётся из fiber_coeff, чтобы step_fiber и fiber_velocity_update
// использовали одно значение, fiber/fiber_friction_projection.md §4).
template <class N, class XArr, class VxArr, class DvArr>
void euler_step_position_velocity(N n, double dt, double dv_coeff, XArr& x, VxArr& vx, const DvArr& dv)
{
    for (std::size_t d = 0; d < dim; ++d)
        for (N i = 0; i < n; ++i) x[d][i] += dt * vx[d][i]; // x^{n+1} = x^n + dt*v^n (старая скорость)
    for (std::size_t d = 0; d < dim; ++d)
        for (N i = 0; i < n; ++i) vx[d][i] += dv_coeff * dv[d][i];
}

template <class N, class XArr, class VxArr, class DvArr>
void euler_step_position_velocity(N n, double dt, XArr& x, VxArr& vx, const DvArr& dv)
{
    euler_step_position_velocity(n, dt, dt, x, vx, dv);
}

// Скалярная величина (Solid: u, rho, s-компоненты): y^{n+1} = y^n + dt*dy.
template <class N, class Arr, class DArr>
void euler_step_scalar(N n, double dt, Arr& y, const DArr& dy)
{
    for (N i = 0; i < n; ++i) y[i] += dt * dy[i];
}

// ============================================================================
// FiberVelocityUpdate -- формула обновления скорости нити на текущем шаге:
//   v^{n+1} = (*base)[d][i] + h * dv[d][i]
// (fiber/fiber_friction_projection.md §4). Значения берутся из ТОГО ЖЕ места,
// что использует advance() (приватный fiber_coeff/fiber_use_base решателя), и
// нужны стадии трения между вычислением dt и advance().
// base указывает на массивы решателя или FiberState (валидны до advance()).
// ============================================================================
struct FiberVelocityUpdate
{
    const std::vector<double>* base[dim]; // база v^base
    double h;                             // коэффициент при dv
};

// ============================================================================
// SolidVelocityUpdate -- формула обновления скорости частицы Solid на текущем
// шаге (fiber/fiber_solid_friction_implementation.md §4):
//   v^{n+1} = ((*base)[d][i] + h * dvx[d][i]) + (*av)[d][i]     (av != nullptr)
//   v^{n+1} =  (*base)[d][i] + h * dvx[d][i]                    (av == nullptr)
// -- ровно арифметика correct_solid / step_solid (PC: pc_finalize_velocity_
// position_with_av; Эйлер: euler_step_position_velocity, без XSPH). При
// xsph_mode = position av == nullptr (поправка идёт только в позицию). Нужна
// стадии трения между вычислением dt и advance(). Указатели валидны до
// advance().
// ============================================================================
struct SolidVelocityUpdate
{
    const std::vector<double>* base[dim];       // база v^base
    double h;                                    // коэффициент при dvx
    const std::vector<std::vector<double>>* av;  // поправка XSPH; nullptr -- нет
};

// ============================================================================
// Буферы производных RHS одного шага (fiber/PLAN.md §3.9: «решатель хранит
// производные прошлого шага f^{n-1}») -- то, что раньше было локальными
// std::vector внутри time_integration (dvx, du, drho, ds*, av, а также
// dx/t/tdsdt/ds, нужные только сигнатуре single_step). Владеются решателем,
// живут между итерациями while (память переиспользуется, как и раньше).
// ============================================================================
struct SolidRhsBuffers
{
    std::vector<std::vector<double>> av, dx, dvx;
    std::vector<double> du, drho, ds, t, tdsdt;
    std::vector<double> dsxx, dsyy, dszz, dsxy, dsyz, dsxz;

    void resize(std::size_t nbuf)
    {
        av.assign(dim, std::vector<double>(nbuf, 0.0));
        dx.assign(dim, std::vector<double>(nbuf, 0.0));
        dvx.assign(dim, std::vector<double>(nbuf, 0.0));
        du.assign(nbuf, 0.0);
        drho.assign(nbuf, 0.0);
        ds.assign(nbuf, 0.0);
        t.assign(nbuf, 0.0);
        tdsdt.assign(nbuf, 0.0);
        dsxx.assign(nbuf, 0.0);
        dsyy.assign(nbuf, 0.0);
        dszz.assign(nbuf, 0.0);
        dsxy.assign(nbuf, 0.0);
        dsyz.assign(nbuf, 0.0);
        dsxz.assign(nbuf, 0.0);
    }
};

struct FiberRhsBuffers
{
    std::vector<double> dv[dim];

    void resize(std::size_t nfiber)
    {
        for (std::size_t d = 0; d < dim; ++d) dv[d].assign(nfiber, 0.0);
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
    // Режим XSPH (cfg.xsph_mode): false = velocity (v += av, x += dt*v; прежнее поведение),
    // true = position (v без av, x += dt*(v + av); классический XSPH Монаганa).
    bool xsph_position = false;
#ifdef SPH_ENERGY_LOG
    double* u_clamp_acc = nullptr; // накопитель энергии, внесённой обрезкой u>=0 (SPH_ENERGY_LOG)
#endif

    // nbuf -- размер Solid RHS-буферов (>= ntotal, покрывает виртуальные
    // частицы, как раньше: max(cfg.maxn, ntotal)); ntotal/nfiber -- число
    // реальных частиц Solid/Fiber.
    void reset(std::size_t nbuf, std::size_t ntotal, std::size_t nfiber)
    {
        solid_rhs.resize(nbuf);
        fiber_rhs.resize(nfiber);
        u_min_.assign(ntotal, 0.0);
        rho_min_.assign(ntotal, 0.0);
        sxx_min_.assign(ntotal, 0.0);
        syy_min_.assign(ntotal, 0.0);
        szz_min_.assign(ntotal, 0.0);
        sxy_min_.assign(ntotal, 0.0);
        syz_min_.assign(ntotal, 0.0);
        sxz_min_.assign(ntotal, 0.0);
        for (std::size_t d = 0; d < dim; ++d) v_min_[d].assign(ntotal, 0.0);
        for (std::size_t d = 0; d < dim; ++d) v_min_fiber_[d].assign(nfiber, 0.0);
        first_step_ = true;
    }

    bool first_step() const { return first_step_; }
    void set_first_step(bool v) { first_step_ = v; } // восстановление из чекпоинта

    // Предиктор Solid на полшага dt (по f^{n-1}, уже лежащей в solid_rhs с
    // прошлого вызова). Нет-оп на первом шаге -- fiber/SPH.md, «Первый шаг»:
    // предиктора нет вовсе.
    void predict_solid(ParticleSystem& sys, std::size_t ntotal, double dt)
    {
        if (first_step_) return;
        for (std::size_t i = 0; i < ntotal; ++i)
        {
            u_min_[i] = sys.u[i];
            double temp_u = 0.0;
            if (dim == 1) temp_u = -nsym * sys.p[i] * sys.vx[0][i] / sys.x[0][i] / sys.rho[i];
            sys.u[i] += (dt / 2.0) * (solid_rhs.du[i] + temp_u);
            if (sys.u[i] < 0) sys.u[i] = 0.0;

            if (!summation_density)
            {
                rho_min_[i] = sys.rho[i];
                double temp_rho = 0.0;
                if (dim == 1) temp_rho = -nsym * sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
                sys.rho[i] += (dt / 2.0) * (solid_rhs.drho[i] + temp_rho);
            }

            sxx_min_[i] = sys.sxx[i];
            sxy_min_[i] = sys.sxy[i];
            syy_min_[i] = sys.syy[i];
            syz_min_[i] = sys.syz[i];
            sxz_min_[i] = sys.sxz[i];
            szz_min_[i] = sys.szz[i];
            sys.sxx[i] += (dt / 2.0) * solid_rhs.dsxx[i];
            sys.syy[i] += (dt / 2.0) * solid_rhs.dsyy[i];
            sys.szz[i] += (dt / 2.0) * solid_rhs.dszz[i];
            sys.sxy[i] += (dt / 2.0) * solid_rhs.dsxy[i];
            sys.syz[i] += (dt / 2.0) * solid_rhs.dsyz[i];
            sys.sxz[i] += (dt / 2.0) * solid_rhs.dsxz[i];
        }
        pc_predict_velocity_half_step(ntotal, dt, sys.vx, solid_rhs.dvx, v_min_);
    }

    // Предиктор Fiber -- только скорость (нет u/rho/s у Fiber).
    void predict_fiber(fiber::FiberState& fb, std::size_t nfiber, double dt)
    {
        if (first_step_) return;
        pc_predict_velocity_half_step(nfiber, dt, fb.vx, fiber_rhs.dv, v_min_fiber_);
    }

    void correct_solid(ParticleSystem& sys, std::size_t ntotal, double dt)
    {
        const bool fs = first_step_;
        for (std::size_t i = 0; i < ntotal; ++i)
        {
            double temp_u = 0.0;
            if (dim == 1) temp_u = -nsym * sys.p[i] * sys.vx[0][i] / sys.x[0][i] / sys.rho[i];

            if (fs)
                sys.u[i] += (dt / 2.0) * (solid_rhs.du[i] + temp_u);
            else
                sys.u[i] = u_min_[i] + dt * (solid_rhs.du[i] + temp_u);
#ifdef SPH_ENERGY_LOG
            if (sys.u[i] < 0)
            {
                if (u_clamp_acc) *u_clamp_acc += energy_mass(sys, i) * (-sys.u[i]);
                sys.u[i] = 0.0;
            }
#else
            if (sys.u[i] < 0) sys.u[i] = 0.0;
#endif

            if (!summation_density)
            {
                double temp_rho = 0.0;
                if (dim == 1) temp_rho = -nsym * sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
                if (fs)
                    sys.rho[i] += (dt / 2.0) * (solid_rhs.drho[i] + temp_rho);
                else
                    sys.rho[i] = rho_min_[i] + dt * (solid_rhs.drho[i] + temp_rho);
            }

            if (fs)
            {
                sys.sxx[i] += (dt / 2.0) * solid_rhs.dsxx[i];
                sys.syy[i] += (dt / 2.0) * solid_rhs.dsyy[i];
                sys.szz[i] += (dt / 2.0) * solid_rhs.dszz[i];
                sys.sxy[i] += (dt / 2.0) * solid_rhs.dsxy[i];
                sys.syz[i] += (dt / 2.0) * solid_rhs.dsyz[i];
                sys.sxz[i] += (dt / 2.0) * solid_rhs.dsxz[i];
            }
            else
            {
                sys.sxx[i] = sxx_min_[i] + dt * solid_rhs.dsxx[i];
                sys.syy[i] = syy_min_[i] + dt * solid_rhs.dsyy[i];
                sys.szz[i] = szz_min_[i] + dt * solid_rhs.dszz[i];
                sys.sxy[i] = sxy_min_[i] + dt * solid_rhs.dsxy[i];
                sys.syz[i] = syz_min_[i] + dt * solid_rhs.dsyz[i];
                sys.sxz[i] = sxz_min_[i] + dt * solid_rhs.dsxz[i];
            }
        }
        // Скорость/координата -- через общий шаблон (axis-ветвление старого
        // кода арифметически совпадало в обеих ветках if/else, поэтому
        // безопасно заменено на безусловный цикл ещё в перегрузке
        // time_integration(fiber::Model&, ...), этап Ф3; здесь -- то же
        // самое, проверено tests/test_fiber_time_integration_equivalence.cpp).
        if (xsph_position)
            pc_finalize_velocity_position_xsph_position(ntotal, dt, solid_coeff(dt), sys.x, sys.vx, solid_rhs.dvx,
                                                         v_min_, /*use_base=*/solid_use_base(), solid_rhs.av);
        else
            pc_finalize_velocity_position_with_av(ntotal, dt, solid_coeff(dt), sys.x, sys.vx, solid_rhs.dvx, v_min_,
                                                   /*use_base=*/solid_use_base(), solid_rhs.av);
    }

    void correct_fiber(fiber::FiberState& fb, std::size_t nfiber, double dt)
    {
        pc_finalize_velocity_position(nfiber, dt, fiber_coeff(dt), fb.x, fb.vx, fiber_rhs.dv, v_min_fiber_,
                                       /*use_base=*/fiber_use_base());
    }

    // База и коэффициент обновления скорости нити ровно как в correct_fiber
    // (fiber/fiber_friction_projection.md §4). Вызывать после predict() и до
    // advance() того же шага. Первый шаг: base = fb.vx, h = dt/2; шаг n>=2:
    // base = v_min_fiber_ (v^n), h = dt.
    FiberVelocityUpdate fiber_velocity_update(const fiber::FiberState& fb, double dt) const
    {
        FiberVelocityUpdate upd{};
        for (std::size_t d = 0; d < dim; ++d) upd.base[d] = fiber_use_base() ? &v_min_fiber_[d] : &fb.vx[d];
        upd.h = fiber_coeff(dt);
        return upd;
    }

    // База, коэффициент и XSPH-поправка обновления скорости частицы Solid ровно
    // как в correct_solid (fiber/fiber_solid_friction_implementation.md §4).
    // Вызывать после predict() и до advance() того же шага. Первый шаг:
    // base = sys.vx, h = dt/2; шаг n>=2: base = v_min_ (v^n), h = dt; av --
    // solid_rhs.av (XSPH).
    SolidVelocityUpdate solid_velocity_update(const ParticleSystem& sys, double dt) const
    {
        SolidVelocityUpdate upd{};
        for (std::size_t d = 0; d < dim; ++d) upd.base[d] = solid_use_base() ? &v_min_[d] : &sys.vx[d];
        upd.h = solid_coeff(dt);
        // xsph_mode = position: av в скорость не входит (только в позицию) -- nullptr,
        // трение Solid-Fiber и лог энергии (D_xsph = 0) av не видят.
        upd.av = xsph_position ? nullptr : &solid_rhs.av;
        return upd;
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
    void predict(ParticleSystem& sys, std::size_t ntotal, fiber::FiberState& fb, std::size_t nfiber, double dt)
    {
        predict_solid(sys, ntotal, dt);
        predict_fiber(fb, nfiber, dt);
    }

    void advance(ParticleSystem& sys, std::size_t ntotal, fiber::FiberState& fb, std::size_t nfiber, double dt)
    {
        correct_solid(sys, ntotal, dt);
        correct_fiber(fb, nfiber, dt);
        mark_step_done();
    }

private:
    // Единственное место, где задаются коэффициент и признак базы обновления
    // скорости нити (используются correct_fiber и fiber_velocity_update).
    double fiber_coeff(double dt) const { return first_step_ ? dt / 2.0 : dt; }
    bool fiber_use_base() const { return !first_step_; }
    // То же для Solid (correct_solid и solid_velocity_update).
    double solid_coeff(double dt) const { return first_step_ ? dt / 2.0 : dt; }
    bool solid_use_base() const { return !first_step_; }

    bool first_step_ = true;
    std::vector<double> u_min_, rho_min_;
    std::vector<double> sxx_min_, syy_min_, szz_min_, sxy_min_, syz_min_, sxz_min_;
    std::vector<double> v_min_[dim];
    std::vector<double> v_min_fiber_[dim];
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
#ifdef SPH_ENERGY_LOG
    double* u_clamp_acc = nullptr; // накопитель энергии, внесённой обрезкой u>=0 (SPH_ENERGY_LOG)
#endif

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

    // Нет-оп: у Эйлера нет признака первого шага (паритет интерфейса).
    void set_first_step(bool) {}
    bool first_step() const { return false; }

    void step_solid(ParticleSystem& sys, std::size_t ntotal, double dt)
    {
        for (std::size_t i = 0; i < ntotal; ++i)
        {
            double temp_u = 0.0;
            if (dim == 1) temp_u = -nsym * sys.p[i] * sys.vx[0][i] / sys.x[0][i] / sys.rho[i];
            sys.u[i] += dt * (solid_rhs.du[i] + temp_u);
#ifdef SPH_ENERGY_LOG
            if (sys.u[i] < 0)
            {
                if (u_clamp_acc) *u_clamp_acc += energy_mass(sys, i) * (-sys.u[i]);
                sys.u[i] = 0.0;
            }
#else
            if (sys.u[i] < 0) sys.u[i] = 0.0;
#endif

            if (!summation_density)
            {
                double temp_rho = 0.0;
                if (dim == 1) temp_rho = -nsym * sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
                sys.rho[i] += dt * (solid_rhs.drho[i] + temp_rho);
            }

            sys.sxx[i] += dt * solid_rhs.dsxx[i];
            sys.syy[i] += dt * solid_rhs.dsyy[i];
            sys.szz[i] += dt * solid_rhs.dszz[i];
            sys.sxy[i] += dt * solid_rhs.dsxy[i];
            sys.syz[i] += dt * solid_rhs.dsyz[i];
            sys.sxz[i] += dt * solid_rhs.dsxz[i];
        }
        euler_step_position_velocity(ntotal, dt, solid_coeff(dt), sys.x, sys.vx, solid_rhs.dvx);
    }

    void step_fiber(fiber::FiberState& fb, std::size_t nfiber, double dt)
    {
        euler_step_position_velocity(nfiber, dt, fiber_coeff(dt), fb.x, fb.vx, fiber_rhs.dv);
    }

    // Эйлер: base = fb.vx (v^n), h = fiber_coeff(dt) = dt -- ровно то, что
    // делает euler_step_position_velocity (vx += dt*dv).
    FiberVelocityUpdate fiber_velocity_update(const fiber::FiberState& fb, double dt) const
    {
        FiberVelocityUpdate upd{};
        for (std::size_t d = 0; d < dim; ++d) upd.base[d] = &fb.vx[d];
        upd.h = fiber_coeff(dt);
        return upd;
    }

    // Эйлер: base = sys.vx (v^n), h = solid_coeff(dt) = dt, XSPH нет
    // (av == nullptr) -- ровно euler_step_position_velocity в step_solid.
    SolidVelocityUpdate solid_velocity_update(const ParticleSystem& sys, double dt) const
    {
        SolidVelocityUpdate upd{};
        for (std::size_t d = 0; d < dim; ++d) upd.base[d] = &sys.vx[d];
        upd.h = solid_coeff(dt);
        upd.av = nullptr;
        return upd;
    }

    // ------------------------------------------------------------------
    // Единообразный интерфейс OdeSolver::step (см. комментарий у
    // PredictorCorrectorSolver::predict/advance выше). У Эйлера предикторной
    // фазы нет -- RHS вычисляется прямо в y^n, поэтому predict() -- нет-оп.
    // ------------------------------------------------------------------
    void predict(ParticleSystem& /*sys*/, std::size_t /*ntotal*/, fiber::FiberState& /*fb*/, std::size_t /*nfiber*/,
                 double /*dt*/)
    {
        // нет-оп: явный Эйлер вычисляет f = RHS(y^n) без промежуточного
        // предиктора (fiber/PLAN.md §3.9, «Явный Эйлер»).
    }

    void advance(ParticleSystem& sys, std::size_t ntotal, fiber::FiberState& fb, std::size_t nfiber, double dt)
    {
        step_solid(sys, ntotal, dt);
        step_fiber(fb, nfiber, dt);
    }

private:
    double fiber_coeff(double dt) const { return dt; }
    double solid_coeff(double dt) const { return dt; }
};

} // namespace ode
