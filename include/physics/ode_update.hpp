#pragma once
// physics/ode_update.hpp
// Поэлементное обновление состояния решателями (Эйлер, PredictorCorrector) --
// fiber/CUDA_PLAN.md §3.2/§4 Г2, файл из списка §3.2: "формулы предиктора/
// корректора/Эйлера для одной частицы". Построчный перенос арифметики
// include/ode/ode_solver.hpp (PredictorCorrectorSolver::predict_solid/
// correct_solid/predict_fiber/correct_fiber, EulerSolver::step_solid/
// step_fiber и общих шаблонов pc_predict_velocity_half_step/
// pc_finalize_velocity_position[_with_av]/euler_step_position_velocity) в
// SPH_HD-функции над скалярами -- то, что решатели ГПУ (этап Г5) будут
// вызывать на каждую (частицу, компоненту).
//
// Данные Solid/Fiber хранятся SoA по компонентам (ParticleSystem::x[dim][i],
// fiber::FiberState::x[dim][i]) -- поэтому "одна частица" здесь означает один
// СКАЛЯРНЫЙ элемент (частица, компонента) или один скаляр состояния (u, rho,
// s*): те же формулы, что применялись поэлементно в двойных циклах "for d {
// for i {...} }" ode_solver.hpp, вынесены без изменения арифметики. Порядок
// операций (deltaвремени/2 вычисляется как dt_prev/real_t(2.0), затем
// умножается -- НЕ dt_prev*real_t(0.5)) не изменён -- требование побитной
// регрессии (план §3.2).
#include "param.hpp"       // dim, nsym
#include "physics/hd.hpp"  // SPH_HD

namespace physics
{

// ----------------------------------------------------------------------
// Скорость/координата (общие Solid и Fiber -- ode_solver.hpp: v_min/vx/x
// одного и того же вида у ParticleSystem и fiber::FiberState).
// ----------------------------------------------------------------------

// Половина предикторного шага скорости (PredictorCorrectorSolver::
// predict_solid/predict_fiber, pc_predict_velocity_half_step): сохранение
// v_min -- побочный эффект вызывающей стороны (не формула), здесь -- только
// v_new = v + (dt_prev/2)*dv.
SPH_HD inline real_t ode_velocity_half_step(real_t v, real_t dt_prev, real_t dv)
{
    return v + (dt_prev / real_t(2.0)) * dv;
}

// Корректор скорости С поправкой XSPH (Solid, pc_finalize_velocity_position_
// with_av): base уже выбран вызывающей стороной (v_base при корректоре не на
// первом шаге, текущая vx иначе) -- см. комментарий у одноимённого шаблона
// ode_solver.hpp.
SPH_HD inline real_t ode_velocity_with_av(real_t base, real_t coeff, real_t dv, real_t av)
{
    return base + coeff * dv + av;
}

// То же БЕЗ поправки XSPH (Fiber, pc_finalize_velocity_position; Эйлер,
// euler_step_position_velocity -- coeff = dt).
SPH_HD inline real_t ode_velocity(real_t base, real_t coeff, real_t dv)
{
    return base + coeff * dv;
}

// Координата: x_new = x + dt*vx (vx -- УЖЕ ОБНОВЛЁННАЯ на этом шаге скорость,
// вызывающая сторона обязана передать новое значение -- порядок "сначала
// скорость, потом координата" не меняется).
SPH_HD inline real_t ode_position(real_t x, real_t dt, real_t vx)
{
    return x + dt * vx;
}

// ----------------------------------------------------------------------
// Скалярные поля Solid (u, rho, sxx..sxz) -- PredictorCorrectorSolver::
// predict_solid/correct_solid, EulerSolver::step_solid. temp -- осевая
// поправка (dim==1), 0 для sxx..sxz (там её нет в оригинале).
// ----------------------------------------------------------------------

// y_new = y_base + coeff*dy, БЕЗ температурной поправки -- sxx..sxz
// (ode_solver.hpp никогда не прибавляет temp к этим полям, в отличие от u/
// rho ниже; отдельная функция, а не ode_solid_update(...,0), НАМЕРЕННО --
// dy + real_t(0.0) добавил бы операцию сложения, которой в оригинале для
// этих полей нет: для большинства dy это не меняет значение, НО для
// dy == -0.0 (достижимо в симметричных конфигурациях BC2/BC6) (-0.0)+0.0 =
// +0.0 -- отличие на знак нуля, недопустимое при побитной регрессии).
SPH_HD inline real_t ode_solid_advance(real_t y_base, real_t coeff, real_t dy)
{
    return y_base + coeff * dy;
}

// y_new = y_base + coeff*(dy + temp), БЕЗ ограничения (u, rho -- temp здесь
// РЕАЛЬНО участвует уже в оригинале: temp_u/temp_rho -- свежая локальная
// real_t, инициализированная 0 (не результат вычисления, поэтому не может
// быть -0.0), и оригинал уже пишет "solid_rhs.du[i] + temp_u" безусловно,
// даже при dim!=1 -- см. ode_solid_energy_axis_term/
// ode_solid_density_axis_term ниже).
SPH_HD inline real_t ode_solid_update(real_t y_base, real_t coeff, real_t dy, real_t temp)
{
    return y_base + coeff * (dy + temp);
}

// Та же формула + ограничение y_new >= 0 (энергия u -- int_force.cpp/
// ode_solver.hpp: "if (sys.u[i] < 0) sys.u[i] = 0").
SPH_HD inline real_t ode_solid_update_nonneg(real_t y_base, real_t coeff, real_t dy, real_t temp)
{
    real_t y_new = y_base + coeff * (dy + temp);
    if (y_new < 0) y_new = 0;
    return y_new;
}

// Осевая поправка энергии (dim==1): -nsym*p*vx0/x0/rho, 0 иначе. dim/nsym --
// глобальные constexpr (param.hpp); "if (dim == 1)" -- буквально как в
// ode_solver.hpp (runtime по форме, компилятор сворачивает при dim!=1).
SPH_HD inline real_t ode_solid_energy_axis_term(real_t p, real_t vx0, real_t x0, real_t rho)
{
    real_t temp_u = 0;
    if (dim == 1) temp_u = -nsym * p * vx0 / x0 / rho;
    return temp_u;
}

// Осевая поправка плотности (dim==1): -nsym*rho*vx0/x0, 0 иначе.
SPH_HD inline real_t ode_solid_density_axis_term(real_t rho, real_t vx0, real_t x0)
{
    real_t temp_rho = 0;
    if (dim == 1) temp_rho = -nsym * rho * vx0 / x0;
    return temp_rho;
}

} // namespace physics
