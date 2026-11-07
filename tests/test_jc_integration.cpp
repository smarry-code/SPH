// tests/test_jc_integration.cpp
// Интеграционные тесты прочности и разрушения Джонсона-Кука I1-I8
// (fiber/jc-implementation-plan.md, §5.3): одна-две частицы, без SPH-соседей,
// через interactions::apply_material_strength (ветка JC -> apply_jc_strength).
//
// Драйвер. Плотность фиксирована (rho = rho0), на каждом шаге задаются
// скорости деформации sys.exx.. (rxy.. = 0, вращения нет). Пробный девиатор
// считается ровно как закон Гука в finalize_solid (solid_sph.hpp):
//   dsxx = 2G (exx + drho/(3 rho)), ..., dsxy = 2G exy, ...
// с drho = -rho * tr(e) (согласованная с кинематикой скорость плотности; для
// isochoric-нагружений drho = 0, и это ровно девиаторная скорость). G берётся
// через effective_shear_modulus (0 для разрушенной частицы, как в решателе).
// Затем s += ds*dt и один вызов apply_material_strength. Давление и энергия
// задаются тестом (p, u), уравнения состояния и энергии не интегрируются.
// I1 -- сборка SPH_DIM=1, I2-I9 -- SPH_DIM=3, I5 -- только при SPH_ENERGY_LOG.
#include "interactions/postprocessing.hpp"

#include "eos.hpp"
#include "materials.hpp"
#include "param.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1e-300, std::fabs(b)); }

const SimulationConfig& cfg()
{
    static const SimulationConfig c;
    return c;
}

const JcMaterial& mat(int itype) { return jc_material(cfg().mats, itype); }

ParticleSystem make_sys(std::size_t n, int itype)
{
    ParticleSystem sys;
    sys.resize(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        sys.inzializ_index(i);
        sys.itype[i] = itype;
        sys.mass[i] = 1.0;
        sys.rho[i] = sys.rho_0[i] = mat(itype).rho0;
        sys.T[i] = mat(itype).s.Tr;
    }
    return sys;
}

// Шаг драйвера: закон Гука без поворота + apply_material_strength.
// dt, e* уже записаны в sys.exx..; drho = -rho*tr(e).
void step(ParticleSystem& sys, double dt, double* acc = nullptr, const SimulationConfig& c = cfg())
{
    const std::size_t n = sys.size();
    std::vector<double> drho(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double G = effective_shear_modulus(c.mats, sys, i);
        if constexpr (dim == 1)
        {
            drho[i] = -sys.rho[i] * sys.exx[i];
            const double dv = drho[i] / (3.0 * sys.rho[i]);
            sys.sxx[i] += dt * 2.0 * G * (sys.exx[i] + dv);
            sys.syy[i] += dt * 2.0 * G * dv;
            sys.szz[i] += dt * 2.0 * G * dv;
        }
        else
        {
            drho[i] = -sys.rho[i] * (sys.exx[i] + sys.eyy[i] + sys.ezz[i]);
            const double dv = drho[i] / (3.0 * sys.rho[i]);
            sys.sxx[i] += dt * 2.0 * G * (sys.exx[i] + dv);
            sys.syy[i] += dt * 2.0 * G * (sys.eyy[i] + dv);
            sys.szz[i] += dt * 2.0 * G * (sys.ezz[i] + dv);
            sys.sxy[i] += dt * 2.0 * G * sys.exy[i];
            sys.syz[i] += dt * 2.0 * G * sys.eyz[i];
            sys.sxz[i] += dt * 2.0 * G * sys.exz[i];
        }
    }
#ifdef SPH_ENERGY_LOG
    interactions::apply_material_strength(c, sys, n, dt, drho, acc);
#else
    (void)acc;
    interactions::apply_material_strength(c, sys, n, dt, drho);
#endif
}

// Эквивалентное напряжение 3D-девиатора, независимо от модуля jc.
double sigma_eq3(const ParticleSystem& s, std::size_t i)
{
    return std::sqrt(1.5 * (s.sxx[i] * s.sxx[i] + s.syy[i] * s.syy[i] + s.szz[i] * s.szz[i]) +
                     3.0 * (s.sxy[i] * s.sxy[i] + s.syz[i] * s.syz[i] + s.sxz[i] * s.sxz[i]));
}

// ---------------------------------------------------------------- I1
void test_i1()
{
    std::printf("I1 одноосная деформация, dim=1\n");
    const JcMaterial& m = mat(12);
    ParticleSystem sys = make_sys(1, 12);
    const double edot = 100.0; // эквивалентная скорость; exx = 1.5*edot
    const double dt = 1.0e-5;
    bool ok_sig = true, ok_ep = true;
    double max_sig = 0.0, max_ep = 0.0;
    int n_plastic = 0;
    for (int k = 1; k <= 400; ++k)
    {
        sys.exx[0] = 1.5 * edot;
        step(sys, dt);
        if (sys.e_pl[0] > 0.0)
        {
            ++n_plastic;
            // |sxx| = (2/3) Y(eps_p, edot, T)
            const double Y = jc::jc_yield(m.s, sys.e_pl[0], edot, m.s.Tr);
            max_sig = std::max(max_sig, std::fabs(std::fabs(sys.sxx[0]) - 2.0 / 3.0 * Y) / Y);
            // Аналитика: eps_p + Y(eps_p)/(3G) = edot*t (бисекция).
            const double rhs = edot * dt * k;
            double lo = 0.0, hi = rhs;
            for (int it = 0; it < 200; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                (mid + jc::jc_yield(m.s, mid, edot, m.s.Tr) / (3.0 * m.s.G) > rhs ? hi : lo) = mid;
            }
            max_ep = std::max(max_ep, std::fabs(sys.e_pl[0] - 0.5 * (lo + hi)) / (0.5 * (lo + hi)));
        }
    }
    ok_sig = max_sig < 1.0e-6;
    ok_ep = max_ep < 1.0e-6;
    check(n_plastic > 300, "пластический режим достигнут");
    check(ok_sig, "|sxx| = (2/3) Y на каждом шаге, отн. ошибка " + std::to_string(max_sig));
    check(ok_ep, "eps_p совпадает с аналитическим интегралом, отн. ошибка " + std::to_string(max_ep));
}

// ---------------------------------------------------------------- I2-I4, I8: чистый сдвиг (3D)
// Сдвиговая скорость exy: эквивалентная edot = exy*sqrt(4/3).
double shear_edot(double exy) { return exy * std::sqrt(4.0 / 3.0); }

void test_i2()
{
    std::printf("I2 чистый сдвиг, 3D\n");
    const JcMaterial& m = mat(12);
    ParticleSystem sys = make_sys(1, 12);
    const double exy = 50.0, dt = 1.0e-5;
    const double edot = shear_edot(exy);
    double max_tau = 0.0, max_tr = 0.0;
    for (int k = 0; k < 300; ++k)
    {
        sys.exy[0] = exy;
        step(sys, dt);
        if (sys.e_pl[0] > 0.0)
        {
            const double Y = jc::jc_yield(m.s, sys.e_pl[0], edot, m.s.Tr);
            max_tau = std::max(max_tau, std::fabs(std::fabs(sys.sxy[0]) - Y / std::sqrt(3.0)) / Y);
        }
        max_tr = std::max(max_tr, std::fabs(sys.sxx[0] + sys.syy[0] + sys.szz[0]) / m.s.A);
    }
    check(sys.e_pl[0] > 0.01, "пластическая деформация накоплена");
    check(max_tau < 1.0e-9, "tau = Y/sqrt3 после текучести, отн. ошибка " + std::to_string(max_tau));
    check(max_tr < 1.0e-12, "девиатор бесследов");
}

void test_i3()
{
    std::printf("I3 скоростная чувствительность\n");
    const JcMaterial& m = mat(12);
    double ratio_src[2] = {0.0, 0.0};
    const double rates[2] = {1.0, 1.0e4};
    for (int c = 0; c < 2; ++c)
    {
        ParticleSystem sys = make_sys(1, 12);
        const double exy = rates[c] / std::sqrt(4.0 / 3.0);
        const double dt = 0.05 / rates[c] / 500.0; // eps_p ~ 0.05 за 500 шагов
        for (int k = 0; k < 500; ++k)
        {
            sys.exy[0] = exy;
            step(sys, dt);
        }
        // sigma_eq / (A + B eps_p^n) = R(edot) при любом eps_p.
        const double hard = m.s.A + m.s.B * std::pow(sys.e_pl[0], m.s.n);
        ratio_src[c] = sigma_eq3(sys, 0) / hard;
        check(sys.e_pl[0] > 0.03, "edot=" + std::to_string(rates[c]) + ": eps_p=" + std::to_string(sys.e_pl[0]));
    }
    // Y при равном eps_p отличается только множителем R: 1 + C ln(1e4) = 1.128945.
    check(near_rel(ratio_src[1] / ratio_src[0], 1.128945, 1.0e-6),
          "Y(1e4)/Y(1) = " + std::to_string(ratio_src[1] / ratio_src[0]));
}

void test_i4()
{
    std::printf("I4 температура и расплав\n");
    const JcMaterial& m = mat(12);
    ParticleSystem sys = make_sys(1, 12);
    const double exy = 10.0, dt = 1.0e-4;
    const double edot = shear_edot(exy);
    sys.u[0] = m.th.cp * 500.0;
    for (int k = 0; k < 100; ++k)
    {
        sys.exy[0] = exy;
        step(sys, dt);
    }
    const double T = m.s.Tr + sys.u[0] / m.th.cp;
    check(sys.T[0] == T, "sys.T = Tr + u/cp (T=" + std::to_string(sys.T[0]) + ")");
    const double Y = jc::jc_yield(m.s, sys.e_pl[0], edot, T);
    check(near_rel(sigma_eq3(sys, 0), Y, 1.0e-9), "Y учитывает T (термическое размягчение)");
    check(Y < jc::jc_yield(m.s, sys.e_pl[0], edot, m.s.Tr), "Y(T) < Y(Tr)");
    // T >= Tm: девиатор обнуляется.
    sys.u[0] = m.th.cp * (m.s.Tm - m.s.Tr) * 1.01;
    sys.exy[0] = exy;
    step(sys, dt);
    check(sys.T[0] >= m.s.Tm, "T >= Tm");
    check(sys.sxx[0] == 0.0 && sys.syy[0] == 0.0 && sys.szz[0] == 0.0 && sys.sxy[0] == 0.0 && sys.syz[0] == 0.0 &&
              sys.sxz[0] == 0.0,
          "при T >= Tm девиатор = 0");
}

// ---------------------------------------------------------------- I9: нагрев Тейлора-Куинни (3D)
// heating = TaylorQuinney: T интегрируется отдельно, u не трогается;
// T - Tr = chi * sum(Y_k * d_eps_p,k) / (rho * cp), Y_k -- напряжение текучести
// после шага (девиатор лежит на поверхности текучести), rho = const.
void test_i9()
{
    std::printf("I9 нагрев Тейлора-Куинни\n");
    const JcMaterial& m = mat(12);
    const double exy = 50.0, dt = 1.0e-5;
    for (int variant = 0; variant < 2; ++variant)
    {
        const double chi = variant == 0 ? 0.9 : 0.0;
        SimulationConfig c;
        for (JcMaterial& jm : c.mats.jc)
        {
            jm.th.heating = jc::Heating::TaylorQuinney;
            jm.th.chi = chi;
        }
        ParticleSystem sys = make_sys(1, 12);
        sys.u[0] = 123.0; // в режиме TQ энергия не влияет на T и не меняется
        const double u0 = sys.u[0];
        double sum_Yd = 0.0, T_prev = sys.T[0];
        bool monotone = true;
        for (int k = 0; k < 300; ++k)
        {
            const double ep_before = sys.e_pl[0];
            sys.exy[0] = exy;
            step(sys, dt, nullptr, c);
            sum_Yd += sigma_eq3(sys, 0) * (sys.e_pl[0] - ep_before);
            if (sys.T[0] < T_prev) monotone = false;
            T_prev = sys.T[0];
        }
        const double dT_ref = chi * sum_Yd / (sys.rho[0] * m.th.cp);
        const std::string tag = "chi=" + std::to_string(chi) + ": ";
        check(sys.e_pl[0] > 0.01, tag + "пластическая деформация накоплена");
        check(sys.u[0] == u0, tag + "u не меняется");
        check(monotone, tag + "T не убывает");
        if (variant == 0)
        {
            check(dT_ref > 0.0 && near_rel(sys.T[0] - m.s.Tr, dT_ref, 1.0e-9),
                  tag + "T - Tr = chi*sum(Y*d_eps_p)/(rho*cp), dT=" + std::to_string(sys.T[0] - m.s.Tr));
        }
        else
        {
            check(sys.T[0] == m.s.Tr, tag + "T == Tr точно");
        }
    }
}

#ifdef SPH_ENERGY_LOG
void test_i5()
{
    std::printf("I5 энергия пластической диссипации\n");
    ParticleSystem sys = make_sys(2, 12);
    sys.mass[0] = 2.5;
    sys.mass[1] = 0.7;
    sys.rho[1] = 1.1 * sys.rho_0[1]; // rho != rho0: проверка деления на rho
    const double exy[2] = {40.0, 120.0}, dt = 1.0e-5;
    double acc = 0.0, ref = 0.0;
    for (int k = 0; k < 300; ++k)
    {
        double ep_before[2] = {sys.e_pl[0], sys.e_pl[1]};
        sys.exy[0] = exy[0];
        sys.exy[1] = exy[1];
        step(sys, dt, &acc);
        for (std::size_t i = 0; i < 2; ++i)
        {
            // Независимо: Y = sigma_eq итогового девиатора (на поверхности текучести).
            ref += sys.mass[i] * sigma_eq3(sys, i) * (sys.e_pl[i] - ep_before[i]) / sys.rho[i];
        }
    }
    check(ref > 0.0, "пластическая работа ненулевая");
    check(near_rel(acc, ref, 1.0e-6), "plastic_acc = sum m*Y*d_eps_p/rho: " + std::to_string(acc) + " vs " + std::to_string(ref));
}
#endif

// ---------------------------------------------------------------- I6, I7: одноосное растяжение, sigma* = 1/3 (3D)
// Несжимаемое одноосное растяжение: exx = e, eyy = ezz = -e/2 (tr = 0, drho = 0),
// edot = e. p задаётся тестом = -sigma_eq/3 (sigma_eq оценивается как Y в конце шага).
void test_i6_i7()
{
    std::printf("I6/I7 разрушение при sigma*=1/3 и разрушенная частица\n");
    const JcMaterial& m = mat(12);
    ParticleSystem sys = make_sys(1, 12);
    const double e = 0.5, dt = 0.01; // edot <= 1: множитель скорости = 1; u = 0: T = Tr
    const double ef = m.d.D1 + m.d.D2 * std::exp(m.d.D3 / 3.0);
    check(near_rel(ef, 1.7469, 1.0e-4), "eps_f(1/3) = " + std::to_string(ef));
    const double dstep = e * dt;
    int k = 0;
    for (; k < 1000 && sys.D[0] < 1.0; ++k)
    {
        sys.exx[0] = e; sys.eyy[0] = -0.5 * e; sys.ezz[0] = -0.5 * e;
        sys.p[0] = -jc::jc_yield(m.s, sys.e_pl[0] + dstep, e, m.s.Tr) / 3.0;
        step(sys, dt);
    }
    check(sys.D[0] == 1.0, "D достигла 1");
    check(std::fabs(sys.e_pl[0] - ef) <= 1.01 * dstep, "eps_p при D=1: " + std::to_string(sys.e_pl[0]) + " (ожидалось " + std::to_string(ef) + " +- " + std::to_string(dstep) + ")");

    // I7: после D = 1.
    check(sys.sxx[0] == 0.0 && sys.syy[0] == 0.0 && sys.szz[0] == 0.0 && sys.sxy[0] == 0.0 && sys.syz[0] == 0.0 &&
              sys.sxz[0] == 0.0,
          "I7: девиатор = 0 сразу при D=1");
    check(effective_shear_modulus(cfg().mats, sys, 0) == 0.0, "I7: effective_shear_modulus = 0");
    check(shear_modulus_by_itype(cfg().mats, 12) == m.s.G, "I7: shear_modulus_by_itype не изменён");
    const double ep_frozen = sys.e_pl[0];
    bool zero = true;
    for (int j = 0; j < 100; ++j)
    {
        sys.exx[0] = e; sys.eyy[0] = -0.5 * e; sys.ezz[0] = -0.5 * e;
        sys.exy[0] = 3.0 * e; // и сдвиг
        step(sys, dt);
        zero = zero && sys.sxx[0] == 0.0 && sys.syy[0] == 0.0 && sys.szz[0] == 0.0 && sys.sxy[0] == 0.0 &&
               sys.syz[0] == 0.0 && sys.sxz[0] == 0.0;
    }
    check(zero, "I7: девиатор остаётся 0 при продолжении нагружения");
    check(sys.e_pl[0] == ep_frozen, "I7: eps_p не растёт");
    check(sys.D[0] == 1.0, "I7: D остаётся 1");
    double p_d = 0.0, p_u = 0.0, c = 0.0;
    const double rho_t = 0.99 * m.rho0; // растяжение
    eos_jc(m, rho_t, m.rho0, 0.0, 1.0, p_d, c);
    eos_jc(m, rho_t, m.rho0, 0.0, 0.0, p_u, c);
    check(p_u < 0.0, "I7: целая частица держит растяжение (p=" + std::to_string(p_u) + ")");
    check(p_d >= 0.0, "I7: разрушенная частица p >= 0");
}

void test_i8()
{
    std::printf("I8 разрушение выключено\n");
    for (int it : {13, 14})
    {
        ParticleSystem sys = make_sys(1, it);
        const double exy = 50.0, dt = 1.0e-4;
        for (int k = 0; k < 400; ++k)
        {
            sys.exy[0] = exy;
            sys.p[0] = -1.0e9;
            step(sys, dt);
        }
        check(sys.e_pl[0] > 1.0, "itype " + std::to_string(it) + ": eps_p=" + std::to_string(sys.e_pl[0]));
        check(sys.D[0] == 0.0, "itype " + std::to_string(it) + ": D = 0");
    }
}

} // namespace

int main()
{
    if constexpr (dim == 1)
    {
        test_i1();
    }
    else if constexpr (dim == 3)
    {
        test_i2();
        test_i3();
        test_i4();
        test_i9();
#ifdef SPH_ENERGY_LOG
        test_i5();
#endif
        test_i6_i7();
        test_i8();
    }
    else
    {
        std::printf("test_jc_integration: dim=%d не поддерживается\n", static_cast<int>(dim));
        return 1;
    }
    std::printf("test_jc_integration: %d failures\n", failures);
    return failures == 0 ? 0 : 1;
}
