#pragma once
// strength/johnson_cook.hpp
// Модель прочности и разрушения Джонсона-Кука (JC): чистые функции без
// ParticleSystem / Materials / param.hpp (notes/johnson-cook/jc-implementation-plan.md,
// разделы 2.1, 2.2, 2.5, 3.1). Единицы СИ: Па, К, с^-1.
// Запись со слайда лекции [1 - T/T0]^m НЕ используется (jc-analysis.md, §7).
#include <algorithm>
#include <cmath>
#include <limits>

namespace jc
{

enum class RateForm { FirstOrder, Exponential };
enum class RateSource { Total, Plastic };
enum class DamageRateForm { JC, Borvik };
enum class Heating { Energy, TaylorQuinney };

// Прочность (§2.1). G нужен для возврата на поверхность текучести.
struct JcStrength
{
    double G = 0.0, A = 0.0, B = 0.0, n = 1.0, C = 0.0, m = 1.0;
    double Tm = 1.0, Tr = 0.0, edot0 = 1.0;
    RateForm rate_form = RateForm::FirstOrder;
    RateSource rate_source = RateSource::Total;
};

// Разрушение (§2.5). Tm_f <= 0 означает "взять Tm из прочности".
struct JcDamage
{
    bool enabled = false;
    double D1 = 0.0, D2 = 0.0, D3 = 0.0, D4 = 0.0, D5 = 0.0;
    double Tm_f = 0.0, edot0_f = 1.0, eps_f_min = 1.0e-6;
    DamageRateForm form = DamageRateForm::JC;
};

struct JcThermal
{
    double cp = 1.0;
    Heating heating = Heating::Energy;
    double chi = 0.9; // доля Тейлора-Куинни
};

// Р6: константа Di >= 1e20 (плейсхолдер Autodyn 1.01e20) = "не задана".
inline bool jc_const_disabled(double Di) { return Di >= 1.0e20; }

// Разрушение включено, только если все D1..D5 заданы (Р6) и флаг enabled.
inline bool jc_damage_active(const JcDamage& d)
{
    return d.enabled && !jc_const_disabled(d.D1) && !jc_const_disabled(d.D2) && !jc_const_disabled(d.D3) &&
           !jc_const_disabled(d.D4) && !jc_const_disabled(d.D5);
}

// Гомологическая температура T* = (T - Tr)/(Tm - Tr) (без ограничений).
inline double jc_hom_temperature(double T, double Tr, double Tm)
{
    const double d = Tm - Tr;
    return d > 0.0 ? (T - Tr) / d : 0.0;
}

// Множитель скорости R(edot); при edot < edot0 ровно 1.
inline double rate_factor(const JcStrength& s, double edot)
{
    const double x = edot / s.edot0;
    if (!(x > 1.0)) return 1.0;
    return s.rate_form == RateForm::FirstOrder ? 1.0 + s.C * std::log(x) : std::pow(x, s.C);
}

// Температурный множитель Theta(T): 1 / 1 - T*^m / 0 (расплав).
inline double thermal_factor(const JcStrength& s, double T)
{
    const double ts = jc_hom_temperature(T, s.Tr, s.Tm);
    if (!(ts > 0.0)) return 1.0;
    if (ts >= 1.0) return 0.0;
    return 1.0 - std::pow(ts, s.m);
}

// Предел текучести Y = [A + B eps_p^n] R Theta.
inline double jc_yield(const JcStrength& s, double eps_p, double edot, double T)
{
    const double ep = std::max(eps_p, 0.0);
    const double hard = s.A + (s.B != 0.0 ? s.B * std::pow(ep, s.n) : 0.0);
    return hard * rate_factor(s, edot) * thermal_factor(s, T);
}

struct JcReturnResult
{
    double scale = 1.0;   // s = s_tr * scale
    double d_eps_p = 0.0; // приращение пластической деформации
    double Y = 0.0;       // новый предел текучести (= scale * sigma_eq_trial)
    int iters = 0;
    bool converged = true;
};

// Возврат на поверхность текучести с упрочнением (§2.2): корень
// f(D) = sig_tr - 3 G D - Y(eps_p + D) на [0, sig_tr/(3G)], безопасный Ньютон
// (при выходе из скобки или нефинитной производной - шаг бисекции).
inline JcReturnResult jc_radial_return(const JcStrength& s, double sigma_eq_trial, double eps_p, double edot, double T)
{
    JcReturnResult r;
    const double Y0 = jc_yield(s, eps_p, edot, T);
    r.Y = sigma_eq_trial;
    if (!(sigma_eq_trial > Y0)) // упруго
    {
        r.Y = Y0;
        return r;
    }
    if (thermal_factor(s, T) <= 0.0) // расплав: девиатор обнуляется
    {
        r.scale = 0.0;
        r.d_eps_p = 0.0;
        r.Y = 0.0;
        return r;
    }
    const double g3 = 3.0 * s.G;
    const double RT = rate_factor(s, edot) * thermal_factor(s, T);
    auto f = [&](double d) { return sigma_eq_trial - g3 * d - jc_yield(s, eps_p + d, edot, T); };
    auto df = [&](double d)
    {
        const double e = eps_p + d;
        const double dY = (s.B != 0.0 && e > 0.0) ? s.B * s.n * std::pow(e, s.n - 1.0) * RT
                                                  : (s.B != 0.0 ? std::numeric_limits<double>::infinity() : 0.0);
        return -g3 - dY;
    };
    double lo = 0.0, hi = sigma_eq_trial / g3; // f(lo) > 0, f(hi) <= 0
    double d = 0.5 * (lo + hi);
    double fd = f(d);
    const double tol = 1.0e-12 * Y0;
    int it = 0;
    bool ok = false;
    for (; it < 50; ++it)
    {
        if (std::fabs(fd) <= tol) { ok = true; break; }
        if (fd > 0.0) lo = d; else hi = d;
        const double dfd = df(d);
        double dn = d - fd / dfd;
        if (!std::isfinite(dfd) || !std::isfinite(dn) || !(dn > lo) || !(dn < hi)) dn = 0.5 * (lo + hi);
        d = dn;
        fd = f(d);
    }
    if (!ok)
    {
        // Откат: чистая бисекция до 1e-14 относительной ширины скобки.
        for (int k = 0; k < 200 && (hi - lo) > 1.0e-14 * hi; ++k)
        {
            d = 0.5 * (lo + hi);
            fd = f(d);
            if (fd > 0.0) lo = d; else hi = d;
            ++it;
        }
        d = 0.5 * (lo + hi);
        ok = (hi - lo) <= 1.0e-14 * hi;
    }
    r.d_eps_p = d;
    r.Y = sigma_eq_trial - g3 * d;
    r.scale = r.Y / sigma_eq_trial;
    r.iters = it;
    r.converged = ok;
    return r;
}

// Деформация разрушения eps_f (§2.5). triax = sigma* = -p/sigma_eq.
// JC-форма: 1 + D4 ln(max(edot/edot0f, 1)); Borvik: (1 + edot/edot0f)^D4.
inline double jc_fracture_strain(const JcDamage& d, const JcStrength& s, double triax, double edot, double T)
{
    const double ed = std::max(edot, 0.0) / d.edot0_f;
    const double rate = d.form == DamageRateForm::Borvik ? std::pow(1.0 + ed, d.D4)
                                                         : 1.0 + d.D4 * std::log(std::max(ed, 1.0));
    const double Tm = d.Tm_f > 0.0 ? d.Tm_f : s.Tm;
    // Ниже T_r температурный множитель не увеличивает и не уменьшает eps_f.
    const double ts = std::max(jc_hom_temperature(T, s.Tr, Tm), 0.0);
    const double ef = (d.D1 + d.D2 * std::exp(d.D3 * triax)) * rate * (1.0 + d.D5 * ts);
    return std::max(ef, d.eps_f_min);
}

// Геометрия решателя (caller сопоставляет dim/axis_sim).
enum class Geom { D1, D2Plane, D2Axi, D3 };

// Эквивалентная скорость деформации sqrt(2/3 e':e'), e' - девиатор.
// D1: одноосная деформация diag(exx,0,0); D2Plane: ezz = 0, сдвиг только exy;
// D2Axi: exx, eyy, ezz (окружная), exy; D3: все шесть. Тензорные сдвиги: 2 e_ij^2.
inline double equivalent_strain_rate(Geom g, double exx, double eyy, double ezz, double exy, double eyz, double exz)
{
    if (g == Geom::D1) { eyy = 0.0; ezz = 0.0; exy = 0.0; }
    if (g == Geom::D2Plane) ezz = 0.0;
    if (g != Geom::D3) { eyz = 0.0; exz = 0.0; }
    const double tr3 = (exx + eyy + ezz) / 3.0;
    const double a = exx - tr3, b = eyy - tr3, c = ezz - tr3;
    const double dd = a * a + b * b + c * c + 2.0 * (exy * exy + eyz * eyz + exz * exz);
    return std::sqrt(2.0 / 3.0 * dd);
}

// J2 девиатора напряжений - в точности алгебра sigma_i (src/ppochnost.cpp),
// sigma_eq = sqrt(3 J2).
inline double j2_of_deviator(Geom g, double sxx, double syy, double sxy, double syz, double sxz)
{
    if (g == Geom::D1) return 0.75 * sxx * sxx;
    if (g == Geom::D2Plane) return sxx * sxx + sxy * sxy;
    if (g == Geom::D2Axi) return sxx * sxx + sxx * syy + syy * syy + sxy * sxy;
    return sxx * sxx + sxx * syy + syy * syy + sxy * sxy + sxz * sxz + syz * syz;
}

} // namespace jc
