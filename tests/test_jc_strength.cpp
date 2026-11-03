// tests/test_jc_strength.cpp
// Unit-тесты чистого модуля Джонсона-Кука (include/strength/johnson_cook.hpp;
// fiber/jc-implementation-plan.md, раздел 5.1, тесты U1-U14). Решатель не
// используется; U10 сверяется с prochnost_mises для геометрии текущей сборки.
#include "strength/johnson_cook.hpp"

#include "materials.hpp"
#include "param.hpp"
#include "ppochnost.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool same(double a, double b) { return !(a < b) && !(a > b); } // точное равенство без -Wfloat-equal

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1e-300, std::fabs(b)); }

// Сталь 4340 (§3.2), СИ.
jc::JcStrength steel4340()
{
    jc::JcStrength s;
    s.G = 81.8e9; s.A = 792e6; s.B = 510e6; s.n = 0.26; s.C = 0.014; s.m = 1.03;
    s.Tm = 1793.0; s.Tr = 300.0; s.edot0 = 1.0;
    return s;
}

jc::JcDamage damage4340()
{
    jc::JcDamage d;
    d.enabled = true; d.D1 = 0.05; d.D2 = 3.44; d.D3 = -2.12; d.D4 = 0.002; d.D5 = 0.61;
    d.Tm_f = 1793.0; d.edot0_f = 1.0;
    return d;
}

jc::JcStrength al5083()
{
    jc::JcStrength s;
    s.G = 26.92e9; s.A = 167e6; s.B = 596e6; s.n = 0.551; s.C = 0.001; s.m = 0.859;
    s.Tm = 893.0; s.Tr = 293.0; s.edot0 = 1.0;
    return s;
}

jc::JcDamage damage5083()
{
    jc::JcDamage d;
    d.enabled = true; d.D1 = 0.178; d.D2 = 0.389; d.D3 = -2.25; d.D4 = 0.147; d.D5 = 0.0;
    d.Tm_f = 893.0; d.edot0_f = 1.0; d.form = jc::DamageRateForm::Borvik;
    return d;
}

void test_yield()
{
    std::printf("U1-U6: предел текучести\n");
    const jc::JcStrength s = steel4340();
    check(near_rel(jc::jc_yield(s, 0.0, s.edot0, s.Tr), 792e6, 1e-14), "U1 Y(0, edot0, Tr) = A");
    const double y2 = jc::jc_yield(s, 0.1, 1.0, 300.0);
    check(near_rel(y2, 792e6 + 510e6 * std::pow(0.1, 0.26), 1e-12), "U2 упрочнение, формула");
    check(near_rel(y2, 1072.266e6, 1e-6), "U2 упрочнение = 1072.266 МПа");
    const double r3 = jc::rate_factor(s, 1e4);
    check(near_rel(r3, 1.0 + 0.014 * std::log(1e4), 1e-12), "U3 скорость, формула");
    check(near_rel(r3, 1.128945, 1e-6), "U3 множитель = 1.128945");
    check(same(jc::rate_factor(s, 0.1), 1.0) && same(jc::rate_factor(s, 0.0), 1.0) && same(jc::rate_factor(s, 1.0), 1.0),
          "U4 edot <= edot0 -> ровно 1");
    const double th = jc::thermal_factor(s, 800.0);
    check(near_rel(th, 1.0 - std::pow(500.0 / 1493.0, 1.03), 1e-12), "U5 Theta(800), формула");
    check(near_rel(th, 0.675916, 1e-6), "U5 Theta(800) = 0.675916");
    check(same(jc::thermal_factor(s, 300.0), 1.0) && same(jc::thermal_factor(s, 200.0), 1.0), "U5 T <= Tr -> 1");
    check(same(jc::thermal_factor(s, 1793.0), 0.0) && same(jc::thermal_factor(s, 2500.0), 0.0), "U5 T >= Tm -> 0");
    jc::JcStrength se = s;
    se.rate_form = jc::RateForm::Exponential;
    check(near_rel(jc::rate_factor(se, 1e4), std::pow(1e4, 0.014), 1e-12), "U6 Exponential = 1e4^C");
    check(same(jc::rate_factor(se, 0.5), 1.0), "U6 Exponential ниже edot0 -> 1");
}

void check_return(const jc::JcStrength& s, double sig_tr, double ep, double edot, double T, const std::string& tag)
{
    const jc::JcReturnResult r = jc::jc_radial_return(s, sig_tr, ep, edot, T);
    const double Y = jc::jc_yield(s, ep + r.d_eps_p, edot, T);
    const double Y0 = jc::jc_yield(s, ep, edot, T);
    const double res = std::fabs(sig_tr - 3.0 * s.G * r.d_eps_p - Y);
    check(r.converged, tag + ": converged");
    check(res <= 1e-10 * Y0 && r.d_eps_p >= 0.0, tag + ": невязка <= 1e-10 Y, dep >= 0");
    check(near_rel(r.scale * sig_tr, Y, 1e-9) && near_rel(r.Y, r.scale * sig_tr, 1e-14), tag + ": scale*sig_tr = Y(ep+dep)");
}

void test_return()
{
    std::printf("U7-U9: radial return\n");
    const jc::JcStrength s = steel4340();
    const jc::JcReturnResult e = jc::jc_radial_return(s, 700e6, 0.0, 1.0, 300.0);
    check(same(e.scale, 1.0) && same(e.d_eps_p, 0.0), "U7 упругий случай: scale = 1, dep = 0");
    check_return(s, 1500e6, 0.0, 1.0, 300.0, "U8 4340 ep=0");
    check_return(s, 1800e6, 0.05, 1e3, 500.0, "U8 4340 ep=0.05, hot");
    check_return(s, 3000e6, 0.0, 1.0, 300.0, "U9 4340 ep=0, большой перенапряг");
    check_return(al5083(), 600e6, 0.0, 1.0, 293.0, "U9 5083 ep=0");
    jc::JcStrength s01 = steel4340();
    s01.n = 0.1;
    check_return(s01, 1500e6, 0.0, 1.0, 300.0, "U9 n = 0.1 ep=0");
    const jc::JcReturnResult m = jc::jc_radial_return(s, 1e9, 0.0, 1.0, 1900.0);
    check(same(m.scale, 0.0), "расплав: scale = 0");
}

// Эталон: prochnost_mises для геометрии текущей сборки.
jc::Geom build_geom()
{
    if constexpr (dim == 1) return jc::Geom::D1;
    else if constexpr (dim == 2) return axis_sim ? jc::Geom::D2Axi : jc::Geom::D2Plane;
    else return jc::Geom::D3;
}

void test_mises()
{
    std::printf("U10: вырожденный JC = Мизес\n");
    jc::JcStrength s = steel4340();
    s.B = 0.0; s.C = 0.0;
    Materials mats;
    const size_t k = 1;
    mats.G[k] = s.G;
    mats.sigy[k] = s.A;
    const size_t nmat = Materials::n_Me + k;
    const jc::Geom g = build_geom();
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    double worst = 0.0, worst_alg = 0.0;
    int n_plastic = 0;
    for (int it = 0; it < 1000; ++it)
    {
        const double amp = 3e9 * (it % 5 == 0 ? 0.1 : 1.0); // часть случаев упругих
        double sxx = amp * u(rng), syy = amp * u(rng), sxy = amp * u(rng), syz = amp * u(rng), sxz = amp * u(rng);
        double szz = -sxx - syy; // девиатор
        if (g == jc::Geom::D1) { syy = -0.5 * sxx; szz = -0.5 * sxx; }
        const double j2 = jc::j2_of_deviator(g, sxx, syy, sxy, syz, sxz);
        double rxx = sxx, ryy = syy, rzz = szz, rxy = sxy, ryz = syz, rxz = sxz;
        prochnost_mises(mats, nmat, rxx, ryy, rzz, rxy, ryz, rxz);
        // sigma_i как эталон алгебры
        const double ref_j2 = sigma_i(sxx, syy, sxy, syz, sxz);
        worst_alg = std::max(worst_alg, std::fabs(j2 - ref_j2) / std::max(1e-300, ref_j2));
        const jc::JcReturnResult r = jc::jc_radial_return(s, std::sqrt(3.0 * j2), 0.0, 0.0, s.Tr);
        n_plastic += r.scale < 1.0 ? 1 : 0;
        const double xs[6] = {sxx, syy, szz, sxy, syz, sxz};
        const double ys[6] = {rxx, ryy, rzz, rxy, ryz, rxz};
        for (int c = 0; c < 6; ++c)
        {
            const double want = ys[c], got = xs[c] * r.scale;
            worst = std::max(worst, std::fabs(got - want) / s.A);
        }
    }
    std::printf("        макс. отн. расхождение компонент (на A): %.3e, J2: %.3e, пластических: %d\n", worst, worst_alg,
                n_plastic);
    check(n_plastic > 300 && n_plastic < 1000, "U10 есть и упругие, и пластические случаи");
    check(worst <= 1e-14, "U10 JC(B=C=0) = prochnost_mises, <= 1e-14");
    check(worst_alg <= 1e-14, "U10 j2_of_deviator = sigma_i, <= 1e-14");
    // Алгебра остальных геометрий против ручных формул
    const double a = 1.3, b = -0.7, c = 0.4, d = 0.9, e = -0.2;
    check(same(jc::j2_of_deviator(jc::Geom::D1, a, b, c, d, e), 0.75 * a * a), "U10 D1 алгебра");
    check(same(jc::j2_of_deviator(jc::Geom::D2Plane, a, b, c, d, e), a * a + c * c), "U10 D2Plane алгебра");
    check(same(jc::j2_of_deviator(jc::Geom::D2Axi, a, b, c, d, e), a * a + a * b + b * b + c * c), "U10 D2Axi алгебра");
    check(same(jc::j2_of_deviator(jc::Geom::D3, a, b, c, d, e), a * a + a * b + b * b + c * c + e * e + d * d),
          "U10 D3 алгебра");
}

void test_fracture()
{
    std::printf("U11-U13: деформация разрушения\n");
    const jc::JcStrength s = steel4340();
    const jc::JcDamage d = damage4340();
    check(near_rel(jc::jc_fracture_strain(d, s, 0.0, 1.0, 300.0), 3.49, 1e-12), "U11 sigma*=0: 3.49");
    const double e3 = jc::jc_fracture_strain(d, s, 1.0 / 3.0, 1.0, 300.0);
    check(near_rel(e3, 0.05 + 3.44 * std::exp(-2.12 / 3.0), 1e-12), "U11 sigma*=1/3, формула");
    check(near_rel(e3, 1.746903, 1e-6), "U11 sigma*=1/3 = 1.746903");
    check(near_rel(jc::jc_fracture_strain(d, s, 0.0, 1.0, 1046.5), 3.49 * (1.0 + 0.61 * 0.5), 1e-12), "U11 T* = 0.5");
    const jc::JcStrength sa = al5083();
    const jc::JcDamage da = damage5083();
    const double eb = jc::jc_fracture_strain(da, sa, 1.0 / 3.0, 100.0, 293.0);
    check(near_rel(eb, (0.178 + 0.389 * std::exp(-0.75)) * std::pow(101.0, 0.147), 1e-12), "U12 Borvik 5083");
    jc::JcDamage dj = da;
    dj.form = jc::DamageRateForm::JC;
    check(near_rel(jc::jc_fracture_strain(dj, sa, 1.0 / 3.0, 100.0, 293.0),
                   (0.178 + 0.389 * std::exp(-0.75)) * (1.0 + 0.147 * std::log(100.0)), 1e-12),
          "U12 JC-форма скорости 5083");
    jc::JcDamage dm = d;
    dm.D1 = -1.0; dm.D2 = 0.0;
    check(same(jc::jc_fracture_strain(dm, s, 5.0, 1.0, 300.0), 1e-6), "U13 ограничение eps_f_min");
    jc::JcDamage dn = d;
    dn.D2 = 1.01e20;
    check(!jc::jc_damage_active(dn) && jc::jc_damage_active(d), "Р6: D >= 1e20 отключает разрушение");
    jc::JcDamage dt = d;
    dt.Tm_f = 0.0;
    check(near_rel(jc::jc_fracture_strain(dt, s, 0.0, 1.0, 1046.5), 3.49 * 1.305, 1e-12), "Tm_f <= 0 -> Tm прочности");
}

void test_strain_rate()
{
    std::printf("U14: эквивалентная скорость деформации\n");
    const double gd = 2.0e3, ed = 1.0e3;
    const jc::Geom all[4] = {jc::Geom::D1, jc::Geom::D2Plane, jc::Geom::D2Axi, jc::Geom::D3};
    const char* nm[4] = {"D1", "D2Plane", "D2Axi", "D3"};
    check(near_rel(jc::equivalent_strain_rate(jc::Geom::D1, ed, 0, 0, 0, 0, 0), 2.0 / 3.0 * ed, 1e-14),
          "D1: одноосная деформация, 2/3 |exx|");
    for (int i = 1; i < 4; ++i)
    {
        check(near_rel(jc::equivalent_strain_rate(all[i], 0, 0, 0, 0.5 * gd, 0, 0), gd / std::sqrt(3.0), 1e-14),
              std::string(nm[i]) + ": чистый сдвиг exy = gd/2 -> gd/sqrt(3)");
    }
    check(near_rel(jc::equivalent_strain_rate(jc::Geom::D3, 0, 0, 0, 0, 0.5 * gd, 0), gd / std::sqrt(3.0), 1e-14) &&
              near_rel(jc::equivalent_strain_rate(jc::Geom::D3, 0, 0, 0, 0, 0, 0.5 * gd), gd / std::sqrt(3.0), 1e-14),
          "D3: сдвиги eyz, exz");
    check(near_rel(jc::equivalent_strain_rate(jc::Geom::D2Axi, ed, -0.5 * ed, -0.5 * ed, 0, 0, 0), ed, 1e-14),
          "D2Axi: одноосное без изменения объёма (eyy = ezz = -exx/2) -> exx");
    check(near_rel(jc::equivalent_strain_rate(jc::Geom::D3, ed, -0.5 * ed, -0.5 * ed, 0, 0, 0), ed, 1e-14),
          "D3: одноосное без изменения объёма -> exx");
    check(near_rel(jc::equivalent_strain_rate(jc::Geom::D2Plane, ed, -ed, 0, 0, 0, 0), 2.0 * ed / std::sqrt(3.0), 1e-14),
          "D2Plane: exx = -eyy = e -> 2e/sqrt(3)");
}

} // namespace

int main()
{
    test_yield();
    test_return();
    test_mises();
    test_fracture();
    test_strain_rate();
    std::printf("%s (%d failures)\n", failures == 0 ? "ALL OK" : "FAILED", failures);
    return failures == 0 ? 0 : 1;
}
