#pragma once
// Искусственное напряжение Монагана против неустойчивости SPH при растяжении
// (tensile instability).
//
// Источники:
//   - J.J. Monaghan, SPH without a tensile instability, J. Comput. Phys. 159 (2000) 290-311;
//   - J.P. Gray, J.J. Monaghan, R.P. Swift, SPH elastic dynamics,
//     Comput. Methods Appl. Mech. Eng. 190 (2001) 6641-6662;
//   - сверено с реализацией PySPH (pysph.sph.solid_mech.basic:
//     MonaghanArtificialStress, MomentumEquationWithStress): eps = 0.3, n = 4.
//
// Схема (GMS 2001, §4):
//   1. Полный тензор напряжений частицы sigma = s - p*I приводится к главным осям.
//   2. Для каждого главного напряжения sigma'_k > 0 (растяжение)
//      R'_k = -eps * sigma'_k, иначе R'_k = 0.
//   3. R поворачивается обратно в исходные оси.
//   4. В парную силу добавляется f_ij^n * (R_i + R_j) (в той же форме, что и
//      напряжения: /(rho_i rho_j) при pa_sph = 1 или R/rho^2 при pa_sph = 2),
//      f_ij = W(r_ij) / W(dp_ij), dp -- начальный шаг частиц.
//
// Отличие от PySPH: здесь R хранится в единицах напряжения (без 1/rho^2) --
// деление на плотность делает ядро в той же форме, что и для s и p. Работа
// искусственной силы записывается в du парно (как работа давления), поэтому
// полная энергия сохраняется; в PySPH энергия этот член не видит.
//
// В осесимметрии (dim = 2, axis_sim) R строится только в плоскости (r, z)
// по (sxx - p, syy - p, sxy); окружная компонента не корректируется
// (осевой член решателя использует только s).
//
// Чистые функции, без ParticleSystem: проверяются в tests/test_art_stress.cpp
// для любой размерности в одной сборке.

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace art_stress
{

// Симметричный тензор 3x3 в нотации xx, yy, zz, xy, yz, xz.
struct Sym3
{
    double xx = 0.0, yy = 0.0, zz = 0.0, xy = 0.0, yz = 0.0, xz = 0.0;
};

// Добавка для одного главного напряжения.
inline double principal_r(double eps, double s) { return s > 0.0 ? -eps * s : 0.0; }

// 1D: тензор -- одно число.
inline double r_1d(double eps, double sxx, double p) { return principal_r(eps, sxx - p); }

// 2D (плоскость x-y): аналитическое разложение 2x2.
// На входе полный тензор (a = sigma_xx, b = sigma_yy, c = sigma_xy).
inline void r_2d(double eps, double a, double b, double c, double& rxx, double& ryy, double& rxy)
{
    const double mean = 0.5 * (a + b);
    const double half = 0.5 * (a - b);
    const double rad = std::sqrt(half * half + c * c);
    const double l1 = mean + rad, l2 = mean - rad; // l1 >= l2
    const double r1 = principal_r(eps, l1), r2 = principal_r(eps, l2);
    if (rad == 0.0)
    {
        // Изотропный тензор: любые оси главные.
        rxx = r1; ryy = r1; rxy = 0.0;
        return;
    }
    // Собственный вектор при l1: угол theta, tan(2 theta) = 2c / (a - b).
    const double cos2 = half / rad, sin2 = c / rad;
    // R = r1 e1 e1^T + r2 e2 e2^T; e1 = (cos t, sin t), e2 = (-sin t, cos t).
    const double rm = 0.5 * (r1 + r2), rh = 0.5 * (r1 - r2);
    rxx = rm + rh * cos2;
    ryy = rm - rh * cos2;
    rxy = rh * sin2;
}

// 3D: разложение методом Якоби (циклический, до сходимости). Матрица 3x3
// симметричная, сходится за несколько проходов; результат R = V diag(r) V^T.
inline Sym3 r_3d(double eps, const Sym3& sig)
{
    double a[3][3] = {{sig.xx, sig.xy, sig.xz}, {sig.xy, sig.yy, sig.yz}, {sig.xz, sig.yz, sig.zz}};
    double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int sweep = 0; sweep < 50; ++sweep)
    {
        const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        const double diag = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2];
        if (off <= 1e-30 * diag || off == 0.0) break;
        for (int p = 0; p < 2; ++p)
        {
            for (int q = p + 1; q < 3; ++q)
            {
                if (a[p][q] == 0.0) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k)
                {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    const double r[3] = {principal_r(eps, a[0][0]), principal_r(eps, a[1][1]), principal_r(eps, a[2][2])};
    auto comp = [&](int i, int j) {
        return r[0] * v[i][0] * v[j][0] + r[1] * v[i][1] * v[j][1] + r[2] * v[i][2] * v[j][2];
    };
    Sym3 out;
    out.xx = comp(0, 0); out.yy = comp(1, 1); out.zz = comp(2, 2);
    out.xy = comp(0, 1); out.yz = comp(1, 2); out.xz = comp(0, 2);
    return out;
}

// Множитель пары f^n, f = W(r) / W(dp). При W(dp) <= 0 -- 0 (выключено).
inline double pair_factor(double w_r, double w_dp, double n)
{
    if (!(w_dp > 0.0)) return 0.0;
    return std::pow(w_r / w_dp, n);
}

} // namespace art_stress
