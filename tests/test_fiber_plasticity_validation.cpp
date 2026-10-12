// tests/test_fiber_plasticity_validation.cpp
// Этап 5 дизайна пластичности нити (fiber/fiber-plasticity-design.md, раздел 10,
// п. 5): валидация на цепочке бусин.
//   Т12  двухволновая картина в стержне: заданная скорость конца v0 выше предела
//        упругости, билинейный закон. Упругий предвестник c_e = sqrt(E/ρ) с
//        амплитудой σ_y0, пластический фронт c_p = sqrt(E_t/ρ), E_t = E·H/(E+H)
//        (теория Кармана-Рахматулина: касательный модуль по ПОЛНОЙ деформации; у
//        закона σ = σ_y0 + H·ε_p модуль H -- по пластической, поэтому в дизайне
//        «c_p = sqrt(H/ρ)» было неточно, см. раздел 8), деформация за фронтом из
//        v0 = c_e·ε_y + c_p·(ε - ε_y), профиль ε_p за фронтом.
//   Т13  независимость от дискретизации: отрезок нити в 10 и 40 бусин,
//        квазистатическое растяжение заданным движением конца; кривая σ-ε
//        одинакова (отдельно во всём диапазоне ε >= 0.01 и в пластической части
//        ε >= 0.03), нет артефактов дискретизации в ε_p; дефект σ_y одной связи
//        на 1 % при H > 0 не локализует ε_p (разность ε_p дефектной и прочих
//        связей = 0.01·σ_y/H); то же с RATE Cowper-Symonds.
//
// Динамика -- локальный драйвер PredictorCorrector (формулы src/time_integration.cpp,
// как в test_fiber_plasticity_integration.cpp, но здесь самодостаточный): порядок
// шага как в боевом цикле (предиктор -> RHS Fiber-Bond -> advance -> проход
// Postprocessing apply_fiber_bond_constitutive). Заданный конец: ускорение
// обнуляется, скорость задаётся, реакция = m·a_связи.
// Собирается с SPH_DIM=3 (CMakeLists.txt).
#include "fiber/fiber_plasticity.hpp"
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/postprocessing.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_plasticity_validation требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

void info(const std::string& what) { std::printf("  info : %s\n", what.c_str()); }

std::string fmt(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.6g", v);
    return buf;
}

struct Chain
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

Chain make_chain(std::size_t N, double l0, const fiber::FiberMaterial& mat)
{
    const double A = interactions::fiber_cross_section_area(mat);
    Chain c;
    c.materials.push_back(mat);
    c.state.resize(N);
    c.topology.resize(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d)
        {
            c.state.x[d][i] = 0.0;
            c.state.vx[d][i] = 0.0;
        }
        c.state.x[0][i] = static_cast<double>(i) * l0;
        c.state.mass[i] = mat.rho * A * l0;
        c.state.object_id[i] = 0;
        c.state.material_id[i] = 0;
        if (i + 1 < N)
        {
            c.topology.next[i] = static_cast<long long>(i + 1);
            c.topology.prev[i + 1] = static_cast<long long>(i);
            c.topology.bonds.push_back(fiber::FiberBond{i, i + 1, l0, l0 * 1.0e9});
        }
    }
    return c;
}

// Драйвер шага; заданные бусины (driven) движутся со скоростью vdrive[i] вдоль x.
struct Driver
{
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
    std::vector<double> dv_prev[dim], v_min[dim];
    std::vector<char> driven;
    std::vector<double> vdrive;
    std::vector<double> react_x; // x-компонента силы связей на заданную бусину, последний шаг [Н]
    bool first = true;

    void ensure(std::size_t n)
    {
        if (dv_prev[0].size() == n) return;
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_prev[d].assign(n, 0.0);
            v_min[d].assign(n, 0.0);
        }
        if (driven.size() != n) driven.assign(n, 0);
        if (vdrive.size() != n) vdrive.assign(n, 0.0);
        react_x.assign(n, 0.0);
    }

    void step(Chain& c, double dt)
    {
        fiber::FiberState& fb = c.state;
        const std::size_t n = fb.size();
        ensure(n);
        if (first) interactions::init_fiber_bond_state(fb, c.topology, c.materials);
        for (std::size_t i = 0; i < n; ++i)
            if (driven[i])
            {
                fb.vx[0][i] = vdrive[i];
                fb.vx[1][i] = 0.0;
                fb.vx[2][i] = 0.0;
            }
        if (!first)
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < n; ++i)
                {
                    v_min[d][i] = fb.vx[d][i];
                    fb.vx[d][i] += 0.5 * dt * dv_prev[d][i];
                }
        const std::vector<fiber::FiberDelta> delta =
            interactions::compute_fiber_bond(exec, fb, c.topology, c.materials, buffers);
        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t d = 0; d < dim; ++d) dv_prev[d][i] = delta[i].dv[d];
            if (driven[i])
            {
                react_x[i] = fb.mass[i] * delta[i].dv[0];
                for (std::size_t d = 0; d < dim; ++d) dv_prev[d][i] = 0.0;
            }
        }
        const double coeff = first ? 0.5 * dt : dt;
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < n; ++i)
            {
                const double base = first ? fb.vx[d][i] : v_min[d][i];
                fb.vx[d][i] = base + coeff * dv_prev[d][i];
                fb.x[d][i] += dt * fb.vx[d][i];
            }
        first = false;
        interactions::apply_fiber_bond_constitutive(exec, fb, c.topology, c.materials, dt);
    }
};

double bond_strain(const Chain& c, std::size_t k)
{
    const fiber::FiberBond& b = c.topology.bonds[k];
    return interactions::fiber_bond_length(c.state, b) / b.l0 - 1.0;
}

// ----------------------------------------------------------------------------
// Т12
// ----------------------------------------------------------------------------
// Положение (в долях длины связи, дробный индекс) внешней границы фронта: первый
// от дальнего конца индекс k, где f[k] >= thr; линейная интерполяция между k и k+1.
double front_position(const std::vector<double>& f, double thr, std::size_t from_far)
{
    for (std::size_t k = from_far; k-- > 1;)
        if (f[k] >= thr)
        {
            const double a = f[k], b = f[k + 1];
            const double frac = (a != b) ? (a - thr) / (a - b) : 0.0; // f[k] >= thr > f[k+1]
            return static_cast<double>(k) + std::clamp(frac, 0.0, 1.0);
        }
    return 0.0;
}

double mean_range(const std::vector<double>& f, std::size_t k0, std::size_t k1)
{
    double s = 0.0;
    for (std::size_t k = k0; k < k1; ++k) s += f[k];
    return s / static_cast<double>(k1 - k0);
}

struct WaveCase
{
    const char* name;
    double E, rho, d, l0, SY, H, v0;
    std::size_t N;
    double tol_ce, tol_cp, tol_eps, tol_epsp, tol_sig, tol_unif; // допуски: ~2-3x над наблюдаемым для случая
    double cs_C = 0.0; // > 0 -- добавить RATE Cowper-Symonds (C, p = 5): аналитики нет, только качественные проверки
};

void run_wave_case(const WaveCase& w)
{
    std::printf("T12 [%s]: волна в стержне (%zu бусин), E=%.3g Па, rho=%.0f, SY=%.3g Па, H=%.3g Па, v0=%.1f м/с:\n", w.name,
                w.N, w.E, w.rho, w.SY, w.H, w.v0);
    fiber::FiberMaterial mat;
    mat.E = w.E;
    mat.nu = 0.3;
    mat.d = w.d;
    mat.rho = w.rho;
    mat.eps_break = 1.0e9;
    mat.plastic.type = fiber::HardeningType::Linear;
    mat.plastic.sy = w.SY;
    mat.plastic.H = w.H;
    if (w.cs_C > 0.0)
    {
        mat.rate.model = fiber::RateModel::CowperSymonds;
        mat.rate.a = w.cs_C;
        mat.rate.b = 5.0;
    }
    Chain c = make_chain(w.N, w.l0, mat);

    const double ce = std::sqrt(w.E / w.rho);
    const double Et = w.E * w.H / (w.E + w.H);
    const double cp = std::sqrt(Et / w.rho);
    const double eps_y = w.SY / w.E;
    const double eps1 = eps_y + (w.v0 - ce * eps_y) / cp; // деформация за пластическим фронтом
    const double epsp1 = (w.E * eps1 - w.SY) / (w.E + w.H);
    const double sig1 = w.SY + w.H * epsp1;
    check(w.v0 > ce * eps_y, "самопроверка: v0 = " + fmt(w.v0) + " > c_e·ε_y = " + fmt(ce * eps_y) + " (выше предела упругости)");
    info("аналитика: c_e = " + fmt(ce) + ", c_p = sqrt(E_t/rho) = " + fmt(cp) + " (sqrt(H/rho) = " + fmt(std::sqrt(w.H / w.rho)) +
         "), eps_y = " + fmt(eps_y) + ", eps1 = " + fmt(eps1) + ", eps_p1 = " + fmt(epsp1) + ", sigma1 = " + fmt(sig1));

    Driver drv;
    drv.ensure(w.N);
    drv.driven[0] = 1;
    drv.vdrive[0] = -w.v0; // конец уходит наружу -> связь (0,1) растягивается
    rhs::CpuExecutor exec_probe(2);
    const double dt = 0.4 * interactions::compute_fiber_bond_dt_limit(exec_probe, c.state, c.topology, c.materials);

    // Два момента: предвестник у связи ~0.25·N и ~0.6·N.
    const double t1 = 0.25 * static_cast<double>(w.N) * w.l0 / ce;
    const double t2 = 0.60 * static_cast<double>(w.N) * w.l0 / ce;
    struct Snap
    {
        double t = 0.0, xe = 0.0, xp = 0.0, sig_pre = 0.0, eps_pl = 0.0, epsp_pl = 0.0, sig_pl = 0.0;
        bool ok = false;
    } s[2];
    double t = 0.0;
    std::size_t taken = 0;
    const std::size_t nb = c.topology.bonds.size();
    while (taken < 2)
    {
        drv.step(c, dt);
        t += dt;
        const double tt = taken == 0 ? t1 : t2;
        if (t < tt) continue;
        std::vector<double> eps(nb), epsp(nb), sig(nb);
        for (std::size_t k = 0; k < nb; ++k)
        {
            eps[k] = bond_strain(c, k);
            epsp[k] = c.topology.bonds[k].st.eps_p;
            sig[k] = w.E * (eps[k] - epsp[k]);
        }
        Snap& sn = s[taken];
        sn.t = t;
        // Позиции: предвестник -- полувысота σ_y/2 (по σ), пластический фронт -- полувысота ε в середине (ε_y, ε1).
        sn.xe = front_position(sig, 0.5 * w.SY, nb - 2);
        const double thr_p = 0.5 * (eps_y + eps1);
        sn.xp = front_position(eps, thr_p, nb - 2);
        // Окна усреднения: за пластическим фронтом (от 30 связей от заданного конца до xp-40),
        // между фронтами (xp+40 .. xe-40).
        const std::size_t kp0 = 30, kp1 = static_cast<std::size_t>(std::max(sn.xp - 40.0, 31.0));
        const std::size_t kq0 = static_cast<std::size_t>(sn.xp + 40.0), kq1 = static_cast<std::size_t>(std::max(sn.xe - 40.0, sn.xp + 41.0));
        sn.eps_pl = mean_range(eps, kp0, kp1);
        sn.epsp_pl = mean_range(epsp, kp0, kp1);
        sn.sig_pl = mean_range(sig, kp0, kp1);
        sn.sig_pre = (kq1 > kq0) ? mean_range(sig, kq0, kq1) : 0.0;
        sn.ok = true;
        if (taken == 1)
        {
            // профиль ε_p за фронтом: однородность на окне (разброс/среднее)
            double mn = 1e300, mx = -1e300;
            for (std::size_t k = kp0; k < kp1; ++k)
            {
                mn = std::min(mn, epsp[k]);
                mx = std::max(mx, epsp[k]);
            }
            info("разброс eps_p на плато за фронтом (max-min)/mean = " + fmt((mx - mn) / sn.epsp_pl));
            // Шероховатость (шум от связи к связи): среднее |ε_p[k] - (ε_p[k-1]+ε_p[k+1])/2| / среднее ε_p. Монотонный
            // дрейф профиля (физический, у RATE -- конец нагружался дольше) её не даёт, пилообразный шум -- даёт.
            double rough = 0.0;
            for (std::size_t k = kp0 + 1; k + 1 < kp1; ++k) rough += std::fabs(epsp[k] - 0.5 * (epsp[k - 1] + epsp[k + 1]));
            rough /= static_cast<double>(kp1 - kp0 - 2) * sn.epsp_pl;
            info("шероховатость eps_p (связь к связи) = " + fmt(rough));
            if (w.cs_C > 0.0)
            {
                char b[160];
                info("профиль eps_p c RATE (каждая 40-я связь от 30 до xp-40):");
                std::string line;
                for (std::size_t k = kp0; k < kp1; k += 40)
                {
                    std::snprintf(b, sizeof b, " [%zu: %.4g]", k, epsp[k]);
                    line += b;
                }
                info(line);
                check(rough < 5.0e-6, "eps_p с RATE: нет пилообразного шума между связями (шероховатость < 5e-6)");
            }
            else
                check((mx - mn) / sn.epsp_pl < w.tol_unif,
                      "eps_p за фронтом однородна (разброс окна < " + fmt(100 * w.tol_unif) + " %, дисперсия цепочки)");
            // за дальним концом предвестника (невозмущённая зона) eps_p == 0 строго
            bool clean = true;
            for (std::size_t k = static_cast<std::size_t>(sn.xe + 20.0); k < nb; ++k) clean = clean && epsp[k] == 0.0;
            check(clean, "перед предвестником eps_p == 0 (пластика только после прихода предвестника)");
        }
        ++taken;
    }
    const double dtm = s[1].t - s[0].t;
    const double ce_m = (s[1].xe - s[0].xe) * w.l0 / dtm;
    const double cp_m = (s[1].xp - s[0].xp) * w.l0 / dtm;
    info("измерено: c_e = " + fmt(ce_m) + " (" + fmt(100.0 * (ce_m / ce - 1.0)) + " %), c_p = " + fmt(cp_m) + " (" +
         fmt(100.0 * (cp_m / cp - 1.0)) + " %)");
    info("плато предвестника σ = " + fmt(s[1].sig_pre) + " (σ_y " + fmt(w.SY) + "), за фронтом ε = " + fmt(s[1].eps_pl) +
         " (ожид. " + fmt(eps1) + "), ε_p = " + fmt(s[1].epsp_pl) + " (ожид. " + fmt(epsp1) + "), σ = " + fmt(s[1].sig_pl) +
         " (ожид. " + fmt(sig1) + ")");
    if (w.cs_C > 0.0)
    {
        // Со скоростной зависимостью замкнутой аналитики нет; сверяемся с rate-independent решением:
        // предел текучести растёт -> пластическая деформация за фронтом меньше, а полная -- не больше.
        // Предвестник остаётся упругим и идёт с c_e.
        check(std::fabs(ce_m / ce - 1.0) < w.tol_ce, "RATE: скорость упругого предвестника c_e в пределах " + fmt(100 * w.tol_ce) + " %");
        check(s[1].epsp_pl < epsp1 && s[1].epsp_pl > 0.0, "RATE: eps_p за фронтом ниже rate-independent значения (" + fmt(epsp1) + ") и положительна");
        check(s[1].sig_pl > sig1, "RATE: напряжение за фронтом выше rate-independent (" + fmt(sig1) + ")");
        check(std::isfinite(s[1].epsp_pl) && std::isfinite(s[1].sig_pl), "RATE: профиль конечен");
        return;
    }
    check(std::fabs(ce_m / ce - 1.0) < w.tol_ce, "скорость упругого предвестника c_e в пределах " + fmt(100 * w.tol_ce) + " %");
    check(std::fabs(cp_m / cp - 1.0) < w.tol_cp, "скорость пластического фронта c_p = sqrt(E_t/rho) в пределах " + fmt(100 * w.tol_cp) + " %");
    check(std::fabs(s[1].sig_pre / w.SY - 1.0) < w.tol_sig, "амплитуда предвестника = σ_y0 в пределах " + fmt(100 * w.tol_sig) + " %");
    check(std::fabs(s[1].eps_pl / eps1 - 1.0) < w.tol_eps, "ε за пластическим фронтом (Карман-Рахматулин) в пределах " + fmt(100 * w.tol_eps) + " %");
    check(std::fabs(s[1].epsp_pl / epsp1 - 1.0) < w.tol_epsp, "ε_p за фронтом в пределах " + fmt(100 * w.tol_epsp) + " %");
}

void test_t12()
{
    // E=2 ГПа, rho=1200: c_e = 1291 м/с. SY=40 МПа (ε_y=2 %), H=0.2E: E_t=E/6... c_p=527 м/с.
    run_wave_case({"H=0.2E", 2.0e9, 1200.0, 1.0e-3, 1.0e-3, 40.0e6, 0.4e9, 40.0, 800, 0.006, 0.01, 0.0015, 0.006, 0.003, 0.02});
    // Более жёсткое упрочнение: H=E, E_t=E/2, c_p = 913 м/с (фронты ближе).
    run_wave_case({"H=E", 2.0e9, 1200.0, 1.0e-3, 1.0e-3, 40.0e6, 2.0e9, 60.0, 800, 0.006, 0.008, 0.003, 0.018, 0.013, 0.09});
    // Та же постановка (H=0.2E) со скоростной зависимостью: ε̇_p во фронте ~ 4e3 1/с, C = 1e4 -> R ~ 1.8. Риск 3 дизайна (шум ε̇).
    run_wave_case({"H=0.2E, CS C=1e4", 2.0e9, 1200.0, 1.0e-3, 1.0e-3, 40.0e6, 0.4e9, 40.0, 800, 0.003, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e4});
}

// ----------------------------------------------------------------------------
// Т13
// ----------------------------------------------------------------------------
struct Curve
{
    std::vector<double> eps, sig;
    std::vector<double> epsp_end; // eps_p по связям в конце
    double eps_end = 0.0;
};

// defect_k < nbeads-1 -- у связи с этим номером предел текучести σ_y уменьшен на 1 % (отдельный материал 1;
// материал связи берётся по material_id первой бусины). Остальные связи -- материал 0.
constexpr std::size_t kNoDefect = static_cast<std::size_t>(-1);
constexpr double kDefectFactor = 0.99;

Curve run_stretch(std::size_t nbeads, bool with_rate, double eta_kelvin, std::size_t defect_k = kNoDefect)
{
    const double L = 0.04, E = 2.0e9, rho = 1200.0, d = 1.0e-3;
    const double edot = 100.0;      // 1/с -- номинальная скорость деформации
    const double eps_final = 0.07;
    const double t_ramp = 2.0e-4;   // плавный разгон конца (6 пробегов волны по отрезку)
    fiber::FiberMaterial mat;
    mat.E = E;
    mat.nu = 0.3;
    mat.d = d;
    mat.rho = rho;
    mat.eps_break = 1.0e9;
    mat.plastic.type = fiber::HardeningType::Linear;
    mat.plastic.sy = 40.0e6;
    mat.plastic.H = 0.4e9;
    if (with_rate)
    {
        mat.rate.model = fiber::RateModel::CowperSymonds;
        mat.rate.a = 100.0; // C [1/с]
        mat.rate.b = 5.0;   // p
    }
    mat.visco.eta = eta_kelvin;
    const double l0 = L / static_cast<double>(nbeads - 1);
    Chain c = make_chain(nbeads, l0, mat);
    if (defect_k != kNoDefect)
    {
        fiber::FiberMaterial md = mat;
        md.plastic.sy = kDefectFactor * mat.plastic.sy;
        c.materials.push_back(md);
        c.state.material_id[defect_k] = 1;
    }
    const double A = interactions::fiber_cross_section_area(mat);

    Driver drv;
    drv.ensure(nbeads);
    drv.driven[0] = 1;
    drv.driven[nbeads - 1] = 1;
    drv.vdrive[0] = 0.0;
    rhs::CpuExecutor exec_probe(2);
    const double dt = 0.8 * interactions::compute_fiber_bond_dt_limit(exec_probe, c.state, c.topology, c.materials);

    // v(t) = vmax·s(t/t_ramp) на разгоне (s -- сглаженная ступень), далее постоянна; vmax = edot·L.
    const double vmax = edot * L;
    // путь за разгон: ∫v = vmax·t_ramp/2 (для s = (1-cos)/2); дальше линейно
    Curve out;
    double t = 0.0, u = 0.0; // u -- смещение конца
    const double t_end = t_ramp / 2.0 + eps_final * L / vmax;
    const double pi_ = 3.14159265358979323846;
    while (u < eps_final * L)
    {
        const double sfun = t < t_ramp ? 0.5 * (1.0 - std::cos(pi_ * t / t_ramp)) : 1.0;
        drv.vdrive[nbeads - 1] = vmax * sfun;
        drv.step(c, dt);
        t += dt;
        u += drv.vdrive[nbeads - 1] * dt;
        const double eps = (c.state.x[0][nbeads - 1] - c.state.x[0][0] - L) / L;
        out.eps.push_back(eps);
        out.sig.push_back(std::fabs(drv.react_x[0]) / A);
        if (t > 2.0 * t_end) break;
    }
    out.eps_end = out.eps.back();
    for (const fiber::FiberBond& b : c.topology.bonds) out.epsp_end.push_back(b.st.eps_p);
    return out;
}

double interp(const Curve& c, double e)
{
    for (std::size_t i = 1; i < c.eps.size(); ++i)
        if (c.eps[i] >= e && c.eps[i - 1] < e)
        {
            const double f = (e - c.eps[i - 1]) / (c.eps[i] - c.eps[i - 1]);
            return c.sig[i - 1] + f * (c.sig[i] - c.sig[i - 1]);
        }
    return std::nan("");
}

void compare_curves(const char* tag, const Curve& a, const Curve& b, double e_from, double tol_rel)
{
    // Сравнение в точках от e_from: 0.01 -- за пределом разгона (σ-ε уже квазистатична), 0.03 -- пластическая часть.
    double worst = 0.0;
    std::string line;
    for (double e = e_from; e <= 0.0651; e += 0.005)
    {
        const double sa = interp(a, e), sb = interp(b, e);
        const double rel = std::fabs(sa - sb) / std::max(sa, sb);
        worst = std::max(worst, rel);
        char buf[96];
        std::snprintf(buf, sizeof buf, " [%.3f: %.4g/%.4g]", e, sa / 1e6, sb / 1e6);
        line += buf;
    }
    info(std::string(tag) + " σ(ε), МПа (10/40 бусин), ε >= " + fmt(e_from) + ":" + line);
    info(std::string(tag) + " максимальное относительное расхождение кривых при ε >= " + fmt(e_from) + " = " + fmt(worst));
    check(worst < tol_rel, std::string(tag) + ": кривая σ-ε 10 и 40 бусин совпадает при ε >= " + fmt(e_from) + " (< " + fmt(100 * tol_rel) + " %)");
}

void check_no_discretization_artifacts(const char* tag, const Curve& c)
{
    double mn = 1e300, mx = -1e300, mean = 0.0;
    for (double e : c.epsp_end)
    {
        mn = std::min(mn, e);
        mx = std::max(mx, e);
        mean += e;
    }
    mean /= static_cast<double>(c.epsp_end.size());
    info(std::string(tag) + ": eps_p по связям: min " + fmt(mn) + ", max " + fmt(mx) + ", среднее " + fmt(mean) + ", max/mean " + fmt(mx / mean));
    check(mx / mean < 1.01 && mn / mean > 0.99, std::string(tag) + ": нет артефактов дискретизации -- eps_p одинакова по связям (±1 %)");
    // Среднее eps_p согласовано с полным ε: суммарное удлинение = упругое + пластическое.
    check(mean > 0.0, std::string(tag) + ": пластическая деформация развилась");
}

// Дефект σ_y одной связи на 1 %: при H > 0 связь упрочняется и «догоняется», локализации нет. Квазистатика:
// у всех связей одно σ = σ_y,i + H·ε_p,i => ε_p,деф - ε_p = (σ_y - σ_y,деф)/H = 0.01·σ_y/H (с RATE R(ε̇_p) разный у связей --
// аналитика приближённая: наблюдается 0.87 от неё, допуск 35 %). Допуск -- ~2-3x над наблюдаемым.
void check_defect_no_localization(const char* tag, const Curve& c, std::size_t defect_k, double tol_rel_diff, double tol_rel_others)
{
    const double SY = 40.0e6, H = 0.4e9;
    const double expected = (1.0 - kDefectFactor) * SY / H;
    double mean = 0.0, mn = 1e300, mx = -1e300;
    std::size_t n = 0;
    for (std::size_t k = 0; k < c.epsp_end.size(); ++k)
    {
        if (k == defect_k) continue;
        mean += c.epsp_end[k];
        mn = std::min(mn, c.epsp_end[k]);
        mx = std::max(mx, c.epsp_end[k]);
        ++n;
    }
    mean /= static_cast<double>(n);
    const double diff = c.epsp_end[defect_k] - mean;
    info(std::string(tag) + ": eps_p дефектной связи " + fmt(c.epsp_end[defect_k]) + ", прочих (среднее) " + fmt(mean) +
         ", разность " + fmt(diff) + ", аналитика 0.01·σ_y/H = " + fmt(expected) + " (отн. отклонение " + fmt(diff / expected - 1.0) +
         "); разброс прочих (max-min)/mean = " + fmt((mx - mn) / mean));
    check(c.epsp_end[defect_k] > mean && c.epsp_end[defect_k] < 1.1 * mean,
          std::string(tag) + ": eps_p дефектной связи выше средней, но не катастрофически (< +10 %), локализации нет");
    check(std::fabs(diff / expected - 1.0) < tol_rel_diff,
          std::string(tag) + ": разность eps_p(дефект) - eps_p(прочие) = 0.01·σ_y/H в пределах " + fmt(100 * tol_rel_diff) + " %");
    check((mx - mn) / mean < tol_rel_others,
          std::string(tag) + ": прочие связи остаются однородными (разброс < " + fmt(100 * tol_rel_others) + " %)");
}

void test_t13()
{
    std::printf("T13: независимость от дискретизации (отрезок 40 мм, 10 и 40 бусин, растяжение до 7 %%, ε̇ = 100 1/с):\n");
    const double eta = 1.0e4; // Па·с -- Кельвин гасит колебания бусин; σ_η = η·ε̇ = 1 МПа (2.5 % от σ_y), тот же для обеих дискретизаций
    const auto w0 = std::chrono::steady_clock::now();
    const Curve a = run_stretch(10, false, eta);
    const Curve b = run_stretch(40, false, eta);
    compare_curves("без RATE", a, b, 0.01, 0.008);
    compare_curves("без RATE", a, b, 0.03, 0.0008);
    check_no_discretization_artifacts("10 бусин", a);
    check_no_discretization_artifacts("40 бусин", b);
    {
        // Сверка с аналитикой статического билинейного закона (эффект η·ε̇ добавляет ~η·ε̇ к σ).
        const double E = 2.0e9, SY = 40.0e6, H = 0.4e9, eta_s = 1.0e4 * 100.0;
        const double e = 0.06;
        const double sig_an = SY + E * H / (E + H) * (e - SY / E) + eta_s;
        const double sa = interp(a, e), sb = interp(b, e);
        info("σ(0.06): 10 бусин " + fmt(sa) + ", 40 бусин " + fmt(sb) + ", аналитика билинейная + η·ε̇ = " + fmt(sig_an));
        check(std::fabs(sa / sig_an - 1.0) < 0.002 && std::fabs(sb / sig_an - 1.0) < 0.002,
              "σ(0.06) обеих дискретизаций совпадает с билинейной аналитикой + η·ε̇ (0.2 %)");
    }
    const Curve ar = run_stretch(10, true, eta);
    const Curve br = run_stretch(40, true, eta);
    compare_curves("с RATE CS", ar, br, 0.01, 0.008);
    compare_curves("с RATE CS", ar, br, 0.03, 0.003);
    check_no_discretization_artifacts("10 бусин (CS)", ar);
    check_no_discretization_artifacts("40 бусин (CS)", br);
    {
        const double e = 0.06;
        check(interp(ar, e) > 1.2 * interp(a, e), "RATE CS повышает напряжение (σ с CS > 1.2 σ без)");
    }
    {
        // Дефект σ_y (-1 %) средней связи: 10 бусин -- связь 4, 40 бусин -- связь 19.
        const Curve ad = run_stretch(10, false, eta, 4);
        const Curve bd = run_stretch(40, false, eta, 19);
        check_defect_no_localization("дефект, 10 бусин", ad, 4, 1.0e-6, 1.0e-8);
        check_defect_no_localization("дефект, 40 бусин", bd, 19, 1.0e-6, 1.0e-8);
        const Curve ard = run_stretch(10, true, eta, 4);
        const Curve brd = run_stretch(40, true, eta, 19);
        check_defect_no_localization("дефект, 10 бусин (CS)", ard, 4, 0.35, 1.5e-5);
        check_defect_no_localization("дефект, 40 бусин (CS)", brd, 19, 0.35, 1.5e-5);
    }
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
    info("T13 wall = " + fmt(wall) + " с");
}

} // namespace

int main()
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    const auto w0 = std::chrono::steady_clock::now();
    test_t12();
    test_t13();
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
    std::printf("\n=== test_fiber_plasticity_validation: %s (провалов: %d, wall %.2f с) ===\n",
                failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", failures, wall);
    return failures == 0 ? 0 : 1;
}
