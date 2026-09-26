// tests/test_solid_contact.cpp
// Unit-тесты Solid-Solid-Contact (fiber/CUDA_PLAN.md, этап Г1; fiber/PLAN.md,
// этап Ф4-A; IMPLEMENTATION_REPORT.md §4.5, «У нового модуля Solid-Solid-
// Contact нет собственного модульного теста»).
//
// interactions/solid_contact.hpp переиспользует ТУ ЖЕ математику
// (contact::penalty_force -> pinball_force/swegle_force/impulse_force), что
// tests/test_contact.cpp уже проверяет с допуском 1e-12 -- здесь формулы НЕ
// передоказываются. Проверяется собственный код модуля:
//   - баланс импульса пары (m_i*dv_i + m_j*dv_j == 0, третий закон Ньютона);
//   - знак силы (отталкивание вдоль линии центров);
//   - фильтр body[i]==body[j] -> контакт не регистрируется (Ф4-A, ОТЛИЧИЕ от
//     старого ext_force::is_contact_material_pair);
//   - dt_limit (SolidContactDiagnostics::min_dt) для типов 1-3 совпадает с
//     penalty_force на тех же геометрии/свойствах пары;
//   - contact_force_type == 0 и отсутствие проникания (p<=0) -> нулевая сила.
//
// Стиль -- по образцу tests/test_fiber_bond.cpp/test_fiber_solid.cpp
// (тот же check()/near_rel()/near_abs(), CpuExecutor/ContributionBuffer
// напрямую). dim по умолчанию (1, param.hpp) -- геометрия контакта dim-
// независима (contact_force.hpp), как и в tests/test_contact.cpp.
//
// fp32-дисциплина (fiber/CUDA_PLAN.md, этап Г1, п.4; в отличие от
// tests/test_fiber_*.cpp этот файл -- НОВЫЙ, поэтому пишется fp32-чистым
// сразу, без исключения из -Wdouble-promotion/-Wimplicit-float-conversion):
// все поля ParticleSystem/contact::PairGeom/PairProps -- real_t, граница с
// double -- явные static_cast<double>/real_t(...) при чтении/записи.
#include "contact_force.hpp"
#include "interactions/solid_contact.hpp"
#include "materials.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b)); }

// fiber/CUDA_PLAN.md, этап Г1: fp32-вариант (SPH_PRECISION=32) считает через
// real_t=float (эпсилон ~1.2e-7). Баланс импульса и dt_limit -- отношения
// f/m_i*m_i и деления в penalty_force -- в double точны до ~1e-12
// (представление f/m_i*m_i в общем случае не бит-в-бит f), во float то же
// вычисление накапливает ошибку на уровне эпсилон float, поэтому допуск
// ослаблен на несколько порядков. Знак силы и фильтр body[] -- логические
// проверки, от точности арифметики не зависят и не ослабляются.
#if SPH_PRECISION == 32
constexpr double kTolBalanceRel = 1e-5;  // было 1e-9 (m_i*dv_i + m_j*dv_j == 0)
constexpr double kTolDtRel = 1e-5;       // было 1e-9 (dt_limit диагностики против penalty_force)
#else
constexpr double kTolBalanceRel = 1e-9;
constexpr double kTolDtRel = 1e-9;
#endif

// Пара Solid-частиц: itype 3 (тантал, G=69e9) / 4 (алюминий, G=26e9) -- те же
// материалы, что tests/test_contact.cpp::make_pair, у обоих есть G в реестре
// (materials.hpp) для contact_force_type==1 (pinball, читает g_eff).
// rho задаётся из требуемого радиуса R = 0.5*(mass/rho)^(1/dim): при dim=1
// rho = mass/h даёт ровно R = 0.5*h (тот же приём, что в test_contact.cpp).
struct PairSetup
{
    ParticleSystem sys;
    std::vector<real_t> c;
    double R0 = 0.0, R1 = 0.0;
};

PairSetup make_pair(double gap_fraction, double v0, double v1, int body0, int body1)
{
    PairSetup ps;
    ParticleSystem& sys = ps.sys;
    sys.resize(2);
    sys.inzializ_index(0);
    sys.inzializ_index(1);

    const double m0 = 2e-6, m1 = 3e-6;
    const double h0 = 1e-3, h1 = 1e-3;
    sys.mass[0] = real_t(m0); sys.mass[1] = real_t(m1);
    sys.hsml[0] = real_t(h0); sys.hsml[1] = real_t(h1);
    // rho = mass/h^dim -> R = 0.5*h (см. шапку файла).
    sys.rho[0] = real_t(m0 / std::pow(h0, static_cast<double>(dim)));
    sys.rho[1] = real_t(m1 / std::pow(h1, static_cast<double>(dim)));
    sys.rho_0[0] = sys.rho[0];
    sys.rho_0[1] = sys.rho[1];
    ps.R0 = 0.5 * h0;
    ps.R1 = 0.5 * h1;

    sys.itype[0] = 3; sys.itype[1] = 4;
    sys.body[0] = body0; sys.body[1] = body1;

    // Пара на оси x: i слева (x=0), j справа на расстоянии
    // gap_fraction*(R0+R1) -- < 1 даёт проникание p>0, >= 1 -- p<=0.
    const double Rsum = ps.R0 + ps.R1;
    sys.x[0][0] = real_t(0.0);
    sys.x[0][1] = real_t(gap_fraction * Rsum);
    sys.vx[0][0] = real_t(v0);
    sys.vx[0][1] = real_t(v1);

    ps.c.assign(2, real_t(500.0)); // скорость звука пары -- как в test_contact.cpp make_pair
    return ps;
}

SimulationConfig make_cfg(int contact_type, double kp = 1.0, double gamma = 1.0)
{
    SimulationConfig cfg;
    cfg.apply_defaults();
    cfg.contact_force_type = contact_type;
    cfg.contact_penalty_scale = kp;
    cfg.contact_impulse_gamma = gamma;
    return cfg;
}

// Прогоняет compute_solid_contact_contributions + reduce, возвращает
// готовые SolidDelta (только dv[0] содержателен при dim по умолчанию).
std::vector<interactions::SolidDelta> run_contributions(const SimulationConfig& cfg, const Materials& mats,
                                                         const PairSetup& ps, double dt_applied)
{
    rhs::CpuExecutor exec(1);
    std::vector<std::pair<rhs::index_t, rhs::index_t>> real_pairs{{0, 1}};
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> buffers;
    interactions::compute_solid_contact_contributions(exec, cfg, mats, ps.sys, ps.c, real_t(dt_applied), real_pairs,
                                                       buffers);
    return rhs::reduce_contributions(exec, buffers, ps.sys.size());
}

interactions::SolidContactDiagnostics run_diagnostics(const SimulationConfig& cfg, const Materials& mats,
                                                       const PairSetup& ps, double dt_applied)
{
    rhs::CpuExecutor exec(1);
    std::vector<std::pair<rhs::index_t, rhs::index_t>> real_pairs{{0, 1}};
    return interactions::compute_solid_contact_diagnostics(exec, cfg, mats, ps.sys, ps.c, real_t(dt_applied),
                                                            real_pairs);
}

// -----------------------------------------------------------------------
// T1: баланс импульса пары и знак силы (отталкивание) -- типы 1-3.
// -----------------------------------------------------------------------
void test_balance_and_sign()
{
    std::printf("test_balance_and_sign:\n");
    const Materials mats; // реестр по умолчанию (G для itype 3, 4 -- см. шапку)
    const double dt_applied = 1e-8;

    for (int type : {contact::Pinball, contact::Swegle, contact::Impulse})
    {
        // gap_fraction=0.9 -> p = 0.1*(R0+R1) > 0; i движется вправо (+v),
        // j движется влево (-v) -- пара сближается (pdot>0), сила ненулевая
        // у любого из трёх типов.
        PairSetup ps = make_pair(0.9, /*v0*/ 2.0, /*v1*/ -2.0, /*body0*/ 1, /*body1*/ 2);
        SimulationConfig cfg = make_cfg(type);
        std::vector<interactions::SolidDelta> delta = run_contributions(cfg, mats, ps, dt_applied);
        const double dv0 = static_cast<double>(delta[0].dv[0]);
        const double dv1 = static_cast<double>(delta[1].dv[0]);

        check(dv0 != 0.0 && dv1 != 0.0,
              "тип " + std::to_string(type) + ": контакт даёт ненулевую силу обеим частицам");
        check(dv0 < 0.0, "тип " + std::to_string(type) + ": частица i (слева) отталкивается в -x (от j)");
        check(dv1 > 0.0, "тип " + std::to_string(type) + ": частица j (справа) отталкивается в +x (от i)");

        const double balance = static_cast<double>(ps.sys.mass[0]) * dv0 + static_cast<double>(ps.sys.mass[1]) * dv1;
        const double scale = std::fabs(static_cast<double>(ps.sys.mass[0]) * dv0);
        check(near_rel(balance, 0.0, kTolBalanceRel * std::max(scale, 1.0)),
              "тип " + std::to_string(type) +
                  ": баланс импульса пары m_i*dv_i + m_j*dv_j == 0 (3-й закон Ньютона), balance=" +
                  std::to_string(balance));
    }
}

// -----------------------------------------------------------------------
// T2: фильтр body[i]==body[j] -- контакт не регистрируется даже при
// проникании и сближении (Ф4-A: ОТЛИЧИЕ от itype[i]!=itype[j]).
// -----------------------------------------------------------------------
void test_body_filter()
{
    std::printf("test_body_filter:\n");
    const Materials mats;
    const double dt_applied = 1e-8;

    // Тот же penetrating/сближающийся кейс, что в T1, но body0==body1.
    PairSetup ps = make_pair(0.9, 2.0, -2.0, /*body0*/ 5, /*body1*/ 5);
    SimulationConfig cfg = make_cfg(contact::Impulse);
    std::vector<interactions::SolidDelta> delta = run_contributions(cfg, mats, ps, dt_applied);
    check(delta[0].dv[0] == real_t(0.0) && delta[1].dv[0] == real_t(0.0),
          "body[i]==body[j]: сила отсутствует, даже если геометрически пара проникает и сближается");

    interactions::SolidContactDiagnostics diag = run_diagnostics(cfg, mats, ps, dt_applied);
    check(diag.pair_count == 0, "body[i]==body[j]: диагностика не регистрирует пару (pair_count==0)");
}

// -----------------------------------------------------------------------
// T3: contact_force_type==0 и отсутствие проникания (p<=0) -> нулевая сила.
// -----------------------------------------------------------------------
void test_off_and_no_penetration()
{
    std::printf("test_off_and_no_penetration:\n");
    const Materials mats;
    const double dt_applied = 1e-8;

    // (a) тип Off, но геометрически пара проникает и сближается.
    {
        PairSetup ps = make_pair(0.9, 2.0, -2.0, 1, 2);
        SimulationConfig cfg = make_cfg(contact::Off);
        std::vector<interactions::SolidDelta> delta = run_contributions(cfg, mats, ps, dt_applied);
        check(delta[0].dv[0] == real_t(0.0) && delta[1].dv[0] == real_t(0.0),
              "contact_force_type==0: сила отсутствует, даже если пара проникает и сближается");
        interactions::SolidContactDiagnostics diag = run_diagnostics(cfg, mats, ps, dt_applied);
        check(diag.pair_count == 0, "contact_force_type==0: диагностика не регистрирует пару");
    }

    // (b) тип включён, но p<=0 (частицы не проникают: gap_fraction=1.5).
    for (int type : {contact::Pinball, contact::Swegle, contact::Impulse})
    {
        PairSetup ps = make_pair(1.5, 2.0, -2.0, 1, 2);
        SimulationConfig cfg = make_cfg(type);
        std::vector<interactions::SolidDelta> delta = run_contributions(cfg, mats, ps, dt_applied);
        check(delta[0].dv[0] == real_t(0.0) && delta[1].dv[0] == real_t(0.0),
              "тип " + std::to_string(type) + ", p<=0: сила отсутствует (нет проникания)");
        interactions::SolidContactDiagnostics diag = run_diagnostics(cfg, mats, ps, dt_applied);
        check(diag.pair_count == 0, "тип " + std::to_string(type) + ", p<=0: диагностика не регистрирует пару");
    }
}

// -----------------------------------------------------------------------
// T4: dt_limit (SolidContactDiagnostics::min_dt) для типов 1-3 совпадает с
// penalty_force(), вызванным на ТЕХ ЖЕ геометрии/свойствах пары, которые
// строит compute_solid_contact_diagnostics -- проверяет обвязку модуля, не
// саму формулу (её проверяет test_contact.cpp M4/M5/M8/M13).
// -----------------------------------------------------------------------
void test_dt_limit()
{
    std::printf("test_dt_limit:\n");
    const Materials mats;
    const double dt_applied = 1e-8;
    const SimulationConfig base = make_cfg(contact::Off); // только для cfg.CFL по умолчанию

    for (int type : {contact::Pinball, contact::Swegle, contact::Impulse})
    {
        PairSetup ps = make_pair(0.9, 2.0, -2.0, 1, 2);
        SimulationConfig cfg = make_cfg(type);

        interactions::SolidContactDiagnostics diag = run_diagnostics(cfg, mats, ps, dt_applied);
        check(diag.pair_count == 1, "тип " + std::to_string(type) + ": ровно одна контактная пара");

        // Ручной пересчёт той же геометрии, что строит compute_solid_contact_*
        // (interactions/solid_contact.hpp): r, p, pdot, v_rel, m_red, c_pair,
        // g_eff, r_eff, dt_ref -- дословно те же формулы. Границы double/real_t
        // -- явные (см. шапку файла).
        const double x0 = static_cast<double>(ps.sys.x[0][0]), x1 = static_cast<double>(ps.sys.x[0][1]);
        const double vx0 = static_cast<double>(ps.sys.vx[0][0]), vx1 = static_cast<double>(ps.sys.vx[0][1]);
        const double r = std::fabs(x1 - x0);
        const double p = (ps.R0 + ps.R1) - r;
        const double dot = (x0 - x1) * (vx0 - vx1);
        contact::PairGeom g;
        g.r = real_t(r); g.p = real_t(p); g.pdot = real_t(-dot / r);
        g.v_rel = real_t(std::fabs(vx0 - vx1));

        const double G0 = static_cast<double>(shear_modulus_by_itype(mats, ps.sys.itype[0]));
        const double G1 = static_cast<double>(shear_modulus_by_itype(mats, ps.sys.itype[1]));
        const double m0 = static_cast<double>(ps.sys.mass[0]), m1 = static_cast<double>(ps.sys.mass[1]);
        contact::PairProps pr;
        pr.m_red = real_t(m0 * m1 / (m0 + m1));
        pr.c_pair = std::min(ps.c[0], ps.c[1]);
        pr.g_eff = real_t((G0 > 0 && G1 > 0) ? G0 * G1 / (G0 + G1) : 0.0);
        pr.r_eff = real_t(ps.R0 * ps.R1 / (ps.R0 + ps.R1));
        auto pair_dt_ref = [&](std::size_t k) {
            const double vk = static_cast<double>(ps.sys.vx[0][k]);
            const double ck = static_cast<double>(ps.c[k]);
            const double hk = static_cast<double>(ps.sys.hsml[k]);
            const double denom = ck + std::sqrt(vk * vk);
            return denom > 0 ? base.CFL * hk / denom : 1e30;
        };
        pr.dt_ref = real_t(std::min(pair_dt_ref(0), pair_dt_ref(1)));

        const contact::Result expected = contact::penalty_force(type, g, pr, real_t(cfg.contact_penalty_scale),
                                                                  real_t(dt_applied), real_t(cfg.contact_impulse_gamma));
        check(near_rel(static_cast<double>(diag.min_dt), static_cast<double>(expected.dt_limit), kTolDtRel),
              "тип " + std::to_string(type) + ": diag.min_dt совпадает с penalty_force().dt_limit (" +
                  std::to_string(static_cast<double>(diag.min_dt)) + " vs " +
                  std::to_string(static_cast<double>(expected.dt_limit)) + ")");
        check(near_rel(static_cast<double>(diag.max_delta), p, kTolDtRel),
              "тип " + std::to_string(type) + ": diag.max_delta совпадает с проникновением p");
    }
}

} // namespace

int main()
{
    std::printf("=== interactions/solid_contact.hpp: модульные тесты (fiber/CUDA_PLAN.md, Г1) ===\n");
    test_balance_and_sign();
    test_body_filter();
    test_off_and_no_penetration();
    test_dt_limit();

    std::printf("\n=== %s (провалов: %d) ===\n", failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", failures);
    return failures == 0 ? 0 : 1;
}
