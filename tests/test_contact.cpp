// test_contact.cpp
// Unit-тест контактных сил «частица-частица» (part_0.md, ADR-007).
// Шаг 2: чистая математика (contact_force.hpp) + реестр модуля сдвига
// (materials.hpp). Аналитика типов 1-3 сверяется с формулами (2.3)-(2.10)
// с относительным допуском 1e-12.
// Шаг 4: тесты оркестратора (ext_force) — пары строятся вручную, без поиска
// соседей (niac/pair_i/pair_j задаются напрямую).
#include "contact_force.hpp"
#include "materials.hpp"
#include "ext_force.hpp"
#include <cmath>
#include <cstdio>

namespace
{

int failures = 0;

bool rel_close(double a, double b, double tol = 1e-12)
{
    const double scale = std::max({std::fabs(a), std::fabs(b), 1.0});
    return std::fabs(a - b) <= tol * scale;
}

void check(bool cond, const char* what)
{
    if (!cond)
    {
        std::printf("FAIL: %s\n", what);
        ++failures;
    }
}

void test_shear_modulus_registry()
{
    Materials m;
    // M1: itype 4, 10, 11 -> 26e9, 199.4e9, 1.165e11
    // Для itype 10 (B4C) значение обновлено вместе с materials.hpp: упругие
    // константы приведены к марке Eagle Picher (SAND96-0551), G = 197e9 ->
    // 199.4e9. Проверка здесь — про диспетчеризацию по itype, а не про физику.
    check(rel_close(shear_modulus_by_itype(m, 4), 26e9), "M1 itype 4 -> 26e9");
    check(rel_close(shear_modulus_by_itype(m, 10), 199.4e9), "M1 itype 10 -> 199.4e9");
    check(rel_close(shear_modulus_by_itype(m, 11), 1.165e11), "M1 itype 11 -> 1.165e11");

    // M2: itype 1, 2, 12, 99, -3 -> 0.0 (без выхода за границы массивов)
    check(shear_modulus_by_itype(m, 1) == 0.0, "M2 itype 1 -> 0");
    check(shear_modulus_by_itype(m, 2) == 0.0, "M2 itype 2 -> 0");
    check(shear_modulus_by_itype(m, 12) == 0.0, "M2 itype 12 -> 0");
    check(shear_modulus_by_itype(m, 99) == 0.0, "M2 itype 99 -> 0");
    check(shear_modulus_by_itype(m, -3) == 0.0, "M2 itype -3 -> 0");
}

void test_type1_pinball()
{
    // M3: p>0, pdot=0 -> fn == 0 (F1 нулевая, min даёт 0)
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-5; g.pdot = 0.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7;
        pr.g_eff = 1e10; pr.r_eff = 1e-5;
        const auto res = contact::pinball_force(g, pr, 1.0);
        check(res.fn == 0.0, "M3 pdot=0 -> fn==0");
    }

    // M4: ветка Герца доминирует (f2 < f1) -> fn == kp*f2
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7;
        pr.g_eff = 1e3; pr.r_eff = 1e-4;
        const double kp = 1.0;
        const double f1 = pr.m_red * g.pdot / pr.dt_ref;
        const double f2 = pr.g_eff * std::sqrt(pr.r_eff) * g.p * std::sqrt(g.p);
        check(f2 < f1, "M4 setup: F2 доминирует (меньше F1)");
        const auto res = contact::pinball_force(g, pr, kp);
        check(rel_close(res.fn, kp * f2), "M4 fn == kp*f2 (Герц)");
        const double k_eff = 1.5 * res.fn / g.p;
        const double expected_dt = contact::kOscSafety / std::sqrt(k_eff / pr.m_red);
        check(rel_close(res.dt_limit, expected_dt), "M4 dt_limit по Герцу");
    }

    // M5: ветка F1 доминирует (f1 < f2) -> fn == kp*f1
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-4; g.pdot = 1e-3; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7;
        pr.g_eff = 1e12; pr.r_eff = 1e-3;
        const double kp = 1.0;
        const double f1 = pr.m_red * g.pdot / pr.dt_ref;
        const double f2 = pr.g_eff * std::sqrt(pr.r_eff) * g.p * std::sqrt(g.p);
        check(f1 < f2, "M5 setup: F1 доминирует (меньше F2)");
        const auto res = contact::pinball_force(g, pr, kp);
        check(rel_close(res.fn, kp * f1), "M5 fn == kp*f1 (вязкая)");
        const double expected_dt = contact::kOscSafety * 2.0 * pr.dt_ref / kp;
        check(rel_close(res.dt_limit, expected_dt), "M5 dt_limit вязкая ветка");
    }

    // M6: g_eff == 0 при pdot>0 -> fn == kp*f1 (не 0 — проверка отката)
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7;
        pr.g_eff = 0.0; pr.r_eff = 1e-5;
        const double kp = 1.0;
        const double f1 = pr.m_red * g.pdot / pr.dt_ref;
        const auto res = contact::pinball_force(g, pr, kp);
        check(res.fn > 0.0, "M6 fn > 0 (не обнулилась)");
        check(rel_close(res.fn, kp * f1), "M6 fn == kp*f1 при неизвестном G");
    }

    // M7 (тип 1): p<=0 -> fn==0, dt_limit==kNoLimit
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 0.0; g.pdot = 1.0; g.v_rel = 1.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.g_eff = 1e10; pr.r_eff = 1e-5;
        const auto res = contact::pinball_force(g, pr, 1.0);
        check(res.fn == 0.0 && res.dt_limit == contact::kNoLimit, "M7 тип1 p<=0");
    }
}

void test_type2_swegle()
{
    // M8: обе ветки по очереди доминируют.
    // pdot подобрано так, чтобы p <= pdot*dt_ref, то есть проникание не больше
    // накопленного сближением за опорный шаг. Это область применимости (2.6),
    // в ней потолок силы неактивен и max(f1,f2) проходит как есть (см. M8c).
    {
        // F2 доминирует: большой c_pair.
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 100.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 5000.0;
        const double kp = 1.0;
        const double f1 = pr.m_red * g.p / (pr.dt_ref * pr.dt_ref);
        const double f2 = 2.0 * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r);
        check(std::max(f1, f2) <= pr.m_red * g.pdot / pr.dt_ref,
              "M8 setup case A: потолок неактивен");
        const auto res = contact::swegle_force(g, pr, kp);
        check(rel_close(res.fn, kp * std::max(f1, f2)), "M8 fn == kp*max(f1,f2) [case A]");
        const double k_eff = res.fn / g.p;
        const double expected_dt = contact::kOscSafety / std::sqrt(k_eff / pr.m_red);
        check(rel_close(res.dt_limit, expected_dt), "M8 dt_limit [case A]");
    }
    {
        // Меняем c_pair так, чтобы доминировала другая ветка.
        contact::PairGeom g; g.r = 1e-2; g.p = 1e-6; g.pdot = 100.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 1.0;
        const double kp = 1.0;
        const double f1 = pr.m_red * g.p / (pr.dt_ref * pr.dt_ref);
        const double f2 = 2.0 * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r);
        check(f2 < f1, "M8 setup case B: F1 доминирует");
        check(std::max(f1, f2) <= pr.m_red * g.pdot / pr.dt_ref,
              "M8 setup case B: потолок неактивен");
        const auto res = contact::swegle_force(g, pr, kp);
        check(rel_close(res.fn, kp * std::max(f1, f2)), "M8 fn == kp*max(f1,f2) [case B]");
    }

    // M8c: потолок силы. Проникание вдесятеро больше накопленного сближением
    // за опорный шаг (p = 10*pdot*dt_ref) — (2.6) вне области применимости.
    // Сила обязана срезаться до m_red*pdot/dt_ref (импульса, гасящего
    // сближение), а ограничитель шага — перейти на вязкую ветку: на потолке
    // F от p не зависит, и k_eff = fn/p перестаёт быть жёсткостью.
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 5000.0;
        const double kp = 1.0;
        const double f_raw = std::max(pr.m_red * g.p / (pr.dt_ref * pr.dt_ref),
                                      2.0 * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r));
        const double f_cap = pr.m_red * g.pdot / pr.dt_ref;
        check(f_raw > f_cap, "M8c setup: без потолка сила превышает гасящую");
        const auto res = contact::swegle_force(g, pr, kp);
        check(rel_close(res.fn, kp * f_cap), "M8c fn срезана до m_red*pdot/dt_ref");
        const double expected_dt = contact::kOscSafety * 2.0 * pr.dt_ref / kp;
        check(rel_close(res.dt_limit, expected_dt), "M8c dt_limit по вязкой ветке");
    }

    // M8d: импульс за шаг не превышает нормального импульса пары. Это то
    // свойство, ради которого введён потолок: dt <= dt_ref по построению
    // pair_dt_ref, значит fn*dt <= m_red*pdot при любом проникании.
    {
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 5000.0;
        bool ok = true;
        for (double p_over : {1.0, 10.0, 1e3, 1e6})
        {
            contact::PairGeom g;
            g.r = 1e-3; g.pdot = 1.0; g.v_rel = 0.0;
            g.p = p_over * g.pdot * pr.dt_ref;
            const auto res = contact::swegle_force(g, pr, 1.0);
            if (res.fn * pr.dt_ref > pr.m_red * g.pdot * (1.0 + 1e-12)) ok = false;
        }
        check(ok, "M8d fn*dt_ref <= m_red*pdot при p до 1e6 накопленных");
    }

    // M9: pdot<0 -> fn==0
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = -1.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 5000.0;
        const auto res = contact::swegle_force(g, pr, 1.0);
        check(res.fn == 0.0, "M9 pdot<0 -> fn==0");
    }

    // M7 (тип 2): p<=0
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 0.0; g.pdot = 1.0; g.v_rel = 1.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 5000.0;
        const auto res = contact::swegle_force(g, pr, 1.0);
        check(res.fn == 0.0 && res.dt_limit == contact::kNoLimit, "M7 тип2 p<=0");
    }
}

void test_type3_impulse()
{
    // M10: pdot>0 -> fn == gamma*m_red*pdot/dt точно.
    // Сила строится на НОРМАЛЬНОЙ составляющей: v_rel задана заведомо больше
    // pdot, и если бы формула брала её, значение отличалось бы в 2.5 раза.
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 2.5;
        contact::PairProps pr; pr.m_red = 1e-3;
        const double dt_applied = 1e-8;
        const auto res = contact::impulse_force(g, pr, dt_applied);
        check(rel_close(res.fn, pr.m_red * g.pdot / dt_applied), "M10 fn == m_red*pdot/dt (gamma=1)");
        check(!rel_close(res.fn, pr.m_red * g.v_rel / dt_applied),
              "M10 fn построена на pdot, а не на полной v_rel");
    }

    // M10a: gamma масштабирует силу линейно.
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 2.5;
        contact::PairProps pr; pr.m_red = 1e-3;
        const double dt_applied = 1e-8;
        const auto full = contact::impulse_force(g, pr, dt_applied, 1.0);
        const auto part = contact::impulse_force(g, pr, dt_applied, 0.25);
        check(rel_close(part.fn, 0.25 * full.fn), "M10a fn линейна по gamma");
        check(contact::impulse_force(g, pr, dt_applied, 0.0).fn == 0.0,
              "M10a gamma<=0 -> контакт выключен");
    }

    // M10b: сила непрерывна в точке включения — при pdot -> +0 она стремится
    // к нулю. Это и есть главное отличие от буквальной (2.10), где при том же
    // предельном переходе сила оставалась конечной (m_red*v_rel/dt).
    {
        contact::PairProps pr; pr.m_red = 1e-3;
        const double dt_applied = 1e-8;
        double prev = 1e300;
        for (double pdot : {1e-2, 1e-4, 1e-6, 1e-8})
        {
            contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = pdot; g.v_rel = 2.5;
            const double fn = contact::impulse_force(g, pr, dt_applied).fn;
            check(fn < prev, "M10b fn убывает вместе с pdot");
            prev = fn;
        }
        check(prev < 1e-6 * (pr.m_red * 2.5 / dt_applied), "M10b fn -> 0 при pdot -> +0");
    }

    // M11: pdot<0 -> fn==0 (контакт не тянет)
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = -1.0; g.v_rel = 2.5;
        contact::PairProps pr; pr.m_red = 1e-3;
        const auto res = contact::impulse_force(g, pr, 1e-8);
        check(res.fn == 0.0, "M11 pdot<0 -> fn==0");
    }

    // M7 (тип 3): p<=0
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 0.0; g.pdot = 1.0; g.v_rel = 1.0;
        contact::PairProps pr; pr.m_red = 1e-3;
        const auto res = contact::impulse_force(g, pr, 1e-8);
        check(res.fn == 0.0 && res.dt_limit == contact::kNoLimit, "M7 тип3 p<=0");
    }

    // M13: ограничитель шага — вязкая ветка dt < 0.8*2*dt_applied/gamma.
    // При gamma <= 1.6 он не строже применяемого шага, то есть не срабатывает.
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 2.5;
        contact::PairProps pr; pr.m_red = 1e-3;
        const double dt_applied = 1e-8;
        for (double gamma : {0.25, 1.0, 1.6})
        {
            const auto res = contact::impulse_force(g, pr, dt_applied, gamma);
            check(rel_close(res.dt_limit, contact::kOscSafety * 2.0 * dt_applied / gamma),
                  "M13 dt_limit == 1.6*dt/gamma");
            check(res.dt_limit >= dt_applied, "M13 при gamma<=1.6 ограничитель не режет шаг");
        }
        const auto tight = contact::impulse_force(g, pr, dt_applied, 4.0);
        check(tight.dt_limit < dt_applied, "M13 при gamma>1.6 ограничитель режет шаг");
    }
}

void test_type0_unknown()
{
    // M12: тип 0 / неизвестный тип -> fn==0, dt_limit==kNoLimit
    contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 1.0;
    contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 5000.0;
    pr.g_eff = 1e10; pr.r_eff = 1e-5;

    const auto res0 = contact::penalty_force(contact::Off, g, pr, 1.0, 1e-8);
    check(res0.fn == 0.0 && res0.dt_limit == contact::kNoLimit, "M12 type Off");

    const auto res_unknown = contact::penalty_force(99, g, pr, 1.0, 1e-8);
    check(res_unknown.fn == 0.0 && res_unknown.dt_limit == contact::kNoLimit, "M12 type unknown");
}

// ------------------------- Шаг 4: тесты оркестратора -------------------------
// dim == 1 (param.hpp) — тесты идут в 1D, формулы контакта dim-независимы.

// Пара частиц: itype[0]=3 (тантал, G=69e9), itype[1]=4 (алюминий, G=26e9) —
// оба материала известны реестру shear_modulus_by_itype (нужно для типа 1).
ParticleSystem make_pair(double x0, double x1, double v0, double v1,
                          double m0, double m1, double h0, double h1,
                          double c0, double c1, int itype0 = 3, int itype1 = 4)
{
    ParticleSystem sys;
    sys.resize(2);
    sys.inzializ_index(0);
    sys.inzializ_index(1);
    sys.x[0][0] = x0;  sys.x[0][1] = x1;
    sys.vx[0][0] = v0; sys.vx[0][1] = v1;
    sys.mass[0] = m0;  sys.mass[1] = m1;
    sys.hsml[0] = h0;  sys.hsml[1] = h1;
    sys.c[0] = c0;     sys.c[1] = c1;
    sys.itype[0] = itype0;
    sys.itype[1] = itype1;
    // Радиус тела оркестратор берёт как 0.5*(mass/rho)^(1/dim) (particle_radius
    // в ext_force.cpp), поэтому rho здесь задаётся из требуемого радиуса, а не
    // произвольно: rho = mass/h^dim даёт ровно R = h/2. Тем самым порог контакта
    // остаётся 0.5*(h_i+h_j) — как в тестах M/O до перехода на радиус по объёму,
    // и все их ожидаемые значения сохраняют силу.
    // Без этой строки rho = 0 (obnul_index), радиус вырождается в нуль и контакт
    // не возникает ни в одном тесте.
    sys.rho[0] = m0 / std::pow(h0, static_cast<double>(dim));
    sys.rho[1] = m1 / std::pow(h1, static_cast<double>(dim));
    return sys;
}

SimulationConfig make_cfg(int contact_type, double kp = 1.0)
{
    SimulationConfig cfg;
    cfg.apply_defaults();
    cfg.contact_force_type = contact_type;
    cfg.contact_penalty_scale = kp;
    return cfg;
}

// Тот же критерий, что pair_dt_ref() в ext_force.cpp (static, недоступна извне) —
// продублирована здесь для сверки ожидаемого результата.
double expected_pair_dt_ref(const SimulationConfig& cfg, const ParticleSystem& sys, size_t k)
{
    double v2 = 0.0;
    for (size_t d = 0; d < dim; ++d) v2 += sys.vx[d][k] * sys.vx[d][k];
    const double denom = sys.c[k] + std::sqrt(v2);
    return denom > 0.0 ? cfg.CFL * sys.hsml[k] / denom : 1e30;
}

bool finite_dvxdt(const std::vector<std::vector<double>>& dvxdt, size_t ntotal)
{
    for (size_t d = 0; d < dim; ++d)
        for (size_t i = 0; i < ntotal; ++i)
            if (!std::isfinite(dvxdt[d][i])) return false;
    return true;
}

void test_orchestrator_impulse_balance()
{
    // O1: тип 3, две частицы навстречу, один вызов ext_force -> точный баланс импульса.
    ParticleSystem sys = make_pair(/*x*/0.0, 1e-4, /*v*/1.0, -1.0,
                                    /*m*/2e-6, 3e-6, /*h*/1e-3, 1e-3, /*c*/500.0, 500.0);
    SimulationConfig cfg = make_cfg(contact::Impulse);
    WorkBuffers wb;
    std::vector<size_t> pair_i{0}, pair_j{1};
    std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
    const double dt = 1e-8;

    ext_force(cfg, sys, dt, 2, 2, 1, pair_i, pair_j, dvxdt, wb);

    const double balance = dvxdt[0][0] * sys.mass[0] + dvxdt[0][1] * sys.mass[1];
    check(std::fabs(balance) <= 1e-12 * std::max(std::fabs(dvxdt[0][0] * sys.mass[0]), 1.0),
          "O1 баланс импульса пары (тип 3)");
}

void test_orchestrator_sign()
{
    // O2: типы 1-3, знак силы — ускорение i направлено от j (отталкивание).
    for (int type : {contact::Pinball, contact::Swegle, contact::Impulse})
    {
        ParticleSystem sys = make_pair(0.0, 1e-4, 1.0, -1.0,
                                        2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0);
        SimulationConfig cfg = make_cfg(type);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);

        // i=0 слева от j=1: отталкивание => a[0] < 0 (влево), a[1] > 0 (вправо).
        check(dvxdt[0][0] < 0.0, "O2 знак ускорения частицы i (отталкивание)");
        check(dvxdt[0][1] > 0.0, "O2 знак ускорения частицы j (отталкивание)");
    }
}

void test_orchestrator_filters()
{
    // O3: одинаковый материал -> контакт не регистрируется.
    {
        ParticleSystem sys = make_pair(0.0, 1e-4, 1.0, -1.0, 2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0,
                                        /*itype0*/3, /*itype1*/3);
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        check(wb.contact_pair_count == 0, "O3 одинаковый материал -> pair_count==0");
        check(dvxdt[0][0] == 0.0 && dvxdt[0][1] == 0.0, "O3 одинаковый материал -> dvxdt==0");
    }
    // O3: виртуальная частица (itype < 0).
    {
        ParticleSystem sys = make_pair(0.0, 1e-4, 1.0, -1.0, 2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0,
                                        /*itype0*/3, /*itype1*/-4);
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        check(wb.contact_pair_count == 0, "O3 виртуальная частица -> pair_count==0");
    }
    // O3: itype == 0.
    {
        ParticleSystem sys = make_pair(0.0, 1e-4, 1.0, -1.0, 2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0,
                                        /*itype0*/3, /*itype1*/0);
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        check(wb.contact_pair_count == 0, "O3 itype==0 -> pair_count==0");
    }
}

void test_orchestrator_ghost_excluded()
{
    // O10: зеркальная (виртуальная) частица не образует контактной пары.
    // Отличие от O3 «виртуальная частица»: там признаком служил itype < 0, а
    // virt_part копирует itype оригинала БЕЗ смены знака, поэтому призрак по
    // материалу неотличим от реальной частицы. Единственный признак — индекс:
    // призраки лежат в [nreal, ntotal).
    //
    // Геометрия обеих веток одинакова и заведомо даёт проникание (p > 0),
    // так что различие в результате вызвано ровно границей nreal.
    const auto make_overlapping = []()
    {
        return make_pair(0.0, 1e-4, 1.0, -1.0, 2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0,
                         /*itype0*/3, /*itype1*/4);
    };

    // Контроль: обе реальные -> контакт есть (иначе тест ниже проходил бы
    // по любой причине, а не из-за фильтра по nreal).
    {
        ParticleSystem sys = make_overlapping();
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, /*ntotal*/2, /*nreal*/2, 1, pair_i, pair_j, dvxdt, wb);
        check(wb.contact_pair_count == 1, "O10 контроль: обе реальные -> pair_count==1");
        check(dvxdt[0][0] != 0.0, "O10 контроль: обе реальные -> сила приложена");
    }

    // Частица 1 объявлена зеркальной (nreal == 1) -> пара отбрасывается целиком:
    // ни силы, ни счётчиков, ни ограничителя шага.
    {
        ParticleSystem sys = make_overlapping();
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, /*ntotal*/2, /*nreal*/1, 1, pair_i, pair_j, dvxdt, wb);
        check(wb.contact_pair_count == 0, "O10 призрак -> pair_count==0");
        check(dvxdt[0][0] == 0.0 && dvxdt[0][1] == 0.0, "O10 призрак -> dvxdt==0");
        check(wb.contact_max_delta == 0.0, "O10 призрак -> contact_max_delta не засорён");
        check(wb.min_dt_contact == 1e30, "O10 призрак -> шаг не ограничивается");
    }

    // Тип 1 ставит ограничитель шага, поэтому засорение min_dt_contact
    // фантомной парой проверяется отдельно именно на нём.
    {
        ParticleSystem sys = make_overlapping();
        SimulationConfig cfg = make_cfg(contact::Pinball);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, /*ntotal*/2, /*nreal*/1, 1, pair_i, pair_j, dvxdt, wb);
        check(wb.min_dt_contact == 1e30, "O10 призрак (тип 1) -> min_dt_contact==kNoLimit");
    }

    // Рассогласование вызывающей стороны диагностируется, а не молчит.
    {
        ParticleSystem sys = make_overlapping();
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        bool threw = false;
        try
        {
            ext_force(cfg, sys, 1e-8, /*ntotal*/2, /*nreal*/3, 1, pair_i, pair_j, dvxdt, wb);
        }
        catch (const std::runtime_error&)
        {
            threw = true;
        }
        check(threw, "O10 nreal > ntotal -> исключение");
    }
}

void test_orchestrator_radius_from_volume()
{
    // O11: радиус тела берётся из объёма частицы (mass/rho), а не из hsml.
    // Проверяется то, ради чего сделан переход: в покое соприкасающиеся частицы
    // НЕ проникают друг в друга, и контакт включается только при сближении
    // ближе суммы радиусов.
    //
    // Построение: задаём радиусы напрямую через rho, а hsml делаем заведомо
    // БОЛЬШИМ (десятикратным). Если бы порог по-прежнему считался от hsml,
    // пара срабатывала бы во всех трёх проверках ниже.
    const double dxp = 2e-4;              // «шаг сетки»: расстояние между центрами
    const double R   = 0.5 * dxp;         // радиус тела -> сумма радиусов = dxp
    const double m0 = 2e-6, m1 = 3e-6;
    const double h_big = 10.0 * dxp;      // hsml намеренно не связан с размером тела

    const auto make_at = [&](double gap)
    {
        ParticleSystem sys = make_pair(0.0, gap, 1.0, -1.0, m0, m1,
                                       h_big, h_big, 500.0, 500.0);
        // rho из требуемого радиуса: R = 0.5*(m/rho)^(1/dim) => rho = m/(2R)^dim
        sys.rho[0] = m0 / std::pow(2.0 * R, static_cast<double>(dim));
        sys.rho[1] = m1 / std::pow(2.0 * R, static_cast<double>(dim));
        return sys;
    };

    const auto pair_count_at = [&](double gap)
    {
        ParticleSystem sys = make_at(gap);
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        return wb.contact_pair_count;
    };

    // Расстояние больше суммы радиусов — тела ещё не касаются.
    check(pair_count_at(1.5 * dxp) == 0, "O11 r > R_i+R_j -> контакта нет");
    // Ровно касание: p = 0, контакт требует строгого p > 0.
    check(pair_count_at(dxp) == 0, "O11 r == R_i+R_j (касание) -> контакта нет");
    // Сблизились — контакт возник.
    check(pair_count_at(0.5 * dxp) == 1, "O11 r < R_i+R_j -> контакт есть");

    // Проникание отсчитывается от суммы радиусов, а не от hsml.
    {
        const double gap = 0.5 * dxp;
        ParticleSystem sys = make_at(gap);
        SimulationConfig cfg = make_cfg(contact::Pinball);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        check(rel_close(wb.contact_max_delta, 2.0 * R - gap),
              "O11 p = (R_i+R_j) - r, а не (h_i+h_j)/2 - r");
    }

    // Радиус следует за сжатием: вдвое большая rho -> вдвое меньший объём.
    // При dim == 1 радиус падает вдвое, поэтому касавшиеся тела расходятся.
    {
        ParticleSystem sys = make_at(dxp * 0.9);   // слегка проникают
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        check(wb.contact_pair_count == 1, "O11 контроль: до сжатия контакт есть");

        ParticleSystem dense = make_at(dxp * 0.9);
        dense.rho[0] *= 8.0;  dense.rho[1] *= 8.0;   // объём /8 при любом dim <= 3
        WorkBuffers wb2;
        std::vector<std::vector<double>> dvxdt2(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, dense, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt2, wb2);
        check(wb2.contact_pair_count == 0, "O11 сжатие уменьшает радиус -> контакт пропал");
    }

    // rho <= 0 — размер не восстановить; пара просто не даёт контакта
    // (падать здесь не наша забота: отрицательную плотность ловит int_force).
    {
        ParticleSystem sys = make_at(0.5 * dxp);
        sys.rho[0] = 0.0;
        SimulationConfig cfg = make_cfg(contact::Impulse);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        check(wb.contact_pair_count == 0, "O11 rho <= 0 -> контакта нет, без падения");
    }
}

void test_orchestrator_no_penetration()
{
    // O4: p <= 0 (частицы разведены) -> pair_count==0, min_dt_contact==1e30.
    ParticleSystem sys = make_pair(0.0, 1.0 /* далеко */, 1.0, -1.0,
                                    2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0);
    SimulationConfig cfg = make_cfg(contact::Pinball);
    WorkBuffers wb;
    std::vector<size_t> pair_i{0}, pair_j{1};
    std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
    ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
    check(wb.contact_pair_count == 0, "O4 p<=0 -> pair_count==0");
    check(wb.min_dt_contact == 1e30, "O4 p<=0 -> min_dt_contact==1e30");
}

void test_orchestrator_type_off()
{
    // O5: contact_force_type=0 -> dvxdt совпадает с расчётом без контакта
    // бит-в-бит (для этой пары контакт — единственный источник dvxdt);
    // pair_count==0.
    ParticleSystem sys = make_pair(0.0, 1e-4, 1.0, -1.0, 2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0);
    SimulationConfig cfg = make_cfg(contact::Off);
    WorkBuffers wb;
    std::vector<size_t> pair_i{0}, pair_j{1};
    std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
    ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
    check(wb.contact_pair_count == 0, "O5 type=0 -> pair_count==0");
    check(dvxdt[0][0] == 0.0 && dvxdt[0][1] == 0.0, "O5 type=0 -> dvxdt==0 бит-в-бит");
}

void test_orchestrator_dt_limit_types_1_2()
{
    // O6: ограничитель шага для типов 1-2 равен аналитическому dt_limit
    // (пересчитанному вручную по той же формуле, что и contact_force.hpp).
    for (int type : {contact::Pinball, contact::Swegle})
    {
        const double x0 = 0.0, x1 = 1e-4;
        const double v0 = 1.0, v1 = -1.0;
        const double m0 = 2e-6, m1 = 3e-6;
        const double h0 = 1e-3, h1 = 1e-3;
        const double c0 = 500.0, c1 = 500.0;
        ParticleSystem sys = make_pair(x0, x1, v0, v1, m0, m1, h0, h1, c0, c1);
        SimulationConfig cfg = make_cfg(type);
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        const double dt_applied = 1e-8;
        ext_force(cfg, sys, dt_applied, 2, 2, 1, pair_i, pair_j, dvxdt, wb);

        // Пересчёт ожидаемых величин по той же геометрии/свойствам, что в ext_force.
        const double r = std::fabs(x0 - x1);
        const double p = 0.5 * (h0 + h1) - r;
        const double dx = x0 - x1;
        const double dv = v0 - v1;
        const double g_pdot = -(dx * dv) / r;
        const double m_red = m0 * m1 / (m0 + m1);
        const double c_pair = std::min(c0, c1);
        const double Gi = shear_modulus_by_itype(cfg.mats, sys.itype[0]);
        const double Gj = shear_modulus_by_itype(cfg.mats, sys.itype[1]);
        const double g_eff = (Gi > 0.0 && Gj > 0.0) ? Gi * Gj / (Gi + Gj) : 0.0;
        const double Ri = 0.5 * h0, Rj = 0.5 * h1;
        const double r_eff = Ri * Rj / (Ri + Rj);
        const double dt_ref = std::min(expected_pair_dt_ref(cfg, sys, 0), expected_pair_dt_ref(cfg, sys, 1));

        contact::PairGeom g; g.r = r; g.p = p; g.pdot = g_pdot; g.v_rel = std::fabs(dv);
        contact::PairProps pr; pr.m_red = m_red; pr.c_pair = c_pair; pr.g_eff = g_eff; pr.r_eff = r_eff; pr.dt_ref = dt_ref;

        const auto expected = contact::penalty_force(type, g, pr, cfg.contact_penalty_scale, dt_applied);
        check(rel_close(wb.min_dt_contact, expected.dt_limit), "O6 min_dt_contact совпадает с аналитикой");
    }
}

void test_orchestrator_dt_limit_type_3()
{
    // O7: тип 3 ставит ограничитель вязкой ветки 0.8*2*dt_applied/gamma.
    // Раньше здесь проверялось обратное (min_dt_contact == kNoLimit): у типа 3
    // не было НИ ОДНОГО контура обратной связи — ни потолка силы, ни предела по
    // шагу. При gamma в рабочем диапазоне ограничитель по-прежнему не режет
    // применяемый шаг, но перестаёт быть тождественно отсутствующим.
    ParticleSystem sys = make_pair(0.0, 1e-4, 1.0, -1.0, 2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0);
    const double dt_applied = 1e-8;

    for (double gamma : {0.25, 1.0})
    {
        SimulationConfig cfg = make_cfg(contact::Impulse);
        cfg.contact_impulse_gamma = gamma;
        WorkBuffers wb;
        std::vector<size_t> pair_i{0}, pair_j{1};
        std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
        ext_force(cfg, sys, dt_applied, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        check(rel_close(wb.min_dt_contact, contact::kOscSafety * 2.0 * dt_applied / gamma),
              "O7 min_dt_contact типа 3 совпадает с аналитикой");
        check(wb.min_dt_contact > dt_applied,
              "O7 в рабочем диапазоне gamma шаг не урезается");
    }
}

void test_orchestrator_impulse_normal_component()
{
    // O12: сила типа 3 строится на нормальной составляющей скорости.
    // Две пары с ОДИНАКОВЫМ сближением вдоль линии центров, но разной полной
    // относительной скоростью, должны давать одинаковое ускорение. При
    // буквальной (2.10) ускорения различались бы.
    //
    // dim == 1, поэтому тангенциальной составляющей в самой геометрии нет:
    // разницу вносим прямо в PairGeom и сверяем чистую функцию.
    contact::PairProps pr; pr.m_red = 1.2e-6;
    const double dt_applied = 1e-8;

    contact::PairGeom head_on;  head_on.r = 1e-3; head_on.p = 1e-6;
    head_on.pdot = 3.0; head_on.v_rel = 3.0;           // чистое сближение

    contact::PairGeom sliding;  sliding.r = 1e-3; sliding.p = 1e-6;
    sliding.pdot = 3.0; sliding.v_rel = 30.0;          // то же сближение + скольжение

    const double f_head = contact::impulse_force(head_on, pr, dt_applied).fn;
    const double f_slide = contact::impulse_force(sliding, pr, dt_applied).fn;
    check(rel_close(f_head, f_slide),
          "O12 скольжение не добавляет нормальной силы (тип 3)");
}

void test_orchestrator_diagnostics()
{
    // O8: contact_max_delta равен фактическому p пары.
    const double x0 = 0.0, x1 = 1e-4, h0 = 1e-3, h1 = 1e-3;
    ParticleSystem sys = make_pair(x0, x1, 1.0, -1.0, 2e-6, 3e-6, h0, h1, 500.0, 500.0);
    SimulationConfig cfg = make_cfg(contact::Pinball);
    WorkBuffers wb;
    std::vector<size_t> pair_i{0}, pair_j{1};
    std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));
    ext_force(cfg, sys, 1e-8, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
    const double expected_p = 0.5 * (h0 + h1) - std::fabs(x0 - x1);
    check(rel_close(wb.contact_max_delta, expected_p), "O8 contact_max_delta == p фактический");
}

void test_orchestrator_dt_collapse_regression()
{
    // O9: регрессия коллапса шага. Цикл x64: ext_force(dt) -> dt = min(CFL-модель, wb.min_dt_contact).
    // Ожидание: dt выходит на плато (не убывает монотонно к нулю), ускорения
    // конечны, NaN/Inf нет.
    ParticleSystem sys = make_pair(0.0, 1e-4, 1.0, -1.0, 2e-6, 3e-6, 1e-3, 1e-3, 500.0, 500.0);
    SimulationConfig cfg = make_cfg(contact::Pinball, /*kp*/1.0);
    WorkBuffers wb;
    std::vector<size_t> pair_i{0}, pair_j{1};
    std::vector<std::vector<double>> dvxdt(dim, std::vector<double>(2, 0.0));

    double dt = 1e-7;   // начальный CFL-подобный шаг
    double prev_dt = dt;
    bool collapsed = false;
    bool all_finite = true;
    for (int iter = 0; iter < 64; ++iter)
    {
        ext_force(cfg, sys, dt, 2, 2, 1, pair_i, pair_j, dvxdt, wb);
        if (!finite_dvxdt(dvxdt, 2)) { all_finite = false; break; }

        const double dt_cfl_model = 1e-7;  // фиксированный «CFL» верхний предел
        double next_dt = std::min(dt_cfl_model, wb.min_dt_contact);
        if (next_dt < prev_dt * 1e-6) { collapsed = true; }
        prev_dt = next_dt;
        dt = next_dt;
    }
    check(all_finite, "O9 ускорения конечны на всех итерациях");
    check(!collapsed, "O9 шаг не коллапсирует (выходит на плато)");
    check(dt > 0.0, "O9 итоговый dt положителен");
}

} // namespace

int main()
{
    test_shear_modulus_registry();
    test_type1_pinball();
    test_type2_swegle();
    test_type3_impulse();
    test_type0_unknown();

    test_orchestrator_impulse_balance();
    test_orchestrator_sign();
    test_orchestrator_filters();
    test_orchestrator_ghost_excluded();
    test_orchestrator_radius_from_volume();
    test_orchestrator_no_penetration();
    test_orchestrator_type_off();
    test_orchestrator_dt_limit_types_1_2();
    test_orchestrator_dt_limit_type_3();
    test_orchestrator_impulse_normal_component();
    test_orchestrator_diagnostics();
    test_orchestrator_dt_collapse_regression();

    if (failures == 0)
    {
        std::printf("contact_test: OK\n");
        return 0;
    }
    std::printf("contact_test: FAILED (%d ошибок)\n", failures);
    return 1;
}
