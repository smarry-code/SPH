// test_contact.cpp
// Unit-тест контактных сил «частица-частица» (part_0.md, ADR-007).
// Шаг 2: чистая математика (contact_force.hpp) + реестр модуля сдвига
// (materials.hpp). Аналитика типов 1-3 сверяется с формулами (2.3)-(2.10)
// с относительным допуском 1e-12. Тесты оркестратора (ext_force) — шаг 4.
#include "contact_force.hpp"
#include "materials.hpp"
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
    // M1: itype 4, 10, 11 -> 26e9, 197e9, 1.165e11
    check(rel_close(shear_modulus_by_itype(m, 4), 26e9), "M1 itype 4 -> 26e9");
    check(rel_close(shear_modulus_by_itype(m, 10), 197e9), "M1 itype 10 -> 197e9");
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
    // M8: обе ветки по очереди доминируют
    {
        // F1 доминирует: большой dt_ref делает F1 малым, подберём так, чтобы F2 доминировала сначала.
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 5000.0;
        const double kp = 1.0;
        const double f1 = pr.m_red * g.p / (pr.dt_ref * pr.dt_ref);
        const double f2 = 2.0 * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r);
        const auto res = contact::swegle_force(g, pr, kp);
        check(rel_close(res.fn, kp * std::max(f1, f2)), "M8 fn == kp*max(f1,f2) [case A]");
        const double k_eff = res.fn / g.p;
        const double expected_dt = contact::kOscSafety / std::sqrt(k_eff / pr.m_red);
        check(rel_close(res.dt_limit, expected_dt), "M8 dt_limit [case A]");
    }
    {
        // Меняем c_pair так, чтобы доминировала другая ветка.
        contact::PairGeom g; g.r = 1e-2; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 0.0;
        contact::PairProps pr; pr.m_red = 1e-3; pr.dt_ref = 1e-7; pr.c_pair = 1.0;
        const double kp = 1.0;
        const double f1 = pr.m_red * g.p / (pr.dt_ref * pr.dt_ref);
        const double f2 = 2.0 * pr.m_red * pr.c_pair * g.p / (pr.dt_ref * g.r);
        check(f2 < f1, "M8 setup case B: F1 доминирует");
        const auto res = contact::swegle_force(g, pr, kp);
        check(rel_close(res.fn, kp * std::max(f1, f2)), "M8 fn == kp*max(f1,f2) [case B]");
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
    // M10: pdot>0 -> fn == m_red*v_rel/dt точно; dt_limit == kNoLimit
    {
        contact::PairGeom g; g.r = 1e-3; g.p = 1e-6; g.pdot = 1.0; g.v_rel = 2.5;
        contact::PairProps pr; pr.m_red = 1e-3;
        const double dt_applied = 1e-8;
        const auto res = contact::impulse_force(g, pr, dt_applied);
        check(rel_close(res.fn, pr.m_red * g.v_rel / dt_applied), "M10 fn точно");
        check(res.dt_limit == contact::kNoLimit, "M10 dt_limit == kNoLimit");
    }

    // M11: pdot<0 -> fn==0
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

} // namespace

int main()
{
    test_shear_modulus_registry();
    test_type1_pinball();
    test_type2_swegle();
    test_type3_impulse();
    test_type0_unknown();

    if (failures == 0)
    {
        std::printf("contact_test: OK\n");
        return 0;
    }
    std::printf("contact_test: FAILED (%d ошибок)\n", failures);
    return 1;
}
