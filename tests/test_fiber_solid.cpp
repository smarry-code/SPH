// tests/test_fiber_solid.cpp
// Unit-тесты Fiber-Solid-Contact (fiber/PLAN.md, этап Ф6; fiber/SPH.md,
// раздел «Fiber-Solid-Contact»):
//   - lambda и clamp (середина сегмента, t<0, t>=1 -- обычный/первый/последний
//     сегмент);
//   - контакт с концом нити (изолированная бусина без соседей -- точечный
//     контакт);
//   - узел нити не получает двойную силу (Solid-частица точно в узле между
//     двумя сегментами);
//   - F_a+F_b+F_s=0 (сохранение импульса на primitive);
//   - E, nu материалов проекта против справочных значений.
//
// Стиль -- по образцу tests/test_fiber_contact.cpp (тот же формат
// check()/near_rel()/near_abs(), тот же CpuExecutor/ContributionBuffer
// напрямую). Собирается с SPH_DIM=3 (CMakeLists.txt).
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "materials.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_solid требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near_rel(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b)); }
bool near_abs(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

fiber::FiberMaterial make_material(double E, double nu, double d, double rho, double eps_break = 1.0)
{
    fiber::FiberMaterial m;
    m.E = E; m.nu = nu; m.d = d; m.rho = rho; m.eps_break = eps_break;
    return m;
}

// Прямая цепочка бусин (одна нить, object_id=0, material_id=0).
struct FiberFixture
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

FiberFixture make_chain(const std::vector<std::array<double, 3>>& beads, double mass_each, double E, double nu,
                         double d, double rho)
{
    FiberFixture fx;
    fx.materials.push_back(make_material(E, nu, d, rho));
    const std::size_t n = beads.size();
    fx.state.resize(n);
    fx.topology.resize(n);
    for (std::size_t k = 0; k < n; ++k)
    {
        for (int d3 = 0; d3 < 3; ++d3) fx.state.x[d3][k] = beads[k][d3];
        for (std::size_t dd = 0; dd < dim; ++dd) fx.state.vx[dd][k] = 0.0;
        fx.state.mass[k] = mass_each;
        fx.state.object_id[k] = 0;
        fx.state.material_id[k] = 0;
        if (k > 0) fx.topology.prev[k] = static_cast<long long>(k - 1);
        if (k + 1 < n) fx.topology.next[k] = static_cast<long long>(k + 1);
    }
    return fx;
}

// Одна Solid-частица в ParticleSystem (itype -- реальный материал проекта,
// materials.hpp). rho == rho_0 -- "неповреждённое" состояние (D=0 по
// умолчанию у inzializ_index, elastic_by_itype не зависит от D в любом
// случае, см. комментарий у elastic_by_itype).
std::size_t add_solid_particle(ParticleSystem& sys, std::size_t i, double x0, double x1, double x2, double mass,
                                double rho, int itype)
{
    sys.inzializ_index(i);
    sys.x[0][i] = x0; sys.x[1][i] = x1; sys.x[2][i] = x2;
    sys.mass[i] = mass;
    sys.rho[i] = rho;
    sys.rho_0[i] = rho;
    sys.itype[i] = itype;
    sys.hsml[i] = 1.0; // не используется Fiber-Solid-Contact (нет SPH-ядра)
    return i;
}

interactions::FiberSolidContactPrimitive make_segment(rhs::index_t a, rhs::index_t b, rhs::index_t s,
                                                        long long a_prev, bool b_is_last)
{
    interactions::FiberSolidContactPrimitive p;
    p.a = a; p.b = b; p.s = s;
    p.has_b = true;
    p.a_prev = a_prev;
    p.b_is_last = b_is_last;
    return p;
}

interactions::FiberSolidContactPrimitive make_point(rhs::index_t a, rhs::index_t s)
{
    interactions::FiberSolidContactPrimitive p;
    p.a = a; p.s = s;
    p.has_b = false;
    // a_prev остаётся -1 по умолчанию -- изолированная бусина без соседей.
    return p;
}

// -----------------------------------------------------------------------
// 1) lambda и clamp: сегмент (a=1,b=2) ВНУТРЕННИЙ (a_prev=0!=-1,
//    b_is_last=false) у 4-бусинной прямой нити 0-1-2-3 вдоль оси x, шаг l0;
//    Solid-частица смещена по z на h < r_f+r_s (гарантированный контакт,
//    если сегмент примет вклад).
// -----------------------------------------------------------------------
void test_lambda_and_clamp()
{
    std::printf("test_lambda_and_clamp:\n");
    const double l0 = 2.0e-3;
    const double E_f = 2.0e9, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    const double h = 0.3e-3; // зазор по z -- гарантирует контакт при попадании в [0,1]

    std::vector<std::array<double, 3>> beads = {{0, 0, 0}, {l0, 0, 0}, {2 * l0, 0, 0}, {3 * l0, 0, 0}};
    FiberFixture fx = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);

    Materials mats; // проектные материалы (materials.hpp)
    ParticleSystem sys;
    sys.resize(8);
    const double rho_s = 2700.0;
    const double mass_s = 1e-7;

    // (a) t=0.5 (середина сегмента 1-2, x=1.5*l0) -- нормальный контакт.
    {
        add_solid_particle(sys, 0, 1.5 * l0, 0.0, h, mass_s, rho_s, 4);
        interactions::FiberSolidContactPrimitive prim = make_segment(1, 2, 0, /*a_prev=*/0, false);
        interactions::FiberSolidContactGeometry geom =
            interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, prim);
        check(geom.contact, "t=0.5 (середина внутреннего сегмента): контакт есть");
        check(near_rel(geom.lambda, 0.5, 1e-9), "t=0.5: lambda == 0.5");
        check(near_rel(geom.D, h, 1e-9), "t=0.5: D == h");
    }

    // (b) t ~ -0.3*l0/l0 = -0.3 (Solid ближе к бусине 0, ЗА пределами сегмента
    //     1-2 со стороны a) -- сегмент НЕ первый (prev(1)=0!=-1) -> пропуск,
    //     ответственность передана сегменту (0,1).
    {
        add_solid_particle(sys, 1, 0.7 * l0, 0.0, h, mass_s, rho_s, 4);
        interactions::FiberSolidContactPrimitive prim = make_segment(1, 2, 1, /*a_prev=*/0, false);
        interactions::FiberSolidContactGeometry geom =
            interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, prim);
        check(!geom.contact, "t<0 у НЕ первого сегмента: primitive пропущен (ответственность у предыдущего сегмента)");
    }

    // (c) t ~ 1.3 (Solid ЗА пределами сегмента 1-2 со стороны b) -- сегмент
    //     НЕ последний (next(2)=3!=-1) -> пропуск.
    {
        add_solid_particle(sys, 2, 2.3 * l0, 0.0, h, mass_s, rho_s, 4);
        interactions::FiberSolidContactPrimitive prim = make_segment(1, 2, 2, /*a_prev=*/0, false);
        interactions::FiberSolidContactGeometry geom =
            interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, prim);
        check(!geom.contact, "t>=1 у НЕ последнего сегмента: primitive пропущен (ответственность у следующего сегмента)");
    }

    // (d) ПЕРВЫЙ сегмент (0,1): a_prev=-1 -- t<0 клэмпится в 0 (законно).
    // Небольшой (не 0.3*l0, как в (b)) отрицательный оффсет по x: клэмп
    // проецирует Solid-частицу на q=x[0] (не на исходную точку под ней), а
    // проникновение delta=r_f+r_s-D считается уже от q -- слишком большой
    // боковой оффсет сделал бы D>r_f+r_s (контакта нет ни у какого сегмента).
    // Здесь цель -- проверить именно факт клэмпа при t<0, а не геометрию D.
    {
        add_solid_particle(sys, 3, -0.05 * l0, 0.0, h, mass_s, rho_s, 4);
        interactions::FiberSolidContactPrimitive prim = make_segment(0, 1, 3, /*a_prev=*/-1, false);
        interactions::FiberSolidContactGeometry geom =
            interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, prim);
        check(geom.contact, "t<0 у ПЕРВОГО сегмента: клэмп в 0 -- контакт есть");
        check(near_abs(geom.lambda, 0.0, 1e-12), "t<0 у первого сегмента: lambda == 0 (клэмп)");
    }

    // (e) ПОСЛЕДНИЙ сегмент (2,3): b_is_last=true -- t>=1 клэмпится в 1
    // (тот же небольшой оффсет, см. комментарий у (d)).
    {
        add_solid_particle(sys, 4, 3.05 * l0, 0.0, h, mass_s, rho_s, 4);
        interactions::FiberSolidContactPrimitive prim = make_segment(2, 3, 4, /*a_prev=*/1, true);
        interactions::FiberSolidContactGeometry geom =
            interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, prim);
        check(geom.contact, "t>=1 у ПОСЛЕДНЕГО сегмента: клэмп в 1 -- контакт есть");
        check(near_abs(geom.lambda, 1.0, 1e-12), "t>=1 у последнего сегмента: lambda == 1 (клэмп)");
    }
}

// -----------------------------------------------------------------------
// 2) Контакт с концом нити: изолированная бусина без соседей (нить из одной
//    бусины, has_b=false) -- точечный контакт, сила целиком на бусину.
// -----------------------------------------------------------------------
void test_end_of_fiber_point_contact()
{
    std::printf("test_end_of_fiber_point_contact:\n");
    const double E_f = 2.0e9, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    const double h = 0.3e-3;

    std::vector<std::array<double, 3>> beads = {{0, 0, 0}}; // одна бусина -- prev=next=-1
    FiberFixture fx = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);
    check(fx.topology.prev[0] == -1 && fx.topology.next[0] == -1, "самопроверка: бусина изолирована (нет соседей)");

    Materials mats;
    ParticleSystem sys;
    sys.resize(1);
    add_solid_particle(sys, 0, 0.0, 0.0, h, 1e-7, 2700.0, 4);

    interactions::FiberSolidContactPrimitive prim = make_point(0, 0);
    interactions::FiberSolidContactGeometry geom =
        interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, prim);
    check(geom.contact, "точечный контакт: есть проникновение");
    check(near_rel(geom.D, h, 1e-9), "точечный контакт: D == h (расстояние до самой бусины)");
    check(near_abs(geom.n[0], 0.0, 1e-9) && near_abs(geom.n[1], 0.0, 1e-9) && near_abs(geom.n[2], 1.0, 1e-9),
          "точечный контакт: n вдоль +z (от бусины к Solid)");

    rhs::CpuExecutor exec(2);
    std::vector<interactions::FiberSolidContactPrimitive> primitives = {prim};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
    interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx.state, fx.materials, primitives, fbuf, sbuf);
    std::vector<fiber::FiberDelta> fd = rhs::reduce_contributions(exec, fbuf, fx.state.size());
    std::vector<interactions::SolidDelta> sd = rhs::reduce_contributions(exec, sbuf, sys.size());

    const double Fa = fx.state.mass[0] * fd[0].dv[2];
    const double Fs = sys.mass[0] * sd[0].dv[2];
    check(Fa < 0.0, "точечный контакт: сила на бусину -- в сторону -z (от Solid)");
    check(Fs > 0.0, "точечный контакт: сила на Solid -- в сторону +z (от бусины)");
    check(near_rel(Fa, -Fs, 1e-9), "точечный контакт: F_a == -F_s (вся сила -- одной бусине, распределения нет)");
}

// -----------------------------------------------------------------------
// 3) Узел нити не получает двойную силу: 3-бусинная нить 0-1-2, Solid-частица
//    ТОЧНО в узле (бусина 1). Оба соседних primitives ((0,1,s) и (1,2,s))
//    строятся, как их построил бы build_fiber_solid_contact_search (по
//    решению 2 в шапке include/interactions/fiber_solid_contact.hpp), и
//    суммарный вклад сравнивается с ОДНИМ эталонным точечным контактом в
//    бусине 1.
// -----------------------------------------------------------------------
void test_node_no_double_force()
{
    std::printf("test_node_no_double_force:\n");
    const double l0 = 2.0e-3;
    const double E_f = 2.0e9, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    const double h = 0.3e-3;

    std::vector<std::array<double, 3>> beads = {{0, 0, 0}, {l0, 0, 0}, {2 * l0, 0, 0}};
    FiberFixture fx = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);

    Materials mats;
    ParticleSystem sys;
    sys.resize(1);
    // Solid-частица точно над бусиной 1 (узлом) -- x=l0, z=h.
    add_solid_particle(sys, 0, l0, 0.0, h, 1e-7, 2700.0, 4);

    // Сегмент (0,1): b=1 -- узел, next(1)=2 -> b_is_last=false. Solid у x=l0
    // даёт t=1 РОВНО -- сегмент НЕ последний -> пропуск (см. решение 3 в
    // шапке fiber_solid_contact.hpp).
    interactions::FiberSolidContactPrimitive seg01 = make_segment(0, 1, 0, /*a_prev=*/-1, /*b_is_last=*/false);
    // Сегмент (1,2): a=1 -- узел, prev(1)=0 -> a_prev=0. Solid у x=l0=x[a]
    // даёт t=0 РОВНО -- правило владения пересчитывает t_prev для (0,1) той
    // же функцией: там та же Solid-частица у x=l0 даёт t_prev=1 РОВНО
    // (симметрично seg01), т.е. t_prev>=1 -- сегмент (1,2) подтверждён
    // владельцем узла и принимает контакт с lambda=0.
    interactions::FiberSolidContactPrimitive seg12 = make_segment(1, 2, 0, /*a_prev=*/0, /*b_is_last=*/false);

    interactions::FiberSolidContactGeometry g01 =
        interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, seg01);
    interactions::FiberSolidContactGeometry g12 =
        interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, seg12);
    check(!g01.contact, "узел: сегмент (0,1) НЕ даёт вклад (t=1 у не-последнего сегмента -- пропуск)");
    check(g12.contact, "узел: сегмент (1,2) даёт вклад (t=0 -- обычный нижний конец диапазона)");
    check(near_abs(g12.lambda, 0.0, 1e-12), "узел: сегмент (1,2) -- lambda==0 (вся сила бусине 1)");

    // Суммарный вклад через kernel по ОБОИМ primitives.
    rhs::CpuExecutor exec(2);
    std::vector<interactions::FiberSolidContactPrimitive> primitives = {seg01, seg12};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
    interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx.state, fx.materials, primitives, fbuf, sbuf);
    std::vector<fiber::FiberDelta> fd = rhs::reduce_contributions(exec, fbuf, fx.state.size());
    std::vector<interactions::SolidDelta> sd = rhs::reduce_contributions(exec, sbuf, sys.size());

    // Эталон: ОДИН точечный контакт бусины 1 с той же Solid-частицей (та же
    // геометрия, что даёт seg12 при lambda=0 -- F_a целиком на бусину-узел).
    interactions::FiberSolidContactPrimitive point1 = make_point(1, 0);
    interactions::FiberSolidContactGeometry gpoint =
        interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, point1);
    check(gpoint.contact, "самопроверка эталона: точечный контакт бусины 1 тоже даёт проникновение");
    check(near_rel(gpoint.F, g12.F, 1e-9), "самопроверка эталона: F совпадает с F сегмента (1,2) при lambda=0 (та же геометрия)");

    const double F_solid_total = sys.mass[0] * sd[0].dv[2];
    const double F_solid_single = gpoint.F; // F^c_fs, приложенная к Solid вдоль +n (n==+z здесь)
    check(near_rel(F_solid_total, F_solid_single, 1e-9),
          "узел: суммарная сила на Solid равна силе ОДНОГО контакта (не удвоена)");

    const double F_bead1 = fx.state.mass[1] * fd[1].dv[2];
    check(near_rel(F_bead1, -F_solid_single, 1e-9),
          "узел: сила на бусину 1 равна -F одного контакта (не удвоена, вся сила -- узлу, не бусине 2)");
    check(near_abs(fx.state.mass[2] * fd[2].dv[2], 0.0, 1e-15),
          "узел: бусина 2 не получает вклада (lambda=0 у сегмента (1,2), сегмент (0,1) пропущен целиком)");
}

// -----------------------------------------------------------------------
// 4) F_a+F_b+F_s=0 -- сохранение импульса на primitive (общий случай,
//    lambda интерьерная, не 0/1).
// -----------------------------------------------------------------------
void test_momentum_conservation()
{
    std::printf("test_momentum_conservation:\n");
    const double l0 = 2.0e-3;
    const double E_f = 2.0e9, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    const double h = 0.3e-3;

    std::vector<std::array<double, 3>> beads = {{0, 0, 0}, {l0, 0, 0}, {2 * l0, 0, 0}, {3 * l0, 0, 0}};
    // Разные массы бусин -- проверка не на F, а на m*a (как test_fiber_contact.cpp).
    FiberFixture fx = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);
    fx.state.mass[1] = 1.3e-6;
    fx.state.mass[2] = 0.8e-6;

    Materials mats;
    ParticleSystem sys;
    sys.resize(1);
    add_solid_particle(sys, 0, 1.35 * l0, 0.2 * l0, h, 1.7e-7, 2700.0, 4); // lambda не 0/1, не по оси x

    interactions::FiberSolidContactPrimitive prim = make_segment(1, 2, 0, /*a_prev=*/0, false);
    interactions::FiberSolidContactGeometry geom =
        interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, prim);
    check(geom.contact, "самопроверка: контакт есть");
    check(geom.lambda > 0.0 && geom.lambda < 1.0, "самопроверка: lambda строго внутри (0,1)");

    rhs::CpuExecutor exec(2);
    std::vector<interactions::FiberSolidContactPrimitive> primitives = {prim};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
    interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx.state, fx.materials, primitives, fbuf, sbuf);
    std::vector<fiber::FiberDelta> fd = rhs::reduce_contributions(exec, fbuf, fx.state.size());
    std::vector<interactions::SolidDelta> sd = rhs::reduce_contributions(exec, sbuf, sys.size());

    double total[dim] = {};
    double scale = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        const double Fa = fx.state.mass[1] * fd[1].dv[d];
        const double Fb = fx.state.mass[2] * fd[2].dv[d];
        const double Fs = sys.mass[0] * sd[0].dv[d];
        total[d] = Fa + Fb + Fs;
        scale = std::max({scale, std::fabs(Fa), std::fabs(Fb), std::fabs(Fs)});
    }
    if (scale < 1e-30) scale = 1.0;
    for (std::size_t d = 0; d < dim; ++d)
        check(std::fabs(total[d]) < 1e-9 * scale,
              "F_a+F_b+F_s == 0 по компоненте " + std::to_string(d) + ": " + std::to_string(total[d]));
}

// -----------------------------------------------------------------------
// 5) E, nu против справочных значений (fiber/PLAN.md, этап Ф6, пункт 1).
//    Источник сравнения: справочные диапазоны для конструкционных
//    материалов (см. комментарии у каждой проверки) -- допуск НЕ калибровка,
//    а sanity-check порядка величины и физичности 0<nu<0.5.
// -----------------------------------------------------------------------
void test_elastic_by_itype_reference_values()
{
    std::printf("test_elastic_by_itype_reference_values:\n");
    Materials mats;

    // Газ/вода (itype 1-2): не определено.
    check(!elastic_by_itype(mats, 1).has_value(), "itype=1 (газ): E/nu не определены");
    check(!elastic_by_itype(mats, 2).has_value(), "itype=2 (вода): E/nu не определены");

    // Тантал (itype=3): справочно E ~ 186 ГПа (напр. ASM Handbook), здесь K=a*z
    // даёт "холодный" модуль по УРС проекта -- допуск щедрый (порядок величины).
    {
        auto e = elastic_by_itype(mats, 3);
        check(e.has_value(), "itype=3 (Ta): определено");
        if (e)
        {
            std::printf("    Ta: E=%.4g nu=%.4g\n", e->E, e->nu);
            check(near_rel(e->E, 186e9, 0.3), "Ta: E в пределах 30% от справочных ~186 ГПа");
            check(e->nu > 0.0 && e->nu < 0.5, "Ta: 0 < nu < 0.5 (физичность)");
        }
    }

    // Алюминий (itype=4): справочно E ~ 68-70 ГПа, nu ~ 0.33 (общеизвестные
    // табличные значения для Al 2024/6061-подобных сплавов); заметка в памяти
    // проекта даёт на ТЕКУЩИХ константах E=70.6 ГПа, nu=0.357.
    {
        auto e = elastic_by_itype(mats, 4);
        check(e.has_value(), "itype=4 (Al): определено");
        if (e)
        {
            std::printf("    Al: E=%.4g nu=%.4g\n", e->E, e->nu);
            check(near_rel(e->E, 70.6e9, 0.05), "Al: E совпадает с заметкой памяти проекта (70.6 ГПа, <=5%)");
            check(near_rel(e->nu, 0.357, 0.05), "Al: nu совпадает с заметкой памяти проекта (0.357, <=5%)");
            check(near_rel(e->E, 69e9, 0.2), "Al: E в пределах 20% от справочных ~68-70 ГПа");
            check(near_rel(e->nu, 0.33, 0.2), "Al: nu в пределах 20% от справочных ~0.33");
        }
    }

    // B4C (itype=10, JH-2): K=K1=230.5e9, G=199.4e9 (materials.hpp) ->
    // справочно E ~ 440-460 ГПа, nu ~ 0.17-0.21 (карбид бора).
    {
        auto e = elastic_by_itype(mats, 10);
        check(e.has_value(), "itype=10 (B4C): определено");
        if (e)
        {
            std::printf("    B4C: E=%.4g nu=%.4g\n", e->E, e->nu);
            check(near_rel(e->E, 464e9, 0.05), "B4C: E совпадает с заметкой памяти проекта (464 ГПа, <=5%)");
            check(near_rel(e->nu, 0.164, 0.05), "B4C: nu совпадает с заметкой памяти проекта (0.164, <=5%)");
            check(e->E > 400e9 && e->E < 500e9, "B4C: E в разумном справочном диапазоне (400-500 ГПа)");
            check(e->nu > 0.1 && e->nu < 0.3, "B4C: nu в разумном справочном диапазоне (0.1-0.3)");
        }
    }
}

// =======================================================================
// Тесты fiber/fiber_solid_contact_fix.md, §3.1-§3.5. Общие правила:
//   - primitives строятся РЕАЛЬНЫМ поиском build_fiber_solid_contact_search
//     (реальный радиус, реальное правило "сегмент принадлежит первой бусине");
//   - все проверки сил -- по суммарным вкладам после reduce_contributions по
//     ВСЕМ primitives сцены;
//   - число контактов = число primitives с compute_fiber_solid_contact_geometry
//     (...).contact.
// =======================================================================

using Vec3 = std::array<double, 3>;

Vec3 v_add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 v_sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 v_mul(const Vec3& a, double k) { return {a[0] * k, a[1] * k, a[2] * k}; }
double v_dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double v_norm(const Vec3& a) { return std::sqrt(v_dot(a, a)); }
Vec3 v_unit(const Vec3& a) { return v_mul(a, 1.0 / v_norm(a)); }

// Общие параметры сцен (как в существующих тестах).
constexpr double kL0 = 2.0e-3;
constexpr double kEf = 2.0e9, kNuF = 0.3, kDf = 1.0e-3, kRhoF = 1200.0;
constexpr double kMassF = 1e-6;
constexpr double kMassS = 1e-7, kRhoS = 2700.0;
constexpr int kItypeS = 4; // Al

// Цепочка с заполненными topology.bonds (номинальная l0, l_break большой):
// нужны build_fiber_solid_contact_search (радиус по ТЕКУЩЕЙ длине связей).
FiberFixture make_bonded_chain(const std::vector<Vec3>& beads, double l0_nominal)
{
    std::vector<std::array<double, 3>> b(beads.begin(), beads.end());
    FiberFixture fx = make_chain(b, kMassF, kEf, kNuF, kDf, kRhoF);
    for (std::size_t k = 0; k + 1 < beads.size(); ++k)
    {
        fiber::FiberBond bond;
        bond.i = k; bond.j = k + 1; bond.l0 = l0_nominal; bond.l_break = 100.0 * l0_nominal;
        fx.topology.bonds.push_back(bond);
    }
    return fx;
}

struct SceneResult
{
    std::vector<interactions::FiberSolidContactPrimitive> primitives;
    std::vector<interactions::FiberSolidContactGeometry> geoms; // по primitives
    std::vector<std::size_t> contact_idx;                        // индексы primitives с contact
    Vec3 F_solid{};                                              // суммарная сила на Solid s
    std::vector<Vec3> F_bead;                                    // суммарная сила на каждую бусину
    std::vector<Vec3> q;                                         // точки q принятых контактов
};

SceneResult run_scene(const Materials& mats, const ParticleSystem& sys, std::size_t s, const FiberFixture& fx)
{
    SceneResult r;
    rhs::CpuExecutor exec(2);
    interactions::FiberSolidContactSearchResult search;
    interactions::build_fiber_solid_contact_search(exec, mats, sys, sys.size(), fx.state, fx.topology, fx.materials, search);
    r.primitives = search.primitives;
    for (std::size_t k = 0; k < r.primitives.size(); ++k)
    {
        auto g = interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, r.primitives[k]);
        r.geoms.push_back(g);
        if (g.contact)
        {
            r.contact_idx.push_back(k);
            const auto& p = r.primitives[k];
            Vec3 qa = {fx.state.x[0][p.a], fx.state.x[1][p.a], fx.state.x[2][p.a]};
            if (p.has_b)
            {
                Vec3 qb = {fx.state.x[0][p.b], fx.state.x[1][p.b], fx.state.x[2][p.b]};
                qa = v_add(qa, v_mul(v_sub(qb, qa), g.lambda));
            }
            r.q.push_back(qa);
        }
    }
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
    interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx.state, fx.materials, r.primitives, fbuf, sbuf);
    std::vector<fiber::FiberDelta> fd = rhs::reduce_contributions(exec, fbuf, fx.state.size());
    std::vector<interactions::SolidDelta> sd = rhs::reduce_contributions(exec, sbuf, sys.size());
    for (std::size_t d = 0; d < 3; ++d) r.F_solid[d] = sys.mass[s] * sd[s].dv[d];
    r.F_bead.assign(fx.state.size(), Vec3{0, 0, 0});
    for (std::size_t k = 0; k < fx.state.size(); ++k)
        for (std::size_t d = 0; d < 3; ++d) r.F_bead[k][d] = fx.state.mass[k] * fd[k].dv[d];
    return r;
}

// Сила точечного контакта бусины a с частицей s (вектор на Solid = F*n).
Vec3 point_contact_force_on_solid(const Materials& mats, const ParticleSystem& sys, const FiberFixture& fx,
                                   rhs::index_t a, rhs::index_t s, bool* has_contact = nullptr)
{
    auto g = interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, make_point(a, s));
    if (has_contact) *has_contact = g.contact;
    return {g.F * g.n[0], g.F * g.n[1], g.F * g.n[2]};
}

double solid_r_f_plus_r_s(const ParticleSystem& sys, std::size_t s)
{
    return 0.5 * kDf + interactions::solid_contact_radius(sys, s);
}

ParticleSystem make_one_solid(const Vec3& x)
{
    ParticleSystem sys;
    sys.resize(1);
    add_solid_particle(sys, 0, x[0], x[1], x[2], kMassS, kRhoS, kItypeS);
    return sys;
}

Vec3 sum_all_forces(const SceneResult& r)
{
    Vec3 tot = r.F_solid;
    for (const auto& fb : r.F_bead) tot = v_add(tot, fb);
    return tot;
}

// -----------------------------------------------------------------------
// 3.1) Вершина V: излом нити, частица во ВНЕШНЕМ секторе под вершиной.
//    Нить из 5 бусин, вершина -- внутренняя бусина 2. Плечи u1, u2 (единичные)
//    образуют угол theta; частица -- на биссектрисе снаружи, на расстоянии
//    0.7*(r_f+r_s) (проникновение 0.3*(r_f+r_s)). Старое правило теряло
//    контакт целиком (оба соседних сегмента отбрасывали его).
// -----------------------------------------------------------------------
struct KinkResult
{
    std::size_t n_contacts;
    double rel_err_solid;   // |F_solid - F_point| / |F_point|
    double rel_err_vertex;  // |F_vertex + F_point| / |F_point|
    double max_other_rel;   // max |F_other_bead| / |F_point|
};

KinkResult run_kink_scene(const std::string& label, const Vec3& u1, const Vec3& u2, const Vec3& w1, const Vec3& w2)
{
    const Vec3 V = {0.0, 0.0, 0.0};
    const Vec3 b1 = v_add(V, v_mul(u1, kL0));
    const Vec3 b3 = v_add(V, v_mul(u2, kL0));
    const Vec3 b0 = v_add(b1, v_mul(w1, kL0));
    const Vec3 b4 = v_add(b3, v_mul(w2, kL0));
    FiberFixture fx = make_bonded_chain({b0, b1, V, b3, b4}, kL0);

    const Vec3 bis = v_unit(v_add(u1, u2));
    ParticleSystem sys = make_one_solid({0, 0, 0});
    const double rsum = solid_r_f_plus_r_s(sys, 0);
    const double Dsep = 0.7 * rsum; // проникновение 0.3*(r_f+r_s)
    const Vec3 xs = v_sub(V, v_mul(bis, Dsep));
    for (int d = 0; d < 3; ++d) sys.x[d][0] = xs[d];

    Materials mats;
    SceneResult r = run_scene(mats, sys, 0, fx);
    bool pt_ok = false;
    const Vec3 Fp = point_contact_force_on_solid(mats, sys, fx, 2, 0, &pt_ok);
    const double Fmag = v_norm(Fp);
    check(pt_ok && Fmag > 0.0, label + ": самопроверка -- точечный контакт вершины с частицей есть");

    KinkResult kr{};
    kr.n_contacts = r.contact_idx.size();
    check(kr.n_contacts == 1, label + ": ровно один контакт (получено " + std::to_string(kr.n_contacts) + ")");
    kr.rel_err_solid = v_norm(v_sub(r.F_solid, Fp)) / Fmag;
    check(kr.rel_err_solid <= 1e-12, label + ": сила на Solid == силе точечного контакта вершины (отн. ошибка " +
                                          std::to_string(kr.rel_err_solid) + ")");
    kr.rel_err_vertex = v_norm(v_add(r.F_bead[2], Fp)) / Fmag;
    check(kr.rel_err_vertex <= 1e-12, label + ": вся сила нити -- на бусине-вершине (отн. ошибка " +
                                           std::to_string(kr.rel_err_vertex) + ")");
    kr.max_other_rel = 0.0;
    for (std::size_t k = 0; k < r.F_bead.size(); ++k)
        if (k != 2) kr.max_other_rel = std::max(kr.max_other_rel, v_norm(r.F_bead[k]) / Fmag);
    check(kr.max_other_rel <= 1e-12, label + ": соседние бусины -- нулевой вклад (max отн. " +
                                          std::to_string(kr.max_other_rel) + ")");
    return kr;
}

void test_kink_vertex()
{
    std::printf("test_kink_vertex (3.1):\n");
    const double pi = std::acos(-1.0);
    struct Row { std::string name; KinkResult kr; };
    std::vector<Row> rows;
    for (double deg : {90.0, 170.0, 120.0, 60.0, 20.0})
    {
        const double h = 0.5 * deg * pi / 180.0;
        const Vec3 u1 = {-std::sin(h), std::cos(h), 0.0}, u2 = {std::sin(h), std::cos(h), 0.0};
        char label[64];
        std::snprintf(label, sizeof label, "излом %.0f град", deg);
        rows.push_back({label, run_kink_scene(label, u1, u2, u1, u2)});
    }
    {
        // Излом ВНЕ плоскости: плечи u1,u2 в разных плоскостях, и внешние
        // бусины 0 и 4 уходят в стороны от плоскостей плеч (нить настоящая 3D).
        const Vec3 u1 = v_unit({-0.5, 0.7, 0.4}), u2 = v_unit({0.6, 0.6, -0.5});
        const Vec3 w1 = v_unit(v_add(u1, {0.3, 0.0, 0.6})), w2 = v_unit(v_add(u2, {0.0, 0.5, 0.4}));
        rows.push_back({"3D (вне плоскости)", run_kink_scene("излом 3D", u1, u2, w1, w2)});
    }
    std::printf("  таблица 3.1: угол -> контактов, отн. ошибка силы на Solid, на вершине, max по соседям\n");
    for (const auto& r : rows)
        std::printf("    %-20s contacts=%zu  err_solid=%.3e  err_vertex=%.3e  err_other=%.3e\n", r.name.c_str(),
                    r.kr.n_contacts, r.kr.rel_err_solid, r.kr.rel_err_vertex, r.kr.max_other_rel);
}

// -----------------------------------------------------------------------
// 3.2) Вогнутый излом: частица ВНУТРИ V (нить изогнута к частице), проецируется
//    внутрь обоих сегментов и касается обоих -> два контакта в разных точках q,
//    суммарный импульс (Solid + нить) ноль.
// -----------------------------------------------------------------------
void test_concave_kink()
{
    std::printf("test_concave_kink (3.2):\n");
    const Vec3 offsets[2] = {{0.0, 0.2 * kL0, 0.0}, {0.05 * kL0, 0.25 * kL0, 0.03 * kL0}}; // симм. и несимм. положение
    for (int c = 0; c < 2; ++c)
    {
        const std::string tag = c == 0 ? "симметричная частица" : "несимметричная частица";
        // V-образная нить, вершина (бусина 2) внизу, частица над вершиной.
        FiberFixture fx = make_bonded_chain({{-2 * kL0, 2 * kL0, 0}, {-kL0, kL0, 0}, {0, 0, 0}, {kL0, kL0, 0}, {2 * kL0, 2 * kL0, 0}},
                                            std::sqrt(2.0) * kL0);
        ParticleSystem sys = make_one_solid(offsets[c]);
        Materials mats;
        SceneResult r = run_scene(mats, sys, 0, fx);
        check(r.contact_idx.size() == 2, tag + ": два контакта (получено " + std::to_string(r.contact_idx.size()) + ")");
        if (r.contact_idx.size() == 2)
        {
            const double sep = v_norm(v_sub(r.q[0], r.q[1]));
            check(sep > 1e-3 * kL0, tag + ": точки q разные (расстояние " + std::to_string(sep) + ")");
            for (std::size_t k : r.contact_idx)
            {
                const auto& g = r.geoms[k];
                check(g.lambda > 0.0 && g.lambda < 1.0, tag + ": lambda внутри (0,1): " + std::to_string(g.lambda));
            }
        }
        const double scale = v_norm(r.F_solid);
        const double res = v_norm(sum_all_forces(r));
        check(scale > 0.0 && res <= 1e-12 * scale,
              tag + ": суммарный импульс ноль (|sum|/|F_s| = " + std::to_string(res / std::max(scale, 1e-300)) + ")");
    }
}

// -----------------------------------------------------------------------
// 3.3) Округление у узла: прямая нить, частица над узлом (бусина 2 из 0..4) со
//    смещением вдоль нити 0, +-1e-13*l0, +-1e-6*l0: ровно один контакт,
//    суммарная сила == силе ОДНОГО контакта.
// -----------------------------------------------------------------------
void test_node_rounding()
{
    std::printf("test_node_rounding (3.3):\n");
    std::vector<Vec3> beads;
    for (int k = 0; k < 5; ++k) beads.push_back({k * kL0, 0.0, 0.0});
    FiberFixture fx = make_bonded_chain(beads, kL0);
    const double h = 0.3e-3;
    for (double off : {0.0, 1e-13, -1e-13, 1e-6, -1e-6})
    {
        ParticleSystem sys = make_one_solid({2 * kL0 + off * kL0, 0.0, h});
        Materials mats;
        SceneResult r = run_scene(mats, sys, 0, fx);
        char tag[64];
        std::snprintf(tag, sizeof tag, "смещение %+.0e*l0", off);
        check(r.contact_idx.size() == 1, std::string(tag) + ": ровно один контакт (получено " +
                                              std::to_string(r.contact_idx.size()) + ")");
        // эталон -- ОДИН контакт (точечный контакт бусины-узла): сравниваем МОДУЛЬ
        // силы (D отличается на ~(off*l0)^2/h -- пренебрежимо); направление у
        // сегмента при off!=0 строго +z (q под частицей), а у точечного эталона
        // наклонено на off*l0/h ~ 7e-6, поэтому векторно эталон не сравнивается.
        const Vec3 Fp = point_contact_force_on_solid(mats, sys, fx, 2, 0);
        const double err = std::fabs(v_norm(r.F_solid) - v_norm(Fp)) / v_norm(Fp);
        check(err <= 1e-9, std::string(tag) + ": |сила на Solid| == силе одного контакта (отн. ошибка " + std::to_string(err) + ")");
        check(std::fabs(r.F_solid[2] - v_norm(r.F_solid)) <= 1e-12 * v_norm(r.F_solid), std::string(tag) + ": сила на Solid направлена вдоль +z (нормаль к нити)");
        check(v_norm(sum_all_forces(r)) <= 1e-12 * v_norm(r.F_solid), std::string(tag) + ": импульс сохраняется");
    }
}

// -----------------------------------------------------------------------
// 3.3a) Случайные ломаные (1000 шт., 6 бусин, шаги ~l0 со случайным
//    направлением; std::mt19937_64 с фиксированным seed) и случайная частица
//    рядом. Для каждой: min D по принятым контактам == глобальному минимуму
//    расстояния от частицы до ломаной (контакт не теряется), нет двух принятых
//    контактов с совпадающей q (нет удвоения).
// -----------------------------------------------------------------------
double dist_point_segment(const Vec3& a, const Vec3& b, const Vec3& p)
{
    const Vec3 e = v_sub(b, a);
    double t = v_dot(v_sub(p, a), e) / v_dot(e, e);
    t = std::min(1.0, std::max(0.0, t));
    return v_norm(v_sub(p, v_add(a, v_mul(e, t))));
}

void test_random_polylines()
{
    std::printf("test_random_polylines (3.3a):\n");
    std::mt19937_64 gen(20240607ULL);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    auto rand_dir = [&]() {
        const double z = 2.0 * U(gen) - 1.0, ph = 2.0 * std::acos(-1.0) * U(gen);
        const double rr = std::sqrt(std::max(0.0, 1.0 - z * z));
        return Vec3{rr * std::cos(ph), rr * std::sin(ph), z};
    };

    Materials mats;
    ParticleSystem sys = make_one_solid({0, 0, 0});
    const double rsum = solid_r_f_plus_r_s(sys, 0);

    const int nscenes = 1000;
    int done = 0, skipped = 0, n1 = 0, n2 = 0, n3p = 0, lost = 0, dups = 0, no_contact = 0, min_mismatch = 0;
    double worst_rel = 0.0;
    while (done < nscenes)
    {
        std::vector<Vec3> beads(1, Vec3{0, 0, 0});
        for (int k = 1; k < 6; ++k)
            beads.push_back(v_add(beads.back(), v_mul(rand_dir(), kL0 * (0.7 + 0.6 * U(gen)))));

        // частица: базовая точка на продолжении случайного сегмента (t в [-0.3,1.3]) + случайный сдвиг
        const std::size_t sg = std::min<std::size_t>(4, static_cast<std::size_t>(U(gen) * 5.0));
        const double tb = -0.3 + 1.6 * U(gen);
        const Vec3 base = v_add(beads[sg], v_mul(v_sub(beads[sg + 1], beads[sg]), tb));
        const Vec3 xs = v_add(base, v_mul(rand_dir(), rsum * (0.05 + 1.25 * U(gen))));

        double dmin = std::numeric_limits<double>::max();
        for (int k = 0; k < 5; ++k) dmin = std::min(dmin, dist_point_segment(beads[k], beads[k + 1], xs));
        if (!(dmin < 0.999 * rsum) || !(dmin > 1e-6 * rsum)) { ++skipped; continue; }

        FiberFixture fx = make_bonded_chain(beads, kL0);
        for (int d = 0; d < 3; ++d) sys.x[d][0] = xs[d];
        SceneResult r = run_scene(mats, sys, 0, fx);
        ++done;

        const std::size_t nc = r.contact_idx.size();
        if (nc == 1) ++n1; else if (nc == 2) ++n2; else if (nc >= 3) ++n3p;

        bool scene_lost = false;
        if (nc == 0) { ++no_contact; scene_lost = true; }
        else
        {
            double dmin_acc = std::numeric_limits<double>::max();
            for (std::size_t k : r.contact_idx) dmin_acc = std::min(dmin_acc, r.geoms[k].D);
            const double rel = std::fabs(dmin_acc - dmin) / dmin;
            worst_rel = std::max(worst_rel, rel);
            if (rel > 1e-12) { ++min_mismatch; scene_lost = true; }
        }
        if (scene_lost) ++lost;

        bool dup = false;
        for (std::size_t i = 0; i < r.q.size(); ++i)
            for (std::size_t j = i + 1; j < r.q.size(); ++j)
                if (v_norm(v_sub(r.q[i], r.q[j])) <= 1e-12 * kL0) dup = true;
        if (dup) ++dups;
    }
    std::printf("    статистика: сцен=%d (пропущено при генерации: %d), контактов 1/2/3+: %d/%d/%d, "
                "потерь=%d (без контакта=%d, min D не совпал=%d), дублей=%d, худшая отн. ошибка min D=%.3e\n",
                done, skipped, n1, n2, n3p, lost, no_contact, min_mismatch, dups, worst_rel);
    check(done == nscenes, "сгенерировано 1000 сцен");
    check(no_contact == 0, "во всех сценах хотя бы один контакт принят");
    check(min_mismatch == 0, "min D среди принятых контактов == глобальному минимуму расстояния до ломаной (отн. 1e-12)");
    check(lost == 0, "потерь контакта: 0");
    check(dups == 0, "дублей (две принятые q совпадают): 0");
}

// -----------------------------------------------------------------------
// 3.4) Концы нити: частица за последней бусиной (t>1 на последнем сегменте) ->
//    контакт, lambda=1; перед первой (t<0) -> контакт, lambda=0; изолированная
//    бусина -- точечный контакт без изменений.
// -----------------------------------------------------------------------
void test_fiber_ends()
{
    std::printf("test_fiber_ends (3.4):\n");
    std::vector<Vec3> beads;
    for (int k = 0; k < 4; ++k) beads.push_back({k * kL0, 0.0, 0.0});
    FiberFixture fx = make_bonded_chain(beads, kL0);
    Materials mats;
    const double h = 0.3e-3;

    { // за последней бусиной
        ParticleSystem sys = make_one_solid({3 * kL0 + 0.05 * kL0, 0.0, h});
        SceneResult r = run_scene(mats, sys, 0, fx);
        check(r.contact_idx.size() == 1, "за последней бусиной: ровно один контакт");
        if (r.contact_idx.size() == 1)
        {
            const auto& g = r.geoms[r.contact_idx[0]];
            const auto& p = r.primitives[r.contact_idx[0]];
            check(p.a == 2 && p.b == 3, "за последней бусиной: контакт у последнего сегмента (2,3)");
            check(near_abs(g.lambda, 1.0, 1e-12), "за последней бусиной: lambda == 1");
        }
        check(v_norm(sum_all_forces(r)) <= 1e-12 * v_norm(r.F_solid), "за последней бусиной: импульс сохраняется");
        check(v_norm(r.F_bead[3]) >= (1.0 - 1e-12) * v_norm(r.F_solid), "за последней бусиной: вся сила на последней бусине");
    }
    { // перед первой бусиной
        ParticleSystem sys = make_one_solid({-0.05 * kL0, 0.0, h});
        SceneResult r = run_scene(mats, sys, 0, fx);
        check(r.contact_idx.size() == 1, "перед первой бусиной: ровно один контакт");
        if (r.contact_idx.size() == 1)
        {
            const auto& g = r.geoms[r.contact_idx[0]];
            const auto& p = r.primitives[r.contact_idx[0]];
            check(p.a == 0 && p.b == 1, "перед первой бусиной: контакт у первого сегмента (0,1)");
            check(near_abs(g.lambda, 0.0, 1e-12), "перед первой бусиной: lambda == 0");
        }
        check(v_norm(sum_all_forces(r)) <= 1e-12 * v_norm(r.F_solid), "перед первой бусиной: импульс сохраняется");
        check(v_norm(r.F_bead[0]) >= (1.0 - 1e-12) * v_norm(r.F_solid), "перед первой бусиной: вся сила на первой бусине");
    }
    { // изолированная бусина: точечный контакт через реальный поиск
        FiberFixture one = make_bonded_chain({{0, 0, 0}}, kL0);
        ParticleSystem sys = make_one_solid({0.0, 0.0, h});
        SceneResult r = run_scene(mats, sys, 0, one);
        check(r.primitives.size() == 1 && !r.primitives[0].has_b, "изолированная бусина: один точечный primitive");
        check(r.contact_idx.size() == 1, "изолированная бусина: один контакт");
        const Vec3 Fp = point_contact_force_on_solid(mats, sys, one, 0, 0);
        check(v_norm(v_sub(r.F_solid, Fp)) <= 1e-12 * v_norm(Fp), "изолированная бусина: сила == силе точечного контакта");
        check(v_norm(v_add(r.F_bead[0], r.F_solid)) <= 1e-12 * v_norm(Fp), "изолированная бусина: вся сила на бусине");
    }
}

// -----------------------------------------------------------------------
// 3.5) Растянутые связи: нить из 5 бусин вдоль x растянута до 1.3*l0
//    (topology.bonds[*].l0 -- номинальная, координаты -- растянутые).
//    Контакт у дальнего конца ПОСЛЕДНЕГО сегмента (3,4): его b=4 -- последняя
//    бусина нити, у неё нет next, поэтому пара (4,s) primitive НЕ строит
//    (сегмент принадлежит первой бусине) -- контакт обязана найти пара (3,s).
//    Частица чуть дальше бусины 4 вдоль оси (t>1, lambda=1) и на высоте h.
//    Геометрия: |x_3 - x_s| = sqrt((1.3*l0+dx)^2+h^2) > r_f+r_s+l0 (старый
//    радиус), т.е. пара (3,s) в СТАРЫЙ радиус не попадала; в старый радиус
//    попадала ТОЛЬКО пара (4,s) (расстояние ~ dx,h), но она primitive не
//    строит -- контакт терялся. С радиусом по текущей длине L_max=1.3*l0 пара
//    (3,s) находится, сегмент (3,4) принимает контакт (t>=1, b_is_last ->
//    lambda=1).
// -----------------------------------------------------------------------
void test_stretched_bonds()
{
    std::printf("test_stretched_bonds (3.5):\n");
    const double L = 1.3 * kL0;
    std::vector<Vec3> beads;
    for (int k = 0; k < 5; ++k) beads.push_back({k * L, 0.0, 0.0});
    FiberFixture fx = make_bonded_chain(beads, kL0); // l0 номинальная, l_break = 100*l0

    const double dx = 0.1e-3, h = 0.3e-3;
    ParticleSystem sys = make_one_solid({4 * L + dx, 0.0, h});
    Materials mats;
    const double rsum = solid_r_f_plus_r_s(sys, 0);
    const double r_old = rsum + kL0; // r_f + r_s + l0 (l_max = номинальная l0)
    const double r_s = interactions::solid_contact_radius(sys, 0);

    rhs::CpuExecutor exec(2);
    const double r_new = interactions::fiber_solid_contact_search_radius(exec, fx.state, fx.topology, fx.materials, r_s);
    const double Lmax = interactions::fiber_max_current_bond_length(exec, fx.state, fx.topology);
    check(near_abs(Lmax, L, 1e-12 * L), "самопроверка: L_max == 1.3*l0 (текущая длина связи), получено/l0 = " + std::to_string(Lmax / kL0));

    const Vec3 xs = {sys.x[0][0], sys.x[1][0], sys.x[2][0]};
    const double d3 = v_norm(v_sub(beads[3], xs)), d4 = v_norm(v_sub(beads[4], xs));
    std::printf("    r_f+r_s=%.6e  старый радиус=%.6e  новый радиус=%.6e\n", rsum, r_old, r_new);
    std::printf("    расстояние от первой бусины сегмента (3) до частицы = %.6e (= %.4f старого радиуса)\n", d3, d3 / r_old);
    std::printf("    расстояние от последней бусины (4) до частицы = %.6e\n", d4);
    check(d3 > r_old, "самопроверка: расстояние от первой бусины сегмента до частицы > старого радиуса r_f+r_s+l0");
    check(d3 <= r_new, "самопроверка: оно же <= нового радиуса (пара (3,s) находится)");
    int within_old = 0;
    for (int k = 0; k < 5; ++k)
        if (v_norm(v_sub(beads[k], xs)) <= r_old) ++within_old;
    check(within_old == 1 && d4 <= r_old, "самопроверка: в старый радиус попала бы только последняя бусина 4 (primitive не строит)");

    SceneResult r = run_scene(mats, sys, 0, fx);
    check(r.contact_idx.size() == 1, "растянутая нить: контакт найден через build_fiber_solid_contact_search (получено " +
                                          std::to_string(r.contact_idx.size()) + ")");
    if (r.contact_idx.size() == 1)
    {
        const auto& p = r.primitives[r.contact_idx[0]];
        const auto& g = r.geoms[r.contact_idx[0]];
        check(p.a == 3 && p.b == 4 && p.b_is_last, "контакт принадлежит последнему сегменту (3,4)");
        check(near_abs(g.lambda, 1.0, 1e-12), "lambda == 1 (t>1 на последнем сегменте)");
        check(std::fabs(g.D - std::sqrt(dx * dx + h * h)) <= 1e-9 * g.D, "D == расстоянию до последней бусины");
    }
    check(v_norm(r.F_solid) > 0.0 && v_norm(sum_all_forces(r)) <= 1e-12 * v_norm(r.F_solid), "импульс сохраняется");
}

// -----------------------------------------------------------------------
// Слоты трения Fiber-Solid (fiber/fiber_solid_friction_implementation.md,
// §2, §5.1): kernel compute_fiber_solid_contact_contributions считает ТОЛЬКО
// нормальную силу и дополнительно заполняет по одному слоту
// FiberSolidFrictionContact на primitive (запись по индексу k).
//   - fiber_solid_friction_coefficient возвращает mu_fs материала нити;
//   - mu_fs = 0: слоты неактивны, has_active == false, слотов столько же, сколько
//     primitives; нормальные вклады побитово те же, что без слотов;
//   - mu_fs > 0 и контакт: активный слот с a, b, s, lambda, n, N, mu из геометрии;
//     точечный контакт: b == -1, lambda == 0; нет контакта: слот неактивен;
//   - коэффициент берётся по материалу бусины a;
//   - вектор слотов переиспользуется между шагами: размер пересчитывается,
//     устаревшие поля не остаются.
// -----------------------------------------------------------------------
void test_friction_slots()
{
    std::printf("test_friction_slots:\n");
    const double l0 = 2.0e-3;
    const double E_f = 2.0e9, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    const double h = 0.3e-3;
    Materials mats;
    rhs::CpuExecutor exec(2);

    std::vector<std::array<double, 3>> beads = {{0, 0, 0}, {l0, 0, 0}, {2 * l0, 0, 0}, {3 * l0, 0, 0}};

    // --- fiber_solid_friction_coefficient --------------------------------
    {
        fiber::FiberMaterial m = make_material(E_f, nu_f, d_f, rho_f);
        check(interactions::fiber_solid_friction_coefficient(m, 4) == 0.0, "mu_fs по умолчанию == 0 -> коэффициент 0");
        m.mu_fs = 0.37;
        m.mu_ff = 0.9;
        check(interactions::fiber_solid_friction_coefficient(m, 4) == 0.37 &&
                  interactions::fiber_solid_friction_coefficient(m, 6) == 0.37,
              "коэффициент == mu_fs материала нити (не mu_ff), не зависит от itype Solid");
    }

    // --- mu_fs = 0 против mu_fs > 0 на одном и том же контакте ---------------
    {
        FiberFixture fx0 = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);
        FiberFixture fx1 = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);
        fx1.materials[0].mu_fs = 0.4;
        ParticleSystem sys;
        sys.resize(1);
        add_solid_particle(sys, 0, 1.35 * l0, 0.2 * l0, h, 1.7e-7, 2700.0, 4);
        const std::vector<interactions::FiberSolidContactPrimitive> prims = {make_segment(1, 2, 0, 0, false)};

        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fb0, fb1, fb2;
        std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sb0, sb1, sb2;
        std::vector<interactions::FiberSolidFrictionContact> slots0, slots1;
        const bool act0 = interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx0.state,
                                                                                   fx0.materials, prims, fb0, sb0, slots0);
        const bool act1 = interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx1.state,
                                                                                   fx1.materials, prims, fb1, sb1, slots1);
        // перегрузка без слотов -- тот же нормальный вклад
        interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx1.state, fx1.materials, prims, fb2, sb2);

        check(!act0 && slots0.size() == 1 && !slots0[0].active, "mu_fs = 0: has_active == false, слот неактивен, слотов 1 (по числу primitives)");
        check(act1 && slots1.size() == 1 && slots1[0].active, "mu_fs = 0.4: has_active == true, слот активен");

        const auto geom = interactions::compute_fiber_solid_contact_geometry(mats, sys, fx1.state, fx1.materials, prims[0]);
        const auto& s = slots1[0];
        check(geom.contact && s.a == 1 && s.b == 2 && s.s == 0, "слот: a = 1, b = 2, s = 0 (индексы primitive)");
        check(s.lambda == geom.lambda && s.lambda > 0.0 && s.lambda < 1.0, "слот: lambda == geom.lambda, строго внутри (0,1)");
        check(s.n[0] == geom.n[0] && s.n[1] == geom.n[1] && s.n[2] == geom.n[2], "слот: n == geom.n (от нити к Solid)");
        check(s.N == geom.F && s.N > 0.0, "слот: N == geom.F (сила Герца, > 0)");
        check(s.mu == 0.4, "слот: mu == mu_fs материала нити");
        // n направлена от нити к Solid: (x_s - q).n > 0, где q -- точка на сегменте
        const double xs[3] = {sys.x[0][0], sys.x[1][0], sys.x[2][0]};
        double q[3];
        for (int d = 0; d < 3; ++d)
            q[d] = fx1.state.x[d][1] + s.lambda * (fx1.state.x[d][2] - fx1.state.x[d][1]);
        const double dxq[3] = {xs[0] - q[0], xs[1] - q[1], xs[2] - q[2]};
        check(dxq[0] * s.n[0] + dxq[1] * s.n[1] + dxq[2] * s.n[2] > 0.0, "n направлена от нити к Solid: (x_s - q).n > 0");

        // нормальные вклады не зависят от наличия слотов и от mu_fs
        const auto f0 = rhs::reduce_contributions(exec, fb0, fx0.state.size());
        const auto f1 = rhs::reduce_contributions(exec, fb1, fx1.state.size());
        const auto f2 = rhs::reduce_contributions(exec, fb2, fx1.state.size());
        const auto s0 = rhs::reduce_contributions(exec, sb0, sys.size());
        const auto s1 = rhs::reduce_contributions(exec, sb1, sys.size());
        const auto s2 = rhs::reduce_contributions(exec, sb2, sys.size());
        bool same = true;
        for (std::size_t i = 0; i < fx0.state.size(); ++i)
            for (std::size_t d = 0; d < dim; ++d)
                same = same && f0[i].dv[d] == f1[i].dv[d] && f1[i].dv[d] == f2[i].dv[d];
        for (std::size_t d = 0; d < dim; ++d) same = same && s0[0].dv[d] == s1[0].dv[d] && s1[0].dv[d] == s2[0].dv[d];
        check(same, "нормальные вклады (нить и Solid) побитово одинаковы при mu_fs = 0, mu_fs > 0 и в перегрузке без слотов");
    }

    // --- точечный контакт, отсутствие контакта, выравнивание слотов по индексу ---
    {
        // изолированная бусина (точечный примитив) + прямая цепочка
        FiberFixture fx = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);
        fx.materials[0].mu_fs = 0.25;
        ParticleSystem sys;
        sys.resize(3);
        add_solid_particle(sys, 0, 1.5 * l0, 0.0, h, 1.7e-7, 2700.0, 4);   // контакт с сегментом (1,2)
        add_solid_particle(sys, 1, 1.5 * l0, 0.0, 10.0 * h, 1.7e-7, 2700.0, 4); // далеко: нет контакта
        add_solid_particle(sys, 2, 0.0, 0.0, h, 1.7e-7, 2700.0, 4);        // над бусиной 0 (точечный примитив)
        std::vector<interactions::FiberSolidContactPrimitive> prims = {make_segment(1, 2, 0, 0, false),
                                                                       make_segment(1, 2, 1, 0, false), make_point(0, 2)};
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
        std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
        std::vector<interactions::FiberSolidFrictionContact> slots;
        const bool act = interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx.state, fx.materials,
                                                                                  prims, fbuf, sbuf, slots);
        check(act && slots.size() == 3, "3 primitives -> 3 слота, has_active == true");
        check(slots[0].active && !slots[1].active && slots[2].active,
              "слоты выровнены по индексу primitive: контакт / нет контакта / точечный контакт");
        check(slots[2].b == -1 && slots[2].lambda == 0.0 && slots[2].a == 0 && slots[2].s == 2,
              "точечный контакт: b == -1, lambda == 0, a = 0, s = 2");
        const auto& z = slots[1];
        check(z.a == 0 && z.b == -1 && z.s == 0 && z.lambda == 0.0 && z.N == 0.0 && z.mu == 0.0 && z.n[0] == 0.0 &&
                  z.n[1] == 0.0 && z.n[2] == 0.0,
              "слот без контакта -- value-initialized (нулевые поля, active == false)");

        // повторный вызов с другим числом primitives на том же векторе: размер и поля пересчитаны
        std::vector<interactions::FiberSolidContactPrimitive> prims2 = {make_segment(1, 2, 1, 0, false)};
        const bool act2 = interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx.state, fx.materials,
                                                                                   prims2, fbuf, sbuf, slots);
        check(!act2 && slots.size() == 1 && !slots[0].active && slots[0].N == 0.0 && slots[0].mu == 0.0,
              "повторный вызов на том же векторе: размер == 1, устаревшие поля не остались, has_active == false");
        std::vector<interactions::FiberSolidContactPrimitive> none;
        const bool act3 = interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx.state, fx.materials,
                                                                                   none, fbuf, sbuf, slots);
        check(!act3 && slots.empty(), "пустой список primitives: слотов нет, has_active == false");
    }

    // --- коэффициент берётся по материалу бусины a ------------------------
    {
        FiberFixture fx = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);
        fx.materials.push_back(make_material(E_f, nu_f, d_f, rho_f)); // материал 1
        fx.materials[0].mu_fs = 0.0;
        fx.materials[1].mu_fs = 0.3;
        ParticleSystem sys;
        sys.resize(1);
        add_solid_particle(sys, 0, 1.5 * l0, 0.0, h, 1.7e-7, 2700.0, 4);
        const std::vector<interactions::FiberSolidContactPrimitive> prims = {make_segment(1, 2, 0, 0, false)};
        for (int variant = 0; variant < 2; ++variant)
        {
            fx.state.material_id[1] = (variant == 0) ? 1 : 0; // a
            fx.state.material_id[2] = (variant == 0) ? 0 : 1; // b
            std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
            std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
            std::vector<interactions::FiberSolidFrictionContact> slots;
            const bool act = interactions::compute_fiber_solid_contact_contributions(
                exec, mats, sys, fx.state, fx.materials, prims, fbuf, sbuf, slots);
            if (variant == 0)
                check(act && slots[0].active && slots[0].mu == 0.3, "бусина a -- материал с mu_fs = 0.3: слот активен, mu == 0.3");
            else
                check(!act && !slots[0].active, "бусина a -- материал с mu_fs = 0 (у b mu_fs = 0.3): слот неактивен (берётся материал a)");
        }
    }

    // --- обёртка compute_fiber_solid_contact: флаг has_active_friction ---------
    {
        FiberFixture fx = make_chain(beads, 1e-6, E_f, nu_f, d_f, rho_f);
        fx.materials[0].mu_fs = 0.4;
        for (std::size_t i = 0; i + 1 < fx.state.size(); ++i) // связи нужны радиусу поиска (L_max)
        {
            fiber::FiberBond bond;
            bond.i = i;
            bond.j = i + 1;
            bond.l0 = l0;
            bond.l_break = l0 * 1.0e9;
            fx.topology.bonds.push_back(bond);
        }
        ParticleSystem sys;
        sys.resize(1);
        add_solid_particle(sys, 0, 1.5 * l0, 0.0, h, 1.7e-7, 2700.0, 4);
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
        std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
        std::vector<interactions::FiberSolidFrictionContact> slots;
        const auto with = interactions::compute_fiber_solid_contact(exec, mats, sys, 1, fx.state, fx.topology, fx.materials,
                                                                     fbuf, sbuf, &slots);
        const auto without = interactions::compute_fiber_solid_contact(exec, mats, sys, 1, fx.state, fx.topology,
                                                                        fx.materials, fbuf, sbuf);
        bool same = true;
        for (std::size_t i = 0; i < fx.state.size(); ++i)
            for (std::size_t d = 0; d < dim; ++d) same = same && with.fiber_delta[i].dv[d] == without.fiber_delta[i].dv[d];
        check(with.has_active_friction && !slots.empty() && !without.has_active_friction,
              "compute_fiber_solid_contact: has_active_friction == true со слотами, false без них");
        check(same, "compute_fiber_solid_contact: нормальный вклад одинаков со слотами и без");
    }
}

// -----------------------------------------------------------------------
// Закон 0: предел шага из ядра (FiberContactStepParams::law0_dt_min_out) побитово совпадает с
// compute_fiber_solid_contact_dt_limit при 1, 2, 7, 8 потоках; сцена без контактов -> max().
// -----------------------------------------------------------------------
bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

void test_law0_dt_from_kernel()
{
    std::printf("test_law0_dt_from_kernel:\n");
    Materials mats;
    std::vector<Vec3> beads;
    for (int k = 0; k < 24; ++k) beads.push_back({k * kL0, 0.0, 0.0});
    FiberFixture fx = make_bonded_chain(beads, kL0);

    for (int scene = 0; scene < 2; ++scene)
    {
        const std::size_t np = 20;
        ParticleSystem sys;
        sys.resize(np);
        for (std::size_t i = 0; i < np; ++i)
        {
            // scene 0: частицы вдоль нити на разной высоте (часть в контакте); scene 1: все далеко.
            const double zoff = (scene == 0) ? (0.2 + 0.05 * static_cast<double>(i % 7)) * 1.0e-3 : 5.0e-2;
            add_solid_particle(sys, i, (0.7 + 2.1 * static_cast<double>(i)) * 1.0e-3, 0.0, zoff, kMassS, kRhoS,
                               kItypeS);
        }
        for (int nt : {1, 2, 7, 8})
        {
            rhs::CpuExecutor exec(nt);
            interactions::FiberSolidContactSearchResult search;
            interactions::build_fiber_solid_contact_search(exec, mats, sys, sys.size(), fx.state, fx.topology,
                                                            fx.materials, search);
            std::size_t ncontact = 0;
            for (const auto& pr : search.primitives)
                if (interactions::compute_fiber_solid_contact_geometry(mats, sys, fx.state, fx.materials, pr).contact)
                    ++ncontact;

            interactions::FiberContactStepParams p;
            double dt_kernel = -1.0;
            p.law0_dt_min_out = &dt_kernel;
            std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> fbuf;
            std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sbuf;
            std::vector<interactions::FiberSolidFrictionContact> slots;
            interactions::compute_fiber_solid_contact_contributions(exec, mats, sys, fx.state, fx.materials,
                                                                     search.primitives, fbuf, sbuf, slots, p);
            const double dt_ref = interactions::compute_fiber_solid_contact_dt_limit(
                exec, mats, sys, fx.state, fx.materials, search.primitives);
            const std::string lab = "сцена " + std::to_string(scene) + ", потоков " + std::to_string(nt);
            check(same_bits(dt_kernel, dt_ref), lab + ": dt ядра == compute_fiber_solid_contact_dt_limit (побитово)");
            if (scene == 0) check(ncontact > 0, lab + ": в сцене есть контакты (предпосылка)");
            else check(ncontact == 0 && dt_kernel == std::numeric_limits<double>::max(), lab + ": без контактов -> max()");
        }
    }
}

} // namespace

int main()
{
    test_lambda_and_clamp();
    test_end_of_fiber_point_contact();
    test_node_no_double_force();
    test_momentum_conservation();
    test_elastic_by_itype_reference_values();
    test_kink_vertex();
    test_concave_kink();
    test_node_rounding();
    test_random_polylines();
    test_fiber_ends();
    test_stretched_bonds();
    test_friction_slots();
    test_law0_dt_from_kernel();

    if (failures == 0)
        std::printf("test_fiber_solid: OK\n");
    else
        std::printf("test_fiber_solid: %d failure(s)\n", failures);

    return failures;
}
