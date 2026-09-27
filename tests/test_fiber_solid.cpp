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
#include <limits>
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

} // namespace

int main()
{
    test_lambda_and_clamp();
    test_end_of_fiber_point_contact();
    test_node_no_double_force();
    test_momentum_conservation();
    test_elastic_by_itype_reference_values();

    if (failures == 0)
        std::printf("test_fiber_solid: OK\n");
    else
        std::printf("test_fiber_solid: %d failure(s)\n", failures);

    return failures;
}
