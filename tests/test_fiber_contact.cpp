// tests/test_fiber_contact.cpp
// Unit-тесты Fiber-Fiber-Contact (fiber/PLAN.md, этап Ф5; fiber/SPH.md,
// раздел «Fiber-Fiber-Contact»):
//   - перпендикулярные нити: D, n, P_ij, w_ij против ручного расчёта;
//   - почти параллельные нити: нет NaN/Inf, срабатывает переход на Surface-End;
//   - Surface-End и End-Surface, End-End: геометрия и сила;
//   - сумма импульсов по всем вкладам одного primitive равна нулю;
//   - диагностика: число пар на одно пересечение.
//
// Стиль -- по образцу tests/test_fiber_bond.cpp (тот же формат check()/
// near_rel()/near_abs(), тот же CpuExecutor/ContributionBuffer напрямую).
// Собирается с SPH_DIM=3 (CMakeLists.txt).
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_contact.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_contact требует сборки с SPH_DIM=3");

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
bool finite3(const double v[3]) { return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]); }

fiber::FiberMaterial make_material(double E, double nu, double d, double rho, double eps_break = 1.0)
{
    fiber::FiberMaterial m;
    m.E = E; m.nu = nu; m.d = d; m.rho = rho; m.eps_break = eps_break;
    return m;
}

// --- Двухнитевая фикстура: нить A (object_id=0), нить B (object_id=1) ------
// beadsA/beadsB -- точки в порядке следования по нити (prev/next строятся
// цепочкой), massesA/massesB -- массы бусин (по одной на точку). Обе нити
// используют ОДИН материал (materials[0]) -- для тестов с разными
// материалами используется отдельный конструктор ниже.
struct TwoFiberFixture
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

void append_chain(fiber::FiberState& state, fiber::FiberTopology& topology, std::size_t base, int object_id,
                   int material_id, const std::vector<std::array<double, 3>>& beads,
                   const std::vector<double>& masses)
{
    const std::size_t n = beads.size();
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t gi = base + k;
        for (int d3 = 0; d3 < 3; ++d3) state.x[d3][gi] = beads[k][d3];
        for (std::size_t dd = 0; dd < dim; ++dd) state.vx[dd][gi] = 0.0;
        state.mass[gi] = masses[k];
        state.object_id[gi] = object_id;
        state.material_id[gi] = material_id;
        if (k > 0) topology.prev[gi] = static_cast<long long>(gi - 1);
        if (k + 1 < n)
        {
            topology.next[gi] = static_cast<long long>(gi + 1);
            fiber::FiberBond b;
            b.i = gi;
            b.j = gi + 1;
            double l2 = 0.0;
            for (int d3 = 0; d3 < 3; ++d3)
            {
                const double dx = beads[k + 1][d3] - beads[k][d3];
                l2 += dx * dx;
            }
            b.l0 = std::sqrt(l2);
            b.l_break = b.l0 * 1e9; // разрыв не в предмете Ф5 (Ф7)
            topology.bonds.push_back(b);
        }
    }
}

TwoFiberFixture make_two_fiber_fixture(const std::vector<std::array<double, 3>>& beadsA,
                                        const std::vector<double>& massesA,
                                        const std::vector<std::array<double, 3>>& beadsB,
                                        const std::vector<double>& massesB, double E, double nu, double d,
                                        double rho)
{
    TwoFiberFixture fx;
    fx.materials.push_back(make_material(E, nu, d, rho));
    const std::size_t nA = beadsA.size(), nB = beadsB.size();
    fx.state.resize(nA + nB);
    fx.topology.resize(nA + nB);
    append_chain(fx.state, fx.topology, 0, 0, 0, beadsA, massesA);
    append_chain(fx.state, fx.topology, nA, 1, 0, beadsB, massesB);
    return fx;
}

// Прямая точечная проверка геометрии для заданного (заранее известного) primitive,
// без поиска -- вызывающая сторона задаёт индексы/case вручную (тесты 1-3).
interactions::FiberContactPrimitive make_primitive(std::size_t i, long long ip, long long in, std::size_t j,
                                                     long long jp, long long jn, interactions::FiberContactCase c)
{
    interactions::FiberContactPrimitive prim;
    prim.i = i; prim.i_prev = ip; prim.i_next = in;
    prim.j = j; prim.j_prev = jp; prim.j_next = jn;
    prim.case_ = c;
    return prim;
}

// -----------------------------------------------------------------------
// 1) Перпендикулярные нити (Surface-Surface, нетривиальные P_ij/P_ji).
//
// Нить A вдоль оси x (через начало координат), нить B вдоль оси y, со
// смещением (x0,y0,h): пересечение осей-прямых -- в точке (x0,0,0) на оси A
// / (x0,0,h) на оси B (расстояние по перпендикуляру -- h, не зависит от
// x0,y0 -- см. вывод в шапке include/interactions/fiber_contact.hpp).
//
// Ручной расчёт (p_A=(-1,0,0) "назад" к A0, p_B=(0,-1,0) "назад" к B0):
//   r = xA1-xB1 = (-x0,-y0,-h),  c = p_A.p_B = 0,  denom = 1
//   P_ij = -dot(r,p_A) = -x0,   P_ji = dot(r,p_B) = y0
//   Q_i = (x0,0,0), Q_j = (x0,0,h), diff = (0,0,-h), D = h, n = (0,0,-1)
//   d_n_i = |A1-A2| = l0 -> P_ij<0 -> w_i = (l0+P_ij)/l0 = (l0-x0)/l0
//   d_p_j = |B1-B0| = l0 -> P_ji>0 -> w_j = (l0-P_ji)/l0 = (l0-y0)/l0
//   сосед i -- A2 (next, т.к. P_ij<=0), сосед j -- B0 (prev, т.к. P_ji>0)
// -----------------------------------------------------------------------
void test_perpendicular_manual()
{
    std::printf("test_perpendicular_manual:\n");
    const double l0 = 0.01;
    const double x0 = 0.4 * l0, y0 = 0.3 * l0, h = 0.5e-3;
    const double E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;

    std::vector<std::array<double, 3>> beadsA = {{-l0, 0, 0}, {0, 0, 0}, {l0, 0, 0}};
    std::vector<std::array<double, 3>> beadsB = {{x0, y0 - l0, h}, {x0, y0, h}, {x0, y0 + l0, h}};
    std::vector<double> mA = {1e-6, 1.1e-6, 1.2e-6};
    std::vector<double> mB = {1.3e-6, 1.4e-6, 1.5e-6};

    TwoFiberFixture fx = make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
    // A1 -- индекс 1, B1 -- индекс 3+1=4
    interactions::FiberContactPrimitive prim =
        make_primitive(1, 0, 2, 4, 3, 5, interactions::FiberContactCase::SurfaceSurface);

    interactions::FiberContactGeometry geom = interactions::compute_fiber_contact_geometry(fx.state, fx.materials, prim);

    check(near_rel(geom.D, h, 1e-9), "D == h (перпендикулярное расстояние осей)");
    check(near_abs(geom.n[0], 0.0, 1e-9) && near_abs(geom.n[1], 0.0, 1e-9) && near_abs(geom.n[2], -1.0, 1e-9),
          "n == (0,0,-1)");
    check(near_rel(geom.P_ij, -x0, 1e-9), "P_ij == -x0");
    check(near_rel(geom.P_ji, y0, 1e-9), "P_ji == y0");

    const double w_i_expected = (l0 - x0) / l0;
    const double w_j_expected = (l0 - y0) / l0;
    check(near_rel(geom.w_i, w_i_expected, 1e-9), "w_i совпадает с ручной формулой (P_ij<=0 -> ветка d_n)");
    check(near_rel(geom.w_j, w_j_expected, 1e-9), "w_j совпадает с ручной формулой (P_ji>0 -> ветка d_p)");
    check(geom.i_has_neighbor && geom.i_neighbor == 2, "сосед i -- A2 (next, P_ij<=0)");
    check(geom.j_has_neighbor && geom.j_neighbor == 3, "сосед j -- B0 (prev, P_ji>0)");

    // Сила Герца: одинаковые материалы -> самопроверка против формулы Meyer.
    const double d_eff = d; // r_i+r_j = d/2+d/2 = d (одинаковые нити)
    check(geom.contact, "D < d_eff -- контакт есть");
    const double delta_expected = d_eff - h;
    const double E_star = E / (2.0 * (1.0 - nu * nu));
    const double R_star = std::sqrt(d / 4.0);
    const double F_expected = (4.0 / 3.0) * E_star * R_star * std::pow(delta_expected, 1.5);
    check(near_rel(geom.F, F_expected, 1e-9), "F совпадает с формулой Meyer (одинаковые нити)");
}

// -----------------------------------------------------------------------
// 2) Почти параллельные нити: |p_i x p_j| < eps -> переход на Surface-End.
// Нить B наклонена относительно A на угол theta=1e-5 рад (sin(theta) << eps
// =1e-3) -- геометрия должна СОВПАСТЬ с прямым ручным расчётом по формуле
// Surface-End (i=A -- поверхность, j=B1 -- точка), а не с (в общем случае
// отличающимся) результатом полной системы Surface-Surface. Проверяется
// также отсутствие NaN/Inf.
// -----------------------------------------------------------------------
void test_near_parallel_fallback()
{
    std::printf("test_near_parallel_fallback:\n");
    const double l0 = 0.01;
    const double x0 = 0.4 * l0, h = 0.5e-3;
    const double theta = 1.0e-5; // << kFiberContactParallelEps = 1e-3
    const double E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;

    const double c = std::cos(theta), s = std::sin(theta);
    std::array<double, 3> B1 = {x0, 0.0, h};
    std::array<double, 3> dirB = {c, s, 0.0}; // почти вдоль +y, чуть наклонена
    std::array<double, 3> B0 = {B1[0] - l0 * dirB[0], B1[1] - l0 * dirB[1], B1[2] - l0 * dirB[2]};
    std::array<double, 3> B2 = {B1[0] + l0 * dirB[0], B1[1] + l0 * dirB[1], B1[2] + l0 * dirB[2]};

    std::vector<std::array<double, 3>> beadsA = {{-l0, 0, 0}, {0, 0, 0}, {l0, 0, 0}};
    std::vector<std::array<double, 3>> beadsB = {B0, B1, B2};
    std::vector<double> mA = {1e-6, 1e-6, 1e-6};
    std::vector<double> mB = {1e-6, 1e-6, 1e-6};

    TwoFiberFixture fx = make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
    interactions::FiberContactPrimitive prim =
        make_primitive(1, 0, 2, 4, 3, 5, interactions::FiberContactCase::SurfaceSurface);
    interactions::FiberContactGeometry geom = interactions::compute_fiber_contact_geometry(fx.state, fx.materials, prim);

    check(std::isfinite(geom.D) && finite3(geom.n) && std::isfinite(geom.F) && std::isfinite(geom.w_i) &&
              std::isfinite(geom.w_j),
          "нет NaN/Inf в геометрии почти параллельных нитей");

    // Ручной расчёт по формуле Surface-End (i=A -- поверхность):
    // p_A = (-1,0,0) (как в тесте 1); P_ij = dot(p_A, xB1-xA1) = -x0;
    // Q_i = (x0,0,0); diff = Q_i - xB1 = (0,0,-h); D=h, n=(0,0,-1).
    check(near_rel(geom.D, h, 1e-6), "D совпадает с формулой Surface-End (не с вырождающейся Surface-Surface)");
    check(near_abs(geom.n[0], 0.0, 1e-6) && near_abs(geom.n[1], 0.0, 1e-6) && near_abs(geom.n[2], -1.0, 1e-6),
          "n совпадает с формулой Surface-End");
    check(near_rel(geom.P_ij, -x0, 1e-6), "P_ij совпадает с формулой Surface-End");
    check(geom.w_j == 1.0, "j (в вырожденной ветке) получает w_j=1 -- без распределения на соседей");
    check(!geom.j_has_neighbor, "j в вырожденной ветке не имеет соседа для распределения");

    // Сумма импульсов ИМЕННО в вырожденной ветке (не только в общем случае
    // Surface-Surface, test_momentum_conservation ниже): w_j=1 означает, что
    // на стороне j всей силой -F*n наделяется одна бусина j (без соседа), а
    // на стороне i -- w_i*F*n самой i плюс (1-w_i)*F*n соседу i_next; сумма
    // обязана быть равна нулю и в этом асимметричном случае (координатор,
    // пункт 3: явно проверить сохранение импульса именно в вырожденной
    // ветке, а не полагаться на то, что near-parallel тест проверяет только
    // геометрию).
    {
        rhs::CpuExecutor exec(2);
        std::vector<interactions::FiberContactPrimitive> primitives = {prim};
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
        interactions::compute_fiber_contact_contributions(exec, fx.state, fx.materials, primitives, buffers);
        std::vector<fiber::FiberDelta> out = rhs::reduce_contributions(exec, buffers, fx.state.size());

        double total[dim] = {};
        double scale = 0.0;
        for (std::size_t idx = 0; idx < fx.state.size(); ++idx)
            for (std::size_t dd = 0; dd < dim; ++dd)
            {
                const double p = fx.state.mass[idx] * out[idx].dv[dd];
                total[dd] += p;
                scale = std::max(scale, std::fabs(p));
            }
        if (scale < 1e-30) scale = 1.0;
        for (std::size_t dd = 0; dd < dim; ++dd)
            check(std::fabs(total[dd]) < 1e-9 * scale,
                  "вырожденная ветка (near-parallel): суммарный импульс по компоненте " + std::to_string(dd) +
                      " равен нулю: " + std::to_string(total[dd]));

        // Асимметрия упрощения: сторона j (w_j=1) отдаёт вклад ровно одной
        // бусине, сторона i (Surface-End-формула) -- двум (i и сосед) --
        // ожидаемо 3 бусины с ненулевым вкладом, а не 4, как в общем случае
        // Surface-Surface (test_momentum_conservation).
        int nonzero = 0;
        for (std::size_t idx = 0; idx < fx.state.size(); ++idx)
        {
            bool any = false;
            for (std::size_t dd = 0; dd < dim; ++dd) any = any || (out[idx].dv[dd] != 0.0);
            if (any) ++nonzero;
        }
        check(nonzero == 3, "вырожденная ветка: 3 бусины получили вклад (i, сосед i, j -- без соседа j)");
    }
}

// -----------------------------------------------------------------------
// 3) Surface-End / End-Surface / End-End -- геометрия и сила.
// Переиспользует ту же сцену теста 1 (перпендикулярные оси, x0,h), но j (или
// i) -- концевая бусина двухбусинной нити (одна связь, один сосед).
// -----------------------------------------------------------------------
void test_surface_end_and_end_end()
{
    std::printf("test_surface_end_and_end_end:\n");
    const double l0 = 0.01;
    const double x0 = 0.4 * l0, h = 0.5e-3;
    const double E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;

    // --- Surface-End: A -- 3 бусины (A1 внутренняя), B -- 2 бусины (B1 конец) --
    {
        std::vector<std::array<double, 3>> beadsA = {{-l0, 0, 0}, {0, 0, 0}, {l0, 0, 0}};
        std::vector<std::array<double, 3>> beadsB = {{x0, -l0, h}, {x0, 0, h}}; // B1 = индекс 1 -- конец (next=-1)
        std::vector<double> mA = {1e-6, 1e-6, 1e-6};
        std::vector<double> mB = {1e-6, 1e-6};
        TwoFiberFixture fx = make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
        // A1 = индекс 1 (prev=0,next=2); B1 = индекс 3+1=4 (prev=3,next=-1)
        interactions::FiberContactPrimitive prim =
            make_primitive(1, 0, 2, 4, 3, -1, interactions::FiberContactCase::SurfaceEnd);
        interactions::FiberContactGeometry geom =
            interactions::compute_fiber_contact_geometry(fx.state, fx.materials, prim);

        check(near_rel(geom.D, h, 1e-9), "Surface-End: D == h");
        check(near_abs(geom.n[2], -1.0, 1e-9), "Surface-End: n вдоль -z");
        check(near_rel(geom.P_ij, -x0, 1e-9), "Surface-End: P_ij == -x0 (та же формула, что и деградировавшая Surface-Surface)");
        check(geom.w_j == 1.0 && !geom.j_has_neighbor, "Surface-End: j (конец) получает w_j=1, соседа нет");
        check(geom.i_has_neighbor && geom.i_neighbor == 2, "Surface-End: сосед i -- A2 (next, P_ij<=0)");
    }

    // --- End-Surface: симметрично, i -- конец, j -- внутренняя -------------
    {
        std::vector<std::array<double, 3>> beadsA = {{x0, -l0, h}, {x0, 0, h}}; // A1 = индекс1 -- конец
        std::vector<std::array<double, 3>> beadsB = {{-l0, 0, 0}, {0, 0, 0}, {l0, 0, 0}};
        std::vector<double> mA = {1e-6, 1e-6};
        std::vector<double> mB = {1e-6, 1e-6, 1e-6};
        TwoFiberFixture fx = make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
        // A1 = индекс 1 (prev=0,next=-1); B1 = индекс 2+1=3 (prev=2,next=4)
        interactions::FiberContactPrimitive prim =
            make_primitive(1, 0, -1, 3, 2, 4, interactions::FiberContactCase::EndSurface);
        interactions::FiberContactGeometry geom =
            interactions::compute_fiber_contact_geometry(fx.state, fx.materials, prim);

        check(near_rel(geom.D, h, 1e-9), "End-Surface: D == h");
        check(geom.w_i == 1.0 && !geom.i_has_neighbor, "End-Surface: i (конец) получает w_i=1, соседа нет");
        check(geom.j_has_neighbor && geom.j_neighbor == 4, "End-Surface: сосед j -- B2 (next, P_ji<=0)");
    }

    // --- End-End: обе бусины -- концы (двухбусинные нити) -------------------
    {
        const double gap = 0.6e-3; // < d -- гарантирован контакт
        std::vector<std::array<double, 3>> beadsA = {{0, 0, 0}, {0, 0, l0}};       // A1=(0,0,l0) -- конец
        std::vector<std::array<double, 3>> beadsB = {{0, 0, l0 + gap}, {0, 0, l0 + gap + l0}}; // B0 -- конец, соседний с A1
        std::vector<double> mA = {1e-6, 1e-6};
        std::vector<double> mB = {1e-6, 1e-6};
        TwoFiberFixture fx = make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
        // A1 = индекс1 (prev=0,next=-1); B0 = индекс2 (prev=-1,next=3)
        interactions::FiberContactPrimitive prim =
            make_primitive(1, 0, -1, 2, -1, 3, interactions::FiberContactCase::EndEnd);
        interactions::FiberContactGeometry geom =
            interactions::compute_fiber_contact_geometry(fx.state, fx.materials, prim);

        check(near_rel(geom.D, gap, 1e-9), "End-End: D == |xi-xj|");
        // n указывает "от j к i" (репульсивная конвенция, см. шапку файла,
        // пункт 1): j=B0 находится в направлении +z от i=A1, поэтому сила на
        // i (вдоль +n) должна толкать i в СТОРОНУ -z, прочь от j.
        check(near_abs(geom.n[2], -1.0, 1e-9) && near_abs(geom.n[0], 0.0, 1e-9) && near_abs(geom.n[1], 0.0, 1e-9),
              "End-End: n вдоль -z (от B к A, репульсивно для i)");
        check(geom.w_i == 1.0 && geom.w_j == 1.0, "End-End: оба веса == 1 (нет распределения)");
        check(!geom.i_has_neighbor && !geom.j_has_neighbor, "End-End: соседей нет ни у одной стороны");
        check(geom.contact && geom.delta > 0.0, "End-End: контакт есть (gap < d)");
    }
}

// -----------------------------------------------------------------------
// 4) Сумма импульсов по всем вкладам ОДНОГО primitive равна нулю (третий
// закон Ньютона на паре нитей, PLAN.md §3.2: "нить i получает Fn, нить j --
// -Fn"). Сцена -- перпендикулярные нити из теста 1 (нетривиальные w_i,w_j,
// оба соседа задействованы -- 4 частицы получают вклад).
//
// Primitive строится ВРУЧНУЮ (не через build_fiber_contact_search): при
// l0 (шаг бусин вдоль нити) того же порядка, что радиус поиска d+l_max,
// поиск у пересечения естественно находит НЕСКОЛЬКО spatial-пар одновременно
// (см. test_diagnostics_pairs_per_intersection и шапку fiber_contact.hpp) --
// это документированное свойство модели, но здесь предмет проверки -- именно
// ОДИН primitive (буквальная формулировка PLAN.md), поэтому его нужно
// зафиксировать явно, а не полагаться на то, сколько пар вернёт поиск.
// -----------------------------------------------------------------------
void test_momentum_conservation()
{
    std::printf("test_momentum_conservation:\n");
    const double l0 = 0.01;
    const double x0 = 0.4 * l0, y0 = 0.3 * l0, h = 0.5e-3;
    const double E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;

    std::vector<std::array<double, 3>> beadsA = {{-l0, 0, 0}, {0, 0, 0}, {l0, 0, 0}};
    std::vector<std::array<double, 3>> beadsB = {{x0, y0 - l0, h}, {x0, y0, h}, {x0, y0 + l0, h}};
    std::vector<double> mA = {1.0e-6, 1.3e-6, 0.8e-6}; // разные массы -- проверка не на F, а на m*a
    std::vector<double> mB = {1.7e-6, 0.9e-6, 1.5e-6};

    TwoFiberFixture fx = make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);

    rhs::CpuExecutor exec(2);
    std::vector<interactions::FiberContactPrimitive> primitives = {
        make_primitive(1, 0, 2, 4, 3, 5, interactions::FiberContactCase::SurfaceSurface)};

    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
    interactions::compute_fiber_contact_contributions(exec, fx.state, fx.materials, primitives, buffers);
    std::vector<fiber::FiberDelta> out = rhs::reduce_contributions(exec, buffers, fx.state.size());

    double total[dim] = {};
    for (std::size_t i = 0; i < fx.state.size(); ++i)
        for (std::size_t dd = 0; dd < dim; ++dd) total[dd] += fx.state.mass[i] * out[i].dv[dd];

    double scale = 0.0;
    for (std::size_t i = 0; i < fx.state.size(); ++i)
        for (std::size_t dd = 0; dd < dim; ++dd) scale = std::max(scale, std::fabs(fx.state.mass[i] * out[i].dv[dd]));
    if (scale < 1e-30) scale = 1.0;

    for (std::size_t dd = 0; dd < dim; ++dd)
        check(std::fabs(total[dd]) < 1e-9 * scale,
              "суммарный импульс по компоненте " + std::to_string(dd) + " равен нулю: " + std::to_string(total[dd]));

    // Дополнительно: ровно 4 бусины получили ненулевой вклад (i, i_next=A2, j, j_prev=B0).
    int nonzero = 0;
    for (std::size_t i = 0; i < fx.state.size(); ++i)
    {
        bool any = false;
        for (std::size_t dd = 0; dd < dim; ++dd) any = any || (out[i].dv[dd] != 0.0);
        if (any) ++nonzero;
    }
    check(nonzero == 4, "ровно 4 бусины получили вклад (i, сосед i, j, сосед j)");
}

// -----------------------------------------------------------------------
// 5) Диагностика: число пар на одно пересечение. Две нити пересекаются
// (перпендикулярно, с проникновением) -- рядом с точкой пересечения
// несколько соседних бусин каждой нити одновременно попадают в радиус
// поиска d+l_max друг у друга, поэтому contact_count в норме БОЛЬШЕ 1 на
// одно геометрическое пересечение, но остаётся вменяемым (небольшой
// множитель), не абсурдно большим.
// -----------------------------------------------------------------------
void test_diagnostics_pairs_per_intersection()
{
    std::printf("test_diagnostics_pairs_per_intersection:\n");
    const double l0 = 1.0e-3; // шаг бусин -- того же порядка, что d -- соседние бусины физически близко
    const double d = 1.0e-3, E = 2.0e9, nu = 0.3, rho = 1200.0;
    const double h = 0.3e-3; // проникновение в точке пересечения

    // Нить A вдоль x, 7 бусин, нить B вдоль y, 7 бусин, пересекаются в (0,0,*).
    std::vector<std::array<double, 3>> beadsA, beadsB;
    std::vector<double> mA, mB;
    const int N = 7;
    for (int k = -N / 2; k <= N / 2; ++k)
    {
        beadsA.push_back({static_cast<double>(k) * l0, 0.0, 0.0});
        beadsB.push_back({0.0, static_cast<double>(k) * l0, h});
        mA.push_back(1e-6);
        mB.push_back(1e-6);
    }

    TwoFiberFixture fx = make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
    rhs::CpuExecutor exec(2);
    interactions::FiberContactSearchResult search;
    interactions::build_fiber_contact_search(exec, fx.state, fx.topology, fx.materials, search);
    interactions::FiberContactDiagnostics diag =
        interactions::compute_fiber_contact_diagnostics(exec, fx.state, fx.materials, search.primitives);

    std::printf("    pair_count=%zu contact_count=%zu max_penetration=%.3e\n", diag.pair_count, diag.contact_count,
                diag.max_penetration);

    check(diag.contact_count > 0, "диагностика: есть хотя бы один реальный контакт");
    check(diag.contact_count < static_cast<std::size_t>(N) * static_cast<std::size_t>(N),
          "диагностика: contact_count не абсурдно велик (< N*N всех возможных пар бусин)");
    // "Одно пересечение" здесь -- одно геометрическое место (окрестность (0,0,h));
    // разумный множитель "пар на пересечение" -- не более desятка бусин с каждой
    // стороны, задействованных в радиусе поиска.
    check(diag.contact_count <= 20, "диагностика: разумное число контактных пар на одно пересечение (<=20)");
    check(diag.pair_count >= diag.contact_count, "pair_count (широкая фаза) >= contact_count (реальные контакты)");
}

} // namespace

int main()
{
    test_perpendicular_manual();
    test_near_parallel_fallback();
    test_surface_end_and_end_end();
    test_momentum_conservation();
    test_diagnostics_pairs_per_intersection();

    if (failures == 0)
        std::printf("test_fiber_contact: OK\n");
    else
        std::printf("test_fiber_contact: %d failure(s)\n", failures);

    return failures;
}
