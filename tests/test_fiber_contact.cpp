// tests/test_fiber_contact.cpp
// Unit-тесты Fiber-Fiber-Contact в схеме «пара сегментов -- одна сила»
// (fiber/fiber_contact_segment_pairs.md, §8): primitive -- пара отрезков осей
// двух нитей, ближайшие точки ищутся с ограничением по отрезкам, правило
// владения гарантирует ровно один контакт на физическое пересечение (кроме
// изогнутых нитей с несколькими локальными минимумами -- §3, §8.3).
//
// Разделы теста соответствуют §8.1-8.8 спецификации:
//   8.1 -- главный тест: ровно один контакт на пересечение прямых нитей,
//          сила равна силе Герца одной пары, для сетки углов x положений;
//   8.2 -- концы нитей (бывшие Surface-End/End-Surface/End-End) и
//          изолированная бусина;
//   8.3 -- изогнутая нить: выпуклый узел (один контакт, s=0) и вогнутое
//          «седло» (два контакта, разные q_A);
//   8.4 -- сохранение импульса на сценах 8.1-8.3;
//   8.5 -- полнота поиска при растянутых связях (радиус по ТЕКУЩЕЙ, не
//          номинальной, длине);
//   8.6 -- почти параллельные и совпадающие оси (нет NaN, импульс ноль);
//   8.7 -- compute_fiber_contact_dt_limit совпадает с ручной формулой (§6);
//   8.8 -- смешение mu_ff и нормальная часть при mu_ff=0 (само трение --
//          отдельная стадия, interactions/fiber_friction.hpp, tests/test_fiber_friction.cpp).
//
// Стиль -- по образцу старого test_fiber_contact.cpp / test_fiber_bond.cpp
// (тот же check()/near_rel()/near_abs(), тот же CpuExecutor/ContributionBuffer
// напрямую). Собирается с SPH_DIM=3 (CMakeLists.txt).
//
// Все проверки сил -- на СУММАРНЫХ вкладах после reduce_contributions по ВСЕМ
// primitives сцены (build_fiber_contact_search), а не по одному вручную
// собранному primitive -- именно в этом месте старая Meyer-подобная схема
// ошибалась (fiber_contact_segment_pairs.md, §8, шапка).
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_contact.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <random>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_contact требует сборки с SPH_DIM=3");

namespace
{

const double kPi = 3.14159265358979323846;

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
    m.E = E;
    m.nu = nu;
    m.d = d;
    m.rho = rho;
    m.eps_break = eps_break;
    return m;
}

// -----------------------------------------------------------------------
// hertz_force -- ТА ЖЕ формула, что kernel (compute_fiber_contact_geometry,
// §4 спецификации): 1/E* = (1-nu_a^2)/E_a+(1-nu_b^2)/E_b, R*=ra*rb/(ra+rb),
// F=(4/3)*E**sqrt(R*)*delta^1.5. Используется, чтобы сравнивать суммарный
// вклад ВСЕХ primitives сцены с "силой одной пары" (§8, шапка).
// -----------------------------------------------------------------------
double hertz_force(double delta, const fiber::FiberMaterial& mat_a, const fiber::FiberMaterial& mat_b)
{
    if (!(delta > 0.0)) return 0.0;
    const double ra = 0.5 * mat_a.d;
    const double rb = 0.5 * mat_b.d;
    const double inv_Estar = (mat_a.E > 0.0 ? (1.0 - mat_a.nu * mat_a.nu) / mat_a.E : 0.0) +
                              (mat_b.E > 0.0 ? (1.0 - mat_b.nu * mat_b.nu) / mat_b.E : 0.0);
    if (!(inv_Estar > 0.0)) return 0.0;
    const double E_star = 1.0 / inv_Estar;
    const double R_star = (ra > 0.0 && rb > 0.0) ? (ra * rb) / (ra + rb) : 0.0;
    return (4.0 / 3.0) * E_star * std::sqrt(R_star) * std::pow(delta, 1.5);
}

// --- Двухнитевая фикстура: нить A (object_id=0, бусины [0,n_a)), нить B
// (object_id=1, бусины [n_a,n_a+n_b)) -- как в старом файле, плюс n_a/n_b
// сохранены явно (нужны для суммирования силы ПО НИТИ, §8: "суммарная сила
// на нить A по всем бусинам"). Обе нити используют ОДИН материал
// (materials[0]).
// -----------------------------------------------------------------------
struct TwoFiberFixture
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    std::size_t n_a = 0;
    std::size_t n_b = 0;
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
            b.l_break = b.l0 * 1e9; // разрыв не в предмете этого теста
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
    fx.n_a = nA;
    fx.n_b = nB;
    return fx;
}

// Задаёт скорость одной бусины (нужна тестам трения §8.8 -- ненулевые, РАЗНЫЕ
// скорости бусин, чего не делают тесты геометрии/нормальной силы).
void set_velocity(fiber::FiberState& state, std::size_t idx, double vx, double vy, double vz)
{
    state.vx[0][idx] = vx;
    state.vx[1][idx] = vy;
    state.vx[2][idx] = vz;
}

// -----------------------------------------------------------------------
// build_straight_pair -- строит ДВЕ ПРЯМЫЕ нити (произвольное число бусин,
// направления eA_hat/eB_hat, шаг step) так, чтобы ближайшие точки
// БЕСКОНЕЧНЫХ прямых, содержащих их оси, находились РОВНО на расстоянии
// D_true, в точках с параметрами u0 (вдоль A от бусины index a_origin) и v0
// (вдоль B от бусины index b_origin). Аналитическая конструкция: n_hat =
// normalized(eA_hat x eB_hat) (перпендикуляр к обеим осям), qA = A_ref +
// u0*eA_hat, qB = qA - D_true*n_hat (лежит на прямой B по построению, т.к.
// n_hat перпендикулярен eB_hat), B_ref = qB - v0*eB_hat. Тогда qA,qB --
// ТОЧНЫЙ (аналитический, не приближённый) минимум расстояния между
// бесконечными прямыми: вектор qA-qB перпендикулярен ОБЕИМ осям, что и есть
// условие минимума для скрещивающихся прямых. Даёт относительную точность
// построения сцены, ограниченную только double (существенно лучше 1e-12,
// нужных §8.1).
// -----------------------------------------------------------------------
struct StraightFiberPair
{
    std::vector<std::array<double, 3>> beadsA, beadsB;
    double n_hat[3] = {}; // от B к A (см. ниже)
};

StraightFiberPair build_straight_pair(int nA, int a_origin, double stepA, int nB, int b_origin, double stepB,
                                       const std::array<double, 3>& eA_hat, const std::array<double, 3>& eB_hat,
                                       double D_true, double u0, double v0)
{
    auto cross = [](const std::array<double, 3>& a, const std::array<double, 3>& b) -> std::array<double, 3> {
        return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    };
    const std::array<double, 3> cr = cross(eA_hat, eB_hat);
    const double crn = std::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
    std::array<double, 3> n_hat = {cr[0] / crn, cr[1] / crn, cr[2] / crn};

    const std::array<double, 3> A_ref = {0.0, 0.0, 0.0};
    const std::array<double, 3> qA = {A_ref[0] + u0 * eA_hat[0], A_ref[1] + u0 * eA_hat[1],
                                       A_ref[2] + u0 * eA_hat[2]};
    const std::array<double, 3> qB = {qA[0] - D_true * n_hat[0], qA[1] - D_true * n_hat[1],
                                       qA[2] - D_true * n_hat[2]};
    const std::array<double, 3> B_ref = {qB[0] - v0 * eB_hat[0], qB[1] - v0 * eB_hat[1], qB[2] - v0 * eB_hat[2]};

    StraightFiberPair out;
    for (int k = 0; k < nA; ++k)
    {
        const double u = static_cast<double>(k - a_origin) * stepA;
        out.beadsA.push_back(
            {A_ref[0] + u * eA_hat[0], A_ref[1] + u * eA_hat[1], A_ref[2] + u * eA_hat[2]});
    }
    for (int k = 0; k < nB; ++k)
    {
        const double v = static_cast<double>(k - b_origin) * stepB;
        out.beadsB.push_back(
            {B_ref[0] + v * eB_hat[0], B_ref[1] + v * eB_hat[1], B_ref[2] + v * eB_hat[2]});
    }
    for (int d3 = 0; d3 < 3; ++d3) out.n_hat[d3] = n_hat[d3];
    return out;
}

// -----------------------------------------------------------------------
// run_contact -- search + contributions + reduce + diagnostics для ЦЕЛОЙ
// сцены (все primitives, не один вручную собранный) -- именно так, как
// требует §8 ("все проверки сил -- на суммарных вкладах... по ВСЕМ
// primitives сцены").
// -----------------------------------------------------------------------
struct ContactRunResult
{
    interactions::FiberContactSearchResult search;
    interactions::FiberContactDiagnostics diag;
    std::vector<fiber::FiberDelta> out;
};

ContactRunResult run_contact(TwoFiberFixture& fx)
{
    rhs::CpuExecutor exec(2);
    ContactRunResult r;
    interactions::build_fiber_contact_search(exec, fx.state, fx.topology, fx.materials, r.search);
    r.diag = interactions::compute_fiber_contact_diagnostics(exec, fx.state, fx.materials, r.search);
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
    interactions::compute_fiber_contact_contributions(exec, fx.state, fx.materials, r.search.primitives, buffers);
    r.out = rhs::reduce_contributions(exec, buffers, fx.state.size());
    return r;
}

// Суммарная сила (Σ m_k*dv_k) на диапазон бусин [first,first+count) -- "сила
// на нить" (§8.1: "суммарная сила на нить A по всем бусинам").
void total_force(const fiber::FiberState& state, const std::vector<fiber::FiberDelta>& out, std::size_t first,
                  std::size_t count, double F[dim])
{
    for (std::size_t d = 0; d < dim; ++d) F[d] = 0.0;
    for (std::size_t k = 0; k < count; ++k)
    {
        const std::size_t i = first + k;
        for (std::size_t d = 0; d < dim; ++d) F[d] += state.mass[i] * out[i].dv[d];
    }
}

// Проверка сохранения импульса (§8.4: Σ m_k*a_k = 0) по ВСЕЙ системе (обе
// нити сразу) -- допуск проекта, относительный к масштабу самих вкладов
// (как в старом файле).
void check_momentum_zero(const std::string& label, const fiber::FiberState& state,
                          const std::vector<fiber::FiberDelta>& out)
{
    double total[dim] = {};
    double scale = 0.0;
    for (std::size_t i = 0; i < state.size(); ++i)
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double p = state.mass[i] * out[i].dv[d];
            total[d] += p;
            scale = std::max(scale, std::fabs(p));
        }
    if (scale < 1e-30) scale = 1.0;
    for (std::size_t d = 0; d < dim; ++d)
        check(std::fabs(total[d]) < 1e-9 * scale,
              label + ": сумма m_k*a_k по компоненте " + std::to_string(d) + " равна нулю (" +
                  std::to_string(total[d]) + ")");
}

// -----------------------------------------------------------------------
// FoundContact / find_all_contacts -- находит ВСЕ принятые (geom.contact)
// primitives сцены через build_fiber_contact_search +
// compute_fiber_contact_geometry по каждому. Используется там, где тесту
// нужны s,t/примитив КОНКРЕТНОГО контакта (8.3 -- кто владелец, 8.7 --
// dt_limit вручную, 8.8 -- трение), а не только суммарная сила.
// -----------------------------------------------------------------------
struct FoundContact
{
    interactions::FiberContactPrimitive prim;
    interactions::FiberContactGeometry geom;
};

std::vector<FoundContact> find_all_contacts(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                             const fiber::FiberTopology& topology,
                                             const std::vector<fiber::FiberMaterial>& materials,
                                             interactions::FiberContactSearchResult* out_search = nullptr)
{
    interactions::FiberContactSearchResult search;
    interactions::build_fiber_contact_search(exec, state, topology, materials, search);
    std::vector<FoundContact> result;
    for (const auto& p : search.primitives)
    {
        interactions::FiberContactGeometry g = interactions::compute_fiber_contact_geometry(state, materials, p);
        if (g.contact) result.push_back({p, g});
    }
    if (out_search) *out_search = search;
    return result;
}

// =========================================================================
// 8.1. Ровно один контакт на пересечение (главный тест).
//
// Две прямые нити по 7 бусин, скрещивание под углом theta, проникновение h
// (D_true=0.9*d). Перебор: theta in {90,60,30,10,3} град x положение
// пересечения in {середина сегмента на обеих нитях; узел на одной нити;
// узел-узел на обеих; 1e-12*l0 от узла}. Геометрия сцены строится АНАЛИТИЧЕСКИ
// (build_straight_pair): точный минимум расстояния между бесконечными
// прямыми задаётся параметрами u0 (на A, от центральной бусины index=3) и v0
// (на B, аналогично), поэтому "D_true" в сцене СОВПАДАЕТ с "D_true" в
// ручной формуле силы Герца с точностью double, а не приближённо.
// =========================================================================
struct CrossBuildResult
{
    TwoFiberFixture fx;
    double n_hat[dim] = {}; // от B к A (см. build_straight_pair)
};

CrossBuildResult build_cross_variant(double theta_deg, double u0_factor, double v0_factor, double l0,
                                      double D_true, double E, double nu, double d, double rho)
{
    const double theta = theta_deg * kPi / 180.0;
    const std::array<double, 3> eA = {1.0, 0.0, 0.0};
    const std::array<double, 3> eB = {std::cos(theta), std::sin(theta), 0.0};
    const int N = 7;
    const int origin = 3; // центральная бусина (индекс 3 из 0..6) -- u=0/v=0

    StraightFiberPair pair =
        build_straight_pair(N, origin, l0, N, origin, l0, eA, eB, D_true, u0_factor * l0, v0_factor * l0);

    const std::vector<double> mA(static_cast<std::size_t>(N), 1.0e-6);
    const std::vector<double> mB(static_cast<std::size_t>(N), 1.0e-6);

    CrossBuildResult res;
    res.fx = make_two_fiber_fixture(pair.beadsA, mA, pair.beadsB, mB, E, nu, d, rho);
    for (int d3 = 0; d3 < 3; ++d3) res.n_hat[d3] = pair.n_hat[d3];
    return res;
}

void test_8_1_single_contact_per_intersection()
{
    std::printf("test_8_1_single_contact_per_intersection:\n");
    const double l0 = 0.01;
    const double E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;
    const double D_true = 0.9 * d; // проникновение h -- D_true < d_eff=d (одинаковые нити)

    struct Pos
    {
        const char* label;
        double u0f, v0f; // в единицах l0, от центральной бусины (index=3)
        double force_tol; // относительный допуск на F_A (см. ниже про near-node)
    };
    // near-node: снапинг s (kFiberContactSnapEta=1e-9) заменяет ИСТИННУЮ
    // ближайшую точку (u0=1e-12*l0 от узла) на сам узел (s=0 точно) -- при
    // этом D (и, значит, |F|) остаётся ТОЧНЫМ (узел и истинная точка
    // отличаются только вдоль оси нити, а D измеряется поперёк), но
    // направление n получает поперечную добавку порядка u0/D_true (точка
    // приложения силы на A сдвинута на u0 вдоль нити относительно точки, под
    // которую построена сцена B) -- это ОЖИДАЕМОЕ, ограниченное и физически
    // объяснимое следствие snapping (проверено отдельно: D совпадает с
    // D_true побитово, отклоняется только n), а не дефект ядра. При
    // u0=1e-12*l0, D_true=0.9*d=9e-4 отклонение ~u0/D_true~1.1e-11 --
    // используем допуск 1e-8 (запас x1000) именно для этого варианта,
    // остальные три варианта держат полные 1e-12.
    const double kSnapForceTol = 1e-8;
    const std::vector<Pos> positions = {
        {"mid-mid (середина сегментов)", -0.5, -0.5, 1e-12},
        {"node-mid (узел на A)", 0.0, -0.5, 1e-12},
        {"node-node (узел-узел)", 0.0, 0.0, 1e-12},
        {"near-node (1e-12*l0 от узла)", 1.0e-12, -0.5, kSnapForceTol},
    };
    const std::vector<double> angles = {90.0, 60.0, 30.0, 10.0, 3.0};

    const double delta = d - D_true; // d_eff=d (одинаковые нити)
    int variant_count = 0;

    std::printf("  таблица вариантов (угол x позиция): contact_count | rel_err(F_A)\n");
    for (double theta_deg : angles)
    {
        for (const Pos& pos : positions)
        {
            ++variant_count;
            CrossBuildResult built = build_cross_variant(theta_deg, pos.u0f, pos.v0f, l0, D_true, E, nu, d, rho);
            ContactRunResult run = run_contact(built.fx);

            const double F_expected = hertz_force(delta, built.fx.materials[0], built.fx.materials[0]);
            double F_A[dim], F_B[dim];
            total_force(built.fx.state, run.out, 0, built.fx.n_a, F_A);
            total_force(built.fx.state, run.out, built.fx.n_a, built.fx.n_b, F_B);

            double err = 0.0;
            for (std::size_t d3 = 0; d3 < dim; ++d3)
                err = std::max(err, std::fabs(F_A[d3] - F_expected * built.n_hat[d3]));
            const double rel_err = err / std::max(std::fabs(F_expected), 1e-30);

            std::printf("  theta=%5.1f  %-32s contact_count=%zu  rel_err=%.3e\n", theta_deg, pos.label,
                        run.diag.contact_count, rel_err);

            const std::string tag = "theta=" + std::to_string(theta_deg) + " pos=" + pos.label;
            check(run.diag.contact_count == 1, tag + ": contact_count == 1");
            check(rel_err < pos.force_tol,
                  tag + ": сила на A == hertz_force(delta)*n (допуск " + std::to_string(pos.force_tol) + ")");

            double sumF_scale = std::max(std::fabs(F_expected), 1e-30);
            bool balanced = true;
            for (std::size_t d3 = 0; d3 < dim; ++d3)
                if (std::fabs(F_A[d3] + F_B[d3]) > 1e-9 * sumF_scale) balanced = false;
            check(balanced, tag + ": сумма сил на обе нити равна нулю");
        }
    }
    std::printf("  всего вариантов: %d (5 углов x 4 положения)\n", variant_count);
}

// =========================================================================
// 8.2. Концы нитей.
// =========================================================================

// "Конец" нити (2 бусины, вторая -- свободный конец) касается СЕРЕДИНЫ
// сегмента "внутренней" нити (2 бусины, один сегмент). end_is_a выбирает,
// какая из двух нитей -- концевая (Surface-End vs End-Surface старой схемы).
TwoFiberFixture build_end_surface_scene(bool end_is_a, double l0, double D_true, double E, double nu, double d,
                                         double rho)
{
    // Конец -- индекс1=(0,0,0); середина сегмента другой нити -- (0,0,D_true).
    const std::vector<std::array<double, 3>> end_beads = {{-l0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    const std::vector<std::array<double, 3>> mid_beads = {{0.0, -0.5 * l0, D_true}, {0.0, 0.5 * l0, D_true}};
    const std::vector<double> m_end = {1e-6, 1e-6}, m_mid = {1e-6, 1e-6};

    if (end_is_a) return make_two_fiber_fixture(end_beads, m_end, mid_beads, m_mid, E, nu, d, rho);
    return make_two_fiber_fixture(mid_beads, m_mid, end_beads, m_end, E, nu, d, rho);
}

TwoFiberFixture build_end_end_scene(double l0, double gap, double E, double nu, double d, double rho)
{
    const std::vector<std::array<double, 3>> beadsA = {{0.0, 0.0, 0.0}, {0.0, 0.0, l0}};
    const std::vector<std::array<double, 3>> beadsB = {{0.0, 0.0, l0 + gap}, {0.0, 0.0, l0 + gap + l0}};
    const std::vector<double> mA = {1e-6, 1e-6}, mB = {1e-6, 1e-6};
    return make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
}

TwoFiberFixture build_isolated_bead_scene(double l0, double D_true, double E, double nu, double d, double rho)
{
    const std::vector<std::array<double, 3>> beadsA = {{0.0, 0.0, 0.0}}; // изолированная бусина
    const std::vector<std::array<double, 3>> beadsB = {{0.0, -0.5 * l0, D_true}, {0.0, 0.5 * l0, D_true}};
    const std::vector<double> mA = {1e-6}, mB = {1e-6, 1e-6};
    return make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
}

void test_8_2_fiber_ends()
{
    std::printf("test_8_2_fiber_ends:\n");
    const double l0 = 0.01, E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;
    const double D_true = 0.9 * d;
    rhs::CpuExecutor exec(2);

    // --- конец A на середине сегмента B (бывший End-Surface) ---------------
    {
        TwoFiberFixture fx = build_end_surface_scene(true, l0, D_true, E, nu, d, rho);
        std::vector<FoundContact> contacts = find_all_contacts(exec, fx.state, fx.topology, fx.materials);
        check(contacts.size() == 1, "конец A на середине сегмента B: ровно один контакт");
        if (contacts.size() == 1)
        {
            const auto& c = contacts[0];
            check(c.geom.s == 0.0 || c.geom.s == 1.0, "конец A: s==0 или 1 (крайний сегмент)");
            check(near_rel(c.geom.t, 0.5, 1e-9), "конец A: t==0.5 (середина сегмента B)");

            ContactRunResult run = run_contact(fx);
            const std::size_t end_bead =
                (c.geom.s >= 1.0) ? static_cast<std::size_t>(c.prim.a1) : static_cast<std::size_t>(c.prim.a0);
            for (std::size_t i = 0; i < fx.n_a; ++i)
            {
                if (i == end_bead) continue;
                bool any = false;
                for (std::size_t d3 = 0; d3 < dim; ++d3) any = any || (run.out[i].dv[d3] != 0.0);
                check(!any, "конец A: остальные бусины A вклада не получили (бусина " + std::to_string(i) + ")");
            }
            bool end_nonzero = false;
            for (std::size_t d3 = 0; d3 < dim; ++d3) end_nonzero = end_nonzero || (run.out[end_bead].dv[d3] != 0.0);
            check(end_nonzero, "конец A: концевая бусина получила вклад (вся сила A на ней)");
        }
    }

    // --- конец B на середине сегмента A (End-Surface, симметрично) ---------
    {
        TwoFiberFixture fx = build_end_surface_scene(false, l0, D_true, E, nu, d, rho);
        std::vector<FoundContact> contacts = find_all_contacts(exec, fx.state, fx.topology, fx.materials);
        check(contacts.size() == 1, "конец B на середине сегмента A: ровно один контакт");
        if (contacts.size() == 1)
        {
            const auto& c = contacts[0];
            check(c.geom.t == 0.0 || c.geom.t == 1.0, "конец B: t==0 или 1 (крайний сегмент)");
            check(near_rel(c.geom.s, 0.5, 1e-9), "конец B: s==0.5 (середина сегмента A)");
        }
    }

    // --- конец к концу (End-End) --------------------------------------------
    {
        const double gap = 0.5 * d; // < d -- гарантирован контакт
        TwoFiberFixture fx = build_end_end_scene(l0, gap, E, nu, d, rho);
        std::vector<FoundContact> contacts = find_all_contacts(exec, fx.state, fx.topology, fx.materials);
        check(contacts.size() == 1, "конец к концу: ровно один контакт");
        if (contacts.size() == 1)
        {
            const auto& c = contacts[0];
            check(c.geom.s == 1.0, "End-End: s==1 (конец A)");
            check(c.geom.t == 0.0, "End-End: t==0 (конец B)");
            check(near_rel(c.geom.D, gap, 1e-9), "End-End: D == gap");
        }
    }

    // --- изолированная бусина у середины сегмента B -------------------------
    {
        TwoFiberFixture fx = build_isolated_bead_scene(l0, D_true, E, nu, d, rho);
        std::vector<FoundContact> contacts = find_all_contacts(exec, fx.state, fx.topology, fx.materials);
        check(contacts.size() == 1, "изолированная бусина: ровно один контакт");
        if (contacts.size() == 1)
        {
            const auto& c = contacts[0];
            check(c.prim.a1 == -1, "изолированная бусина: точечный сегмент (a1==-1)");
            check(c.geom.s == 0.0, "изолированная бусина: s==0");
            check(near_rel(c.geom.t, 0.5, 1e-9), "изолированная бусина: t==0.5 (середина сегмента B)");
        }
    }
}

// =========================================================================
// 8.3. Изогнутая нить.
// =========================================================================

// Выпуклый изгиб: A0=(-l0,0,delta), A1=(0,0,0) [апекс, ближайшая точка к B],
// A2=(l0,0,delta) -- вершина "смотрит" на B (обе "руки" приподняты дальше от
// B). Производная расстояния-до-B вдоль каждого плеча в апексе положительна
// (доказано в отчёте разработки) -- апекс есть НАСТОЯЩИЙ локальный минимум
// ОБОИХ соседних сегментов, что и требуется для честной проверки правила
// владения (не "случайно совпавший узел прямой нити", как в 8.1).
TwoFiberFixture build_convex_curved_scene(std::size_t* apex_index)
{
    const double l0 = 1.5e-3, delta_z = 0.5e-3;
    const double d = 1.0e-3, E = 2.0e9, nu = 0.3, rho = 1200.0;
    const double D_true = 0.9 * d;

    const std::vector<std::array<double, 3>> beadsA = {
        {-l0, 0.0, delta_z}, {0.0, 0.0, 0.0}, {l0, 0.0, delta_z}};
    const std::vector<std::array<double, 3>> beadsB = {{0.0, -0.5 * l0, -D_true}, {0.0, 0.5 * l0, -D_true}};
    const std::vector<double> mA = {1e-6, 1e-6, 1e-6}, mB = {1e-6, 1e-6};

    if (apex_index) *apex_index = 1;
    return make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
}

// Вогнутое "седло": A0/A2 -- "плечи", каждое по отдельности проходит близко к
// B в СВОЕЙ середине (интерьерный минимум 0<s<1), апекс A1 -- ДАЛЕКО от B
// (вершина "смотрит" ОТ B). Числа подобраны так (см. отчёт разработки), что
// интерьерный минимум каждого плеча приходится примерно на s=0.5 с запасом
// до d_eff. B -- ОДИН сегмент (2 бусины), чтобы исключить побочные пары с
// другими сегментами B (упрощает проверку "ровно два контакта").
TwoFiberFixture build_concave_curved_scene(std::size_t* wing0_owner, std::size_t* wing1_owner)
{
    const double l0 = 1.5e-3;
    const double c = 1.6e-3;       // высота апекса над линией B (апекс НЕ в контакте)
    const double dz = -1.0432236e-3; // решение Delta^2+2c*Delta+l0^2=0 (см. отчёт)
    const double d = 1.4e-3, E = 2.0e9, nu = 0.3, rho = 1200.0;

    const double z_wing = c + dz;
    const std::vector<std::array<double, 3>> beadsA = {
        {-l0, 0.0, z_wing}, {0.0, 0.0, c}, {l0, 0.0, z_wing}};
    const std::vector<std::array<double, 3>> beadsB = {{0.0, -0.5e-3, 0.0}, {0.0, 0.5e-3, 0.0}};
    const std::vector<double> mA = {1e-6, 1e-6, 1e-6}, mB = {1e-6, 1e-6};

    if (wing0_owner) *wing0_owner = 0; // сегмент (A0,A1), владелец -- бусина 0
    if (wing1_owner) *wing1_owner = 1; // сегмент (A1,A2), владелец -- бусина 1
    return make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);
}

void test_8_3_curved_fiber()
{
    std::printf("test_8_3_curved_fiber:\n");
    rhs::CpuExecutor exec(2);

    // --- выпуклый изгиб: ровно один контакт, принят сегментом S_a0 с s=0 ---
    {
        std::size_t apex_index = 0;
        TwoFiberFixture fx = build_convex_curved_scene(&apex_index);
        std::vector<FoundContact> contacts = find_all_contacts(exec, fx.state, fx.topology, fx.materials);
        check(contacts.size() == 1, "выпуклый изгиб: ровно один контакт");
        if (contacts.size() == 1)
        {
            check(static_cast<std::size_t>(contacts[0].prim.a0) == apex_index,
                  "выпуклый изгиб: контакт принят сегментом, начинающимся в апексе (S_a0 == апекс)");
            check(contacts[0].geom.s == 0.0, "выпуклый изгиб: s==0 (апекс -- начало принятого сегмента)");
        }
    }

    // --- вогнутое седло: два контакта с разными q_A -------------------------
    {
        std::size_t w0 = 0, w1 = 0;
        TwoFiberFixture fx = build_concave_curved_scene(&w0, &w1);
        interactions::FiberContactSearchResult search;
        std::vector<FoundContact> contacts = find_all_contacts(exec, fx.state, fx.topology, fx.materials, &search);
        std::printf("    сегмент-пар построено: %zu, принятых контактов: %zu\n", search.primitives.size(),
                    contacts.size());
        check(contacts.size() == 2, "вогнутое седло: ровно два контакта (оба плеча)");
        if (contacts.size() == 2)
        {
            auto q_a = [&](const FoundContact& c) -> std::array<double, 3> {
                const std::array<double, 3> xa0 = {fx.state.x[0][c.prim.a0], fx.state.x[1][c.prim.a0],
                                                    fx.state.x[2][c.prim.a0]};
                const std::array<double, 3> xa1 = (c.prim.a1 == -1)
                                                       ? xa0
                                                       : std::array<double, 3>{fx.state.x[0][c.prim.a1],
                                                                                fx.state.x[1][c.prim.a1],
                                                                                fx.state.x[2][c.prim.a1]};
                return {xa0[0] + c.geom.s * (xa1[0] - xa0[0]), xa0[1] + c.geom.s * (xa1[1] - xa0[1]),
                        xa0[2] + c.geom.s * (xa1[2] - xa0[2])};
            };
            const std::array<double, 3> q0 = q_a(contacts[0]);
            const std::array<double, 3> q1 = q_a(contacts[1]);
            double dist2 = 0.0;
            for (int d3 = 0; d3 < 3; ++d3) dist2 += (q0[d3] - q1[d3]) * (q0[d3] - q1[d3]);
            check(std::sqrt(dist2) > 1e-6, "вогнутое седло: q_A двух контактов РАЗНЫЕ (не дубликат)");

            // Владельцы -- разные сегменты (оба плеча), не один и тот же.
            check(contacts[0].prim.a0 != contacts[1].prim.a0 ||
                      (contacts[0].geom.s != contacts[1].geom.s),
                  "вогнутое седло: контакты принадлежат разным сегментам/точкам");
        }
    }
}

// =========================================================================
// 8.4. Сохранение импульса на сценах 8.1-8.3.
// =========================================================================
void test_8_4_momentum_conservation()
{
    std::printf("test_8_4_momentum_conservation:\n");
    const double l0 = 0.01, E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;
    const double D_true = 0.9 * d;

    // Представительная сцена 8.1 (90 град, середина сегментов).
    {
        CrossBuildResult built = build_cross_variant(90.0, -0.5, -0.5, l0, D_true, E, nu, d, rho);
        ContactRunResult run = run_contact(built.fx);
        check_momentum_zero("8.1 (90 град, mid-mid)", built.fx.state, run.out);
    }
    // Мелкий угол -- несколько кандидатов-primitives (не только истинный минимум).
    {
        CrossBuildResult built = build_cross_variant(3.0, 0.0, 0.0, l0, D_true, E, nu, d, rho);
        ContactRunResult run = run_contact(built.fx);
        check_momentum_zero("8.1 (3 град, node-node)", built.fx.state, run.out);
    }

    // Сцены 8.2.
    {
        TwoFiberFixture fx = build_end_surface_scene(true, l0, D_true, E, nu, d, rho);
        ContactRunResult run = run_contact(fx);
        check_momentum_zero("8.2 (конец A на середине B)", fx.state, run.out);
    }
    {
        TwoFiberFixture fx = build_end_end_scene(l0, 0.5 * d, E, nu, d, rho);
        ContactRunResult run = run_contact(fx);
        check_momentum_zero("8.2 (конец к концу)", fx.state, run.out);
    }
    {
        TwoFiberFixture fx = build_isolated_bead_scene(l0, D_true, E, nu, d, rho);
        ContactRunResult run = run_contact(fx);
        check_momentum_zero("8.2 (изолированная бусина)", fx.state, run.out);
    }

    // Сцены 8.3.
    {
        TwoFiberFixture fx = build_convex_curved_scene(nullptr);
        ContactRunResult run = run_contact(fx);
        check_momentum_zero("8.3 (выпуклый изгиб)", fx.state, run.out);
    }
    {
        TwoFiberFixture fx = build_concave_curved_scene(nullptr, nullptr);
        ContactRunResult run = run_contact(fx);
        check_momentum_zero("8.3 (вогнутое седло)", fx.state, run.out);
    }
}

// =========================================================================
// 8.5. Полнота поиска при растянутых связях.
//
// Обе нити растянуты до l=1.3*l0 (bond.l0 в topology -- НОМИНАЛЬНАЯ, т.е.
// l0_nominal, координаты -- фактически растянутые: append_chain вычисляет
// l0 из фактических координат, поэтому bond.l0 переписывается вручную ПОСЛЕ
// построения). Нити -- 2 бусины каждая (один сегмент = весь их "далёкий
// конец"), почти встречные (угол между e_A и -e_B ~5 град), контакт у
// параметров s=t=0.95 (не прижаты к 1) на обеих. Построение -- той же
// аналитической схемой build_straight_pair (u0=v0=0.95*L, L=1.3*l0_nominal).
// =========================================================================
void test_8_5_stretched_bonds_search_completeness()
{
    std::printf("test_8_5_stretched_bonds_search_completeness:\n");
    const double l0_nominal = 0.01;
    const double stretch = 1.3;
    const double L = stretch * l0_nominal; // ФАКТИЧЕСКАЯ длина сегмента
    const double d = 1.0e-3, E = 2.0e9, nu = 0.3, rho = 1200.0;
    const double D_true = 0.95 * d;
    const double theta5 = 5.0 * kPi / 180.0;
    const double s0 = 0.95, t0 = 0.95;

    const std::array<double, 3> eA = {1.0, 0.0, 0.0};
    // Угол между eA и (-eB) == 5 град ("почти встречные"): -eB = (cos5,sin5,0).
    const std::array<double, 3> eB = {-std::cos(theta5), std::sin(theta5), 0.0};

    StraightFiberPair pair = build_straight_pair(2, 0, L, 2, 0, L, eA, eB, D_true, s0 * L, t0 * L);
    const std::vector<double> mA = {1e-6, 1e-6}, mB = {1e-6, 1e-6};
    TwoFiberFixture fx = make_two_fiber_fixture(pair.beadsA, mA, pair.beadsB, mB, E, nu, d, rho);

    // bond.l0 -- НОМИНАЛЬНОЕ значение (append_chain вычислило ФАКТИЧЕСКОЕ
    // растянутое l0=L -- переписываем на номинальное l0_nominal, как в
    // реальном расчёте: l0 фиксируется при построении сетки, а координаты
    // потом растягиваются под нагрузкой).
    for (auto& b : fx.topology.bonds)
    {
        b.l0 = l0_nominal;
        b.l_break = l0_nominal * 1e9;
    }

    rhs::CpuExecutor exec(2);
    interactions::FiberContactSearchResult search;
    interactions::build_fiber_contact_search(exec, fx.state, fx.topology, fx.materials, search);

    // Самопроверка сцены: расстояние между "якорными" бусинами a0,b0 (=A0,B0,
    // единственные бусины с next!=-1 в 2-бусинных нитях) БОЛЬШЕ номинального
    // радиуса d+2*l0_nominal -- старая (по номинальному l0) оценка радиуса
    // пропустила бы эту пару.
    double dist2 = 0.0;
    for (int d3 = 0; d3 < 3; ++d3)
    {
        const double dx = fx.state.x[d3][0] - fx.state.x[d3][fx.n_a];
        dist2 += dx * dx;
    }
    const double dist_a0_b0 = std::sqrt(dist2);
    const double r_nominal = d + 2.0 * l0_nominal;
    const double r_actual = d + 2.0 * L;

    std::printf("    s=t=%.2f, D_true/d=%.3f, |x_a0-x_b0|=%.6e m\n", s0, D_true / d, dist_a0_b0);
    std::printf("    d+2*l0(номинал)=%.6e m,  d+2*L(факт.)=%.6e m\n", r_nominal, r_actual);

    check(dist_a0_b0 > r_nominal,
          "самопроверка сцены: |x_a0-x_b0| > d+2*l0_nominal (номинальная оценка пропустила бы пару)");
    check(dist_a0_b0 < r_actual, "самопроверка сцены: |x_a0-x_b0| < d+2*L (новая оценка покрывает пару)");

    check(!search.primitives.empty(), "новый радиус (d+2*L_max по фактической длине): пара сегментов найдена");

    std::vector<FoundContact> contacts = find_all_contacts(exec, fx.state, fx.topology, fx.materials);
    check(!contacts.empty(), "контакт найден с новым радиусом поиска");
    if (!contacts.empty())
    {
        std::printf("    найден контакт: s=%.6f t=%.6f D/d=%.6f\n", contacts[0].geom.s, contacts[0].geom.t,
                    contacts[0].geom.D / d);
        check(near_rel(contacts[0].geom.s, s0, 1e-6), "8.5: s ~= 0.95 (не прижато к 1)");
        check(near_rel(contacts[0].geom.t, t0, 1e-6), "8.5: t ~= 0.95 (не прижато к 1)");
    }
}

// =========================================================================
// 8.6. Почти параллельные и совпадающие нити.
// =========================================================================
void test_8_6_near_parallel_and_parallel()
{
    std::printf("test_8_6_near_parallel_and_parallel:\n");
    const double l0 = 0.01, E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;
    rhs::CpuExecutor exec(2);

    // --- почти параллельные: угол 1e-4 рад, смещение по нормали < d --------
    {
        const double D_true = 0.5 * d; // < d -- смещение по нормали меньше диаметра
        CrossBuildResult built = build_cross_variant(1.0e-4 * 180.0 / kPi, -0.5, -0.5, l0, D_true, E, nu, d, rho);
        ContactRunResult run = run_contact(built.fx);

        std::printf("    contact_count=%zu (nA-1=nB-1=6 сегментов на нить)\n", run.diag.contact_count);
        check(run.diag.contact_count > 0, "почти параллельные: есть хотя бы один контакт");
        check(run.diag.contact_count <= 6, "почти параллельные: contact_count <= числа сегментов перекрытия (6)");

        bool all_finite = true;
        for (const auto& d3 : run.out)
            for (std::size_t k = 0; k < dim; ++k)
                if (!std::isfinite(d3.dv[k])) all_finite = false;
        check(all_finite, "почти параллельные: нет NaN/Inf в результирующих dv");

        check_momentum_zero("почти параллельные", built.fx.state, run.out);
    }

    // --- совпадающие оси: та же прямая, частичное перекрытие ----------------
    {
        const int N = 7;
        std::vector<std::array<double, 3>> beadsA, beadsB;
        for (int k = -3; k <= 3; ++k) beadsA.push_back({static_cast<double>(k) * l0, 0.0, 0.0});
        for (int k = -3; k <= 3; ++k) beadsB.push_back({(static_cast<double>(k) + 0.5) * l0, 0.0, 0.0});
        const std::vector<double> mA(static_cast<std::size_t>(N), 1e-6), mB(static_cast<std::size_t>(N), 1e-6);
        TwoFiberFixture fx = make_two_fiber_fixture(beadsA, mA, beadsB, mB, E, nu, d, rho);

        interactions::FiberContactSearchResult search;
        interactions::FiberContactDiagnostics diag;
        interactions::build_fiber_contact_search(exec, fx.state, fx.topology, fx.materials, search);
        diag = interactions::compute_fiber_contact_diagnostics(exec, fx.state, fx.materials, search);

        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
        interactions::compute_fiber_contact_contributions(exec, fx.state, fx.materials, search.primitives, buffers);
        std::vector<fiber::FiberDelta> out = rhs::reduce_contributions(exec, buffers, fx.state.size());

        std::printf("    совпадающие оси: degenerate_count=%zu, contact_count=%zu\n", diag.degenerate_count,
                    diag.contact_count);

        bool all_finite = true;
        for (const auto& d3 : out)
            for (std::size_t k = 0; k < dim; ++k)
                if (!std::isfinite(d3.dv[k])) all_finite = false;
        check(all_finite, "совпадающие оси: нет NaN/Inf");
        check(diag.degenerate_count > 0, "совпадающие оси: degenerate_count > 0");

        bool all_zero = true;
        for (const auto& d3 : out)
            for (std::size_t k = 0; k < dim; ++k)
                if (d3.dv[k] != 0.0) all_zero = false;
        check(all_zero, "совпадающие оси: силы нулевые (контакт не сформирован)");
    }
}

// =========================================================================
// 8.7. dt_limit.
// =========================================================================
void test_8_7_dt_limit()
{
    std::printf("test_8_7_dt_limit:\n");
    const double l0 = 0.01, E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;
    const double D_true = 0.9 * d;
    rhs::CpuExecutor exec(2);

    CrossBuildResult built = build_cross_variant(90.0, -0.5, -0.5, l0, D_true, E, nu, d, rho);
    interactions::FiberContactSearchResult search;
    std::vector<FoundContact> contacts =
        find_all_contacts(exec, built.fx.state, built.fx.topology, built.fx.materials, &search);
    check(contacts.size() == 1, "8.7: сцена 8.1 (90 град, mid-mid) даёт ровно один контакт (предпосылка)");
    if (contacts.empty()) return;

    const auto& c = contacts[0];
    const double m_a0 = built.fx.state.mass[c.prim.a0];
    const double m_a1 = built.fx.state.mass[static_cast<rhs::index_t>(c.prim.a1)];
    const double m_b0 = built.fx.state.mass[c.prim.b0];
    const double m_b1 = built.fx.state.mass[static_cast<rhs::index_t>(c.prim.b1)];

    const double s = c.geom.s, t = c.geom.t;
    const double inv_mA = (1.0 - s) * (1.0 - s) / m_a0 + s * s / m_a1;
    const double inv_mB = (1.0 - t) * (1.0 - t) / m_b0 + t * t / m_b1;
    const double m_A = 1.0 / inv_mA;
    const double m_B = 1.0 / inv_mB;
    const double m_red = (m_A * m_B) / (m_A + m_B);
    const double k_eff = 1.5 * c.geom.F / c.geom.delta;
    const double dt_expected = 0.8 / std::sqrt(k_eff / m_red);

    const double dt_actual = interactions::compute_fiber_contact_dt_limit(exec, built.fx.state, built.fx.materials,
                                                                            search.primitives);

    std::printf("    dt_expected=%.9e s, dt_actual=%.9e s\n", dt_expected, dt_actual);
    check(near_rel(dt_actual, dt_expected, 1e-12), "8.7: compute_fiber_contact_dt_limit == ручная формула (§6)");

    // fiber_contact_side_inv_mass должна давать то же самое 1/m_A/1/m_B.
    const double inv_mA_fn = interactions::fiber_contact_side_inv_mass(built.fx.state, c.prim.a0, c.prim.a1, s);
    const double inv_mB_fn = interactions::fiber_contact_side_inv_mass(built.fx.state, c.prim.b0, c.prim.b1, t);
    check(inv_mA_fn == inv_mA, "8.7: fiber_contact_side_inv_mass(A) совпадает с ручной формулой (побитово)");
    check(inv_mB_fn == inv_mB, "8.7: fiber_contact_side_inv_mass(B) совпадает с ручной формулой (побитово)");
}

// =========================================================================
// 8.8. Трение на новой геометрии.
// =========================================================================

void test_8_8_friction_coefficient_mixing()
{
    std::printf("test_8_8_friction_coefficient_mixing:\n");
    fiber::FiberMaterial a = make_material(2.0e9, 0.3, 1.0e-3, 1200.0);
    fiber::FiberMaterial b = make_material(2.0e9, 0.3, 1.0e-3, 1200.0);

    a.mu_ff = 0.4;
    b.mu_ff = 0.4;
    check(interactions::fiber_friction_coefficient(a, b) == 0.4, "равные mu_ff -> точное ==");

    a.mu_ff = 0.2;
    b.mu_ff = 0.8;
    check(near_rel(interactions::fiber_friction_coefficient(a, b), 0.4, 1e-12), "(0.2,0.8) -> sqrt(0.16)=0.4");

    a.mu_ff = 0.0;
    b.mu_ff = 0.5;
    check(interactions::fiber_friction_coefficient(a, b) == 0.0, "(0,0.5) -> 0");
}

// Набор представительных сцен для трения -- переиспользует билдеры 8.1/8.2
// (Surface-Surface ~ 8.1 90/mid-mid; Surface-End/End-Surface ~ 8.2;
// End-End ~ 8.2).
struct FrictionScene
{
    const char* label;
    TwoFiberFixture fx;
};

std::vector<FrictionScene> build_friction_scenes()
{
    const double l0 = 0.01, E = 2.0e9, nu = 0.3, d = 1.0e-3, rho = 1200.0;
    const double D_true = 0.9 * d;

    std::vector<FrictionScene> scenes;
    scenes.push_back(
        {"Surface-Surface", build_cross_variant(90.0, -0.5, -0.5, l0, D_true, E, nu, d, rho).fx});
    scenes.push_back({"Surface-End (конец A)", build_end_surface_scene(true, l0, D_true, E, nu, d, rho)});
    scenes.push_back({"End-Surface (конец B)", build_end_surface_scene(false, l0, D_true, E, nu, d, rho)});
    scenes.push_back({"End-End", build_end_end_scene(l0, 0.5 * d, E, nu, d, rho)});
    return scenes;
}

// Задаёт РАЗНЫЕ ненулевые скорости всем бусинам сцены (детерминированно, по
// индексу) -- нужны тестам трения (скольжение).
void set_distinct_velocities(fiber::FiberState& state)
{
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        const double k = static_cast<double>(i) + 1.0;
        set_velocity(state, i, 0.3 * std::sin(0.7 * k), -0.2 * std::cos(0.5 * k), 0.15 * std::sin(0.3 * k + 1.0));
    }
}

// --- test_friction_zero_mu_regression: при mu_ff=0 dv СОВПАДАЕТ ПОБИТОВО
// (==) с формулой ЧИСТО нормального контакта, посчитанной в тесте ИЗ geom
// (§8.8: "сравнивается ПОБИТОВО с формулой нового нормального контакта").
// Формула -- буквально как в compute_fiber_contact_contributions:
// d.dv[dd] = w*geom.F*geom.n[dd]/m, тот же порядок операций.
void test_8_8_friction_zero_mu_regression()
{
    std::printf("test_8_8_friction_zero_mu_regression:\n");
    rhs::CpuExecutor exec(2);

    for (auto& scene : build_friction_scenes())
    {
        scene.fx.materials[0].mu_ff = 0.0;
        set_distinct_velocities(scene.fx.state);

        std::vector<FoundContact> contacts = find_all_contacts(exec, scene.fx.state, scene.fx.topology,
                                                                 scene.fx.materials);
        check(contacts.size() == 1, std::string(scene.label) + ": ровно один контакт (предпосылка)");
        if (contacts.size() != 1) continue;
        const auto& c = contacts[0];

        ContactRunResult run = run_contact(scene.fx);

        auto expect_bead = [&](rhs::index_t idx, double w, double sign) {
            const double m = scene.fx.state.mass[idx];
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double expected = sign * w * c.geom.F * c.geom.n[d] / m;
                check(run.out[idx].dv[d] == expected, std::string(scene.label) + ": dv[" + std::to_string(idx) +
                                                            "][" + std::to_string(d) + "] == формула нормального контакта (mu=0)");
            }
        };

        if (c.geom.s < 1.0) expect_bead(c.prim.a0, 1.0 - c.geom.s, 1.0);
        if (c.prim.a1 != -1 && c.geom.s > 0.0) expect_bead(static_cast<rhs::index_t>(c.prim.a1), c.geom.s, 1.0);
        if (c.geom.t < 1.0) expect_bead(c.prim.b0, 1.0 - c.geom.t, -1.0);
        if (c.prim.b1 != -1 && c.geom.t > 0.0) expect_bead(static_cast<rhs::index_t>(c.prim.b1), c.geom.t, -1.0);
    }
}

// =========================================================================
// Консервативный отсев fiber_contact_surely_apart и предел шага закона 0 из ядра.
// =========================================================================

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

using V3 = std::array<double, 3>;

// Четыре бусины (a0,a1,b0,b1) в состоянии; материалы 0 и 1 с разным d. a1/b1 == -1 -- точечный сегмент.
struct PairCase
{
    fiber::FiberState state;
    std::vector<fiber::FiberMaterial> materials;
    interactions::FiberContactPrimitive prim;
};

PairCase make_pair_case(const V3& a0, const V3& a1, bool a_point, const V3& b0, const V3& b1, bool b_point,
                        double d_a, double d_b)
{
    PairCase c;
    c.state.resize(4);
    const V3* pts[4] = {&a0, &a1, &b0, &b1};
    for (std::size_t k = 0; k < 4; ++k)
    {
        for (std::size_t d = 0; d < 3; ++d)
        {
            c.state.x[d][k] = (*pts[k])[d];
            c.state.vx[d][k] = 0.0;
        }
        c.state.mass[k] = 1.0e-6;
        c.state.object_id[k] = k < 2 ? 0 : 1;
        c.state.material_id[k] = k < 2 ? 0 : 1;
    }
    c.materials.push_back(make_material(2.0e9, 0.3, d_a, 1200.0));
    c.materials.push_back(make_material(2.0e9, 0.3, d_b, 1200.0));
    c.prim.a0 = 0;
    c.prim.a1 = a_point ? -1 : 1;
    c.prim.b0 = 2;
    c.prim.b1 = b_point ? -1 : 3;
    return c;
}

V3 vadd(const V3& a, const V3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
V3 vmul(const V3& a, double k) { return {a[0] * k, a[1] * k, a[2] * k}; }
V3 vsub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double vdot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V3 vcross(const V3& a, const V3& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
V3 vunit(const V3& a) { return vmul(a, 1.0 / std::sqrt(vdot(a, a))); }

void test_surely_apart_conservative()
{
    std::printf("test_surely_apart_conservative:\n");
    std::mt19937_64 gen(20261112ULL);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    auto rand_dir = [&]() {
        const double z = 2.0 * U(gen) - 1.0, ph = 2.0 * kPi * U(gen);
        const double rr = std::sqrt(std::max(0.0, 1.0 - z * z));
        return V3{rr * std::cos(ph), rr * std::sin(ph), z};
    };

    std::size_t violations = 0, n_apart = 0, n_contact = 0, n_total = 0, n_far_apart = 0, n_far = 0;
    auto test_case = [&](const PairCase& c, bool far_expected) {
        const bool apart = interactions::fiber_contact_surely_apart(c.state, c.materials, c.prim);
        const auto geom = interactions::compute_fiber_contact_geometry(c.state, c.materials, c.prim);
        ++n_total;
        if (apart) ++n_apart;
        if (geom.contact) ++n_contact;
        if (apart && geom.contact) ++violations;
        if (far_expected)
        {
            ++n_far;
            if (apart) ++n_far_apart;
        }
    };

    // 1) Случайные пары в кубе (в т.ч. точечные сегменты, пересекающиеся и почти параллельные).
    for (int it = 0; it < 60000; ++it)
    {
        const double scale = (it % 3 == 0) ? 0.3 : 1.0;
        const double d_a = 0.01 + 0.2 * U(gen), d_b = 0.01 + 0.2 * U(gen);
        const V3 a0 = {scale * U(gen), scale * U(gen), scale * U(gen)};
        const V3 b0 = {scale * U(gen), scale * U(gen), scale * U(gen)};
        const bool a_point = (U(gen) < 0.15), b_point = (U(gen) < 0.15);
        const V3 ea = rand_dir();
        // Часть пар -- почти параллельные оси.
        V3 eb = rand_dir();
        if (it % 5 == 0) eb = vunit(vadd(ea, vmul(rand_dir(), 1.0e-6)));
        const V3 a1 = vadd(a0, vmul(ea, 0.5 * U(gen)));
        const V3 b1 = vadd(b0, vmul(eb, 0.5 * U(gen)));
        test_case(make_pair_case(a0, a1, a_point, b0, b1, b_point, d_a, d_b), false);
    }

    // 2) Адверсариальные: истинное расстояние d_eff*(1 +- eps). Нить A -- отрезок, B -- перпендикулярный
    //    отрезок/точка над внутренней точкой A (расстояние ровно D) либо сдвиг за конец A.
    const double eps_list[] = {1e-10, -1e-10, 1e-6, -1e-6, 1e-12, -1e-12, 1e-9, -1e-9, 0.0};
    for (int it = 0; it < 20000; ++it)
    {
        const double d_a = 0.01 + 0.1 * U(gen), d_b = 0.01 + 0.1 * U(gen);
        const double d_eff = 0.5 * (d_a + d_b);
        const double eps = eps_list[it % 9];
        const double dist = d_eff * (1.0 + eps);
        const V3 e1 = rand_dir();
        V3 e2 = rand_dir();
        e2 = vunit(vsub(e2, vmul(e1, vdot(e2, e1))));
        const V3 e3 = vcross(e1, e2);
        const double La = (U(gen) < 0.2) ? 0.0 : 0.05 + 0.5 * U(gen);
        const V3 p0 = {U(gen), U(gen), U(gen)};
        const double off = (it % 4 == 0) ? 100.0 : 0.0; // большие координаты -- проверка округления
        const V3 pbase = {p0[0] + off, p0[1] + off, p0[2] + off};
        const V3 a1 = vadd(pbase, vmul(e1, La));
        V3 qb;
        if (it % 3 == 0 && La > 0.0)
        {
            // Ближайшая точка -- внутри A.
            qb = vadd(vadd(pbase, vmul(e1, La * (0.1 + 0.8 * U(gen)))), vmul(e2, dist));
        }
        else
        {
            // Ближайшая точка -- за концом A: расстояние dist по диагонали (долю даёт угол).
            const double ang = 0.2 + 1.2 * U(gen);
            qb = vadd(vadd(a1, vmul(e1, dist * std::sin(ang))), vmul(e2, dist * std::cos(ang)));
        }
        const double Lb = (U(gen) < 0.3) ? 0.0 : 0.05 + 0.5 * U(gen);
        const V3 b0 = vsub(qb, vmul(e3, Lb * U(gen)));
        const V3 b1 = vadd(b0, vmul(e3, Lb));
        test_case(make_pair_case(pbase, a1, La == 0.0, b0, b1, Lb == 0.0, d_a, d_b), false);
    }

    // 3) Заведомо далёкие пары: фильтр должен срабатывать.
    for (int it = 0; it < 5000; ++it)
    {
        const double d_a = 0.01 + 0.1 * U(gen), d_b = 0.01 + 0.1 * U(gen);
        const V3 a0 = {U(gen), U(gen), U(gen)};
        const V3 a1 = vadd(a0, vmul(rand_dir(), 0.3 * U(gen)));
        const V3 b0 = {U(gen) + 5.0, U(gen), U(gen)};
        const V3 b1 = vadd(b0, vmul(rand_dir(), 0.3 * U(gen)));
        test_case(make_pair_case(a0, a1, (it % 7 == 0), b0, b1, (it % 11 == 0), d_a, d_b), true);
    }

    // 4) Большие абсолютные координаты (ревью Codex): округление середин и длин при |x| >> длины.
    //    Точный контрпример без абсолютного запаса: точка A в base-8, отрезок B [base-5, base+8] по оси x,
    //    d_eff = 3.25, base = 2^52: истинное расстояние 3 < d_eff -- контакт есть, отсев недопустим.
    {
        const double base = std::ldexp(1.0, 52);
        const PairCase ce = make_pair_case({base - 8.0, 0, 0}, {base - 8.0, 0, 0}, true, {base - 5.0, 0, 0},
                                           {base + 8.0, 0, 0}, false, 3.25, 3.25);
        const auto g = interactions::compute_fiber_contact_geometry(ce.state, ce.materials, ce.prim);
        check(g.contact, "контрпример 2^52: полная геометрия видит контакт (кейс содержателен)");
        check(!interactions::fiber_contact_surely_apart(ce.state, ce.materials, ce.prim),
              "контрпример 2^52: отсева нет");
        test_case(ce, false);
    }
    //    Семейство: адверсариальные пары у порога, сдвинутые на 2^20 .. 2^60 по всем осям.
    for (int it = 0; it < 8000; ++it)
    {
        const double off = std::ldexp(1.0, 20 + (it % 41));
        const double d_a = 0.5 + U(gen), d_b = 0.5 + U(gen);
        const double d_eff = 0.5 * (d_a + d_b);
        const double dist = d_eff * (1.0 + eps_list[it % 9]);
        const V3 e1 = rand_dir();
        V3 e2 = rand_dir();
        e2 = vunit(vsub(e2, vmul(e1, vdot(e2, e1))));
        const V3 e3 = vcross(e1, e2);
        const double La = (U(gen) < 0.2) ? 0.0 : 1.0 + 20.0 * U(gen);
        const V3 pbase = {off + U(gen), off + U(gen), off + U(gen)};
        const V3 a1 = vadd(pbase, vmul(e1, La));
        const V3 qb = vadd(vadd(pbase, vmul(e1, La * U(gen))), vmul(e2, dist));
        const double Lb = (U(gen) < 0.3) ? 0.0 : 1.0 + 20.0 * U(gen);
        const V3 b0 = vsub(qb, vmul(e3, Lb * U(gen)));
        const V3 b1 = vadd(b0, vmul(e3, Lb));
        test_case(make_pair_case(pbase, a1, La == 0.0, b0, b1, Lb == 0.0, d_a, d_b), false);
    }

    std::printf("    пар=%zu, контактов=%zu, отсеяно=%zu, далёких отсеяно %zu из %zu, нарушений=%zu\n", n_total,
                n_contact, n_apart, n_far_apart, n_far, violations);
    check(violations == 0, "surely_apart => нет контакта (нарушений нет)");
    check(n_contact > 1000, "в наборе достаточно настоящих контактов (тест не пустой)");
    check(n_far_apart == n_far, "заведомо далёкие пары отсеиваются фильтром");
    check(n_apart > n_total / 10, "фильтр отсеивает заметную долю пар");

    // Некорректный материал -> false (решает полная геометрия).
    PairCase bad = make_pair_case({0, 0, 0}, {1, 0, 0}, false, {0, 10, 0}, {1, 10, 0}, false, 0.01, 0.01);
    check(interactions::fiber_contact_surely_apart(bad.state, bad.materials, bad.prim),
          "корректная далёкая пара: отсев");
    bad.state.material_id[0] = 7;
    check(!interactions::fiber_contact_surely_apart(bad.state, bad.materials, bad.prim),
          "индекс материала вне диапазона -> false");
}

// Сцена с множеством контактов: nx нитей вдоль x и ny вдоль y, пересекающихся на разной глубине.
struct GridScene
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

GridScene make_grid_scene(int nf, double z_gap_scale)
{
    GridScene g;
    const double d = 1.0e-3, l0 = 0.01;
    g.materials.push_back(make_material(2.0e9, 0.3, d, 1200.0));
    const int nb = 12;
    const std::size_t total = static_cast<std::size_t>(2 * nf * nb);
    g.state.resize(total);
    g.topology.resize(total);
    std::size_t base = 0;
    int oid = 0;
    for (int j = 0; j < nf; ++j)
    {
        std::vector<std::array<double, 3>> beads;
        std::vector<double> masses;
        for (int i = 0; i < nb; ++i)
        {
            beads.push_back({(i - 5.5) * l0, 0.5 * l0 + j * 0.9 * l0, 1.0e-4 * z_gap_scale * j});
            masses.push_back(1.0e-6);
        }
        append_chain(g.state, g.topology, base, oid++, 0, beads, masses);
        base += static_cast<std::size_t>(nb);
    }
    for (int k = 0; k < nf; ++k)
    {
        std::vector<std::array<double, 3>> beads;
        std::vector<double> masses;
        for (int i = 0; i < nb; ++i)
        {
            beads.push_back({0.5 * l0 + k * 0.9 * l0 - 2.5 * l0, (i - 2.0) * l0, 0.95e-3 - 0.7e-4 * z_gap_scale * k});
            masses.push_back(1.0e-6 * (1.0 + 0.1 * k));
        }
        append_chain(g.state, g.topology, base, oid++, 0, beads, masses);
        base += static_cast<std::size_t>(nb);
    }
    return g;
}

void test_law0_dt_from_kernel()
{
    std::printf("test_law0_dt_from_kernel:\n");
    for (int scene = 0; scene < 3; ++scene)
    {
        // 0: плотная сетка с контактами; 1: контактов нет; 2: одна y-нить против x-нитей.
        GridScene g = (scene == 0) ? make_grid_scene(6, 1.0) : make_grid_scene(6, scene == 1 ? 40.0 : 1.0);
        if (scene == 1)
            for (std::size_t i = 0; i < g.state.size(); ++i) g.state.x[2][i] += (g.state.object_id[i] >= 6) ? 0.05 : 0.0;
        if (scene == 2)
            for (std::size_t i = 0; i < g.state.size(); ++i)
                if (g.state.object_id[i] >= 7) g.state.x[2][i] += 1.0; // в контакте остаётся одна y-нить (6) с x-нитями
        for (int nt : {1, 2, 7, 8})
        {
            rhs::CpuExecutor exec(nt);
            interactions::FiberContactSearchResult search;
            interactions::build_fiber_contact_search(exec, g.state, g.topology, g.materials, search);
            std::size_t ncontact = 0;
            for (const auto& pr : search.primitives)
                if (interactions::compute_fiber_contact_geometry(g.state, g.materials, pr).contact) ++ncontact;

            interactions::FiberContactStepParams p;
            double dt_kernel = -1.0;
            p.law0_dt_min_out = &dt_kernel;
            std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
            std::vector<interactions::FiberFrictionContact> slots;
            interactions::compute_fiber_contact_contributions(exec, g.state, g.materials, search.primitives, buffers,
                                                              slots, p);
            const double dt_ref = interactions::compute_fiber_contact_dt_limit(exec, g.state, g.materials,
                                                                                search.primitives);
            const std::string lab = "сцена " + std::to_string(scene) + ", потоков " + std::to_string(nt);
            std::printf("    %s: primitives=%zu, контактов=%zu, dt=%.9e\n", lab.c_str(), search.primitives.size(),
                        ncontact, dt_ref);
            check(same_bits(dt_kernel, dt_ref), lab + ": dt ядра == compute_fiber_contact_dt_limit (побитово)");
            // Независимый эталон без отсева fiber_contact_surely_apart (ревью Codex): прежняя формула
            // parallel_reduce(min) по полной геометрии всех primitive.
            const auto& prims = search.primitives;
            const double dt_nofilter = exec.parallel_reduce(
                prims.size(), std::numeric_limits<double>::max(),
                [&](std::size_t k) {
                    return interactions::fiber_contact_pair_dt_limit(
                        g.state, prims[k], interactions::compute_fiber_contact_geometry(g.state, g.materials, prims[k]));
                },
                [](double a, double b) { return std::min(a, b); });
            check(same_bits(dt_kernel, dt_nofilter), lab + ": dt ядра == эталон без отсева (побитово)");
            if (scene == 0) check(ncontact > 5, lab + ": в сцене много контактов (предпосылка)");
            if (scene == 1)
                check(ncontact == 0 && dt_kernel == std::numeric_limits<double>::max(),
                      lab + ": без контактов -> max()");
        }
    }
    // n == 0: пустой список primitive.
    {
        rhs::CpuExecutor exec(2);
        GridScene g = make_grid_scene(2, 1.0);
        std::vector<interactions::FiberContactPrimitive> none;
        interactions::FiberContactStepParams p;
        double dt_kernel = -1.0;
        p.law0_dt_min_out = &dt_kernel;
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
        std::vector<interactions::FiberFrictionContact> slots;
        interactions::compute_fiber_contact_contributions(exec, g.state, g.materials, none, buffers, slots, p);
        check(dt_kernel == std::numeric_limits<double>::max(), "n == 0: dt ядра = max()");
    }
}

// =========================================================================
// Поиск с отсечкой (FiberContactBuildOptions::cull_apart) и ядра по виду FiberContactPrimView.
// Эталон сборки -- последовательные find_pairs + make_fiber_contact_prim (без for_each_pair_segmented).
//  * полный режим: те же primitives (в том же порядке) и pair_count, что у эталона;
//  * режим отсечки: prim_k строго возрастает, primitives[j] == эталон[prim_k[j]], n_full == размеру эталона,
//    {k} == {k : !surely_apart(эталон[k])}, все k с контактом присутствуют, seg_bound верен;
//  * ядро/предел шага по culled-виду против полного вектора для законов 0..4 (в т. ч. режим проекции и
//    режим списка слотов трения): вклады, предел шага, active, слоты, pair_dt, нормальные слоты -- побитово.
// Сцены: сетка нитей, сетка с разными материалами, растянутые связи, одна пара, пусто, поиск без выживших.
// 1/2/7/8 потоков; объекты поиска переиспользуются между сценами (рабочая память растёт и сжимается).
// =========================================================================
struct CScene
{
    std::string name;
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

CScene scene_of(const std::string& name, const GridScene& g)
{
    CScene c;
    c.name = name;
    c.state = g.state;
    c.topology = g.topology;
    c.materials = g.materials;
    return c;
}

void give_velocities(CScene& c)
{
    for (std::size_t i = 0; i < c.state.size(); ++i)
        for (std::size_t d = 0; d < dim; ++d)
            c.state.vx[d][i] = 0.4 * std::sin(0.7 * static_cast<double>(i + 1) + static_cast<double>(d)) +
                                0.1 * static_cast<double>(d);
}

std::vector<CScene> make_culled_scenes()
{
    std::vector<CScene> v;
    {
        CScene c = scene_of("сетка 6+6", make_grid_scene(6, 1.0));
        c.materials[0].mu_ff = 0.3;
        give_velocities(c);
        v.push_back(c);
    }
    {
        CScene c = scene_of("сетка 4+4, два материала", make_grid_scene(4, 1.0));
        c.materials[0].mu_ff = 0.3;
        c.materials.push_back(make_material(7.0e8, 0.35, 0.8e-3, 900.0));
        c.materials[1].mu_ff = 0.5;
        for (std::size_t i = 0; i < c.state.size(); ++i) c.state.material_id[i] = c.state.object_id[i] % 2;
        give_velocities(c);
        v.push_back(c);
    }
    {
        CScene c = scene_of("сетка 6+6, зазор z x40", make_grid_scene(6, 40.0));
        c.materials[0].mu_ff = 0.3;
        give_velocities(c);
        v.push_back(c);
    }
    {
        // «Лестница» по z: нить с object_id = o сдвинута на 0.018*o. Соседние нити (|Δo| = 1) -- бусины на
        // расстоянии ~0.020 < радиуса поиска 0.021 (пары есть, n_full > 0), но середины сегментов не ближе
        // 0.018 > R = d + L = 0.011 -> отсечены все; |Δo| >= 2 -- вне радиуса поиска. Итог: m = 0.
        // (Сдвиг только одного семейства не годится: параллельные нити семейства в 0.009 < R остаются.)
        CScene c = scene_of("сетка 6+6, поиск без выживших", make_grid_scene(6, 1.0));
        for (std::size_t i = 0; i < c.state.size(); ++i)
            c.state.x[2][i] += 0.018 * static_cast<double>(c.state.object_id[i]);
        give_velocities(c);
        v.push_back(c);
    }
    {
        CScene c = scene_of("сетка 1+1 (n_full < 64)", make_grid_scene(1, 1.0));
        c.materials[0].mu_ff = 0.3;
        give_velocities(c);
        v.push_back(c);
    }
    {
        // Растянутые связи (8.5): радиус по текущей длине; одна пара.
        const double l0_nominal = 0.01, L = 1.3 * l0_nominal, d = 1.0e-3;
        const double th = 5.0 * kPi / 180.0;
        StraightFiberPair pair = build_straight_pair(2, 0, L, 2, 0, L, {1.0, 0.0, 0.0}, {-std::cos(th), std::sin(th), 0.0},
                                                      0.95 * d, 0.95 * L, 0.95 * L);
        TwoFiberFixture fx = make_two_fiber_fixture(pair.beadsA, {1e-6, 1e-6}, pair.beadsB, {1e-6, 1e-6}, 2.0e9, 0.3, d,
                                                    1200.0);
        for (auto& b : fx.topology.bonds)
        {
            b.l0 = l0_nominal;
            b.l_break = l0_nominal * 1e9;
        }
        CScene c;
        c.name = "растянутые связи, одна пара";
        c.state = fx.state;
        c.topology = fx.topology;
        c.materials = fx.materials;
        c.materials[0].mu_ff = 0.3;
        give_velocities(c);
        v.push_back(c);
    }
    {
        // Одна пара: два креста из 2 бусин.
        const double d = 1.0e-3, step = 0.01;
        StraightFiberPair pair = build_straight_pair(2, 0, step, 2, 0, step, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                                      0.9 * d, 0.5 * step, 0.5 * step);
        TwoFiberFixture fx = make_two_fiber_fixture(pair.beadsA, {1e-6, 1e-6}, pair.beadsB, {1e-6, 1e-6}, 2.0e9, 0.3, d,
                                                    1200.0);
        CScene c;
        c.name = "одна пара";
        c.state = fx.state;
        c.topology = fx.topology;
        c.materials = fx.materials;
        c.materials[0].mu_ff = 0.3;
        give_velocities(c);
        v.push_back(c);
    }
    {
        CScene c = scene_of("пусто (n = 0)", make_grid_scene(1, 1.0));
        c.state.resize(0);
        c.topology = fiber::FiberTopology{};
        c.topology.resize(0);
        v.push_back(c);
    }
    return v;
}

bool prim_eq(const interactions::FiberContactPrimitive& x, const interactions::FiberContactPrimitive& y)
{
    return x.a0 == y.a0 && x.a1 == y.a1 && x.a_prev == y.a_prev && x.a_has_next_segment == y.a_has_next_segment &&
           x.b0 == y.b0 && x.b1 == y.b1 && x.b_prev == y.b_prev && x.b_has_next_segment == y.b_has_next_segment;
}

bool prims_eq(const std::vector<interactions::FiberContactPrimitive>& a,
              const std::vector<interactions::FiberContactPrimitive>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!prim_eq(a[i], b[i])) return false;
    return true;
}

bool fslot_eq(const interactions::FiberFrictionContact& a, const interactions::FiberFrictionContact& b)
{
    bool eq = a.active == b.active && a.a0 == b.a0 && a.a1 == b.a1 && a.b0 == b.b0 && a.b1 == b.b1 &&
              same_bits(a.s, b.s) && same_bits(a.t, b.t) && same_bits(a.N, b.N) && same_bits(a.mu, b.mu);
    for (std::size_t d = 0; d < dim; ++d) eq = eq && same_bits(a.n[d], b.n[d]);
    return eq;
}

bool nslot_eq(const interactions::FiberNormalContact& a, const interactions::FiberNormalContact& b)
{
    bool eq = a.active == b.active && a.a0 == b.a0 && a.a1 == b.a1 && a.b0 == b.b0 && a.b1 == b.b1 &&
              a.prim == b.prim && a.seg == b.seg && a.key == b.key && same_bits(a.s, b.s) && same_bits(a.t, b.t) &&
              same_bits(a.delta, b.delta) && same_bits(a.mu, b.mu);
    for (std::size_t d = 0; d < dim; ++d) eq = eq && same_bits(a.n[d], b.n[d]);
    return eq;
}

bool deltas_eq(const std::vector<fiber::FiberDelta>& a, const std::vector<fiber::FiberDelta>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t d = 0; d < dim; ++d)
            if (!same_bits(a[i].dv[d], b[i].dv[d])) return false;
    return true;
}

// Эталонная сборка: последовательные find_pairs + make_fiber_contact_prim.
std::vector<interactions::FiberContactPrimitive> reference_primitives(const CScene& sc, std::size_t& pair_count)
{
    std::vector<interactions::FiberContactPrimitive> prims;
    pair_count = 0;
    const std::size_t n = sc.state.size();
    if (n == 0) return prims;
    rhs::CpuExecutor exec1(1);
    const double radius = interactions::fiber_contact_search_radius(exec1, sc.state, sc.topology, sc.materials);
    std::array<const double*, dim> coords{};
    for (std::size_t d = 0; d < dim; ++d) coords[d] = sc.state.x[d].data();
    rhs::SpatialIndex<dim> index;
    index.build(coords, n, radius > 0.0 ? radius : 1.0);
    const auto pairs = index.find_pairs(radius, [&](rhs::index_t a, rhs::index_t b) {
        return sc.state.object_id[a] != sc.state.object_id[b];
    });
    pair_count = pairs.size();
    for (const auto& pr : pairs)
    {
        interactions::FiberContactPrimitive prim;
        if (interactions::make_fiber_contact_prim(sc.topology, pr.first, pr.second, prim)) prims.push_back(prim);
    }
    return prims;
}

struct KernelRun
{
    bool active = false;
    std::vector<fiber::FiberDelta> dv;
    std::vector<interactions::FiberFrictionContact> slots;
    std::vector<double> pair_dt;
    double law0 = -1.0, dt_slots = -1.0, dt_recalc = -1.0;
    std::vector<interactions::FiberNormalContact> normal;
    std::vector<interactions::FrictionActiveRef> act;
};

KernelRun run_view_kernel(const rhs::CpuExecutor& exec, const CScene& sc, const interactions::FiberContactPrimView& view,
                          int law, bool proj, bool list, const std::vector<double>& bead_dt)
{
    KernelRun r;
    interactions::FiberContactStepParams p;
    p.law.type = law;
    p.dt_prev = 1.0e-7;
    p.bead_dt_ref = &bead_dt;
    std::vector<rhs::CpuExecutor::SegPart<interactions::FiberNormalContact>> nparts;
    std::vector<rhs::CpuExecutor::SegPart<interactions::FrictionActiveRef>> aparts;
    if (law != 0) p.pair_dt_limit = &r.pair_dt;
    if (law == 0) p.law0_dt_min_out = &r.law0;
    if (proj)
    {
        p.normal_ff = &r.normal;
        p.normal_ff_parts = &nparts;
    }
    if (list)
    {
        p.fric_active = &r.act;
        p.fric_active_parts = &aparts;
    }
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buf;
    r.active = interactions::compute_fiber_contact_contributions(exec, sc.state, sc.materials, view, buf, r.slots, p);
    r.dv = rhs::reduce_contributions(exec, buf, sc.state.size());
    r.dt_slots = interactions::compute_fiber_contact_dt_limit(exec, sc.state, sc.materials, view, p);
    p.pair_dt_limit = nullptr;
    r.dt_recalc = interactions::compute_fiber_contact_dt_limit(exec, sc.state, sc.materials, view, p);
    return r;
}

// Сравнение прогона по culled-виду (c) с прогоном по полному вектору (f): j-й выживший -- k = prim_k[j].
bool kernel_runs_match(const KernelRun& c, const KernelRun& f, const std::vector<std::size_t>& prim_k,
                       std::size_t n_full, int law, bool proj, bool list)
{
    const std::size_t m = prim_k.size();
    bool ok = c.active == f.active && deltas_eq(c.dv, f.dv) && same_bits(c.law0, f.law0) &&
              same_bits(c.dt_slots, f.dt_slots) && same_bits(c.dt_recalc, f.dt_recalc);
    ok = ok && c.slots.size() == m && f.slots.size() == n_full;
    std::vector<long long> j_of_k(n_full, -1);
    for (std::size_t j = 0; j < m; ++j) j_of_k[prim_k[j]] = static_cast<long long>(j);
    for (std::size_t k = 0; ok && k < n_full; ++k)
    {
        const interactions::FiberFrictionContact none{};
        const auto& cs = j_of_k[k] >= 0 ? c.slots[static_cast<std::size_t>(j_of_k[k])] : none;
        ok = ok && fslot_eq(cs, f.slots[k]);
    }
    if (law != 0)
    {
        ok = ok && c.pair_dt.size() == m && f.pair_dt.size() == n_full;
        const double kInf = std::numeric_limits<double>::max();
        for (std::size_t k = 0; ok && k < n_full; ++k)
            ok = ok && same_bits(j_of_k[k] >= 0 ? c.pair_dt[static_cast<std::size_t>(j_of_k[k])] : kInf, f.pair_dt[k]);
    }
    if (proj)
    {
        ok = ok && c.normal.size() == f.normal.size();
        for (std::size_t q = 0; ok && q < c.normal.size(); ++q)
        {
            interactions::FiberNormalContact cn = c.normal[q];
            ok = ok && cn.prim < m;
            if (!ok) break;
            cn.prim = prim_k[cn.prim];
            ok = ok && nslot_eq(cn, f.normal[q]);
        }
    }
    if (list)
    {
        ok = ok && c.act.size() == f.act.size();
        for (std::size_t q = 0; ok && q < c.act.size(); ++q)
            ok = ok && c.act[q].k < m && prim_k[c.act[q].k] == f.act[q].k && c.act[q].seg == f.act[q].seg;
    }
    return ok;
}

void test_culled_build()
{
    std::printf("test_culled_build:\n");
    constexpr std::size_t S = rhs::CpuExecutor::kDetSegments;
    bool saw_small = false, saw_m0 = false, saw_n0 = false, saw_cull = false, saw_contacts = false;
    interactions::FiberContactSearchResult full_s, cull_s; // переиспользуются между сценами и числами потоков
    const std::vector<CScene> scenes = make_culled_scenes();
    for (const CScene& sc : scenes)
    {
        std::size_t pair_ref = 0;
        const auto ref = reference_primitives(sc, pair_ref);
        const std::size_t n_full = ref.size();
        std::vector<std::size_t> exp_k, contact_k;
        for (std::size_t k = 0; k < n_full; ++k)
        {
            if (!interactions::fiber_contact_surely_apart(sc.state, sc.materials, ref[k])) exp_k.push_back(k);
            if (interactions::compute_fiber_contact_geometry(sc.state, sc.materials, ref[k]).contact)
                contact_k.push_back(k);
        }
        if (n_full > 0 && n_full < 64) saw_small = true;
        if (n_full > 0 && exp_k.empty()) saw_m0 = true;
        if (n_full == 0) saw_n0 = true;
        if (exp_k.size() < n_full) saw_cull = true;
        if (!contact_k.empty()) saw_contacts = true;

        // Полный прогон 1 поток -- опорный для вкладов (бит-в-бит между числами потоков).
        std::vector<double> bead_dt;
        {
            rhs::CpuExecutor exec1(1);
            if (sc.state.size() > 0) interactions::compute_fiber_bead_dt_ref(exec1, sc.state, sc.topology, sc.materials, bead_dt);
        }
        std::vector<KernelRun> ref1(5 * 4);
        if (sc.state.size() > 0)
        {
            rhs::CpuExecutor exec1(1);
            for (int law = 0; law <= 4; ++law)
                for (int mode = 0; mode < 4; ++mode)
                {
                    const bool proj = (mode & 1) != 0, list = (mode & 2) != 0;
                    if (proj && (law < 1 || law > 3)) continue;
                    ref1[static_cast<std::size_t>(law * 4 + mode)] =
                        run_view_kernel(exec1, sc, interactions::fiber_contact_view(ref), law, proj, list, bead_dt);
                }
        }

        for (int nt : {1, 2, 7, 8})
        {
            rhs::CpuExecutor exec(nt);
            const std::string lab = sc.name + ", потоков " + std::to_string(nt) + ": ";

            // --- полный режим ---
            interactions::build_fiber_contact_search(exec, sc.state, sc.topology, sc.materials, full_s);
            bool ok = prims_eq(full_s.primitives, ref) && full_s.pair_count == pair_ref && !full_s.culled &&
                      full_s.prim_k.empty() && full_s.n_full == n_full;
            for (std::size_t s = 0; s < S; ++s) ok = ok && full_s.seg_bound[s] == n_full * s / S;
            ok = ok && full_s.seg_bound[S] == n_full;
            check(ok, lab + "полный режим: primitives, pair_count, seg_bound == эталон (последовательная сборка)");

            // --- режим отсечки (L_max задан вызывающей стороной, кроме пустого состояния) ---
            interactions::FiberContactBuildOptions opt;
            opt.cull_apart = true;
            double lmax = 0.0;
            if (sc.state.size() > 0)
            {
                lmax = interactions::fiber_max_current_bond_length(exec, sc.state, sc.topology);
                opt.l_max = &lmax;
            }
            interactions::build_fiber_contact_search(exec, sc.state, sc.topology, sc.materials, cull_s, opt);
            const std::size_t m = cull_s.primitives.size();
            bool okc = cull_s.culled && cull_s.n_full == n_full && cull_s.pair_count == pair_ref &&
                       cull_s.prim_k == exp_k && m == exp_k.size();
            for (std::size_t j = 0; okc && j < m; ++j)
                okc = okc && (j == 0 || cull_s.prim_k[j - 1] < cull_s.prim_k[j]) &&
                      prim_eq(cull_s.primitives[j], ref[cull_s.prim_k[j]]);
            for (std::size_t k : contact_k)
                okc = okc && std::binary_search(cull_s.prim_k.begin(), cull_s.prim_k.end(), k);
            check(okc, lab + "отсечка: prim_k строго возрастает, primitives[j] == эталон[prim_k[j]], множество k, все контакты на месте");
            bool okb = cull_s.seg_bound[0] == 0 && cull_s.seg_bound[S] == m;
            for (std::size_t s = 0; s < S; ++s)
            {
                okb = okb && cull_s.seg_bound[s] <= cull_s.seg_bound[s + 1];
                for (std::size_t j = cull_s.seg_bound[s]; okb && j < cull_s.seg_bound[s + 1]; ++j)
                    okb = okb && rhs::CpuExecutor::seg_of_index(cull_s.prim_k[j], n_full) == s;
            }
            check(okb, lab + "отсечка: seg_bound -- границы отрезков полного индекса k");

            // Диагностика не определена для отсечённого поиска.
            bool threw = false;
            try
            {
                (void)interactions::compute_fiber_contact_diagnostics(exec, sc.state, sc.materials, cull_s);
            }
            catch (const std::logic_error&)
            {
                threw = true;
            }
            check(threw, lab + "compute_fiber_contact_diagnostics на culled-поиске -> logic_error");

            // --- ядра по виду ---
            if (sc.state.size() == 0) continue;
            const auto cview = interactions::fiber_contact_view(cull_s);
            const auto fview = interactions::fiber_contact_view(full_s);
            bool okk = true, ok_thr = true, saw_active = false;
            for (int law = 0; law <= 4; ++law)
                for (int mode = 0; mode < 4; ++mode)
                {
                    const bool proj = (mode & 1) != 0, list = (mode & 2) != 0;
                    if (proj && (law < 1 || law > 3)) continue;
                    const KernelRun c = run_view_kernel(exec, sc, cview, law, proj, list, bead_dt);
                    const KernelRun f = run_view_kernel(exec, sc, fview, law, proj, list, bead_dt);
                    okk = okk && kernel_runs_match(c, f, cull_s.prim_k, n_full, law, proj, list);
                    // полный вид == прежний вызов по вектору; опорный прогон 1 потока
                    const KernelRun& r1 = ref1[static_cast<std::size_t>(law * 4 + mode)];
                    ok_thr = ok_thr && deltas_eq(c.dv, r1.dv) && deltas_eq(f.dv, r1.dv) && same_bits(c.dt_slots, r1.dt_slots);
                    if (c.active) saw_active = true;
                }
            check(okk, lab + "ядро/предел шага по culled-виду == полный вектор, законы 0..4 (вклады, active, слоты, pair_dt, нормальные слоты, список)");
            check(ok_thr, lab + "вклады и dt_limit не зависят от числа потоков (== полный прогон на 1 потоке)");
            if (sc.name == "сетка 6+6") check(saw_active, lab + "в сетке есть активные контакты (предпосылка)");
        }
    }
    check(saw_small, "сцена с 0 < n_full < 64 (предпосылка)");
    check(saw_m0, "сцена с n_full > 0 и m = 0 (предпосылка)");
    check(saw_n0, "сцена с n_full = 0 (предпосылка)");
    check(saw_cull, "отсечка убирает primitive (предпосылка)");
    check(saw_contacts, "в сценах есть контакты (предпосылка)");
}

} // namespace

int main()
{
    test_8_1_single_contact_per_intersection();
    test_8_2_fiber_ends();
    test_8_3_curved_fiber();
    test_8_4_momentum_conservation();
    test_8_5_stretched_bonds_search_completeness();
    test_8_6_near_parallel_and_parallel();
    test_8_7_dt_limit();
    test_8_8_friction_coefficient_mixing();
    test_8_8_friction_zero_mu_regression();
    test_surely_apart_conservative();
    test_law0_dt_from_kernel();
    test_culled_build();

    if (failures == 0)
        std::printf("test_fiber_contact: OK\n");
    else
        std::printf("test_fiber_contact: %d failure(s)\n", failures);

    return failures;
}
