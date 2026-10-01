// tests/test_clamp.cpp
// Unit-тесты кинематического закрепления (interactions/clamp.hpp;
// fiber/SPH.md, «Clamp (кинематическое закрепление)»): U1 список, U2 проекция и детерминизм
// реакции, U3 блочная сумма, U4 счётчики трения, U5 mark_clamp_box,
// U6 директива CLAMP_BOX в файле нитей, плюс validate_clamps. Параллельные
// варианты -- явные rhs::CpuExecutor(1), (2), (7). Double сравниваются
// побитово. Собирается с SPH_DIM=3.
#include "fiber/fiber_input.hpp"
#include "fiber/fiber_system.hpp"
#include "fiber/model.hpp"
#include "interactions/clamp.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/cpu_executor.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

static_assert(dim == 3, "test_clamp требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool same_bits(double a, double b)
{
    std::uint64_t ua, ub;
    std::memcpy(&ua, &a, sizeof a);
    std::memcpy(&ub, &b, sizeof b);
    return ua == ub;
}

bool same_vec(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!same_bits(a[i], b[i])) return false;
    return true;
}

template <class F>
bool throws(F&& f)
{
    try { f(); }
    catch (const std::runtime_error&) { return true; }
    return false;
}

using interactions::ClampList;
using interactions::ClampSet;
using interactions::ClampWorkspace;
using interactions::kClampAll;
using interactions::kClampReduceBlock;
using interactions::kClampSerialThreshold;

// Случайные маски: доля clamped_frac частиц закреплена, маска 1..kClampAll.
std::vector<std::uint8_t> random_mask(std::mt19937_64& rng, std::size_t n, std::size_t n_clamped)
{
    std::vector<std::uint8_t> m(n, 0);
    std::uniform_int_distribution<std::size_t> pick(0, n - 1);
    std::uniform_int_distribution<int> mk(1, kClampAll);
    std::size_t placed = 0;
    while (placed < n_clamped)
    {
        const std::size_t i = pick(rng);
        if (m[i] != 0) continue;
        m[i] = static_cast<std::uint8_t>(mk(rng));
        ++placed;
    }
    return m;
}

std::vector<double> random_vec(std::mt19937_64& rng, std::size_t n, double lo, double hi)
{
    std::uniform_real_distribution<double> u(lo, hi);
    std::vector<double> v(n);
    for (auto& e : v) e = u(rng);
    return v;
}

// Эталон суммы по той же блочной схеме, последовательно.
double reference_block_sum(const std::vector<double>& items, std::size_t n)
{
    double total = 0.0;
    for (std::size_t k0 = 0; k0 < n; k0 += kClampReduceBlock)
    {
        double s = 0.0;
        for (std::size_t k = k0; k < std::min(n, k0 + kClampReduceBlock); ++k) s += items[k];
        total += s;
    }
    return total;
}

// --------------------------------------------------------------------------
void test_u1()
{
    std::printf("U1. build_clamp_list\n");
    std::mt19937_64 rng(101);
    const std::size_t n = 1000;
    auto mask = random_mask(rng, n, 300);
    for (std::size_t i = 0; i < 50; ++i) mask[i] = kClampAll; // гарантируем full
    ClampList l;
    interactions::build_clamp_list(mask, n, l);
    bool incr = true, masks_ok = true;
    for (std::size_t k = 0; k < l.size(); ++k)
    {
        if (k > 0 && l.index[k] <= l.index[k - 1]) incr = false;
        if (l.mask[k] == 0 || l.mask[k] != mask[l.index[k]]) masks_ok = false;
    }
    check(incr, "index строго возрастает");
    check(masks_ok, "маски списка совпадают с плотной и != 0");
    std::size_t nz = 0, nfull = 0;
    for (auto m : mask) { nz += (m != 0); nfull += (m == kClampAll); }
    check(l.size() == nz, "размер списка == числу ненулевых масок");
    bool full_ok = l.full.size() == nfull;
    for (std::size_t k = 0; k < l.full.size(); ++k)
    {
        if (mask[l.full[k]] != kClampAll) full_ok = false;
        if (k > 0 && l.full[k] <= l.full[k - 1]) full_ok = false;
    }
    check(full_ok, "full -- ровно kClampAll, по возрастанию");

    ClampList e;
    interactions::build_clamp_list({}, 10, e);
    check(e.size() == 0 && e.full.empty(), "пустая маска -> пустой список");
    interactions::build_clamp_list(std::vector<std::uint8_t>(10, 0), 10, e);
    check(e.size() == 0, "нулевая плотная маска -> пустой список");

    check(throws([&] { interactions::build_clamp_list(std::vector<std::uint8_t>(5, 1), 6, e); }),
          "ошибка: размер маски != n");
    std::vector<std::uint8_t> bad(10, 0);
    bad[3] = static_cast<std::uint8_t>(1u << dim);
    check(throws([&] { interactions::build_clamp_list(bad, 10, e); }), "ошибка: бит >= dim");
}

// --------------------------------------------------------------------------
// Solid-стиль: vector<vector<double>> + av; Fiber-стиль: vector<double>[dim], av=nullptr.
struct SolidData
{
    std::vector<std::vector<double>> dv, av;
    std::vector<double> mass;
};
struct FiberData
{
    std::vector<double> dv[dim];
    std::vector<double> mass;
};

SolidData make_solid(std::mt19937_64& rng, std::size_t n)
{
    SolidData s;
    s.mass = random_vec(rng, n, 0.1, 10.0);
    for (std::size_t d = 0; d < dim; ++d)
    {
        s.dv.push_back(random_vec(rng, n, -1e3, 1e3));
        s.av.push_back(random_vec(rng, n, -1.0, 1.0));
    }
    return s;
}

FiberData make_fiber(std::mt19937_64& rng, std::size_t n)
{
    FiberData f;
    f.mass = random_vec(rng, n, 0.1, 10.0);
    for (std::size_t d = 0; d < dim; ++d) f.dv[d] = random_vec(rng, n, -1e3, 1e3);
    return f;
}

void test_u2_size(std::size_t n, std::size_t nclamped)
{
    std::printf("U2. apply_clamp_projection: n=%zu, зажато %zu\n", n, nclamped);
    std::mt19937_64 rng(200 + nclamped);
    const auto mask = random_mask(rng, n, nclamped);
    ClampList list;
    interactions::build_clamp_list(mask, n, list);
    check(list.size() == nclamped, "размер списка");

    const SolidData s0 = make_solid(rng, n);
    const FiberData f0 = make_fiber(rng, n);

    // эталон: поэлементно, затем блочная схема
    std::vector<double> ref_items[dim];
    SolidData sref = s0;
    for (std::size_t d = 0; d < dim; ++d) ref_items[d].assign(list.size(), 0.0);
    for (std::size_t k = 0; k < list.size(); ++k)
    {
        const std::size_t i = list.index[k];
        for (std::size_t d = 0; d < dim; ++d)
        {
            if ((list.mask[k] >> d) & 1u)
            {
                ref_items[d][k] = -s0.mass[i] * s0.dv[d][i];
                sref.dv[d][i] = 0.0;
                sref.av[d][i] = 0.0;
            }
        }
    }
    double ref_force[dim];
    for (std::size_t d = 0; d < dim; ++d) ref_force[d] = reference_block_sum(ref_items[d], list.size());
    // эталон реакции Fiber: те же маски, свои массы и dv
    std::vector<double> ref_items_f[dim];
    for (std::size_t d = 0; d < dim; ++d) ref_items_f[d].assign(list.size(), 0.0);
    for (std::size_t k = 0; k < list.size(); ++k)
        for (std::size_t d = 0; d < dim; ++d)
            if ((list.mask[k] >> d) & 1u) ref_items_f[d][k] = -f0.mass[list.index[k]] * f0.dv[d][list.index[k]];
    double ref_force_f[dim];
    for (std::size_t d = 0; d < dim; ++d) ref_force_f[d] = reference_block_sum(ref_items_f[d], list.size());

    const int threads[] = {1, 2, 7};
    double first_force_s[dim] = {}, first_force_f[dim] = {};
    bool have_first = false;
    for (int nt : threads)
    {
        rhs::CpuExecutor exec(nt);
        const std::string tag = " (потоков " + std::to_string(nt) + ")";

        SolidData s = s0;
        ClampWorkspace ws;
        interactions::ClampSet cs;
        cs.solid_mask = mask;
        interactions::prepare_clamp_workspace(cs, n, 0, ws);
        double fs[dim] = {};
        interactions::apply_clamp_projection(exec, ws.solid, s.mass, s.dv, &s.av, ws, fs);
        bool dv_ok = true, av_ok = true;
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv_ok = dv_ok && same_vec(s.dv[d], sref.dv[d]);
            av_ok = av_ok && same_vec(s.av[d], sref.av[d]);
        }
        check(dv_ok, "Solid: dv: закреплённые == 0, остальные побитово прежние" + tag);
        check(av_ok, "Solid: av: закреплённые == 0, остальные побитово прежние" + tag);
        bool force_ref = true;
        for (std::size_t d = 0; d < dim; ++d) force_ref = force_ref && same_bits(fs[d], ref_force[d]);
        check(force_ref, "Solid: реакция == последовательной блочной сумме" + tag);

        // второй вызов прибавляет (force_out += ...): при нулевых dv добавка 0
        double fs2[dim] = {1.0, 2.0, 3.0};
        interactions::apply_clamp_projection(exec, ws.solid, s.mass, s.dv, &s.av, ws, fs2);
        check(same_bits(fs2[0], 1.0) && same_bits(fs2[1], 2.0) && same_bits(fs2[2], 3.0),
              "повторная проекция уже нулевых компонент не меняет накопитель" + tag);

        FiberData f = f0;
        interactions::ClampSet cf;
        cf.fiber_mask = mask;
        ClampWorkspace wf;
        interactions::prepare_clamp_workspace(cf, 0, n, wf);
        double ff[dim] = {};
        interactions::apply_clamp_projection(exec, wf.fiber, f.mass, f.dv,
                                             static_cast<std::vector<std::vector<double>>*>(nullptr), wf, ff);
        bool fdv_ok = true, free_ok = true;
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                const bool clamped = (mask[i] >> d) & 1u;
                if (clamped) fdv_ok = fdv_ok && same_bits(f.dv[d][i], 0.0);
                else free_ok = free_ok && same_bits(f.dv[d][i], f0.dv[d][i]);
            }
        }
        check(fdv_ok, "Fiber: закреплённые компоненты == 0" + tag);
        check(free_ok, "Fiber: свободные компоненты побитово прежние" + tag);
        bool fforce_ref = true;
        for (std::size_t d = 0; d < dim; ++d) fforce_ref = fforce_ref && same_bits(ff[d], ref_force_f[d]);
        check(fforce_ref, "Fiber: реакция == последовательной блочной сумме" + tag);
        if (!have_first)
        {
            for (std::size_t d = 0; d < dim; ++d) { first_force_s[d] = fs[d]; first_force_f[d] = ff[d]; }
            have_first = true;
        }
        else
        {
            bool eq = true;
            for (std::size_t d = 0; d < dim; ++d)
                eq = eq && same_bits(fs[d], first_force_s[d]) && same_bits(ff[d], first_force_f[d]);
            check(eq, "реакция Solid и Fiber побитово совпадает с результатом 1 потока" + tag);
        }
    }
}

void test_u2()
{
    test_u2_size(500, 100);
    test_u2_size(20000, 10000);
    test_u2_size(kClampSerialThreshold + 50, kClampSerialThreshold);     // ровно на пороге
    test_u2_size(kClampSerialThreshold, kClampSerialThreshold - 1);      // чуть ниже порога
    std::printf("U2. clamp_project_item (av == nullptr и av с nullptr-элементом)\n");
    double a[dim], b[dim], fv[dim];
    for (std::size_t d = 0; d < dim; ++d) { a[d] = 1.5 + d; b[d] = 7.0; }
    double* dvp[dim];
    double* avp[dim];
    for (std::size_t d = 0; d < dim; ++d) { dvp[d] = &a[d]; avp[d] = &b[d]; }
    avp[0] = nullptr;
    interactions::clamp_project_item(0x3, 2.0, dvp, avp, fv);
    check(same_bits(fv[0], -3.0) && same_bits(a[0], 0.0), "бит 0: сила и обнуление");
    check(same_bits(b[1], 0.0) && same_bits(b[2], 7.0), "av обнулён только у закреплённых с указателем");
    check(same_bits(fv[2], 0.0) && same_bits(a[2], 3.5), "свободная компонента не тронута, сила 0");
    interactions::clamp_project_item(kClampAll, 1.0, dvp, nullptr, fv);
    check(same_bits(a[1], 0.0) && same_bits(fv[2], -3.5), "av == nullptr допустим");
}

// --------------------------------------------------------------------------
void test_u3()
{
    std::printf("U3. clamp_fixed_block_sum\n");
    std::mt19937_64 rng(303);
    const std::size_t lengths[] = {0, 1, 255, 256, 257, 100000};
    for (std::size_t n : lengths)
    {
        const auto items = random_vec(rng, n == 0 ? 1 : n, -1e6, 1e6);
        std::vector<double> bp;
        double r1 = 0.0;
        bool ok = true;
        bool ref_ok = true;
        const int threads[] = {1, 2, 7};
        for (int nt : threads)
        {
            rhs::CpuExecutor exec(nt);
            std::vector<double> bpt((n + kClampReduceBlock - 1) / kClampReduceBlock, 0.0);
            const double r = interactions::clamp_fixed_block_sum(exec, items, n, bpt);
            if (nt == 1) r1 = r;
            else ok = ok && same_bits(r, r1);
            ref_ok = ref_ok && same_bits(r, reference_block_sum(items, n));
        }
        check(ok, "длина " + std::to_string(n) + ": побитово одинаково при 1, 2, 7 потоках");
        check(ref_ok, "длина " + std::to_string(n) + ": совпадает с последовательной блочной схемой");
    }
    std::vector<double> small_bp;
    rhs::CpuExecutor exec(2);
    const auto items = random_vec(rng, 1000, -1.0, 1.0);
    const double r = interactions::clamp_fixed_block_sum(exec, items, 1000, small_bp);
    check(small_bp.size() >= 4 && same_bits(r, reference_block_sum(items, 1000)),
          "block_partial растёт, если мал");
}

// --------------------------------------------------------------------------
void test_u4()
{
    std::printf("U4. apply_clamp_friction_counts\n");
    for (std::size_t n : {std::size_t(300), std::size_t(12000)})
    {
        std::mt19937_64 rng(400 + n);
        ClampSet cs;
        cs.solid_mask = random_mask(rng, n, n / 3);
        cs.fiber_mask = random_mask(rng, n / 2, n / 6);
        const auto cs0 = random_vec(rng, n, 1.0, 5.0);
        const auto cf0 = random_vec(rng, n / 2, 1.0, 5.0);
        const int threads[] = {1, 2, 7};
        for (int nt : threads)
        {
            rhs::CpuExecutor exec(nt);
            ClampWorkspace ws;
            interactions::prepare_clamp_workspace(cs, n, n / 2, ws);
            auto c_s = cs0;
            auto c_f = cf0;
            interactions::apply_clamp_friction_counts(exec, ws, c_f, c_s);
            bool ok = true;
            for (std::size_t i = 0; i < n; ++i)
                ok = ok && same_bits(c_s[i], cs.solid_mask[i] == kClampAll ? 0.0 : cs0[i]);
            for (std::size_t k = 0; k < n / 2; ++k)
                ok = ok && same_bits(c_f[k], cs.fiber_mask[k] == kClampAll ? 0.0 : cf0[k]);
            check(ok, "n=" + std::to_string(n) + ", потоков " + std::to_string(nt) +
                          ": обнулены ровно полностью закреплённые");
        }
    }
    // Параллельный путь: списки full длиннее kClampSerialThreshold.
    {
        const std::size_t n = 3 * kClampSerialThreshold;
        std::mt19937_64 rng(404);
        ClampSet cs;
        cs.solid_mask.assign(n, 0);
        cs.fiber_mask.assign(n, 0);
        for (std::size_t i = 0; i < n; i += 2) cs.solid_mask[i] = kClampAll;
        for (std::size_t i = 1; i < n; i += 2) cs.fiber_mask[i] = kClampAll;
        const auto cs0 = random_vec(rng, n, 1.0, 5.0);
        const auto cf0 = random_vec(rng, n, 1.0, 5.0);
        const int threads[] = {1, 2, 7};
        for (int nt : threads)
        {
            rhs::CpuExecutor exec(nt);
            ClampWorkspace ws;
            interactions::prepare_clamp_workspace(cs, n, n, ws);
            check(ws.solid.full.size() >= kClampSerialThreshold && ws.fiber.full.size() >= kClampSerialThreshold,
                  "списки full выше порога последовательного пути");
            auto c_s = cs0;
            auto c_f = cf0;
            interactions::apply_clamp_friction_counts(exec, ws, c_f, c_s);
            bool ok = true;
            for (std::size_t i = 0; i < n; ++i)
            {
                ok = ok && same_bits(c_s[i], cs.solid_mask[i] == kClampAll ? 0.0 : cs0[i]);
                ok = ok && same_bits(c_f[i], cs.fiber_mask[i] == kClampAll ? 0.0 : cf0[i]);
            }
            check(ok, "параллельный путь, потоков " + std::to_string(nt) + ": обнулены ровно полностью закреплённые");
        }
    }
    // пустой набор: ничего не меняется
    rhs::CpuExecutor exec(2);
    ClampWorkspace ws;
    std::vector<double> a(5, 2.0), b(5, 3.0);
    interactions::apply_clamp_friction_counts(exec, ws, a, b);
    check(same_bits(a[0], 2.0) && same_bits(b[4], 3.0), "пустой workspace не меняет счётчики");
}

// --------------------------------------------------------------------------
void test_u5()
{
    std::printf("U5. mark_clamp_box\n");
    // сетка 10x10x10, шаг 1
    const std::size_t g = 10, n = g * g * g;
    std::vector<double> x[dim];
    for (std::size_t d = 0; d < dim; ++d) x[d].resize(n);
    for (std::size_t i = 0; i < g; ++i)
        for (std::size_t j = 0; j < g; ++j)
            for (std::size_t k = 0; k < g; ++k)
            {
                const std::size_t p = (i * g + j) * g + k;
                x[0][p] = double(i); x[1][p] = double(j); x[2][p] = double(k);
            }
    const double lo1[dim] = {0.0, 0.0, 0.0}, hi1[dim] = {2.0, 9.0, 9.0};
    const double lo2[dim] = {2.0, 0.0, 0.0}, hi2[dim] = {4.0, 9.0, 9.0};

    std::vector<std::uint8_t> m1, m7;
    rhs::CpuExecutor e1(1), e7(7);
    interactions::mark_clamp_box(e1, x, n, lo1, hi1, 0x1, m1);
    check(m1.size() == n, "out_mask пустой -> размер n");
    bool edge = true;
    for (std::size_t p = 0; p < n; ++p)
    {
        const bool in = x[0][p] <= 2.0;
        edge = edge && (m1[p] == (in ? 1 : 0));
    }
    check(edge, "замкнутый бокс: частицы на грани x==2 закреплены, x==3 нет");
    interactions::mark_clamp_box(e1, x, n, lo2, hi2, 0x2, m1);
    bool orr = true;
    for (std::size_t p = 0; p < n; ++p)
    {
        std::uint8_t want = 0;
        if (x[0][p] <= 2.0) want |= 1;
        if (x[0][p] >= 2.0 && x[0][p] <= 4.0) want |= 2;
        orr = orr && m1[p] == want;
    }
    check(orr, "две области объединяются по «или» (пересечение x==2 -> 3)");
    interactions::mark_clamp_box(e7, x, n, lo1, hi1, 0x1, m7);
    interactions::mark_clamp_box(e7, x, n, lo2, hi2, 0x2, m7);
    check(m7 == m1, "результат при 1 и 7 потоках одинаков");
    std::vector<std::uint8_t> small(3, 0);
    interactions::mark_clamp_box(e7, x, n, lo1, hi1, 0x4, small);
    check(small.size() == n, "out_mask короче n дорастает до n");
}

// --------------------------------------------------------------------------
void test_validate()
{
    std::printf("validate_clamps\n");
    const std::size_t nt = 10, nf = 6;
    ParticleSystem sys;
    sys.resize(nt + 5); // запас под виртуальные частицы
    for (std::size_t i = 0; i < sys.size(); ++i)
    {
        for (std::size_t d = 0; d < dim; ++d) sys.vx[d][i] = 0.0;
        sys.hsml[i] = 0.01 * double(i + 1);
    }
    fiber::FiberState fb;
    fb.resize(nf);
    for (std::size_t d = 0; d < dim; ++d)
        for (std::size_t k = 0; k < nf; ++k) fb.vx[d][k] = 0.0;

    ClampSet ok;
    ok.solid_mask.assign(nt, 0);
    ok.solid_mask[2] = kClampAll;
    ok.fiber_mask.assign(nf, 0);
    ok.fiber_mask[0] = 1;
    check(!throws([&] { interactions::validate_clamps(ok, sys, nt, fb); }), "корректный набор проходит");

    ClampSet empty;
    check(!throws([&] { interactions::validate_clamps(empty, sys, nt, fb); }), "пустой набор проходит");

    ClampSet bad = ok;
    bad.solid_mask.assign(sys.size(), 0); // размер sys.size(), а не ntotal
    check(throws([&] { interactions::validate_clamps(bad, sys, nt, fb); }), "ошибка: solid_mask != ntotal");
    bad = ok;
    bad.fiber_mask.assign(nf + 1, 0);
    check(throws([&] { interactions::validate_clamps(bad, sys, nt, fb); }), "ошибка: fiber_mask != nfiber");
    bad = ok;
    bad.solid_mask[3] = static_cast<std::uint8_t>(1u << dim);
    check(throws([&] { interactions::validate_clamps(bad, sys, nt, fb); }), "ошибка: бит >= dim (Solid)");
    bad = ok;
    bad.fiber_mask[1] = static_cast<std::uint8_t>(1u << dim);
    check(throws([&] { interactions::validate_clamps(bad, sys, nt, fb); }), "ошибка: бит >= dim (Fiber)");

    sys.vx[1][2] = 0.5;
    check(throws([&] { interactions::validate_clamps(ok, sys, nt, fb); }),
          "ошибка: ненулевая скорость закреплённой компоненты Solid");
    sys.vx[1][2] = 0.0;
    ClampSet partial = ok;
    partial.solid_mask[2] = 0x1; // закреплена только x
    sys.vx[1][2] = 0.5;          // скорость по свободной компоненте допустима
    check(!throws([&] { interactions::validate_clamps(partial, sys, nt, fb); }),
          "ненулевая скорость по свободной компоненте допустима");
    sys.vx[1][2] = 0.0;
    fb.vx[0][0] = -1.0;
    check(throws([&] { interactions::validate_clamps(ok, sys, nt, fb); }),
          "ошибка: ненулевая скорость закреплённой бусины");
    fb.vx[0][0] = 0.0;

    ClampWorkspace ws;
    interactions::prepare_clamp_workspace(ok, nt, nf, ws);
    interactions::ClampReport rep;
    interactions::clamp_report_counts(ws, sys, rep);
    check(rep.solid_clamped == 1 && rep.fiber_clamped == 1, "clamp_report_counts: число зажатых");
    check(same_bits(rep.solid_max_h, sys.hsml[2]), "clamp_report_counts: solid_max_h");
}

// --------------------------------------------------------------------------
// Временный файл нитей: имя уникально для процесса (три регистрации теста
// идут параллельно при ctest -j), удаляется и при исключении (RAII).
struct TempFiberFile
{
    std::filesystem::path path;
    explicit TempFiberFile(const std::string& content)
    {
        static int counter = 0;
        path = std::filesystem::temp_directory_path() /
               ("test_clamp_box_" + std::to_string(std::random_device{}()) + "_" + std::to_string(counter++) + ".fib");
        std::ofstream out(path);
        out << content;
    }
    ~TempFiberFile()
    {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
};

// Читает файл нитей из строки и возвращает модель.
fiber::Model load_from_content(const std::string& content)
{
    TempFiberFile f(content);
    fiber::Model model;
    fiber::load_fiber_file(f.path.string(), model);
    return model;
}

// Текст std::runtime_error разбора или "" (исключения не было).
std::string load_error(const std::string& content)
{
    try
    {
        load_from_content(content);
    }
    catch (const std::runtime_error& e)
    {
        return e.what();
    }
    return "";
}

// --------------------------------------------------------------------------
void test_u6_clamp_box()
{
    std::printf("U6. CLAMP_BOX в файле нитей\n");

    // Тест 1: директива с маской и без маски
    {
        std::ostringstream content;
        content << "MATERIAL 0 2.0e9 0.35 1.0e-4 1200.0 0.05\n";
        content << "FIBER 0 3\n";
        content << "0.0 0.0 0.0\n";
        content << "1.0 0.0 0.0\n";
        content << "2.0 0.0 0.0\n";
        content << "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0 1.0\n"; // закрепляет бусину 0 полностью
        content << "CLAMP_BOX 1.5 -1.0 -1.0 2.5 1.0 1.0 xy\n"; // закрепляет бусину 2 по x,y

        auto model = load_from_content(content.str());
        check(!model.clamps.fiber_mask.empty(), "директива CLAMP_BOX: fiber_mask не пустой");
        check(model.clamps.fiber_mask.size() == 3, "директива CLAMP_BOX: размер fiber_mask = nfiber");
        check(model.clamps.fiber_mask[0] == interactions::kClampAll,
              "директива CLAMP_BOX: бусина 0 закреплена полностью");
        check(model.clamps.fiber_mask[2] == 0x3, "директива CLAMP_BOX: бусина 2 закреплена по xy (маска 0x3)");
        check(model.clamps.fiber_mask[1] == 0, "директива CLAMP_BOX: бусина 1 свободна");
    }

    // Тест 2: несколько директив объединяются по «или»
    {
        std::ostringstream content;
        content << "MATERIAL 0 2.0e9 0.35 1.0e-4 1200.0 0.05\n";
        content << "FIBER 0 2\n";
        content << "0.0 0.0 0.0\n";
        content << "1.0 0.0 0.0\n";
        content << "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0 1.0 x\n"; // бусина 0 по x
        content << "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0 1.0 y\n"; // бусина 0 по y

        auto model = load_from_content(content.str());
        check(model.clamps.fiber_mask[0] == 0x3, "две директивы на одну бусину: маска объединена");
    }

    // Тест 3: директива до и после FIBER блоков
    {
        std::ostringstream content;
        content << "MATERIAL 0 2.0e9 0.35 1.0e-4 1200.0 0.05\n";
        content << "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0 1.0\n"; // до FIBER
        content << "FIBER 0 2\n";
        content << "0.0 0.0 0.0\n";
        content << "1.0 0.0 0.0\n";
        content << "CLAMP_BOX 0.5 -1.0 -1.0 1.5 1.0 1.0\n"; // после FIBER

        auto model = load_from_content(content.str());
        check(model.clamps.fiber_mask[0] == interactions::kClampAll, "директива до FIBER работает");
        check(model.clamps.fiber_mask[1] == interactions::kClampAll, "директива после FIBER работает");
    }

    // Тест 4: бусина на грани бокса закреплена (замкнутый бокс)
    {
        std::ostringstream content;
        content << "MATERIAL 0 2.0e9 0.35 1.0e-4 1200.0 0.05\n";
        content << "FIBER 0 3\n";
        content << "0.0 0.0 0.0\n";
        content << "1.0 0.0 0.0\n";
        content << "2.0 0.0 0.0\n";
        content << "CLAMP_BOX 0.0 -1.0 -1.0 2.0 1.0 1.0\n"; // все три бусины на границе включены

        auto model = load_from_content(content.str());
        check(model.clamps.fiber_mask[0] == interactions::kClampAll, "бусина на lo закреплена");
        check(model.clamps.fiber_mask[1] == interactions::kClampAll, "бусина внутри закреплена");
        check(model.clamps.fiber_mask[2] == interactions::kClampAll, "бусина на hi закреплена");
    }

    // Тест 5: без директив fiber_mask пуст
    {
        std::ostringstream content;
        content << "MATERIAL 0 2.0e9 0.35 1.0e-4 1200.0 0.05\n";
        content << "FIBER 0 2\n";
        content << "0.0 0.0 0.0\n";
        content << "1.0 0.0 0.0\n";

        auto model = load_from_content(content.str());
        check(model.clamps.fiber_mask.empty(), "без директив fiber_mask пуст (не заполнен нулями)");
    }

    // Тест 6: ошибки формата -- std::runtime_error с номером строки директивы.
    {
        const std::string head = "MATERIAL 0 2.0e9 0.35 1.0e-4 1200.0 0.05\n"; // строка 1
        const std::string fib = "FIBER 0 2\n0.0 0.0 0.0\n1.0 0.0 0.0\n";       // строки 2-4
        struct Case
        {
            const char* what;
            std::string text;
            int line;
        };
        const Case cases[] = {
            {"не число в lo", head + fib + "CLAMP_BOX abc -1.0 -1.0 0.5 1.0 1.0\n", 5},
            {"число с мусором", head + fib + "CLAMP_BOX 0.0x -1.0 -1.0 0.5 1.0 1.0\n", 5},
            {"незаконченная экспонента", head + fib + "CLAMP_BOX 1e -1.0 -1.0 0.5 1.0 1.0\n", 5},
            {"nan", head + fib + "CLAMP_BOX nan -1.0 -1.0 0.5 1.0 1.0\n", 5},
            {"inf", head + fib + "CLAMP_BOX 0.0 -1.0 -1.0 inf 1.0 1.0\n", 5},
            {"не хватает hi", head + fib + "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0\n", 5},
            {"xmin > xmax", head + fib + "CLAMP_BOX 1.0 -1.0 -1.0 0.5 1.0 1.0\n", 5},
            {"zmin > zmax", head + fib + "CLAMP_BOX 0.0 -1.0 2.0 0.5 1.0 1.0\n", 5},
            {"недопустимая буква", head + fib + "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0 1.0 xyw\n", 5},
            {"заглавная буква", head + fib + "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0 1.0 X\n", 5},
            {"повтор буквы", head + fib + "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0 1.0 xyx\n", 5},
            {"лишний токен", head + fib + "CLAMP_BOX 0.0 -1.0 -1.0 0.5 1.0 1.0 xy extra\n", 5},
            // директива до FIBER и после комментариев/пустых строк: номер строки файла
            {"ошибка до FIBER", head + "# комментарий\n\nCLAMP_BOX 1.0 0 0 0.0 1 1\n" + fib, 4},
        };
        for (const Case& c : cases)
        {
            const std::string msg = load_error(c.text);
            const std::string want = "строке " + std::to_string(c.line);
            check(!msg.empty() && msg.find(want) != std::string::npos &&
                      msg.find("строке " + std::to_string(c.line) + "0") == std::string::npos,
                  std::string("ошибка CLAMP_BOX (") + c.what + "): runtime_error с '" + want + "'");
        }
    }
}

} // namespace

#include "clamp_integration_solid.hpp" // S1-S4 (используют check/failures выше)
#include "clamp_integration_fiber.hpp" // F1-F5 (используют check/failures выше)

int main()
{
    test_u1();
    test_u2();
    test_u3();
    test_u4();
    test_u5();
    test_validate();
    test_u6_clamp_box();
    // F1-F5, S1-S4: интеграционные тесты через time_integration -- добавляются отдельно.
    run_clamp_fiber_integration_tests(); // F1-F5
    run_clamp_solid_integration_tests(); // S1-S4

    std::printf("\n%s (failures = %d)\n", failures == 0 ? "ALL OK" : "FAILED", failures);
    return failures == 0 ? 0 : 1;
}
