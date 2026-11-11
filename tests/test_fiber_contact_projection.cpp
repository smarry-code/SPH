// tests/test_fiber_contact_projection.cpp
// Unit-тесты стадии нормального контакта как ПРОЕКЦИИ СКОРОСТИ
// (interactions/fiber_normal_projection.hpp; fiber/LAW3_VARIANT1_VELOCITY_PROJECTION.md §7,
// fiber/TZ_LAW3_PROJECTION.md §2 A5).
//
// Реального поиска контактов нет: сцены крошечные, геометрия слота считается
// простыми помощниками (точечные бусины a1 = b1 = -1 либо прямые сегменты с
// ЗАДАННЫМИ параметрами s, t; нормаль и проникание -- аналитически). Шаг сцены
// повторяет арифметику решателя:
//     v^{n+1} = v^n + dt * a,    x^{n+1} = x^n + dt * v^{n+1},
// где a = внешнее ускорение + Δv_стадии/dt; ode::FiberVelocityUpdate{base = v^n, h = dt}.
//
// Разделы:
//   U1..U6  -- формула одного контакта, обход слота и счётчики (+ зажимы), одиночный
//              контакт на сегментах, нить-Solid с XSPH av, тёплый старт, слоты трения;
//   Т1      -- статика под нагрузкой, 1e5 шагов (+ закон 3 для сравнения: дрейф);
//   Т2      -- лобовой удар (+ сравнение с законом 3 при gamma = 1);
//   Т3      -- связанные контакты (узел из 3 нитей, Solid на трёх бусинах, цепочка
//              Solid-нить-нить): таблица остатка для K = 1,2,4,8 с тёплым стартом и без;
//   Т4/Т5   -- энергия стадии и импульс в Т1-Т3 и на случайной сцене;
//   Т6      -- сходимость по шагу (dt_scale 1; 0,5; 0,25): конечные значения и ряды delta(t), v(t) на 20 моментах;
//   Т7      -- побитовая одинаковость на 1 и 8 потоках.
// Собирается с SPH_DIM=3 (CMakeLists.txt).
#include "contact_force.hpp"
#include "fiber/fiber_system.hpp"
#include "interactions/clamp.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_friction.hpp"
#include "interactions/fiber_normal_projection.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "ode/ode_solver.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_contact_projection требует сборки с SPH_DIM=3");

namespace
{

using namespace interactions;

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

double dot3(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm3(const double* a) { return std::sqrt(dot3(a, a)); }

struct Rng
{
    std::uint64_t s;
    explicit Rng(std::uint64_t seed) : s(seed) {}
    std::uint64_t next()
    {
        std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uni() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }
    double range(double a, double b) { return a + (b - a) * uni(); }
    std::size_t index(std::size_t n) { return static_cast<std::size_t>(next() % n); }
    void unit(double (&v)[dim])
    {
        double l;
        do
        {
            for (std::size_t d = 0; d < dim; ++d) v[d] = range(-1.0, 1.0);
            l = norm3(v);
        } while (l < 0.1 || l > 1.0);
        for (std::size_t d = 0; d < dim; ++d) v[d] /= l;
    }
};

bool bits_equal(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

// Физические масштабы сцен (порядок как на сцене Guo 1T): диаметр d, шаг 1e-8 с.
constexpr double kD = 1e-3;     // r_a + r_b, м
constexpr double kM = 1e-6;     // масса бусины, кг
constexpr double kG = 1e8;      // нагрузка (ускорение), м/с^2: g*dt = 1 м/с при dt = 1e-8
constexpr double kDt0 = 1e-8;   // базовый шаг, с
constexpr double kPi = 3.14159265358979323846;

// ============================================================================
// Мини-решатель сцены: бусины (FiberState), частицы Solid (ParticleSystem),
// внешние ускорения, закреплённые элементы (c = 0 через ClampWorkspace), определения
// контактов (параметры s,t,lambda заданы, нормаль и проникание считаются по x^n),
// слоты, тёплый старт, стадия.
// ============================================================================
struct FFDef
{
    rhs::index_t a0;
    long long a1;
    rhs::index_t b0;
    long long b1;
    double s, t;
    double dsum; // r_a + r_b
    double mu;
};

struct FSDef
{
    rhs::index_t a;
    long long b;
    rhs::index_t s;
    double lambda;
    double dsum;
    double mu;
};

struct Sim
{
    rhs::CpuExecutor exec;
    fiber::FiberState fb;
    ParticleSystem sys;
    std::size_t nf = 0, ns = 0;
    std::vector<double> ext_f[dim], ext_s[dim]; // внешние ускорения
    std::vector<std::vector<double>> av;        // XSPH-поправка Solid ([dim][ns]); пусто -- нет
    std::vector<char> fixed_f, fixed_s;
    std::vector<FFDef> ffdef;
    std::vector<FSDef> fsdef;
    double dt = kDt0;
    int K = 8;
    bool warm_on = true;
    bool check_momentum = true;

    std::vector<double> dv_f[dim];
    std::vector<std::vector<double>> dvx_s;
    std::vector<FiberNormalContact> ff;
    std::vector<FiberSolidNormalContact> fs;
    NormalWarmStart warm;
    NormalProjectionWorkspace ws;
    ClampWorkspace cw;
    NormalProjectionStats st;

    // накопители по шагам (Т4, Т5)
    double Wpos_max = 0.0;    // max_n max(0, W_stage)
    double E_ref = 0.0;       // max_n KE пробных скоростей
    double mom_rel_max = 0.0; // max_n |sum m dv| / sum m |dv|
    std::size_t nsteps = 0;

    explicit Sim(int threads) : exec(threads) {}

    void init(std::size_t nfiber, std::size_t nsolid)
    {
        nf = nfiber;
        ns = nsolid;
        fb.resize(nf);
        sys.resize(ns);
        for (std::size_t d = 0; d < dim; ++d)
        {
            ext_f[d].assign(nf, 0.0);
            ext_s[d].assign(ns, 0.0);
            dv_f[d].assign(nf, 0.0);
            for (std::size_t i = 0; i < nf; ++i)
            {
                fb.x[d][i] = 0.0;
                fb.vx[d][i] = 0.0;
            }
            for (std::size_t i = 0; i < ns; ++i)
            {
                sys.x[d][i] = 0.0;
                sys.vx[d][i] = 0.0;
            }
        }
        for (std::size_t i = 0; i < nf; ++i) fb.mass[i] = kM;
        for (std::size_t i = 0; i < ns; ++i) sys.mass[i] = kM;
        dvx_s.assign(dim, std::vector<double>(ns, 0.0));
        fixed_f.assign(nf, 0);
        fixed_s.assign(ns, 0);
    }

    // Закреплённые элементы (c = 0): ClampWorkspace.full по возрастанию.
    void finalize_fixed()
    {
        cw.fiber.full.clear();
        cw.solid.full.clear();
        for (std::size_t i = 0; i < nf; ++i)
            if (fixed_f[i]) cw.fiber.full.push_back(i);
        for (std::size_t i = 0; i < ns; ++i)
            if (fixed_s[i]) cw.solid.full.push_back(i);
    }

    void fiber_point(rhs::index_t k0, long long k1, double w1, double (&p)[dim]) const
    {
        for (std::size_t d = 0; d < dim; ++d)
        {
            p[d] = fb.x[d][k0];
            if (k1 != -1) p[d] = (1.0 - w1) * fb.x[d][k0] + w1 * fb.x[d][static_cast<std::size_t>(k1)];
        }
    }

    // Расстояние между точками контакта нить-нить.
    double ff_distance(const FFDef& f) const
    {
        double pa[dim], pb[dim], df[dim];
        fiber_point(f.a0, f.a1, f.s, pa);
        fiber_point(f.b0, f.b1, f.t, pb);
        for (std::size_t d = 0; d < dim; ++d) df[d] = pa[d] - pb[d];
        return norm3(df);
    }
    double ff_delta(std::size_t i) const { return ffdef[i].dsum - ff_distance(ffdef[i]); }

    double fs_distance(const FSDef& f) const
    {
        double pa[dim], df[dim];
        fiber_point(f.a, f.b, f.lambda, pa);
        for (std::size_t d = 0; d < dim; ++d) df[d] = sys.x[d][f.s] - pa[d];
        return norm3(df);
    }
    double fs_delta(std::size_t i) const { return fsdef[i].dsum - fs_distance(fsdef[i]); }

    // Слоты по x^n: контакт активен при delta > 0.
    void build_slots()
    {
        ff.assign(ffdef.size(), FiberNormalContact{});
        for (std::size_t k = 0; k < ffdef.size(); ++k)
        {
            const FFDef& f = ffdef[k];
            FiberNormalContact& c = ff[k];
            double pa[dim], pb[dim], df[dim];
            fiber_point(f.a0, f.a1, f.s, pa);
            fiber_point(f.b0, f.b1, f.t, pb);
            for (std::size_t d = 0; d < dim; ++d) df[d] = pa[d] - pb[d];
            const double D = norm3(df);
            c.a0 = f.a0;
            c.a1 = f.a1;
            c.b0 = f.b0;
            c.b1 = f.b1;
            c.s = f.s;
            c.t = f.t;
            c.mu = f.mu;
            c.delta = f.dsum - D;
            c.active = (D > 0.0 && c.delta > 0.0);
            if (D > 0.0)
                for (std::size_t d = 0; d < dim; ++d) c.n[d] = df[d] / D;
            c.key = make_ff_pair_key(static_cast<long long>(f.a0), f.a1, static_cast<long long>(f.b0), f.b1);
        }
        fs.assign(fsdef.size(), FiberSolidNormalContact{});
        for (std::size_t k = 0; k < fsdef.size(); ++k)
        {
            const FSDef& f = fsdef[k];
            FiberSolidNormalContact& c = fs[k];
            double pa[dim], df[dim];
            fiber_point(f.a, f.b, f.lambda, pa);
            for (std::size_t d = 0; d < dim; ++d) df[d] = sys.x[d][f.s] - pa[d];
            const double D = norm3(df);
            c.a = f.a;
            c.b = f.b;
            c.s = f.s;
            c.lambda = f.lambda;
            c.mu = f.mu;
            c.delta = f.dsum - D;
            c.active = (D > 0.0 && c.delta > 0.0);
            if (D > 0.0)
                for (std::size_t d = 0; d < dim; ++d) c.n[d] = df[d] / D;
            c.key = make_fs_pair_key(static_cast<long long>(f.a), f.b, static_cast<long long>(f.s));
        }
    }

    // Внешнее ускорение шага без нормального контакта; закреплённые -- 0.
    void prepare_dv()
    {
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < nf; ++i) dv_f[d][i] = fixed_f[i] ? 0.0 : ext_f[d][i];
            for (std::size_t i = 0; i < ns; ++i) dvx_s[d][i] = fixed_s[i] ? 0.0 : ext_s[d][i];
        }
    }

    // Стадия: тёплый старт -> счётчики (+ зажимы) -> проекция -> таблица.
    void stage()
    {
        ode::FiberVelocityUpdate fupd{};
        ode::SolidVelocityUpdate supd{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            fupd.base[d] = &fb.vx[d];
            supd.base[d] = &sys.vx[d];
        }
        fupd.h = dt;
        supd.h = dt;
        supd.av = av.empty() ? nullptr : &av;

        warm.apply(exec, ff, fs, dt, warm_on);
        compute_normal_counts(exec, nf, ns, ff, fs, ws);
        apply_clamp_friction_counts(exec, cw, ws.bead_count, ws.solid_count);

        // КЭ пробных скоростей (масштаб E для допусков энергии)
        double E = 0.0;
        for (std::size_t i = 0; i < nf; ++i)
        {
            if (fixed_f[i]) continue;
            double v2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double v = fb.vx[d][i] + dt * dv_f[d][i];
                v2 += v * v;
            }
            E += 0.5 * fb.mass[i] * v2;
        }
        for (std::size_t i = 0; i < ns; ++i)
        {
            if (fixed_s[i]) continue;
            double v2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                double v = sys.vx[d][i] + dt * dvx_s[d][i];
                if (!av.empty()) v += av[d][i];
                v2 += v * v;
            }
            E += 0.5 * sys.mass[i] * v2;
        }
        E_ref = std::max(E_ref, E);

        run_normal_projection(exec, fb, fupd, dv_f, ff, sys, supd, dvx_s, fs, K, ws, st);
        warm.rebuild(ff, fs, dt);

        Wpos_max = std::max(Wpos_max, st.W_stage);

        // импульс: sum m Δv стадии
        double p[dim] = {0.0, 0.0, 0.0};
        double den = 0.0;
        for (std::size_t i = 0; i < nf; ++i)
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double a = ws.dvel_fiber[d][i];
                p[d] += fb.mass[i] * a;
                den += fb.mass[i] * std::fabs(a);
            }
        for (std::size_t i = 0; i < ns; ++i)
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double a = ws.dvel_solid[d][i];
                p[d] += sys.mass[i] * a;
                den += sys.mass[i] * std::fabs(a);
            }
        if (den > 0.0) mom_rel_max = std::max(mom_rel_max, norm3(p) / den);
        ++nsteps;
    }

    void advance()
    {
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < nf; ++i)
            {
                if (fixed_f[i]) continue;
                fb.vx[d][i] += dt * dv_f[d][i];
                fb.x[d][i] += dt * fb.vx[d][i];
            }
            for (std::size_t i = 0; i < ns; ++i)
            {
                if (fixed_s[i]) continue;
                sys.vx[d][i] += dt * dvx_s[d][i];
                if (!av.empty()) sys.vx[d][i] += av[d][i];
                sys.x[d][i] += dt * sys.vx[d][i];
            }
        }
    }

    void step()
    {
        build_slots();
        prepare_dv();
        stage();
        advance();
    }

    // Внешние ускорения, уравновешивающие заданные N_c нагрузки: на каждую бусину/частицу
    // -(сумма сил контакта на неё), где сила контакта на элемент = sigma w N n (A получает +N n).
    void set_equilibrium_loads(const std::vector<double>& Nff, const std::vector<double>& Nfs)
    {
        build_slots();
        for (std::size_t k = 0; k < ff.size(); ++k)
        {
            const FiberNormalContact& c = ff[k];
            if (!c.active) continue;
            fiber_normal_for_each_bead(c, [&](rhs::index_t b, double w, double sg, std::uint16_t) {
                for (std::size_t d = 0; d < dim; ++d) ext_f[d][b] -= sg * w * Nff[k] * c.n[d] / fb.mass[b];
            });
        }
        for (std::size_t k = 0; k < fs.size(); ++k)
        {
            const FiberSolidNormalContact& c = fs[k];
            if (!c.active) continue;
            fiber_solid_normal_for_each_bead(c, [&](rhs::index_t b, double w, double sg, std::uint16_t) {
                for (std::size_t d = 0; d < dim; ++d) ext_f[d][b] -= sg * w * Nfs[k] * c.n[d] / fb.mass[b];
            });
            for (std::size_t d = 0; d < dim; ++d) ext_s[d][c.s] -= Nfs[k] * c.n[d] / sys.mass[c.s];
        }
    }
};

// ============================================================================
// Сцены Т1 и Т3.
// ============================================================================

// Т1: бусина (variant 0) либо сегмент из двух бусин (variant 1, контакт в s = 0.5) лежит на
// закреплённой бусине под нагрузкой kG. Начальное проникание delta0 = 0.1 d.
void make_static_scene(Sim& s, int variant, double delta0)
{
    if (variant == 0)
    {
        s.init(2, 0);
        s.fixed_f[0] = 1;
        s.fb.x[2][1] = kD - delta0;
        s.ext_f[2][1] = -kG;
        s.ffdef.push_back(FFDef{1, -1, 0, -1, 0.0, 0.0, kD, 0.3});
    }
    else
    {
        s.init(3, 0);
        s.fixed_f[0] = 1;
        const double L = 2.0 * kD;
        s.fb.x[0][1] = -L;
        s.fb.x[0][2] = L;
        s.fb.x[2][1] = kD - delta0;
        s.fb.x[2][2] = kD - delta0;
        s.ext_f[2][1] = -kG;
        s.ext_f[2][2] = -kG;
        s.ffdef.push_back(FFDef{1, 2, 0, -1, 0.5, 0.0, kD, 0.3});
    }
    s.finalize_fixed();
}

// Узел из трёх нитей: центральная бусина 0, соседи 1..3 вокруг неё (3D, не в плоскости),
// контакты A = сосед, B = центр. Нагрузки N_c = N0*(1,2,3) (разные).
void make_knot(Sim& s, double N0)
{
    s.init(4, 0);
    const double delta0 = 0.05 * kD;
    const double ang[3] = {0.5 * kPi, 0.5 * kPi + 2.0 * kPi / 3.0, 0.5 * kPi + 4.0 * kPi / 3.0};
    for (int i = 0; i < 3; ++i)
    {
        double u[dim] = {std::cos(ang[i]), std::sin(ang[i]), 0.3};
        const double l = norm3(u);
        for (std::size_t d = 0; d < dim; ++d) s.fb.x[d][static_cast<std::size_t>(i + 1)] = (kD - delta0) * u[d] / l;
        s.ffdef.push_back(FFDef{static_cast<rhs::index_t>(i + 1), -1, 0, -1, 0.0, 0.0, kD, 0.3});
    }
    s.set_equilibrium_loads({N0, 2.0 * N0, 3.0 * N0}, {});
}

// Solid-частица (масса 4 m) на трёх бусинах кольца: три контакта нить-Solid (точечные), N_c = N0*(1,2,3).
void make_tripod(Sim& s, double N0)
{
    s.init(3, 1);
    s.sys.mass[0] = 4.0 * kM;
    const double delta0 = 0.05 * kD;
    const double rho = 0.8 * kD;
    const double H = std::sqrt((kD - delta0) * (kD - delta0) - rho * rho);
    s.sys.x[2][0] = H;
    for (int i = 0; i < 3; ++i)
    {
        const double a = 2.0 * kPi * i / 3.0;
        s.fb.x[0][static_cast<std::size_t>(i)] = rho * std::cos(a);
        s.fb.x[1][static_cast<std::size_t>(i)] = rho * std::sin(a);
        s.fb.x[2][static_cast<std::size_t>(i)] = 0.0;
        s.fsdef.push_back(FSDef{static_cast<rhs::index_t>(i), -1, 0, 0.0, kD, 0.3});
    }
    s.set_equilibrium_loads({}, {N0, 2.0 * N0, 3.0 * N0});
}

// Цепочка: Solid S (масса 4 m) -- сегмент (бусины 0 и 2, точка контакта lambda = s = 0.3) -- бусина 1.
// Контакты: нить-Solid (сегмент -- S) и нить-нить (сегмент -- бусина 1). Бусины сегмента в двух
// слотах каждая (c = 2). N_fs = N0, N_ff = 2 N0.
void make_chain(Sim& s, double N0)
{
    s.init(3, 1);
    s.sys.mass[0] = 4.0 * kM;
    const double delta0 = 0.05 * kD;
    const double z1 = kD - delta0;
    const double L = 2.0 * kD;
    const double xc = -L + 0.3 * 2.0 * L; // точка контакта на сегменте
    s.fb.x[0][0] = -L;
    s.fb.x[2][0] = z1;
    s.fb.x[0][2] = L;
    s.fb.x[2][2] = z1;
    s.fb.x[0][1] = xc;
    s.fb.x[2][1] = 0.0;
    s.sys.x[0][0] = xc;
    s.sys.x[2][0] = 2.0 * z1;
    s.ffdef.push_back(FFDef{0, 2, 1, -1, 0.3, 0.0, kD, 0.3});
    s.fsdef.push_back(FSDef{0, 2, 0, 0.3, kD, 0.3});
    s.set_equilibrium_loads({2.0 * N0}, {N0});
}

// ============================================================================
// Закон 3 (нынешний) для сравнения: одномерные сцены.
// F = gamma m_red pdot_pred / dt (contact::impulse_force), pdot_pred -- по скорости
// предиктора v~ = v^n + (dt/2) a^{n-1}.
// ============================================================================
struct Law3Static
{
    double delta_end = 0.0, N_end = 0.0;
};

// Бусина массы m над закреплённой бусиной (ось z), нагрузка g вниз.
Law3Static law3_static(double dt, std::size_t steps, double delta0)
{
    double z = kD - delta0, v = 0.0, a_prev = -kG, N = 0.0;
    for (std::size_t i = 0; i < steps; ++i)
    {
        const double vp = v + 0.5 * dt * a_prev;
        contact::PairGeom pg;
        pg.r = z;
        pg.p = kD - z;
        pg.pdot = -vp;
        pg.v_rel = std::fabs(vp);
        contact::PairProps pr;
        pr.m_red = kM;
        const contact::Result res = contact::impulse_force(pg, pr, dt, 1.0);
        const double a = -kG + res.fn / kM;
        v += dt * a;
        z += dt * v;
        a_prev = a;
        N = res.fn;
    }
    return {kD - z, N};
}

// Лобовой удар двух бусин по закону 3 (gamma = 1); ось x, A слева. Возвращает скорости после steps шагов.
void law3_impact(double mA, double mB, double vA0, double vB0, double dt, std::size_t steps, double delta0,
                 double& vA, double& vB)
{
    double xA = -0.5 * (kD - delta0), xB = 0.5 * (kD - delta0);
    vA = vA0;
    vB = vB0;
    double aA_prev = 0.0, aB_prev = 0.0;
    for (std::size_t i = 0; i < steps; ++i)
    {
        const double vAp = vA + 0.5 * dt * aA_prev, vBp = vB + 0.5 * dt * aB_prev;
        contact::PairGeom pg;
        pg.r = xB - xA;
        pg.p = kD - pg.r;
        pg.pdot = (vAp - vBp); // A слева: сближение при vA > vB
        pg.v_rel = std::fabs(vAp - vBp);
        contact::PairProps pr;
        pr.m_red = mA * mB / (mA + mB);
        const contact::Result res = contact::impulse_force(pg, pr, dt, 1.0);
        const double aA = -res.fn / mA, aB = res.fn / mB;
        vA += dt * aA;
        vB += dt * aB;
        xA += dt * vA;
        xB += dt * vB;
        aA_prev = aA;
        aB_prev = aB;
    }
}

// ============================================================================
// U1. Формула одного контакта.
// ============================================================================
void test_impulse_formula()
{
    std::printf("[U1] normal_projection_impulse\n");
    double Pn, dP;
    normal_projection_impulse(2.0, 4.0, 1.0, Pn, dP);
    check(Pn == 1.5 && dP == 0.5, "P_new = P_old + pdot/W (1.5), dP = 0.5");
    normal_projection_impulse(-8.0, 4.0, 1.0, Pn, dP);
    check(Pn == 0.0 && dP == -1.0, "накопленный импульс зажимается в 0, приращение -- нет (dP = -P_old)");
    normal_projection_impulse(-1.0, 4.0, 0.0, Pn, dP);
    check(Pn == 0.0 && dP == 0.0, "нет сближения и P_old = 0 -- импульса нет");
    normal_projection_impulse(5.0, 0.0, 0.7, Pn, dP);
    check(Pn == 0.7 && dP == 0.0, "W = 0 (бесконечная масса) -- P не меняется");
    normal_projection_impulse(5.0, -1.0, 0.7, Pn, dP);
    check(Pn == 0.7 && dP == 0.0, "W < 0 -- P не меняется");
}

// ============================================================================
// U2. Обход слота и счётчики c_k, c_s (+ зажимы).
// ============================================================================
void test_counts_and_traversal()
{
    std::printf("[U2] обход слота, нулевые веса, счётчики c_k/c_s, зажимы\n");
    FiberNormalContact c;
    c.active = true;
    c.a0 = 0;
    c.a1 = 1;
    c.s = 1.0; // вес a0 = 0
    c.b0 = 2;
    c.b1 = 3;
    c.t = 0.0; // вес b1 = 0
    int n = 0;
    bool saw0 = false, saw3 = false;
    fiber_normal_for_each_bead(c, [&](rhs::index_t b, double w, double sg, std::uint16_t) {
        ++n;
        if (b == 0) saw0 = true;
        if (b == 3) saw3 = true;
        if (b == 1) check(w == 1.0 && sg == 1.0, "a1: вес s = 1, знак +1");
        if (b == 2) check(w == 1.0 && sg == -1.0, "b0: вес 1 - t = 1, знак -1");
    });
    check(n == 2 && !saw0 && !saw3, "бусины с нулевым весом не обходятся");

    FiberSolidNormalContact fc;
    fc.active = true;
    fc.a = 4;
    fc.b = 5;
    fc.lambda = 0.25;
    double wsum = 0.0;
    int nn = 0;
    fiber_solid_normal_for_each_bead(fc, [&](rhs::index_t, double w, double sg, std::uint16_t) {
        wsum += w;
        ++nn;
        check(sg == -1.0, "нить-Solid: знак нити -1 (Solid -- сторона A)");
    });
    check(nn == 2 && std::fabs(wsum - 1.0) < 1e-15, "нить-Solid: веса (1-lambda, lambda), сумма 1");

    // сцена: 6 бусин, 3 частицы Solid
    Sim s(1);
    s.init(6, 3);
    auto mk_ff = [](rhs::index_t a0, long long a1, double sa, rhs::index_t b0, long long b1, double tb, bool act) {
        FiberNormalContact x;
        x.active = act;
        x.a0 = a0;
        x.a1 = a1;
        x.s = sa;
        x.b0 = b0;
        x.b1 = b1;
        x.t = tb;
        return x;
    };
    auto mk_fs = [](rhs::index_t a, long long b, double lam, rhs::index_t sp) {
        FiberSolidNormalContact x;
        x.active = true;
        x.a = a;
        x.b = b;
        x.lambda = lam;
        x.s = sp;
        return x;
    };
    s.ff = {mk_ff(0, 1, 0.5, 2, -1, 0.0, true), mk_ff(2, 3, 1.0, 4, 5, 0.3, true), mk_ff(0, 1, 0.5, 2, -1, 0.0, false)};
    s.fs = {mk_fs(1, 2, 0.5, 0), mk_fs(2, -1, 0.0, 0), mk_fs(5, 4, 0.0, 1)};
    compute_normal_counts(s.exec, 6, 3, s.ff, s.fs, s.ws);
    const double expect_b[6] = {1, 2, 3, 1, 1, 2};
    bool ok = true;
    for (int i = 0; i < 6; ++i) ok = ok && (s.ws.bead_count[static_cast<std::size_t>(i)] == expect_b[i]);
    check(ok, "c_k: [1,2,3,1,1,2] (нулевые веса и неактивный слот не считаются)");
    check(s.ws.solid_count[0] == 2.0 && s.ws.solid_count[1] == 1.0 && s.ws.solid_count[2] == 0.0,
          "c_s: [2,1,0]");

    s.fixed_f[1] = 1;
    s.fixed_s[0] = 1;
    s.finalize_fixed();
    apply_clamp_friction_counts(s.exec, s.cw, s.ws.bead_count, s.ws.solid_count);
    check(s.ws.bead_count[1] == 0.0 && s.ws.solid_count[0] == 0.0 && s.ws.bead_count[2] == 3.0,
          "зажим: c = 0 у закреплённых бусины и частицы, остальные не тронуты");
}

// ============================================================================
// U3. Одиночный контакт нить-нить на сегментах: точное решение за один проход.
// ============================================================================
void test_single_ff()
{
    std::printf("[U3] одиночный контакт нить-нить (сегменты), K = 1\n");
    for (int flip = 0; flip < 2; ++flip)
    {
        Rng rng(77 + static_cast<std::uint64_t>(flip));
        Sim s(1);
        s.init(4, 0);
        s.warm_on = false;
        s.K = 1;
        for (std::size_t i = 0; i < 4; ++i) s.fb.mass[i] = rng.range(0.5, 3.0) * kM;
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < 4; ++i)
            {
                s.fb.vx[d][i] = rng.range(-5.0, 5.0);
                s.ext_f[d][i] = rng.range(-1.0, 1.0) * kG;
            }
        FiberNormalContact c;
        c.active = true;
        c.a0 = 0;
        c.a1 = 1;
        c.s = 0.37;
        c.b0 = 2;
        c.b1 = 3;
        c.t = 0.62;
        rng.unit(c.n);
        s.ff = {c};
        s.prepare_dv();
        double vtrial[dim][4];
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < 4; ++i) vtrial[d][i] = s.fb.vx[d][i] + s.dt * s.dv_f[d][i];
        auto pdot_of = [&](auto vfun) {
            double rel = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double vA = (1.0 - c.s) * vfun(d, 0) + c.s * vfun(d, 1);
                const double vB = (1.0 - c.t) * vfun(d, 2) + c.t * vfun(d, 3);
                rel += c.n[d] * (vA - vB);
            }
            return -rel;
        };
        // сделать сближение (flip = 0) либо разлёт (flip = 1)
        double pd = pdot_of([&](std::size_t d, std::size_t i) { return vtrial[d][i]; });
        if ((pd > 0.0) == (flip == 1))
        {
            for (std::size_t d = 0; d < dim; ++d) s.ff[0].n[d] = -s.ff[0].n[d];
            c = s.ff[0];
            pd = -pd;
        }
        s.stage();
        const double W = (1.0 - c.s) * (1.0 - c.s) / s.fb.mass[0] + c.s * c.s / s.fb.mass[1] +
                         (1.0 - c.t) * (1.0 - c.t) / s.fb.mass[2] + c.t * c.t / s.fb.mass[3];
        const double pd_after = pdot_of([&](std::size_t d, std::size_t i) { return s.fb.vx[d][i] + s.dt * s.dv_f[d][i]; });
        if (flip == 0)
        {
            check(pd > 0.0, "контакт сближается");
            check(std::fabs(s.ff[0].P / (pd / W) - 1.0) < 1e-12, "P = pdot/W с точностью округления");
            check(std::fabs(pd_after) < 1e-11, "после шага скорость сближения = 0 (< 1e-11 м/с)");
            check(s.st.resid < 1e-11 && s.st.n_slots_ff == 1, "stats: остаток ~0, один слот");
        }
        else
        {
            check(pd < 0.0, "контакт разлетается");
            bool unchanged = true;
            for (std::size_t d = 0; d < dim; ++d)
                for (std::size_t i = 0; i < 4; ++i) unchanged = unchanged && (s.dv_f[d][i] == s.ext_f[d][i]);
            check(s.ff[0].P == 0.0 && unchanged, "нет сближения: P = 0, ускорения не тронуты (побитово)");
        }
    }

    // нулевой вес: s = 1 -- бусина a0 не участвует
    Sim s(1);
    s.init(3, 0);
    s.warm_on = false;
    s.K = 1;
    s.fb.vx[0][1] = 3.0;
    s.fb.vx[0][2] = -3.0;
    FiberNormalContact c;
    c.active = true;
    c.a0 = 0;
    c.a1 = 1;
    c.s = 1.0;
    c.b0 = 2;
    c.b1 = -1;
    c.t = 0.0;
    c.n[0] = -1.0; // B -> A: A справа от... сближаются: vA - vB = 6 вдоль x, n = -x => pdot = +6
    s.ff = {c};
    s.prepare_dv();
    s.stage();
    check(s.dv_f[0][0] == 0.0 && s.dv_f[1][0] == 0.0 && s.dv_f[2][0] == 0.0, "бусина с нулевым весом не получает Δv");
    check(s.ws.bead_count[0] == 0.0 && s.ws.bead_count[1] == 1.0, "её счётчик c = 0");
    check(std::fabs(s.fb.vx[0][1] + s.dt * s.dv_f[0][1]) < 1e-12 && std::fabs(s.fb.vx[0][2] + s.dt * s.dv_f[0][2]) < 1e-12,
          "равные массы: обе скорости -> 0 (e = 0)");
}

// ============================================================================
// U4. Нить-Solid с XSPH av: пробная скорость частицы включает av.
// ============================================================================
void test_fs_av()
{
    std::printf("[U4] нить-Solid: пробная скорость частицы с XSPH av\n");
    Sim s(1);
    s.init(1, 1);
    s.warm_on = false;
    s.K = 1;
    s.sys.mass[0] = 3.0 * kM;
    s.sys.x[2][0] = 0.9 * kD;
    s.av.assign(dim, std::vector<double>(1, 0.0));
    s.av[0][0] = 0.3;
    s.av[2][0] = -2.0; // частица сближается с бусиной за счёт av
    s.fb.vx[0][0] = 0.1;
    FiberSolidNormalContact c;
    c.active = true;
    c.a = 0;
    c.b = -1;
    c.s = 0;
    c.lambda = 0.0;
    c.n[2] = 1.0; // от нити к Solid
    s.fs = {c};
    s.prepare_dv();
    s.stage();
    const double W = 1.0 / s.fb.mass[0] + 1.0 / s.sys.mass[0];
    const double pd = 2.0; // -n.(v_s - v_q) = 2 - 0
    check(std::fabs(s.fs[0].P / (pd / W) - 1.0) < 1e-12, "P = pdot/W, pdot учитывает av = -2");
    const double vs = s.sys.vx[2][0] + s.dt * s.dvx_s[2][0] + s.av[2][0];
    const double vq = s.fb.vx[2][0] + s.dt * s.dv_f[2][0];
    check(std::fabs(vs - vq) < 1e-12, "после шага (с av) нормальные скорости совпадают");
    check(s.dvx_s[2][0] > 0.0 && s.dv_f[2][0] < 0.0, "Solid отталкивается вдоль n, нить -- против");
    check(s.dv_f[0][0] == 0.0 && s.dvx_s[0][0] == 0.0, "касательных компонент нет");
}

// ============================================================================
// U5. Тёплый старт: таблица, ключи, масштаб h/h_prev, дубликаты, экспорт/импорт.
// ============================================================================
void test_warm_start()
{
    std::printf("[U5] NormalWarmStart\n");
    rhs::CpuExecutor exec(1);
    std::vector<FiberNormalContact> ff(5);
    std::vector<FiberSolidNormalContact> fs(2);
    auto setff = [&](std::size_t i, long long a0, long long b0, double P, bool act) {
        ff[i].active = act;
        ff[i].a0 = static_cast<rhs::index_t>(a0);
        ff[i].a1 = a0 + 1;
        ff[i].b0 = static_cast<rhs::index_t>(b0);
        ff[i].b1 = b0 + 1;
        ff[i].key = make_ff_pair_key(a0, a0 + 1, b0, b0 + 1);
        ff[i].P = P;
    };
    setff(0, 10, 4, 5.0, true);  // a0 > b0 -- ключ канонизируется (стороны переставляются)
    setff(1, 2, 8, 7.0, true);
    setff(2, 2, 8, 9.0, true);   // повтор ключа -- остаётся первая запись (7.0)
    setff(3, 20, 30, 0.0, true); // P = 0 -- не попадает в таблицу
    setff(4, 40, 50, 3.0, false); // неактивный -- не попадает
    fs[0].active = true;
    fs[0].key = make_fs_pair_key(3, 4, 11);
    fs[0].P = 2.0;
    fs[1].active = true;
    fs[1].key = make_fs_pair_key(1, -1, 7);
    fs[1].P = 4.0;
    NormalWarmStart w;
    w.rebuild(ff, fs, 1e-8);
    check(w.size() == 4 && w.h_prev() == 1e-8, "в таблице 4 записи (без P=0, неактивных и повторов)");
    bool sorted = true;
    for (std::size_t i = 1; i < w.entries().size(); ++i) sorted = sorted && (w.entries()[i - 1].key < w.entries()[i].key);
    check(sorted, "таблица строго отсортирована по ключу");

    // новые слоты: те же пары (одна со стороны, переставленной местами)
    std::vector<FiberNormalContact> ff2(4);
    std::vector<FiberSolidNormalContact> fs2(3);
    auto setkey = [&](std::size_t i, long long a0, long long b0) {
        ff2[i].active = true;
        ff2[i].key = make_ff_pair_key(a0, a0 + 1, b0, b0 + 1);
        ff2[i].P = 123.0; // должен быть перезаписан
    };
    setkey(0, 4, 10); // тот же ключ, что (10,4)
    setkey(1, 2, 8);
    setkey(2, 99, 100); // нет в таблице
    setkey(3, 2, 8);
    ff2[3].active = false; // неактивный: P = 0
    fs2[0].active = true;
    fs2[0].key = make_fs_pair_key(3, 4, 11);
    fs2[1].active = true;
    fs2[1].key = make_fs_pair_key(3, 4, 12); // другой s
    fs2[2].active = true;
    fs2[2].key = make_fs_pair_key(1, -1, 7);
    w.apply(exec, ff2, fs2, 0.5e-8, true);
    check(ff2[0].warm && ff2[0].P == 5.0 * 0.5, "канонический ключ найден, P = P_prev*h/h_prev");
    check(ff2[1].warm && ff2[1].P == 7.0 * 0.5, "дубликат: взята первая запись");
    check(!ff2[2].warm && ff2[2].P == 0.0, "пара не найдена -- холодный старт");
    check(!ff2[3].warm && ff2[3].P == 0.0, "неактивный слот обнуляется");
    check(fs2[0].warm && fs2[0].P == 1.0 && !fs2[1].warm && fs2[1].P == 0.0 && fs2[2].warm && fs2[2].P == 2.0,
          "нить-Solid: ключи различают частицу");
    w.apply(exec, ff2, fs2, 0.5e-8, false);
    check(!ff2[0].warm && ff2[0].P == 0.0 && !fs2[0].warm && fs2[0].P == 0.0, "enabled = false: всё холодное");

    // экспорт / импорт
    NormalWarmStart w2;
    std::vector<NormalWarmEntry> exported = w.entries();
    std::reverse(exported.begin(), exported.end());
    w2.set(exported, w.h_prev());
    bool same = (w2.size() == w.size()) && (w2.h_prev() == w.h_prev());
    for (std::size_t i = 0; same && i < w.size(); ++i)
        same = (w2.entries()[i].key == w.entries()[i].key) && (w2.entries()[i].P == w.entries()[i].P);
    check(same, "импорт переупорядоченного массива даёт ту же таблицу");
    w.clear();
    check(w.size() == 0 && w.h_prev() == 0.0, "clear");
}

// ============================================================================
// U6. Слоты трения из нормальных.
// ============================================================================
void test_friction_fill()
{
    std::printf("[U6] слоты трения из нормальных\n");
    rhs::CpuExecutor exec(1);
    std::vector<FiberNormalContact> ff(4);
    std::vector<FiberSolidNormalContact> fs(3);
    for (std::size_t k = 0; k < 4; ++k)
    {
        ff[k].active = true;
        ff[k].prim = k;
        ff[k].a0 = static_cast<rhs::index_t>(k);
        ff[k].b0 = static_cast<rhs::index_t>(k + 10);
        ff[k].a1 = -1;
        ff[k].b1 = -1;
        ff[k].n[1] = 1.0;
        ff[k].P = 2e-6;
        ff[k].mu = 0.27;
    }
    ff[1].P = 0.0;      // нет импульса
    ff[2].mu = 0.0;     // нет трения
    ff[3].active = false;
    for (std::size_t k = 0; k < 3; ++k)
    {
        fs[k].active = true;
        fs[k].prim = k;
        fs[k].a = static_cast<rhs::index_t>(k);
        fs[k].b = -1;
        fs[k].s = static_cast<rhs::index_t>(k + 5);
        fs[k].n[2] = -1.0;
        fs[k].P = 3e-6;
        fs[k].mu = 0.27;
        fs[k].lambda = 0.0;
    }
    fs[2].P = 0.0;
    std::vector<FiberFrictionContact> fr;
    std::vector<FiberSolidFrictionContact> sfr;
    const bool any = fill_friction_slots_from_normal(exec, 1e-8, ff, fs, fr, sfr);
    check(any && fr.size() == 4 && sfr.size() == 3, "есть активные слоты; размеры совпадают с нормальными");
    check(fr[0].active && !fr[1].active && !fr[2].active && !fr[3].active,
          "нить-нить активен только при active && P>0 && mu>0");
    check(std::fabs(fr[0].N - 200.0) < 1e-9 && fr[0].mu == 0.27 && fr[0].a0 == 0 && fr[0].b0 == 10 && fr[0].n[1] == 1.0,
          "N = P/h, геометрия и mu скопированы");
    check(sfr[0].active && sfr[1].active && !sfr[2].active && std::fabs(sfr[1].N - 300.0) < 1e-9 && sfr[1].s == 6 &&
              sfr[1].n[2] == -1.0,
          "нить-Solid: то же");
    // Слоты трения на шаге сбрасывает ядро контакта (все неактивны, по primitive) -- так и здесь.
    fr.assign(4, FiberFrictionContact{});
    sfr.assign(3, FiberSolidFrictionContact{});
    for (auto& f : ff) f.P = 0.0;
    for (auto& f : fs) f.P = 0.0;
    check(!fill_friction_slots_from_normal(exec, 1e-8, ff, fs, fr, sfr), "нет нагрузки -- возвращает false");
    check(!fr[0].active && !sfr[0].active && fr.size() == 4 && sfr.size() == 3,
          "без нагрузки слоты не записываются, размеры не меняются");
    // Сжатый массив нормальных слотов: слот пишется на место своего primitive (prim), остальные не трогаются.
    std::vector<FiberNormalContact> ffc(1, ff[0]);
    ffc[0].prim = 6;
    ffc[0].P = 2e-6;
    fr.assign(4, FiberFrictionContact{});
    check(fill_friction_slots_from_normal_ff(exec, 1e-8, ffc, fr) && fr.size() == 7 && fr[6].active &&
              std::fabs(fr[6].N - 200.0) < 1e-9 && !fr[0].active && !fr[5].active,
          "сжатый массив: слот трения на месте prim, массив дополнен до prim + 1");
}

// ============================================================================
// Т1. Статика под постоянной нагрузкой.
// ============================================================================
struct T1Result
{
    double delta0 = 0.0, delta_end = 0.0, N_end = 0.0, P_load = 0.0;
    double Wpos_rel = 0.0, mom_rel = 0.0;
};

T1Result run_t1(int variant, double dt, std::size_t steps, int K, bool warm)
{
    Sim s(1);
    s.dt = dt;
    s.K = K;
    s.warm_on = warm;
    s.check_momentum = false; // опора закреплена -- реакция не учитывается
    make_static_scene(s, variant, 0.1 * kD);
    T1Result r;
    r.delta0 = s.ff_delta(0);
    r.P_load = (variant == 0 ? 1.0 : 2.0) * kM * kG;
    for (std::size_t i = 0; i < steps; ++i) s.step();
    r.delta_end = s.ff_delta(0);
    r.N_end = s.ff[0].P / dt;
    r.Wpos_rel = s.E_ref > 0.0 ? s.Wpos_max / s.E_ref : 0.0;
    return r;
}

void test_t1()
{
    std::printf("[Т1] статика: бусина/сегмент на закреплённой бусине, нагрузка kG, 1e5 шагов\n");
    const std::size_t steps = 100000;
    for (int variant = 0; variant < 2; ++variant)
    {
        const T1Result r = run_t1(variant, kDt0, steps, 8, true);
        std::printf("    variant %d: delta0/d = %.6f, delta_end/d = %.9f, N/P_load = %.9f, W_pos/E = %.3e\n", variant,
                    r.delta0 / kD, r.delta_end / kD, r.N_end / r.P_load, r.Wpos_rel);
        check(std::fabs(r.delta_end - r.delta0) < 1e-3 * kD, std::string("проникание не растёт: |Δδ| < 1e-3 d (variant ") +
                                                              std::to_string(variant) + ")");
        check(std::fabs(r.N_end / r.P_load - 1.0) < 1e-3,
              std::string("|N/P_load - 1| < 1e-3 (variant ") + std::to_string(variant) + ")");
        check(r.Wpos_rel <= 1e-12, std::string("Т4: W_stage <= 1e-12 E на каждом шаге (variant ") +
                                       std::to_string(variant) + ")");
    }

    // закон 3 для сравнения: тот же критерий должен его ловить
    const Law3Static l3 = law3_static(kDt0, steps, 0.1 * kD);
    const double drift = (l3.delta_end - 0.1 * kD) / kD;
    std::printf("    закон 3 (gamma = 1): delta_end/d = %.4f (рост %.4f d), N/P_load = %.4f\n", l3.delta_end / kD, drift,
                l3.N_end / (kM * kG));
    check(std::fabs(l3.delta_end - 0.1 * kD) > 1e-3 * kD, "критерий Т1 ловит дрейф закона 3");
    check(drift > 0.5, "закон 3: проникание растёт более чем на 0.5 d за 1e5 шагов");
}

// ============================================================================
// Т2. Лобовой удар двух бусин.
// ============================================================================
void test_t2()
{
    std::printf("[Т2] лобовой удар двух бусин (e = 0)\n");
    const double v0 = 79.0;
    for (int heavy = 0; heavy < 2; ++heavy)
    {
        Sim s(1);
        s.init(2, 0);
        const double mA = kM * (heavy ? 5.0 : 1.0), mB = kM;
        s.fb.mass[0] = mA;
        s.fb.mass[1] = mB;
        const double delta0 = 1e-3 * kD;
        s.fb.x[0][0] = -0.5 * (kD - delta0);
        s.fb.x[0][1] = 0.5 * (kD - delta0);
        s.fb.vx[0][0] = 0.5 * v0;
        s.fb.vx[0][1] = -0.5 * v0;
        s.ffdef.push_back(FFDef{0, -1, 1, -1, 0.0, 0.0, kD, 0.3});
        s.K = 8;
        const double p_before = mA * s.fb.vx[0][0] + mB * s.fb.vx[0][1];
        s.step();
        const double rel1 = (s.fb.vx[0][0] - s.fb.vx[0][1]) / v0;
        const double p_after = mA * s.fb.vx[0][0] + mB * s.fb.vx[0][1];
        double vA3, vB3;
        law3_impact(mA, mB, 0.5 * v0, -0.5 * v0, kDt0, 1, delta0, vA3, vB3);
        for (int i = 1; i < 10; ++i) s.step();
        const double rel10 = (s.fb.vx[0][0] - s.fb.vx[0][1]) / v0;
        double vA3_10, vB3_10;
        law3_impact(mA, mB, 0.5 * v0, -0.5 * v0, kDt0, 10, delta0, vA3_10, vB3_10);
        std::printf("    mA/mB = %.0f: rel после 1 шага = %.3e v0, после 10 = %.3e v0; vs закон 3: |ΔvA| = %.3e, |ΔvB| = %.3e (после 10 шагов)\n",
                    mA / mB, rel1, rel10, std::fabs(s.fb.vx[0][0] - vA3_10), std::fabs(s.fb.vx[0][1] - vB3_10));
        check(std::fabs(rel1) < 1e-10 && std::fabs(rel10) < 1e-10,
              std::string("относительная нормальная скорость после удара < 1e-10 v0 (mA/mB = ") +
                  std::to_string(static_cast<int>(mA / mB)) + ")");
        check(std::fabs(p_after - p_before) < 1e-12 * (mA + mB) * v0, "импульс пары сохраняется");
        // «До округления»: прежний допуск 1e-9 v0 был на 5 порядков шире фактического отличия (7e-15 м/с,
        // ревью Codex); 1e-12 v0 -- около 4500 ulp от v0 = 79 м/с, с запасом на 10 шагов накопления.
        check(std::fabs(s.fb.vx[0][0] - vA3_10) < 1e-12 * v0 && std::fabs(s.fb.vx[0][1] - vB3_10) < 1e-12 * v0,
              "скорости после удара совпадают с законом 3 (gamma = 1) до округления (< 1e-12 v0)");
        (void)vA3;
        (void)vB3;
        check(s.Wpos_max <= 1e-12 * s.E_ref, "Т4: W_stage <= 1e-12 E");
        check(s.mom_rel_max < 1e-12, "Т5: |sum m a_stage| на уровне округления");
        // энергия: удар неупругий, потери = 1/2 m_red v0^2
        const double Ek0 = 0.5 * (mA * 0.25 + mB * 0.25) * v0 * v0;
        const double Ek1 = 0.5 * mA * s.fb.vx[0][0] * s.fb.vx[0][0] + 0.5 * mB * s.fb.vx[0][1] * s.fb.vx[0][1];
        const double mred = mA * mB / (mA + mB);
        check(std::fabs((Ek0 - Ek1) - 0.5 * mred * v0 * v0) < 1e-9 * Ek0, "потеря кинетической энергии = m_red v0^2 / 2");
    }

    // удар из зазора: контакт включается, когда delta > 0 на x^n
    Sim s(1);
    s.init(2, 0);
    s.fb.x[0][0] = -0.5 * (kD + 3.3 * v0 * kDt0);
    s.fb.x[0][1] = 0.5 * (kD + 3.3 * v0 * kDt0);
    s.fb.vx[0][0] = 0.5 * v0;
    s.fb.vx[0][1] = -0.5 * v0;
    s.ffdef.push_back(FFDef{0, -1, 1, -1, 0.0, 0.0, kD, 0.3});
    for (int i = 0; i < 30; ++i) s.step();
    check(std::fabs(s.fb.vx[0][0]) < 1e-9 * v0 && std::fabs(s.fb.vx[0][1]) < 1e-9 * v0 && s.ff_delta(0) > 0.0,
          "удар из зазора: оба остановились, проникание = перелёт <= v0*dt");
}

// ============================================================================
// Т3. Связанные контакты. Остаток после 100 шагов, таблица по K и тёплому старту.
// ============================================================================
struct SceneInfo
{
    const char* name;
    void (*make)(Sim&, double);
    std::size_t nff, nfs; // число контактов
};

const SceneInfo kScenes[3] = {{"узел из 3 нитей", &make_knot, 3, 0},
                              {"Solid на 3 бусинах", &make_tripod, 0, 3},
                              {"цепочка Solid-сегмент-бусина", &make_chain, 1, 1}};

struct T3Run
{
    NormalProjectionStats st;
    double Wpos_rel = 0.0;
    double mom_rel = 0.0;
    double maxerr_N = 0.0; // max |N_c/N_exp - 1| в конце
    std::size_t nff = 0;
};

T3Run run_t3(const SceneInfo& sc, int threads, int K, bool warm, std::size_t steps, Sim* out = nullptr)
{
    const double N0 = kM * kG;
    Sim local(threads);
    Sim& s = out ? *out : local;
    s.K = K;
    s.warm_on = warm;
    sc.make(s, N0);
    for (std::size_t i = 0; i < steps; ++i) s.step();
    T3Run r;
    r.st = s.st;
    r.Wpos_rel = s.E_ref > 0.0 ? s.Wpos_max / s.E_ref : 0.0;
    r.mom_rel = s.mom_rel_max;
    // ожидаемые N: узел и Solid -- N0*(1,2,3); цепочка -- ff 2*N0, fs N0
    for (std::size_t k = 0; k < s.ff.size(); ++k)
    {
        const double expect = (sc.make == &make_chain) ? 2.0 * N0 : N0 * static_cast<double>(k + 1);
        r.maxerr_N = std::max(r.maxerr_N, std::fabs(s.ff[k].P / s.dt / expect - 1.0));
    }
    for (std::size_t k = 0; k < s.fs.size(); ++k)
    {
        const double expect = (sc.make == &make_chain) ? N0 : N0 * static_cast<double>(k + 1);
        r.maxerr_N = std::max(r.maxerr_N, std::fabs(s.fs[k].P / s.dt / expect - 1.0));
    }
    r.nff = s.ff.size() + s.fs.size();
    return r;
}

void test_t3()
{
    std::printf("[Т3] связанные контакты, 100 шагов: остаток resid_rel = max pdot^+/(W_c P_c), W_c = sum w^2/m без c_k "
                "(доли h N_c/m_red, ТЗ A4); resid_load = sum(pdot^+/W~) / sum P (доля полной нагрузки)\n");
    const int Ks[4] = {1, 2, 4, 8};
    for (const SceneInfo& sc : kScenes)
    {
        std::printf("  сцена: %s\n", sc.name);
        std::printf("    %-3s | %-16s | %-16s | %-16s | %-16s\n", "K", "resid_rel хол.", "resid_rel тёпл.",
                    "resid_load хол.", "resid_load тёпл.");
        for (int K : Ks)
        {
            const T3Run cold = run_t3(sc, 1, K, false, 100);
            const T3Run wrm = run_t3(sc, 1, K, true, 100);
            std::printf("    %-3d | %-16.3e | %-16.3e | %-16.3e | %-16.3e   (resid, м/с: %.2e | %.2e)\n", K,
                        cold.st.resid_rel, wrm.st.resid_rel, cold.st.resid_load, wrm.st.resid_load, cold.st.resid,
                        wrm.st.resid);
        }
        const T3Run r = run_t3(sc, 1, 8, true, 100);
        check(r.st.n_slots_ff == sc.nff && r.st.n_slots_fs == sc.nfs, "число активных слотов совпадает со сценой");
        check(r.st.resid_rel < 0.01, std::string("K = 8 + тёплый старт: остаток < 1% нагрузки (") + sc.name + ")");
        check(r.maxerr_N < 1e-2, std::string("N_c в пределах 1% от заданной нагрузки (") + sc.name + ")");
        check(r.st.warm_frac == 1.0, "доля слотов с тёплым стартом = 1");
        check(r.Wpos_rel <= 1e-12, std::string("Т4: W_stage <= 1e-12 E (") + sc.name + ")");
        check(r.mom_rel < 1e-12, std::string("Т5: |sum m a_stage| на уровне округления (") + sc.name + ")");
        std::printf("    stats: max c_k = %.0f, max c_s = %.0f, W_stage(последний) = %.3e Дж, max N_err = %.2e, "
                    "resid_rel = %.3e, resid_load = %.3e (K = 8, тёплый старт)\n",
                    r.st.max_ck, r.st.max_cs, r.st.W_stage, r.maxerr_N, r.st.resid_rel, r.st.resid_load);
    }
}

// ============================================================================
// Т5/Т7. Случайная сцена на уровне стадии: импульс, детерминизм на 1 и 8 потоках.
// ============================================================================
void make_random(Sim& s, std::uint64_t seed, bool with_fixed)
{
    Rng rng(seed);
    const std::size_t nf = 150, ns = 40, nff = 260, nfs = 90;
    s.init(nf, ns);
    s.check_momentum = !with_fixed;
    for (std::size_t i = 0; i < nf; ++i) s.fb.mass[i] = rng.range(0.5, 5.0) * kM;
    for (std::size_t i = 0; i < ns; ++i) s.sys.mass[i] = rng.range(1.0, 20.0) * kM;
    for (std::size_t d = 0; d < dim; ++d)
    {
        for (std::size_t i = 0; i < nf; ++i)
        {
            s.fb.vx[d][i] = rng.range(-10.0, 10.0);
            s.ext_f[d][i] = rng.range(-1.0, 1.0) * kG;
        }
        for (std::size_t i = 0; i < ns; ++i)
        {
            s.sys.vx[d][i] = rng.range(-10.0, 10.0);
            s.ext_s[d][i] = rng.range(-1.0, 1.0) * kG;
        }
    }
    if (with_fixed)
    {
        for (std::size_t i = 0; i < nf; ++i) s.fixed_f[i] = (rng.uni() < 0.1) ? 1 : 0;
        for (std::size_t i = 0; i < ns; ++i) s.fixed_s[i] = (rng.uni() < 0.1) ? 1 : 0;
        // Закреплённые -- НЕПОДВИЖНЫЕ опоры, как зажимы сцены (CLAMP_BOX): у движущейся опоры бесконечной массы
        // положительная работа над свободными частицами законна (опора толкает), и критерий Т4 к ней не
        // применим. Прежде здесь оставались случайные скорости -- отсюда «+3e-3 Дж за шаг» (ревью Codex).
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < nf; ++i)
                if (s.fixed_f[i]) s.fb.vx[d][i] = 0.0;
            for (std::size_t i = 0; i < ns; ++i)
                if (s.fixed_s[i]) s.sys.vx[d][i] = 0.0;
        }
    }
    s.finalize_fixed();
    auto weight = [&](double& w) {
        const double u = rng.uni();
        w = (u < 0.15) ? 0.0 : (u < 0.3 ? 1.0 : rng.uni());
    };
    s.ff.assign(nff, FiberNormalContact{});
    for (std::size_t k = 0; k < nff; ++k)
    {
        FiberNormalContact& c = s.ff[k];
        std::size_t idx[4];
        for (int j = 0; j < 4; ++j)
        {
            bool again;
            do
            {
                idx[j] = rng.index(nf);
                again = false;
                for (int q = 0; q < j; ++q) again = again || (idx[q] == idx[j]);
            } while (again);
        }
        c.active = rng.uni() < 0.95;
        c.a0 = idx[0];
        c.a1 = (rng.uni() < 0.2) ? -1 : static_cast<long long>(idx[1]);
        c.b0 = idx[2];
        c.b1 = (rng.uni() < 0.2) ? -1 : static_cast<long long>(idx[3]);
        weight(c.s);
        weight(c.t);
        // Точечный сегмент: геометрия даёт параметр 0 (fiber_segment_closest_params,
        // вырожденный отрезок), вес второй бусины не задан -- иначе сцена нефизична.
        if (c.a1 == -1) c.s = 0.0;
        if (c.b1 == -1) c.t = 0.0;
        rng.unit(c.n);
        c.mu = 0.3;
        c.key = make_ff_pair_key(static_cast<long long>(c.a0), c.a1, static_cast<long long>(c.b0), c.b1);
    }
    s.fs.assign(nfs, FiberSolidNormalContact{});
    for (std::size_t k = 0; k < nfs; ++k)
    {
        FiberSolidNormalContact& c = s.fs[k];
        const std::size_t a = rng.index(nf);
        std::size_t b = rng.index(nf);
        if (b == a) b = (a + 1) % nf;
        c.active = rng.uni() < 0.95;
        c.a = a;
        c.b = (rng.uni() < 0.2) ? -1 : static_cast<long long>(b);
        c.s = rng.index(ns);
        weight(c.lambda);
        if (c.b == -1) c.lambda = 0.0; // точечный контакт: lambda = 0 (fiber_solid_contact.hpp)
        rng.unit(c.n);
        c.mu = 0.3;
        c.key = make_fs_pair_key(static_cast<long long>(c.a), c.b, static_cast<long long>(c.s));
    }
}

// Три стадии подряд (геометрия слотов неизменна, скорости обновляются, тёплый старт включён).
void run_random_steps(Sim& s)
{
    s.K = 8;
    s.warm_on = true;
    for (int it = 0; it < 3; ++it)
    {
        // слоты: геометрия та же, P/warm выставляет warm.apply
        s.prepare_dv();
        s.stage();
        s.advance();
    }
}

void test_random_scene()
{
    std::printf("[Т5/Т7] случайная сцена: 150 бусин, 40 частиц Solid, 260 + 90 контактов (сегменты, нулевые веса)\n");
    {
        Sim s(1);
        make_random(s, 2026, false);
        run_random_steps(s);
        std::printf("    импульс: max |sum m Δv| / sum m|Δv| = %.3e; max c_k = %.0f, max c_s = %.0f; resid_rel = %.3e; "
                    "resid_load = %.3e\n",
                    s.mom_rel_max, s.st.max_ck, s.st.max_cs, s.st.resid_rel, s.st.resid_load);
        // Т4 на случайной (связной, c_k >= 3) сцене: положительная работа стадии <= 1e-12 E на каждом шаге.
        const double wpos_rel = s.E_ref > 0.0 ? s.Wpos_max / s.E_ref : 0.0;
        std::printf("T4 random: W_pos/E = %.3e (порог 1e-12)\n", wpos_rel);
        check(wpos_rel <= 1e-12, "Т4: W_stage <= 1e-12 E на каждом шаге (случайная сцена)");
        check(s.mom_rel_max < 1e-12, "Т5: суммарный импульс стадии = 0 с точностью округления (все свободны)");
        check(s.st.max_ck >= 3.0, "сцена действительно связная: max c_k >= 3");
    }
    for (int with_fixed = 0; with_fixed < 2; ++with_fixed)
    {
        Sim a(1), b(8);
        make_random(a, 31 + static_cast<std::uint64_t>(with_fixed), with_fixed != 0);
        make_random(b, 31 + static_cast<std::uint64_t>(with_fixed), with_fixed != 0);
        run_random_steps(a);
        run_random_steps(b);
        bool eq = true;
        for (std::size_t d = 0; d < dim; ++d)
        {
            eq = eq && bits_equal(a.dv_f[d], b.dv_f[d]) && bits_equal(a.dvx_s[d], b.dvx_s[d]) &&
                 bits_equal(a.fb.vx[d], b.fb.vx[d]) && bits_equal(a.sys.vx[d], b.sys.vx[d]) &&
                 bits_equal(a.ws.dvel_fiber[d], b.ws.dvel_fiber[d]) && bits_equal(a.ws.dvel_solid[d], b.ws.dvel_solid[d]);
        }
        bool pe = (a.ff.size() == b.ff.size()) && (a.fs.size() == b.fs.size());
        for (std::size_t k = 0; pe && k < a.ff.size(); ++k) pe = std::memcmp(&a.ff[k].P, &b.ff[k].P, sizeof(double)) == 0;
        for (std::size_t k = 0; pe && k < a.fs.size(); ++k) pe = std::memcmp(&a.fs[k].P, &b.fs[k].P, sizeof(double)) == 0;
        bool st = std::memcmp(&a.st.resid, &b.st.resid, sizeof(double)) == 0 &&
                  std::memcmp(&a.st.resid_rel, &b.st.resid_rel, sizeof(double)) == 0 &&
                  std::memcmp(&a.st.resid_load, &b.st.resid_load, sizeof(double)) == 0 &&
                  std::memcmp(&a.st.W_stage, &b.st.W_stage, sizeof(double)) == 0 && a.st.max_ck == b.st.max_ck;
        bool wt = a.warm.size() == b.warm.size();
        for (std::size_t k = 0; wt && k < a.warm.size(); ++k)
            wt = std::memcmp(&a.warm.entries()[k].P, &b.warm.entries()[k].P, sizeof(double)) == 0 &&
                 a.warm.entries()[k].key == b.warm.entries()[k].key;
        check(eq && pe && st && wt,
              std::string("Т7: стадия на 1 и 8 потоках побитово одинакова (со ") + (with_fixed ? "" : "без ") +
                  "закреплённых), W_stage = " + std::to_string(a.st.W_stage));
        check(std::memcmp(&a.st.warm_scale, &b.st.warm_scale, sizeof(double)) == 0 &&
                  std::memcmp(&a.st.W_warm, &b.st.W_warm, sizeof(double)) == 0,
              "Т7: множитель тёплого старта и его работа одинаковы на 1 и 8 потоках (побитово)");
        // Т4 на той же сцене (ревью: со закреплёнными бусинами прежде было +3e-3 Дж за шаг).
        const double wpos_rel = a.E_ref > 0.0 ? a.Wpos_max / a.E_ref : 0.0;
        std::printf("    Т4 random (%s закреплённых): W_pos/E = %.3e (порог 1e-12); последний шаг: W_warm = %.3e Дж, "
                    "warm_scale = %.17g, W_stage = %.3e Дж\n",
                    with_fixed ? "со" : "без", wpos_rel, a.st.W_warm, a.st.warm_scale, a.st.W_stage);
        check(wpos_rel <= 1e-12, std::string("Т4: W_stage <= 1e-12 E на каждом шаге (случайная сцена ") +
                                     (with_fixed ? "со" : "без") + " закреплённых)");
    }
}

// ============================================================================
// Т4а. Тёплый старт на расходящейся паре (ревью Codex): импульс прошлого шага у пары,
// которая уже разлетается, дал бы положительную работу. Энергетически безопасный тёплый
// старт уменьшает P^(0) так, что его работа <= 0 (здесь -- до нуля): стадия пару не разгоняет.
// ============================================================================
void test_warm_separating()
{
    std::printf("[Т4а] тёплый старт на расходящейся паре\n");
    Sim s(1);
    s.init(2, 0);
    const double delta0 = 0.05 * kD;
    s.fb.x[0][0] = -0.5 * (kD - delta0); // A слева, n -- от B к A, то есть -x
    s.fb.x[0][1] = 0.5 * (kD - delta0);
    s.fb.vx[0][0] = -1.0; // пара разлетается: pdot = -n.(vA - vB) = -2 м/с
    s.fb.vx[0][1] = 1.0;
    s.ffdef.push_back(FFDef{0, -1, 1, -1, 0.0, 0.0, kD, 0.3});
    // Таблица тёплого старта: большой импульс прошлого шага у этой пары, h_prev = h.
    std::vector<NormalWarmEntry> tab(1);
    tab[0].key = make_ff_pair_key(0, -1, 1, -1);
    tab[0].P = 10.0 * kM; // без отбрасывания добавил бы паре 10 м/с на бусину
    s.warm.set(tab, s.dt);
    const double vA0 = s.fb.vx[0][0], vB0 = s.fb.vx[0][1];
    s.step();
    std::printf("    расходящаяся: W_warm = %.3e Дж, warm_scale = %.17g, P = %.3e, W_stage = %.3e Дж\n", s.st.W_warm,
                s.st.warm_scale, s.ff[0].P, s.st.W_stage);
    check(s.st.n_slots_ff == 1 && s.st.W_warm > 0.0, "слот активен, работа тёплого старта положительна");
    check(s.st.warm_scale == 0.0, "тёплый старт обнулён (warm_scale = 0: он только разгонял пару)");
    check(s.ff[0].P == 0.0, "импульс шага P = 0");
    check(s.fb.vx[0][0] == vA0 && s.fb.vx[0][1] == vB0, "скорости пары не изменились (побитово)");
    check(s.st.W_stage == 0.0, "работа стадии 0");
    check(s.warm.size() == 0, "в таблицу тёплого старта пара не попала (P = 0)");

    // Та же пара, но сближается: тёплый старт сохраняется (отбрасывается только у расходящихся).
    Sim c(1);
    c.init(2, 0);
    c.fb.x[0][0] = -0.5 * (kD - delta0);
    c.fb.x[0][1] = 0.5 * (kD - delta0);
    c.fb.vx[0][0] = 1.0;
    c.fb.vx[0][1] = -1.0;
    c.ffdef.push_back(FFDef{0, -1, 1, -1, 0.0, 0.0, kD, 0.3});
    c.warm.set(tab, c.dt);
    c.step();
    std::printf("    сближающаяся: W_warm = %.3e Дж, warm_scale = %.17g, P/(m v) = %.6f\n", c.st.W_warm, c.st.warm_scale,
                c.ff[0].P / kM);
    // P^(0) = 10 m при нужном 1 m: работа тёплого старта > 0, lambda = -2A0/B0 = 0.2, затем итерации доводят до 1 m.
    check(c.ff[0].warm && c.st.warm_scale > 0.0 && c.st.warm_scale < 1.0,
          "сближающаяся пара с завышенным P^(0): тёплый старт уменьшен, но не отброшен");
    check(std::fabs(c.ff[0].P / kM - 1.0) < 1e-12, "сближающаяся пара: итоговый импульс = m_red * 2 м/с = 1 m");
    check(std::fabs(c.fb.vx[0][0] - c.fb.vx[0][1]) < 1e-12, "сближающаяся пара: относительная скорость погашена");
    check(c.Wpos_max <= 1e-12 * c.E_ref, "сближающаяся пара: W_stage <= 1e-12 E");
}

void test_t7_scenes()
{
    std::printf("[Т7] сцены Т3 на 1 и 8 потоках (100 шагов, K = 8, тёплый старт)\n");
    for (const SceneInfo& sc : kScenes)
    {
        Sim a(1), b(8);
        run_t3(sc, 1, 8, true, 100, &a);
        run_t3(sc, 8, 8, true, 100, &b);
        bool eq = true;
        for (std::size_t d = 0; d < dim; ++d)
        {
            eq = eq && bits_equal(a.fb.x[d], b.fb.x[d]) && bits_equal(a.fb.vx[d], b.fb.vx[d]) &&
                 bits_equal(a.sys.x[d], b.sys.x[d]) && bits_equal(a.sys.vx[d], b.sys.vx[d]) &&
                 bits_equal(a.dv_f[d], b.dv_f[d]) && bits_equal(a.dvx_s[d], b.dvx_s[d]);
        }
        for (std::size_t k = 0; eq && k < a.ff.size(); ++k) eq = std::memcmp(&a.ff[k].P, &b.ff[k].P, sizeof(double)) == 0;
        for (std::size_t k = 0; eq && k < a.fs.size(); ++k) eq = std::memcmp(&a.fs[k].P, &b.fs[k].P, sizeof(double)) == 0;
        check(eq, std::string("побитово одинаковы P, ускорения и траектории (") + sc.name + ")");
    }
}

// ============================================================================
// Т6. Сходимость по шагу: dt_scale = 1; 0,5; 0,25 на одном физическом времени.
// ============================================================================

// Удар из зазора на закреплённую бусину под нагрузкой (одна и та же сцена для проверки по концу и по рядам).
constexpr double kT6ImpactV = 79.0; // скорость удара, м/с (характерная скорость сцены удара)
void setup_t6_impact(Sim& s, double dt)
{
    s.dt = dt;
    s.init(2, 0);
    s.fixed_f[0] = 1;
    s.finalize_fixed();
    s.fb.x[2][1] = kD + 20.0 * kT6ImpactV * kDt0; // фиксированный зазор
    s.fb.vx[2][1] = -kT6ImpactV;
    s.ext_f[2][1] = -kG;
    s.ffdef.push_back(FFDef{1, -1, 0, -1, 0.0, 0.0, kD, 0.3});
}

// Ряд delta(t) и v(t) подвижной бусины 1 в npts одинаковых физических моментах t_j = t_first + j*t_step.
// Моменты кратны 0,25*kDt0 (наименьший шаг серии), поэтому для dt_scale 0,5 и 0,25 номер шага
// n_j = t_j/dt целый -- шаг точен, интерполяция не нужна (для dt_scale 1 в серии удара полуцелые t_j
// округляются до ближайшего шага; эта пара только печатается). Снимается состояние (x^n, v^n) после n_j шагов;
// скорость v^n сдвинута от t_j на dt/2 (полунеявная схема), разница между шагами <= g*dt/2 = 0,125 м/с.
struct T6Series
{
    std::vector<double> delta, v;
};

T6Series sample_t6_series(Sim& s, double t_first, double t_step, int npts)
{
    T6Series r;
    std::size_t done = 0;
    for (int j = 0; j < npts; ++j)
    {
        const double tj = t_first + static_cast<double>(j) * t_step;
        const std::size_t target = static_cast<std::size_t>(std::llround(tj / s.dt));
        while (done < target)
        {
            s.step();
            ++done;
        }
        r.delta.push_back(s.ff_delta(0));
        r.v.push_back(s.fb.vx[2][1]);
    }
    return r;
}

// max_j |a_j - b_j| / ref для delta и v ряда.
void series_diff(const T6Series& a, const T6Series& b, double delta_ref, double v_ref, double& ed, double& ev)
{
    ed = 0.0;
    ev = 0.0;
    for (std::size_t j = 0; j < a.delta.size() && j < b.delta.size(); ++j)
    {
        ed = std::max(ed, std::fabs(a.delta[j] - b.delta[j]) / delta_ref);
        ev = std::max(ev, std::fabs(a.v[j] - b.v[j]) / v_ref);
    }
}

void test_t6()
{
    std::printf("[Т6] сходимость по шагу: dt_scale 1; 0.5; 0.25\n");
    const double scales[3] = {1.0, 0.5, 0.25};
    const double T = 2000.0 * kDt0; // физическое время

    // (a) статика Т1 (variant 0): delta(T)
    double dl[3], d3[3];
    for (int i = 0; i < 3; ++i)
    {
        const double dt = kDt0 * scales[i];
        const std::size_t steps = static_cast<std::size_t>(std::llround(T / dt));
        dl[i] = run_t1(0, dt, steps, 8, true).delta_end;
        d3[i] = law3_static(dt, steps, 0.1 * kD).delta_end;
        std::printf("    scale %.2f: статика  delta/d: проекция = %.9f, закон 3 = %.9f\n", scales[i], dl[i] / kD, d3[i] / kD);
    }
    check(std::fabs(dl[1] - dl[2]) < 0.01 * 0.1 * kD, "статика: delta(T) между шагами 0.5 и 0.25 отличается < 1% от delta0");
    check(std::fabs(dl[0] - dl[1]) < 0.01 * 0.1 * kD, "статика: то же между 1 и 0.5");

    // (b) удар из зазора на закреплённую бусину под нагрузкой: delta_end и скорость
    double de[3], ve[3];
    for (int i = 0; i < 3; ++i)
    {
        const double dt = kDt0 * scales[i];
        Sim s(1);
        setup_t6_impact(s, dt);
        const std::size_t steps = static_cast<std::size_t>(std::llround(60.0 * kDt0 / dt));
        for (std::size_t k = 0; k < steps; ++k) s.step();
        de[i] = s.ff_delta(0);
        ve[i] = s.fb.vx[2][1];
        std::printf("    scale %.2f: удар     delta_end/d = %.6f, v_end = %.3e м/с, N/P = %.6f\n", scales[i], de[i] / kD, ve[i],
                    s.ff[0].P / dt / (kM * kG));
    }
    check(std::fabs(de[1] - de[2]) < 0.01 * kD && std::fabs(de[0] - de[1]) < 0.01 * kD,
          "удар: delta_end сходится (разница < 1% d)");
    check(std::fabs(ve[0]) < 1e-9 && std::fabs(ve[1]) < 1e-9 && std::fabs(ve[2]) < 1e-9,
          "удар: скорость после удара 0 при всех шагах (e = 0, нагрузка удержана)");

    // (c) РЯДЫ по времени: 20 одинаковых физических моментов, допуск 1 % (ТЗ A5) между dt_scale 0,5 и 0,25
    //     (между 1 и 0,5 -- только печать). Масштабы отнесения:
    //     - статика: delta(t) -- к delta0 = 0,1 d; скорость -- к g*dt0 = 1 м/с (скорость, которую бусина
    //       набрала бы за базовый шаг без контакта; сама скорость в статике ~0, так что делить на неё нельзя);
    //     - удар: delta(t) -- к d (как проверка по концу выше); скорость -- к скорости удара 79 м/с
    //       (мгновенная скорость после остановки ~0, делить на неё нельзя).
    //     Моменты удара сдвинуты на 1,5*kDt0: остановка бусины происходит при t_c ~ 18..19 kDt0, а сетка
    //     3*j*kDt0 попала бы прямо в t_c, где ряды разных шагов различались бы на скачок скорости.
    {
        const int npts = 20;
        T6Series st_ser[3], im_ser[3];
        for (int i = 0; i < 3; ++i)
        {
            const double dt = kDt0 * scales[i];
            {
                Sim s(1);
                s.dt = dt;
                s.K = 8;
                s.warm_on = true;
                s.check_momentum = false;
                make_static_scene(s, 0, 0.1 * kD);
                st_ser[i] = sample_t6_series(s, T / npts, T / npts, npts);
            }
            {
                Sim s(1);
                setup_t6_impact(s, dt);
                im_ser[i] = sample_t6_series(s, 1.5 * kDt0, 3.0 * kDt0, npts);
            }
        }
        double ed12, ev12, ed23, ev23;
        series_diff(st_ser[0], st_ser[1], 0.1 * kD, kG * kDt0, ed12, ev12);
        series_diff(st_ser[1], st_ser[2], 0.1 * kD, kG * kDt0, ed23, ev23);
        std::printf("    ряды, статика (20 точек): max |d delta|/delta0 = %.3e (1 vs 0.5), %.3e (0.5 vs 0.25); "
                    "max |d v|/(g dt0) = %.3e, %.3e\n", ed12, ed23, ev12, ev23);
        check(ed23 < 0.01 && ev23 < 0.01, "статика: ряды delta(t), v(t) на 20 моментах между шагами 0.5 и 0.25 отличаются < 1%");
        series_diff(im_ser[0], im_ser[1], kD, kT6ImpactV, ed12, ev12);
        series_diff(im_ser[1], im_ser[2], kD, kT6ImpactV, ed23, ev23);
        std::printf("    ряды, удар (20 точек): max |d delta|/d = %.3e (1 vs 0.5), %.3e (0.5 vs 0.25); "
                    "max |d v|/V0 = %.3e, %.3e\n", ed12, ed23, ev12, ev23);
        check(ed23 < 0.01 && ev23 < 0.01, "удар: ряды delta(t), v(t) на 20 моментах между шагами 0.5 и 0.25 отличаются < 1%");
    }
}

} // namespace

int main()
{
    test_impulse_formula();
    test_counts_and_traversal();
    test_single_ff();
    test_fs_av();
    test_warm_start();
    test_friction_fill();
    test_t2();
    test_t3();
    test_random_scene();
    test_warm_separating();
    test_t7_scenes();
    test_t6();
    test_t1();
    if (failures == 0)
        std::printf("fiber_contact_projection_test: ALL OK\n");
    else
        std::printf("fiber_contact_projection_test: %d FAILED\n", failures);
    return failures == 0 ? 0 : 1;
}
