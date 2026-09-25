// test_rhs_equivalence.cpp
// Стенд эквивалентности старого и нового путей RHS для Solid (fiber/PLAN.md,
// этап Ф4, ЧАСТЬ A: "стенд test_rhs_equivalence (SPH_DIM=2): старый и новый
// RHS на одном состоянии, все производные совпадают с относительной точностью
// ~1e-12").
//
// СХЕМА:
//   1) строится ОДНО физически осмысленное состояние Solid: два тела разных
//      материалов (B4C itype=10, JH-2 EOS; PMMA itype=9, УА-Тэта EOS) в
//      контакте у оси симметрии (r мало -> часть частиц получает зеркальные
//      частицы, virt_part.cpp/GhostSet), с ненулевыми радиальной и осевой
//      скоростями (упражняет плотность, кинематику, искусственную вязкость,
//      осесимметричные члены) и ненулевым начальным девиатором напряжений
//      (упражняет производную Яуманна и девиаторную работу в энергии);
//   2) СТАРЫЙ путь: virt_part -> direct_find -> con_density -> int_force ->
//      art_visc -> ext_force, В ТОЧНОСТИ как их вызывает src/single_step.cpp
//      (тот же порядок и сигнатуры, что в tests/test_axisym.cpp);
//   3) НОВЫЙ путь: interactions/solid_sph.hpp + interactions/solid_contact.hpp
//      (Solid derived -> GhostSet -> SolidExtended -> поиск -> kernel'ы SPH и
//      контакта -> rhs::reduce_contributions -> finalize_solid);
//   4) сравниваются ВСЕ производные (drho, du, dvxdt по всем dim, dsxx..dsxz)
//      частица за частицей по относительной норме.
//
// Прогон повторяется для нескольких комбинаций (pa_sph, density_method,
// contact_force_type), чтобы упражнять обе ветви pa_sph, все три
// density_method и все типы контактной силы -- не только тот единственный
// набор (pa_sph=2, density_method=1), что реально используют BC2/BC6
// (tests/data/bc2_initial_conditions.txt, tests/data/bc6_initial_conditions.txt).
#include "contact_force.hpp"
#include "density.hpp"
#include "direct_find.hpp"
#include "ext_force.hpp"
#include "int_force.hpp"
#include "art_visc.hpp"
#include "interactions/solid_contact.hpp"
#include "interactions/solid_sph.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"
#include "virt_part.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

static_assert(dim == 2, "test_rhs_equivalence требует сборки с SPH_DIM=2");

namespace
{

int g_failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++g_failures;
}

// --- Геометрия испытательного состояния -------------------------------------
// Два тела в контакте у оси: A (B4C, itype=10) движется вдоль +z навстречу
// покоящемуся B (PMMA, itype=9). Оба тела заходят вглубь друг друга на
// полшага решётки (сознательный небольшой overlap: контакт СРАБАТЫВАЕТ уже
// на этом, единственном, вычислении RHS, а не только на следующем шаге).
// Радиальная скорость -- сходящаяся в A (v_r = -k*r, упражняет искусственную
// вязкость и осесимметричные члены), слегка расходящаяся в B.
constexpr std::size_t kNr = 5;
constexpr std::size_t kNzA = 5;
constexpr std::size_t kNzB = 5;
constexpr double kDx = 1.0e-4;
constexpr int kItypeA = 10;   // B4C, JH-2
constexpr int kItypeB = 9;    // PMMA, УА-Тэта
constexpr double kRhoA = 2517.0;
constexpr double kRhoB = 1186.0;

std::size_t idxA(std::size_t ir, std::size_t iz) { return ir * kNzA + iz; }
std::size_t idxB(std::size_t ir, std::size_t iz) { return kNr * kNzA + ir * kNzB + iz; }

void build_state(ParticleSystem& sys, std::size_t& nreal)
{
    const std::size_t nA = kNr * kNzA;
    const std::size_t nB = kNr * kNzB;
    nreal = nA + nB;
    sys.resize(nreal * 2); // запас под зеркальные частицы (см. tests/test_axisym.cpp)

    const double ZB0 = (static_cast<double>(kNzA) - 0.5) * kDx; // сознательный overlap 0.5*dx

    for (std::size_t ir = 0; ir < kNr; ++ir)
    {
        const double r = (static_cast<double>(ir) + 0.5) * kDx;
        for (std::size_t iz = 0; iz < kNzA; ++iz)
        {
            const std::size_t k = idxA(ir, iz);
            const double z = static_cast<double>(iz) * kDx;
            sys.x[0][k] = r; sys.x[1][k] = z;
            sys.vx[0][k] = -50.0 * r;                 // сходящееся радиальное поле
            sys.vx[1][k] = 300.0 + 20.0 * r;           // движение к телу B + сдвиг
            sys.itype[k] = kItypeA; sys.body[k] = kItypeA;
            // rho слегка ВЫШЕ rho_0 (сжатие) -- НЕ равно rho_0: иначе mu =
            // rho/rho_0 - 1 == 0 и обе ветви JH-2/УА-Тэта (mu>=0 при D==0)
            // дают p == 0 тождественно, из-за чего энергетический член "he"
            // (int_force.cpp, часть от ДАВЛЕНИЯ) был бы нулевым в ЛЮБОМ
            // сценарии -- давление никогда бы не упражнялось, а pa_sph=1
            // против pa_sph=2 в "he" различались бы формулой, но не числом
            // (см. диагностику ниже, в main()).
            sys.rho[k] = kRhoA * 1.02; sys.rho_0[k] = kRhoA;
            sys.mass[k] = kRhoA * kDx * kDx;
            sys.hsml[k] = 1.3 * kDx;                   // ОДНОРОДНАЯ hsml -- см. solid_sph.hpp
            sys.u[k] = 1.0e3; sys.dp[k] = 0.0; sys.D[k] = 0.0; sys.e_pl[k] = 0.0;
            sys.p[k] = 0.0;
            sys.sxx[k] = 1.0e7; sys.syy[k] = -0.4e7; sys.szz[k] = -0.6e7;
            sys.sxy[k] = 0.35e7; sys.syz[k] = 0.0; sys.sxz[k] = 0.0;
            sys.axis[k] = (ir == 0) ? 1 : 0;
        }
    }
    for (std::size_t ir = 0; ir < kNr; ++ir)
    {
        const double r = (static_cast<double>(ir) + 0.5) * kDx;
        for (std::size_t iz = 0; iz < kNzB; ++iz)
        {
            const std::size_t k = idxB(ir, iz);
            const double z = ZB0 + static_cast<double>(iz) * kDx;
            sys.x[0][k] = r; sys.x[1][k] = z;
            sys.vx[0][k] = 30.0 * r;                   // слегка расходящееся радиальное поле
            sys.vx[1][k] = 10.0 * (r - 2.5 * kDx);      // небольшой сдвиг, тело в целом покоится
            sys.itype[k] = kItypeB; sys.body[k] = kItypeB;
            // rho слегка НИЖЕ rho_0 (разрежение) -- противоположный знак mu
            // относительно тела A, чтобы обе ветви (mu>=0 / mu<0) обеих EOS
            // (JH-2 и УА-Тэта) были упражнены хотя бы одной частицей теста.
            sys.rho[k] = kRhoB * 0.98; sys.rho_0[k] = kRhoB;
            sys.mass[k] = kRhoB * kDx * kDx;
            sys.hsml[k] = 1.3 * kDx;
            sys.u[k] = 5.0e2; sys.dp[k] = 0.0; sys.D[k] = 0.0; sys.e_pl[k] = 0.0;
            sys.p[k] = 0.0;
            sys.sxx[k] = -0.6e7; sys.syy[k] = 0.5e7; sys.szz[k] = 0.1e7;
            sys.sxy[k] = -0.2e7; sys.syz[k] = 0.0; sys.sxz[k] = 0.0;
            sys.axis[k] = (ir == 0) ? 1 : 0;
        }
    }
}

// --- Старый путь --------------------------------------------------------------
struct OldResult
{
    std::vector<double> drho, du;
    std::vector<double> dv[dim];
    std::vector<double> dsxx, dsyy, dszz, dsxy, dsyz, dsxz;
};

OldResult run_old_path(const SimulationConfig& cfg, ParticleSystem sys, std::size_t nreal, double dt_applied)
{
    std::size_t nvirt = 0;
    virt_part(sys, nreal, nvirt, sys.size());
    const std::size_t nall = nreal + nvirt;

    std::size_t niac = 0;
    std::vector<std::size_t> pair_i(cfg.max_interaction), pair_j(cfg.max_interaction);
    std::vector<double> w(cfg.max_interaction, 0.0);
    std::vector<std::vector<double>> dwdx(dim, std::vector<double>(cfg.max_interaction, 0.0));
    std::vector<int> countiac(nall, 0);
    direct_find(cfg, sys, nall, niac, pair_i, pair_j, w, dwdx, countiac);

    PairGraph graph;
    graph.num_colors = 0; // atomic-ветвь; с omp_set_num_threads(1) даёт детерминированную сумму

    std::vector<double> drho(nall, 0.0);
    con_density(cfg, sys, nall, niac, pair_i, pair_j, dwdx, drho, graph);

    std::vector<double> eta(nall, 0.0), t(nall, 0.0);
    std::vector<double> tdsdt(nall, 0.0), dedt(nall, 0.0);
    std::vector<double> dsxx(nall, 0.0), dsyy(nall, 0.0), dszz(nall, 0.0);
    std::vector<double> dsxy(nall, 0.0), dsyz(nall, 0.0), dsxz(nall, 0.0);
    std::vector<std::vector<double>> indvxdt(dim, std::vector<double>(nall, 0.0));
    WorkBuffers wb;
    int_force(cfg, sys, nall, niac, eta, pair_i, pair_j, dwdx, w, t, indvxdt,
              tdsdt, dedt, dsxx, dsyy, dszz, dsxy, dsyz, dsxz, drho, wb, graph);

    std::vector<std::vector<double>> ardvxdt(dim, std::vector<double>(nall, 0.0));
    std::vector<double> avdudt(nall, 0.0);
    art_visc(cfg, sys, nall, niac, pair_i, pair_j, dwdx, ardvxdt, avdudt, graph);

    std::vector<std::vector<double>> exdvxdt(dim, std::vector<double>(nall, 0.0));
    ext_force(cfg, sys, dt_applied, nall, nreal, niac, pair_i, pair_j, exdvxdt, wb);

    OldResult out;
    out.drho.assign(drho.begin(), drho.begin() + static_cast<std::ptrdiff_t>(nreal));
    out.du.resize(nreal);
    out.dsxx.assign(dsxx.begin(), dsxx.begin() + static_cast<std::ptrdiff_t>(nreal));
    out.dsyy.assign(dsyy.begin(), dsyy.begin() + static_cast<std::ptrdiff_t>(nreal));
    out.dszz.assign(dszz.begin(), dszz.begin() + static_cast<std::ptrdiff_t>(nreal));
    out.dsxy.assign(dsxy.begin(), dsxy.begin() + static_cast<std::ptrdiff_t>(nreal));
    out.dsyz.assign(dsyz.begin(), dsyz.begin() + static_cast<std::ptrdiff_t>(nreal));
    out.dsxz.assign(dsxz.begin(), dsxz.begin() + static_cast<std::ptrdiff_t>(nreal));
    for (std::size_t d = 0; d < dim; ++d) out.dv[d].resize(nreal);
    for (std::size_t i = 0; i < nreal; ++i)
    {
        out.du[i] = dedt[i] + avdudt[i]; // ahdudt == 0 (heat_artificial == false)
        for (std::size_t d = 0; d < dim; ++d)
            out.dv[d][i] = indvxdt[d][i] + exdvxdt[d][i] + ardvxdt[d][i];
    }
    return out;
}

// --- Новый путь ----------------------------------------------------------------
std::vector<interactions::SolidFinal> run_new_path(const SimulationConfig& cfg, const ParticleSystem& sys,
                                                     std::size_t nreal, double dt_applied)
{
    rhs::CpuExecutor exec(1); // один поток -- побитная воспроизводимость сравнения

    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec, cfg.mats, sys, nreal, derived);

    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, nreal, ghosts);

    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, nreal, derived, ghosts, ext);

    interactions::SolidSphSearchResult search;
    interactions::build_solid_sph_search(exec, cfg, ext, search);

    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sph_buffers, contact_buffers;
    interactions::compute_solid_sph_contributions(exec, cfg, ext, search.primitives, sph_buffers);
    interactions::compute_solid_contact_contributions(exec, cfg, cfg.mats, sys, derived.c, dt_applied,
                                                        search.real_pairs, contact_buffers);

    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> all_buffers = sph_buffers;
    all_buffers.insert(all_buffers.end(), contact_buffers.begin(), contact_buffers.end());

    const std::vector<interactions::SolidDelta> reduced = rhs::reduce_contributions(exec, all_buffers, nreal);

    std::vector<interactions::SolidFinal> out;
    interactions::finalize_solid(exec, cfg.mats, sys, nreal, reduced, out);
    return out;
}

// --- Сравнение -------------------------------------------------------------
struct WorstDiff
{
    double rel = 0.0;
    std::size_t particle = 0;
    std::string field;
    double old_val = 0.0, new_val = 0.0;
};

double rel_diff(double a, double b, double abs_floor)
{
    const double scale = std::max({std::fabs(a), std::fabs(b), abs_floor});
    return std::fabs(a - b) / scale;
}

void track(WorstDiff& worst, double a, double b, std::size_t i, const char* field, double abs_floor)
{
    const double r = rel_diff(a, b, abs_floor);
    if (r > worst.rel)
    {
        worst.rel = r;
        worst.particle = i;
        worst.field = field;
        worst.old_val = a;
        worst.new_val = b;
    }
}

// abs_floor -- масштаб, ниже которого величины считаются "практически нулевыми"
// и сравниваются в абсолютной, а не относительной норме (иначе 1e-20 против
// -1e-20 давали бы relative diff = 2 при физически нулевых значениях).
bool compare(const OldResult& oldr, const std::vector<interactions::SolidFinal>& newr,
             std::size_t nreal, const char* label, double tol)
{
    WorstDiff worst;
    const double floor_rho = 1.0;      // кг/(м^3 с) -- на фоне |drho| ~ 1e6..1e9
    const double floor_du  = 1.0;      // Дж/(кг с)  -- на фоне |du|   ~ 1e3..1e9
    const double floor_dv  = 1.0e-3;   // м/с^2      -- на фоне |dv|  ~ 1e6..1e12
    const double floor_ds  = 1.0;      // Па/с       -- на фоне |ds|  ~ 1e13..1e18

    for (std::size_t i = 0; i < nreal; ++i)
    {
        track(worst, oldr.drho[i], newr[i].drho, i, "drho", floor_rho);
        track(worst, oldr.du[i],   newr[i].du,   i, "du",   floor_du);
        for (std::size_t d = 0; d < dim; ++d)
        {
            const std::string fld = "dv[" + std::to_string(d) + "]";
            track(worst, oldr.dv[d][i], newr[i].dv[d], i, fld.c_str(), floor_dv);
        }
        track(worst, oldr.dsxx[i], newr[i].dsxx, i, "dsxx", floor_ds);
        track(worst, oldr.dsyy[i], newr[i].dsyy, i, "dsyy", floor_ds);
        track(worst, oldr.dszz[i], newr[i].dszz, i, "dszz", floor_ds);
        track(worst, oldr.dsxy[i], newr[i].dsxy, i, "dsxy", floor_ds);
        track(worst, oldr.dsyz[i], newr[i].dsyz, i, "dsyz", floor_ds);
        track(worst, oldr.dsxz[i], newr[i].dsxz, i, "dsxz", floor_ds);
    }

    const bool ok = worst.rel <= tol;
    std::printf("  %-5s: %s -- худшее относительное расхождение %.3e (допуск %.1e)\n",
                ok ? "ok" : "FAIL", label, worst.rel, tol);
    std::printf("        частица #%zu, поле %s: старое=%.10e новое=%.10e\n",
                worst.particle, worst.field.c_str(), worst.old_val, worst.new_val);
    if (!ok) ++g_failures;
    return ok;
}

// --- Один сценарий -----------------------------------------------------------
// Возвращает результат НОВОГО пути -- нужен вызывающей стороне (main) для
// межсценарных проверок: убедиться, что варьируемые параметры (pa_sph,
// density_method, contact_force_type) реально дают РАЗНЫЕ числа между собой
// в НОВОМ коде, а не только совпадают со старым путём при каждом фиксированном
// наборе параметров (см. "ГЛАВНЫЙ ВОПРОС 1" в отчёте ревью части A: во всех
// 7 сценариях худшее старое/новое расхождение приходится на одну и ту же
// частицу №25, поле dsxy -- т.к. dsxy там определяется исключительно
// деформацией/поворотом и state.sxx/syy/sxz/syz, которые НЕ зависят от
// density_method/pa_sph/contact_force_type; чтобы убедиться, что ветки
// реально задействованы, а не просто "прогоняются вхолостую", сравниваются
// поля, которые ПО ФОРМУЛЕ ДОЛЖНЫ отличаться -- drho между density_method,
// du/dv между pa_sph, dv между contact_force_type).
std::vector<interactions::SolidFinal> run_scenario(int pa_sph, int density_method, int contact_type, const char* label)
{
    std::printf("\n[сценарий] %s (pa_sph=%d, density_method=%d, contact_force_type=%d)\n",
                label, pa_sph, density_method, contact_type);

    SimulationConfig cfg;
    cfg.apply_defaults();
    cfg.pa_sph = pa_sph;
    cfg.density_method = density_method;
    cfg.skf = 1;
    cfg.contact_force_type = contact_type;
    cfg.contact_penalty_scale = 2.5;
    cfg.contact_impulse_gamma = 0.5;
    cfg.CFL = 0.1;
    cfg.alpha = 1.0; cfg.beta = 1.0; cfg.etq = 0.01;
    cfg.xi = 0.001; cfg.epsilon_h = 0.01;

    ParticleSystem sys;
    std::size_t nreal = 0;
    build_state(sys, nreal);

    cfg.maxn = sys.size();
    cfg.max_interaction = nreal * 400;

    // density_method == 3 (density.cpp::con_density) читает sys.c ДО того,
    // как urs() внутри int_force его обновит на ЭТОМ шаге -- то есть видит
    // значение звуковой скорости С ПРЕДЫДУЩЕГО шага (в реальном многошаговом
    // расчёте оно уже есть; в тесте с "нулевого" состояния sys.c был бы 0).
    // Новый путь считает derived.c ДО SPH-стадии (fiber/PLAN.md §3.6, "Solid
    // derived"), то есть видит СВЕЖЕЕ значение. Чтобы сравнение было честным
    // (а не артефактом свежепостроенного sys.c == 0), предзаполняем sys.c тем
    // же значением, что и derived.c -- как это фактически было бы после
    // первого шага расчёта.
    for (std::size_t i = 0; i < nreal; ++i)
    {
        interactions::solid_eos(cfg.mats, sys, i, sys.p[i], sys.c[i]);
    }

    const double dt_applied = 1.0e-9;

    const OldResult oldr = run_old_path(cfg, sys, nreal, dt_applied);
    std::vector<interactions::SolidFinal> newr = run_new_path(cfg, sys, nreal, dt_applied);

    compare(oldr, newr, nreal, label, 1.0e-9);
    return newr;
}

// --- Диагностика "ветки реально разные" (ГЛАВНЫЙ ВОПРОС 1) ------------------
// max относительная разница ОДНОГО поля между двумя прогонами НОВОГО пути по
// всем реальным частицам. Используется, чтобы явно показать: варьируемые
// параметры меняют результат заметно (на уровне физики, >> порога roundoff
// 1e-9 из compare()), а не совпадают тождественно, как это происходит у
// dsxy/rxy-подобных полей, не зависящих от этих параметров.
template <class Field>
double max_rel_diff_field(const std::vector<interactions::SolidFinal>& a,
                           const std::vector<interactions::SolidFinal>& b,
                           Field field, double abs_floor)
{
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i)
        worst = std::max(worst, rel_diff(field(a[i]), field(b[i]), abs_floor));
    return worst;
}

} // namespace

int main()
{
#ifdef _OPENMP
    omp_set_num_threads(1); // детерминированное суммирование старого пути (atomic-ветвь)
#endif
    std::printf("=== Стенд эквивалентности RHS Solid: старый путь vs новый (fiber/PLAN.md, Ф4-A) ===\n");

    // Основной сценарий -- РЕАЛЬНАЯ конфигурация BC2/BC6
    // (tests/data/bc2_initial_conditions.txt, tests/data/bc6_initial_conditions.txt):
    // pa_sph=2, density_method=1 (умолчание SimulationConfig).
    const auto base = run_scenario(2, 1, 0, "BC2/BC6: без контакта");
    const auto contact1 = run_scenario(2, 1, 1, "BC2/BC6: контакт тип 1 (pinball)");
    run_scenario(2, 1, 2, "BC2/BC6: контакт тип 2 (Swegle)");
    run_scenario(2, 1, 3, "BC2/BC6: контакт тип 3 (импульсный)");

    // Дополнительные комбинации -- для полноты переноса (pa_sph=1,
    // density_method 2 и 3), не используемые реальными кейсами, но
    // предусмотренные SPH.md/PLAN.md.
    const auto pasph1 = run_scenario(1, 1, 0, "pa_sph=1, без контакта");
    const auto density2 = run_scenario(2, 2, 0, "density_method=2, без контакта");
    const auto density3 = run_scenario(2, 3, 0, "density_method=3 (диффузия), без контакта");

    // --- ГЛАВНЫЙ ВОПРОС 1 (отчёт ревью части A) --------------------------
    // Во всех 7 сценариях выше compare() сообщает ОДНО и то же худшее
    // расхождение старое/новое (4.122e-12, частица #25, поле dsxy) -- это НЕ
    // значит, что density_method/pa_sph/contact_force_type не влияют на
    // расчёт: dsxy = 2*G*exy + rxy*(syy-sxx) + sxz*ryz + syz*rxz зависит
    // ТОЛЬКО от деформации/поворота (exy, rxy, rxz, ryz) и state.s** -- ни
    // одно из варьируемых полей там не участвует, поэтому и roundoff-разница
    // старое/новое там одна и та же во всех сценариях. Ниже -- прямая
    // проверка, что сами варьируемые параметры дают в НОВОМ пути РАЗНЫЕ
    // (не совпадающие в пределах roundoff) числа, то есть ветки
    // density_method 1/2/3, pa_sph 1/2 и contact_force_type реально
    // задействованы, а не молча пропускаются/схлопываются в одну и ту же
    // формулу.
    std::printf("\n[диагностика] ветки density_method/pa_sph/contact_force_type действительно расходятся между собой:\n");
    {
        // drho -- единственное поле, формула которого меняется по
        // density_method (SPH.md, "Неразрывность"): 1 vs 2 vs 3 должны
        // отличаться на уровне физики, а не roundoff.
        auto get_drho = [](const interactions::SolidFinal& f) { return f.drho; };
        const double d12 = max_rel_diff_field(base, density2, get_drho, 1.0);
        const double d13 = max_rel_diff_field(base, density3, get_drho, 1.0);
        const double d23 = max_rel_diff_field(density2, density3, get_drho, 1.0);
        std::printf("  drho: density_method 1 vs 2 -- rel %.3e; 1 vs 3 -- rel %.3e; 2 vs 3 -- rel %.3e\n", d12, d13, d23);
        check(d12 > 1.0e-3, "density_method=1 и 2 дают РАЗНЫЙ drho (ветка density_method реально задействована)");
        check(d13 > 1.0e-3, "density_method=1 и 3 дают РАЗНЫЙ drho (диффузионный член реально задействован)");
        check(d23 > 1.0e-3, "density_method=2 и 3 дают РАЗНЫЙ drho");

        // du и dv -- формулы pa_sph=1 и pa_sph=2 в SPH.md явно разные
        // (симметризация давления/напряжений: (p_i+p_j)/(rho_i*rho_j) против
        // p_i/rho_i^2+p_j/rho_j^2). dv сравнивается по ОБЕИМ компонентам
        // (max по d): контактная сила и сходящееся радиальное поле body A
        // направлены преимущественно по-разному (dv[1] -- ось z, между
        // телами; dv[0] -- радиус r), фиксированная компонента могла бы
        // случайно оказаться нечувствительной к варьируемому параметру.
        auto get_du = [](const interactions::SolidFinal& f) { return f.du; };
        const double pdu = max_rel_diff_field(base, pasph1, get_du, 1.0);
        double pdv = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto get_dvd = [d](const interactions::SolidFinal& f) { return f.dv[d]; };
            pdv = std::max(pdv, max_rel_diff_field(base, pasph1, get_dvd, 1.0e-3));
        }
        std::printf("  du: pa_sph=2 vs 1 -- rel %.3e; dv: pa_sph=2 vs 1 -- rel %.3e\n", pdu, pdv);
        check(pdu > 1.0e-3, "pa_sph=1 и 2 дают РАЗНЫЙ du (ветка pa_sph реально задействована)");
        check(pdv > 1.0e-3, "pa_sph=1 и 2 дают РАЗНЫЙ dv (ветка pa_sph реально задействована)");

        // dv -- контакт добавляет ускорение ТОЛЬКО частицам проникающей пары
        // (Solid-Solid-Contact); "без контакта" и "контакт тип 1" должны
        // разойтись хотя бы в норме по dv, раз состояние построено с
        // сознательным overlap (см. build_state). Геометрия контакта в этом
        // тесте -- преимущественно ВДОЛЬ оси z (dv[1]), поэтому сравниваются
        // обе компоненты, а не только dv[0]. Порог ниже, чем у
        // density_method/pa_sph (1e-3): контактная сила -- малая добавка
        // поверх уже большого внутреннего SPH-ускорения (девиатор ~1e7 Па),
        // поэтому её ОТНОСИТЕЛЬНЫЙ вклад в dv скромный; порог 1e-6 всё равно
        // на ~5 порядков выше чистого roundoff старое/новое (~1e-11 из
        // compare() выше) -- отличить "контакт реально приложен" от
        // "контакт молча не сработал" он позволяет уверенно.
        double cdv = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto get_dvd = [d](const interactions::SolidFinal& f) { return f.dv[d]; };
            cdv = std::max(cdv, max_rel_diff_field(base, contact1, get_dvd, 1.0e-3));
        }
        std::printf("  dv: без контакта vs контакт тип 1 -- rel %.3e\n", cdv);
        check(cdv > 1.0e-6, "contact_force_type=0 и 1 дают РАЗНЫЙ dv (контактная сила реально применяется)");
    }

    std::printf("\n=== %s (провалов: %d) ===\n", g_failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", g_failures);
    return g_failures == 0 ? 0 : 1;
}
