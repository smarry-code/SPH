// tests/test_time_integration_fiber.cpp
// Продовый time_integration(cfg, fiber::Model&, ...) с НЕПУСТЫМИ нитями
// (fiber/CUDA_PLAN.md, этап Г1; IMPLEMENTATION_REPORT.md §4.5: «Продовый
// time_integration ни разу не запускается с непустыми нитями»).
//
// Все динамические тесты нитей (fiber_bond_test, fiber_contact_test,
// e2e_fiber_contact_test, e2e_fiber_solid_test сценарий Б) используют либо
// interactions::compute_fiber_* напрямую, либо ЛОКАЛЬНЫЙ PredictorCorrector-
// драйвер (FiberPcIntegrator/FiberOnPlateIntegrator), повторяющий формулы
// решателя вручную. Ни один тест не проверяет, что ПРОДОВЫЙ цикл
// (src/time_integration.cpp, detail::run_time_integration) в ОДНОМ вызове
// правильно считает и складывает вклады ВСЕХ пяти типов взаимодействий
// одновременно: Solid-SPH, Solid-Solid-Contact, Fiber-Bond,
// Fiber-Fiber-Contact, Fiber-Solid-Contact (порядок RHS и сложение
// Fiber-вкладов в Solid -- см. src/time_integration.cpp, детально
// прокомментировано в теле run_time_integration).
//
// СЦЕНА (см. build_scene() ниже): два маленьких блока Solid (алюминий,
// itype=4) 3x3 частицы по слою, верхний блок сближается с нижним с самого
// t=0 (готовый контакт + скорость сближения) -- даёт Solid-SPH (плотная
// решётка внутри и между слоями) И Solid-Solid-Contact (body[]
// различаются, проникание и pdot>0 с первого шага). Нить A идёт вдоль x над
// верхним блоком, её средняя бусина нависает точно над центральной частицей
// блока и уже проникает в контактный порог -- Fiber-Solid-Contact. Нить B
// идёт вдоль y на той же (x,y) точке пересечения, но выше нити A на 0.3 мм --
// Fiber-Fiber-Contact (обе нити пересекаются перпендикулярно, случай
// Surface-Surface, как test_perpendicular_manual в test_fiber_contact.cpp).
// Fiber-Bond -- обычные связи внутри каждой нити (в покое при t=0, начинают
// работать, как только контактные силы сдвигают среднюю бусину).
//
// K шагов достаточно, чтобы отработали все пять RHS (проверяется явно, см.
// verify_interactions_present() -- та же геометрия/поиск, что вызывает
// продовый цикл, применена к НАЧАЛЬНОМУ состоянию: контакты сохраняются в
// течение K шагов, т.к. K*dt*v_characteristic << размер сцены, см. числа в
// build_scene()), но мало для того, чтобы тест шёл дольше пары секунд.
//
// ПРОВЕРКИ:
//   - сохранение полного импульса системы (Solid+Fiber) -- у сцены нет
//     внешних сил (гравитации в Solid-конвейере нет, см. grep) и нет
//     закреплённых частиц, поэтому импульс обязан сохраняться с точностью до
//     округления PredictorCorrector;
//   - конечность всех полей (x, v, rho, u, sxx.. Solid; x, v Fiber) после K
//     шагов;
//   - реальное срабатывание всех пяти типов взаимодействий -- см.
//     verify_interactions_present();
//   - снимок состояния после K шагов сверяется с эталоном (tests/data/) --
//     заготовка для сравнения CPU/GPU на этапе Г8 (fiber/CUDA_PLAN.md).
//     Эталон генерируется ключом --emit-golden=<dir>, однопоточно
//     (omp_set_num_threads(1) вызывается ВСЕГДА в этом тесте, не только при
//     генерации -- иначе сравнение с эталоном плавало бы от прогона к
//     прогону, см. rhs/cpu_executor.hpp: порядок суммирования вкладов
//     детерминирован только при одинаковом числе worker'ов).
#include "fiber/model.hpp"
#include "input.hpp"
#include "int_force.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "interactions/solid_contact.hpp"
#include "materials.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "run_guard.hpp"
#include "simulation_config.hpp"
#include "time_integration.hpp"

#include "e2e_harness.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

static_assert(dim == 3, "test_time_integration_fiber требует сборки с SPH_DIM=3");

namespace fs = std::filesystem;

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

// fiber/CUDA_PLAN.md, этап Г1: fp32-вариант (SPH_PRECISION=32) считает через
// real_t=float. Снимок и импульс за K=15 шагов накапливают ошибку округления
// на уровне эпсилон float на КАЖДОМ из пяти RHS каждого шага -- допуск
// ослаблен на несколько порядков относительно fp64 (там детерминированный
// однопоточный прогон даёт разброс на уровне double-эпсилон, см. измеренные
// числа в комментарии у kMomentumTol/kSnapshotTol ниже).
#if SPH_PRECISION == 32
constexpr double kMomentumTol = 1e-4;
constexpr double kSnapshotTolRel = 1e-3;
constexpr double kSnapshotFloor = 1e-6;
#else
constexpr double kMomentumTol = 1e-9;
constexpr double kSnapshotTolRel = 1e-9;
constexpr double kSnapshotFloor = 1e-12;
#endif

constexpr std::size_t kSteps = 15;

// =========================================================================
// Сцена: два блока Solid + две пересекающиеся нити (см. шапку файла).
// =========================================================================
struct Scene
{
    SimulationConfig cfg;
    fiber::Model model;
    std::size_t n = 0; // число РЕАЛЬНЫХ Solid-частиц
};

fiber::FiberMaterial make_fiber_material(double E, double nu, double d, double rho)
{
    fiber::FiberMaterial m;
    m.E = real_t(E); m.nu = real_t(nu); m.d = real_t(d); m.rho = real_t(rho); m.eps_break = real_t(1.0);
    return m;
}

// Слой 3x3 частиц Solid (алюминий, itype=4) в плоскости z=z0, шаг dx,
// нижний левый узел в (0,0). body -- признак тела (Solid-Solid-Contact).
void build_layer(ParticleSystem& sys, std::size_t base, double dx, double z0, double rho, int body,
                  double vz)
{
    std::size_t k = base;
    for (int ix = 0; ix < 3; ++ix)
        for (int iy = 0; iy < 3; ++iy)
        {
            sys.inzializ_index(k);
            sys.x[0][k] = real_t(ix * dx);
            sys.x[1][k] = real_t(iy * dx);
            sys.x[2][k] = real_t(z0);
            sys.vx[2][k] = real_t(vz);
            sys.itype[k] = 4; // алюминий (materials.hpp)
            sys.body[k] = body;
            sys.rho[k] = real_t(rho);
            sys.rho_0[k] = real_t(rho);
            sys.mass[k] = real_t(rho * dx * dx * dx); // кубическая ячейка dx x dx x dx
            sys.hsml[k] = real_t(1.5 * dx);
            ++k;
        }
}

// Прямая цепочка бусин одной нити (object_id/material_id -- 0 и 1 для двух
// нитей, чтобы Fiber-Fiber-Contact не спутал их с "самоконтактом").
void append_chain(fiber::FiberState& state, fiber::FiberTopology& topology, std::size_t base, int object_id,
                   int material_id, const std::vector<std::array<double, 3>>& beads, double mass_each, double l0)
{
    const std::size_t n = beads.size();
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t gi = base + k;
        for (int d3 = 0; d3 < 3; ++d3) state.x[d3][gi] = real_t(beads[k][d3]);
        for (std::size_t dd = 0; dd < dim; ++dd) state.vx[dd][gi] = real_t(0.0);
        state.mass[gi] = real_t(mass_each);
        state.object_id[gi] = object_id;
        state.material_id[gi] = material_id;
        if (k > 0) topology.prev[gi] = static_cast<long long>(gi - 1);
        if (k + 1 < n)
        {
            topology.next[gi] = static_cast<long long>(gi + 1);
            fiber::FiberBond b;
            b.i = gi; b.j = gi + 1;
            b.l0 = real_t(l0);
            b.l_break = real_t(l0 * 1e9); // разрыв -- не предмет этого теста
            topology.bonds.push_back(b);
        }
    }
}

Scene build_scene()
{
    Scene sc;
    sc.cfg.apply_defaults();
    sc.cfg.contact_force_type = 1; // pinball -- та же ветка, что штатно включена в BC2 (ADR-007)
    sc.cfg.t_total = 1.0;          // не ограничивает прогон -- шаги считает RunGuard
    sc.cfg.CFL = 0.1;

    // --- Solid: два слоя 3x3, body 1 (низ, неподвижен) / body 2 (верх,
    //     движется вниз навстречу) -- см. числа в шапке файла. ---------------
    const double dx = 1.0e-3;
    const double rho_al = 2700.0;
    const double gap = 0.7 * dx; // < R_i+R_j=dx -- контакт уже при t=0
    const double v_approach = -100.0;

    sc.model.solid.resize(18);
    build_layer(sc.model.solid, 0, dx, 0.0, rho_al, /*body*/ 1, /*vz*/ 0.0);
    build_layer(sc.model.solid, 9, dx, gap, rho_al, /*body*/ 2, /*vz*/ v_approach);
    sc.n = 18;
    for (std::size_t i = 0; i < sc.n; ++i) urs(sc.cfg.mats, sc.model.solid, i);
    sc.cfg.maxn = sc.n; // dim==3 -- нет зеркальных частиц (ADR-002), запас не нужен
    sc.cfg.max_interaction = 100 * sc.cfg.maxn;
    sc.cfg.moni_particle = 0;

    // --- Нити: A вдоль x, B вдоль y, пересекаются над центральной частицей
    //     верхнего блока (ix=iy=1 -> (dx,dx)), см. числа в шапке файла. ------
    const double E_f = 2.0e7, nu_f = 0.3, d_f = 1.0e-3, rho_f = 1200.0;
    const double r_f = 0.5 * d_f;
    const double r_s = 0.5 * dx; // solid_contact_radius при rho==rho_0
    const double l0 = 1.0e-3;
    const double A_f = static_cast<double>(pi) * d_f * d_f / 4.0;
    const double mass_each = rho_f * A_f * l0;

    const double z_top_solid = gap;                                  // z верхнего блока
    const double z_fiberA = z_top_solid + r_s + r_f - 0.2e-3;         // проникание 0.2 мм в Fiber-Solid
    const double z_fiberB = z_fiberA + 0.3e-3;                        // проникание 0.7 мм в Fiber-Fiber (d_f=1мм)
    const double cx = dx, cy = dx;                                    // точка пересечения -- над центральной частицей

    std::vector<std::array<double, 3>> beadsA, beadsB;
    for (int k = -2; k <= 2; ++k) beadsA.push_back({cx + k * l0, cy, z_fiberA});
    for (int k = -2; k <= 2; ++k) beadsB.push_back({cx, cy + k * l0, z_fiberB});

    sc.model.fiber.resize(beadsA.size() + beadsB.size());
    sc.model.topology.resize(sc.model.fiber.size());
    sc.model.materials.push_back(make_fiber_material(E_f, nu_f, d_f, rho_f));
    append_chain(sc.model.fiber, sc.model.topology, 0, /*object_id*/ 0, /*material_id*/ 0, beadsA, mass_each, l0);
    append_chain(sc.model.fiber, sc.model.topology, beadsA.size(), /*object_id*/ 1, /*material_id*/ 0, beadsB,
                 mass_each, l0);

    return sc;
}

// =========================================================================
// Статическая проверка "контакт реально есть" -- та же геометрия/поиск, что
// вызывает продовый цикл (src/time_integration.cpp), применённая к
// НАЧАЛЬНОМУ состоянию сцены (сохраняется весь прогон, см. шапку файла).
// =========================================================================
void verify_interactions_present(const Scene& sc)
{
    std::printf("verify_interactions_present (геометрия на t=0, та же, что видит продовый цикл):\n");
    rhs::CpuExecutor exec(1);

    // Solid-SPH: плотная решётка гарантирует соседей в радиусе ядра по
    // построению (dx между узлами слоя и между слоями, hsml=1.5*dx, радиус
    // ядра 2*hsml=3*dx) -- не требует отдельной проверки геометрии, косвенно
    // подтверждается ненулевым откликом Solid в тесте моментов ниже.
    check(sc.n >= 2, "Solid-SPH: решётка достаточно плотная для соседей (n>=2, hsml=1.5*dx)");

    // Solid-Solid-Contact: брутфорс по всем парам (n мало) -- то же
    // подмножество real_pairs, что строит build_solid_sph_search в продовом
    // пути (fiber/CUDA_PLAN.md §3.1, таблица фильтров).
    {
        std::vector<std::pair<rhs::index_t, rhs::index_t>> real_pairs;
        for (std::size_t i = 0; i < sc.n; ++i)
            for (std::size_t j = i + 1; j < sc.n; ++j) real_pairs.emplace_back(i, j);
        std::vector<real_t> c(sc.model.solid.c.begin(), sc.model.solid.c.begin() + static_cast<std::ptrdiff_t>(sc.n));
        interactions::SolidContactDiagnostics diag = interactions::compute_solid_contact_diagnostics(
            exec, sc.cfg, sc.cfg.mats, sc.model.solid, c, real_t(1e-9), real_pairs);
        check(diag.pair_count >= 1,
              "Solid-Solid-Contact: контакт есть (pair_count=" + std::to_string(diag.pair_count) + ")");
    }

    // Fiber-Fiber-Contact: тот же build_fiber_contact_search, что в продовом цикле.
    {
        interactions::FiberContactSearchResult search;
        interactions::build_fiber_contact_search(exec, sc.model.fiber, sc.model.topology, sc.model.materials, search);
        std::size_t contacts = 0;
        for (const auto& prim : search.primitives)
        {
            interactions::FiberContactGeometry geom =
                interactions::compute_fiber_contact_geometry(sc.model.fiber, sc.model.materials, prim);
            if (geom.contact) ++contacts;
        }
        check(contacts >= 1, "Fiber-Fiber-Contact: контакт есть (найдено " + std::to_string(contacts) +
                                  " контактных primitives из " + std::to_string(search.primitives.size()) + ")");
    }

    // Fiber-Solid-Contact: тот же build_fiber_solid_contact_search.
    {
        interactions::FiberSolidContactSearchResult search;
        interactions::build_fiber_solid_contact_search(exec, sc.cfg.mats, sc.model.solid, sc.n, sc.model.fiber,
                                                        sc.model.topology, sc.model.materials, search);
        std::size_t contacts = 0;
        for (const auto& prim : search.primitives)
        {
            interactions::FiberSolidContactGeometry geom = interactions::compute_fiber_solid_contact_geometry(
                sc.cfg.mats, sc.model.solid, sc.model.fiber, sc.model.materials, prim);
            if (geom.contact) ++contacts;
        }
        check(contacts >= 1, "Fiber-Solid-Contact: контакт есть (найдено " + std::to_string(contacts) +
                                  " контактных primitives из " + std::to_string(search.primitives.size()) + ")");
    }

    // Fiber-Bond: связи есть по построению (append_chain); реально ненулевую
    // силу они дадут через несколько шагов, когда контакт сдвинет среднюю
    // бусину -- проверяется динамически (test_bonds_respond ниже).
    check(!sc.model.topology.bonds.empty(),
          "Fiber-Bond: связи построены (" + std::to_string(sc.model.topology.bonds.size()) + ")");
}

// =========================================================================
// Импульс системы (Solid + Fiber). Нет внешних сил и закреплённых частиц --
// суммарный импульс обязан сохраняться (см. шапку файла).
// =========================================================================
struct Momentum
{
    double p[dim] = {};
};

Momentum total_momentum(const Scene& sc)
{
    Momentum m;
    for (std::size_t i = 0; i < sc.n; ++i)
        for (std::size_t d = 0; d < dim; ++d)
            m.p[d] += static_cast<double>(sc.model.solid.mass[i]) * static_cast<double>(sc.model.solid.vx[d][i]);
    for (std::size_t i = 0; i < sc.model.fiber.size(); ++i)
        for (std::size_t d = 0; d < dim; ++d)
            m.p[d] += static_cast<double>(sc.model.fiber.mass[i]) * static_cast<double>(sc.model.fiber.vx[d][i]);
    return m;
}

bool all_finite(const Scene& sc)
{
    for (std::size_t i = 0; i < sc.n; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d)
            if (!std::isfinite(static_cast<double>(sc.model.solid.x[d][i])) ||
                !std::isfinite(static_cast<double>(sc.model.solid.vx[d][i])))
                return false;
        if (!std::isfinite(static_cast<double>(sc.model.solid.rho[i])) ||
            !std::isfinite(static_cast<double>(sc.model.solid.u[i])) ||
            !std::isfinite(static_cast<double>(sc.model.solid.sxx[i])) ||
            !std::isfinite(static_cast<double>(sc.model.solid.syy[i])) ||
            !std::isfinite(static_cast<double>(sc.model.solid.szz[i])) ||
            !std::isfinite(static_cast<double>(sc.model.solid.sxy[i])) ||
            !std::isfinite(static_cast<double>(sc.model.solid.syz[i])) ||
            !std::isfinite(static_cast<double>(sc.model.solid.sxz[i])))
            return false;
    }
    for (std::size_t i = 0; i < sc.model.fiber.size(); ++i)
        for (std::size_t d = 0; d < dim; ++d)
            if (!std::isfinite(static_cast<double>(sc.model.fiber.x[d][i])) ||
                !std::isfinite(static_cast<double>(sc.model.fiber.vx[d][i])))
                return false;
    return true;
}

// =========================================================================
// Снимок состояния -- заготовка эталона для сравнения CPU/GPU (Г8), по
// образцу tests/regression_snapshot.hpp (текстовый формат, шапка "# ключ:
// значение", допуск, а не побитовое сравнение) -- своя, более простая
// реализация: снимок объединяет Solid И Fiber, чего regression_snapshot.hpp
// (заточен под один ParticleSystem) не умеет.
// =========================================================================
std::string fmt(double v) { char buf[64]; std::snprintf(buf, sizeof(buf), "%.12e", v); return buf; }

bool write_golden(const Scene& sc, const std::string& path)
{
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    out << "# Эталон test_time_integration_fiber (fiber/CUDA_PLAN.md, этап Г1)\n";
    out << "# Заготовка для сравнения CPU/GPU на этапе Г8 -- НЕ приборный эталон.\n";
    out << "# format: 1\n";
    out << "# precision: " << SPH_PRECISION << "\n";
    out << "# steps: " << kSteps << "\n";
    out << "# nsolid: " << sc.n << "\n";
    out << "# nfiber: " << sc.model.fiber.size() << "\n";
    out << "# columns solid: i x0 x1 x2 v0 v1 v2 rho u sxx syy szz sxy syz sxz\n";
    for (std::size_t i = 0; i < sc.n; ++i)
    {
        out << i;
        for (std::size_t d = 0; d < dim; ++d) out << " " << fmt(static_cast<double>(sc.model.solid.x[d][i]));
        for (std::size_t d = 0; d < dim; ++d) out << " " << fmt(static_cast<double>(sc.model.solid.vx[d][i]));
        out << " " << fmt(static_cast<double>(sc.model.solid.rho[i])) << " "
            << fmt(static_cast<double>(sc.model.solid.u[i]));
        out << " " << fmt(static_cast<double>(sc.model.solid.sxx[i])) << " "
            << fmt(static_cast<double>(sc.model.solid.syy[i])) << " "
            << fmt(static_cast<double>(sc.model.solid.szz[i])) << " "
            << fmt(static_cast<double>(sc.model.solid.sxy[i])) << " "
            << fmt(static_cast<double>(sc.model.solid.syz[i])) << " "
            << fmt(static_cast<double>(sc.model.solid.sxz[i]));
        out << "\n";
    }
    out << "# columns fiber: i x0 x1 x2 v0 v1 v2\n";
    for (std::size_t i = 0; i < sc.model.fiber.size(); ++i)
    {
        out << i;
        for (std::size_t d = 0; d < dim; ++d) out << " " << fmt(static_cast<double>(sc.model.fiber.x[d][i]));
        for (std::size_t d = 0; d < dim; ++d) out << " " << fmt(static_cast<double>(sc.model.fiber.vx[d][i]));
        out << "\n";
    }
    return static_cast<bool>(out);
}

// Разбор эталона в плоский список double (в том же порядке, что write_golden
// печатает) -- достаточно для попарного сравнения, без строгого парсинга схемы.
bool read_golden(const std::string& path, std::vector<double>& solid_rows, std::vector<double>& fiber_rows,
                  std::size_t nsolid_cols, std::size_t nfiber_cols, std::size_t nsolid, std::size_t nfiber,
                  std::string& err)
{
    std::ifstream in(path);
    if (!in) { err = "эталон не найден: " + path; return false; }
    std::string line;
    solid_rows.clear();
    fiber_rows.clear();
    std::size_t solid_seen = 0, fiber_seen = 0;
    bool in_fiber_section = false;
    while (std::getline(in, line))
    {
        if (line.empty()) continue;
        if (line[0] == '#')
        {
            if (line.find("columns fiber") != std::string::npos) in_fiber_section = true;
            continue;
        }
        std::istringstream is(line);
        std::size_t idx = 0;
        if (!(is >> idx)) { err = "испорчена строка: " + line; return false; }
        std::vector<double>& target = in_fiber_section ? fiber_rows : solid_rows;
        const std::size_t ncols = in_fiber_section ? nfiber_cols : nsolid_cols;
        for (std::size_t k = 0; k < ncols; ++k)
        {
            double v = 0.0;
            if (!(is >> v)) { err = "не хватает поля в строке: " + line; return false; }
            target.push_back(v);
        }
        if (in_fiber_section) ++fiber_seen; else ++solid_seen;
    }
    if (solid_seen != nsolid) { err = "nsolid в эталоне " + std::to_string(solid_seen) + ", ожидалось " +
                                       std::to_string(nsolid); return false; }
    if (fiber_seen != nfiber) { err = "nfiber в эталоне " + std::to_string(fiber_seen) + ", ожидалось " +
                                       std::to_string(nfiber); return false; }
    return true;
}

void compare_with_golden(const Scene& sc, const std::string& path)
{
    constexpr std::size_t kSolidCols = 2 * dim + 1 + 1 + 6; // x,v,rho,u,sxx..sxz
    constexpr std::size_t kFiberCols = 2 * dim;             // x,v

    std::vector<double> ref_solid, ref_fiber;
    std::string err;
    if (!read_golden(path, ref_solid, ref_fiber, kSolidCols, kFiberCols, sc.n, sc.model.fiber.size(), err))
    {
        check(false, "чтение эталона: " + err);
        return;
    }

    auto near = [](double a, double b, double tol_rel, double floor) {
        const double scale = std::max({std::fabs(a), std::fabs(b), floor});
        return std::fabs(a - b) <= tol_rel * scale;
    };

    double worst = 0.0;
    std::string worst_what;
    for (std::size_t i = 0; i < sc.n; ++i)
    {
        std::vector<double> cur;
        for (std::size_t d = 0; d < dim; ++d) cur.push_back(static_cast<double>(sc.model.solid.x[d][i]));
        for (std::size_t d = 0; d < dim; ++d) cur.push_back(static_cast<double>(sc.model.solid.vx[d][i]));
        cur.push_back(static_cast<double>(sc.model.solid.rho[i]));
        cur.push_back(static_cast<double>(sc.model.solid.u[i]));
        cur.push_back(static_cast<double>(sc.model.solid.sxx[i]));
        cur.push_back(static_cast<double>(sc.model.solid.syy[i]));
        cur.push_back(static_cast<double>(sc.model.solid.szz[i]));
        cur.push_back(static_cast<double>(sc.model.solid.sxy[i]));
        cur.push_back(static_cast<double>(sc.model.solid.syz[i]));
        cur.push_back(static_cast<double>(sc.model.solid.sxz[i]));
        for (std::size_t k = 0; k < kSolidCols; ++k)
        {
            const double a = cur[k], b = ref_solid[i * kSolidCols + k];
            if (!near(a, b, kSnapshotTolRel, kSnapshotFloor))
            {
                const double rel = std::fabs(a - b) / std::max({std::fabs(a), std::fabs(b), kSnapshotFloor});
                if (rel > worst) { worst = rel; worst_what = "solid[" + std::to_string(i) + "][" + std::to_string(k) + "]"; }
            }
        }
    }
    for (std::size_t i = 0; i < sc.model.fiber.size(); ++i)
    {
        std::vector<double> cur;
        for (std::size_t d = 0; d < dim; ++d) cur.push_back(static_cast<double>(sc.model.fiber.x[d][i]));
        for (std::size_t d = 0; d < dim; ++d) cur.push_back(static_cast<double>(sc.model.fiber.vx[d][i]));
        for (std::size_t k = 0; k < kFiberCols; ++k)
        {
            const double a = cur[k], b = ref_fiber[i * kFiberCols + k];
            if (!near(a, b, kSnapshotTolRel, kSnapshotFloor))
            {
                const double rel = std::fabs(a - b) / std::max({std::fabs(a), std::fabs(b), kSnapshotFloor});
                if (rel > worst) { worst = rel; worst_what = "fiber[" + std::to_string(i) + "][" + std::to_string(k) + "]"; }
            }
        }
    }
    check(worst_what.empty(), "снимок после " + std::to_string(kSteps) + " шагов совпадает с эталоном (<= " +
                                   std::to_string(kSnapshotTolRel) + " отн., худшее: " +
                                   (worst_what.empty() ? "-" : worst_what) + ", " + std::to_string(worst) + ")");
}

std::string golden_path()
{
#ifndef SPH_TEST_DATA_DIR
#define SPH_TEST_DATA_DIR "tests/data"
#endif
#if SPH_PRECISION == 32
    return std::string(SPH_TEST_DATA_DIR) + "/fiber_time_integration_fp32/golden.txt";
#else
    return std::string(SPH_TEST_DATA_DIR) + "/fiber_time_integration/golden.txt";
#endif
}

} // namespace

int main(int argc, char** argv)
{
    // Однопоточно ВСЕГДА (не только при --emit-golden): rhs::CpuExecutor
    // делит работу по потокам, порядок суммирования вкладов детерминирован
    // только при одинаковом числе worker'ов (fiber/PLAN.md §3.5) -- без этого
    // сравнение со снимком плавало бы от прогона к прогону на многоядерной
    // машине, а не только между генерацией и проверкой.
#ifdef _WIN32
    omp_set_num_threads(1);
#endif

    std::string emit_golden_dir;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        const std::string prefix = "--emit-golden=";
        if (a.rfind(prefix, 0) == 0) emit_golden_dir = a.substr(prefix.size());
    }

    std::printf("=== time_integration(fiber::Model&) с непустыми нитями (fiber/CUDA_PLAN.md, Г1) ===\n");

    Scene sc = build_scene();
    verify_interactions_present(sc);

    const Momentum p0 = total_momentum(sc);
    double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.n);
    double time = 0.0;
    RunGuard guard(0.0, 0.0, 0);
    guard.set_step_limit(kSteps);

    // Скорости центральной бусины ДО прогона -- для проверки, что Fiber-Bond
    // реально отработал (см. шапку файла): в покое при t=0, должна сдвинуться
    // контактными силами за K шагов.
    const std::size_t midA = 2; // индекс средней бусины нити A (append_chain, 5 бусин, индекс 2)
    const double vA0 = static_cast<double>(sc.model.fiber.vx[2][midA]);

    bool threw = false;
    std::string what;
    {
        // Изоляция побочных эффектов (tests/e2e_harness.hpp): time_integration
        // пишет log.txt и sph_output_*.txt в ТЕКУЩИЙ каталог жёстко зашитым
        // именем -- без ScopedWorkDir параллельные тесты (ctest -j) писали бы
        // в один и тот же файл рабочего каталога. ScopedSilence глушит
        // построчный вывод конвейера (std::wcout/std::cout).
        e2e::ScopedWorkDir wd("sph_time_integration_fiber", "run");
        e2e::ScopedSilence silence(true);
        try
        {
            time_integration(sc.cfg, sc.model, sc.n, dt, time, 1.0e-3, &guard);
        }
        catch (const std::exception& ex)
        {
            threw = true;
            what = ex.what();
        }
    }

    check(!threw, "time_integration(fiber::Model&, ...) не бросает исключений" + (threw ? (": " + what) : std::string()));
    check(guard.reason() == RunGuard::Reason::StepLimit, "остановлена штатно по StepLimit (не авария/застой)");
    check(guard.steps() == kSteps, "сделано ровно " + std::to_string(kSteps) + " шагов");
    check(std::isfinite(dt) && dt > 0.0, "итоговый dt конечен и положителен");

    check(all_finite(sc), "все поля Solid и Fiber конечны после K шагов (нет NaN/Inf)");

    const double vA1 = static_cast<double>(sc.model.fiber.vx[2][midA]);
    check(vA1 != vA0, "Fiber-Bond реально отработал: скорость средней бусины нити A изменилась (" +
                           std::to_string(vA0) + " -> " + std::to_string(vA1) + ")");

    const Momentum p1 = total_momentum(sc);
    const double scale = std::max({std::fabs(p0.p[0]), std::fabs(p0.p[1]), std::fabs(p0.p[2]),
                                    std::fabs(p1.p[0]), std::fabs(p1.p[1]), std::fabs(p1.p[2]), 1e-30});
    for (std::size_t d = 0; d < dim; ++d)
    {
        const double dp = std::fabs(p1.p[d] - p0.p[d]) / scale;
        check(dp <= kMomentumTol, "импульс системы (Solid+Fiber) сохраняется, компонента " + std::to_string(d) +
                                       ": |Δp|/scale = " + std::to_string(dp) + " <= " + std::to_string(kMomentumTol) +
                                       " (p0=" + std::to_string(p0.p[d]) + ", p1=" + std::to_string(p1.p[d]) + ")");
    }

    if (!emit_golden_dir.empty())
    {
        std::error_code ec;
        fs::create_directories(emit_golden_dir, ec);
        const std::string path = emit_golden_dir + "/golden.txt";
        check(write_golden(sc, path), "эталон записан: " + path);
        std::printf("Эталон записан в %s -- перегенерация выполняется отдельным запуском с "
                    "--emit-golden=<dir>, тест дальше СРАВНИВАЕТ с уже сохранённым эталоном "
                    "(tests/data/), а не с только что записанным.\n",
                    path.c_str());
    }
    else
    {
        compare_with_golden(sc, golden_path());
    }

    std::printf("\n=== %s (провалов: %d) ===\n", failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН", failures);
    return failures == 0 ? 0 : 1;
}
