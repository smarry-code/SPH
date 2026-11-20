#pragma once
// fiber/fiber_system.hpp
// Модель данных нити (notes/fiber-model/fiber-plan.md, этап Ф2, раздел «3.1. Данные»):
// FiberState (переменные состояния бусин), FiberTopology (связи и соседство
// по нити), FiberMaterial (свойства материала нити).
//
// Все поля — double, хранение — Structure of Arrays, частицы всех нитей
// одного типа лежат в общих массивах (notes/fiber-model/fiber-problem-statement.md, «Организация хранения
// частиц»), по образцу ParticleSystem (particle_system.hpp).
//
// Компилируется при ЛЮБОМ SPH_DIM (заголовок — часть общего конвейера
// src/*.cpp, см. CMakeLists.txt: SPH_PIPELINE_SOURCES собирается один раз на
// каждую цель со своим SPH_DIM). Физически нить требует dim == 3
// (notes/fiber-model/fiber-problem-statement.md: «Нити поддерживаются только в трёхмерной постановке»,
// notes/fiber-model/fiber-plan.md Р9) — эту границу проверяет check_fiber_requires_dim3()
// НИЖЕ, а не сам факт компиляции структур: FiberState::x[dim]/vx[dim]
// зеркалирует ParticleSystem::x[dim]/vx[dim] и безопасен при любом dim,
// пока код обращается только к индексам [0, dim).
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#include "fiber/fiber_plasticity.hpp"
#include "param.hpp"

namespace fiber
{

// Проверка границы Р9. Бросает исключение, если код, реально использующий
// модель нитей (разбор файла нитей, построение FiberState/FiberTopology,
// вывод), вызван в сборке с SPH_DIM, отличным от 3. Вызывается ДО любого
// обращения к данным нити — как из parse_sph_config (как только в конфиге
// обнаружен путь к файлу нитей), так и из load_fiber_file/fiber_output
// повторно (defence in depth), чтобы проверка срабатывала «на этапе разбора
// конфигурации/подготовки расчёта, а не потом в середине шага»
// (notes/fiber-model/fiber-plan.md, этап Ф2, пункт 5).
//
// dim — compile-time constexpr (param.hpp, ADR-002), поэтому ветка
// dim != 3 недостижима в бинарнике, собранном с SPH_DIM=3: её проверяет
// отдельная лёгкая цель test_fiber_dim_guard, собранная БЕЗ SPH_DIM
// (dim = 1 по умолчанию) — см. CMakeLists.txt и
// tests/test_fiber_dim_guard.cpp.
inline void check_fiber_requires_dim3()
{
    if constexpr (dim != 3)
    {
        throw std::runtime_error(
            "Нити (fiber_file) заданы в конфигурации, но сборка выполнена с "
            "SPH_DIM=" + std::to_string(dim) +
            " -- нити поддерживаются только при SPH_DIM=3 (notes/fiber-model/fiber-problem-statement.md, "
            "раздел \"Fiber\"; notes/fiber-model/fiber-plan.md, решение Р9).");
    }
}

// --- FiberState -------------------------------------------------------------
// Переменные состояния бусин нити: координаты, скорости, постоянная масса
// (из ненагруженной конфигурации, notes/fiber-model/fiber-problem-statement.md), номер нити (меняется при
// разрыве, этап Ф7) и индекс материала (не меняется при разрыве).
struct FiberState
{
    std::vector<double> x[dim];
    std::vector<double> vx[dim];
    std::vector<double> mass;
    std::vector<int> object_id;
    std::vector<int> material_id;

    std::size_t n = 0;

    void resize(std::size_t new_n)
    {
        n = new_n;
        for (std::size_t d = 0; d < dim; ++d)
        {
            x[d].resize(n);
            vx[d].resize(n);
        }
        mass.resize(n);
        object_id.resize(n);
        material_id.resize(n);
    }

    std::size_t size() const { return n; }
};

// --- FiberTopology -----------------------------------------------------------
// Связь топологических соседей (notes/fiber-model/fiber-problem-statement.md, Fiber-Bond, «Постановка
// проекта: только растяжение»): пара бусин i<j по цепочке нити, длина в
// ненагруженной конфигурации l0 и порог разрыва l_break = (1+eps_break)*l0.
struct FiberBond
{
    std::size_t i = 0;
    std::size_t j = 0;
    double l0 = 0.0;
    double l_break = 0.0;
    // Внутренние переменные определяющего закона (notes/fiber-model/fiber-plasticity-design.md,
    // раздел 2): eps_p, q, eps_prev, w_p, w_v. ЛЕЖАТ В КОНЦЕ структуры, чтобы
    // агрегатная инициализация FiberBond{i, j, l0, l_break} продолжала
    // компилироваться. Хранятся внутри связи, а не в отдельном массиве:
    // apply_fiber_break_events пересобирает topology.bonds копированием связей,
    // и состояние оставшихся связей переезжает вместе с ними. RHS-ядро
    // (interactions/fiber_bond.hpp) состояние только читает; пишет его проход
    // Postprocessing interactions::apply_fiber_bond_constitutive. Энергии w_p,
    // w_v -- Дж/м^3 (энергия связи = w·A·l0). Нулевое состояние -- упругая
    // ненагруженная связь.
    BondState st{};
    // Индекс материала связи. -1 (по умолчанию) -- материал бусины i, побитово
    // прежнее поведение обычных связей цепочки. >= 0 -- собственный материал
    // связи-стежка (notes/fabrics/plan-quilt-fabric.md, раздел 10): связь между бусинами
    // разных нитей, prev/next она не меняет. Лежит ПОСЛЕ st по той же причине
    // (агрегатная инициализация FiberBond{i, j, l0, l_break}).
    int mat = -1;
};

// prev[]/next[] — соседи по нити, -1 на конце (частица — конец нити, если
// prev == -1 или next == -1). bonds — каждая связь один раз. version —
// topology_version (notes/fiber-model/fiber-plan.md §3.1): растёт при изменении топологии
// (разрыв связи, этап Ф7); используется, чтобы инвалидировать кэшированные
// от старой топологии структуры (Postprocessing, notes/fiber-model/fiber-problem-statement.md).
struct FiberTopology
{
    std::vector<long long> prev;
    std::vector<long long> next;
    std::vector<FiberBond> bonds;
    unsigned long long version = 0;

    void resize(std::size_t new_n)
    {
        prev.assign(new_n, -1);
        next.assign(new_n, -1);
    }

    std::size_t size() const { return prev.size(); }
};

// Материал связи: собственный (стежок, mat >= 0) или материал бусины i.
inline int bond_material_id(const FiberBond& b, const FiberState& s)
{
    return b.mat >= 0 ? b.mat : s.material_id[b.i];
}

// Связь входит в цепочку нити (next[i]==j и prev[j]==i); иначе это связь-стежок.
inline bool is_chain_bond(const FiberBond& b, const FiberTopology& t)
{
    return t.next[b.i] == static_cast<long long>(b.j) && t.prev[b.j] == static_cast<long long>(b.i);
}

// --- FiberDelta ---------------------------------------------------------------
// Изменение состояния Fiber-частицы, суммируемое на стадии Reduce
// (notes/fiber-model/fiber-plan.md §3.3: «FiberDelta: dv[dim]»; notes/fiber-model/fiber-problem-statement.md, «Основные
// требования к interaction kernel», п.6: тип Delta определяется
// compute-структурой target частицы). Единственное поле — ускорение;
// interaction kernels (Fiber-Bond, позже Fiber-Fiber/Fiber-Solid) пишут в
// него через rhs::ContributionWriter<FiberDelta>, суммирование —
// rhs::reduce_contributions (operator+= ниже).
struct FiberDelta
{
    double dv[dim] = {};

    FiberDelta& operator+=(const FiberDelta& other)
    {
        for (std::size_t d = 0; d < dim; ++d) dv[d] += other.dv[d];
        return *this;
    }
};

// --- FiberMaterial -----------------------------------------------------------
// Свойства материала нити (notes/fiber-model/fiber-plan.md §3.1; notes/fiber-model/fiber-problem-statement.md, Fiber-Bond и
// Fiber-Solid-Contact): модуль Юнга, коэффициент Пуассона, диаметр,
// плотность, критическое удлинение разрыва, коэффициент трения нить-нить
// mu_ff (notes/fiber-contact/fiber-friction-projection.md, §2, п. 2.3). Коэффициент
// одновременно работает и как статический, и как кинематический: прилипание
// определяется достаточностью mu_ff*N для остановки касательного скольжения
// под действием других сил шага. mu_ff = 0 по умолчанию — поведение без
// трения между нитями, как раньше. Поле mu_ff — ПОСЛЕДНЕЕ в структуре:
// агрегатная инициализация вида FiberMaterial{E, nu, d, rho, eps_break}
// (используется в тестах) должна продолжать компилироваться и после
// добавления этого поля.
//
// mu_fs (notes/fiber-contact/fiber-solid-friction-implementation.md, §2) -- коэффициент
// трения этой нити о ЛЮБОЕ Solid-тело; статический и кинематический
// одновременно (та же проекция на уровне скоростей, что и у mu_ff). Берётся по
// материалу бусины a сегмента (как K_fs; правило выбора --
// interactions::fiber_solid_friction_coefficient). mu_fs = 0 по умолчанию --
// трения нить-Solid нет. mu_fs -- ПОСЛЕДНЕЕ поле (после mu_ff), чтобы не
// ломать агрегатную инициализацию вида FiberMaterial{E, nu, d, rho, eps_break,
// mu_ff}.
struct FiberMaterial
{
    double E = 0.0;
    double nu = 0.0;
    double d = 0.0;
    double rho = 0.0;
    double eps_break = 0.0;
    double mu_ff = 0.0;
    double mu_fs = 0.0; // трение этой нити о любое Solid-тело; статический и кинематический одновременно
    // Неупругость связи (notes/fiber-model/fiber-plasticity-design.md, раздел 6): пластика,
    // скоростная зависимость, вязкоупругость, критерии разрыва. ПОСЛЕ mu_fs --
    // агрегатная инициализация FiberMaterial{E, nu, d, rho, eps_break, mu_ff,
    // mu_fs} не ломается. По умолчанию всё выключено (упругий материал).
    PlasticParams plastic{};
    RateParams rate{};
    ViscoParams visco{};
    FailureParams failure{};
    // Повреждённость с разупрочнением (notes/fiber-model/fiber-damage-plan.md). По умолчанию выключена.
    DamageParams damage{};

    // Включена ли у материала хоть одна неупругая составляющая связи (пластика,
    // ветвь Максвелла, демпфер Кельвина через eta или xi, повреждённость). false =>
    // RHS-ядро и dt-лимит идут по прежней упругой формуле побитово, проход
    // Postprocessing apply_fiber_bond_constitutive связь не трогает.
    bool inelastic() const
    {
        return plastic.enabled() || (visco.E1 > 0.0 && visco.tau > 0.0) || visco.eta != 0.0 || visco.xi != 0.0 ||
               damage.enabled();
    }

    // Включён ли хоть один критерий разрыва FAILURE (EPS_P > 0 или SIGMA_U > 0).
    // Независим от inelastic(): FAILURE без PLASTIC/VISCO оставляет ядро и
    // dt-лимит на упругой быстрой ветке (st нулевой, sigma_ep = E*eps), а
    // проход разрыва проверяет только SIGMA_U (EPS_P без пластики не сработает:
    // eps_p == 0).
    bool has_failure() const
    {
        return failure.eps_p > 0.0 || failure.sigma_u > 0.0;
    }
};

} // namespace fiber
