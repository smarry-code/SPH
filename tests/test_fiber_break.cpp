// tests/test_fiber_break.cpp
// Unit-тесты Postprocessing: разрыв связей Fiber-Bond (fiber/PLAN.md, этап
// Ф7; fiber/SPH.md, раздел «Fiber-Bond», «Разрыв связи», и «Postprocessing и
// изменение topology»):
//   - связь рвётся РОВНО при пороге l > l_break (дискретная проверка по
//     шагам -- не раньше, не позже, с точностью до шага стептинга);
//   - после разрыва две половины нити получают РАЗНЫЕ object_id и
//     обнаруживаются Fiber-Fiber-Contact (чего не было при одинаковом
//     object_id до разрыва);
//   - после разрыва Fiber-Bond RHS не создаёт силу между бывшими соседями,
//     даже если их сблизить обратно (порвавшаяся связь удалена из
//     topology.bonds -- kernel строит primitives заново на каждом вызове).
//
// Стиль -- по образцу tests/test_fiber_bond.cpp / tests/test_fiber_contact.cpp
// (тот же check()/near_rel()/near_abs(), тот же CpuExecutor напрямую).
// Собирается с SPH_DIM=3 (CMakeLists.txt).
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_contact_law.hpp"
#include "interactions/postprocessing.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_break требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near_abs(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

fiber::FiberMaterial make_material(double E, double nu, double d, double rho, double eps_break)
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
// 1) Связь рвётся РОВНО при пороге l > l_break, не раньше и не позже --
//    бусина 1 отодвигается от бусины 0 с постоянной скоростью (кинематически,
//    l(t) = l0_geom + v*t известно аналитически -- проверяем логику самого
//    порога find_fiber_break_events, без зависимости от динамики
//    Fiber-Bond), find_fiber_break_events вызывается на каждом шаге.
//    Проверяется: на шаге ДО первого срабатывания l <= l_break, на шаге
//    срабатывания l > l_break -- переход происходит ровно на границе с
//    точностью до шага дискретизации.
// -----------------------------------------------------------------------
void test_break_threshold()
{
    std::printf("test_break_threshold:\n");
    const double l0 = 0.01;
    const double eps_b = 0.05;
    const double l_break = (1.0 + eps_b) * l0;
    const double v = 0.2;        // м/с, скорость разведения бусин
    const double dt = 1.0e-4;    // шаг дискретизации Postprocessing

    fiber::FiberState state;
    fiber::FiberTopology topology;
    state.resize(2);
    topology.resize(2);
    for (std::size_t d = 0; d < dim; ++d)
    {
        state.x[d][0] = 0.0;
        state.x[d][1] = 0.0;
    }
    state.x[0][1] = l0;
    state.mass[0] = state.mass[1] = 1.0e-6;
    state.object_id[0] = state.object_id[1] = 0;
    state.material_id[0] = state.material_id[1] = 0;
    topology.next[0] = 1;
    topology.prev[1] = 0;
    fiber::FiberBond b;
    b.i = 0;
    b.j = 1;
    b.l0 = l0;
    b.l_break = l_break;
    topology.bonds.push_back(b);

    rhs::CpuExecutor exec(2);
    std::vector<std::vector<interactions::FiberBreakEvent>> buffers;

    double l_prev = l0;
    bool fired = false;
    double l_at_fire = 0.0;
    for (int step = 1; step <= 2000 && !fired; ++step)
    {
        state.x[0][1] = l0 + v * dt * static_cast<double>(step);
        const double l_now = state.x[0][1] - state.x[0][0];

        std::vector<interactions::FiberBreakEvent> events =
            interactions::find_fiber_break_events(exec, state, topology, buffers);

        if (events.empty())
        {
            check(l_now <= l_break, "нет события -> l ещё не превысило порог (шаг " + std::to_string(step) + ")");
            l_prev = l_now;
        }
        else
        {
            check(events.size() == 1, "ровно одно событие при срабатывании");
            check(events[0].bond_index == 0, "индекс связи в событии совпадает с единственной связью");
            check(l_now > l_break, "событие -> l превысило порог");
            check(l_prev <= l_break, "l на ПРЕДЫДУЩЕМ шаге ещё не превышало порог (не сработало раньше времени)");
            l_at_fire = l_now;
            fired = true;
        }
    }

    check(fired, "разрыв произошёл за отведённое число шагов");
    if (fired)
    {
        std::printf("    l_break = %.6f, l_prev = %.6f (<= l_break), l_at_fire = %.6f (> l_break)\n", l_break,
                     l_prev, l_at_fire);
    }

    // Ровно на пороге (l == l_break): события НЕТ (условие строгое l > l_break).
    {
        fiber::FiberState s2 = state;
        for (std::size_t d = 0; d < dim; ++d) s2.x[d][1] = s2.x[d][0];
        s2.x[0][1] = s2.x[0][0] + l_break;
        std::vector<interactions::FiberBreakEvent> events =
            interactions::find_fiber_break_events(exec, s2, topology, buffers);
        check(events.empty(), "l == l_break ровно -> события нет (строгое неравенство)");
    }
}

// -----------------------------------------------------------------------
// Фикстура: одна нить из 4 бусин (0-1-2-3), связь (1,2) -- та, что рвётся.
// Геометрия подобрана так, чтобы:
//   - |x1-x0| == l0 бусины (0,1) не в напряжении;
//   - |x2-x1| >> l_break(1,2) -- связь (1,2) рвётся немедленно;
//   - |x3-x2| == l0 бусины (2,3) не в напряжении;
//   - |x3-x0| мало (в пределах радиуса поиска Fiber-Fiber-Contact) -- после
//     разрыва бусины 0 и 3 разных нитей и физически близки.
// -----------------------------------------------------------------------
struct BrokenChainFixture
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

BrokenChainFixture make_chain_fixture()
{
    BrokenChainFixture fx;
    const double l0 = 0.01;
    const double d = 0.002; // диаметр нити
    fx.materials.push_back(make_material(2.0e9, 0.3, d, 1200.0, 0.05));

    fx.state.resize(4);
    fx.topology.resize(4);

    const std::array<double, 3> x0{0.0, 0.0, 0.0};
    const std::array<double, 3> x1{l0, 0.0, 0.0};
    const std::array<double, 3> x2{l0 + 0.05, 0.0, 0.0}; // |x2-x1| = 0.05 >> l_break(1,2)
    const std::array<double, 3> x3{0.001, 0.0, 0.0};     // близко к x0

    const std::array<std::array<double, 3>, 4> pts{x0, x1, x2, x3};
    for (std::size_t k = 0; k < 4; ++k)
    {
        for (std::size_t d3 = 0; d3 < dim; ++d3) fx.state.x[d3][k] = pts[k][d3];
        for (std::size_t d3 = 0; d3 < dim; ++d3) fx.state.vx[d3][k] = 0.0;
        fx.state.mass[k] = 1.0e-6;
        fx.state.object_id[k] = 0; // одна нить, все бусины -- один object_id
        fx.state.material_id[k] = 0;
    }
    fx.topology.next[0] = 1;
    fx.topology.prev[1] = 0;
    fx.topology.next[1] = 2;
    fx.topology.prev[2] = 1;
    fx.topology.next[2] = 3;
    fx.topology.prev[3] = 2;

    auto make_bond = [&](std::size_t i, std::size_t j, double l0_bond, double l_break_bond) {
        fiber::FiberBond b;
        b.i = i;
        b.j = j;
        b.l0 = l0_bond;
        b.l_break = l_break_bond;
        return b;
    };
    // (0,1): не в напряжении, разрыв заведомо не сработает.
    fx.topology.bonds.push_back(make_bond(0, 1, l0, l0 * 1e9));
    // (1,2): l0 номинально мал, но фактическое расстояние 0.05 -- триггер разрыва.
    fx.topology.bonds.push_back(make_bond(1, 2, l0, (1.0 + 0.05) * l0));
    // (2,3): не в напряжении, разрыв заведомо не сработает.
    fx.topology.bonds.push_back(make_bond(2, 3, l0, l0 * 1e9));

    return fx;
}

// -----------------------------------------------------------------------
// 2) После разрыва -- разные object_id у двух половин, и пара, которая была
//    отфильтрована по одинаковому object_id ДО разрыва, находится
//    Fiber-Fiber-Contact ПОСЛЕ разрыва.
// -----------------------------------------------------------------------
void test_break_new_object_id_and_contact()
{
    std::printf("test_break_new_object_id_and_contact:\n");
    BrokenChainFixture fx = make_chain_fixture();
    rhs::CpuExecutor exec(2);

    const unsigned long long version_before = fx.topology.version;

    // --- ДО разрыва: одна нить, все object_id одинаковы -> Fiber-Fiber-Contact
    //     не находит пару (0,3), даже физически близкую (self-contact одной
    //     нити не формулируется, fiber/SPH.md).
    {
        interactions::FiberContactSearchResult search;
        interactions::build_fiber_contact_search(exec, fx.state, fx.topology, fx.materials, search);
        check(search.primitives.empty(),
              "до разрыва: Fiber-Fiber-Contact не находит пар (общий object_id фильтрует самоконтакт)");
    }

    // --- Разрыв: поиск + применение события.
    std::vector<std::vector<interactions::FiberBreakEvent>> break_buffers;
    std::vector<interactions::FiberBreakEvent> events =
        interactions::find_fiber_break_events(exec, fx.state, fx.topology, break_buffers);
    check(events.size() == 1, "ровно одно событие разрыва (связь (1,2))");
    check(!events.empty() && events[0].bond_index == 1, "индекс события -- связь (1,2) (bond_index == 1)");

    interactions::apply_fiber_break_events(fx.state, fx.topology, events);

    // --- Топология: next[1]==-1, prev[2]==-1, связь удалена, version вырос.
    check(fx.topology.next[1] == -1, "next[1] обнулён после разрыва");
    check(fx.topology.prev[2] == -1, "prev[2] обнулён после разрыва");
    check(fx.topology.bonds.size() == 2, "порвавшаяся связь удалена из topology.bonds (осталось 2 из 3)");
    check(fx.topology.version == version_before + 1, "topology.version увеличился ровно на 1 за проход");

    // --- object_id: {0,1} -- старый (неизменный), {2,3} -- новый, общий для
    //     обеих бусин j-стороны; отличается от старого.
    check(fx.state.object_id[0] == fx.state.object_id[1], "бусины 0,1 (i-сторона) сохранили общий object_id");
    check(fx.state.object_id[2] == fx.state.object_id[3], "бусины 2,3 (j-сторона) получили ОДИН новый object_id");
    check(fx.state.object_id[0] != fx.state.object_id[2], "object_id двух половин РАЗНЫЕ после разрыва");
    check(fx.state.object_id[0] == 0, "i-сторона сохранила исходный object_id (0)");

    // --- ПОСЛЕ разрыва: Fiber-Fiber-Contact находит пару сегментов (0,1)-(2,3)
    //     (разные object_id, физически близко). В новой схеме «пара
    //     сегментов — одна сила» (fiber/fiber_contact_segment_pairs.md)
    //     primitive хранит сегменты, а не отдельные бусины: бусина 3 --
    //     последняя в своей половине (prev[3]==2 != -1), поэтому spatial-пара
    //     (0,3) отбрасывается целиком правилом «сегмент принадлежит первой
    //     бусине» -- эту зону покрывает primitive, построенный из пары (0,2),
    //     где a0=0 -- первая бусина сегмента (0,1), b0=2 -- первая бусина
    //     сегмента (2,3).
    {
        interactions::FiberContactSearchResult search;
        interactions::build_fiber_contact_search(exec, fx.state, fx.topology, fx.materials, search);
        bool found_segment_pair = false;
        for (const auto& p : search.primitives)
        {
            if ((p.a0 == 0 && p.b0 == 2) || (p.a0 == 2 && p.b0 == 0)) found_segment_pair = true;
        }
        check(!search.primitives.empty(), "после разрыва: Fiber-Fiber-Contact находит хотя бы одну пару");
        check(found_segment_pair,
              "после разрыва: пара сегментов (0,1)-(2,3) ({a0,b0}={0,2}) найдена");
    }
}

// -----------------------------------------------------------------------
// 3) После разрыва Fiber-Bond RHS не создаёт силу между бывшими соседями (1,2),
//    даже если их снова сблизить (в смысле "растянуть, как раньше делала бы
//    порвавшаяся связь") -- связь удалена из topology.bonds, kernel строит
//    primitives заново на каждом вызове (fiber/PLAN.md §3.8).
// -----------------------------------------------------------------------
void test_break_removes_bond_from_rhs()
{
    std::printf("test_break_removes_bond_from_rhs:\n");
    BrokenChainFixture fx = make_chain_fixture();
    rhs::CpuExecutor exec(2);

    std::vector<std::vector<interactions::FiberBreakEvent>> break_buffers;
    std::vector<interactions::FiberBreakEvent> events =
        interactions::find_fiber_break_events(exec, fx.state, fx.topology, break_buffers);
    interactions::apply_fiber_break_events(fx.state, fx.topology, events);
    check(fx.topology.bonds.size() == 2, "после разрыва осталось 2 связи");

    // Переставляем бусины так, чтобы ОБЕ оставшиеся связи (0,1) и (2,3) были
    // РОВНО в покое (l == l0 -> нулевая сила), а расстояние (1,2) -- такое,
    // что БЫВШАЯ связь (strain 100%), если бы она ещё существовала, дала бы
    // ОГРОМНУЮ силу растяжения. Это делает проверку "нет фантомной силы от
    // удалённой связи" содержательной, а не тривиальной.
    // Компрессия (l < l0), а не точное равенство l == l0, — гарантированно
    // нулевая сила (условие kernel'а `l > b.l0`) независимо от округления
    // double: маленький запас (0.5*l0) исключает случайное превышение из-за
    // погрешности представления при вычитании координат.
    const double l0 = 0.01;
    fx.state.x[0][0] = 0.0;
    fx.state.x[0][1] = l0;              // |x1-x0| == l0 -> связь (0,1) в покое (точно, тот же double)
    fx.state.x[0][2] = l0 + 0.02;       // |x2-x1| == 2*l0 -> strain 100% ЕСЛИ БЫ связь (1,2) существовала
    fx.state.x[0][3] = fx.state.x[0][2] + 0.5 * l0; // |x3-x2| == 0.5*l0 < l0 -> связь (2,3) в компрессии (0 силы)
    for (std::size_t k = 0; k < 4; ++k)
        for (std::size_t d = 1; d < dim; ++d) fx.state.x[d][k] = 0.0;

    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> bond_buffers;
    std::vector<fiber::FiberDelta> delta =
        interactions::compute_fiber_bond(exec, fx.state, fx.topology, fx.materials, bond_buffers);

    for (std::size_t i = 0; i < 4; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d)
        {
            check(near_abs(delta[i].dv[d], 0.0, 1e-12),
                  "delta[" + std::to_string(i) + "].dv[" + std::to_string(d) +
                      "] == 0 (нет силы -- оставшиеся связи в покое, разорванная не участвует)");
        }
    }
}

// -----------------------------------------------------------------------
// 4) FiberDtRefCache: ensure побитово равен прямым compute_fiber_bead_dt_ref /
//    compute_fiber_bond_dt_limit; пустой набор событий не меняет version (кэш
//    не пересчитывается); непустой -- version растёт, кэш пересчитывается и
//    отличается от старого у бусины, потерявшей связь.
// -----------------------------------------------------------------------
bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

void test_dt_ref_cache()
{
    std::printf("test_dt_ref_cache:\n");
    BrokenChainFixture fx = make_chain_fixture();
    fx.topology.bonds[1].l0 = 0.005; // другая l0 у связи (1,2): минимум бусин 1, 2 определяется ею
    for (int nt : {1, 2, 7, 8})
    {
        rhs::CpuExecutor exec(nt);
        interactions::FiberDtRefCache cache;
        check(!cache.valid, "кэш изначально invalid");
        cache.ensure(exec, fx.state, fx.topology, fx.materials, true);
        std::vector<double> ref;
        interactions::compute_fiber_bead_dt_ref(exec, fx.state, fx.topology, fx.materials, ref);
        const double lim = interactions::compute_fiber_bond_dt_limit(exec, fx.state, fx.topology, fx.materials);
        check(same_bits(cache.bead_dt_ref, ref) && same_bits(cache.bond_dt_limit, lim),
              "ensure == прямые вызовы побитово (потоков " + std::to_string(nt) + ")");
    }

    rhs::CpuExecutor exec(2);
    interactions::FiberDtRefCache cache;
    cache.ensure(exec, fx.state, fx.topology, fx.materials, false);
    check(cache.valid && !cache.have_bead && cache.bead_dt_ref.empty(), "need_bead=false: bead_dt_ref не считается");
    {
        const double lim = interactions::compute_fiber_bond_dt_limit(exec, fx.state, fx.topology, fx.materials);
        check(same_bits(cache.bond_dt_limit, lim), "need_bead=false: bond_dt_limit == прямой вызов");
    }
    cache.ensure(exec, fx.state, fx.topology, fx.materials, true);
    check(cache.have_bead && cache.bead_dt_ref.size() == fx.state.size(), "need_bead=true после false: bead_dt_ref досчитан");

    // Пустой набор событий: version не меняется, кэш остаётся валидным и тем же.
    const unsigned long long v0 = fx.topology.version;
    const std::vector<double> old_ref = cache.bead_dt_ref;
    const double old_lim = cache.bond_dt_limit;
    interactions::apply_fiber_break_events(fx.state, fx.topology, {});
    check(fx.topology.version == v0, "пустой набор событий: topology.version не меняется");
    cache.ensure(exec, fx.state, fx.topology, fx.materials, true);
    check(same_bits(cache.bead_dt_ref, old_ref) && same_bits(cache.bond_dt_limit, old_lim),
          "пустой набор: кэш прежний");

    // Разрыв связи (1,2).
    std::vector<std::vector<interactions::FiberBreakEvent>> bufs;
    std::vector<interactions::FiberBreakEvent> events = interactions::find_fiber_break_events(exec, fx.state, fx.topology, bufs);
    check(events.size() == 1, "разрыв: одно событие");
    interactions::apply_fiber_break_events(fx.state, fx.topology, events);
    check(fx.topology.version == v0 + 1, "непустой набор: version выросла на 1");
    cache.ensure(exec, fx.state, fx.topology, fx.materials, true);
    std::vector<double> ref;
    interactions::compute_fiber_bead_dt_ref(exec, fx.state, fx.topology, fx.materials, ref);
    const double lim = interactions::compute_fiber_bond_dt_limit(exec, fx.state, fx.topology, fx.materials);
    check(same_bits(cache.bead_dt_ref, ref) && same_bits(cache.bond_dt_limit, lim),
          "после разрыва: кэш пересчитан и == прямому вызову побитово");
    check(cache.version == fx.topology.version && cache.nbonds == fx.topology.bonds.size(), "ключ кэша обновлён");
    check(!same_bits(cache.bead_dt_ref[1], old_ref[1]) && !same_bits(cache.bead_dt_ref[2], old_ref[2]),
          "после разрыва dt_ref бусин 1 и 2 отличается от прежнего (потеряли короткую связь)");

    cache.invalidate();
    check(!cache.valid, "invalidate(): кэш invalid");
}

} // namespace

int main()
{
    test_break_threshold();
    test_break_new_object_id_and_contact();
    test_break_removes_bond_from_rhs();
    test_dt_ref_cache();

    if (failures == 0)
        std::printf("test_fiber_break: OK\n");
    else
        std::printf("test_fiber_break: %d failure(s)\n", failures);

    return failures;
}
