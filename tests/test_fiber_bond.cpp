// tests/test_fiber_bond.cpp
// Unit-тесты Fiber-Bond (notes/fiber-model/fiber-plan.md, этап Ф3; notes/fiber-model/fiber-problem-statement.md, раздел
// «Fiber-Bond», «Постановка проекта: только растяжение»):
//   - сила двух бусин против формулы вручную (растяжение/сжатие);
//   - dt_limit связи против формулы, не зависит от знака деформации;
//   - частота колебаний в растянутой фазе (точный одномерный гармонический
//     осциллятор в тензии — см. вывод ниже);
//   - волна растяжения в прямой нити со скоростью sqrt(E/rho_f);
//   - сохранение суммарного импульса при интегрировании нескольких шагов.
//
// Тесты работают НАПРЯМУЮ с interactions::compute_fiber_bond /
// compute_fiber_bond_dt_limit (include/interactions/fiber_bond.hpp) — это и
// есть предмет проверки этапа Ф3. Для динамических тестов (колебания, волна,
// импульс) используется локальный PredictorCorrector-драйвер
// (FiberPcIntegrator ниже), который буквально воспроизводит формулы Solid/
// Fiber-цикла из src/time_integration.cpp (первый шаг — половина dt по
// текущей скорости, шаги n>=2 — предиктор на dt_prev/2 по производной
// прошлого шага, затем корректор на полный dt от сохранённого v^n; см.
// notes/fiber-model/fiber-problem-statement.md, «ODESolver», «Нити интегрируются той же схемой»), но не
// требует полноценной Solid-подсистемы (у реального time_integration(cfg,
// fiber::Model&, ...) поиск соседей Solid при ntotal == 0 — путь, который в
// проекте никогда не упражнялся при пустом Solid и не относится к
// Fiber-Bond; полноценная Solid+Fiber e2e связка — материал будущих этапов).
//
// Собирается с SPH_DIM=3 (CMakeLists.txt), по образцу test_fiber_model.cpp.
#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

static_assert(dim == 3, "test_fiber_bond требует сборки с SPH_DIM=3");

namespace
{

int failures = 0;

void check(bool cond, const std::string& what)
{
    std::printf("  %-5s: %s\n", cond ? "ok" : "FAIL", what.c_str());
    if (!cond) ++failures;
}

bool near_rel(double a, double b, double tol)
{
    return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

bool near_abs(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

fiber::FiberMaterial make_material(double E, double d, double rho, double eps_break = 1.0)
{
    fiber::FiberMaterial m;
    m.E = E;
    m.nu = 0.3; // не используется Fiber-Bond -- нарочно отлично от E, чтобы
                // перепутанное E<->nu в kernel сразу дало неверный результат.
    m.d = d;
    m.rho = rho;
    m.eps_break = eps_break;
    return m;
}

double cross_section_area(double d) { return pi * d * d / 4.0; }

// --- Общая двухбусинная фикстура (тесты 1 и 2) ------------------------------
struct TwoBeadSystem
{
    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
};

// Бусина 0 в начале координат, бусина 1 на расстоянии l вдоль оси x (0,0,0)/(l,0,0).
TwoBeadSystem make_two_bead_system(double l0, double l, double mi, double mj, double E, double d, double rho)
{
    TwoBeadSystem sys;
    sys.materials.push_back(make_material(E, d, rho));
    sys.state.resize(2);
    for (std::size_t dcomp = 0; dcomp < dim; ++dcomp)
    {
        sys.state.x[dcomp][0] = 0.0;
        sys.state.x[dcomp][1] = 0.0;
        sys.state.vx[dcomp][0] = 0.0;
        sys.state.vx[dcomp][1] = 0.0;
    }
    sys.state.x[0][1] = l;
    sys.state.mass[0] = mi;
    sys.state.mass[1] = mj;
    sys.state.object_id[0] = sys.state.object_id[1] = 0;
    sys.state.material_id[0] = sys.state.material_id[1] = 0;

    sys.topology.resize(2);
    sys.topology.next[0] = 1;
    sys.topology.prev[1] = 0;
    fiber::FiberBond b;
    b.i = 0;
    b.j = 1;
    b.l0 = l0;
    b.l_break = l0 * 1e9; // не используется Ф3 (разрыв -- Ф7); заведомо не сработает
    sys.topology.bonds.push_back(b);
    return sys;
}

// -----------------------------------------------------------------------
// 1) Сила двух бусин против формулы вручную; при сжатии сила равна нулю.
// -----------------------------------------------------------------------
void test_force_formula()
{
    std::printf("test_force_formula:\n");
    const double l0 = 0.01, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double mi = 2.0e-3, mj = 3.0e-3; // разные массы -- третий закон Ньютона на a=F/m, не на F
    const double A = cross_section_area(d);

    // --- Растяжение: l = 1.02*l0 -----------------------------------------
    {
        const double l = 1.02 * l0;
        TwoBeadSystem sys = make_two_bead_system(l0, l, mi, mj, E, d, rho);

        rhs::CpuExecutor exec(2);
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
        std::vector<fiber::FiberDelta> out =
            interactions::compute_fiber_bond(exec, sys.state, sys.topology, sys.materials, buffers);

        const double strain = l / l0 - 1.0;
        const double F = E * A * strain; // F_i<-j вдоль x (xij = (l,0,0))

        check(near_rel(out[0].dv[0], F / mi, 1e-9), "растяжение: dv_x бусины 0 == F/m0");
        check(near_rel(out[1].dv[0], -F / mj, 1e-9), "растяжение: dv_x бусины 1 == -F/m1");
        check(near_abs(out[0].dv[1], 0.0, 1e-12) && near_abs(out[0].dv[2], 0.0, 1e-12),
              "растяжение: поперечные компоненты dv бусины 0 нулевые");
        check(near_abs(out[1].dv[1], 0.0, 1e-12) && near_abs(out[1].dv[2], 0.0, 1e-12),
              "растяжение: поперечные компоненты dv бусины 1 нулевые");
        check(mi * out[0].dv[0] + mj * out[1].dv[0] < 1e-6, "растяжение: сумма m*dv по x равна нулю (3-й закон Ньютона)");
    }

    // --- Сжатие: l = 0.98*l0 -> сила ровно ноль ----------------------------
    {
        const double l = 0.98 * l0;
        TwoBeadSystem sys = make_two_bead_system(l0, l, mi, mj, E, d, rho);

        rhs::CpuExecutor exec(2);
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
        std::vector<fiber::FiberDelta> out =
            interactions::compute_fiber_bond(exec, sys.state, sys.topology, sys.materials, buffers);

        for (std::size_t dcomp = 0; dcomp < dim; ++dcomp)
        {
            check(out[0].dv[dcomp] == 0.0, "сжатие: dv бусины 0 == 0 (компонента " + std::to_string(dcomp) + ")");
            check(out[1].dv[dcomp] == 0.0, "сжатие: dv бусины 1 == 0 (компонента " + std::to_string(dcomp) + ")");
        }
    }

    // --- Ровно на пороге l == l0 -> сила тоже ноль (условие l > l0, не >=) --
    {
        TwoBeadSystem sys = make_two_bead_system(l0, l0, mi, mj, E, d, rho);
        rhs::CpuExecutor exec(1);
        std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
        std::vector<fiber::FiberDelta> out =
            interactions::compute_fiber_bond(exec, sys.state, sys.topology, sys.materials, buffers);
        check(out[0].dv[0] == 0.0 && out[1].dv[0] == 0.0, "l == l0: сила равна нулю (только строгое растяжение)");
    }
}

// -----------------------------------------------------------------------
// dt_limit связи: 0.8/sqrt(k_eff/m_red), k_eff = EA/l0, m_red = mi*mj/(mi+mj);
// не зависит от текущей деформации (растяжение/сжатие дают одно и то же).
// -----------------------------------------------------------------------
void test_dt_limit_formula()
{
    std::printf("test_dt_limit_formula:\n");
    const double l0 = 0.01, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double mi = 2.0e-3, mj = 3.0e-3;
    const double A = cross_section_area(d);
    const double k_eff = E * A / l0;
    const double m_red = (mi * mj) / (mi + mj);
    const double expected = 0.8 / std::sqrt(k_eff / m_red);

    rhs::CpuExecutor exec(2);

    TwoBeadSystem stretched = make_two_bead_system(l0, 1.05 * l0, mi, mj, E, d, rho);
    double dt_stretched = interactions::compute_fiber_bond_dt_limit(exec, stretched.state, stretched.topology, stretched.materials);
    check(near_rel(dt_stretched, expected, 1e-9), "dt_limit (растянутая связь) совпадает с формулой");

    TwoBeadSystem compressed = make_two_bead_system(l0, 0.90 * l0, mi, mj, E, d, rho);
    double dt_compressed = interactions::compute_fiber_bond_dt_limit(exec, compressed.state, compressed.topology, compressed.materials);
    check(near_rel(dt_compressed, expected, 1e-9), "dt_limit (сжатая связь) совпадает с формулой");
    check(dt_stretched == dt_compressed, "dt_limit не зависит от знака деформации");

    fiber::FiberTopology empty_topology;
    empty_topology.resize(2);
    double dt_empty = interactions::compute_fiber_bond_dt_limit(exec, stretched.state, empty_topology, stretched.materials);
    check(dt_empty == std::numeric_limits<double>::max(), "dt_limit без связей == +inf");
}

// -----------------------------------------------------------------------
// Локальный PredictorCorrector-драйвер для fiber::FiberState -- см. заголовок
// файла: буквально формулы src/time_integration.cpp (Solid/Fiber-цикл, Ф3),
// применённые к тестовой системе без Solid. dt фиксирован по шагам (тестам
// удобно), что не противоречит формулам -- production-цикл лишь пересчитывает
// dt каждый шаг заново тем же способом (min CFL/contact/fiber_dt_limit).
// -----------------------------------------------------------------------
struct FiberPcIntegrator
{
    rhs::CpuExecutor exec{1};
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
    std::vector<double> dv_prev[dim];
    std::vector<double> v_min[dim];
    bool first = true;

    void step(fiber::FiberState& fb, const fiber::FiberTopology& topo, const std::vector<fiber::FiberMaterial>& mats,
              double dt)
    {
        const std::size_t n = fb.size();
        if (dv_prev[0].size() != n)
        {
            for (std::size_t d = 0; d < dim; ++d)
            {
                dv_prev[d].assign(n, 0.0);
                v_min[d].assign(n, 0.0);
            }
        }

        if (!first)
        {
            for (std::size_t d = 0; d < dim; ++d)
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    v_min[d][i] = fb.vx[d][i];
                    fb.vx[d][i] += 0.5 * dt * dv_prev[d][i];
                }
            }
        }

        std::vector<fiber::FiberDelta> delta = interactions::compute_fiber_bond(exec, fb, topo, mats, buffers);
        for (std::size_t d = 0; d < dim; ++d)
            for (std::size_t i = 0; i < n; ++i) dv_prev[d][i] = delta[i].dv[d];

        const double coeff = first ? dt * 0.5 : dt;
        for (std::size_t d = 0; d < dim; ++d)
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                const double base = first ? fb.vx[d][i] : v_min[d][i];
                fb.vx[d][i] = base + coeff * dv_prev[d][i];
                fb.x[d][i] += dt * fb.vx[d][i];
            }
        }
        first = false;
    }
};

// -----------------------------------------------------------------------
// 2) Частота колебаний в растянутой фазе.
//
// Вывод (одномерный случай, движение строго вдоль оси связи): при l > l0
// сила F = EA*(l/l0 - 1) = EA*s/l0, где s = l - l0 -- РОВНО линейна по s (l в
// коэффициенте направления x_ij/l сокращается с l в числителе x_ij), поэтому
// относительное удлинение s(t) точно удовлетворяет
//     s'' = -(EA/(l0*m_red)) * s = -omega^2 * s,   omega^2 = k/m_red,
// без приближения малых s -- ПОКА l(t) > l0 (иначе сила обнуляется и формула
// перестаёт быть корректной). Освобождение из состояния s=s0, v=0 даёт
// s(t) = s0*cos(omega t), корректно вплоть до t = T/4 (s обращается в 0).
// Берём t1 = T/8 << T/4 с большим запасом.
// -----------------------------------------------------------------------
void test_oscillation_frequency()
{
    std::printf("test_oscillation_frequency:\n");
    const double l0 = 0.01, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double m = 1.0e-6; // равные массы -> m_red = m/2
    const double A = cross_section_area(d);
    const double k = E * A / l0;
    const double m_red = m / 2.0;
    const double omega = std::sqrt(k / m_red);
    const double T = 2.0 * pi / omega;

    const double s0 = l0 * 1.0e-3; // малое начальное растяжение (не принципиально -- формула точна при любом s0>0)
    TwoBeadSystem sys = make_two_bead_system(l0, l0 + s0, m, m, E, d, rho);

    const double dt = T / 2000.0;
    const std::size_t nsteps = static_cast<std::size_t>(std::round((T / 8.0) / dt));

    FiberPcIntegrator integ;
    for (std::size_t k_step = 0; k_step < nsteps; ++k_step)
        integ.step(sys.state, sys.topology, sys.materials, dt);

    const double t1 = static_cast<double>(nsteps) * dt;
    const double l_sim = sys.state.x[0][1] - sys.state.x[0][0];
    const double s_sim = l_sim - l0;
    const double v_rel_sim = sys.state.vx[0][1] - sys.state.vx[0][0];

    const double s_exact = s0 * std::cos(omega * t1);
    const double v_exact = -s0 * omega * std::sin(omega * t1);

    check(s_exact > 0.0, "самопроверка: t1 < T/4, s_exact ещё положительно (растяжение)");
    check(near_rel(s_sim, s_exact, 5e-4), "удлинение s(t1) совпадает с аналитикой гармонического осциллятора");
    // Схема (notes/fiber-model/fiber-problem-statement.md, «ODESolver») вычисляет силу по координате ДО
    // текущего шага (сила от скорости не зависит), т.е. эквивалентна
    // симплектическому Эйлеру для Fiber-only системы -- скорость несёт
    // небольшой фазовый сдвиг относительно координаты того же порядка dt,
    // поэтому допуск для неё шире, чем для s(t1).
    check(near_rel(v_rel_sim, v_exact, 3e-3), "относительная скорость v(t1) близка к аналитике (допуск шире -- фазовый сдвиг схемы)");
}

// -----------------------------------------------------------------------
// 3) Волна растяжения в прямой нити идёт со скоростью sqrt(E/rho_f).
//
// Прямая нить из Nbeads бусин с шагом l0, все бусины одинаковой массы
// m = rho*A*l0. Бусина 0 получает импульс СКОРОСТИ, направленный ПРОЧЬ от
// нити (растягивает связь (0,1)) -- "резкое растяжение с одного конца"
// (notes/fiber-model/fiber-problem-statement.md). Фронт возмущения отслеживается по времени первого
// превышения порога скорости на двух зондовых бусинах, средняя скорость
// фронта между зондами сравнивается с c = sqrt(E/rho).
// -----------------------------------------------------------------------
void test_wave_speed()
{
    std::printf("test_wave_speed:\n");
    const std::size_t Nbeads = 300;
    const double l0 = 1.0e-3, E = 2.0e9, d = 1.0e-3, rho = 1200.0;
    const double A = cross_section_area(d);
    const double m = rho * A * l0;
    const double c = std::sqrt(E / rho);

    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    materials.push_back(make_material(E, d, rho));

    state.resize(Nbeads);
    topology.resize(Nbeads);
    for (std::size_t i = 0; i < Nbeads; ++i)
    {
        for (std::size_t dcomp = 0; dcomp < dim; ++dcomp)
        {
            state.x[dcomp][i] = 0.0;
            state.vx[dcomp][i] = 0.0;
        }
        state.x[0][i] = static_cast<double>(i) * l0;
        state.mass[i] = m;
        state.object_id[i] = 0;
        state.material_id[i] = 0;
        if (i + 1 < Nbeads)
        {
            topology.next[i] = static_cast<long long>(i + 1);
            topology.prev[i + 1] = static_cast<long long>(i);
            fiber::FiberBond b;
            b.i = i;
            b.j = i + 1;
            b.l0 = l0;
            b.l_break = l0 * 1e9;
            topology.bonds.push_back(b);
        }
    }

    // Импульс: бусина 0 уходит в отрицательном направлении x -> связь (0,1)
    // растягивается (l = x1-x0 растёт).
    const double v_pulse = 1.0; // м/с, много меньше c -- линейный режим
    state.vx[0][0] = -v_pulse;

    rhs::CpuExecutor exec_probe(2);
    const double dt_limit = interactions::compute_fiber_bond_dt_limit(exec_probe, state, topology, materials);
    const double dt = 0.25 * dt_limit; // с запасом ниже предела устойчивости

    const std::size_t probeA = 40;
    const std::size_t probeB = 220;
    const double threshold = 0.02 * v_pulse; // порог первого заметного отклика

    const double target_time = 1.3 * (static_cast<double>(probeB) * l0) / c; // время с запасом, чтобы фронт дошёл до probeB
    const std::size_t max_steps = static_cast<std::size_t>(target_time / dt) + 10;

    FiberPcIntegrator integ;
    double t = 0.0;
    double tA = -1.0, tB = -1.0;
    double prev_vA = 0.0, prev_vB = 0.0, prev_t = 0.0;
    for (std::size_t s = 0; s < max_steps && (tA < 0.0 || tB < 0.0); ++s)
    {
        integ.step(state, topology, materials, dt);
        t += dt;
        const double vA = std::fabs(state.vx[0][probeA]);
        const double vB = std::fabs(state.vx[0][probeB]);
        if (tA < 0.0 && vA >= threshold)
        {
            // линейная интерполяция момента пересечения порога между (prev_t, prev_vA) и (t, vA)
            const double frac = (vA > prev_vA) ? (threshold - prev_vA) / (vA - prev_vA) : 0.0;
            tA = prev_t + frac * (t - prev_t);
        }
        if (tB < 0.0 && vB >= threshold)
        {
            const double frac = (vB > prev_vB) ? (threshold - prev_vB) / (vB - prev_vB) : 0.0;
            tB = prev_t + frac * (t - prev_t);
        }
        prev_vA = vA;
        prev_vB = vB;
        prev_t = t;
    }

    check(tA >= 0.0, "фронт достиг зондовой бусины A за отведённое время");
    check(tB >= 0.0, "фронт достиг зондовой бусины B за отведённое время");
    if (tA >= 0.0 && tB >= 0.0 && tB > tA)
    {
        const double distance = static_cast<double>(probeB - probeA) * l0;
        const double c_measured = distance / (tB - tA);
        std::printf("    c_measured = %.3f m/s, c = sqrt(E/rho) = %.3f m/s\n", c_measured, c);
        check(near_rel(c_measured, c, 0.20), "измеренная скорость фронта близка к sqrt(E/rho) (20%)");
    }
}

// -----------------------------------------------------------------------
// 4) Сумма импульсов сохраняется на нескольких шагах интегрирования.
// -----------------------------------------------------------------------
void test_momentum_conservation()
{
    std::printf("test_momentum_conservation:\n");
    const double l0 = 0.02, E = 1.5e9, d = 8.0e-4, rho = 1000.0;
    const std::size_t N = 5;
    const std::vector<double> masses = {1.0e-3, 1.6e-3, 0.9e-3, 2.2e-3, 1.3e-3}; // разные массы

    fiber::FiberState state;
    fiber::FiberTopology topology;
    std::vector<fiber::FiberMaterial> materials;
    materials.push_back(make_material(E, d, rho));

    state.resize(N);
    topology.resize(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t dcomp = 0; dcomp < dim; ++dcomp)
        {
            state.x[dcomp][i] = 0.0;
            state.vx[dcomp][i] = 0.0;
        }
        // Небольшое неравномерное растяжение каждого сегмента (все > l0),
        // чтобы связи были в тензии и создавали нетривиальные силы.
        const double stretch = 1.0 + 0.01 * static_cast<double>(i + 1) / static_cast<double>(N);
        state.x[0][i] = (i == 0) ? 0.0 : (state.x[0][i - 1] + l0 * stretch);
        state.mass[i] = masses[i];
        state.object_id[i] = 0;
        state.material_id[i] = 0;
        // Ненулевые начальные скорости в разных направлениях -- не только
        // вдоль оси связи, чтобы проверка не была тривиальной по симметрии.
        state.vx[0][i] = 0.1 * static_cast<double>(i) - 0.2;
        state.vx[1][i] = (i % 2 == 0) ? 0.05 : -0.03;
        if (dim > 2) state.vx[2][i] = 0.02 * static_cast<double>(i);

        if (i + 1 < N)
        {
            topology.next[i] = static_cast<long long>(i + 1);
            topology.prev[i + 1] = static_cast<long long>(i);
            fiber::FiberBond b;
            b.i = i;
            b.j = i + 1;
            b.l0 = l0;
            b.l_break = l0 * 1e9;
            topology.bonds.push_back(b);
        }
    }

    auto total_momentum = [&](double p[dim]) {
        for (std::size_t d = 0; d < dim; ++d) p[d] = 0.0;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t d = 0; d < dim; ++d) p[d] += state.mass[i] * state.vx[d][i];
    };

    double p0[dim];
    total_momentum(p0);

    rhs::CpuExecutor exec_probe(1);
    const double dt_limit = interactions::compute_fiber_bond_dt_limit(exec_probe, state, topology, materials);
    const double dt = 0.5 * dt_limit;

    FiberPcIntegrator integ;
    const int nsteps = 25;
    for (int s = 0; s < nsteps; ++s) integ.step(state, topology, materials, dt);

    double p1[dim];
    total_momentum(p1);

    double p0_scale = 0.0;
    for (std::size_t d = 0; d < dim; ++d) p0_scale = std::max(p0_scale, std::fabs(p0[d]));
    if (p0_scale < 1e-12) p0_scale = 1.0;

    for (std::size_t d = 0; d < dim; ++d)
    {
        check(std::fabs(p1[d] - p0[d]) < 1e-9 * p0_scale,
              "импульс сохраняется (компонента " + std::to_string(d) + "): p0=" + std::to_string(p0[d]) +
                  " p_after=" + std::to_string(p1[d]));
    }
}

} // namespace

int main()
{
    test_force_formula();
    test_dt_limit_formula();
    test_oscillation_frequency();
    test_wave_speed();
    test_momentum_conservation();

    if (failures == 0)
        std::printf("test_fiber_bond: OK\n");
    else
        std::printf("test_fiber_bond: %d failure(s)\n", failures);

    return failures;
}
