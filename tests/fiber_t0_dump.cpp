// tests/fiber_t0_dump.cpp
// Вспомогательная программа Т0 (notes/fiber-model/fiber-plasticity-design.md, раздел 8):
// печатает с %.17g состояние нитей после фиксированного числа шагов, чтобы
// побитово сравнить результат ДО и ПОСЛЕ внедрения пластичности (этап 2).
// ИНСТРУМЕНТ, не тест: не регистрируется в ctest и собирается только по
// запросу (EXCLUDE_FROM_ALL: cmake --build <dir> --target fiber_t0_dump).
// Побитовое сравнение: собрать эту же программу (без аргумента "inelastic") на
// старом дереве (до этапа 2) и на новом, запускать с OMP_NUM_THREADS=1 и
// сравнить выводы. Для этого упругая часть использует только API, существовавший
// ДО этапа 2 (FiberBond{i,j,l0,l_break}, FiberMaterial{E,nu,d,rho,eps_break},
// find/apply_fiber_break_events); поля st/plastic/visco трогает только режим
// "inelastic" (аргумент 2), который старым деревом не собирается.
//   fiber_t0_dump [nsteps=150] [inelastic]
//   (A) цепочка нити под Fiber-Bond с разрывом (прямой драйвер PC + find/apply break);
//   (B) боевой time_integration на сцене Guo 2019, конфигурация 1T, N шагов.
#include "guo2019_monitor.hpp"
#include "guo2019_scene.hpp"

#include "e2e_harness.hpp"

#include "interactions/fiber_bond.hpp"
#include "interactions/postprocessing.hpp"
#include "run_guard.hpp"
#include "time_integration.hpp"

#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

namespace
{

void dump_fiber(const char* tag, const fiber::FiberState& fb, const fiber::FiberTopology& topo)
{
    std::printf("[%s] n=%zu bonds=%zu version=%llu\n", tag, fb.size(), topo.bonds.size(), topo.version);
    for (std::size_t i = 0; i < fb.size(); ++i)
        std::printf("%zu %d %.17g %.17g %.17g %.17g %.17g %.17g\n", i, fb.object_id[i], fb.x[0][i], fb.x[1][i],
                    fb.x[2][i], fb.vx[0][i], fb.vx[1][i], fb.vx[2][i]);
}

void part_chain()
{
    const std::size_t N = 12;
    const double l0 = 1.0e-3, E = 2.0e9, d = 1.0e-4, rho = 1140.0;
    fiber::FiberMaterial mat{E, 0.35, d, rho, 0.05};
    std::vector<fiber::FiberMaterial> mats{mat};
    fiber::FiberState st;
    fiber::FiberTopology topo;
    st.resize(N);
    topo.resize(N);
    const double A = pi * d * d / 4.0;
    for (std::size_t i = 0; i < N; ++i)
    {
        st.x[0][i] = static_cast<double>(i) * l0 * (1.0 + 0.004 * static_cast<double>(i % 3));
        st.x[1][i] = 1e-5 * static_cast<double>(i % 4);
        st.x[2][i] = 0.0;
        st.vx[0][i] = 3.0 * static_cast<double>(i) - 10.0;
        st.vx[1][i] = (i % 2) ? 1.0 : -2.0;
        st.vx[2][i] = 0.5 * static_cast<double>(i % 5);
        st.mass[i] = rho * A * l0 * (1.0 + 0.1 * static_cast<double>(i % 3));
        st.object_id[i] = 0;
        st.material_id[i] = 0;
        if (i + 1 < N)
        {
            topo.next[i] = static_cast<long long>(i + 1);
            topo.prev[i + 1] = static_cast<long long>(i);
            topo.bonds.push_back(fiber::FiberBond{i, i + 1, l0, l0 * 1.05});
        }
    }
    rhs::CpuExecutor exec(2);
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>> buffers;
    std::vector<std::vector<interactions::FiberBreakEvent>> brk;
    std::vector<double> dvp[dim];
    for (std::size_t dd = 0; dd < dim; ++dd) dvp[dd].assign(N, 0.0);
    double dt_lim_sum = 0.0;
    for (int s = 0; s < 400; ++s)
    {
        const double dt = 0.4 * interactions::compute_fiber_bond_dt_limit(exec, st, topo, mats);
        dt_lim_sum += dt;
        std::vector<fiber::FiberDelta> dl = interactions::compute_fiber_bond(exec, st, topo, mats, buffers);
        for (std::size_t dd = 0; dd < dim; ++dd)
            for (std::size_t i = 0; i < N; ++i)
            {
                st.vx[dd][i] += dt * dl[i].dv[dd];
                st.x[dd][i] += dt * st.vx[dd][i];
            }
        // Поиск + применение разрыва напрямую (API, существовавший до этапа 2):
        // та же последовательность, что делала прежняя обёртка.
        interactions::apply_fiber_break_events(st, topo, interactions::find_fiber_break_events(exec, st, topo, brk));
    }
    std::printf("dt_sum=%.17g\n", dt_lim_sum);
    dump_fiber("chain", st, topo);
}

void part_guo(std::size_t nsteps, bool inelastic)
{
    namespace fs = std::filesystem;
    e2e::ScopedWorkDir wd("sph_t0_dump", "1T");
    guo2019::GuoScene sc;
    {
        e2e::ScopedSilence silence(true);
        sc = guo2019::build_guo_scene(guo2019::Config::T1, guo2019::Level::Smoke, guo2019::kSmokeV0, wd.dir(),
                                      guo2019::kMuDefault, true);
    }
    sc.cfg.t_total = guo2019::kSmokeTTotal;
    sc.cfg.save_t = guo2019::kSmokeSaveT;
    sc.cfg.print_step = 1000000000;
    if (inelastic) // только для проверки прохода в боевом цикле (с ним побитовое сравнение не делается)
        for (auto& m : sc.model.materials)
        {
            m.plastic.sy = 1.0e6;
            m.plastic.H = 0.1e9;
            m.visco.xi = 0.05;
            m.visco.E1 = 0.2e9;
            m.visco.tau = 1.0e-5;
        }
    double dt = compute_cfl_timestep(sc.cfg, sc.model.solid, sc.ntotal), time = 0.0;
    RunGuard guard(100.0, 2.0e6, 200);
    guard.set_step_limit(nsteps);
    {
        e2e::ScopedSilence silence(true);
        time_integration(sc.cfg, sc.model, sc.ntotal, dt, time, sc.spec.dx_s, &guard);
    }
    std::printf("guo steps=%zu dt=%.17g time=%.17g\n", guard.steps(), dt, time);
    dump_fiber("guo_fiber", sc.model.fiber, sc.model.topology);
    if (inelastic)
    {
        double ep = 0.0, wp = 0.0, wv = 0.0, epmax = 0.0;
        for (const auto& b : sc.model.topology.bonds)
        {
            ep += b.st.eps_p;
            epmax = std::max(epmax, b.st.eps_p);
            wp += b.st.w_p;
            wv += b.st.w_v;
        }
        std::printf("inelastic: sum eps_p=%.9g max eps_p=%.9g sum w_p=%.9g sum w_v=%.9g bonds=%zu\n", ep, epmax, wp, wv,
                    sc.model.topology.bonds.size());
    }
    double sx = 0.0, sv = 0.0;
    for (std::size_t i = 0; i < sc.ntotal; ++i)
        for (std::size_t dd = 0; dd < dim; ++dd)
        {
            sx += sc.model.solid.x[dd][i] * static_cast<double>(i + 1);
            sv += sc.model.solid.vx[dd][i] * static_cast<double>(i + 1);
        }
    std::printf("solid_checksum x=%.17g v=%.17g\n", sx, sv);
}

} // namespace

int main(int argc, char** argv)
{
    const std::size_t nsteps = argc > 1 ? static_cast<std::size_t>(std::atoi(argv[1])) : 150;
    const bool inelastic = argc > 2 && std::string(argv[2]) == "inelastic";
    if (!inelastic) part_chain();
    part_guo(nsteps, inelastic);
    return 0;
}
