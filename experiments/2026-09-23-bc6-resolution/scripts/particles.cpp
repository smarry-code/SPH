// Число частиц по телам для каждого шага решётки лестницы — той же функцией
// bc6::plan_geometry, которой пользуется test_e2e_bc6. Прогона не требует;
// нужен потому, что оборванный прогон (0,10 мм) отчёт о геометрии не печатает.
//
// Сборка и запуск из корня репозитория (см. REPRO.md):
//   clang++ -std=c++20 -O1 -DSPH_DIM=2 -Iinclude -Itests \
//       experiments/2026-09-23-bc6-resolution/scripts/particles.cpp -o particles.exe
//   ./particles.exe > experiments/2026-09-23-bc6-resolution/data/particles.txt
#include <cstdio>
#include "bc6_geometry.hpp"

int main()
{
    const double dxs[] = {4.0e-4, 2.8e-4, 2.0e-4, 1.4e-4, 1.0e-4, 0.5e-4};
    const char* names[] = {"пена", "ударник", "образец", "окно"};
    std::printf("# Число частиц BC6 по bc6::plan_geometry (tests/bc6_geometry.hpp)\n");
    std::printf("# Δx,мм  всего   | по телам: n_r × n_z = частиц | монитор: индекс (с нуля), r, мм\n");
    for (double dx : dxs)
    {
        const bc6::Geometry g = bc6::plan_geometry(dx);
        std::printf("%6.2f %8zu  |", dx * 1e3, g.ntotal);
        for (size_t b = 0; b < 4; ++b)
            std::printf(" %s %zu×%zu=%zu;", names[b], g.n_r[b], g.n_z[b], g.count[b]);
        std::printf(" | монитор #%zu, r = %.3f\n", g.moni_axis, g.moni_axis_r * 1e3);
    }
    return 0;
}
