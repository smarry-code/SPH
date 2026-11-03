// tests/test_jc_urs_consistency.cpp
// urs() (src/int_force.cpp) и solid_eos() (interactions/solid_sph.hpp) обязаны
// давать побитово одинаковые p, c для itype 12-15 (риск 2 в
// fiber/jc-implementation-plan.md): оба зовут eos_jc. Сетка (rho, u, D).
#include "int_force.hpp"
#include "interactions/solid_sph.hpp"

#include <cstdio>

using interactions::solid_eos;
using interactions::solid_itype_known;

int main()
{
    const Materials mats;
    int failures = 0, n = 0;
    for (int it = 12; it <= 15; ++it)
    {
        const JcMaterial& m = jc_material(mats, it);
        if (!solid_itype_known(it)) { std::printf("FAIL solid_itype_known(%d)\n", it); ++failures; }
        ParticleSystem sys;
        sys.resize(1);
        sys.inzializ_index(0);
        sys.itype[0] = it;
        sys.rho_0[0] = m.rho0;
        for (double mu : {-0.1, -0.02, -0.001, 0.0, 0.001, 0.02, 0.1, 0.3})
            for (double u : {-1.0e4, 0.0, 1.0e5, 1.0e6})
                for (double D : {0.0, 0.5, 1.0})
                {
                    sys.rho[0] = m.rho0 * (1.0 + mu);
                    sys.u[0] = u;
                    sys.D[0] = D;
                    urs(mats, sys, 0);
                    double p = 0.0, c = 0.0;
                    solid_eos(mats, sys, 0, p, c);
                    ++n;
                    if (!(p == sys.p[0]) || !(c == sys.c[0]))
                    {
                        std::printf("FAIL itype=%d mu=%g u=%g D=%g: urs p=%.17g c=%.17g, solid_eos p=%.17g c=%.17g\n", it,
                                    mu, u, D, sys.p[0], sys.c[0], p, c);
                        ++failures;
                    }
                }
    }
    std::printf("jc_urs_consistency: %d points, %d failures\n", n, failures);
    return failures == 0 ? 0 : 1;
}
