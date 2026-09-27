// tests/cuda/physics_compile_probe.cu
// Пробная CUDA-цель этапа Г2 (fiber/CUDA_PLAN.md §3.2, §4 Г2): проверяет,
// что ВСЕ заголовки include/physics/* компилируются и хост-компилятором
// nvcc (cl на Windows -- решение Р-Г1), и Clang (обычные .cpp проекта, как
// прочие CPU-тесты), и что функции реально инстанцируются и вызываются из
// __global__-функции (не только парсятся под __CUDACC__) в обеих точностях
// (SPH_PRECISION=64/32, как test_cuda_toolchain).
//
// Это НЕ тест физики Solid-SPH (для него в Г4 будет cuda_rhs_equivalence_3d
// на реальных сценах): здесь -- один синтетический вызов sph_kernel_w_dwdx
// и compute_solid_sph_pair на GPU против ТЕХ ЖЕ функций на CPU (SPH_HD
// разворачивается в обычный inline и на хосте), с допуском на ULP/FMA
// (план §5, "Одна SPH_HD-функция: CPU vs GPU на случайных входах"), а не на
// побитное совпадение -- побитная регрессия CPU-пути уже проверена отдельно
// (снимки regression/regression_fp32 до и после перевода Solid-SPH на
// include/physics/*, без участия GPU).
//
// SPH_DIM=3 (CMakeLists.txt): целевой первый путь Solid на GPU -- 3D
// (решение Р-Г7), поэтому здесь дополнительно упражняются ветки dim==3
// (скорости деформаций/поворота, импульс+энергия) -- самые объёмные из
// SolidSphPairInput/Output.
//
// SKIP_RETURN_CODE 77 (см. CMakeLists.txt) -- как test_cuda_toolchain: если
// cudaGetDeviceCount не находит рабочий GPU, тест пропущен, а не провален.

#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "precision.hpp"
#include "param.hpp"
#include "physics/hd.hpp"
#include "physics/sph_kernel.hpp"
#include "physics/solid_sph_pair.hpp"

namespace {

constexpr int kSkipReturnCode = 77;

}  // namespace

#define SPH_CUDA_CHECK(expr)                                                 \
    do {                                                                     \
        cudaError_t sph_cuda_err_ = (expr);                                  \
        if (sph_cuda_err_ != cudaSuccess) {                                  \
            std::fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, \
                         cudaGetErrorString(sph_cuda_err_));                 \
            std::exit(1);                                                    \
        }                                                                    \
    } while (0)

namespace {

// Допуск сверки GPU/CPU для одной SPH_HD-функции (план §5): "≤ 4 ULP без
// FMA; ≤ 1e-14 отн. с FMA" в fp64 -- берём с запасом, как в
// test_cuda_toolchain::kMathRelTol.
#if SPH_PRECISION == 64
constexpr real_t kRelTol = real_t(1e-12);
#else
constexpr real_t kRelTol = real_t(1e-5);
#endif

// Результат обоих вызовов -- W/∇W и полная пара di/dj -- в одном буфере,
// чтобы обойтись одним запуском ядра и одной парой cudaMemcpy.
struct ProbeResult
{
    real_t w;
    real_t dwdx[dim];
    physics::SolidSphPairOutput pair_out;
};

// Синтетические, но физически осмысленные (r < 2h, ненулевые скорости и
// напряжения) входные данные -- одинаковые на host и device, чтобы сверка
// проверяла именно арифметику SPH_HD-функций, а не различие входов.
SPH_HD inline void fill_sample_pair(physics::SolidSphPairInput& in)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        in.x_i[d] = real_t(0.0);
        in.x_j[d] = real_t(0.0);
        in.vx_i[d] = real_t(0.0);
        in.vx_j[d] = real_t(0.0);
    }
    in.x_j[0] = real_t(0.02);
    if constexpr (dim >= 2) in.x_j[1] = real_t(0.01);
    if constexpr (dim >= 3) in.x_j[2] = real_t(-0.005);

    in.vx_i[0] = real_t(-1.5);
    in.vx_j[0] = real_t(2.0);
    if constexpr (dim >= 2) { in.vx_i[1] = real_t(0.5); in.vx_j[1] = real_t(-0.3); }
    if constexpr (dim >= 3) { in.vx_i[2] = real_t(0.1); in.vx_j[2] = real_t(-0.2); }

    in.mass_i = real_t(0.001); in.mass_j = real_t(0.0012);
    in.rho_i = real_t(2700.0); in.rho_j = real_t(2650.0);
    in.p_i = real_t(1.0e6);    in.p_j = real_t(0.9e6);
    in.c_i = real_t(5200.0);   in.c_j = real_t(5100.0);
    in.hsml_i = real_t(0.03);  in.hsml_j = real_t(0.031);

    in.sxx_i = real_t(1.0e5); in.syy_i = real_t(-0.5e5); in.szz_i = real_t(0.2e5);
    in.sxy_i = real_t(0.3e5); in.syz_i = real_t(0.1e5);  in.sxz_i = real_t(-0.15e5);
    in.sxx_j = real_t(0.8e5); in.syy_j = real_t(-0.4e5); in.szz_j = real_t(0.15e5);
    in.sxy_j = real_t(0.25e5); in.syz_j = real_t(0.05e5); in.sxz_j = real_t(-0.1e5);

    real_t r2 = 0;
    real_t dx[dim];
    for (std::size_t d = 0; d < dim; ++d) { dx[d] = in.x_i[d] - in.x_j[d]; r2 += dx[d] * dx[d]; }
    const real_t r = std::sqrt(r2);
    const real_t mh = real_t(0.5) * (in.hsml_i + in.hsml_j);

    // W, ∇W -- через ту же SPH_HD-функцию, что и build_solid_sph_search на
    // CPU (kernel.cpp делегирует туда же, см. Г2): один вызов покрывает обе
    // физические функции, добавленные в include/physics/ этой частью Г2.
    physics::sph_kernel_w_dwdx(/*skf=*/1, r, dx, mh, in.w, in.dwdx);

    in.density_method = 1;
    in.pa_sph = 1;
    in.xi = real_t(0.0);
    in.epsilon_h = real_t(0.01);
    in.alpha = real_t(1.0);
    in.beta = real_t(1.0);
    in.etq = real_t(0.1);
}

__global__ void physics_probe_kernel(ProbeResult* out)
{
    if (threadIdx.x != 0 || blockIdx.x != 0) return;

    physics::SolidSphPairInput in{};
    fill_sample_pair(in);

    out->w = in.w;
    for (std::size_t d = 0; d < dim; ++d) out->dwdx[d] = in.dwdx[d];

    physics::SolidSphPairOutput pair_out{};
    physics::compute_solid_sph_pair(in, pair_out);
    out->pair_out = pair_out;
}

real_t max_rel_diff(real_t a, real_t b)
{
    const real_t denom = std::fabs(a) > real_t(1e-30) ? std::fabs(a) : real_t(1);
    return std::fabs(a - b) / denom;
}

}  // namespace

int main()
{
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0)
    {
        std::printf("CUDA-устройство не найдено (%s) -- тест пропущен.\n", cudaGetErrorString(err));
        return kSkipReturnCode;
    }

    cudaDeviceProp prop{};
    SPH_CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::printf("Устройство: %s, sm_%d%d, SPH_PRECISION=%d, dim=%zu\n", prop.name, prop.major, prop.minor,
                SPH_PRECISION, dim);

    // Эталон -- те же SPH_HD-функции, вызванные прямо на хосте (без
    // __CUDACC__ они -- обычные inline-функции, компилируются Clang/cl).
    physics::SolidSphPairInput host_in{};
    fill_sample_pair(host_in);
    physics::SolidSphPairOutput host_out{};
    physics::compute_solid_sph_pair(host_in, host_out);

    ProbeResult* d_out = nullptr;
    SPH_CUDA_CHECK(cudaMalloc(&d_out, sizeof(ProbeResult)));
    physics_probe_kernel<<<1, 1>>>(d_out);
    SPH_CUDA_CHECK(cudaGetLastError());
    SPH_CUDA_CHECK(cudaDeviceSynchronize());

    ProbeResult gpu_out{};
    SPH_CUDA_CHECK(cudaMemcpy(&gpu_out, d_out, sizeof(ProbeResult), cudaMemcpyDeviceToHost));
    SPH_CUDA_CHECK(cudaFree(d_out));

    bool ok = true;
    real_t max_rel = real_t(0);

    auto check = [&](const char* name, real_t host_v, real_t gpu_v) {
        const real_t rel = max_rel_diff(host_v, gpu_v);
        if (rel > max_rel) max_rel = rel;
        if (!(rel <= kRelTol))
        {
            std::printf("  %s: host=%.9e gpu=%.9e rel=%.3e (допуск %.3e)\n", name, static_cast<double>(host_v),
                        static_cast<double>(gpu_v), static_cast<double>(rel), static_cast<double>(kRelTol));
            ok = false;
        }
    };

    check("w", host_in.w, gpu_out.w);
    for (std::size_t d = 0; d < dim; ++d) check("dwdx[d]", host_in.dwdx[d], gpu_out.dwdx[d]);

    check("drho_i", host_out.drho_i, gpu_out.pair_out.drho_i);
    check("drho_j", host_out.drho_j, gpu_out.pair_out.drho_j);
    check("du_i", host_out.du_i, gpu_out.pair_out.du_i);
    check("du_j", host_out.du_j, gpu_out.pair_out.du_j);
    for (std::size_t d = 0; d < dim; ++d)
    {
        check("dv_i[d]", host_out.dv_i[d], gpu_out.pair_out.dv_i[d]);
        check("dv_j[d]", host_out.dv_j[d], gpu_out.pair_out.dv_j[d]);
        check("xsph_i[d]", host_out.xsph_i[d], gpu_out.pair_out.xsph_i[d]);
        check("xsph_j[d]", host_out.xsph_j[d], gpu_out.pair_out.xsph_j[d]);
    }
    check("dexx_i", host_out.dexx_i, gpu_out.pair_out.dexx_i);
    check("deyy_i", host_out.deyy_i, gpu_out.pair_out.deyy_i);
    check("dezz_i", host_out.dezz_i, gpu_out.pair_out.dezz_i);
    check("dexy_i", host_out.dexy_i, gpu_out.pair_out.dexy_i);
    check("dexz_i", host_out.dexz_i, gpu_out.pair_out.dexz_i);
    check("deyz_i", host_out.deyz_i, gpu_out.pair_out.deyz_i);
    check("drxy_i", host_out.drxy_i, gpu_out.pair_out.drxy_i);
    check("drxz_i", host_out.drxz_i, gpu_out.pair_out.drxz_i);
    check("dryz_i", host_out.dryz_i, gpu_out.pair_out.dryz_i);

    std::printf("  макс. отн. расхождение host/GPU: %.3e (допуск %.3e)\n", static_cast<double>(max_rel),
                static_cast<double>(kRelTol));

    if (!ok)
    {
        std::fprintf(stderr, "physics_compile_probe: ПРОВАЛ\n");
        return 1;
    }
    std::printf("physics_compile_probe: OK\n");
    return 0;
}
