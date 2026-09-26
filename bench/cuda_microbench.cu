// bench/cuda_microbench.cu
// Микробенчмарк тулчейна (fiber/CUDA_PLAN.md, этап Г0): пропускная
// способность памяти и пиковая арифметика FP32/FP64 на этой карте, чтобы
// заранее знать потолок ускорения каждого режима точности (Р-Г6: рабочий
// режим GPU -- fp32, fp64 -- только эталон корректности). НЕ тест ctest:
// отдельная цель, запускается вручную, числа идут в ADR-009 и в план §1.1.
//
// Запуск (Windows, нужно окружение MSVC для сборки -- см. tools/cuda_env.bat):
//     tools\cuda_env.bat cmake --build build-cuda --target cuda_microbench
//     build-cuda\cuda_microbench.exe

#include <cuda_runtime.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

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

__global__ void copy_kernel(const float* in, float* out, size_t n) {
    size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) out[i] = in[i];
}

// Цепочка FMA без зависимости от входа компилятора наперёд (a и b меняются
// местами каждую итерацию) -- не даёт оптимизатору свернуть цикл в константу
// и не создаёт цепочки более длинной задержки, чем одна FMA.
template <typename T>
__global__ void fma_kernel(T* data, int iters) {
    size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    T a = data[i];
    T b = static_cast<T>(1.0000001);
    const T c = static_cast<T>(0.0000001);
#pragma unroll 4
    for (int k = 0; k < iters; ++k) {
        a = a * b + c;
        b = b * a + c;
    }
    data[i] = a + b;
}

double measure_bandwidth() {
    const size_t n = size_t(256) * 1024 * 1024;  // 256M float = 1 ГиБ на буфер
    float *d_in = nullptr, *d_out = nullptr;
    SPH_CUDA_CHECK(cudaMalloc(&d_in, n * sizeof(float)));
    SPH_CUDA_CHECK(cudaMalloc(&d_out, n * sizeof(float)));
    SPH_CUDA_CHECK(cudaMemset(d_in, 1, n * sizeof(float)));

    const int threads = 256;
    const int blocks = static_cast<int>((n + threads - 1) / threads);

    copy_kernel<<<blocks, threads>>>(d_in, d_out, n);  // прогрев
    SPH_CUDA_CHECK(cudaDeviceSynchronize());

    const int reps = 20;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < reps; ++r) {
        copy_kernel<<<blocks, threads>>>(d_in, d_out, n);
    }
    SPH_CUDA_CHECK(cudaDeviceSynchronize());
    auto t1 = std::chrono::high_resolution_clock::now();

    double seconds = std::chrono::duration<double>(t1 - t0).count() / reps;
    double bytes = 2.0 * static_cast<double>(n) * sizeof(float);  // чтение + запись
    double gbps = bytes / seconds / 1e9;

    SPH_CUDA_CHECK(cudaFree(d_in));
    SPH_CUDA_CHECK(cudaFree(d_out));
    return gbps;
}

template <typename T>
double measure_flops(const char* label) {
    const size_t n = size_t(4) * 1024 * 1024;
    std::vector<T> host(n, static_cast<T>(1.0));
    T* d = nullptr;
    SPH_CUDA_CHECK(cudaMalloc(&d, n * sizeof(T)));
    SPH_CUDA_CHECK(cudaMemcpy(d, host.data(), n * sizeof(T), cudaMemcpyHostToDevice));

    const int threads = 256;
    const int blocks = static_cast<int>((n + threads - 1) / threads);
    const int iters = 4000;

    fma_kernel<T><<<blocks, threads>>>(d, iters);  // прогрев
    SPH_CUDA_CHECK(cudaDeviceSynchronize());

    auto t0 = std::chrono::high_resolution_clock::now();
    fma_kernel<T><<<blocks, threads>>>(d, iters);
    SPH_CUDA_CHECK(cudaDeviceSynchronize());
    auto t1 = std::chrono::high_resolution_clock::now();
    double seconds = std::chrono::duration<double>(t1 - t0).count();

    // 2 FMA на итерацию тела цикла, 2 flop на FMA.
    double flops = 2.0 * 2.0 * static_cast<double>(iters) * static_cast<double>(n);
    double gflops = flops / seconds / 1e9;
    std::printf("  %s: %.1f GFLOP/s (n=%zu, iters=%d, %.4f с)\n", label, gflops, n, iters, seconds);

    SPH_CUDA_CHECK(cudaFree(d));
    return gflops;
}

}  // namespace

int main() {
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count == 0) {
        std::printf("CUDA-устройство не найдено -- бенчмарк пропущен.\n");
        return 0;
    }

    cudaDeviceProp prop{};
    SPH_CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::printf("Устройство: %s, sm_%d%d, %.1f ГБ\n", prop.name, prop.major, prop.minor,
                static_cast<double>(prop.totalGlobalMem) / (1024.0 * 1024.0 * 1024.0));

    double gbps = measure_bandwidth();
    std::printf("Пропускная способность памяти (copy, чтение+запись): %.1f ГБ/с\n", gbps);

    measure_flops<float>("FP32 FMA");
    measure_flops<double>("FP64 FMA");

    return 0;
}
