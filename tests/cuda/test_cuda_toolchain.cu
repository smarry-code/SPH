// tests/cuda/test_cuda_toolchain.cu
// Первый CUDA-тест проекта (fiber/CUDA_PLAN.md, этап Г0). Ловит поломку
// ТУЛЧЕЙНА (обновление CUDA Toolkit, драйвера, MSVC, Clang) -- не физику:
// физики здесь ещё нет, она появится в Г2/Г4 (include/physics/*).
//
// Повторяет проверку внешнего спайка 2026-09-27 (план §1.1), теперь внутри
// сборки проекта:
//   - device-лямбда в ядре (--extended-lambda, Р-Г1);
//   - cub::DeviceRadixSort::SortPairs (со стабильностью),
//     cub::DeviceScan::ExclusiveSum, cub::DeviceReduce::Min;
//   - pow/exp/sqrt/log в real_t на GPU против std:: на CPU;
//   - "precision.hpp" и "param.hpp" (заголовок конвейера) компилируются и
//     Clang (обычные .cpp проекта), и MSVC cl (хост-часть .cu под nvcc на
//     Windows, решение Р-Г1) -- сам факт компиляции этого файла cl уже
//     проверяет требование.
//
// SKIP_RETURN_CODE 77 (см. CMakeLists.txt): если cudaGetDeviceCount не
// находит работоспособный GPU, тест возвращает 77 и CTest показывает его
// как пропущенный, а не проваленный -- на машине без CUDA-карты.
//
// Цели test_cuda_toolchain (SPH_PRECISION=64) и test_cuda_toolchain_fp32
// (SPH_PRECISION=32) собирают ровно этот файл с разным макросом (по
// образцу test_regression_fp32, CMakeLists.txt).

#include <cuda_runtime.h>
#include <cub/cub.cuh>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "precision.hpp"
#include "param.hpp"  // заголовок конвейера -- проверка компилируемости под cl

namespace {

constexpr int kSkipReturnCode = 77;

}  // namespace

// Не в анонимном namespace: макрос, используется во всём файле.
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

// Допуск сверки GPU/CPU для pow/exp/sqrt/log (план §1.1): проверка
// тулчейна на машине разработки дала 4e-16 отн. (fp64) и 2e-7 отн. (fp32);
// берём с запасом на другую карту/драйвер/версию CUDA (Р-Г1а).
#if SPH_PRECISION == 64
constexpr real_t kMathRelTol = real_t(1e-15);
#else
constexpr real_t kMathRelTol = real_t(1e-6);
#endif

constexpr int kMathN = 1 << 14;

__global__ void math_kernel(const real_t* x, real_t* out, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    // device-лямбда общего вида, вызывается прямо в ядре (--extended-lambda,
    // Р-Г1): в физике (Г2, include/physics/*) SPH_HD-функции строятся так же.
    auto combine = [] __device__(real_t a, real_t b, real_t c, real_t d) {
        return a + b + c + d;
    };
    real_t v = x[i];
    real_t p = pow(v, real_t(2.5));
    real_t e = exp(v * real_t(0.1));
    real_t s = sqrt(v + real_t(10.0));
    real_t l = log(v + real_t(10.0));
    out[i] = combine(p, e, s, l);
}

bool check_math() {
    std::vector<real_t> host_x(kMathN), host_out(kMathN), device_out(kMathN);
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist(0.01, 5.0);
    for (int i = 0; i < kMathN; ++i) {
        host_x[i] = static_cast<real_t>(dist(rng));
    }
    for (int i = 0; i < kMathN; ++i) {
        real_t v = host_x[i];
        real_t p = std::pow(v, real_t(2.5));
        real_t e = std::exp(v * real_t(0.1));
        real_t s = std::sqrt(v + real_t(10.0));
        real_t l = std::log(v + real_t(10.0));
        host_out[i] = p + e + s + l;
    }

    real_t* d_x = nullptr;
    real_t* d_out = nullptr;
    SPH_CUDA_CHECK(cudaMalloc(&d_x, kMathN * sizeof(real_t)));
    SPH_CUDA_CHECK(cudaMalloc(&d_out, kMathN * sizeof(real_t)));
    SPH_CUDA_CHECK(cudaMemcpy(d_x, host_x.data(), kMathN * sizeof(real_t), cudaMemcpyHostToDevice));

    int threads = 256;
    int blocks = (kMathN + threads - 1) / threads;
    math_kernel<<<blocks, threads>>>(d_x, d_out, kMathN);
    SPH_CUDA_CHECK(cudaGetLastError());
    SPH_CUDA_CHECK(cudaDeviceSynchronize());

    SPH_CUDA_CHECK(cudaMemcpy(device_out.data(), d_out, kMathN * sizeof(real_t), cudaMemcpyDeviceToHost));
    SPH_CUDA_CHECK(cudaFree(d_x));
    SPH_CUDA_CHECK(cudaFree(d_out));

    real_t max_rel = real_t(0);
    int worst = -1;
    for (int i = 0; i < kMathN; ++i) {
        real_t denom = std::fabs(host_out[i]) > real_t(1e-30) ? std::fabs(host_out[i]) : real_t(1);
        real_t rel = std::fabs(device_out[i] - host_out[i]) / denom;
        if (rel > max_rel) {
            max_rel = rel;
            worst = i;
        }
    }
    std::printf("  pow/exp/sqrt/log GPU vs CPU: макс. отн. расхождение %.3e (допуск %.3e), частица %d\n",
                static_cast<double>(max_rel), static_cast<double>(kMathRelTol), worst);
    return max_rel <= kMathRelTol;
}

bool check_cub_radix_sort() {
    constexpr int n = 1 << 16;
    std::vector<unsigned int> host_keys(n), host_vals(n);
    std::mt19937 rng(777);
    // Узкий диапазон ключей -> много повторов -> видна стабильность сортировки.
    std::uniform_int_distribution<unsigned int> dist(0u, 999u);
    for (int i = 0; i < n; ++i) {
        host_keys[i] = dist(rng);
        host_vals[i] = static_cast<unsigned int>(i);
    }

    unsigned int *d_keys_in = nullptr, *d_keys_out = nullptr;
    unsigned int *d_vals_in = nullptr, *d_vals_out = nullptr;
    SPH_CUDA_CHECK(cudaMalloc(&d_keys_in, n * sizeof(unsigned int)));
    SPH_CUDA_CHECK(cudaMalloc(&d_keys_out, n * sizeof(unsigned int)));
    SPH_CUDA_CHECK(cudaMalloc(&d_vals_in, n * sizeof(unsigned int)));
    SPH_CUDA_CHECK(cudaMalloc(&d_vals_out, n * sizeof(unsigned int)));
    SPH_CUDA_CHECK(cudaMemcpy(d_keys_in, host_keys.data(), n * sizeof(unsigned int), cudaMemcpyHostToDevice));
    SPH_CUDA_CHECK(cudaMemcpy(d_vals_in, host_vals.data(), n * sizeof(unsigned int), cudaMemcpyHostToDevice));

    void* d_temp = nullptr;
    size_t temp_bytes = 0;
    cub::DeviceRadixSort::SortPairs(d_temp, temp_bytes, d_keys_in, d_keys_out, d_vals_in, d_vals_out, n);
    SPH_CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
    cub::DeviceRadixSort::SortPairs(d_temp, temp_bytes, d_keys_in, d_keys_out, d_vals_in, d_vals_out, n);
    SPH_CUDA_CHECK(cudaGetLastError());
    SPH_CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<unsigned int> gpu_keys(n), gpu_vals(n);
    SPH_CUDA_CHECK(cudaMemcpy(gpu_keys.data(), d_keys_out, n * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    SPH_CUDA_CHECK(cudaMemcpy(gpu_vals.data(), d_vals_out, n * sizeof(unsigned int), cudaMemcpyDeviceToHost));

    SPH_CUDA_CHECK(cudaFree(d_temp));
    SPH_CUDA_CHECK(cudaFree(d_keys_in));
    SPH_CUDA_CHECK(cudaFree(d_keys_out));
    SPH_CUDA_CHECK(cudaFree(d_vals_in));
    SPH_CUDA_CHECK(cudaFree(d_vals_out));

    // Эталон на CPU: stable_sort по ключу; значение -- исходный индекс,
    // поэтому у стабильной сортировки значения при равных ключах идут по
    // возрастанию.
    std::vector<unsigned int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = static_cast<unsigned int>(i);
    std::stable_sort(idx.begin(), idx.end(),
                      [&](unsigned int a, unsigned int b) { return host_keys[a] < host_keys[b]; });

    bool ok = true;
    for (int i = 0; i < n; ++i) {
        unsigned int expect_key = host_keys[idx[i]];
        unsigned int expect_val = host_vals[idx[i]];
        if (gpu_keys[i] != expect_key || gpu_vals[i] != expect_val) {
            std::printf("  cub::DeviceRadixSort::SortPairs: расхождение на позиции %d (key %u/%u, value %u/%u)\n",
                        i, gpu_keys[i], expect_key, gpu_vals[i], expect_val);
            ok = false;
            break;
        }
    }
    std::printf("  cub::DeviceRadixSort::SortPairs (стабильность, n=%d): %s\n", n, ok ? "OK" : "ПРОВАЛ");
    return ok;
}

bool check_cub_scan() {
    constexpr int n = 1 << 14;
    std::vector<int> host_in(n);
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(0, 9);
    for (int i = 0; i < n; ++i) host_in[i] = dist(rng);

    std::vector<long long> host_out(n);
    long long running = 0;
    for (int i = 0; i < n; ++i) {
        host_out[i] = running;
        running += host_in[i];
    }

    int* d_in = nullptr;
    long long* d_out = nullptr;
    SPH_CUDA_CHECK(cudaMalloc(&d_in, n * sizeof(int)));
    SPH_CUDA_CHECK(cudaMalloc(&d_out, n * sizeof(long long)));
    SPH_CUDA_CHECK(cudaMemcpy(d_in, host_in.data(), n * sizeof(int), cudaMemcpyHostToDevice));

    void* d_temp = nullptr;
    size_t temp_bytes = 0;
    cub::DeviceScan::ExclusiveSum(d_temp, temp_bytes, d_in, d_out, n);
    SPH_CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
    cub::DeviceScan::ExclusiveSum(d_temp, temp_bytes, d_in, d_out, n);
    SPH_CUDA_CHECK(cudaGetLastError());
    SPH_CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<long long> gpu_out(n);
    SPH_CUDA_CHECK(cudaMemcpy(gpu_out.data(), d_out, n * sizeof(long long), cudaMemcpyDeviceToHost));
    SPH_CUDA_CHECK(cudaFree(d_temp));
    SPH_CUDA_CHECK(cudaFree(d_in));
    SPH_CUDA_CHECK(cudaFree(d_out));

    bool ok = (host_out == gpu_out);
    std::printf("  cub::DeviceScan::ExclusiveSum (n=%d): %s\n", n, ok ? "OK" : "ПРОВАЛ");
    return ok;
}

bool check_cub_reduce() {
    constexpr int n = 1 << 14;
    std::vector<real_t> host_in(n);
    std::mt19937 rng(9001);
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
    for (int i = 0; i < n; ++i) host_in[i] = static_cast<real_t>(dist(rng));
    real_t host_min = *std::min_element(host_in.begin(), host_in.end());

    real_t* d_in = nullptr;
    real_t* d_out = nullptr;
    SPH_CUDA_CHECK(cudaMalloc(&d_in, n * sizeof(real_t)));
    SPH_CUDA_CHECK(cudaMalloc(&d_out, sizeof(real_t)));
    SPH_CUDA_CHECK(cudaMemcpy(d_in, host_in.data(), n * sizeof(real_t), cudaMemcpyHostToDevice));

    void* d_temp = nullptr;
    size_t temp_bytes = 0;
    cub::DeviceReduce::Min(d_temp, temp_bytes, d_in, d_out, n);
    SPH_CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
    cub::DeviceReduce::Min(d_temp, temp_bytes, d_in, d_out, n);
    SPH_CUDA_CHECK(cudaGetLastError());
    SPH_CUDA_CHECK(cudaDeviceSynchronize());

    real_t gpu_min = real_t(0);
    SPH_CUDA_CHECK(cudaMemcpy(&gpu_min, d_out, sizeof(real_t), cudaMemcpyDeviceToHost));
    SPH_CUDA_CHECK(cudaFree(d_temp));
    SPH_CUDA_CHECK(cudaFree(d_in));
    SPH_CUDA_CHECK(cudaFree(d_out));

    bool ok = (gpu_min == host_min);
    std::printf("  cub::DeviceReduce::Min (n=%d): %s (host=%.6f, gpu=%.6f)\n", n, ok ? "OK" : "ПРОВАЛ",
                static_cast<double>(host_min), static_cast<double>(gpu_min));
    return ok;
}

}  // namespace

int main() {
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        std::printf("CUDA-устройство не найдено (%s) -- тест пропущен.\n", cudaGetErrorString(err));
        return kSkipReturnCode;
    }

    cudaDeviceProp prop{};
    SPH_CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::printf("Устройство: %s, sm_%d%d, %.1f ГБ, драйвер поддерживает CUDA %d.%d, SPH_PRECISION=%d\n", prop.name,
                prop.major, prop.minor, static_cast<double>(prop.totalGlobalMem) / (1024.0 * 1024.0 * 1024.0),
                CUDART_VERSION / 1000, (CUDART_VERSION % 100) / 10, SPH_PRECISION);

    bool ok = true;
    ok &= check_math();
    ok &= check_cub_radix_sort();
    ok &= check_cub_scan();
    ok &= check_cub_reduce();

    // Повторный прогон радикс-сортировки: GPU-путь обязан быть
    // детерминирован от запуска к запуску (Р-Г5, план §5) -- дешёвая
    // проверка этого свойства уже на первом CUDA-тесте.
    ok &= check_cub_radix_sort();

    if (!ok) {
        std::fprintf(stderr, "test_cuda_toolchain: ПРОВАЛ\n");
        return 1;
    }
    std::printf("test_cuda_toolchain: OK\n");
    return 0;
}
