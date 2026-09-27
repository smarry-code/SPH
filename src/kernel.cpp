// kernel.cpp
// Определение функции вычисления ядра сглаживания и его производной.
// Объявление — в kernel.hpp.
//
// Математика перенесена в include/physics/sph_kernel.hpp как SPH_HD-функцию
// (fiber/CUDA_PLAN.md §3.2, этап Г2: "функции ядра сглаживания (W, ∇W)");
// здесь остаётся только тонкая обёртка под старую сигнатуру (std::span) —
// её используют существующие вызывающие стороны (interactions/solid_sph.hpp)
// и, потенциально, будущие части Г2. Арифметика не изменена ни на бит:
// sph_kernel_w_dwdx — дословная копия тела этой функции.
#include "kernel.hpp"
#include "physics/sph_kernel.hpp"

void kernel(int skf, real_t r, std::span<const real_t> dx, real_t mhsml, real_t& w, real_t* dwdx)
{
    physics::sph_kernel_w_dwdx(skf, r, dx.data(), mhsml, w, dwdx);
}
