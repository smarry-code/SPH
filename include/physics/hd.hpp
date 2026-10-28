#pragma once
// physics/hd.hpp
// Инфраструктура SPH_HD (fiber/CUDA_PLAN.md §3.2, этап Г2): макрос для
// функций, которые обязаны собираться и для CPU (host), и для GPU (device)
// без изменения кода. Под nvcc (__CUDACC__ определён только при компиляции
// .cu -- см. решение Р-Г1: .cu компилирует nvcc, весь остальной проект,
// включая эти же заголовки внутри обычных .cpp, -- Clang) разворачивается в
// __host__ __device__; при обычной CPU-сборке -- пусто, т.е. обычная inline-
// функция.
//
// Всё, что попадает в include/physics/*, -- inline-функции без выделения
// памяти, исключений, std::vector/std::string/iostream и виртуальных
// вызовов (план §3.2, "Правила кода"): это код, который должен быть
// одинаково компилируем и хостовым компилятором (Clang/cl), и device-частью
// nvcc.
#ifdef __CUDACC__
#define SPH_HD __host__ __device__
#else
#define SPH_HD
#endif

// sph_min/sph_max -- замена std::min/std::max для SPH_HD-кода. Под Windows
// windows.h определяет макросы min/max, которые ломают CUB в CUDA-целях
// (план §1.1/§3.2: CUDA-цели собираются с -DNOMINMAX); использовать
// std::min/std::max в SPH_HD-функциях запрещено правилами кода (план §3.2)
// -- вместо них везде, где нужен минимум/максимум пары real_t (или другого
// сравнимого типа), используются эти две функции.
template <typename T>
SPH_HD inline T sph_min(T a, T b) { return (a < b) ? a : b; }

template <typename T>
SPH_HD inline T sph_max(T a, T b) { return (a > b) ? a : b; }

// Код ошибки SPH_HD-функции (план §3.2 перечисляет "коды ошибок" как часть
// hd.hpp; §3.7 -- "на GPU пишут код ошибки частицы", редукция max по шагу).
// kUnknownItype добавлен этапом Г2д вместе с переводом диспетчера EOS
// (physics::solid_eos_dispatch, physics/eos.hpp) на SPH_HD: там, где старый
// диспетчер (interactions::solid_sph.hpp::solid_eos, src/int_force.cpp::
// urs()) бросал исключение на неизвестном itype в default-ветке switch,
// SPH_HD-функции исключения запрещены (план §3.2) -- вместо этого код
// ошибки, а исключение с ТЕМ ЖЕ текстом бросает CPU-обёртка, получившая
// этот код (план §3.7). Остальные проверки (rho<=0/NaN, r=0 у оси, NaN в h)
// пока остаются последовательными CPU-проверками ДО параллельного участка
// (compute_solid_derived, finalize_solid, h_upgrade) -- перечисление кодов
// расширяется по мере того, как конкретная часть Г2/Г4 переводит
// соответствующую функцию на SPH_HD.
enum class SphErrorCode : int
{
    kOk = 0,
    kUnknownItype = 1, // неизвестный itype в диспетчере EOS (physics::solid_eos_dispatch)
};
