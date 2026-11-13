#pragma once
// interactions/clamp.hpp
// Кинематическое закрепление (Clamp) частиц Solid и бусин Fiber: условие
// Дирихле x_{i,d}(t) = x_{i,d}(0), v_{i,d}(t) = 0 для отмеченных компонент d.
// Реализуется проекцией ускорений (у Solid -- и XSPH-поправки) и нулевой
// обратной массой в стадии трения; advance() не меняется: при a = 0, av = 0 и
// v(0) = 0 и предиктор-корректор, и Эйлер дают v = 0 и x += dt*0 побитово.
// Постановка и обоснование: fiber/SPH.md, раздел «Clamp (кинематическое
// закрепление)».
//
// МАСКА. Для частицы Solid i in [0,ntotal) и бусины k in [0,nfiber) --
// std::uint8_t, бит d == закреплена компонента d; 0 == свободна;
// kClampAll == закреплена целиком. Биты >= dim запрещены.
//
// ТРИ ТОЧКИ ВСТАВКИ В ШАГ (src/time_integration.cpp), все под !clamps.empty():
//   [1] после RHS (Solid-SPH, контакты, Fiber-Bond/Fiber-Fiber/Fiber-Solid
//       полностью сложены в dvx/dv), до dt: apply_clamp_projection для Solid
//       (dvx и av) и для Fiber (dv); накопление реакции R = -sum m*a;
//   [2] сразу после compute_friction_counts, до compute_friction_contributions:
//       apply_clamp_friction_counts -- c_k = c_s = 0 у полностью закреплённых
//       (бесконечная масса в проекции трения);
//   [3] после сложения вкладов трения, до advance: повторная проекция dvx и dv
//       БЕЗ av; реакция прибавляется к накопителю [1].
//
// ПАРАЛЛЕЛИЗМ И ДЕТЕРМИНИЗМ (ТЗ §6). Это стадия constraint над уже сведённым
// RHS, не interaction kernel. Все функции вызываются с главного потока между
// стадиями, параллелизм -- только через rhs::CpuExecutor (без #pragma omp,
// mutex и atomics). Элемент k списка пишет только в ячейки своей частицы
// (index строго возрастает => уникален) и в свой слот item_force[d][k], поэтому
// гонок нет. Реакция -- НЕ parallel_reduce (результат зависит от числа
// потоков), а двухуровневая сумма с фиксированным разбиением на блоки по
// kClampReduceBlock элементов (clamp_fixed_block_sum): внутри блока сумма по
// возрастанию k, затем по возрастанию b на главном потоке; результат побитово
// одинаков при любом числе потоков. Короткие списки (< kClampSerialThreshold)
// обходятся последовательно на главном потоке -- результат тот же. Память
// выделяется один раз до цикла (prepare_clamp_workspace); в шаге ни
// аллокаций, ни исключений, ни ввода-вывода.
//
// GPU-задел: clamp_project_item -- операция над одним элементом, inline, без
// executor, std::vector и исключений (кандидат на __host__ __device__).
// ClampList -- SoA плоских массивов.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"

namespace interactions
{

inline constexpr std::uint8_t kClampAll = static_cast<std::uint8_t>((1u << dim) - 1u);

// Размер блока детерминированной суммы (не зависит от числа потоков).
inline constexpr std::size_t kClampReduceBlock = 256;
// Списки короче этого порога обходятся последовательно на главном потоке.
// Замер (8 потоков, g++ -O3): открытие OpenMP-региона ~30 мкс, проекция
// последовательно ~4.7 нс/частицу; последовательный путь быстрее до ~6-7 тыс.
// частиц (4096: 19 против 33 мкс; 8192: 43 против 33 мкс).
inline constexpr std::size_t kClampSerialThreshold = 8192;

// Входные маски (плотные). Пустой вектор == зажимов этого типа нет.
struct ClampSet
{
    std::vector<std::uint8_t> solid_mask; // размер ntotal или 0
    std::vector<std::uint8_t> fiber_mask; // размер nfiber или 0
    bool empty() const { return solid_mask.empty() && fiber_mask.empty(); }
};

// Компактный список для цикла (SoA). Инвариант: index строго возрастает.
struct ClampList
{
    std::vector<rhs::index_t> index;  // зажатые частицы
    std::vector<std::uint8_t> mask;   // их маски (!= 0)
    std::vector<rhs::index_t> full;   // подсписок с mask == kClampAll, по возрастанию
    std::size_t size() const { return index.size(); }
};

// Реакция опоры за последний шаг (точки [1] + [3]); пишет ядро, читает сцена/хук.
struct ClampReport
{
    double solid_force[dim] = {};
    double fiber_force[dim] = {};
    std::size_t solid_clamped = 0, fiber_clamped = 0;
    double solid_max_h = 0.0; // max hsml среди зажатых Solid (печать при старте)
    // Диагностика шага (решение владельца 2026-11-11): по полной реакции R_k = [1] + [3] каждой
    // зажатой частицы -- сумма модулей sum_k |R_k| и max_k |R_k|, Н. Векторная сумма solid_force /
    // fiber_force гасит вклады противоположных концов нити, сумма модулей -- нет. Пишется ядром
    // каждый шаг (как force), в чекпоинт НЕ входит: формат чекпоинта прежний.
    double solid_force_abs = 0.0, fiber_force_abs = 0.0;
    double solid_force_max = 0.0, fiber_force_max = 0.0;
};

// Рабочая память стадии: списки + буфер для детерминированной суммы.
// Выделяется до цикла по времени; в шаге память только переиспользуется.
struct ClampWorkspace
{
    ClampList solid, fiber;
    std::vector<double> item_force[dim];    // размер >= max(solid.size(), fiber.size())
    std::vector<double> block_partial[dim]; // размер >= ceil(n / kClampReduceBlock)
    // Полная реакция шага по элементам списка ([1] пишет, [3] прибавляет): размер solid.size() /
    // fiber.size(); item_abs -- |R_k| для суммы модулей (размер >= max из двух).
    std::vector<double> acc_solid[dim], acc_fiber[dim];
    std::vector<double> item_abs;
};

// ----------------------------------------------------------------------------
// Операция над ОДНИМ элементом списка (GPU-friendly).
//   mask    -- маска частицы;
//   m       -- масса;
//   dv_comp -- dv_comp[d] указывает на ячейку dv[d][i] (обязательно, dim указателей);
//   av_comp -- nullptr (XSPH не трогать) или массив из dim указателей на av[d][i];
//   item_force[d] -- выход: -m*a_d для закреплённых d до обнуления, 0 для свободных.
// ----------------------------------------------------------------------------
inline void clamp_project_item(std::uint8_t mask, double m, double* const* dv_comp,
                               double* const* av_comp, double (&item_force)[dim])
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        if ((mask >> d) & 1u)
        {
            item_force[d] = -m * (*dv_comp[d]);
            *dv_comp[d] = 0.0;
            if (av_comp != nullptr && av_comp[d] != nullptr) *av_comp[d] = 0.0;
        }
        else
        {
            item_force[d] = 0.0;
        }
    }
}

// Детерминированная двухуровневая сумма items[0..n): блоки по kClampReduceBlock,
// сумма внутри блока по возрастанию k, затем по возрастанию b на вызывающем
// потоке. block_partial растёт только если мал (в шаге заранее размерен).
inline double clamp_fixed_block_sum(const rhs::CpuExecutor& exec, const std::vector<double>& items,
                                    std::size_t n, std::vector<double>& block_partial)
{
    if (n == 0) return 0.0;
    const std::size_t nblocks = (n + kClampReduceBlock - 1) / kClampReduceBlock;
    if (block_partial.size() < nblocks) block_partial.resize(nblocks);
    const double* src = items.data();
    double* part = block_partial.data();
    auto block_fn = [src, part, n](std::size_t b)
    {
        const std::size_t k0 = b * kClampReduceBlock;
        const std::size_t k1 = std::min(n, k0 + kClampReduceBlock);
        double s = 0.0;
        for (std::size_t k = k0; k < k1; ++k) s += src[k];
        part[b] = s;
    };
    if (n < kClampSerialThreshold)
        for (std::size_t b = 0; b < nblocks; ++b) block_fn(b);
    else
        exec.parallel_for(nblocks, 1, block_fn);
    double total = 0.0;
    for (std::size_t b = 0; b < nblocks; ++b) total += part[b];
    return total;
}

// Построение списка из плотной маски (один проход по возрастанию индекса).
inline void build_clamp_list(const std::vector<std::uint8_t>& mask, std::size_t n, ClampList& out)
{
    out.index.clear();
    out.mask.clear();
    out.full.clear();
    if (mask.empty()) return;
    if (mask.size() != n)
        throw std::runtime_error("clamp: размер маски (" + std::to_string(mask.size()) +
                                 ") не равен числу частиц (" + std::to_string(n) + ")");
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::uint8_t mk = mask[i];
        if (mk == 0) continue;
        if ((static_cast<unsigned>(mk) >> dim) != 0u)
            throw std::runtime_error("clamp: у частицы " + std::to_string(i) +
                                     " заданы биты маски >= dim (" + std::to_string(dim) + ")");
        out.index.push_back(i);
        out.mask.push_back(mk);
        if (mk == kClampAll) out.full.push_back(i);
    }
}

// Проверки при старте. Бросает std::runtime_error с индексом частицы.
inline void validate_clamps(const ClampSet& clamps, const ParticleSystem& sys, std::size_t ntotal,
                            const fiber::FiberState& fb)
{
    if (!clamps.solid_mask.empty() && clamps.solid_mask.size() != ntotal)
        throw std::runtime_error("clamp: размер solid_mask (" + std::to_string(clamps.solid_mask.size()) +
                                 ") не равен ntotal (" + std::to_string(ntotal) + ")");
    if (!clamps.fiber_mask.empty() && clamps.fiber_mask.size() != fb.size())
        throw std::runtime_error("clamp: размер fiber_mask (" + std::to_string(clamps.fiber_mask.size()) +
                                 ") не равен числу бусин (" + std::to_string(fb.size()) + ")");

    std::size_t n_solid = 0, n_solid_full = 0, n_bead = 0, n_bead_full = 0;
    double max_h = 0.0;
    for (std::size_t i = 0; i < clamps.solid_mask.size(); ++i)
    {
        const std::uint8_t mk = clamps.solid_mask[i];
        if (mk == 0) continue;
        if ((static_cast<unsigned>(mk) >> dim) != 0u)
            throw std::runtime_error("clamp: у частицы Solid " + std::to_string(i) +
                                     " заданы биты маски >= dim");
        for (std::size_t d = 0; d < dim; ++d)
            if (((mk >> d) & 1u) && sys.vx[d][i] != 0.0)
                throw std::runtime_error("clamp: у закреплённой частицы Solid " + std::to_string(i) +
                                         " компонента " + std::to_string(d) +
                                         " имеет ненулевую начальную скорость " +
                                         std::to_string(sys.vx[d][i]));
        ++n_solid;
        if (mk == kClampAll) ++n_solid_full;
        max_h = std::max(max_h, sys.hsml[i]);
    }
    for (std::size_t k = 0; k < clamps.fiber_mask.size(); ++k)
    {
        const std::uint8_t mk = clamps.fiber_mask[k];
        if (mk == 0) continue;
        if ((static_cast<unsigned>(mk) >> dim) != 0u)
            throw std::runtime_error("clamp: у бусины " + std::to_string(k) +
                                     " заданы биты маски >= dim");
        for (std::size_t d = 0; d < dim; ++d)
            if (((mk >> d) & 1u) && fb.vx[d][k] != 0.0)
                throw std::runtime_error("clamp: у закреплённой бусины " + std::to_string(k) +
                                         " компонента " + std::to_string(d) +
                                         " имеет ненулевую начальную скорость " +
                                         std::to_string(fb.vx[d][k]));
        ++n_bead;
        if (mk == kClampAll) ++n_bead_full;
    }
    if (!clamps.empty())
    {
        std::wcout << L"Clamp: закреплено частиц Solid " << n_solid << L" (полностью " << n_solid_full
                   << L"), бусин " << n_bead << L" (полностью " << n_bead_full
                   << L"), solid_max_h = " << max_h << std::endl;
    }
}

// Проекция P + реакция. DvArr: std::vector<std::vector<double>> (Solid dvx)
// или std::vector<double>[dim] (Fiber dv) -- достаточно operator[][].
// av == nullptr -- XSPH не трогать. force_out += -sum m a (добавляется).
// ws.item_force[d] должен быть размерен заранее (prepare_clamp_workspace).
// acc != nullptr -- полная реакция по элементам: acc[d][k] = f (acc_add = false, точка [1]) или
// acc[d][k] += f (acc_add = true, точка [3]); размер acc[d] >= list.size() (prepare_clamp_workspace).
template <class DvArr, class AvArr>
void apply_clamp_projection(const rhs::CpuExecutor& exec, const ClampList& list,
                            const std::vector<double>& mass, DvArr& dv, AvArr* av,
                            ClampWorkspace& ws, double (&force_out)[dim],
                            std::vector<double> (*acc)[dim] = nullptr, bool acc_add = false)
{
    const std::size_t n = list.size();
    if (n == 0) return;
    const rhs::index_t* idx = list.index.data();
    const std::uint8_t* msk = list.mask.data();
    double* items[dim];
    for (std::size_t d = 0; d < dim; ++d) items[d] = ws.item_force[d].data();
    double* accp[dim] = {};
    if (acc != nullptr)
        for (std::size_t d = 0; d < dim; ++d) accp[d] = (*acc)[d].data();
    const bool use_acc = acc != nullptr;

    auto item_fn = [&, idx, msk](std::size_t k)
    {
        const rhs::index_t i = idx[k];
        double* dvp[dim];
        double* avp[dim];
        for (std::size_t d = 0; d < dim; ++d)
        {
            dvp[d] = &dv[d][i];
            avp[d] = (av != nullptr) ? &(*av)[d][i] : nullptr;
        }
        double f[dim];
        clamp_project_item(msk[k], mass[i], dvp, (av != nullptr) ? avp : nullptr, f);
        for (std::size_t d = 0; d < dim; ++d) items[d][k] = f[d];
        if (use_acc)
            for (std::size_t d = 0; d < dim; ++d) accp[d][k] = acc_add ? accp[d][k] + f[d] : f[d];
    };

    if (n < kClampSerialThreshold)
        for (std::size_t k = 0; k < n; ++k) item_fn(k);
    else
        exec.parallel_for(n, item_fn);

    for (std::size_t d = 0; d < dim; ++d)
        force_out[d] += clamp_fixed_block_sum(exec, ws.item_force[d], n, ws.block_partial[d]);
}

// Сумма модулей и максимум полной реакции по элементам списка (acc после [1] и [3]).
// Сумма -- той же фиксированной двухуровневой схемой, что и force (побитово одинакова при любом
// числе потоков); максимум от порядка не зависит.
inline void clamp_abs_summary(const rhs::CpuExecutor& exec, const std::vector<double> (&acc)[dim],
                              std::size_t n, ClampWorkspace& ws, double& sum_abs, double& max_abs)
{
    sum_abs = 0.0;
    max_abs = 0.0;
    if (n == 0) return;
    if (ws.item_abs.size() < n) ws.item_abs.resize(n, 0.0);
    double* out = ws.item_abs.data();
    const double* a[dim];
    for (std::size_t d = 0; d < dim; ++d) a[d] = acc[d].data();
    auto fn = [out, &a](std::size_t k)
    {
        double s = 0.0;
        for (std::size_t d = 0; d < dim; ++d) s += a[d][k] * a[d][k];
        out[k] = std::sqrt(s);
    };
    if (n < kClampSerialThreshold)
        for (std::size_t k = 0; k < n; ++k) fn(k);
    else
        exec.parallel_for(n, fn);
    sum_abs = clamp_fixed_block_sum(exec, ws.item_abs, n, ws.block_partial[0]);
    for (std::size_t k = 0; k < n; ++k) max_abs = std::max(max_abs, out[k]);
}

// c_k = 0, c_s = 0 для полностью закреплённых (списки full).
inline void apply_clamp_friction_counts(const rhs::CpuExecutor& exec, const ClampWorkspace& ws,
                                        std::vector<double>& bead_count,
                                        std::vector<double>& solid_count)
{
    auto zero_list = [&exec](const std::vector<rhs::index_t>& full, std::vector<double>& count)
    {
        const std::size_t n = full.size();
        if (n == 0) return;
        const rhs::index_t* idx = full.data();
        double* c = count.data();
        if (n < kClampSerialThreshold)
            for (std::size_t k = 0; k < n; ++k) c[idx[k]] = 0.0;
        else
            exec.parallel_for(n, [idx, c](std::size_t k) { c[idx[k]] = 0.0; });
    };
    zero_list(ws.fiber.full, bead_count);
    zero_list(ws.solid.full, solid_count);
}

// Разметка области: out_mask[i] |= mask для всех i с lo <= x_i <= hi (замкнутый бокс).
inline void mark_clamp_box(const rhs::CpuExecutor& exec, const std::vector<double> (&x)[dim],
                           std::size_t n, const double (&lo)[dim], const double (&hi)[dim],
                           std::uint8_t mask, std::vector<std::uint8_t>& out_mask)
{
    if (out_mask.size() < n) out_mask.resize(n, 0);
    std::uint8_t* om = out_mask.data();
    exec.parallel_for(n, [&x, &lo, &hi, mask, om](std::size_t i)
    {
        bool inside = true;
        for (std::size_t d = 0; d < dim; ++d)
            if (x[d][i] < lo[d] || x[d][i] > hi[d]) { inside = false; break; }
        if (inside) om[i] = static_cast<std::uint8_t>(om[i] | mask);
    });
}

// Подготовка рабочей памяти до цикла по времени: списки из масок и размеры
// буферов (в шаге аллокаций нет).
inline void prepare_clamp_workspace(const ClampSet& clamps, std::size_t ntotal, std::size_t nfiber,
                                    ClampWorkspace& ws)
{
    build_clamp_list(clamps.solid_mask, ntotal, ws.solid);
    build_clamp_list(clamps.fiber_mask, nfiber, ws.fiber);
    const std::size_t nmax = std::max(ws.solid.size(), ws.fiber.size());
    const std::size_t nblocks = (nmax + kClampReduceBlock - 1) / kClampReduceBlock;
    for (std::size_t d = 0; d < dim; ++d)
    {
        if (ws.item_force[d].size() < nmax) ws.item_force[d].resize(nmax, 0.0);
        if (ws.block_partial[d].size() < nblocks) ws.block_partial[d].resize(nblocks, 0.0);
        if (ws.acc_solid[d].size() < ws.solid.size()) ws.acc_solid[d].resize(ws.solid.size(), 0.0);
        if (ws.acc_fiber[d].size() < ws.fiber.size()) ws.acc_fiber[d].resize(ws.fiber.size(), 0.0);
    }
    if (ws.item_abs.size() < nmax) ws.item_abs.resize(nmax, 0.0);
}

// Заполнение счётчиков и solid_max_h отчёта (реакцию пишет шаг).
inline void clamp_report_counts(const ClampWorkspace& ws, const ParticleSystem& sys,
                                ClampReport& report)
{
    report.solid_clamped = ws.solid.size();
    report.fiber_clamped = ws.fiber.size();
    double h = 0.0;
    for (std::size_t k = 0; k < ws.solid.size(); ++k) h = std::max(h, sys.hsml[ws.solid.index[k]]);
    report.solid_max_h = h;
}

} // namespace interactions
