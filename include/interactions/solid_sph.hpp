#pragma once
// interactions/solid_sph.hpp
// Solid-SPH: новая схема RHS для сплошной среды (fiber/SPH.md, разделы
// "Solid -- Solid" / "Solid-SPH", "Общие уравнения динамики" -> "Solid";
// fiber/PLAN.md, этап Ф4, часть A).
//
// ЭТО ЧАСТЬ A этапа Ф4: новый путь строится РЯДОМ со старым (src/density.cpp,
// src/int_force.cpp, src/art_visc.cpp, src/av_vel.cpp, src/virt_part.cpp,
// src/eos.cpp) и переиспользует их МАТЕМАТИКУ построчно, не меняя старые
// файлы и не переключая продовые вызовы (это будет часть Б, после ревью).
//
// Архитектурный контракт (SPH.md, "Основные требования к interaction kernel"):
//   - Search и primitive creation — только readonly доступ к состоянию;
//   - Compute interaction kernel не пишет напрямую в общий RHS, а формирует
//     typed Contribution<SolidDelta> через ContributionWriter;
//   - Reduce (rhs::reduce_contributions) суммирует contributions по частицам;
//   - Finalize — отдельный parallel_for по частицам, вне кернела.
//
// Что перестаёт писаться в sys (PLAN.md §3.6):
//   - p, c (EOS) -> SolidDerived, стадия "Solid derived";
//   - exx..eyz, rxy..ryz -> накопители SolidDelta, затем finalize;
//   - виртуальные (зеркальные) частицы -> GhostSet + SolidExtended, а не
//     дозапись в sys за ntotal (virt_part.cpp).
//
// Мёртвый код старого пути, СОЗНАТЕЛЬНО НЕ перенесённый (см. итоговый отчёт
// части A):
//   - физическая вязкость (int_force.cpp, `if (visc)`) -- visc == false
//     constexpr (param.hpp), ветка недостижима при текущей конфигурации;
//   - отталкивание "реальная <-> виртуальная" в ext_force.cpp
//     (is_real_virtual_pair, itype[j] < 0) -- виртуальные частицы всегда
//     наследуют ПОЛОЖИТЕЛЬНЫЙ itype (particle_system.hpp::
//     inzializ_2d_symm_index не меняет знак itype), поэтому это условие
//     тождественно ложно на всём текущем пути построения зеркальных частиц.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "eos.hpp"
#include "kernel.hpp"
#include "materials.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "physics/solid_sph_pair.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/spatial_index.hpp"
#include "rhs/types.hpp"
#include "simulation_config.hpp"

namespace interactions
{

// Метка interaction_type для ContributionMeta/детерминированного reduce
// (rhs/contributions.hpp, "Детерминизм") -- уникальна среди типов,
// вносящих вклад в Solid. kSolidContactInteractionType -- в
// interactions/solid_contact.hpp (второй interaction type того же reduce).
constexpr std::uint32_t kSolidSphInteractionType = 10;

// ----------------------------------------------------------------------
// SolidDelta -- то, что суммируется по primitives обоих interaction types,
// формирующих вклад в Solid (fiber/PLAN.md §3.3):
//     dv[dim]      ускорение (ФИНАЛЬНОЕ -- не требует постобработки:
//                  int_force momentum + art_visc momentum + contact
//                  складываются напрямую, как в single_step.cpp:
//                  "dvx[d][i] = indvxdt[d][i] + exdvxdt[d][i] + ardvxdt[d][i]");
//     drho         СЫРАЯ сумма по парам (density.cpp, con_density body).
//                  Осевой член (axis_sim && dim==2) добавляется в finalize
//                  (density.cpp: "drhodt[i] -= rho_i*v_r_i/r_i", по частице);
//     du           СЫРАЯ (нескорректированная) сумма энергии от давления
//                  (int_force.cpp body2 "he") И от искусственной вязкости
//                  (art_visc.cpp), сложенные в ОДИН накопитель. В finalize
//                  умножается на 0.5 (см. int_force.cpp "dedt[i] = tdsdt[i]
//                  + 0.5*dedt[i] + ..." и art_visc.cpp "dedt[i] = 0.5*dedt[i]" --
//                  раздельные 0.5 над раздельными суммами равны 0.5 над их
//                  суммой) и к ней прибавляется девиаторная работа s:eps/rho;
//     eps[6] (dexx..deyz), rot[3] (drxy..dryz) -- накопители скорости
//                  деформаций/поворота (int_force.cpp body1), СЫРЫЕ: при
//                  dim==2 компонента dezz заменяется в finalize на v_r/r
//                  (окружная, не приходит из парных сумм -- см. int_force.cpp,
//                  комментарий у "Окружная скорость деформации");
//     xsph[dim]    СЫРАЯ сумма XSPH-поправки (av_vel.cpp); финальное значение
//                  требует умножения на cfg.epsilon -- часть A этот
//                  постпроцессинг не выполняет (поле не входит в сравниваемый
//                  набор производных, av_vel не производная времени);
//     dh           зарезервировано под sle == 2 (PLAN.md §3.3); в части A не
//                  заполняется -- BC2/BC6 используют sle == 1, а sle == 1
//                  физически не производит производной (h пересчитывается
//                  напрямую из массы/плотности, переезжает в Postprocessing
//                  этапа Ф7, PLAN.md §3.6). Поле оставлено для интерфейсной
//                  полноты по SPH.md/PLAN.md §3.3.
// ----------------------------------------------------------------------
struct SolidDelta
{
    real_t dv[dim] = {};
    real_t drho = 0;
    real_t du = 0;
    real_t dexx = 0, deyy = 0, dezz = 0, dexy = 0, dexz = 0, deyz = 0;
    real_t drxy = 0, drxz = 0, dryz = 0;
    real_t xsph[dim] = {};
    real_t dh = 0;

    SolidDelta& operator+=(const SolidDelta& o)
    {
        for (std::size_t d = 0; d < dim; ++d)
        {
            dv[d] += o.dv[d];
            xsph[d] += o.xsph[d];
        }
        drho += o.drho;
        du += o.du;
        dexx += o.dexx; deyy += o.deyy; dezz += o.dezz;
        dexy += o.dexy; dexz += o.dexz; deyz += o.deyz;
        drxy += o.drxy; drxz += o.drxz; dryz += o.dryz;
        dh += o.dh;
        return *this;
    }
};

// Итог finalize -- поля, что сравниваются test_rhs_equivalence со старым
// путём (drho, du, dv[dim], dsxx..dsxz), ПЛЮС "сырые" скорости деформаций/
// поворота (exx..eyz, rxy..ryz), которые test_rhs_equivalence не сравнивает,
// но которые нужны fiber/PLAN.md, этап Ф4-B1: временный шим совместимости
// пишет их обратно в sys для prochnost_cem_JH2 (src/ppochnost.cpp читает
// sys.exx.. и sys.rxy.. напрямую).
//
// ОБНОВЛЕНИЕ (этап Ф7, interactions/postprocessing.hpp): перенос
// prochnost_cem_JH2/prochnost_mises в Postprocessing выполнен, но координатор
// Ф7 СОЗНАТЕЛЬНО ОСТАВИЛ этот шим -- сигнатуру prochnost_cem_JH2 не меняли,
// она по-прежнему читает sys.exx.. / sys.rxy.. напрямую (см. комментарий у
// single_step() в src/single_step.cpp и у apply_material_strength в
// interactions/postprocessing.hpp). Шим НЕ временный до Ф7 -- он остаётся
// как минимум до этапа Ф8 (унификация ODESolver, fiber/PLAN.md §3.9), где
// его наконец стоит убрать вместе с переводом prochnost_cem_JH2 на явные
// параметры. НЕ удалять этот комментарий и запись в sys без одновременного
// обновления prochnost_cem_JH2 -- иначе Ф8 рискует молча оставить мёртвую
// запись в sys или, наоборот, сломать JH-2, читающий уже неактуальные поля.
struct SolidFinal
{
    real_t drho = 0;
    real_t du = 0;
    real_t dv[dim] = {};
    real_t dsxx = 0, dsyy = 0, dszz = 0, dsxy = 0, dsyz = 0, dsxz = 0;
    real_t exx = 0, eyy = 0, ezz = 0, exy = 0, exz = 0, eyz = 0;
    real_t rxy = 0, rxz = 0, ryz = 0;
};

// ----------------------------------------------------------------------
// Solid derived: EOS (p, c) -- fiber/PLAN.md §3.6, "стадия Solid derived:
// EOS в RhsWorkspace". Реализация -- ТОЧНАЯ копия диспетчера urs()
// (src/int_force.cpp), но пишет в переданные p/c по ссылке, а не в sys.p/
// sys.c: старое urs() мутирует ParticleSystem и поэтому не годится для
// readonly-стадии RHS.
// ----------------------------------------------------------------------
// Известен ли itype диспетчеру solid_eos -- вынесено отдельно, чтобы
// compute_solid_derived могла проверить это ДО входа в параллельный участок
// (см. там же: бросок исключения из тела OpenMP `parallel for` -- UB).
inline bool solid_itype_known(int itype)
{
    switch (itype)
    {
        case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
            return true;
        default:
            return false;
    }
}

inline void solid_eos(const Materials& mats, const ParticleSystem& sys, std::size_t i, real_t& p, real_t& c)
{
    switch (sys.itype[i])
    {
        case 1: // идеальный газ
            p_gas(sys.rho[i], sys.u[i], p, c);
            break;
        case 2: // вода
            p_art_water(sys.rho[i], p, c);
            break;
        case 3: case 4: // тантал, алюминий
            eos_mi_grunaizen(mats, static_cast<std::size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.u[i], p, c);
            break;
        case 5: case 6: case 7: case 8: case 9: // LiF, полиуретаны, PMMA
            eos_ua_teta(mats, static_cast<std::size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], p, c);
            break;
        case 10: case 11: // B4C, Al2O3 (JH-2)
        {
            // urs_cem_JH2 принимает D по неконстантной ссылке, но не изменяет
            // её (см. src/eos.cpp -- существующая особенность сигнатуры,
            // не наша). Локальная копия достаточна и не мутирует sys.
            real_t D_local = sys.D[i];
            urs_cem_JH2(mats, static_cast<std::size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.dp[i], c, p, D_local);
            break;
        }
        default:
            throw std::runtime_error("solid_eos: неизвестный itype=" + std::to_string(sys.itype[i]));
    }
}

struct SolidDerived
{
    std::vector<real_t> p, c;
};

inline void compute_solid_derived(const rhs::CpuExecutor& exec, const Materials& mats,
                                   const ParticleSystem& sys, std::size_t nreal, SolidDerived& out)
{
    out.p.assign(nreal, real_t(0.0));
    out.c.assign(nreal, real_t(0.0));
    // Защитная проверка -- БУКВАЛЬНО как в int_force.cpp (тот же диапазон:
    // ghost'ы -- точные копии rho/mass реальных источников, поэтому проверка
    // только реальных частиц достаточна). Нулевая/отрицательная/NaN плотность
    // или масса здесь немедленно дают NaN в EOS и во всех производных дальше
    // по конвейеру -- лучше упасть с понятной диагностикой на этом шаге.
    // Намеренно ПОСЛЕДОВАТЕЛЬНЫЙ цикл (не exec.parallel_for), как и в
    // int_force.cpp: бросок исключения из тела OpenMP `parallel for` --
    // неопределённое поведение (исключение обязано быть поймано в том же
    // потоке, где брошено); валидация поэтому отделена от параллельного
    // вычисления EOS ниже.
    for (std::size_t i = 0; i < nreal; ++i)
    {
        if (!std::isfinite(sys.rho[i]) || sys.rho[i] <= 0 ||
            !std::isfinite(sys.mass[i]) || sys.mass[i] <= 0)
        {
            throw std::runtime_error(
                "compute_solid_derived: некорректная частица " + std::to_string(i) +
                " (rho=" + std::to_string(sys.rho[i]) +
                ", mass=" + std::to_string(sys.mass[i]) +
                ", itype=" + std::to_string(sys.itype[i]) + ")");
        }
        if (!solid_itype_known(sys.itype[i]))
        {
            // Тот же самый случай, что default: в solid_eos -- проверен здесь
            // заранее и последовательно по той же причине (UB-бросок из
            // параллельного региона), а не потому что solid_eos избыточен:
            // switch внутри solid_eos остаётся защитой на случай прямого
            // вызова в обход этой проверки.
            throw std::runtime_error(
                "compute_solid_derived: неизвестный itype=" + std::to_string(sys.itype[i]) +
                " у частицы " + std::to_string(i));
        }
    }
    exec.parallel_for(nreal, [&](std::size_t i) { solid_eos(mats, sys, i, out.p[i], out.c[i]); });
}

// ----------------------------------------------------------------------
// GhostSet -- виртуальные (зеркальные) частицы у оси при axis_sim && dim==2
// (fiber/PLAN.md §3.6: "GhostSet в RhsWorkspace; вкладов к виртуальным
// частицам нет"). Условие отбора -- БУКВАЛЬНО как в src/virt_part.cpp
// ("sys.x[0][i] < 2 * sys.hsml[i]"). Хранит только индекс исходной реальной
// частицы: координаты/скорости/масса/напряжения зеркала выводятся из него в
// build_solid_extended (см. там же обоснование p_ghost == p_source).
// ----------------------------------------------------------------------
struct GhostSet
{
    std::vector<std::size_t> source_index;
};

inline void build_ghost_set(const ParticleSystem& sys, std::size_t nreal, GhostSet& out)
{
    out.source_index.clear();
    if constexpr (axis_sim && dim == 2)
    {
        for (std::size_t i = 0; i < nreal; ++i)
        {
            if (sys.x[0][i] < real_t(2.0) * sys.hsml[i]) out.source_index.push_back(i);
        }
    }
}

// ----------------------------------------------------------------------
// SolidExtended -- readonly "расширенное" представление (реальные + ghost),
// по которому строится поиск и по которому кернел читает состояние пары.
// Заполняется ОДИН раз в начале RHS из sys (реальные) и из ghost-правила
// (зеркальные), дальше НЕ мутируется -- отличие от старого virt_part.cpp,
// который дописывал зеркальные частицы прямо в sys.
//
// Почему ext.p[ghost] = derived.p[source], БЕЗ отдельного вызова EOS на
// зеркальной частице: в старом коде virt_part запускается ДО поиска/
// density/int_force и копирует (rho, u, itype, rho_0, dp, D) зеркала из
// ТЕКУЩЕГО (входного для этого RHS) состояния источника; следующий за этим
// urs() внутри int_force затем вызывается для ВСЕГО диапазона [0, ntotal+
// nvirt), в т.ч. для зеркала, с ТЕМИ ЖЕ аргументами, что и для источника.
// EOS -- чистая функция своих аргументов, поэтому p/c зеркала совпадают с
// p/c источника побитно, и отдельный вызов solid_eos на ghost избыточен.
// ----------------------------------------------------------------------
struct SolidExtended
{
    std::vector<real_t> x[dim];
    std::vector<real_t> vx[dim];
    std::vector<real_t> mass, rho, p, c;
    std::vector<real_t> sxx, syy, szz, sxy, syz, sxz;
    std::vector<real_t> hsml;
    std::size_t nreal = 0;
    std::size_t next = 0; // nreal + nghost
};

inline void build_solid_extended(const ParticleSystem& sys, std::size_t nreal, const SolidDerived& derived,
                                  const GhostSet& ghosts, SolidExtended& ext)
{
    const std::size_t nghost = ghosts.source_index.size();
    ext.nreal = nreal;
    ext.next = nreal + nghost;

    for (std::size_t d = 0; d < dim; ++d)
    {
        ext.x[d].resize(ext.next);
        ext.vx[d].resize(ext.next);
    }
    ext.mass.resize(ext.next); ext.rho.resize(ext.next);
    ext.p.resize(ext.next);    ext.c.resize(ext.next);
    ext.sxx.resize(ext.next);  ext.syy.resize(ext.next);  ext.szz.resize(ext.next);
    ext.sxy.resize(ext.next);  ext.syz.resize(ext.next);  ext.sxz.resize(ext.next);
    ext.hsml.resize(ext.next);

    for (std::size_t i = 0; i < nreal; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d) { ext.x[d][i] = sys.x[d][i]; ext.vx[d][i] = sys.vx[d][i]; }
        ext.mass[i] = sys.mass[i]; ext.rho[i] = sys.rho[i];
        ext.p[i] = derived.p[i];   ext.c[i] = derived.c[i];
        ext.sxx[i] = sys.sxx[i];   ext.syy[i] = sys.syy[i];   ext.szz[i] = sys.szz[i];
        ext.sxy[i] = sys.sxy[i];   ext.syz[i] = sys.syz[i];   ext.sxz[i] = sys.sxz[i];
        ext.hsml[i] = sys.hsml[i];
    }
    for (std::size_t g = 0; g < nghost; ++g)
    {
        const std::size_t dst = nreal + g;
        const std::size_t src = ghosts.source_index[g];
        // Зеркалирование -- ровно как ParticleSystem::inzializ_2d_symm_index:
        // компонента 0 (r) меняет знак у координаты и скорости, остальные --
        // без изменений; из тензора напряжений знак меняет только sxy.
        for (std::size_t d = 0; d < dim; ++d)
        {
            if (d == 0) { ext.x[d][dst] = -sys.x[d][src]; ext.vx[d][dst] = -sys.vx[d][src]; }
            else        { ext.x[d][dst] =  sys.x[d][src]; ext.vx[d][dst] =  sys.vx[d][src]; }
        }
        ext.mass[dst] = sys.mass[src]; ext.rho[dst] = sys.rho[src];
        ext.p[dst] = derived.p[src];   ext.c[dst] = derived.c[src];
        ext.sxx[dst] = sys.sxx[src];   ext.syy[dst] = sys.syy[src];   ext.szz[dst] = sys.szz[src];
        ext.sxy[dst] = -sys.sxy[src];  ext.syz[dst] = sys.syz[src];   ext.sxz[dst] = sys.sxz[src];
        ext.hsml[dst] = sys.hsml[src];
    }
}

// ----------------------------------------------------------------------
// SolidSphPrimitive = (i, j, W_ij, grad_i W_ij) -- fiber/SPH.md, "Solid-SPH".
// i -- всегда реальная частица (i < ext.nreal); j -- реальная ИЛИ ghost.
// Метаданные (W, grad W) вычисляются при построении primitive из readonly
// состояния -- SPH.md: "вычисляются... и хранятся в нём как metadata".
// ----------------------------------------------------------------------
struct SolidSphPrimitive
{
    rhs::index_t i = 0, j = 0;
    real_t w = 0;
    real_t dwdx[dim] = {};
};

struct SolidSphSearchResult
{
    std::vector<SolidSphPrimitive> primitives;                          // (real, real|ghost)
    std::vector<std::pair<rhs::index_t, rhs::index_t>> real_pairs;      // подмножество (real, real)
};

// Поиск пар -- та же логика, что был у прежнего прямого поиска соседей
// (direct_find, удалён на fiber/PLAN.md, этапе Ф4-B2): канонический порядок
// i<j, критерий r < scale_k*(h_i+h_j)/2 (scale_k = 2 для кубического сплайна,
// 3 для гаусса/квинтического). SpatialIndex даёт БРОУЗ-фазу с фиксированным
// радиусом (rhs/spatial_index.hpp поддерживает только его), точный per-pair
// критерий проверяется в accept() -- воспроизводит прежний поиск для ЛЮБОГО
// (не обязательно однородного) распределения hsml.
// ghost-ghost пары отбрасываются в accept() (i >= nreal): вклад в Solid они
// дать не могут (обе стороны -- ghost), а старый путь их тоже фактически
// не учитывает (обе стороны отбрасываются downstream, см. single_step.cpp).
inline void build_solid_sph_search(const rhs::CpuExecutor& exec, const SimulationConfig& cfg,
                                    const SolidExtended& ext, SolidSphSearchResult& out)
{
    out.primitives.clear();
    out.real_pairs.clear();
    if (ext.next == 0) return;

    real_t max_h = 0;
    for (std::size_t k = 0; k < ext.next; ++k) max_h = std::max(max_h, ext.hsml[k]);
    const real_t scale_k = (cfg.skf == 1) ? real_t(2.0) : real_t(3.0);
    const real_t radius = scale_k * max_h;

    rhs::SpatialIndex<dim> index;
    std::array<const real_t*, dim> coords{};
    for (std::size_t d = 0; d < dim; ++d) coords[d] = ext.x[d].data();
    index.build(coords, ext.next, radius > 0 ? radius : real_t(1.0));

    const std::size_t nreal = ext.nreal;
    auto accept = [&](rhs::index_t i, rhs::index_t j) -> bool {
        if (i >= nreal) return false;
        real_t r2 = 0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const real_t dd = ext.x[d][i] - ext.x[d][j];
            r2 += dd * dd;
        }
        const real_t mh = real_t(0.5) * (ext.hsml[i] + ext.hsml[j]);
        return std::sqrt(r2) < scale_k * mh;
    };

    const auto pairs = index.find_pairs_parallel(exec, radius, accept);

    out.primitives.reserve(pairs.size());
    for (const auto& pr : pairs)
    {
        const rhs::index_t i = pr.first, j = pr.second;
        real_t dx[dim];
        real_t r2 = 0;
        for (std::size_t d = 0; d < dim; ++d) { dx[d] = ext.x[d][i] - ext.x[d][j]; r2 += dx[d] * dx[d]; }
        const real_t r = std::sqrt(r2);
        const real_t mh = real_t(0.5) * (ext.hsml[i] + ext.hsml[j]);

        SolidSphPrimitive prim;
        prim.i = i; prim.j = j;
        real_t w = 0;
        real_t dwdx_buf[dim];
        kernel(cfg.skf, r, std::span<const real_t>(dx, dim), mh, w, dwdx_buf);
        prim.w = w;
        for (std::size_t d = 0; d < dim; ++d) prim.dwdx[d] = dwdx_buf[d];
        out.primitives.push_back(prim);
        if (j < nreal) out.real_pairs.emplace_back(i, j);
    }
}

// ----------------------------------------------------------------------
// Compute interaction kernel -- по одному primitive формирует contributions
// к i и (если j реальна) к j. Математика -- построчная копия:
//   - неразрывность: src/density.cpp, con_density (все три density_method);
//   - скорости деформаций/поворота: src/int_force.cpp, body1 (dim 1..3);
//   - импульс+энергия (pa_sph 1/2) и осесимметричные члены:
//     src/int_force.cpp, body2;
//   - искусственная вязкость: src/art_visc.cpp;
//   - XSPH: src/av_vel.cpp.
// Физическая вязкость (visc) НЕ перенесена -- см. комментарий вверху файла.
// ----------------------------------------------------------------------
inline void compute_solid_sph_contributions(
    const rhs::CpuExecutor& exec, const SimulationConfig& cfg, const SolidExtended& ext,
    const std::vector<SolidSphPrimitive>& primitives,
    std::vector<rhs::ContributionBuffer<SolidDelta>>& thread_buffers)
{
    const std::size_t nworkers = static_cast<std::size_t>(std::max(1, exec.num_workers()));
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const std::size_t n = primitives.size();
    if (n == 0) return;

    const int density_method = cfg.density_method;
    const int pa_sph = cfg.pa_sph;
    // Параметры конфига -- double (SimulationConfig), арифметика -- real_t.
    const real_t xi = static_cast<real_t>(cfg.xi);
    const real_t epsilon_h = static_cast<real_t>(cfg.epsilon_h);
    const real_t alpha = static_cast<real_t>(cfg.alpha), beta = static_cast<real_t>(cfg.beta),
                 etq = static_cast<real_t>(cfg.etq);

    exec.parallel_for(n, [&](std::size_t k) {
        const SolidSphPrimitive& prim = primitives[k];
        const rhs::index_t i = prim.i, j = prim.j;
        const real_t w = prim.w;
        const real_t* dwdx = prim.dwdx;
        const bool j_real = (j < ext.nreal);

        SolidDelta di{}, dj{};

        // Математика пары (неразрывность, скорости деформаций/поворота,
        // импульс+энергия, искусственная вязкость, XSPH) вынесена в
        // physics::compute_solid_sph_pair (fiber/CUDA_PLAN.md §3.2, этап
        // Г2) -- построчная копия прежнего тела этого лямбда-кернела,
        // порядок операций не изменён (побитная регрессия). Здесь -- только
        // сборка входа из SolidExtended/primitive и раскладка выхода в di/dj.
        physics::SolidSphPairInput pair_in{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            pair_in.x_i[d] = ext.x[d][i];   pair_in.x_j[d] = ext.x[d][j];
            pair_in.vx_i[d] = ext.vx[d][i]; pair_in.vx_j[d] = ext.vx[d][j];
            pair_in.dwdx[d] = dwdx[d];
        }
        pair_in.mass_i = ext.mass[i]; pair_in.mass_j = ext.mass[j];
        pair_in.rho_i = ext.rho[i];   pair_in.rho_j = ext.rho[j];
        pair_in.p_i = ext.p[i];       pair_in.p_j = ext.p[j];
        pair_in.c_i = ext.c[i];       pair_in.c_j = ext.c[j];
        pair_in.hsml_i = ext.hsml[i]; pair_in.hsml_j = ext.hsml[j];
        pair_in.sxx_i = ext.sxx[i]; pair_in.syy_i = ext.syy[i]; pair_in.szz_i = ext.szz[i];
        pair_in.sxy_i = ext.sxy[i]; pair_in.syz_i = ext.syz[i]; pair_in.sxz_i = ext.sxz[i];
        pair_in.sxx_j = ext.sxx[j]; pair_in.syy_j = ext.syy[j]; pair_in.szz_j = ext.szz[j];
        pair_in.sxy_j = ext.sxy[j]; pair_in.syz_j = ext.syz[j]; pair_in.sxz_j = ext.sxz[j];
        pair_in.w = w;
        pair_in.density_method = density_method;
        pair_in.pa_sph = pa_sph;
        pair_in.xi = xi; pair_in.epsilon_h = epsilon_h;
        pair_in.alpha = alpha; pair_in.beta = beta; pair_in.etq = etq;

        physics::SolidSphPairOutput pair_out;
        physics::compute_solid_sph_pair(pair_in, pair_out);

        for (std::size_t d = 0; d < dim; ++d)
        {
            di.dv[d] = pair_out.dv_i[d];       dj.dv[d] = pair_out.dv_j[d];
            di.xsph[d] = pair_out.xsph_i[d];   dj.xsph[d] = pair_out.xsph_j[d];
        }
        di.drho = pair_out.drho_i; dj.drho = pair_out.drho_j;
        di.du = pair_out.du_i;     dj.du = pair_out.du_j;
        di.dexx = pair_out.dexx_i; di.deyy = pair_out.deyy_i; di.dezz = pair_out.dezz_i;
        di.dexy = pair_out.dexy_i; di.dexz = pair_out.dexz_i; di.deyz = pair_out.deyz_i;
        dj.dexx = pair_out.dexx_j; dj.deyy = pair_out.deyy_j; dj.dezz = pair_out.dezz_j;
        dj.dexy = pair_out.dexy_j; dj.dexz = pair_out.dexz_j; dj.deyz = pair_out.deyz_j;
        di.drxy = pair_out.drxy_i; di.drxz = pair_out.drxz_i; di.dryz = pair_out.dryz_i;
        dj.drxy = pair_out.drxy_j; dj.drxz = pair_out.drxz_j; dj.dryz = pair_out.dryz_j;

        int w_id = rhs::CpuExecutor::worker_id();
        if (w_id < 0 || static_cast<std::size_t>(w_id) >= thread_buffers.size()) w_id = 0;
        rhs::ContributionWriter<SolidDelta> writer(thread_buffers[static_cast<std::size_t>(w_id)], kSolidSphInteractionType);
        writer.emit(i, di, static_cast<std::uint64_t>(k), 0);
        if (j_real) writer.emit(j, dj, static_cast<std::uint64_t>(k), 1);
    });
}

// ----------------------------------------------------------------------
// Finalize -- parallel_for по РЕАЛЬНЫМ частицам (fiber/PLAN.md §3.3):
// осевые члены (dim==2), закон Гука с производной Яуманна (ИСПРАВЛЕННЫЙ
// знак у syz*ryz, как в текущем src/int_force.cpp), девиаторная работа
// s:eps/rho в du, домножение накопленной энергии на 0.5.
// ----------------------------------------------------------------------
inline void finalize_solid(const rhs::CpuExecutor& exec, const Materials& mats,
                            const ParticleSystem& sys, std::size_t nreal,
                            const std::vector<SolidDelta>& reduced, std::vector<SolidFinal>& out)
{
    out.assign(nreal, SolidFinal{});

    // Защита от деления на ноль -- БУКВАЛЬНО как в density.cpp, con_density
    // (осевая поправка к уравнению непрерывности): клампы, удерживавшие
    // осевые частицы от r == 0, закомментированы в time_integration.cpp,
    // поэтому r == 0 достижимо. Без этой проверки ezz/drho молча получили бы
    // Inf/NaN, которые разошлись бы по dszz/dsxy/... и энергии значительно
    // позже и далеко от причины -- падаем в точке возникновения, как старый
    // путь. Намеренно ПОСЛЕДОВАТЕЛЬНЫЙ цикл (не exec.parallel_for), как и в
    // compute_solid_derived: бросок исключения из тела OpenMP `parallel for`
    // -- неопределённое поведение, поэтому валидация отделена от
    // параллельного вычисления ниже.
    if constexpr (axis_sim && dim == 2)
    {
        for (std::size_t i = 0; i < nreal; ++i)
        {
            if (sys.x[0][i] == 0)
            {
                throw std::runtime_error(
                    "finalize_solid: частица " + std::to_string(i) + " с r = "
                    + std::to_string(sys.x[0][i]) + " == 0 (itype "
                    + std::to_string(sys.itype[i])
                    + "): осевой член неразрывности не определён");
            }
        }
    }

    exec.parallel_for(nreal, [&](std::size_t i) {
        SolidFinal& f = out[i];
        const SolidDelta& d = reduced[i];

        real_t exx = d.dexx, eyy = d.deyy, ezz = d.dezz;
        real_t exy = d.dexy, exz = d.dexz, eyz = d.deyz;
        const real_t rxy = d.drxy, rxz = d.drxz, ryz = d.dryz;

        real_t drho = d.drho;
        if constexpr (axis_sim && dim == 2)
        {
            // Окружная скорость деформации -- ПОЧАСТИЧНАЯ формула, не парная
            // сумма (int_force.cpp, "Окружная скорость деформации").
            ezz = sys.vx[0][i] / sys.x[0][i];
            // Осевой член неразрывности (density.cpp, "Осесимметричная
            // поправка к уравнению непрерывности").
            drho -= sys.rho[i] * sys.vx[0][i] / sys.x[0][i];
        }
        f.drho = drho;

        for (std::size_t dd = 0; dd < dim; ++dd) f.dv[dd] = d.dv[dd];

        // Шим совместимости (fiber/PLAN.md, этап Ф4-B1) -- ровно те значения,
        // что старый int_force.cpp пишет в sys.exx.. и sys.rxy.. (те же
        // локальные exx/eyy/ezz(с осевым переопределением)/exy/exz/eyz/rxy/
        // rxz/ryz, что участвуют в формулах dsxx.. ниже).
        f.exx = exx; f.eyy = eyy; f.ezz = ezz;
        f.exy = exy; f.exz = exz; f.eyz = eyz;
        f.rxy = rxy; f.rxz = rxz; f.ryz = ryz;

        // Энергия: 0.5*(сырая сумма давление+вязкость) + девиаторная работа.
        const real_t dedt_half = real_t(0.5) * d.du;
        real_t deviatoric;
        if constexpr (axis_sim && dim == 2)
        {
            deviatoric = (sys.sxx[i] * exx + sys.syy[i] * eyy + sys.szz[i] * ezz + real_t(2.0) * sys.sxy[i] * exy) / sys.rho[i];
        }
        else
        {
            deviatoric = (sys.sxx[i] * exx + sys.syy[i] * eyy + sys.szz[i] * ezz
                         + real_t(2.0) * sys.sxy[i] * exy + real_t(2.0) * sys.syz[i] * eyz + real_t(2.0) * sys.sxz[i] * exz) / sys.rho[i];
        }
        f.du = dedt_half + deviatoric;

        // Закон Гука + производная Яуманна (int_force.cpp, финальный блок).
        if (sys.itype[i] > 2)
        {
            const real_t G_ = shear_modulus_by_itype(mats, sys.itype[i]);
            if constexpr (dim > 1)
            {
                f.dsxx = real_t(2.0) * G_ * (exx + drho / (real_t(3.0) * sys.rho[i])) + real_t(2.0) * (sys.sxy[i] * rxy + sys.sxz[i] * rxz);
                f.dsyy = real_t(2.0) * G_ * (eyy + drho / (real_t(3.0) * sys.rho[i])) + real_t(2.0) * (sys.syz[i] * ryz - sys.sxy[i] * rxy);
                // Исправленный знак у syz*ryz (см. tests/test_kinematics_3d.cpp,
                // проверки C2-C4: прежний знак «− syz·ryz» делал след
                // поворотной части dS/dt равным 4·syz·ryz, и при вращении
                // вокруг x девиатор переставал быть бесследовым).
                f.dszz = real_t(2.0) * G_ * (ezz + drho / (real_t(3.0) * sys.rho[i])) - real_t(2.0) * (sys.sxz[i] * rxz + sys.syz[i] * ryz);
                f.dsxy = real_t(2.0) * G_ * exy + rxy * (sys.syy[i] - sys.sxx[i]) + sys.sxz[i] * ryz + sys.syz[i] * rxz;
                f.dsyz = real_t(2.0) * G_ * eyz - sys.sxy[i] * rxz - sys.sxz[i] * rxy + ryz * (sys.szz[i] - sys.syy[i]);
                f.dsxz = real_t(2.0) * G_ * exz + rxz * (sys.szz[i] - sys.sxx[i]) - sys.sxy[i] * ryz + sys.syz[i] * rxy;
            }
            else
            {
                f.dsxx = real_t(2.0) * G_ * (exx + drho / (real_t(3.0) * sys.rho[i]));
                f.dsyy = real_t(2.0) * G_ * (drho / (real_t(3.0) * sys.rho[i]));
                f.dszz = real_t(2.0) * G_ * (drho / (real_t(3.0) * sys.rho[i]));
            }
        }
    });
}

} // namespace interactions
