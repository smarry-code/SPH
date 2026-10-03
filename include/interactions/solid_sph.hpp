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
    double dv[dim] = {};
    double drho = 0.0;
    double du = 0.0;
    double dexx = 0.0, deyy = 0.0, dezz = 0.0, dexy = 0.0, dexz = 0.0, deyz = 0.0;
    double drxy = 0.0, drxz = 0.0, dryz = 0.0;
    double xsph[dim] = {};
    double dh = 0.0;

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
    double drho = 0.0;
    double du = 0.0;
    double dv[dim] = {};
    double dsxx = 0.0, dsyy = 0.0, dszz = 0.0, dsxy = 0.0, dsyz = 0.0, dsxz = 0.0;
    double exx = 0.0, eyy = 0.0, ezz = 0.0, exy = 0.0, exz = 0.0, eyz = 0.0;
    double rxy = 0.0, rxz = 0.0, ryz = 0.0;
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

inline void solid_eos(const Materials& mats, const ParticleSystem& sys, std::size_t i, double& p, double& c)
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
            double D_local = sys.D[i];
            urs_cem_JH2(mats, static_cast<std::size_t>(sys.itype[i]), sys.rho[i], sys.rho_0[i], sys.dp[i], c, p, D_local);
            break;
        }
        default:
            throw std::runtime_error("solid_eos: неизвестный itype=" + std::to_string(sys.itype[i]));
    }
}

struct SolidDerived
{
    std::vector<double> p, c;
};

inline void compute_solid_derived(const rhs::CpuExecutor& exec, const Materials& mats,
                                   const ParticleSystem& sys, std::size_t nreal, SolidDerived& out)
{
    out.p.assign(nreal, 0.0);
    out.c.assign(nreal, 0.0);
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
        if (!std::isfinite(sys.rho[i]) || sys.rho[i] <= 0.0 ||
            !std::isfinite(sys.mass[i]) || sys.mass[i] <= 0.0)
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
            if (sys.x[0][i] < 2.0 * sys.hsml[i]) out.source_index.push_back(i);
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
    std::vector<double> x[dim];
    std::vector<double> vx[dim];
    std::vector<double> mass, rho, p, c;
    std::vector<double> sxx, syy, szz, sxy, syz, sxz;
    std::vector<double> hsml;
    std::vector<int> body; // номер тела (sys.body); ghost наследует body источника
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
    ext.body.resize(ext.next);

    for (std::size_t i = 0; i < nreal; ++i)
    {
        for (std::size_t d = 0; d < dim; ++d) { ext.x[d][i] = sys.x[d][i]; ext.vx[d][i] = sys.vx[d][i]; }
        ext.mass[i] = sys.mass[i]; ext.rho[i] = sys.rho[i];
        ext.p[i] = derived.p[i];   ext.c[i] = derived.c[i];
        ext.sxx[i] = sys.sxx[i];   ext.syy[i] = sys.syy[i];   ext.szz[i] = sys.szz[i];
        ext.sxy[i] = sys.sxy[i];   ext.syz[i] = sys.syz[i];   ext.sxz[i] = sys.sxz[i];
        ext.hsml[i] = sys.hsml[i];
        ext.body[i] = sys.body[i];
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
        ext.body[dst] = sys.body[src]; // как inzializ_2d_symm_index
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
    double w = 0.0;
    double dwdx[dim] = {};
};

// Межтельная пара real-real из real_pairs (body[i] != body[j]) -- список для
// Solid-Solid-Contact. k -- исходный индекс пары в real_pairs: он передаётся
// как primitive_id в meta вкладов, чтобы meta (и детерминированный порядок
// reduce) не зависели от того, что контакт идёт по сокращённому списку.
struct SolidContactPair
{
    rhs::index_t i, j;
    std::uint64_t k;
};

// Счётчики одного блока сборки примитивов (проход 1 -> смещения для прохода 2).
// alignas(64): запись счётчиков каждым потоком -- в своей кэш-линии (false sharing).
struct alignas(64) SolidBlockCount
{
    std::size_t prim = 0, real = 0, cont = 0;
};

struct SolidSphSearchResult
{
    std::vector<SolidSphPrimitive> primitives;                          // (real, real|ghost)
    std::vector<std::pair<rhs::index_t, rhs::index_t>> real_pairs;      // подмножество (real, real)
    // Подмножество real_pairs с body[i] != body[j], порядок real_pairs сохранён.
    // По нему идут контактное ядро и диагностика контакта: межтельных пар на
    // порядки меньше, чем всех real_pairs.
    std::vector<SolidContactPair> contact_pairs;

    // Рабочая память поиска. Лежит здесь, а не в локальных переменных
    // build_solid_sph_search, чтобы при постоянном SolidSphSearchResult
    // (SolidStepWorkspace) сетка и буферы пар переиспользовались между
    // шагами, а не выделялись заново на каждом (glibc mmap/munmap).
    rhs::SpatialIndex<dim> index;
    std::vector<std::pair<rhs::index_t, rhs::index_t>> pairs;
    std::vector<SolidBlockCount> blocks;
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
    // Изоляция SPH-сумм по телам (ТЗ эксперименты Guo 2019, §4.2): без
    // Solid-Solid-Contact разные тела не взаимодействовали бы вовсе.
    if (cfg.sph_same_body_only && cfg.contact_force_type == 0)
    {
        throw std::runtime_error(
            "sph_same_body_only = 1 требует contact_force_type != 0: иначе разные тела не взаимодействуют вовсе");
    }

    // clear() в начале нет намеренно: primitives/real_pairs/contact_pairs
    // получают размер resize() ниже и целиком перезаписываются по индексам, а
    // clear()+resize() каждый шаг заново обнуляли бы десятки МБ на главном потоке.
    // Буферы живут в SolidStepWorkspace и переиспользуются; только на ранних
    // выходах массивы обнуляются явно.
    if (ext.next == 0)
    {
        out.primitives.resize(0);
        out.real_pairs.resize(0);
        out.contact_pairs.resize(0);
        return;
    }

    double max_h = 0.0;
    for (std::size_t k = 0; k < ext.next; ++k) max_h = std::max(max_h, ext.hsml[k]);
    const double scale_k = (cfg.skf == 1) ? 2.0 : 3.0;
    const double radius = scale_k * max_h;

    rhs::SpatialIndex<dim>& index = out.index;
    std::array<const double*, dim> coords{};
    for (std::size_t d = 0; d < dim; ++d) coords[d] = ext.x[d].data();
    index.build(coords, ext.next, radius > 0.0 ? radius : 1.0);

    const std::size_t nreal = ext.nreal;
    const bool same_body_only = cfg.sph_same_body_only;
    auto accept = [&](rhs::index_t i, rhs::index_t j) -> bool {
        if (i >= nreal) return false;
        double r2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double dd = ext.x[d][i] - ext.x[d][j];
            r2 += dd * dd;
        }
        const double mh = 0.5 * (ext.hsml[i] + ext.hsml[j]);
        return std::sqrt(r2) < scale_k * mh;
    };

    const auto& pairs = out.pairs;
    index.find_pairs_parallel(exec, radius, accept, out.pairs);

    const std::size_t n = pairs.size();
    if (n == 0)
    {
        out.primitives.resize(0);
        out.real_pairs.resize(0);
        out.contact_pairs.resize(0);
        return;
    }

    // Сборка по индексам: счёт -> префиксная сумма -> запись. Пары делятся на nb
    // непрерывных блоков; проход 1 считает у каждого блока число примитивов,
    // real_pairs и contact_pairs, исключающая сумма даёт смещения блоков, проход 2
    // пишет элементы по своим смещениям. Позиция элемента определяется номером
    // пары, а не блоком или потоком: массивы совпадают с прежней последовательной
    // сборкой (push_back в порядке pairs) при любом числе потоков. Блок --
    // отдельная кэш-линия SolidBlockCount.
    const std::size_t nb = std::clamp<std::size_t>(n / 2048, 1, 8 * static_cast<std::size_t>(exec.num_workers()));
    auto block_begin = [&](std::size_t t) { return n * t / nb; };
    out.blocks.resize(nb + 1);
    SolidBlockCount* const blocks = out.blocks.data();

    // Пара отсечена изоляцией по телам (ТЗ эксперименты Guo 2019, §4.2): real-real
    // разных тел не даёт SPH-primitive, но остаётся в real_pairs -- по этому же
    // списку работает Solid-Solid-Contact (solid_contact.hpp). Пары с ghost
    // (j >= nreal) не затрагиваются.
    exec.parallel_for(nb, 1, [&](std::size_t t) {
        std::size_t prim = 0, real = 0, cont = 0;
        for (std::size_t q = block_begin(t); q < block_begin(t + 1); ++q)
        {
            const rhs::index_t i = pairs[q].first, j = pairs[q].second;
            const bool cross = j < nreal && ext.body[i] != ext.body[j];
            if (!(same_body_only && cross)) ++prim;
            if (j < nreal) ++real;
            if (cross) ++cont;
        }
        blocks[t].prim = prim; blocks[t].real = real; blocks[t].cont = cont;
    });
    SolidBlockCount run;
    for (std::size_t t = 0; t < nb; ++t)
    {
        const SolidBlockCount c = blocks[t];
        blocks[t] = run;
        run.prim += c.prim; run.real += c.real; run.cont += c.cont;
    }
    blocks[nb] = run;

    out.primitives.resize(run.prim);
    out.real_pairs.resize(run.real);
    out.contact_pairs.resize(run.cont);
    SolidSphPrimitive* const primitives = out.primitives.data();
    std::pair<rhs::index_t, rhs::index_t>* const real_pairs = out.real_pairs.data();
    SolidContactPair* const contact_pairs = out.contact_pairs.data();

    // Тело прохода 2; ядро -- параметр (kernel_fn), чтобы skf выбирался один
    // раз вне цикла, а вызов ядра встраивался (у каждой лямбды свой тип -- своя
    // инстанциация build_loop; указатель на функцию дал бы непрямой вызов).
    // Межтельная пара попадает и в contact_pairs (k -- её индекс в real_pairs).
    auto build_loop = [&](auto kernel_fn) {
        exec.parallel_for(nb, 1, [&](std::size_t t) {
            std::size_t p = blocks[t].prim, r = blocks[t].real, c = blocks[t].cont;
            for (std::size_t q = block_begin(t); q < block_begin(t + 1); ++q)
            {
                const rhs::index_t i = pairs[q].first, j = pairs[q].second;
                if (same_body_only && j < nreal && ext.body[i] != ext.body[j])
                {
                    contact_pairs[c++] = {i, j, static_cast<std::uint64_t>(r)};
                    real_pairs[r++] = {i, j};
                    continue;
                }
                double dx[dim];
                double r2 = 0.0;
                for (std::size_t d = 0; d < dim; ++d) { dx[d] = ext.x[d][i] - ext.x[d][j]; r2 += dx[d] * dx[d]; }
                const double rr = std::sqrt(r2);
                const double mh = 0.5 * (ext.hsml[i] + ext.hsml[j]);

                SolidSphPrimitive prim;
                prim.i = i; prim.j = j;
                double w = 0.0;
                double dwdx_buf[dim];
                kernel_fn(rr, std::span<const double>(dx, dim), mh, w, dwdx_buf);
                prim.w = w;
                for (std::size_t d = 0; d < dim; ++d) prim.dwdx[d] = dwdx_buf[d];
                primitives[p++] = prim;
                if (j < nreal)
                {
                    if (ext.body[i] != ext.body[j])
                        contact_pairs[c++] = {i, j, static_cast<std::uint64_t>(r)};
                    real_pairs[r++] = {i, j};
                }
            }
        });
    };
    switch (cfg.skf)
    {
    case 1:
        build_loop([](double r, std::span<const double> dx, double mh, double& w, double* dwdx) {
            kernel_cubic(r, dx, mh, w, dwdx);
        });
        break;
    case 2:
        build_loop([](double r, std::span<const double> dx, double mh, double& w, double* dwdx) {
            kernel_gauss(r, dx, mh, w, dwdx);
        });
        break;
    case 3:
        build_loop([](double r, std::span<const double> dx, double mh, double& w, double* dwdx) {
            kernel_quintic(r, dx, mh, w, dwdx);
        });
        break;
    // Неизвестный skf: kernel() обнуляет w и dwdx.
    default:
        build_loop([&](double r, std::span<const double> dx, double mh, double& w, double* dwdx) {
            kernel(cfg.skf, r, dx, mh, w, dwdx);
        });
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
    rhs::IndexedContributions<SolidDelta>& out)
{
    // Число вкладов известно заранее: примитив k даёт вклад для i и, если j
    // реальна, для j. Позиции -- исключающая префиксная сумма; порядок вкладов
    // i0, j0, i1, j1, ... канонический (не зависит от числа потоков), общего
    // буфера с заголовком, в который пишут все потоки, нет.
    const std::size_t n = primitives.size();
    if (n == 0)
    {
        out.prepare(0);
        return;
    }
    out.counts.resize(n);
    out.offsets.resize(n + 1);
    std::size_t* const cnt = out.counts.data();
    std::size_t* const offset = out.offsets.data();
    exec.parallel_for(n, [&](std::size_t k) { cnt[k] = 1 + (primitives[k].j < ext.nreal ? 1 : 0); });
    exec.parallel_scan(cnt, offset, n);
    out.prepare(offset[n]);
    const auto w_out = out.writer();

    const int density_method = cfg.density_method;
    const int pa_sph = cfg.pa_sph;
    const double xi = cfg.xi;
    const double epsilon_h = cfg.epsilon_h;
    const double alpha = cfg.alpha, beta = cfg.beta, etq = cfg.etq;

    exec.parallel_for(n, [&](std::size_t k) {
        const SolidSphPrimitive& prim = primitives[k];
        const rhs::index_t i = prim.i, j = prim.j;
        const double w = prim.w;
        const double* dwdx = prim.dwdx;
        const bool j_real = (j < ext.nreal);

        SolidDelta di{}, dj{};

        // --- неразрывность (density.cpp::con_density, тело body) -----------
        {
            double dvx_ij[dim]; // v_i - v_j, как в density.cpp
            for (std::size_t d = 0; d < dim; ++d) dvx_ij[d] = ext.vx[d][i] - ext.vx[d][j];
            double vcc = 0.0;
            for (std::size_t d = 0; d < dim; ++d) vcc += dvx_ij[d] * dwdx[d];

            if (density_method == 1)
            {
                di.drho = ext.mass[j] * vcc;
                dj.drho = ext.mass[i] * vcc;
            }
            else if (density_method == 2)
            {
                di.drho = (ext.rho[i] / ext.rho[j]) * ext.mass[j] * vcc;
                dj.drho = (ext.rho[j] / ext.rho[i]) * ext.mass[i] * vcc;
            }
            else // density_method == 3 (искусственная диффузия, Bonet-Lok/Molteni-Colagrossi)
            {
                const double q_i = ext.mass[j] * vcc;
                const double q_j = ext.mass[i] * vcc;
                const double h_ij = 0.5 * (ext.hsml[i] + ext.hsml[j]);
                double rr = 0.0;
                double dx[dim];
                for (std::size_t d = 0; d < dim; ++d) { dx[d] = ext.x[d][i] - ext.x[d][j]; rr += dx[d] * dx[d]; }
                const double vi = 1.0 / ext.rho[i], vj = 1.0 / ext.rho[j];
                const double denom = rr + epsilon_h * h_ij * h_ij;
                const double coef_i = 2.0 * (vi / vj - 1.0) / denom;
                const double coef_j = 2.0 * (vj / vi - 1.0) / denom;
                double psi_i = 0.0, psi_j = 0.0;
                for (std::size_t d = 0; d < dim; ++d)
                {
                    psi_i += (coef_i * dx[d]) * dwdx[d];
                    psi_j += (coef_j * (-dx[d])) * dwdx[d];
                }
                const double dVi = ext.mass[i] * vi, dVj = ext.mass[j] * vj;
                const double c0 = 0.5 * (ext.c[i] + ext.c[j]);
                const double D_i = xi * h_ij * c0 * psi_i * dVj;
                const double D_j = xi * h_ij * c0 * psi_j * dVi;
                di.drho = q_i - ext.rho[i] * D_i;
                dj.drho = q_j - ext.rho[j] * D_j;
            }
        }

        // --- скорости деформаций/поворота (int_force.cpp body1) ------------
        {
            double dvx_ji[dim]; // v_j - v_i, как в int_force.cpp body1
            for (std::size_t d = 0; d < dim; ++d) dvx_ji[d] = ext.vx[d][j] - ext.vx[d][i];

            if constexpr (dim == 1)
            {
                di.dexx = (ext.mass[j] * dvx_ji[0] * dwdx[0]) / ext.rho[j];
                dj.dexx = (ext.mass[i] * dvx_ji[0] * dwdx[0]) / ext.rho[i];
            }
            else if constexpr (dim == 2)
            {
                di.dexx = (ext.mass[j] * dvx_ji[0] * dwdx[0]) / ext.rho[j];
                di.deyy = (ext.mass[j] * dvx_ji[1] * dwdx[1]) / ext.rho[j];
                di.dexy = 0.5 * (((ext.mass[j] * dvx_ji[0] * dwdx[1]) / ext.rho[j])
                                + ((ext.mass[j] * dvx_ji[1] * dwdx[0]) / ext.rho[j]));
                dj.dexx = (ext.mass[i] * dvx_ji[0] * dwdx[0]) / ext.rho[i];
                dj.deyy = (ext.mass[i] * dvx_ji[1] * dwdx[1]) / ext.rho[i];
                dj.dexy = 0.5 * (((ext.mass[i] * dvx_ji[0] * dwdx[1]) / ext.rho[i])
                                + ((ext.mass[i] * dvx_ji[1] * dwdx[0]) / ext.rho[i]));
                di.drxy = 0.5 * (((ext.mass[j] * dvx_ji[0] * dwdx[1]) / ext.rho[j])
                                - ((ext.mass[j] * dvx_ji[1] * dwdx[0]) / ext.rho[j]));
                dj.drxy = 0.5 * (((ext.mass[i] * dvx_ji[0] * dwdx[1]) / ext.rho[i])
                                - ((ext.mass[i] * dvx_ji[1] * dwdx[0]) / ext.rho[i]));
            }
            else if constexpr (dim == 3)
            {
                di.dexx = (ext.mass[j] * dvx_ji[0] * dwdx[0]) / ext.rho[j];
                di.deyy = (ext.mass[j] * dvx_ji[1] * dwdx[1]) / ext.rho[j];
                di.dezz = (ext.mass[j] * dvx_ji[2] * dwdx[2]) / ext.rho[j];
                di.dexy = 0.5 * (((ext.mass[j] * dvx_ji[0] * dwdx[1]) / ext.rho[j]) + ((ext.mass[j] * dvx_ji[1] * dwdx[0]) / ext.rho[j]));
                di.deyz = 0.5 * (((ext.mass[j] * dvx_ji[1] * dwdx[2]) / ext.rho[j]) + ((ext.mass[j] * dvx_ji[2] * dwdx[1]) / ext.rho[j]));
                di.dexz = 0.5 * (((ext.mass[j] * dvx_ji[0] * dwdx[2]) / ext.rho[j]) + ((ext.mass[j] * dvx_ji[2] * dwdx[0]) / ext.rho[j]));
                dj.dexx = (ext.mass[i] * dvx_ji[0] * dwdx[0]) / ext.rho[i];
                dj.deyy = (ext.mass[i] * dvx_ji[1] * dwdx[1]) / ext.rho[i];
                dj.dezz = (ext.mass[i] * dvx_ji[2] * dwdx[2]) / ext.rho[i];
                dj.dexy = 0.5 * (((ext.mass[i] * dvx_ji[0] * dwdx[1]) / ext.rho[i]) + ((ext.mass[i] * dvx_ji[1] * dwdx[0]) / ext.rho[i]));
                dj.deyz = 0.5 * (((ext.mass[i] * dvx_ji[1] * dwdx[2]) / ext.rho[i]) + ((ext.mass[i] * dvx_ji[2] * dwdx[1]) / ext.rho[i]));
                dj.dexz = 0.5 * (((ext.mass[i] * dvx_ji[0] * dwdx[2]) / ext.rho[i]) + ((ext.mass[i] * dvx_ji[2] * dwdx[0]) / ext.rho[i]));
                di.drxy = 0.5 * (((ext.mass[j] * dvx_ji[0] * dwdx[1]) / ext.rho[j]) - ((ext.mass[j] * dvx_ji[1] * dwdx[0]) / ext.rho[j]));
                di.dryz = 0.5 * (((ext.mass[j] * dvx_ji[1] * dwdx[2]) / ext.rho[j]) - ((ext.mass[j] * dvx_ji[2] * dwdx[1]) / ext.rho[j]));
                di.drxz = 0.5 * (((ext.mass[j] * dvx_ji[0] * dwdx[2]) / ext.rho[j]) - ((ext.mass[j] * dvx_ji[2] * dwdx[0]) / ext.rho[j]));
                dj.drxy = 0.5 * (((ext.mass[i] * dvx_ji[0] * dwdx[1]) / ext.rho[i]) - ((ext.mass[i] * dvx_ji[1] * dwdx[0]) / ext.rho[i]));
                dj.dryz = 0.5 * (((ext.mass[i] * dvx_ji[1] * dwdx[2]) / ext.rho[i]) - ((ext.mass[i] * dvx_ji[2] * dwdx[1]) / ext.rho[i]));
                dj.drxz = 0.5 * (((ext.mass[i] * dvx_ji[0] * dwdx[2]) / ext.rho[i]) - ((ext.mass[i] * dvx_ji[2] * dwdx[0]) / ext.rho[i]));
            }
        }

        // --- импульс + энергия (int_force.cpp body2, pa_sph 1/2) -----------
        {
            double dvx_ji[dim];
            for (std::size_t d = 0; d < dim; ++d) dvx_ji[d] = ext.vx[d][j] - ext.vx[d][i];
            const double rhoij = 1.0 / (ext.rho[i] * ext.rho[j]);
            double he = 0.0;

            if (pa_sph == 1)
            {
                for (std::size_t d = 0; d < dim; ++d)
                {
                    double h = -(ext.p[i] + ext.p[j]) * dwdx[d];
                    he += dvx_ji[d] * h; // he использует h ДО добавления члена напряжений (int_force.cpp)

                    if (d == 0)
                    {
                        h += (ext.sxx[i] + ext.sxx[j]) * dwdx[0];
                        if constexpr (dim >= 2) h += (ext.sxy[i] + ext.sxy[j]) * dwdx[1];
                        if constexpr (dim == 3) h += (ext.sxz[i] + ext.sxz[j]) * dwdx[2];
                    }
                    if constexpr (dim >= 2)
                    {
                        if (d == 1)
                        {
                            h += (ext.sxy[i] + ext.sxy[j]) * dwdx[0] + (ext.syy[i] + ext.syy[j]) * dwdx[1];
                            if constexpr (dim == 3) h += (ext.syz[i] + ext.syz[j]) * dwdx[2];
                        }
                    }
                    if constexpr (dim == 3)
                    {
                        if (d == 2)
                        {
                            h += (ext.sxz[i] + ext.sxz[j]) * dwdx[0] + (ext.syz[i] + ext.syz[j]) * dwdx[1]
                               + (ext.szz[i] + ext.szz[j]) * dwdx[2];
                        }
                    }

                    // int_force.cpp (pa_sph==1): "h *= rhoij;" применяется к
                    // ПОЛНОМУ h (давление+напряжения) ПЕРЕД domножением на
                    // массу -- в отличие от he, который использует h ДО этого
                    // домножения (см. "he += ... * h;" выше, до if(d==0)/...).
                    h *= rhoij;
                    const double dvxdt_i = ext.mass[j] * h;
                    const double dvxdt_j = ext.mass[i] * h;
                    di.dv[d] += dvxdt_i;
                    dj.dv[d] += -dvxdt_j;

                    if constexpr (axis_sim && dim == 2)
                    {
                        if (d == 0)
                        {
                            di.dv[0] += rhoij * ext.mass[j] * ((ext.sxx[j] - ext.szz[j]) / ext.x[0][j]) * w;
                            dj.dv[0] += rhoij * ext.mass[i] * ((ext.sxx[i] - ext.szz[i]) / ext.x[0][i]) * w;
                        }
                        else if (d == 1)
                        {
                            di.dv[1] += rhoij * ext.mass[j] * (ext.sxy[j] / ext.x[0][j]) * w;
                            dj.dv[1] += rhoij * ext.mass[i] * (ext.sxy[i] / ext.x[0][i]) * w;
                        }
                    }
                }
                he *= rhoij;
                di.du += ext.mass[j] * he;
                dj.du += ext.mass[i] * he;
            }
            else // pa_sph == 2
            {
                for (std::size_t d = 0; d < dim; ++d)
                {
                    double h = -(ext.p[i] / (ext.rho[i] * ext.rho[i]) + ext.p[j] / (ext.rho[j] * ext.rho[j])) * dwdx[d];
                    he += dvx_ji[d] * h;

                    if (d == 0)
                    {
                        h += (ext.sxx[i] / (ext.rho[i] * ext.rho[i]) + ext.sxx[j] / (ext.rho[j] * ext.rho[j])) * dwdx[0];
                        if constexpr (dim >= 2)
                            h += (ext.sxy[i] / (ext.rho[i] * ext.rho[i]) + ext.sxy[j] / (ext.rho[j] * ext.rho[j])) * dwdx[1];
                        if constexpr (dim == 3)
                            h += (ext.sxz[i] / (ext.rho[i] * ext.rho[i]) + ext.sxz[j] / (ext.rho[j] * ext.rho[j])) * dwdx[2];
                    }
                    if constexpr (dim >= 2)
                    {
                        if (d == 1)
                        {
                            h += (ext.sxy[i] / (ext.rho[i] * ext.rho[i]) + ext.sxy[j] / (ext.rho[j] * ext.rho[j])) * dwdx[0]
                               + (ext.syy[i] / (ext.rho[i] * ext.rho[i]) + ext.syy[j] / (ext.rho[j] * ext.rho[j])) * dwdx[1];
                            if constexpr (dim == 3)
                                h += (ext.syz[i] / (ext.rho[i] * ext.rho[i]) + ext.syz[j] / (ext.rho[j] * ext.rho[j])) * dwdx[2];
                        }
                    }
                    if constexpr (dim == 3)
                    {
                        if (d == 2)
                        {
                            h += (ext.sxz[i] / (ext.rho[i] * ext.rho[i]) + ext.sxz[j] / (ext.rho[j] * ext.rho[j])) * dwdx[0]
                               + (ext.syz[i] / (ext.rho[i] * ext.rho[i]) + ext.syz[j] / (ext.rho[j] * ext.rho[j])) * dwdx[1]
                               + (ext.szz[i] / (ext.rho[i] * ext.rho[i]) + ext.szz[j] / (ext.rho[j] * ext.rho[j])) * dwdx[2];
                        }
                    }

                    const double dvxdt_i = ext.mass[j] * h;
                    const double dvxdt_j = ext.mass[i] * h;
                    di.dv[d] += dvxdt_i;
                    dj.dv[d] += -dvxdt_j;

                    if constexpr (axis_sim && dim == 2)
                    {
                        if (d == 0)
                        {
                            di.dv[0] += rhoij * ext.mass[j] * ((ext.sxx[j] - ext.szz[j]) / ext.x[0][j]) * w;
                            dj.dv[0] += rhoij * ext.mass[i] * ((ext.sxx[i] - ext.szz[i]) / ext.x[0][i]) * w;
                        }
                        else if (d == 1)
                        {
                            di.dv[1] += rhoij * ext.mass[j] * (ext.sxy[j] / ext.x[0][j]) * w;
                            dj.dv[1] += rhoij * ext.mass[i] * (ext.sxy[i] / ext.x[0][i]) * w;
                        }
                    }
                }
                di.du += ext.mass[j] * he;
                dj.du += ext.mass[i] * he;
            }
        }

        // --- искусственная вязкость (art_visc.cpp) --------------------------
        {
            const double mhsml = 0.5 * (ext.hsml[i] + ext.hsml[j]);
            double dvx_iminusj[dim], dx[dim];
            double vr = 0.0, rr = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                dvx_iminusj[d] = ext.vx[d][i] - ext.vx[d][j];
                dx[d] = ext.x[d][i] - ext.x[d][j];
                vr += dvx_iminusj[d] * dx[d];
                rr += dx[d] * dx[d];
            }
            if (vr < 0.0)
            {
                const double muv = mhsml * vr / (rr + mhsml * mhsml * etq * etq);
                const double mc = 0.5 * (ext.c[i] + ext.c[j]);
                const double mrho = 0.5 * (ext.rho[i] + ext.rho[j]);
                const double piv = (beta * muv - alpha * mc) * muv / mrho;
                for (std::size_t d = 0; d < dim; ++d)
                {
                    const double h = -piv * dwdx[d];
                    const double dvxdt_i = ext.mass[j] * h;
                    const double dvxdt_j = ext.mass[i] * h;
                    di.dv[d] += dvxdt_i;
                    dj.dv[d] += -dvxdt_j;
                    const double dedt_i = ext.mass[j] * dvx_iminusj[d] * h;
                    const double dedt_j = ext.mass[i] * dvx_iminusj[d] * h;
                    di.du += -dedt_i; // Acc::add(dedt[i], -(dedt_i)) в art_visc.cpp
                    dj.du += -dedt_j;
                }
            }
        }

        // --- XSPH (av_vel.cpp) ----------------------------------------------
        {
            const double rho_sum = ext.rho[i] + ext.rho[j];
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double dvx_ij = ext.vx[d][i] - ext.vx[d][j];
                const double av_i = 2.0 * ext.mass[j] * dvx_ij / rho_sum * w;
                const double av_j = 2.0 * ext.mass[i] * dvx_ij / rho_sum * w;
                di.xsph[d] += -av_i;
                dj.xsph[d] += av_j;
            }
        }

        const std::size_t pos = offset[k];
        w_out.set(pos, rhs::Contribution<SolidDelta>{
            i, di, rhs::ContributionMeta{kSolidSphInteractionType, static_cast<std::uint64_t>(k), 0}});
        if (j_real)
        {
            w_out.set(pos + 1, rhs::Contribution<SolidDelta>{
                j, dj, rhs::ContributionMeta{kSolidSphInteractionType, static_cast<std::uint64_t>(k), 1}});
        }
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
            if (sys.x[0][i] == 0.0)
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

        double exx = d.dexx, eyy = d.deyy, ezz = d.dezz;
        double exy = d.dexy, exz = d.dexz, eyz = d.deyz;
        const double rxy = d.drxy, rxz = d.drxz, ryz = d.dryz;

        double drho = d.drho;
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
        const double dedt_half = 0.5 * d.du;
        double deviatoric;
        if constexpr (axis_sim && dim == 2)
        {
            deviatoric = (sys.sxx[i] * exx + sys.syy[i] * eyy + sys.szz[i] * ezz + 2.0 * sys.sxy[i] * exy) / sys.rho[i];
        }
        else
        {
            deviatoric = (sys.sxx[i] * exx + sys.syy[i] * eyy + sys.szz[i] * ezz
                         + 2.0 * sys.sxy[i] * exy + 2.0 * sys.syz[i] * eyz + 2.0 * sys.sxz[i] * exz) / sys.rho[i];
        }
        f.du = dedt_half + deviatoric;

        // Закон Гука + производная Яуманна (int_force.cpp, финальный блок).
        if (sys.itype[i] > 2)
        {
            const double G_ = shear_modulus_by_itype(mats, sys.itype[i]);
            if constexpr (dim > 1)
            {
                f.dsxx = 2.0 * G_ * (exx + drho / (3.0 * sys.rho[i])) + 2.0 * (sys.sxy[i] * rxy + sys.sxz[i] * rxz);
                f.dsyy = 2.0 * G_ * (eyy + drho / (3.0 * sys.rho[i])) + 2.0 * (sys.syz[i] * ryz - sys.sxy[i] * rxy);
                // Исправленный знак у syz*ryz (см. tests/test_kinematics_3d.cpp,
                // проверки C2-C4: прежний знак «− syz·ryz» делал след
                // поворотной части dS/dt равным 4·syz·ryz, и при вращении
                // вокруг x девиатор переставал быть бесследовым).
                f.dszz = 2.0 * G_ * (ezz + drho / (3.0 * sys.rho[i])) - 2.0 * (sys.sxz[i] * rxz + sys.syz[i] * ryz);
                f.dsxy = 2.0 * G_ * exy + rxy * (sys.syy[i] - sys.sxx[i]) + sys.sxz[i] * ryz + sys.syz[i] * rxz;
                f.dsyz = 2.0 * G_ * eyz - sys.sxy[i] * rxz - sys.sxz[i] * rxy + ryz * (sys.szz[i] - sys.syy[i]);
                f.dsxz = 2.0 * G_ * exz + rxz * (sys.szz[i] - sys.sxx[i]) - sys.sxy[i] * ryz + sys.syz[i] * rxy;
            }
            else
            {
                f.dsxx = 2.0 * G_ * (exx + drho / (3.0 * sys.rho[i]));
                f.dsyy = 2.0 * G_ * (drho / (3.0 * sys.rho[i]));
                f.dszz = 2.0 * G_ * (drho / (3.0 * sys.rho[i]));
            }
        }
    });
}

} // namespace interactions
