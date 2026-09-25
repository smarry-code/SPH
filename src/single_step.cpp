// single_step.cpp
// Определения time_print и оркестратора одного шага симуляции (single_step).
// Объявления — в single_step.hpp.
//
// fiber/PLAN.md, этап Ф4-B1: single_step() ниже — НОВЫЙ путь RHS для Solid
// (interactions/solid_sph.hpp, interactions/solid_contact.hpp), подключённый
// к продовому циклу. Обе перегрузки time_integration (ParticleSystem&,
// fiber::Model&, src/time_integration.cpp) вызывают именно эту функцию,
// поэтому переключение её тела переключает оба цикла одновременно — это
// единственная точка врезки, минимально затрагивающая остальной код.
//
// Старый путь (direct_find/link_list -> con_density -> int_force -> art_visc
// -> ext_force -> h_upgrade -> av_vel) сохранён НИЖЕ под именем
// single_step_legacy — компилируется, но не вызывается (fiber/PLAN.md, этап
// Ф4-B1: "сохранить старый single_step... физически доступным... не
// удалять"). Удаление — часть Ф4-B2, после независимого ревью части B1.
#include "single_step.hpp"
#include "direct_find.hpp"
#include "linked_list.hpp"
#include "density.hpp"
#include "int_force.hpp"
#include "art_visc.hpp"
#include "virt_part.hpp"
#include "av_vel.hpp"
#include "ext_force.hpp"
#include "h_upgrade.hpp"

#include "interactions/solid_contact.hpp"
#include "interactions/solid_sph.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <chrono>
#include <ctime>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <omp.h>

#ifdef _WIN32
#include <windows.h>
#endif

void time_print()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    // Получаем текущее системное время
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    // Создаем буфер для форматированного времени
    char buffer[80];

    // Преобразуем в структуру tm
#ifdef _WIN32
    struct tm timeinfo;
    localtime_s(&timeinfo, &now_c);
#else
    struct tm timeinfo;
    localtime_r(&now_c, &timeinfo);
#endif

    // Форматируем время в буфер с помощью strftime
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", &timeinfo);

    // Выводим отформатированное время
    std::wcout << "Текущее время: " << buffer << std::endl;
}

// Заглушки для недостающих функций (используются только single_step_legacy)
static void tree_search()
{
    std::cout << "tree_search (stub) called" << std::endl;
}

static void viscosity()
{
    std::cout << "viscosity (stub) called" << std::endl;
}

static void art_heat()
{
    std::cout << "art_heat (stub) called" << std::endl;
}

// ============================================================================
// single_step — НОВЫЙ путь (fiber/PLAN.md, этап Ф4-B1).
//
// Граф вычисления (fiber/PLAN.md §3.4): Solid derived (EOS) -> GhostSet ->
// SolidExtended -> поиск пар (SpatialIndex) -> Solid-SPH kernel + Solid-
// Solid-Contact kernel -> reduce_contributions -> finalize_solid -> запись
// в выходные параметры + шим совместимости в sys.
//
// ШИМ СОВМЕСТИМОСТИ (решение координатора, НЕ часть Ф7): p, c (EOS) и
// exx..eyz, rxy..ryz (скорости деформаций/поворота) по-прежнему пишутся в
// sys, хотя физически они теперь считаются в RhsWorkspace-эквиваленте
// (SolidDerived/SolidFinal). Причина — src/ppochnost.cpp (prochnost_mises,
// prochnost_cem_JH2) и compute_cfl_timestep (src/time_integration.cpp)
// читают их из sys НАПРЯМУЮ и вызываются time_integration.cpp ПОСЛЕ
// single_step на каждом шаге; переносить их читателей на новую схему —
// содержание этапа Ф7 (fiber/PLAN.md §3.6, раздел «### Ф7. Postprocessing»:
// «перенос прочности (Мизес, JH-2)»), не Ф4. Шим временный и будет убран
// вместе с этим переносом.
//
// h_upgrade (sle) НЕ входит в новый RHS-путь по плану (переезжает в
// Postprocessing на Ф7) — вызывается как раньше, старым кодом. sle == 1
// (BC2/BC6 и все текущие конфиги) не использует пары соседей (h == fac*
// (mass/rho)^(1/dim), src/h_upgrade.cpp) — вызывается с пустым списком пар.
// sle == 2 пары ИСПОЛЬЗУЕТ (сумма div v по niac), а новый путь их в формате
// старого поиска не строит — обнаружено явно исключением, а не молча
// неверным результатом (kNoLimit-подобная защита).
// ============================================================================
void single_step(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, double dt, size_t ntotal, std::vector<double>& t, std::vector<double>& tdsdt,
                 std::vector<std::vector<double>>& dx, std::vector<std::vector<double>>& dvx, std::vector<double>& du, std::vector<double>& ds,
                 std::vector<double>& drho, std::vector<std::vector<double>>& av, std::vector<double>& dsxx, std::vector<double>& dsyy,
                 std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, WorkBuffers& wb,
                 PairGraph& graph)
{
    // itimestep — нужен только link_list (старый поиск соседей с кэшем по
    // шагам, nnps == 2); новый путь строит SpatialIndex заново на каждом
    // вызове. graph (PairGraph/раскраска) — инфраструктура atomic-ветви
    // старого int_force/art_visc/av_vel; rhs::reduce_contributions её не
    // использует (reduce по шардам owner-computes, fiber/PLAN.md §3.5).
    // t, ds — унаследованные от старого single_step мёртвые выходные
    // параметры (ни старый int_force, ни time_integration.cpp их не
    // читают после вызова single_step) — сохранены нетронутыми только ради
    // неизменности сигнатуры (её меняют оба вызова в time_integration.cpp).
    (void)itimestep;
    (void)graph;
    (void)t;
    (void)ds;

    if (cfg.sle == 2)
    {
        throw std::runtime_error(
            "single_step: sle=2 требует сумму по парам соседей (div v), которую "
            "новый путь RHS (fiber/PLAN.md, этап Ф4-B1) не строит в формате "
            "старого поиска для h_upgrade; ни один текущий конфиг (BC2/BC6 -> "
            "sle=1) sle=2 не использует. Поддержать sle=2 в новом пути — часть "
            "переноса h_upgrade в Postprocessing (fiber/PLAN.md, этап Ф7).");
    }

    // CpuExecutor — один на вызов (num_threads=0 => omp_get_max_threads(),
    // тот же счётчик потоков, что main.cpp выставляет один раз через
    // omp_set_num_threads(cfg.num_threads), rhs/cpu_executor.hpp). Часть A
    // (test_rhs_equivalence) строит executor так же.
    rhs::CpuExecutor exec;

    // --- Solid derived: EOS (p, c) -- fiber/PLAN.md §3.6 -------------------
    interactions::SolidDerived derived;
    interactions::compute_solid_derived(exec, cfg.mats, sys, ntotal, derived);

    // --- GhostSet + SolidExtended -- замена virt_part.cpp -------------------
    interactions::GhostSet ghosts;
    interactions::build_ghost_set(sys, ntotal, ghosts);

    interactions::SolidExtended ext;
    interactions::build_solid_extended(sys, ntotal, derived, ghosts, ext);

    // --- Поиск пар Solid-SPH (общий для Solid-SPH и Solid-Solid-Contact,
    // fiber/SPH.md: "Одна и та же построенная spatial structure может
    // использоваться несколькими interaction types") --------------------
    interactions::SolidSphSearchResult search;
    interactions::build_solid_sph_search(exec, cfg, ext, search);

    // --- Compute interaction kernels -> contributions -----------------------
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> sph_buffers, contact_buffers;
    interactions::compute_solid_sph_contributions(exec, cfg, ext, search.primitives, sph_buffers);
    interactions::compute_solid_contact_contributions(exec, cfg, cfg.mats, sys, derived.c, dt,
                                                        search.real_pairs, contact_buffers);

    // --- Reduce (owner-computes по шардам, fiber/PLAN.md §3.5) -------------
    std::vector<rhs::ContributionBuffer<interactions::SolidDelta>> all_buffers = sph_buffers;
    all_buffers.insert(all_buffers.end(), contact_buffers.begin(), contact_buffers.end());
    const std::vector<interactions::SolidDelta> reduced = rhs::reduce_contributions(exec, all_buffers, ntotal);

    // --- Finalize (parallel_for по реальным частицам) -----------------------
    std::vector<interactions::SolidFinal> finalres;
    interactions::finalize_solid(exec, cfg.mats, sys, ntotal, reduced, finalres);

    // --- Диагностика контакта -- заполняет wb (ADR-007), как старый
    // ext_force: min_dt_contact ограничивает шаг в time_integration.cpp,
    // contact_pair_count/contact_max_delta -- вывод прогресса. Вызывается
    // каждый шаг (не "по запросу", в отличие от комментария в
    // interactions/solid_contact.hpp у compute_solid_contact_diagnostics) --
    // ADR-007 требует min_dt_contact СЕЙЧАС же ограничивающим ПРОДОВЫЙ шаг,
    // а не по отдельному запросу диагностики.
    wb.reset_contact();
    if (cfg.contact_force_type != contact::Off)
    {
        const interactions::SolidContactDiagnostics diag = interactions::compute_solid_contact_diagnostics(
            exec, cfg, cfg.mats, sys, derived.c, dt, search.real_pairs);
        wb.contact_pair_count = diag.pair_count;
        wb.contact_max_delta = diag.max_delta;
        wb.min_dt_contact = diag.min_dt;
    }

    // --- Запись в выходные параметры + шим совместимости в sys -------------
    for (size_t i = 0; i < ntotal; ++i)
    {
        const interactions::SolidFinal& f = finalres[i];
        drho[i] = f.drho;
        du[i] = f.du;
        for (size_t d = 0; d < dim; ++d)
        {
            dvx[d][i] = f.dv[d];
            dx[d][i] = sys.vx[d][i]; // dx/dt = vx, как в старом пути
            // XSPH: СЫРАЯ сумма (SolidDelta::xsph) -> *cfg.epsilon, как
            // старый av_vel.cpp ("Масштабирование средней скорости").
            av[d][i] = reduced[i].xsph[d] * cfg.epsilon;
        }
        dsxx[i] = f.dsxx; dsyy[i] = f.dsyy; dszz[i] = f.dszz;
        dsxy[i] = f.dsxy; dsyz[i] = f.dsyz; dsxz[i] = f.dsxz;

        // Шим совместимости (см. комментарий у функции): p, c -- EOS этого
        // RHS; exx..eyz, rxy..ryz -- сырые скорости деформаций/поворота
        // (с осевым переопределением ezz при axis_sim && dim == 2, уже
        // применённым в finalize_solid). Читатели: src/ppochnost.cpp
        // (prochnost_cem_JH2), compute_cfl_timestep (sys.c).
        sys.p[i] = derived.p[i];
        sys.c[i] = derived.c[i];
        sys.exx[i] = f.exx; sys.eyy[i] = f.eyy; sys.ezz[i] = f.ezz;
        sys.exy[i] = f.exy; sys.exz[i] = f.exz; sys.eyz[i] = f.eyz;
        sys.rxy[i] = f.rxy; sys.rxz[i] = f.rxz; sys.ryz[i] = f.ryz;

        tdsdt[i] = 0.0; // физическая вязкость (visc) выключена constexpr -- как в старом int_force
    }

    // h_upgrade -- НЕ часть нового RHS-пути (fiber/PLAN.md §3.6: "sle=1 --
    // Postprocessing" на этапе Ф7; сейчас остаётся старым кодом). sle == 1
    // не читает pair_i/pair_j/dwdx (src/h_upgrade.cpp) -- пустые буферы
    // корректны; sle == 0 -- no-op; sle == 2 отклонён исключением выше.
    if (cfg.sle == 1)
    {
        static const std::vector<size_t> kEmptyPairIdx;
        static const std::vector<std::vector<double>> kEmptyDwdx(dim);
        h_upgrade(cfg, sys, dt, ntotal, 0, kEmptyPairIdx, kEmptyPairIdx, kEmptyDwdx);
    }

    // Единая диагностика состояния (подэтап 5.8, перенесена из старого
    // single_step): read-only проверка конечности p/rho/c/u/hsml после шага.
    Diagnostics::validate(sys, ntotal, "post-solid-rhs");

    // Старый путь дозаписывал зеркальные частицы прямо в sys (virt_part.cpp)
    // и поэтому обнулял [ntotal, sys.size()) в конце каждого шага, чтобы не
    // накапливать мусор. Новый путь НИКОГДА не пишет в sys за пределами
    // [0, ntotal) (зеркала живут только в interactions::SolidExtended,
    // fiber/PLAN.md §3.6: "виртуальные частицы -> GhostSet + SolidExtended,
    // а не дозапись в sys"), поэтому очищать там нечего.
}

// ============================================================================
// single_step_legacy — старый путь, побитная копия прежнего тела single_step
// (fiber/PLAN.md, этап Ф4-B1: "сохранить старый single_step физически
// доступным... для дальнейшего сравнения в части B2/на случай отката").
// НЕ вызывается ни из одной перегрузки time_integration.
// TODO(Ф4-B2): удалить вместе со старым путём (direct_find, link_list,
// int_force, art_visc, av_vel, ext_force, virt_part) после независимого
// ревью части B1.
// ============================================================================
void single_step_legacy(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, double dt, size_t ntotal, std::vector<double>& t, std::vector<double>& tdsdt,
                 std::vector<std::vector<double>>& dx, std::vector<std::vector<double>>& dvx, std::vector<double>& du, std::vector<double>& ds,
                 std::vector<double>& drho, std::vector<std::vector<double>>& av, std::vector<double>& dsxx, std::vector<double>& dsyy,
                 std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, WorkBuffers& wb,
                 PairGraph& graph)
{
    // Локальные псевдонимы параметров (тело функции не изменяется):
    // max_interaction, nnps, sle, coloring_niac_threshold, print_step ранее — глобалы param.hpp.
    const size_t max_interaction = cfg.max_interaction;
    const int nnps = cfg.nnps;
    const int sle = cfg.sle;
    const int coloring_niac_threshold = cfg.coloring_niac_threshold;
    const size_t print_step = cfg.print_step;

    size_t niac = 0;
    size_t nvirt = 0;
    double s1, s2 = 0.0;

    // Генерация виртуальных (зеркальных) частиц у оси симметрии — ДО выделения
    // рабочих массивов, потому что их размер зависит от nvirt.
    // При axis_sim && dim == 2 приосевой столбец иначе считается с неполным
    // носителем ядра: соседей со стороны r < 0 нет, и осевой член rho*v_r/r
    // в уравнении неразрывности (density.cpp) ничем не уравновешен.
    // Отбор и зеркалирование — внутри virt_part (условие axis_sim && dim == 2
    // проверяется там же, поэтому при dim == 1 вызов даёт nvirt = 0).
    if (virtual_part)
    {
        virt_part(sys, ntotal, nvirt, cfg.maxn);
    }

    // Рабочие массивы охватывают реальные И виртуальные частицы: все узлы
    // конвейера ниже вызываются с ntotal + nvirt и индексируют массивы этим
    // диапазоном. Пока virt_part был отключён, nvirt был тождественно нулём и
    // размера ntotal хватало; с его восстановлением ntotal привёл бы к выходу
    // за границы во всех потребителях сразу.
    const size_t nall = ntotal + nvirt;
    std::vector<double> avdudt(nall, 0.0);
    std::vector<double> ahdudt(nall, 0.0);
    std::vector<std::vector<double>> indvxdt(dim, std::vector<double>(nall, 0.0));
    std::vector<std::vector<double>> ardvxdt(dim, std::vector<double>(nall, 0.0));
    std::vector<std::vector<double>> exdvxdt(dim, std::vector<double>(nall, 0.0));

    std::vector<size_t> pair_i(max_interaction);
    std::vector<size_t> pair_j(max_interaction);
    std::vector<double> w(max_interaction);
    std::vector<std::vector<double>> dwdx(dim, std::vector<double>(max_interaction, 0.0));
    std::vector<int> ns(nall, 0);
    std::vector<double> eta(nall, 0.0);
    // Поиск соседей
    if (nnps == 1)
    {
        time_elapsed(s1);
        direct_find(cfg, sys, ntotal + nvirt, niac, pair_i, pair_j, w, dwdx, ns);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на direct_find = " << s2 - s1 << " секунд" << std::endl;
    }
    else if (nnps == 2)
    {
        time_elapsed(s1);
        link_list(cfg, sys, itimestep, ntotal + nvirt, niac, pair_i, pair_j, w, dwdx);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на link_list = " << s2 - s1 << " секунд" << std::endl;
        //std::wcout << "Количество виртуальных частиц = " << nvirt << std::endl;
    }
    else if (nnps == 3)
    {
        tree_search();
    }

    // Построение графа пар + жадная окраска (один раз за шаг, после поиска соседей).
    // Переиспользуется всеми потребителями пар: пары одного цвета не разделяют
    // частиц => конкурентной записи нет => #pragma omp atomic устранены.
    // Условный путь: graph coloring включается только при niac >= порога.
    // При меньших niac (текущие кейсы N<=5000) работает atomic-ветвь потребителей
    // (graph.num_colors = 0 — сигнал «без окраски»).
    if (niac >= static_cast<size_t>(coloring_niac_threshold))
    {
#ifdef SPH_PROFILE
        double _gb0 = omp_get_wtime();
#endif
        graph.build(pair_i, pair_j, niac, ntotal + nvirt);
#ifdef SPH_PROFILE
        double _gb1 = omp_get_wtime();
        static double g_build_time = 0.0;
        static size_t g_build_calls = 0;
        g_build_time += (_gb1 - _gb0);
        ++g_build_calls;
        if (g_build_calls % print_step == 0)
        {
            std::ofstream _bf("force_timing.csv", std::ios::app);
            if (_bf)
            {
                _bf << "build," << g_build_calls << "," << graph.num_colors << ","
                    << (g_build_time / static_cast<double>(g_build_calls)) * 1e6 << "\n";
            }
        }
#endif
    }
    else
    {
        graph.num_colors = 0;  // сигнал потребителям: использовать atomic-ветвь
    }

    // Аппроксимация плотности
    if (summation_density)
    {
        sum_density(cfg, sys, ntotal + nvirt, niac, pair_i, pair_j, w);
    }
    else
    {
        time_elapsed(s1);
        con_density(cfg, sys, ntotal + nvirt, niac, pair_i, pair_j, dwdx, drho, graph);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на con_density = " << s2 - s1 << " секунд" << std::endl;
    }

    // Единая диагностика состояния (подэтап 5.8): read-only проверка конечности
    // ключевых термодинамических полей после этапа плотности. Печать в stderr
    // только при NaN/Inf; численно нейтральна (поля sys не изменяет).
    Diagnostics::validate(sys, ntotal, "post-density");

    // Расчет физической вязкости
    if (visc)
    {
        viscosity();
    }

    // Внутренние силы
    time_elapsed(s1);
    int_force(cfg, sys, ntotal + nvirt, niac, eta, pair_i, pair_j, dwdx, w, t, indvxdt, tdsdt, du, dsxx, dsyy, dszz, dsxy, dsyz, dsxz, drho, wb, graph);
    time_elapsed(s2);
    //std::wcout << "Затраченное время на int_force = " << s2 - s1 << " секунд" << std::endl;
    // Единая диагностика (подэтап 5.8) после сил внутреннего взаимодействия.
    Diagnostics::validate(sys, ntotal, "post-int_force");
    // Искусственная вязкость
    if (visc_artificial)
    {
        time_elapsed(s1);
        art_visc(cfg, sys, ntotal + nvirt, niac, pair_i, pair_j, dwdx, ardvxdt, avdudt, graph);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на art_visc = " << s2 - s1 << " секунд" << std::endl;
    }
    // Внешние силы
    if (ex_force)
    {
        // Полный диапазон — ntotal + nvirt (по нему размерен exdvxdt), реальных
        // из них — ntotal: контакт не должен видеть зеркальные частицы.
        ext_force(cfg, sys, dt, ntotal + nvirt, ntotal, niac, pair_i, pair_j, exdvxdt, wb);
    }
    // Обновление длины сглаживания
    if (sle != 0)
    {
        time_elapsed(s1);
        h_upgrade(cfg, sys, dt, ntotal + nvirt, niac, pair_i, pair_j, dwdx);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на h_upgrade = " << s2 - s1 << " секунд" << std::endl;
    }
    // Искусственный нагрев
    if (heat_artificial)
    {
        art_heat();
    }
    // Средняя скорость для предотвращения проникания
    if (average_velocity)
    {
        time_elapsed(s1);
        av_vel(cfg, sys, ntotal + nvirt, niac, pair_i, pair_j, w, av, graph);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на av_vel = " << s2 - s1 << " секунд" << std::endl;
    }
    // Комбинирование всех сил
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            dvx[d][i] = indvxdt[d][i] + exdvxdt[d][i] + ardvxdt[d][i];
            dx[d][i] = sys.vx[d][i]; // dx/dt = vx
        }
        du[i] = du[i] + avdudt[i] + ahdudt[i];
    }
    for (size_t i = ntotal; i < sys.size(); ++i)
    {
        sys.obnul_index(i);
    }
    //sys.resize(ntotal);

}
