// single_step.cpp
// Определения time_print и оркестратора одного шага симуляции (single_step).
// Объявления — в single_step.hpp.
//
// fiber/PLAN.md, этап Ф4: single_step() ниже — путь RHS для Solid
// (interactions/solid_sph.hpp, interactions/solid_contact.hpp), подключённый
// к продовому циклу. Обе перегрузки time_integration (ParticleSystem&,
// fiber::Model&, src/time_integration.cpp) вызывают именно эту функцию,
// поэтому переключение её тела переключает оба цикла одновременно — это
// единственная точка врезки, минимально затрагивающая остальной код.
//
// Старый путь (direct_find/link_list -> con_density -> int_force -> art_visc
// -> ext_force -> h_upgrade -> av_vel) удалён на этапе Ф4-B2 после
// независимого ревью части Ф4-B1 (переключение продового расчёта) и переноса
// tests/test_axisym.cpp / tests/test_kinematics_3d.cpp на этот путь.
#include "single_step.hpp"
#include "ext_force.hpp"

#include "interactions/solid_contact.hpp"
#include "interactions/solid_sph.hpp"
#include "interactions/solid_workspace.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"

#include <chrono>
#include <ctime>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <memory>
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

// ============================================================================
// single_step — путь RHS для Solid (fiber/PLAN.md, этап Ф4).
//
// Граф вычисления (fiber/PLAN.md §3.4): Solid derived (EOS) -> GhostSet ->
// SolidExtended -> поиск пар (SpatialIndex) -> Solid-SPH kernel + Solid-
// Solid-Contact kernel -> reduce_contributions -> finalize_solid -> запись
// в выходные параметры + шим совместимости в sys.
//
// ШИМ СОВМЕСТИМОСТИ (решение координатора этапа Ф7 — оставлен сознательно,
// НЕ убран вместе с переносом прочности): p, c (EOS) и exx..eyz, rxy..ryz
// (скорости деформаций/поворота) по-прежнему пишутся в sys, хотя физически
// они теперь считаются в RhsWorkspace-эквиваленте (SolidDerived/SolidFinal).
// prochnost_cem_JH2 (src/ppochnost.cpp) по-прежнему читает 9 величин
// деформации/поворота из sys напрямую (сигнатуру не меняли — риск регрессии
// от смены источника данных не оправдан архитектурной чистотой на этом
// этапе, fiber/PLAN.md §Ф7 явно допускает оба варианта); compute_cfl_timestep
// (src/time_integration.cpp) читает sys.c. Оба вызываются из
// interactions::postprocessing (interactions/postprocessing.hpp) ПОСЛЕ
// single_step на каждом шаге.
//
// h_upgrade (sle) — этап Ф7 (fiber/PLAN.md §3.8: «sle=1 — Postprocessing")
// вынес вызов ИЗ single_step в interactions::update_smoothing_length
// (interactions/postprocessing.hpp), вызываемую из ОБЕИХ перегрузок
// time_integration.cpp СРАЗУ после single_step() — см. там же подробное
// обоснование, ПОЧЕМУ именно это место (а не буквально "после corrector'а")
// обязательно для нулевого изменения численного результата.
// ============================================================================
void single_step(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, double dt, size_t ntotal, std::vector<double>& t, std::vector<double>& tdsdt,
                 std::vector<std::vector<double>>& dx, std::vector<std::vector<double>>& dvx, std::vector<double>& du, std::vector<double>& ds,
                 std::vector<double>& drho, std::vector<std::vector<double>>& av, std::vector<double>& dsxx, std::vector<double>& dsyy,
                 std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz, WorkBuffers& wb)
{
    // itimestep — нужен был только link_list (старый поиск соседей с кэшем
    // по шагам, nnps == 2, удалён на Ф4-B2); новый путь строит SpatialIndex
    // заново на каждом вызове. t, ds — унаследованные от старого single_step
    // мёртвые выходные параметры (их не читал ни старый int_force, ни
    // time_integration.cpp) — сохранены нетронутыми только ради неизменности
    // сигнатуры (её меняют оба вызова в time_integration.cpp).
    (void)itimestep;
    (void)t;
    (void)ds;

    // CpuExecutor — один на вызов (num_threads=0 => omp_get_max_threads(),
    // тот же счётчик потоков, что main.cpp выставляет один раз через
    // omp_set_num_threads(cfg.num_threads), rhs/cpu_executor.hpp).
    rhs::CpuExecutor exec;

    // Рабочая память шага живёт в wb между вызовами (interactions/
    // solid_workspace.hpp) -- создаётся при первом вызове, дальше только
    // переиспользуется.
    if (!wb.solid_ws) wb.solid_ws = std::make_shared<interactions::SolidStepWorkspace>();
    interactions::SolidStepWorkspace& ws = *static_cast<interactions::SolidStepWorkspace*>(wb.solid_ws.get());

    // --- Solid derived: EOS (p, c) -- fiber/PLAN.md §3.6 -------------------
    interactions::SolidDerived& derived = ws.derived;
    interactions::compute_solid_derived(exec, cfg.mats, sys, ntotal, derived);

    // --- GhostSet + SolidExtended -- замена virt_part.cpp -------------------
    interactions::GhostSet& ghosts = ws.ghosts;
    interactions::build_ghost_set(sys, ntotal, ghosts);

    interactions::SolidExtended& ext = ws.ext;
    interactions::build_solid_extended(sys, ntotal, derived, ghosts, ext);

    // --- Поиск пар Solid-SPH (общий для Solid-SPH и Solid-Solid-Contact,
    // fiber/SPH.md: "Одна и та же построенная spatial structure может
    // использоваться несколькими interaction types") --------------------
    interactions::SolidSphSearchResult& search = ws.search;
    interactions::build_solid_sph_search(exec, cfg, ext, search);

    // --- Compute interaction kernels -> contributions -----------------------
    interactions::compute_solid_sph_contributions(exec, cfg, ext, search.primitives, ws.sph_buffers);
    interactions::compute_solid_contact_contributions(exec, cfg, cfg.mats, sys, derived.c, dt,
                                                        search.real_pairs, ws.contact_buffers);

    // --- Reduce (owner-computes по шардам, fiber/PLAN.md §3.5) -------------
    // Наборы буферов склеиваются логически (без копирования), порядок прежний:
    // сначала Solid-SPH, затем контакт.
    rhs::reduce_contributions_into(exec, {&ws.sph_buffers, &ws.contact_buffers}, ntotal, ws.reduced,
                                   ws.reduce_scratch);
    const std::vector<interactions::SolidDelta>& reduced = ws.reduced;

    // --- Finalize (parallel_for по реальным частицам) -----------------------
    std::vector<interactions::SolidFinal>& finalres = ws.finalres;
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

    // h_upgrade (sle) -- этап Ф7: вызов вынесен ИЗ single_step в
    // interactions::update_smoothing_length (interactions/postprocessing.hpp),
    // вызываемую из time_integration.cpp сразу ПОСЛЕ single_step() (см.
    // комментарий у update_smoothing_length -- почему именно там, а не
    // буквально в Postprocessing-стадии после corrector'а). Diagnostics::
    // validate ниже поэтому проверяет sys.hsml СО ЗНАЧЕНИЕМ ПРЕДЫДУЩЕГО шага
    // (ещё не обновлённым на этом шаге) -- чисто диагностическая проверка на
    // конечность, на численный результат не влияет.

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

