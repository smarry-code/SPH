// single_step.cpp
// Определения time_print и оркестратора одного шага симуляции (single_step).
// Объявления — в single_step.hpp.
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

#include <chrono>
#include <ctime>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
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

// Заглушки для недостающих функций
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

// Основная функция single_step
void single_step(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, double dt, size_t ntotal, std::vector<double>& t, std::vector<double>& tdsdt,
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
        ext_force(cfg, sys, dt, ntotal + nvirt, niac, pair_i, pair_j, exdvxdt, wb);
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
