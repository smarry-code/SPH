// single_step.cpp
#include "direct_find.hpp"
#include "linked_list.hpp"
#include "density.hpp"
#include "int_force.hpp"
#include "art_visc.hpp"
#include "virt_part.hpp"
#include "av_vel.hpp"
#include "ext_force.hpp"
#include "h_upgrade.hpp"
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
void tree_search()
{
    std::cout << "tree_search (stub) called" << std::endl;
}

void viscosity()
{
    std::cout << "viscosity (stub) called" << std::endl;
}

void art_heat()
{
    std::cout << "art_heat (stub) called" << std::endl;
}

// Основная функция single_step
void single_step(std::vector<Particle>& Elem, size_t itimestep, double dt, size_t ntotal, std::vector<double>& t, std::vector<double>& tdsdt,
                 std::vector<std::vector<double>>& dx, std::vector<std::vector<double>>& dvx, std::vector<double>& du, std::vector<double>& ds,
                 std::vector<double>& drho, std::vector<std::vector<double>>& av, std::vector<double>& dsxx, std::vector<double>& dsyy,
                 std::vector<double>& dszz, std::vector<double>& dsxy, std::vector<double>& dsyz, std::vector<double>& dsxz)
{

    // Инициализация массивов без виртуальных частиц
    std::vector<double> avdudt(ntotal, 0.0);
    std::vector<double> ahdudt(ntotal, 0.0);
    std::vector<std::vector<double>> indvxdt(dim, std::vector<double>(ntotal, 0.0));
    std::vector<std::vector<double>> ardvxdt(dim, std::vector<double>(ntotal, 0.0));
    std::vector<std::vector<double>> exdvxdt(dim, std::vector<double>(ntotal, 0.0));

    std::vector<size_t> pair_i(max_interaction);
    std::vector<size_t> pair_j(max_interaction);
    std::vector<double> w(max_interaction);
    std::vector<std::vector<double>> dwdx(dim, std::vector<double>(max_interaction, 0.0));
    std::vector<int> ns(ntotal, 0);
    std::vector<double> eta(ntotal, 0.0);
    size_t niac = 0;
    size_t nvirt = 0;
    double s1, s2 = 0.0;
    // Генерация виртуальных частиц
    /*if (axis_sim == true && dim == 2)
    {
        virt_part(Elem, itimestep, ntotal, nvirt);
    }*/
    //virt_part(Elem, ntotal, nvirt);
    // Поиск соседей
    if (nnps == 1)
    {
        time_elapsed(s1);
        direct_find(Elem, ntotal + nvirt, niac, pair_i, pair_j, w, dwdx, ns);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на direct_find = " << s2 - s1 << " секунд" << std::endl;
    }
    else if (nnps == 2)
    {
        time_elapsed(s1);
        link_list(Elem, itimestep, ntotal + nvirt, niac, pair_i, pair_j, w, dwdx);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на link_list = " << s2 - s1 << " секунд" << std::endl;
        //std::wcout << "Количество виртуальных частиц = " << nvirt << std::endl;
    }
    else if (nnps == 3)
    {
        tree_search();
    }

    // Аппроксимация плотности
    if (summation_density)
    {
        sum_density(Elem, ntotal + nvirt, niac, pair_i, pair_j, w);
    }
    else
    {
        time_elapsed(s1);
        con_density(Elem, ntotal + nvirt, niac, pair_i, pair_j, dwdx, drho);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на con_density = " << s2 - s1 << " секунд" << std::endl;
    }

    // Расчет физической вязкости
    if (visc)
    {
        viscosity();
    }

    // Внутренние силы
    time_elapsed(s1);
    int_force(Elem, ntotal + nvirt, niac, eta, pair_i, pair_j, dwdx, w, t, indvxdt, tdsdt, du, dsxx, dsyy, dszz, dsxy, dsyz, dsxz, drho);
    time_elapsed(s2);
    //std::wcout << "Затраченное время на int_force = " << s2 - s1 << " секунд" << std::endl;
    // Искусственная вязкость
    if (visc_artificial)
    {
        time_elapsed(s1);
        art_visc(Elem, ntotal + nvirt, niac, pair_i, pair_j, dwdx, ardvxdt, avdudt);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на art_visc = " << s2 - s1 << " секунд" << std::endl;
    }
    // Внешние силы
    if (ex_force)
    {
        ext_force(Elem, ntotal + nvirt, niac, pair_i, pair_j, exdvxdt);
    }
    // Обновление длины сглаживания
    if (sle != 0)
    {
        time_elapsed(s1);
        h_upgrade(Elem, dt, ntotal + nvirt, niac, pair_i, pair_j, dwdx);
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
        av_vel(Elem, ntotal + nvirt, niac, pair_i, pair_j, w, av);
        time_elapsed(s2);
        //std::wcout << "Затраченное время на av_vel = " << s2 - s1 << " секунд" << std::endl;
    }
    // Комбинирование всех сил
    for (size_t i = 0; i < ntotal; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            dvx[d][i] = indvxdt[d][i] + exdvxdt[d][i] + ardvxdt[d][i];
            dx[d][i] = Elem[i].vx[d]; // dx/dt = vx
        }
        du[i] = du[i] + avdudt[i] + ahdudt[i];
    }
    for (size_t i = ntotal; i < Elem.size(); ++i)
    {
        Elem[i].obnul();
    }
    //Elem.resize(ntotal);

}
