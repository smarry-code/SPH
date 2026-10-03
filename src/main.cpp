// sph.cpp
#include <vector>
#include <omp.h>
#include "param.hpp"
#include "particle_system.hpp"
#include "simulation_config.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <stdexcept>
#include <locale>       // Для setlocale
#include <codecvt>      // Для UTF-8 в Windows
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX       // иначе макросы min/max из windows.h ломают std::min/std::max в заголовках
#endif
#include <windows.h>   // Для SetConsoleOutputCP
#endif
#include "input.hpp"
#include "output.hpp"
#include "time_integration.hpp"

namespace
{
// Замер времени выполнения программы. Раньше жил в linked_list.cpp вместе с
// link_list() (старый поиск соседей); при удалении старого пути
// (fiber/PLAN.md, этап Ф4-B2) единственный оставшийся потребитель —
// main.cpp, поэтому функция перенесена сюда без изменений.
const auto program_start = std::chrono::high_resolution_clock::now();

void time_elapsed(double& elapsed)
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = now - program_start;
    elapsed = duration.count();
}
} // namespace

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    // SimulationConfig: инкапсуляция runtime-параметров 
    // Создаётся здесь, наполняется в input()/parse_sph_config() и передаётся по
    // const& во весь конвейер. Ранее эти параметры были неконстантными глобалами param.hpp.
    SimulationConfig cfg;
    cfg.apply_defaults();
    omp_set_num_threads(cfg.num_threads);   //число потоков из конфига (ранее хардкод)
    std::ofstream file("log.txt", std::ios::out);
    if (file.tellp() == 0)
    {
        file << "\xEF\xBB\xBF"; // BOM для UTF-8
    }
    file.close();
    // Контейнер частиц создаётся пустым: фактический размер устанавливается
    // в input() после определения ntotal. Раньше контейнер аллоцировался
    // сразу под maxn (по умолчанию 5.7 млн) с инициализацией всех элементов —
    // это приводило к избыточной аллокации ~1.7 ГБ и миллионам вызовов
    // inzializ() ещё до загрузки реальных частиц (Этап 3.4).
    // Первичное хранилище частиц — Structure-of-Arrays (Этап 4.3.3). AoS
    // vector<Particle> упразднён: все поля в плоских непрерывных массивах sys.
    ParticleSystem sys;
    // Инициализация динамических массивов
    //std::vector<std::vector<double>> x(dim, std::vector<double>(maxn, 0.0));
    //std::vector<std::vector<double>> vx(dim, std::vector<double>(maxn, 0.0));
    //std::vector<double> mass(maxn, 0.0);
    //std::vector<double> rho(maxn, 0.0);
    //std::vector<double> rho_0(maxn, 0.0);
    //std::vector<double> p(maxn, 0.0);
    //std::vector<double> u(maxn, 0.0);
    //std::vector<double> c(maxn, 0.0);
    //std::vector<double> s(maxn, 0.0);
    //std::vector<double> e(maxn, 0.0);
    //std::vector<double> sxx(maxn, 0.0);
    //std::vector<double> syy(maxn, 0.0);
    //std::vector<double> szz(maxn, 0.0);
    //std::vector<double> sxy(maxn, 0.0);
    //std::vector<double> syz(maxn, 0.0);
    //std::vector<double> sxz(maxn, 0.0);
    //std::vector<double> hsml(maxn, 0.0);
    //std::vector<int> itype(maxn, 0);

    size_t ntotal = 0;
    size_t maxtimestep = 3000;
    double dt = 0.0;
    double time = 0.0;
    double dx = 0.0;
    // Установка шага по времени
    if (shocktube)
    {
        dt = 0.005;
    }
    else if (shearcavity)
    {
        dt = 5.0e-5;
    }

    // Замер времени
    double s1, s2;
    time_print();
    time_elapsed(s1);

    input(cfg, sys, ntotal, dx);
    cfg.maxn = static_cast<size_t>(ntotal * 1.4);
    cfg.max_interaction = 100 * cfg.maxn;
    // Подгоняем размер системы под актуальное число частиц с запасом maxn
    // (запас нужен под возможные виртуальные частицы). Реальные частицы [0, ntotal)
    // уже заполнены в input(); диапазон [ntotal, maxn) инициализируется нулями.
    {
        size_t cap = cfg.maxn;
        if (sys.size() < cap)
        {
            size_t prev = sys.size();
            sys.resize(cap);
            for (size_t i = prev; i < sys.size(); ++i)
            {
                sys.inzializ_index(i);
            }
        }
    }
    time_integration(cfg, sys, ntotal, dt, time, dx);
    output(sys, ntotal, time, maxtimestep, dt);

    time_print();
    time_elapsed(s2);
    std::wcout << "Затраченное время на расчет = " << s2 - s1 << " секунд" << std::endl;


    return 0;
}
