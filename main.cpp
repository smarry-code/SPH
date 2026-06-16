// sph.cpp
#include <vector>
#include <omp.h>
#include "param.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <stdexcept>
#include <locale>       // Для setlocale
#include <codecvt>      // Для UTF-8 в Windows
#include <windows.h>   // Для SetConsoleOutputCP
#include "input.hpp"
#include "output.hpp"
#include "time_integration.hpp"

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    omp_set_num_threads(7);
    std::ofstream file("log.txt", std::ios::out);
    if (file.tellp() == 0)
    {
        file << "\xEF\xBB\xBF"; // BOM для UTF-8
    }
    file.close();
    std::vector <Particle> Elem(maxn);
    for (size_t i = 0; i < Elem.size(); i++)
    {
        Elem[i].inzializ();
    }
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
    int maxtimestep = 3000;
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

    input(Elem, ntotal, dx);
    maxn = ntotal * 1.4;
    max_interaction = 100 * maxn;
    time_integration(Elem, ntotal, dt, time, dx);
    output(Elem, ntotal, time, maxtimestep, dt);

    time_print();
    time_elapsed(s2);
    std::wcout << "Затраченное время на расчет = " << s2 - s1 << " секунд" << std::endl;


    return 0;
}
