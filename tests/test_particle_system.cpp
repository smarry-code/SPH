// test_particle_system.cpp
// Unit-тест ParticleSystem (SoA) — первичного хранилища частиц (Этап 4.3.3).
// Мост AoS<->SoA упразднён; тест проверяет актуальные методы sys как primary:
// resize, inzializ_index (обнуление), vuvod_index (вывод), inzializ_2d_symm_index
// (осесимметричное зеркалирование), корректность плоского доступа sys.field[i].
#include "particle_system.hpp"
#include <cmath>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

int main()
{
    int failures = 0;
    ParticleSystem sys;

    // --- resize + size ---
    const size_t N = 50;
    sys.resize(N);
    if (sys.size() != N) { std::printf("FAIL size: %zu != %zu\n", sys.size(), N); ++failures; }

    // --- inzializ_index: все поля должны быть 0 ---
    for (size_t i = 0; i < N; ++i)
        sys.inzializ_index(i);
    for (size_t i = 0; i < N; ++i)
    {
        if (sys.mass[i] != 0.0 || sys.rho[i] != 0.0 || sys.p[i] != 0.0 ||
            sys.sxx[i] != 0.0 || sys.exx[i] != 0.0 || sys.hsml[i] != 0.0)
        { std::printf("FAIL inzializ_index: ненулевое поле у частицы %zu\n", i); ++failures; break; }
        for (size_t d = 0; d < dim; ++d)
            if (sys.x[d][i] != 0.0 || sys.vx[d][i] != 0.0)
            { std::printf("FAIL inzializ_index: ненулевая kinematics\n"); ++failures; break; }
    }

    // --- плоский доступ: запись и чтение совпадают ---
    for (size_t i = 0; i < N; ++i)
    {
        sys.rho[i]  = 1000.0 + static_cast<double>(i);
        sys.mass[i] = 1e-6 * static_cast<double>(i + 1);
        sys.itype[i] = 3;
        for (size_t d = 0; d < dim; ++d)
            sys.x[d][i] = static_cast<double>(i) + static_cast<double>(d) * 0.1;
    }
    for (size_t i = 0; i < N; ++i)
    {
        if (sys.rho[i] != 1000.0 + static_cast<double>(i)) { std::printf("FAIL плоский доступ rho[%zu]\n", i); ++failures; break; }
        if (sys.x[0][i] != static_cast<double>(i)) { std::printf("FAIL плоский доступ x[0][%zu]\n", i); ++failures; break; }
    }

    // --- vuvod_index: вывод частицы в файл без исключения ---
    {
        // vuvod_index принимает std::ofstream; проверяем через реальный файл.
        std::ofstream tmp("__ps_test_out.txt");
        sys.vuvod_index(0, tmp);
        tmp.close();
        std::ifstream chk("__ps_test_out.txt");
        std::string content;
        std::getline(chk, content);
        chk.close();
        if (content.empty()) { std::printf("FAIL vuvod_index: пустой вывод\n"); ++failures; }
    }

    // --- поля Джонсона-Кука T, edot: размер, обнуление, вывод в конец строки ---
    {
        sys.resize(N);
        if (sys.T.size() != N || sys.edot.size() != N) { std::printf("FAIL resize T/edot\n"); ++failures; }
        sys.T[3] = 123.0; sys.edot[3] = 4.5;
        sys.inzializ_index(3);
        if (sys.T[3] != 0.0 || sys.edot[3] != 0.0) { std::printf("FAIL obnul T/edot\n"); ++failures; }
        sys.T[0] = 300.0; sys.edot[0] = 7.0;
        std::ofstream tmp("__ps_test_out2.txt");
        sys.vuvod_index(0, tmp);
        tmp.close();
        std::ifstream chk("__ps_test_out2.txt");
        std::string line;
        std::getline(chk, line);
        std::istringstream ss(line);
        std::vector<double> v;
        double x;
        while (ss >> x) v.push_back(x);
        // itype mass x[dim] vx[dim] + 8 скаляров + hsml + 6 компонент s + T edot
        const size_t expect = 2 + 2 * dim + 5 + 3 + 6 + 2;
        if (v.size() != expect) { std::printf("FAIL vuvod столбцов %zu != %zu\n", v.size(), expect); ++failures; }
        else if (v[expect - 2] != 300.0 || v[expect - 1] != 7.0) { std::printf("FAIL vuvod T/edot в конце\n"); ++failures; }
        std::remove("__ps_test_out2.txt");
    }

    // --- inzializ_2d_symm_index: корректность зеркалирования (проверяем поля, общие для любого dim) ---
    // Для dim>=1: x[0] отражается, x[1..] копируется. itype/mass/rho копируются.
    if (dim >= 1)
    {
        sys.resize(2);
        sys.inzializ_index(0);
        sys.inzializ_index(1);
        sys.x[0][0] = 5.0;   // src
        sys.mass[0] = 2.5;
        sys.itype[0] = 4;
        sys.T[0] = 321.0; sys.edot[0] = 9.0;
        sys.inzializ_2d_symm_index(1, 0);  // dst=1 = образ src=0
        // x[0] отражается: -5.0
        if (std::fabs(sys.x[0][1] - (-5.0)) > 1e-12) { std::printf("FAIL inzializ_2d_symm x[0] mirror: %g\n", sys.x[0][1]); ++failures; }
        if (std::fabs(sys.mass[1] - 2.5) > 1e-12)    { std::printf("FAIL inzializ_2d_symm mass copy\n"); ++failures; }
        if (sys.T[1] != 321.0 || sys.edot[1] != 9.0) { std::printf("FAIL inzializ_2d_symm T/edot copy\n"); ++failures; }
        if (sys.itype[1] != 4)                       { std::printf("FAIL inzializ_2d_symm itype copy\n"); ++failures; }
    }

    if (failures == 0)
    {
        std::printf("particle_system_test: OK (size=%zu, dim=%zu)\n", N, dim);
        return 0;
    }
    std::printf("particle_system_test: FAILED (%d ошибок)\n", failures);
    return 1;
}
