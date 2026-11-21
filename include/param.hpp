#ifndef PARAMS_HPP
#define PARAMS_HPP
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstddef>
// Размерность задачи (1, 2 или 3).
// Значение задаётся при сборке макросом SPH_DIM (решение E1,
// notes/contact/04-e2e-test-brief.md §2); без него — 1, как раньше. Семантика
// не меняется: dim остаётся constexpr и определяет размер C-массивов
// x[dim]/vx[dim] на этапе компиляции (ADR-002).
#ifndef SPH_DIM
#define SPH_DIM 1
#endif
constexpr size_t dim = SPH_DIM;
// Общие параметры
struct Particle
{
    Particle();
    std::vector<double> x;
    std::vector<double> vx;
    int itype;
    double mass;
    double rho;
    double rho_0;
    double p;
    double u;
    double c;
    double s;
    double e;
    //--Тензор напряжений
    double sxx;
    double syy;
    double szz;
    double sxy;
    double syz;
    double sxz;
    //--интенсивность пластических деформаций
    double e_pl;
    //--Тензор скоростей деформаций
    double exx;
    double eyy;
    double ezz;
    double exy;
    double exz;
    double eyz;
    //-- Тензор скоростей поворота
    double rxy;
    double rxz;
    double ryz;
    //--Параметр поврежденности
    double D;
    // Для модели JH-2
    double dp;
    double hsml;
    void inzializ();
    void inzializ_2d_symm(Particle buf);
    void obnul();
    void vuvod(std::ofstream& filename);
    //--
    bool axis;
};
// определения методов Particle помечены `inline`, чтобы param.hpp
// оставался header-only и безопасно включался в несколько единиц трансляции
// (param.hpp — компиляционный инвариант: dim, флаги, тип Particle).
inline Particle::Particle()
{
    for (size_t i = 0; i < dim; i++)
    {
        x.push_back(0.0);
        vx.push_back(0.0);
    }
    itype = 0;
    mass = 0.0;
    rho = 0.0;
    rho_0 = 0.0;
    p = 0.0;
    u = 0.0;
    c = 0.0;
    s = 0.0;
    e = 0.0;
    sxx = 0.0;
    syy = 0.0;
    szz = 0.0;
    sxy = 0.0;
    syz = 0.0;
    sxz = 0.0;
    exx = 0.0;
    eyy = 0.0;
    ezz = 0.0;
    exy = 0.0;
    eyz = 0.0;
    exz = 0.0;
    rxy = 0.0;
    rxz = 0.0;
    ryz = 0.0;
    hsml = 0.0;
    D = 0.0;
    dp = 0.0;
    e_pl = 0.0;
    axis = false;
}
inline void Particle::inzializ()
{
    for (size_t i = 0; i < dim; i++)
    {
        x.push_back(0.0);
        vx.push_back(0.0);
    }
    itype = 0;
    mass = 0.0;
    rho = 0.0;
    rho_0 = 0.0;
    p = 0.0;
    u = 0.0;
    c = 0.0;
    s = 0.0;
    e = 0.0;
    sxx = 0.0;
    syy = 0.0;
    szz = 0.0;
    sxy = 0.0;
    syz = 0.0;
    sxz = 0.0;
    exx = 0.0;
    eyy = 0.0;
    ezz = 0.0;
    exy = 0.0;
    eyz = 0.0;
    exz = 0.0;
    rxy = 0.0;
    rxz = 0.0;
    ryz = 0.0;
    hsml = 0.0;
    D = 0.0;
    dp = 0.0;
    e_pl = 0.0;
    axis = false;
}
inline void Particle::inzializ_2d_symm(Particle buf)
{
    for (size_t i = 0; i < dim; i++)
    {
        if (i == 0)
        {
            x[i] = -buf.x[i];
            vx[i] = -buf.vx[i];
        }
        else
        {
            x[i] = buf.x[i];
            vx[i] = buf.vx[i];
        }
    }
    sxy = -buf.sxy;
    exy = buf.exy;
    rxy = buf.rxy;
    itype = buf.itype;
    mass = buf.mass;
    rho = buf.rho;
    rho_0 = buf.rho_0;
    p = buf.p;
    u = buf.u;
    c = buf.c;
    s = buf.s;
    e = buf.e;
    sxx = buf.sxx;
    syy = buf.syy;
    szz = buf.szz;
    syz = buf.syz;
    sxz = buf.sxz;

    exx = buf.exx;
    eyy = buf.eyy;
    ezz = buf.ezz;
    
    eyz = buf.eyz;
    exz = buf.exz;

    
    rxz = buf.rxz;
    ryz = buf.ryz;
    D = buf.D;
    dp = buf.dp;
    hsml = buf.hsml;
    e_pl = buf.e_pl;
}
inline void Particle::obnul()
{
    for (size_t i = 0; i < dim; i++)
    {
        x[i] = 0.0;
        vx[i] = 0.0;
    }
    itype = 0;
    mass = 0.0;
    rho = 0.0;
    rho_0 = 0.0;
    p = 0.0;
    u = 0.0;
    c = 0.0;
    s = 0.0;
    e = 0.0;
    sxx = 0.0;
    syy = 0.0;
    szz = 0.0;
    sxy = 0.0;
    syz = 0.0;
    sxz = 0.0;
    exx = 0.0;
    eyy = 0.0;
    ezz = 0.0;
    exy = 0.0;
    eyz = 0.0;
    exz = 0.0;
    rxy = 0.0;
    rxz = 0.0;
    ryz = 0.0;
    hsml = 0.0;
    D = 0.0;
    dp = 0.0;
    e_pl = 0.0;
}
inline void Particle::vuvod(std::ofstream& filename)
{
    filename << " " << itype;
    filename << " " << mass;

    // Координаты
    for (size_t d = 0; d < dim; ++d)
    {
        filename << " " << x[d];
    }
    // Скорости
    for (size_t d = 0; d < dim; ++d)
    {
        filename << " " << vx[d];
    }
    // Остальные параметры
    double sigma_i = sqrt((sxx - syy) * (sxx - syy) + (sxx - szz) * (sxx - szz) + (syy - szz) * (syy - szz) + 6 * (sxy * sxy + syz * syz + sxz * sxz))/sqrt(2.0);
    filename << " " << rho << " " << p << " " << u << " " << c << " " << sigma_i << " " << e_pl << " " << D << " " << hsml 
             << " " << sxx << " " << syy << " " << szz << " " << sxy << " " << syz << " " << sxz <<"\n";
}
// Параметры геометрии
constexpr double x_maxgeom = 10.0;
constexpr double x_mingeom = -10.0;
constexpr double y_maxgeom = 10.0;
constexpr double y_mingeom = -10.0;
constexpr double z_maxgeom = 10.0;
constexpr double z_mingeom = -10.0;

// Параметры SPH-алгоритма (pa_sph, nnps, sle, skf, density_method,
// coloring_niac_threshold, CFL, epsilon, alpha, beta, etq, xi, epsilon_h)
// перенесены в struct SimulationConfig (simulation_config.hpp), подэтап 5.7.
// Логические флаги
constexpr bool summation_density = false; // Использовать суммирование плотности
constexpr bool average_velocity = true;   // Использовать усреднение скорости
constexpr bool config_input = false;      // Загрузка начальной конфигурации
constexpr bool virtual_part = true;       // Использовать виртуальные частицы
constexpr bool vp_input = false;          // Загрузка виртуальных частиц
constexpr bool visc = false;              // Физическая вязкость
constexpr bool ex_force = true;           // Внешние силы
constexpr bool visc_artificial = true;    // Искусственная вязкость
constexpr bool heat_artificial = false;   // Искусственный нагрев
constexpr bool self_gravity = false;      // Самогравитация
constexpr bool nor_density = false;       // Нормализация плотности
constexpr bool axis_sim = true;           // осесимметричные расчеты для 2D

constexpr bool is_contact_pair = true;
// Параметры симметрии
constexpr int nsym = 0; // Тип симметрии (0: нет, 1: осевая, 2: центральная) для одномерных расчетов

// Параметры вывода
constexpr bool int_stat = true; // Вывод статистики взаимодействия
constexpr int save_step = 4;    // Шаги для сохранения на диск
// print_step, save_t, t_total, moni_particle перенесены в SimulationConfig.

// Математические константы
constexpr double pi = 3.14159265358979323846; // Число Пи

// Режимы работы
constexpr bool shocktube = true;    // Тест "ударная трубка"
constexpr bool shearcavity = false; // Тест "сдвиговая каверна"
// Материальные константы (G, sigy, aa, zz, gamm, const_JH_2, n_Me) перенесены
// в struct Materials (materials.hpp) — подобъект SimulationConfig (ADR-3 этапа 6.3).

#endif // PARAMS_HPP
