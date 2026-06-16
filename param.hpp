#include <iostream>
#include <fstream>
#include <sstream>
#ifndef PARAMS_HPP
#define PARAMS_HPP
constexpr size_t dim = 1; // Размерность задачи (1, 2 или 3)
// Общие параметры
struct Particle
{
    Particle::Particle();
    std::vector<double> x;
    std::vector<double> vx;
    size_t itype;
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
Particle::Particle()
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
void Particle::inzializ()
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
void Particle::inzializ_2d_symm(Particle buf)
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
void Particle::obnul()
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
void Particle::vuvod(std::ofstream& filename)
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
int maxn = 5700000;              // Максимальное число частиц
int max_interaction = 100 * maxn; // Максимальное число взаимодействующих пар

// Параметры геометрии
constexpr double x_maxgeom = 10.0;
constexpr double x_mingeom = -10.0;
constexpr double y_maxgeom = 10.0;
constexpr double y_mingeom = -10.0;
constexpr double z_maxgeom = 10.0;
constexpr double z_mingeom = -10.0;

// Параметры SPH-алгоритма
int pa_sph = 2;       // Метод аппроксимации давления
int nnps = 1;         // Метод поиска соседних частиц
int sle = 1;          // Алгоритм обновления длины сглаживания
int skf = 1;          // Тип ядра сглаживания
double CFL = 0.1;     // Число Куранта
double epsilon = 0.3; // Константа для усреднения скоростей (Monaghan 1992)
// Псевдовязкость
double alpha = 1.0;
double beta = 1.0;
double etq = 0.01;
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
int print_step = 50;   // Шаги для вывода на экран
constexpr int save_step = 4;    // Шаги для сохранения на диск
double save_t = 2.5e-8;
double t_total = 3e-6;
size_t moni_particle = 2403366; // Номер мониторимой частицы

// Математические константы
constexpr double pi = 3.14159265358979323846; // Число Пи

// Режимы работы
constexpr bool shocktube = true;    // Тест "ударная трубка"
constexpr bool shearcavity = false; // Тест "сдвиговая каверна"
// Формальная нумерация материалов
// 3 - тантал, 4 - алюминий, 5 - LiF
// 6 - полиуретан 160, 7 - полиуретан 320, 8 - полиуретан 626
// 9 - PMMA, 10 - B4C, 11 - Al2O3
constexpr size_t n_Me = 3;
// Модуль сдвига
std::vector<double> G{69e9, 26e9, 5.5e+10, 5.27e7, 1.836e8, 4e+08, 3e+09};
// Предел текучести для модели Мизеса
std::vector<double> sigy{450e6, 240e6, 2e+07, 1e+06, 3.0e+6, 4e+06, 7.7e+07};
// Константы для УА Тэта и УРС Ми-Грюнайзена
std::vector<double> aa{1.4506e11, 1.97e+10, 1.6659e+10, 2.1953e8, 3.595e8,3.0227e8, 1.8797e9};
std::vector<double> zz{1.987, 4.2, 4.272, 1.245, 1.736, 3.5, 4.668};
std::vector<double> gamm{2.6, 1.97};
// Наборы констант для модели керамики JH-2; 0 - B4C, 1 - Al2O3
// 0 - плотность, 1 - модуль сдвига, 2 - P_hel, 3 - s_hel, 4 - T_n, 5 - A, 6 - N, 7 - B, 8 - M
// 9 - s_f_max, 10 - C, 11 - D1, 12 - D2, 13 - K1, 14 - K2, 15 - K3, 16 - T
std::vector<std::vector<double>> const_JH_2{
    {2510, 197e9, 8.71e9, 15.44e9, 0.03, 0.927, 0.67, 0.70, 0.85, 0.2, 0.005, 0.001, 0.5, 233e9, -593e9, 2800e9, 0.26e9},
    {3555, 1.165e11, 2.626e9, 3.561e9, 0.137, 0.930, 0.600, 0.110, 0.600, 1.0, 0.0, 0.005, 1.000, 1.692e11, 0.0, 0.0, 0.2e9}};

#endif // PARAMS_HPP
