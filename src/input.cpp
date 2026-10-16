// input.cpp
// Определения функций ввода/инициализации данных.
// Объявления — в input.hpp. Вспомогательные функции разбора (shock_tube,
// shear_cavity) — внутренние (internal linkage / static) этого TU.
#include "input.hpp"
#include "fiber/fiber_system.hpp"
#include "fiber/fiber_input.hpp"
#include "k_input/k_geometry.hpp"
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <unordered_map>
#include <array>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

// Вспомогательная функция для форматированного вывода
static void format_output(std::ofstream& file, size_t particle_id, const std::vector<double>& values)
{
    file << std::setw(5) << particle_id;
    for (double val : values)
    {
        file << std::setw(8) << std::scientific << val;
    }
    file << "\n";
}

// Вспомогательная функция для форматированного вывода целых чисел
static void format_output_other(std::ofstream& file, size_t particle_id, int itype_val, double hsml_val)
{
    file << std::setw(5) << particle_id << "  "
         << std::setw(2) << itype_val << "  "
         << std::setw(8) << std::scientific << hsml_val << "\n";
}
static void trim_spaces(std::string& s)
{
    // Приведение к unsigned char обязательно. std::isspace принимает int, но
    // допустимы только значения, представимые в unsigned char, и EOF; передача
    // отрицательного — неопределённое поведение. На char со знаком (x86) любой
    // байт UTF-8 со старшим битом даёт отрицательное значение, а в
    // Initial_Conditions.txt комментарии русские, и строка, кончающаяся
    // кириллицей, приводит сюда байт вроде 0xB9 = -71.
    // На Windows с UCRT это срабатывание внутреннего assert и немедленный
    // abort: процесс умирает молча, до единой строки вывода, ещё на разборе
    // конфигурации. Так падали ОБА сквозных теста, BC2 и BC6.
    auto not_space = [](unsigned char ch){ return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
}
// Положительное конечное число для ключей k_*; иначе исключение с именем ключа.
static double parse_k_positive(const char* key, const std::string& value)
{
    std::size_t pos = 0;
    double v = 0.0;
    try
    {
        v = std::stod(value, &pos);
    }
    catch (const std::exception&)
    {
        pos = 0;
    }
    if (pos == 0 || pos != value.size() || !std::isfinite(v) || !(v > 0.0))
    {
        throw std::runtime_error(std::string(key) + ": ожидается конечное число > 0, получено '" + value + "'");
    }
    return v;
}
// Закон нормальной силы контакта нитей (fiber_*_contact_type): 0..3,
// 0 = Герц (НЕ «выключено», в отличие от contact_force_type). Вне диапазона
// -- сообщение и откат к 0.
static int parse_fiber_contact_type(const std::string& key, const std::string& value)
{
    const int v = std::stoi(value);
    if (v < 0 || v > 3)
    {
        std::cerr << key << " вне диапазона 0..3: " << v << ", откат к 0 (Герц)\n";
        return 0;
    }
    return v;
}

// gamma типа 3 контакта нитей: те же границы и сообщения, что у
// contact_impulse_gamma.
static double parse_fiber_contact_gamma(const std::string& key, const std::string& value)
{
    const double v = std::stod(value);
    if (!(v > 0.0))
    {
        std::cerr << key << " = " << v << " <= 0: тип 3 фактически выключен\n";
    }
    else if (v > 1.6)
    {
        std::cerr << key << " = " << v
                  << " > 1.6: включается ограничитель шага контакта,"
                     " возможен коллапс dt\n";
    }
    return v;
}

void parse_sph_config(SimulationConfig& cfg, const std::string& filename)
{
    std::ifstream in(filename);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(in, line))
    {
        ++line_number;
        trim_spaces(line);
        if (line.empty())
        {
            continue;
        }
        // Строки-комментарии начинаются с '*'
        if (line[0] == '*')
        {
            continue;
        }
        std::size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos)
        {
            continue;
        }        
        std::string key = line.substr(0, eq_pos);
        std::string value = line.substr(eq_pos + 1);
        trim_spaces(key);
        trim_spaces(value);
        if (key.empty() || value.empty())
        {
            continue;            
        }
        if (key == "pa_sph")
        {
            cfg.pa_sph = std::stoi(value);
        }
        else if (key == "nnps")
        {
            cfg.nnps = std::stoi(value);
        }
        else if (key == "sle")
        {
            cfg.sle = std::stoi(value);
        }
        else if (key == "skf")
        {
            cfg.skf = std::stoi(value);
        }
        else if (key == "CFL")
        {
            cfg.CFL = std::stod(value);
        }
        else if (key == "epsilon")
        {
            cfg.epsilon = std::stod(value);
        }
        else if (key == "alpha")
        {
            cfg.alpha = std::stod(value);
        }
        else if (key == "beta")
        {
            cfg.beta = std::stod(value);
        }
        else if (key == "etq")
        {
            cfg.etq = std::stod(value);
        }
        else if (key == "density_method")
        {
            cfg.density_method = std::stoi(value);
        }
        else if (key == "coloring_niac_threshold")
        {
            cfg.coloring_niac_threshold = std::stoi(value);
        }
        else if (key == "t_total")
        {
            cfg.t_total = std::stod(value);
        }
        else if (key == "save_t")
        {
            cfg.save_t = std::stod(value);
        }
        else if (key == "print_step")
        {
            cfg.print_step = static_cast<size_t>(std::stoull(value));
        }
        else if (key == "contact_force_type")
        {
            const int v = std::stoi(value);
            if (v < 0 || v > 3)
            {
                std::cerr << "contact_force_type вне диапазона 0..3: " << v
                          << ", откат к 0 (контакт выключен)\n";
                cfg.contact_force_type = 0;
            }
            else
            {
                cfg.contact_force_type = v;
            }
        }
        else if (key == "contact_penalty_scale")
        {
            cfg.contact_penalty_scale = std::stod(value);
        }
        else if (key == "contact_impulse_gamma")
        {
            const double v = std::stod(value);
            // Границы: <= 0 равносильно выключенному контакту типа 3 (ранний
            // выход в impulse_force), > 1.6 включает ограничитель шага, который
            // сам зависит от урезанного шага (петля «жёсткость -> шаг»).
            // Значение принимается как есть, но молча это не проходит.
            if (!(v > 0.0))
            {
                std::cerr << "contact_impulse_gamma = " << v
                          << " <= 0: тип 3 фактически выключен\n";
            }
            else if (v > 1.6)
            {
                std::cerr << "contact_impulse_gamma = " << v
                          << " > 1.6: включается ограничитель шага контакта,"
                             " возможен коллапс dt\n";
            }
            cfg.contact_impulse_gamma = v;
        }
        else if (key == "fiber_solid_contact_type")
        {
            cfg.fiber_solid_contact_type = parse_fiber_contact_type(key, value);
        }
        else if (key == "fiber_solid_contact_penalty_scale")
        {
            cfg.fiber_solid_contact_penalty_scale = std::stod(value);
        }
        else if (key == "fiber_solid_contact_impulse_gamma")
        {
            cfg.fiber_solid_contact_impulse_gamma = parse_fiber_contact_gamma(key, value);
        }
        else if (key == "fiber_fiber_contact_type")
        {
            cfg.fiber_fiber_contact_type = parse_fiber_contact_type(key, value);
        }
        else if (key == "fiber_fiber_contact_penalty_scale")
        {
            cfg.fiber_fiber_contact_penalty_scale = std::stod(value);
        }
        else if (key == "fiber_fiber_contact_impulse_gamma")
        {
            cfg.fiber_fiber_contact_impulse_gamma = parse_fiber_contact_gamma(key, value);
        }
        else if (key == "ode_solver")
        {
            const int v = std::stoi(value);
            if (v != 0 && v != 1)
            {
                std::cerr << "ode_solver вне диапазона 0..1: " << v
                          << ", откат к 0 (PredictorCorrector)\n";
                cfg.ode_solver = 0;
            }
            else
            {
                cfg.ode_solver = v;
            }
        }
        else if (key == "sph_same_body_only")
        {
            // Принимаются 0/1 и true/false (без учёта регистра); иное -- откат к false.
            std::string low = value;
            for (char& ch : low) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
            if (low == "1" || low == "true")
            {
                cfg.sph_same_body_only = true;
            }
            else if (low == "0" || low == "false")
            {
                cfg.sph_same_body_only = false;
            }
            else
            {
                std::cerr << "sph_same_body_only: недопустимое значение '" << value
                          << "' (ожидается 0/1/true/false), откат к 0\n";
                cfg.sph_same_body_only = false;
            }
        }
        else if (key == "fiber_file")
        {
            cfg.fiber_file = value;
            // fiber/PLAN.md, этап Ф2, пункт 5: если нити реально заданы, а
            // dim != 3, это ошибка, и она должна выявляться здесь, на этапе
            // разбора конфигурации, а не позже в середине расчёта.
            fiber::check_fiber_requires_dim3();
        }
        else if (key == "k_file")
        {
            cfg.k_file = value;
        }
        else if (key == "k_mode")
        {
            // Геометрия важнее удобства: неверное значение -- ошибка, не откат.
            if (value != "solid" && value != "solid_fiber")
            {
                throw std::runtime_error("k_mode: допустимы solid и solid_fiber, получено '" + value + "'");
            }
            cfg.k_mode = value;
            if (value == "solid_fiber")
            {
                fiber::check_fiber_requires_dim3();
            }
        }
        else if (key == "k_parts_file")
        {
            cfg.k_parts_file = value;
        }
        else if (key == "k_length_scale")
        {
            cfg.k_length_scale = parse_k_positive("k_length_scale", value);
        }
        else if (key == "k_mass_scale")
        {
            cfg.k_mass_scale = parse_k_positive("k_mass_scale", value);
        }
        else if (key == "k_monitor_nid")
        {
            std::size_t pos = 0;
            long long v = 0;
            try
            {
                v = std::stoll(value, &pos);
            }
            catch (const std::exception&)
            {
                pos = 0;
            }
            if (pos == 0 || pos != value.size())
            {
                throw std::runtime_error("k_monitor_nid: ожидается целое число, получено '" + value + "'");
            }
            cfg.k_monitor_nid = v;
        }
        else if (key == "k_solid_elements")
        {
            if (value != "skip" && value != "particles")
            {
                throw std::runtime_error("k_solid_elements: допустимы skip и particles, получено '" + value + "'");
            }
            cfg.k_solid_elements = value;
        }
    }
    in.close();
}
// Функция инициализации частиц для теста "ударная трубка"
void shock_tube(SimulationConfig& cfg, ParticleSystem& sys, size_t& ntotal, double& delta_x)
{
    /*ntotal = 3200;
    constexpr size_t n_l_i = 2560;
    size_t n_p_i = ntotal - n_l_i;
    double delta_x = 0.6 / double (n_p_i);

    // Инициализация всех частиц
    for (size_t i = 0; i < ntotal; ++i)
    {
        mass[i] = 0.75 / double(ntotal);
        hsml[i] = 2 * delta_x;
        itype[i] = 1;  // Идеальный газ

        for (size_t d = 0; d < dim; ++d)
        {
            x[d][i] = 0.0;
            vx[d][i] = 0.0;
        }
    }
    double dx_l = delta_x / 4.0;
    // Установка координат для левой части
    for (size_t i = 0; i < n_l_i; ++i)
    {
        x[0][i] = -0.6 + dx_l * double(i);
    }

    // Установка координат для правой части
    for (size_t i = n_l_i; i < ntotal; ++i)
    {
        x[0][i] = delta_x * (double(i) - (n_l_i - 1));
    }

    // Установка начальных параметров
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (x[0][i] <= 1e-8)
        {
            // Левое состояние
            u[i] = 2.5;
            rho[i] = 1.0;
            p[i] = 1.0;
        }
        else
        {
            // Правое состояние
            u[i] = 1.795;
            rho[i] = 0.25;
            p[i] = 0.1795;
        }
    } */
    ntotal = 5000;
    // Обеспечиваем размер контейнера (при старте он может быть пустым)
    if (sys.size() < ntotal)
    {
        sys.resize(ntotal);
    }
    //double R_ud = 0.05;
    double Z1 = 0.006;
    double Z2 = Z1 + 0.004831;
    double Z3 = Z2 + 0.010322;
    double Z4 = Z3 + 0.0254;
    double L = Z4;
    double vx0 = 1633;
    delta_x = L / double (ntotal);
    bool flag_moni_particle = false;
    for (size_t i = 0; i < ntotal; ++i)
    {
        sys.hsml[i] = 2 * delta_x;
        sys.u[i] = 0.0;
        sys.p[i] = 0.0;
        sys.x[0][i] = delta_x * (double(i) + 1);
        if (sys.x[0][i] <= Z1)
        {
            sys.itype[i] = 7;
            sys.rho[i] = 320;
            sys.rho_0[i] = sys.rho[i];
            sys.mass[i] = sys.rho[i] * delta_x;
            sys.vx[0][i] = vx0;

        }
        else if (sys.x[0][i] <= Z2)
        {
            sys.itype[i] = 10;
            sys.rho[i] = 2510;
            sys.rho_0[i] = sys.rho[i];
            sys.mass[i] = sys.rho[i] * delta_x;
            sys.vx[0][i] = vx0;
        }
        else if (sys.x[0][i] <= Z3)
        {
            sys.itype[i] = 10;
            sys.rho[i] = 2510;
            sys.rho_0[i] = sys.rho[i];
            sys.mass[i] = sys.rho[i] * delta_x;
            sys.vx[0][i] = 0.0;
            sys.x[0][i] += delta_x; //=======
        }
        else // последняя область: любой x > Z3 (в т.ч. x > Z4 из-за округления) → LiF
             // ВАЖНО: здесь именно `else`, а не `else if (x <= Z4)`. Из-за неточности
             // плавающей арифметики delta_x*ntotal (= x последней частицы) может
             // оказаться чуть больше Z4, и тогда частица "выпадала" из всех веток,
             // оставаясь с rho=0, mass=0, itype=0. Это порождало NaN/Inf в dvxdt,
             // dedt и hsml (деление 0/0 и на нулевую плотность) для последних частиц.
        {
            if (flag_moni_particle == false)
            {
                flag_moni_particle = true;
                cfg.moni_particle = i - 1;
            }
            sys.itype[i] = 5;
            sys.rho[i] = 2638;
            sys.rho_0[i] = sys.rho[i];
            sys.mass[i] = sys.rho[i] * delta_x;
            sys.vx[0][i] = 0.0;
            sys.x[0][i] += delta_x; //=======
        }
    }

    // body[] = itype[] (fiber/PLAN.md §3.1): встроенные генераторы не знают
    // отдельного номера тела, поэтому body совпадает с материалом — это
    // сохраняет нынешний критерий Solid-Solid-Contact «разные материалы»
    // (fiber/SPH.md, раздел «Solid -- Solid») без изменений на BC2/BC6.
    for (size_t i = 0; i < ntotal; ++i)
    {
        sys.body[i] = sys.itype[i];
    }

    // Контроль корректности инициализации: каждая частица обязана получить
    // материал (itype > 0) и положительные плотность/массу. Это страховка от
    // ошибок граничных сравнений (см. комментарий в ветке `else` выше): если
    // хоть одна частица останется с itype=0/rho=0, она даст деление на ноль и
    // каскад NaN по всей области взаимодействия. Лучше упасть здесь с понятным
    // сообщением, чем ловить -nan(ind) глубоко в силовом расчёте.
    for (size_t i = 0; i < ntotal; ++i)
    {
        if (sys.itype[i] == 0 || sys.rho[i] <= 0.0 || sys.mass[i] <= 0.0)
        {
            throw std::runtime_error(
                "shock_tube: частица " + std::to_string(i) +
                " не получила корректный материал (itype=" + std::to_string(sys.itype[i]) +
                ", rho=" + std::to_string(sys.rho[i]) +
                ", mass=" + std::to_string(sys.mass[i]) + ")");
        }
    }
    /*ntotal = 0;
    delta_x = 0.0002;
    double Z1 = 0.006;
    double Z2 = Z1 + 0.003917;
    double Z3 = Z2 + 0.009033;
    double Z4 = Z3 + 0.0254;
    //for (double r = -0.04375; r <= 0.04375; r += delta_x)
    for (double r = 0.5 * delta_x; r <= 0.04375; r += delta_x)
    {
        for (double z = 0.0; z <= 0.006; z += delta_x)
        {
            if (fabs (r) <= 0.5 * delta_x)
            {
                sys.axis[ntotal] = true;
            }
            sys.hsml[ntotal] = 2.0 * delta_x;
            sys.u[ntotal] = 0.0;
            sys.p[ntotal] = 0.0;
            sys.x[0][ntotal] = r;
            sys.x[1][ntotal] = z;
            sys.itype[ntotal] = 8;
            sys.rho[ntotal] = 626;
            sys.rho_0[ntotal] = sys.rho[ntotal];
            sys.vx[0][ntotal] = 0.0;
            sys.vx[1][ntotal] = 2210;
            sys.mass[ntotal] = sys.rho[ntotal] * delta_x * delta_x;
            // sys.mass[ntotal] = sys.rho[ntotal] * 2 * pi * pi * 0.25 * delta_x * delta_x * sys.x[0][ntotal];
            // sys.mass[ntotal] = 2 * pi * sys.x[0][ntotal] * delta_x * delta_x * sys.rho[ntotal];
            ntotal++;
        }
    }
    //for (double r = -0.0346; r <= 0.0346; r += delta_x)
    for (double r = 0.5 * delta_x; r <= 0.0346; r += delta_x)
    {
        for (double z = Z1 + delta_x; z <= Z2; z += delta_x)
        {
            if (fabs (r) <= 0.5 * delta_x)
            {
                sys.axis[ntotal] = true;
            }
            sys.hsml[ntotal] = 2.0 * delta_x;
            sys.u[ntotal] = 0.0;
            sys.p[ntotal] = 0.0;
            sys.x[0][ntotal] = r;
            sys.x[1][ntotal] = z;
            sys.itype[ntotal] = 10;
            sys.rho[ntotal] = 2510;
            sys.rho_0[ntotal] = sys.rho[ntotal];
            sys.vx[0][ntotal] = 0.0;
            sys.vx[1][ntotal] = 2210;
            sys.mass[ntotal] = sys.rho[ntotal] * delta_x * delta_x;
            // sys.mass[ntotal] = sys.rho[ntotal] * 2 * pi * pi * 0.25 * delta_x * delta_x * sys.x[0][ntotal];
            // sys.mass[ntotal] = 2 * pi * sys.x[0][ntotal] * delta_x * delta_x * sys.rho[ntotal];
            ntotal++;
        }
    }
    //for (double r = -0.0346; r <= 0.0346; r += delta_x)
    for (double r = 0.5 * delta_x; r <= 0.0346; r += delta_x)
    {
        for (double z = Z2 + 4 * delta_x; z <= Z3 + 4 * delta_x; z += delta_x)
        {
            if (fabs (r) <= 0.5 * delta_x)
            {
                sys.axis[ntotal] = true;
            }
            sys.hsml[ntotal] = 2.0 * delta_x;
            sys.u[ntotal] = 0.0;
            sys.p[ntotal] = 0.0;
            sys.x[0][ntotal] = r;
            sys.x[1][ntotal] = z;
            sys.itype[ntotal] = 10;
            sys.rho[ntotal] = 2510;
            sys.rho_0[ntotal] = sys.rho[ntotal];
            sys.vx[0][ntotal] = 0.0;
            sys.vx[1][ntotal] = 0.0;
            sys.mass[ntotal] = sys.rho[ntotal] * delta_x * delta_x;
            // sys.mass[ntotal] = sys.rho[ntotal] * 2 * pi * pi * 0.25 * delta_x * delta_x * sys.x[0][ntotal];
            // sys.mass[ntotal] = 2 * pi * sys.x[0][ntotal] * delta_x * delta_x * sys.rho[ntotal];
            ntotal++;
        }
    }
    moni_particle = 12973;
    //for (double r = -0.0254; r <= 0.0254; r += delta_x)
    for (double r = 0.5 * delta_x; r <= 0.0254; r += delta_x)
    {
        for (double z = Z3 + 5 * delta_x; z <= Z4 + 4 * delta_x; z += delta_x)
        {
            if (fabs (r) <= 0.5 * delta_x)
            {
                sys.axis[ntotal] = true;
            }
            sys.hsml[ntotal] = 2.0 * delta_x;
            sys.u[ntotal] = 0.0;
            sys.p[ntotal] = 0.0;
            sys.x[0][ntotal] = r;
            sys.x[1][ntotal] = z;
            sys.itype[ntotal] = 5;
            sys.rho[ntotal] = 2638;
            sys.rho_0[ntotal] = sys.rho[ntotal];
            sys.vx[0][ntotal] = 0.0;
            sys.vx[1][ntotal] = 0.0;
            sys.mass[ntotal] = sys.rho[ntotal] * delta_x * delta_x;
            // sys.mass[ntotal] = sys.rho[ntotal] * 2 * pi * pi * 0.25 * delta_x * delta_x * sys.x[0][ntotal];
            // sys.mass[ntotal] = 2 * pi * sys.x[0][ntotal] * delta_x * delta_x * sys.rho[ntotal];
            ntotal++;
        }
    }*/
}

// Функция инициализации частиц для теста "сдвиговая каверна"
void shear_cavity(ParticleSystem& sys, size_t& ntotal)
{
    const size_t m = 41;
    const size_t n = 41;
    const size_t mp = m - 1;
    const size_t np = n - 1;
    ntotal = mp * np;
    // Обеспечиваем размер контейнера (при старте он может быть пустым)
    if (sys.size() < ntotal)
    {
        sys.resize(ntotal);
    }

    const double xl = 1e-3;
    const double yl = 1e-3;
    const double dx = xl / mp;
    const double dy = yl / np;

    // Инициализация координат
    for (size_t i = 0; i < mp; ++i)
    {
        for (size_t j = 0; j < np; ++j)
        {
            size_t k = j + i * np;
            sys.x[0][k] = (double(i) * dx) + dx / 2.0;
            sys.x[1][k] = (double(j) * dy) + dy / 2.0;
        }
    }

    // Инициализация остальных параметров
    for (size_t i = 0; i < ntotal; ++i)
    {
        sys.vx[0][i] = 0.0;
        sys.vx[1][i] = 0.0;
        sys.rho[i] = 1000.0;  // Плотность воды
        sys.mass[i] = dx * dy * sys.rho[i];
        sys.p[i] = 0.0;
        sys.u[i] = 357.1;     // Внутренняя энергия
        sys.itype[i] = 2;     // Водяная частица
        sys.hsml[i] = dx;     // Длина сглаживания = начальный шаг
        sys.body[i] = sys.itype[i];   // body[] = itype[] (fiber/PLAN.md §3.1)
    }
}

// Встроенные источники Solid: particles.dat / shock_tube / shear_cavity.
static void input_builtin(SimulationConfig& cfg, ParticleSystem& sys, size_t& ntotal, double& dx)
{
    // Если загрузка из файла включена
    if (config_input)
    {
        /*std::ifstream xv_file("../data/f_xv.dat");
        std::ifstream state_file("../data/f_state.dat");
        std::ifstream other_file("../data/f_other.dat");

        if (!xv_file.is_open() || !state_file.is_open() || !other_file.is_open())
        {
            throw std::runtime_error("Не удалось открыть файлы для загрузки начальной конфигурации");
        }

        xv_file >> ntotal;
        std::wcout << " ******************************\n";
        std::wcout << " Загрузка начальной конфигурации...\n";
        std::wcout << " Общее количество частиц " << ntotal << "\n";
        std::wcout << " ******************************\n";

        for (size_t i = 0; i < ntotal; ++i)
        {
            int im;
            xv_file >> im;
            for (size_t d = 0; d < dim; ++d) xv_file >> x[d][i];
            for (size_t d = 0; d < dim; ++d) xv_file >> vx[d][i];

            state_file >> im >> mass[i] >> rho[i] >> p[i] >> u[i];
            other_file >> im >> itype[i] >> hsml[i];
        }

        xv_file.close();
        state_file.close();
        other_file.close();
        */
        std::ifstream file("../data/particles.dat");
        if (!file.is_open())
        {
            throw std::runtime_error("Не удалось открыть файл ../data/particles.dat");
        }

        // Чтение общего количества частиц
        file >> ntotal;
        file.ignore((std::numeric_limits<std::streamsize>::max)(), '\n'); // Переход на новую строку

        // Вывод информации о загрузке
        std::cout << " ******************************\n";
        std::cout << " Загрузка начальной конфигурации...\n";
        std::cout << " Общее количество частиц " << ntotal << "\n";
        std::cout << " ******************************\n";

        // Обеспечиваем размер контейнера под загружаемое число частиц
        // (контейнер создаётся пустым в main.cpp и наполняется здесь)
        if (sys.size() < ntotal)
        {
            sys.resize(ntotal);
        }
        // Чтение данных для каждой частицы
        for (size_t i = 0; i < ntotal; ++i)
        {
            std::string line;
            std::getline(file, line);
            if (line.empty()) continue;

            std::istringstream iss(line);
            int particle_id;
            iss >> particle_id;

            // Проверка идентификатора частицы
            if (particle_id != i + 1)
            {
                std::cerr << "Предупреждение: ожидалась частица " << i + 1
                          << ", найдена " << particle_id << "\n";
            }

            // Чтение координат
            for (size_t d = 0; d < dim; ++d)
            {
                iss >> sys.x[d][i];
            }

            // Чтение скоростей
            for (size_t d = 0; d < dim; ++d)
            {
                iss >> sys.vx[d][i];
            }

            // Чтение массы
            iss >> sys.mass[i];

            // Чтение плотности
            iss >> sys.rho[i];

            // Чтение давления
            iss >> sys.p[i];

            // Чтение внутренней энергии
            iss >> sys.u[i];

            // Чтение типа частицы
            iss >> sys.itype[i];

            // Чтение длины сглаживания
            iss >> sys.hsml[i];

            // body[] = itype[] (fiber/PLAN.md §3.1): формат particles.dat не
            // содержит отдельного номера тела.
            sys.body[i] = sys.itype[i];

            // Проверка на ошибки чтения
            if (iss.fail())
            {
                throw std::runtime_error("Ошибка чтения данных для частицы " + std::to_string(i + 1));
            }
        }

        file.close();
    }
    // Иначе генерируем данные
    else
    {
        if (shocktube)
        {
            shock_tube(cfg, sys, ntotal, dx);

        }
        else if (shearcavity)
        {
            shear_cavity(sys, ntotal);
        }
        std::wcout << " ******************************\n";
        std::wcout << " Начальная конфигурация сгенерирована\n";
        std::wcout << " Общее количество частиц " << ntotal << "\n";
        std::wcout << " ******************************\n";
    }
}

// Сводка импорта .k (readk/k-import-plan.md, раздел 4.4).
static void print_k_summary(const SimulationConfig& cfg, const k_input::KImportSummary& sm)
{
    std::cout << " ******************************\n";
    std::cout << " Импорт LS-DYNA .k: " << cfg.k_file << " (режим " << cfg.k_mode << ")\n";
    std::cout << " Узлов в колоде " << sm.nodes << "\n";
    std::cout << " Частиц: SPH " << sm.sph_particles << ", из солидов " << sm.solid_element_particles << "\n";
    std::cout << " Габарит SPH, м: x [" << sm.bbox_min[0] << ", " << sm.bbox_max[0] << "]"
              << ", y [" << sm.bbox_min[1] << ", " << sm.bbox_max[1] << "]"
              << ", z [" << sm.bbox_min[2] << ", " << sm.bbox_max[2] << "]\n";
    std::cout << " dx_median, м: " << sm.dx_median << "\n";
    if (sm.fibers > 0 || sm.beads > 0)
    {
        std::cout << " Нитей " << sm.fibers << ", бусин " << sm.beads << ", связей " << sm.bonds << "\n";
        std::cout << " Сегмент l0, м: [" << sm.seg_min << ", " << sm.seg_max << "], l/d: ["
                  << sm.ld_min << ", " << sm.ld_max << "]\n";
    }
    if (sm.clamped_particles > 0 || sm.clamped_beads > 0)
    {
        std::cout << " Зажато частиц " << sm.clamped_particles << ", бусин " << sm.clamped_beads << "\n";
    }
    for (const std::string& w : sm.warnings)
    {
        std::cout << " Предупреждение: " << w << "\n";
    }
    for (const auto& kv : sm.skipped)
    {
        std::cout << " Пропущено " << kv.first << ": " << kv.second << "\n";
    }
    std::cout << " ******************************\n";
}

// Основная функция ввода данных
void input(SimulationConfig& cfg, fiber::Model& model, size_t& ntotal, double& dx)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    parse_sph_config(cfg, "Initial_Conditions.txt");

    const bool solid_fiber = (cfg.k_mode == "solid_fiber");
    // Р10: два источника нитей одновременно не допускаются.
    if (solid_fiber && !cfg.fiber_file.empty())
    {
        throw std::runtime_error("k_mode=solid_fiber вместе с fiber_file: два источника нитей");
    }

    if (!cfg.k_file.empty())
    {
        if (config_input)
        {
            throw std::runtime_error("k_file несовместим с config_input (загрузка particles.dat)");
        }
        k_input::KImportOptions opts;
        opts.k_file = cfg.k_file;
        opts.mode = solid_fiber ? k_input::KMode::SolidFiber : k_input::KMode::Solid;
        opts.parts_file = cfg.k_parts_file;
        opts.length_scale = cfg.k_length_scale;
        opts.mass_scale = cfg.k_mass_scale;
        opts.monitor_nid = cfg.k_monitor_nid;
        opts.solid_elements = (cfg.k_solid_elements == "particles") ? k_input::KSolidElements::Particles
                                                                    : k_input::KSolidElements::Skip;
        const k_input::KImportSummary summary = k_input::load_k_geometry(opts, cfg.mats, model, cfg.moni_particle);
        ntotal = model.solid.size();
        // Для печати/совместимости: интегратор delta_x не использует.
        dx = summary.dx_median;
        print_k_summary(cfg, summary);
    }
    else
    {
        input_builtin(cfg, model.solid, ntotal, dx);
    }

    if (!cfg.fiber_file.empty())
    {
        fiber::load_fiber_file(cfg.fiber_file, model);
        std::size_t fibers = 0;
        for (std::size_t k = 0; k < model.fiber.size(); ++k)
        {
            fibers = std::max<std::size_t>(fibers, static_cast<std::size_t>(model.fiber.object_id[k]) + 1);
        }
        std::cout << " Нити из " << cfg.fiber_file << ": бусин " << model.fiber.size()
                  << ", связей " << model.topology.bonds.size() << ", нитей " << fibers << "\n";
    }
}
