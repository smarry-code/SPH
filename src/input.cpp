// input.cpp
// Определения функций ввода/инициализации данных.
// Объявления — в input.hpp. Вспомогательные функции разбора (parse_k_file,
// shock_tube, shear_cavity) — внутренние (internal linkage / static) этого TU.
#include "input.hpp"
#include "fiber/fiber_system.hpp"
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

struct RawSph
{
    long long nid;           // ID узла/частицы
    int pid;                 // PID / signed material type
    double mass;             // MASS (>0 по допущению)
};
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
enum class BlockState
{
    None,
    Node,
    ElementSPH
};
static void trim_comment(std::string &line)
{
    std::size_t pos = line.find('$');
    if (pos != std::string::npos)
    {
        line.erase(pos);
    }        
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
static std::vector<std::string> split_by_commas_preserve_empty(const std::string& line)
{
    std::vector<std::string> fields;
    std::string current;

    for (char c : line)
    {
        if (c == ',')
        {
            // конец поля
            trim_spaces(current);
            fields.push_back(current);
            current.clear();
        }
        else
        {
            current.push_back(c);
        }
    }
    // последнее поле
    trim_spaces(current);
    fields.push_back(current);

    return fields;
}
static long long to_int_or_zero(const std::string& s)
{
    if (s.empty())
    {
        return 0;
    }      
    return std::stoll(s);
}

// безопасное чтение double: пустая строка => 0.0
static double to_double_or_zero(const std::string& s)
{
    if (s.empty())
    {
        return 0.0;
    }      
    return std::stod(s);
}

void parse_k_file(SimulationConfig& cfg, const std::string& filename, ParticleSystem& sys)
{
    std::ifstream in(filename);
    // узлы: NID -> (x,y,z)
    std::unordered_map<long long, std::array<double,3>> nodes;
    nodes.reserve(1024);

    // сырые SPH-элементы
    std::vector<RawSph> sph_elements;
    sph_elements.reserve(1024);

    BlockState state = BlockState::None;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }
        // заголовки ключевых слов
        if (!line.empty() && line[0] == '*')
        {
            std::string header = line;
            trim_comment(header);
            for (char& c : header)
            {
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
            if (header.find("*NODE") == 0)
            {
                state = BlockState::Node;
            }
            else if (header.find("*ELEMENT_SPH") == 0)
            {
                state = BlockState::ElementSPH;
            }
            else
            {
                state = BlockState::None;
            }
            continue;
        }
        // обрезаем комментарий и проверяем на пустоту
        trim_comment(line);
        bool only_space = true;
        for (char c : line)
        {
            if (!std::isspace(static_cast<unsigned char>(c)))
            {
                only_space = false;
                break;
            }
        }
        if (only_space)
        {
            continue;
        }
        // free format: делим по запятым
        auto fields = split_by_commas_preserve_empty(line);
        if (fields.empty())
        {
            continue;
        }  
        if (state == BlockState::Node)
        {
            // NID, x, y, z (недостающие поля -> 0)
            if (fields.size() < 4)
                fields.resize(4);

            long long nid = to_int_or_zero(fields[0]);
            double x = to_double_or_zero(fields[1]);
            double y = to_double_or_zero(fields[2]);
            double z = to_double_or_zero(fields[3]);
            if (nid < 0)
            {
                continue;
            }            
            nodes[nid] = {x, y, z};
        }
        else if (state == BlockState::ElementSPH)
        {
            // NID, PID, MASS (недостающие поля -> 0)
            if (fields.size() < 3)
            {
                fields.resize(3);
            }                
            long long nid = to_int_or_zero(fields[0]);
            int pid = fields[1].empty() ? 0 : std::stoi(fields[1]);
            double mass = to_double_or_zero(fields[2]);
            if (nid < 0)
            {
                continue;
            }                
            RawSph e{nid, pid, mass};
            sph_elements.push_back(e);
        }                  
    }
    // объединение: по каждому RawSph ищем координаты узла.
    // Сначала подсчитываем число валидных элементов (с найденным узлом),
    // затем resize sys и заполняем плоские массивы.
    size_t valid = 0;
    for (const auto& e : sph_elements)
    {
        if (nodes.find(e.nid) != nodes.end()) ++valid;
    }
    sys.resize(valid);
    size_t idx = 0;
    for (const auto& e : sph_elements)
    {
        auto it = nodes.find(e.nid);
        if (it == nodes.end())
        {
            // узел не найден — пропускаем
            continue;
        }
        sys.inzializ_index(idx);  // обнулить все поля (эквивалент default Particle)
        sys.x[0][idx] = it->second[0];
        sys.x[1][idx] = it->second[1];
        sys.x[2][idx] = it->second[2];
        sys.mass[idx]  = e.mass;
        sys.itype[idx] = e.pid;
        // body[] = PID (fiber/PLAN.md §3.1): "для k-файла это PID". Совпадает
        // с itype здесь только потому, что itype тоже читается из PID выше —
        // но body[] хранит именно номер тела (object_id), а не материал, и
        // остаётся отдельным полем для будущего Solid-Solid-Contact (Ф4).
        sys.body[idx]  = e.pid;
        if (sys.x[0][idx] > 0.0131 && sys.itype[idx] == 3 &&
            sys.x[1][idx] > -1.9e-4 && sys.x[1][idx] < 1.9e-4 &&
            sys.x[2][idx] > -1.9e-4 && sys.x[2][idx] < 1.9e-4)
        {
            cfg.moni_particle = idx;
            std::ofstream file("log.txt", std::ios::app);
            file << "Частица номер " << cfg.moni_particle << " с координатами " << sys.x[0][idx] << ", " << sys.x[1][idx]
                 << ", " << sys.x[2][idx] << std::endl;
            file.close();
        }
        ++idx;
    }
    in.close();
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
        else if (key == "fiber_file")
        {
            cfg.fiber_file = value;
            // fiber/PLAN.md, этап Ф2, пункт 5: если нити реально заданы, а
            // dim != 3, это ошибка, и она должна выявляться здесь, на этапе
            // разбора конфигурации, а не позже в середине расчёта.
            fiber::check_fiber_requires_dim3();
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

// Основная функция ввода данных
void input(SimulationConfig& cfg, ParticleSystem& sys, size_t& ntotal, double& dx)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    parse_sph_config(cfg, "Initial_Conditions.txt");
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
        /*parse_k_file("1.k", Elem);
        double delta_x = 0.0003;
        ntotal = sys.size();
        std::ofstream file("test.txt");
        file << "N" << " " << "Nmat" << " " << "M" << " ";
        if (dim == 1)
        {
            file << "x" << " ";
        }
        else if (dim == 2)
        {
            file << "x" << " " << "y" << " ";
        }
        else if (dim == 3)
        {
            file << "x" << " " << "y" << " " << "z" << " ";
        }
        if (dim == 1)
        {
            file << "vx" << " ";
        }
        else if (dim == 2)
        {
            file << "vx" << " " << "vy" << " ";
        }
        else if (dim == 3)
        {
            file << "vx" << " " << "vy" << " " << "vz" << " ";
        }
        file << "rho" << " " << "p" << " " << "u" << " " << "c" << " " << "sigma_i" << " " << "e_pl" << " " << "D" << " " << "hsml" << " "
             << "sxx" << " " << "syy" << " " << "szz" << " " << "sxy" << " " << "syz" << " " << "sxz" << "\n";
        for (size_t i = 0; i < ntotal; ++i)
        {
            sys.hsml[i] = 2 * delta_x;

            if (sys.x[0][i] < 0.0045)
            {
                sys.vx[0][i] = 1546;
            }
            if (sys.itype[i] == 1)
            {
                sys.itype[i] = 7;
                sys.rho[i] = 320;
                sys.rho_0[i] = sys.rho[i];
                //sys.vx[0][i] = 1546;
            }
            else if (sys.itype[i] == 2)
            {
                sys.itype[i] = 10;
                sys.rho[i] = 2510;
                sys.rho_0[i] = sys.rho[i];
            }
            else if (sys.itype[i] == 3)
            {
                sys.itype[i] = 10;
                sys.rho[i] = 2510;
                sys.rho_0[i] = sys.rho[i];
            }
            else if (sys.itype[i] == 4)
            {
                sys.itype[i] = 5;
                sys.rho[i] = 2638;
                sys.rho_0[i] = sys.rho[i];
            }
            // Запись основных данных
            file << i + 1; // ID частицы (в Fortran индексы начинаются с 1)
            sys.vuvod_index(i, file);
        }
        file.close();*/
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
