// input.cpp
#include <cmath>
#include <stdexcept>
#include <limits>
#include <unordered_map>
#include <array>
#include <algorithm>
struct RawSph
{
    long long nid;           // ID узла/частицы
    std::size_t pid;         // PID
    double mass;             // MASS (>0 по допущению)
};
// Вспомогательная функция для форматированного вывода
void format_output(std::ofstream& file, size_t particle_id, const std::vector<double>& values)
{
    file << std::setw(5) << particle_id;
    for (double val : values)
    {
        file << std::setw(8) << std::scientific << val;
    }
    file << "\n";
}

// Вспомогательная функция для форматированного вывода целых чисел
void format_output_other(std::ofstream& file, size_t particle_id, int itype_val, double hsml_val)
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
    auto not_space = [](int ch){ return !std::isspace(ch); };
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

void parse_k_file(const std::string& filename, std::vector<Particle>& result)
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
            std::size_t pid = static_cast<std::size_t>(to_int_or_zero(fields[1]));
            double mass = to_double_or_zero(fields[2]);
            if (nid < 0)
            {
                continue;
            }                
            RawSph e{nid, pid, mass};
            sph_elements.push_back(e);
        }                  
    }
    // объединение: по каждому RawSph ищем координаты узла
    result.clear();
    //result.reserve(sph_elements.size());
    for (const auto& e : sph_elements)
    {
        auto it = nodes.find(e.nid);
        if (it == nodes.end())
        {
            // узел не найден — можно пропустить или кинуть исключение
            // здесь просто пропускаем
            continue;
        }
        Particle p;
        p.x = {it->second[0], it->second[1], it->second[2]};
        p.mass = e.mass;
        p.itype = e.pid;
        result.push_back(std::move(p));
        if (result[result.size() - 1].x[0] > 0.0131 && result[result.size() - 1].itype == 3 && result[result.size() - 1].x[1] > -1.9e-4 && result[result.size() - 1].x[1] < 1.9e-4 && 
            result[result.size() - 1].x[2] > -1.9e-4 && result[result.size() - 1].x[2] < 1.9e-4)
        {
            moni_particle = result.size() - 1;
            std::ofstream file("log.txt", std::ios::app);
            file << "Частица номер " << moni_particle << " с координатами " << result[moni_particle].x[0] << ", " << result[moni_particle].x[1] 
                 << ", " << result[moni_particle].x[2] << std::endl;
            file.close();
        }
        //std::wcout << result.back().x[0] << " " << result.back().x[1] << " " << result.back().x[2] <<"\n";
    }
    in.close();
}
void parse_sph_config(const std::string& filename)
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
            pa_sph = std::stoi(value);
        }
        else if (key == "nnps")
        {
            nnps = std::stoi(value);
        }
        else if (key == "sle")
        {
            sle = std::stoi(value);
        }
        else if (key == "skf")
        {
            skf = std::stoi(value);
        }
        else if (key == "CFL")
        {
            CFL = std::stod(value);
        }
        else if (key == "epsilon")
        {
            epsilon = std::stod(value);
        }
        else if (key == "alpha")
        {
            alpha = std::stod(value);
        }
        else if (key == "beta")
        {
            beta = std::stod(value);
        }
        else if (key == "etq")
        {
            etq = std::stod(value);
        }
        else if (key == "t_total")
        {
            t_total = std::stod(value);
        }
        else if (key == "save_t")
        {
            save_t = std::stod(value);
        }
        else if (key == "print_step")
        {
            print_step = std::stoi(value);
        }
    }
    in.close();    
}
// Функция инициализации частиц для теста "ударная трубка"
void shock_tube(std::vector<Particle>& Elem, size_t& ntotal, double& delta_x)
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
    ntotal = 1000;
    //double R_ud = 0.05;
    double Z1 = 0.005;
    double Z2 = Z1 + 0.005;
    double Z3 = Z2 + 0.005;
    double Z4 = Z3 + 0.005;
    double L = Z4;
    double vx0 = 1000;
    delta_x = L / double (ntotal);
    bool flag_moni_particle = false;
    for (size_t i = 0; i < ntotal; ++i)
    {
        Elem[i].hsml = 2 * delta_x;
        Elem[i].u = 0.0;
        Elem[i].p = 0.0;
        Elem[i].x[0] = delta_x * (double(i) + 1);
        if (Elem[i].x[0] <= Z1)
        {
            Elem[i].itype = 4;
            Elem[i].rho = 2750;
            Elem[i].rho_0 = Elem[i].rho;
            Elem[i].mass = Elem[i].rho * delta_x;
            Elem[i].vx[0] = vx0;

        }
        else if (Elem[i].x[0] <= Z2)
        {
            Elem[i].itype = 4;
            Elem[i].rho = 2750;
            Elem[i].rho_0 = Elem[i].rho;
            Elem[i].mass = Elem[i].rho * delta_x;
            Elem[i].vx[0] = vx0;
        }
        else if (Elem[i].x[0] <= Z3)
        {
            Elem[i].itype = 4;
            Elem[i].rho = 2750;
            Elem[i].rho_0 = Elem[i].rho;
            Elem[i].mass = Elem[i].rho * delta_x;
            Elem[i].vx[0] = 0.0;
            Elem[i].x[0] += delta_x; //=======
        }
        else if (Elem[i].x[0] <= Z4)
        {
            if (flag_moni_particle == false)
            {
                flag_moni_particle = true;
                moni_particle = i - 1;
            }
            Elem[i].itype = 4;
            Elem[i].rho = 2750;
            Elem[i].rho_0 = Elem[i].rho;
            Elem[i].mass = Elem[i].rho * delta_x;
            Elem[i].vx[0] = 0.0;
            Elem[i].x[0] += delta_x; //=======
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
                Elem[ntotal].axis = true;
            }
            Elem[ntotal].hsml = 2.0 * delta_x;
            Elem[ntotal].u = 0.0;
            Elem[ntotal].p = 0.0;
            Elem[ntotal].x[0] = r;
            Elem[ntotal].x[1] = z;
            Elem[ntotal].itype = 8;
            Elem[ntotal].rho = 626;
            Elem[ntotal].rho_0 = Elem[ntotal].rho;
            Elem[ntotal].vx[0] = 0.0;
            Elem[ntotal].vx[1] = 2210;
            Elem[ntotal].mass = Elem[ntotal].rho * delta_x * delta_x;
            // Elem[ntotal].mass = Elem[ntotal].rho * 2 * pi * pi * 0.25 * delta_x * delta_x * Elem[ntotal].x[0];
            // Elem[ntotal].mass = 2 * pi * Elem[ntotal].x[0] * delta_x * delta_x * Elem[ntotal].rho;
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
                Elem[ntotal].axis = true;
            }
            Elem[ntotal].hsml = 2.0 * delta_x;
            Elem[ntotal].u = 0.0;
            Elem[ntotal].p = 0.0;
            Elem[ntotal].x[0] = r;
            Elem[ntotal].x[1] = z;
            Elem[ntotal].itype = 10;
            Elem[ntotal].rho = 2510;
            Elem[ntotal].rho_0 = Elem[ntotal].rho;
            Elem[ntotal].vx[0] = 0.0;
            Elem[ntotal].vx[1] = 2210;
            Elem[ntotal].mass = Elem[ntotal].rho * delta_x * delta_x;
            // Elem[ntotal].mass = Elem[ntotal].rho * 2 * pi * pi * 0.25 * delta_x * delta_x * Elem[ntotal].x[0];
            // Elem[ntotal].mass = 2 * pi * Elem[ntotal].x[0] * delta_x * delta_x * Elem[ntotal].rho;
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
                Elem[ntotal].axis = true;
            }
            Elem[ntotal].hsml = 2.0 * delta_x;
            Elem[ntotal].u = 0.0;
            Elem[ntotal].p = 0.0;
            Elem[ntotal].x[0] = r;
            Elem[ntotal].x[1] = z;
            Elem[ntotal].itype = 10;
            Elem[ntotal].rho = 2510;
            Elem[ntotal].rho_0 = Elem[ntotal].rho;
            Elem[ntotal].vx[0] = 0.0;
            Elem[ntotal].vx[1] = 0.0;
            Elem[ntotal].mass = Elem[ntotal].rho * delta_x * delta_x;
            // Elem[ntotal].mass = Elem[ntotal].rho * 2 * pi * pi * 0.25 * delta_x * delta_x * Elem[ntotal].x[0];
            // Elem[ntotal].mass = 2 * pi * Elem[ntotal].x[0] * delta_x * delta_x * Elem[ntotal].rho;
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
                Elem[ntotal].axis = true;
            }
            Elem[ntotal].hsml = 2.0 * delta_x;
            Elem[ntotal].u = 0.0;
            Elem[ntotal].p = 0.0;
            Elem[ntotal].x[0] = r;
            Elem[ntotal].x[1] = z;
            Elem[ntotal].itype = 5;
            Elem[ntotal].rho = 2638;
            Elem[ntotal].rho_0 = Elem[ntotal].rho;
            Elem[ntotal].vx[0] = 0.0;
            Elem[ntotal].vx[1] = 0.0;
            Elem[ntotal].mass = Elem[ntotal].rho * delta_x * delta_x;
            // Elem[ntotal].mass = Elem[ntotal].rho * 2 * pi * pi * 0.25 * delta_x * delta_x * Elem[ntotal].x[0];
            // Elem[ntotal].mass = 2 * pi * Elem[ntotal].x[0] * delta_x * delta_x * Elem[ntotal].rho;
            ntotal++;
        }
    }*/
}

// Функция инициализации частиц для теста "сдвиговая каверна"
void shear_cavity(std::vector<Particle>& Elem, size_t& ntotal)
{
    const size_t m = 41;
    const size_t n = 41;
    const size_t mp = m - 1;
    const size_t np = n - 1;
    ntotal = mp * np;

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
            Elem[k].x[0] = (double(i) * dx) + dx / 2.0;
            Elem[k].x[1] = (double(j) * dy) + dy / 2.0;
        }
    }

    // Инициализация остальных параметров
    for (size_t i = 0; i < ntotal; ++i)
    {
        Elem[i].vx[0] = 0.0;
        Elem[i].vx[1] = 0.0;
        Elem[i].rho = 1000.0;  // Плотность воды
        Elem[i].mass = dx * dy * Elem[i].rho;
        Elem[i].p = 0.0;
        Elem[i].u = 357.1;     // Внутренняя энергия
        Elem[i].itype = 2;     // Водяная частица
        Elem[i].hsml = dx;     // Длина сглаживания = начальный шаг
    }
}

// Основная функция ввода данных
void input(std::vector<Particle>& Elem, size_t& ntotal, double& dx)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    parse_sph_config("Initial_Conditions.txt");
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

        // Проверка размеров векторов
        for (size_t d = 0; d < dim; ++d)
        {
            if (Elem.size() < ntotal)
            {
                throw std::runtime_error("Недостаточный размер массивов x и vx");
            }
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
                iss >> Elem[i].x[d];
            }

            // Чтение скоростей
            for (size_t d = 0; d < dim; ++d)
            {
                iss >> Elem[i].vx[d];
            }

            // Чтение массы
            iss >> Elem[i].mass;

            // Чтение плотности
            iss >> Elem[i].rho;

            // Чтение давления
            iss >> Elem[i].p;

            // Чтение внутренней энергии
            iss >> Elem[i].u;

            // Чтение типа частицы
            iss >> Elem[i].itype;

            // Чтение длины сглаживания
            iss >> Elem[i].hsml;

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
        ntotal = Elem.size();
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
            Elem[i].hsml = 2 * delta_x;

            if (Elem[i].x[0] < 0.0045)
            {
                Elem[i].vx[0] = 1546;
            }
            if (Elem[i].itype == 1)
            {
                Elem[i].itype = 7;
                Elem[i].rho = 320;
                Elem[i].rho_0 = Elem[i].rho;
                //Elem[i].vx[0] = 1546;
            }
            else if (Elem[i].itype == 2)
            {
                Elem[i].itype = 10;
                Elem[i].rho = 2510;
                Elem[i].rho_0 = Elem[i].rho;
            }
            else if (Elem[i].itype == 3)
            {
                Elem[i].itype = 10;
                Elem[i].rho = 2510;
                Elem[i].rho_0 = Elem[i].rho;
            }
            else if (Elem[i].itype == 4)
            {
                Elem[i].itype = 5;
                Elem[i].rho = 2638;
                Elem[i].rho_0 = Elem[i].rho;
            }
            // Запись основных данных
            file << i + 1; // ID частицы (в Fortran индексы начинаются с 1)
            Elem[i].vuvod(file);
        }
        file.close();*/
        if (shocktube)
        {
            shock_tube(Elem, ntotal, dx);

        }
        else if (shearcavity)
        {
            shear_cavity(Elem, ntotal);
        }
        std::wcout << " ******************************\n";
        std::wcout << " Начальная конфигурация сгенерирована\n";
        std::wcout << " Общее количество частиц " << ntotal << "\n";
        std::wcout << " ******************************\n";
    }
}
