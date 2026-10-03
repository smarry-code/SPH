// test_k_legacy_equiv.cpp
// Эквивалентность старого parse_k_file (копия из git, эталон) и нового
// k_input::load_k_geometry без таблицы деталей на корневом 1.k (этап K4,
// readk/k-import-plan.md, разд. 5): x, mass, itype, body - побитово.
// Из эталона вырезан хак мониторной частицы с записью в log.txt (к геометрии
// не относится).
#include "k_input/k_geometry.hpp"
#include "param.hpp"
#include "simulation_config.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef SPH_SOURCE_DIR
#error "SPH_SOURCE_DIR должен быть задан CMake"
#endif

namespace
{

// ----------------------------- эталон: старый parse_k_file -----------------
struct RawSph
{
    long long nid;           // ID узла/частицы
    int pid;                 // PID / signed material type
    double mass;             // MASS (>0 по допущению)
};
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

void legacy_parse_k_file(const std::string& filename, ParticleSystem& sys)
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
        ++idx;
    }
    in.close();
}

bool bit_equal(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

} // namespace

int main()
{
    const std::string path = std::string(SPH_SOURCE_DIR) + "/1.k";
    ParticleSystem legacy;
    legacy_parse_k_file(path, legacy);

    SimulationConfig cfg;
    cfg.apply_defaults();
    fiber::Model model;
    k_input::KImportOptions opts;
    opts.k_file = path;
    std::size_t moni = 0;
    const auto t0 = std::chrono::steady_clock::now();
    const k_input::KImportSummary sum = k_input::load_k_geometry(opts, cfg.mats, model, moni);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("load_k_geometry(1.k): %zu частиц, %.3f с\n", model.solid.size(), sec);
    for (const auto& w : sum.warnings) std::printf("  warning: %s\n", w.c_str());

    const ParticleSystem& now = model.solid;
    int failures = 0;
    auto check = [&](bool c, const char* what) {
        std::printf("  %-5s: %s\n", c ? "ok" : "FAIL", what);
        if (!c) ++failures;
    };
    check(legacy.size() > 0 && legacy.size() == now.size(), "одинаковое число частиц");
    if (legacy.size() == now.size())
    {
        bool bx = true, bm = true, bi = true, bb = true;
        for (std::size_t i = 0; i < legacy.size(); ++i)
        {
            for (std::size_t d = 0; d < dim; ++d)
                if (!bit_equal(legacy.x[d][i], now.x[d][i])) bx = false;
            if (!bit_equal(legacy.mass[i], now.mass[i])) bm = false;
            if (legacy.itype[i] != now.itype[i]) bi = false;
            if (legacy.body[i] != now.body[i]) bb = false;
        }
        check(bx, "x побитово");
        check(bm, "mass побитово");
        check(bi, "itype");
        check(bb, "body");
    }
    if (failures)
    {
        std::printf("\nFAILED: %d\n", failures);
        return 1;
    }
    std::printf("\nвсе проверки пройдены\n");
    return 0;
}
