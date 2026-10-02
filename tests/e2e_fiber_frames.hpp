#pragma once
// e2e_fiber_frames.hpp
// Разбор кадров нитей fiber_output_<itimestep>.txt (src/fiber_output.cpp):
// строгое чтение одного кадра, перечисление кадров каталога по возрастанию
// номера шага, проверка живучести (NaN/Inf) и поиск бусины по id.
// Парный к tests/e2e_frames.hpp (тот читает sph_output_*.txt Solid-части).
//
// ФОРМАТ (проверено по src/fiber_output.cpp; значения печатаются потоком по
// умолчанию, то есть с 6 значащими цифрами — сравнивать с записанным надо с
// относительным допуском порядка 1e-5, а не побитово):
//   строка 1: "TIME <t> STEP <it> N_PARTICLES <n> N_BONDS <nb>"
//   "# id object_id material_id x y z vx vy vz"   (комментарий)
//   "PARTICLES"
//   n строк: "<id> <object_id> <material_id> <x> <y> <z> <vx> <vy> <vz>"
//   "# i j l0 l l_break elongation"               (комментарий)
//   "BONDS"
//   nb строк: "<i> <j> <l0> <l> <l_break> <elongation>"
//
// id — индекс бусины в FiberState (0-based) на момент записи; он устойчив
// после разрыва связи (меняется object_id, а не id), поэтому бусину между
// кадрами надо искать именно по id. Читатель не требует, чтобы id шли подряд:
// поиск fiber_frame_bead() линейный по id, а при плотной нумерации — прямой.
//
// Читатель ЧИСЛА разбирает через strtod, поэтому "nan"/"inf" в данных не
// ломают разбор формата: они читаются как есть, а их учёт — дело
// scan_fiber_frames(). Нарушение структуры (счётчики, число полей, лишний
// текст) — ошибка чтения.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace e2e
{

namespace fs = std::filesystem;

struct FiberFrame
{
    struct Bond
    {
        std::size_t i = 0, j = 0;
        double l0 = 0.0, l = 0.0, l_break = 0.0, elong = 0.0;
    };

    double time = 0.0;
    std::size_t step = 0;
    std::vector<std::size_t> id;
    std::vector<int> object_id, material_id;
    std::vector<std::array<double, 3>> x, v;
    std::vector<Bond> bonds;
    // Счётчики из заголовка (проверяются на равенство фактическому числу строк).
    std::size_t n_particles_header = 0, n_bonds_header = 0;
};

namespace fiber_frames_detail
{

inline bool fail(std::string* err, const std::string& msg)
{
    if (err) *err = msg;
    return false;
}

// Токен целиком -> double (strtod: принимает nan/inf, в отличие от operator>>).
inline bool parse_double(const std::string& tok, double& out)
{
    if (tok.empty()) return false;
    char* end = nullptr;
    out = std::strtod(tok.c_str(), &end);
    return end == tok.c_str() + tok.size();
}

inline bool parse_size(const std::string& tok, std::size_t& out)
{
    if (tok.empty()) return false;
    for (char ch : tok)
        if (ch < '0' || ch > '9') return false;
    char* end = nullptr;
    out = static_cast<std::size_t>(std::strtoull(tok.c_str(), &end, 10));
    return end == tok.c_str() + tok.size();
}

inline bool parse_int(const std::string& tok, int& out)
{
    if (tok.empty()) return false;
    char* end = nullptr;
    const long v = std::strtol(tok.c_str(), &end, 10);
    if (end != tok.c_str() + tok.size()) return false;
    out = static_cast<int>(v);
    return true;
}

inline std::vector<std::string> split(const std::string& line)
{
    std::vector<std::string> t;
    std::istringstream is(line);
    std::string s;
    while (is >> s) t.push_back(s);
    return t;
}

// Следующая значащая строка (не пустая, не комментарий '#'); false — конец файла.
inline bool next_line(std::ifstream& in, std::string& line)
{
    while (std::getline(in, line))
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::vector<std::string> t = split(line);
        if (t.empty() || t[0][0] == '#') continue;
        return true;
    }
    return false;
}

} // namespace fiber_frames_detail

// Строгий разбор кадра. false + *err при любом отклонении от формата.
inline bool read_fiber_frame(const fs::path& path, FiberFrame& fr, std::string* err = nullptr)
{
    using namespace fiber_frames_detail;
    fr = FiberFrame{};
    const std::string name = path.filename().string();
    std::ifstream in(path);
    if (!in) return fail(err, "не удалось открыть " + name);

    std::string line;
    if (!next_line(in, line)) return fail(err, name + ": пустой файл");
    {
        const auto t = split(line);
        if (t.size() != 8 || t[0] != "TIME" || t[2] != "STEP" || t[4] != "N_PARTICLES" || t[6] != "N_BONDS")
            return fail(err, name + ": неверная строка заголовка");
        if (!parse_double(t[1], fr.time) || !parse_size(t[3], fr.step) ||
            !parse_size(t[5], fr.n_particles_header) || !parse_size(t[7], fr.n_bonds_header))
            return fail(err, name + ": не разобраны числа в заголовке");
    }

    if (!next_line(in, line) || split(line) != std::vector<std::string>{"PARTICLES"})
        return fail(err, name + ": нет слова PARTICLES после заголовка");
    for (std::size_t k = 0; k < fr.n_particles_header; ++k)
    {
        if (!next_line(in, line))
            return fail(err, name + ": строк PARTICLES меньше, чем N_PARTICLES (" + std::to_string(k) + " из " +
                                 std::to_string(fr.n_particles_header) + ")");
        const auto t = split(line);
        if (t.size() != 9)
            return fail(err, name + ": строка частицы " + std::to_string(k) + ": полей " + std::to_string(t.size()) +
                                 " вместо 9 (либо BONDS раньше времени)");
        std::size_t id = 0;
        int oid = 0, mid = 0;
        std::array<double, 3> x{}, v{};
        bool ok = parse_size(t[0], id) && parse_int(t[1], oid) && parse_int(t[2], mid);
        for (std::size_t d = 0; ok && d < 3; ++d) ok = parse_double(t[3 + d], x[d]) && parse_double(t[6 + d], v[d]);
        if (!ok) return fail(err, name + ": не разобрана строка частицы " + std::to_string(k));
        fr.id.push_back(id);
        fr.object_id.push_back(oid);
        fr.material_id.push_back(mid);
        fr.x.push_back(x);
        fr.v.push_back(v);
    }

    if (!next_line(in, line) || split(line) != std::vector<std::string>{"BONDS"})
        return fail(err, name + ": нет слова BONDS (строк PARTICLES больше N_PARTICLES или блок оборван)");
    for (std::size_t k = 0; k < fr.n_bonds_header; ++k)
    {
        if (!next_line(in, line))
            return fail(err, name + ": строк BONDS меньше, чем N_BONDS (" + std::to_string(k) + " из " +
                                 std::to_string(fr.n_bonds_header) + ")");
        const auto t = split(line);
        if (t.size() != 6)
            return fail(err, name + ": строка связи " + std::to_string(k) + ": полей " + std::to_string(t.size()) +
                                 " вместо 6");
        FiberFrame::Bond b;
        if (!(parse_size(t[0], b.i) && parse_size(t[1], b.j) && parse_double(t[2], b.l0) && parse_double(t[3], b.l) &&
              parse_double(t[4], b.l_break) && parse_double(t[5], b.elong)))
            return fail(err, name + ": не разобрана строка связи " + std::to_string(k));
        fr.bonds.push_back(b);
    }
    if (next_line(in, line))
        return fail(err, name + ": лишние строки после N_BONDS связей (BONDS больше N_BONDS)");
    return true;
}

// Номер шага из имени fiber_output_<ddddd>.txt; false, если имя не подходит.
inline bool fiber_frame_timestep(const std::string& name, std::size_t& out)
{
    static const std::string prefix = "fiber_output_";
    if (name.rfind(prefix, 0) != 0) return false;
    const std::size_t dot = name.rfind('.');
    if (dot == std::string::npos || name.substr(dot) != ".txt" || dot <= prefix.size()) return false;
    const std::string digits = name.substr(prefix.size(), dot - prefix.size());
    for (char ch : digits)
        if (ch < '0' || ch > '9') return false;
    try { out = static_cast<std::size_t>(std::stoull(digits)); }
    catch (const std::exception&) { return false; }
    return true;
}

// Кадры каталога, по возрастанию номера шага (числовому, не лексикографическому).
inline std::vector<fs::path> list_fiber_frames(const fs::path& dir)
{
    std::vector<std::pair<std::size_t, fs::path>> files;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(dir, ec))
    {
        std::size_t ts = 0;
        if (fiber_frame_timestep(e.path().filename().string(), ts)) files.emplace_back(ts, e.path());
    }
    std::sort(files.begin(), files.end());
    std::vector<fs::path> out;
    out.reserve(files.size());
    for (auto& f : files) out.push_back(std::move(f.second));
    return out;
}

struct FiberFrameScan
{
    std::size_t frames = 0;
    bool health = true;                 // все кадры разобраны и без NaN/Inf
    std::string first_problem;
    std::size_t nonfinite = 0;          // число не-конечных значений в x, v, l0, l, l_break, elong
    std::vector<std::size_t> n_particles, n_bonds;   // по кадрам (из разобранных строк)
    std::vector<double> time;                        // по кадрам (метка TIME)
};

inline FiberFrameScan scan_fiber_frames(const fs::path& dir)
{
    FiberFrameScan scan;
    const std::vector<fs::path> files = list_fiber_frames(dir);
    scan.frames = files.size();
    for (const fs::path& p : files)
    {
        FiberFrame fr;
        std::string err;
        if (!read_fiber_frame(p, fr, &err))
        {
            scan.health = false;
            if (scan.first_problem.empty()) scan.first_problem = err;
            scan.n_particles.push_back(0);
            scan.n_bonds.push_back(0);
            scan.time.push_back(std::nan(""));
            continue;
        }
        scan.n_particles.push_back(fr.id.size());
        scan.n_bonds.push_back(fr.bonds.size());
        scan.time.push_back(fr.time);

        std::size_t bad = 0;
        if (!std::isfinite(fr.time)) ++bad;
        for (std::size_t k = 0; k < fr.id.size(); ++k)
            for (std::size_t d = 0; d < 3; ++d)
            {
                if (!std::isfinite(fr.x[k][d])) ++bad;
                if (!std::isfinite(fr.v[k][d])) ++bad;
            }
        for (const auto& b : fr.bonds)
        {
            if (!std::isfinite(b.l0)) ++bad;
            if (!std::isfinite(b.l)) ++bad;
            if (!std::isfinite(b.l_break)) ++bad;
            if (!std::isfinite(b.elong)) ++bad;
        }
        if (bad > 0)
        {
            scan.nonfinite += bad;
            scan.health = false;
            if (scan.first_problem.empty())
                scan.first_problem = "не-конечное значение в " + p.filename().string();
        }
    }
    return scan;
}

// Бусина по id (id устойчив после разрыва связи). false, если такой нет.
inline bool fiber_frame_bead(const FiberFrame& fr, std::size_t id, std::array<double, 3>& x, std::array<double, 3>& v)
{
    std::size_t k = id;
    if (k >= fr.id.size() || fr.id[k] != id)   // быстрый путь: плотная нумерация
    {
        const auto it = std::find(fr.id.begin(), fr.id.end(), id);
        if (it == fr.id.end()) return false;
        k = static_cast<std::size_t>(it - fr.id.begin());
    }
    x = fr.x[k];
    v = fr.v[k];
    return true;
}

} // namespace e2e
