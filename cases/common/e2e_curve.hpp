#pragma once
// e2e_curve.hpp
// Кривая монитора, её чтение, прореживание и обмен между процессами.
//
// Происхождение — tests/test_e2e_bc2.cpp (раздел 2, начало, и раздел 8).
// Сам BC2 на этот заголовок не переключён, см. пояснение в e2e_harness.hpp.
//
// Что параметризовано против оригинала: префикс имени файла выгрузки (был
// литерал "bc2_curve_") и название кейса в шапке выгрузки.

#include <array>
#include <algorithm>
#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace e2e
{

namespace fs = std::filesystem;

struct Curve
{
    std::vector<double> t;
    std::vector<double> v;
    std::vector<int>    itype;
    bool empty() const { return t.empty(); }
    size_t size() const { return t.size(); }
};

// Чтение файла формата log.txt / эталона: три поля «время v itype», строки
// с '#' и пустые — комментарии. Возможен UTF-8 BOM в начале (его пишет
// main.cpp в реальный log.txt).
inline Curve read_curve(const fs::path& p, std::string& err)
{
    Curve c;
    std::ifstream f(p);
    if (!f)
    {
        err = "не удалось открыть " + p.string();
        return c;
    }
    std::string line;
    while (std::getline(f, line))
    {
        size_t b = 0;
        if (line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF &&
            static_cast<unsigned char>(line[1]) == 0xBB && static_cast<unsigned char>(line[2]) == 0xBF)
        {
            b = 3;
        }
        while (b < line.size() && (line[b] == ' ' || line[b] == '\t' || line[b] == '\r')) ++b;
        if (b >= line.size() || line[b] == '#') continue;
        std::istringstream is(line.substr(b));
        double t = 0.0, v = 0.0;
        int it = 0;
        if (!(is >> t >> v >> it)) continue;
        c.t.push_back(t);
        c.v.push_back(v);
        c.itype.push_back(it);
    }
    return c;
}

// Прореживание: на N интервалов, в каждом сохраняются первая точка и точки
// минимума/максимума. Полка и пики выживают, файл не раздувается.
inline Curve thin_curve(const Curve& c, size_t buckets = 500)
{
    if (c.t.size() <= buckets) return c;
    Curve out;
    const double t0 = c.t.front(), t1 = c.t.back();
    const double step = (t1 - t0) / static_cast<double>(buckets);
    if (!(step > 0.0)) return c;

    size_t k = 0;
    double edge = t0 + step;
    while (k < c.t.size())
    {
        size_t lo = k, hi = k;
        const size_t first = k;
        while (k < c.t.size() && c.t[k] <= edge)
        {
            if (c.v[k] < c.v[lo]) lo = k;
            if (c.v[k] > c.v[hi]) hi = k;
            ++k;
        }
        if (k == first)   // ни одна точка не попала в интервал — двигаем границу
        {
            edge += step;
            continue;
        }
        std::array<size_t, 3> pick{first, lo, hi};
        std::sort(pick.begin(), pick.end());
        size_t prev = static_cast<size_t>(-1);
        for (size_t idx : pick)
        {
            if (idx == prev) continue;
            prev = idx;
            out.t.push_back(c.t[idx]);
            out.v.push_back(c.v[idx]);
            out.itype.push_back(c.itype[idx]);
        }
        edge += step;
    }
    // Последняя точка кривой должна попасть на график всегда.
    if (!out.t.empty() && out.t.back() != c.t.back())
    {
        out.t.push_back(c.t.back());
        out.v.push_back(c.v.back());
        out.itype.push_back(c.itype.back());
    }
    return out;
}

// ===========================================================================
// Обмен кривыми между процессами
// ===========================================================================
// Сценарии контакта разнесены по ОТДЕЛЬНЫМ прогонам теста (по одному на
// contact_force_type), поэтому в одном процессе все четыре кривые не живут.
// Каждый прогон выгружает свою кривую (--emit-curve), а рисование — отдельный
// дешёвый шаг без расчёта (--plot-from). Связь «график не расходится с
// метриками отчёта» при этом сохраняется: в файл пишется ровно тот объект
// Curve, который проверил компаратор, без пересчёта и без повторного чтения
// log.txt.
struct PlotSeries
{
    Curve  curve;
    int    type    = 0;      // contact_force_type
    bool   aborted = false;  // прогон не дошёл до t_total
    size_t steps   = 0;
};

// Имя файла выгрузки внутри каталога обмена. Фиксировано, потому что по нему
// --plot-from собирает серии, не получая списка файлов.
inline fs::path curve_dump_path(const fs::path& dir, const char* case_prefix, int type)
{
    return dir / (std::string(case_prefix) + "_curve_" + std::to_string(type) + ".txt");
}

// Формат — тот же, что у эталона и log.txt: <время> <скорость> <itype>.
// Признаки прогона идут заголовочными комментариями, которые read_curve
// пропускает, а load_series читает.
inline bool write_curve_dump(const PlotSeries& s, const fs::path& path,
                             const char* case_name, std::string& err)
{
    std::error_code ec;
    if (path.has_parent_path()) fs::create_directories(path.parent_path(), ec);
    std::ofstream f(path, std::ios::trunc);
    if (!f)
    {
        err = "не удалось открыть " + path.string();
        return false;
    }
    f << "# Кривая " << case_name << " для сводного графика (см. --plot-from)\n"
      << "# contact_force_type = " << s.type << "\n"
      << "# steps = " << s.steps << "\n"
      << "# aborted = " << (s.aborted ? 1 : 0) << "\n"
      << "# Формат: <время, с>  <v_z монитора, м/с>  <itype монитора>\n";
    f.precision(12);
    for (size_t i = 0; i < s.curve.t.size(); ++i)
    {
        f << s.curve.t[i] << ' ' << s.curve.v[i] << ' ' << s.curve.itype[i] << '\n';
    }
    f.close();
    if (!f)
    {
        err = "ошибка записи " + path.string();
        return false;
    }
    return true;
}

// Разбор одного числа из строки вида "# ключ = значение".
inline bool header_value(const std::string& line, const char* key, long long& out)
{
    const std::string pat = std::string("# ") + key + " = ";
    if (line.rfind(pat, 0) != 0) return false;
    try { out = std::stoll(line.substr(pat.size())); }
    catch (const std::exception&) { return false; }
    return true;
}

inline bool load_series(const fs::path& path, PlotSeries& out, std::string& err)
{
    out.curve = read_curve(path, err);
    if (!err.empty()) return false;
    if (out.curve.empty())
    {
        err = "в " + path.string() + " нет точек кривой";
        return false;
    }
    std::ifstream f(path);
    std::string line;
    long long v = 0;
    while (std::getline(f, line))
    {
        if (!line.empty() && line[0] != '#') break;
        if (header_value(line, "contact_force_type", v)) out.type = static_cast<int>(v);
        else if (header_value(line, "steps", v))         out.steps = static_cast<size_t>(v);
        else if (header_value(line, "aborted", v))       out.aborted = (v != 0);
    }
    return true;
}

// Названия типов для легенды. Короткие: легенда идёт в одну строку.
inline const char* contact_type_name(int type)
{
    switch (type)
    {
        case 0:  return "выключен";
        case 1:  return "Пинбол";
        case 2:  return "Свигл";
        case 3:  return "импульс";
        default: return "?";
    }
}

} // namespace e2e
