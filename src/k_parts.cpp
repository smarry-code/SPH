// k_parts.cpp
// Чтение таблицы деталей для импорта LS-DYNA .k (этап K3,
// readk/k-import-plan.md, Р2). Формат - в k_input/k_parts.hpp.
#include "k_input/k_parts.hpp"

#include <cmath>
#include <fstream>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>

#include "fiber/fiber_build.hpp"

namespace k_input
{

namespace
{

[[noreturn]] void fail(const std::string& file, std::size_t line, const std::string& msg)
{
    throw std::runtime_error(file + ":" + std::to_string(line) + ": " + msg);
}

// Токен целиком как double (локаль classic), только конечные значения.
bool parse_double(const std::string& tok, double& out)
{
    std::istringstream s(tok);
    s.imbue(std::locale::classic());
    double v = 0.0;
    if (!(s >> v)) return false;
    if (s.peek() != std::char_traits<char>::eof()) return false;
    if (!std::isfinite(v)) return false;
    out = v;
    return true;
}

bool parse_ll(const std::string& tok, long long& out)
{
    std::istringstream s(tok);
    s.imbue(std::locale::classic());
    long long v = 0;
    if (!(s >> v)) return false;
    if (s.peek() != std::char_traits<char>::eof()) return false;
    out = v;
    return true;
}

} // namespace

const SolidPart* KPartsTable::find_solid(long long pid) const
{
    auto it = solid_index.find(pid);
    return it == solid_index.end() ? nullptr : &solid[it->second];
}

const FiberPart* KPartsTable::find_fiber(long long pid) const
{
    auto it = fiber_index.find(pid);
    return it == fiber_index.end() ? nullptr : &fiber[it->second];
}

KPartsTable read_k_parts_from_string(const std::string& text, const std::string& virtual_name)
{
    KPartsTable t;
    t.file = virtual_name;

    std::istringstream in(text);
    std::string raw;
    std::size_t line_no = 0;
    while (std::getline(in, raw))
    {
        ++line_no;
        const auto hash = raw.find('#');
        if (hash != std::string::npos) raw.erase(hash);
        std::istringstream iss(raw);
        std::string word;
        if (!(iss >> word)) continue;

        if (word == "SOLID_PART")
        {
            SolidPart p;
            p.line = line_no;
            std::string tp, ti, tr;
            if (!(iss >> tp >> ti >> tr))
                fail(t.file, line_no, "SOLID_PART: ожидается <pid> <itype> <rho>");
            long long itype = 0;
            if (!parse_ll(tp, p.pid)) fail(t.file, line_no, "SOLID_PART: PID не целое: '" + tp + "'");
            if (!parse_ll(ti, itype)) fail(t.file, line_no, "SOLID_PART: itype не целое: '" + ti + "'");
            if (!parse_double(tr, p.rho))
                fail(t.file, line_no, "SOLID_PART: rho не конечное число: '" + tr + "'");
            if (itype < 1 || itype > 11)
                fail(t.file, line_no, "SOLID_PART: itype=" + ti + " вне диапазона 1..11");
            p.itype = static_cast<int>(itype);
            if (!(p.rho > 0.0)) fail(t.file, line_no, "SOLID_PART: rho должно быть > 0");

            // Необязательные именованные параметры: hsml=<h>, v=<vx> [vy [vz]].
            std::string tok;
            bool got = static_cast<bool>(iss >> tok);
            while (got)
            {
                if (tok.rfind("hsml=", 0) == 0)
                {
                    if (p.has_hsml) fail(t.file, line_no, "SOLID_PART: hsml= повторяется");
                    const std::string val = tok.substr(5);
                    if (!parse_double(val, p.hsml))
                        fail(t.file, line_no, "SOLID_PART: hsml= не конечное число: '" + val + "'");
                    if (!(p.hsml > 0.0)) fail(t.file, line_no, "SOLID_PART: hsml должно быть > 0");
                    p.has_hsml = true;
                    got = static_cast<bool>(iss >> tok);
                }
                else if (tok.rfind("v=", 0) == 0)
                {
                    if (p.has_v) fail(t.file, line_no, "SOLID_PART: v= повторяется");
                    p.has_v = true;
                    std::string val = tok.substr(2);
                    if (val.empty())
                    {
                        if (!(iss >> val)) fail(t.file, line_no, "SOLID_PART: v= без значений");
                    }
                    // Первое число уже в val, затем до двух следующих токенов-чисел.
                    for (;;)
                    {
                        double c = 0.0;
                        if (!parse_double(val, c))
                            fail(t.file, line_no, "SOLID_PART: v= нечисловое значение: '" + val + "'");
                        p.v[static_cast<std::size_t>(p.nv)] = c;
                        ++p.nv;
                        got = static_cast<bool>(iss >> tok);
                        if (!got || p.nv == 3) break;
                        // Следующий токен - именованный параметр или лишний мусор.
                        if (tok.find('=') != std::string::npos) break;
                        val = tok;
                    }
                }
                else
                {
                    fail(t.file, line_no, "SOLID_PART: лишний или неизвестный токен '" + tok + "'");
                }
            }

            if (t.solid_index.count(p.pid) || t.fiber_index.count(p.pid))
                fail(t.file, line_no, "повтор PID " + std::to_string(p.pid));
            t.solid_index[p.pid] = t.solid.size();
            t.solid.push_back(p);
        }
        else if (word == "MATERIAL")
        {
            // Сообщение общего разбора MATERIAL дополняем префиксом "<файл>:<строка>:".
            try
            {
                fiber::parse_material_directive(iss, line_no, "k_parts", t.materials,
                                                t.material_defined);
            }
            catch (const std::runtime_error& e)
            {
                fail(t.file, line_no, e.what());
            }
        }
        else if (word == "FIBER_PART")
        {
            FiberPart p;
            p.line = line_no;
            std::string tp, tm, extra;
            if (!(iss >> tp >> tm)) fail(t.file, line_no, "FIBER_PART: ожидается <pid> <material_id>");
            long long mid = 0;
            if (!parse_ll(tp, p.pid)) fail(t.file, line_no, "FIBER_PART: PID не целое: '" + tp + "'");
            if (!parse_ll(tm, mid) || mid < 0 || mid > 1000000)
                fail(t.file, line_no, "FIBER_PART: неверный material_id: '" + tm + "'");
            if (iss >> extra) fail(t.file, line_no, "FIBER_PART: лишний токен '" + extra + "'");
            p.material_id = static_cast<int>(mid);
            if (t.solid_index.count(p.pid) || t.fiber_index.count(p.pid))
                fail(t.file, line_no, "повтор PID " + std::to_string(p.pid));
            t.fiber_index[p.pid] = t.fiber.size();
            t.fiber.push_back(p);
        }
        else
        {
            fail(t.file, line_no, "неизвестная директива '" + word + "'");
        }
    }

    // MATERIAL может идти после FIBER_PART - проверяем по концу файла.
    for (const auto& fp : t.fiber)
    {
        const auto id = static_cast<std::size_t>(fp.material_id);
        if (id >= t.material_defined.size() || !t.material_defined[id])
            fail(t.file, fp.line,
                 "FIBER_PART ссылается на необъявленный MATERIAL " + std::to_string(fp.material_id));
    }
    return t;
}

KPartsTable read_k_parts_file(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error(path + ":0: не удалось открыть файл таблицы деталей");
    std::ostringstream ss;
    ss << f.rdbuf();
    return read_k_parts_from_string(ss.str(), path);
}

void check_parts_against_deck(const KPartsTable& table, const KDeck& deck, bool fiber_mode,
                              std::vector<std::string>& warnings)
{
    for (const auto& e : deck.sph)
        if (!table.find_solid(e.pid))
            throw std::runtime_error(table.file + ": PID " + std::to_string(e.pid) + " из *ELEMENT_SPH ("
                                     + deck.where(e.src) + ") не описан в SOLID_PART");
    if (fiber_mode)
        for (const auto& b : deck.beams)
            if (!table.find_fiber(b.pid))
                throw std::runtime_error(table.file + ": PID " + std::to_string(b.pid)
                                         + " из *ELEMENT_BEAM (" + deck.where(b.src)
                                         + ") не описан в FIBER_PART");

    std::set<long long> known;
    for (const auto& e : deck.sph) known.insert(e.pid);
    for (const auto& e : deck.beams) known.insert(e.pid);
    for (const auto& e : deck.solids) known.insert(e.pid);
    for (const auto& e : deck.parts) known.insert(e.pid);

    auto warn = [&](long long pid, std::size_t line, const char* kind) {
        warnings.push_back(table.file + ":" + std::to_string(line) + ": " + kind + " PID "
                           + std::to_string(pid) + " отсутствует в .k");
    };
    for (const auto& p : table.solid)
        if (!known.count(p.pid)) warn(p.pid, p.line, "SOLID_PART");
    for (const auto& p : table.fiber)
        if (!known.count(p.pid)) warn(p.pid, p.line, "FIBER_PART");
}

} // namespace k_input
