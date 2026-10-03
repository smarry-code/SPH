// k_reader.cpp
// Чтение LS-DYNA keyword (.k): ключевые слова, комментарии, free/fixed/long,
// *INCLUDE, сборка KDeck и проверка ссылок (этап K2, readk/k-import-plan.md).
// Объявления - в include/k_input/k_reader.hpp, ширины карточек - в k_cards.hpp.
#include "k_input/k_reader.hpp"

#include "k_input/k_cards.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace fs = std::filesystem;

namespace k_input
{

std::string KDeck::where(const SrcLoc& s) const
{
    const std::string& f = s.file < files.size() ? files[s.file] : std::string("?");
    return f + ":" + std::to_string(s.line);
}

const KNode* KDeck::find_node(long long nid) const
{
    auto it = node_index.find(nid);
    return it == node_index.end() ? nullptr : &nodes[it->second];
}

namespace
{

// ---- Мелкие строковые утилиты (без аллокаций) -------------------------------

std::string_view trim(std::string_view s)
{
    std::size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::string upper(std::string_view s)
{
    std::string r(s);
    for (char& c : r) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return r;
}

bool starts_with(const std::string& s, const char* p)
{
    return s.rfind(p, 0) == 0;
}

// Значение одного поля карточки.
struct Value
{
    long long        i = 0;
    double           r = 0.0;
    std::string_view s;
    bool             set = false;   // поле было непустым
};

struct Values
{
    std::array<Value, kMaxCardFields> v;
};

// ---- Состояние блока ---------------------------------------------------------

enum class Kind
{
    None, Skip, Keyword, Title, Include, IncludePath, IncludePathRel,
    Node, Sph, Beam, Solid, Solid10, Part, SetList, SpcNode, SpcSet
};

struct Block
{
    Kind        kind = Kind::None;
    bool        is_long = false;
    bool        volume = false;       // *ELEMENT_SPH_VOLUME
    bool        generate = false;     // *SET_NODE_LIST_GENERATE
    bool        set_title = false;    // *SET_NODE_LIST_TITLE
    int         card = 0;             // номер карточки внутри блока
    bool        title_done = false;   // *TITLE
    bool        inc_cont = false;     // *INCLUDE: ожидается продолжение имени
    std::string inc_acc;              // накопленное имя файла
    SrcLoc      inc_loc;
    SrcLoc      start;                // строка ключевого слова
    std::string name;                 // для сообщений
    // накопители
    long long   eid = 0, pid = 0;     // Solid10, карточка 1
    SrcLoc      card1_src;
    std::string heading;              // *PART
    long long   cur_sid = 0;          // текущий набор
};

// ---- Читатель ------------------------------------------------------------------

class Reader
{
public:
    explicit Reader(KDeck& d) : d_(d) {}

    void read_file(const fs::path& path, const SrcLoc* includer)
    {
        std::ifstream in(path, std::ios::binary);
        if (!in) fail_at(includer, "не удалось открыть файл " + path.string());
        std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        read_text(text, path);
    }

    void read_text(std::string_view text, const fs::path& path)
    {
        std::error_code ec;
        fs::path canon = fs::weakly_canonical(path, ec);
        if (ec) canon = path;
        if (stack_.empty()) main_dir_ = path.has_parent_path() ? path.parent_path() : fs::path(".");
        stack_.push_back(canon);

        const std::uint32_t fidx = static_cast<std::uint32_t>(d_.files.size());
        d_.files.push_back(path.string());
        const fs::path parent = path.has_parent_path() ? path.parent_path() : fs::path(".");

        Block blk;
        std::uint32_t lineno = 0;
        std::size_t pos = 0;
        bool ended = false;
        SrcLoc loc{fidx, 0};
        while (pos < text.size() && !ended)
        {
            std::size_t nl = text.find('\n', pos);
            std::size_t end = (nl == std::string_view::npos) ? text.size() : nl;
            std::string_view line = text.substr(pos, end - pos);
            pos = (nl == std::string_view::npos) ? text.size() : nl + 1;
            ++lineno;
            if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
            loc.line = lineno;

            if (!line.empty() && line[0] == '$') continue;       // комментарий
            if (!line.empty() && line[0] == '*')
            {
                finish_block(blk, loc);
                if (start_block(line, loc, blk)) ended = true;   // *END
                continue;
            }
            // пустая строка: значима только как заголовок *PART
            bool blank = trim(line).empty();
            if (blank && !(blk.kind == Kind::Part && blk.card == 0)) continue;
            data_line(blk, line, loc, parent);
        }
        if (!ended)
        {
            loc.line = lineno;
            finish_block(blk, loc);
        }
        else
        {
            // при *END незавершённость проверяется в start_block (finish_block вызван до него)
        }
        stack_.pop_back();
    }

private:
    KDeck&                 d_;
    bool                   global_long_ = false;
    std::vector<fs::path>  inc_paths_;
    fs::path               main_dir_;
    std::vector<fs::path>  stack_;

    // ---- Ошибки ----
    [[noreturn]] void fail(const SrcLoc& l, const std::string& msg) const
    {
        throw std::runtime_error(d_.where(l) + ": " + msg);
    }
    [[noreturn]] void fail_at(const SrcLoc* l, const std::string& msg) const
    {
        if (l) fail(*l, msg);
        throw std::runtime_error(msg);
    }

    // ---- Числа ----
    long long parse_int(std::string_view s, const char* fname, const SrcLoc& l) const
    {
        if (!s.empty() && s[0] == '&') fail(l, "*PARAMETER не поддерживается (поле " + std::string(fname) + ": '" + std::string(s) + "')");
        std::string_view t = s;
        if (!t.empty() && t[0] == '+') t.remove_prefix(1);
        long long v = 0;
        auto r = std::from_chars(t.data(), t.data() + t.size(), v);
        if (t.empty() || r.ec != std::errc() || r.ptr != t.data() + t.size())
            fail(l, "поле " + std::string(fname) + ": ожидалось целое, получено '" + std::string(s) + "'");
        return v;
    }

    double parse_real(std::string_view s, const char* fname, const SrcLoc& l) const
    {
        if (!s.empty() && s[0] == '&') fail(l, "*PARAMETER не поддерживается (поле " + std::string(fname) + ": '" + std::string(s) + "')");
        char buf[64];
        if (s.empty() || s.size() >= sizeof(buf))
            fail(l, "поле " + std::string(fname) + ": ожидалось число, получено '" + std::string(s) + "'");
        std::size_t n = 0;
        for (char c : s)
        {
            if (c == 'D' || c == 'd') c = 'E';
            if (!((c >= '0' && c <= '9') || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E'))
                fail(l, "поле " + std::string(fname) + ": ожидалось число, получено '" + std::string(s) + "'");
            buf[n++] = c;
        }
        const char* b = buf;
        const char* e = buf + n;
        if (*b == '+') ++b;
        double v = 0.0;
        auto r = std::from_chars(b, e, v);
        if (b == e || r.ec != std::errc() || r.ptr != e)
            fail(l, "поле " + std::string(fname) + ": ожидалось число, получено '" + std::string(s) + "'");
        return v;
    }

    // ---- Разбор строки данных по карточке ----
    // free, если есть запятая; иначе fixed по ширинам (long -> max(w,20)).
    void parse_card(std::string_view line, const CardDef& card, bool is_long,
                    const SrcLoc& l, Values& out) const
    {
        std::string_view toks[kMaxCardFields];
        const bool free_fmt = line.find(',') != std::string_view::npos;
        if (free_fmt)
        {
            std::size_t p = 0;
            for (int i = 0; i < card.count; ++i)
            {
                if (p > line.size()) { toks[i] = {}; continue; }
                std::size_t c = line.find(',', p);
                std::size_t e = (c == std::string_view::npos) ? line.size() : c;
                toks[i] = trim(line.substr(p, e - p));
                p = (c == std::string_view::npos) ? line.size() + 1 : c + 1;
            }
        }
        else
        {
            std::size_t off = 0;
            for (int i = 0; i < card.count; ++i)
            {
                std::size_t w = static_cast<std::size_t>(field_width(card.fields[i], is_long));
                if (off < line.size()) toks[i] = trim(line.substr(off, w));
                else toks[i] = {};
                off += w;
            }
        }
        for (int i = 0; i < card.count; ++i)
        {
            const FieldDef& f = card.fields[i];
            Value& v = out.v[static_cast<std::size_t>(i)];
            v = Value{};
            if (toks[i].empty())
            {
                if (f.required)
                    fail(l, std::string("обязательное поле ") + f.name + " пусто");
                v.i = static_cast<long long>(f.def);
                v.r = f.def;
                continue;
            }
            v.set = true;
            v.s = toks[i];
            if (f.type == FieldType::I)      { v.i = parse_int(toks[i], f.name, l); v.r = static_cast<double>(v.i); }
            else if (f.type == FieldType::R) { v.r = parse_real(toks[i], f.name, l); }
        }
    }

    // ---- Ключевое слово ----
    // Возвращает true при *END.
    bool start_block(std::string_view line, const SrcLoc& l, Block& blk)
    {
        blk = Block{};
        blk.start = l;
        std::string_view body = line.substr(1);
        std::size_t e = 0;
        while (e < body.size() && !std::isspace(static_cast<unsigned char>(body[e])) && body[e] != '(') ++e;
        std::string name = upper(body.substr(0, e));
        std::string rest = upper(trim(body.substr(e)));

        bool force_long = false, force_std = false;
        if (!name.empty() && name.back() == '+') { force_long = true; name.pop_back(); }
        else if (!name.empty() && name.back() == '-') { force_std = true; name.pop_back(); }
        blk.name = name;
        blk.is_long = force_long ? true : (force_std ? false : global_long_);

        if (name == "END") return true;
        if (name == "KEYWORD")
        {
            std::size_t p = rest.find("LONG=");
            if (p != std::string::npos && p + 5 < rest.size())
                global_long_ = (rest[p + 5] == 'Y');
            blk.kind = Kind::Keyword;
            return false;
        }
        if (name == "TITLE") { blk.kind = Kind::Title; return false; }
        if (name == "INCLUDE") { blk.kind = Kind::Include; return false; }
        if (name == "INCLUDE_PATH") { blk.kind = Kind::IncludePath; return false; }
        if (name == "INCLUDE_PATH_RELATIVE") { blk.kind = Kind::IncludePathRel; return false; }
        if (name == "INCLUDE_TRANSFORM")
            fail(l, "*INCLUDE_TRANSFORM не поддерживается (смещения ID и масштабы не применяются)");
        if (name == "NODE") { blk.kind = Kind::Node; return false; }
        if (name == "ELEMENT_SPH") { blk.kind = Kind::Sph; return false; }
        if (name == "ELEMENT_SPH_VOLUME") { blk.kind = Kind::Sph; blk.volume = true; return false; }
        if (starts_with(name, "ELEMENT_SPH_"))
            fail(l, "*" + name + ": опции *ELEMENT_SPH (кроме _VOLUME) не поддерживаются");
        if (name == "ELEMENT_BEAM") { blk.kind = Kind::Beam; return false; }
        if (starts_with(name, "ELEMENT_BEAM_"))
            fail(l, "*" + name + ": опции *ELEMENT_BEAM добавляют карточки и не поддерживаются");
        if (name == "ELEMENT_SOLID")
        {
            blk.kind = (rest.find("TEN NODES") != std::string::npos) ? Kind::Solid10 : Kind::Solid;
            return false;
        }
        if (starts_with(name, "ELEMENT_SOLID_"))
            fail(l, "*" + name + ": опции *ELEMENT_SOLID не поддерживаются");
        if (name == "PART") { blk.kind = Kind::Part; return false; }
        if (name == "BOUNDARY_SPC_NODE") { blk.kind = Kind::SpcNode; return false; }
        if (name == "BOUNDARY_SPC_SET") { blk.kind = Kind::SpcSet; return false; }
        if (starts_with(name, "SET_NODE_LIST"))
        {
            // допустимые опции: _TITLE, _GENERATE (в любом порядке)
            std::string opt = name.substr(13);
            bool ok = true, title = false, gen = false;
            std::size_t p = 0;
            while (ok && p < opt.size())
            {
                if (opt.compare(p, 6, "_TITLE") == 0) { title = true; p += 6; }
                else if (opt.compare(p, 9, "_GENERATE") == 0) { gen = true; p += 9; }
                else ok = false;
            }
            if (ok)
            {
                blk.kind = Kind::SetList;
                blk.set_title = title;
                blk.generate = gen;
                return false;
            }
        }
        // прочие ключевые слова - пропуск со сводкой
        blk.kind = Kind::Skip;
        ++d_.skipped[name];
        return false;
    }

    // Проверка завершённости блока при смене/конце файла.
    void finish_block(Block& blk, const SrcLoc& l) const
    {
        switch (blk.kind)
        {
        case Kind::Include:
            if (blk.inc_cont) fail(blk.inc_loc, "*INCLUDE: имя файла оборвано продолжением ' +'");
            break;
        case Kind::Part:
            if (blk.card == 1) fail(blk.start, "*PART: после заголовка нет карточки с PID");
            break;
        case Kind::Solid10:
            if (blk.card == 1) fail(blk.card1_src, "*ELEMENT_SOLID (ten nodes format): нет второй карточки с узлами");
            break;
        case Kind::SetList:
            if (blk.card == 0 && blk.set_title) fail(blk.start, "*" + blk.name + ": нет заголовка набора");
            break;
        default: break;
        }
        (void)l;
    }

    // ---- Строка данных ----
    void data_line(Block& blk, std::string_view line, const SrcLoc& l, const fs::path& parent)
    {
        Values vs;
        auto I = [&](int i) { return vs.v[static_cast<std::size_t>(i)].i; };
        auto R = [&](int i) { return vs.v[static_cast<std::size_t>(i)].r; };
        switch (blk.kind)
        {
        case Kind::None: case Kind::Skip: case Kind::Keyword:
            return;
        case Kind::Title:
            if (!blk.title_done)
            {
                d_.title = std::string(trim(line));
                blk.title_done = true;
            }
            return;
        case Kind::IncludePath: case Kind::IncludePathRel:
        {
            std::string_view t = trim(line);
            std::size_t p = 0;
            while (p <= t.size())
            {
                std::size_t c = t.find(';', p);
                std::size_t e = (c == std::string_view::npos) ? t.size() : c;
                std::string_view part = trim(t.substr(p, e - p));
                if (!part.empty())
                {
                    fs::path dir{std::string(part)};
                    if (blk.kind == Kind::IncludePathRel && dir.is_relative()) dir = main_dir_ / dir;
                    inc_paths_.push_back(dir);
                }
                if (c == std::string_view::npos) break;
                p = c + 1;
            }
            return;
        }
        case Kind::Include:
        {
            std::string_view t = trim(line);
            bool cont = t.size() >= 1 && t.back() == '+' && (t.size() == 1 || t[t.size() - 2] == ' ');
            if (cont) t = t.substr(0, t.size() - 1), t = t.substr(0, t.size() >= 1 ? t.size() - 1 : 0);
            if (!blk.inc_cont) { blk.inc_acc.clear(); blk.inc_loc = l; }
            blk.inc_acc.append(blk.inc_cont ? trim(t) : t);
            blk.inc_cont = cont;
            if (!cont) do_include(blk.inc_acc, blk.inc_loc, parent);
            return;
        }
        case Kind::Node:
        {
            parse_card(line, kNodeCard, blk.is_long, l, vs);
            d_.nodes.push_back(KNode{I(0), R(1), R(2), R(3), l});
            d_.node_index.emplace(I(0), d_.nodes.size() - 1);
            return;
        }
        case Kind::Sph:
        {
            parse_card(line, kSphCard, blk.is_long, l, vs);
            long long nid = I(0), nend = I(3);
            if (nend != 0 && nend < nid)
                fail(l, "*ELEMENT_SPH: NEND (" + std::to_string(nend) + ") < NID (" + std::to_string(nid) + ")");
            if (nend <= nid) nend = nid;
            for (long long n = nid; n <= nend; ++n)
                d_.sph.push_back(KSph{n, I(1), R(2), blk.volume, l});
            return;
        }
        case Kind::Beam:
        {
            parse_card(line, kBeamCard, blk.is_long, l, vs);
            d_.beams.push_back(KBeam{I(0), I(1), I(2), I(3), I(4), l});
            return;
        }
        case Kind::Solid:
        {
            parse_card(line, kSolidCard, blk.is_long, l, vs);
            KSolid s;
            s.eid = I(0); s.pid = I(1); s.src = l;
            bool ext = I(6) != 0 || I(7) != 0 || I(8) != 0 || I(9) != 0;
            s.nn = ext ? 8 : 4;
            for (int k = 0; k < s.nn; ++k) s.n[static_cast<std::size_t>(k)] = I(2 + k);
            d_.solids.push_back(s);
            return;
        }
        case Kind::Solid10:
        {
            if (blk.card == 0)
            {
                parse_card(line, kSolid10Card1, blk.is_long, l, vs);
                blk.eid = I(0); blk.pid = I(1); blk.card1_src = l; blk.card = 1;
            }
            else
            {
                parse_card(line, kSolid10Card2, blk.is_long, l, vs);
                KSolid s;
                s.eid = blk.eid; s.pid = blk.pid; s.nn = 10; s.src = blk.card1_src;
                for (int k = 0; k < 10; ++k) s.n[static_cast<std::size_t>(k)] = I(k);
                d_.solids.push_back(s);
                blk.card = 0;
            }
            return;
        }
        case Kind::Part:
        {
            if (blk.card == 0) { blk.heading = std::string(trim(line)); blk.card = 1; return; }
            if (blk.card == 1)
            {
                parse_card(line, kPartCard, blk.is_long, l, vs);
                d_.parts.push_back(KPart{I(0), blk.heading, l});
                blk.card = 2;
            }
            return;   // остальные карточки *PART пропускаются
        }
        case Kind::SpcNode:
            parse_card(line, kSpcNodeCard, blk.is_long, l, vs);
            d_.spc_nodes.push_back(KSpcNode{I(0), static_cast<int>(I(2)), static_cast<int>(I(3)), static_cast<int>(I(4)), l});
            return;
        case Kind::SpcSet:
            parse_card(line, kSpcSetCard, blk.is_long, l, vs);
            d_.spc_sets.push_back(KSpcSet{I(0), static_cast<int>(I(2)), static_cast<int>(I(3)), static_cast<int>(I(4)), l});
            return;
        case Kind::SetList:
        {
            if (blk.set_title && blk.card == 0) { blk.card = 1; return; }   // строка заголовка
            if (blk.card <= (blk.set_title ? 1 : 0))
            {
                parse_card(line, kSetHeadCard, blk.is_long, l, vs);
                blk.cur_sid = I(0);
                auto ins = d_.node_sets.emplace(blk.cur_sid, KNodeSet{{}, l});
                (void)ins;
                blk.card = blk.set_title ? 2 : 1;
                return;
            }
            std::vector<long long>& dst = d_.node_sets[blk.cur_sid].nids;
            if (blk.generate)
            {
                parse_card(line, kSetGenCard, blk.is_long, l, vs);
                for (int k = 0; k < 8; k += 2)
                {
                    long long b = I(k), e2 = I(k + 1);
                    if (b == 0 && e2 == 0) continue;
                    if (e2 < b) fail(l, "диапазон набора: конец (" + std::to_string(e2) + ") < начала (" + std::to_string(b) + ")");
                    for (long long n = b; n <= e2; ++n) dst.push_back(n);
                }
            }
            else
            {
                parse_card(line, kSetListCard, blk.is_long, l, vs);
                for (int k = 0; k < 8; ++k)
                    if (I(k) != 0) dst.push_back(I(k));
            }
            return;
        }
        }
    }

    // ---- *INCLUDE ----
    void do_include(const std::string& name, const SrcLoc& l, const fs::path& parent)
    {
        std::string_view nm = trim(name);
        if (nm.empty()) fail(l, "*INCLUDE: пустое имя файла");
        fs::path want{std::string(nm)};
        std::error_code ec;
        fs::path found;
        std::vector<fs::path> tried;
        auto try_path = [&](const fs::path& p) {
            tried.push_back(p);
            if (found.empty() && fs::is_regular_file(p, ec)) found = p;
        };
        if (want.is_absolute()) try_path(want);
        else
        {
            try_path(parent / want);
            for (const fs::path& dir : inc_paths_)
            {
                try_path(dir / want);
                if (dir.is_relative()) try_path(main_dir_ / dir / want);
            }
        }
        if (found.empty())
        {
            std::string msg = "*INCLUDE: файл '" + std::string(nm) + "' не найден (искали:";
            for (const auto& p : tried) msg += " " + p.string();
            fail(l, msg + ")");
        }
        fs::path canon = fs::weakly_canonical(found, ec);
        if (ec) canon = found;
        for (const fs::path& s : stack_)
            if (s == canon) fail(l, "*INCLUDE: цикл включений, файл '" + found.string() + "' уже читается");
        read_file(found, &l);
    }
};

} // namespace

KDeck read_k_deck_raw(const std::string& path)
{
    KDeck d;
    Reader r(d);
    r.read_file(fs::path(path), nullptr);
    return d;
}

KDeck read_k_deck(const std::string& path)
{
    KDeck d = read_k_deck_raw(path);
    validate_k_deck(d);
    return d;
}

KDeck read_k_deck_from_string(const std::string& text, const std::string& virtual_name, bool validate)
{
    KDeck d;
    Reader r(d);
    r.read_text(text, fs::path(virtual_name));
    if (validate) validate_k_deck(d);
    return d;
}

void validate_k_deck(const KDeck& d)
{
    auto fail = [&](const SrcLoc& l, const std::string& msg) {
        throw std::runtime_error(d.where(l) + ": " + msg);
    };

    // дубликаты NID в *NODE: индекс хранит первый
    for (std::size_t i = 0; i < d.nodes.size(); ++i)
    {
        auto it = d.node_index.find(d.nodes[i].nid);
        if (it != d.node_index.end() && it->second != i)
            fail(d.nodes[i].src, "*NODE: дубликат NID " + std::to_string(d.nodes[i].nid) +
                 " (первый раз: " + d.where(d.nodes[it->second].src) + ")");
    }

    std::unordered_map<long long, SrcLoc> parts;
    parts.reserve(d.parts.size());
    for (const KPart& p : d.parts)
    {
        auto ins = parts.emplace(p.pid, p.src);
        if (!ins.second)
            fail(p.src, "*PART: дубликат PID " + std::to_string(p.pid) + " (первый раз: " + d.where(ins.first->second) + ")");
    }

    std::unordered_map<long long, SrcLoc> sph_ids;
    sph_ids.reserve(d.sph.size() * 2);
    for (const KSph& s : d.sph)
    {
        auto ins = sph_ids.emplace(s.nid, s.src);
        if (!ins.second)
            fail(s.src, "*ELEMENT_SPH: дубликат NID " + std::to_string(s.nid) + " (первый раз: " + d.where(ins.first->second) + ")");
        if (!d.find_node(s.nid))
            fail(s.src, "*ELEMENT_SPH: узел " + std::to_string(s.nid) + " не найден в *NODE");
        if (!parts.count(s.pid))
            fail(s.src, "*ELEMENT_SPH: PID " + std::to_string(s.pid) + " не описан в *PART");
    }

    std::unordered_map<long long, SrcLoc> beam_ids;
    for (const KBeam& b : d.beams)
    {
        auto ins = beam_ids.emplace(b.eid, b.src);
        if (!ins.second)
            fail(b.src, "*ELEMENT_BEAM: дубликат EID " + std::to_string(b.eid) + " (первый раз: " + d.where(ins.first->second) + ")");
        for (long long n : {b.n1, b.n2})
        {
            if (!d.find_node(n))
                fail(b.src, "*ELEMENT_BEAM: узел " + std::to_string(n) + " не найден в *NODE");
            auto it = sph_ids.find(n);
            if (it != sph_ids.end())
                fail(b.src, "узел " + std::to_string(n) + " одновременно в *ELEMENT_BEAM и в *ELEMENT_SPH (" + d.where(it->second) + ")");
        }
        if (b.n3 != 0 && !d.find_node(b.n3))
            fail(b.src, "*ELEMENT_BEAM: узел ориентации " + std::to_string(b.n3) + " не найден в *NODE");
        if (!parts.count(b.pid))
            fail(b.src, "*ELEMENT_BEAM: PID " + std::to_string(b.pid) + " не описан в *PART");
    }

    std::unordered_map<long long, SrcLoc> solid_ids;
    for (const KSolid& s : d.solids)
    {
        auto ins = solid_ids.emplace(s.eid, s.src);
        if (!ins.second)
            fail(s.src, "*ELEMENT_SOLID: дубликат EID " + std::to_string(s.eid) + " (первый раз: " + d.where(ins.first->second) + ")");
        for (int k = 0; k < s.nn; ++k)
        {
            long long n = s.n[static_cast<std::size_t>(k)];
            if (!d.find_node(n))
                fail(s.src, "*ELEMENT_SOLID: узел " + std::to_string(n) + " не найден в *NODE");
        }
        if (!parts.count(s.pid))
            fail(s.src, "*ELEMENT_SOLID: PID " + std::to_string(s.pid) + " не описан в *PART");
    }

    for (const KSpcNode& c : d.spc_nodes)
        if (!d.find_node(c.nid))
            fail(c.src, "*BOUNDARY_SPC_NODE: узел " + std::to_string(c.nid) + " не найден в *NODE");
    for (const KSpcSet& c : d.spc_sets)
        if (!d.node_sets.count(c.nsid))
            fail(c.src, "*BOUNDARY_SPC_SET: набор узлов " + std::to_string(c.nsid) + " не найден (*SET_NODE_LIST)");
}

} // namespace k_input
