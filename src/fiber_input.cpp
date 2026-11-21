// fiber_input.cpp
// Разбор файла ломаных нитей в fiber::Model (notes/fiber-model/fiber-plan.md, этап Ф2, пункт 3).
// Объявление load_fiber_file — в include/fiber/fiber_input.hpp.
//
// ФОРМАТ ФАЙЛА (текстовый, свободный формат, разделители — пробелы):
//   - строки-комментарии и хвостовые комментарии начинаются с '#';
//   - пустые строки игнорируются;
//   - "MATERIAL <id> <E> <nu> <d> <rho> <eps_break> [mu_ff [mu_fs]]"
//       определяет материал нити с индексом <id> (целое >= 0; индексы не
//       обязаны быть подряд, но должны быть объявлены ДО первой ссылки на
//       них из FIBER). E, nu, d, rho, eps_break — double (модуль Юнга,
//       коэффициент Пуассона, диаметр, плотность, критическое удлинение
//       разрыва; см. notes/fiber-model/fiber-problem-statement.md, разделы Fiber-Bond и Fiber-Solid-Contact).
//       mu_ff — НЕОБЯЗАТЕЛЬНОЕ седьмое поле, коэффициент динамического
//       трения нить-нить (notes/fiber-contact/fiber-friction-implementation.md, §8);
//       если не задан, mu_ff = 0 (старое поведение без трения). Если задан,
//       должен быть числом >= 0, иначе — ошибка формата.
//       mu_fs — НЕОБЯЗАТЕЛЬНОЕ восьмое поле (только после mu_ff),
//       коэффициент трения нить-Solid (fiber/fiber_solid_friction_
//       implementation.md, §2); по умолчанию 0; правила те же, что у mu_ff.
//   - "FIBER <material_id> <point_count>"
//       открывает одну ломаную (= одну нить). За строкой заголовка следуют
//       РОВНО <point_count> строк "x y z" — вершины ломаной в порядке вдоль
//       нити, в ненагруженной (недеформированной) конфигурации. <point_count>
//       >= 2. Каждая ломаная становится отдельной Fiber с собственным
//       object_id (0, 1, 2, ... в порядке следования блоков FIBER в файле).
//   - "BOND <bead_i> <bead_j> <material_id>"
//       НЕОБЯЗАТЕЛЬНАЯ директива: связь-стежок между бусинами разных нитей
//       (notes/fabrics/plan-quilt-fabric.md, раздел 10). <bead_i>, <bead_j> -- глобальные
//       индексы бусин (с нуля, в порядке блоков FIBER и точек внутри них),
//       i != j; <material_id> -- материал стежка (E, d, eps_break), объявленный
//       директивой MATERIAL (в любом месте файла, до или после BOND). Связь
//       только растягивающая, l0 = расстояние между бусинами в начальной
//       конфигурации (> 0), l_break = (1 + eps_break) * l0. Стежок не меняет
//       prev/next и object_id. Масса бусин не меняется: масса швейной нити не
//       учитывается. Ошибка формата: индекс вне диапазона, i == j, неопределённый
//       материал, l0 == 0, лишние токены.
//   - "CLAMP_BOX <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> [<mask>]"
//       НЕОБЯЗАТЕЛЬНАЯ директива кинематического закрепления бусин в
//       прямоугольной области (interactions/clamp.hpp, notes/fiber-model/fiber-problem-statement.md,
//       раздел «Clamp (кинематическое закрепление)»). Координат столько,
//       сколько dim (в 3D — шесть чисел). <mask> — непустая строка из букв
//       x, y, z без повторов; по умолчанию все компоненты закреплены. Бусины,
//       попавшие в замкнутый бокс, закрепляются по указанным компонентам.
//       Несколько директив объединяются по «или». Может появиться в файле
//       в любом месте (до/после блоков FIBER).
//
//   - Неупругость связи (notes/fiber-model/fiber-plasticity-design.md, раздел 6; все
//     директивы необязательны, ключи КЛЮЧ=значение в любом порядке, регистр
//     ключей не важен; могут стоять до или после MATERIAL/FIBER, материал с
//     этим id должен быть объявлен где-нибудь в файле):
//       PLASTIC <id> SY=.. [H=..] | SY=.. K=.. N=..
//       PLASTIC_TABLE <id> <N>  + N строк "eps_p sigma_s" (eps_p от 0, строго
//                               возрастает; sigma_s > 0, не убывает)
//       RATE <id> CS C=.. P=.. | JC C=.. EDOT0=.. | PERZYNA ETA=.. M=..
//       VISCO <id> [ETA=..] [XI=..] [E1=.. TAU=..]
//       FAILURE <id> [EPS_P=..] [SIGMA_U=..]
//     Правила и диапазоны -- docs/03-reference/file-formats/fib-file.md.
//
// Пример:
//   # Нить из капрона, две прямые ломаные по 3 и 2 точки, концы закреплены
//   MATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.02
//   FIBER 0 3
//   0.0 0.0 0.0
//   0.01 0.0 0.0
//   0.02 0.0 0.0
//   CLAMP_BOX 0.0 -0.01 -0.01 0.0 0.01 0.01
//   FIBER 0 2
//   0.0 0.01 0.0
//   0.0 0.02 0.0
//   CLAMP_BOX 0.0 0.01 -0.01 0.0 0.02 0.01 xy
//
// Вычисления при разборе (notes/fiber-model/fiber-problem-statement.md, раздел «Fiber»):
//   l0 связи (k, k+1) = расстояние между точками файла (как заданы, без
//   деформации);
//   l_break = (1 + eps_break) * l0 (материал связи; notes/fiber-model/fiber-problem-statement.md,
//   Fiber-Bond, «Постановка проекта: только растяжение»);
//   масса бусины k: m_k = rho_f * A * (l_{k-1,k} + l_{k,k+1}) / 2,
//   A = pi*d^2/4, для концевой бусины берётся только один прилежащий
//   сегмент (второе слагаемое отсутствует, а не делится пополам дважды).
#include "fiber/fiber_input.hpp"

#include <array>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "fiber/fiber_build.hpp"
#include "interactions/clamp.hpp"
#include "param.hpp"
#include "rhs/cpu_executor.hpp"

namespace fiber
{

namespace
{

// Обрезает комментарий, начинающийся с '#'. Возвращает true, если после
// обрезки строка непустая (содержит не-пробельные символы).
bool strip_comment_nonblank(std::string& line)
{
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    const std::size_t b = line.find_first_not_of(" \t\r\n");
    return b != std::string::npos;
}

// ----------------------------------------------------------------------
// Директивы неупругости связи (notes/fiber-model/fiber-plasticity-design.md, раздел 6):
// PLASTIC, PLASTIC_TABLE, RATE, VISCO, FAILURE. Разбираются в промежуточную
// структуру по id материала и применяются к материалам ПОСЛЕ чтения всего
// файла (директива может стоять и до, и после MATERIAL/FIBER; повторный
// MATERIAL с тем же id не затирает неупругие параметры).
// ----------------------------------------------------------------------
struct PendingInelastic
{
    std::size_t line_plastic = 0; // номер строки директивы; 0 = не было
    std::size_t line_table = 0;
    std::size_t line_rate = 0;
    std::size_t line_visco = 0;
    std::size_t line_failure = 0;
    std::size_t line_damage = 0;
    std::size_t first_line = 0; // строка первой директивы (для ошибки "нет материала")
    PlasticParams plastic{};
    RateParams rate{};
    ViscoParams visco{};
    FailureParams failure{};
    DamageParams damage{};
};

[[noreturn]] void inelastic_fail(const char* dir, std::size_t line_no, const std::string& msg)
{
    throw std::runtime_error(std::string("fiber::load_fiber_file: ") + dir + ": " + msg + " в строке " +
                             std::to_string(line_no));
}

std::string upper_copy(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// Целое число из токена целиком.
bool parse_whole_ll(const std::string& tok, long long& v)
{
    std::istringstream t(tok);
    return static_cast<bool>(t >> v) && (t >> std::ws).eof();
}

// Конечное double из токена целиком (поток с классической локалью).
bool parse_whole_double(const std::string& tok, double& v)
{
    std::istringstream t(tok);
    return static_cast<bool>(t >> v) && (t >> std::ws).eof() && std::isfinite(v);
}

// Читает id материала (целое >= 0) из потока.
std::size_t read_material_id(std::istringstream& iss, const char* dir, std::size_t line_no)
{
    std::string tok;
    long long id = -1;
    if (!(iss >> tok) || !parse_whole_ll(tok, id) || id < 0)
        inelastic_fail(dir, line_no, "ожидался id материала (целое >= 0)");
    return static_cast<std::size_t>(id);
}

using KeyValues = std::map<std::string, double>;

// Читает до конца строки токены КЛЮЧ=значение (ключ без учёта регистра).
// Неизвестный ключ, повтор, токен без '=', нечисловое/неконечное значение --
// std::runtime_error.
KeyValues read_key_values(std::istringstream& iss, const char* dir, std::size_t line_no,
                          const std::vector<std::string>& allowed)
{
    KeyValues kv;
    std::string tok;
    while (iss >> tok)
    {
        const std::size_t eq = tok.find('=');
        if (eq == std::string::npos || eq == 0 || eq + 1 == tok.size())
            inelastic_fail(dir, line_no, "ожидалось КЛЮЧ=значение, получено '" + tok + "'");
        const std::string key = upper_copy(tok.substr(0, eq));
        bool ok = false;
        for (const auto& a : allowed) ok = ok || (a == key);
        if (!ok) inelastic_fail(dir, line_no, "неизвестный ключ '" + key + "'");
        double v = 0.0;
        if (!parse_whole_double(tok.substr(eq + 1), v))
            inelastic_fail(dir, line_no,
                           "значение ключа " + key + " не конечное число: '" + tok.substr(eq + 1) + "'");
        if (!kv.emplace(key, v).second) inelastic_fail(dir, line_no, "ключ " + key + " повторяется");
    }
    return kv;
}

double require_key(const KeyValues& kv, const char* key, const char* dir, std::size_t line_no)
{
    const auto it = kv.find(key);
    if (it == kv.end()) inelastic_fail(dir, line_no, std::string("не задан обязательный ключ ") + key);
    return it->second;
}

bool opt_key(const KeyValues& kv, const char* key, double& v)
{
    const auto it = kv.find(key);
    if (it == kv.end()) return false;
    v = it->second;
    return true;
}

std::string repeat_msg(std::size_t first_line)
{
    return "повторная директива для этого материала (первая в строке " + std::to_string(first_line) + ")";
}

void parse_plastic_directive(std::istringstream& iss, std::size_t line_no, PendingInelastic& pi)
{
    static const char* dir = "PLASTIC";
    const KeyValues kv = read_key_values(iss, dir, line_no, {"SY", "H", "K", "N"});
    if (pi.line_plastic) inelastic_fail(dir, line_no, repeat_msg(pi.line_plastic));
    if (pi.line_table) inelastic_fail(dir, line_no, "PLASTIC и PLASTIC_TABLE для одного материала несовместимы");
    const double sy = require_key(kv, "SY", dir, line_no);
    if (!(sy > 0.0)) inelastic_fail(dir, line_no, "SY должен быть > 0");
    double H = 0.0, K = 0.0, n = 1.0;
    const bool hasH = opt_key(kv, "H", H);
    const bool hasK = opt_key(kv, "K", K);
    const bool hasN = opt_key(kv, "N", n);
    if (hasH && (hasK || hasN)) inelastic_fail(dir, line_no, "линейное (H) и степенное (K, N) упрочнение несовместимы");
    if (hasK != hasN) inelastic_fail(dir, line_no, "степенное упрочнение требует оба ключа K и N");
    if (hasH && !(H >= 0.0)) inelastic_fail(dir, line_no, "H должен быть >= 0");
    PlasticParams p;
    p.sy = sy;
    if (hasK)
    {
        if (!(K >= 0.0)) inelastic_fail(dir, line_no, "K должен быть >= 0");
        if (!(n > 0.0 && n <= 1.0)) inelastic_fail(dir, line_no, "N должен быть в (0, 1]");
        p.type = HardeningType::Power;
        p.K = K;
        p.n = n;
    }
    else
    {
        p.type = HardeningType::Linear;
        p.H = H;
    }
    pi.plastic = p;
    pi.line_plastic = line_no;
}

void parse_plastic_table_directive(std::istringstream& iss, std::ifstream& in, std::size_t& line_no,
                                   PendingInelastic& pi)
{
    static const char* dir = "PLASTIC_TABLE";
    const std::size_t head_line = line_no;
    std::string tok;
    long long cnt = 0;
    if (!(iss >> tok) || !parse_whole_ll(tok, cnt)) inelastic_fail(dir, head_line, "ожидалось число строк таблицы N");
    if (cnt < 2) inelastic_fail(dir, head_line, "число строк таблицы должно быть >= 2");
    if (iss >> tok) inelastic_fail(dir, head_line, "лишний токен '" + tok + "' после числа строк");
    if (pi.line_table) inelastic_fail(dir, head_line, repeat_msg(pi.line_table));
    if (pi.line_plastic) inelastic_fail(dir, head_line, "PLASTIC и PLASTIC_TABLE для одного материала несовместимы");

    PlasticParams p;
    p.type = HardeningType::Table;
    std::string line;
    for (long long r = 0; r < cnt; ++r)
    {
        bool got = false;
        while (std::getline(in, line))
        {
            ++line_no;
            if (strip_comment_nonblank(line))
            {
                got = true;
                break;
            }
        }
        if (!got)
            inelastic_fail(dir, line_no,
                           "неожиданный конец файла: ожидалось " + std::to_string(cnt) + " строк таблицы, прочитано " +
                               std::to_string(r));
        std::istringstream row(line);
        std::string a, b, extra;
        double e = 0.0, s = 0.0;
        if (!(row >> a >> b) || !parse_whole_double(a, e) || !parse_whole_double(b, s))
            inelastic_fail(dir, line_no,
                           "строка таблицы " + std::to_string(r + 1) + " из " + std::to_string(cnt) +
                               ": ожидалось 'eps_p sigma_s'");
        if (row >> extra) inelastic_fail(dir, line_no, "лишний токен '" + extra + "' в строке таблицы");
        if (r == 0 && e != 0.0) inelastic_fail(dir, line_no, "первая точка таблицы должна иметь eps_p = 0");
        if (!(s > 0.0)) inelastic_fail(dir, line_no, "sigma_s должен быть > 0");
        if (r > 0)
        {
            if (!(e > p.table_eps.back())) inelastic_fail(dir, line_no, "eps_p должен строго возрастать");
            if (!(s >= p.table_sigma.back())) inelastic_fail(dir, line_no, "sigma_s не должен убывать");
        }
        p.table_eps.push_back(e);
        p.table_sigma.push_back(s);
    }
    p.sy = p.table_sigma.front();
    pi.plastic = p;
    pi.line_table = head_line;
}

void parse_rate_directive(std::istringstream& iss, std::size_t line_no, PendingInelastic& pi)
{
    static const char* dir = "RATE";
    std::string model_tok;
    if (!(iss >> model_tok)) inelastic_fail(dir, line_no, "ожидалась модель CS, JC или PERZYNA");
    const std::string model = upper_copy(model_tok);
    RateParams r;
    if (model == "CS")
    {
        const KeyValues kv = read_key_values(iss, dir, line_no, {"C", "P"});
        r.model = RateModel::CowperSymonds;
        r.a = require_key(kv, "C", dir, line_no);
        r.b = require_key(kv, "P", dir, line_no);
        if (!(r.a > 0.0)) inelastic_fail(dir, line_no, "C должен быть > 0");
        if (!(r.b >= 1.0)) inelastic_fail(dir, line_no, "P должен быть >= 1");
    }
    else if (model == "JC")
    {
        const KeyValues kv = read_key_values(iss, dir, line_no, {"C", "EDOT0"});
        r.model = RateModel::JohnsonCook;
        r.a = require_key(kv, "C", dir, line_no);
        r.b = require_key(kv, "EDOT0", dir, line_no);
        if (!(r.a > 0.0)) inelastic_fail(dir, line_no, "C должен быть > 0");
        if (!(r.b > 0.0)) inelastic_fail(dir, line_no, "EDOT0 должен быть > 0");
    }
    else if (model == "PERZYNA")
    {
        const KeyValues kv = read_key_values(iss, dir, line_no, {"ETA", "M"});
        r.model = RateModel::Perzyna;
        r.a = require_key(kv, "ETA", dir, line_no);
        r.b = require_key(kv, "M", dir, line_no);
        if (!(r.a > 0.0)) inelastic_fail(dir, line_no, "ETA должен быть > 0");
        if (!(r.b > 0.0)) inelastic_fail(dir, line_no, "M должен быть > 0");
    }
    else
    {
        inelastic_fail(dir, line_no, "неизвестная модель '" + model_tok + "' (допустимы CS, JC, PERZYNA)");
    }
    if (pi.line_rate) inelastic_fail(dir, line_no, repeat_msg(pi.line_rate));
    pi.rate = r;
    pi.line_rate = line_no;
}

void parse_visco_directive(std::istringstream& iss, std::size_t line_no, PendingInelastic& pi)
{
    static const char* dir = "VISCO";
    const KeyValues kv = read_key_values(iss, dir, line_no, {"ETA", "XI", "E1", "TAU"});
    if (kv.empty()) inelastic_fail(dir, line_no, "не задан ни один ключ (ETA, XI, E1, TAU)");
    if (pi.line_visco) inelastic_fail(dir, line_no, repeat_msg(pi.line_visco));
    ViscoParams v;
    const bool hasEta = opt_key(kv, "ETA", v.eta);
    const bool hasXi = opt_key(kv, "XI", v.xi);
    const bool hasE1 = opt_key(kv, "E1", v.E1);
    const bool hasTau = opt_key(kv, "TAU", v.tau);
    if (hasEta && !(v.eta >= 0.0)) inelastic_fail(dir, line_no, "ETA должен быть >= 0");
    if (hasXi && !(v.xi >= 0.0)) inelastic_fail(dir, line_no, "XI должен быть >= 0");
    if (hasE1 && !(v.E1 >= 0.0)) inelastic_fail(dir, line_no, "E1 должен быть >= 0");
    if (hasTau && !(v.tau > 0.0)) inelastic_fail(dir, line_no, "TAU должен быть > 0");
    if (v.eta > 0.0 && v.xi > 0.0) inelastic_fail(dir, line_no, "ETA и XI одновременно задавать нельзя (выберите одно)");
    if (v.E1 > 0.0 && !hasTau) inelastic_fail(dir, line_no, "E1 > 0 требует TAU");
    if (hasTau && !(v.E1 > 0.0)) inelastic_fail(dir, line_no, "TAU задан без E1 > 0");
    pi.visco = v;
    pi.line_visco = line_no;
}

void parse_failure_directive(std::istringstream& iss, std::size_t line_no, PendingInelastic& pi)
{
    static const char* dir = "FAILURE";
    const KeyValues kv = read_key_values(iss, dir, line_no, {"EPS_P", "SIGMA_U"});
    if (kv.empty()) inelastic_fail(dir, line_no, "не задан ни один ключ (EPS_P, SIGMA_U)");
    if (pi.line_failure) inelastic_fail(dir, line_no, repeat_msg(pi.line_failure));
    FailureParams f;
    if (opt_key(kv, "EPS_P", f.eps_p) && !(f.eps_p > 0.0)) inelastic_fail(dir, line_no, "EPS_P должен быть > 0");
    if (opt_key(kv, "SIGMA_U", f.sigma_u) && !(f.sigma_u > 0.0))
        inelastic_fail(dir, line_no, "SIGMA_U должен быть > 0");
    pi.failure = f;
    pi.line_failure = line_no;
}

// DAMAGE <id> DEPS_F=<x> | GF=<Дж/м²> -- повреждённость с разупрочнением после σ_u
// (notes/fiber-model/fiber-damage-plan.md). Ровно один из ключей; требует FAILURE SIGMA_U.
// THETA_C (угловой критерий) зарезервирован, пока не реализован -- ошибка, а не молчание.
void parse_damage_directive(std::istringstream& iss, std::size_t line_no, PendingInelastic& pi)
{
    static const char* dir = "DAMAGE";
    const KeyValues kv = read_key_values(iss, dir, line_no, {"DEPS_F", "GF", "THETA_C"});
    if (pi.line_damage) inelastic_fail(dir, line_no, repeat_msg(pi.line_damage));
    if (kv.count("THETA_C")) inelastic_fail(dir, line_no, "THETA_C (угловой критерий) пока не реализован");
    DamageParams d;
    const bool hasDeps = opt_key(kv, "DEPS_F", d.deps_f);
    const bool hasGf = opt_key(kv, "GF", d.gf);
    if (hasDeps == hasGf) inelastic_fail(dir, line_no, "нужен ровно один ключ: DEPS_F или GF");
    if (hasDeps && !(d.deps_f > 0.0)) inelastic_fail(dir, line_no, "DEPS_F должен быть > 0");
    if (hasGf && !(d.gf > 0.0)) inelastic_fail(dir, line_no, "GF должен быть > 0");
    pi.damage = d;
    pi.line_damage = line_no;
}

// Проверки "после всего файла" и перенос параметров в материалы.
void apply_inelastic(const std::map<std::size_t, PendingInelastic>& pending, std::vector<FiberMaterial>& materials,
                     const std::vector<char>& material_defined)
{
    for (const auto& entry : pending)
    {
        const std::size_t id = entry.first;
        const PendingInelastic& pi = entry.second;
        if (id >= materials.size() || !material_defined[id])
        {
            throw std::runtime_error("fiber::load_fiber_file: директива неупругости (строка " +
                                     std::to_string(pi.first_line) + ") ссылается на неопределённый материал " +
                                     std::to_string(id));
        }
        if (pi.line_rate && !pi.line_plastic && !pi.line_table)
        {
            throw std::runtime_error("fiber::load_fiber_file: RATE (строка " + std::to_string(pi.line_rate) +
                                     ") для материала " + std::to_string(id) + " требует PLASTIC или PLASTIC_TABLE");
        }
        // FAILURE EPS_P проверяет eps_p связи; без PLASTIC/PLASTIC_TABLE eps_p == 0
        // и критерий никогда не сработает. Это не ошибка -- только предупреждение.
        if (pi.line_failure && pi.failure.eps_p > 0.0 && !pi.line_plastic && !pi.line_table)
        {
            std::cerr << "Предупреждение: FAILURE EPS_P (строка " << pi.line_failure << ") для материала "
                      << id << " задан без PLASTIC/PLASTIC_TABLE: критерий EPS_P никогда не сработает\n";
        }
        if (pi.line_damage && !(pi.line_failure && pi.failure.sigma_u > 0.0))
        {
            throw std::runtime_error("fiber::load_fiber_file: DAMAGE (строка " + std::to_string(pi.line_damage) +
                                     ") для материала " + std::to_string(id) +
                                     " требует FAILURE SIGMA_U (напряжение начала повреждения)");
        }
        FiberMaterial& m = materials[id];
        if (pi.line_plastic || pi.line_table) m.plastic = pi.plastic;
        if (pi.line_rate) m.rate = pi.rate;
        if (pi.line_visco) m.visco = pi.visco;
        if (pi.line_failure) m.failure = pi.failure;
        if (pi.line_damage) m.damage = pi.damage;
    }
}

} // namespace

void parse_material_directive(std::istringstream& iss, std::size_t line_no, const std::string& who,
                              std::vector<FiberMaterial>& materials, std::vector<char>& material_defined)
{
    long long id = -1;
    double E = 0.0, nu = 0.0, d = 0.0, rho = 0.0, eps_break = 0.0;
    if (!(iss >> id >> E >> nu >> d >> rho >> eps_break) || id < 0)
    {
        throw std::runtime_error(
            who + ": неверный формат MATERIAL в строке " + std::to_string(line_no));
    }

    // mu_ff (коэффициент динамического трения нить-нить,
    // notes/fiber-contact/fiber-friction-implementation.md, §8) — НЕОБЯЗАТЕЛЬНОЕ
    // седьмое поле строки MATERIAL. Если после шести обязательных
    // полей в строке больше нет токена, оставляем mu_ff = 0.0
    // (старое поведение без трения, ничего не меняющее для уже
    // существующих файлов нитей). Если токен есть, он обязан
    // разбираться как double ЦЕЛИКОМ (иначе — ошибка формата, как и
    // для остальных полей MATERIAL выше) и быть неотрицательным.
    // Разбор — потоком, как у остальных полей (локаль потока
    // "classic"), а не std::stod, зависящим от глобальной C-локали.
    // mu_fs (коэффициент трения нить-Solid, fiber/
    // fiber-solid-friction-implementation.md, §2) -- НЕОБЯЗАТЕЛЬНОЕ
    // восьмое поле (только после mu_ff), те же правила разбора и
    // проверок; по умолчанию 0.0.
    auto read_optional_mu = [&](const char* name, double& value) {
        std::string tok_str;
        if (!(iss >> tok_str)) return;
        std::istringstream tok(tok_str);
        const bool ok = static_cast<bool>(tok >> value) && (tok >> std::ws).eof();
        if (!ok)
        {
            throw std::runtime_error(who + ": неверный формат MATERIAL (" + name +
                                     ") в строке " + std::to_string(line_no));
        }
        if (value < 0.0)
        {
            throw std::runtime_error(who + ": " + std::string(name) +
                                     " должен быть >= 0 в MATERIAL, строка " + std::to_string(line_no));
        }
    };
    double mu_ff = 0.0;
    double mu_fs = 0.0;
    read_optional_mu("mu_ff", mu_ff);
    read_optional_mu("mu_fs", mu_fs);

    const std::size_t uid = static_cast<std::size_t>(id);
    if (uid >= materials.size())
    {
        materials.resize(uid + 1);
        material_defined.resize(uid + 1, 0);
    }
    materials[uid] = FiberMaterial{E, nu, d, rho, eps_break, mu_ff, mu_fs};
    material_defined[uid] = 1;
}

void build_model_from_polylines(const std::vector<PolylineSpec>& polylines,
                                const std::vector<FiberMaterial>& materials,
                                const std::vector<ClampBoxSpec>& boxes, Model& model, const std::string& who)
{
    for (std::size_t fid = 0; fid < polylines.size(); ++fid)
    {
        const PolylineSpec& f = polylines[fid];
        if (f.material_id < 0 || static_cast<std::size_t>(f.material_id) >= materials.size())
        {
            throw std::runtime_error("fiber::build_model_from_polylines: нить " + std::to_string(fid) +
                                     " ссылается на неопределённый материал " + std::to_string(f.material_id));
        }
        if (f.points.size() < 2)
        {
            throw std::runtime_error("fiber::build_model_from_polylines: нить " + std::to_string(fid) +
                                     " содержит меньше двух точек");
        }
    }

    std::size_t total_n = 0;
    for (const auto& f : polylines) total_n += f.points.size();

    model.materials = materials;
    model.fiber.resize(total_n);
    model.topology.resize(total_n);
    model.topology.bonds.clear();
    model.topology.bonds.reserve(total_n);
    ++model.topology.version;

    std::size_t base = 0;
    for (std::size_t fid = 0; fid < polylines.size(); ++fid)
    {
        const PolylineSpec& f = polylines[fid];
        const std::size_t n = f.points.size();
        const FiberMaterial& mat = model.materials[static_cast<std::size_t>(f.material_id)];
        const double A = pi * mat.d * mat.d / 4.0;

        std::vector<double> seg(n - 1);
        for (std::size_t k = 0; k + 1 < n; ++k)
        {
            double d2 = 0.0;
            for (std::size_t c = 0; c < 3; ++c)
            {
                const double dc = f.points[k + 1][c] - f.points[k][c];
                d2 += dc * dc;
            }
            seg[k] = std::sqrt(d2);
            if (!(seg[k] > 0.0))
            {
                throw std::runtime_error(
                    who + ": нулевая или отрицательная длина сегмента в нити " +
                    std::to_string(fid) + " (точки " + std::to_string(k) + ", " + std::to_string(k + 1) + ")");
            }
        }

        for (std::size_t k = 0; k < n; ++k)
        {
            const std::size_t gi = base + k;
            for (std::size_t d = 0; d < dim; ++d)
            {
                model.fiber.x[d][gi] = f.points[k][d];
                model.fiber.vx[d][gi] = 0.0;
            }
            const double l_prev = (k > 0) ? seg[k - 1] : 0.0;
            const double l_next = (k + 1 < n) ? seg[k] : 0.0;
            model.fiber.mass[gi] = mat.rho * A * (l_prev + l_next) / 2.0;
            model.fiber.object_id[gi] = static_cast<int>(fid);
            model.fiber.material_id[gi] = f.material_id;

            model.topology.prev[gi] = (k > 0) ? static_cast<long long>(gi - 1) : -1;
            model.topology.next[gi] = (k + 1 < n) ? static_cast<long long>(gi + 1) : -1;
        }

        for (std::size_t k = 0; k + 1 < n; ++k)
        {
            FiberBond b;
            b.i = base + k;
            b.j = base + k + 1;
            b.l0 = seg[k];
            b.l_break = (1.0 + mat.eps_break) * seg[k];
            model.topology.bonds.push_back(b);
        }

        base += n;
    }

    // Применяем все собранные boxes к clamps.fiber_mask
    if (!boxes.empty())
    {
        rhs::CpuExecutor exec(1);
        model.clamps.fiber_mask.assign(total_n, 0);
        for (const auto& box : boxes)
        {
            interactions::mark_clamp_box(exec, model.fiber.x, total_n, box.lo, box.hi, box.mask,
                                        model.clamps.fiber_mask);
        }
    }
    else
    {
        // Нет директив CLAMP_BOX -> оставляем fiber_mask пустым (не заполненным нулями)
        model.clamps.fiber_mask.clear();
    }
}

// BOND <i> <j> <mat>: связь-стежок, строка файла хранится для сообщений об ошибках.
struct StitchSpec
{
    long long i = 0;
    long long j = 0;
    long long mat = 0;
    std::size_t line = 0;
};

void load_fiber_file(const std::string& filename, Model& model)
{
    // Границу Р9 проверяем прежде всего: дальше код пишет ровно dim
    // компонент в FiberState::x[dim]/vx[dim] из троек (x,y,z) файла, и это
    // корректно только при dim == 3 (иначе часть заданных координат молча
    // терялась бы).
    check_fiber_requires_dim3();

    std::ifstream in(filename);
    if (!in.is_open())
    {
        throw std::runtime_error("fiber::load_fiber_file: не удалось открыть файл '" + filename + "'");
    }

    std::vector<FiberMaterial> materials;
    std::vector<char> material_defined;
    std::vector<PolylineSpec> fibers;
    std::vector<ClampBoxSpec> boxes;
    std::vector<StitchSpec> stitches;
    std::map<std::size_t, PendingInelastic> pending; // PLASTIC/RATE/VISCO/FAILURE по id

    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line))
    {
        ++line_no;
        if (!strip_comment_nonblank(line)) continue;

        std::istringstream iss(line);
        std::string kw;
        iss >> kw;

        if (kw == "MATERIAL")
        {
            parse_material_directive(iss, line_no, "fiber::load_fiber_file", materials, material_defined);
        }
        else if (kw == "FIBER")
        {
            long long material_id = -1;
            long long point_count = 0;
            if (!(iss >> material_id >> point_count) || material_id < 0 || point_count < 2)
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: неверный заголовок FIBER (материал/число точек) в строке " +
                    std::to_string(line_no));
            }
            const std::size_t umat = static_cast<std::size_t>(material_id);
            if (umat >= materials.size() || !material_defined[umat])
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: FIBER ссылается на неопределённый материал " +
                    std::to_string(material_id) + " в строке " + std::to_string(line_no));
            }

            PolylineSpec rf;
            rf.material_id = material_id;
            rf.points.reserve(static_cast<std::size_t>(point_count));
            for (long long p = 0; p < point_count; ++p)
            {
                if (!std::getline(in, line))
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: неожиданный конец файла внутри FIBER (после строки " +
                        std::to_string(line_no) + ")");
                }
                ++line_no;
                strip_comment_nonblank(line);
                std::istringstream piss(line);
                std::array<double, 3> pt{};
                if (!(piss >> pt[0] >> pt[1] >> pt[2]))
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: неверный формат точки в строке " + std::to_string(line_no));
                }
                rf.points.push_back(pt);
            }
            fibers.push_back(std::move(rf));
        }
        else if (kw == "PLASTIC" || kw == "PLASTIC_TABLE" || kw == "RATE" || kw == "VISCO" || kw == "FAILURE" ||
                 kw == "DAMAGE")
        {
            const std::size_t mid = read_material_id(iss, kw.c_str(), line_no);
            PendingInelastic& pi = pending[mid];
            if (pi.first_line == 0) pi.first_line = line_no;
            if (kw == "PLASTIC") parse_plastic_directive(iss, line_no, pi);
            else if (kw == "PLASTIC_TABLE") parse_plastic_table_directive(iss, in, line_no, pi);
            else if (kw == "RATE") parse_rate_directive(iss, line_no, pi);
            else if (kw == "VISCO") parse_visco_directive(iss, line_no, pi);
            else if (kw == "DAMAGE") parse_damage_directive(iss, line_no, pi);
            else parse_failure_directive(iss, line_no, pi);
        }
        else if (kw == "BOND")
        {
            // BOND <bead_i> <bead_j> <material_id>: связь-стежок между бусинами
            // разных нитей (диапазоны индексов проверяются после разбора всех FIBER).
            StitchSpec st;
            st.line = line_no;
            if (!(iss >> st.i >> st.j >> st.mat))
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: неверный формат BOND (ожидается: bead_i bead_j material_id) в строке " +
                    std::to_string(line_no));
            }
            std::string extra;
            if (iss >> extra)
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: лишние токены после BOND в строке " + std::to_string(line_no));
            }
            stitches.push_back(st);
        }
        else if (kw == "CLAMP_BOX")
        {
            ClampBoxSpec box;
            for (std::size_t d = 0; d < dim; ++d)
            {
                if (!(iss >> box.lo[d]))
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: не удалось прочитать минимальную координату CLAMP_BOX в строке " +
                        std::to_string(line_no));
                }
            }
            for (std::size_t d = 0; d < dim; ++d)
            {
                if (!(iss >> box.hi[d]))
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: не удалось прочитать максимальную координату CLAMP_BOX в строке " +
                        std::to_string(line_no));
                }
            }

            // Проверяем, что координаты конечны (NaN прошёл бы проверку min <= max
            // и закрепил бы всё) и xmin <= xmax для всех компонент
            for (std::size_t d = 0; d < dim; ++d)
            {
                if (!std::isfinite(box.lo[d]) || !std::isfinite(box.hi[d]))
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: неконечная координата CLAMP_BOX в строке " + std::to_string(line_no));
                }
                if (box.lo[d] > box.hi[d])
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: в CLAMP_BOX координата min[" + std::to_string(d) +
                        "] > max[" + std::to_string(d) + "] в строке " + std::to_string(line_no));
                }
            }

            // Читаем опциональную маску
            std::string mask_str;
            if (iss >> mask_str)
            {
                // Маска не пустая и содержит только x, y, z без повторов
                if (mask_str.empty())
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: маска CLAMP_BOX не должна быть пустой в строке " +
                        std::to_string(line_no));
                }

                box.mask = 0;
                std::vector<bool> seen(dim, false);
                for (char c : mask_str)
                {
                    std::size_t component = 0;
                    if (c == 'x') component = 0;
                    else if (c == 'y') component = 1;
                    else if (c == 'z') component = 2;
                    else
                    {
                        throw std::runtime_error(
                            "fiber::load_fiber_file: недопустимая буква '" + std::string(1, c) +
                            "' в маске CLAMP_BOX в строке " + std::to_string(line_no));
                    }

                    if (component >= dim)
                    {
                        throw std::runtime_error(
                            "fiber::load_fiber_file: буква '" + std::string(1, c) + "' недопустима для dim=" +
                            std::to_string(dim) + " в строке " + std::to_string(line_no));
                    }

                    if (seen[component])
                    {
                        throw std::runtime_error(
                            "fiber::load_fiber_file: повтор буквы '" + std::string(1, c) +
                            "' в маске CLAMP_BOX в строке " + std::to_string(line_no));
                    }
                    seen[component] = true;
                    box.mask |= static_cast<std::uint8_t>(1u << component);
                }
            }

            // Проверяем, что нет лишних токенов после маски
            std::string extra;
            if (iss >> extra)
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: лишние токены после CLAMP_BOX в строке " + std::to_string(line_no));
            }

            boxes.push_back(box);
        }
        else
        {
            throw std::runtime_error(
                "fiber::load_fiber_file: неизвестное ключевое слово '" + kw + "' в строке " + std::to_string(line_no));
        }
    }

    if (fibers.empty())
    {
        throw std::runtime_error("fiber::load_fiber_file: файл '" + filename + "' не содержит ни одной нити (FIBER)");
    }

    apply_inelastic(pending, materials, material_defined);

    build_model_from_polylines(fibers, materials, boxes, model);

    // Связи-стежки BOND: добавляются после цепочечных связей, prev/next не меняют.
    // l0 -- расстояние между бусинами в начальной конфигурации; l_break по eps_break
    // материала стежка. Масса бусин не меняется (масса швейной нити не учитывается).
    const std::size_t total_beads = model.fiber.size();
    for (const StitchSpec& st : stitches)
    {
        const std::string at = " в строке " + std::to_string(st.line);
        if (st.i < 0 || st.j < 0 || static_cast<std::size_t>(st.i) >= total_beads ||
            static_cast<std::size_t>(st.j) >= total_beads)
        {
            throw std::runtime_error("fiber::load_fiber_file: индекс бусины BOND вне диапазона [0, " +
                                     std::to_string(total_beads) + ")" + at);
        }
        if (st.i == st.j)
        {
            throw std::runtime_error("fiber::load_fiber_file: BOND соединяет бусину саму с собой" + at);
        }
        if (st.mat < 0 || static_cast<std::size_t>(st.mat) >= materials.size() ||
            !material_defined[static_cast<std::size_t>(st.mat)])
        {
            throw std::runtime_error("fiber::load_fiber_file: BOND ссылается на неопределённый материал " +
                                     std::to_string(st.mat) + at);
        }
        const std::size_t bi = static_cast<std::size_t>(st.i);
        const std::size_t bj = static_cast<std::size_t>(st.j);
        double d2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double dc = model.fiber.x[d][bj] - model.fiber.x[d][bi];
            d2 += dc * dc;
        }
        const double l0 = std::sqrt(d2);
        if (!(l0 > 0.0))
        {
            throw std::runtime_error("fiber::load_fiber_file: BOND между совпадающими бусинами (l0 = 0)" + at);
        }
        FiberBond b;
        b.i = bi;
        b.j = bj;
        b.l0 = l0;
        b.l_break = (1.0 + materials[static_cast<std::size_t>(st.mat)].eps_break) * l0;
        b.mat = static_cast<int>(st.mat);
        model.topology.bonds.push_back(b);
    }
}

} // namespace fiber
