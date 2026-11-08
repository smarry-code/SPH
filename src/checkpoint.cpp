// checkpoint.cpp
// Бинарный чекпоинт для побитового рестарта (include/checkpoint.hpp).
#include "checkpoint.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <type_traits>

#ifdef _WIN32
#include <io.h> // _commit, _fileno
#else
#include <unistd.h> // fsync, fileno
#endif

namespace checkpoint
{

namespace
{

static_assert(sizeof(double) == 8, "чекпоинт рассчитан на 8-байтовый double");
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "чекпоинт хранится в little-endian");
#endif

constexpr char kMagic[8] = {'S', 'P', 'H', 'C', 'K', 'P', 'T', '\0'};
constexpr std::uint32_t kFormatVersion = 3; // v2: в конец double-полей добавлены T, edot (Джонсон-Кук); v3: r0 (вес энергии в осесимметрии)

// Флаги сборки, влияющие на состав сохраняемых данных.
constexpr std::uint32_t build_flags()
{
    std::uint32_t f = 0;
#ifdef SPH_ENERGY_LOG
    f |= 1u;
#endif
#ifdef SPH_CONTACT_STATS
    f |= 2u;
#endif
    return f;
}

constexpr std::uint32_t make_tag(char a, char b, char c, char d)
{
    return static_cast<std::uint32_t>(static_cast<unsigned char>(a)) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(b)) << 8) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(c)) << 16) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(d)) << 24);
}

constexpr std::uint32_t kTagConf = make_tag('C', 'O', 'N', 'F');
constexpr std::uint32_t kTagPart = make_tag('P', 'A', 'R', 'T');
constexpr std::uint32_t kTagFibr = make_tag('F', 'I', 'B', 'R');
constexpr std::uint32_t kTagTopo = make_tag('T', 'O', 'P', 'O');
constexpr std::uint32_t kTagMats = make_tag('M', 'A', 'T', 'S');
constexpr std::uint32_t kTagClmp = make_tag('C', 'L', 'M', 'P');
constexpr std::uint32_t kTagLoop = make_tag('L', 'O', 'O', 'P');
constexpr std::uint32_t kTagHist = make_tag('H', 'I', 'S', 'T');
constexpr std::uint32_t kTagEnrg = make_tag('E', 'N', 'R', 'G');

[[noreturn]] void fail(const std::string& msg) { throw std::runtime_error("чекпоинт: " + msg); }

// --- CRC32 (IEEE 802.3, полином 0xEDB88320) ---------------------------------
const std::uint32_t* crc_table()
{
    static const auto table = []
    {
        struct T { std::uint32_t v[256]; } t{};
        for (std::uint32_t i = 0; i < 256; ++i)
        {
            std::uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            t.v[i] = c;
        }
        return t;
    }();
    return table.v;
}

std::uint32_t crc32(const unsigned char* p, size_t n)
{
    const std::uint32_t* t = crc_table();
    std::uint32_t c = 0xFFFFFFFFu;
    for (size_t i = 0; i < n; ++i) c = t[(c ^ p[i]) & 0xFFu] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

// --- Writer: накопление в памяти --------------------------------------------
class Writer
{
public:
    void raw(const void* p, size_t n)
    {
        const auto* b = static_cast<const unsigned char*>(p);
        buf_.insert(buf_.end(), b, b + n);
    }
    void u8(std::uint8_t v) { raw(&v, 1); }
    void u32(std::uint32_t v) { raw(&v, 4); }
    void u64(std::uint64_t v) { raw(&v, 8); }
    void i32(std::int32_t v) { raw(&v, 4); }
    void i64(std::int64_t v) { raw(&v, 8); }
    void f64(double v) { raw(&v, 8); }
    void str(const std::string& s)
    {
        u64(s.size());
        raw(s.data(), s.size());
    }
    // Вектор с префиксом-длиной.
    void f64vec(const std::vector<double>& v)
    {
        u64(v.size());
        raw(v.data(), v.size() * 8);
    }
    // Первые n элементов без префикса (n известно из контекста).
    void f64range(const std::vector<double>& v, size_t n) { raw(v.data(), n * 8); }

    size_t begin_block(std::uint32_t tag)
    {
        u32(tag);
        const size_t pos = buf_.size();
        u64(0);
        return pos;
    }
    void end_block(size_t pos)
    {
        const std::uint64_t len = buf_.size() - pos - 8;
        std::memcpy(buf_.data() + pos, &len, 8);
    }

    std::vector<unsigned char>& buf() { return buf_; }

private:
    std::vector<unsigned char> buf_;
};

// --- Reader: чтение с контролем границ --------------------------------------
class Reader
{
public:
    Reader(const unsigned char* p, size_t n) : p_(p), e_(p + n) {}

    size_t remaining() const { return static_cast<size_t>(e_ - p_); }
    void raw(void* dst, size_t n)
    {
        if (n > remaining()) fail("файл обрезан или повреждён (неожиданный конец данных)");
        std::memcpy(dst, p_, n);
        p_ += n;
    }
    std::uint8_t u8() { std::uint8_t v; raw(&v, 1); return v; }
    std::uint32_t u32() { std::uint32_t v; raw(&v, 4); return v; }
    std::uint64_t u64() { std::uint64_t v; raw(&v, 8); return v; }
    std::int32_t i32() { std::int32_t v; raw(&v, 4); return v; }
    std::int64_t i64() { std::int64_t v; raw(&v, 8); return v; }
    double f64() { double v; raw(&v, 8); return v; }
    bool boolean() { return u8() != 0; }
    std::string str()
    {
        const std::uint64_t n = u64();
        if (n > remaining()) fail("файл обрезан или повреждён (длина строки)");
        std::string s(reinterpret_cast<const char*>(p_), static_cast<size_t>(n));
        p_ += n;
        return s;
    }
    // Длина-префикс элементов по elem байт: проверка до выделения памяти.
    size_t count(size_t elem)
    {
        const std::uint64_t n = u64();
        if (elem != 0 && n > remaining() / elem) fail("файл обрезан или повреждён (длина массива)");
        return static_cast<size_t>(n);
    }
    std::vector<double> f64vec()
    {
        std::vector<double> v(count(8));
        raw(v.data(), v.size() * 8);
        return v;
    }
    std::vector<double> f64range(size_t n)
    {
        if (n > remaining() / 8) fail("файл обрезан или повреждён (длина массива)");
        std::vector<double> v(n);
        raw(v.data(), n * 8);
        return v;
    }

    // Вход в блок: проверка тега, подчитыватель на содержимое блока.
    Reader block(std::uint32_t tag, const char* name)
    {
        const std::uint32_t t = u32();
        const std::uint64_t len = u64();
        if (t != tag) fail(std::string("ожидался блок ") + name + ", в файле другой тег");
        if (len > remaining()) fail(std::string("блок ") + name + " обрезан");
        Reader sub(p_, static_cast<size_t>(len));
        p_ += len;
        return sub;
    }
    void expect_end(const char* name) const
    {
        if (p_ != e_) fail(std::string("блок ") + name + ": лишние данные (несовпадение формата)");
    }

private:
    const unsigned char* p_;
    const unsigned char* e_;
};

// --- Типизированные put/get/text для полей SimulationConfig -----------------
void put(Writer& w, int v) { w.i32(v); }
void put(Writer& w, long long v) { w.i64(v); }
void put(Writer& w, size_t v) { w.u64(v); }
void put(Writer& w, double v) { w.f64(v); }
void put(Writer& w, bool v) { w.u8(v ? 1 : 0); }
void put(Writer& w, const std::string& v) { w.str(v); }
void put(Writer& w, const std::vector<double>& v) { w.f64vec(v); }
void put(Writer& w, const std::vector<std::vector<double>>& v)
{
    w.u64(v.size());
    for (const auto& r : v) w.f64vec(r);
}

void get(Reader& r, int& v) { v = r.i32(); }
void get(Reader& r, long long& v) { v = r.i64(); }
void get(Reader& r, size_t& v) { v = static_cast<size_t>(r.u64()); }
void get(Reader& r, double& v) { v = r.f64(); }
void get(Reader& r, bool& v) { v = r.boolean(); }
void get(Reader& r, std::string& v) { v = r.str(); }
void get(Reader& r, std::vector<double>& v) { v = r.f64vec(); }
void get(Reader& r, std::vector<std::vector<double>>& v)
{
    v.resize(r.count(8));
    for (auto& row : v) row = r.f64vec();
}

std::string dtext(double v)
{
    char b[40];
    std::snprintf(b, sizeof b, "%.17g", v);
    return b;
}
std::string text(int v) { return std::to_string(v); }
std::string text(long long v) { return std::to_string(v); }
std::string text(size_t v) { return std::to_string(v); }
std::string text(double v) { return dtext(v); }
std::string text(bool v) { return v ? "true" : "false"; }
std::string text(const std::string& v) { return "\"" + v + "\""; }
std::string text(const std::vector<double>& v)
{
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + dtext(v[i]);
    return s + "]";
}
std::string text(const std::vector<std::vector<double>>& v)
{
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ";" : "") + text(v[i]);
    return s + "]";
}

// Перечень полей SimulationConfig: f(имя, поле, свободно_ли_менять). Порядок
// фиксирован версией формата (меняя его или состав полей, поднять kFormatVersion).
// Поля checkpoint_* / restart_file -- настройки запуска, не сериализуются.
template <class C, class F>
void visit_cfg(C& c, F&& f)
{
    f("pa_sph", c.pa_sph, false);
    f("nnps", c.nnps, false);
    f("sle", c.sle, false);
    f("skf", c.skf, false);
    f("density_method", c.density_method, false);
    f("coloring_niac_threshold", c.coloring_niac_threshold, false);
    f("ode_solver", c.ode_solver, false);
    f("contact_force_type", c.contact_force_type, false);
    f("contact_penalty_scale", c.contact_penalty_scale, false);
    f("contact_impulse_gamma", c.contact_impulse_gamma, false);
    f("fiber_solid_contact_type", c.fiber_solid_contact_type, false);
    f("fiber_fiber_contact_type", c.fiber_fiber_contact_type, false);
    f("fiber_solid_contact_penalty_scale", c.fiber_solid_contact_penalty_scale, false);
    f("fiber_solid_contact_impulse_gamma", c.fiber_solid_contact_impulse_gamma, false);
    f("fiber_fiber_contact_penalty_scale", c.fiber_fiber_contact_penalty_scale, false);
    f("fiber_fiber_contact_impulse_gamma", c.fiber_fiber_contact_impulse_gamma, false);
    f("fiber_solid_contact_impulse_push_tau", c.fiber_solid_contact_impulse_push_tau, false);
    f("fiber_fiber_contact_impulse_push_tau", c.fiber_fiber_contact_impulse_push_tau, false);
    f("fiber_contact_impulse_push_depth", c.fiber_contact_impulse_push_depth, false);
    f("fiber_contact_impulse_push_vmax", c.fiber_contact_impulse_push_vmax, false);
    f("dt_scale", c.dt_scale, false);
    f("sph_same_body_only", c.sph_same_body_only, false);
    f("CFL", c.CFL, false);
    f("t_total", c.t_total, true);
    f("save_t", c.save_t, true);
    f("print_step", c.print_step, true);
    f("alpha", c.alpha, false);
    f("beta", c.beta, false);
    f("etq", c.etq, false);
    f("epsilon", c.epsilon, false);
    f("xi", c.xi, false);
    f("epsilon_h", c.epsilon_h, false);
    f("num_threads", c.num_threads, true);
    f("maxn", c.maxn, false);
    f("max_interaction", c.max_interaction, false);
    f("moni_particle", c.moni_particle, false);
    f("mats.G", c.mats.G, false);
    f("mats.sigy", c.mats.sigy, false);
    f("mats.aa", c.mats.aa, false);
    f("mats.zz", c.mats.zz, false);
    f("mats.gamm", c.mats.gamm, false);
    f("mats.const_JH_2", c.mats.const_JH_2, false);
    f("fiber_file", c.fiber_file, false);
    f("k_file", c.k_file, false);
    f("k_mode", c.k_mode, false);
    f("k_parts_file", c.k_parts_file, false);
    f("k_length_scale", c.k_length_scale, false);
    f("k_mass_scale", c.k_mass_scale, false);
    f("k_monitor_nid", c.k_monitor_nid, false);
    f("k_solid_elements", c.k_solid_elements, false);
    f("jc_heating", c.jc_heating, false);
    f("jc_chi", c.jc_chi, false);
    f("xsph_mode", c.xsph_mode, false);
#ifdef SPH_ENERGY_LOG
    f("energy_log_every", c.energy_log_every, true);
    f("energy_log_file", c.energy_log_file, true);
#endif
#ifdef SPH_CONTACT_STATS
    f("contact_stats_every", c.contact_stats_every, true);
    f("contact_stats_file", c.contact_stats_file, true);
    f("contact_stats_summary", c.contact_stats_summary, true);
    f("contact_hold_n", c.contact_hold_n, true);
#endif
}

struct CfgField
{
    std::string name;
    bool free_to_change;
    std::vector<unsigned char> bytes; // побитовое представление для сравнения
    std::string txt;
};

std::vector<CfgField> collect_cfg(const SimulationConfig& c)
{
    std::vector<CfgField> out;
    visit_cfg(c, [&out](const char* name, const auto& v, bool fr)
              {
                  Writer w;
                  put(w, v);
                  out.push_back({name, fr, std::move(w.buf()), text(v)});
              });
    return out;
}

void write_cfg(Writer& w, const SimulationConfig& c)
{
    std::uint32_t cnt = 0;
    visit_cfg(c, [&cnt](const char*, const auto&, bool) { ++cnt; });
    w.u32(cnt);
    visit_cfg(c, [&w](const char*, const auto& v, bool) { put(w, v); });
}

SimulationConfig read_cfg(Reader& r)
{
    SimulationConfig c;
    std::uint32_t cnt = 0;
    visit_cfg(c, [&cnt](const char*, auto&, bool) { ++cnt; });
    if (r.u32() != cnt) fail("число полей SimulationConfig не совпадает (другой состав настроек)");
    visit_cfg(c, [&r](const char*, auto& v, bool) { get(r, v); });
    return c;
}

// --- Поля ParticleSystem -----------------------------------------------------
// Все double-массивы по порядку: x[d], vx[d], далее скаляры и тензоры.
template <class S>
auto double_fields(S& s)
{
    using Ptr = std::conditional_t<std::is_const_v<S>, const std::vector<double>*, std::vector<double>*>;
    std::vector<Ptr> v;
    for (size_t d = 0; d < dim; ++d) v.push_back(&s.x[d]);
    for (size_t d = 0; d < dim; ++d) v.push_back(&s.vx[d]);
    for (auto* p : {&s.mass, &s.rho, &s.rho_0, &s.p, &s.u, &s.c, &s.s, &s.e, &s.hsml, &s.dp, &s.e_pl, &s.D,
                    &s.sxx, &s.syy, &s.szz, &s.sxy, &s.syz, &s.sxz, &s.exx, &s.eyy, &s.ezz, &s.exy, &s.exz,
                    &s.eyz, &s.rxy, &s.rxz, &s.ryz, &s.T, &s.edot, &s.r0})
        v.push_back(p);
    return v;
}

// --- Запись/чтение блоков ----------------------------------------------------
void write_solid(Writer& w, const ParticleSystem& s, size_t n)
{
    const size_t b = w.begin_block(kTagPart);
    w.u64(n);
    for (const auto* f : double_fields(s)) w.f64range(*f, n);
    for (size_t i = 0; i < n; ++i) w.i32(s.itype[i]);
    for (size_t i = 0; i < n; ++i) w.u8(static_cast<std::uint8_t>(s.axis[i]));
    for (size_t i = 0; i < n; ++i) w.i32(s.body[i]);
    w.end_block(b);
}

void write_fiber(Writer& w, const fiber::FiberState& f)
{
    const size_t b = w.begin_block(kTagFibr);
    const size_t n = f.size();
    w.u64(n);
    for (size_t d = 0; d < dim; ++d) w.f64range(f.x[d], n);
    for (size_t d = 0; d < dim; ++d) w.f64range(f.vx[d], n);
    w.f64range(f.mass, n);
    for (size_t i = 0; i < n; ++i) w.i32(f.object_id[i]);
    for (size_t i = 0; i < n; ++i) w.i32(f.material_id[i]);
    w.end_block(b);
}

void write_topology(Writer& w, const fiber::FiberTopology& t)
{
    const size_t b = w.begin_block(kTagTopo);
    w.u64(t.prev.size());
    for (long long v : t.prev) w.i64(v);
    w.u64(t.next.size());
    for (long long v : t.next) w.i64(v);
    w.u64(t.version);
    w.u64(t.bonds.size());
    for (const auto& bd : t.bonds)
    {
        w.u64(bd.i);
        w.u64(bd.j);
        w.f64(bd.l0);
        w.f64(bd.l_break);
        const auto& s = bd.st;
        w.f64(s.eps_p);
        w.f64(s.q);
        w.f64(s.eps_prev);
        w.f64(s.w_p);
        w.f64(s.w_v);
        w.f64(s.d);
        w.f64(s.kappa);
        w.f64(s.eps0);
        w.f64(s.w_d);
    }
    w.end_block(b);
}

void write_materials(Writer& w, const std::vector<fiber::FiberMaterial>& ms)
{
    const size_t b = w.begin_block(kTagMats);
    w.u64(ms.size());
    for (const auto& m : ms)
    {
        w.f64(m.E);
        w.f64(m.nu);
        w.f64(m.d);
        w.f64(m.rho);
        w.f64(m.eps_break);
        w.f64(m.mu_ff);
        w.f64(m.mu_fs);
        w.i64(static_cast<std::int64_t>(m.plastic.type));
        w.f64(m.plastic.sy);
        w.f64(m.plastic.H);
        w.f64(m.plastic.K);
        w.f64(m.plastic.n);
        w.f64vec(m.plastic.table_eps);
        w.f64vec(m.plastic.table_sigma);
        w.i64(static_cast<std::int64_t>(m.rate.model));
        w.f64(m.rate.a);
        w.f64(m.rate.b);
        w.f64(m.visco.eta);
        w.f64(m.visco.xi);
        w.f64(m.visco.E1);
        w.f64(m.visco.tau);
        w.f64(m.failure.eps_p);
        w.f64(m.failure.sigma_u);
        w.f64(m.damage.deps_f);
        w.f64(m.damage.gf);
    }
    w.end_block(b);
}

void write_u8vec(Writer& w, const std::vector<std::uint8_t>& v)
{
    w.u64(v.size());
    w.raw(v.data(), v.size());
}

std::vector<std::uint8_t> read_u8vec(Reader& r)
{
    std::vector<std::uint8_t> v(r.count(1));
    r.raw(v.data(), v.size());
    return v;
}

void write_clamps(Writer& w, const interactions::ClampSet& c, const interactions::ClampReport& rep)
{
    const size_t b = w.begin_block(kTagClmp);
    write_u8vec(w, c.solid_mask);
    write_u8vec(w, c.fiber_mask);
    for (size_t d = 0; d < dim; ++d) w.f64(rep.solid_force[d]);
    for (size_t d = 0; d < dim; ++d) w.f64(rep.fiber_force[d]);
    w.u64(rep.solid_clamped);
    w.u64(rep.fiber_clamped);
    w.f64(rep.solid_max_h);
    w.end_block(b);
}

void write_loop(Writer& w, const LoopState& l)
{
    const size_t b = w.begin_block(kTagLoop);
    w.f64(l.time);
    w.f64(l.dt);
    w.f64(l.next_output_time);
    w.u64(l.step);
    w.u64(l.guard_steps);
    w.f64(l.guard_dt_first);
    w.f64(l.guard_dt_min);
    w.end_block(b);
}

void write_hist(Writer& w, const SolverHistory& h)
{
    const size_t b = w.begin_block(kTagHist);
    w.u8(h.first_step ? 1 : 0);
    for (const auto* v : {&h.du, &h.drho, &h.dsxx, &h.dsyy, &h.dszz, &h.dsxy, &h.dsyz, &h.dsxz}) w.f64vec(*v);
    for (size_t d = 0; d < dim; ++d) w.f64vec(h.dvx[d]);
    for (size_t d = 0; d < dim; ++d) w.f64vec(h.fiber_dv[d]);
    w.end_block(b);
}

void write_energy(Writer& w, const Payload& p)
{
    const size_t b = w.begin_block(kTagEnrg);
    w.u8(p.has_energy ? 1 : 0);
    for (double c : p.energy.chan) w.f64(c);
    w.f64(p.energy.e_mech0);
    w.u8(p.energy.have_e0 ? 1 : 0);
    w.end_block(b);
}

// Проверка длины вектора истории: 0 (нет истории) или n.
void check_hist_len(const std::vector<double>& v, size_t n, const char* name)
{
    if (!v.empty() && v.size() != n)
        fail(std::string("история решателя: длина ") + name + " (" + std::to_string(v.size()) +
             ") не равна " + std::to_string(n));
}

std::string temp_path(const std::string& path) { return path + ".tmp"; }

} // namespace

// ============================================================================
void write(const std::string& path, const SimulationConfig& cfg, const fiber::Model& model, size_t ntotal,
           const Payload& payload)
{
    if (model.solid.size() < ntotal) fail("запись: размер Solid меньше ntotal");
    const size_t nf = model.fiber.size();
    if (!model.topology.prev.empty() && model.topology.size() != nf)
        fail("запись: размер топологии не равен числу бусин нитей");
    {
        const auto& h = payload.hist;
        for (const auto* v : {&h.du, &h.drho, &h.dsxx, &h.dsyy, &h.dszz, &h.dsxy, &h.dsyz, &h.dsxz})
            check_hist_len(*v, ntotal, "solid-истории");
        for (size_t d = 0; d < dim; ++d)
        {
            check_hist_len(h.dvx[d], ntotal, "dvx");
            check_hist_len(h.fiber_dv[d], nf, "fiber_dv");
        }
    }

    Writer w;
    w.raw(kMagic, 8);
    w.u32(kFormatVersion);
    w.u32(static_cast<std::uint32_t>(dim));
    w.u32(build_flags());
    w.u32(static_cast<std::uint32_t>(sizeof(double)));

    {
        const size_t b = w.begin_block(kTagConf);
        write_cfg(w, cfg);
        w.end_block(b);
    }
    write_solid(w, model.solid, ntotal);
    write_fiber(w, model.fiber);
    write_topology(w, model.topology);
    write_materials(w, model.materials);
    write_clamps(w, model.clamps, model.clamp_report);
    write_loop(w, payload.loop);
    write_hist(w, payload.hist);
    write_energy(w, payload);
    w.u32(crc32(w.buf().data(), w.buf().size()));

    const std::string tmp = temp_path(path);
    std::FILE* f = std::fopen(tmp.c_str(), "wb");
    if (!f) fail("не удалось создать временный файл " + tmp);
    const auto& buf = w.buf();
    bool ok = std::fwrite(buf.data(), 1, buf.size(), f) == buf.size();
    ok = ok && std::fflush(f) == 0;
#ifdef _WIN32
    ok = ok && _commit(_fileno(f)) == 0;
#else
    ok = ok && fsync(fileno(f)) == 0;
#endif
    ok = (std::fclose(f) == 0) && ok;
    if (!ok)
    {
        std::remove(tmp.c_str());
        fail("ошибка записи " + tmp + " (диск заполнен?)");
    }

    // Ротация: предыдущий чекпоинт -> path.prev, затем атомарная подмена.
    if (std::FILE* ex = std::fopen(path.c_str(), "rb"))
    {
        std::fclose(ex);
#ifdef _WIN32
        // На Windows rename не перезаписывает существующий файл — убираем старый .prev.
        // На POSIX rename заменяет .prev атомарно, удаление не нужно.
        std::remove((path + ".prev").c_str());
#endif
        if (std::rename(path.c_str(), (path + ".prev").c_str()) != 0)
        {
            std::remove(tmp.c_str());
            fail("не удалось переименовать " + path + " в .prev");
        }
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0)
    {
        std::remove(tmp.c_str());
        fail("не удалось переименовать " + tmp + " в " + path);
    }
}

// ============================================================================
void read(const std::string& path, SimulationConfig& cfg, fiber::Model& model, size_t& ntotal, Payload& payload,
          const RestoreOptions& opt, ConfigDiff* diff)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) fail("не удалось открыть файл " + path);
    std::vector<unsigned char> data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    constexpr size_t kHeader = 8 + 4 * 4;
    if (data.size() < kHeader + 4) fail("файл слишком короткий (обрезан): " + path);
    if (std::memcmp(data.data(), kMagic, 8) != 0) fail("неверная сигнатура файла (не чекпоинт SPH): " + path);

    Reader hr(data.data() + 8, kHeader - 8);
    const std::uint32_t ver = hr.u32(), fdim = hr.u32(), flags = hr.u32(), dsz = hr.u32();
    if (ver != kFormatVersion)
        fail("версия формата " + std::to_string(ver) + " не поддерживается (ожидается " +
             std::to_string(kFormatVersion) + ")");
    if (fdim != dim)
        fail("SPH_DIM чекпоинта (" + std::to_string(fdim) + ") не равен SPH_DIM сборки (" + std::to_string(dim) + ")");
    if (flags != build_flags())
        fail("флаги сборки чекпоинта (" + std::to_string(flags) + ") не совпадают с текущей сборкой (" +
             std::to_string(build_flags()) + "; бит 1 = SPH_ENERGY_LOG, бит 2 = SPH_CONTACT_STATS)");
    if (dsz != sizeof(double)) fail("размер double в чекпоинте не совпадает с текущей сборкой");

    const size_t body = data.size() - 4;
    std::uint32_t stored_crc;
    std::memcpy(&stored_crc, data.data() + body, 4);
    if (crc32(data.data(), body) != stored_crc) fail("CRC не совпадает (файл повреждён или обрезан): " + path);

    Reader top(data.data() + kHeader, body - kHeader);

    // --- разбор в промежуточные объекты; модель меняется только после проверок ---
    Reader rc = top.block(kTagConf, "CONF");
    const SimulationConfig saved = read_cfg(rc);
    rc.expect_end("CONF");

    Reader rp = top.block(kTagPart, "PART");
    const size_t n_saved = static_cast<size_t>(rp.u64());
    if (n_saved != ntotal)
        fail("число частиц Solid в чекпоинте (" + std::to_string(n_saved) + ") не равно ntotal сцены (" +
             std::to_string(ntotal) + ")");
    if (model.solid.size() < ntotal) fail("размер Solid сцены меньше ntotal");
    ParticleSystem ps;
    ps.resize(ntotal);
    for (auto* f : double_fields(ps)) *f = rp.f64range(ntotal);
    for (size_t i = 0; i < ntotal; ++i) ps.itype[i] = rp.i32();
    for (size_t i = 0; i < ntotal; ++i) ps.axis[i] = static_cast<char>(rp.u8());
    for (size_t i = 0; i < ntotal; ++i) ps.body[i] = rp.i32();
    rp.expect_end("PART");

    Reader rf = top.block(kTagFibr, "FIBR");
    const size_t nf = static_cast<size_t>(rf.u64());
    if (nf != model.fiber.size())
        fail("число бусин нитей в чекпоинте (" + std::to_string(nf) + ") не равно сцене (" +
             std::to_string(model.fiber.size()) + ")");
    fiber::FiberState fs;
    fs.resize(nf);
    for (size_t d = 0; d < dim; ++d) fs.x[d] = rf.f64range(nf);
    for (size_t d = 0; d < dim; ++d) fs.vx[d] = rf.f64range(nf);
    fs.mass = rf.f64range(nf);
    for (size_t i = 0; i < nf; ++i) fs.object_id[i] = rf.i32();
    for (size_t i = 0; i < nf; ++i) fs.material_id[i] = rf.i32();
    rf.expect_end("FIBR");

    Reader rt = top.block(kTagTopo, "TOPO");
    fiber::FiberTopology topo;
    topo.prev.resize(rt.count(8));
    for (auto& v : topo.prev) v = rt.i64();
    topo.next.resize(rt.count(8));
    for (auto& v : topo.next) v = rt.i64();
    topo.version = rt.u64();
    topo.bonds.resize(rt.count(8 * 13));
    for (auto& bd : topo.bonds)
    {
        bd.i = static_cast<size_t>(rt.u64());
        bd.j = static_cast<size_t>(rt.u64());
        bd.l0 = rt.f64();
        bd.l_break = rt.f64();
        auto& s = bd.st;
        s.eps_p = rt.f64();
        s.q = rt.f64();
        s.eps_prev = rt.f64();
        s.w_p = rt.f64();
        s.w_v = rt.f64();
        s.d = rt.f64();
        s.kappa = rt.f64();
        s.eps0 = rt.f64();
        s.w_d = rt.f64();
        if (bd.i >= nf || bd.j >= nf) fail("индекс связи нити выходит за число бусин");
    }
    rt.expect_end("TOPO");
    if (topo.prev.size() != model.topology.size() || topo.next.size() != model.topology.size())
        fail("размер топологии нитей не совпадает со сценой");

    Reader rm = top.block(kTagMats, "MATS");
    std::vector<fiber::FiberMaterial> mats(rm.count(8 * 20));
    for (auto& m : mats)
    {
        m.E = rm.f64();
        m.nu = rm.f64();
        m.d = rm.f64();
        m.rho = rm.f64();
        m.eps_break = rm.f64();
        m.mu_ff = rm.f64();
        m.mu_fs = rm.f64();
        m.plastic.type = static_cast<fiber::HardeningType>(rm.i64());
        m.plastic.sy = rm.f64();
        m.plastic.H = rm.f64();
        m.plastic.K = rm.f64();
        m.plastic.n = rm.f64();
        m.plastic.table_eps = rm.f64vec();
        m.plastic.table_sigma = rm.f64vec();
        m.rate.model = static_cast<fiber::RateModel>(rm.i64());
        m.rate.a = rm.f64();
        m.rate.b = rm.f64();
        m.visco.eta = rm.f64();
        m.visco.xi = rm.f64();
        m.visco.E1 = rm.f64();
        m.visco.tau = rm.f64();
        m.failure.eps_p = rm.f64();
        m.failure.sigma_u = rm.f64();
        m.damage.deps_f = rm.f64();
        m.damage.gf = rm.f64();
    }
    rm.expect_end("MATS");

    Reader rk = top.block(kTagClmp, "CLMP");
    interactions::ClampSet clamps;
    interactions::ClampReport rep;
    clamps.solid_mask = read_u8vec(rk);
    clamps.fiber_mask = read_u8vec(rk);
    for (size_t d = 0; d < dim; ++d) rep.solid_force[d] = rk.f64();
    for (size_t d = 0; d < dim; ++d) rep.fiber_force[d] = rk.f64();
    rep.solid_clamped = static_cast<size_t>(rk.u64());
    rep.fiber_clamped = static_cast<size_t>(rk.u64());
    rep.solid_max_h = rk.f64();
    rk.expect_end("CLMP");
    if (!clamps.solid_mask.empty() && clamps.solid_mask.size() != ntotal)
        fail("размер solid_mask зажимов не равен ntotal");
    if (!clamps.fiber_mask.empty() && clamps.fiber_mask.size() != nf)
        fail("размер fiber_mask зажимов не равен числу бусин");

    Reader rl = top.block(kTagLoop, "LOOP");
    Payload pl;
    pl.loop.time = rl.f64();
    pl.loop.dt = rl.f64();
    pl.loop.next_output_time = rl.f64();
    pl.loop.step = static_cast<size_t>(rl.u64());
    pl.loop.guard_steps = static_cast<size_t>(rl.u64());
    pl.loop.guard_dt_first = rl.f64();
    pl.loop.guard_dt_min = rl.f64();
    rl.expect_end("LOOP");

    Reader rh = top.block(kTagHist, "HIST");
    auto& h = pl.hist;
    h.first_step = rh.boolean();
    for (auto* v : {&h.du, &h.drho, &h.dsxx, &h.dsyy, &h.dszz, &h.dsxy, &h.dsyz, &h.dsxz})
    {
        *v = rh.f64vec();
        check_hist_len(*v, ntotal, "solid-истории");
    }
    for (size_t d = 0; d < dim; ++d)
    {
        h.dvx[d] = rh.f64vec();
        check_hist_len(h.dvx[d], ntotal, "dvx");
    }
    for (size_t d = 0; d < dim; ++d)
    {
        h.fiber_dv[d] = rh.f64vec();
        check_hist_len(h.fiber_dv[d], nf, "fiber_dv");
    }
    rh.expect_end("HIST");

    Reader re = top.block(kTagEnrg, "ENRG");
    pl.has_energy = re.boolean();
    for (double& c : pl.energy.chan) c = re.f64();
    pl.energy.e_mech0 = re.f64();
    pl.energy.have_e0 = re.boolean();
    re.expect_end("ENRG");
    top.expect_end("файл");

    // --- сравнение настроек ---
    ConfigDiff d;
    {
        const auto a = collect_cfg(saved);
        const auto b = collect_cfg(cfg);
        for (size_t k = 0; k < a.size(); ++k)
            if (a[k].bytes != b[k].bytes)
            {
                const std::string line = a[k].name + ": " + a[k].txt + " -> " + b[k].txt;
                (a[k].free_to_change ? d.allowed : d.forbidden).push_back(line);
            }
    }
    if (diff) *diff = d;
    if (!d.forbidden.empty())
    {
        std::string list;
        for (const auto& s : d.forbidden) list += "\n  " + s;
        if (!opt.force)
            fail("настройки запуска отличаются от сохранённых (допускается только с --force):" + list);
        std::cerr << "ПРЕДУПРЕЖДЕНИЕ: рестарт с --force, изменены настройки:" << list << "\n"
                  << force_risks_text() << std::endl;
    }

    // --- фиксация: проверки пройдены, перезаписываем динамическое состояние ---
    ParticleSystem& s = model.solid;
    {
        auto dst = double_fields(s);
        auto src = double_fields(static_cast<const ParticleSystem&>(ps));
        for (size_t k = 0; k < dst.size(); ++k) std::copy_n(src[k]->begin(), ntotal, dst[k]->begin());
    }
    std::copy_n(ps.itype.begin(), ntotal, s.itype.begin());
    std::copy_n(ps.axis.begin(), ntotal, s.axis.begin());
    std::copy_n(ps.body.begin(), ntotal, s.body.begin());
    model.fiber = std::move(fs);
    model.topology = std::move(topo);
    model.materials = std::move(mats);
    model.clamps = std::move(clamps);
    model.clamp_report = rep;
    payload = std::move(pl);
    // cfg остаётся новым: свободные и (при force) изменённые поля берутся из него
}

std::string force_risks_text()
{
    return "Риски --force (рестарт с изменёнными настройками):\n"
           "  1. Траектория перестаёт быть продолжением исходного расчёта; результат нельзя сравнивать с непрерывным прогоном.\n"
           "  2. Смена CFL/dt_scale/ode_solver ломает историю предиктора (f^{n-1} посчитана другим шагом/схемой): возможен скачок и потеря устойчивости.\n"
           "  3. Смена контактов/материалов даёт разрыв по энергии: энергобаланс EnergyLogger теряет смысл, база E_mech0 взята от старых настроек.\n"
           "  4. Смена sle/skf/nnps/alpha/beta: hsml, напряжения и искусственная вязкость несогласованы с новыми параметрами.\n"
           "  5. Смена num_threads нарушает побитовую воспроизводимость (порядок редукции OpenMP); без --force это допускается, предупреждение остаётся.\n"
           "  6. Серия с разными настройками в одном выстреле нарушает вердикт Guo: окно пробития подбиралось при других настройках.";
}

} // namespace checkpoint
