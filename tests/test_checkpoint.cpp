// tests/test_checkpoint.cpp
// Регрессия чекпоинта/рестарта (ветка fiber_checkpoint_restart).
//
// Проверяется, что «N шагов -> чекпоинт -> новая сцена + рестарт -> ещё N шагов»
// даёт ПОБИТОВО то же состояние, что непрерывный прогон 2N шагов (OMP_NUM_THREADS=1),
// а также целостность файла, политика настроек, периодичность и ротация записи.
// Сцена -- tests/checkpoint_scene.hpp. Собирается с SPH_DIM=3 (CMakeLists.txt);
// вторая цель test_checkpoint_energy отличается флагом SPH_ENERGY_LOG и добавляет
// проверку продолжения CSV-лога энергии.

#include "checkpoint.hpp"
#include "checkpoint_scene.hpp"
#include "e2e_harness.hpp"
#include "run_guard.hpp"
#include "time_integration.hpp"

#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace
{

namespace fs = std::filesystem;

constexpr std::size_t N = 20; // половина полного прогона, шагов

// Мини-харнесс: каждый кейс печатает PASS/FAIL, итог -- ненулевой код при провале.
int g_failed_cases = 0;

void run_case(const char* name, const std::function<void()>& body)
{
    std::printf("[ case ] %s\n", name);
    int before = e2e::failures;
    try
    {
        body();
    }
    catch (const std::exception& e)
    {
        e2e::check(false, std::string("исключение: ") + e.what());
    }
    const bool ok = e2e::failures == before;
    std::printf("%s: %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++g_failed_cases;
}

// Прогон до step_limit шагов (общих, с учётом восстановленных счётчиков сторожа).
// Возвращает остановку по StepLimit; time и dt -- по ссылке, как в time_integration.
bool run_steps(CkptScene& sc, std::size_t step_limit, double& time)
{
    RunGuard guard(0.0, 0.0, 0);
    guard.set_step_limit(step_limit);
    {
        e2e::ScopedSilence silence(true);
        time_integration(sc.cfg, sc.model, sc.ntotal, sc.dt, time, sc.delta_x, &guard);
    }
    return guard.reason() == RunGuard::Reason::StepLimit && guard.steps() == step_limit;
}

// ---------------------------------------------------------------------------
// Побитовый снимок состояния: именованные блоки байтов.
// ---------------------------------------------------------------------------
struct Blob
{
    std::string name;
    std::vector<unsigned char> bytes;
};

template <class T>
void add_vec(std::vector<Blob>& out, const std::string& name, const std::vector<T>& v, std::size_t count)
{
    Blob b;
    b.name = name;
    b.bytes.resize(count * sizeof(T));
    if (count > 0) std::memcpy(b.bytes.data(), v.data(), count * sizeof(T));
    out.push_back(std::move(b));
}

template <class T>
void add_val(std::vector<Blob>& out, const std::string& name, const T& v)
{
    Blob b;
    b.name = name;
    b.bytes.resize(sizeof(T));
    std::memcpy(b.bytes.data(), &v, sizeof(T));
    out.push_back(std::move(b));
}

std::vector<Blob> snapshot(const CkptScene& sc, double time)
{
    std::vector<Blob> o;
    const std::size_t n = sc.ntotal;
    const ParticleSystem& s = sc.model.solid;
    for (std::size_t d = 0; d < dim; ++d)
    {
        add_vec(o, "solid.x[" + std::to_string(d) + "]", s.x[d], n);
        add_vec(o, "solid.vx[" + std::to_string(d) + "]", s.vx[d], n);
    }
    add_vec(o, "solid.rho", s.rho, n);
    add_vec(o, "solid.u", s.u, n);
    add_vec(o, "solid.p", s.p, n);
    add_vec(o, "solid.hsml", s.hsml, n);
    add_vec(o, "solid.sxx", s.sxx, n);
    add_vec(o, "solid.syy", s.syy, n);
    add_vec(o, "solid.szz", s.szz, n);
    add_vec(o, "solid.sxy", s.sxy, n);
    add_vec(o, "solid.syz", s.syz, n);
    add_vec(o, "solid.sxz", s.sxz, n);
    add_vec(o, "solid.e_pl", s.e_pl, n);
    add_vec(o, "solid.D", s.D, n);
    add_vec(o, "solid.dp", s.dp, n);

    const fiber::FiberState& f = sc.model.fiber;
    for (std::size_t d = 0; d < dim; ++d)
    {
        add_vec(o, "fiber.x[" + std::to_string(d) + "]", f.x[d], f.size());
        add_vec(o, "fiber.vx[" + std::to_string(d) + "]", f.vx[d], f.size());
    }
    add_vec(o, "fiber.object_id", f.object_id, f.size());

    const fiber::FiberTopology& t = sc.model.topology;
    add_vec(o, "topology.prev", t.prev, t.prev.size());
    add_vec(o, "topology.next", t.next, t.next.size());
    add_val(o, "topology.version", t.version);
    add_val(o, "topology.bonds.size", t.bonds.size());
    for (std::size_t k = 0; k < t.bonds.size(); ++k)
    {
        const fiber::FiberBond& b = t.bonds[k];
        const std::string p = "bond[" + std::to_string(k) + "].";
        add_val(o, p + "i", b.i);
        add_val(o, p + "j", b.j);
        add_val(o, p + "l0", b.l0);
        add_val(o, p + "l_break", b.l_break);
        // Каждое поле BondState отдельно: отступы структуры в сравнение не попадают.
        add_val(o, p + "eps_p", b.st.eps_p);
        add_val(o, p + "q", b.st.q);
        add_val(o, p + "eps_prev", b.st.eps_prev);
        add_val(o, p + "w_p", b.st.w_p);
        add_val(o, p + "w_v", b.st.w_v);
        add_val(o, p + "d", b.st.d);
        add_val(o, p + "kappa", b.st.kappa);
        add_val(o, p + "eps0", b.st.eps0);
        add_val(o, p + "w_d", b.st.w_d);
    }
    add_val(o, "time", time);
    add_val(o, "dt", sc.dt);
    return o;
}

// Сравнивает снимки, печатает первое отличие (поле и индекс элемента).
bool compare_snapshots(const std::vector<Blob>& a, const std::vector<Blob>& b)
{
    if (a.size() != b.size())
    {
        std::printf("         число блоков: %zu против %zu\n", a.size(), b.size());
        return false;
    }
    int reported = 0;
    for (std::size_t k = 0; k < a.size(); ++k)
    {
        if (a[k].name != b[k].name || a[k].bytes != b[k].bytes)
        {
            if (reported++ < 5)
            {
                std::size_t first = 0;
                while (first < a[k].bytes.size() && first < b[k].bytes.size() && a[k].bytes[first] == b[k].bytes[first])
                    ++first;
                std::printf("         отличие: поле %s, байт %zu (элемент double %zu); размеры %zu / %zu\n",
                            a[k].name.c_str(), first, first / 8, a[k].bytes.size(), b[k].bytes.size());
            }
        }
    }
    return reported == 0;
}

fs::path scratch_root(const char* tag)
{
    return fs::temp_directory_path() / (std::string("sph_ckpt_") + tag);
}

// Общая часть кейсов 1 и 2: A -- 2N шагов подряд; B -- N шагов с чекпоинтом, затем рестарт.
// Возвращает число связей в A на шаге N (до рестарта оно же) для проверки разрыва.
void roundtrip(const char* tag, double eps_break, int ode_solver, bool expect_break)
{
    e2e::ScopedWorkDir wd("sph_ckpt", tag);
    const fs::path root = wd.dir();
    const std::string ck = (root / "state.ckpt").string();

    // Прогон A: непрерывный 2N. Заодно запомним число связей на шаге N отдельным прогоном A_N.
    std::size_t bonds_initial = 0, bonds_at_N = 0;
    {
        CkptScene sc = build_checkpoint_scene((root / "probe").string(), eps_break, ode_solver);
        bonds_initial = sc.model.topology.bonds.size();
        double time = 0.0;
        e2e::check(run_steps(sc, N, time), "зонд: N шагов выполнены");
        bonds_at_N = sc.model.topology.bonds.size();
    }
    std::ostringstream os;
    os << "связей: исходно " << bonds_initial << ", на шаге " << N << " -- " << bonds_at_N;
    e2e::note(os.str());
    if (expect_break)
        e2e::check(bonds_at_N < bonds_initial, "связи рвутся ДО чекпоинта (сценарий осмысленный)");
    else
        e2e::check(bonds_at_N == bonds_initial, "связи не рвутся (базовый сценарий)");

    std::vector<Blob> snap_a;
    {
        CkptScene sc = build_checkpoint_scene((root / "A").string(), eps_break, ode_solver);
        double time = 0.0;
        e2e::check(run_steps(sc, 2 * N, time), "A: 2N шагов выполнены");
        snap_a = snapshot(sc, time);
    }

    // Прогон B, часть 1: N шагов, чекпоинт на шаге N (периодический = финальный).
    {
        CkptScene sc = build_checkpoint_scene((root / "B1").string(), eps_break, ode_solver);
        sc.cfg.checkpoint_file = ck;
        sc.cfg.checkpoint_every_steps = N;
        double time = 0.0;
        e2e::check(run_steps(sc, N, time), "B1: N шагов выполнены");
    }
    e2e::check(fs::exists(ck), "чекпоинт записан");

    // Часть 2: новая сцена в отдельном каталоге, рестарт, лимит 2N (ПОЛНЫЙ).
    std::vector<Blob> snap_b;
    {
        CkptScene sc = build_checkpoint_scene((root / "B2").string(), eps_break, ode_solver);
        sc.cfg.restart_file = ck;
        double time = 0.0;
        e2e::check(run_steps(sc, 2 * N, time), "B2: рестарт дошёл до 2N шагов");
        snap_b = snapshot(sc, time);
    }
    e2e::check(compare_snapshots(snap_a, snap_b), "состояние после рестарта побитово равно непрерывному прогону");
}

// Файл целиком в память / из памяти.
std::vector<char> slurp(const std::string& p)
{
    std::ifstream in(p, std::ios::binary);
    return std::vector<char>((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}
void spit(const std::string& p, const std::vector<char>& v)
{
    std::ofstream out(p, std::ios::binary | std::ios::trunc);
    out.write(v.data(), static_cast<std::streamsize>(v.size()));
}

// Смещения начал блоков верхнего уровня (тег u32 + длина u64 + тело) после заголовка 24 байта; CRC -- последние 4.
std::vector<std::size_t> block_starts(const std::vector<char>& raw)
{
    std::vector<std::size_t> v;
    std::size_t pos = 24;
    while (pos + 12 + 4 <= raw.size())
    {
        std::uint64_t len;
        std::memcpy(&len, raw.data() + pos + 4, 8);
        v.push_back(pos);
        pos += 12 + static_cast<std::size_t>(len);
    }
    return v;
}

// CRC32 (IEEE, как в src/checkpoint.cpp) по всему буферу.
std::uint32_t crc32_of(const std::vector<char>& raw)
{
    std::uint32_t crc = 0xFFFFFFFFu;
    for (char ch : raw)
    {
        crc ^= static_cast<unsigned char>(ch);
        for (int b = 0; b < 8; ++b) crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
    }
    return ~crc;
}

// Пишет чекпоинт после N шагов в path; сцена остаётся в cwd workdir.
void make_checkpoint(const std::string& workdir, const std::string& path, std::size_t steps, std::size_t every)
{
    CkptScene sc = build_checkpoint_scene(workdir, 0.5, 0);
    sc.cfg.checkpoint_file = path;
    sc.cfg.checkpoint_every_steps = every;
    double time = 0.0;
    run_steps(sc, steps, time);
}

// Читает path новой сценой с заданной правкой настроек; возвращает текст ошибки ("" -- успех).
std::string try_read(const std::string& workdir, const std::string& path, bool force,
                     const std::function<void(SimulationConfig&)>& tweak, checkpoint::Payload* out = nullptr,
                     checkpoint::ConfigDiff* diff = nullptr)
{
    CkptScene sc = build_checkpoint_scene(workdir, 0.5, 0);
    tweak(sc.cfg);
    checkpoint::Payload pl;
    checkpoint::RestoreOptions opt;
    opt.force = force;
    std::size_t nt = sc.ntotal;
    try
    {
        checkpoint::read(path, sc.cfg, sc.model, nt, out ? *out : pl, opt, diff);
    }
    catch (const std::exception& e)
    {
        return e.what();
    }
    return "";
}

#ifdef SPH_ENERGY_LOG
std::vector<std::string> read_lines(const std::string& p)
{
    std::ifstream in(p);
    std::vector<std::string> v;
    for (std::string l; std::getline(in, l);) v.push_back(l);
    return v;
}
#endif

} // namespace

int main()
{
    // 1. Побитовый round-trip, предиктор-корректор и Эйлер
    run_case("1a. round-trip PC (ode_solver=0), связи целы", [] { roundtrip("rt_pc", 0.5, 0, false); });
    run_case("1b. round-trip Euler (ode_solver=1), связи целы", [] { roundtrip("rt_eu", 0.5, 1, false); });

    // 2. Разрыв связей до чекпоинта (eps_break=0.002 рвёт 4 из 6 связей за 40 шагов -- проверяется зондом)
    run_case("2a. round-trip PC, связи рвутся до шага N", [] { roundtrip("brk_pc", 0.002, 0, true); });
    run_case("2b. round-trip Euler, связи рвутся до шага N", [] { roundtrip("brk_eu", 0.002, 1, true); });

    // 3. Целостность
    run_case("3. целостность: порча байта, обрезка, неверный magic", [] {
        e2e::ScopedWorkDir wd("sph_ckpt", "integrity");
        const std::string root = wd.dir().string();
        const std::string ck = root + "/ok.ckpt";
        make_checkpoint(root + "/w0", ck, 5, 5);
        e2e::check(try_read(root + "/w1", ck, false, [](SimulationConfig&) {}).empty(), "исходный файл читается");

        std::vector<char> good = slurp(ck);
        e2e::check(good.size() > 100, "файл непустой");

        std::vector<char> flipped = good;
        flipped[flipped.size() / 2] ^= 0x01;
        spit(root + "/flip.ckpt", flipped);
        const std::string e1 = try_read(root + "/w2", root + "/flip.ckpt", false, [](SimulationConfig&) {});
        e2e::check(!e1.empty(), "порча байта -> ошибка");
        e2e::note("сообщение: " + e1);

        std::vector<char> trunc(good.begin(), good.begin() + static_cast<std::ptrdiff_t>(good.size() / 2));
        spit(root + "/trunc.ckpt", trunc);
        const std::string e2 = try_read(root + "/w3", root + "/trunc.ckpt", false, [](SimulationConfig&) {});
        e2e::check(!e2.empty(), "обрезанный файл -> ошибка");
        e2e::note("сообщение: " + e2);

        std::vector<char> magic = good;
        magic[0] = 'X';
        spit(root + "/magic.ckpt", magic);
        const std::string e3 = try_read(root + "/w4", root + "/magic.ckpt", false, [](SimulationConfig&) {});
        e2e::check(!e3.empty(), "неверный magic -> ошибка");
        e2e::note("сообщение: " + e3);

        const std::string e4 = try_read(root + "/w5", root + "/net_takogo.ckpt", false, [](SimulationConfig&) {});
        e2e::check(!e4.empty(), "несуществующий файл -> ошибка");
    });

    // 4. Политика настроек
    run_case("4. политика настроек: t_total/save_t -- да, CFL -- только с force", [] {
        e2e::ScopedWorkDir wd("sph_ckpt", "policy");
        const std::string root = wd.dir().string();
        const std::string ck = root + "/ok.ckpt";
        make_checkpoint(root + "/w0", ck, 5, 5);

        e2e::check(try_read(root + "/w1", ck, false,
                            [](SimulationConfig& c) {
                                c.t_total = 2.0;
                                c.save_t = 1.0e-7;
                            })
                       .empty(),
                   "смена t_total/save_t без force проходит");

        const std::string e = try_read(root + "/w2", ck, false, [](SimulationConfig& c) { c.CFL = 0.05; });
        e2e::check(!e.empty() && e.find("CFL") != std::string::npos, "смена CFL без force -> ошибка со словом CFL");
        e2e::note("сообщение: " + e);

        const std::string e2 = try_read(root + "/w3", ck, false, [](SimulationConfig& c) { c.dt_scale = 0.5; });
        e2e::check(!e2.empty() && e2.find("dt_scale") != std::string::npos, "смена dt_scale без force -> ошибка со словом dt_scale");

        checkpoint::ConfigDiff diff;
        double cfl_after = 0.0;
        {
            CkptScene sc = build_checkpoint_scene(root + "/w4", 0.5, 0);
            sc.cfg.CFL = 0.05;
            checkpoint::Payload pl;
            checkpoint::RestoreOptions opt;
            opt.force = true;
            std::size_t nt = sc.ntotal;
            checkpoint::read(ck, sc.cfg, sc.model, nt, pl, opt, &diff);
            cfl_after = sc.cfg.CFL;
        }
        e2e::check(cfl_after == 0.05, "с force действует новое значение CFL");
        e2e::check(!diff.forbidden.empty(), "с force ConfigDiff.forbidden непуст");
    });

    // 4a. Периодичность
    run_case("4a. периодичность checkpoint_every_steps=7", [] {
        e2e::ScopedWorkDir wd("sph_ckpt", "period");
        const std::string root = wd.dir().string();
        auto last_step = [&](const std::string& tag, std::size_t limit, std::size_t every) -> long long {
            const std::string ck = root + "/" + tag + ".ckpt";
            make_checkpoint(root + "/w_" + tag, ck, limit, every);
            if (!fs::exists(ck)) return -1;
            checkpoint::Payload pl;
            const std::string err = try_read(root + "/r_" + tag, ck, false, [](SimulationConfig&) {}, &pl);
            if (!err.empty()) return -2;
            return static_cast<long long>(pl.loop.step);
        };
        // Лимит 15 при периоде 7: периодические на 7 и 14, финальный на 15 (останов сторожа).
        e2e::check(last_step("p15", 15, 7) == 15, "лимит 15: последняя запись на шаге 15 (финальная)");
        // Лимит 14 при периоде 7: 14 одновременно периодический и финальный.
        e2e::check(last_step("p14", 14, 7) == 14, "лимит 14: последняя запись на шаге 14");
        // Лимит 6 при периоде 7: периодических нет, есть только финальная на шаге 6.
        e2e::check(last_step("p6", 6, 7) == 6, "лимит 6: периодических нет, финальная на шаге 6");
        // Без периода (0) файл не пишется вообще.
        e2e::check(last_step("p0", 10, 0) == -1, "период 0: файл не пишется");

        // Периодическая запись действительно была на шаге 14: после лимита 15 в .prev
        // лежит предыдущая запись (шаг 14), а в основном файле -- 15.
        const std::string ck = root + "/p15.ckpt";
        e2e::check(fs::exists(ck + ".prev"), "лимит 15: есть .prev (запись на 14 была)");
        checkpoint::Payload pl;
        e2e::check(try_read(root + "/r_prev", ck + ".prev", false, [](SimulationConfig&) {}, &pl).empty() &&
                       pl.loop.step == 14,
                   ".prev содержит шаг 14");
    });

    // 5. Ротация
    run_case("5. ротация: .ckpt и .ckpt.prev, без временного файла", [] {
        e2e::ScopedWorkDir wd("sph_ckpt", "rot");
        const std::string root = wd.dir().string();
        const std::string ck = root + "/state.ckpt";
        make_checkpoint(root + "/w", ck, 9, 3); // записи: 3, 6, 9 (9 -- и периодическая, и финальная)
        e2e::check(fs::exists(ck), "основной файл есть");
        e2e::check(fs::exists(ck + ".prev"), ".prev есть");
        int tmp = 0;
        for (const auto& e : fs::directory_iterator(root))
        {
            const std::string n = e.path().filename().string();
            if (n.find(".tmp") != std::string::npos) ++tmp;
        }
        e2e::check(tmp == 0, "временных .tmp файлов нет");
        checkpoint::Payload a, b;
        e2e::check(try_read(root + "/r1", ck, false, [](SimulationConfig&) {}, &a).empty() && a.loop.step == 9,
                   "основной файл -- шаг 9");
        e2e::check(try_read(root + "/r2", ck + ".prev", false, [](SimulationConfig&) {}, &b).empty() && b.loop.step == 6,
                   ".prev -- шаг 6");
    });

#ifdef SPH_ENERGY_LOG
    // 6. Продолжение CSV-лога энергии при рестарте
    run_case("6. лог энергии: рестарт продолжает CSV без второго заголовка и дублей", [] {
        e2e::ScopedWorkDir wd("sph_ckpt", "elog");
        const std::string root = wd.dir().string();

        // Непрерывный прогон 2N, лог каждый шаг.
        std::vector<std::string> cont;
        {
            CkptScene sc = build_checkpoint_scene(root + "/A", 0.002, 0);
            sc.cfg.energy_log_every = 1;
            sc.cfg.energy_log_file = "energy.csv";
            double time = 0.0;
            e2e::check(run_steps(sc, 2 * N, time), "A: 2N шагов выполнены");
            cont = read_lines(root + "/A/energy.csv");
        }

        // Часть 1 и рестарт -- в ОДНОМ каталоге, как при реальном продолжении расчёта.
        std::vector<std::string> split;
        {
            const std::string dir = root + "/B";
            const std::string ck = dir + "/state.ckpt";
            {
                CkptScene sc = build_checkpoint_scene(dir, 0.002, 0);
                sc.cfg.energy_log_every = 1;
                sc.cfg.energy_log_file = "energy.csv";
                sc.cfg.checkpoint_file = ck;
                sc.cfg.checkpoint_every_steps = N;
                double time = 0.0;
                e2e::check(run_steps(sc, N, time), "B1: N шагов выполнены");
            }
            const std::size_t rows_after_1 = read_lines(dir + "/energy.csv").size();
            {
                // build_checkpoint_scene переписывает файлы сцены, но energy.csv не трогает.
                CkptScene sc = build_checkpoint_scene(dir, 0.002, 0);
                sc.cfg.energy_log_every = 1;
                sc.cfg.energy_log_file = "energy.csv";
                sc.cfg.restart_file = ck;
                double time = 0.0;
                e2e::check(run_steps(sc, 2 * N, time), "B2: рестарт дошёл до 2N шагов");
            }
            split = read_lines(dir + "/energy.csv");
            e2e::check(rows_after_1 == N + 2, "после части 1: заголовок + строки шагов 0..N");
        }

        std::size_t headers = 0;
        for (const auto& l : split)
            if (l.rfind("step,", 0) == 0) ++headers;
        e2e::check(headers == 1, "ровно один заголовок");

        // Номера шагов строго возрастают без повторов: 0,1,...,2N.
        bool seq = true;
        for (std::size_t k = 1; k < split.size(); ++k)
            if (std::stoull(split[k].substr(0, split[k].find(','))) != k - 1) seq = false;
        e2e::check(seq, "шаги в CSV идут 0,1,...,2N без дублей строк 0 и N");

        bool same = cont.size() == split.size();
        if (!same) std::printf("         строк: непрерывный %zu, рестарт %zu\n", cont.size(), split.size());
        int diffs = 0;
        for (std::size_t k = 0; k < cont.size() && k < split.size(); ++k)
            if (cont[k] != split[k] && diffs++ < 5)
                std::printf("         строка %zu:\n           A: %s\n           B: %s\n", k, cont[k].c_str(), split[k].c_str());
        e2e::check(same && diffs == 0, "CSV после рестарта совпадает с непрерывным построчно (в т.ч. residual)");
    });
#endif

    // 7. Проекция контакта (формат v2): таблица тёплого старта, поля конфига, чтение файла v1
    run_case("7. v2: proj_warm/proj_h_prev и поля проекции round-trip; файл v1 читается", [] {
        e2e::ScopedWorkDir wd("sph_ckpt", "proj");
        const std::string root = wd.dir().string();
        const std::string ck = root + "/proj.ckpt";
        const auto set_proj = [](SimulationConfig& c) {
            c.fiber_fiber_contact_type = 3;
            c.fiber_fiber_contact_predictor = 1;
            c.fiber_contact_projection_iters = 17;
            c.fiber_contact_projection_warm = 0;
        };
        const auto make_entries = [] {
            std::vector<interactions::NormalWarmEntry> v(2);
            v[0].key = interactions::make_ff_pair_key(7, 8, 3, -1);
            v[0].P = 1.25e-7;
            v[1].key = interactions::make_fs_pair_key(4, 5, 99);
            v[1].P = -0.0;
            return v;
        };
        // Запись без интегрирования: состояние сцены сразу после построения.
        const auto write_ck = [&](const std::string& wdir, const std::string& path, bool empty) {
            CkptScene sc = build_checkpoint_scene(wdir, 0.5, 0);
            set_proj(sc.cfg);
            checkpoint::Payload pl;
            if (!empty) pl.proj_warm = make_entries();
            pl.proj_h_prev = empty ? 0.0 : 3.5e-9;
            checkpoint::write(path, sc.cfg, sc.model, sc.ntotal, pl);
        };

        write_ck(root + "/w0", ck, false);
        checkpoint::Payload out;
        e2e::check(try_read(root + "/w1", ck, false, set_proj, &out).empty(), "v2 читается");
        const auto ref = make_entries();
        bool same = out.proj_warm.size() == ref.size() && out.proj_h_prev == 3.5e-9;
        for (std::size_t k = 0; same && k < ref.size(); ++k)
            same = out.proj_warm[k].key == ref[k].key && out.proj_warm[k].P == ref[k].P;
        e2e::check(same, "proj_warm и proj_h_prev совпадают побитово");

        // Пустая таблица
        write_ck(root + "/w2", root + "/empty.ckpt", true);
        checkpoint::Payload out_e;
        out_e.proj_warm = make_entries();
        e2e::check(try_read(root + "/w3", root + "/empty.ckpt", false, set_proj, &out_e).empty() &&
                       out_e.proj_warm.empty() && out_e.proj_h_prev == 0.0,
                   "пустая таблица round-trip");

        // Смена поля проекции без --force запрещена
        const std::string e = try_read(root + "/w4", ck, false, [&](SimulationConfig& c) {
            set_proj(c);
            c.fiber_contact_projection_iters = 8;
        });
        e2e::check(e.find("fiber_contact_projection_iters") != std::string::npos, "смена iters без force -> ошибка");

        // Файл v1: версия 1, из CONF убраны 4 последних int, блоки PROJ и USER удалены, CRC пересчитан.
        std::vector<char> raw = slurp(ck);
        auto rd64 = [&](std::size_t off) { std::uint64_t v; std::memcpy(&v, raw.data() + off, 8); return v; };
        std::size_t pos = 24;
        const std::uint64_t conf_len = rd64(pos + 4);
        const std::vector<std::size_t> starts = block_starts(raw);
        const std::size_t proj_pos = starts[starts.size() - 2]; // предпоследний блок -- PROJ, последний -- USER
        raw.resize(proj_pos);      // отрезали PROJ, USER и CRC
        const std::uint64_t nl = conf_len - 16;
        std::memcpy(raw.data() + pos + 4, &nl, 8);
        std::uint32_t cnt;
        std::memcpy(&cnt, raw.data() + pos + 12, 4);
        cnt -= 4;
        std::memcpy(raw.data() + pos + 12, &cnt, 4);
        raw.erase(raw.begin() + static_cast<std::ptrdiff_t>(pos + 12 + conf_len - 16),
                  raw.begin() + static_cast<std::ptrdiff_t>(pos + 12 + conf_len));
        const std::uint32_t one = 1;
        std::memcpy(raw.data() + 8, &one, 4);
        std::uint32_t crc = 0xFFFFFFFFu;
        for (unsigned char c : raw)
        {
            crc ^= c;
            for (int b = 0; b < 8; ++b) crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
        }
        crc = ~crc;
        raw.resize(raw.size() + 4);
        std::memcpy(raw.data() + raw.size() - 4, &crc, 4);
        spit(root + "/v1.ckpt", raw);

        checkpoint::Payload out1;
        checkpoint::ConfigDiff diff1;
        const std::string e1 = try_read(root + "/w5", root + "/v1.ckpt", false, set_proj, &out1, &diff1);
        e2e::check(e1.empty(), "файл v1 читается" + (e1.empty() ? std::string() : ": " + e1));
        e2e::check(out1.proj_warm.empty() && diff1.forbidden.empty(), "v1: таблицы нет, различий новых полей нет");
    });

    // 8. Блок USER (формат v3): состояние вызывающего round-trip; файл v2 (без USER) читается с пустым состоянием
    run_case("8. v3: user_state round-trip (в т.ч. нулевые байты); файл v2 читается", [] {
        e2e::ScopedWorkDir wd("sph_ckpt", "user");
        const std::string root = wd.dir().string();
        const std::string ck = root + "/user.ckpt";
        std::string blob("MONx\0\x01\xff", 7);
        blob += std::string(1000, '\x7f');
        {
            CkptScene sc = build_checkpoint_scene(root + "/w0", 0.5, 0);
            checkpoint::Payload pl;
            pl.user_state = blob;
            checkpoint::write(ck, sc.cfg, sc.model, sc.ntotal, pl);
        }
        checkpoint::Payload out;
        const std::string e = try_read(root + "/w1", ck, false, [](SimulationConfig&) {}, &out);
        e2e::check(e.empty(), "v3 читается" + (e.empty() ? std::string() : ": " + e));
        e2e::check(out.user_state == blob, "user_state совпадает побайтно");

        // Файл v2: версия 2, блок USER удалён, CRC пересчитан.
        std::vector<char> raw = slurp(ck);
        const std::vector<std::size_t> starts = block_starts(raw);
        raw.resize(starts.back());
        const std::uint32_t two = 2;
        std::memcpy(raw.data() + 8, &two, 4);
        const std::uint32_t crc = crc32_of(raw);
        raw.resize(raw.size() + 4);
        std::memcpy(raw.data() + raw.size() - 4, &crc, 4);
        spit(root + "/v2.ckpt", raw);
        checkpoint::Payload out2;
        out2.user_state = "мусор";
        const std::string e2 = try_read(root + "/w2", root + "/v2.ckpt", false, [](SimulationConfig&) {}, &out2);
        e2e::check(e2.empty(), "файл v2 читается" + (e2.empty() ? std::string() : ": " + e2));
        e2e::check(out2.user_state.empty(), "v2: user_state пуст");
    });

    std::printf("\nпровалено кейсов: %d\n", g_failed_cases);
    return (g_failed_cases == 0 && e2e::failures == 0) ? 0 : 1;
}
