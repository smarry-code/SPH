// tests/test_e2e_clamp.cpp
// Сквозные (e2e) тесты кинематического закрепления Solid и Fiber против
// аналитических решений (постановка Clamp/CLAMP_TESTS.md; механизм --
// include/interactions/clamp.hpp, точки [1]-[3] в src/time_integration.cpp).
// Дополняют модульные и интеграционные тесты tests/test_clamp.cpp (U1-U6,
// F1-F5, S1-S4): там проверяется «частицы стоят», здесь -- что зажим ведёт
// себя как опора (реакция, энергия, отражение волн).
//
// Все сцены идут ТОЛЬКО через боевой time_integration(cfg, fiber::Model&, ...)
// одним непрерывным прогоном; монитор §1.2 снимается наблюдателем шага
// (RunGuard::set_step_observer, tests/clamp_e2e_common.hpp).
//
// Один запуск -- один вариант:
//   test_e2e_clamp --case=ES1 [--solver=pc|euler] [--full] [--variant=...]
//                  [--time-limit=СЕК] [--keep-workdir] [--monitor=ФАЙЛ.csv]
//                  [--frames-dir=DIR [--frames=N]]   (кадры для tools/clamp_animate.py)
// Режим по умолчанию -- smoke (ctest, метка e2e); --full -- полный (метка slow).
// Код возврата: 0 -- все проверки прошли, 1 -- есть провалы, 2 -- ошибка ключей.

#include "clamp_e2e_common.hpp"
#include "clamp_e2e_fiber.hpp"
#include "clamp_e2e_solid.hpp"

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <map>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

using CaseFn = void (*)(const clamp_e2e::Options&);

const std::map<std::string, CaseFn>& cases()
{
    using namespace clamp_e2e;
    static const std::map<std::string, CaseFn> m = {
        {"ES1", &run_ES1}, {"ES2", &run_ES2}, {"ES3", &run_ES3}, {"ES4", &run_ES4},
        {"EF1", &run_EF1}, {"EF2", &run_EF2}, {"EF3", &run_EF3}, {"EF4", &run_EF4},
        {"EF5", &run_EF5}, {"EF6", &run_EF6}, {"EK1", &run_EK1},
    };
    return m;
}

void usage()
{
    std::printf(
        "Использование: test_e2e_clamp --case=ИМЯ [ключи]\n"
        "  --case=ИМЯ        ES1 ES2 ES3 ES4 EF1 EF2 EF3 EF4 EF5 EF6 EK1\n"
        "  --solver=pc|euler решатель (по умолчанию pc)\n"
        "  --full            полный режим (по умолчанию smoke)\n"
        "  --variant=ИМЯ     вариант сцены (см. шапку кейса)\n"
        "  --time-limit=СЕК  бюджет на один прогон (по умолчанию 600)\n"
        "  --keep-workdir    не удалять временный каталог прогона\n"
        "  --monitor=ФАЙЛ    записать монитор §1.2 в CSV\n"
        "  --frames-dir=DIR  сохранить кадры для tools/clamp_animate.py в DIR/<прогон>\n"
        "  --frames=N        число кадров на прогон (по умолчанию 80)\n");
}

} // namespace

int main(int argc, char** argv)
{
    clamp_e2e::Options opt;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto val = [&](const char* key) { return a.substr(std::strlen(key)); };
        if (a.rfind("--case=", 0) == 0) opt.case_name = val("--case=");
        else if (a == "--solver=pc") opt.solver = 0;
        else if (a == "--solver=euler") opt.solver = 1;
        else if (a == "--full") opt.full = true;
        else if (a.rfind("--variant=", 0) == 0) opt.variant = val("--variant=");
        else if (a == "--keep-workdir") opt.keep_workdir = true;
        else if (a.rfind("--monitor=", 0) == 0) opt.monitor_csv = val("--monitor=");
        else if (a.rfind("--frames-dir=", 0) == 0)
            // Абсолютный путь: run_case копирует кадры, находясь во временном каталоге.
            opt.frames_dir = std::filesystem::absolute(val("--frames-dir=")).string();
        else if (a.rfind("--frames=", 0) == 0)
        {
            try { opt.frames = static_cast<std::size_t>(std::stoul(val("--frames="))); }
            catch (...) { std::fprintf(stderr, "Неверное значение %s\n", a.c_str()); return 2; }
        }
        else if (a.rfind("--time-limit=", 0) == 0)
        {
            try { opt.time_limit = std::stod(val("--time-limit=")); }
            catch (...) { std::fprintf(stderr, "Неверное значение %s\n", a.c_str()); return 2; }
        }
        else if (a == "--help" || a == "-h") { usage(); return 0; }
        else { std::fprintf(stderr, "Неизвестный ключ: %s\n", a.c_str()); usage(); return 2; }
    }
    const auto it = cases().find(opt.case_name);
    if (it == cases().end())
    {
        std::fprintf(stderr, "Неизвестный или пустой --case=%s\n", opt.case_name.c_str());
        usage();
        return 2;
    }
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    std::printf("=== e2e закрепления: %s, решатель %s, режим %s%s, потоков %d ===\n", opt.case_name.c_str(),
                clamp_e2e::solver_name(opt.solver), opt.full ? "full" : "smoke",
                opt.variant.empty() ? "" : (", вариант " + opt.variant).c_str(), threads);
    try
    {
        it->second(opt);
    }
    catch (const std::exception& ex)
    {
        std::printf("  FAIL : исключение вне прогона: %s\n", ex.what());
        ++e2e::failures;
    }
    std::printf("\n=== %s: %s (провалов %d) ===\n", opt.case_name.c_str(), e2e::failures == 0 ? "PASS" : "FAIL",
                e2e::failures);
    return e2e::failures == 0 ? 0 : 1;
}
