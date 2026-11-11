// Генератор геометрии полотняного переплетения в формате .fib
// Использование:
//   plain_weave_gen --out=<файл.fib> [--layout=<файл>] [--preset=<имя>]
//                   [--d=<м>] [--pitch=<м>] [--amp=<м>]
//                   [--nx=<чис>] [--ny=<чис>] [--seg=<чис>]
//                   [--overhang=<м>] [--overhang-x=<м>] [--overhang-y=<м>]
//                   [--profile=cos|arcs] [--center=cx,cy,cz]
//                   [--clamp=none|x|y|xy] [--clamp-edge=<м>]
//                   [--E=<Па>] [--nu=] [--rho=<кг/м³>] [--eps-break=]
//                   [--mu-ff=] [--mu-fs=] [--allow-penetration]
//                   [--patch=NXxNY --tile=TXxTY [--allow-odd-patch]] [--help]
//
// Пример:
//   plain_weave_gen --preset=test --out=weave.fib --layout=weave.layout.txt
//   plain_weave_gen --preset=twaron_ct717 --out=twaron.fib --clamp-edge=0.06
//   plain_weave_gen --preset=twaron_ct717 --patch=34x34 --tile=3x3 --out=quilt.fib

#include <fiber/plain_weave.hpp>
#include <fiber/plain_weave_quilt.hpp>
#include <chrono>
#include <cerrno>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

namespace
{

// Строгий парсинг double: есть ли в конце строки мусор?
bool parse_double(const char* str, double& out)
{
    if (!str || *str == '\0') return false;
    char* endptr = nullptr;
    errno = 0;
    const double v = std::strtod(str, &endptr);
    // Мусор в конце, переполнение/потеря точности (ERANGE), inf и nan -- ошибка.
    if (endptr == nullptr || *endptr != '\0' || endptr == str || errno == ERANGE || !std::isfinite(v))
        return false;
    out = v;
    return true;
}

// Строгий парсинг int
bool parse_int(const char* str, int& out)
{
    if (!str || *str == '\0') return false;
    char* endptr = nullptr;
    errno = 0;
    const long val = std::strtol(str, &endptr, 10);
    if (endptr == nullptr || *endptr != '\0' || endptr == str || errno == ERANGE) return false;
    if (val < std::numeric_limits<int>::min() || val > std::numeric_limits<int>::max()) return false;
    out = static_cast<int>(val);
    return true;
}

// Парсинг center=cx,cy,cz
bool parse_center(const char* str, double& cx, double& cy, double& cz)
{
    std::string s(str);
    size_t pos1 = s.find(',');
    if (pos1 == std::string::npos) return false;
    size_t pos2 = s.find(',', pos1 + 1);
    if (pos2 == std::string::npos) return false;
    return parse_double(s.substr(0, pos1).c_str(), cx) &&
           parse_double(s.substr(pos1 + 1, pos2 - pos1 - 1).c_str(), cy) &&
           parse_double(s.substr(pos2 + 1).c_str(), cz);
}

// Парсинг "AxB" (два целых через 'x').
bool parse_pair(const char* str, int& a, int& b)
{
    const std::string s(str);
    const std::size_t pos = s.find('x');
    if (pos == std::string::npos) return false;
    return parse_int(s.substr(0, pos).c_str(), a) && parse_int(s.substr(pos + 1).c_str(), b);
}

void print_help()
{
    fprintf(stderr,
            "Генератор геометрии полотняного переплетения (plain weave)\n\n"
            "Использование:\n"
            "  plain_weave_gen --out=<файл.fib> [параметры]\n\n"
            "Обязательные:\n"
            "  --out=<файл.fib>         Выходной файл нитей\n\n"
            "Опциональные:\n"
            "  --preset=<имя>           Пресет: test, circle_draw, twaron_ct717\n"
            "  --layout=<файл>          Раскладка (по умолчанию <out>.layout.txt)\n"
            "  --d=<м>                  Диаметр нити\n"
            "  --pitch=<м>              Шаг между нитями\n"
            "  --amp=<м>                Амплитуда волнистости\n"
            "  --nx=<n>                 Нитей основы\n"
            "  --ny=<n>                 Нитей утка\n"
            "  --seg=<k>                Сегментов на полуволну (>=2)\n"
            "  --overhang=<м>           Выступ (оба направления)\n"
            "  --overhang-x=<м>         Выступ основы\n"
            "  --overhang-y=<м>         Выступ утка\n"
            "  --profile=cos|arcs       Профиль нити\n"
            "  --center=cx,cy,cz        Центр полотна\n"
            "  --clamp=none|x|y|xy      Зажим\n"
            "  --clamp-edge=<м>         Граница зажима c_e\n"
            "  --E=<Па>                 Модуль Юнга\n"
            "  --nu=<>                  Коэффициент Пуассона\n"
            "  --rho=<кг/м³>            Плотность\n"
            "  --eps-break=<>           Деформация разрыва\n"
            "  --mu-ff=<>               Трение нить–нить\n"
            "  --mu-fs=<>               Трение нить–Solid\n"
            "  --allow-penetration      Разрешить проникание\n"
            "  --patch=NXxNY            Одеяло: патч NX нитей основы x NY нитей утка (например 34x34)\n"
            "  --tile=TXxTY             Одеяло: патчей вдоль x и вдоль y (например 3x3); нужны оба ключа,\n"
            "                           несовместимо с --nx/--ny; пишется <out>.patches.txt\n"
            "  --allow-odd-patch        Разрешить нечётный патч при tile > 1 (с предупреждением)\n"
            "  --help                   Эта справка\n");
}

} // namespace

int main(int argc, char** argv)
{
    auto start_time = std::chrono::steady_clock::now();

    fiber::PlainWeaveParams q;
    std::string out_file, layout_file, preset_name;
    int patch_nx = 0, patch_ny = 0, tile_x = 0, tile_y = 0;
    bool patch_set = false, tile_set = false, allow_odd_patch = false, nx_set = false, ny_set = false;

    // Пресет применяется ДО явных параметров (явный параметр сильнее пресета,
    // ТЗ §4.2), поэтому сначала ищем --preset среди всех аргументов.
    for (int i = 1; i < argc; ++i)
        if (std::strncmp(argv[i], "--preset=", 9) == 0) preset_name = argv[i] + 9;
    if (!preset_name.empty())
    {
        try
        {
            fiber::apply_plain_weave_preset(preset_name, q);
        }
        catch (const std::exception& e)
        {
            fprintf(stderr, "%s\n", e.what());
            return 2;
        }
    }

    // --overhang= задаёт оба выступа, --overhang-x/-y сильнее него при любом
    // порядке аргументов: общий выступ применяется после разбора.
    double overhang_both = 0.0, overhang_x = 0.0, overhang_y = 0.0;
    bool overhang_both_set = false, overhang_x_set = false, overhang_y_set = false;

    // Таблица числовых ключей: имя -> куда писать.
    struct DoubleKey { const char* name; double* dst; };
    const DoubleKey dkeys[] = {
        {"--d=", &q.d},           {"--pitch=", &q.pitch},       {"--amp=", &q.amp},
        {"--overhang=", &overhang_both}, {"--overhang-x=", &overhang_x}, {"--overhang-y=", &overhang_y},
        {"--clamp-edge=", &q.clamp_edge},
        {"--E=", &q.E},           {"--nu=", &q.nu},             {"--rho=", &q.rho},
        {"--eps-break=", &q.eps_break}, {"--mu-ff=", &q.mu_ff}, {"--mu-fs=", &q.mu_fs},
    };
    struct IntKey { const char* name; int* dst; };
    const IntKey ikeys[] = {{"--nx=", &q.nx}, {"--ny=", &q.ny}, {"--seg=", &q.seg}};

    int failures = 0;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        auto value_of = [&arg](const char* key) -> const char* {
            const std::size_t n = std::strlen(key);
            return arg.compare(0, n, key) == 0 ? arg.c_str() + n : nullptr;
        };

        if (arg == "--help")
        {
            print_help();
            return 0;
        }
        if (arg == "--allow-penetration")
        {
            q.allow_penetration = true;
            continue;
        }
        if (arg == "--allow-odd-patch")
        {
            allow_odd_patch = true;
            continue;
        }
        if (const char* v = value_of("--patch="))
        {
            patch_set = true;
            if (!parse_pair(v, patch_nx, patch_ny) || patch_nx < 1 || patch_ny < 1)
            {
                fprintf(stderr, "Ошибка: --patch=%s (ожидается NXxNY, целые >= 1)\n", v);
                ++failures;
            }
            continue;
        }
        if (const char* v = value_of("--tile="))
        {
            tile_set = true;
            if (!parse_pair(v, tile_x, tile_y) || tile_x < 1 || tile_y < 1)
            {
                fprintf(stderr, "Ошибка: --tile=%s (ожидается TXxTY, целые >= 1)\n", v);
                ++failures;
            }
            continue;
        }
        if (value_of("--nx=")) nx_set = true;
        if (value_of("--ny=")) ny_set = true;
        if (const char* v = value_of("--out=")) { out_file = v; continue; }
        if (const char* v = value_of("--layout=")) { layout_file = v; continue; }
        if (value_of("--preset=")) continue; // уже применён

        bool handled = false;
        for (const auto& k : dkeys)
            if (const char* v = value_of(k.name))
            {
                handled = true;
                if (!parse_double(v, *k.dst))
                {
                    fprintf(stderr, "Ошибка: не число в %s\n", arg.c_str());
                    ++failures;
                }
                if (k.dst == &overhang_both) overhang_both_set = true;
                if (k.dst == &overhang_x) overhang_x_set = true;
                if (k.dst == &overhang_y) overhang_y_set = true;
                if (k.dst == &q.clamp_edge) q.clamp_edge_set = true;
                if (k.dst == &q.mu_fs) q.mu_fs_set = true;
                break;
            }
        if (!handled)
            for (const auto& k : ikeys)
                if (const char* v = value_of(k.name))
                {
                    handled = true;
                    if (!parse_int(v, *k.dst))
                    {
                        fprintf(stderr, "Ошибка: не целое в %s\n", arg.c_str());
                        ++failures;
                    }
                    break;
                }
        if (handled) continue;

        if (const char* v = value_of("--profile="))
        {
            const std::string prof = v;
            if (prof == "cos") q.profile = fiber::PlainWeaveProfile::Cos;
            else if (prof == "arcs") q.profile = fiber::PlainWeaveProfile::Arcs;
            else
            {
                fprintf(stderr, "Неизвестный профиль: %s (cos или arcs)\n", v);
                ++failures;
            }
            continue;
        }
        if (const char* v = value_of("--center="))
        {
            if (!parse_center(v, q.cx, q.cy, q.cz))
            {
                fprintf(stderr, "Ошибка: --center=%s (ожидается cx,cy,cz)\n", v);
                ++failures;
            }
            continue;
        }
        if (const char* v = value_of("--clamp="))
        {
            const std::string c = v;
            if (c == "none") q.clamp = fiber::PlainWeaveClamp::None;
            else if (c == "x") q.clamp = fiber::PlainWeaveClamp::X;
            else if (c == "y") q.clamp = fiber::PlainWeaveClamp::Y;
            else if (c == "xy") q.clamp = fiber::PlainWeaveClamp::XY;
            else
            {
                fprintf(stderr, "Неизвестный режим зажима: %s (none, x, y или xy)\n", v);
                ++failures;
            }
            continue;
        }

        fprintf(stderr, "Неизвестный ключ: %s\n", arg.c_str());
        ++failures;
    }

    if (out_file.empty())
    {
        fprintf(stderr, "Ошибка: необходимо --out=<файл.fib>\n");
        ++failures;
    }
    const bool quilt = patch_set || tile_set;
    if (quilt && !(patch_set && tile_set))
    {
        fprintf(stderr, "Ошибка: --patch и --tile задаются вместе\n");
        ++failures;
    }
    if (quilt && (nx_set || ny_set))
    {
        fprintf(stderr, "Ошибка: --patch/--tile несовместимы с --nx/--ny\n");
        ++failures;
    }
    if (allow_odd_patch && !quilt)
    {
        fprintf(stderr, "Ошибка: --allow-odd-patch требует --patch и --tile\n");
        ++failures;
    }
    // Отрицательный выступ в PlainWeaveParams означает «по умолчанию p/2»,
    // поэтому явный отрицательный ввод отсекается здесь.
    for (const auto& oh : {std::make_pair(overhang_both_set, overhang_both),
                           std::make_pair(overhang_x_set, overhang_x),
                           std::make_pair(overhang_y_set, overhang_y)})
        if (oh.first && oh.second < 0.0)
        {
            fprintf(stderr, "Ошибка: выступ должен быть >= 0\n");
            ++failures;
            break;
        }
    if (failures > 0)
    {
        print_help();
        return 2;
    }
    if (overhang_both_set) q.overhang_x = q.overhang_y = overhang_both;
    if (overhang_x_set) q.overhang_x = overhang_x;
    if (overhang_y_set) q.overhang_y = overhang_y;

    // Построить полотно
    fiber::PlainWeave W;
    fiber::PlainWeaveQuilt Qlt;
    try
    {
        if (quilt)
        {
            fiber::QuiltParams QP;
            QP.base = q;
            QP.patch_nx = patch_nx;
            QP.patch_ny = patch_ny;
            QP.tile_x = tile_x;
            QP.tile_y = tile_y;
            QP.allow_odd_patch = allow_odd_patch;
            Qlt = fiber::build_plain_weave_quilt(QP);
            W = Qlt.W;
        }
        else
            W = fiber::build_plain_weave(q);
    }
    catch (const std::exception& e)
    {
        fprintf(stderr, "%s\n", e.what());
        return 1;
    }

    // Собрать строку команды
    std::ostringstream cmd_line;
    for (int i = 0; i < argc; ++i)
    {
        if (i > 0) cmd_line << " ";
        cmd_line << argv[i];
    }

    // Написать .fib
    try
    {
        std::ofstream fib_stream(out_file);
        if (!fib_stream)
        {
            fprintf(stderr, "Не удалось открыть файл для записи: %s\n", out_file.c_str());
            return 1;
        }
        fiber::write_plain_weave_fib(fib_stream, W, cmd_line.str());
        fib_stream.close();
    }
    catch (const std::exception& e)
    {
        fprintf(stderr, "%s\n", e.what());
        return 1;
    }

    // Написать раскладку
    if (layout_file.empty())
    {
        layout_file = out_file + ".layout.txt";
    }
    try
    {
        std::ofstream layout_stream(layout_file);
        if (!layout_stream)
        {
            fprintf(stderr, "Не удалось открыть файл раскладки для записи: %s\n", layout_file.c_str());
            return 1;
        }
        fiber::write_plain_weave_layout(layout_stream, W);
        layout_stream.close();
    }
    catch (const std::exception& e)
    {
        fprintf(stderr, "%s\n", e.what());
        return 1;
    }

    // Раскладка по патчам
    std::string patches_file;
    if (quilt)
    {
        patches_file = out_file + ".patches.txt";
        std::ofstream ps(patches_file);
        if (!ps)
        {
            fprintf(stderr, "Не удалось открыть файл патчей для записи: %s\n", patches_file.c_str());
            return 1;
        }
        fiber::write_plain_weave_quilt_patches(ps, Qlt);
    }

    // Напечатать сводку и время
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    for (const auto& w : W.warnings) printf("%s\n", w.c_str());
    fiber::print_plain_weave_summary(stdout, W);
    if (quilt)
    {
        const double p = W.params.pitch;
        const double Lpx = patch_ny * p, Lpy = patch_nx * p; // размер патча по x (patch_ny нитей утка) и по y
        double Lx = 0.0, Ly = 0.0;
        fiber::plain_weave_window(W.params, Lx, Ly);
        printf("патчей %dx%d по %dx%d нитей, размер патча %.9g x %.9g м, окно ткани %.9g x %.9g м\n", tile_x, tile_y,
               patch_nx, patch_ny, Lpx, Lpy, Lx, Ly);
        printf("патчи: %s\n", patches_file.c_str());
    }
    printf("файл: %s\n", out_file.c_str());
    printf("раскладка: %s\n", layout_file.c_str());
    printf("время генерации: %lld мс\n", (long long)duration.count());

    return 0;
}
