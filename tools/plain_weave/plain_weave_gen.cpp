// Генератор геометрии полотняного переплетения в формате .fib
// Использование:
//   plain_weave_gen --out=<файл.fib> [--layout=<файл>] [--preset=<имя>]
//                   [--d=<м>] [--pitch=<м>] [--amp=<м>]
//                   [--nx=<чис>] [--ny=<чис>] [--seg=<чис>]
//                   [--overhang=<м>] [--overhang-x=<м>] [--overhang-y=<м>]
//                   [--profile=cos|arcs] [--center=cx,cy,cz]
//                   [--clamp=none|x|y|xy] [--clamp-edge=<м>]
//                   [--E=<Па>] [--nu=] [--rho=<кг/м³>] [--eps-break=]
//                   [--mu-ff=] [--mu-fs=] [--allow-penetration] [--help]
//
// Пример:
//   plain_weave_gen --preset=test --out=weave.fib --layout=weave.layout.txt
//   plain_weave_gen --preset=twaron_ct717 --out=twaron.fib --clamp-edge=0.06

#include <fiber/plain_weave.hpp>
#include <chrono>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace
{

// Строгий парсинг double: есть ли в конце строки мусор?
bool parse_double(const char* str, double& out)
{
    if (!str || *str == '\0') return false;
    char* endptr = nullptr;
    out = std::strtod(str, &endptr);
    // Если endptr не указывает на конец строки, значит был мусор
    return endptr != nullptr && *endptr == '\0';
}

// Строгий парсинг int
bool parse_int(const char* str, int& out)
{
    if (!str || *str == '\0') return false;
    char* endptr = nullptr;
    long val = std::strtol(str, &endptr, 10);
    out = static_cast<int>(val);
    return endptr != nullptr && *endptr == '\0';
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
            "  --help                   Эта справка\n");
}

} // namespace

int main(int argc, char** argv)
{
    auto start_time = std::chrono::steady_clock::now();

    fiber::PlainWeaveParams q;
    std::string out_file, layout_file, preset_name;

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
    const double unset = std::numeric_limits<double>::quiet_NaN();
    double overhang_both = unset, overhang_x = unset, overhang_y = unset;

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
    // Отрицательный выступ в PlainWeaveParams означает «по умолчанию p/2»,
    // поэтому явный отрицательный ввод отсекается здесь.
    for (double e : {overhang_both, overhang_x, overhang_y})
        if (e < 0.0) // NaN (не задан) сюда не проходит
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
    if (!std::isnan(overhang_both)) q.overhang_x = q.overhang_y = overhang_both;
    if (!std::isnan(overhang_x)) q.overhang_x = overhang_x;
    if (!std::isnan(overhang_y)) q.overhang_y = overhang_y;

    // Построить полотно
    fiber::PlainWeave W;
    try
    {
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

    // Напечатать сводку и время
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    for (const auto& w : W.warnings) printf("%s\n", w.c_str());
    fiber::print_plain_weave_summary(stdout, W);
    printf("файл: %s\n", out_file.c_str());
    printf("раскладка: %s\n", layout_file.c_str());
    printf("время генерации: %lld мс\n", (long long)duration.count());

    return 0;
}
