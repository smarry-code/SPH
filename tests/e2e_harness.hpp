#pragma once
// e2e_harness.hpp
// Обвязка сквозных тестов, не зависящая от кейса: счётчик проверок, изоляция
// побочных эффектов, глушение вывода конвейера, русское форматирование чисел.
//
// Происхождение. Это перенос разделов «3. Изоляция побочных эффектов» и части
// раздела 8 из tests/test_e2e_bc2.cpp — того же кода, проверенного на кейсе
// BC2. Сам test_e2e_bc2.cpp НЕ переключён на этот заголовок и продолжает жить
// со своими копиями: цель второго кейса — не переписывать работающий тест.
// Когда BC6 устоится, BC2 можно будет перевести сюда отдельной задачей.
// Из BC2 пока вынесена только геометрия (tests/bc2_geometry.hpp) — её
// переиспользует регрессионный тест tests/test_regression.cpp.
//
// Что параметризовано против оригинала: префикс имени временного каталога
// (был литерал "sph_e2e_bc2_").

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace e2e
{

namespace fs = std::filesystem;

// ===========================================================================
// Счётчик проверок
// ===========================================================================
// Единственное изменяемое глобальное состояние обвязки. inline-переменная
// (C++17), поэтому заголовок можно включать из нескольких единиц трансляции.
inline int failures = 0;

// Возвращает cond — чтобы вызывающий мог доложить подробности провала сразу
// после строки FAIL.
inline bool check(bool cond, const std::string& what)
{
    if (cond)
    {
        std::printf("  ok   : %s\n", what.c_str());
    }
    else
    {
        std::printf("  FAIL : %s\n", what.c_str());
        ++failures;
    }
    return cond;
}

inline void note(const std::string& what)
{
    std::printf("  info : %s\n", what.c_str());
}

// Предупреждение: не провал, но и не рядовая справка. Выделено отдельно,
// потому что у BC6 есть исходы, которые постановка прямо называет ожидаемым
// диагностическим результатом, а не сбоем (§5.4: отсутствие провалов на
// кривой). Такой исход обязан быть заметен в выводе и не обязан ронять тест.
inline void warn(const std::string& what)
{
    std::printf("  WARN : %s\n", what.c_str());
}

// ===========================================================================
// Изоляция побочных эффектов
// ===========================================================================
// time_integration() пишет log.txt в ТЕКУЩИЙ каталог жёстко зашитым именем,
// output() — sph_output_*.txt и sph_gp.csv туда же. Прогон в каталоге
// разработчика затёр бы лог реального расчёта. Поэтому каждый прогон уходит
// в собственный временный каталог, а после — проверка, что рабочий каталог
// не изменился (check_no_side_effects).
class ScopedWorkDir
{
public:
    ScopedWorkDir(const std::string& prefix, const std::string& tag)
        : prev_(fs::current_path())
    {
        std::ostringstream os;
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        os << prefix << "_" << tag << "_" << static_cast<long long>(stamp)
           << "_" << counter_++;
        dir_ = fs::temp_directory_path() / os.str();
        fs::remove_all(dir_);
        fs::create_directories(dir_);
        fs::current_path(dir_);
    }
    ~ScopedWorkDir()
    {
        std::error_code ec;
        fs::current_path(prev_, ec);
        if (!keep_) fs::remove_all(dir_, ec);
    }
    ScopedWorkDir(const ScopedWorkDir&) = delete;
    ScopedWorkDir& operator=(const ScopedWorkDir&) = delete;

    const fs::path& dir() const { return dir_; }
    void keep() { keep_ = true; }

private:
    fs::path prev_;
    fs::path dir_;
    bool keep_ = false;
    inline static int counter_ = 0;
};

// Глушение вывода конвейера: single_step/output/time_integration печатают
// в std::wcout/std::cout десятки строк на кадр. Без этого отчёт теста
// нечитаем.
class ScopedSilence
{
public:
    explicit ScopedSilence(bool active) : active_(active)
    {
        if (!active_) return;
        old_w_ = std::wcout.rdbuf(sink_w_.rdbuf());
        old_n_ = std::cout.rdbuf(sink_n_.rdbuf());
    }
    ~ScopedSilence()
    {
        if (!active_) return;
        std::wcout.rdbuf(old_w_);
        std::cout.rdbuf(old_n_);
    }
    ScopedSilence(const ScopedSilence&) = delete;
    ScopedSilence& operator=(const ScopedSilence&) = delete;

private:
    bool active_;
    std::wostringstream sink_w_;
    std::ostringstream  sink_n_;
    std::wstreambuf* old_w_ = nullptr;
    std::streambuf*  old_n_ = nullptr;
};

// Снимок берётся ТОЛЬКО по именам, которые пишет конвейер. Сравнивать весь
// каталог нельзя: при `ctest -j` соседние тесты создают в том же каталоге свои
// файлы, и проверка стала бы флакающей.
inline bool is_pipeline_artifact(const std::string& n)
{
    return n == "log.txt" || n == "sph_gp.csv" || n == "force_timing.csv" ||
           (n.rfind("sph_output_", 0) == 0);
}

struct DirSnapshot
{
    // имя -> размер; отсутствующий файл в отображении просто отсутствует
    std::vector<std::pair<std::string, std::uintmax_t>> artifacts;
};

inline DirSnapshot snapshot_dir(const fs::path& d)
{
    DirSnapshot s;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(d, ec))
    {
        const std::string n = e.path().filename().string();
        if (!is_pipeline_artifact(n)) continue;
        std::error_code ec2;
        s.artifacts.emplace_back(n, fs::file_size(e.path(), ec2));
    }
    std::sort(s.artifacts.begin(), s.artifacts.end());
    return s;
}

inline void check_no_side_effects(const fs::path& cwd, const DirSnapshot& before)
{
    std::printf("\n[§7] Побочные эффекты в рабочем каталоге\n");
    const DirSnapshot after = snapshot_dir(cwd);
    {
        std::ostringstream os;
        os << "артефактов конвейера в рабочем каталоге: было " << before.artifacts.size()
           << ", стало " << after.artifacts.size();
        note(os.str());
    }
    check(before.artifacts == after.artifacts,
          "log.txt / sph_output_*.txt / sph_gp.csv в рабочем каталоге не созданы и не изменены");
    if (before.artifacts != after.artifacts)
    {
        for (const auto& a : after.artifacts)
        {
            const auto it = std::find_if(before.artifacts.begin(), before.artifacts.end(),
                                         [&](const auto& b) { return b == a; });
            if (it == before.artifacts.end())
            {
                std::printf("         появился/изменился: %s (%llu байт)\n",
                            a.first.c_str(), static_cast<unsigned long long>(a.second));
            }
        }
    }
}

// ===========================================================================
// Русское форматирование чисел
// ===========================================================================
// Запятая как десятичный разделитель, неразрывные пробелы по три разряда.
// График и отчёт русскоязычные, и «0.2 мм» рядом с «0,2 мм» в тексте выглядит
// как разные величины.
inline std::string ru_num(double value, int precision)
{
    std::ostringstream s;
    s.setf(std::ios::fixed);
    s.precision(precision);
    s << value;
    std::string out = s.str();

    const size_t dot = out.find('.');
    std::string head = (dot == std::string::npos) ? out : out.substr(0, dot);
    const std::string tail = (dot == std::string::npos) ? std::string() : out.substr(dot + 1);

    const bool negative = !head.empty() && head[0] == '-';
    if (negative) head.erase(0, 1);
    for (size_t pos = head.size(); pos > 3;)
    {
        pos -= 3;
        head.insert(pos, " ");   // неразрывный пробел
    }
    if (negative) head.insert(0, "-");

    return tail.empty() ? head : head + "," + tail;
}

// Согласование существительного с числом: 1 частица, 2 частицы, 5 частиц.
// Нужно потому, что ntotal и число шагов меняются от прогона к прогону, и
// зашитая форма при другом delta_x станет неверной.
inline const char* ru_plural(size_t n, const char* one, const char* few, const char* many)
{
    const size_t mod100 = n % 100;
    if (mod100 >= 11 && mod100 <= 14) return many;
    switch (n % 10)
    {
        case 1:          return one;
        case 2: case 3:
        case 4:          return few;
        default:         return many;
    }
}

} // namespace e2e
