#pragma once
// regression_snapshot.hpp
// Снимок состояния частиц и его сверка с эталоном — основа регрессионного
// теста tests/test_regression.cpp.
//
// ЧЕМ ЭТО ОТЛИЧАЕТСЯ ОТ E2E. Сквозные тесты BC2/BC6 сверяют расчёт с ПРИБОРНОЙ
// кривой по форме и с грубыми воротами: они отвечают на вопрос «похоже ли на
// физику». Регрессия отвечает на другой вопрос — «не изменился ли расчёт» — и
// сверяет ВСЕ поля КАЖДОЙ частицы с эталоном, снятым этим же кодом.
//
// ПОЧЕМУ СНИМОК БЕРЁТСЯ ИЗ ПАМЯТИ. Кадры sph_output_*.txt и log.txt конвейер
// пишет с точностью потока по умолчанию — 6 значащих цифр. Для регрессии этого
// мало: допуск по полям лежит на уровне 1e-9..1e-6 относительных.
//
// ПОЧЕМУ СВЕРКА С ДОПУСКОМ, А НЕ ПОБИТОВАЯ. При числе потоков OpenMP больше
// одного результат плавает на уровне округления: порядок суммирования
// вкладов в частицу (rhs::reduce_contributions, шарды owner-computes,
// fiber/PLAN.md §3.5) меняется от прогона к прогону. Побитовое совпадение к
// тому же сломалось бы на первой же правке, меняющей порядок суммирования.
//
// НОРМА. Для каждого поля f по частицам [0, ntotal):
//     scale_f = max(max_i |f_ref_i|, floor_f)
//     Linf_f  = max_i |f_i − f_ref_i| / scale_f
//     L2_f    = ‖f − f_ref‖₂ / max(‖f_ref‖₂, floor_f·√n)
// floor_f — физический пол поля: без него почти нулевые поля (D, e_pl, sxy на
// первом шаге) раздували бы относительную ошибку до бесконечности.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "param.hpp"
#include "particle_system.hpp"

static_assert(dim == 2 || dim == 3,
              "regression_snapshot.hpp рассчитан на dim = 2 (BC2/BC6) и dim = 3 (уменьшенный BC2)");

namespace reg
{

// Версия формата файла эталона. Меняется при любом изменении состава колонок
// ПРИ ДАННОЙ РАЗМЕРНОСТИ; сама размерность пишется в шапку отдельно (ключ dim),
// и эталон другой размерности не читается.
inline constexpr int kFormatVersion = 1;

// Поля снимка. Порядок совпадает с порядком колонок в файле эталона:
//   x0..x{dim-1}, v0..v{dim-1}, rho p u c hsml e_pl D, sxx syy szz sxy
//   и при dim = 3 ещё syz sxz.
// mass в эталон не входит: в расчёте она не меняется. При dim = 2 не входят и
// syz, sxz — там они тождественно нулевые; это проверяет сам тест.
inline constexpr size_t kX0 = 0, kX1 = 1;
inline constexpr size_t kV0 = dim, kV1 = dim + 1;
inline constexpr size_t kRho = 2 * dim, kP = kRho + 1, kU = kRho + 2, kC = kRho + 3;
inline constexpr size_t kHsml = kRho + 4, kEpl = kRho + 5, kD = kRho + 6;
inline constexpr size_t kSxx = kRho + 7, kSyy = kSxx + 1, kSzz = kSxx + 2, kSxy = kSxx + 3;
inline constexpr size_t kSyz = kSxx + 4, kSxz = kSxx + 5;   // только dim = 3
inline constexpr size_t kFieldCount = (dim == 3) ? kSxz + 1 : kSxy + 1;

namespace detail
{
inline constexpr std::array<const char*, 19> kName3 = {
    "x0", "x1", "x2", "v0", "v1", "v2", "rho", "p", "u", "c", "hsml", "e_pl", "D",
    "sxx", "syy", "szz", "sxy", "syz", "sxz"};
inline constexpr std::array<const char*, 15> kName2 = {
    "x0", "x1", "v0", "v1", "rho", "p", "u", "c", "hsml", "e_pl", "D",
    "sxx", "syy", "szz", "sxy"};
} // namespace detail

// Копированием, а не выбором инициализатора: список на 19 имён в массив на 15
// не компилируется даже в отброшенной ветке (лямбда здесь не шаблон).
inline constexpr auto kFieldName = [] {
    std::array<const char*, kFieldCount> a{};
    for (size_t k = 0; k < kFieldCount; ++k)
        a[k] = (dim == 3) ? detail::kName3[k] : detail::kName2[k];
    return a;
}();

// Физический пол поля (см. шапку). Единицы СИ.
//   x    — 1 мм: координаты тел порядка сантиметров, пол не срабатывает;
//   v    — 1 м/с; rho — 1 кг/м³; p и девиатор — 1 МПа; u — 1 Дж/кг;
//   c    — 1 м/с; hsml — 1 мкм;
//   e_pl — 0.01, D — 1: безразмерные, пол задаёт абсолютную чувствительность
//          там, где поле почти везде ноль (повреждение до прихода разгрузки).
inline constexpr auto kFieldFloor = [] {
    std::array<double, kFieldCount> f{};
    for (size_t d = 0; d < dim; ++d)
    {
        f[kX0 + d] = 1e-3;
        f[kV0 + d] = 1.0;
    }
    f[kRho] = 1.0; f[kP] = 1e6; f[kU] = 1.0; f[kC] = 1.0;
    f[kHsml] = 1e-6; f[kEpl] = 1e-2; f[kD] = 1.0;
    for (size_t k = kSxx; k < kFieldCount; ++k) f[k] = 1e6;
    return f;
}();

struct Snapshot
{
    // --- метаданные ---
    std::string case_name;          // "bc2" | "bc6" | "bc2mini3d"
    int    contact_force_type = 0;
    double delta_x = 0.0;
    size_t steps = 0;               // сколько шагов сделано до снимка
    double time = 0.0;              // модельное время после steps шагов
    double dt = 0.0;                // последний принятый шаг
    int    threads = 0;             // omp_get_max_threads() при снятии
    std::map<std::string, std::string> extra;   // прочие строки шапки

    // --- частицы [0, ntotal) ---
    std::vector<int> itype;
    std::array<std::vector<double>, kFieldCount> f;

    size_t size() const { return itype.size(); }
};

inline Snapshot capture(const ParticleSystem& sys, size_t ntotal)
{
    Snapshot s;
    s.itype.assign(sys.itype.begin(), sys.itype.begin() + static_cast<std::ptrdiff_t>(ntotal));
    // Поля конвейера — real_t (include/precision.hpp); снимок всегда double:
    // float -> double без потерь, формат эталона от точности не зависит.
    auto take = [&](size_t k, const std::vector<real_t>& src) {
        s.f[k].assign(src.begin(), src.begin() + static_cast<std::ptrdiff_t>(ntotal));
    };
    for (size_t d = 0; d < dim; ++d)
    {
        take(kX0 + d, sys.x[d]);
        take(kV0 + d, sys.vx[d]);
    }
    take(kRho, sys.rho);   take(kP, sys.p);     take(kU, sys.u);   take(kC, sys.c);
    take(kHsml, sys.hsml); take(kEpl, sys.e_pl); take(kD, sys.D);
    take(kSxx, sys.sxx);   take(kSyy, sys.syy); take(kSzz, sys.szz); take(kSxy, sys.sxy);
    if constexpr (dim == 3)
    {
        take(kSyz, sys.syz);
        take(kSxz, sys.sxz);
    }
    return s;
}

// ===========================================================================
// Запись и чтение
// ===========================================================================
// Шапка — строки «# ключ: значение»; данные — одна строка на частицу:
// «i itype» и kFieldCount полей в %.12e. 17 цифр не нужны: допуск по полям не
// бывает тоньше 1e-11, а ошибка представления %.12e — около 5e-13.

// printf-форматирование в строку: std::ofstream сам %.12e не умеет, а
// fopen на MSVC-рантайме помечен устаревшим.
template <typename... Args>
inline void appendf(std::string& out, const char* fmt, Args... args)
{
    char buf[256];
    const int n = std::snprintf(buf, sizeof(buf), fmt, args...);
    if (n > 0) out.append(buf, std::min(static_cast<size_t>(n), sizeof(buf) - 1));
}

inline bool write_snapshot(const Snapshot& s, const std::string& path, std::string& err)
{
    std::string t;
    t += "# Регрессионный эталон SPH (tests/test_regression.cpp)\n";
    t += "#\n";
    t += "# СГЕНЕРИРОВАН ТЕКУЩИМ КОДОМ, а не прибором: фиксирует поведение\n";
    t += "# расчёта, а не его физическую верность. Перегенерировать только при\n";
    t += "# осознанном изменении физики — с причиной в сообщении коммита\n";
    t += "# (docs/02-how-to/01-run-tests.md, раздел «Регрессия»).\n";
    t += "#\n";
    appendf(t, "# format: %d\n", kFormatVersion);
    appendf(t, "# dim: %zu\n", dim);
    t += "# case: " + s.case_name + "\n";
    appendf(t, "# contact_force_type: %d\n", s.contact_force_type);
    appendf(t, "# delta_x: %.17g\n", s.delta_x);
    appendf(t, "# ntotal: %zu\n", s.size());
    appendf(t, "# steps: %zu\n", s.steps);
    appendf(t, "# time: %.17g\n", s.time);
    appendf(t, "# dt: %.17g\n", s.dt);
    appendf(t, "# threads: %d\n", s.threads);
    for (const auto& [k, v] : s.extra) t += "# " + k + ": " + v + "\n";
    t += "#\n# columns: i itype";
    for (const char* n : kFieldName) t += std::string(" ") + n;
    t += "\n";
    for (size_t i = 0; i < s.size(); ++i)
    {
        appendf(t, "%zu %d", i, s.itype[i]);
        for (size_t k = 0; k < kFieldCount; ++k) appendf(t, " %.12e", s.f[k][i]);
        t += "\n";
    }

    // Двоичный режим: переводы строк — LF на любой платформе, чтобы эталон,
    // перегенерированный на другой ОС, не менялся целиком.
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out)
    {
        err = "не удалось открыть на запись: " + path;
        return false;
    }
    out.write(t.data(), static_cast<std::streamsize>(t.size()));
    out.close();
    if (!out)
    {
        err = "ошибка записи: " + path;
        return false;
    }
    return true;
}

inline bool read_snapshot(const std::string& path, Snapshot& s, std::string& err)
{
    std::ifstream in(path);
    if (!in)
    {
        err = "эталон не найден: " + path;
        return false;
    }
    s = Snapshot{};
    std::map<std::string, std::string> head;
    std::string line;
    size_t lineno = 0;
    while (std::getline(in, line))
    {
        ++lineno;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#')
        {
            const size_t colon = line.find(": ");
            if (colon == std::string::npos || colon < 2) continue;
            std::string key = line.substr(1, colon - 1);
            key.erase(0, key.find_first_not_of(' '));
            head[key] = line.substr(colon + 2);
            continue;
        }
        std::istringstream is(line);
        size_t idx = 0;
        int it = 0;
        if (!(is >> idx >> it) || idx != s.itype.size())
        {
            err = path + ":" + std::to_string(lineno) + ": испорчена строка частицы";
            return false;
        }
        s.itype.push_back(it);
        for (size_t k = 0; k < kFieldCount; ++k)
        {
            double v = 0.0;
            if (!(is >> v))
            {
                err = path + ":" + std::to_string(lineno) + ": не хватает поля " + kFieldName[k];
                return false;
            }
            s.f[k].push_back(v);
        }
    }

    auto need = [&](const char* key) -> const std::string* {
        const auto p = head.find(key);
        if (p == head.end())
        {
            err = path + ": в шапке нет ключа «" + key + "»";
            return nullptr;
        }
        return &p->second;
    };
    const std::string* v = nullptr;
    try
    {
        if (!(v = need("format"))) return false;
        if (std::stoi(*v) != kFormatVersion)
        {
            err = path + ": формат " + *v + ", ожидался " + std::to_string(kFormatVersion) +
                  " — эталон нужно перегенерировать";
            return false;
        }
        // Ключ dim появился вместе с 3D-эталонами; файл без него — 2D.
        const auto dim_it = head.find("dim");
        const size_t file_dim =
            dim_it == head.end() ? 2u : static_cast<size_t>(std::stoul(dim_it->second));
        if (file_dim != dim)
        {
            err = path + ": эталон для dim = " + std::to_string(file_dim) +
                  ", а тест собран с dim = " + std::to_string(dim);
            return false;
        }
        if (!(v = need("case"))) return false;
        s.case_name = *v;
        if (!(v = need("contact_force_type"))) return false;
        s.contact_force_type = std::stoi(*v);
        if (!(v = need("delta_x"))) return false;
        s.delta_x = std::stod(*v);
        if (!(v = need("steps"))) return false;
        s.steps = static_cast<size_t>(std::stoull(*v));
        if (!(v = need("time"))) return false;
        s.time = std::stod(*v);
        if (!(v = need("dt"))) return false;
        s.dt = std::stod(*v);
        if (!(v = need("threads"))) return false;
        s.threads = std::stoi(*v);
        if (!(v = need("ntotal"))) return false;
        if (static_cast<size_t>(std::stoull(*v)) != s.size())
        {
            err = path + ": в шапке ntotal = " + *v + ", строк частиц " + std::to_string(s.size());
            return false;
        }
    }
    catch (const std::exception&)
    {
        err = path + ": не разобрано значение в шапке: " + (v ? *v : std::string("?"));
        return false;
    }
    for (const char* key : {"format", "dim", "case", "contact_force_type", "delta_x", "ntotal",
                            "steps", "time", "dt", "threads", "columns"})
    {
        head.erase(key);
    }
    s.extra = std::move(head);
    return true;
}

// ===========================================================================
// Сверка
// ===========================================================================

struct Tolerance
{
    double inf = 0.0;   // предел Linf
    double l2  = 0.0;   // предел L2
};

struct Tolerances
{
    std::array<Tolerance, kFieldCount> field{};
    double time_rel = 0.0;   // |t − t_ref| / t_ref
    double dt_rel   = 0.0;   // |dt − dt_ref| / dt_ref
};

struct FieldResult
{
    double linf = 0.0;
    double l2   = 0.0;
    size_t worst = 0;        // частица с наибольшим |Δ|
    bool   ok = true;
};

struct Report
{
    bool comparable = true;          // размеры и itype совпали — поля сравнимы
    std::vector<std::string> fatal;  // почему несравнимо / жёсткие расхождения
    std::array<FieldResult, kFieldCount> field{};
    double time_rel = 0.0;
    double dt_rel   = 0.0;
    bool   time_ok = true;
    bool   dt_ok   = true;
    size_t itype_mismatch = 0;
    size_t first_itype_mismatch = 0;

    bool ok() const
    {
        if (!comparable || !fatal.empty() || !time_ok || !dt_ok) return false;
        for (const auto& r : field) if (!r.ok) return false;
        return true;
    }
};

inline double rel_diff(double a, double ref)
{
    const double d = std::fabs(a - ref);
    const double s = std::fabs(ref);
    return s > 0.0 ? d / s : d;
}

inline Report compare(const Snapshot& ref, const Snapshot& cur, const Tolerances& tol)
{
    Report r;
    if (ref.size() != cur.size())
    {
        r.comparable = false;
        r.fatal.push_back("ntotal: эталон " + std::to_string(ref.size()) + ", расчёт " +
                          std::to_string(cur.size()));
        return r;
    }
    if (ref.steps != cur.steps)
    {
        r.fatal.push_back("число шагов: эталон " + std::to_string(ref.steps) + ", расчёт " +
                          std::to_string(cur.steps));
    }
    if (ref.contact_force_type != cur.contact_force_type || ref.case_name != cur.case_name)
    {
        r.fatal.push_back("эталон снят с другого кейса: " + ref.case_name + "/contact" +
                          std::to_string(ref.contact_force_type));
    }
    for (size_t i = 0; i < ref.size(); ++i)
    {
        if (ref.itype[i] != cur.itype[i])
        {
            if (r.itype_mismatch == 0) r.first_itype_mismatch = i;
            ++r.itype_mismatch;
        }
    }
    if (r.itype_mismatch > 0)
    {
        r.comparable = false;
        r.fatal.push_back("itype разошёлся у " + std::to_string(r.itype_mismatch) +
                          " частиц, первая — " + std::to_string(r.first_itype_mismatch));
        return r;
    }

    r.time_rel = rel_diff(cur.time, ref.time);
    r.dt_rel   = rel_diff(cur.dt, ref.dt);
    r.time_ok  = r.time_rel <= tol.time_rel;
    r.dt_ok    = r.dt_rel <= tol.dt_rel;

    const size_t n = ref.size();
    for (size_t k = 0; k < kFieldCount; ++k)
    {
        const auto& a = cur.f[k];
        const auto& b = ref.f[k];
        double max_ref = 0.0, ref2 = 0.0, d2 = 0.0, dmax = -1.0;
        size_t worst = 0;
        bool finite = true;
        for (size_t i = 0; i < n; ++i)
        {
            const double d = std::fabs(a[i] - b[i]);
            if (!std::isfinite(a[i])) finite = false;
            max_ref = std::max(max_ref, std::fabs(b[i]));
            ref2 += b[i] * b[i];
            d2 += d * d;
            // !(d <= dmax) ловит и NaN: такая частица сразу становится худшей.
            if (!(d <= dmax)) { dmax = d; worst = i; }
        }
        const double floor = kFieldFloor[k];
        const double scale = std::max(max_ref, floor);
        const double norm2 = std::max(std::sqrt(ref2), floor * std::sqrt(static_cast<double>(n)));
        FieldResult& fr = r.field[k];
        fr.linf  = finite ? dmax / scale : INFINITY;
        fr.l2    = finite ? std::sqrt(d2) / norm2 : INFINITY;
        fr.worst = worst;
        fr.ok    = finite && fr.linf <= tol.field[k].inf && fr.l2 <= tol.field[k].l2;
    }
    return r;
}

// Отчёт печатается ВСЕГДА — и при успехе: по нему калибруются допуски
// (docs/02-how-to/01-run-tests.md, «Регрессия»), и по нему видно, насколько
// близко к порогу прошёл расчёт.
inline void print_report(const Report& r, const Snapshot& ref, const Snapshot& cur,
                         const Tolerances& tol)
{
    for (const auto& m : r.fatal) std::printf("  FAIL : %s\n", m.c_str());
    if (!r.comparable) return;

    std::printf("  info : t:  эталон %.9e, расчёт %.9e, отн. %.2e (допуск %.1e) %s\n",
                ref.time, cur.time, r.time_rel, tol.time_rel, r.time_ok ? "ok" : "FAIL");
    std::printf("  info : dt: эталон %.9e, расчёт %.9e, отн. %.2e (допуск %.1e) %s\n",
                ref.dt, cur.dt, r.dt_rel, tol.dt_rel, r.dt_ok ? "ok" : "FAIL");
    // Заголовок выровнен вручную: printf считает ширину в байтах, а
    // кириллица в UTF-8 занимает по два байта на букву.
    // Координаты худшей частицы: (r z) в осесимметричном 2D, (x y z) в 3D.
    std::printf("       поле        Linf   допуск         L2   допуск  "
                "%s           эталон           расчёт\n",
                dim == 3 ? "худшая частица (itype x y z)         "
                         : "худшая частица (itype r z)");
    for (size_t k = 0; k < kFieldCount; ++k)
    {
        const FieldResult& fr = r.field[k];
        const size_t i = fr.worst;
        char where[80];
        if constexpr (dim == 3)
            std::snprintf(where, sizeof(where), "%zu (%d %.2e %.2e %.2e)", i, ref.itype[i],
                          ref.f[kX0][i], ref.f[kX1][i], ref.f[kX0 + 2][i]);
        else
            std::snprintf(where, sizeof(where), "%zu (%d %.2e %.2e)", i, ref.itype[i],
                          ref.f[kX0][i], ref.f[kX1][i]);
        std::printf("  %-5s%-5s %10.2e %8.1e %10.2e %8.1e  %-*s %16.9e %16.9e\n",
                    fr.ok ? "ok" : "FAIL", kFieldName[k], fr.linf, tol.field[k].inf, fr.l2,
                    tol.field[k].l2, dim == 3 ? 37 : 26, where, ref.f[k][i], cur.f[k][i]);
    }
}

} // namespace reg
