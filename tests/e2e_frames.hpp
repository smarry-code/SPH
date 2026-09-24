#pragma once
// e2e_frames.hpp
// Разбор кадров, которые конвейер пишет сам (output.cpp -> sph_output_<ts>.txt):
// проверка живучести (NaN/Inf, rho > 0, mass > 0), снятие осевой скорости
// произвольных частиц и аудит повреждения JH-2 по материалу.
//
// Происхождение проверки живучести — tests/test_e2e_bc2.cpp, раздел 4.
// Расширения (зонды по частицам, аудит D, привязка кадров ко времени) —
// новые: их требуют проверки §6 постановки BC6, которых у BC2 не было.
// Сам BC2 на этот заголовок не переключён, см. пояснение в e2e_harness.hpp.
//
// ЗАЧЕМ ЧИТАТЬ КАДРЫ, А НЕ ПРАВИТЬ КОНВЕЙЕР. Боевой time_integration пишет в
// log.txt РОВНО ОДИН монитор (cfg.moni_particle). Постановке BC6 нужны два
// монитора (на оси и на r = 12.5 мм, §2.4) и повреждение D по всему образцу
// (§6). И то, и другое уже лежит в кадрах, которые конвейер пишет и так —
// значит правка боевого кода не нужна вовсе.
//
// ЦЕНА. Кадр — это ntotal строк по 21 числу. Полный прогон BC6: 86 кадров на
// 44 638 частиц ≈ 0,9 ГБ текста во временном каталоге прогона. Каталог
// удаляется ScopedWorkDir, но место на диске под него нужно.
//
// ФОРМАТ СТРОКИ (output.cpp + ParticleSystem::vuvod_index, dim = 2):
//   0:N 1:itype 2:mass 3:r 4:z 5:v_r 6:v_z 7:rho 8:p 9:u 10:c
//   11:sigma_i 12:e_pl 13:D 14:hsml 15..20:sxx syy szz sxy syz sxz
// Индексы ниже выражены через dim, а не зашиты, и проверяются по числу полей.
//
// ПРИВЯЗКА КАДРА КО ВРЕМЕНИ. Разбирается по порядку операций в цикле
// time_integration (src/time_integration.cpp): счётчик шагов увеличивается в
// начале итерации, затем частицы обновляются на полный шаг и проходят модель
// прочности, затем пишется кадр (строка 300), затем время увеличивается
// (строка 303) и пишется строка log.txt с УЖЕ увеличенным временем (строка 319).
//
// Отсюда: кадр sph_output_<n> содержит состояние ПОСЛЕ n шагов, и это ровно то
// состояние, которое попало в n-ю строку лога. Значит кадр n отвечает строке
// лога n, то есть элементу n-1 при нумерации с нуля.
//
// Легко ошибиться на единицу в другую сторону: кадр помечен значением time ДО
// инкремента, то есть его собственная метка времени отстаёт от его же
// содержимого на один шаг. Это особенность боевого вывода; здесь берётся
// время, отвечающее СОСТОЯНИЮ, а не метка. Проверяется это прямо в тесте:
// осевой монитор, снятый с кадров, обязан совпасть с ним же из log.txt
// побитово, и при сдвиге на единицу совпадения не будет.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "param.hpp"
#include "e2e_curve.hpp"

namespace e2e
{

namespace fs = std::filesystem;

// Индексы колонок кадра, выраженные через dim.
inline constexpr size_t kIdxItype = 1;
inline constexpr size_t kIdxMass  = 2;
inline constexpr size_t kIdxX     = 3;                              // r
// Осевая скорость — последняя компонента (v_z при dim = 2 и dim = 3), как и
// в log.txt (time_integration: moni_axis = dim − 1).
inline constexpr size_t kIdxVAxis = 3 + dim + (dim - 1);
inline constexpr size_t kIdxRho   = 3 + 2 * dim;
inline constexpr size_t kIdxD     = kIdxRho + 6;                     // rho p u c sigma_i e_pl D
inline constexpr size_t kMinCols  = kIdxD + 1;

struct FrameCheck
{
    size_t frames = 0;
    size_t rows = 0;
    size_t bad_finite = 0;
    size_t bad_rho = 0;
    size_t bad_mass = 0;
    size_t short_rows = 0;      // строк с недостаточным числом полей
    std::string first_problem;
};

// Что снимать с кадров сверх проверки живучести.
struct FrameProbe
{
    // Сквозные индексы частиц (0-based), чью осевую скорость надо собрать.
    std::vector<size_t> particles;
    // Материал, по которому ведётся аудит повреждения D; -1 — не вести.
    int damage_itype = -1;
};

struct FrameScan
{
    FrameCheck health;
    std::vector<size_t> timestep;                 // номер из имени файла, по возрастанию
    std::vector<std::vector<double>> vz;          // [зонд][кадр]
    std::vector<double> max_damage_by_frame;      // по кадрам, для damage_itype
    double max_damage = 0.0;
    size_t max_damage_ts = 0;
    // Где стоит самая повреждённая частица. Без этого «D = 1» неразличимо
    // между настоящим отколом внутри образца и краевым артефактом на
    // свободной поверхности, где у SPH недобор соседей.
    double max_damage_r = 0.0;
    double max_damage_z = 0.0;
    size_t max_damage_index = 0;
    bool   probes_resolved = true;                // все ли зонды найдены в каждом кадре
};

// Номер шага из имени sph_output_<ts>.txt; npos, если имя не подходит.
inline bool frame_timestep(const std::string& name, size_t& out)
{
    if (name.rfind("sph_output_", 0) != 0) return false;
    const size_t dot = name.rfind('.');
    if (dot == std::string::npos) return false;
    const std::string digits = name.substr(11, dot - 11);
    if (digits.empty()) return false;
    for (char ch : digits) if (ch < '0' || ch > '9') return false;
    try { out = static_cast<size_t>(std::stoull(digits)); }
    catch (const std::exception&) { return false; }
    return true;
}

inline FrameScan scan_frames(const fs::path& dir, const FrameProbe& probe)
{
    FrameScan scan;
    std::vector<std::pair<size_t, fs::path>> files;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(dir, ec))
    {
        const std::string n = e.path().filename().string();
        size_t ts = 0;
        if (e.path().extension() == ".txt" && frame_timestep(n, ts))
        {
            files.emplace_back(ts, e.path());
        }
    }
    std::sort(files.begin(), files.end());
    scan.health.frames = files.size();
    scan.vz.assign(probe.particles.size(), {});

    std::vector<double> vals;
    for (const auto& [ts, path] : files)
    {
        scan.timestep.push_back(ts);
        for (auto& series : scan.vz) series.push_back(std::nan(""));
        double frame_damage = 0.0;
        double frame_damage_r = 0.0, frame_damage_z = 0.0;
        size_t frame_damage_index = 0;
        bool   resolved_here = probe.particles.empty();
        size_t resolved_count = 0;

        std::ifstream in(path);
        std::string line;
        size_t row = 0;   // порядковый номер частицы в кадре = её сквозной индекс
        while (std::getline(in, line))
        {
            if (line.empty()) continue;
            // Строка заголовка начинается с 'N'; данные — с цифры.
            const char c0 = line[0];
            if (!(c0 >= '0' && c0 <= '9')) continue;
            vals.clear();
            std::istringstream is(line);
            double v;
            while (is >> v) vals.push_back(v);
            if (vals.size() < kMinCols)
            {
                ++scan.health.short_rows;
                continue;
            }
            ++scan.health.rows;

            for (double x : vals)
            {
                if (!std::isfinite(x))
                {
                    ++scan.health.bad_finite;
                    if (scan.health.first_problem.empty())
                    {
                        scan.health.first_problem =
                            "не-конечное значение в " + path.filename().string();
                    }
                    break;
                }
            }
            if (!(vals[kIdxMass] > 0.0))
            {
                ++scan.health.bad_mass;
                if (scan.health.first_problem.empty())
                {
                    scan.health.first_problem = "mass <= 0 в " + path.filename().string();
                }
            }
            if (!(vals[kIdxRho] > 0.0))
            {
                ++scan.health.bad_rho;
                if (scan.health.first_problem.empty())
                {
                    scan.health.first_problem = "rho <= 0 в " + path.filename().string();
                }
            }

            if (probe.damage_itype >= 0 &&
                static_cast<int>(vals[kIdxItype]) == probe.damage_itype)
            {
                const double d = vals[kIdxD];
                if (std::isfinite(d) && d > frame_damage)
                {
                    frame_damage = d;
                    // r — расстояние от оси (в 3D ось — z), z — осевая координата.
                    frame_damage_r = (dim == 3) ? std::hypot(vals[kIdxX], vals[kIdxX + 1])
                                                : vals[kIdxX];
                    frame_damage_z = vals[kIdxX + dim - 1];
                    frame_damage_index = row;
                }
            }

            for (size_t k = 0; k < probe.particles.size(); ++k)
            {
                if (probe.particles[k] == row)
                {
                    scan.vz[k].back() = vals[kIdxVAxis];
                    ++resolved_count;
                }
            }
            ++row;
        }

        if (!resolved_here && resolved_count < probe.particles.size())
        {
            scan.probes_resolved = false;
        }
        scan.max_damage_by_frame.push_back(frame_damage);
        if (frame_damage > scan.max_damage)
        {
            scan.max_damage = frame_damage;
            scan.max_damage_ts = ts;
            scan.max_damage_r = frame_damage_r;
            scan.max_damage_z = frame_damage_z;
            scan.max_damage_index = frame_damage_index;
        }
    }
    return scan;
}

// Совместимость с прежним вызовом BC2: только проверка живучести.
inline FrameCheck check_frames(const fs::path& dir)
{
    return scan_frames(dir, FrameProbe{}).health;
}

// Время состояния, записанного в кадре с номером шага ts: это время n-й строки
// лога (нумерация строк с единицы). Возвращает false, если номер выходит за
// пределы лога — например, лог оборван вместе с прогоном.
inline bool frame_time(const Curve& log_curve, size_t ts, double& t_out)
{
    if (ts == 0) return false;
    const size_t idx = ts - 1;
    if (idx >= log_curve.t.size()) return false;
    t_out = log_curve.t[idx];
    return true;
}

// Кривая по одному зонду: время берётся из лога, значение — из кадров.
// Точки, для которых время или значение не определены, пропускаются.
inline Curve probe_curve(const FrameScan& scan, size_t probe, const Curve& log_curve,
                         int itype_hint)
{
    Curve c;
    if (probe >= scan.vz.size()) return c;
    for (size_t k = 0; k < scan.timestep.size(); ++k)
    {
        double t = 0.0;
        if (!frame_time(log_curve, scan.timestep[k], t)) continue;
        const double v = scan.vz[probe][k];
        if (!std::isfinite(v)) continue;
        c.t.push_back(t);
        c.v.push_back(v);
        c.itype.push_back(itype_hint);
    }
    return c;
}

} // namespace e2e
