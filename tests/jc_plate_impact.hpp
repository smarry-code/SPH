#pragma once
// jc_plate_impact.hpp
// Геометрия и разбор кадров для сквозных тестов плоского удара E2 (HEL) и
// E3 (откол) модели Джонсона-Кука (fiber/jc-implementation-plan.md §5.3).
//
// ПОСТАНОВКА (SPH_DIM=1, одноосная деформация). Ударник занимает x in [-Lf, 0),
// мишень x in [0, Lt]; ось x направлена в мишень, тыльная (свободная)
// поверхность мишени — x = Lt. Частицы стоят в центрах ячеек: x = (k + 1/2)·dx,
// масса m = rho0·dx (длина ячейки dx, единичное сечение).
//
// ВЗАИМОДЕЙСТВИЕ УДАРНИКА С МИШЕНЬЮ. Материал один (STEEL 4340, itype 12), тела
// касаются решёткой без зазора (шаг через границу тоже dx). Поэтому они — один
// SPH-континуум: передача импульса идёт обычным градиентом давления, а разрыв
// начальной скорости на плоскости удара гасит псевдовязкость. Контакт
// (contact_force_type != 0) не нужен и выключен; sys.body у всех частиц
// одинаков (= itype), что и для контакта Solid-Solid давало бы «одно тело».
// Начальная скорость задаётся поштучно, отдельные тела для этого не нужны.
// Зазор между телами не делается — он вносил бы нефизичный разгон (см.
// bc6_case.hpp, «СБОРКА ВСТЫК»).
//
// ЗНАКИ. Конвейер: сжатие даёт p > 0, девиатор s = sxx. Продольное напряжение
// (растяжение положительно) sigma_x = -p + sxx; сжимающее напряжение
// (сжатие положительно) sigma_c = p - sxx. Везде ниже используется sigma_c.
//
// КАДРЫ (dim = 1): колонки строки 0:N 1:itype 2:mass 3:x 4:v 5:rho 6:p 7:u
// 8:c 9:sigma_i 10:e_pl 11:D 12:hsml 13:sxx 14:syy 15:szz 16:sxy 17:syz
// 18:sxz 19:T 20:edot. Время кадра берётся по e2e::frame_time из log.txt.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "param.hpp"
#include "particle_system.hpp"
#include "materials.hpp"

#include "e2e_curve.hpp"
#include "e2e_frames.hpp"

static_assert(dim == 1, "jc_plate_impact.hpp требует сборки с -DSPH_DIM=1");

namespace jcplate
{

namespace fs = std::filesystem;

inline constexpr int kItype = 12;   // STEEL 4340

struct Case
{
    const char* name = "hel";
    double v_impact = 400.0;    // м/с
    double L_flyer  = 5.0e-3;   // м
    double L_target = 20.0e-3;  // м
    double dx       = 1.25e-5;  // м
    double t_total  = 4.0e-6;   // с
    double save_t   = 1.0e-8;   // с
    bool   damage   = true;
    // Модель границы ударник-мишень: false = единый континуум (сварены), true = разные тела
    // (sph_same_body_only + контакт Solid-Solid, передаёт только сжатие).
    bool   contact  = false;
    int    contact_type = 3;        // 1 пинбол, 2 Свигл, 3 импульсный
    double contact_gamma = 0.5;     // доля гашения (тип 3), см. bc6_initial_conditions.txt
    std::vector<double> gauge_x;   // начальные координаты манометров, м
};

struct Geometry
{
    size_t nf = 0, nt = 0, ntotal = 0;
    double dx = 0.0;
    std::vector<double> x0;               // начальные координаты (лагранжевы)
    std::vector<size_t> gauge;            // индексы частиц-манометров
    std::vector<double> gauge_x0;         // их фактические координаты
    size_t index_of(double x) const
    {
        size_t best = 0;
        for (size_t i = 0; i < x0.size(); ++i)
            if (std::fabs(x0[i] - x) < std::fabs(x0[best] - x)) best = i;
        return best;
    }
};

inline Geometry build_geometry(ParticleSystem& sys, const Case& c)
{
    Geometry g;
    g.dx = c.dx;
    g.nf = static_cast<size_t>(std::llround(c.L_flyer / c.dx));
    g.nt = static_cast<size_t>(std::llround(c.L_target / c.dx));
    g.ntotal = g.nf + g.nt;
    sys.resize(g.ntotal);
    g.x0.resize(g.ntotal);

    const double rho0 = 7830.0;   // материал itype 12; проверяется в тесте
    for (size_t k = 0; k < g.ntotal; ++k)
    {
        const bool flyer = k < g.nf;
        const double x = flyer ? -(static_cast<double>(g.nf - k) - 0.5) * c.dx
                               :  (static_cast<double>(k - g.nf) + 0.5) * c.dx;
        g.x0[k] = x;
        sys.inzializ_index(k);
        sys.x[0][k]  = x;
        sys.vx[0][k] = flyer ? c.v_impact : 0.0;
        sys.hsml[k]  = 2.0 * c.dx;       // на 1-м шаге h_upgrade (sle=1) → 1.3·dx
        sys.itype[k] = kItype;
        sys.body[k]  = (c.contact && flyer) ? kItype + 1 : kItype;
        sys.rho[k]   = rho0;
        sys.rho_0[k] = rho0;
        sys.mass[k]  = rho0 * c.dx;
    }
    for (double xg : c.gauge_x)
    {
        const size_t i = g.index_of(xg);
        g.gauge.push_back(i);
        g.gauge_x0.push_back(g.x0[i]);
    }
    return g;
}

// ---------------------------------------------------------------------------
// Разбор кадров: полные колонки по всем частицам.
// ---------------------------------------------------------------------------
struct Frame
{
    size_t timestep = 0;
    double t = 0.0;               // время СОСТОЯНИЯ кадра, с
    std::vector<double> x, v, p, sxx, D, edot, epl;
};

inline bool read_frame(const fs::path& path, Frame& f)
{
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    std::vector<double> vals;
    while (std::getline(in, line))
    {
        if (line.empty() || !(line[0] >= '0' && line[0] <= '9')) continue;
        vals.clear();
        std::istringstream is(line);
        double d;
        while (is >> d) vals.push_back(d);
        if (vals.size() < 21) return false;
        f.x.push_back(vals[3]);
        f.v.push_back(vals[4]);
        f.p.push_back(vals[6]);
        f.epl.push_back(vals[10]);
        f.D.push_back(vals[11]);
        f.sxx.push_back(vals[13]);
        f.edot.push_back(vals[20]);
    }
    return true;
}

// Все кадры каталога, привязанные ко времени по log.txt.
inline std::vector<Frame> read_frames(const fs::path& dir, const e2e::Curve& log)
{
    std::vector<std::pair<size_t, fs::path>> files;
    for (const auto& e : fs::directory_iterator(dir))
    {
        size_t ts = 0;
        if (e.path().extension() == ".txt" &&
            e2e::frame_timestep(e.path().filename().string(), ts))
            files.emplace_back(ts, e.path());
    }
    std::sort(files.begin(), files.end());
    std::vector<Frame> out;
    for (const auto& [ts, p] : files)
    {
        Frame f;
        f.timestep = ts;
        if (!e2e::frame_time(log, ts, f.t)) continue;
        if (!read_frame(p, f)) continue;
        out.push_back(std::move(f));
    }
    return out;
}

// Сжимающее напряжение частицы i кадра f, Па.
inline double sigma_c(const Frame& f, size_t i) { return f.p[i] - f.sxx[i]; }

// Время первого пересечения порога снизу вверх рядом (линейная интерполяция
// между кадрами); NaN, если не пересекло.
template <class Getter>
inline double first_crossing(const std::vector<Frame>& fr, Getter get, double thr)
{
    for (size_t k = 1; k < fr.size(); ++k)
    {
        const double a = get(fr[k - 1]), b = get(fr[k]);
        if (a < thr && b >= thr)
            return fr[k - 1].t + (thr - a) / (b - a) * (fr[k].t - fr[k - 1].t);
    }
    return std::nan("");
}

} // namespace jcplate
