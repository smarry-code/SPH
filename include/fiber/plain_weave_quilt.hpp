// plain_weave_quilt.hpp
// Ткань-«одеяло» из патчей поверх генератора полотняного переплетения
// (include/fiber/plain_weave.hpp) и сетка стёжки слоёв (этап 2б).
//
// Header-only. НЕ зависит от dim: только plain_weave.hpp и стандартная библиотека.
//
// ИДЕЯ. Прямоугольная сетка tile_x x tile_y одинаковых патчей по patch_ny x patch_nx
// нитей. Нити НЕПРЕРЫВНЫ через стыки, поэтому геометрия «одеяла» -- это одна
// обычная ткань с nx = tile_y*patch_nx нитей основы и ny = tile_x*patch_ny нитей
// утка; build_plain_weave вызывается без изменений, API и формат .fib прежние.
// Здесь добавлены: метка принадлежности бусин патчам, координаты стыков и запись
// раскладки по патчам. Фаза волны на стыках сходится сама, если patch_nx и
// patch_ny чётны (иначе патч (1,0) зеркален патчу (0,0), а не сдвинут).
//
// СТЫКИ. Стык между патчами t-1 и t лежит на середине полуволны между вершинами
// t*P-1 и t*P (P = patch_ny вдоль x, patch_nx вдоль y):
//   seam = c + (t*P - 1/2 - (n-1)/2)*p,  n = ny (вдоль x) или nx (вдоль y).
// При чётном seg на стыке есть бусина каждой поперечной нити (z = cz).
//
// ПРИНАДЛЕЖНОСТЬ БУСИНЫ. id патча = ty*tile_x + tx. Ряд вдоль оси нити определяется
// координатой бусины относительно стыков (бусина на стыке с допуском 1e-9*p
// относится к патчу с меньшим индексом; выступы -- крайним патчам), ряд поперёк
// нити -- номером нити: основа j -> ty = j / patch_nx, уток i -> tx = i / patch_ny.
//
// СТЁЖКА (этап 2б). Стежок = пара бусин одной вертикали в соседних слоях пакета на
// линии шва (стык патчей). Номера бусин -- ГЛОБАЛЬНЫЕ: сквозной счёт по слоям,
// нитям и бусинам в том же порядке, в каком .fib пишется слой за слоем.
#pragma once

#include "fiber/plain_weave.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fiber
{

struct QuiltParams
{
    PlainWeaveParams base;   // d, pitch, amp, seg, overhang_*, profile, c*, clamp, материал;
                             // поля base.nx / base.ny игнорируются и вычисляются
    int patch_nx = 34;       // нитей основы (вдоль x) в патче
    int patch_ny = 34;       // нитей утка  (вдоль y) в патче
    int tile_x = 3;          // патчей вдоль x (патчи по утку: ny = tile_x*patch_ny)
    int tile_y = 3;          // патчей вдоль y (nx = tile_y*patch_nx)
    bool allow_odd_patch = false; // разрешить нечётные patch_* при tile > 1 (с предупреждением)
};

struct PlainWeaveQuilt
{
    PlainWeave W;                              // обычная ткань: nx = tile_y*patch_nx, ny = tile_x*patch_ny
    QuiltParams params;
    std::vector<std::vector<int>> bead_patch;  // [нить][бусина] -> id патча = ty*tile_x + tx
    std::vector<double> seam_x, seam_y;        // стыки вдоль x и y (tile_x-1 и tile_y-1 значений)
};

// Координаты стыков по полной ткани full (используются full.pitch, full.cx, full.cy,
// full.nx, full.ny). Не требует чётных патчей: годится и для стёжки.
// sx -- линии x = const (tile_x-1 штук), sy -- линии y = const (tile_y-1 штук).
inline void quilt_seams(const PlainWeaveParams& full, int patch_nx, int patch_ny, int tile_x, int tile_y,
                        std::vector<double>& sx, std::vector<double>& sy)
{
    if (patch_nx < 1 || patch_ny < 1 || tile_x < 1 || tile_y < 1)
        throw std::invalid_argument("plain_weave_quilt: patch_nx, patch_ny, tile_x, tile_y должны быть >= 1");
    sx.clear();
    sy.clear();
    const double p = full.pitch;
    for (int t = 1; t < tile_x; ++t)
        sx.push_back(full.cx + (static_cast<double>(t) * patch_ny - 0.5 - 0.5 * (full.ny - 1)) * p);
    for (int t = 1; t < tile_y; ++t)
        sy.push_back(full.cy + (static_cast<double>(t) * patch_nx - 0.5 - 0.5 * (full.nx - 1)) * p);
}

namespace plain_weave_quilt_detail
{
// Номер ряда по координате вдоль нити: число стыков строго левее v (с допуском).
inline int row_of(double v, const std::vector<double>& seams, double tol)
{
    int r = 0;
    for (double s : seams)
        if (v > s + tol) ++r;
    return r;
}
} // namespace plain_weave_quilt_detail

inline PlainWeaveQuilt build_plain_weave_quilt(const QuiltParams& Q)
{
    if (Q.patch_nx < 1 || Q.patch_ny < 1 || Q.tile_x < 1 || Q.tile_y < 1)
        throw std::invalid_argument("plain_weave_quilt: patch_nx, patch_ny, tile_x, tile_y должны быть >= 1");
    const long long nxl = static_cast<long long>(Q.tile_y) * Q.patch_nx;
    const long long nyl = static_cast<long long>(Q.tile_x) * Q.patch_ny;
    if (nxl > 1000000 || nyl > 1000000)
        throw std::invalid_argument("plain_weave_quilt: число нитей tile*patch слишком велико");

    std::vector<std::string> odd_warn;
    if (Q.patch_ny % 2 != 0 && Q.tile_x > 1)
        odd_warn.push_back("WARN plain_weave_quilt: patch_ny нечётно при tile_x > 1: патчи по x зеркальны, а не сдвинуты");
    if (Q.patch_nx % 2 != 0 && Q.tile_y > 1)
        odd_warn.push_back("WARN plain_weave_quilt: patch_nx нечётно при tile_y > 1: патчи по y зеркальны, а не сдвинуты");
    if (!odd_warn.empty() && !Q.allow_odd_patch)
        throw std::invalid_argument(odd_warn.front().substr(5) + " (разрешить: --allow-odd-patch)");

    PlainWeaveQuilt R;
    R.params = Q;
    PlainWeaveParams q = Q.base;
    q.nx = static_cast<int>(nxl);
    q.ny = static_cast<int>(nyl);
    R.W = build_plain_weave(q);
    for (const auto& w : odd_warn)
    {
        R.W.warnings.push_back(w);
        std::fprintf(stderr, "%s\n", w.c_str());
    }
    quilt_seams(q, Q.patch_nx, Q.patch_ny, Q.tile_x, Q.tile_y, R.seam_x, R.seam_y);

    const double tol = 1e-9 * q.pitch;
    R.bead_patch.resize(R.W.fibers.size());
    for (std::size_t f = 0; f < R.W.fibers.size(); ++f)
    {
        const auto& F = R.W.fibers[f];
        auto& bp = R.bead_patch[f];
        bp.resize(F.pts.size());
        if (F.direction == 'x')
        {
            const int ty = F.index / Q.patch_nx;
            for (std::size_t b = 0; b < F.pts.size(); ++b)
                bp[b] = ty * Q.tile_x + plain_weave_quilt_detail::row_of(F.pts[b][0], R.seam_x, tol);
        }
        else
        {
            const int tx = F.index / Q.patch_ny;
            for (std::size_t b = 0; b < F.pts.size(); ++b)
                bp[b] = plain_weave_quilt_detail::row_of(F.pts[b][1], R.seam_y, tol) * Q.tile_x + tx;
        }
    }
    return R;
}

// <out>.patches.txt: строка на нить. Колонки: object_id direction index first_bead
// last_bead n_ranges, затем n_ranges троек "patch:first:last" (индексы бусин
// глобальные, как в .layout.txt; диапазоны идут вдоль нити подряд).
inline void write_plain_weave_quilt_patches(std::ostream& os, const PlainWeaveQuilt& R)
{
    const auto& Q = R.params;
    os << "# патчей " << Q.tile_x << "x" << Q.tile_y << " по " << Q.patch_nx << "x" << Q.patch_ny
       << " нитей (patch_nx x patch_ny); id патча = ty*tile_x + tx\n";
    os << "# object_id direction index first_bead last_bead n_ranges {patch:first:last}...\n";
    std::size_t first = 0;
    for (std::size_t f = 0; f < R.W.fibers.size(); ++f)
    {
        const auto& F = R.W.fibers[f];
        const auto& bp = R.bead_patch[f];
        const std::size_t last = first + F.pts.size() - 1;
        std::vector<std::array<std::size_t, 3>> rg; // patch, first, last (глобальные)
        for (std::size_t b = 0; b < bp.size(); ++b)
        {
            if (rg.empty() || static_cast<int>(rg.back()[0]) != bp[b])
                rg.push_back({static_cast<std::size_t>(bp[b]), first + b, first + b});
            else
                rg.back()[2] = first + b;
        }
        os << F.object_id << " " << F.direction << " " << F.index << " " << first << " " << last << " " << rg.size();
        for (const auto& r : rg) os << " " << r[0] << ":" << r[1] << ":" << r[2];
        os << "\n";
        first = last + 1;
    }
}

// ---------------------------------------------------------------------------
// Стёжка слоёв (этап 2б).
// ---------------------------------------------------------------------------
struct QuiltStitchParams
{
    int every = 1;        // прокол на каждой every-й поперечной нити (index % every == 0)
    bool lines_x = true;  // линии шва x = seam_x[t] (проколы на нитях основы)
    bool lines_y = true;  // линии шва y = seam_y[t] (проколы на нитях утка)
};

// Пары ГЛОБАЛЬНЫХ индексов бусин соседних слоёв (k, k+1) в точках прокола.
// Глобальный индекс: сквозной счёт по слоям, нитям и бусинам (порядок записи .fib
// слой за слоем). Слои обязаны совпадать по составу нитей и бусин и лежать
// строго друг над другом (в x-y). Нет бусины на стыке (нечётный seg) -- исключение.
inline std::vector<std::pair<std::size_t, std::size_t>>
build_quilt_stitches(const std::vector<PlainWeave>& plies, const std::vector<double>& seam_x,
                     const std::vector<double>& seam_y, const QuiltStitchParams& sp)
{
    std::vector<std::pair<std::size_t, std::size_t>> out;
    if (sp.every < 1) throw std::invalid_argument("plain_weave_quilt: stitch every должно быть >= 1");
    if (plies.size() < 2) return out;

    const PlainWeave& P0 = plies.front();
    const double p = P0.params.pitch;
    const double tol_seam = 1e-9 * p, tol_xy = 1e-12 * p;

    // Начальный глобальный индекс каждого слоя и проверка согласованности.
    std::vector<std::size_t> base(plies.size(), 0);
    std::size_t run = 0;
    for (std::size_t k = 0; k < plies.size(); ++k)
    {
        base[k] = run;
        const PlainWeave& P = plies[k];
        if (P.fibers.size() != P0.fibers.size())
            throw std::invalid_argument("plain_weave_quilt: слои пакета различаются числом нитей");
        for (std::size_t f = 0; f < P.fibers.size(); ++f)
        {
            const auto &F = P.fibers[f], &F0 = P0.fibers[f];
            if (F.pts.size() != F0.pts.size() || F.direction != F0.direction || F.index != F0.index)
                throw std::invalid_argument("plain_weave_quilt: слои пакета различаются составом бусин");
            for (std::size_t b = 0; b < F.pts.size(); ++b)
                if (std::fabs(F.pts[b][0] - F0.pts[b][0]) > tol_xy || std::fabs(F.pts[b][1] - F0.pts[b][1]) > tol_xy)
                    throw std::invalid_argument("plain_weave_quilt: слои пакета не лежат строго друг над другом");
            run += F.pts.size();
        }
    }

    // Начало нити f внутри слоя.
    std::vector<std::size_t> fstart(P0.fibers.size(), 0);
    {
        std::size_t s = 0;
        for (std::size_t f = 0; f < P0.fibers.size(); ++f)
        {
            fstart[f] = s;
            s += P0.fibers[f].pts.size();
        }
    }

    auto add_line = [&](char dir, double seam) {
        const int axis = dir == 'x' ? 0 : 1; // у основы бусина на x = seam, у утка на y = seam
        for (std::size_t f = 0; f < P0.fibers.size(); ++f)
        {
            const auto& F = P0.fibers[f];
            if (F.direction != dir || F.index % sp.every != 0) continue;
            std::size_t hit = F.pts.size();
            for (std::size_t b = 0; b < F.pts.size(); ++b)
                if (std::fabs(F.pts[b][axis] - seam) < tol_seam)
                {
                    hit = b;
                    break;
                }
            if (hit == F.pts.size())
            {
                std::ostringstream s;
                s << "plain_weave_quilt: на линии шва " << (dir == 'x' ? "x" : "y") << " = " << seam << " у нити "
                  << dir << F.index << " нет бусины (нечётный seg?)";
                throw std::invalid_argument(s.str());
            }
            for (std::size_t k = 0; k + 1 < plies.size(); ++k)
                out.emplace_back(base[k] + fstart[f] + hit, base[k + 1] + fstart[f] + hit);
        }
    };
    if (sp.lines_x)
        for (double s : seam_x) add_line('x', s);
    if (sp.lines_y)
        for (double s : seam_y) add_line('y', s);
    return out;
}

// Шов-строчка (в плоскости слоя): на линии шва каждая поперечная нить (основа на x = seam_x[t], уток на
// y = seam_y[t]) имеет одну бусину-прокол; проколы сортируются вдоль линии (по y у основы, по x у утка), и
// соседние проколы (индексы нитей j и j + every при every-фильтре index % every == 0) связываются стежком.
// Каждый слой прошивается независимо; число стежков на линию = (проколов - 1) на слой. Индексы -- глобальные
// (тот же порядок, что у build_quilt_stitches): слои по порядку, внутри слоя линии x, затем линии y.
inline std::vector<std::pair<std::size_t, std::size_t>>
build_quilt_seam_stitches(const std::vector<PlainWeave>& plies, const std::vector<double>& seam_x,
                          const std::vector<double>& seam_y, const QuiltStitchParams& sp)
{
    std::vector<std::pair<std::size_t, std::size_t>> out;
    if (sp.every < 1) throw std::invalid_argument("plain_weave_quilt: stitch every должно быть >= 1");
    if (plies.empty()) return out;

    const double p = plies.front().params.pitch;
    const double tol_seam = 1e-9 * p;

    std::size_t base = 0;
    for (const PlainWeave& P : plies)
    {
        // Начало каждой нити внутри слоя.
        std::vector<std::size_t> fstart(P.fibers.size(), 0);
        std::size_t s = 0;
        for (std::size_t f = 0; f < P.fibers.size(); ++f)
        {
            fstart[f] = s;
            s += P.fibers[f].pts.size();
        }
        auto add_line = [&](char dir, double seam) {
            const int axis = dir == 'x' ? 0 : 1;  // у основы бусина на x = seam, у утка на y = seam
            const int along = 1 - axis;           // координата вдоль линии шва
            std::vector<std::pair<double, std::size_t>> punct; // (координата вдоль линии, глобальный индекс бусины)
            for (std::size_t f = 0; f < P.fibers.size(); ++f)
            {
                const auto& F = P.fibers[f];
                if (F.direction != dir || F.index % sp.every != 0) continue;
                std::size_t hit = F.pts.size();
                for (std::size_t b = 0; b < F.pts.size(); ++b)
                    if (std::fabs(F.pts[b][axis] - seam) < tol_seam)
                    {
                        hit = b;
                        break;
                    }
                if (hit == F.pts.size())
                {
                    std::ostringstream es;
                    es << "plain_weave_quilt: на линии шва " << (dir == 'x' ? "x" : "y") << " = " << seam
                       << " у нити " << dir << F.index << " нет бусины (нечётный seg?)";
                    throw std::invalid_argument(es.str());
                }
                punct.emplace_back(F.pts[hit][along], base + fstart[f] + hit);
            }
            std::sort(punct.begin(), punct.end());
            for (std::size_t k = 0; k + 1 < punct.size(); ++k) out.emplace_back(punct[k].second, punct[k + 1].second);
        };
        if (sp.lines_x)
            for (double sx : seam_x) add_line('x', sx);
        if (sp.lines_y)
            for (double sy : seam_y) add_line('y', sy);
        base += s;
    }
    return out;
}

} // namespace fiber
