// k_geometry.cpp
// Сборка Solid из LS-DYNA .k (этап K4, readk/k-import-plan.md, Р3-Р6, Р9).
#include "k_input/k_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

#include "int_force.hpp"
#include "param.hpp"

namespace k_input
{

namespace
{

[[noreturn]] void fail(const std::string& where, const std::string& msg)
{
    throw std::runtime_error(where + ": " + msg);
}

std::string join_pids(const std::set<long long>& pids)
{
    std::string l;
    for (long long p : pids) l += (l.empty() ? "" : ",") + std::to_string(p);
    return l;
}

} // namespace

void build_solid_from_deck(const KDeck& deck, const KPartsTable* table, const KImportOptions& opts,
                           const Materials& mats, ParticleSystem& sys, std::size_t& moni_particle,
                           KImportSummary& summary)
{
    const std::size_t n = deck.sph.size();
    summary.nodes = deck.nodes.size();
    summary.sph_particles = n;
    summary.skipped = deck.skipped;

    // Роли, не входящие в режим, пропускаем с предупреждением.
    if (opts.mode == KMode::Solid && !deck.beams.empty())
    {
        std::set<long long> pids;
        for (const auto& b : deck.beams) pids.insert(b.pid);
        summary.warnings.push_back("*ELEMENT_BEAM пропущены (режим solid): " + std::to_string(deck.beams.size())
                                   + " элементов, PID: " + join_pids(pids));
    }
    if (!deck.solids.empty())
    {
        std::set<long long> pids;
        for (const auto& s : deck.solids) pids.insert(s.pid);
        summary.warnings.push_back("*ELEMENT_SOLID пропущены (пока не поддерживаются): "
                                   + std::to_string(deck.solids.size()) + " элементов, PID: " + join_pids(pids));
    }

    // Координаты в СИ (Р5) и габарит по SPH-узлам.
    std::vector<std::array<double, 3>> pos(n);
    std::array<double, 3> lo{}, hi{};
    for (std::size_t i = 0; i < n; ++i)
    {
        const KSph& e = deck.sph[i];
        const KNode* nd = deck.find_node(e.nid);
        if (!nd) fail(deck.where(e.src), "узел " + std::to_string(e.nid) + " не найден");
        pos[i] = {nd->x * opts.length_scale, nd->y * opts.length_scale, nd->z * opts.length_scale};
        for (int c = 0; c < 3; ++c)
        {
            if (i == 0) { lo[c] = hi[c] = pos[i][c]; }
            else { lo[c] = std::min(lo[c], pos[i][c]); hi[c] = std::max(hi[c], pos[i][c]); }
        }
    }
    summary.bbox_min = lo;
    summary.bbox_max = hi;

    // Р6: лишние координаты должны быть нулевыми (допуск относительный от габарита).
    if (n > 0 && dim < 3)
    {
        double ext = 0.0;
        for (int c = 0; c < 3; ++c) ext = std::max(ext, hi[c] - lo[c]);
        const double tol = 1e-9 * std::max(ext, 1e-300);
        for (std::size_t i = 0; i < n; ++i)
            for (int c = static_cast<int>(dim); c < 3; ++c)
                if (std::fabs(pos[i][c]) > tol)
                    fail(deck.where(deck.sph[i].src),
                         "узел NID " + std::to_string(deck.sph[i].nid) + ": |" + "xyz"[c] + "|="
                             + std::to_string(std::fabs(pos[i][c])) + " > допуска "
                             + std::to_string(tol) + " при dim=" + std::to_string(dim));
    }

    sys.resize(n);
    std::vector<double> dxs;
    for (std::size_t i = 0; i < n; ++i)
    {
        const KSph& e = deck.sph[i];
        sys.inzializ_index(i);
        for (std::size_t d = 0; d < dim; ++d) sys.x[d][i] = pos[i][d];
        if (e.mass == 0.0) fail(deck.where(e.src), "MASS = 0 у узла " + std::to_string(e.nid));

        if (!table)
        {
            // Режим без таблицы: как старый parse_k_file.
            if (e.mass < 0.0 || e.is_volume)
                fail(deck.where(e.src), "MASS-объём (<0 или _VOLUME) требует таблицы деталей (нужна плотность)");
            sys.mass[i] = e.mass * opts.mass_scale;
            sys.itype[i] = static_cast<int>(e.pid);
            sys.body[i] = static_cast<int>(e.pid);
            continue;
        }

        const SolidPart* sp = table->find_solid(e.pid);
        if (!sp)
            fail(table->file, "PID " + std::to_string(e.pid) + " из " + deck.where(e.src)
                                  + " не описан в SOLID_PART");
        sys.itype[i] = sp->itype;
        sys.body[i] = static_cast<int>(e.pid);
        sys.rho[i] = sp->rho;
        sys.rho_0[i] = sp->rho;
        // MASS > 0 - масса, масштабируется mass_scale. MASS < 0 или _VOLUME -
        // "объём": масса = |MASS|*length_scale^dim*rho. При dim < 3 это объём в
        // смысле конвенции проекта (площадь/длина), а не геометрический объём.
        if (e.mass > 0.0 && !e.is_volume)
            sys.mass[i] = e.mass * opts.mass_scale;
        else
            sys.mass[i] = std::fabs(e.mass) * std::pow(opts.length_scale, static_cast<double>(dim)) * sp->rho;
        const double dx = std::pow(sys.mass[i] / sp->rho, 1.0 / static_cast<double>(dim));
        dxs.push_back(dx);
        sys.hsml[i] = sp->has_hsml ? sp->hsml : 2.0 * dx;
        for (int c = 0; c < sp->nv; ++c)
        {
            if (c < static_cast<int>(dim)) sys.vx[c][i] = sp->v[c];
            else if (sp->v[c] != 0.0)
                fail(table->file + ":" + std::to_string(sp->line),
                     "v= задаёт ненулевую компоненту " + std::to_string(c + 1) + " при dim="
                         + std::to_string(dim));
        }
    }

    if (table)
    {
        for (std::size_t i = 0; i < n; ++i) urs(mats, sys, i);
        if (!dxs.empty())
        {
            const std::size_t m = dxs.size() / 2;
            std::nth_element(dxs.begin(), dxs.begin() + static_cast<std::ptrdiff_t>(m), dxs.end());
            summary.dx_median = dxs[m];
        }
    }

    if (opts.monitor_nid >= 0)
    {
        bool found = false;
        for (std::size_t i = 0; i < n; ++i)
            if (deck.sph[i].nid == opts.monitor_nid) { moni_particle = i; found = true; break; }
        if (!found)
            fail(deck.files.empty() ? opts.k_file : deck.files[0],
                 "k_monitor_nid=" + std::to_string(opts.monitor_nid) + " не среди SPH-частиц");
    }
}

KImportSummary load_k_geometry(const KImportOptions& opts, const Materials& mats, fiber::Model& model,
                               std::size_t& moni_particle)
{
    if (opts.mode == KMode::SolidFiber)
    {
        if (opts.parts_file.empty())
            throw std::runtime_error("k_mode=solid_fiber требует таблицу деталей (k_parts_file)");
        throw std::runtime_error("k_mode=solid_fiber: сборка нитей из .k ещё не реализована (K5)");
    }
    const KDeck deck = read_k_deck(opts.k_file);
    KImportSummary summary;
    KPartsTable table;
    const KPartsTable* tp = nullptr;
    if (!opts.parts_file.empty())
    {
        table = read_k_parts_file(opts.parts_file);
        check_parts_against_deck(table, deck, false, summary.warnings);
        tp = &table;
    }
    build_solid_from_deck(deck, tp, opts, mats, model.solid, moni_particle, summary);
    return summary;
}

} // namespace k_input
