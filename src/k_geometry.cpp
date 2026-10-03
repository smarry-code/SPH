// k_geometry.cpp
// Сборка Solid из LS-DYNA .k (этап K4, readk/k-import-plan.md, Р3-Р6, Р9).
#include "k_input/k_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
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

using Vec3 = std::array<double, 3>;

// Знаковый объём тетраэдра: det(b-a, c-a, d-a)/6.
double tet_volume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d)
{
    const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const double w[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0])) / 6.0;
}

// Разбиение гексаэдра N1..N8 на 6 тетраэдров по диагонали N1-N7 (при нумерации
// LS-DYNA все объёмы положительны; вырожденные формы дают нулевые тетраэдры).
constexpr int kHexTets[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};

struct SolidGeom
{
    double volume = 0.0;
    Vec3 center{0.0, 0.0, 0.0};
    bool mixed_sign = false;
    double scale3 = 0.0;   // (макс. размах узлов)^3 - для порога нулевого объёма
};

SolidGeom solid_geometry(const std::array<Vec3, 10>& p, int nn)
{
    SolidGeom g;
    std::vector<std::array<int, 4>> tets;
    if (nn == 8)
        for (const auto& t : kHexTets) tets.push_back({t[0], t[1], t[2], t[3]});
    else
        tets.push_back({0, 1, 2, 3});   // 4 узла или углы 10-узлового
    const int used = nn == 8 ? 8 : 4;
    Vec3 lo = p[0], hi = p[0];
    for (int k = 1; k < used; ++k)
        for (int c = 0; c < 3; ++c)
        {
            lo[c] = std::min(lo[c], p[static_cast<std::size_t>(k)][c]);
            hi[c] = std::max(hi[c], p[static_cast<std::size_t>(k)][c]);
        }
    const double ext = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
    g.scale3 = ext * ext * ext;

    double sum = 0.0, maxabs = 0.0;
    Vec3 mom{0.0, 0.0, 0.0};
    std::vector<double> vols;
    for (const auto& t : tets)
    {
        const Vec3 &a = p[static_cast<std::size_t>(t[0])], &b = p[static_cast<std::size_t>(t[1])],
                   &c = p[static_cast<std::size_t>(t[2])], &d = p[static_cast<std::size_t>(t[3])];
        const double v = tet_volume(a, b, c, d);
        vols.push_back(v);
        sum += v;
        maxabs = std::max(maxabs, std::fabs(v));
        for (int k = 0; k < 3; ++k) mom[k] += v * 0.25 * (a[k] + b[k] + c[k] + d[k]);
    }
    bool pos = false, neg = false;
    for (double v : vols)
    {
        if (v > 1e-12 * maxabs) pos = true;
        if (v < -1e-12 * maxabs) neg = true;
    }
    g.mixed_sign = pos && neg;
    g.volume = std::fabs(sum);
    if (sum != 0.0)
        for (int k = 0; k < 3; ++k) g.center[k] = mom[k] / sum;
    return g;
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
    const bool solid_parts = opts.solid_elements == KSolidElements::Particles;
    if (solid_parts && !deck.solids.empty())
    {
        if (dim != 3) fail(opts.k_file, "*ELEMENT_SOLID → частицы требует SPH_DIM=3");
        if (!table) fail(opts.k_file, "k_solid_elements=particles требует таблицу деталей (k_parts_file)");
    }
    if (!solid_parts && !deck.solids.empty())
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

    // K9: частицы из *ELEMENT_SOLID (после всех SPH, в порядке файла, Р9).
    struct SolidParticle { Vec3 pos; double volume; const SolidPart* sp; const KSolid* el; };
    std::vector<SolidParticle> sol;
    if (solid_parts && !deck.solids.empty())
    {
        std::set<long long> sph_nodes, beam_nodes;
        for (const KSph& e : deck.sph) sph_nodes.insert(e.nid);
        for (const KBeam& b : deck.beams) { beam_nodes.insert(b.n1); beam_nodes.insert(b.n2); }
        for (const KSolid& el : deck.solids)
        {
            const SolidPart* sp = table->find_solid(el.pid);
            if (!sp)
                fail(table->file, "PID " + std::to_string(el.pid) + " из " + deck.where(el.src)
                                      + " не описан в SOLID_PART");
            std::array<Vec3, 10> pts{};
            for (int k = 0; k < el.nn; ++k)
            {
                const long long nid = el.n[static_cast<std::size_t>(k)];
                const KNode* nd = deck.find_node(nid);
                if (!nd) fail(deck.where(el.src), "узел " + std::to_string(nid) + " не найден");
                if (sph_nodes.count(nid))
                    fail(deck.where(el.src), "узел " + std::to_string(nid) + " одновременно в *ELEMENT_SOLID и в *ELEMENT_SPH");
                if (beam_nodes.count(nid))
                    fail(deck.where(el.src), "узел " + std::to_string(nid) + " одновременно в *ELEMENT_SOLID и в *ELEMENT_BEAM");
                pts[static_cast<std::size_t>(k)] = {nd->x * opts.length_scale, nd->y * opts.length_scale,
                                                    nd->z * opts.length_scale};
            }
            const SolidGeom g = solid_geometry(pts, el.nn);
            if (!std::isfinite(g.volume) || g.volume <= 1e-12 * g.scale3)
                fail(deck.where(el.src), "*ELEMENT_SOLID EID " + std::to_string(el.eid)
                                             + ": нулевой или неконечный объём V=" + std::to_string(g.volume));
            if (g.mixed_sign)
                summary.warnings.push_back(deck.where(el.src) + ": *ELEMENT_SOLID EID " + std::to_string(el.eid)
                                           + ": тетраэдры разбиения разной ориентации (вывернутый элемент?)");
            sol.push_back({g.center, g.volume, sp, &el});
            for (int c = 0; c < 3; ++c)
            {
                if (n == 0 && sol.size() == 1) { lo[c] = hi[c] = g.center[c]; }
                else { lo[c] = std::min(lo[c], g.center[c]); hi[c] = std::max(hi[c], g.center[c]); }
            }
        }
    }
    summary.solid_element_particles = sol.size();
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

    sys.resize(n + sol.size());
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

    for (std::size_t k = 0; k < sol.size(); ++k)
    {
        const std::size_t i = n + k;
        const SolidParticle& s = sol[k];
        const SolidPart* sp = s.sp;
        sys.inzializ_index(i);
        for (std::size_t d = 0; d < dim; ++d) sys.x[d][i] = s.pos[d];
        sys.itype[i] = sp->itype;
        sys.body[i] = static_cast<int>(s.el->pid);
        sys.rho[i] = sp->rho;
        sys.rho_0[i] = sp->rho;
        sys.mass[i] = sp->rho * s.volume;
        const double dx = std::cbrt(s.volume);
        dxs.push_back(dx);
        sys.hsml[i] = sp->has_hsml ? sp->hsml : 2.0 * dx;
        for (int c = 0; c < sp->nv; ++c) sys.vx[c][i] = sp->v[c];
    }

    if (table)
    {
        for (std::size_t i = 0; i < n + sol.size(); ++i) urs(mats, sys, i);
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

void build_fibers_from_deck(const KDeck& deck, const KPartsTable& table, const KImportOptions& opts,
                            fiber::Model& model, KImportSummary& summary, std::vector<long long>* bead_nids)
{
    fiber::check_fiber_requires_dim3();
    if (deck.beams.empty())
        fail(deck.files.empty() ? opts.k_file : deck.files[0],
             "в режиме solid_fiber нет ни одной нити (нет *ELEMENT_BEAM)");

    const std::size_t nb = deck.beams.size();
    // Граф "узел - балка" по N1/N2 (N3 - узел ориентации - игнорируется).
    std::unordered_map<long long, std::vector<std::size_t>> adj;
    std::map<std::pair<long long, long long>, std::size_t> edges;
    for (std::size_t b = 0; b < nb; ++b)
    {
        const KBeam& e = deck.beams[b];
        if (e.n1 == e.n2)
            fail(deck.where(e.src), "*ELEMENT_BEAM EID " + std::to_string(e.eid) + ": N1 == N2 (NID "
                                        + std::to_string(e.n1) + ")");
        const auto key = std::make_pair(std::min(e.n1, e.n2), std::max(e.n1, e.n2));
        const auto ins = edges.emplace(key, b);
        if (!ins.second)
            fail(deck.where(e.src), "дубликат ребра (узлы " + std::to_string(key.first) + ", "
                                        + std::to_string(key.second) + "): уже в балке EID "
                                        + std::to_string(deck.beams[ins.first->second].eid) + " ("
                                        + deck.where(deck.beams[ins.first->second].src) + ")");
        for (long long nid : {e.n1, e.n2})
        {
            auto& v = adj[nid];
            v.push_back(b);
            if (v.size() > 2)
                fail(deck.where(e.src), "ветвление нити: узел NID " + std::to_string(nid) + " принадлежит "
                                            + std::to_string(v.size()) + " балкам (допустимо не более 2)");
        }
    }

    const double ls = opts.length_scale;
    std::vector<char> seen(nb, 0);
    std::vector<fiber::PolylineSpec> lines;
    std::vector<long long> nids_order;   // NID бусин в порядке сборки (бусина k <-> nids_order[k])
    for (std::size_t b0 = 0; b0 < nb; ++b0)
    {
        if (seen[b0]) continue;
        // Компонента: обход балок от b0 по смежным узлам.
        std::vector<std::size_t> comp;
        std::vector<std::size_t> stack{b0};
        seen[b0] = 1;
        while (!stack.empty())
        {
            const std::size_t b = stack.back();
            stack.pop_back();
            comp.push_back(b);
            for (long long nid : {deck.beams[b].n1, deck.beams[b].n2})
                for (std::size_t o : adj[nid])
                    if (!seen[o]) { seen[o] = 1; stack.push_back(o); }
        }
        // Один PID на нить.
        const KBeam& first = deck.beams[b0];
        for (std::size_t b : comp)
            if (deck.beams[b].pid != first.pid)
                fail(deck.where(deck.beams[b].src),
                     "нить из балок с разными PID: " + std::to_string(deck.beams[b].pid) + " и "
                         + std::to_string(first.pid) + " (первая балка: " + deck.where(first.src) + ")");
        const FiberPart* fp = table.find_fiber(first.pid);
        if (!fp)
            fail(table.file, "PID " + std::to_string(first.pid) + " из " + deck.where(first.src)
                                 + " не описан в FIBER_PART");
        // Концы: узлы степени 1; без концов - кольцо.
        long long start = 0;
        bool have = false;
        long long ring_nid = first.n1;
        for (std::size_t b : comp)
            for (long long nid : {deck.beams[b].n1, deck.beams[b].n2})
            {
                ring_nid = std::min(ring_nid, nid);
                if (adj[nid].size() == 1 && (!have || nid < start)) { start = nid; have = true; }
            }
        if (!have)
            fail(deck.where(first.src), "замкнутое кольцо из балок (узел NID " + std::to_string(ring_nid)
                                            + "), нить должна быть незамкнутой цепочкой");
        // Обход от start.
        fiber::PolylineSpec ps;
        ps.material_id = fp->material_id;
        long long cur = start;
        std::size_t prev_b = nb;
        for (;;)
        {
            const KNode* nd = deck.find_node(cur);
            if (!nd) fail(deck.where(first.src), "узел " + std::to_string(cur) + " не найден");
            ps.points.push_back({nd->x * ls, nd->y * ls, nd->z * ls});
            nids_order.push_back(cur);
            std::size_t nxt = nb;
            for (std::size_t o : adj[cur])
                if (o != prev_b) nxt = o;
            if (nxt == nb) break;
            prev_b = nxt;
            cur = deck.beams[nxt].n1 == cur ? deck.beams[nxt].n2 : deck.beams[nxt].n1;
        }
        lines.push_back(std::move(ps));
    }

    fiber::build_model_from_polylines(lines, table.materials, {}, model, opts.k_file);

    if (nids_order.size() != model.fiber.size())
        fail(opts.k_file, "внутренняя ошибка: число бусин не совпало с числом узлов нитей");
    if (bead_nids) *bead_nids = std::move(nids_order);

    summary.fibers = lines.size();
    summary.beads = model.fiber.size();
    summary.bonds = model.topology.bonds.size();
    bool first_bond = true;
    for (const auto& bd : model.topology.bonds)
    {
        const double d = model.materials[static_cast<std::size_t>(model.fiber.material_id[bd.i])].d;
        const double ld = d > 0.0 ? bd.l0 / d : 0.0;
        if (first_bond)
        {
            summary.seg_min = summary.seg_max = bd.l0;
            summary.ld_min = summary.ld_max = ld;
            first_bond = false;
        }
        else
        {
            summary.seg_min = std::min(summary.seg_min, bd.l0);
            summary.seg_max = std::max(summary.seg_max, bd.l0);
            summary.ld_min = std::min(summary.ld_min, ld);
            summary.ld_max = std::max(summary.ld_max, ld);
        }
    }
}

void apply_spc_clamps(const KDeck& deck, const std::vector<long long>& bead_nids, fiber::Model& model,
                      KImportSummary& summary)
{
    if (deck.spc_nodes.empty() && deck.spc_sets.empty()) return;

    std::unordered_map<long long, std::size_t> solid_of, bead_of;
    for (std::size_t i = 0; i < deck.sph.size(); ++i) solid_of.emplace(deck.sph[i].nid, i);
    for (std::size_t k = 0; k < bead_nids.size(); ++k) bead_of.emplace(bead_nids[k], k);

    const std::size_t n_solid = model.solid.size();
    const std::size_t n_bead = model.fiber.size();
    std::vector<std::uint8_t> sm(n_solid, 0), fm(n_bead, 0);
    if (model.clamps.fiber_mask.size() == n_bead) fm = model.clamps.fiber_mask;   // CLAMP_BOX: OR
    if (model.clamps.solid_mask.size() == n_solid) sm = model.clamps.solid_mask;

    std::set<long long> no_role;
    std::set<std::string> high_dof;   // предупреждения о компонентах >= dim (без повторов)

    auto apply = [&](long long nid, int cid, const int (&dof)[3], const SrcLoc& src, const char* kw)
    {
        if (cid != 0)
            fail(deck.where(src), std::string(kw) + ": CID=" + std::to_string(cid)
                                      + " (локальная система координат) не поддерживается");
        unsigned bits = 0;
        for (int c = 0; c < 3; ++c)
        {
            if (dof[c] == 0) continue;
            if (c < static_cast<int>(dim)) bits |= 1u << c;
            else high_dof.insert(std::string(kw) + ": DOF" + "XYZ"[c] + " закреплён, но dim=" + std::to_string(dim)
                                 + " - компонента отсутствует в модели, пропущено");
        }
        const auto si = solid_of.find(nid);
        const auto bi = bead_of.find(nid);
        if (si == solid_of.end() && bi == bead_of.end()) { no_role.insert(nid); return; }
        if (si != solid_of.end() && si->second < n_solid) sm[si->second] |= static_cast<std::uint8_t>(bits);
        if (bi != bead_of.end()) fm[bi->second] |= static_cast<std::uint8_t>(bits);
    };

    for (const KSpcNode& c : deck.spc_nodes)
    {
        const int dof[3] = {c.dofx, c.dofy, c.dofz};
        apply(c.nid, c.cid, dof, c.src, "*BOUNDARY_SPC_NODE");
    }
    for (const KSpcSet& c : deck.spc_sets)
    {
        const int dof[3] = {c.dofx, c.dofy, c.dofz};
        const auto it = deck.node_sets.find(c.nsid);
        if (it == deck.node_sets.end())
            fail(deck.where(c.src), "*BOUNDARY_SPC_SET: набор узлов " + std::to_string(c.nsid) + " не найден");
        for (long long nid : it->second.nids) apply(nid, c.cid, dof, c.src, "*BOUNDARY_SPC_SET");
    }

    for (const std::string& w : high_dof) summary.warnings.push_back(w);
    if (!no_role.empty())
        summary.warnings.push_back("*BOUNDARY_SPC_*: " + std::to_string(no_role.size())
                                   + " узлов не являются ни SPH-частицей, ни бусиной нити - закрепление пропущено");

    auto count = [](const std::vector<std::uint8_t>& m)
    {
        std::size_t k = 0;
        for (std::uint8_t v : m) k += v != 0;
        return k;
    };
    summary.clamped_particles = count(sm);
    summary.clamped_beads = count(fm);
    // Пустая маска = "зажимов нет" (как в .fib).
    model.clamps.solid_mask = summary.clamped_particles ? std::move(sm) : std::vector<std::uint8_t>{};
    model.clamps.fiber_mask = summary.clamped_beads ? std::move(fm) : std::vector<std::uint8_t>{};
}

KImportSummary load_k_geometry(const KImportOptions& opts, const Materials& mats, fiber::Model& model,
                               std::size_t& moni_particle)
{
    const bool fib = opts.mode == KMode::SolidFiber;
    if (fib && opts.parts_file.empty())
        throw std::runtime_error("k_mode=solid_fiber требует таблицу деталей (k_parts_file)");
    if (fib) fiber::check_fiber_requires_dim3();
    const KDeck deck = read_k_deck(opts.k_file);
    KImportSummary summary;
    KPartsTable table;
    const KPartsTable* tp = nullptr;
    if (!opts.parts_file.empty())
    {
        table = read_k_parts_file(opts.parts_file);
        check_parts_against_deck(table, deck, fib, summary.warnings, opts.solid_elements == KSolidElements::Particles);
        tp = &table;
    }
    build_solid_from_deck(deck, tp, opts, mats, model.solid, moni_particle, summary);
    std::vector<long long> bead_nids;
    if (fib) build_fibers_from_deck(deck, table, opts, model, summary, &bead_nids);
    apply_spc_clamps(deck, bead_nids, model, summary);
    return summary;
}

} // namespace k_input
