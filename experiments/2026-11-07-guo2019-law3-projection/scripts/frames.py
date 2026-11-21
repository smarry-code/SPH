#!/usr/bin/env python3
"""frames.py <rundir> [--out CSV]
Ряды по кадрам fiber_output_*/sph_output_* одного прогона.

Колонки CSV:
  step, t, z_ball          — шаг, время, средняя z частиц шара (м);
  dff_max                  — max δ_ff/d по парам сегментов РАЗНЫХ нитей, δ = d - dist(отрезок, отрезок);
  dff_n_contact            — число пар с δ > 0;
  dff_n_gt09               — число пар с δ_ff/d > 0,9 («проходы», критерий б);
  dfs_max                  — max δ_fs/d = (r_f + r_s - dist)/d  (как в мониторе: delta/(2 r_f));
  n_tunnel                 — число бусин, у которых сегмент ближе r_f + r_s - 0,5 r_f к частице шара;
  n_cross                  — число перекрестий основа x уток (102x102);
  n_cross_flip             — перекрестья, у которых знак (z_x-нить - z_y-нить) изменился относительно кадра 0 (критерий а);
  n_cross_lost             — перекрестья, у которых в кадре не нашлось пересечения проекций нитей в плане (см. Crossings);
  n_bonds                  — число связей в кадре.
Критерии повторяют tests/guo2019_monitor.hpp::diagnostics() (стр. ~1053-1175): сегмент = связь топологии (BONDS кадра),
r_f = d/2, r_s = 0,5*(m/max(rho, rho_0))^(1/3) (interactions::solid_contact_radius, dim=3), rho_0 берётся из
кадра шара самого первого шага (sph_output_00001; rho_0 в кадрах не пишется).
"""
import argparse, os, sys
import numpy as np
from scipy.spatial import cKDTree
from common import (read_fiber_frame, read_sph_frame, read_layout, frame_steps, fiber_diameter, read_shot)


def seg_seg(p1, q1, p2, q2):
    """Расстояние между отрезками [p1,q1] и [p2,q2] (по Эриксону, векторно); возвращает dist, s, t."""
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a = (d1 * d1).sum(1); e = (d2 * d2).sum(1)
    f = (d2 * r).sum(1); c = (d1 * r).sum(1); b = (d1 * d2).sum(1)
    tiny = 1e-30
    a = np.maximum(a, tiny); e = np.maximum(e, tiny)
    den = a * e - b * b
    ok = den > 1e-9 * a * e
    s = np.where(ok, np.clip((b * f - c * e) / np.where(ok, den, 1.0), 0, 1), 0.0)
    t = (b * s + f) / e
    lo = t < 0; hi = t > 1
    t = np.clip(t, 0, 1)
    s = np.where(lo, np.clip(-c / a, 0, 1), np.where(hi, np.clip((b - c) / a, 0, 1), s))
    c1 = p1 + s[:, None] * d1
    c2 = p2 + t[:, None] * d2
    return np.linalg.norm(c1 - c2, axis=1), s, t


def point_seg(p, a, b):
    """Расстояние от точек p до отрезков [a,b] (попарно)."""
    ab = b - a
    L2 = np.maximum((ab * ab).sum(1), 1e-30)
    u = np.clip(((p - a) * ab).sum(1) / L2, 0, 1)
    return np.linalg.norm(p - (a + u[:, None] * ab), axis=1)


def ff_pairs(A, B, d):
    """Кандидатные пары сегментов (i<j) с расстоянием между серединами <= (L_i+L_j)/2 + d.
    «Короткие» сегменты (L <= 1,5 медианы) ищутся одним kd-деревом с общим радиусом; редкие длинные
    (оборванные, улетевшие бусины) — отдельным запросом на каждый."""
    mid = 0.5 * (A + B)
    L = np.linalg.norm(B - A, axis=1)
    Lq = 1.5 * float(np.median(L))
    short = np.nonzero(L <= Lq)[0]
    long_ = np.nonzero(L > Lq)[0]
    pr = []
    ts = cKDTree(mid[short])
    q = ts.query_pairs(r=Lq + d, output_type="ndarray")
    if len(q):
        pr.append(np.stack([short[q[:, 0]], short[q[:, 1]]], 1))
    for i in long_:
        near = ts.query_ball_point(mid[i], r=0.5 * L[i] + 0.5 * Lq + d)
        if near:
            pr.append(np.stack([np.full(len(near), i), short[np.asarray(near)]], 1))
        for j in long_:
            if j > i and np.linalg.norm(mid[i] - mid[j]) <= 0.5 * (L[i] + L[j]) + d:
                pr.append(np.array([[i, j]]))
    return np.concatenate(pr) if pr else np.zeros((0, 2), dtype=np.int64)


def ff_metrics(F, d):
    """δ_ff: пары сегментов разных нитей. Возвращает (max δ/d, n δ>0, n δ/d>0.9)."""
    X = F["X"]; bi, bj = F["bi"], F["bj"]
    A, B = X[bi], X[bj]
    pr = ff_pairs(A, B, d)
    if len(pr) == 0:
        return 0.0, 0, 0
    i, j = pr[:, 0], pr[:, 1]
    m = F["obj"][bi[i]] != F["obj"][bi[j]]
    i, j = i[m], j[m]
    dist = np.empty(len(i))
    for k in range(0, len(i), 400000):          # блоками, чтобы не раздувать память
        sl = slice(k, k + 400000)
        dist[sl] = seg_seg(A[i[sl]], B[i[sl]], A[j[sl]], B[j[sl]])[0]
    delta = (d - dist) / d
    return float(max(delta.max(), 0.0)), int((delta > 0).sum()), int((delta > 0.9).sum())


def fs_metrics(F, S, rho0, d):
    """δ_fs и n_tunnel — как в мониторе: сегмент против частиц ШАРА (все частицы кадра шара)."""
    X = F["X"]; bi, bj = F["bi"], F["bj"]
    A, B = X[bi], X[bj]
    rf = 0.5 * d
    rs = 0.5 * (S["m"] / np.maximum(S["rho"], rho0)) ** (1.0 / 3.0)
    P = S["X"]
    mid = 0.5 * (A + B)
    half = 0.5 * np.linalg.norm(B - A, axis=1)
    tree = cKDTree(P)
    reach = half + rf + rs.max()   # радиус поиска свой для каждого сегмента
    lists = tree.query_ball_point(mid, r=reach)
    cnt = np.fromiter((len(l) for l in lists), dtype=np.int64, count=len(lists))
    if cnt.sum() == 0:
        return 0.0, 0
    sid = np.repeat(np.arange(len(lists)), cnt)
    pid = np.concatenate([np.asarray(l, dtype=np.int64) for l in lists if len(l)])
    dist = point_seg(P[pid], A[sid], B[sid])
    delta = rf + rs[pid] - dist
    dmax = float(max(delta.max(), 0.0) / (2.0 * rf))
    deep = dist < rf + rs[pid] - 0.5 * rf
    tun = np.zeros(len(F["X"]), dtype=bool)
    tun[bi[sid[deep]]] = True
    tun[bj[sid[deep]]] = True
    return dmax, int(tun.sum())


class Crossings:
    """Перекрестья основа x уток, без окон поиска. Нити 'x' (вдоль x, off = y) и 'y' (вдоль y, off = x) — из
    раскладки; сегменты нити — соседние бусины по порядку в раскладке (геометрически, независимо от разрыва связи).
    В каждом кадре: все пары сегментов (x-нить, y-нить), пересекающиеся В ПЛАНЕ xy (расстояние между проекциями
    отрезков < tol*d; кандидаты — kd-дерево по серединам). Для пары нитей (k,l) в кадре 0 выбирается пересечение,
    ближайшее в плане к (off_l, off_k); в следующих кадрах — пересечение, ближайшее по номерам сегментов вдоль обеих
    нитей к найденному в предыдущем кадре (отслеживание материальной точки; концы нитей после втягивания могут давать
    лишние пересечения в плане); знак s = sign(z_x - z_y) в точках пересечения.
    n_cross_flip — число пар (k,l), для которых пересечение найдено в обоих кадрах (0 и текущем) и s изменился.
    n_cross_lost — число пар, у которых в текущем кадре пересечения в плане не нашлось (длинные/оборванные сегменты)."""
    TOL = 0.02          # допуск пересечения в плане, доли d

    def __init__(self, lay):
        self.fx = [f for f in lay if f["dir"] == "x"]
        self.fy = [f for f in lay if f["dir"] == "y"]
        self.n = len(self.fx) * len(self.fy)
        def segs(fibs):
            A, K = [], []
            for k, f in enumerate(fibs):
                i = np.arange(f["first"], f["last"])
                A.append(i); K.append(np.full(len(i), k))
            return np.concatenate(A), np.concatenate(K)
        self.ax, self.kx = segs(self.fx)
        self.ay, self.ky = segs(self.fy)
        self.firstx = np.array([f["first"] for f in self.fx])
        self.firsty = np.array([f["first"] for f in self.fy])
        self.P0 = np.stack(np.broadcast_arrays(np.array([f["off"] for f in self.fy])[None, :],
                                               np.array([f["off"] for f in self.fx])[:, None]), -1)  # (nx, ny, 2)
        self.sign0 = None
        self.ref = None

    def state(self, X, d):
        nx, ny = len(self.fx), len(self.fy)
        a1, b1 = X[self.ax], X[self.ax + 1]
        a2, b2 = X[self.ay], X[self.ay + 1]
        m1, m2 = 0.5 * (a1 + b1)[:, :2], 0.5 * (a2 + b2)[:, :2]
        L = np.concatenate([np.linalg.norm((b1 - a1)[:, :2], axis=1), np.linalg.norm((b2 - a2)[:, :2], axis=1)])
        r = 1.5 * float(np.median(L))
        lists = cKDTree(m1).query_ball_tree(cKDTree(m2), r)
        cnt = np.fromiter((len(l) for l in lists), dtype=np.int64, count=len(lists))
        i1 = np.repeat(np.arange(len(lists)), cnt)
        i2 = np.concatenate([np.asarray(l, dtype=np.int64) for l in lists if len(l)])
        p1, q1, p2, q2 = a1[i1].copy(), b1[i1].copy(), a2[i2].copy(), b2[i2].copy()
        z1 = np.stack([p1[:, 2], q1[:, 2]]); z2 = np.stack([p2[:, 2], q2[:, 2]])
        for q in (p1, q1, p2, q2):
            q[:, 2] = 0.0                                  # расстояние в плане
        dist, s, t = seg_seg(p1, q1, p2, q2)
        ok = dist < self.TOL * d
        i1, i2, s, t = i1[ok], i2[ok], s[ok], t[ok]
        z1, z2 = z1[:, ok], z2[:, ok]
        dz = (z1[0] + s * (z1[1] - z1[0])) - (z2[0] + t * (z2[1] - z2[0]))
        xy = (p1[ok][:, :2] + s[:, None] * (q1[ok] - p1[ok])[:, :2])
        k, l = self.kx[i1], self.ky[i2]
        pid = k * ny + l
        loc1 = self.ax[i1] - self.firstx[k]          # номер сегмента вдоль своей нити (материальная координата)
        loc2 = self.ay[i2] - self.firsty[l]
        if self.ref is None:                          # кадр 0: ближайшее в плане к (off_l, off_k)
            dd = np.linalg.norm(xy - self.P0[k, l], axis=1)
        else:                                         # далее: ближайшее по материальным координатам к предыдущему кадру
            dd = np.abs(loc1 - self.ref[0][pid]) + np.abs(loc2 - self.ref[1][pid])
        o = np.lexsort((dd, pid))
        pid_s = pid[o]
        first = np.concatenate([[True], pid_s[1:] != pid_s[:-1]])
        sel = o[first]
        sign = np.zeros(nx * ny)
        found = np.zeros(nx * ny, dtype=bool)
        sign[pid[sel]] = np.sign(dz[sel])
        found[pid[sel]] = True
        if self.ref is None:
            self.ref = (np.zeros(nx * ny), np.zeros(nx * ny))
        self.ref[0][pid[sel]] = loc1[sel]            # опорные индексы обновляются каждый кадр (цепочка кадров)
        self.ref[1][pid[sel]] = loc2[sel]
        return sign, found

    def count(self, X, d):
        sg, found = self.state(X, d)
        if self.sign0 is None:
            self.sign0, self.found0 = sg, found
        flip = found & self.found0 & (sg != self.sign0) & (self.sign0 != 0)
        return int(flip.sum()), int((~found).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = a.rundir.rstrip("/")
    d = fiber_diameter(run)
    lay = read_layout(os.path.join(run, "guo2019_fabric_layout.txt"))
    fr = frame_steps(run)
    sph = {s: os.path.join(run, "sph_output_%05d.txt" % s) for s in fr}
    steps = list(fr)
    rho0 = read_sph_frame(sph[steps[0]])["rho"]
    cr = None
    rows = []
    for s in steps:
        F = read_fiber_frame(fr[s])
        S = read_sph_frame(sph[s])
        if cr is None:
            cr = Crossings(lay)
        dff, nc, n09 = ff_metrics(F, d)
        dfs, ntun = fs_metrics(F, S, rho0, d)
        nflip, nlost = cr.count(F["X"], d)
        rows.append((s, F["t"], float(S["X"][:, 2].mean()), dff, nc, n09, dfs, ntun, cr.n, nflip, nlost, len(F["bi"])))
        print(os.path.basename(run), rows[-1], file=sys.stderr, flush=True)
    head = "step,t,z_ball,dff_max,dff_n_contact,dff_n_gt09,dfs_max,n_tunnel,n_cross,n_cross_flip,n_cross_lost,n_bonds"
    out = a.out or os.path.join(run, "frames.csv")
    with open(out, "w") as fh:
        fh.write(head + "\n")
        for r in rows:
            fh.write(",".join(repr(x) if isinstance(x, float) else str(x) for x in r) + "\n")
    sh = read_shot(os.path.join(run, "shot.txt"))
    A = np.array(rows, dtype=float)
    print(f"# {os.path.basename(run)}: по кадрам max dff={A[:,3].max():.4g} (shot {sh.get('max_delta_ff/d')}), "
          f"max dfs={A[:,6].max():.4g} (shot {sh.get('max_delta_fs/d')}), max n_tunnel={A[:,7].max():.0f} (shot {sh.get('n_tunnel')}); "
          f"последний кадр: dff={A[-1,3]:.4g} dfs={A[-1,6]:.4g} n_tun={A[-1,7]:.0f}; flip макс {A[:,9].max():.0f}; "
          f"-> {out}")


if __name__ == "__main__":
    main()
