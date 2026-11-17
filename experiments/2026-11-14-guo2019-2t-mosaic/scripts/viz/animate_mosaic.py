#!/usr/bin/env python3
"""Мозаичная ткань (--random-cut): копия experiments/2026-11-09-fabric-kT-ball1g-v550/scripts/viz/animate_ply.py с режимом
--color=frag (по умолчанию): нити окрашены по фрагментам раскройки (guo2019_fabric.cut.json, диапазоны бусин; слой 2
темнее слоя 1), стежки -- тонкие чёрные отрезки (в слоях и сквозные), порванные стежки -- красные крестики в середине
стежка по текущим позициям его бусин; в заголовке -- число порванных стежков по кадру. Без cut.json --color=frag
недоступен.

Копия experiments/2026-10-26-guo2019-damage-summary/scripts/viz/animate_ply.py для пакетов kT с большим числом слоёв
(24T и далее): при n_ply > 6 слои окрашены непрерывной шкалой viridis с цветовой полосой вместо легенды, а счётчики
разрывов в заголовке сжаты (всего, слои с разрывами, целые слои). При n_ply <= 6 поведение как у оригинала.

Анимация и снимки прогона Guo 2019, конфигурации 1T/2T/4T (многослойная ткань), на основе animate_1t.py;
в стиле
experiments/2026-09-27-fiber-layers-demo/scripts/animate.py и tools/clamp_animate.py.

Каталог прогона (--run-dir) содержит
    sph_output_<NNNNN>.txt     -- частицы Solid (здесь только шарик),
    fiber_output_<NNNNN>.txt   -- бусины и связи ткани (формат src/fiber_output.cpp),
    monitor.csv                -- построчный монитор (step,t,dt,z_c,v_z,...,w_f,...),
    guo2019_fabric_layout.txt  -- раскладка нитей (направление, диапазон бусин); необязательно.

Рисуется (только реальные кадры, без интерполяции):
  * слева -- 3D: нити ломаными по связям, цвет -- eps/eps_break (фиксированная шкала
    0..1); порванные связи в кадре отсутствуют. Вертикальные нити (вдоль y, зажаты) -- толще,
    горизонтальные (вдоль x, свободны) -- тоньше. Частицы шарика -- серые точки;
  * справа -- боковое сечение (проекция y-z или x-z бусин из полосы около центральной
    плоскости) с контуром шарика: виден прогиб ткани.
Время на кадре -- время монитора для STEP кадра (кадр шага k пишется после advance, но до
time += dt, поэтому TIME в заголовке fiber_output может быть сдвинут на шаг; расхождение
печатается при загрузке). Оси фиксированы по всем кадрам.

Многослойность: раскладка guo2019_fabric_layout.txt при N>=2 имеет первой колонкой ply (слой 1 -- лицевой,
по нему бьёт шарик); при 1T колонки ply нет (один слой). --color=ply -- нити и бусины по слоям (палитра
из 4 цветов + легенда), --color=strain -- нити по eps/eps_break, бусины в сечении всё равно по слоям.
В заголовке: время, прогиб и число порванных нитей по слоям (n_broken_v_pK / n_broken_h_pK из monitor.csv).
--zoom auto (по умолчанию) и --zlim auto подбирают окно по данным; --zoom 0 -- вся ткань.

Usage:
    python animate_ply.py --run-dir out/4T_scaled_v269.0_r0 --out anim_4t.gif --max-frames 40 --snapshots auto
    python animate_ply.py --run-dir out/... --out anim.gif --color strain --zoom 8e-3 --zlim -30,6
    python animate_ply.py --run-dir out/... --out anim.gif --no-gif --snapshots 0,10,20
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import time
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (регистрирует проекцию 3d)
from mpl_toolkits.mplot3d.art3d import Line3DCollection

FRAME_RE = re.compile(r"(?:sph|fiber)_output_(\d+)\.txt$")


# ---------------------------------------------------------------------------
# Чтение
# ---------------------------------------------------------------------------

def discover_frame_ids(run_dir: str) -> list[str]:
    """Номера кадров, для которых есть ОБА файла (sph и fiber)."""
    sph = {FRAME_RE.search(os.path.basename(f)).group(1)
           for f in glob.glob(os.path.join(run_dir, "sph_output_*.txt"))}
    fib = {FRAME_RE.search(os.path.basename(f)).group(1)
           for f in glob.glob(os.path.join(run_dir, "fiber_output_*.txt"))}
    return sorted(sph & fib, key=int)


def parse_sph_frame(path: str) -> np.ndarray:
    """Позиции частиц (N,3). Колонки берём по заголовку."""
    with open(path) as f:
        cols = f.readline().lstrip("# ").split()
        rows = [ln.split() for ln in f if ln.strip()]
    if not rows:
        return np.zeros((0, 3))
    ci = {n: i for i, n in enumerate(cols)}
    d = np.array(rows, dtype=float)
    return d[:, [ci["x"], ci["y"], ci["z"]]]


class FiberFrame:
    __slots__ = ("time", "step", "pos", "ij", "l0", "lbreak", "elong")


def parse_fiber_frame(path: str) -> FiberFrame:
    """Формат src/fiber_output.cpp: TIME/STEP/N_PARTICLES/N_BONDS, PARTICLES, BONDS."""
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]
    h = lines[0].split()
    hd = {h[i]: h[i + 1] for i in range(0, len(h) - 1, 2)}
    n, nb = int(hd["N_PARTICLES"]), int(hd["N_BONDS"])
    idx = 1
    while idx < len(lines) and lines[idx].lstrip().startswith("#"):
        idx += 1
    if lines[idx].strip() == "PARTICLES":
        idx += 1
    pr = np.array([ln.split()[:6] for ln in lines[idx:idx + n]], dtype=float)
    pos = np.zeros((n, 3))
    pos[pr[:, 0].astype(int)] = pr[:, 3:6]
    idx += n
    while idx < len(lines) and lines[idx].lstrip().startswith("#"):
        idx += 1
    if idx < len(lines) and lines[idx].strip() == "BONDS":
        idx += 1
    rows = np.array([ln.split()[:6] for ln in lines[idx:idx + nb]], dtype=float).reshape(-1, 6)
    ff = FiberFrame()
    ff.time, ff.step = float(hd["TIME"]), int(float(hd["STEP"]))
    ff.pos = pos
    ff.ij = rows[:, :2].astype(int)
    ff.l0, ff.lbreak, ff.elong = rows[:, 2], rows[:, 4], rows[:, 5]
    return ff


def read_monitor_file(path: str) -> dict[str, np.ndarray]:
    """monitor.csv; недописанная последняя строка (прогон ещё идёт) отбрасывается."""
    with open(path) as f:
        names = f.readline().strip().split(",")
        rows = [ln.split(",") for ln in f]
    rows = [r for r in rows if len(r) == len(names)]
    d = np.array(rows, dtype=float).reshape(-1, len(names))
    return {n: d[:, i] for i, n in enumerate(names)}


def read_monitor(run_dir: str) -> dict[str, np.ndarray] | None:
    path = os.path.join(run_dir, "monitor.csv")
    return read_monitor_file(path) if os.path.exists(path) else None


def read_layout(run_dir: str, n_beads: int):
    """Раскладка нитей. Возвращает (dirs, yidx, ply, pitch):
    dirs -- 0: x (горизонтальная), 1: y (вертикальная, зажатая), -1 неизвестно;
    yidx -- номер нити в своём направлении (-1 неизвестно); ply -- слой 0-based (-1 неизвестно;
    при 1T без колонки ply все бусины -- слой 0); pitch -- шаг нитей, м, либо 0.
    Формат: '# [ply] object_id direction index offset clamped first_bead last_bead ...' (ply -- только при N>=2)."""
    dirs = -np.ones(n_beads, dtype=int)
    yidx = -np.ones(n_beads, dtype=int)
    ply = -np.ones(n_beads, dtype=int)
    path = os.path.join(run_dir, "guo2019_fabric_layout.txt")
    offs: dict[int, list[float]] = {}
    if os.path.exists(path):
        has_ply = False
        with open(path) as f:
            for ln in f:
                p = ln.split()
                if not p:
                    continue
                if p[0].startswith("#"):
                    has_ply = len(p) > 1 and p[1] == "ply"
                    continue
                o = 1 if has_ply else 0
                if len(p) < 7 + o:
                    continue
                pl = int(p[0]) - 1 if has_ply else 0
                a, b = int(p[o + 5]), int(p[o + 6])
                dirs[a:b + 1] = 0 if p[o + 1] == "x" else 1
                yidx[a:b + 1] = int(p[o + 2])
                ply[a:b + 1] = pl
                if p[o + 1] == "x":
                    offs.setdefault(pl, []).append(float(p[o + 3]))
    ds = [np.median(np.diff(sorted(v))) for v in offs.values() if len(v) > 1]
    pitch = float(np.median(ds)) if ds else 0.0
    return dirs, yidx, ply, pitch


def read_cut(run_dir: str, n_beads: int):
    """guo2019_fabric.cut.json: (frag (n_beads,) -- номер фрагмента бусины или -1, stitches (M,2) -- пары бусин стежков,
    through (M,) -- сквозной ли стежок) либо None, если файла нет."""
    import json
    path = os.path.join(run_dir, "guo2019_fabric.cut.json")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    frag = -np.ones(n_beads, dtype=int)
    for r in d["ranges"]:
        frag[r["bead_first"]:r["bead_first"] + r["bead_count"]] = r["frag"]
    cols = d["stitch_columns"]
    st = np.array(d["stitches"], dtype=np.int64).reshape(-1, len(cols))
    return frag, st[:, [cols.index("a"), cols.index("b")]], st[:, cols.index("through")].astype(bool), float(d["p"])


def frag_palette(n_frag: int, n_ply: int) -> np.ndarray:
    """(n_ply, n_frag) цветов: фрагменты -- tab20 + tab20b вперемешку, слой p > 0 -- затемнён."""
    from matplotlib.colors import to_rgb, to_hex
    base = [plt.get_cmap("tab20")(i) for i in range(20)] + [plt.get_cmap("tab20b")(i) for i in range(20)]
    rng = np.random.default_rng(7)
    order = rng.permutation(len(base))
    out = np.empty((max(n_ply, 1), max(n_frag, 1)), dtype=object)
    for p in range(out.shape[0]):
        k = 1.0 - 0.35 * p / max(out.shape[0] - 1, 1)
        for i in range(out.shape[1]):
            out[p, i] = to_hex(np.array(to_rgb(base[order[i % len(base)]])) * k)
    return out


PLY_COLORS = ["#D55E00","#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]  # Okabe-Ito; слой 1 -- вермильон


def ratio_of(ff: FiberFrame) -> np.ndarray:
    """eps / eps_break в [0,1]; сжатие (elong<0) даёт 0."""
    eb = np.maximum(ff.lbreak / ff.l0 - 1.0, 1e-12)
    return np.clip(ff.elong / eb, 0.0, 1.0)


# ---------------------------------------------------------------------------

def main() -> None:
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(encoding="utf-8")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="каталог прогона (sph_output_*, fiber_output_*, monitor.csv)")
    ap.add_argument("--out", required=True, help="выходной .gif; снимки -- рядом: <имя>_frame<k>.png")
    ap.add_argument("--color", choices=["frag", "ply", "strain"], default="frag",
                    help="цвет нитей: frag -- по фрагментам раскройки, ply -- по слоям, strain -- eps/eps_break "
                         "(бусины сечения: frag -- по фрагментам, иначе по слоям)")
    ap.add_argument("--no-stitches", action="store_true", help="не рисовать целые стежки (порванные -- рисуются)")
    ap.add_argument("--fps", type=float, default=6.0)
    ap.add_argument("--dpi", type=int, default=75)
    ap.add_argument("--snap-dpi-factor", type=float, default=1.5, help="dpi снимков = dpi * множитель")
    ap.add_argument("--max-frames", type=int, default=0,
                    help="если кадров больше -- равномерно проредить до стольких (0 = все)")
    ap.add_argument("--snapshots", default=None,
                    help="индексы кадров через запятую (в выбранной последовательности) или 'auto' "
                         "(первый, первый разрыв, максимальный прогиб, последний)")
    ap.add_argument("--zoom", default="auto",
                    help="полуразмер окна вокруг места удара по x,y, м; 'auto' (по умолчанию) -- по области "
                         "прогиба и размеру шарика; 0 -- вся ткань")
    ap.add_argument("--section", choices=["yz", "xz"], default="yz",
                    help="проекция бокового сечения (по умолчанию yz: вертикальные нити вдоль y)")
    ap.add_argument("--band", type=float, default=0.0,
                    help="полуширина полосы сечения по поперечной оси, м (0 = 1.5 шага нитей)")
    ap.add_argument("--zlim", default="auto",
                    help="пределы по z в мм 'lo,hi' или 'auto' (по бусинам окна и шарику по всем кадрам)")
    ap.add_argument("--yarn-stride", type=int, default=1,
                    help="в 3D рисовать только каждую N-ю нить в каждом направлении (сечение показывает все бусины полосы)")
    ap.add_argument("--no-clip", action="store_true",
                    help="при окне не отбрасывать в 3D связи вне окна (по умолчанию рисуются только связи с обоими концами в окне)")
    ap.add_argument("--lw", type=float, default=1.0, help="множитель толщины линий 3D (вертикальные 1.6, горизонтальные 0.9)")
    ap.add_argument("--vmax", type=float, default=1.0, help="(--color strain) верх шкалы eps/eps_break")
    ap.add_argument("--cmap", default="turbo")
    ap.add_argument("--ball-alpha", type=float, default=0.35, help="прозрачность точек шарика в 3D")
    ap.add_argument("--elev", type=float, default=24.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--no-gif", action="store_true", help="только снимки")
    args = ap.parse_args()
    by_frag = args.color == "frag"
    by_ply = args.color in ("ply", "frag")   # цвет задаётся массивом цветов, а не шкалой

    ids = discover_frame_ids(args.run_dir)
    if not ids:
        print(f"{args.run_dir}: нет пар sph_output/fiber_output", file=sys.stderr)
        sys.exit(1)
    if args.max_frames and len(ids) > args.max_frames:
        sel = np.linspace(0, len(ids) - 1, args.max_frames).round().astype(int)
        ids = [ids[i] for i in sorted(set(sel.tolist()))]
    nf = len(ids)
    print(f"Загрузка {nf} кадров из {args.run_dir} ...")
    t_load = time.time()
    ball = [parse_sph_frame(os.path.join(args.run_dir, f"sph_output_{i}.txt")) for i in ids]
    fib = [parse_fiber_frame(os.path.join(args.run_dir, f"fiber_output_{i}.txt")) for i in ids]
    print(f"  чтение кадров: {time.time() - t_load:.1f} с ({(time.time() - t_load) / nf:.2f} с/кадр)")
    mon = read_monitor(args.run_dir)
    dirs, yidx, ply, pitch = read_layout(args.run_dir, len(fib[0].pos))
    have_ply = (ply >= 0).any()
    n_ply = int(ply.max()) + 1 if have_ply else 1
    if not have_ply:
        ply = np.zeros(len(fib[0].pos), dtype=int)
    ply = np.where(ply < 0, 0, ply)
    band = args.band or (1.5 * pitch if pitch > 0 else 1.5e-3)
    many = n_ply > len(PLY_COLORS)
    if many:
        from matplotlib.colors import to_hex
        pcol = np.array([to_hex(plt.get_cmap("viridis")(p / max(n_ply - 1, 1))) for p in range(n_ply)])
    else:
        pcol = np.array([PLY_COLORS[p % len(PLY_COLORS)] for p in range(max(n_ply, 1))])

    # --- мозаика: фрагменты и стежки ---
    nB = len(fib[0].pos)
    cut = read_cut(args.run_dir, nB)
    if by_frag and cut is None:
        print("--color=frag: нет guo2019_fabric.cut.json", file=sys.stderr)
        sys.exit(1)
    if cut is not None:
        frag, st_ab, st_thr, pitch = cut   # шаг нитей -- из раскройки (смещения раскладки у фрагментов свои)
        band = args.band or 1.5 * pitch
        n_frag = int(frag.max()) + 1
        fpal = frag_palette(n_frag, n_ply)
        bcol = np.array([fpal[ply[i], max(frag[i], 0)] for i in range(nB)])    # цвет бусины: (слой, фрагмент)
        st_key = np.minimum(st_ab[:, 0], st_ab[:, 1]) * nB + np.maximum(st_ab[:, 0], st_ab[:, 1])
        print(f"Мозаика: фрагментов {n_frag}, стежков {len(st_key)} (сквозных {int(st_thr.sum())})")
    else:
        st_ab, st_key = np.zeros((0, 2), dtype=np.int64), np.zeros(0, dtype=np.int64)
        bcol = pcol[ply]

    def bond_keys(f):
        a, b = f.ij[:, 0].astype(np.int64), f.ij[:, 1].astype(np.int64)
        return np.minimum(a, b) * nB + np.maximum(a, b)

    is_st = [np.isin(bond_keys(f), st_key) for f in fib]                     # связь -- стежок
    broken_st = [~np.isin(st_key, bond_keys(f)) for f in fib]                # стежок порван
    pcol_b = bcol if by_frag else pcol[ply]                                  # цвет бусин и нитей по бусине

    # --- время и величины монитора для каждого кадра (по STEP) ---
    hdr_t = np.array([f.time for f in fib])
    times = hdr_t.copy()
    ms = None
    pk_v = [f"n_broken_v_p{p + 1}" for p in range(n_ply)]
    pk_h = [f"n_broken_h_p{p + 1}" for p in range(n_ply)]
    if mon is not None:
        ms = {}
        steps = mon["step"].astype(int)
        order = {s: r for r, s in enumerate(steps)}
        rows, exact = [], []
        for f in fib:
            r = order.get(f.step)
            exact.append(r is not None)
            rows.append(r if r is not None else int(np.argmin(np.abs(steps - f.step))))
        exact = np.array(exact)
        for key in ["t", "v_z", "w_f", "n_broken_v", "n_broken_h", "z_c"] + pk_v + pk_h:
            if key in mon:
                ms[key] = mon[key][rows]
        if "t" in ms and exact.any():
            times = np.where(exact, ms["t"], hdr_t)
            diff = np.abs(ms["t"][exact] - hdr_t[exact]).max()
            print(f"Время: монитор для STEP кадра; max|t_монитор - TIME_заголовка| = {diff * 1e6:.4f} мкс "
                  f"(кадров без строки монитора: {int((~exact).sum())}, для них берётся TIME заголовка)")
    else:
        print("monitor.csv нет: время берётся из заголовка fiber_output", file=sys.stderr)
    per_ply_cnt = n_ply > 1 and ms is not None and all(k in ms for k in pk_v)

    # --- окно и границы осей ---
    c0 = fib[0]
    ctr_xy = 0.5 * (c0.pos.min(axis=0) + c0.pos.max(axis=0))[:2]
    r_ball = 0.5 * float(np.ptp(ball[0][:, 2])) if len(ball[0]) else 3e-3
    half_fabric = 0.5 * float((c0.pos.max(axis=0) - c0.pos.min(axis=0))[:2].max())
    if str(args.zoom).lower() == "auto":
        dz = np.max([np.abs(f.pos[:, 2] - c0.pos[:, 2]) for f in fib], axis=0)   # макс. смещение бусин по z
        moved = dz > 0.05 * max(dz.max(), 1e-12)
        rr = float(np.max(np.hypot(c0.pos[moved, 0] - ctr_xy[0], c0.pos[moved, 1] - ctr_xy[1]))) if moved.any() else 0.0
        zoom = float(min(max(1.25 * rr, 4.0 * r_ball), half_fabric))
        if zoom > 0.9 * half_fabric:
            zoom = 0.0
    else:
        zoom = float(args.zoom)
    clip_on = zoom > 0 and not args.no_clip

    def in_win(pos):
        if zoom <= 0:
            return np.ones(len(pos), dtype=bool)
        return (np.abs(pos[:, 0] - ctr_xy[0]) <= zoom) & (np.abs(pos[:, 1] - ctr_xy[1]) <= zoom)

    lo = np.zeros(3)
    hi = np.zeros(3)
    if zoom > 0:
        lo[:2], hi[:2] = ctr_xy - zoom, ctr_xy + zoom
    else:
        allxy = np.concatenate([f.pos for f in fib], axis=0)
        lo[:2], hi[:2] = allxy.min(axis=0)[:2], allxy.max(axis=0)[:2]
    if args.zlim != "auto":
        lo[2], hi[2] = (float(v) * 1e-3 for v in args.zlim.split(","))
    else:
        zs = np.concatenate([f.pos[in_win(f.pos), 2] for f in fib] + [b[in_win(b), 2] for b in ball if len(b)])
        lo[2], hi[2] = zs.min(), zs.max()
    pad = 0.06 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    span = hi - lo
    aspect = np.maximum(span / span[:2].max(), 0.25)
    print(f"Окно: zoom = {zoom * 1e3:.2f} мм ({'вся ткань' if zoom == 0 else 'полуразмер'}); слоёв: {n_ply}; "
          f"пределы (мм): x[{lo[0]*1e3:.2f},{hi[0]*1e3:.2f}] y[{lo[1]*1e3:.2f},{hi[1]*1e3:.2f}] "
          f"z[{lo[2]*1e3:.2f},{hi[2]*1e3:.2f}]; шаг нитей {pitch*1e3:.3f} мм; полоса сечения ±{band*1e3:.2f} мм")

    # --- фигура ---
    fig = plt.figure(figsize=(13, 6.8))
    ax = fig.add_axes((0.0, 0.03, 0.55, 0.84), projection="3d")
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(tuple(aspect))
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_xlabel("x, м"); ax.set_ylabel("y, м"); ax.set_zlabel("z, м")
    ax.tick_params(labelsize=7)
    for a, axis in zip(aspect, (ax.xaxis, ax.yaxis, ax.zaxis)):
        axis.set_major_locator(MaxNLocator(nbins=5 if a > 0.5 else 3))

    stride = max(args.yarn_stride, 1)
    have_layout = len(yidx) == len(c0.pos) and (yidx >= 0).any()
    if stride > 1 and not have_layout:
        warnings.warn("--yarn-stride без guo2019_fabric_layout.txt игнорируется")

    lc3 = Line3DCollection(np.zeros((0, 2, 3)), cmap=args.cmap, clim=(0.0, args.vmax), linewidths=1.2)
    ax.add_collection3d(lc3, autolim=False)
    sc3 = ax.scatter(*ball[0].T, c="0.35", s=3, alpha=args.ball_alpha, depthshade=True)
    br3 = ax.scatter([], [], [], c="red", marker="x", s=22, linewidths=1.2, depthshade=False)
    handles = [Line2D([], [], color="0.35", marker=".", ls="", label="шарик")]
    if len(st_key):
        handles += [Line2D([], [], color="k", lw=0.8, label="стежок"),
                    Line2D([], [], color="red", marker="x", ls="", label="порванный стежок")]
    if by_frag:
        handles = [Line2D([], [], color=fpal[p, 0], lw=2.5,
                          label=f"слой {p + 1}: цвет -- фрагмент" + (" (темнее)" if p > 0 else ""))
                   for p in range(n_ply)] + handles
    elif by_ply and not many:
        handles = [Line2D([], [], color=pcol[p], lw=2.5,
                          label=f"слой {p + 1}" + (" (лицевой)" if p == 0 and n_ply > 1 else ""))
                   for p in range(n_ply)] + handles
    ax.legend(handles=handles, loc="upper left", fontsize=8)
    if not by_ply:
        cb = fig.colorbar(lc3, ax=ax, shrink=0.5, pad=0.0, location="left")
        cb.set_label(r"$\varepsilon / \varepsilon_{разр}$ связи" + (f" (шкала до {args.vmax:g})" if args.vmax != 1.0 else ""))

    # боковое сечение
    ia = 1 if args.section == "yz" else 0     # горизонтальная ось сечения
    ic = 0 if args.section == "yz" else 1     # поперечная ось полосы
    axS = fig.add_axes((0.62, 0.10, 0.36, 0.74))
    axS.set_xlim(lo[ia], hi[ia]); axS.set_ylim(lo[2], hi[2])
    axS.set_aspect("equal", adjustable="box")
    axS.set_xlabel(f"{'xyz'[ia]}, м"); axS.set_ylabel("z, м")
    axS.set_title(f"Сечение {args.section}: бусины с |{'xyz'[ic]}| < {band*1e3:.2f} мм (цвет точек -- слой)", fontsize=10)
    axS.grid(alpha=0.3)
    lcS = LineCollection([], cmap=args.cmap, clim=(0.0, args.vmax), linewidths=1.0)
    axS.add_collection(lcS)
    dotsS = axS.scatter([], [], s=9, zorder=3)
    ballS, = axS.plot([], [], ".", ms=2, color="0.6")
    circ = Circle((0, 0), r_ball, fill=False, ec="0.2", lw=1.0, ls="--")
    axS.add_patch(circ)
    if many:
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import BoundaryNorm, ListedColormap
        smp = ScalarMappable(norm=BoundaryNorm(np.arange(0.5, n_ply + 1.0), n_ply), cmap=ListedColormap(list(pcol)))
        cbp = fig.colorbar(smp, ax=axS, shrink=0.8, pad=0.02, ticks=[1] + list(range(4, n_ply + 1, 4)))
        cbp.set_label("слой (1 -- лицевой)")
    brS = axS.scatter([], [], c="red", marker="x", s=30, linewidths=1.2, zorder=5)
    if by_frag:
        axS.set_title(f"Сечение {args.section}: бусины с |{'xyz'[ic]}| < {band*1e3:.2f} мм (цвет -- фрагмент)",
                      fontsize=10)
    elif n_ply > 1:
        axS.legend(handles=[Line2D([], [], color=pcol[p], marker="o", ls="-", ms=4, label=f"слой {p + 1}")
                            for p in range(n_ply)], loc="upper right", fontsize=8)
    title = fig.text(0.02, 0.975, "", fontsize=11, va="top", ha="left", linespacing=1.5)
    ply_ok = len(ply) == len(c0.pos)

    def update(k: int):
        f = fib[k]
        seg = np.stack([f.pos[f.ij[:, 0]], f.pos[f.ij[:, 1]]], axis=1)
        r = np.minimum(ratio_of(f), args.vmax)

        if stride > 1 and have_layout:
            ya, yb = yidx[f.ij[:, 0]], yidx[f.ij[:, 1]]
            k3 = ((ya % stride == 0) | (ya < 0)) & ((yb % stride == 0) | (yb < 0))
        else:
            k3 = np.ones(len(f.ij), dtype=bool)
        if clip_on:
            win = in_win(f.pos)
            k3 = k3 & win[f.ij[:, 0]] & win[f.ij[:, 1]]
        stb = is_st[k]
        if args.no_stitches:
            k3 = k3 & ~stb
        lc3.set_segments(seg[k3])
        if by_ply and ply_ok:
            lc3.set_array(None)
            cols = pcol_b[f.ij[k3, 0]].astype(object)
            cols[stb[k3]] = "k"
            lc3.set_color(cols)
        else:
            lc3.set_array(r[k3])
        if len(dirs) == len(f.pos) and (dirs >= 0).any():
            lw = np.where(dirs[f.ij[k3, 0]] == 1, 1.6, 0.9)
            lc3.set_linewidths(args.lw * np.where(stb[k3], 0.5, lw))
        bst = broken_st[k]
        mid = 0.5 * (f.pos[st_ab[bst, 0]] + f.pos[st_ab[bst, 1]]) if bst.any() else np.zeros((0, 3))
        if clip_on and len(mid):
            mid = mid[in_win(mid)]
        br3._offsets3d = (mid[:, 0], mid[:, 1], mid[:, 2])
        sc3._offsets3d = tuple(ball[k].T)
        inb = np.abs(f.pos[:, ic] - ctr_xy[ic]) < band
        bi = inb[f.ij[:, 0]] & inb[f.ij[:, 1]] & ~stb
        lcS.set_segments(seg[bi][:, :, [ia, 2]])
        if by_ply and ply_ok:
            lcS.set_array(None)
            lcS.set_color(pcol_b[f.ij[bi, 0]])
        else:
            lcS.set_array(r[bi])
        dotsS.set_offsets(np.column_stack([f.pos[inb, ia], f.pos[inb, 2]]))
        dotsS.set_facecolor(pcol_b[inb] if ply_ok else "k")
        ms_ = mid[np.abs(mid[:, ic] - ctr_xy[ic]) < 3 * band] if len(mid) else mid
        brS.set_offsets(np.column_stack([ms_[:, ia], ms_[:, 2]]) if len(ms_) else np.zeros((0, 2)))
        dotsS.set_edgecolor("none")
        bb = ball[k]
        if len(bb):
            bsel = np.abs(bb[:, ic] - ctr_xy[ic]) < band
            ballS.set_data(bb[bsel, ia], bb[bsel, 2])
            circ.center = (float(bb[:, ia].mean()), float(0.5 * (bb[:, 2].min() + bb[:, 2].max())))
        s = f"t = {times[k] * 1e6:.2f} мкс   шаг = {f.step}"
        if ms is not None:
            if "v_z" in ms:
                s += f"   v_z шарика = {ms['v_z'][k]:.1f} м/с"
            if "w_f" in ms:
                s += f"   прогиб w_f = {ms['w_f'][k] * 1e3:.2f} мм"
        s2 = ""
        if per_ply_cnt and many:
            hit = [p + 1 for p in range(n_ply) if ms[pk_v[p]][k] + (ms[pk_h[p]][k] if pk_h[p] in ms else 0) > 0]
            s2 = f"порвано нитей верт./гориз. = {int(ms['n_broken_v'][k])}/{int(ms['n_broken_h'][k])}"
            s2 += (f";  слои с разрывами: {min(hit)}-{max(hit)} ({len(hit)} из {n_ply})" if hit else ";  разрывов нет")
        elif per_ply_cnt:
            s2 = "порвано верт. по слоям: " + " ".join(f"{p + 1}:{int(ms[pk_v[p]][k])}" for p in range(n_ply))
            if all(h in ms for h in pk_h):
                s2 += "   гориз.: " + " ".join(f"{p + 1}:{int(ms[pk_h[p]][k])}" for p in range(n_ply))
        elif ms is not None and "n_broken_v" in ms:
            s2 = f"порвано нитей верт./гориз. = {int(ms['n_broken_v'][k])}/{int(ms['n_broken_h'][k])}"
        if len(st_key):
            s2 += (";   " if s2 else "") + f"порвано стежков {int(bst.sum())} из {len(st_key)}"
            if bst.any() and cut is not None:
                s2 += f" (сквозных {int((bst & st_thr).sum())})"
        title.set_text(s + ("\n" + s2 if s2 else ""))

    # --- снимки ---
    base = os.path.splitext(args.out)[0]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    snaps: list[int] = []
    if args.snapshots == "auto":
        snaps = [0]
        if ms is not None and "n_broken_v" in ms:
            hit = np.where(ms["n_broken_v"] + ms["n_broken_h"] > 0)[0]
            if len(hit):
                snaps.append(int(hit[0]))
        if ms is not None and "w_f" in ms:
            snaps.append(int(np.argmax(ms["w_f"])))
        snaps.append(nf - 1)
        snaps = sorted(set(snaps))
    elif args.snapshots:
        snaps = [int(s) for s in args.snapshots.split(",") if s.strip()]
    for si in snaps:
        if 0 <= si < nf:
            t0 = time.time()
            update(si)
            p = f"{base}_frame{si}.png"
            fig.savefig(p, dpi=int(args.dpi * args.snap_dpi_factor))
            print(f"снимок: {p} ({time.time() - t0:.1f} с)")
        else:
            warnings.warn(f"--snapshots {si} вне [0, {nf - 1}]")

    if not args.no_gif:
        anim = FuncAnimation(fig, update, frames=nf, interval=1000.0 / args.fps, blit=False)
        t0 = time.time()
        anim.save(args.out, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
        print(f"GIF: рендер {time.time() - t0:.1f} с ({(time.time() - t0) / nf:.2f} с/кадр)")
        print(f"GIF: {args.out} ({os.path.getsize(args.out) / 1e6:.2f} МБ, {nf} кадров)")


if __name__ == "__main__":
    main()
