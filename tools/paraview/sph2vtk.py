#!/usr/bin/env python3
"""Перевод кадров расчёта в формат VTK XML для ParaView.

Каталог прогона (--run-dir) может содержать кадры одного из трёх видов:

  * pair     -- sph_output_<N>.txt (частицы Solid) + fiber_output_<N>.txt (бусины и связи
                нитей, формат src/fiber_output.cpp); так пишут guo2019_shot и демо нитей;
  * combined -- только sph_output_<N>.txt, в котором нити записаны вместе с частицами:
                колонки ... T k fiber i j eps eps_pl tension alive; строка бусины несёт
                связь i -> j (i, j -- номера строк с 0), alive = 1, если связь есть;
  * solid    -- только sph_output_<N>.txt без нитей (src/output.cpp).

Вид определяется сам. Необязательные файлы рядом с кадрами:
    monitor.csv                -- построчный монитор (step,t,...): время кадра и величины
                                  монитора берутся из строки с номером шага кадра;
    guo2019_fabric_layout.txt  -- раскладка нитей: направление (x/y) и зажатые бусины.

Результат (--out-dir):
    solid_<N>.vtp,  solid.pvd   -- частицы Solid: точки + все колонки кадра как PointData,
                                   velocity (вектор), speed;
    fibers_<N>.vtp, fibers.pvd  -- бусины (точки) и связи (отрезки-ячейки); PointData бусин,
                                   CellData связей (eps, eps_ratio = eps/eps_разр, tension, ...).
    breaks_<N>.vtp, breaks.pvd  -- точки разорванных связей, накопленные к кадру N;
    frames.csv                  -- журнал: строка на кадр (время, сводка по Solid и нитям:
                                   центр масс, скорость, KE, прогиб, порванные связи, eps, ...,
                                   и строка monitor.csv на шаге кадра с приставкой mon_);
    breaks.csv                  -- журнал разрывов: строка на связь (когда, какая нить, где);
    columns.txt                 -- что значит каждая колонка обоих журналов.
В FieldData каждого файла -- TIME, STEP и (если есть monitor.csv) величины монитора для шага
кадра: их видит Python Annotation в ParaView. Время в .pvd -- физическое (TIME заголовка
fiber_output или t монитора); если времени в кадрах нет (combined/solid без монитора), в .pvd
пишется номер шага, а не время.

Пример:
    python tools/paraview/sph2vtk.py --run-dir out/1T_scaled_v76.0_r0 --out-dir out/vtk
    python tools/paraview/sph2vtk.py --run-dir ParaviewEx --out-dir ParaviewEx/vtk --jobs 3

Нужен только numpy (pandas -- ускоряет чтение, необязателен); VTK и ParaView не нужны.
"""
from __future__ import annotations

import argparse
import base64
import glob
import os
import re
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor

import numpy as np

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None

FRAME_RE = re.compile(r"(sph|fiber)_output_(\d+)\.txt$")
COMBINED_COLS = ("fiber", "i", "j", "alive")


# ---------------------------------------------------------------------------
# Чтение
# ---------------------------------------------------------------------------

def read_table(path: str) -> tuple[list[str], np.ndarray]:
    """Текстовая таблица: строка заголовка + числа через пробелы."""
    with open(path) as f:
        cols = f.readline().lstrip("# ").split()
    if pd is not None:
        d = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None, dtype=np.float64,
                        engine="c").to_numpy()
    else:
        d = np.loadtxt(path, skiprows=1, ndmin=2)
    if d.size == 0:
        d = np.zeros((0, len(cols)))
    if d.shape[1] != len(cols):
        raise ValueError(f"{path}: в заголовке {len(cols)} колонок, в данных {d.shape[1]}")
    return cols, d


def xyz_of(cols: list[str], d: np.ndarray, names=("x", "y", "z")) -> np.ndarray:
    """(N,3) из колонок; отсутствующие (dim < 3) -- нули."""
    out = np.zeros((len(d), 3))
    for k, n in enumerate(names):
        if n in cols:
            out[:, k] = d[:, cols.index(n)]
    return out


def read_fiber_frame(path: str) -> dict:
    """Формат src/fiber_output.cpp: TIME/STEP/N_PARTICLES/N_BONDS, PARTICLES, BONDS."""
    with open(path) as f:
        lines = f.read().splitlines()
    h = lines[0].split()
    hd = {h[k]: h[k + 1] for k in range(0, len(h) - 1, 2)}
    n, nb = int(hd["N_PARTICLES"]), int(hd["N_BONDS"])
    idx = 1
    while lines[idx].lstrip().startswith("#"):
        idx += 1
    if lines[idx].strip() == "PARTICLES":
        idx += 1
    p = np.loadtxt(lines[idx:idx + n], ndmin=2) if n else np.zeros((0, 9))
    idx += n
    while idx < len(lines) and lines[idx].lstrip().startswith("#"):
        idx += 1
    if idx < len(lines) and lines[idx].strip() == "BONDS":
        idx += 1
    b = np.loadtxt(lines[idx:idx + nb], ndmin=2) if nb else np.zeros((0, 6))
    dim = (p.shape[1] - 3) // 2
    order = np.argsort(p[:, 0].astype(np.int64))      # бусины по id
    p = p[order]
    pos = np.zeros((n, 3)); vel = np.zeros((n, 3))
    pos[:, :dim] = p[:, 3:3 + dim]
    vel[:, :dim] = p[:, 3 + dim:3 + 2 * dim]
    return {
        "time": float(hd["TIME"]), "step": int(float(hd["STEP"])),
        "pos": pos, "vel": vel,
        "object_id": p[:, 1].astype(np.int32), "material_id": p[:, 2].astype(np.int32),
        "ij": b[:, :2].astype(np.int64), "l0": b[:, 2], "l": b[:, 3],
        "l_break": b[:, 4], "elongation": b[:, 5],
    }


def read_monitor(path: str) -> dict[str, np.ndarray] | None:
    """monitor.csv; недописанная последняя строка (прогон идёт) отбрасывается."""
    if not os.path.exists(path):
        return None
    with open(path) as f:
        names = f.readline().strip().split(",")
        rows = [ln.split(",") for ln in f]
    rows = [r for r in rows if len(r) == len(names)]
    d = np.array(rows, dtype=float).reshape(-1, len(names))
    return {n: d[:, k] for k, n in enumerate(names)}


def read_layout(path: str, n_beads: int) -> dict[str, np.ndarray] | None:
    """guo2019_fabric_layout.txt -> исходный номер нити, направление (0 -- x, 1 -- y), зажатие бусины.
    Номер нити берётся отсюда: object_id в fiber_output после разрыва получает оторванный кусок."""
    if not os.path.exists(path):
        return None
    fiber = -np.ones(n_beads, dtype=np.int32)
    direction = -np.ones(n_beads, dtype=np.int32)
    clamped = np.zeros(n_beads, dtype=np.int32)
    with open(path) as f:
        for ln in f:
            p = ln.split()
            if len(p) < 7 or p[0].startswith("#"):
                continue
            a, b = int(p[5]), int(p[6])
            fiber[a:b + 1] = int(p[0])
            direction[a:b + 1] = 0 if p[1] == "x" else 1
            if len(p) >= 9:
                lo, hi = int(p[7]), int(p[8])
                clamped[a:a + lo] = 1
                if hi:
                    clamped[b + 1 - hi:b + 1] = 1
    return {"fiber": fiber, "direction": direction, "clamped": clamped}


# ---------------------------------------------------------------------------
# Запись VTK XML (PolyData, binary base64 + zlib) и PVD
# ---------------------------------------------------------------------------

VTK_TYPES = {np.dtype(np.float32): "Float32", np.dtype(np.float64): "Float64",
             np.dtype(np.int32): "Int32", np.dtype(np.int64): "Int64",
             np.dtype(np.uint8): "UInt8"}


def _encode(a: np.ndarray) -> str:
    """Сжатый блок VTK (vtkZLibDataCompressor, заголовок UInt64), base64."""
    raw = np.ascontiguousarray(a).tobytes()
    if not raw:
        return base64.b64encode(np.array([0, 0, 0], dtype="<u8").tobytes()).decode()
    comp = zlib.compress(raw, 6)
    hdr = np.array([1, len(raw), len(raw), len(comp)], dtype="<u8").tobytes()
    return base64.b64encode(hdr).decode() + base64.b64encode(comp).decode()


def _data_array(name: str, a: np.ndarray, ncomp: int = 1) -> str:
    t = VTK_TYPES[a.dtype]
    nc = f' NumberOfComponents="{ncomp}"' if ncomp > 1 else ""
    return (f'<DataArray type="{t}" Name="{name}"{nc} format="binary">'
            f"{_encode(a)}</DataArray>\n")


def _arr(a: np.ndarray, double: bool) -> np.ndarray:
    a = np.asarray(a)
    if a.dtype.kind == "f":
        return a.astype(np.float64 if double else np.float32)
    if a.dtype.kind in "iub":
        return a.astype(np.int64 if a.dtype.itemsize > 4 else np.int32)
    raise TypeError(a.dtype)


def write_vtp(path: str, points: np.ndarray, *, verts: bool = False,
              lines: np.ndarray | None = None, point_data: dict | None = None,
              cell_data: dict | None = None, field_data: dict | None = None,
              double: bool = False) -> None:
    """PolyData: точки, либо вершины (verts=True), либо отрезки lines (M,2)."""
    n = len(points)
    nl = 0 if lines is None else len(lines)
    nv = n if verts else 0
    out = ['<?xml version="1.0"?>\n',
           '<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" '
           'header_type="UInt64" compressor="vtkZLibDataCompressor">\n<PolyData>\n']
    if field_data:
        out.append("<FieldData>\n")
        for k, v in field_data.items():
            a = np.atleast_1d(np.asarray(v))
            a = a.astype(np.int64) if a.dtype.kind in "iub" else a.astype(np.float64)
            out.append(_data_array(k, a).replace(">", f' NumberOfTuples="{len(a)}">', 1))
        out.append("</FieldData>\n")
    out.append(f'<Piece NumberOfPoints="{n}" NumberOfVerts="{nv}" NumberOfLines="{nl}" '
               f'NumberOfStrips="0" NumberOfPolys="0">\n')
    if point_data:
        out.append("<PointData>\n")
        for k, v in point_data.items():
            v = _arr(v, double)
            out.append(_data_array(k, v, v.shape[1] if v.ndim == 2 else 1))
        out.append("</PointData>\n")
    if cell_data:
        # порядок ячеек VTK: Verts, Lines, ... -- данные связей только для Lines (verts нет)
        assert not verts
        out.append("<CellData>\n")
        for k, v in cell_data.items():
            v = _arr(v, double)
            out.append(_data_array(k, v, v.shape[1] if v.ndim == 2 else 1))
        out.append("</CellData>\n")
    out.append("<Points>\n" + _data_array("Points", _arr(points, double), 3) + "</Points>\n")
    if verts:
        out.append("<Verts>\n" + _data_array("connectivity", np.arange(n, dtype=np.int64))
                   + _data_array("offsets", np.arange(1, n + 1, dtype=np.int64)) + "</Verts>\n")
    if nl:
        out.append("<Lines>\n" + _data_array("connectivity", lines.astype(np.int64).ravel())
                   + _data_array("offsets", np.arange(2, 2 * nl + 1, 2, dtype=np.int64))
                   + "</Lines>\n")
    out.append("</Piece>\n</PolyData>\n</VTKFile>\n")
    tmp = path + ".tmp"
    with open(tmp, "w", newline="\n") as f:
        f.write("".join(out))
    os.replace(tmp, path)


def write_pvd(path: str, entries: list[tuple[float, str]]) -> None:
    with open(path, "w", newline="\n") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
                'byte_order="LittleEndian">\n<Collection>\n')
        for t, name in entries:
            f.write(f'  <DataSet timestep="{t:.12g}" part="0" file="{name}"/>\n')
        f.write("</Collection>\n</VTKFile>\n")


# ---------------------------------------------------------------------------
# Кадр -> файлы
# ---------------------------------------------------------------------------

def solid_point_data(cols: list[str], d: np.ndarray, skip=()) -> dict:
    pd_ = {}
    for k, n in enumerate(cols):
        if n in ("x", "y", "z") or n in skip:
            continue
        v = d[:, k]
        if n in ("N", "Nmat", "alive", "fiber", "k", "i", "j"):
            v = v.astype(np.int32)
        pd_[n] = v
    vel = xyz_of(cols, d, ("vx", "vy", "vz"))
    pd_["velocity"] = vel
    pd_["speed"] = np.linalg.norm(vel, axis=1)
    return pd_


def solid_stats(cols: list[str], d: np.ndarray) -> dict:
    """Сводка по частицам Solid, отдельно для каждого материала Nmat: solid<Nmat>_<величина>."""
    out = {}
    if not len(d) or "Nmat" not in cols:
        return out
    ci = {n: k for k, n in enumerate(cols)}
    for m in np.unique(d[:, ci["Nmat"]]).astype(int):
        s = d[d[:, ci["Nmat"]] == m]
        pre = f"solid{m}_"
        pos, vel = xyz_of(cols, s), xyz_of(cols, s, ("vx", "vy", "vz"))
        w = s[:, ci["M"]] if "M" in ci else np.ones(len(s))
        mt = float(w.sum())
        out[pre + "n"] = len(s)
        out[pre + "mass"] = mt
        for k, a in enumerate("xyz"):
            out[pre + a + "_c"] = float(w @ pos[:, k] / mt)
        for k, a in enumerate("xyz"):
            out[pre + "v" + a + "_c"] = float(w @ vel[:, k] / mt)
        out[pre + "KE"] = float(0.5 * w @ (vel ** 2).sum(axis=1))
        out[pre + "z_min"] = float(pos[:, 2].min())
        out[pre + "z_max"] = float(pos[:, 2].max())
        out[pre + "speed_max"] = float(np.linalg.norm(vel, axis=1).max())
        for n in ("p", "sigma_i", "e_pl", "D"):
            if n in ci:
                if n == "p":
                    out[pre + "p_min"] = float(s[:, ci[n]].min())
                out[pre + n + "_max"] = float(s[:, ci[n]].max())
    return out


def fiber_stats(pos: np.ndarray, vel: np.ndarray, cdata: dict, mass: np.ndarray | None) -> dict:
    """Сводка по нитям кадра (без величин, которым нужен первый кадр: их считает main)."""
    out = {"fiber_n_beads": len(pos), "fiber_n_bonds": len(next(iter(cdata.values())))}
    sp = np.linalg.norm(vel, axis=1)
    out["fiber_speed_max"] = float(sp.max()) if len(sp) else 0.0
    out["fiber_z_min"] = float(pos[:, 2].min()) if len(pos) else 0.0
    if mass is not None:
        out["fiber_KE"] = float(0.5 * mass @ sp ** 2)
    if "eps" in cdata and len(cdata["eps"]):
        e = cdata["eps"]
        out["fiber_eps_max"] = float(e.max())
        out["fiber_n_tensile"] = int((e > 0).sum())
    if "eps_ratio" in cdata and len(cdata["eps_ratio"]):
        r = cdata["eps_ratio"]
        out["fiber_eps_ratio_max"] = float(r.max())
        out["fiber_n_ratio_gt_0.9"] = int((r > 0.9).sum())
    if "tension" in cdata and len(cdata["tension"]):
        out["fiber_tension_max"] = float(cdata["tension"].max())
        out["fiber_tension_sum"] = float(cdata["tension"].sum())
    return out


def convert_frame(job: dict) -> dict:
    """Один кадр: пишет .vtp и возвращает время, шаг, сводку и связи (для журнала разрывов)."""
    fid, mode, run_dir, out_dir = job["fid"], job["mode"], job["run_dir"], job["out_dir"]
    double, mon, layout_path = job["double"], job["monitor"], job["layout"]
    step = int(fid)
    time = None
    sph_path = os.path.join(run_dir, f"sph_output_{fid}.txt")
    summary = []
    stats: dict = {}

    fiber = None
    if mode == "pair":
        fr = read_fiber_frame(os.path.join(run_dir, f"fiber_output_{fid}.txt"))
        step, time = fr["step"], fr["time"]
        ij = fr["ij"]
        eb = fr["l_break"] / fr["l0"] - 1.0
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(eb > 0, fr["elongation"] / eb, 0.0)
        pdata = {"velocity": fr["vel"], "speed": np.linalg.norm(fr["vel"], axis=1),
                 "object_id": fr["object_id"], "material_id": fr["material_id"]}
        cdata = {"eps": fr["elongation"], "eps_ratio": ratio, "l0": fr["l0"], "l": fr["l"],
                 "l_break": fr["l_break"], "object_id": fr["object_id"][ij[:, 0]]}
        lay = read_layout(layout_path, len(fr["pos"])) if layout_path else None
        if lay is not None:
            pdata.update(lay)
            cdata["direction"] = lay["direction"][ij[:, 0]]
            cdata["fiber"] = lay["fiber"][ij[:, 0]]
        fiber = (fr["pos"], ij, pdata, cdata)
        stats.update(fiber_stats(fr["pos"], fr["vel"], cdata, None))
        summary.append(f"бусин {len(fr['pos'])}, связей {len(ij)}")

    solid = None
    if os.path.exists(sph_path):
        cols, d = read_table(sph_path)
        if mode == "combined":
            ci = {n: k for k, n in enumerate(cols)}
            alive = d[:, ci["alive"]] > 0.5
            fmat = np.unique(d[alive, ci["Nmat"]])            # материалы, у которых есть связи
            isf = np.isin(d[:, ci["Nmat"]], fmat)
            gidx = np.flatnonzero(isf)
            loc = -np.ones(len(d), dtype=np.int64); loc[gidx] = np.arange(len(gidx))
            bead = d[isf]
            brow = np.flatnonzero(alive)                       # строки, несущие связь
            gi = d[brow, ci["i"]].astype(np.int64); gj = d[brow, ci["j"]].astype(np.int64)
            ok = (gi >= 0) & (gi < len(d)) & (gj >= 0) & (gj < len(d))
            ok &= (loc[np.clip(gi, 0, len(d) - 1)] >= 0) & (loc[np.clip(gj, 0, len(d) - 1)] >= 0)
            if not ok.all():
                print(f"  кадр {fid}: {int((~ok).sum())} связей с i/j вне бусин -- пропущены",
                      file=sys.stderr)
            brow, gi, gj = brow[ok], gi[ok], gj[ok]
            ij = np.stack([loc[gi], loc[gj]], axis=1)
            pdata = solid_point_data(cols, bead, skip=("k", "i", "j", "eps", "eps_pl", "tension"))
            # у последней бусины нити колонка fiber = 0 (связи нет): берём номер нити соседа
            fib = pdata["fiber"].copy()
            fib[ij[:, 1]] = fib[ij[:, 0]]
            pdata["fiber"] = fib
            cdata = {n: d[brow, ci[n]] for n in ("eps", "eps_pl", "tension") if n in ci}
            cdata["k"] = d[brow, ci["k"]].astype(np.int32)
            if job["eps_break"] > 0 and "eps" in cdata:
                cdata["eps_ratio"] = cdata["eps"] / job["eps_break"]
            cdata["fiber"] = fib[ij[:, 0]]
            fiber = (xyz_of(cols, bead), ij, pdata, cdata)
            stats.update(fiber_stats(fiber[0], pdata["velocity"], cdata,
                                     bead[:, ci["M"]] if "M" in ci else None))
            d = d[~isf]
            summary.append(f"бусин {len(bead)}, связей {len(ij)}")
        skip = ("T", "k", "fiber", "i", "j", "eps", "eps_pl", "tension", "alive") \
            if mode == "combined" else ()
        solid = (xyz_of(cols, d), solid_point_data(cols, d, skip=skip))
        stats.update(solid_stats(cols, d))
        summary.append(f"частиц Solid {len(d)}")

    fdata = {"STEP": step}
    if mon is not None:
        r = np.flatnonzero(mon["step"].astype(np.int64) == step)
        if len(r):
            for k, v in mon.items():
                fdata[k] = v[r[0]]
            if "t" in mon:
                time = float(mon["t"][r[0]])
    pvd_t = time if time is not None else float(step)
    fdata["TIME"] = pvd_t

    if solid is not None:
        write_vtp(os.path.join(out_dir, f"solid_{fid}.vtp"), solid[0], verts=True,
                  point_data=solid[1], field_data=fdata, double=double)
    if fiber is not None:
        write_vtp(os.path.join(out_dir, f"fibers_{fid}.vtp"), fiber[0], lines=fiber[1],
                  point_data=fiber[2], cell_data=fiber[3], field_data=fdata, double=double)
    res = {"fid": fid, "time": pvd_t, "step": step, "summary": ", ".join(summary),
           "stats": stats, "monitor_row": {k: v for k, v in fdata.items()
                                           if k not in ("STEP", "TIME")}}
    if fiber is not None:
        pos, ij, pdata, cdata = fiber
        res["pos"] = pos.astype(np.float32)
        res["ij"] = ij
        res["bond"] = {k: cdata[k] for k in ("eps", "eps_ratio", "tension") if k in cdata}
        bead_fiber = pdata.get("fiber", pdata.get("object_id"))   # исходная нить, иначе кусок
        if bead_fiber is not None:
            res["bead_fiber"] = np.asarray(bead_fiber)
        if "direction" in pdata:
            res["bead_dir"] = np.asarray(pdata["direction"])
    return res


# ---------------------------------------------------------------------------
# Журналы для графиков: frames.csv, breaks.csv, breaks_<N>.vtp, columns.txt
# ---------------------------------------------------------------------------

COLUMN_DOC = {
    "frame": "порядковый номер кадра в серии (с 0)",
    "step": "номер шага расчёта",
    "time": "время кадра, с (если времени в кадрах нет -- номер шага)",
    "time_us": "время кадра, мкс",
    "n": "число частиц",
    "mass": "масса, кг",
    "x_c": "центр масс по x, м", "y_c": "центр масс по y, м", "z_c": "центр масс по z, м",
    "vx_c": "скорость центра масс по x, м/с", "vy_c": "скорость центра масс по y, м/с",
    "vz_c": "скорость центра масс по z, м/с",
    "KE": "кинетическая энергия, Дж",
    "z_min": "наименьшая координата z, м", "z_max": "наибольшая координата z, м",
    "speed_max": "наибольшая скорость частицы, м/с",
    "p_min": "наименьшее давление, Па", "p_max": "наибольшее давление, Па",
    "sigma_i_max": "наибольшая интенсивность напряжений, Па",
    "e_pl_max": "наибольшая интенсивность пластических деформаций",
    "D_max": "наибольшая повреждённость",
    "fiber_n_beads": "число бусин",
    "fiber_n_bonds": "число целых связей",
    "fiber_n_broken": "порвано связей с первого кадра серии (накопительно)",
    "fiber_n_broken_x": "из них в нитях вдоль x (по раскладке)",
    "fiber_n_broken_y": "из них в нитях вдоль y (по раскладке)",
    "fiber_w_max": "наибольший прогиб бусины вниз от её положения в первом кадре, м",
    "fiber_z_min": "наименьшая координата z бусин, м",
    "fiber_speed_max": "наибольшая скорость бусины, м/с",
    "fiber_KE": "кинетическая энергия бусин, Дж (если в кадре есть масса)",
    "fiber_eps_max": "наибольшее удлинение связи eps",
    "fiber_n_tensile": "число растянутых связей (eps > 0)",
    "fiber_eps_ratio_max": "наибольшее eps/eps_разр",
    "fiber_n_ratio_gt_0.9": "число связей с eps/eps_разр > 0,9",
    "fiber_tension_max": "наибольшее натяжение связи, Н",
    "fiber_tension_sum": "сумма натяжений связей, Н",
}

BREAK_DOC = {
    "frame": "кадр, в котором связи уже нет",
    "step": "шаг этого кадра",
    "time": "время этого кадра, с (или шаг); разрыв -- между time_before и time",
    "time_before": "время предыдущего кадра, где связь ещё была",
    "time_us": "time, мкс",
    "time_before_us": "time_before, мкс",
    "bead_i": "бусина-начало связи (id с 0, как в fibers_<N>.vtp)",
    "bead_j": "бусина-конец связи",
    "fiber": "номер нити (object_id или fiber)",
    "direction": "направление нити: 0 -- x, 1 -- y, -1 -- неизвестно",
    "eps_last": "eps связи в последнем кадре до разрыва",
    "eps_ratio_last": "eps/eps_разр в последнем кадре до разрыва",
    "tension_last": "натяжение в последнем кадре до разрыва, Н",
    "x": "середина связи в кадре разрыва, м",
    "y": "середина связи в кадре разрыва, м",
    "z": "середина связи в кадре разрыва, м",
}


def _doc(name: str) -> str:
    if name in COLUMN_DOC:
        return COLUMN_DOC[name]
    if name.startswith("mon_"):
        return f"monitor.csv: {name[4:]} на шаге кадра"
    m = re.match(r"solid(\d+)_(.+)$", name)
    if m and m.group(2) in COLUMN_DOC:
        return f"Solid, материал Nmat={m.group(1)}: {COLUMN_DOC[m.group(2)]}"
    return ""


def _fmt(v) -> str:
    if v is None:
        return "nan"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    v = float(v)
    return f"{v:.9g}" if np.isfinite(v) else "nan"


def write_csv(path: str, cols: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(_fmt(r.get(c)) for c in cols) + "\n")


def write_logs(out_dir: str, res: list[dict], timed: bool, double: bool) -> None:
    """frames.csv (строка на кадр), breaks.csv (строка на разрыв), breaks_<N>.vtp + breaks.pvd
    (точки разрывов, накопленные к кадру) и columns.txt (описание колонок).
    Разрыв -- связь (i, j), которая была в предыдущем кадре серии и пропала в текущем."""
    rows, breaks = [], []
    z0 = res[0]["pos"][:, 2] if "pos" in res[0] else None
    have_dir = any("bead_dir" in r for r in res)
    n_br = n_br_x = n_br_y = 0
    prev = None
    for k, r in enumerate(res):
        row = {"frame": k, "step": r["step"], "time": r["time"]}
        if timed:
            row["time_us"] = r["time"] * 1e6
        if "pos" in r:
            nb = len(r["pos"])
            key = r["ij"][:, 0] * nb + r["ij"][:, 1]
            if prev is not None and len(prev["pos"]) == nb:
                pkey = prev["ij"][:, 0] * nb + prev["ij"][:, 1]
                for b in np.flatnonzero(~np.isin(pkey, key)):
                    i, j = (int(v) for v in prev["ij"][b])
                    mid = 0.5 * (r["pos"][i].astype(float) + r["pos"][j].astype(float))
                    ev = {"frame": k, "step": r["step"], "time": r["time"],
                          "time_before": prev["time"], "bead_i": i, "bead_j": j,
                          "fiber": int(r["bead_fiber"][i]) if "bead_fiber" in r else -1,
                          "direction": int(r["bead_dir"][i]) if "bead_dir" in r else -1,
                          "x": mid[0], "y": mid[1], "z": mid[2]}
                    if timed:
                        ev["time_us"], ev["time_before_us"] = r["time"] * 1e6, prev["time"] * 1e6
                    for n, a in prev["bond"].items():
                        ev[n + "_last"] = float(a[b])
                    breaks.append(ev)
                    n_br += 1
                    n_br_x += ev["direction"] == 0
                    n_br_y += ev["direction"] == 1
            row["fiber_n_broken"] = n_br
            if have_dir:
                row["fiber_n_broken_x"], row["fiber_n_broken_y"] = n_br_x, n_br_y
            if z0 is not None and len(z0) == nb:
                row["fiber_w_max"] = max(0.0, float((z0 - r["pos"][:, 2]).max()))
            prev = r
        row.update(r["stats"])
        row.update({"mon_" + n: v for n, v in r["monitor_row"].items()})
        rows.append(row)

    cols: list[str] = []
    for row in rows:
        cols += [c for c in row if c not in cols]
    write_csv(os.path.join(out_dir, "frames.csv"), cols, rows)
    print(f"frames.csv: {len(rows)} строк, {len(cols)} колонок")

    bcols = ["frame", "step", "time", "time_before"] \
        + (["time_us", "time_before_us"] if timed else []) \
        + ["bead_i", "bead_j", "fiber", "direction", "eps_last", "eps_ratio_last", "tension_last",
           "x", "y", "z"]
    present = set().union(*(e.keys() for e in breaks)) if breaks else None
    if present is not None:
        bcols = [c for c in bcols if c in present]
    write_csv(os.path.join(out_dir, "breaks.csv"), bcols, breaks)
    print(f"breaks.csv: {len(breaks)} разрывов")

    if "pos" in res[0]:
        ent = []
        for k, r in enumerate(res):
            evs = [e for e in breaks if e["frame"] <= k]
            pts = np.array([[e["x"], e["y"], e["z"]] for e in evs], dtype=float).reshape(-1, 3)
            pdat = {"break_time": np.array([e["time"] for e in evs], dtype=float),
                    "break_frame": np.array([e["frame"] for e in evs], dtype=np.int32),
                    "fiber": np.array([e["fiber"] for e in evs], dtype=np.int32),
                    "direction": np.array([e["direction"] for e in evs], dtype=np.int32)}
            if timed:
                pdat["break_time_us"] = pdat["break_time"] * 1e6
            name = f"breaks_{r['fid']}.vtp"
            write_vtp(os.path.join(out_dir, name), pts, verts=True, point_data=pdat,
                      field_data={"TIME": r["time"], "STEP": r["step"]}, double=double)
            ent.append((r["time"], name))
        write_pvd(os.path.join(out_dir, "breaks.pvd"), ent)
        print(f"breaks.pvd: {len(ent)} кадров")

    with open(os.path.join(out_dir, "columns.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write("frames.csv -- строка на кадр; разделитель запятая, имена колонок в первой строке.\n")
        for c in cols:
            f.write(f"  {c:28s} {_doc(c)}\n")
        f.write("\nbreaks.csv -- строка на разорванную связь (разрыв найден между соседними кадрами"
                " серии: при --every > 1 время разрыва известно грубее).\n")
        for c in bcols:
            f.write(f"  {c:28s} {BREAK_DOC.get(c, '')}\n")
        f.write("\nbreaks_<N>.vtp, breaks.pvd -- точки разрывов, случившихся к кадру N;\n"
                "  PointData: break_time, break_time_us, break_frame, fiber, direction.\n")
        if not timed:
            f.write("\nВремени в кадрах нет: колонки time -- номер шага.\n")


# ---------------------------------------------------------------------------

def discover(run_dir: str) -> tuple[str, list[str]]:
    ids = {"sph": set(), "fiber": set()}
    for f in glob.glob(os.path.join(run_dir, "*_output_*.txt")):
        m = FRAME_RE.search(os.path.basename(f))
        if m:
            ids[m.group(1)].add(m.group(2))
    if ids["fiber"]:
        return "pair", sorted(ids["fiber"], key=int)
    if not ids["sph"]:
        return "", []
    fids = sorted(ids["sph"], key=int)
    with open(os.path.join(run_dir, f"sph_output_{fids[0]}.txt")) as f:
        cols = f.readline().split()
    return ("combined" if all(c in cols for c in COMBINED_COLS) else "solid"), fids


def main() -> None:
    for st in (sys.stdout, sys.stderr):   # Windows: иначе вывод в pipe идёт в cp1251
        if hasattr(st, "reconfigure"):
            st.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="каталог с кадрами *_output_<N>.txt")
    ap.add_argument("--out-dir", required=True, help="куда писать .vtp и .pvd")
    ap.add_argument("--every", type=int, default=1, help="брать каждый k-й кадр")
    ap.add_argument("--double", action="store_true", help="Float64 вместо Float32")
    ap.add_argument("--eps-break", type=float, default=0.0,
                    help="для кадров combined: деформация разрыва, даёт eps_ratio = eps/eps_break "
                         "(в pair берётся из l_break каждой связи)")
    ap.add_argument("--jobs", type=int, default=max(1, min(8, (os.cpu_count() or 2) // 2)),
                    help="процессов для перевода кадров")
    args = ap.parse_args()

    mode, fids = discover(args.run_dir)
    if not fids:
        sys.exit(f"{args.run_dir}: нет кадров sph_output_*/fiber_output_*")
    fids = fids[::max(1, args.every)]
    os.makedirs(args.out_dir, exist_ok=True)
    mon = read_monitor(os.path.join(args.run_dir, "monitor.csv"))
    lay = os.path.join(args.run_dir, "guo2019_fabric_layout.txt")
    print(f"{args.run_dir}: вид кадров '{mode}', кадров {len(fids)}; "
          f"monitor.csv {'есть' if mon is not None else 'нет'}; "
          f"раскладка нитей {'есть' if os.path.exists(lay) else 'нет'}")
    jobs = [{"fid": f, "mode": mode, "run_dir": args.run_dir, "out_dir": args.out_dir,
             "double": args.double, "monitor": mon, "eps_break": args.eps_break,
             "layout": lay if os.path.exists(lay) else None} for f in fids]

    if args.jobs > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            res = list(ex.map(convert_frame, jobs))
    else:
        res = [convert_frame(j) for j in jobs]

    timed = mode == "pair" or mon is not None
    for r in res[:1] + res[-1:]:
        print(f"  кадр {r['fid']}: шаг {r['step']}, "
              f"{'t = %.4g с' % r['time'] if timed else 'время = шаг'}; {r['summary']}")
    for kind in ("solid", "fibers"):
        ent = [(r["time"], f"{kind}_{r['fid']}.vtp") for r in res
               if os.path.exists(os.path.join(args.out_dir, f"{kind}_{r['fid']}.vtp"))]
        if ent:
            write_pvd(os.path.join(args.out_dir, f"{kind}.pvd"), ent)
            print(f"{kind}.pvd: {len(ent)} кадров")
    write_logs(args.out_dir, res, timed, args.double)
    if not timed:
        print("Времени в кадрах нет: в .pvd и журналах вместо времени номер шага.")


if __name__ == "__main__":
    main()
