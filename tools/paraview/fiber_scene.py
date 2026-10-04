"""Сцена ParaView для кадров, переведённых tools/paraview/sph2vtk.py (solid.pvd + fibers.pvd).

Запуск -- через pvpython из комплекта ParaView (проверено на 6.2):

    pvpython tools/paraview/fiber_scene.py --vtk-dir out/vtk --state out/vtk/scene.pvsm
    pvpython tools/paraview/fiber_scene.py --vtk-dir out/vtk --out anim.mp4 --zoom 9e-3
    pvpython tools/paraview/fiber_scene.py --vtk-dir out/vtk --snapshots 0,-1 --out snap.png

Что строится (как в experiments/2026-10-01-guo2019-1t-scaled/scripts/animate_1t.py):
  * слева -- 3D: нити трубками (Tube по отрезкам-связям), цвет -- eps_ratio = eps/eps_разр
    (шкала 0..1); если eps_ratio нет -- tension, иначе eps. Порванные связи в кадре отсутствуют.
    Частицы Solid -- Point Gaussian, серые или по массиву --solid-array;
  * справа (--section yz|xz, по умолчанию yz) -- сечение: связи и частицы из полосы
    |x - x_c| < band (для yz), вид сбоку в параллельной проекции;
  * текст -- время и величины монитора из FieldData кадра (v_z, w_f, порванные нити).
--state сохраняет состояние .pvsm: его можно открыть в GUI (File > Load State) и крутить руками.
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
from paraview import servermanager as sm
from paraview.simple import (ColorBy, ExtractSurface, CreateRenderView, CreateLayout, ExtractCellsByRegion,
                             GetAnimationScene, GetColorTransferFunction, GetOpacityTransferFunction, GetScalarBar, PVDReader,
                             PythonAnnotation, SaveAnimation, SaveScreenshot, SaveState, Show,
                             Tube, AssignViewToLayout, Render)
from vtkmodules.util.numpy_support import vtk_to_numpy


def fetch(src, t):
    src.UpdatePipeline(t)
    return sm.Fetch(src)


def cell_arrays(o) -> list[str]:
    cd = o.GetCellData()
    return [cd.GetArrayName(k) for k in range(cd.GetNumberOfArrays())]


def field_arrays(o) -> list[str]:
    fd = o.GetFieldData()
    return [fd.GetArrayName(k) for k in range(fd.GetNumberOfArrays())]


def median_bond_length(o) -> float:
    pts = vtk_to_numpy(o.GetPoints().GetData())
    conn = vtk_to_numpy(o.GetLines().GetConnectivityArray()).reshape(-1, 2)
    return float(np.median(np.linalg.norm(pts[conn[:, 0]] - pts[conn[:, 1]], axis=1)))


def annotation_expr(fields: list[str]) -> str:
    """Выражение PythonAnnotation по тем полям FieldData, что есть в кадре.
    Поле из одного значения приходит в выражение то массивом, то числом -- отсюда v(...);
    lambda не годится: внутри неё не видны имена полей."""
    def v(n: str) -> str:
        return f"float(getattr({n}, 'flat', [{n}])[0])"
    timed = "t" in fields or "dt" in fields
    head = [f"'t = %.2f мкс' % ({v('TIME')}*1e6)" if timed else f"'шаг %d' % {v('STEP')}"]
    if timed:
        head.append(f"'шаг %d' % {v('STEP')}")
    parts = []
    if "v_z" in fields:
        parts.append(f"'v_z = %.1f м/с' % {v('v_z')}")
    if "w_f" in fields:
        parts.append(f"'прогиб w_f = %.2f мм' % ({v('w_f')}*1e3)")
    if "n_broken_v" in fields and "n_broken_h" in fields:
        parts.append(f"'порвано нитей верт./гориз. = %d/%d' % ({v('n_broken_v')}, {v('n_broken_h')})")
    def line(xs: list[str]) -> str:
        return "'   '.join([" + ", ".join(xs) + "])"
    return line(head) + (" + chr(10) + " + line(parts) if parts else "")


FONT_CANDIDATES = (
    r"C:/Windows/Fonts/arial.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
    "/Library/Fonts/Arial Unicode.ttf",
)


def find_font(user: str) -> str | None:
    """Шрифт с кириллицей: встроенный в ParaView Arial её не содержит (русские буквы пропадают)."""
    for f in ((user,) if user else FONT_CANDIDATES):
        if f and os.path.exists(f):
            return f
    return None


def use_font(obj, font: str | None, *prefixes: str) -> None:
    """FontFamily='File' + FontFile для свойств <prefix>FontFamily/<prefix>FontFile."""
    if not font:
        return
    for p in prefixes:
        setattr(obj, p + "FontFamily", "File")
        setattr(obj, p + "FontFile", font)


def set_camera(view, focal, dist, elev, azim):
    """Камера как view_init(elev, azim) у matplotlib: ось z вверх."""
    e, a = math.radians(elev), math.radians(azim)
    d = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
    view.CameraFocalPoint = list(focal)
    view.CameraPosition = list(np.asarray(focal) + dist * d)
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraViewAngle = 30.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vtk-dir", required=True, help="каталог с solid.pvd / fibers.pvd")
    ap.add_argument("--out", default=None,
                    help="анимация (.mp4/.avi/.ogv) или серия .png; со --snapshots -- база имени снимков")
    ap.add_argument("--snapshots", default=None,
                    help="номера кадров через запятую (отрицательные -- с конца): только снимки")
    ap.add_argument("--state", default=None, help="сохранить состояние .pvsm для GUI")
    ap.add_argument("--zoom", type=float, default=0.0,
                    help="полуразмер видимого окна вокруг центра ткани, м (0 -- вся ткань)")
    ap.add_argument("--section", choices=["yz", "xz", "none"], default="yz")
    ap.add_argument("--band", type=float, default=1.5e-3, help="полуширина полосы сечения, м")
    ap.add_argument("--tube-radius", type=float, default=0.0,
                    help="радиус трубок нитей, м (0 -- 0,3 медианной длины связи)")
    ap.add_argument("--solid-array", default="", help="массив PointData для цвета Solid (по умолчанию серый)")
    ap.add_argument("--color-range", default="0,1", help="шкала eps_ratio 'lo,hi'")
    ap.add_argument("--elev", type=float, default=24.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--size", default="1600x800", help="размер кадра WxH")
    ap.add_argument("--fps", type=int, default=6)
    ap.add_argument("--font", default="", help="TTF-шрифт с кириллицей (по умолчанию ищется Arial/DejaVu)")
    args = ap.parse_args()
    font = find_font(args.font)
    if font is None:
        print("шрифт с кириллицей не найден (--font): русские подписи могут пропасть", file=sys.stderr)

    fib_pvd = os.path.join(args.vtk_dir, "fibers.pvd")
    sol_pvd = os.path.join(args.vtk_dir, "solid.pvd")
    have_f, have_s = os.path.exists(fib_pvd), os.path.exists(sol_pvd)
    if not (have_f or have_s):
        sys.exit(f"{args.vtk_dir}: нет fibers.pvd / solid.pvd (сначала tools/paraview/sph2vtk.py)")
    W, H = (int(v) for v in args.size.lower().split("x"))

    fib = PVDReader(registrationName="fibers", FileName=os.path.abspath(fib_pvd)) if have_f else None
    sol = PVDReader(registrationName="solid", FileName=os.path.abspath(sol_pvd)) if have_s else None
    base = fib or sol
    times = list(base.TimestepValues) or [0.0]
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()

    o0 = fetch(base, times[0])
    b = o0.GetBounds()
    ctr = np.array([0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5])])
    half = args.zoom if args.zoom > 0 else 0.5 * max(b[1] - b[0], b[3] - b[2])

    # ---- 3D-вид ----
    layout = CreateLayout("scene")
    v3 = CreateRenderView()
    AssignViewToLayout(view=v3, layout=layout, hint=0)
    v3.OrientationAxesVisibility = 1
    v3.Background = [1, 1, 1]
    v3.UseColorPaletteForBackground = 0

    carr, r_tube, c_range = None, 0.0, None
    if fib is not None:
        arrs = cell_arrays(o0)
        carr = next((a for a in ("eps_ratio", "tension", "eps") if a in arrs), None)
        r_tube = args.tube_radius or 0.3 * median_bond_length(o0)
        tube = Tube(registrationName="fiber_tubes", Input=fib)
        tube.Radius = r_tube
        tube.NumberofSides = 8
        tube.Capping = 0
        td = Show(tube, v3)
        if carr:
            ColorBy(td, ("CELLS", carr))
            lut = GetColorTransferFunction(carr)
            lut.ApplyPreset("Turbo", True)
            lut.AutomaticRescaleRangeMode = "Never"   # иначе шкала подстраивается под кадр
            if carr == "eps_ratio":
                lo, hi = (float(x) for x in args.color_range.split(","))
            else:
                rg = fetch(fib, times[-1]).GetCellData().GetArray(carr).GetRange()
                lo, hi = (0.0 if carr == "tension" else rg[0]), rg[1]
            c_range = (lo, hi)
            td.SetScalarBarVisibility(v3, True)
            bar = GetScalarBar(lut, v3)
            bar.Title = {"eps_ratio": "ε / ε разр", "tension": "натяжение, Н", "eps": "ε"}[carr]
            bar.ComponentTitle = ""
            bar.TitleColor = bar.LabelColor = [0, 0, 0]
            bar.TitleFontSize, bar.LabelFontSize = 18, 14
            use_font(bar, font, "Title", "Label")
            bar.Orientation = "Horizontal"
            bar.WindowLocation = "Lower Right Corner"
            bar.ScalarBarLength = 0.4
            bar.AddRangeLabels = 0
            bar.DrawBackground = 1
            bar.BackgroundColor = [1, 1, 1, 0.85]
        print(f"нити: трубки r = {r_tube * 1e3:.3f} мм, цвет -- {carr or 'однотонный'}")

    def show_solid(src, view):
        sd = Show(src, view)
        sd.SetRepresentationType("Point Gaussian")
        so = fetch(src, times[0])
        hs = so.GetPointData().GetArray("hsml")
        sd.GaussianRadius = 0.4 * (hs.GetRange()[0] if hs is not None else 1e-3)
        sd.ShaderPreset = "Sphere"
        if args.solid_array:
            ColorBy(sd, ("POINTS", args.solid_array))
            GetColorTransferFunction(args.solid_array).ApplyPreset("Cool to Warm", True)
            sd.SetScalarBarVisibility(view, True)
        else:
            sd.ColorArrayName = ["POINTS", ""]
            sd.AmbientColor = sd.DiffuseColor = [0.55, 0.55, 0.55]
        return sd

    if sol is not None:
        show_solid(sol, v3)

    ann = PythonAnnotation(registrationName="caption", Input=base)
    ann.ArrayAssociation = "Field Data"
    ann.Expression = annotation_expr(field_arrays(o0))
    ad = Show(ann, v3)
    ad.FontSize = 18
    ad.Color = [0, 0, 0]
    ad.WindowLocation = "Upper Left Corner"
    ad.Justification = "Left"
    use_font(ad, font, "")

    # ход по z за весь расчёт (первый и последний кадры): камеры смотрят на его середину
    zs = []
    for t in (times[0], times[-1]):
        for src in (fib, sol):
            if src is not None:
                bb = fetch(src, t).GetBounds()
                zs += [bb[4], bb[5]]
    zc = 0.5 * (min(zs) + max(zs))
    focal = np.array([ctr[0], ctr[1], zc])
    reach = max(half, 0.5 * (max(zs) - min(zs)))

    # камеры ставятся после первой отрисовки: при ней ParaView сбрасывает камеру по данным
    cameras = [lambda: set_camera(v3, focal, reach / math.tan(math.radians(15.0)) * 2.0,
                                  args.elev, args.azim)]

    # ---- сечение ----
    vs = None
    if args.section != "none":
        ic = 0 if args.section == "yz" else 1       # поперечная ось полосы
        big = 10.0
        pos = [-big / 2] * 3; ln = [big] * 3
        pos[ic], ln[ic] = ctr[ic] - args.band, 2.0 * args.band
        vs = CreateRenderView()
        AssignViewToLayout(view=vs, layout=layout, hint=2)
        vs.Background = [1, 1, 1]
        vs.UseColorPaletteForBackground = 0
        vs.OrientationAxesVisibility = 0
        vs.CameraParallelProjection = 1
        vs.AxesGrid.Visibility = 1
        vs.AxesGrid.XTitle, vs.AxesGrid.YTitle, vs.AxesGrid.ZTitle = "x, м", "y, м", "z, м"
        for p in ("XTitleColor", "YTitleColor", "ZTitleColor", "XLabelColor", "YLabelColor",
                  "ZLabelColor", "GridColor"):
            setattr(vs.AxesGrid, p, [0, 0, 0])
        # подписи только у осей плоскости сечения: грани min/max этой плоскости, оси min
        vs.AxesGrid.FacesToRender = 0x01 | 0x08 if args.section == "yz" else 0x02 | 0x10
        vs.AxesGrid.AxesToLabel = 0x02 | 0x04 if args.section == "yz" else 0x01 | 0x04
        use_font(vs.AxesGrid, font, *(a + k for a in "XYZ" for k in ("Title", "Label")))
        for src, name in ((fib, "fiber_section"), (sol, "solid_section")):
            if src is None:
                continue
            ex = ExtractCellsByRegion(registrationName=name, Input=src)
            ex.IntersectWith = "Box"
            ex.IntersectWith.UseReferenceBounds = 0
            ex.IntersectWith.Position = pos
            ex.IntersectWith.Length = ln
            ex.Extractintersected = 1
            if src is fib:
                surf = ExtractSurface(registrationName="fiber_section_poly", Input=ex)
                st = Tube(registrationName="fiber_section_tubes", Input=surf)
                st.Radius = r_tube
                st.NumberofSides = 8
                sd = Show(st, vs)
                if carr:
                    ColorBy(sd, ("CELLS", carr))
                pts = Show(ex, vs)  # бусины -- точками поверх отрезков
                pts.SetRepresentationType("Points")
                pts.ColorArrayName = ["POINTS", ""]
                pts.AmbientColor = pts.DiffuseColor = [0, 0, 0]
                pts.PointSize = 3
            else:
                show_solid(ex, vs)
        # смотрим вдоль поперечной оси: для yz камера на +x, y вправо; для xz -- на -y, x вправо
        d = np.zeros(3); d[ic] = 1.0 if args.section == "yz" else -1.0
        # окно сечения: по горизонтали ±half, по z -- весь ход по кадрам
        scale = 1.25 * max(half * H / (W / 2), 0.5 * (max(zs) - min(zs)))

        def section_camera():
            vs.CameraFocalPoint = [ctr[0], ctr[1], zc]
            vs.CameraPosition = [ctr[0] + d[0], ctr[1] + d[1], zc]
            vs.CameraViewUp = [0, 0, 1]
            vs.CameraParallelScale = scale
        cameras.append(section_camera)
        # сетка осей -- по видимому окну, иначе ребро с подписями z уходит за кадр (ткань шире окна)
        gb = [ctr[0] - half, ctr[0] + half, ctr[1] - half, ctr[1] + half, min(zs), max(zs)]
        gb[2 * ic], gb[2 * ic + 1] = ctr[ic] - args.band, ctr[ic] + args.band
        vs.AxesGrid.UseCustomBounds = 1
        vs.AxesGrid.CustomBounds = gb
        tt = PythonAnnotation(registrationName="section_title", Input=base)
        tt.ArrayAssociation = "Field Data"
        tt.Expression = (f"'сечение {args.section}: полоса {'xyz'[ic]} = "
                         f"{ctr[ic]*1e3:.2f} ± {args.band*1e3:.2f} мм'")
        td2 = Show(tt, vs)
        td2.FontSize = 16
        td2.Color = [0, 0, 0]
        td2.WindowLocation = "Upper Center"
        use_font(td2, font, "")

    if c_range is not None:
        # после всех ColorBy: ColorBy для сечения успевает перенастроить диапазон шкалы
        lut = GetColorTransferFunction(carr)
        lut.RescaleTransferFunction(*c_range)
        GetOpacityTransferFunction(carr).RescaleTransferFunction(*c_range)
    layout.SetSize(W, H)
    for v in (v3, vs):
        if v is not None:
            Render(v)
    for c in cameras:
        c()

    if args.state:
        SaveState(os.path.abspath(args.state))
        print(f"состояние: {args.state}")

    if args.snapshots:
        stem, ext = os.path.splitext(args.out or os.path.join(args.vtk_dir, "snap.png"))
        for s in args.snapshots.split(","):
            k = int(s)
            k = k + len(times) if k < 0 else k
            scene.AnimationTime = times[k]
            p = f"{stem}_frame{k}{ext or '.png'}"
            SaveScreenshot(p, layout, ImageResolution=[W, H])
            print(f"снимок: {p}")
    elif args.out:
        SaveAnimation(args.out, layout, ImageResolution=[W, H], FrameRate=max(1, round(args.fps)))
        print(f"анимация: {args.out} ({len(times)} кадров)")


if __name__ == "__main__":
    main()
