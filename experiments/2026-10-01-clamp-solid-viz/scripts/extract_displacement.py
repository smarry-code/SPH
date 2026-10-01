#!/usr/bin/env python3
"""Смещение конца консоли (ES2) и середины шарнирной балки (ES3) по кадрам --frames-dir.

Кадры (sph_output_*.txt) в репозиторий не кладутся, поэтому скрипт запускается один раз
на машине, где они есть, и пишет data/<прогон>/disp.csv (t, z) -- по одной точке на кадр.
Номер файла кадра = номер шага; время шага берётся из reaction.csv того же прогона.

    python scripts/extract_displacement.py <каталог --frames-dir>

Индекс частицы N-1 = (ix*ny + iy)*nz + iz (tests/clamp_e2e_solid.hpp), число слоёв и сечение
берутся из clamp_case.hpp::solid: ES2: 5 слоёв заделки + 10*4 = 40 слоёв, сечение 4x4, конец --
последний слой; ES3: 10*5 = 50 слоёв, сечение 5x5, середина -- слои 24 и 25. Начальное среднее z
сечения равно нулю (сечение симметрично), поэтому смещение = среднее z.
"""
import glob
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

RUNS = {
    # прогон: (nx, ny, nz, слои для усреднения)
    "ES2_pc_ES2": (45, 4, 4, [44]),
    "ES3_pc_ES3": (50, 5, 5, [24, 25]),
}


def main(frames_root):
    for run, (nx, ny, nz, cols) in RUNS.items():
        rdir = os.path.join(frames_root, run)
        rc = np.genfromtxt(os.path.join(rdir, "reaction.csv"), delimiter=",", names=True)
        tstep = dict(zip(rc["step"].astype(int), rc["t"]))
        idx = np.array([(ix * ny + iy) * nz + iz for ix in cols for iy in range(ny) for iz in range(nz)])
        rows = []
        for f in sorted(glob.glob(os.path.join(rdir, "sph_output_*.txt"))):
            step = int(re.search(r"(\d{5})\.txt$", f).group(1))
            a = np.loadtxt(f, skiprows=1)
            rows.append((tstep[step], a[idx, 5].mean()))   # колонка 5 = z
        np.savetxt(os.path.join(DATA, run, "disp.csv"), np.array(rows), delimiter=",",
                   header="t,z_mean", comments="", fmt="%.10e")
        print(run, len(rows), "кадров")


if __name__ == "__main__":
    main(sys.argv[1])
