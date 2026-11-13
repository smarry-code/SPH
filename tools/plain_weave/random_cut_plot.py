#!/usr/bin/env python3
"""Графики случайной раскройки по выходу генератора: <out>.fib + <out>.cut.json.

Скрипт НЕ строит геометрию заново: бусины и связи нитей читаются из .fib (FIBER-блоки,
BOND-стежки), раскройки, слои, швы, стежки и контакты -- из .cut.json (формат
random_cut_v1, см. write_random_cut_json в include/fiber/plain_weave_random_cut.hpp).
Тип нити (основа/уток) берётся из <fib>.layout.txt (колонка direction: 'x' -- основа,
'y' -- уток), если файла нет -- бусины рисуются без различия маркеров.

Выход (в --out-dir):
  cut_layout.png     -- по панели на раскройку: фрагменты, направления нитей, точки стежков швов;
  beads_contacts.png -- слой --ply: весь слой и увеличение у самого длинного шва раскройки;
  layers_zoom.png    -- тот же участок увеличения для каждого слоя (если слоёв > 1).

Все длины во входе в метрах, на графиках -- в мм.
Зависимости: numpy, matplotlib (backend Agg).

Запуск:
  python random_cut_plot.py --fib out.fib [--cut out.cut.json] [--out-dir DIR] [--ply N] [--zoom 6]
"""
import argparse
import json
import math
import os
import sys
from types import SimpleNamespace

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPolygon

M2MM = 1e3
TAB = matplotlib.colormaps["tab10"]
FRAG_RGB = np.array([TAB(i)[:3] for i in range(10)])
# цвета раскроек для стежков между слоями (без красного: красный -- стежки в слое)
CUT_COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple",
              "tab:brown", "tab:pink", "tab:olive", "tab:cyan"]

# столбцы таблицы стежков (stitch_columns в .cut.json)
SA, SB, STH, SPA, SPB, SFA, SFB, SCUT, SSEAM = range(9)


def frag_rgb(ids):
    return FRAG_RGB[np.asarray(ids, int) % 10]


def cut_color(g):
    return CUT_COLORS[g % len(CUT_COLORS)]


def inbox(P, box):
    """Маска точек P (n,2), попадающих в окно box = (x0, x1, y0, y1)."""
    x0, x1, y0, y1 = box
    return (P[:, 0] >= x0) & (P[:, 0] <= x1) & (P[:, 1] >= y0) & (P[:, 1] <= y1)


# ---------------------------------------------------------------- чтение входа
def read_fib(path):
    """Координаты бусин (м, (N,3)), размеры FIBER-блоков (нитей) и пары BOND (стежки)."""
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    blocks, bonds = [], []
    i = 0
    while i < len(lines):
        tok = lines[i].split()
        i += 1
        if not tok or tok[0].startswith("#"):
            continue
        if tok[0] == "FIBER":
            n = int(tok[2])
            blocks.append(np.array([ln.split() for ln in lines[i:i + n]], float).reshape(n, 3))
            i += n
        elif tok[0] == "BOND":
            bonds.append((int(tok[1]), int(tok[2])))
    sizes = np.array([len(b) for b in blocks], int)
    return np.vstack(blocks), sizes, np.array(bonds, int).reshape(-1, 2)


def read_layout(path, nfib):
    """Тип нити по строкам .layout.txt: 0 -- основа ('x'), 1 -- уток ('y'); None если нет."""
    if not os.path.exists(path):
        print(f"нет {path}: бусины рисуются без различия основы и утка")
        return None
    header, rows = None, []
    with open(path, encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#"):
                cols = s.lstrip("#").split()
                if header is None and "direction" in cols:
                    header = cols
                continue
            rows.append(s.split())
    if header is None or len(rows) != nfib:
        print(f"предупреждение: {path} не разобран (строк {len(rows)}, нитей {nfib}); без различия маркеров")
        return None
    di = header.index("direction")
    dirs = [r[di] for r in rows]
    if not set(dirs) <= {"x", "y"}:
        print(f"предупреждение: в {path} неизвестные direction; без различия маркеров")
        return None
    return np.array([0 if d == "x" else 1 for d in dirs], int)


def load(fib_path, cut_path):
    """Сводит .fib и .cut.json в один набор массивов; геометрия переведена в мм."""
    xyz, sizes, bonds = read_fib(fib_path)
    N, F = len(xyz), len(sizes)
    with open(cut_path, encoding="utf-8") as f:
        J = json.load(f)
    if J.get("format") != "random_cut_v1":
        raise ValueError(f"неизвестный формат {J.get('format')!r} в {cut_path}")
    if J["n_beads"] != N or J["n_fibers"] != F:
        raise ValueError(f".fib: бусин {N}, нитей {F}; .cut.json: бусин {J['n_beads']}, нитей {J['n_fibers']}")

    # слой и фрагмент каждой бусины -- по ranges (бусины и нити идут подряд по блокам)
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    bply = np.full(N, -1)
    bfrag = np.full(N, -1)
    for r in J["ranges"]:
        b0, nb, f0, nf = r["bead_first"], r["bead_count"], r["fiber_first"], r["fiber_count"]
        if starts[f0] != b0 or sizes[f0:f0 + nf].sum() != nb:
            raise ValueError(f"ranges не согласованы с FIBER-блоками: нить {f0}, бусина {b0}")
        bply[b0:b0 + nb] = r["ply"]
        bfrag[b0:b0 + nb] = r["frag"]
    if (bply < 0).any():
        raise ValueError("ranges не покрывают все бусины")
    bfiber = np.repeat(np.arange(F), sizes)

    # связи нити: соседние бусины внутри одного FIBER-блока
    ia = np.arange(N - 1)
    keep = bfiber[ia] == bfiber[ia + 1]
    chains = np.column_stack([ia[keep], ia[keep] + 1])

    S = np.array(J["stitches"], int).reshape(-1, 9)
    C = np.array(J["contacts"], int).reshape(-1, 2)
    if len(bonds) != len(S):
        print(f"предупреждение: BOND в .fib {len(bonds)} != стежков в .cut.json {len(S)}")
    if (bply[S[:, SA]] != S[:, SPA]).any() or (bply[S[:, SB]] != S[:, SPB]).any():
        print("предупреждение: слой концов стежка не совпадает со слоем бусины")

    kind_f = read_layout(fib_path + ".layout.txt", F)
    has_kind = kind_f is not None
    bkind = np.repeat(kind_f, sizes) if has_kind else np.zeros(N, int)

    cuts = []
    for c in J["cuts"]:
        frags = [dict(id=fr["id"], center=np.array(fr["center"]) * M2MM, theta=fr["theta"],
                      polygon=np.array(fr["polygon"]).reshape(-1, 2) * M2MM)
                 for fr in c["fragments"]]
        seams = []
        for s in c["seams"]:
            a = np.array(s["a"]) * M2MM
            b = np.array(s["b"]) * M2MM
            seams.append(dict(i=s["i"], j=s["j"], a=a, b=b, length=float(np.hypot(*(b - a))),
                              points=np.array(s["points"]).reshape(-1, 2) * M2MM))
        cuts.append(dict(frags=frags, seams=seams))

    return SimpleNamespace(
        xy=xyz[:, :2] * M2MM, bply=bply, bfrag=bfrag, bkind=bkind, has_kind=has_kind,
        chains=chains, S=S, C=C, cuts=cuts,
        plies=[dict(cut=p["cut"], cz=p["cz"] * M2MM) for p in J["plies"]],
        layers=len(J["plies"]), layers_per_cut=J["layers_per_cut"], n_frag=J["n_frag"],
        gap=J["gap"] * M2MM, W=J["size"][0] * M2MM, H=J["size"][1] * M2MM,
        ox=J["origin"][0] * M2MM, oy=J["origin"][1] * M2MM,
        n_beads=N, n_bonds=len(bonds),
    )


# ---------------------------------------------------------------- хелперы рисования
def zoom_window(D, l, zoom):
    """Окно увеличения: середина самого длинного шва раскройки слоя l (мм)."""
    cut = D.cuts[D.plies[l]["cut"]]
    if cut["seams"]:
        s = max(cut["seams"], key=lambda q: q["length"])
        cx, cy = (s["a"] + s["b"]) / 2.0
        title = f"шов {s['i']}–{s['j']}"
    else:
        cx, cy = D.ox + D.W / 2, D.oy + D.H / 2
        title = "центр одеяла"
    return (cx - zoom, cx + zoom, cy - zoom, cy + zoom), title


def draw_polys(ax, cut, fill_alpha=0.0, lw=0.6):
    """Границы рабочих многоугольников слоя (ec черный), заливка по фрагменту."""
    for fr in cut["frags"]:
        rgb = FRAG_RGB[fr["id"] % 10]
        ax.add_patch(MplPolygon(fr["polygon"], closed=True, fc=(*rgb, fill_alpha), ec="k",
                                lw=lw, zorder=3))


def draw_layer(ax, D, l, box=None, big=False, through_mode="blue", fill_alpha=0.0):
    """Слой l: нити, бусины, контакты, стежки в слое и концы сквозных стежков.

    box -- окно (x0, x1, y0, y1) мм или None (весь слой); big -- крупный масштаб.
    through_mode: 'blue' -- концы сквозных стежков синими треугольниками;
                  'cut'  -- полые треугольники цветом раскройки шва и отрезок до второго конца.
    """
    xy, bply, bfrag, bkind = D.xy, D.bply, D.bfrag, D.bkind
    ca, cb = D.chains[:, 0], D.chains[:, 1]

    # нити слоя: связи между соседними бусинами одного FIBER-блока
    sel = bply[ca] == l
    if box is not None:
        sel &= inbox(xy[ca], box) | inbox(xy[cb], box)
    ax.add_collection(LineCollection(np.stack([xy[ca[sel]], xy[cb[sel]]], axis=1),
                                     colors=frag_rgb(bfrag[ca[sel]]),
                                     lw=1.0 if big else 0.35, alpha=1.0 if big else 0.6, zorder=1))

    # бусины слоя: основа 'o', уток 's' (без layout -- все 'o')
    bsel = bply == l
    if box is not None:
        bsel &= inbox(xy, box)
    for kid, mk in ((0, "o"), (1, "s")):
        s = bsel & (bkind == kid)
        if s.any():
            ax.scatter(xy[s, 0], xy[s, 1], s=14 if big else 3, c=frag_rgb(bfrag[s]), marker=mk,
                       edgecolors="none", zorder=2)

    draw_polys(ax, D.cuts[D.plies[l]["cut"]], fill_alpha=fill_alpha, lw=0.9 if big else 0.6)

    # контакты: только пары бусин этого слоя
    C = D.C
    cs = (bply[C[:, 0]] == l) & (bply[C[:, 1]] == l)
    if box is not None:
        cs &= inbox(xy[C[:, 0]], box) | inbox(xy[C[:, 1]], box)
    ax.add_collection(LineCollection(np.stack([xy[C[cs, 0]], xy[C[cs, 1]]], axis=1),
                                     colors="0.25", lw=1.0 if big else 0.5, zorder=4))

    # стежки в слое (оба конца в слое l)
    S = D.S
    inl = (S[:, STH] == 0) & (S[:, SPA] == l) & (S[:, SPB] == l)
    ax.add_collection(LineCollection(np.stack([xy[S[inl, SA]], xy[S[inl, SB]]], axis=1),
                                     colors="red", lw=2.2 if big else 1.6, zorder=5))

    # сквозные стежки, касающиеся слоя l: конец в слое l
    T = S[(S[:, STH] == 1) & ((S[:, SPA] == l) | (S[:, SPB] == l))]
    here_a = T[:, SPA] == l
    end = np.where(here_a, T[:, SA], T[:, SB])
    other = np.where(here_a, T[:, SB], T[:, SA])
    oply = np.where(here_a, T[:, SPB], T[:, SPA])
    cut_g = T[:, SCUT]
    if box is not None:
        k = inbox(xy[end], box)
        end, other, oply, cut_g = end[k], other[k], oply[k], cut_g[k]
    down = oply > l  # второй конец в более глубоком слое -- стежок «вниз»
    size = 70 if big else 40
    for mask, mk in ((down, "^"), (~down, "v")):
        if not mask.any():
            continue
        if through_mode == "blue":
            ax.scatter(xy[end[mask], 0], xy[end[mask], 1], s=size, marker=mk, color="tab:blue",
                       edgecolors="white", linewidths=0.5, zorder=6)
        else:
            cols = [cut_color(g) for g in cut_g[mask]]
            ax.add_collection(LineCollection(np.stack([xy[end[mask]], xy[other[mask]]], axis=1),
                                             colors=cols, lw=1.0, zorder=5))
            ax.scatter(xy[end[mask], 0], xy[end[mask], 1], s=size, marker=mk, facecolors="none",
                       edgecolors=cols, linewidths=1.4, zorder=6)

    ax.set_aspect("equal")
    ax.set_xlabel("x, мм")
    ax.set_ylabel("y, мм")


def set_view(ax, xl, yl):
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)


def full_view(D):
    return (D.ox - 1, D.ox + D.W + 1), (D.oy - 1, D.oy + D.H + 1)


def rect_outline(ax, D):
    ax.plot([D.ox, D.ox + D.W, D.ox + D.W, D.ox, D.ox], [D.oy, D.oy, D.oy + D.H, D.oy + D.H, D.oy],
            "k-", lw=2)


# ---------------------------------------------------------------- графики
def plot_cut_layout(D, path):
    nc = len(D.cuts)
    W = D.W
    fig, axs = plt.subplots(1, nc, figsize=(7.2 * nc, 7.8), squeeze=False)
    for g, ax in enumerate(axs[0]):
        cut = D.cuts[g]
        for fr in cut["frags"]:
            rgb = FRAG_RGB[fr["id"] % 10]
            ax.add_patch(MplPolygon(fr["polygon"], closed=True, fc=(*rgb, 0.45), ec="k", lw=1.2))
        # швы других раскроек -- пунктиром
        for g2, other in enumerate(D.cuts):
            if g2 == g:
                continue
            for s in other["seams"]:
                ax.plot([s["a"][0], s["b"][0]], [s["a"][1], s["b"][1]], "--", color=cut_color(g2), lw=1.0)
        # точки стежков своих швов
        for s in cut["seams"]:
            ax.plot(s["points"][:, 0], s["points"][:, 1], "o", color="red", ms=3, zorder=5)
        r = 0.07 * W  # полудлина стрелки
        for fr in cut["frags"]:
            c, th = fr["center"], fr["theta"]
            u = np.array([math.cos(th), math.sin(th)])
            v = np.array([-u[1], u[0]])
            ax.annotate("", xy=tuple(c + r * u), xytext=tuple(c - r * u),
                        arrowprops=dict(arrowstyle="-|>", color="k", lw=1.0), zorder=6)
            ax.annotate("", xy=tuple(c + r * v), xytext=tuple(c - r * v),
                        arrowprops=dict(arrowstyle="-|>", color="0.45", lw=1.0), zorder=6)
            ax.text(c[0], c[1], f"{fr['id']}\n{math.degrees(th):+.0f}°", ha="center", va="center",
                    fontsize=9, weight="bold", zorder=7,
                    bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.2))
        rect_outline(ax, D)
        set_view(ax, *full_view(D))
        ax.set_aspect("equal")
        ax.set_xlabel("x, мм")
        ax.set_ylabel("y, мм")
        lay = ", ".join(str(l + 1) for l, p in enumerate(D.plies) if p["cut"] == g)
        ax.set_title(f"раскройка {g + 1}: слои {lay}\n"
                     "красные точки -- стежки своих швов, пунктир -- швы другой раскройки", fontsize=10)
    fig.suptitle(f"Раскройки пакета: N = {D.n_frag} фрагментов, зазор {D.gap:.2f} мм, "
                 f"слоёв {D.layers}, k = {D.layers_per_cut} (слоёв на раскройку)\n"
                 "чёрная стрелка -- основа, серая -- уток (угол theta)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_beads_contacts(D, l, zoom, path):
    g = D.plies[l]["cut"]
    box, seam_name = zoom_window(D, l, zoom)
    bply, S, C = D.bply, D.S, D.C
    nb = int((bply == l).sum())
    nc = int(((bply[C[:, 0]] == l) & (bply[C[:, 1]] == l)).sum())
    ns = int(((S[:, STH] == 0) & (S[:, SPA] == l) & (S[:, SPB] == l)).sum())
    nt = int(((S[:, STH] == 1) & ((S[:, SPA] == l) | (S[:, SPB] == l))).sum())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8.6))
    draw_layer(ax1, D, l, None, big=False, through_mode="blue")
    rect_outline(ax1, D)
    set_view(ax1, *full_view(D))
    ax1.plot([box[0], box[1], box[1], box[0], box[0]], [box[2], box[2], box[3], box[3], box[2]],
             "k--", lw=1)
    ax1.set_title(f"Слой {l + 1} целиком: бусин {nb}, контактов {nc}, стежков в слое {ns}", fontsize=10)

    draw_layer(ax2, D, l, box, big=True, through_mode="blue")
    set_view(ax2, box[0:2], box[2:4])
    ax2.set_title(f"Увеличение у {seam_name} раскройки {g + 1}, окно ±{zoom:g} мм", fontsize=10)

    handles = [Line2D([], [], marker="o", ls="none", color="0.35", ms=5,
                      label="бусины основы" if D.has_kind else "бусины нитей"),
               Line2D([], [], color="0.35", lw=1, label="связь нити (соседние бусины)"),
               Line2D([], [], color="0.25", lw=1, label=f"контакт (пара бусин слоя, < rc)"),
               Line2D([], [], color="red", lw=2, label="стежок в слое (через шов)"),
               Line2D([], [], marker="^", ls="none", color="tab:blue", ms=7,
                      label="конец сквозного стежка к слою ниже"),
               Line2D([], [], marker="v", ls="none", color="tab:blue", ms=7,
                      label="конец сквозного стежка к слою выше"),
               Line2D([], [], color="k", lw=0.9, label="граница фрагмента")]
    if D.has_kind:
        handles.insert(1, Line2D([], [], marker="s", ls="none", color="0.35", ms=5, label="бусины утка"))
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9, frameon=False)
    fig.suptitle(f"Слой {l + 1} (раскройка {g + 1}, z = {D.plies[l]['cz']:.3f} мм): бусин {nb}, "
                 f"контактов {nc}, стежков {ns} в слое и {nt} сквозных; цвет бусин -- фрагмент",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0.1, 1, 0.94))
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_layers_zoom(D, l0, zoom, path):
    L = D.layers
    box, seam_name = zoom_window(D, l0, zoom)
    ncol = math.ceil(L / 2)
    fig, axs = plt.subplots(2, ncol, figsize=(6.8 * ncol, 13.2), squeeze=False)
    axs = axs.ravel()
    for l in range(L):
        ax = axs[l]
        draw_layer(ax, D, l, box, big=True, through_mode="cut", fill_alpha=0.12)
        set_view(ax, box[0:2], box[2:4])
        ax.set_title(f"слой {l + 1} (раскройка {D.plies[l]['cut'] + 1}), "
                     f"z = {D.plies[l]['cz']:.3f} мм", fontsize=11)
    for ax in axs[L:]:
        ax.axis("off")

    handles = [Line2D([], [], color="0.4", lw=1, label="нить: связи соседних бусин (цвет -- фрагмент)"),
               Line2D([], [], color="red", lw=2.2, label="стежок в слое (через шов)"),
               Line2D([], [], marker="^", ls="none", mfc="none", color="tab:blue", ms=7,
                      label="сквозной стежок к слою ниже (цвет -- раскройка шва)"),
               Line2D([], [], marker="v", ls="none", mfc="none", color="tab:orange", ms=7,
                      label="сквозной стежок к слою выше (цвет -- раскройка шва)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=10, frameon=False)
    fig.suptitle(f"Участок {2 * zoom:g}×{2 * zoom:g} мм у {seam_name} раскройки "
                 f"{D.plies[l0]['cut'] + 1} (окно по слою {l0 + 1}); все слои пакета",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------- главное
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fib", required=True, help="файл пакета .fib")
    ap.add_argument("--cut", default=None, help="json раскройки (по умолчанию <fib>.cut.json)")
    ap.add_argument("--out-dir", default=None, help="каталог графиков (по умолчанию каталог .fib)")
    ap.add_argument("--ply", type=int, default=0, help="слой для beads_contacts.png (с нуля)")
    ap.add_argument("--zoom", type=float, default=6.0, help="полуширина окна увеличения, мм")
    a = ap.parse_args()
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")  # сообщения по-русски в любой консоли

    cut_path = a.cut or a.fib + ".cut.json"
    out_dir = a.out_dir or os.path.dirname(os.path.abspath(a.fib))
    os.makedirs(out_dir, exist_ok=True)

    D = load(a.fib, cut_path)
    if not 0 <= a.ply < D.layers:
        raise SystemExit(f"--ply {a.ply} вне диапазона 0..{D.layers - 1}")

    print(f"бусин {D.n_beads}, связей BOND {D.n_bonds}, слоёв {D.layers}, раскроек {len(D.cuts)}")
    created = [
        plot_cut_layout(D, os.path.join(out_dir, "cut_layout.png")),
        plot_beads_contacts(D, a.ply, a.zoom, os.path.join(out_dir, "beads_contacts.png")),
    ]
    if D.layers > 1:
        created.append(plot_layers_zoom(D, a.ply, a.zoom, os.path.join(out_dir, "layers_zoom.png")))
    for p in created:
        print(p)


if __name__ == "__main__":
    main()
