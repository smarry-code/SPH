#!/usr/bin/env python3
# Схема осесимметричной геометрии кейса BC2 (меридиональное сечение).
# Масштаб по обеим осям одинаковый — пропорции сохранены.
#
# Запуск:  python3 notes/contact/e2e/mkgeom.py
# Результат: notes/contact/e2e/bc2_geometry.svg
#
# Источник чисел — bc2_layers() в tests/test_e2e_bc2.cpp. Схема иллюстративная
# и строится по НОМИНАЛЬНЫМ границам слоёв; узлы в расчёте ставятся по центрам
# ячеек со своим шагом в каждом слое (cell_count там же), поэтому фактические
# крайние узлы отстоят от границ на полшага.
#
# Скрипт лежит рядом с картинкой намеренно: прежняя версия жила во временном
# каталоге, из-за чего схема разошлась с постановкой и полгода показывала
# зазоры, которых уже не было.

OUT = "/Users/vladimirbogachev/PycharmProjects/SPH/notes/contact/e2e/bc2_geometry.svg"

# Опорные плоскости (мм), подтверждены статьёй Grady & Moody, SAND96-0551,
# опыт BC2 = ориг. CE6: толщины 6,000 / 3,917 / 9,033 / 25,400 мм.
Z1, Z2, Z3, Z4 = 6.0, 9.917, 18.950, 44.350

# Слои задаются номинальными границами: тела в КОНТАКТЕ, искусственных
# зазоров нет (прежние +dx / +4dx / +5dx убраны — статья их не подтверждает,
# а сами они зависели от разрешения).
layers = [
    (0.0, Z1, 43.75, "1", "var(--m-pu)",  2210),
    (Z1,  Z2, 34.60, "2", "var(--m-b4c)", 2210),
    (Z2,  Z3, 34.60, "3", "var(--m-b4c)", 0),
    (Z3,  Z4, 25.40, "4", "var(--m-lif)", 0),
]

# Частица-монитор: ось симметрии, тыл слоя 3 (граница с окном).
# В расчёте это последний узел по z внутреннего радиального столбца.
MONI_Z = Z3 - (Z3 - Z2) / 45 / 2     # центр последней ячейки слоя 3
MONI_R = 0.1                          # 0.5 * dr внутреннего столбца

S = 9.0
ML, MT, MR, MB = 78, 104, 116, 162
ZMAX, RMAX = Z4, 43.75
PW, PH = ZMAX * S, RMAX * S
W, H = int(ML + PW + MR), int(MT + PH + MB)


def X(z): return ML + z * S
def Y(r): return MT + PH - r * S


def num(v):
    return f"{v:.3f}".rstrip("0").rstrip(".").replace(".", ",")


p = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
     f'role="img" aria-label="Меридиональное сечение осесимметричной сборки BC2: '
     f'ударник из полиуретана, две пластины карбида бора и окно из фторида лития">']
p.append("""<style>
svg { --surface-1:#fcfcfb; --text-primary:#0b0b0b; --text-secondary:#52514e;
      --axis:#b9b8b2; --dim:#8a8984; --moni:#c2255c;
      --m-pu:#eb6834; --m-b4c:#2a78d6; --m-lif:#1baf7a;
      font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }
@media (prefers-color-scheme: dark) {
  svg { --surface-1:#1a1a19; --text-primary:#ffffff; --text-secondary:#c3c2b7;
        --axis:#5a5955; --dim:#6f6e69; --moni:#f06595;
        --m-pu:#d95926; --m-b4c:#3987e5; --m-lif:#199e70; }
}
.ttl { fill:var(--text-primary); font-size:15px; font-weight:600; }
.sub { fill:var(--text-secondary); font-size:11.5px; }
.lb  { fill:#ffffff; font-size:14px; font-weight:700; }
.dm  { fill:var(--dim); font-size:10px; }
.ax  { fill:var(--text-secondary); font-size:10.5px; }
.lg  { fill:var(--text-primary); font-size:11px; }
.vl  { font-size:11px; font-weight:600; }
.mn  { fill:var(--moni); font-size:11px; font-weight:600; }
</style>""")
p.append('<defs><marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
         'orient="auto"><path d="M0,1 L9,5 L0,9 z" fill="context-stroke"/></marker></defs>')
p.append(f'<rect width="{W}" height="{H}" fill="var(--surface-1)"/>')
p.append(f'<text class="ttl" x="{ML}" y="26">Геометрия кейса BC2 — меридиональное сечение</text>')
p.append(f'<text class="sub" x="{ML}" y="44">Осесимметричная постановка, показана половина r &gt; 0.</text>')
p.append(f'<text class="sub" x="{ML}" y="60">Масштаб по обеим осям одинаковый; размеры в мм.</text>')

# тела
for z0, z1, rm, lab, col, v in layers:
    x0, x1, y0, y1 = X(z0), X(z1), Y(rm), Y(0)
    p.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.1f}" height="{y1-y0:.1f}" '
             f'fill="{col}" fill-opacity="0.82" stroke="{col}" stroke-width="1.5"/>')
    cx = (x0 + x1) / 2
    p.append(f'<text class="lb" x="{cx:.1f}" y="{(y0+y1)/2+5:.1f}" text-anchor="middle">{lab}</text>')
    p.append(f'<text class="dm" x="{cx:.1f}" y="{Y(0)+16:.1f}" text-anchor="middle">{num(z1-z0)}</text>')
    if v:
        ay = y0 - 14
        p.append(f'<line x1="{cx-18:.1f}" y1="{ay:.1f}" x2="{cx+18:.1f}" y2="{ay:.1f}" '
                 f'stroke="{col}" stroke-width="2" marker-end="url(#ar)"/>')
        p.append(f'<text class="vl" x="{cx:.1f}" y="{ay-7:.1f}" text-anchor="middle" fill="{col}">{v} м/с</text>')

# ось симметрии
p.append(f'<line x1="{ML-20}" y1="{Y(0):.1f}" x2="{X(ZMAX)+20:.1f}" y2="{Y(0):.1f}" '
         f'stroke="var(--axis)" stroke-width="1.2" stroke-dasharray="10 4 2 4"/>')
p.append(f'<text class="ax" x="{ML-24}" y="{Y(0)+4:.1f}" text-anchor="end">r = 0</text>')
p.append(f'<text class="ax" x="{X(ZMAX)+24:.1f}" y="{Y(0)+4:.1f}">z</text>')
p.append(f'<text class="ax" x="{ML-52}" y="{Y(RMAX/2):.1f}" text-anchor="middle" '
         f'transform="rotate(-90 {ML-52} {Y(RMAX/2):.1f})">r, мм</text>')

# частица-монитор: на оси, тыл слоя 3 (граница с окном)
mx, my = X(MONI_Z), Y(MONI_R)
p.append(f'<circle cx="{mx:.1f}" cy="{my:.1f}" r="4.5" fill="var(--moni)" '
         f'stroke="var(--surface-1)" stroke-width="1.5"/>')
p.append(f'<line x1="{mx:.1f}" y1="{my-8:.1f}" x2="{mx:.1f}" y2="{my-46:.1f}" '
         f'stroke="var(--moni)" stroke-width="1.2"/>')
p.append(f'<text class="mn" x="{mx:.1f}" y="{my-52:.1f}" text-anchor="middle">монитор</text>')
p.append(f'<text class="mn" x="{mx:.1f}" y="{my-38:.1f}" text-anchor="middle" '
         f'style="font-weight:400;font-size:9.5px">на оси</text>')

# засечки по r слева
for rv in (10, 20, 30, 43.75):
    y = Y(rv)
    p.append(f'<line x1="{ML-8}" y1="{y:.1f}" x2="{ML-2}" y2="{y:.1f}" stroke="var(--dim)" stroke-width="1"/>')
    p.append(f'<text class="dm" x="{ML-12}" y="{y+3.4:.1f}" text-anchor="end">{num(rv)}</text>')

# радиусы тел — справа
for rm, who in ((43.75, "тело 1"), (34.60, "тела 2, 3"), (25.40, "тело 4")):
    y = Y(rm)
    p.append(f'<line x1="{X(ZMAX)+4:.1f}" y1="{y:.1f}" x2="{X(ZMAX)+14:.1f}" y2="{y:.1f}" '
             f'stroke="var(--dim)" stroke-width="1"/>')
    p.append(f'<text class="dm" x="{X(ZMAX)+18:.1f}" y="{y+3.4:.1f}">R = {num(rm)} ({who})</text>')

# легенда
ly = MT + PH + 46
items = [("1", "полиуретановая подложка, itype 8, ρ₀ = 640 кг/м³", "var(--m-pu)"),
         ("2, 3", "карбид бора B4C (Eagle Picher), itype 10, ρ₀ = 2517 кг/м³", "var(--m-b4c)"),
         ("4", "окно из фторида лития LiF, itype 5, ρ₀ = 2640 кг/м³", "var(--m-lif)")]
lx = ML
for lab, txt, col in items:
    p.append(f'<rect x="{lx}" y="{ly-9}" width="12" height="12" rx="2" fill="{col}" fill-opacity="0.82"/>')
    p.append(f'<text class="lg" x="{lx+18}" y="{ly+1}">{lab} — {txt}</text>')
    ly += 20

p.append(f'<text class="dm" x="{ML}" y="{H-32}">Числа под осью — толщина тела, мм. '
         f'Тела 1 и 2 летят как единый снаряд; тела 3 и 4 покоятся. Удар вдоль z слева направо.</text>')
p.append(f'<text class="dm" x="{ML}" y="{H-16}">Тела в контакте: искусственных зазоров нет. '
         f'Монитор — на оси, у тыльной поверхности тела 3 (граница с окном), как в опыте.</text>')
p.append("</svg>")

open(OUT, "w", encoding="utf-8").write("\n".join(p))
print("записано:", OUT, f"({W}x{H})")
