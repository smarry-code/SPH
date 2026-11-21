# -*- coding: utf-8 -*-
# Сборка фигуры bc2_scheme_study.svg по кривым из ../data.
# Запуск: python scripts/fig.py — пути считаются от положения самого скрипта.
import matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
EXP  = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(EXP, "..", ".."))

S    = os.path.join(EXP, "data")
REF  = os.path.join(ROOT, "contact", "e2e", "gold_standart.txt")
OUT  = os.path.join(EXP, "bc2_scheme_study.svg")

def rd(p, tscale):
    T=[];V=[]
    for L in open(p, encoding="utf-8", errors="replace"):
        L=L.strip().lstrip("\ufeff")
        if not L or L[0] in "#*": continue
        a=L.split()
        try: t=float(a[0]); v=float(a[1])
        except ValueError: continue
        T.append(t*tscale); V.append(v)
    return T,V

def t50(T,V):
    pk=max(V); ip=V.index(pk); lv=0.5*pk
    for i in range(ip,0,-1):
        if V[i-1]<=lv<=V[i]:
            d=V[i]-V[i-1]
            return (T[i-1]+((lv-V[i-1])/d)*(T[i]-T[i-1])) if d>0 else T[i], pk
    return T[ip], pk

RT,RV = rd(REF,1.0)                 # эталон уже в мкс
r50,rpk = t50(RT,RV)

def series(ax, fname, label, **kw):
    T,V = rd(os.path.join(S, fname), 1e6)
    s,pk = t50(T,V)
    ax.plot([t-s+r50 for t in T], V, label=f"{label} — пик {pk:.0f}", **kw)
    return pk

fig,axes = plt.subplots(1,2, figsize=(15,6.2), sharey=True)
for ax in axes:
    ax.plot(RT,RV,color="black",lw=2.4,label=f"эталон (приборный) — пик {rpk:.0f}")
    ax.axhline(1105,color="gray",ls=":",lw=1.2)
    ax.text(3.35,1120,"$u_p$ = 1105 м/с (симметричный удар)",fontsize=8.5,color="gray",ha="right")
    ax.set_xlim(1.6,3.6); ax.set_ylim(0,1450)
    ax.set_xlabel("время, мкс  (кривые совмещены по $t_{50}$ фронта)")
    ax.grid(alpha=0.3)

ax=axes[0]
series(ax,"curve0.txt",      "Δx = 0,20 мм", color="#1f77b4", lw=1.6)
series(ax,"curve_dx14e4.txt","Δx = 0,14 мм", color="#2ca02c", lw=1.6)
series(ax,"curve_dx10e4.txt","Δx = 0,10 мм", color="#d62728", lw=1.6)
ax.set_title("Сгущение сетки: амплитуда не меняется")
ax.set_ylabel("осевая скорость монитора $v_z$, м/с")
ax.legend(loc="upper left",fontsize=9)

ax=axes[1]
series(ax,"curve0.txt",      "ε = 0,25, CFL = 0,1 (база)", color="#1f77b4", lw=1.6)
series(ax,"curve_eps005.txt","ε = 0,05",                   color="#ff7f0e", lw=1.6)
series(ax,"curve_cfl03.txt", "CFL = 0,3",                  color="#9467bd", lw=1.6)
series(ax,"curve_xsph.txt",  "XSPH на переносе (обрыв)",   color="#8c564b", lw=1.3, ls="--", alpha=0.85)
ax.set_title("Параметры схемы: амплитуда не меняется")
ax.legend(loc="upper left",fontsize=9)

fig.suptitle("BC2: расчёт против приборного эталона — граница B4C/LiF на оси", fontsize=13)
fig.tight_layout(rect=[0,0.03,1,0.96])
fig.text(0.5,0.005,"Эталон: Grady & Moody, SAND96-0551, опыт BC2 (CE6). Штриховая кривая обрывается на 2,28 мкс "
                   "(отрицательная плотность); её хвост — численная пила, не волна.",ha="center",fontsize=8.5)
fig.savefig(OUT)
print("записано:",OUT)
