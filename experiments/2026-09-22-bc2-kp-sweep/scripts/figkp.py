# -*- coding: utf-8 -*-
# Сборка фигуры bc2_kp_sweep.svg по кривым из ../data.
# Запуск:  python scripts/figkp.py     (из каталога эксперимента или любого другого —
#          пути считаются от положения самого скрипта)
import matplotlib, os, bisect, sys
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
EXP  = os.path.dirname(HERE)                       # каталог эксперимента
ROOT = os.path.abspath(os.path.join(EXP, "..", "..")) # корень репозитория

D    = os.path.join(EXP, "data")
REF  = os.path.join(ROOT, "contact", "e2e", "gold_standart.txt")
OUT  = os.path.join(EXP, "bc2_kp_sweep.svg")
KPS  = [0.5, 1, 2, 4, 7, 10]

def rd(p,tscale):
    T=[];V=[];st=0
    for L in open(p,encoding="utf-8",errors="replace"):
        L=L.strip().lstrip("\ufeff")
        if L.startswith("#"):
            if "steps" in L: st=int(L.split("=")[1])
            continue
        if not L or L[0]=="*": continue
        a=L.split()
        try: t=float(a[0]); v=float(a[1])
        except ValueError: continue
        T.append(t*tscale); V.append(v)
    return T,V,st

def feat(T,V):
    pk=max(V); ip=V.index(pk)
    def cr(lv):
        for i in range(ip,0,-1):
            if V[i-1]<=lv<=V[i]:
                d=V[i]-V[i-1]
                return T[i-1]+((lv-V[i-1])/d)*(T[i]-T[i-1]) if d>0 else T[i]
        return T[0]
    t10,t50,t90=cr(.1*pk),cr(.5*pk),cr(.9*pk)
    fall=None
    for i in range(ip,len(V)-1):
        if V[i]>=.7*pk>=V[i+1]:
            d=V[i]-V[i+1]
            fall=(T[i]+((V[i]-.7*pk)/d)*(T[i+1]-T[i]))-T[ip] if d>0 else T[i]-T[ip]
            break
    return pk,T[ip],t10,t50,t90,fall

RT,RV,_ = rd(REF,1.0)
rpk,rtp,r10,r50,r90,rfall = feat(RT,RV)
print("эталон: пик=%.1f фронт=%.4f спад=%.4f"%(rpk,r90-r10,rfall))
print()
print("тип   Kp   шагов     пик   доля  фронт10-90  спад70  шоулдер")
rows={}
fig,axes=plt.subplots(1,2,figsize=(15.5,6.4),sharey=True)
cmap=plt.get_cmap("viridis")
for ax,ty,name in ((axes[0],1,"тип 1 — Пинбол"),(axes[1],2,"тип 2 — Свигл")):
    ax.plot(RT,RV,color="black",lw=2.6,label="эталон (прибор) — пик 1318",zorder=5)
    for n,kp in enumerate(KPS):
        tag="t%d_kp%s"%(ty,str(kp).replace(".",""))
        f=os.path.join(D,"curve_%s.txt"%tag)
        T,V,st=rd(f,1e6)
        pk,tp,t10,t50,t90,fall=feat(T,V)
        i=bisect.bisect_left(T,t10+0.15); sh=V[min(i,len(V)-1)]
        rows[(ty,kp)]=(st,pk,t90-t10,fall,sh)
        print("%3d %5.1f %7d %8.1f %6.3f %10.4f %8s %8.0f"%(
            ty,kp,st,pk,pk/rpk,t90-t10,("%.4f"%fall) if fall else "—",sh))
        ax.plot([t-t50+r50 for t in T],V,color=cmap(n/(len(KPS)-1)),lw=1.5,
                label="Kp = %g — пик %.0f, шагов %d"%(kp,pk,st))
    ax.axhline(1105,color="gray",ls=":",lw=1.1)
    ax.set_xlim(1.6,3.6); ax.set_ylim(0,1450); ax.grid(alpha=0.3)
    ax.set_xlabel("время, мкс  (кривые совмещены по $t_{50}$ фронта)")
    ax.set_title(name)
    ax.legend(loc="upper left",fontsize=8.5)
axes[0].set_ylabel("осевая скорость монитора $v_z$, м/с")
axes[0].text(3.5,1120,"$u_p$ = 1105 м/с",fontsize=8.5,color="gray",ha="right")
fig.suptitle("BC2: развёртка по масштабу штрафной силы Kp, контакт типов 1 и 2",fontsize=13)
fig.tight_layout(rect=[0,0.035,1,0.955])
fig.text(0.5,0.006,"Эталон: Grady & Moody, SAND96-0551, опыт BC2 (CE6). Все двенадцать прогонов дошли до t_total. "
                   "Рост числа шагов у типа 2 — работа ограничителя dt_limit = 1,6·dt_ref/Kp.",ha="center",fontsize=8.5)
fig.savefig(OUT)
print("\nзаписано:",OUT)
