import sys, math
f=sys.argv[1]; ce=0.025
mats={}; L=0; Lw=0; nb=0; zmin=1; zmax=-1; nf=0
with open(f) as fh:
    lines=fh.read().split('\n')
i=0
while i<len(lines):
    t=lines[i].split()
    if t and t[0]=='MATERIAL': mats[int(t[1])]=(float(t[4]),float(t[5]))
    if t and t[0]=='FIBER':
        m=int(t[1]); n=int(t[2]); pts=[tuple(map(float,lines[i+1+k].split()[:3])) for k in range(n)]
        nf+=1; nb+=n
        for a,b in zip(pts,pts[1:]):
            l=math.dist(a,b); L+=l
            if abs(a[1])<ce and abs(b[1])<ce: Lw+=l
        for p in pts: zmin=min(zmin,p[2]); zmax=max(zmax,p[2])
        i+=n
    i+=1
d,rho=mats[0]; A=math.pi*d*d/4
print(f"fibers={nf} beads={nb} L={L:.4f} m mass={rho*A*L*1e3:.4f} g window_mass={rho*A*Lw*1e3:.4f} g thick={1e3*(zmax-zmin)+d*1e3:.3f} mm m_per_bead={rho*A*L/nb*1e6:.4f} mg")
