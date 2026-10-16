import csv,glob,statistics as st
for fn in sorted(glob.glob('*_v0.05*.csv')):
    rows=[r for r in csv.DictReader(open(fn)) if int(r['n_geom'])>0]
    sg=[]; 
    for r in rows:
        x=float(r['imp_q'])
        sg.append(0 if x==0 else (1 if x>0 else -1))
    runs=[];cur=None;L=0
    for s in sg:
        if s==0: continue
        if s==cur: L+=1
        else:
            if cur is not None: runs.append(L)
            cur=s;L=1
    runs.append(L)
    short=sum(1 for r in runs if r<=2)
    # period of vy_rel sign changes in steps
    print(f"{fn:28s} runs={len(runs):4d} median_run={st.median(runs):6.1f} short(<=2)={short:4d} ({100*short/max(1,len(runs)):.0f}%) max_run={max(runs)}")
