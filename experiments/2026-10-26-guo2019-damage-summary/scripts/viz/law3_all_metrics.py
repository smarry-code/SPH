"""Сводные показатели всех расчётов закона 3 (label|dir|threads). Выход: таблица markdown + json."""
import csv, glob, json, os, sys

SP = sys.argv[1]
runs = [
    ("L0: 1T 73, E78, шаг x1, без выталк. (ноутбук)", "local1T/impulse_mu0.1_v73_E78", 8),
    ("A1: 1T 73, tau 0,5, шаг x1", "law3_srv/A_tau5e-07", 4),
    ("A2: 1T 73, tau 2, шаг x1", "law3_srv/A_tau2e-06", 4),
    ("B: 1T 73, tau 0,5, шаг x0,5", "law3_srv/B_tau5e-07_dt0.5", 4),
    ("P71: 1T 71, tau 0,5, x0,5", "law3_dt05_srv/v71", 4),
    ("P73: 1T 73, tau 0,5, x0,5", "law3_dt05_srv/v73", 8),
    ("P76: 1T 76, tau 0,5, x0,5", "law3_dt05_srv/v76", 4),
    ("P78: 1T 78, tau 0,5, x0,5", "law3_dt05_srv/v78", 8),
    ("V73: 1T 73, tau 0,5, предел 20 м/с, x0,5", "law3_vmax20/v73", 4),
    ("N71: 1T 71, без выталк., x0,5", "law3_nopush/v71", 2),
    ("N80: 1T 80, без выталк., x0,5", "law3_nopush/v80", 4),
    ("T140: 2T 140, tau 0,5, x0,5", "local2T/law3_dt05_v140", 8),
    ("T142: 2T 142, tau 0,5, x0,5", "local2T/law3_dt05_v142", 8),
    ("T144: 2T 144, tau 0,5, x0,5", "local2T/law3_dt05_v144", 8),
]
out = []
for lab, d, th in runs:
    f = glob.glob(os.path.join(SP, d, "**", "shot.txt"), recursive=True)
    if not f:
        print("НЕТ", lab, d); continue
    rd = os.path.dirname(f[0])
    kv = dict(t.split("=", 1) for t in open(f[0], encoding="utf-8").read().split() if "=" in t)
    e = list(csv.DictReader(open(os.path.join(rd, "energy_log.csv"))))
    E0 = float(e[0]["E_mech"])
    t = [float(q["time"]) * 1e6 for q in e]
    c = [100 * float(q["D_contact_fiber"]) / E0 for q in e]
    worst, j = (0.0, 0.0), 0
    for i in range(len(e)):
        while t[i] - t[j] > 10:
            j += 1
        dd = c[i] - max(c[j:i + 1])
        if dd < worst[0]:
            worst = (dd, t[i])
    m = [q for q in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if q.get("w_p")]
    tl, prev = [], 0
    for q in m:
        n = int(float(q["n_broken_v"])) + int(float(q["n_broken_h"]))
        if n > prev:
            tl.append(round(float(q["t"]) * 1e6)); prev = n
    tfb = float(kv["t_first_break"]) if kv.get("t_first_break") not in (None, "nan") else float("nan")
    efb = min(e, key=lambda q: abs(float(q["time"]) - tfb)) if tfb == tfb else e[0]
    last = e[-1]
    r = dict(label=lab, dir=d, threads=th, V0=float(kv["V0"]), E=float(kv["fiber_E"]) / 1e9,
             tau=kv.get("impulse_push_tau", "0"), vmax=kv.get("impulse_push_vmax", "0"), dt_scale=kv.get("dt_scale", "1"),
             dt_mean=float(kv["dt_mean"]), outcome=kv["outcome"], V_R=float(kv["V_R"]) if kv["V_R"] != "nan" else 0.0,
             nv=int(float(kv["n_broken_v"])), nh=int(float(kv["n_broken_h"])), t_fb=round(tfb * 1e6),
             w_fb=round(1e3 * float(kv["w_f_at_first_break"]), 1), dff=float(kv["max_delta_ff/d"]),
             dfs=float(kv["max_delta_fs/d"]), ntun=int(float(kv["n_tunnel"])), t_end=round(1e6 * float(kv["t_end"])),
             ball_fb=round(100 * float(efb["KE_solid"]) / E0, 1), fric=round(100 * float(last["D_friction"]) / E0, 1),
             cont=round(100 * float(last["D_contact_fiber"]) / E0, 1), wp=round(100 * float(last["W_p_fiber"]) / E0, 2),
             ball_end=round(100 * float(last["KE_solid"]) / E0, 1), res=round(100 * float(last["residual"]) / E0, 2),
             burst=round(worst[0], 2), burst_t=round(worst[1]), timeline=tl, wall_h=round(float(kv["wall_s"]) / 3600, 1))
    out.append(r)
    print("%-48s V0 %3.0f E %2.0f tau %-6s vmax %-3s dt %-4s dt_mean %.3g thr %d | %-15s V_R %5.1f %2d/%d t_fb %4d w_fb %4.1f | dff %.2f dfs %.2f tun %3d | ball_fb %4.1f fric %4.1f cont %5.1f wp %.2f res %.2f | burst %6.2f @%4d | wall %.1f h | breaks %s" % (
        lab, r["V0"], r["E"], r["tau"], r["vmax"], r["dt_scale"], r["dt_mean"], th, r["outcome"], r["V_R"], r["nv"], r["nh"],
        r["t_fb"], r["w_fb"], r["dff"], r["dfs"], r["ntun"], r["ball_fb"], r["fric"], r["cont"], r["wp"], r["res"],
        r["burst"], r["burst_t"], r["wall_h"], r["timeline"]))
json.dump(out, open(os.path.join(SP, "law3_all_metrics.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
