"""Как рвутся нити: по monitor.csv (запись каждые 10 шагов, ~0,1-0,2 мкс) -- моменты приращения числа порванных
нитей, «залпы» (приращение > 1 за одну запись монитора = одновременный разрыв), кластеры (разрывы ближе --gap мкс),
разброс времени; по energy_log.csv -- D_break (упругая энергия, потерянная в момент удаления связи) и W_p_fiber.
Usage: python break_timing.py out.json label=run_dir ...   (run_dir -- каталог с monitor.csv)
"""
import csv, json, os, sys

GAP_US = 2.0


def mon(rd):
    rows = [r for r in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if r.get("w_p")]
    return rows


def energy_end(rd):
    p = os.path.join(rd, "energy_log.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p)) if r.get("D_contact_fiber")]
    if not rows:
        return None
    first, last = rows[0], rows[-1]
    e0 = float(first["E_mech"])
    g = lambda k: float(last[k]) if last.get(k) not in (None, "") else float("nan")
    return {"E0": e0, "D_break": g("D_break"), "W_p_fiber": g("W_p_fiber"), "D_break_frac_E0": g("D_break") / e0,
            "W_p_fiber_frac_E0": g("W_p_fiber") / e0}


out = {}
for a in sys.argv[2:]:
    lab, rd = a.split("=", 1)
    m = mon(rd)
    ev = []  # (t_us, dn, dv, dh, w_f_mm)
    prev = None
    for r in m:
        nv, nh = int(float(r["n_broken_v"])), int(float(r["n_broken_h"]))
        if prev is not None and nv + nh > prev[0] + prev[1]:
            ev.append((float(r["t"]) * 1e6, nv + nh - prev[0] - prev[1], nv - prev[0], nh - prev[1],
                       float(r["w_f"]) * 1e3))
        prev = (nv, nh)
    n = sum(e[1] for e in ev)
    clusters = []
    for e in ev:
        if clusters and e[0] - clusters[-1][-1][0] < GAP_US:
            clusters[-1].append(e)
        else:
            clusters.append([e])
    dt_rec = (float(m[-1]["t"]) - float(m[0]["t"])) / max(len(m) - 1, 1) * 1e6
    res = {
        "n_broken": n, "n_events": len(ev), "rec_us": dt_rec,
        "t_first": ev[0][0] if ev else None, "t_last": ev[-1][0] if ev else None,
        "span_us": ev[-1][0] - ev[0][0] if ev else None,
        "max_simult": max((e[1] for e in ev), default=0),
        "n_in_multi": sum(e[1] for e in ev if e[1] > 1),
        "n_clusters": len(clusters), "max_cluster": max((sum(x[1] for x in c) for c in clusters), default=0),
        "first3_span_us": None,
        "events": [[round(e[0], 2), e[1], e[2], e[3], round(e[4], 1)] for e in ev],
        "energy": energy_end(rd),
    }
    ts = [e[0] for e in ev for _ in range(e[1])]
    if len(ts) >= 3:
        res["first3_span_us"] = ts[2] - ts[0]
    if len(ts) >= 2:
        gaps = sorted(b - a for a, b in zip(ts, ts[1:]))
        res["median_gap_us"] = gaps[len(gaps) // 2]
    out[lab] = res
    print("%-22s n=%2d событий=%2d запись=%.2f мкс  t1=%6.1f  разброс=%6.1f  макс.одновр.=%d  в залпах=%2d  "
          "кластеров(<%g мкс)=%2d макс.кластер=%2d  медиана интервала=%s" % (
              lab, n, len(ev), dt_rec, res["t_first"] or -1, res["span_us"] or -1, res["max_simult"],
              res["n_in_multi"], GAP_US, len(clusters), res["max_cluster"],
              "%.1f" % res["median_gap_us"] if "median_gap_us" in res else "-"))
    if res["energy"]:
        e = res["energy"]
        print("    энергия: D_break = %.4f Дж (%.2f %% E0), W_p_fiber = %.4f Дж (%.2f %% E0)" % (
            e["D_break"], 100 * e["D_break_frac_E0"], e["W_p_fiber"], 100 * e["W_p_fiber_frac_E0"]))
json.dump(out, open(sys.argv[1], "w"), ensure_ascii=False, indent=1)
