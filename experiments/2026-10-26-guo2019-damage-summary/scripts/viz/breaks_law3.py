"""Разбор разрывов по выстрелам: label=run_dir (каталог с shot.txt).
Хронология разрывов (monitor), какие нити и где рвутся (кадры), напряжения перед первым разрывом,
энергия шара/ткани к первому разрыву и в конце, импульс в закрепление, пик реакции закрепления."""
import csv, glob, json, math, os, sys

SIG_U = 2.5e9


def read_frame(path, want_pos):
    t, pos, bonds, mode = None, {}, {}, None
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("TIME"):
            t = float(line.split()[1]); continue
        if line.startswith("PARTICLES"):
            mode = "p"; continue
        if line.startswith("BONDS"):
            mode = "b"; continue
        if not line.strip() or line[0] == "#":
            continue
        s = line.split()
        if mode == "p" and want_pos:
            pos[int(s[0])] = (float(s[3]), float(s[4]), float(s[5]))
        elif mode == "b":
            bonds[(int(s[0]), int(s[1]))] = (float(s[5]), float(s[6]), float(s[7]))
    return t, pos, bonds


def layout(rd):
    ys = []
    for line in open(os.path.join(rd, "guo2019_fabric_layout.txt"), encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        s = line.split()
        ys.append((s[1], float(s[3]), int(s[5]), int(s[6]), int(s[2]) if s[2].lstrip("-").isdigit() else 0))
    return ys


def yarn_of(ys, i):
    for k, y in enumerate(ys):
        if y[2] <= i <= y[3]:
            return k
    return -1


out = {}
for arg in sys.argv[1:]:
    label, root = arg.split("=", 1)
    rd = os.path.dirname(glob.glob(os.path.join(root, "**", "shot.txt"), recursive=True)[0])
    kv = dict(t.split("=", 1) for t in open(os.path.join(rd, "shot.txt"), encoding="utf-8").read().split() if "=" in t)
    V0 = float(kv["V0"])
    mon = [q for q in csv.DictReader(open(os.path.join(rd, "monitor.csv"))) if q.get("w_p")]
    en = list(csv.DictReader(open(os.path.join(rd, "energy_log.csv"))))
    E0 = float(en[0]["E_mech"])
    m_ball = 2 * float(en[0]["KE_solid"]) / V0 ** 2
    tfb = float(kv["t_first_break"])
    tl, prev = [], 0
    for q in mon:
        n = int(float(q["n_broken_v"])) + int(float(q["n_broken_h"]))
        if n > prev:
            tl.append((round(float(q["t"]) * 1e6, 1), n)); prev = n
    efb = min(en, key=lambda q: abs(float(q["time"]) - tfb))
    J, fmax, tfmax = 0.0, 0.0, 0.0
    for a, b in zip(mon, mon[1:]):
        J += 0.5 * (float(a["Ff_z"]) + float(b["Ff_z"])) * (float(b["t"]) - float(a["t"]))
        if abs(float(b["Ff_z"])) > abs(fmax):
            fmax, tfmax = float(b["Ff_z"]), float(b["t"])
    vfb = float(min(mon, key=lambda q: abs(float(q["t"]) - tfb))["v_z"])
    r = {"V0": V0, "outcome": kv["outcome"], "V_R": float(kv["V_R"]), "nv": int(float(kv["n_broken_v"])),
         "nh": int(float(kv["n_broken_h"])), "t_fb_us": round(tfb * 1e6, 1),
         "w_fb_mm": round(float(kv["w_f_at_first_break"]) * 1e3, 2), "v_fb": round(vfb, 1),
         "KEball_fb_pct": round(100 * float(efb["KE_solid"]) / E0, 1),
         "KEfab_fb_pct": round(100 * float(efb["KE_fiber"]) / E0, 1),
         "Ubond_fb_pct": round(100 * float(efb["U_bond"]) / E0, 1),
         "fric_fb_pct": round(100 * float(efb["D_friction"]) / E0, 1),
         "fric_end_pct": round(100 * float(en[-1]["D_friction"]) / E0, 1),
         "cont_end_pct": round(100 * float(en[-1]["D_contact_fiber"]) / E0, 1),
         "Wp_end_pct": round(100 * float(en[-1]["W_p_fiber"]) / E0, 2),
         "KEfab_max_pct": round(100 * max(float(q["KE_fiber"]) for q in en) / E0, 1),
         "Ubond_max_pct": round(100 * max(float(q["U_bond"]) for q in en) / E0, 1),
         "F_clamp_max_N": round(fmax), "t_F_clamp_max_us": round(tfmax * 1e6), "J_clamp_Ns": round(J, 4),
         "p_ball_lost_Ns": round(m_ball * (V0 - abs(float(mon[-1]["v_z"]))), 4), "timeline": tl}
    frames = sorted(glob.glob(os.path.join(rd, "fiber_output_*.txt")), key=lambda p: int(os.path.basename(p)[13:-4]))
    if frames:
        ys = layout(rd)
        t0, pos0, b0 = read_frame(frames[0], True)
        prev_b, prev_t, before, br = b0, t0, None, []
        for fp in frames[1:]:
            t, _, bb = read_frame(fp, False)
            gone = [k for k in prev_b if k not in bb]
            if gone and before is None:
                before = (prev_t, prev_b)
            for (i, j) in gone:
                y = yarn_of(ys, i)
                dr, off = ys[y][0], ys[y][1]
                x, yy, _ = pos0[i]
                br.append({"t_us": round(t * 1e6), "yarn": y, "dir": dr, "off_mm": round(off * 1e3, 2),
                           "along_mm": round((yy if dr == "y" else x) * 1e3, 1)})
            prev_b, prev_t = bb, t
        # нити, порванные хоть раз: смещение от центра
        yarns = {}
        for b in br:
            yarns.setdefault((b["dir"], b["off_mm"]), b["t_us"])
        r["broken_yarns"] = sorted(((d, o, t) for (d, o), t in yarns.items()), key=lambda z: z[2])
        if before:
            tb, bb = before
            sig = [v[2] for v in bb.values()]
            r["pre_t_us"] = round(tb * 1e6)
            r["pre_n_sig_gt_0.8su"] = sum(s > 0.8 * SIG_U for s in sig)
            r["pre_n_sig_gt_0.9su"] = sum(s > 0.9 * SIG_U for s in sig)
            # нити (вертикальные), у которых max sigma > 0.9 sigma_u
            hot = {}
            for (i, j), v in bb.items():
                if v[2] > 0.9 * SIG_U:
                    y = yarn_of(ys, i)
                    hot[(ys[y][0], round(ys[y][1] * 1e3, 2))] = max(hot.get((ys[y][0], round(ys[y][1] * 1e3, 2)), 0), v[2])
            r["pre_hot_yarns"] = sorted((d, o, round(s / 1e9, 3)) for (d, o), s in hot.items())
    out[label] = r
    print("==", label, json.dumps({k: v for k, v in r.items() if k not in ("timeline", "broken_yarns", "pre_hot_yarns")},
                                 ensure_ascii=False))
    print("   хронология (мкс, число разрывов):", r["timeline"])
    print("   порванные нити (dir, смещение мм, t мкс):", r.get("broken_yarns"))
    print("   нити с sigma > 0.9 sigma_u перед первым разрывом:", r.get("pre_hot_yarns"))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "breaks_law3.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
