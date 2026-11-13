#!/usr/bin/env python3
"""metrics.py <root> [--out DIR]
Таблицы по всем прогонам серии (root/<RUN>/{shot.txt,monitor.csv,energy_log.csv}) ->
metrics.json и tables.md. Критерии ТЗ — столбцы pass/fail."""
import argparse, json, os
import numpy as np
from common import find_runs, read_shot, read_csv_dedup

E_JUMP = 1e-3      # скачок residual между строками energy_log, доля E0
DCF_DROP = -1e-6   # допустимое убывание D_contact_fiber между строками, доля E0
RES_MARGIN = 0.005  # допуск residual над C3


def fnum(x, f="{:.4g}"):
    if x is None:
        return "-"
    if isinstance(x, str):
        return x
    if isinstance(x, (bool, np.bool_)):
        return "pass" if x else "FAIL"
    if isinstance(x, float) and not np.isfinite(x):
        return "nan"
    return f.format(x)


def run_metrics(root, run):
    d = os.path.join(root, run)
    shot = read_shot(os.path.join(d, "shot.txt"))
    M = {"shot": {k: shot.get(k) for k in [
        "outcome", "V_R", "t_end", "steps", "n_tunnel", "max_delta_ff/d", "max_delta_fs/d", "wall_s", "wall_per_step_s",
        "t_first_break", "n_broken_v", "n_broken_h", "proj_resid_max", "proj_resid_rel_max", "W_proj_pos_max",
        "proj_slots_max", "contact_proj_iters", "contact_proj_warm", "E_tot_end/E_tot_0", "dt_mean"]}}
    # --- monitor ---
    _, mo = read_csv_dedup(os.path.join(d, "monitor.csv"))
    step, t = mo["step"], mo["t"]
    nb = np.nan_to_num(mo["n_broken_v"]) + np.nan_to_num(mo["n_broken_h"])
    M["monitor_rows"] = int(len(step))
    M["last_step"] = int(step[-1])
    ib = np.nonzero(nb > 0)[0]
    M["t_first_break_monitor"] = float(t[ib[0]]) if len(ib) else None
    M["step_first_break_monitor"] = int(step[ib[0]]) if len(ib) else None
    M["n_broken_v_end"] = int(mo["n_broken_v"][-1])
    M["n_broken_h_end"] = int(mo["n_broken_h"][-1])
    if "proj_resid" in mo:
        act = np.nan_to_num(mo["proj_n_ff"]) + np.nan_to_num(mo["proj_n_fs"]) > 0
        for key in ("proj_resid", "proj_resid_rel"):
            v = mo[key][act & np.isfinite(mo[key])]
            M[key + "_max"] = float(v.max()) if len(v) else None
            M[key + "_median"] = float(np.median(v)) if len(v) else None
        M["proj_active_rows"] = int(act.sum())
        wf = mo["proj_warm_frac"]
        M["proj_warm_frac_min"] = float(np.nanmin(wf)) if np.isfinite(wf).any() else None
        w = mo["W_proj_pos_max"]
        wv = np.nan_to_num(w)
        cm = np.maximum.accumulate(wv)
        M["W_proj_col_max"] = float(wv.max())
        M["W_proj_col_last"] = float(wv[-1])
        M["W_proj_monotone"] = bool(np.all(np.diff(wv) >= 0))
        M["W_proj_nonzero_rows"] = int((wv > 0).sum())
        if wv.max() > 0:
            ir = int(np.nonzero(wv >= wv.max() * (1 - 1e-9))[0][0])
            M["W_proj_first_reach_step"] = int(step[ir])
            M["W_proj_first_reach_t"] = float(t[ir])
        steps_up = np.nonzero(np.diff(cm, prepend=0.0) > 0)[0]
        M["W_proj_cummax_steps"] = [dict(step=int(step[i]), t=float(t[i]), cummax=float(cm[i])) for i in steps_up[:10]]
        M["W_proj_cummax_n_steps"] = int(len(steps_up))
    # --- энергия ---
    _, en = read_csv_dedup(os.path.join(d, "energy_log.csv"))
    E0 = en["E_mech"][0]
    M["E0"] = float(E0)
    M["energy_rows"] = int(len(en["step"]))
    dcf = en["D_contact_fiber"] / E0
    res = en["residual"] / E0
    M["Dcf_end"] = float(dcf[-1])
    M["Dcf_max"] = float(np.nanmax(dcf))
    M["Dcf_min_inc"] = float(np.nanmin(np.diff(dcf))) if len(dcf) > 1 else 0.0
    M["res_end"] = float(res[-1])
    M["res_max"] = float(np.nanmax(res))
    M["res_min"] = float(np.nanmin(res))
    M["res_absmax"] = float(np.nanmax(np.abs(res)))
    M["res_max_jump"] = float(np.nanmax(np.abs(np.diff(res)))) if len(res) > 1 else 0.0
    M["Dfric_end"] = float(en["D_friction"][-1] / E0)
    M["KEs_end"] = float(en["KE_solid"][-1] / E0)
    M["KEf_end"] = float(en["KE_fiber"][-1] / E0)
    M["Emech_end"] = float(en["E_mech"][-1] / E0)
    M["Emech_max"] = float(np.nanmax(en["E_mech"]) / E0)
    # --- цена ---
    wall = shot.get("wall_s")
    M["sec_per_step"] = float(wall) / M["last_step"] if wall else None
    return M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or a.root
    os.makedirs(out, exist_ok=True)
    runs = find_runs(a.root)
    R = {r: run_metrics(a.root, r) for r in runs}
    # критерии ТЗ
    res3 = R["C3"]["res_absmax"] if "C3" in R else None
    for r, M in R.items():
        s = M["shot"]
        c = {
            "dff<=0.3": s["max_delta_ff/d"] <= 0.3,
            "dfs<=0.3": s["max_delta_fs/d"] <= 0.3,
            "n_tunnel<=5": s["n_tunnel"] <= 5,
            "Dcf_nondecr": M["Dcf_min_inc"] >= DCF_DROP,
            "res<=res(C3)+0.005": (M["res_absmax"] <= res3 + RES_MARGIN) if res3 is not None else None,
            "res_jump<=1e-3": M["res_max_jump"] <= E_JUMP,
        }
        M["criteria"] = {k: (bool(v) if v is not None else None) for k, v in c.items()}
    if "C1" in R and "C2" in R and R["C1"]["sec_per_step"] and R["C2"]["sec_per_step"]:
        cost = R["C1"]["sec_per_step"] / R["C2"]["sec_per_step"]
    else:
        cost = None
    with open(os.path.join(out, "metrics.json"), "w") as fh:
        json.dump(dict(runs=R, cost_ratio_C1_over_C2=cost), fh, indent=1, ensure_ascii=False)

    L = []
    def table(title, cols, rows):
        L.append(f"\n### {title}\n")
        L.append("| прогон | " + " | ".join(cols) + " |")
        L.append("|" + "---|" * (len(cols) + 1))
        for r, vals in rows:
            L.append(f"| {r} | " + " | ".join(vals) + " |")
    table("Итог из shot.txt",
          ["outcome", "V_R", "t_end, мс", "шагов", "n_tunnel", "δ_ff/d max", "δ_fs/d max", "wall, ч", "t_first_break, мс"],
          [(r, [fnum(M["shot"]["outcome"]), fnum(M["shot"]["V_R"]), fnum(M["shot"]["t_end"] * 1e3),
                fnum(M["shot"]["steps"], "{:.0f}"), fnum(M["shot"]["n_tunnel"], "{:.0f}"),
                fnum(M["shot"]["max_delta_ff/d"]), fnum(M["shot"]["max_delta_fs/d"]),
                fnum(M["shot"]["wall_s"] / 3600 if M["shot"]["wall_s"] else None),
                fnum(M["shot"]["t_first_break"] * 1e3 if M["shot"]["t_first_break"] is not None else None)])
           for r, M in R.items()])
    table("Разрывы (monitor.csv)", ["t первого разрыва (monitor), мс", "шаг", "n_broken_v (конец)", "n_broken_h (конец)"],
          [(r, [fnum(M["t_first_break_monitor"] * 1e3 if M["t_first_break_monitor"] is not None else None),
                fnum(M["step_first_break_monitor"], "{:.0f}"), str(M["n_broken_v_end"]), str(M["n_broken_h_end"])])
           for r, M in R.items()])
    pr = {r: M for r, M in R.items() if "proj_resid_max" in M}
    table("Проекция: невязка и W_proj_pos_max",
          ["proj_resid max", "proj_resid медиана", "proj_resid_rel max", "proj_resid_rel медиана", "warm_frac min",
           "W_proj max (колонка)", "W_proj монотонна", "строк с W>0", "первый шаг достижения", "t, мс", "shot: W_proj_pos_max"],
          [(r, [fnum(M["proj_resid_max"]), fnum(M["proj_resid_median"]), fnum(M["proj_resid_rel_max"]),
                fnum(M["proj_resid_rel_median"]), fnum(M["proj_warm_frac_min"]), fnum(M["W_proj_col_max"]),
                fnum(M["W_proj_monotone"]).replace("pass", "да").replace("FAIL", "нет"), str(M["W_proj_nonzero_rows"]),
                fnum(M.get("W_proj_first_reach_step"), "{:.0f}"),
                fnum(M["W_proj_first_reach_t"] * 1e3 if "W_proj_first_reach_t" in M else None),
                fnum(M["shot"]["W_proj_pos_max"])]) for r, M in pr.items()])
    L.append("\n### Ступени накопленного максимума W_proj_pos_max (первые 10)\n")
    for r, M in pr.items():
        st = M["W_proj_cummax_steps"]
        L.append(f"- {r} (всего ступеней {M['W_proj_cummax_n_steps']}): " +
                 "; ".join(f"шаг {x['step']}, t={x['t']*1e3:.4f} мс, {x['cummax']:.3e}" for x in st))
    table("Энергия (доли E0)",
          ["E0, Дж", "E_mech конец", "E_mech max", "D_cf конец", "D_cf max", "min прирост D_cf", "residual конец",
           "residual max", "residual min", "max |Δresidual|", "D_friction конец", "KE_solid конец", "KE_fiber конец"],
          [(r, [fnum(M["E0"]), fnum(M["Emech_end"]), fnum(M["Emech_max"]), fnum(M["Dcf_end"]), fnum(M["Dcf_max"]),
                fnum(M["Dcf_min_inc"]), fnum(M["res_end"]), fnum(M["res_max"]), fnum(M["res_min"]),
                fnum(M["res_max_jump"]), fnum(M["Dfric_end"]), fnum(M["KEs_end"]), fnum(M["KEf_end"])])
           for r, M in R.items()])
    table("Цена",
          ["wall, с", "последний шаг", "с/шаг", "dt_mean, нс"],
          [(r, [fnum(M["shot"]["wall_s"]), str(M["last_step"]), fnum(M["sec_per_step"]),
                fnum(M["shot"]["dt_mean"] * 1e9 if M["shot"]["dt_mean"] else None)]) for r, M in R.items()])
    L.append(f"\nОтношение с/шаг C1/C2 = {fnum(cost)} (прогоны шли параллельно по 8 на 8 ядрах, 1 поток каждый; "
             "рестарты из чекпоинтов: wall_s из shot.txt может учитывать не всё время — см. оговорку в отчёте).")
    ck = list(next(iter(R.values()))["criteria"])
    table("Критерии ТЗ (pass/FAIL)", ck,
          [(r, [fnum(M["criteria"][k]) if M["criteria"][k] is not None else "-" for k in ck]) for r, M in R.items()])
    L.append(f"\nПороги: δ_ff/d ≤ 0,3; δ_fs/d ≤ 0,3; n_tunnel ≤ 5; прирост D_contact_fiber между строками ≥ {DCF_DROP}·E0; "
             f"max|residual|/E0 ≤ max|residual|(C3)/E0 + {RES_MARGIN}"
             + (f" (C3: {res3:.4g})" if res3 is not None else "") + f"; max|Δresidual| ≤ {E_JUMP}·E0.")
    open(os.path.join(out, "tables.md"), "w").write("# Метрики серии закон 3 / проекция\n" + "\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
