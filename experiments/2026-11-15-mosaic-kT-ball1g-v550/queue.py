#!/usr/bin/env python3
"""Очередь выстрелов kT на мозаичной ткани, локально (решение владельца 2026-11-15, серверы заняты).

Постановка -- серия kT (experiments/2026-11-09-fabric-kT-ball1g-v550): шар 1 г, D = 6,253 мм, dx = 0,3 мм, V0 = 550 м/с,
закон 3 как проекция скорости (K = 8, тёплый старт), без выталкивания, mu = 0,27, dt-scale 0,5, нить
twaron_s25_dmg0.00415.fib (E = 8e10, eps_b = 5 %). Мозаика: N = 7, seed = 1, поле 50 x 70 мм (свободная часть
50 x 50, зажим |y| >= 25 мм), стежки -- материал стёжки «одеяла» (умолчания guo2019_shot).

Этап 1: 16T вариант A (одна раскройка на весь пакет, кусочек -- стопка из всех слоёв), затем 16T вариант B
(своя раскройка у каждого слоя, k = 1; сквозные стежки по швам всех раскроек пакета).
Выбор варианта -- тот, что держит лучше:
  не пробит лучше пробитого; оба пробиты -- меньше V_R; оба не пробиты -- меньше наибольшая глубина шара.
Этап 2 -- только выбранным вариантом, поиск границы по числу слоёв n делением отрезка пополам:
  нет непробитого n -- следующее n + 8 (16 -> 24 -> 32 ...); нет пробитого -- n - 8; иначе середина (lo + hi) // 2,
  пока hi - lo > 1 (lo -- наибольшее пробитое, hi -- наименьшее непробитое).
Вердикт по шару: perforated -> пробит; not_perforated -> пробит, если V_R > 0, иначе не пробит. undetermined к t_max --
выстрел продолжается из чекпоинта с t_max x 1,5 (до 1,2 мс); если и тогда не определён: |v_z| < 5 % V0 -- «не пробит
(оценка)», иначе «пробит (оценка)».

Состояние -- queue_state.json, журнал -- queue.log. Повторный запуск продолжает: готовые выстрелы не считаются заново,
незаконченный -- продолжается из чекпоинта (--restart).
Запуск: python queue.py   (из каталога эксперимента; в фоне -- run_queue.ps1)
"""
import csv
import glob
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin", "guo2019_shot.exe")
FIB = os.path.join(HERE, "twaron_s25_dmg0.00415.fib")
RUNS = os.path.join(HERE, "runs")
STATE = os.path.join(HERE, "queue_state.json")
LOG = os.path.join(HERE, "queue.log")
V0 = 550.0
THREADS = 6
T_MAX0 = 4e-4
T_MAX_CAP = 1.2e-3
STOP_VZ_FRAC = 0.05
EXTRA = [17]  # доп. выстрелы выбранным вариантом после завершения поиска (решение владельца 2026-11-15: 17T после 15T)
VARIANTS = {"A": [], "B": ["--layers-per-cut=1"]}  # A: --layers-per-cut не задан -- одна раскройка на пакет
COMMON = ["--level=base", "--v0=%g" % V0, "--frames=off", "--save-t=1e-5", "--threads=%d" % THREADS, "--n-mon=10",
          "--n-diag=100", "--repeat=1", "--no-bracket",
          "--random-cut=7", "--cut-seed=1", "--cut-size=0.05x0.07", "--cut-clamp-edge=0.025",
          "--mu=0.27", "--ff-contact=3", "--fs-contact=3", "--ff-contact-predictor=1", "--fs-contact-predictor=1",
          "--dt-scale=0.5", "--ball-d=6.253e-3", "--ball-mass=1e-3", "--ball-dx=3e-4",
          "--fiber-inelastic=" + FIB, "--fiber-E=8e10", "--eps-break=0.05",
          "--energy-log=100", "--checkpoint-every=500"]


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def load_state():
    if os.path.exists(STATE):
        with open(STATE, encoding="utf-8") as f:
            return json.load(f)
    return {"results": {}, "variant": None, "why_variant": None, "done": False, "boundary": None}


def save_state(st):
    tmp = STATE + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(st, f, ensure_ascii=False, indent=1)
    os.replace(tmp, STATE)


def shot_dir(out):
    d = [p for p in glob.glob(os.path.join(out, "*T_base_*_v550.0_r0")) if os.path.isdir(p)]
    return d[0] if d else None


def read_row(out):
    p = os.path.join(out, "shots.csv")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return rows[-1] if rows else None


def last_monitor(sd):
    with open(os.path.join(sd, "monitor.csv"), encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r.get("gap_mean")]
    zmin = min(float(r["z_c"]) for r in rows)
    return rows[-1], zmin


def run_once(n, var, t_max, restart):
    out = os.path.join(RUNS, "%dT_%s" % (n, var))
    os.makedirs(out, exist_ok=True)
    cmd = [BIN, "--config=%dT" % n, "--out=" + out, "--t-max=%g" % t_max] + COMMON + VARIANTS[var]
    if restart:
        cmd.append("--restart")
    env = dict(os.environ, OMP_NUM_THREADS=str(THREADS))
    log("%s %dT вариант %s, t_max = %g: %s" % ("РЕСТАРТ" if restart else "старт", n, var, t_max, " ".join(cmd[1:])))
    t0 = time.time()
    with open(os.path.join(out, "driver.out"), "a", encoding="utf-8") as lg:
        rc = subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT, env=env, cwd=HERE)
    log("  код возврата %d, %.2f ч" % (rc, (time.time() - t0) / 3600.0))
    return out, rc


def shoot(n, var):
    """Один выстрел до определённого исхода. Возвращает словарь результата."""
    out = os.path.join(RUNS, "%dT_%s" % (n, var))
    t_max = T_MAX0
    sd = shot_dir(out)
    restart = bool(sd and os.path.exists(os.path.join(sd, "checkpoint.ckpt")))
    if restart and os.path.exists(os.path.join(out, "t_max.txt")):
        t_max = float(open(os.path.join(out, "t_max.txt")).read())
    fails = 0
    while True:
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "t_max.txt"), "w") as f: # t_max текущего прогона -- для продолжения после сбоя
            f.write("%g" % t_max)
        out, rc = run_once(n, var, t_max, restart)
        sd = shot_dir(out)
        row = read_row(out)
        if rc != 0 or row is None:
            fails += 1
            if fails >= 2 or not (sd and os.path.exists(os.path.join(sd, "checkpoint.ckpt"))):
                raise RuntimeError("%dT %s: выстрел упал (код %d)" % (n, var, rc))
            restart = True
            continue
        outcome = row["outcome"]
        if outcome != "undetermined" or t_max >= T_MAX_CAP:
            break
        t_max = min(T_MAX_CAP, 1.5 * t_max)
        restart = True
        log("  %dT %s: не определён к t_end = %s, продолжаю до t_max = %g" % (n, var, row["t_end"], t_max))
    last, zmin = last_monitor(sd)
    vr = float(row["V_R"]) if row["V_R"] not in ("", "nan", "-") else float("nan")
    vz = float(last["v_z"])
    if outcome == "perforated":
        perf, why = True, "пробит, V_R = %.2f м/с" % vr
    elif outcome == "not_perforated":
        perf = vr > 0.0
        why = ("V_R = %.2f > 0" % vr) if perf else "не пробит: шар остановлен"
    else:
        perf = abs(vz) >= STOP_VZ_FRAC * V0 and vz < 0
        why = "%s (оценка): к t_max v_z = %.2f м/с" % ("пробит" if perf else "не пробит", vz)
    res = {"n": n, "variant": var, "perforated": perf, "why": why, "outcome": outcome, "V_R": vr, "e": row.get("e"),
           "t_end": row.get("t_end"), "z_c_min": zmin, "v_z_end": vz,
           "n_broken_v": row.get("n_broken_v"), "n_broken_h": row.get("n_broken_h"),
           "n_stitch": row.get("n_stitch"), "n_stitch_through": row.get("n_stitch_through"),
           "n_stitch_broken": row.get("n_stitch_broken_end"), "w_f_max": row.get("w_f_max"),
           "max_delta_ff/d": row.get("max_delta_ff/d"), "max_delta_fs/d": row.get("max_delta_fs/d"),
           "E_end/E0": row.get("E_tot_end/E_tot_0"), "wall_s": row.get("wall_s"), "dir": sd}
    log("ИТОГ %dT %s: %s; разрывы нитей %s/%s, стежков %s из %s; z_c min = %.1f мм" % (
        n, var, why, res["n_broken_v"], res["n_broken_h"], res["n_stitch_broken"], res["n_stitch"], zmin * 1e3))
    return res


def get(st, n, var):
    key = "%dT_%s" % (n, var)
    if key not in st["results"]:
        st["results"][key] = shoot(n, var)
        save_state(st)
    return st["results"][key]


def better(a, b):
    """True, если a держит лучше b."""
    if a["perforated"] != b["perforated"]:
        return not a["perforated"]
    if a["perforated"]:
        return a["V_R"] < b["V_R"]
    return a["z_c_min"] > b["z_c_min"]  # z_c < 0: меньше глубина -- больше z_c_min


def extra(st):
    """Дополнительные выстрелы EXTRA выбранным вариантом (на границу не влияют)."""
    v = st["variant"]
    for n in EXTRA:
        if "%dT_%s" % (n, v) not in st["results"]:
            log("доп. выстрел: %dT (вариант %s)" % (n, v))
            get(st, n, v)
    log("ДОП ГОТОВО: %s" % ", ".join("%dT" % n for n in EXTRA))
    return 0


def main():
    st = load_state()
    try:
        if st.get("done"):
            log("поиск уже завершён: %s" % st.get("boundary"))
            return extra(st)
        a = get(st, 16, "A")
        b = get(st, 16, "B")
        if st["variant"] is None:
            v = "A" if better(a, b) else "B"
            st["variant"] = v
            st["why_variant"] = "16T: A -- %s; B -- %s" % (a["why"], b["why"])
            save_state(st)
            log("ВЫБРАН вариант %s (%s)" % (v, st["why_variant"]))
        v = st["variant"]
        while True:
            res = {int(k.split("T_")[0]): r["perforated"] for k, r in st["results"].items() if k.endswith("_" + v)}
            perf = [n for n, p in res.items() if p]
            stop = [n for n, p in res.items() if not p]
            lo = max(perf) if perf else None
            hi = min(stop) if stop else None
            if lo is not None and hi is not None and hi < lo:
                log("ВНИМАНИЕ: немонотонность: пробит %d, не пробит %d -- поиск остановлен" % (lo, hi))
                break
            if lo is not None and hi is not None and hi - lo <= 1:
                break
            if hi is None:
                n = lo + 8
            elif lo is None:
                n = max(1, hi - 8)
                if n == hi:
                    break
            else:
                n = (lo + hi) // 2
            if n > 64:
                log("n > 64 -- поиск остановлен")
                break
            log("следующий выстрел: %dT (вариант %s; пробит до %s, не пробит от %s)" % (n, v, lo, hi))
            get(st, n, v)
        st["done"] = True
        st["boundary"] = {"variant": v, "max_perforated": lo, "min_not_perforated": hi}
        save_state(st)
        log("ГОТОВО: вариант %s, наибольшее пробитое n = %s, наименьшее непробитое n = %s" % (v, lo, hi))
        return extra(st)
    except Exception as exc:
        log("ОСТАНОВЛЕНО: %s" % exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
