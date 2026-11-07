#!/usr/bin/env python3
"""Закон 3 (импульсный) с выталкиванием: проверка п. 1 и п. 4, затем подбор окна непробитие / пробитие.

Постановка (решение владельца, 2026-10-20): Guo 2019 1T base, закон 3 на обоих контактах нитей,
mu = 0.1, E = 80 ГПа, нить twaron_plastic_s25 (без RATE), скорости 70..81 м/с
(статья: при 73 не пробит, V_BL = 76 +3/-3). Кадры каждые 20 мкс и лог энергии включены.

Этап A (п. 1, выталкивание): 73 м/с, два выстрела параллельно (по --threads) с tau из --taus.
  Проходит tau, при котором n_tunnel = 0 и max delta/d (ff и fs) <= --delta-max. Берётся
  наибольший проходящий tau (самое мягкое выталкивание). Не прошёл ни один -> стоп.
Этап B (п. 4, сходимость по шагу): тот же выстрел с --dt-scale=0.5 на --threads-single.
  Проходит, если исход тот же, |V_R - V_R(A)| <= --dvr-max и число разрывов отличается <= 2.
  Не проходит -> стоп: окно зависело бы от шага.
Этап C (окно): window_bisect.py с выбранным tau; выстрел 73 м/с из этапа A переиспользуется.
Итоги этапов -- в <work>/program.md; лог -- <work>/program.log.
"""
import argparse, concurrent.futures, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import window_pinball as wp  # noqa: E402


def num(r, k):
    return wp.num(r.get(k))


def summary(r):
    return ("исход %s, V_R %.2f, разрывы %s/%s, первый разрыв %.0f мкс, прогиб %.1f мм, max delta/d ff %.3f fs %.3f, "
            "n_tunnel %s, W_p %.3f Дж, dt_mean %.3g, счёт %.1f ч" % (
                r.get("outcome"), num(r, "V_R"), r.get("n_broken_v"), r.get("n_broken_h"),
                1e6 * num(r, "t_first_break"), 1e3 * num(r, "w_f_at_first_break"), num(r, "max_delta_ff/d"),
                num(r, "max_delta_fs/d"), r.get("n_tunnel"), num(r, "W_p"), num(r, "dt_mean"),
                num(r, "wall_s") / 3600.0))


def shot_args(args, work, extra):
    a = argparse.Namespace(**vars(args))
    a.work = work
    a.extra = extra
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--fib", required=True)
    ap.add_argument("--fiber-E", type=float, default=80e9)
    ap.add_argument("--eps-break", type=float, default=0.05)
    ap.add_argument("--mu", type=float, default=0.1)
    ap.add_argument("--v-check", type=float, default=73.0)
    ap.add_argument("--taus", default="2e-6,5e-7")
    ap.add_argument("--depth", type=float, default=0.1)
    ap.add_argument("--delta-max", type=float, default=0.3)
    ap.add_argument("--dvr-max", type=float, default=3.0)
    ap.add_argument("--vmin", type=float, default=70.0)
    ap.add_argument("--vmax", type=float, default=81.0)
    ap.add_argument("--v-first", type=float, default=73.0)
    ap.add_argument("--tol", type=float, default=2.0)
    ap.add_argument("--max-rounds", type=int, default=4)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--threads-single", type=int, default=8)
    ap.add_argument("--t-max", type=float, default=1.5e-3)
    ap.add_argument("--save-t", type=float, default=20e-6)
    ap.add_argument("--energy-log", type=int, default=100)
    args = ap.parse_args()
    args.contact = 3
    args.level = "base"
    args.t_max_long = args.t_max
    W = args.work
    os.makedirs(W, exist_ok=True)
    md = []

    def note(line):
        md.append(line)
        wp.log(W, line)
        with open(os.path.join(W, "program.md"), "w", encoding="utf-8") as f:
            f.write("# Закон 3 + выталкивание: п. 1, п. 4, окно\n\n" + "\n".join("- " + x for x in md) + "\n")

    note("параметры: mu=%g E=%g V_check=%g taus=%s depth=%g; окно [%g; %g], первая пара %g/%g, tol=%g, раундов <= %d" % (
        args.mu, args.fiber_E, args.v_check, args.taus, args.depth, args.vmin, args.vmax, args.vmin, args.v_first,
        args.tol, args.max_rounds))

    # ---- Этап A
    taus = [float(x) for x in args.taus.split(",")]
    resA = {}

    def runA(tau):
        a = shot_args(args, os.path.join(W, "A_tau%g" % tau),
                      "--impulse-push-tau=%g --impulse-push-depth=%g" % (tau, args.depth))
        os.makedirs(a.work, exist_ok=True)
        return wp.run_shot(a, args.mu, args.v_check, args.threads, args.t_max)

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(taus)) as ex:
        fut = {t: ex.submit(runA, t) for t in taus}
        for t, f in fut.items():
            resA[t] = f.result()
    ok = []
    for t in taus:
        r = resA[t]
        good = (wp.num(r.get("n_tunnel"), 1) == 0 and num(r, "max_delta_ff/d") <= args.delta_max
                and num(r, "max_delta_fs/d") <= args.delta_max)
        note("A, tau = %g с: %s -> %s" % (t, summary(r), "ПРОХОДИТ" if good else "не проходит"))
        if good:
            ok.append(t)
    if not ok:
        note("ИТОГ: п. 1 не выполнен ни при одном tau -> серия остановлена, нужен разбор")
        return 2
    tau = max(ok)
    rA = resA[tau]
    note("выбран tau = %g с" % tau)
    extra = "--impulse-push-tau=%g --impulse-push-depth=%g" % (tau, args.depth)

    # ---- Этап B
    aB = shot_args(args, os.path.join(W, "B_tau%g_dt0.5" % tau), extra + " --dt-scale=0.5")
    os.makedirs(aB.work, exist_ok=True)
    rB = wp.run_shot(aB, args.mu, args.v_check, args.threads_single, args.t_max)
    nA = wp.num(rA.get("n_broken_v"), 0) + wp.num(rA.get("n_broken_h"), 0)
    nB = wp.num(rB.get("n_broken_v"), 0) + wp.num(rB.get("n_broken_h"), 0)
    dvr = abs(num(rB, "V_R") - num(rA, "V_R"))
    goodB = rB.get("outcome") == rA.get("outcome") and dvr <= args.dvr_max and abs(nB - nA) <= 2
    note("B, dt x 0.5: %s; |dV_R| = %.2f м/с, разрывов %d против %d -> %s" % (
        summary(rB), dvr, nB, nA, "ПРОХОДИТ" if goodB else "НЕ проходит"))
    if not goodB:
        note("ИТОГ: результат зависит от шага -> подбор окна не запускается, нужен разбор")
        return 3

    # ---- Этап C
    wC = os.path.join(W, "C_window_tau%g" % tau)
    os.makedirs(wC, exist_ok=True)
    dst = wp.shot_dir(wC, args.mu, args.v_check)
    if not os.path.exists(dst):
        shutil.copytree(wp.shot_dir(os.path.join(W, "A_tau%g" % tau), args.mu, args.v_check), dst)
        src_log = wp.shot_dir(os.path.join(W, "A_tau%g" % tau), args.mu, args.v_check) + ".log"
        if os.path.exists(src_log):
            shutil.copy(src_log, dst + ".log")
    note("C: окно, выстрел %g м/с взят из этапа A" % args.v_check)
    cmd = [sys.executable, os.path.join(HERE, "window_bisect.py"), "--work=" + wC, "--bin=" + args.bin,
           "--fib=" + args.fib, "--fiber-E=%g" % args.fiber_E, "--eps-break=%g" % args.eps_break, "--contact=3",
           "--mu=%g" % args.mu, "--vmin=%g" % args.vmin, "--vmax=%g" % args.vmax, "--v-first=%g" % args.v_first,
           "--tol=%g" % args.tol, "--max-rounds=%d" % args.max_rounds, "--no-long", "--adaptive", "--cores=8",
           "--threads=%d" % args.threads, "--threads-single=%d" % args.threads_single, "--t-max=%g" % args.t_max,
           "--save-t=%g" % args.save_t, "--energy-log=%d" % args.energy_log, "--extra=" + extra]
    rc = subprocess.run(cmd).returncode
    final = ""
    p = os.path.join(wC, "window.md")
    if os.path.exists(p):
        final = open(p, encoding="utf-8").read().strip().splitlines()[-1]
    note("ИТОГ C (rc=%d): %s; подробности %s" % (rc, final, p))
    return rc


if __name__ == "__main__":
    sys.exit(main())
