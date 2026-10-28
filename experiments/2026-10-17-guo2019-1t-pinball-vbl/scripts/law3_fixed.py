#!/usr/bin/env python3
"""Закон 3 + выталкивание на фиксированных скоростях (решение владельца 2026-10-21).

Условия как у этапа B (не пробит при 73 м/с): закон 3 на нить-нить и нить-шар, tau = 0,5 мкс,
порог 0,1 r, шаг x0,5, E = 80 ГПа, mu = 0,1, нить twaron_plastic_s25 (без RATE), Solid-Solid как в
сцене (pinball). Кадры каждые 20 мкс, лог энергии.
Порядок: --pair параллельно по --threads, затем --single на --threads-single. Готовые выстрелы
(есть shots.csv) не пересчитываются. Сводка -- <work>/summary.md.
"""
import argparse, concurrent.futures, glob, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import window_pinball as wp  # noqa: E402

PAPER = "статья (1T): при 73 м/с не пробит, V_BL = 76 (+3/-3) м/с; 3 верт. / 0 гориз. разрывов; прогиб к первому разрыву 25-27 мм"


def line(v, r):
    n = wp.num
    return "| %g | %s | %.2f | %s / %s | %.0f мкс / %.1f мм | %.3f / %.3f | %s | %.3f | %.1f ч |" % (
        v, r.get("outcome"), n(r.get("V_R"), 0.0), r.get("n_broken_v"), r.get("n_broken_h"),
        1e6 * n(r.get("t_first_break")), 1e3 * n(r.get("w_f_at_first_break")), n(r.get("max_delta_ff/d")),
        n(r.get("max_delta_fs/d")), r.get("n_tunnel"), n(r.get("W_p")), n(r.get("wall_s")) / 3600.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--fib", required=True)
    ap.add_argument("--fiber-E", type=float, default=80e9)
    ap.add_argument("--eps-break", type=float, default=0.05)
    ap.add_argument("--mu", type=float, default=0.1)
    ap.add_argument("--pair", default="71,76")
    ap.add_argument("--single", default="78")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--threads-single", type=int, default=8)
    ap.add_argument("--t-max", type=float, default=1.5e-3)
    ap.add_argument("--save-t", type=float, default=20e-6)
    ap.add_argument("--energy-log", type=int, default=100)
    # Без выталкивания по умолчанию (решение владельца 2026-10-21: --impulse-push-tau -- источник энергии).
    ap.add_argument("--extra", default="--dt-scale=0.5")
    args = ap.parse_args()
    args.contact = 3
    args.level = "base"
    os.makedirs(args.work, exist_ok=True)
    res = {}

    def save():
        # В сводку -- все готовые выстрелы в каталоге (несколько экземпляров драйвера в одной папке).
        allres = {}
        for sd in glob.glob(os.path.join(args.work, "mu%g_v*" % args.mu)):
            if not os.path.isdir(sd):
                continue
            r = wp.read_result(sd)
            if r is not None:
                allres[float(os.path.basename(sd).split("_v")[1])] = r
        allres.update(res)
        rows = "\n".join(line(v, allres[v]) for v in sorted(allres))
        with open(os.path.join(args.work, "summary.md"), "w", encoding="utf-8") as f:
            f.write("# 1T, закон 3, E = %g, mu = %g (ключи расчёта -- ниже)\n\n%s\n\n"
                    "ключи: %s\n\n| V0, м/с | исход | V_R, м/с | разрывы v/h | первый разрыв / прогиб | max delta/d ff/fs "
                    "| n_tunnel | W_p, Дж | счёт |\n|---|---|---|---|---|---|---|---|---|\n%s\n" % (
                        args.fiber_E, args.mu, PAPER, args.extra, rows))

    wp.log(args.work, "серия: пара %s, затем %s; %s" % (args.pair, args.single, args.extra))
    pair = [float(x) for x in args.pair.split(",") if x]
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(len(pair), 1)) as ex:
        fut = {v: ex.submit(wp.run_shot, args, args.mu, v, args.threads, args.t_max) for v in pair}
        for v, f in fut.items():
            res[v] = f.result()
            save()
    for v in [float(x) for x in args.single.split(",") if x]:
        res[v] = wp.run_shot(args, args.mu, v, args.threads_single, args.t_max)
        save()
    wp.log(args.work, "серия готова")
    return 0


if __name__ == "__main__":
    sys.exit(main())
