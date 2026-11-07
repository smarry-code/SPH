#!/usr/bin/env python3
"""Подбор окна непробитие / пробитие Guo 2019, 1T base, делением диапазона [vmin; vmax] на трети.

Постановка (решение владельца, 2026-10-18): закон контакта нитей 4 (Герц + демпфер), E = 78 ГПа,
mu = 0.1, скорости 71..81 м/с (статья: при 73 не пробит, V_BL = 76 +3/-3).
Критерий «не пробит» и запуск выстрела -- из window_pinball.py (V_R = 0, верт. <= 8, гориз. = 0;
«не определён» -> повтор до --t-max-long; готовые выстрелы не пересчитываются, идущие -- дожидаются).

Алгоритм (lo -- наибольшая «НП» ниже hi, hi -- наименьшая «П»):
  1. пара vmin / --v-first (71 / 76) параллельно;
  2. vmin «П" -> непробития в диапазоне нет -> конец;
     пробития ещё нет (hi нет) -> пара: середина [lo; vmax] и vmax; vmax «НП» -> пробития в диапазоне нет -> конец;
     есть lo и hi: пока hi - lo > --tol -> две точки на третях [lo; hi] (шаг 0.5 м/с) параллельно,
     одна новая точка -> один выстрел на --threads-single.
  Немонотонность («НП» выше «П») записывается; окно берётся по наименьшей «П» и наибольшей «НП» под ней.
"""
import argparse, concurrent.futures, glob, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import window_pinball as wp  # noqa: E402


def r05(x):
    return round(x * 2.0) / 2.0


def bounds(res):
    p = [v for v, q in res.items() if q is True]
    hi = min(p) if p else None
    np_ = [v for v, q in res.items() if q is False and (hi is None or v < hi)]
    lo = max(np_) if np_ else None
    nonmono = sorted(v for v, q in res.items() if q is False and hi is not None and v > hi)
    return lo, hi, nonmono


def threads_in_use():
    """Сумма --threads= у идущих guo2019_shot (свои и чужие драйверы, ручные запуски)."""
    total = 0
    for p in glob.glob("/proc/[0-9]*/cmdline"):
        try:
            with open(p, "rb") as f:
                argv = f.read().split(b"\0")
        except OSError:
            continue
        if not argv or not argv[0].endswith(b"guo2019_shot"):
            continue
        for a in argv:
            if a.startswith(b"--threads="):
                try:
                    total += int(a[len(b"--threads="):])
                except ValueError:
                    pass
    return total


def pick_threads(args, n_new):
    """Потоков на новый выстрел: свободные ядра (--cores минус занятые) поровну на n_new, не меньше --threads."""
    if not args.adaptive:
        return args.threads if n_new > 1 else args.threads_single
    free = max(args.cores - threads_in_use(), 0)
    return max(args.threads, min(args.cores, free // max(n_new, 1)))


def shoot_once(args, mu, v0, threads, rec):
    """Выстрел без повтора до t_max_long (--no-long): «не определён» возвращается как None."""
    if not args.no_long:
        return wp.shoot(args, mu, v0, threads, rec)
    r = wp.run_shot(args, mu, v0, threads, args.t_max)
    p, why = wp.verdict(r)
    with wp.LOCK:
        rec["shots"].append({"mu": mu, "V0": v0, "run": "base", "verdict": p, "why": why, "res": r})
    return p


def run_points(args, mu, pts, rec, res):
    pts = [v for v in pts if v not in res]
    if not pts:
        return False
    th = pick_threads(args, len(pts))
    if len(pts) == 1:
        res[pts[0]] = shoot_once(args, mu, pts[0], th, rec)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(pts)) as ex:
            fut = {v: ex.submit(shoot_once, args, mu, v, th, rec) for v in pts}
            for v, f in fut.items():
                res[v] = f.result()
    return True


def process(args, mu, rec, save):
    res = {}
    run_points(args, mu, [args.vmin, args.v_first], rec, res)
    save()
    rounds = 1
    while True:
        lo, hi, nonmono = bounds(res)
        if nonmono:
            rec["nonmonotonic"] = "НП выше наименьшей П (%g): %s" % (hi, ", ".join("%g" % v for v in nonmono))
        if res.get(args.vmin) is True:
            rec["result"] = "пробит уже при %g -> непробития в [%g; %g] нет" % (args.vmin, args.vmin, args.vmax)
            return "no_np"
        if lo is None:
            rec["result"] = "непробитие не установлено (исходы: %s) -- вывод невозможен" % (
                ", ".join("%g: %s" % (v, wp.word(q)) for v, q in sorted(res.items())))
            return "undetermined"
        if hi is None:
            if res.get(args.vmax) is False:
                rec["result"] = "не пробит до %g включительно -> пробития в [%g; %g] нет" % (
                    args.vmax, args.vmin, args.vmax)
                return "no_p"
            pts = [args.vmax] if args.vmax - lo <= args.tol else [r05(lo + (args.vmax - lo) / 2.0), args.vmax]
        else:
            if hi - lo <= args.tol:
                rec["window"] = [lo, hi]
                rec["result"] = "окно [%g; %g]" % (lo, hi)
                return "window"
            if args.halves:
                # половинное деление: одна точка в середине [lo; hi] (выстрел на --threads-single)
                pts = sorted({r05(0.5 * (lo + hi))} - {lo, hi})
            else:
                pts = sorted({r05(lo + (hi - lo) / 3.0), r05(lo + 2.0 * (hi - lo) / 3.0)} - {lo, hi})
        if rounds >= args.max_rounds:
            rec["window"] = [lo, hi]
            rec["result"] = "достигнут предел раундов (%d): окно [%s; %s] м/с" % (args.max_rounds, lo, hi)
            return "window" if hi is not None else "undetermined"
        rounds += 1
        if not run_points(args, mu, pts, rec, res):
            rec["window"] = [lo, hi] if hi is not None else [lo, None]
            rec["result"] = "новых точек нет, окно [%s; %s]" % (lo, hi)
            return "window" if hi is not None else "undetermined"
        save()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--fib", required=True)
    ap.add_argument("--fiber-E", type=float, default=78e9)
    ap.add_argument("--eps-break", type=float, default=0.05)
    ap.add_argument("--contact", type=int, default=4)
    ap.add_argument("--config", default="1T", choices=["1T", "2T"], help="2T: вердикт по шару, без допуска на разрывы")
    ap.add_argument("--level", default="base")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--threads-single", type=int, default=8)
    ap.add_argument("--energy-log", type=int, default=100)
    ap.add_argument("--save-t", type=float, default=20e-6)
    ap.add_argument("--t-max", type=float, default=1.5e-3)
    ap.add_argument("--t-max-long", type=float, default=3.0e-3)
    ap.add_argument("--mu", type=float, default=0.1)
    ap.add_argument("--vmin", type=float, default=71.0)
    ap.add_argument("--vmax", type=float, default=81.0)
    ap.add_argument("--v-first", type=float, default=76.0)
    ap.add_argument("--tol", type=float, default=2.0)
    ap.add_argument("--max-rounds", type=int, default=99, help="предел раундов (пара или одиночный выстрел = 1 раунд)")
    ap.add_argument("--no-long", action="store_true", help="без повтора до --t-max-long; «не определён» -> точка пропускается")
    ap.add_argument("--adaptive", action="store_true", help="потоки нового выстрела = свободные ядра / число выстрелов")
    ap.add_argument("--cores", type=int, default=8)
    ap.add_argument("--halves", action="store_true", help="половинное деление окна (одна точка в середине), а не трети")
    ap.add_argument("--checkpoint-every", type=int, default=0,
                    help="чекпоинт каждые N шагов и рестарт оборванных выстрелов из него (0 -- без чекпоинтов)")
    ap.add_argument("--extra", default="", help="доп. ключи guo2019_shot через пробел (например --impulse-push-tau=2e-6)")
    ap.add_argument("--verdict-ball", action="store_true", help="1T: вердикт по шару, без допуска на разрывы (как 2T)")
    args = ap.parse_args()
    wp.BALL_ONLY = args.verdict_ball
    os.makedirs(args.work, exist_ok=True)
    for path in (args.bin, args.fib):
        if not os.path.exists(path):
            print("нет файла: " + path, file=sys.stderr)
            return 1
    wp.log(args.work, "серия (деление): mu=%g V0=[%g; %g] первая пара %g/%g tol=%g E=%g contact=%d fib=%s" % (
        args.mu, args.vmin, args.vmax, args.vmin, args.v_first, args.tol, args.fiber_E, args.contact, args.fib))
    if args.extra:
        wp.log(args.work, "доп. ключи: " + args.extra)
    rec = {"mu": args.mu, "shots": []}
    final = ""
    save = lambda: wp.write_window(args, [rec], final)
    try:
        res = process(args, args.mu, rec, save)
        final = "mu = %g: %s (%s)" % (args.mu, rec.get("result"), res)
    except Exception as ex:
        final = "Серия ОСТАНОВЛЕНА: %s" % ex
        wp.log(args.work, final)
        save()
        return 1
    wp.log(args.work, final)
    save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
