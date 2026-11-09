#!/usr/bin/env python3
"""Подбор числа слоёв ткани k (конфигурация kT) при шаре 1 г и V0 = 550 м/с (Guo 2019, закон 3 + выталкивание).

Решение владельца (2026-11-09). P -- наибольшее k «пробит», N -- наименьшее k «не пробит» (монотонность предполагается,
нарушения записываются). Ответ: k = P (kT пробивается, (k+1)T -- нет).
  1. пара --k-first (4T / 8T) параллельно;
  2. N известно, P нет: один выстрел на max(1, N // 2); НП -> один выстрел на 1T; 1T «НП» -> «k < 1», конец;
     пробит -> интервал [P; N];
  3. P известно, N нет (все пробиты): пары (P + k_step, P + 2 k_step) параллельно, пока не появится «НП»;
     выход за --k-max -> ошибка;
  4. пока N - P > 1: две точки на третях [P; N] (округление, без дублей и концов) параллельно; одна новая -- один выстрел.
«Не определён» после продолжения до --t-max-long -> серия останавливается (rc 1), угадывания нет.
Результаты: <work>/search.md и search.json после каждого раунда. Выстрелы -- window_pinball.py (очередь, чекпоинты).
--dry-run --fake-kstar=K: выстрелов нет, исход k <= K -> «пробит» (проверка алгоритма).
"""
import argparse, concurrent.futures, copy, json, math, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import window_pinball as wp  # noqa: E402


def rnd(x):
    return int(math.floor(x + 0.5))


def bounds(res):
    """(P, N, немонотонные): N -- наименьшее НП, P -- наибольшее П ниже N."""
    fl = [k for k, q in res.items() if q is False]
    n = min(fl) if fl else None
    pl = [k for k, q in res.items() if q is True and (n is None or k < n)]
    p = max(pl) if pl else None
    nonmono = sorted(k for k, q in res.items() if q is True and n is not None and k > n)
    return p, n, nonmono


def next_points(res, args):
    """Следующий раунд. ('go', [k...]) | ('done', ответ_k) | ('below', None) | ('error', текст)."""
    p, n, _ = bounds(res)
    if p is None and n is None:
        return "go", sorted(int(x) for x in args.k_first.split(","))
    if p is None:
        if n <= 1:
            return "below", None
        for k in (max(1, n // 2), 1):
            if k not in res:
                return "go", [k]
        return "below", None
    if n is None:
        if p >= args.k_max:
            return "error", "все k до --k-max = %d пробиты, N не найдено" % args.k_max
        pts = [k for k in (p + args.k_step, p + 2 * args.k_step) if k <= args.k_max]
        pts = [k for k in pts if k not in res]
        return "go", pts
    if n - p <= 1:
        return "done", p
    a, b = rnd(p + (n - p) / 3.0), rnd(p + 2.0 * (n - p) / 3.0)
    return "go", sorted({a, b} - {p, n} - set(res))


def real_shot(args, k):
    """Выстрел kT: возвращает (вердикт, запись). Копия args с config = kT."""
    a = copy.copy(args)
    a.config = "%dT" % k
    rec = {"shots": []}
    t0 = time.time()
    if args.no_long:
        r = wp.run_shot(a, a.mu, a.v0, a.threads, a.t_max)
        p, why = wp.verdict(r)
        rec["shots"].append({"run": "base", "verdict": p, "why": why, "res": r})
    else:
        wp.shoot(a, a.mu, a.v0, a.threads, rec)
    e = rec["shots"][-1]
    return e["verdict"], {"k": k, "verdict": e["verdict"], "why": e["why"], "run": e["run"], "res": e["res"],
                          "wall_driver_s": round(time.time() - t0)}


def fake_shot(args, k):
    p = k <= args.fake_kstar
    r = {"V_R": "100" if p else "0", "t_end": "4e-4", "wall_s": "0", "n_broken_v": "0", "n_broken_h": "0"}
    return p, {"k": k, "verdict": p, "why": "фиктивный исход (kstar=%d)" % args.fake_kstar, "run": "fake", "res": r,
               "wall_driver_s": 0}


def write_report(args, entries, rounds, final, notes):
    res = {e["k"]: e["verdict"] for e in entries}
    p, n, nonmono = bounds(res)
    keys = sorted({q for e in entries for q in e["res"] if q.startswith("n_broken")})
    L = ["# Подбор числа слоёв k: шар 1 г, V0 = %g м/с, mu = %g, закон контакта %d" % (
        args.v0, args.mu, args.contact), "",
        "Нить %s, E = %g Па, eps_break = %g; t_max = %g (повтор %g). Вердикт по шару (V_R > 0 -- пробит)." % (
            os.path.basename(args.fib or ""), args.fiber_E, args.eps_break, args.t_max, args.t_max_long), "",
        "| k | вердикт | причина | V_R | t_end, мкс | %s | wall, с | прогон |" % " | ".join(keys or ["n_broken"]),
        "|---|---|---|---|---|%s---|---|" % ("---|" * max(len(keys), 1))]
    for e in sorted(entries, key=lambda x: x["k"]):
        r = e["res"]
        L.append("| %d | %s | %s | %.2f | %.0f | %s | %.0f | %s |" % (
            e["k"], wp.word(e["verdict"]), e["why"], wp.num(r.get("V_R")), wp.num(r.get("t_end"), 0.0) * 1e6,
            " | ".join(str(r.get(q, "")) for q in keys) if keys else "", wp.num(r.get("wall_s")), e["run"]))
    L += ["", "- P (наибольшее пробитое) = %s, N (наименьшее непробитое) = %s" % (p, n),
          "- раунды: " + "; ".join("+".join("%dT" % k for k in rd) for rd in rounds)]
    if nonmono:
        L.append("- немонотонно: пробиты при k выше наименьшего НП (%s): %s" % (n, nonmono))
    for t in notes:
        L.append("- " + t)
    L += ["", "## Итог", "", final or "в работе"]
    with open(os.path.join(args.work, "search.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    with open(os.path.join(args.work, "search.json"), "w", encoding="utf-8") as f:
        json.dump({"P": p, "N": n, "nonmonotonic": nonmono, "rounds": rounds, "final": final, "notes": notes,
                   "shots": entries}, f, ensure_ascii=False, indent=1, default=str)


def run_search(args):
    """Возвращает (rc, итоговая строка, ответ k или None, раунды)."""
    shot = fake_shot if args.dry_run else real_shot
    entries, rounds, notes, res = [], [], [], {}
    final = ""
    save = lambda: write_report(args, entries, rounds, final, notes)
    while True:
        kind, val = next_points(res, args)
        if kind == "done":
            final = "Ответ: k = %d (%dT пробивается, %dT -- нет)." % (val, val, val + 1)
            save()
            return 0, final, val, rounds
        if kind == "below":
            final = "Ответ: k < 1 (1T не пробивается)."
            save()
            return 0, final, 0, rounds
        if kind == "error":
            final = "Серия ОСТАНОВЛЕНА: " + val
            save()
            return 1, final, None, rounds
        pts = val
        if not pts:
            final = "Серия ОСТАНОВЛЕНА: новых точек нет, а границы не сошлись (res=%s)" % res
            save()
            return 1, final, None, rounds
        rounds.append(pts)
        wp.log(args.work, "раунд %d: %s" % (len(rounds), ", ".join("%dT" % k for k in pts)))
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(pts)) as ex:
            futs = {k: ex.submit(shot, args, k) for k in pts}
            outs = {}
            for k, f in futs.items():  # дождаться всех, даже если один упал
                try:
                    outs[k] = f.result()
                except Exception as exc:
                    outs[k] = exc
        for k in pts:
            o = outs[k]
            if isinstance(o, Exception):
                final = "Серия ОСТАНОВЛЕНА: выстрел %dT: %s" % (k, o)
                wp.log(args.work, final)
                save()
                return 1, final, None, rounds
            res[k] = o[0]
            entries.append(o[1])
            wp.log(args.work, "%dT: %s (%s)" % (k, wp.word(o[0]), o[1]["why"]))
        p, n, nonmono = bounds(res)
        if nonmono:
            note = "немонотонно: пробиты %s при N = %s" % (nonmono, n)
            if note not in notes:
                notes.append(note)
        und = [k for k in pts if res[k] is None]
        if und:
            final = "не определён при %s -- серия остановлена, вывод невозможен" % ", ".join("%dT" % k for k in und)
            notes.append(final)
            wp.log(args.work, final)
            save()
            return 1, final, None, rounds
        save()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--bin", default="")
    ap.add_argument("--fib", default="")
    ap.add_argument("--fiber-E", type=float, default=8e10)
    ap.add_argument("--eps-break", type=float, default=0.05)
    ap.add_argument("--contact", type=int, default=3)
    ap.add_argument("--level", default="base")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--energy-log", type=int, default=100)
    ap.add_argument("--save-t", type=float, default=20e-6)
    ap.add_argument("--t-max", type=float, default=4e-4)
    ap.add_argument("--t-max-long", type=float, default=1e-3)
    ap.add_argument("--mu", type=float, default=0.27)
    ap.add_argument("--checkpoint-every", type=int, default=0)
    ap.add_argument("--extra", default="")
    ap.add_argument("--verdict-ball", action="store_true")
    ap.add_argument("--stop-vz-frac", type=float, default=0.05,
                    help="«не определён» -> «не пробит», если к t_max |v_z| < доля*V0 и z_c < 0; 0 -- всегда продолжение")
    ap.add_argument("--no-long", action="store_true", help="без продолжения до --t-max-long")
    ap.add_argument("--slot-script", default="~/proj_runs/threads_free.sh")
    ap.add_argument("--proj-status", default="~/proj_runs/proj_status.sh")
    ap.add_argument("--min-uptime", type=float, default=600)
    ap.add_argument("--slot-poll", type=float, default=300)
    ap.add_argument("--slot-verify-s", type=float, default=120)
    ap.add_argument("--v0", type=float, default=550.0)
    ap.add_argument("--k-first", default="4,8")
    ap.add_argument("--k-step", type=int, default=2)
    ap.add_argument("--k-max", type=int, default=40)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--fake-kstar", type=int, default=0)
    args = ap.parse_args()
    wp.BALL_ONLY = args.verdict_ball
    wp.STOP_VZ_FRAC = args.stop_vz_frac
    os.makedirs(args.work, exist_ok=True)
    if not args.dry_run:
        for path in (args.bin, args.fib):
            if not path or not os.path.exists(path):
                print("нет файла: " + path, file=sys.stderr)
                return 1
    wp.log(args.work, "поиск k: V0=%g mu=%g k-first=%s шаг=%d k-max=%d E=%g contact=%d fib=%s%s" % (
        args.v0, args.mu, args.k_first, args.k_step, args.k_max, args.fiber_E, args.contact, args.fib,
        " [DRY-RUN kstar=%d]" % args.fake_kstar if args.dry_run else ""))
    if args.extra:
        wp.log(args.work, "доп. ключи: " + args.extra)
    rc, final, ans, rounds = run_search(args)
    wp.log(args.work, final)
    if args.dry_run:
        print("RESULT kstar=%d раунды: %s -> %s" % (
            args.fake_kstar, " | ".join("/".join("%d" % k for k in rd) for rd in rounds),
            ans if ans is not None else final))
    return rc


if __name__ == "__main__":
    sys.exit(main())
