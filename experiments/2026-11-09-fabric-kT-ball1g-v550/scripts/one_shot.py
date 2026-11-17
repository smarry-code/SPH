#!/usr/bin/env python3
"""Один выстрел kT без подбора k (claude-compute-8t): решение владельца 2026-11-09, 24T на законе контакта 0 (Герц).

Допуск, проверка занятых потоков и продолжение из чекпоинта после перезагрузки -- window_pinball.run_shot.
Выстрел ждёт FREE >= --threads, поэтому стартует сам, когда освободится сервер. Новых точек не выбирает.
«Не определён» к --t-max не продолжается: финальный чекпоинт остаётся, решение за владельцем.
Итог: <work>/result.txt. Код возврата 0 -- выстрел досчитан (при любом вердикте), 1 -- выстрел упал дважды.
"""
import argparse, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import window_pinball as wp  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--fib", required=True)
    ap.add_argument("--config", required=True, help="например 24T")
    ap.add_argument("--contact", type=int, required=True)
    ap.add_argument("--v0", type=float, required=True)
    ap.add_argument("--mu", type=float, required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--t-max", type=float, required=True)
    ap.add_argument("--fiber-E", type=float, default=8e10)
    ap.add_argument("--eps-break", type=float, default=0.05)
    ap.add_argument("--level", default="base")
    ap.add_argument("--extra", default="")
    ap.add_argument("--save-t", type=float, default=7.5e-6)
    ap.add_argument("--energy-log", type=int, default=100)
    ap.add_argument("--checkpoint-every", type=int, default=500)
    ap.add_argument("--slot-script", default=os.path.expanduser("~/proj_runs/threads_free.sh"))
    ap.add_argument("--proj-status", default=os.path.expanduser("~/proj_runs/proj_status.sh"))
    ap.add_argument("--min-uptime", type=float, default=600)
    ap.add_argument("--slot-poll", type=float, default=60)
    ap.add_argument("--slot-verify-s", type=float, default=120)
    args = ap.parse_args()
    wp.BALL_ONLY = True
    wp.STOP_VZ_FRAC = 0.0
    os.makedirs(args.work, exist_ok=True)
    wp.log(args.work, "один выстрел: %s contact=%d V0=%g mu=%g, %d потоков, t_max=%g; доп. ключи: %s" % (
        args.config, args.contact, args.v0, args.mu, args.threads, args.t_max, args.extra))
    try:
        r = wp.run_shot(args, args.mu, args.v0, args.threads, args.t_max)
    except Exception as exc:
        wp.log(args.work, "ОСТАНОВЛЕНО: %s" % exc)
        return 1
    p, why = wp.verdict(r)
    line = "%s contact=%d: %s (%s), V_R=%s, t_end=%s" % (
        args.config, args.contact, wp.word(p), why, r.get("V_R"), r.get("t_end"))
    wp.log(args.work, "итог: " + line)
    with open(os.path.join(args.work, "result.txt"), "w", encoding="utf-8") as f:
        f.write(line + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
