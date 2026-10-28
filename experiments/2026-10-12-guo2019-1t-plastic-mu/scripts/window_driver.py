#!/usr/bin/env python3
"""Подбор окна скоростей Guo 2019 (1T base) с пластично-вязкой нитью: последовательно по mu.

Использует функции plastic_driver.py (запуск guo2019_shot, разбор shot.txt/shots.csv/monitor.csv, строгий критерий,
блоки eval.md / shots.csv / eval.jsonl, дозапись при перезапуске). Критерий ВСЕГДА строгий:
  «не пробит» (НП) = исход not_perforated (в т.ч. по правилу остановки шара) И V_R = 0 И порвано верт. <= 5 И гориз. = 0;
  иначе «пробит» (П).

Алгоритм для текущего mu (V_np -- наибольшая НП, V_p -- наименьшая П):
  1. ОДНОВРЕМЕННО два выстрела 73 и 79 м/с (два процесса guo2019_shot, по --threads=4 потока, OMP_NUM_THREADS=4).
  2. 73 НП и 79 П            -> окно [73; 79].
     73 П и 79 П             -> выстрел 67 (одиночный, --threads-single потоков, по умолчанию 8):
                                67 НП -> окно [67; 73]; 67 П -> для mu окна нет.
     73 НП и 79 НП           -> выстрел 85: 85 П -> окно [79; 85]; 85 НП -> окна нет.
     73 П и 79 НП            -> «немонотонно», окна для mu нет (без дальнейших выстрелов), отмечается в отчёте.
     любой исход «не определён» (после повтора с t_max_long) -> для mu вывод невозможен, окна нет.
  3. Окно допустимо, если V_p - V_np <= 6 и 66 <= V_np, V_p <= 87.
  4. Окно найдено -> драйвер ОСТАНАВЛИВАЕТСЯ (следующие mu не считаются). Иначе -- следующий mu.

Падение выстрела (код возврата guo2019_shot != 0 или нет shots.csv) -- это НЕ «не определён» и НЕ «окна нет»:
  - упавший выстрел перезапускается один раз автоматически (лог первой попытки сохраняется как *.fail1.log);
  - если повтор тоже упал -- серия ОСТАНАВЛИВАЕТСЯ с сообщением в queue.log и window.md (exit 1), следующие mu не считаются.
    Перезапуск той же командой: готовые выстрелы не пересчитываются, упавший запускается заново.
  - при пары 73/79 падение одного выстрела не отбрасывает второй: готовый результат записывается, повторяется только упавший.
Перед стартом проверяются --bin и --fib (нет файла -- сообщение и выход 1).

Каталог выстрела: <work>/mu<mu>_v<V0>/. Итог: <work>/window.md (переписывается после каждого шага).
Перезапуск безопасен: готовые выстрелы не пересчитываются (логика окна восстанавливается по готовым результатам).
"""
import argparse, concurrent.futures, copy, os, sys, threading

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plastic_driver as pd  # noqa: E402

MUS = [0.1, 0.3, 0.5, 0.7, 0.9]
V_LO, V_HI = 73.0, 79.0     # первая пара
V_BELOW, V_ABOVE = 67.0, 85.0
MAX_WIDTH, V_MIN, V_MAX = 6.0, 66.0, 87.0
ATTEMPTS = 2                # выстрел запускается не более двух раз (исходный + один повтор при падении)

_LOCK = threading.Lock()


def _locked(fn):
    def wrapper(*a, **k):
        with _LOCK:
            return fn(*a, **k)
    return wrapper


# eval.md / shots.csv / eval.jsonl пишутся из двух потоков -- сериализуем.
pd.write_outputs = _locked(pd.write_outputs)


class ShotFailed(Exception):
    """Выстрел упал и повтор тоже упал: серия останавливается (окно не определено)."""


def window_ok(v_np, v_p):
    return (v_p - v_np) <= MAX_WIDTH + 1e-9 and v_np >= V_MIN - 1e-9 and v_p <= V_MAX + 1e-9


def run_checked(args, mu, v0, tmax, suffix, run):
    """Запуск выстрела с одним повтором при падении. Падение = код возврата != 0 или нет shots.csv.
    Возвращает результат (dict); если повтор тоже упал -- бросает ShotFailed."""
    sdir = pd.shot_dir(args.work, mu, v0, suffix)
    why = ""
    for attempt in range(1, ATTEMPTS + 1):
        if attempt > 1:
            if os.path.exists(sdir + ".log"):
                os.replace(sdir + ".log", sdir + ".fail1.log")
            pd.log("ПАДЕНИЕ mu=%.3f V0=%.1f (%s): %s -- повтор выстрела" % (mu, v0, run, why))
        pd.log("ЗАПУСК mu=%.3f V0=%.1f (%s) t_max=%g, попытка %d" % (mu, v0, run, tmax, attempt))
        res = pd.run_shot(args, mu, v0, tmax, suffix)
        rc = pd.RC.get(sdir)
        if res is not None and rc == 0:
            return res
        why = "код возврата %s%s" % (rc, "" if res is not None else ", shots.csv не записан")
    raise ShotFailed("mu=%g V0=%g (%s): выстрел упал дважды (%s); лог: %s.log" % (mu, v0, run, why, sdir))


def get_shot(args, mu, v0):
    """Как pd.get, но падение запуска -> повтор (run_checked); повторное падение -> ShotFailed.
    Исход «не определён» (без падения) -> повтор с t_max_long, как раньше."""
    sdir = pd.shot_dir(args.work, mu, v0)
    res = pd.read_result(sdir, v0, args.t_max)
    if res is None:
        res = run_checked(args, mu, v0, args.t_max, "", "base")
        pd.evaluate(args.work, mu, v0, "base", res)
    else:
        pd.ensure_eval(args.work, mu, v0, "base", res)
    if pd.effective(res) == "undetermined":
        rl = pd.read_result(pd.shot_dir(args.work, mu, v0, "_long"), v0, args.t_max_long)
        if rl is None:
            pd.log("исход не определён -> повтор mu=%.3f V0=%.1f с t_max=%g" % (mu, v0, args.t_max_long))
            rl = run_checked(args, mu, v0, args.t_max_long, "_long", "long")
            pd.evaluate(args.work, mu, v0, "long", rl)
        else:
            pd.ensure_eval(args.work, mu, v0, "long", rl)
        if pd.effective(rl) != "undetermined":
            return rl
    return res


def shoot(args, mu, v0, threads):
    """Один выстрел (или чтение готового) -> запись {'V0','verdict': True=П/False=НП/None, 'failed', 'why', ...}."""
    a = copy.copy(args)
    a.threads = threads
    rec = {"V0": v0, "verdict": None, "failed": None, "why": "", "outcome": None, "V_R": None, "nbv": None, "nbh": None,
           "eps_p_max": None, "run": "base"}
    try:
        res = get_shot(a, mu, v0)
        base = pd.read_result(pd.shot_dir(args.work, mu, v0), v0, args.t_max)
        if base is not None and pd.effective(base) == "undetermined":
            rec["run"] = "long"
    except ShotFailed as e:
        rec["failed"] = str(e)
        rec["why"] = "ВЫСТРЕЛ УПАЛ: " + str(e)
        return rec
    except Exception:
        pd.report_error("mu=%.3f V0=%.1f" % (mu, v0))
        rec["failed"] = "ошибка драйвера mu=%g V0=%g (см. queue.log)" % (mu, v0)
        rec["why"] = rec["failed"]
        return rec
    perf, why = pd.strict_verdict(res)
    rec.update({"verdict": perf, "why": why, "outcome": res.get("outcome"), "V_R": res.get("V_R"),
                "nbv": res.get("nbv"), "nbh": res.get("nbh"), "eps_p_max": res.get("eps_p_max")})
    return rec


def shoot_pair(args, mu, v1, v2):
    """Два выстрела одновременно (по args.threads потоков каждый). Падение одного не отменяет второй."""
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        f1 = ex.submit(shoot, args, mu, v1, args.threads)
        f2 = ex.submit(shoot, args, mu, v2, args.threads)
        return f1.result(), f2.result()


def shoot_single(args, mu, v0):
    os.environ["OMP_NUM_THREADS"] = str(args.threads_single)
    return shoot(args, mu, v0, args.threads_single)


def word(r):
    if r.get("failed"):
        return "УПАЛ"
    return {True: "ПРОБИТ", False: "НЕ ПРОБИТ", None: "НЕ ОПРЕДЕЛЁН"}[r["verdict"]]


def check_failed(rec_mu, recs):
    bad = [r for r in recs if r.get("failed")]
    if bad:
        raise ShotFailed(bad[0]["failed"])


def process_mu(args, mu, rec_mu, on_step):
    """Логика окна для одного mu. rec_mu = {'mu','shots','status','window','note'}; on_step() -- запись window.md.
    Падение выстрела -> ShotFailed (записывается в shots и window.md до выхода)."""
    r73, r79 = shoot_pair(args, mu, V_LO, V_HI)
    rec_mu["shots"] += [r73, r79]
    on_step()
    check_failed(rec_mu, (r73, r79))
    p73, p79 = r73["verdict"], r79["verdict"]
    if p73 is None or p79 is None:
        rec_mu["status"] = "no_window"
        rec_mu["note"] = "исход не определён для %s: вывод невозможен" % ", ".join(
            "%g м/с" % r["V0"] for r in (r73, r79) if r["verdict"] is None)
        return
    if p73 is False and p79 is True:
        v_np, v_p = V_LO, V_HI
    elif p73 is True and p79 is False:
        rec_mu["status"] = "nonmonotonic"
        rec_mu["note"] = ("НЕМОНОТОННО: %g м/с пробит, а %g м/с не пробит (для mu окна нет, дальнейшие выстрелы не делались)"
                          % (V_LO, V_HI))
        return
    elif p73 is True and p79 is True:
        r = shoot_single(args, mu, V_BELOW)
        rec_mu["shots"].append(r)
        on_step()
        check_failed(rec_mu, (r,))
        if r["verdict"] is None:
            rec_mu["status"], rec_mu["note"] = "no_window", "исход %g м/с не определён" % V_BELOW
            return
        if r["verdict"] is True:
            rec_mu["status"], rec_mu["note"] = "no_window", "%g м/с уже пробит: окно ниже диапазона" % V_BELOW
            return
        v_np, v_p = V_BELOW, V_LO
    else:  # оба НП
        r = shoot_single(args, mu, V_ABOVE)
        rec_mu["shots"].append(r)
        on_step()
        check_failed(rec_mu, (r,))
        if r["verdict"] is None:
            rec_mu["status"], rec_mu["note"] = "no_window", "исход %g м/с не определён" % V_ABOVE
            return
        if r["verdict"] is False:
            rec_mu["status"], rec_mu["note"] = "no_window", "%g м/с всё ещё не пробит: окно выше диапазона" % V_ABOVE
            return
        v_np, v_p = V_HI, V_ABOVE
    if window_ok(v_np, v_p):
        rec_mu["status"], rec_mu["window"] = "window", (v_np, v_p)
        rec_mu["note"] = "окно [%g; %g] м/с (ширина %g)" % (v_np, v_p, v_p - v_np)
    else:
        rec_mu["status"], rec_mu["note"] = "no_window", "окно [%g; %g] недопустимо (ширина/границы)" % (v_np, v_p)


STATUS_RU = {"window": "ОКНО НАЙДЕНО", "no_window": "окна нет", "nonmonotonic": "немонотонно, окна нет",
             "running": "в процессе", "stopped": "СЕРИЯ ОСТАНОВЛЕНА (падение выстрела)"}


def write_window(args, records, finished, found, stop=None):
    L = ["# Подбор окна скоростей: 1T base, пластично-вязкая нить", "",
         "- fib: `%s`; eps_break = %g; критерий строгий (не пробит = not_perforated, V_R = 0, верт. <= %d, гориз. = %d)"
         % (os.path.basename(args.fib), args.eps_break, pd.STRICT_MAX_V, pd.STRICT_MAX_H),
         "- допустимое окно: V_p - V_np <= %g м/с, %g <= V_np, V_p <= %g м/с" % (MAX_WIDTH, V_MIN, V_MAX), ""]
    for rec in records:
        L.append("## mu = %g: %s" % (rec["mu"], STATUS_RU.get(rec["status"], rec["status"])))
        L.append("")
        L.append("| V0, м/с | вердикт | исход | V_R, м/с | верт. | гориз. | eps_p max | причина |")
        L.append("|---|---|---|---|---|---|---|---|")
        for r in sorted(rec["shots"], key=lambda r: r["V0"]):
            L.append("| %g | %s | %s%s | %s | %s | %s | %s | %s |" % (
                r["V0"], word(r), r["outcome"] or "нет данных", " (повтор t_max_long)" if r["run"] == "long" else "",
                pd.f2(r["V_R"]), pd.f2(r["nbv"], "%.0f"), pd.f2(r["nbh"], "%.0f"), pd.f2(r.get("eps_p_max"), "%.5f"),
                r["why"]))
        if rec.get("note"):
            L += ["", "Итог mu: " + rec["note"]]
        L.append("")
    L.append("## Итог серии")
    L.append("")
    if stop is not None:
        L.append("СЕРИЯ ОСТАНОВЛЕНА: %s. Окно не определено, следующие mu не считались. Перезапуск той же командой: "
                 "готовые выстрелы не пересчитываются, упавший выстрел запускается заново." % stop)
    elif found is not None:
        L.append("ОКНО НАЙДЕНО: mu = %g, [%g; %g] м/с (V_np = %g, V_p = %g). Серия остановлена." % (
            found["mu"], found["window"][0], found["window"][1], found["window"][0], found["window"][1]))
    elif finished:
        L.append("Окно не найдено ни для одного mu из %s. Дальше: менять параметры нити (см. README, порядок варьирования)." %
                 ", ".join("%g" % r["mu"] for r in records))
    else:
        L.append("Серия идёт (или прервана): см. выше состояние по mu.")
    tmp = os.path.join(args.work, "window.md.tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    os.replace(tmp, os.path.join(args.work, "window.md"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.path.expanduser("~/window1T"))
    ap.add_argument("--bin", default=os.path.expanduser("~/SPH_plastic/build-elog/guo2019_shot"))
    ap.add_argument("--fib", default=os.path.join(os.path.expanduser("~/window1T"), "twaron_plastic_visco_rate.fib"))
    ap.add_argument("--eps-break", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=4, help="потоков на каждый из двух параллельных выстрелов")
    ap.add_argument("--threads-single", type=int, default=8, help="потоков на одиночный выстрел (67/85 м/с)")
    ap.add_argument("--energy-log", type=int, default=100, metavar="N")
    ap.add_argument("--t-max", type=float, default=1.5e-3)
    ap.add_argument("--t-max-long", type=float, default=3.0e-3)
    ap.add_argument("--mus", default=",".join(str(m) for m in MUS))
    args = ap.parse_args()
    for label, path in (("--bin", args.bin), ("--fib", args.fib)):
        if not os.path.isfile(path):
            print("ОШИБКА: файл %s не найден: %s -- серия не запущена" % (label, path), flush=True)
            sys.exit(1)
    pd.CRITERION = "strict"
    os.makedirs(args.work, exist_ok=True)
    mus = [float(x) for x in args.mus.split(",")]
    pd.log("серия окна 1T base: mu=%s fib=%s eps_break=%g threads=%d/%d" % (
        mus, args.fib, args.eps_break, args.threads, args.threads_single))
    if args.energy_log > 0:
        pd.check_energy_build(args.bin)
    pd.recover(args, mus, [V_BELOW, V_LO, V_HI, V_ABOVE])
    records, found, stop = [], None, None
    for mu in mus:
        rec = {"mu": mu, "shots": [], "status": "running", "window": None, "note": ""}
        records.append(rec)

        def on_step():
            write_window(args, records, False, None)

        try:
            process_mu(args, mu, rec, on_step)
        except ShotFailed as e:
            stop = str(e)
            rec["status"], rec["note"] = "stopped", "серия остановлена: " + str(e)
            pd.log("ОСТАНОВ СЕРИИ (падение выстрела), следующие mu не считаются: %s" % stop)
            break
        except Exception:
            pd.report_error("mu=%.3f" % mu)
            stop = "mu=%g: ошибка драйвера (см. queue.log)" % mu
            rec["status"], rec["note"] = "stopped", "серия остановлена: ошибка драйвера (см. queue.log)"
            pd.log("ОСТАНОВ СЕРИИ (ошибка драйвера), следующие mu не считаются: %s" % stop)
            break
        pd.log("mu=%g: %s -- %s" % (mu, STATUS_RU.get(rec["status"], rec["status"]), rec["note"]))
        if rec["status"] == "window":
            found = rec
            break
        write_window(args, records, False, None)
    write_window(args, records, stop is None and found is None, found, stop)
    if stop is not None:
        pd.log("СЕРИЯ ОСТАНОВЛЕНА: %s" % stop)
        sys.exit(1)
    pd.log("СЕРИЯ ЗАВЕРШЕНА: %s" % ("окно найдено, mu=%g" % found["mu"] if found else "окна нет"))


if __name__ == "__main__":
    main()
