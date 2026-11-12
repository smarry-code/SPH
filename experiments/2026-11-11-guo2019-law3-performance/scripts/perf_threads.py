#!/usr/bin/env python3
"""Разбивка многопоточного профиля perf (record --call-graph fp) по потокам и параллельным регионам OpenMP.

Каждый сэмпл относится к одному из классов (по верхнему кадру вне ядра ОС и по наличию в стеке кадра libgomp, вызвавшего код проекта):
  par_work  -- работа внутри параллельного региона (код проекта, вызванный из libgomp);
  par_wait  -- ожидание libgomp внутри региона (барьер в конце цикла: дисбаланс нагрузки);
  idle_wait -- ожидание libgomp вне региона (рабочие потоки ждут, пока главный поток считает последовательно);
  serial    -- работа вне параллельного региона (последовательный код; должен быть только у главного потока).
Окна -- доли длительности записи (--win имя:a-b, a и b в [0, 1]); окна по шагам задаёт bench_tables.py.

Запуск: python3 perf_threads.py <perf.data> [--win pre:0.03-0.30 --win imp:0.40-1.0] [--top 25]
"""
import argparse
import collections
import re
import subprocess

HDR = re.compile(r"^\s*(\d+)\s+([\d.]+):")
FRAME = re.compile(r"^\s+([0-9a-f]+)\s+(.*?)\s+\((.*)\)\s*$")
LAMBDA_OWNER = re.compile(r"([A-Za-z_][\w:]*)\((?:[^()]|\([^()]*\))*\)(?: const)?::\{lambda")


def region_key(sym):
    """Короткое имя региона: цепочка функций, в которых определены лямбды, переданные в parallel_for."""
    owners = [o.split("::")[-1] if o.startswith("rhs::CpuExecutor") else o for o in LAMBDA_OWNER.findall(sym)]
    owners = [o for o in owners if o not in ("parallel_for", "parallel_for_segments", "parallel_reduce")]
    seen = []
    for o in owners:
        if o not in seen:
            seen.append(o)
    if seen:
        return " > ".join(seen)
    m = re.search(r"parallel_\w+<([A-Za-z_][\w:]*)", sym)  # parallel_for<rhs::CpuExecutor::concat_segments<...
    if m:
        return m.group(1).split("::")[-1] + " [libc]"
    m = re.search(r"parallel_\w+", sym)
    return (m.group(0) if m else sym[:60]) + " [?]"


def is_kernel(dso):
    return dso.startswith("[kernel") or dso == "[unknown]"


def resolve_inlined(fr):
    """perf >= 6 разворачивает inline-кадры с dso "inlined": им приписывается dso ближайшего настоящего кадра ниже."""
    out, real = [], "?"
    for sym, dso in reversed(fr):
        if dso == "inlined":
            dso = real
        else:
            real = dso
        out.append((sym, dso))
    out.reverse()
    return out


def samples(path):
    """Поток (tid, time, frames[(sym, dso)]) из perf script."""
    p = subprocess.Popen(["perf", "script", "-i", path, "-F", "tid,time,ip,sym,dso"],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, errors="replace", bufsize=1 << 20)
    tid, t, fr = None, None, []
    for line in p.stdout:
        if not line.strip():
            if tid is not None:
                yield tid, t, resolve_inlined(fr)
            tid, t, fr = None, None, []
            continue
        m = HDR.match(line)
        if m and tid is None:
            tid, t = int(m.group(1)), float(m.group(2))
            continue
        m = FRAME.match(line)
        if m and tid is not None:
            fr.append((m.group(2), m.group(3).rsplit("/", 1)[-1]))
    if tid is not None:
        yield tid, t, resolve_inlined(fr)
    p.wait()


MEMFN = re.compile(r"^(__GI_|__)?(mem(cpy|move|set|cmp)|str\w*|wmem\w*)")  # суффиксы ifunc (_evex_/_avx2_/_sse2_...) -- \w*


def is_libc(dso):
    return dso.startswith("libc.") or dso.startswith("libc-")


def classify(fr, binary="guo2019_shot"):
    """Регион -- код проекта, вызванный из libgomp (GOMP_parallel у главного потока, gomp_thread_start у рабочих).
    Имя вынесенной функции берётся из кадра проекта сразу над libgomp: при inline-развёртке perf печатает имя
    из DWARF без суффикса ._omp_fn, поэтому по суффиксу регион не определяется."""
    fr = [(s, d) for s, d in fr if not is_kernel(d)]
    top = fr[0] if fr else ("?", "?")
    omp = None
    for i in range(1, len(fr)):
        if "libgomp" in fr[i][1] and fr[i - 1][1] == binary:
            omp = fr[i - 1][0]
            break
    wait = "libgomp" in top[1]
    if omp is not None:
        return ("par_wait" if wait else "par_work"), omp, top
    if not wait and is_libc(top[1]):
        # Листовая функция libc без frame pointer: при --call-graph fp кадр её вызывающего (вынесенной функции
        # региона) теряется, и над libgomp оказывается кадр libc, а не проекта. Это работа региона ТОЛЬКО для
        # memory/string-процедур (memmove/memset/strlen...); любой другой листовой вызов libc над libgomp
        # (futex, sched_yield, clock_nanosleep, syscall) -- ожидание внутри libgomp, а не работа.
        # Имя региона -- первый кадр проекта ниже libgomp (у главного потока это parallel_for<...>,
        # вызвавший GOMP_parallel), иначе "? (<функция libc>)".
        for i in range(1, len(fr)):
            if "libgomp" in fr[i][1]:
                if MEMFN.match(top[0]):
                    below = next((s for s, d in fr[i + 1:] if d == binary), None)
                    return "par_work", below or f"? ({top[0]})", top
                return "idle_wait", None, top
    return ("idle_wait" if wait else "serial"), None, top


GENERIC = re.compile(r"^(operator|pair|construct|set|get|lookup|push_back|emplace|size|begin|end|_M_|__|std::|"
                     r"allocate|max|min|data|at|swap|move|forward|invoke|~)")


def prog_frames(fr, binary):
    """Кадры проекта без служебных (operator() лямбд, методы контейнеров): первый -- содержательная функция."""
    pf = [s for s, d in fr if d == binary]
    named = [s for s in pf if not GENERIC.match(s)]
    return named or pf


def short(sym, n=110):
    sym = re.sub(r"\(.*", "", sym) if "{lambda" not in sym else sym
    return sym if len(sym) <= n else sym[:n] + "…"


def _selftest():
    """Синтетические стеки (верхний кадр первым): классификация без perf."""
    B, G, C = "guo2019_shot", "libgomp.so.1", "libc.so.6"
    cases = [
        ("главный поток: memmove в регионе", [("__memmove_evex_unaligned_erms", C), ("rhs::parallel_for<X>", B),
                                              ("GOMP_parallel", G), ("rhs::parallel_for<X>", B), ("main", B)], "par_work"),
        ("рабочий поток: memmove в регионе", [("__memset_avx2_unaligned_erms", C), ("gomp_thread_start", G),
                                              ("start_thread", C)], "par_work"),
        ("рабочий поток: syscall над libgomp (ожидание)", [("syscall", C), ("gomp_barrier_wait_end", G),
                                                           ("gomp_thread_start", G), ("start_thread", C)], "idle_wait"),
        ("последовательный memmove (нет libgomp)", [("__memmove_avx_unaligned", C), ("step()", B), ("main", B)], "serial"),
        ("код проекта над libgomp", [("foo()", B), ("foo() [clone ._omp_fn.0]", B), ("GOMP_parallel", G), ("main", B)],
         "par_work"),
        ("ожидание libgomp внутри региона", [("gomp_team_barrier_wait_end", G), ("foo._omp_fn.0", B),
                                             ("GOMP_parallel", G)], "par_wait"),
    ]
    bad = 0
    for name, fr, want in cases:
        got = classify(fr, B)[0]
        ok = got == want
        bad += not ok
        print(f"{'OK ' if ok else 'ОШИБКА'} {name}: {got} (ожидалось {want})")
    print("selftest:", "пройден" if not bad else f"ошибок {bad}")
    raise SystemExit(1 if bad else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("perf_data", nargs="?")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--win", action="append", default=[])
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--binary", default="guo2019_shot")
    a = ap.parse_args()
    if a.selftest:
        _selftest()
    if not a.perf_data:
        ap.error("нужен perf_data (или --selftest)")
    wins = []
    for w in a.win or ["all:0-1"]:
        name, rng = w.split(":")
        lo, hi = map(float, rng.split("-"))
        wins.append((name, lo, hi))

    # проход 1 (без стеков): границы записи и главный поток
    t0 = t1 = None
    tids = set()
    q = subprocess.Popen(["perf", "script", "-i", a.perf_data, "-F", "tid,time", "--hide-call-graph"],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, errors="replace")
    for line in q.stdout:
        m = HDR.match(line)
        if m:
            t = float(m.group(2))
            tids.add(int(m.group(1)))
            t0 = t if t0 is None else min(t0, t)
            t1 = t if t1 is None else max(t1, t)
    q.wait()
    if t0 is None:
        print("нет сэмплов")
        return
    main_tid = min(tids)
    print(f"потоков {len(tids)}, запись {t1 - t0:.1f} с, главный tid {main_tid}")

    W = []
    for name, lo, hi in wins:
        W.append(dict(name=name, lo=lo, hi=hi, a0=t0 + lo * (t1 - t0), a1=t0 + hi * (t1 - t0), n=0,
                      cls=collections.Counter(), per_tid=collections.defaultdict(collections.Counter),
                      reg=collections.defaultdict(collections.Counter), ser=collections.Counter(),
                      ser_ctx=collections.Counter(), par_top=collections.Counter()))
    # проход 2: стеки, счётчики по окнам на лету
    for tid, t, fr in samples(a.perf_data):
        hit = [w for w in W if w["a0"] <= t <= w["a1"]]
        if not hit:
            continue
        c, omp, top = classify(fr, a.binary)
        rk = region_key(omp) if omp is not None else None
        pt = None
        if c == "par_work":
            pf = prog_frames(fr, a.binary)
            pt = short(pf[0]) if pf else f"{short(top[0])} ({top[1]})"
        sk = ck = None
        if c == "serial":
            pf = prog_frames(fr, a.binary)
            sk = short(pf[0]) if pf else f"{short(top[0])} ({top[1]})"
            ck = " < ".join(short(s, 50) for s in pf[:3]) if pf else "?"
        for w in hit:
            w["n"] += 1
            w["cls"][c] += 1
            w["per_tid"][tid][c] += 1
            if rk is not None:
                w["reg"][rk][c] += 1
            if pt is not None:
                w["par_top"][pt] += 1
            if sk is not None:
                w["ser"][sk] += 1
                w["ser_ctx"][ck] += 1

    for w in W:
        name, lo, hi, n = w["name"], w["lo"], w["hi"], w["n"]
        cls, per_tid, reg, ser, ser_ctx, par_top = (w[k] for k in ("cls", "per_tid", "reg", "ser", "ser_ctx", "par_top"))
        if n == 0:
            continue
        nt = len(per_tid)
        print(f"\n=== окно {name} [{lo:.2f}, {hi:.2f}]: сэмплов {n}, потоков {nt}")
        print("классы (доля всех сэмплов окна): " + "  ".join(f"{k}={100 * v / n:.1f}%" for k, v in sorted(cls.items())))
        mt = per_tid[main_tid]
        mtot = sum(mt.values())
        if mtot:
            print(f"главный поток: serial={100 * mt['serial'] / mtot:.1f}% par_work={100 * mt['par_work'] / mtot:.1f}% "
                  f"par_wait={100 * mt['par_wait'] / mtot:.1f}% idle_wait={100 * mt['idle_wait'] / mtot:.1f}%")
        eff = cls["par_work"] + cls["serial"]
        print(f"полезная работа (par_work+serial) {100 * eff / n:.1f}% от {nt} x время; "
              f"доля последовательных сэмплов главного потока (не f_s по Амдалу) = {mt['serial'] / max(mtot, 1):.3f}")
        print("по потокам: tid  par_work par_wait idle_wait serial")
        for tid in sorted(per_tid):
            c = per_tid[tid]
            s = sum(c.values())
            print(f"  {tid:7d} {100 * c['par_work'] / s:7.1f}% {100 * c['par_wait'] / s:7.1f}% "
                  f"{100 * c['idle_wait'] / s:8.1f}% {100 * c['serial'] / s:6.1f}%")
        print("параллельные регионы (доля всех сэмплов окна; wait/(work+wait) -- дисбаланс):")
        for k, c in sorted(reg.items(), key=lambda kv: -sum(kv[1].values()))[:a.top]:
            wk, ww = c["par_work"], c["par_wait"]
            print(f"  {100 * (wk + ww) / n:5.1f}%  work {100 * wk / n:5.1f}%  wait {100 * ww / n:5.1f}%  "
                  f"imb {ww / max(wk + ww, 1):.2f}  {k}")
        print("последовательный код главного потока: верхняя функция проекта (доля всех сэмплов окна, доля serial):")
        st = max(cls["serial"], 1)
        for k, v in ser.most_common(a.top):
            print(f"  {100 * v / n:5.2f}% {100 * v / st:5.1f}%  {k}")
        print("последовательный код: стек проекта (3 кадра):")
        for k, v in ser_ctx.most_common(a.top):
            print(f"  {100 * v / n:5.2f}% {100 * v / st:5.1f}%  {k}")
        print("работа в регионах: верхний кадр:")
        for k, v in par_top.most_common(a.top):
            print(f"  {100 * v / n:5.2f}%  {k}")


if __name__ == "__main__":
    main()
