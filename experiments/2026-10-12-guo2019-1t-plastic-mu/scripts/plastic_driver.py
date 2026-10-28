#!/usr/bin/env python3
"""Серия выстрелов Guo 2019, 1T, уровень base, ПЛАСТИЧНАЯ нить, варьирование трения mu.

Образец -- ~/mu_driver.py (Guo2T, подбор mu): тот же вызов guo2019_shot (--n-mon=10 --n-diag=100 --save-t=20e-6 --frames=off),
тот же разбор shot.txt/shots.csv/monitor.csv, то же правило «undetermined -> один повтор с большим t_max».
Отличия: сетка фиксирована (не подбор); добавлены --fiber-inelastic/--eps-break; в оценку добавлены max eps_p, W_p и n_bonds_plastic.

Порядок: сначала все V0 = 73 по возрастанию mu, затем все V0 = 79. Последовательно, один выстрел на 8 потоков.
После КАЖДОГО выстрела: блок оценки дописывается в <work>/eval.md, строка -- в <work>/shots.csv (и eval.jsonl).
Каталог выстрела: <work>/mu<mu>_v<V0>/ (там 1T_base_v<V0>_r0/shot.txt, monitor.csv, кадры; рядом .log).
Перезапуск безопасен: готовые выстрелы (есть shots.csv в каталоге) не пересчитываются. При старте для них проверяется,
что блоки записаны в eval.md, shots.csv и eval.jsonl (ключ: mu, V0, run = base|long); недостающие части дописываются,
уже записанные не дублируются. Ошибка в одном выстреле печатает traceback (stdout = queue.log) и серия идёт дальше.
Отсутствие бинаря guo2019_shot -- сообщение и выход.
Опции: --energy-log=N (пробросить в guo2019_shot, нужна сборка с -DSPH_ENERGY_LOG=ON; energy_log.csv в подкаталоге выстрела);
--criterion=strict: «не пробит» = исход not_perforated И V_R = 0 И порвано вертикальных <= 5 И горизонтальных = 0, иначе «пробит»
(в т.ч. остановка шара при большом разрушении; причина пишется в eval.md). В eval.md исход пишется по критерию, без «цели».
--config=1T|2T|4T (по умолчанию 1T): мишень, передаётся в guo2019_shot. Для 1T всё прежнее (имена каталогов mu<mu>_v<V0>,
столбцы shots.csv, критерий). Для 2T/4T: каталоги <config>_mu<mu>_v<V0>; в eval.md/shots.csv/eval.jsonl разрывы печатаются
ПОСЛОЙНО (вертикальные = вдоль y, зажатые; горизонтальные = вдоль x, свободные; время первого разрыва, прогиб слоя при нём,
втягивание нитей слоя); критерий по умолчанию strict-слойный: «не пробит» = шар остановлен (V_R = 0) И в КАЖДОМ слое
<= 7 вертикальных и 0 горизонтальных порванных нитей (--criterion=default для 2T/4T -- только исход программы).
"""
import argparse, csv, glob, json, os, subprocess, sys, time, traceback

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = {"broken_v": 3, "broken_h": 0, "wf_lo": 25.0, "wf_hi": 27.0}
MUS = [0.3, 0.5, 0.7, 0.9, 1.0]
V0S = [73.0, 79.0]
CSV_COLS = ["mu", "V0", "run", "outcome", "effective", "V_R", "e", "n_broken_v", "n_broken_h", "t_first_break_us",
            "wf_first_break_mm", "eps_p_max", "W_p", "E_tot_max/E_tot_0", "E_tot_end/E_tot_0", "n_tunnel", "t_end_us",
            "t_max_us", "wall_s", "strict_perforated"]
MD_PREFIX = "### mu = %.3f, V0 = %.1f м/с"
STRICT_MAX_V = 5  # strict: не более стольких порванных вертикальных нитей
STRICT_MAX_H = 0
STRICT_MAX_V_PLY = 7  # strict для многослойного пакета: в каждом слое не более стольких порванных вертикальных нитей
CRITERION = "default"  # задаётся --criterion
CONFIG = "1T"  # задаётся --config; для 1T все имена/столбцы прежние
MU_FF = None  # --mu-ff: трение нить-нить отдельно от --mu (None -- не задано, прежнее поведение)
MU_FS = None  # --mu-fs: трение нить-Solid отдельно от --mu
PLIES = {"1T": 1, "2T": 2, "4T": 4}
RC = {}  # каталог выстрела -> код возврата guo2019_shot (заполняет run_shot; для драйвера окна)


def multi():
    """Многослойный пакет (2T/4T): включает послойный вывод и послойный критерий."""
    return PLIES.get(CONFIG, 1) >= 2


def ply_csv_cols(n):
    c = ["n_plies"]
    for k in range(1, n + 1):
        c += ["n_broken_v_p%d" % k, "n_broken_h_p%d" % k, "t_first_break_us_p%d" % k, "wf_first_break_mm_p%d" % k]
    return c


def csv_cols():
    return CSV_COLS + (ply_csv_cols(PLIES.get(CONFIG, 1)) if multi() else [])


def md_prefix(mu, v0):
    return (MD_PREFIX if CONFIG == "1T" else "### " + CONFIG + ": mu = %.3f, V0 = %.1f м/с") % (mu, v0)


def log(m):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), m, flush=True)


def report_error(where):
    """Traceback в stdout (при запуске по инструкции stdout = queue.log); серия продолжается."""
    log("ОШИБКА (%s), выстрел пропущен, серия продолжается:\n%s" % (where, traceback.format_exc()))


def kv(path):
    d = {}
    with open(path, encoding="utf-8") as f:
        for t in f.read().split():
            if "=" in t:
                k, v = t.split("=", 1)
                d[k] = v
    return d


def num(x):
    """float или None (отсутствует, не число, nan)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def scale(x, k):
    return None if x is None else x * k


def f2(x, fmt="%.2f"):
    return "нет данных" if x is None else fmt % x


def shot_dir(work, mu, v0, suffix=""):
    pre = "" if CONFIG == "1T" else CONFIG + "_"  # 1T: прежние имена каталогов
    sep = ("_ff%g" % MU_FF if MU_FF is not None else "") + ("_fs%g" % MU_FS if MU_FS is not None else "")
    return os.path.join(work, "%smu%g%s_v%g%s" % (pre, mu, sep, v0, suffix))  # без --mu-ff/--mu-fs: прежние имена


def read_result(sdir, v0=None, tmax=None):
    """Результат выстрела из <sdir>/shots.csv и <sdir>/*/shot.txt. Отсутствующие поля -- None (не KeyError)."""
    p = os.path.join(sdir, "shots.csv")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    r = rows[0]
    res = {"outcome": r.get("outcome") or None, "V_R": num(r.get("V_R")), "e": num(r.get("e")),
           "nbv": None, "nbh": None, "wf_first": None, "t_first": None, "etot_max": None, "etot_end": None,
           "tunnel": None, "t_end": None, "wall_s": None, "t_max": None, "eps_p_max": None, "W_p": None,
           "n_bonds_plastic": None, "plies": None}
    for st in sorted(glob.glob(os.path.join(sdir, "*", "shot.txt"))):
        d = kv(st)
        res.update({
            "nbv": num(d.get("n_broken_v")), "nbh": num(d.get("n_broken_h")),
            "wf_first": scale(num(d.get("w_f_at_first_break")), 1e3), "t_first": scale(num(d.get("t_first_break")), 1e6),
            "etot_max": num(d.get("E_tot_max/E_tot_0")), "etot_end": num(d.get("E_tot_end/E_tot_0")),
            "tunnel": num(d.get("n_tunnel")), "t_end": scale(num(d.get("t_end")), 1e6), "wall_s": num(d.get("wall_s")),
            "t_max": scale(num(d.get("t_max")), 1e6), "eps_p_max": num(d.get("eps_p_max")), "W_p": num(d.get("W_p")),
            "n_bonds_plastic": num(d.get("n_bonds_plastic")),
        })
        npl = int(num(d.get("n_plies")) or 0)  # послойные ключи пишутся только для пакета из >= 2 слоёв
        if npl >= 2:
            res["plies"] = [{
                "nbv": num(d.get("n_broken_v_p%d" % k)), "nbh": num(d.get("n_broken_h_p%d" % k)),
                "t_first": scale(num(d.get("t_first_break_p%d" % k)), 1e6),
                "wf_first": scale(num(d.get("w_f_at_first_break_p%d" % k)), 1e3),
                "pullin_first": scale(num(d.get("pullin_at_first_break_p%d" % k)), 1e3),
                "pullin_end": scale(num(d.get("pullin_end_p%d" % k)), 1e3)} for k in range(1, npl + 1)]
    if res["t_max"] is None and tmax is not None:
        res["t_max"] = tmax * 1e6
    if res["outcome"] in (None, "undetermined") and v0 is not None:  # шар остановлен к t_max -> «не пробит (правило)»
        for mon in glob.glob(os.path.join(sdir, "*", "monitor.csv")):
            with open(mon, encoding="utf-8") as f:
                head = f.readline().strip().split(",")
                last = None
                for line in f:
                    last = line
            if last and last.strip():
                m = dict(zip(head, last.strip().split(",")))
                vz, zc = num(m.get("v_z")), num(m.get("z_c"))
                if vz is not None:
                    res["vz_end"] = vz
                if zc is not None:
                    res["zc_end"] = zc * 1e3
                if vz is not None and zc is not None and abs(vz) < 0.05 * v0 and zc < 0:
                    res["rule_not_perf"] = True
    return res


def effective(res):
    if res is None:
        return None
    if res["outcome"] in ("perforated", "not_perforated"):
        return res["outcome"]
    return "not_perforated" if res.get("rule_not_perf") else "undetermined"


def strict_verdict(res):
    """Строгий критерий: (пробит: bool|None, причина). None -- данных недостаточно для вывода."""
    eff = effective(res)
    if eff == "undetermined":
        return None, "исход не определён"
    vr, nbv, nbh = res.get("V_R"), res.get("nbv"), res.get("nbh")
    if vr is None and res.get("rule_not_perf"):
        vr = 0.0  # исход not_perforated по правилу остановки шара: отскока/выбивания нет, V_R = 0 (nan в shots.csv)
    if multi():
        return strict_verdict_plies(res, eff, vr)
    reasons = []
    if eff == "perforated":
        reasons.append("исход perforated")
    if vr is not None and vr != 0.0:
        reasons.append("V_R = %s м/с != 0" % f2(vr))
    if nbv is not None and nbv > STRICT_MAX_V:
        reasons.append("порвано верт. %.0f > %d" % (nbv, STRICT_MAX_V))
    if nbh is not None and nbh > STRICT_MAX_H:
        reasons.append("порвано гориз. %.0f > %d" % (nbh, STRICT_MAX_H))
    if reasons:
        return True, "; ".join(reasons)
    if vr is None or nbv is None or nbh is None:
        return None, "нет данных (V_R/разрывы) для строгого критерия"
    return False, "not_perforated, V_R = 0, верт. %.0f <= %d, гориз. %.0f = %d" % (nbv, STRICT_MAX_V, nbh, STRICT_MAX_H)


def strict_verdict_plies(res, eff, vr):
    """Строгий критерий для 2T/4T: не пробит = шар остановлен (V_R = 0) И в каждом слое <= 7 верт. и 0 гориз. разрывов."""
    plies = res.get("plies")
    reasons = []
    if eff == "perforated":
        reasons.append("исход perforated")
    if vr is not None and vr != 0.0:
        reasons.append("V_R = %s м/с != 0" % f2(vr))
    if plies:
        for k, p in enumerate(plies, 1):
            if p["nbv"] is not None and p["nbv"] > STRICT_MAX_V_PLY:
                reasons.append("слой %d: порвано верт. %.0f > %d" % (k, p["nbv"], STRICT_MAX_V_PLY))
            if p["nbh"] is not None and p["nbh"] > STRICT_MAX_H:
                reasons.append("слой %d: порвано гориз. %.0f > %d" % (k, p["nbh"], STRICT_MAX_H))
    if reasons:
        return True, "; ".join(reasons)
    if vr is None or not plies or any(p["nbv"] is None or p["nbh"] is None for p in plies):
        return None, "нет данных (V_R/послойные разрывы) для строгого критерия"
    return False, "not_perforated, V_R = 0, в каждом из %d слоёв верт. <= %d и гориз. = %d" % (
        len(plies), STRICT_MAX_V_PLY, STRICT_MAX_H)


def ply_lines(res):
    """Строки послойной сводки (только 2T/4T): слой 1 -- лицевой (ближе к шару)."""
    plies = res.get("plies")
    if not plies:
        return ["- послойные разрывы: нет данных (в shot.txt нет n_plies / *_pK)"]
    L = ["- разрывы по слоям (слой 1 -- лицевой; верт. = вдоль y, зажатые; гориз. = вдоль x, свободные):"]
    for k, p in enumerate(plies, 1):
        L.append("  - слой %d: верт. %s, гориз. %s; первый разрыв %s мкс; прогиб слоя при нём %s мм; втягивание при нём %s мм, "
                 "в конце %s мм" % (k, f2(p["nbv"], "%.0f"), f2(p["nbh"], "%.0f"), f2(p["t_first"], "%.0f"),
                                    f2(p["wf_first"], "%.1f"), f2(p["pullin_first"], "%.2f"), f2(p["pullin_end"], "%.2f")))
    return L


def block_lines(mu, v0, run, res):
    """Строки блока оценки для eval.md."""
    eff = effective(res)
    title = md_prefix(mu, v0)
    if run == "long":
        title += " (повтор, t_max = %s мкс)" % f2(res.get("t_max"), "%.0f")
    title += ": %s" % (res.get("outcome") or "нет данных")
    if res.get("rule_not_perf"):
        title += " (правило: шар остановлен к t_max, vz_end = %s м/с, z_c = %s мм)" % (
            f2(res.get("vz_end")), f2(res.get("zc_end"), "%.1f"))
    L = [title]
    L.append("- V_R = %s м/с; e = %s Дж; t_end = %s мкс; wall = %s с" % (
        f2(res.get("V_R")), f2(res.get("e")), f2(res.get("t_end"), "%.0f"), f2(res.get("wall_s"), "%.0f")))
    L.append("- порвано нитей: верт. %s, гориз. %s; прогиб при первом разрыве %s мм; первый разрыв %s мкс" % (
        f2(res.get("nbv"), "%.0f"), f2(res.get("nbh"), "%.0f"), f2(res.get("wf_first"), "%.1f"),
        f2(res.get("t_first"), "%.0f")))
    if multi():
        L.append("  (выше -- суммы по всем слоям; по слоям -- ниже)")
        L.extend(ply_lines(res))
    em = res.get("etot_max")
    L.append("- E_tot_max/E_tot_0 = %s%s; E_tot_end/E_tot_0 = %s; n_tunnel = %s" % (
        f2(em, "%.4f"), "  **ПРЕВЫШЕНИЕ 1.005: возможная генерация энергии**" if em is not None and em > 1.005 else "",
        f2(res.get("etot_end"), "%.3f"), f2(res.get("tunnel"), "%.0f")))
    L.append("  (пояснение: E_tot_end/E_tot_0 пластичной нити меньше 1 ожидаемо -- на величину необратимой работы W_p "
             "(в E_tot не входит) и работы трения; это не потеря энергии схемой. Контроль схемы -- E_tot_max/E_tot_0.)")
    eps = res.get("eps_p_max")
    L.append("- пластика нити: max eps_p = %s (%s %%); W_p = %s Дж (пластическая работа всех связей, включая порванные); "
             "n_bonds_plastic = %s (только целые, непорванные связи с eps_p > 0 на конец расчёта; порванные сюда не входят)" % (
                 f2(eps, "%.5f"), f2(None if eps is None else 100 * eps), f2(res.get("W_p"), "%.4f"),
                 f2(res.get("n_bonds_plastic"), "%.0f")))
    if CRITERION == "strict":
        perf, why = strict_verdict(res)
        cond = ("не пробит = шар остановлен (V_R = 0), в каждом слое верт. <= %d и гориз. = %d" % (STRICT_MAX_V_PLY, STRICT_MAX_H)
                if multi() else "не пробит = not_perforated, V_R = 0, верт. <= %d, гориз. = %d" % (STRICT_MAX_V, STRICT_MAX_H))
        L.append("- исход по строгому критерию (%s): %s (%s)" % (
            cond, "исход не определён" if perf is None else "ПРОБИТ" if perf else "НЕ ПРОБИТ", why))
    else:
        L.append("- исход: %s (без сравнения с целью; V_R = %s м/с, порвано верт./гориз. = %s/%s)" % (
            {"perforated": "пробит", "not_perforated": "не пробит"}.get(eff, "не определён"), f2(res.get("V_R")),
            f2(res.get("nbv"), "%.0f"), f2(res.get("nbh"), "%.0f")))
    return L


def _key(mu, v0, run):
    return (round(float(mu), 6), round(float(v0), 3), str(run))


def done_csv(work):
    keys = set()
    p = os.path.join(work, "shots.csv")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                try:
                    keys.add(_key(r["mu"], r["V0"], r["run"]))
                except (KeyError, TypeError, ValueError):
                    pass
    return keys


def done_jsonl(work):
    keys = set()
    p = os.path.join(work, "eval.jsonl")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            for line in f:
                try:
                    o = json.loads(line)
                    keys.add(_key(o["mu"], o["V0"], o["run"]))
                except (ValueError, KeyError, TypeError):
                    pass
    return keys


def md_has(work, mu, v0, run):
    """Есть ли в eval.md заголовок блока (mu, V0) нужного рода: повтор (long) или базовый (base)."""
    p = os.path.join(work, "eval.md")
    if not os.path.exists(p):
        return False
    prefix = md_prefix(mu, v0)
    with open(p, encoding="utf-8") as f:
        for line in f:
            if line.startswith(prefix) and ("(повтор" in line) == (run == "long"):
                return True
    return False


def strict_cell(res):
    if CRITERION != "strict":
        return ""
    perf, _ = strict_verdict(res)
    return "" if perf is None else int(perf)


def write_outputs(work, mu, v0, run, res, md=True, csv_=True, jsonl=True):
    eff = effective(res)
    L = block_lines(mu, v0, run, res)
    if md:
        open(os.path.join(work, "eval.md"), "a", encoding="utf-8").write("\n".join(L) + "\n\n")
        print("\n".join(L), flush=True)
    row = {"mu": mu, "V0": v0, "run": run, "outcome": res.get("outcome"), "effective": eff, "V_R": res.get("V_R"),
           "e": res.get("e"), "n_broken_v": res.get("nbv"), "n_broken_h": res.get("nbh"),
           "t_first_break_us": res.get("t_first"), "wf_first_break_mm": res.get("wf_first"),
           "eps_p_max": res.get("eps_p_max"), "W_p": res.get("W_p"), "E_tot_max/E_tot_0": res.get("etot_max"),
           "E_tot_end/E_tot_0": res.get("etot_end"), "n_tunnel": res.get("tunnel"), "t_end_us": res.get("t_end"),
           "t_max_us": res.get("t_max"), "wall_s": res.get("wall_s"),
           "strict_perforated": strict_cell(res)}
    if multi():
        row["n_plies"] = len(res["plies"]) if res.get("plies") else None
        for k, p in enumerate(res.get("plies") or [], 1):
            row.update({"n_broken_v_p%d" % k: p["nbv"], "n_broken_h_p%d" % k: p["nbh"],
                        "t_first_break_us_p%d" % k: p["t_first"], "wf_first_break_mm_p%d" % k: p["wf_first"]})
    if csv_:
        pcsv = os.path.join(work, "shots.csv")
        new = not os.path.exists(pcsv) or os.path.getsize(pcsv) == 0
        with open(pcsv, "a", newline="", encoding="utf-8") as f:
            fields = csv_cols()
            if not new:  # дозапись в файл прежнего формата: берём его заголовок
                with open(pcsv, encoding="utf-8") as fh:
                    fields = next(csv.reader(fh), fields)
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerow(row)
    if jsonl:
        open(os.path.join(work, "eval.jsonl"), "a", encoding="utf-8").write(
            json.dumps({**row, "eff": eff, "n_bonds_plastic": res.get("n_bonds_plastic"), "plies": res.get("plies")},
                   ensure_ascii=False) + "\n")


def evaluate(work, mu, v0, run, res):
    """Блок оценки для НОВОГО выстрела: eval.md + строка в shots.csv + eval.jsonl (сразу после выстрела)."""
    write_outputs(work, mu, v0, run, res)


def ensure_eval(work, mu, v0, run, res):
    """Для уже готового выстрела: дописать только недостающие части (eval.md / shots.csv / eval.jsonl)."""
    k = _key(mu, v0, run)
    need_md = not md_has(work, mu, v0, run)
    need_csv = k not in done_csv(work)
    need_jsonl = k not in done_jsonl(work)
    if not (need_md or need_csv or need_jsonl):
        return False
    missing = [n for n, f in (("eval.md", need_md), ("shots.csv", need_csv), ("eval.jsonl", need_jsonl)) if f]
    log("ДОЗАПИСЬ mu=%.3f V0=%.1f %s: нет в %s" % (mu, v0, run, ", ".join(missing)))
    write_outputs(work, mu, v0, run, res, need_md, need_csv, need_jsonl)
    return True


def check_energy_build(binary):
    """Предпроверка --energy-log: сборка без SPH_ENERGY_LOG выходит с кодом 2 и сообщением уже при разборе аргументов
    (до расчёта). Вызов только с --energy-log=1: без --config сборка с логом тоже выходит сразу, расчёта не будет."""
    try:
        p = subprocess.run([binary, "--energy-log=1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
    except FileNotFoundError:
        log("ОШИБКА: бинарь не найден: %s -- серия остановлена" % binary)
        sys.exit(1)
    # Сборка с логом принимает --energy-log и останавливается дальше, на отсутствии --config («Нужен --config»).
    # Сборка без лога (или без поддержки ключа) падает раньше, на самом --energy-log. Справка (usage) упоминает
    # SPH_ENERGY_LOG в любой сборке, поэтому по этой строке судить нельзя.
    if "Нужен --config".encode("utf-8") not in (p.stdout or b"") + (p.stderr or b""):
        log("ОШИБКА: %s собран без SPH_ENERGY_LOG, а --energy-log=N передан. Пересоберите с -DSPH_ENERGY_LOG=ON "
            "или уберите --energy-log. Серия не запущена." % binary)
        sys.exit(1)
    log("предпроверка --energy-log: бинарь с SPH_ENERGY_LOG")


def run_shot(args, mu, v0, tmax, suffix=""):
    out = shot_dir(args.work, mu, v0, suffix)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = [args.bin, "--config=" + getattr(args, "config", "1T"), "--level=base", "--v0=%g" % v0, "--mu=%g" % mu, "--out=" + out, "--frames=off",
           "--threads=%d" % args.threads, "--n-mon=10", "--n-diag=100", "--repeat=1", "--no-bracket", "--save-t=20e-6",
           "--t-max=%g" % tmax, "--fiber-inelastic=" + args.fib, "--eps-break=%g" % args.eps_break]
    if MU_FF is not None:
        cmd.append("--mu-ff=%g" % MU_FF)
    if MU_FS is not None:
        cmd.append("--mu-fs=%g" % MU_FS)
    if args.energy_log > 0:
        cmd.append("--energy-log=%d" % args.energy_log)
    with open(out + ".log", "w") as lg:
        try:
            RC[out] = subprocess.run(cmd, stdout=lg, stderr=subprocess.STDOUT).returncode
        except FileNotFoundError:
            log("ОШИБКА: бинарь не найден: %s -- серия остановлена" % args.bin)
            sys.exit(1)
    return read_result(out, v0, tmax)


def get(args, mu, v0):
    sdir = shot_dir(args.work, mu, v0)
    res = read_result(sdir, v0, args.t_max)
    if res is None:
        log("ЗАПУСК mu=%.3f V0=%.1f t_max=%g" % (mu, v0, args.t_max))
        res = run_shot(args, mu, v0, args.t_max)
        if res is None:
            log("ОШИБКА: нет результата mu=%.3f V0=%.1f (см. %s.log)" % (mu, v0, sdir))
            return None
        evaluate(args.work, mu, v0, "base", res)
    else:
        ensure_eval(args.work, mu, v0, "base", res)
    if effective(res) == "undetermined":
        ldir = shot_dir(args.work, mu, v0, "_long")
        rl = read_result(ldir, v0, args.t_max_long)
        if rl is None:
            log("исход не определён -> повтор mu=%.3f V0=%.1f с t_max=%g" % (mu, v0, args.t_max_long))
            rl = run_shot(args, mu, v0, args.t_max_long, "_long")
            if rl is not None:
                evaluate(args.work, mu, v0, "long", rl)
        else:
            ensure_eval(args.work, mu, v0, "long", rl)
        if rl is not None and effective(rl) != "undetermined":
            return rl
    return res


def recover(args, mus, v0s):
    """При старте: все готовые выстрелы (есть shots.csv) -- проверить и дозаписать недостающие блоки оценки."""
    n = 0
    for v0 in v0s:
        for mu in mus:
            for run, suffix, tmax in (("base", "", args.t_max), ("long", "_long", args.t_max_long)):
                sdir = shot_dir(args.work, mu, v0, suffix)
                if not os.path.exists(os.path.join(sdir, "shots.csv")):
                    continue
                try:
                    res = read_result(sdir, v0, tmax)
                    if res is not None and ensure_eval(args.work, mu, v0, run, res):
                        n += 1
                except Exception:
                    report_error("восстановление mu=%.3f V0=%.1f %s" % (mu, v0, run))
    log("проверка готовых выстрелов завершена, дописано блоков: %d" % n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.path.expanduser("~/plastic1T"))
    ap.add_argument("--bin", default=os.path.expanduser("~/SPH_plastic/build/guo2019_shot"))
    ap.add_argument("--fib", default=os.path.join(os.path.expanduser("~/plastic1T"), "twaron_plastic.fib"))
    ap.add_argument("--eps-break", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--energy-log", type=int, default=0, metavar="N",
                    help="energy_log.csv каждые N шагов (бинарь с -DSPH_ENERGY_LOG=ON); 0 -- выкл.")
    ap.add_argument("--config", choices=("1T", "2T", "4T"), default="1T",
                    help="мишень guo2019_shot; 2T/4T -- послойный вывод и послойный строгий критерий")
    ap.add_argument("--criterion", choices=("default", "strict"), default=None,
                    help="strict: не пробит = not_perforated И V_R = 0 И верт. <= 5 И гориз. = 0 (1T); для 2T/4T: V_R = 0 И в "
                         "каждом слое верт. <= 7 И гориз. = 0. По умолчанию: default для 1T, strict для 2T/4T")
    ap.add_argument("--mu-ff", type=float, default=None, metavar="X",
                    help="трение нить-нить (переопределяет --mu/значение из --mus только для mu_ff); имя каталога: mu<mu>_ff<X>_...")
    ap.add_argument("--mu-fs", type=float, default=None, metavar="X",
                    help="трение нить-Solid (переопределяет только mu_fs); имя каталога: mu<mu>_fs<X>_...")
    ap.add_argument("--t-max", type=float, default=1.5e-3)
    ap.add_argument("--t-max-long", type=float, default=3.0e-3)
    ap.add_argument("--mus", default=",".join(str(m) for m in MUS))
    ap.add_argument("--v0s", default=",".join(str(v) for v in V0S))
    args = ap.parse_args()
    global CRITERION, CONFIG, MU_FF, MU_FS
    CONFIG = args.config
    MU_FF, MU_FS = args.mu_ff, args.mu_fs
    CRITERION = args.criterion or ("strict" if multi() else "default")
    os.makedirs(args.work, exist_ok=True)
    mus = [float(x) for x in args.mus.split(",")]
    v0s = [float(x) for x in args.v0s.split(",")]
    log("серия %s base пластика: mu=%s V0=%s fib=%s eps_break=%g criterion=%s" % (
        CONFIG, mus, v0s, args.fib, args.eps_break, CRITERION))
    if args.energy_log > 0:
        check_energy_build(args.bin)
    recover(args, mus, v0s)
    for v0 in v0s:
        for mu in mus:
            try:
                get(args, mu, v0)
            except Exception:
                report_error("mu=%.3f V0=%.1f" % (mu, v0))
    log("СЕРИЯ ЗАВЕРШЕНА")


if __name__ == "__main__":
    main()
