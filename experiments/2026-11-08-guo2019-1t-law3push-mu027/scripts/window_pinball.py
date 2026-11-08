#!/usr/bin/env python3
"""Подбор скорости непробития / пробития Guo 2019, 1T base, закон контакта нитей pinball.

Постановка (решения владельца, 2026-10-17):
  - fs = ff = 1 (pinball на обоих контактах нитей);
  - нить упругопластичная без вязкости и скоростной зависимости, sigma_u = 2.5 ГПа, E = 97 ГПа;
  - mu по порядку 0.1, 0.3, 0.5, 0.7, 0.9; скорости только из [70; 82];
  - «не пробит» = outcome not_perforated И V_R = 0 И порвано верт. <= 8 И гориз. = 0; иначе «пробит».

Алгоритм на одном mu (два выстрела параллельно по --threads потоков):
  1. пара 70 / 76:
       70 П              -> непробития на этом mu нет -> следующий mu;
       70 НП, 76 П       -> окно [70; 76] -> уточнение 73 (один выстрел, --threads-single) -> конец серии;
       70 НП, 76 НП      -> шаг 2;
       70 П,  76 НП      -> немонотонно (записывается), непробитие на 76 -> шаг 2.
  2. пара 79 / 82:
       79 П              -> окно [76; 79] -> конец;
       79 НП, 82 П       -> окно [79; 82] -> конец;
       79 НП, 82 НП      -> пробитие в [70; 82] не найдено -> конец (нужен анализ причин).
  Исход «не определён» -> повтор с --t-max-long; если и он не определён -> mu «вывод невозможен».
  Если ни при каком mu непробития нет -> конец (нужен анализ причин).

Перезапуск той же командой безопасен: выстрел с готовым shots.csv не пересчитывается,
а идущий выстрел с тем же --out (запущенный вручную) драйвер ждёт, не стирая каталог.
Упавший выстрел (rc != 0 или нет shots.csv) повторяется один раз; повторное падение останавливает серию.
"""
import argparse, concurrent.futures, csv, glob, json, os, re, shutil, subprocess, sys, threading, time

MUS = [0.1, 0.3, 0.5, 0.7, 0.9]
MAX_BROKEN_V = 8  # только 1T; для 2T вердикт по шару, разрывы -- отдельная метрика (владелец, 2026-10-10)
BALL_ONLY = False  # True: и для 1T вердикт по шару, без допуска на разрывы (владелец, 2026-10-25; ключ --verdict-ball)
LOCK = threading.Lock()


def log(work, msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    with LOCK:
        print(line, flush=True)
        with open(os.path.join(work, "queue.log"), "a", encoding="utf-8") as f:
            f.write(line + "\n")


def kv(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for tok in f.read().split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                out[k] = v
    return out


def num(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def shot_dir(work, mu, v0, suffix=""):
    return os.path.join(work, "mu%g_v%g%s" % (mu, v0, suffix))


def read_result(sdir):
    files = glob.glob(os.path.join(sdir, "*", "shot.txt"))
    if not files or not os.path.exists(os.path.join(sdir, "shots.csv")):
        return None
    r = kv(files[0])
    r["_dir"] = os.path.dirname(files[0])
    return r


def verdict(r):
    """True -- «пробит», False -- «не пробит», None -- не определён. Вторая величина -- причина."""
    out = r.get("outcome", "")
    vr = num(r.get("V_R"), 0.0)
    nv = num(r.get("n_broken_v"), 0.0)
    nh = num(r.get("n_broken_h"), 0.0)
    if out == "perforated":
        return True, "пробит шаром (V_R = %.2f)" % vr
    if out == "not_perforated":
        if vr > 0.0:
            return True, "V_R = %.2f > 0" % vr
        if r.get("config", "1T") == "1T" and not BALL_ONLY and (nv > MAX_BROKEN_V or nh > 0):
            return True, "шар остановлен, но порвано верт. %d / гориз. %d (допуск <= %d / 0)" % (nv, nh, MAX_BROKEN_V)
        return False, "не пробит: V_R = 0, порвано %d / %d" % (nv, nh)
    if BALL_ONLY and out == "undetermined" and r.get("_dir"):
        # к t_max шар практически стоит за плоскостью ткани: |v_z| < 5 % V0 и z_c < 0 -> «не пробит» (оценка),
        # без повтора до t_max_long (правило прежних драйверов; владелец, 2026-10-25)
        try:
            with open(os.path.join(r["_dir"], "monitor.csv")) as f:
                last = [q for q in csv.DictReader(f) if q.get("w_p")][-1]
            vz, zc = float(last["v_z"]), float(last["z_c"])
            if abs(vz) < 0.05 * num(r.get("V0"), 0.0) and zc < 0.0:
                return False, "не пробит (оценка): к t_max |v_z| = %.2f м/с < 5 %% V0, z_c = %.1f мм; порвано %d / %d" % (
                    abs(vz), zc * 1e3, nv, nh)
        except (OSError, IndexError, KeyError, ValueError):
            pass
    return None, "исход не определён (%s)" % out


def has_checkpoint(sdir):
    """Есть ли в каталоге выстрела sdir чекпоинт guo2019_shot (checkpoint.ckpt или .prev)."""
    return bool(glob.glob(os.path.join(sdir, "*", "checkpoint.ckpt")) + glob.glob(os.path.join(sdir, "*", "checkpoint.ckpt.prev")))


def drop_checkpoint(sdir):
    """После определённого результата чекпоинты не нужны (сотни МБ на выстрел)."""
    for f in glob.glob(os.path.join(sdir, "*", "checkpoint.ckpt*")):
        try:
            os.remove(f)
        except OSError:
            pass


# ---------------------------------------------------------------------------------------------------------------
# Рестарт без состояния монитора (fiber_new ec2439f: GuoShotMonitor в чекпоинт не входит; contact3_fix не влит).
# guo2019_shot при --restart ПЕРЕЗАПИСЫВАЕТ monitor.csv с шага рестарта, а energy_log.csv ДОПИСЫВАЕТ (хвост после
# последнего чекпоинта повторится). Драйвер перед каждым --restart откладывает monitor.csv в monitor.partK.csv, после
# выстрела склеивает части по времени и пересчитывает в shot.txt величины, которые монитор копит с начала выстрела
# (первый разрыв, прогиб к нему, максимумы прогиба). n_tunnel, max δ/d, E_tot_max, реакции зажимов -- с последнего
# рестарта; в shot.txt ключ restarts=K, исходный файл -- shot.raw.txt.
# ---------------------------------------------------------------------------------------------------------------

def _read_csv(path):
    """Строки CSV без обрезанного хвоста: после kill -9 последняя строка может быть недописана (файл не кончается
    переводом строки) -- она отбрасывается, как и строки с другим числом полей."""
    with open(path, encoding="utf-8", newline="") as f:
        text = f.read()
    if text and not text.endswith(chr(10)):
        text = text[:text.rfind(chr(10)) + 1]
    rows = list(csv.reader(text.splitlines()))
    if not rows:
        return None, []
    return rows[0], [r for r in rows[1:] if len(r) == len(rows[0])]


def splice_rows(header, chunks, tcol):
    """Склейка по времени: строка с t <= последнего принятого t выбрасывает из принятых все строки с t >= её t
    (хвост после последнего чекпоинта, посчитанный заново после рестарта)."""
    ti = header.index(tcol)
    out = []
    gaps = []
    for rows in chunks:
        # Монитор пишет monitor.csv через буфер stdio и при чекпоинте его не сбрасывает: после kill -9 строки между
        # последней записанной и шагом чекпоинта теряются. Дыра -- скачок t на стыке частей больше 3 медианных шагов.
        ts = []
        for row in rows:
            try:
                ts.append(float(row[ti]))
            except (IndexError, ValueError):
                pass
        if out and ts:
            dts = sorted(b - a for a, b in zip(ts, ts[1:]) if b > a)
            step = dts[len(dts) // 2] if dts else 0.0
            t_prev = max((float(r[ti]) for r in out if float(r[ti]) < ts[0]), default=None)
            if t_prev is not None and step > 0.0 and ts[0] - t_prev > 3.0 * step:
                gaps.append((t_prev, ts[0]))
        for row in rows:
            try:
                t = float(row[ti])
            except (IndexError, ValueError):
                continue
            while out and float(out[-1][ti]) >= t:
                out.pop()
            out.append(row)
    splice_rows.gaps = gaps
    return out


def _write_csv(path, header, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    os.replace(tmp, path)


def monitor_parts(sub):
    return sorted(glob.glob(os.path.join(sub, "monitor.part*.csv")),
                  key=lambda q: int(re.search(r"part(\d+)\.csv$", q).group(1)))


def stash_monitor(sdir):
    """Перед --restart: monitor.csv -> monitor.partK.csv (иначе выстрел его перезапишет)."""
    for sub in glob.glob(os.path.join(sdir, "*", "")):
        m = os.path.join(sub, "monitor.csv")
        if os.path.exists(m):
            os.replace(m, os.path.join(sub, "monitor.part%d.csv" % (len(monitor_parts(sub)) + 1)))


def splice_shot(sdir):
    """После выстрела с рестартами: склеить monitor.csv и energy_log.csv, поправить shot.txt. Возвращает K."""
    k = 0
    for sub in glob.glob(os.path.join(sdir, "*", "")):
        parts = monitor_parts(sub)
        m = os.path.join(sub, "monitor.csv")
        if not parts or not os.path.exists(m):
            continue
        k = max(k, len(parts))
        header, last = _read_csv(m)
        rows = splice_rows(header, [_read_csv(q)[1] for q in parts] + [last], "t")
        gaps = splice_rows.gaps
        _write_csv(m, header, rows)
        for q in parts:
            os.remove(q)
        el = os.path.join(sub, "energy_log.csv")
        if os.path.exists(el):
            eh, er = _read_csv(el)
            _write_csv(el, eh, splice_rows(eh, [er], "time"))
        patch_shot_txt(sub, header, rows, len(parts), gaps)
    return k


def patch_shot_txt(sub, header, rows, k, gaps=()):
    st = os.path.join(sub, "shot.txt")
    if not os.path.exists(st) or not rows:
        return
    col = {n: i for i, n in enumerate(header)}
    fix = {"restarts": str(k)}
    nb0 = float(rows[0][col["n_bonds"]])
    fb = next((r for r in rows if float(r[col["n_bonds"]]) < nb0), None)
    fix["t_first_break"] = "%.9g" % float(fb[col["t"]]) if fb else "nan"
    fix["w_f_at_first_break"] = "%.9g" % float(fb[col["w_f"]]) if fb else "nan"
    # дыры в мониторе от буфера stdio; если первый разрыв -- первая строка после дыры, он мог случиться раньше
    fix["monitor_gaps"] = str(len(gaps))
    fix["monitor_gap_max"] = "%.9g" % max((b - a for a, b in gaps), default=0.0)
    if fb:
        tfb = float(fb[col["t"]])
        lo = next((a for a, b in gaps if a < tfb <= b), None)
        if lo is not None:
            fix["t_first_break_lo"] = "%.9g" % lo
    for key, c in (("w_f_max", "w_f"), ("w_p_max", "w_p")):
        vals = []
        for r in rows:
            try:
                vals.append(float(r[col[c]]))
            except (IndexError, ValueError):
                pass
        vals = [v for v in vals if v == v]
        if vals:
            fix[key] = "%.9g" % max(vals)
    raw = os.path.join(sub, "shot.raw.txt")
    if not os.path.exists(raw):
        shutil.copy(st, raw)
    with open(raw, encoding="utf-8") as f:
        text = f.read()
    # Точный шаг первого разрыва монитор знает сам (строки monitor.csv -- раз в n_mon шагов). Его значение верно,
    # если разрыв случился после последнего рестарта: тогда оно лежит между строкой перед первым разрывом и строкой
    # первого разрыва. Иначе (разрыв до рестарта) монитор ставит момент рестарта -- берём значение по строкам.
    if fb:
        i = rows.index(fb)
        m_t = re.search(r"(?<!\S)t_first_break=(\S+)", text)
        m_w = re.search(r"(?<!\S)w_f_at_first_break=(\S+)", text)
        try:
            t_raw = float(m_t.group(1))
            if i > 0 and float(rows[i - 1][col["t"]]) < t_raw <= float(fb[col["t"]]):
                fix["t_first_break"] = m_t.group(1)
                if m_w:
                    fix["w_f_at_first_break"] = m_w.group(1)
        except (AttributeError, ValueError):
            pass
    for key, val in fix.items():
        text, n = re.subn(r"(?<!\S)%s=\S*" % re.escape(key), "%s=%s" % (key, val), text, count=1)
        if n == 0:
            text = text.rstrip("\n") + " %s=%s\n" % (key, val)
    with open(st + ".tmp", "w", encoding="utf-8") as f:
        f.write(text)
    os.replace(st + ".tmp", st)


FRAME_RE = re.compile(r"(?:sph|fiber)_output_\d+\.txt$")


def _link_frame_or_copy(src, dst):
    """Кадры анимации (sph_output_*, fiber_output_* -- ~5 МБ, раз в save_t) после записи не меняются: в копию
    продолжения кладётся жёсткая ссылка, а не вторая копия. Остальное (energy_log.csv дописывается, monitor.csv
    переименовывается) -- копией, чтобы не испортить базовый прогон."""
    if FRAME_RE.search(os.path.basename(src)):
        try:
            os.link(src, dst)
            return dst
        except OSError:
            pass
    return shutil.copy2(src, dst)


def seed_long_from_base(base, sdir):
    """«Не определён» к t_max: продолжить до t_max_long из финального чекпоинта базового прогона (t_total при
    рестарте менять можно). Копия каталога без результатов базового прогона. False -- чекпоинта нет."""
    if not has_checkpoint(base):
        return False
    if os.path.isdir(sdir):
        shutil.rmtree(sdir)
    shutil.copytree(base, sdir, copy_function=_link_frame_or_copy)
    for pat in ("shots.csv", os.path.join("*", "shot.txt"), os.path.join("*", "shot.raw.txt")):
        for f in glob.glob(os.path.join(sdir, pat)):
            os.remove(f)
    if os.path.exists(base + ".log"):
        shutil.copy(base + ".log", sdir + ".log")
    return True


def shot_running(sdir):
    """Идёт ли уже guo2019_shot с --out=sdir (запущен вручную или прежним экземпляром драйвера)."""
    tag = ("--out=" + sdir).encode()
    for p in glob.glob("/proc/[0-9]*/cmdline"):
        try:
            with open(p, "rb") as f:
                if tag in f.read().split(b"\0"):
                    return True
        except OSError:
            pass
    return False


# --- Допуск выстрела к запуску по свободным потокам (правила claude-compute-8t, 2026-11-08):
#     новый выстрел на N потоков ставится, только если uptime > 600 с, proj_status.sh != 2 и
#     threads_free.sh показывает FREE >= N. Через --slot-verify-s после старта проверяется, что выстрел
#     виден в threads_free.sh со своими N потоками (или BUSY вырос на N); иначе выстрел снимается, серия стоп.
#     GATE: проверка и старт -- по одному выстрелу, второй выстрел видит потоки первого.
GATE = threading.Lock()
GATE_FAILED = []  # непусто -- проверка потоков уже не прошла, новых выстрелов не ставить


def _shot_threads():
    """Сумма --threads= у идущих guo2019_shot (запасной счёт, если нет threads_free.sh)."""
    total = 0
    for p in glob.glob("/proc/[0-9]*/cmdline"):
        try:
            with open(p, "rb") as f:
                argv = f.read().split(b"\0")
        except OSError:
            continue
        if argv and argv[0].endswith(b"guo2019_shot"):
            for a in argv:
                if a.startswith(b"--threads="):
                    try:
                        total += int(a[len(b"--threads="):])
                    except ValueError:
                        pass
    return total


def threads_free(args):
    """(BUSY, FREE, вывод) по threads_free.sh; без скрипта -- по --threads= идущих guo2019_shot из 8.
    None -- скрипт не дал итоговой строки."""
    sc = os.path.expanduser(getattr(args, "slot_script", "~/proj_runs/threads_free.sh"))
    if not os.path.exists(sc):
        b = _shot_threads()
        return b, 8 - b, ""
    try:
        out = subprocess.run([sc], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                             timeout=120).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(r"BUSY=(\d+)\s+FREE=(-?\d+)", out)
    return (int(m.group(1)), int(m.group(2)), out) if m else None


def proj_rc(args):
    sc = os.path.expanduser(getattr(args, "proj_status", "~/proj_runs/proj_status.sh"))
    if not os.path.exists(sc):
        return 0  # серия proj_* убрана с сервера
    try:
        return subprocess.run([sc], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120).returncode
    except (OSError, subprocess.SubprocessError):
        return None


def uptime_s():
    try:
        with open("/proc/uptime") as f:
            return float(f.read().split()[0])
    except (OSError, ValueError, IndexError):
        return 1e9


def start_gated(args, cmd, env, lg, threads, label):
    """Дождаться свободных потоков, запустить выстрел, проверить, что он занял потоки. Возвращает Popen."""
    min_up = getattr(args, "min_uptime", 600)
    poll = getattr(args, "slot_poll", 300)
    verify = getattr(args, "slot_verify_s", 120)
    with GATE:
        if GATE_FAILED:
            raise RuntimeError("%s не запущен: %s" % (label, GATE_FAILED[0]))
        last, t_last = None, 0.0
        while True:
            up, rc, tf = uptime_s(), proj_rc(args), threads_free(args)
            if up >= min_up and rc in (0, 1) and tf is not None and tf[1] >= threads:
                break
            state = "uptime=%.0f с, proj rc=%s, %s" % (
                up, rc, "BUSY=%d FREE=%d" % tf[:2] if tf else "threads_free.sh: нет итога")
            if rc == 2:
                state += " -- proj: ОСТАНОВЛЕНО С ОШИБКОЙ, новые выстрелы не ставлю, жду решения владельца"
            if state != last or time.time() - t_last > 3600:
                log(args.work, "жду %d свободных потоков для %s: %s" % (threads, label, state))
                last, t_last = state, time.time()
            time.sleep(poll)
        busy0 = tf[0]
        p = subprocess.Popen(cmd, stdout=lg, stderr=subprocess.STDOUT, env=env)
        log(args.work, "запущен %s: PID %d, %d потоков (было BUSY=%d FREE=%d)" % (label, p.pid, threads, busy0, tf[1]))
        t_end = time.time() + verify
        while time.time() < t_end and p.poll() is None:
            time.sleep(5)
        if p.poll() is None:
            tf2 = threads_free(args)
            seen = tf2 is not None and re.search(r"^%d \S+ threads=%d\b" % (p.pid, threads), tf2[2], re.M)
            if not (seen or (tf2 is not None and tf2[0] >= busy0 + threads)):
                p.kill()
                p.wait()
                msg = "ОШИБКА: %s (PID %d) не занял %d потоков через %d с (BUSY было %d, стало %s) -- выстрел снят" % (
                    label, p.pid, threads, verify, busy0, tf2[0] if tf2 else "?")
                log(args.work, msg)
                GATE_FAILED.append(msg)
                raise RuntimeError(msg)
            log(args.work, "проверка %s: занял %d потоков (BUSY=%s FREE=%s)" % (
                label, threads, tf2[0], tf2[1]))
        return p


def run_shot(args, mu, v0, threads, tmax, suffix="", seed=None):
    """seed -- каталог базового прогона: без своего чекпоинта продолжить из его финального."""
    sdir = shot_dir(args.work, mu, v0, suffix)
    r = read_result(sdir)
    if r is not None:
        return r
    # Идущий выстрел в том же каталоге не стирать: дождаться и взять его результат.
    if shot_running(sdir):
        log(args.work, "жду идущий выстрел mu=%g V0=%g%s" % (mu, v0, suffix))
        while shot_running(sdir):
            time.sleep(60)
        r = read_result(sdir)
        if r is not None:
            log(args.work, "готово (внешний) mu=%g V0=%g%s: %s" % (mu, v0, suffix, verdict(r)[1]))
            return r
    ck_every = getattr(args, "checkpoint_every", 0)
    for attempt in (1, 2):
        # Чекпоинт от оборванного прогона (перезагрузка сервера, kill) -> продолжить, а не считать заново.
        # Вторая попытка (первая упала) всегда с нуля: битый чекпоинт не должен зациклить выстрел.
        if attempt == 1 and seed and ck_every > 0 and not has_checkpoint(sdir) and seed_long_from_base(seed, sdir):
            log(args.work, "продолжение mu=%g V0=%g%s из финального чекпоинта %s" % (
                mu, v0, suffix, os.path.basename(seed)))
        restart = attempt == 1 and ck_every > 0 and has_checkpoint(sdir)
        if restart:
            stash_monitor(sdir)
        if not restart:
            if os.path.isdir(sdir):
                shutil.rmtree(sdir)
            os.makedirs(sdir)
        cmd = [args.bin, "--config=" + getattr(args, "config", "1T"), "--level=" + args.level, "--v0=%g" % v0, "--mu=%g" % mu, "--out=" + sdir,
               "--frames=off", "--save-t=%g" % args.save_t, "--threads=%d" % threads, "--n-mon=10", "--n-diag=100",
               "--repeat=1", "--no-bracket", "--t-max=%g" % tmax, "--fiber-inelastic=" + args.fib,
               "--eps-break=%g" % args.eps_break, "--fiber-E=%g" % args.fiber_E,
               "--fs-contact=%d" % args.contact, "--ff-contact=%d" % args.contact]
        if args.energy_log > 0:
            cmd.append("--energy-log=%d" % args.energy_log)
        cmd += getattr(args, "extra", "").split()
        if ck_every > 0:
            cmd.append("--checkpoint-every=%d" % ck_every)
        if restart:
            cmd.append("--restart")
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        log(args.work, "%s mu=%g V0=%g%s (%d потоков, t_max=%g, попытка %d)" % (
            "РЕСТАРТ из чекпоинта" if restart else "старт", mu, v0, suffix, threads, tmax, attempt))
        t0 = time.time()
        with open(sdir + ".log", "a" if restart else "w", encoding="utf-8") as lg:
            p = start_gated(args, cmd, env, lg, threads, "mu=%g V0=%g%s" % (mu, v0, suffix))
            rc = p.wait()
        r = read_result(sdir)
        if rc == 0 and r is not None:
            try:
                k = splice_shot(sdir)
            except Exception as ex:  # постобработка не должна останавливать серию
                k = 0
                log(args.work, "ПРЕДУПРЕЖДЕНИЕ: склейка монитора mu=%g V0=%g%s не удалась: %r" % (mu, v0, suffix, ex))
            r = read_result(sdir)
            p, why = verdict(r)
            if k:
                why += " [рестартов: %d; n_tunnel, max δ/d, E_tot_max -- с последнего рестарта]" % k
            log(args.work, "готово mu=%g V0=%g%s: %s, %.0f с" % (mu, v0, suffix, why, time.time() - t0))
            if p is not None:
                # «не определён» -- финальный чекпоинт нужен для продолжения до t_max_long
                drop_checkpoint(sdir)
                if seed:
                    drop_checkpoint(seed)
            return r
        log(args.work, "ПАДЕНИЕ mu=%g V0=%g%s: rc=%s, лог %s.log" % (mu, v0, suffix, rc, sdir))
        if attempt == 1 and os.path.exists(sdir + ".log"):
            shutil.copy(sdir + ".log", sdir + ".fail1.log")
    raise RuntimeError("выстрел mu=%g V0=%g%s упал дважды, лог %s.log" % (mu, v0, suffix, sdir))


def shoot(args, mu, v0, threads, rec):
    """Выстрел с повтором при «не определён». Возвращает вердикт (True/False/None)."""
    r = run_shot(args, mu, v0, threads, args.t_max)
    p, why = verdict(r)
    entry = {"mu": mu, "V0": v0, "run": "base", "verdict": p, "why": why, "res": r}
    if p is None:
        r2 = run_shot(args, mu, v0, threads, args.t_max_long, "_long", seed=shot_dir(args.work, mu, v0))
        p, why = verdict(r2)
        entry = {"mu": mu, "V0": v0, "run": "long", "verdict": p, "why": why, "res": r2}
    with LOCK:
        rec["shots"].append(entry)
    return p


def shoot_pair(args, mu, v1, v2, rec):
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        f1 = ex.submit(shoot, args, mu, v1, args.threads, rec)
        f2 = ex.submit(shoot, args, mu, v2, args.threads, rec)
        return f1.result(), f2.result()


def word(p):
    return {True: "П", False: "НП", None: "не опр."}[p]


def process_mu(args, mu, rec, save):
    """Возвращает 'window' | 'no_np' | 'no_p' | 'undetermined'."""
    p70, p76 = shoot_pair(args, mu, 70.0, 76.0, rec)
    save()
    if p70 is None or p76 is None:
        if p70 is True:
            rec["result"] = "на 70 м/с пробит; 76 не определён -> непробития нет"
            return "no_np"
        rec["result"] = "исход не определён (70: %s, 76: %s) -- вывод невозможен" % (word(p70), word(p76))
        return "undetermined"
    if p70 is True and p76 is True:
        rec["result"] = "пробит при 70 и 76 -> непробития в [70; 82] нет"
        return "no_np"
    if p70 is False and p76 is True:
        p73 = shoot(args, mu, 73.0, args.threads_single, rec)
        save()
        if p73 is False:
            rec["window"] = [73.0, 76.0]
        elif p73 is True:
            rec["window"] = [70.0, 73.0]
        else:
            rec["window"] = [70.0, 76.0]
        rec["result"] = "окно [%g; %g]" % tuple(rec["window"])
        return "window"
    if p70 is True and p76 is False:
        rec["nonmonotonic"] = "70 П, 76 НП"
    # непробитие на 76 -> ищем пробитие выше
    p79, p82 = shoot_pair(args, mu, 79.0, 82.0, rec)
    save()
    if p79 is True:
        rec["window"] = [76.0, 79.0]
    elif p79 is False and p82 is True:
        rec["window"] = [79.0, 82.0]
    elif p79 is False and p82 is False:
        rec["result"] = "не пробит при 76, 79, 82 -> пробития в [70; 82] нет"
        return "no_p"
    else:
        rec["result"] = "исход не определён (79: %s, 82: %s) -- вывод невозможен" % (word(p79), word(p82))
        return "undetermined"
    rec["result"] = "окно [%g; %g]" % tuple(rec["window"])
    return "window"


def write_window(args, records, final):
    law = {0: "Герц", 1: "pinball", 2: "Swegle", 3: "импульсный", 4: "Герц + демпфер"}.get(args.contact, "?")
    cfg = getattr(args, "config", "1T")
    lines = ["# Подбор окна %s base, закон контакта нитей %s (fs = ff = %d)" % (cfg, law, args.contact), "",
             "Нить: %s, E = %g Па, eps_break = %g; уровень %s; t_max = %g (повтор %g)." % (
                 os.path.basename(args.fib), args.fiber_E, args.eps_break, args.level, args.t_max, args.t_max_long),
             ("Критерий «не пробит»: not_perforated, V_R = 0, верт. <= %d, гориз. = 0." % MAX_BROKEN_V
              if cfg == "1T" and not BALL_ONLY else
              "Критерий «не пробит»: not_perforated, V_R = 0 (по шару; разрывы нитей -- отдельная метрика)."), ""]
    for rec in records:
        lines.append("## mu = %g" % rec["mu"])
        lines.append("")
        lines.append("| V0 | прогон | вердикт | причина | V_R | порвано v/h | первый разрыв, мкс | прогиб, мм | "
                     "max δff/d | max δfs/d | n_tunnel | E_end/E0 | W_p | t_end, мкс | wall, с |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for e in sorted(rec["shots"], key=lambda x: (x["V0"], x["run"])):
            r = e["res"]
            lines.append("| %g | %s | %s | %s | %.2f | %s / %s | %.0f | %.1f | %.3f | %.3f | %s | %.3f | %.3f | %.0f | %.0f |" % (
                e["V0"], e["run"], word(e["verdict"]), e["why"], num(r.get("V_R")), r.get("n_broken_v"),
                r.get("n_broken_h"), num(r.get("t_first_break")) * 1e6, num(r.get("w_f_at_first_break")) * 1e3,
                num(r.get("max_delta_ff/d")), num(r.get("max_delta_fs/d")), r.get("n_tunnel"),
                num(r.get("E_tot_end/E_tot_0")), num(r.get("W_p")), num(r.get("t_end"), 0.0) * 1e6,
                num(r.get("wall_s"))))
        lines.append("")
        if rec.get("nonmonotonic"):
            lines.append("- немонотонно: " + rec["nonmonotonic"])
        lines.append("- итог: " + rec.get("result", "в работе"))
        lines.append("")
    lines.append("## Итог серии")
    lines.append("")
    lines.append(final or "в работе")
    with open(os.path.join(args.work, "window.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(args.work, "window.json"), "w", encoding="utf-8") as f:
        json.dump({"records": [{k: v for k, v in rec.items()} for rec in records], "final": final}, f,
                  ensure_ascii=False, indent=1, default=str)


def wait_for(args):
    if not args.wait_file:
        return
    log(args.work, "жду %s и окончания других guo2019_shot" % args.wait_file)
    while True:
        busy = subprocess.run(["pgrep", "-f", "guo2019_shot"], stdout=subprocess.PIPE).stdout.strip()
        if os.path.exists(args.wait_file) and not busy:
            return
        time.sleep(60)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.path.expanduser("~/vbl1T_pinball"))
    ap.add_argument("--bin", default=os.path.expanduser("~/SPH_vbl/build-elog/guo2019_shot"))
    ap.add_argument("--fib", default=os.path.expanduser("~/vbl1T_pinball/twaron_plastic_s25.fib"))
    ap.add_argument("--fiber-E", type=float, default=97e9)
    ap.add_argument("--eps-break", type=float, default=0.05)
    ap.add_argument("--contact", type=int, default=1)
    ap.add_argument("--level", default="base")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--threads-single", type=int, default=8)
    ap.add_argument("--energy-log", type=int, default=100)
    ap.add_argument("--save-t", type=float, default=20e-6)
    ap.add_argument("--t-max", type=float, default=1.5e-3)
    ap.add_argument("--t-max-long", type=float, default=3.0e-3)
    ap.add_argument("--mus", default=",".join("%g" % m for m in MUS))
    ap.add_argument("--wait-file", default="")
    args = ap.parse_args()
    os.makedirs(args.work, exist_ok=True)
    for path in (args.bin, args.fib):
        if not os.path.exists(path):
            print("нет файла: " + path, file=sys.stderr)
            return 1
    wait_for(args)
    mus = [float(x) for x in args.mus.split(",")]
    log(args.work, "серия: mu=%s E=%g fib=%s contact=%d" % (mus, args.fiber_E, args.fib, args.contact))
    records, final = [], ""
    save = lambda: write_window(args, records, final)
    try:
        for mu in mus:
            rec = {"mu": mu, "shots": []}
            records.append(rec)
            res = process_mu(args, mu, rec, save)
            log(args.work, "mu=%g: %s" % (mu, rec.get("result")))
            if res == "window":
                final = "Окно найдено при mu = %g: [%g; %g] м/с." % (mu, rec["window"][0], rec["window"][1])
                break
            if res == "no_p":
                final = ("При mu = %g непробитие есть, но пробития в [70; 82] нет. Нужен анализ причин "
                         "(v_z(t), energy_log.csv, разрывы, max δ/d, кадры)." % mu)
                break
        else:
            final = "Окно не найдено ни при одном mu из %s. Нужен анализ причин." % mus
    except Exception as ex:
        final = "Серия ОСТАНОВЛЕНА: %s" % ex
        log(args.work, final)
        save()
        return 1
    log(args.work, final)
    save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
