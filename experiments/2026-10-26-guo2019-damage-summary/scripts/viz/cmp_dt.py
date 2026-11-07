"""П. 1: какой канал зависит от шага. A (шаг x1) против B (шаг x0.5), закон 3 + выталкивание tau = 0,5 мкс, 73 м/с.
Выход: проверка честности сравнения, каналы энергии (% E0) по времени и разность B - A,
импульс шара и реакция закрепления (интеграл Ff_z), график."""
import csv, glob, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dA, dB, out = sys.argv[1], sys.argv[2], sys.argv[3]


def rd(d):
    r = glob.glob(os.path.join(d, "*", "shot.txt"))[0]
    run = os.path.dirname(r)
    kv = dict(t.split("=", 1) for t in open(r, encoding="utf-8").read().split() if "=" in t)
    e = list(csv.DictReader(open(os.path.join(run, "energy_log.csv"))))
    m = [q for q in csv.DictReader(open(os.path.join(run, "monitor.csv"))) if q.get("w_p")]
    return kv, e, m


A, B = rd(dA), rd(dB)
print("== Честность сравнения")
for name, (kv, e, m) in (("A", A), ("B", B)):
    E0 = float(e[0]["E_mech"])
    print("  %s: dt_mean %s dt_min %s, шагов %s, t_end %.0f мкс, исход %s, V_R %s, разрывы %s/%s, первый разрыв %.0f мкс, "
          "max delta/d ff %s fs %s, n_tunnel %s, невязка в конце %.2f %% E0" % (
              name, kv["dt_mean"], kv["dt_min"], kv["steps"], 1e6 * float(kv["t_end"]), kv["outcome"], kv["V_R"],
              kv["n_broken_v"], kv["n_broken_h"], 1e6 * float(kv["t_first_break"]), kv["max_delta_ff/d"],
              kv["max_delta_fs/d"], kv["n_tunnel"], 100 * float(e[-1]["residual"]) / E0))
print("  отношение dt_mean B/A = %.3f" % (float(B[0]["dt_mean"]) / float(A[0]["dt_mean"])))

ch = ["KE_solid", "KE_fiber", "U_bond", "D_friction", "D_contact_fiber", "W_p_fiber", "D_break", "D_xsph", "residual"]
E0 = float(A[1][0]["E_mech"])
tend = min(float(A[1][-1]["time"]), float(B[1][-1]["time"]))


def at(e, t):
    return min(e, key=lambda q: abs(float(q["time"]) - t))


print("\n== Каналы энергии, %% E0 (E0 = %.2f Дж); в скобках B - A" % E0)
print("  t, мкс | " + " | ".join(ch))
for t in [x * 1e-6 for x in range(100, 1501, 100)] + [tend]:
    if t > tend + 1e-9:
        continue
    a, b = at(A[1], t), at(B[1], t)
    print("  %5.0f  | " % (t * 1e6) + " | ".join(
        "%.1f (%+.1f)" % (100 * float(b[k]) / E0, 100 * (float(b[k]) - float(a[k])) / E0) for k in ch))

# Импульс: шар (m v_z) и реакция закрепления ткани (интеграл Ff_z dt).
m_ball = 2 * float(A[1][0]["KE_solid"]) / 73.0 ** 2


def series(m):
    t = [float(q["t"]) for q in m]
    vz = [float(q["v_z"]) for q in m]
    fz = [float(q["Ff_z"]) for q in m]
    J = [0.0]
    for i in range(1, len(t)):
        J.append(J[-1] + 0.5 * (fz[i] + fz[i - 1]) * (t[i] - t[i - 1]))
    return t, vz, fz, J


sA, sB = series(A[2]), series(B[2])


def val(s, t, k):
    i = min(range(len(s[0])), key=lambda j: abs(s[0][j] - t))
    return s[k][i]


print("\n== Импульс (m шара = %.4f кг): потеря импульса шара и импульс, ушедший в закрепление" % m_ball)
print("  t, мкс | m(v_z + V0) A, B (Н с) | интеграл Ff_z A, B (Н с) | доля закрепления в разнице B - A")
for t in [x * 1e-6 for x in range(300, 1501, 100)]:
    if t > tend:
        continue
    pa = m_ball * (val(sA, t, 1) + 73.0)
    pb = m_ball * (val(sB, t, 1) + 73.0)
    ja, jb = val(sA, t, 3), val(sB, t, 3)
    share = (jb - ja) / (pb - pa) if abs(pb - pa) > 1e-6 else float("nan")
    print("  %5.0f  | %.4f  %.4f | %.4f  %.4f | %.2f" % (t * 1e6, pa, pb, ja, jb, share))

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
for (lab, s, e), c in ((("A: шаг x1", sA, A[1]), "#2a78d6"), (("B: шаг x0,5", sB, B[1]), "#d6452a")):
    ax[0].plot([x * 1e6 for x in s[0]], s[1], color=c, label=lab)
    ax[1].plot([x * 1e6 for x in s[0]], s[2], color=c, label=lab)
    te = [float(q["time"]) * 1e6 for q in e]
    ax[2].plot(te, [100 * float(q["D_friction"]) / E0 for q in e], color=c, label=lab + ", трение")
    ax[2].plot(te, [100 * float(q["D_contact_fiber"]) / E0 for q in e], color=c, ls="--", label=lab + ", контакт")
ax[0].set(title="v_z шара", xlabel="t, мкс", ylabel="м/с")
ax[1].set(title="Реакция закрепления ткани Ff_z", xlabel="t, мкс", ylabel="Н")
ax[2].set(title="Потери: трение и контакт нитей", xlabel="t, мкс", ylabel="% E0")
for a in ax:
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
fig.suptitle("1T, закон 3 + выталкивание τ = 0,5 мкс, E = 80 ГПа, μ = 0,1, 73 м/с: шаг x1 против x0,5")
fig.tight_layout()
fig.savefig(out, dpi=110)
print("\nграфик:", out)
